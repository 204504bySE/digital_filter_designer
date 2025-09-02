#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Copyright (c) 2025 ArqAlice 

Released under the MIT license
https://opensource.org/licenses/mit-license.php

PySide6 を用いたデジタルフィルタ設計ツール。
- FIR (カイザー窓) と IIR (Chebyshev II, SOS表現) を切替
- ローパス / ハイパス / バンドパス を選択
- サンプリング周波数 [Hz]、通過域端周波数、阻止域端(カットオフ)周波数、減衰量[dB] を入力
- 出力：ゲイン線図[dB]、位相線図[deg]、インパルス応答（表示時間指定）、フィルタ係数
- 係数は C 言語の配列として出力（FIR taps と SOS 配列、CMSIS-DSP 互換形式も併記）

依存関係（例）:
    pip install numpy scipy matplotlib PySide6

実行:
    python filter_designer.py
"""
from __future__ import annotations
import sys
import math
import numpy as np

from PySide6.QtCore import Qt
from PySide6.QtGui import QClipboard, QCursor
from PySide6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QFormLayout, QLabel, QDoubleSpinBox, QSpinBox, QComboBox,
    QPushButton, QTextEdit, QGroupBox, QGridLayout, QMessageBox, QCheckBox, QToolTip
)

from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

from scipy.signal import (
    cheb2ord, cheby2, sosfreqz, lfilter, freqz, firwin
)

# -------------------------------
# ユーティリティ
# -------------------------------

def kaiser_beta_for_atten(atten_db: float) -> float:
    """カイザー窓のβを減衰量[dB]から求める。
    典型式：
      A <= 21: beta = 0
      21 < A <= 50: beta = 0.5842*(A-21)^0.4 + 0.07886*(A-21)
      A > 50: beta = 0.1102*(A-8.7)
    """
    A = float(atten_db)
    if A <= 21:
        return 0.0
    elif A <= 50:
        return 0.5842 * (A - 21) ** 0.4 + 0.07886 * (A - 21)
    else:
        return 0.1102 * (A - 8.7)


def kaiser_order(atten_db: float, trans_width_hz: float, fs: float) -> int:
    """カイザー窓 FIR の近似次数（タップ数 N-1）を求める。
    公式: N ≈ (A - 8) / (2.285 * Δω) + 1, ここで Δω = 2π * (Δf / fs)
    返すのは tap 数 N （次数+1）。
    """
    if trans_width_hz <= 0 or fs <= 0:
        return 31
    A = max(atten_db, 10.0)
    delta_omega = 2.0 * math.pi * (trans_width_hz / fs)
    if delta_omega <= 1e-12:
        return 8191  # 異常に狭い遷移域のときは上限に張り付く
    N_est = (A - 8.0) / (2.285 * delta_omega) + 1.0
    N = int(math.ceil(N_est))
    # 奇数タップにして線形位相の群遅延特性を扱いやすく
    if N % 2 == 0:
        N += 1
    return max(11, min(N, 8191))


def hz_to_nyq_norm(freq_hz: float | np.ndarray, fs: float) -> float | np.ndarray:
    """周波数[Hz]をナイキスト正規化(0..1)に変換 (scipyのWn規格)。"""
    return np.asarray(freq_hz) / (fs / 2.0)

# 単位コンボユーティリティ
UNIT_SPECS = [("Hz", 1.0), ("kHz", 1e3), ("MHz", 1e6)]

def make_unit_cb(default: str = "Hz") -> QComboBox:
    cb = QComboBox()
    for name, mul in UNIT_SPECS:
        cb.addItem(name, mul)
    try:
        idx = [name for name, _ in UNIT_SPECS].index(default)
        cb.setCurrentIndex(idx)
    except Exception:
        pass
    cb.setFixedWidth(72)
    return cb


def unit_mul(cb: QComboBox) -> float:
    data = cb.currentData()
    try:
        return float(data)
    except Exception:
        return 1.0

# 表示用の単位整形
def nice_hz_str(f: float) -> str:
    f = float(abs(f))
    if f >= 9.995e5:
        return f"{f/1e6:.3g} MHz"
    if f >= 9.995e2:
        return f"{f/1e3:.3g} kHz"
    return f"{f:.3g} Hz"


def nice_sec_str(t: float) -> str:
    t = float(abs(t))
    if t >= 1.0:
        return f"{t:.3g} s"
    if t >= 1e-3:
        return f"{t*1e3:.3g} ms"
    if t >= 1e-6:
        return f"{t*1e6:.3g} µs"
    return f"{t*1e9:.3g} ns"


# -------------------------------
# メインUI
# -------------------------------
class FilterDesigner(QMainWindow):
    def __init__(self) -> None:
        super().__init__()
        self.setWindowTitle("デジタルフィルタ設計ツール (FIR / IIR-SOS)")
        self.resize(1300, 820)

        root = QWidget()
        self.setCentralWidget(root)
        root_layout = QHBoxLayout(root)

        # 左側：パラメータ
        left = QWidget()
        left_layout = QVBoxLayout(left)
        left_layout.setSpacing(8)
        root_layout.addWidget(left, 0)

        # 右側：グラフ
        right = QWidget()
        right_layout = QVBoxLayout(right)
        root_layout.addWidget(right, 1)

        # --- パラメータフォーム ---
        self.family_cb = QComboBox()
        self.family_cb.addItems(["FIR (Kaiser)", "IIR (Chebyshev II, SOS)"])
        self.family_cb.currentIndexChanged.connect(self._on_family_changed)

        self.btype_cb = QComboBox()
        self.btype_cb.addItems(["ローパス", "ハイパス", "バンドパス"])
        self.btype_cb.currentIndexChanged.connect(self._on_btype_changed)

        self.fs_sb = QDoubleSpinBox()
        self.fs_sb.setDecimals(1)
        self.fs_sb.setRange(1.0, 1_000_000_000.0)
        self.fs_sb.setValue(48_000.0)
        self.fs_sb.setSuffix("")

        # 通過域端/阻止域端（ローパス・ハイパス用）
        self.fp_sb = QDoubleSpinBox(); self.fp_sb.setRange(0.0, 1_000_000_000.0); self.fp_sb.setDecimals(1); self.fp_sb.setValue(18_000.0); self.fp_sb.setSuffix("")
        self.fs_edge_sb = QDoubleSpinBox(); self.fs_edge_sb.setRange(0.0, 1_000_000_000.0); self.fs_edge_sb.setDecimals(1); self.fs_edge_sb.setValue(20_000.0); self.fs_edge_sb.setSuffix("")

        # バンドパス用
        self.fp1_sb = QDoubleSpinBox(); self.fp1_sb.setRange(0.0, 1_000_000_000.0); self.fp1_sb.setDecimals(1); self.fp1_sb.setValue(500.0); self.fp1_sb.setSuffix("")
        self.fp2_sb = QDoubleSpinBox(); self.fp2_sb.setRange(0.0, 1_000_000_000.0); self.fp2_sb.setDecimals(1); self.fp2_sb.setValue(5_000.0); self.fp2_sb.setSuffix("")
        self.fs1_sb = QDoubleSpinBox(); self.fs1_sb.setRange(0.0, 1_000_000_000.0); self.fs1_sb.setDecimals(1); self.fs1_sb.setValue(300.0); self.fs1_sb.setSuffix("")
        self.fs2_sb = QDoubleSpinBox(); self.fs2_sb.setRange(0.0, 1_000_000_000.0); self.fs2_sb.setDecimals(1); self.fs2_sb.setValue(7_000.0); self.fs2_sb.setSuffix("")

        # 単位選択コンボ
        self.fs_unit = make_unit_cb("Hz")
        self.fp_unit = make_unit_cb("Hz")
        self.fs_edge_unit = make_unit_cb("Hz")
        self.fp1_unit = make_unit_cb("Hz")
        self.fp2_unit = make_unit_cb("Hz")
        self.fs1_unit = make_unit_cb("Hz")
        self.fs2_unit = make_unit_cb("Hz")

        # 減衰量/リプル
        self.stop_atten_db_sb = QDoubleSpinBox(); self.stop_atten_db_sb.setRange(10.0, 300.0); self.stop_atten_db_sb.setDecimals(2); self.stop_atten_db_sb.setValue(80.0); self.stop_atten_db_sb.setSuffix(" dB")
        self.pass_ripple_db_sb = QDoubleSpinBox(); self.pass_ripple_db_sb.setRange(0.01, 5.0); self.pass_ripple_db_sb.setDecimals(2); self.pass_ripple_db_sb.setValue(0.5); self.pass_ripple_db_sb.setSuffix(" dB")

        # FIR 設計制限
        self.max_fir_order_sb = QSpinBox(); self.max_fir_order_sb.setRange(11, 16385); self.max_fir_order_sb.setValue(2049); self.max_fir_order_sb.setSingleStep(2)
        self.force_odd_taps_cb = QCheckBox("FIRタップ数を奇数に強制（推奨）"); self.force_odd_taps_cb.setChecked(True)

        # インパルス応答 表示時間
        self.imp_ms_sb = QDoubleSpinBox(); self.imp_ms_sb.setRange(0.00001, 5000.0); self.imp_ms_sb.setDecimals(5); self.imp_ms_sb.setValue(5.0); self.imp_ms_sb.setSuffix(" ms")
        self.logx_cb = QCheckBox("周波数軸を対数表示 (Bode)"); self.logx_cb.setChecked(True)
        self.unwrap_phase_cb = QCheckBox("位相をアンラップ表示"); self.unwrap_phase_cb.setChecked(True)

        # ボタン
        self.design_btn = QPushButton("設計する")
        self.design_btn.clicked.connect(self.design_filter)
        self.copy_btn = QPushButton("C配列をコピー")
        self.copy_btn.clicked.connect(self.copy_c_arrays_to_clipboard)

        # 出力テキスト
        self.param_text = QTextEdit(); self.param_text.setReadOnly(False)

        # フォーム配置
        spec_group = QGroupBox("仕様")
        form = QFormLayout()
        form.addRow("ファミリ", self.family_cb)
        form.addRow("タイプ", self.btype_cb)
        self.row_fs_with_unit = QWidget(); _layout_fs = QHBoxLayout(self.row_fs_with_unit); _layout_fs.setContentsMargins(0,0,0,0); _layout_fs.addWidget(self.fs_sb); _layout_fs.addWidget(self.fs_unit)
        form.addRow("サンプリング周波数", self.row_fs_with_unit)
        # 単一境界（LP/HP）
        self.row_fp = QWidget(); row_fp_layout = QHBoxLayout(self.row_fp); row_fp_layout.setContentsMargins(0,0,0,0); row_fp_layout.addWidget(self.fp_sb); row_fp_layout.addWidget(self.fp_unit)
        self.row_fs = QWidget(); row_fs_layout = QHBoxLayout(self.row_fs); row_fs_layout.setContentsMargins(0,0,0,0); row_fs_layout.addWidget(self.fs_edge_sb); row_fs_layout.addWidget(self.fs_edge_unit)
        form.addRow("通過域端 f_p", self.row_fp)
        form.addRow("阻止域端 f_s (カットオフ)", self.row_fs)
        # バンドパス（2境界）
        self.row_fp_bp = QWidget(); fp_bp_layout = QGridLayout(self.row_fp_bp); fp_bp_layout.setContentsMargins(0,0,0,0)
        fp_bp_layout.addWidget(QLabel("f_p1"), 0,0); fp_bp_layout.addWidget(self.fp1_sb, 0,1); fp_bp_layout.addWidget(self.fp1_unit, 0,2)
        fp_bp_layout.addWidget(QLabel("f_p2"), 0,3); fp_bp_layout.addWidget(self.fp2_sb, 0,4); fp_bp_layout.addWidget(self.fp2_unit, 0,5)
        self.row_fs_bp = QWidget(); fs_bp_layout = QGridLayout(self.row_fs_bp); fs_bp_layout.setContentsMargins(0,0,0,0)
        fs_bp_layout.addWidget(QLabel("f_s1"), 0,0); fs_bp_layout.addWidget(self.fs1_sb, 0,1); fs_bp_layout.addWidget(self.fs1_unit, 0,2)
        fs_bp_layout.addWidget(QLabel("f_s2"), 0,3); fs_bp_layout.addWidget(self.fs2_sb, 0,4); fs_bp_layout.addWidget(self.fs2_unit, 0,5)
        form.addRow("通過域端 (BP)", self.row_fp_bp)
        form.addRow("阻止域端 (BP)", self.row_fs_bp)
        form.addRow("阻止域減衰量 A_s", self.stop_atten_db_sb)
        form.addRow("通過域リプル A_p (IIR用)", self.pass_ripple_db_sb)
        form.addRow("インパルス表示時間", self.imp_ms_sb)
        form.addRow("周波数軸(ゲイン/位相)", self.logx_cb)
        form.addRow("位相表示", self.unwrap_phase_cb)

        # FIR詳細
        fir_group = QGroupBox("FIR設計オプション (Kaiser)")
        fir_form = QFormLayout()
        fir_form.addRow("最大タップ数", self.max_fir_order_sb)
        fir_form.addRow(self.force_odd_taps_cb)
        fir_group.setLayout(fir_form)

        # ボタン行
        btn_row = QWidget(); btn_l = QHBoxLayout(btn_row); btn_l.setContentsMargins(0,0,0,0)
        btn_l.addWidget(self.design_btn); btn_l.addWidget(self.copy_btn)

        # 出力グループ
        out_group = QGroupBox("フィルタ係数 (C言語配列)／設計サマリ")
        out_layout = QVBoxLayout(out_group)
        out_layout.addWidget(self.param_text)

        spec_group.setLayout(form)
        left_layout.addWidget(spec_group)
        left_layout.addWidget(fir_group)
        left_layout.addWidget(btn_row)
        left_layout.addWidget(out_group)
        left_layout.addStretch(1)

        # --- プロット領域 ---
        self.fig = Figure(figsize=(6, 6), tight_layout=True)
        self.canvas = FigureCanvas(self.fig)
        right_layout.addWidget(self.canvas)
        self.ax_mag = self.fig.add_subplot(3,1,1)
        self.ax_phase = self.fig.add_subplot(3,1,2)
        self.ax_imp = self.fig.add_subplot(3,1,3)
        self.ax_mag.set_ylabel("Gain [dB]")
        self.ax_phase.set_ylabel("Phase [deg]")
        self.ax_imp.set_ylabel("Amplitude")
        self.ax_imp.set_xlabel("Time [s]")
        # ツールチップ用マウストラッキング
        self.canvas.mpl_connect('motion_notify_event', self._on_motion)

        self._on_btype_changed()
        self._on_family_changed()

        # 内部状態
        self.last_coeffs_fir: np.ndarray | None = None
        self.last_sos_iir: np.ndarray | None = None

    # ---------- UI 切替 ----------
    def _on_family_changed(self) -> None:
        is_fir = (self.family_cb.currentIndex() == 0)
        self.pass_ripple_db_sb.setEnabled(not is_fir)  # IIR のみ有効
        # FIRオプションの有効/無効
        self.force_odd_taps_cb.setEnabled(is_fir)
        self.max_fir_order_sb.setEnabled(is_fir)

    def _on_btype_changed(self) -> None:
        btype = self.btype_cb.currentText()
        is_bp = (btype == "バンドパス")
        # LP/HP 用：単一境界表示
        self.row_fp.setVisible(not is_bp)
        self.row_fs.setVisible(not is_bp)
        # BP 用：2境界表示
        self.row_fp_bp.setVisible(is_bp)
        self.row_fs_bp.setVisible(is_bp)

    # ---------- 設計メイン ----------
    def design_filter(self) -> None:
        try:
            fs = float(self.fs_sb.value()) * unit_mul(self.fs_unit)
            btype = self.btype_cb.currentText()
            family_is_fir = (self.family_cb.currentIndex() == 0)
            Astop = float(self.stop_atten_db_sb.value())
            Apass = float(self.pass_ripple_db_sb.value())

            if btype != "バンドパス":
                fp = float(self.fp_sb.value()) * unit_mul(self.fp_unit)
                fs_edge = float(self.fs_edge_sb.value()) * unit_mul(self.fs_edge_unit)
                if not self._validate_edges_lp_hp(fp, fs_edge, fs, btype):
                    return
            else:
                fp1 = float(self.fp1_sb.value()) * unit_mul(self.fp1_unit); fp2 = float(self.fp2_sb.value()) * unit_mul(self.fp2_unit)
                fs1 = float(self.fs1_sb.value()) * unit_mul(self.fs1_unit); fs2 = float(self.fs2_sb.value()) * unit_mul(self.fs2_unit)
                if not self._validate_edges_bp(fp1, fp2, fs1, fs2, fs):
                    return

            # 設計
            if family_is_fir:
                if btype != "バンドパス":
                    coeffs = self._design_fir_lp_hp(fs, btype, fp, fs_edge, Astop)
                else:
                    coeffs = self._design_fir_bp(fs, fp1, fp2, fs1, fs2, Astop)
                self.last_coeffs_fir = coeffs
                self.last_sos_iir = None
                self._update_plots_fir(coeffs, fs)
                self._update_text_fir(coeffs, fs, btype, Astop)
            else:
                if btype != "バンドパス":
                    sos = self._design_iir_lp_hp(fs, btype, fp, fs_edge, Apass, Astop)
                else:
                    sos = self._design_iir_bp(fs, fp1, fp2, fs1, fs2, Apass, Astop)
                self.last_coeffs_fir = None
                self.last_sos_iir = sos
                self._update_plots_iir(sos, fs)
                self._update_text_iir(sos, fs, btype, Apass, Astop)
        except Exception as e:
            QMessageBox.critical(self, "設計エラー", f"{type(e).__name__}: {e}")

    # ---------- 入力検証 ----------
    def _validate_edges_lp_hp(self, fp: float, fs_edge: float, fs: float, btype_text: str) -> bool:
        if fs <= 0:
            QMessageBox.warning(self, "入力エラー", "サンプリング周波数は正である必要があります。")
            return False
        nyq = fs/2.0
        if not (0 < fp < nyq) or not (0 < fs_edge < nyq):
            QMessageBox.warning(self, "入力エラー", "f_p と f_s は (0, Nyquist) に収めてください。")
            return False
        if btype_text == "ローパス":
            if not (fs_edge > fp):
                QMessageBox.warning(self, "入力エラー", "ローパスでは f_s > f_p となるようにしてください。")
                return False
        elif btype_text == "ハイパス":
            if not (fs_edge < fp):
                QMessageBox.warning(self, "入力エラー", "ハイパスでは f_s < f_p となるようにしてください。")
                return False
        return True

    def _validate_edges_bp(self, fp1: float, fp2: float, fs1: float, fs2: float, fs: float) -> bool:
        if fs <= 0:
            QMessageBox.warning(self, "入力エラー", "サンプリング周波数は正である必要があります。")
            return False
        nyq = fs/2.0
        vals = [fp1, fp2, fs1, fs2]
        if any(v <= 0 or v >= nyq for v in vals):
            QMessageBox.warning(self, "入力エラー", "周波数は (0, Nyquist) に収めてください。")
            return False
        if not (fp1 < fp2 and fs1 < fp1 and fs2 > fp2):
            QMessageBox.warning(self, "入力エラー", "バンドパスは fs1 < fp1 < fp2 < fs2 となるようにしてください。")
            return False
        return True

    # ---------- FIR 設計 ----------
    def _design_fir_lp_hp(self, fs: float, btype_text: str, fp: float, fs_edge: float, Astop: float) -> np.ndarray:
        # 遷移域と中心カットオフを決定
        trans = abs(fs_edge - fp)
        fc = 0.5 * (fp + fs_edge)
        N = kaiser_order(Astop, trans, fs)
        # 上限・奇数化
        N = min(N, int(self.max_fir_order_sb.value()))
        if self.force_odd_taps_cb.isChecked() and N % 2 == 0:
            N = max(11, N-1)
        beta = kaiser_beta_for_atten(Astop)
        # firwin を使う（fs 指定でHz単位）
        pass_zero = (btype_text == "ローパス")
        taps = firwin(numtaps=N, cutoff=fc, window=("kaiser", beta), pass_zero=pass_zero, fs=fs)
        return taps.astype(np.float64)

    def _design_fir_bp(self, fs: float, fp1: float, fp2: float, fs1: float, fs2: float, Astop: float) -> np.ndarray:
        trans1 = abs(fp1 - fs1)
        trans2 = abs(fs2 - fp2)
        trans = max(min(trans1, trans2), 1e-6)
        fc1 = 0.5 * (fp1 + fs1)
        fc2 = 0.5 * (fp2 + fs2)
        N = kaiser_order(Astop, trans, fs)
        N = min(N, int(self.max_fir_order_sb.value()))
        if self.force_odd_taps_cb.isChecked() and N % 2 == 0:
            N = max(11, N-1)
        beta = kaiser_beta_for_atten(Astop)
        taps = firwin(numtaps=N, cutoff=[fc1, fc2], window=("kaiser", beta), pass_zero=False, fs=fs)
        return taps.astype(np.float64)

    # ---------- IIR 設計 (Chebyshev II) ----------
    def _design_iir_lp_hp(self, fs: float, btype_text: str, fp: float, fs_edge: float, Apass: float, Astop: float) -> np.ndarray:
        wp = hz_to_nyq_norm(fp, fs)
        ws = hz_to_nyq_norm(fs_edge, fs)
        btype = 'lowpass' if btype_text == "ローパス" else 'highpass'
        N, wn = cheb2ord(wp, ws, gpass=Apass, gstop=Astop)
        sos = cheby2(N, rs=Astop, Wn=wn, btype=btype, output='sos')
        return sos.astype(np.float64)

    def _design_iir_bp(self, fs: float, fp1: float, fp2: float, fs1: float, fs2: float, Apass: float, Astop: float) -> np.ndarray:
        wp = hz_to_nyq_norm([fp1, fp2], fs)
        ws = hz_to_nyq_norm([fs1, fs2], fs)
        N, wn = cheb2ord(wp, ws, gpass=Apass, gstop=Astop)
        sos = cheby2(N, rs=Astop, Wn=wn, btype='bandpass', output='sos')
        return sos.astype(np.float64)

    # ---------- プロット更新 ----------
    def _update_plots_fir(self, taps: np.ndarray, fs: float) -> None:
        w, h = freqz(taps, worN=4096, fs=fs)
        self._update_plots_common(w, h, fs, lambda n: lfilter(taps, [1.0], n))

    def _update_plots_iir(self, sos: np.ndarray, fs: float) -> None:
        w, h = sosfreqz(sos, worN=4096, fs=fs)
        from scipy.signal import sosfilt
        self._update_plots_common(w, h, fs, lambda n: sosfilt(sos, n))

    def _update_plots_common(self, w: np.ndarray, H: np.ndarray, fs: float, filt_fn) -> None:
        # ゲイン・位相（対数軸オプション対応）
        log_x = bool(self.logx_cb.isChecked())
        fmin = max(1e-3, fs / 1e6)  # 0Hzはlog不可のためクリップ
        mask = w >= (fmin if log_x else 0.0)
        w_plot = w[mask]
        H_plot = H[mask]
        mag_db = 20.0 * np.log10(np.maximum(np.abs(H_plot), 1e-15))
        phi = np.angle(H_plot)
        if self.unwrap_phase_cb.isChecked():
            phi = np.unwrap(phi)
        phase_deg = np.degrees(phi)

        self.ax_mag.clear(); self.ax_phase.clear(); self.ax_imp.clear()
        self.ax_mag.set_xscale('log' if log_x else 'linear')
        self.ax_phase.set_xscale('log' if log_x else 'linear')
        self.ax_mag.plot(w_plot, mag_db)
        self.ax_mag.set_ylabel("Gain [dB]")
        self.ax_mag.grid(True, which='both', linestyle=':')

        self.ax_phase.plot(w_plot, phase_deg)
        self.ax_phase.set_ylabel("Phase [deg]")
        self.ax_phase.grid(True, which='both', linestyle=':')

        # インパルス応答
        Tms = float(self.imp_ms_sb.value())
        Nimp = int(max(2, min(1_000_000, round((Tms/1000.0)*fs))))
        x = np.zeros(Nimp, dtype=np.float64); x[0] = 1.0
        y = filt_fn(x)
        t = np.arange(Nimp) / fs
        self.ax_imp.plot(t, y)
        self.ax_imp.set_ylabel("Amplitude")
        self.ax_imp.set_xlabel("Time [s]")
        self.ax_imp.grid(True, linestyle=':')

        # カーソルラインを各Axesに用意
        self._ensure_cursor(self.ax_mag)
        self._ensure_cursor(self.ax_phase)
        self._ensure_cursor(self.ax_imp)

        self.fig.tight_layout()
        self.canvas.draw_idle()

    # ---------- 出力テキスト生成 ----------
    def _update_text_fir(self, taps: np.ndarray, fs: float, btype_text: str, Astop: float) -> None:
        N = len(taps)
        imp_ms = float(self.imp_ms_sb.value())
        beta = kaiser_beta_for_atten(Astop)
        c_array = self._format_c_array_fir(taps)
        summary = (
            f"//[FIR (Kaiser)]\n"
            f"//タイプ: {btype_text}\n"
            f"//Fs = {fs:.3f} Hz, タップ数 N = {N}, beta = {beta:.3f}, A_s = {Astop:.2f} dB\n"
        )
        self.param_text.setPlainText(summary + c_array)

    def _update_text_iir(self, sos: np.ndarray, fs: float, btype_text: str, Apass: float, Astop: float) -> None:
        nsec = sos.shape[0]
        imp_ms = float(self.imp_ms_sb.value())
        c_sos = self._format_c_array_sos(sos)
        c_cmsis = self._format_c_array_cmsis_biquad(sos)
        summary = (
            f"//[IIR (Chebyshev II, SOS)]\n"
            f"//タイプ: {btype_text}\n"
            f"//Fs = {fs:.3f} Hz, セクション数 = {nsec}, A_p = {Apass:.2f} dB, A_s = {Astop:.2f} dB\n"
        )
        self.param_text.setPlainText(summary + c_sos + "\n" + c_cmsis)

    # ---------- C配列生成 ----------
    def _format_c_array_fir(self, taps: np.ndarray, varname: str = "fir_taps") -> str:
        values = ",\n    ".join(f"{v:.9e}f" for v in taps.astype(np.float32))
        return (
            f"// FIR taps (float32), length = {len(taps)}\n"
            f"static const float {varname}[{len(taps)}] = {{\n    {values}\n}};\n"
        )

    def _format_c_array_sos(self, sos: np.ndarray, varname: str = "iir_sos") -> str:
        lines = []
        for s in sos:
            b0,b1,b2,a0,a1,a2 = s
            # 念のため a0 を 1 に正規化
            if a0 == 0:
                a0 = 1.0
            b0,b1,b2,a1,a2 = b0/a0, b1/a0, b2/a0, a1/a0, a2/a0
            lines.append(f"    {{ {b0:.9e}f, {b1:.9e}f, {b2:.9e}f, 1.000000000e+00f, {a1:.9e}f, {a2:.9e}f }}")
        body = ",\n".join(lines)
        return (
            f"// IIR SOS (float32) as rows: {{b0,b1,b2,a0,a1,a2}} per section\n"
            f"static const float {varname}[{sos.shape[0]}][6] = {{\n{body}\n}};\n"
        )

    def _format_c_array_cmsis_biquad(self, sos: np.ndarray, varname: str = "biquad_coeffs") -> str:
        """CMSIS-DSP arm_biquad_cascade_df1_f32 互換の係数列
        各セクション: {b0, b1, b2, -a1, -a2}
        ゲインは b0 に乗る（SciPyのSOSは全体ゲインを先頭段に含む）。
        """
        coeffs = []
        for s in sos:
            b0,b1,b2,a0,a1,a2 = s
            if a0 == 0:
                a0 = 1.0
            b0,b1,b2,a1,a2 = b0/a0, b1/a0, b2/a0, a1/a0, a2/a0
            coeffs.extend([b0, b1, b2, -a1, -a2])
        values = ",\n    ".join(f"{v:.9e}f" for v in np.asarray(coeffs, dtype=np.float32))
        return (
            f"// CMSIS-DSP biquad coeffs (float32), 5 * {sos.shape[0]} elements\n"
            f"static const float {varname}[{5*sos.shape[0]}] = {{\n    {values}\n}};\n"
        )

    # ---------- マウストラッキング / ツールチップ ----------
    def _ensure_cursor(self, ax):
        if not hasattr(self, "_cursor_lines"):
            self._cursor_lines = {}
        pair = self._cursor_lines.get(ax)
        valid = False
        if pair is not None:
            v, h = pair
            valid = (getattr(v, 'axes', None) is ax) and (getattr(h, 'axes', None) is ax)
        if not valid:
            v = ax.axvline(np.nan, linestyle='--', linewidth=1.0, alpha=0.6)
            h = ax.axhline(np.nan, linestyle='--', linewidth=1.0, alpha=0.6)
            self._cursor_lines[ax] = (v, h)

    def _hide_all_cursors(self):
        if not hasattr(self, "_cursor_lines"):
            return
        for ax, (v, h) in self._cursor_lines.items():
            try:
                v.set_xdata([np.nan, np.nan])
                h.set_ydata([np.nan, np.nan])
            except Exception:
                pass
        self.canvas.draw_idle()

    def _on_motion(self, event) -> None:
        try:
            if event.inaxes is None or event.xdata is None:
                QToolTip.hideText();
                self._hide_all_cursors();
                return
            ax = event.inaxes
            if not ax.lines:
                QToolTip.hideText(); return
            line = ax.lines[0]
            xdata = np.asarray(line.get_xdata())
            ydata = np.asarray(line.get_ydata())
            if xdata.size < 2:
                QToolTip.hideText(); return
            x = float(event.xdata)
            xmin, xmax = float(np.min(xdata)), float(np.max(xdata))
            x = float(np.clip(x, xmin, xmax))
            y = float(np.interp(x, xdata, ydata))

            if ax is self.ax_mag:
                text = f"f = {nice_hz_str(x)}|H| = {y:.2f} dB"
            elif ax is self.ax_phase:
                text = f"f = {nice_hz_str(x)}∠H = {y:.1f}°"
            elif ax is self.ax_imp:
                text = f"t = {nice_sec_str(x)}h = {y:.6g}"
            else:
                QToolTip.hideText(); return
            QToolTip.showText(QCursor.pos(), text, self.canvas)

            # カーソルライン更新
            self._ensure_cursor(ax)
            if hasattr(self, "_cursor_lines"):
                # 他軸は消す
                for a, (vl, hl) in list(self._cursor_lines.items()):
                    if a is not ax:
                        try:
                            vl.set_xdata([np.nan, np.nan])
                            hl.set_ydata([np.nan, np.nan])
                        except Exception:
                            pass
                vl, hl = self._cursor_lines[ax]
                try:
                    vl.set_xdata([x, x])
                    hl.set_ydata([y, y])
                except Exception:
                    pass
                self.canvas.draw_idle()
        except Exception:
            # 失敗時は静かに無視
            pass

    # ---------- クリップボード ----------
    def copy_c_arrays_to_clipboard(self) -> None:
        text = self.param_text.toPlainText()
        QApplication.clipboard().setText(text, mode=QClipboard.Clipboard)
        QMessageBox.information(self, "コピー完了", "係数とサマリをクリップボードにコピーしました。")


# -------------------------------
# エントリポイント
# -------------------------------
if __name__ == '__main__':
    app = QApplication(sys.argv)
    w = FilterDesigner()
    w.show()
    sys.exit(app.exec())
