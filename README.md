# PySide6 デジタルフィルタ設計ツール (FIR / IIR-SOS)

PySide6 + SciPy 製の**デジタルフィルタ設計 GUI**です。FIR(カイザー窓) と IIR(Chebyshev II, SOS) を切り替え、ローパス/ハイパス/バンドパスのいずれかを設計できます。設計結果として **ゲイン/位相線図、インパルス応答、C言語配列の係数** を出力します。

> **ライセンス**: MIT（本リポジトリ内の `LICENSE` 参照）

---

## 特長

* **FIR (Kaiser)** / **IIR (Chebyshev II, SOS)** を切替
* ローパス / ハイパス / バンドパス対応
* **単位選択**: Fs, f<sub>p</sub>/f<sub>s</sub> を **Hz / kHz / MHz** から選択可能
* **Bode 表示**: ゲイン/位相の周波数軸を **log/lin 切替**
* **位相表示**: ラップ (±180°) / アンラップ切替
* **インパルス応答**: 表示時間(ms)を任意指定
* **カーソル＆ツールチップ**: マウス位置の **f(または t)** と **|H| / ∠H / h** を即時表示。縦線/横線のカーソルラインを描画
* **C言語配列出力**:

  * FIR: `static const float fir_taps[N]`  (float32)
  * IIR(SOS): `static const float iir_sos[sec][6]`  (各行 `{b0,b1,b2,1,a1,a2}`)
  * CMSIS-DSP 互換: `static const float biquad_coeffs[5*sec] = {b0,b1,b2,-a1,-a2,...}`
* **クリップボード**にサマリと配列を一括コピー

---

## 動作環境

* Python 3.10+ 推奨
* 主要依存: `numpy`, `scipy`, `matplotlib`, `PySide6`

```bash
pip install numpy scipy matplotlib PySide6
```

---

## 使い方

1. リポジトリを取得

   ```bash
   git clone https://github.com/<your-account>/<your-repo>.git
   cd <your-repo>
   ```
2. 依存をインストール（上記参照）
3. アプリ起動

   ```bash
   python filter_designer.py
   ```
4. 左ペインで仕様を設定し、**\[設計する]** を押す
5. 右ペインに **ゲイン\[dB] / 位相\[deg] / インパルス応答** が表示されます
6. 下部テキストに **C配列** と設計サマリが出力されます（**\[C配列をコピー]** でクリップボードへ）

---

## UI 詳説

### 設計ファミリ

* **FIR (Kaiser)**: 線形位相。阻止域減衰 A<sub>s</sub> と遷移帯域幅からタップ数を推定します
* **IIR (Chebyshev II, SOS)**: Stopband に等リプル。SOS 表現で安定性良好

### フィルタタイプ

* ローパス / ハイパス: f<sub>p</sub>（通過域端）, f<sub>s</sub>（阻止域端）
* バンドパス: f<sub>p1</sub>, f<sub>p2</sub>, f<sub>s1</sub>, f<sub>s2</sub>
  **関係**:  f<sub>s1</sub> < f<sub>p1</sub> < f<sub>p2</sub> < f<sub>s2</sub>

### 単位選択 (Hz/kHz/MHz)

* Fs および各境界周波数の右側で選択可能。内部では自動的に **Hz** に換算して設計します

### スペック

* **阻止域減衰 A<sub>s</sub> \[dB]**: FIR の β／タップ数推定、IIR の `cheb2ord`/`cheby2` に利用
* **通過域リプル A<sub>p</sub> \[dB]**: IIR のみ（`cheb2ord` の gpass）
* **FIRオプション**: 最大タップ数、**奇数タップ強制**（推奨）
* **インパルス表示時間**: ms 指定
* **周波数軸**: log/linear 切替
* **位相表示**: ラップ or アンラップ

### グラフ操作

* マウス移動で **ツールチップ**（f/t と Y 値）表示
* 同時に **縦線/横線**のカーソルラインを描画（グラフ外に出ると自動で非表示）

---

## 出力される C 配列

### FIR (Kaiser)

```c
// FIR taps (float32), length = N
static const float fir_taps[N] = {
    /* 係数 ... */
};
```

> `firwin` により線形位相の係数を生成。既定で **奇数タップ**。

### IIR (SOS)

```c
// IIR SOS (float32) SはIIRフィルタの段数
static const float iir_sos[S][6] = {
    // {b0, b1, b2, 1.0, a1, a2}  (a0正規化済み)
};
```

### CMSIS-DSP 互換 (biquad df1)

```c
// arm_biquad_cascade_df1_f32 用
static const float biquad_coeffs[5*S] = {
    // per section: b0, b1, b2, -a1, -a2
};
// ワークバッファ例: float biquad_state[4*S] = {0};
```

---

## 設計ノート

### FIR タップ数推定

* 近似式:
  $N \approx \frac{A_s - 8}{2.285\,\Delta\omega} + 1$  ($\Delta\omega = 2\pi\,\Delta f/\mathrm{Fs}$)
* β 推定: 21 dB 以下 → 0、50 dB 以下 → `0.5842(A-21)^0.4 + 0.07886(A-21)`、それ以外 → `0.1102(A-8.7)`

### IIR (Chebyshev II)

* `cheb2ord(wp, ws, gpass=A_p, gstop=A_s)` で次数・境界を推定し、`cheby2(..., output='sos')` で実数 SOS に展開
* SciPy の規約で周波数は **ナイキスト正規化** (Hz→`/ (Fs/2)`) を使用

---

## 既知の制限

* **極端に狭い遷移帯域**では FIR タップ数が非常に大きくなります（上限: 16385）。
* IIR は **Stopband 等リプル**で、**通過域は単調**だが **位相は非線形**です。
* 係数は `float32` で C 配列化。固定小数点(Qフォーマット)は未実装。
* グラフの画像保存や係数の `.h`/`.csv` への**ファイル保存**は未実装（今後対応予定）。

---

## トラブルシューティング

* `ModuleNotFoundError`: 依存をインストールしてください

  ```bash
  pip install numpy scipy matplotlib PySide6
  ```

---

## 開発

* 単一ファイル: `filter_designer.py`
* 主要依存 API:

  * SciPy: `firwin`, `cheb2ord`, `cheby2`, `freqz`, `sosfreqz`, `sosfilt`
  * Matplotlib QtAgg: 埋め込みキャンバス
  * PySide6: Qt Widgets

PR / Issue 歓迎です！

---

## ライセンス

MIT License

```
MIT License

Copyright (c) 2025 <Your Name>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```

---

## 謝辞

* SciPy / NumPy / Matplotlib / Qt の各コミュニティに感謝します。
