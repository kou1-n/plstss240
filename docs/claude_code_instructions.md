
# Claude Code 実行用 指示ファイル — DP 非関連塑性の一貫接線を**再検討**し**検証**する

最終更新: 2025-09-30 (JST)

---

## 目的
Drucker–Prager（DP）**非関連**塑性モデルにおける、あなたの実装（`stress_dp_rm` 系）の**一貫接線（consistent tangent）**と**塑性乗数 Δγ**の式が正しいかを**再考**し、再現可能なコード検証で確認する。特に以下を重点確認する：

1. 一貫接線の一般式  
   \[ \boxed{ \mathbf D_{\mathrm{ep}} = \mathbf D_e \;-\; \dfrac{(\mathbf D_e:\mathbf B)\;\bar\otimes\;(\mathbf A:\mathbf D_e)}{\mathbf A:\mathbf D_e:\mathbf B \;+\; H} } \]
   ただし \(\mathbf A=\partial f/\partial\boldsymbol\sigma\), \(\mathbf B=\partial g/\partial\boldsymbol\sigma\), \(H=-\dfrac{\partial f}{\partial \alpha}\dfrac{\partial \alpha}{\partial \gamma}\)。

2. ふぁい
   \[ \mathbf A:\mathbf D_e:\mathbf B = G + K\,\eta\,\bar\eta,\quad
      H = \xi \, \dfrac{\partial K}{\partial \alpha} \, \dfrac{\partial \alpha}{\partial \gamma},\quad
      \dfrac{\partial \alpha}{\partial \gamma} = \dfrac{1}{\sqrt{3}}. \]
   **誤りが生じやすい点**：硬化項を \(\xi^2\) としてしまうミス（**正しくは一次**）。

3. **非対称接線**であること（\(\mathbf D_{\mathrm{ep}} \neq \mathbf D_{\mathrm{ep}}^\mathsf T\)）の確認。

4. **局所戻り写像**の Newton 反復が、上式の \(\mathbf D_{\mathrm{ep}}\) 採用時に**二次収束**することの定量確認。

---

## 前提と記号（本指示内で統一）
- 応力分解：\(\boldsymbol\sigma=\mathbf s + p\mathbf I\)、\(p=\tfrac{1}{3}\mathrm{tr}(\boldsymbol\sigma)\)、\(\mathbf s\) は偏差応力。
- 等方弾性：\(\mathbf D_e = 2G\,\mathbf P_{\mathrm{dev}} + K\,\mathbf I\otimes\mathbf I\)。
- DP の降伏・塑性ポテンシャル（同一規約を採用）
  \[ f=\frac{1}{\sqrt{2}}\|\mathbf s\| + \eta\,p - \xi\big(\sigma_Y + K(\alpha)\big),\quad
     g=\frac{1}{\sqrt{2}}\|\mathbf s\| + \bar\eta\,p. \]
- 勾配：\(\mathbf A=\frac{1}{\sqrt2}\mathbf n+\frac{\eta}{3}\mathbf I\), \(\mathbf B=\frac{1}{\sqrt2}\mathbf n+\frac{\bar\eta}{3}\mathbf I\)（\(\mathbf n=\mathbf s/\|\mathbf s\|\)）。
- 等価塑性ひずみの規約：\(\Delta\alpha=\Delta\gamma/\sqrt{3}\Rightarrow \partial\alpha/\partial\gamma=1/\sqrt{3}\)。

---

## **再考の結論（要約）**
- **分母の硬化項**は \(\boxed{H=\xi\,(\partial K/\partial\alpha)\,(1/\sqrt{3})}\)。  
  以前の \(\xi^2(\partial K/\partial\alpha)\) という形は**誤り**で、次元・収束性の両面で不整合を招く。
- 非関連（\(\eta\neq\bar\eta\)）では \(\mathbf D_{\mathrm{ep}}\) は**非対称**。グローバル Newton で二次収束を得るには**非対称ソルバ**が必要。
- 放射戻り（\(\mathbf n\) 凍結、円錐側面）を前提にした**アルゴリズム接線**として上式を使えば、局所 Newton は**二次収束**する（頂点は別分岐）。

---

## Claude Code への**作業依頼**（自動化タスク）

> **実行順**に上から処理してください。エラー時はログとともに継続不能箇所を報告。

### 1) リポジトリ準備
- 新規ローカルプロジェクトを作成（Python 3.10+、NumPy 必須）。仮想環境可。
- ルート構成：
  ```text
  dp_consistent_tangent/
    src/materials/elastic.py
    src/materials/dp_nonassoc.py
    src/tests/test_consistent_tangent.py
    src/tests/test_local_return_convergence.py
    scripts/run_all.py
    pyproject.toml  (or requirements.txt: numpy)
    README.md
  ```

### 2) 等方弾性ユーティリティ `src/materials/elastic.py`
- 目的：4階テンソル操作の最小限関数群を提供。
- 仕様（NumPy, フル 3D）：
  - `voigt6(sym_tensor_3x3) -> R6`
  - `de_voigt6(R6) -> sym_tensor_3x3`
  - `elastic_stiffness(E, nu) -> De(6x6)` （Voigt 6x6）
  - `proj_dev_6() -> Pdev (6x6)`
  - `split_dev_hydro(sigma) -> (s_dev(3x3), p)`
  - `unit_dev_dir(s_dev) -> n(3x3)`（‖n‖=1、s≠0 でのみ使用）

### 3) DP 非関連材料ルーチン `src/materials/dp_nonassoc.py`
- **関数**：
  - `yield_f(sigma, alpha, params) -> f`
  - `grad_f(sigma, params) -> A (6,)`（Voigt）
  - `grad_g(sigma, params) -> B (6,)`（Voigt）
  - `trial_state(sigma_n, deps, De) -> sigma_tr`
  - `return_map_nonassoc(sigma_tr, alpha_n, params) -> (sigma, alpha, dgamma, iter, converged)`  
    - Radial return（側面）を想定。\(\mathbf n\) は試行方向で凍結（s≠0）
    - \(\Delta\gamma = f^{\mathrm{tr}} / (G + K\eta\bar\eta + \xi(\partial K/\partial\alpha)/\sqrt{3})\)
    - 応力更新：\(\sigma = \sigma^{\mathrm{tr}} - \Delta\gamma \, ( \mathbf D_e:B)\)
    - \(\alpha \leftarrow \alpha + \Delta\gamma/\sqrt{3}\)
  - `consistent_tangent_nonassoc(sigma, params, De) -> Dep(6x6)`  
    \[ \mathbf D_{\mathrm{ep}} = \mathbf D_e - \dfrac{(\mathbf D_e B)\otimes(A^\top \mathbf D_e)}{A^\top \mathbf D_e B + H} \]
    - \(H=\xi(\partial K/\partial\alpha)/\sqrt{3}\) とする。
- **パラメータ**：`params = dict(E, nu, G, K, eta, bar_eta, xi, sigma_y, K_of_alpha, dK_dalpha)`  
  - `G,K` は `E,nu` から計算してもよい。
  - `K_of_alpha(alpha)` と `dK_dalpha(alpha)` はコールバック（例：線形/指数硬化）。

### 4) 数値差分との一致検証 `src/tests/test_consistent_tangent.py`
- **目的**：\(\mathbf D_{\mathrm{ep}}\) が**解析解**と**数値差分**で一致するか検証。
- **手順**：
  1. ランダムな弾性定数と DP 定数（例：E=70e9, nu=0.3, eta=0.2, bar_eta=0.1, xi=1.0, 等）を設定。
  2. ランダムな \(\sigma_n\) と \(\Delta\varepsilon\) を生成し試行状態へ。
  3. 塑性が発生するようスケール調整（`f_tr>0`）。
  4. 解析 \(\mathbf D_{\mathrm{ep}}\) を取得。
  5. **数値ヤコビアン**：\(\partial \sigma/\partial \varepsilon\) を中心差分（擾乱 \(h=10^{-8}\sim10^{-6}\)）で計算。
  6. Frobenius ノルム相対誤差 \(\|\mathrm{Dep}_{\mathrm{ana}}-\mathrm{Dep}_{\mathrm{num}}\|/\|\mathrm{Dep}_{\mathrm{num}}\|\) を評価。
- **合格基準**：`< 1e-6`（倍精度想定、問題により `1e-5` まで許容）。
- **併せて確認**：
  - `np.allclose(Dep, Dep.T, rtol=..., atol=...)` が **偽**（=非対称）であること。
  - **関連化**テスト：\(\eta=\bar\eta\) とし、同条件で対称性が向上（`‖Dep-Dep.T‖` が数値的ゼロ付近）すること。

### 5) 局所 Newton 収束率の検証 `src/tests/test_local_return_convergence.py`
- **目的**：一貫接線の**正しさ**を収束次数で裏付け。
- **手順**：
  1. 同じ条件で `return_map_nonassoc` を Newton（残差：`[r_sigma; r_f]`）で実装（解析接線を使用）。
  2. 収束履歴 \(\|R_k\|\) をログ取りし、反復毎の比 \(\|R_{k+1}\|/\|R_k\|^2\) を評価。
  3. **比較実験**：硬化項を**誤って** \(\xi^2(\partial K/\partial\alpha)\) にした場合も走らせ、収束が**線形化/悪化**することを示す。
- **合格基準**：正しい接線で**二次収束**（最終 2–3 反復で二乗則）。誤接線では**明確に劣化**。

### 6) 端部（頂点）近傍の健全性チェック（オプション）
- \(\|\mathbf s^{\mathrm{tr}}\|\to 0\) のケースで `n` が不定。閾値 \(\epsilon_n\) を設けて**体積塑性戻り**に切替える分岐を実装し、数値的安定を確認。

### 7) 圧縮・引張の符号規約の整合
- \(p=\tfrac13\mathrm{tr}(\sigma)\) の定義で、外力の符号を切替えても `f` の**一貫性**が保たれることを、2–3 ケースで検証。

### 8) `scripts/run_all.py`
- 上記テストを一括実行し、結果を標準出力に要約（誤差、対称性、収束次数）。
- すべて**合格**なら `ALL TESTS PASSED` と表示。

### 9) `README.md`
- 数式（本ファイルの要点）と、実装・テストの実行方法を記載：
  ```bash
  python -m pip install -r requirements.txt   # or uv/pdm/poetry
  python scripts/run_all.py
  ```

---

## 受け入れ判定（Acceptance Criteria）
1. `test_consistent_tangent.py` の相対誤差 < `1e-6` を複数乱数シードで達成。  
2. 非関連時に \(\|\mathrm{Dep}-\mathrm{Dep}^\mathsf T\|\) が十分な非ゼロ。関連化でほぼゼロ化。  
3. 局所 Newton の最終段で**二次収束**。誤接線では二次にならない。  
4. 頂点近傍・符号規約のテストがクラッシュせず所期の分岐で安定。

---

## 付記：よくある実装ミスのチェックポイント
- \(H\) を \(\xi^2\) にしてしまう（**本検証の主眼**）。
- \(\partial\alpha/\partial\gamma\) の規約と `K_of_alpha, dK_dalpha` の次元不整合。
- \(\mathbf D_e B\) と \(A^\top \mathbf D_e\) の**順序**ミス（外積の左右に注意）。
- Voigt 6 成分の**せん断係数**（工学せん断 vs 物理せん断）混同。すべて **工学 Voigt** で統一。
- 近似的に対称化してしまう（非関連では二次収束を失う）。

---

## 参考式（再掲）
- \(\mathbf D_e:\mathbf B=\sqrt2 G\,\mathbf n + K\,\bar\eta\,\mathbf I\)（Voigt では適切に射影）。  
- \(\mathbf A:\mathbf D_e=\sqrt2 G\,\mathbf n + K\,\eta\,\mathbf I\)。  
- \(\mathbf A:\mathbf D_e:\mathbf B=G+K\eta\bar\eta\)。  
- \(H=\xi(\partial K/\partial\alpha)/\sqrt3\)。  
- \(\mathbf D_{\mathrm{ep}}=\mathbf D_e-(\mathbf D_e B)\otimes(A^\top \mathbf D_e)/(A^\top \mathbf D_e B+H)\)。

---

## 実行ログの提出
- `scripts/run_all.py` の全文ログを保存して提出。
- 主要ケースの `Dep_ana` と `Dep_num` の差のヒートマップ（任意）も保存可。

以上。これに従ってコードを生成・実行し、**すべての Acceptance Criteria を満たすこと**。
