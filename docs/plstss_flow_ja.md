# PLSTss ワークフロー概説

このドキュメントでは，PLSTss プログラムが入力ファイルを読み込んでから結果を書き出すまでの流れを日本語で解説します。各段階に関連する主要サブルーチンとデータ構造を併記し，数式を交えて説明します。

## 1. 初期化 (Initialization)

1. `main.f` から解析が開始され，`analys.f` が全体の制御を行います。
2. `redcml.f` により CML 形式の入力ファイルが読み込まれます。ここで座標 `xyz`，要素接続 `conn`，材料定数 `prop` などの配列が初期化されます。【F:redcml.f†L1-L10】
3. `initia.f` では解析パラメータと収束許容値が設定され，出力ファイルが開かれます。

## 2. 組み立て (Assembly)

1. 各要素ルーチン (`phexa8.f`，`pquad4.f` など) が要素剛性行列 `K_e` と内部力ベクトル `f_e` を計算します。要素剛性は次の式に基づきます。

\[
K_e = \int_\Omega B^T\,C\,B\, d\Omega
\]

2. `assemb.f` は上式で得られた `K_e` を大域剛性行列 `K` に組み込みます。【F:assemb.f†L1-L10】

## 3. 解法 (Solving)

1. `solut.f` で求解器が選択され，以下の大域平衡方程式を解きます。

\[
K\,\Delta u = F_{\text{ext}} - F_{\text{int}}
\]

2. Newton–Raphson 反復では残差 \(R\) が \(\|R\| < \varepsilon\) となるまで更新を続けます。線形方程式の解法には Skyline，PARDISO，CG などが利用可能です。

## 4. ポスト処理 (Post-processing)

1. `postpr.f` が要素応力を計算し，`stress.f` や `stress_dp.f` を通じて Drucker–Prager 収束条件を確認します。降伏関数は次式で表されます。

\[
\sqrt{J_2} + \alpha I_1 - k = 0
\]

2. `output.f` により変位 `RES_*` と応力 `STS_*` のファイルが書き出されます。

## 5. 一連の流れ

1. `redcml` で CML ファイルを読み込む
2. `constr` で境界条件を適用
3. `assemb` で大域剛性行列を構築
4. `solut` で平衡方程式を解く
5. `update` で変位および履歴変数を更新
6. `postpr` / `output` で結果ファイルを出力

以上が PLSTss による解析の基本的な手順です。
