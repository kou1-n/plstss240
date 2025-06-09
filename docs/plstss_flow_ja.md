# PLSTss 解析フロー

PLSTss は塑性体解析を行う Fortran 製 FEM プログラムです。ここでは CML 入力の読み込みから結果書き出しまでの一連の処理を日本語で説明します。

## 主要ステージ

### 1. 初期化
- `main.f` で配列を確保した後、`flopen` で各種ファイルを開き、`chkary` により作業領域の大きさを確認します。その後 `analys` を呼び出します【F:main.f†L170-L219】。
- `analys.f` 冒頭では `redcml` を使用して CML ファイルを読み込みます【F:analys.f†L90-L101】。必要に応じて `rowupr` や `presky`、`inform` でインデックスを整備します【F:analys.f†L120-L151】。

### 2. 組み立て
- `initia` で要素の積分点情報などを準備し、ロードステップごとの Newton-Raphson 反復が開始されます【F:analys.f†L180-L205】。
- 反復内では `assemb.f` が呼ばれ、要素剛性行列がグローバル行列 `sk` に組み込まれます【F:analys.f†L229-L238】。要素計算には `quad4a.f` や `hexa8a.f` などが利用されます【F:assemb.f†L22-L56】【F:assemb.f†L80-L117】。

### 3. 解法
- `constr` で境界条件を強制変位に変換後、選択されたソルバー（`skylin`、`PARDISO`、`RCI CG`）で線形方程式を解きます【F:analys.f†L249-L322】。
- `update` が解を再配置し、現在の変位量を更新します【F:update.f†L23-L34】。

### 4. ポスト処理
- `postpr` で内部力や応力を算定し、履歴変数を更新します【F:analys.f†L349-L359】【F:postpr.f†L1-L30】。
- 収束後は `stored` で履歴を保存し、`output` により `RES_*.cml` や `STS_*.txt` などの結果ファイルが生成されます【F:analys.f†L437-L457】【F:output.f†L192-L251】。
- 解析が完了すると `pars99` でソルバーのメモリを解放し、最後に `main.f` が全ファイルを閉じて終了します【F:analys.f†L467-L473】【F:main.f†L221-L237】。

## データの流れ

1. 実行時に入力ベース名を指定すると `redcml.f` が `<name>.cml` を読み込み、節点座標や要素接続、材料定数を取り込みます。
2. `solut.f` で `/SOLUT/` ブロックを作成しておくと、計算時に読み取られる荷重増分 `finc` などの値が定義されます【F:solut.f†L1-L20】。
3. 計算終了後、`output.f` により `RES_<name>.cml`、`DIS_<name>.txt`、`STS_<name>.txt` などが実行ディレクトリに出力されます。
4. `scripts/` 以下の Python ツールを使うと、これらの結果を可視化できます。

PLSTss の詳細なルーチンや変数については [`docs/FORTRAN_FILES.md`](FORTRAN_FILES.md) を参照してください。
