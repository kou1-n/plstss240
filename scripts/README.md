# PLSTss Analysis Scripts - ファイル構成

## 📁 ディレクトリ構造

```
scripts/
├── README.md              # このファイル
│
├── analysis/              # 解析・データ抽出スクリプト
│   ├── parse_iterations_to_csv.py      # 反復回数データ抽出
│   ├── parse_residual_convergence.py   # 残差収束データ抽出
│   └── extract2csv.py                  # 応力-ひずみデータ抽出
│
└── plotting/              # グラフ作成スクリプト
    └── plot_cumulative_iterations.py   # 累積反復回数プロット

output/
├── csv/                   # CSV出力ファイル
│   ├── *_iterations.csv              # 反復回数データ
│   └── *_iterations_summary.txt      # 反復回数サマリー
│
├── convergence/           # 収束解析結果
│   ├── *_convergence_*.csv           # 収束データ
│   └── *_convergence_*_analysis.txt  # 収束解析レポート
│
├── plots/                 # グラフ出力（Excel/matplotlib）
│
└── logs/                  # 解析ログファイル
    └── *.log             # PLSTss実行ログ
```

## 📊 スクリプト使用方法

### 1. 反復回数データの抽出（Yamamoto et al. Figure 4のRMプロット用）

```bash
# ログファイルから反復回数をCSV化
python scripts/analysis/parse_iterations_to_csv.py output/logs/1elem_f_test_current.log

# 出力：
#   output/csv/1elem_f_test_current_iterations.csv
#   output/csv/1elem_f_test_current_iterations_summary.txt
```

**CSVカラム説明：**
- `Step`: ロードステップ番号
- `Time_Normalized`: 正規化時間 (0-1)
- `Iterations`: 各ステップの反復回数
- `Cumulative_Iterations`: 累積反復回数

### 2. 残差収束データの抽出（2次収束の確認用）

```bash
# 全ステップの収束データ
python scripts/analysis/parse_residual_convergence.py output/logs/1elem_f_test_current.log

# 特定ステップ（例：ステップ13）のみ
python scripts/analysis/parse_residual_convergence.py output/logs/1elem_f_test_current.log -s 13

# 出力：
#   output/convergence/1elem_f_test_current_convergence_*.csv
#   output/convergence/1elem_f_test_current_convergence_*_analysis.txt
```

**CSVカラム説明：**
- `Step`: ロードステップ番号
- `Iteration`: 反復番号
- `Rnorm`: 残差ノルム
- `Fnorm`: 力ノルム
- `Residual`: 正規化残差
- `Log10_Residual`: 残差の対数（セミログプロット用）

### 3. 応力-ひずみ曲線データの抽出

```bash
# RES_*.cmlファイルから応力-ひずみデータを抽出
python scripts/analysis/extract2csv.py RES_1elem_f.cml

# 出力ファイル名指定
python scripts/analysis/extract2csv.py RES_1elem_f.cml output/csv/stress_strain.csv

# 出力：
#   タイムスタンプ付きCSV or 指定ファイル名
```

**CSVカラム説明：**
- `step`: ロードステップ
- `element`: 要素番号
- `von`: von Mises応力（相当応力）
- `eps`: 相当塑性ひずみ
- `yystress`: y方向垂直応力
- `yystrain`: y方向垂直ひずみ

## 📈 Excelでのグラフ作成

### 累積反復回数プロット（Figure 4 RMタイプ）
1. `output/csv/*_iterations.csv`を開く
2. X軸：`Time_Normalized`列
3. Y軸：`Cumulative_Iterations`列
4. グラフタイプ：折れ線グラフ

### 残差収束プロット（2次収束確認）
1. `output/convergence/*_convergence_*.csv`を開く
2. X軸：`Iteration`列
3. Y軸：`Log10_Residual`列（セミログ用）または`Residual`列
4. グラフタイプ：散布図＋対数軸

### 応力-ひずみ曲線
1. `extract2csv.py`で生成したCSVを開く
2. X軸：`yystrain`列（またはひずみデータ）
3. Y軸：`yystress`列（またはvon Mises応力）
4. グラフタイプ：折れ線グラフ

## 🔧 Python環境要件

```bash
# 必要なパッケージ
pip install numpy matplotlib pathlib
```

## 📝 使用例

### 完全な解析フロー

```bash
# 1. PLSTss解析実行（リモート）
ssh nagasaku@mie "cd WORK && echo '1elem_f' | ~/bin/plstss2393_nagasaku > output/logs/1elem_f_new.log"

# 2. ログファイルをローカルにコピー
scp nagasaku@mie:WORK/output/logs/1elem_f_new.log output/logs/

# 3. データ抽出
python scripts/analysis/parse_iterations_to_csv.py output/logs/1elem_f_new.log
python scripts/analysis/parse_residual_convergence.py output/logs/1elem_f_new.log

# 4. RESファイルから応力-ひずみ抽出（RESファイルが生成されている場合）
scp nagasaku@mie:WORK/RES_1elem_f.cml ./
python scripts/analysis/extract2csv.py RES_1elem_f.cml output/csv/stress_strain_1elem_f.csv
```

## 📌 注意事項

- ログファイルのエンコーディング問題を避けるため、スクリプトは`utf-8`で読み込み、エラーを無視します
- CSVファイルはExcelで直接開けます
- 大規模な解析の場合、ログファイルが大きくなる可能性があります