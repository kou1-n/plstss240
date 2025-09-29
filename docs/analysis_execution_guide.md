# PLSTss解析実行ガイド

最終更新: 2024-09-30

## 解析実行方法

### 1. 単一ファイルの解析実行

#### 最も簡単な方法
```bash
./run_analysis ファイル名（拡張子なし）
```

**例:**
```bash
./run_analysis 1elem_f
```

このコマンドは以下を自動実行：
- `input_files/1elem_f.cml`を読み込み
- 解析を実行
- ターミナル出力を`output/logs/1elem_f_YYYYMMDD_HHMMSS.log`に保存
- 結果ファイルを`output/`に生成

#### 詳細なコマンド（直接実行）
```bash
echo "input_files/1elem_f" | ~/bin/plstss2393_nagasaku 2>&1 | tee output/logs/analysis.log
```

### 2. 複数ファイルの一括解析

#### すべての入力ファイルを解析
```bash
./scripts/batch_analysis.sh
```

#### 特定パターンのファイルのみ解析
```bash
./scripts/batch_analysis.sh "1elem_*"
```
例：`1elem_`で始まるすべてのファイルを解析

### 3. その他の実行スクリプト

| スクリプト | 説明 | 使用例 |
|-----------|------|--------|
| `./run_analysis` | 単一ファイルの簡単な実行（ルートディレクトリ） | `./run_analysis 1elem_f` |
| `scripts/batch_analysis.sh` | 複数ファイルの一括実行 | `./scripts/batch_analysis.sh` |
| `scripts/run_plstss_with_log.sh` | ログ付き実行（旧版） | `./scripts/run_plstss_with_log.sh input_files/1elem_f` |
| `scripts/run_plstss.py` | Python版実行スクリプト | `python scripts/run_plstss.py 1elem_f` |

## 出力ファイルの場所

### 解析結果ファイル
解析実行後、以下のファイルが生成されます：

| ファイル | 内容 | 場所 |
|---------|------|------|
| `LOG_*.log` | ターミナル出力の完全ログ（デバッグ情報含む） | `output/logs/` |
| `STS_*.txt` | 応力・ひずみデータ | `output/` または `output/results/latest/` |
| `DIS_*.txt` | 節点変位・反力データ | `output/` または `output/results/latest/` |
| `NOR_*.txt` | 収束ノルム・反復情報 | `output/` または `output/results/latest/` |
| `ENE_*.txt` | エネルギー履歴 | `output/` または `output/results/latest/` |
| `TMP_*.txt` | 温度データ | `output/` または `output/results/latest/` |
| `RES_*.cml` | 結果CMLファイル | `output/` または `output/results/latest/` |

### ログファイルの特徴
- **LOG_*.log**: 塑性修正の詳細情報（ftreg値、Newton-Raphson反復、consistent tangent診断など）を含む完全なターミナル出力

## ディレクトリ構造

```
WORK/
├── run_analysis              # 簡単実行スクリプト（単一ファイル用）
├── input_files/             # 入力ファイル（.cml）の置き場所
│   ├── 1elem_f.cml
│   ├── 1elem_phi10psi5.cml
│   └── ...
├── output/                  # 出力ディレクトリ
│   ├── logs/               # ターミナル出力ログ
│   │   └── batch_YYYYMMDD_HHMMSS/  # 一括実行時のログフォルダ
│   ├── results/
│   │   ├── latest/         # 最新の解析結果
│   │   └── archive/        # 過去の結果アーカイブ
│   └── [直接出力される結果ファイル]
└── scripts/                # 各種実行スクリプト
    ├── batch_analysis.sh   # 一括解析用
    ├── run_plstss_with_log.sh
    └── run_plstss.py
```

## 使用例

### 例1: 単一ファイルの解析
```bash
# 1elem_fを解析
./run_analysis 1elem_f

# 出力例：
# 解析開始: 1elem_f
# ログファイル: output/logs/1elem_f_20240930_143022.log
# =====================================
# [解析実行中のターミナル出力]
# =====================================
# 解析完了！
# ログ: output/logs/1elem_f_20240930_143022.log
# 結果ファイル:
# STS_1elem_f.txt
# DIS_1elem_f.txt
# NOR_1elem_f.txt
# ...
```

### 例2: パターン指定での一括解析
```bash
# phi10を含むファイルをすべて解析
./scripts/batch_analysis.sh "*phi10*"

# 出力例：
# ==========================================
# PLSTss Batch Analysis Tool
# ==========================================
# Found 2 input file(s) to analyze:
# 1elem_phi10psi5
# 1elem_phi10psi5_step11
#
# [1/2] Analyzing: 1elem_phi10psi5
# ✓ Success: Analysis completed
# [2/2] Analyzing: 1elem_phi10psi5_step11
# ✓ Success: Analysis completed
# ==========================================
# Batch Analysis Complete!
# Total files: 2
# Successful: 2
# Failed: 0
# Time elapsed: 3 seconds
```

## トラブルシューティング

### 入力ファイルが見つからない場合
```bash
./run_analysis
# 利用可能な入力ファイル一覧が表示される
```

### 解析が失敗する場合
1. ログファイルを確認：`output/logs/`内の該当ログファイル
2. エラーメッセージで「forrtl: 致命的なエラー (64): 入力変換エラー」が出る場合：
   - 入力ファイルのフォーマットを確認（特に科学記数法の表記）
   - 拡張子が`.cml`であることを確認

### SSHで実行する場合
```bash
ssh nagasaku@mie "cd WORK && ./run_analysis 1elem_f"
```

## 注意事項

1. **入力ファイルの場所**: 必ず`input_files/`フォルダに配置
2. **ファイル名**: 拡張子`.cml`は自動付与されるため、実行時は省略
3. **ログファイル**: タイムスタンプ付きで自動保存されるため、過去の実行履歴も残る
4. **一括実行時**: 大量のファイルを処理する場合は時間がかかることがある

## 関連ドキュメント

- `README.md`: プロジェクト全体の説明
- `CLAUDE.md`: プロジェクト固有の設定とガイドライン
- `docs/claude_code_instructions.md`: コード修正に関する技術詳細