# 出力ファイル一覧 - Output File List

## 📊 生成済みデータファイル

### CSV データ (`output/csv/`)
| ファイル名 | 内容 | 用途 |
|-----------|------|------|
| `1elem_f_test_current_iterations.csv` | 反復回数データ（100ステップ） | 累積反復回数グラフ用 |
| `1elem_f_test_current_iterations_summary.txt` | 反復統計サマリー | 全体傾向の把握 |

**データ概要：**
- 総ステップ数：100
- 総反復回数：101
- 平均反復回数/ステップ：1.01
- 複数反復ステップ：13のみ（2回反復）

### 収束データ (`output/convergence/`)
| ファイル名 | 内容 | 用途 |
|-----------|------|------|
| `1elem_f_test_current_convergence_all.csv` | ステップ13の収束データ | 2次収束確認用 |
| `1elem_f_test_current_convergence_all_analysis.txt` | 収束率解析レポート | 収束特性の確認 |

**収束特性（ステップ13）：**
- 初期残差：1.371e-01
- 最終残差：3.419e-15
- 削減率：2.494e-14（2反復で14桁の改善）
- 収束タイプ：2次収束（Newton-Raphson特性）

### ログファイル (`output/logs/`)
| ファイル名 | 内容 |
|-----------|------|
| `1elem_f_test_current.log` | 100ステップ解析の完全ログ |

## 🔍 データの見方

### 1. 累積反復回数（Excelグラフ用）
```
Step, Time_Normalized, Iterations, Cumulative_Iterations
1,    0.01,           1,          1
2,    0.02,           1,          2
...
13,   0.13,           2,          14  ← 唯一の複数反復
...
100,  1.00,           1,          101
```

### 2. 残差収束データ（対数プロット用）
```
Step, Iteration, Residual,     Log10_Residual
13,   1,        1.371e-01,    -0.8631
13,   2,        3.419e-15,    -14.4661
```
→ 約14桁の改善 = 強い2次収束

## 📝 次のステップ

1. **Excelでグラフ作成**
   - CSVファイルを開いて、上記のカラムを使用してグラフ作成

2. **他の入力ファイルでの解析**
   - 非関連流れ則（φ≠ψ）のケースでの比較
   - 異なる材料パラメータでの検証

3. **応力-ひずみ曲線**
   - RES_*.cmlファイルが生成されたら`extract2csv.py`で抽出

## 📂 スクリプト配置

```
scripts/
├── analysis/
│   ├── parse_iterations_to_csv.py      # 反復データ抽出
│   ├── parse_residual_convergence.py   # 収束データ抽出
│   └── extract2csv.py                  # 応力-ひずみ抽出
│
└── plotting/
    └── plot_cumulative_iterations.py   # Pythonプロット（オプション）
```