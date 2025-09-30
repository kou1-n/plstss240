# PLSTss解析結果サマリー - Analysis Summary

## 🔍 重要な発見

### 1. **降伏応力の20%誤差の原因が判明**

入力ファイル`1elem_f.cml`の境界条件を詳細に解析した結果：

| 問題点 | 詳細 |
|--------|------|
| **載荷条件** | 純粋な単軸圧縮ではなく「部分拘束圧縮」 |
| **境界条件** | 各節点で異なる拘束（非対称） |
| **横方向拘束** | X方向に部分的な拘束あり |
| **応力状態** | 多軸応力状態（単軸ではない） |

**結論**: 理論降伏応力372.37 MPaに対し、実際446.84 MPa（20%高い）は、横方向拘束による閉じ込め効果が原因。

### 2. **改善提案**

#### 真の単軸圧縮試験のための境界条件：
```
底面：Z方向のみ拘束（横方向自由）
上面：Y方向に均一変位または圧力
側面：完全自由（拘束なし）
剛体運動防止：1節点のみX,Y拘束
```

## 📁 整理されたファイル構造

```
WORK/
├── src/               # Fortranソースコード
├── input_files/       # 入力ファイル（.cml）
├── output/
│   ├── csv/          # 全CSVデータ（統合済み）
│   ├── convergence/  # 収束解析結果
│   └── logs/         # 実行ログ
├── scripts/
│   ├── analysis/     # データ抽出・検証スクリプト
│   │   ├── parse_iterations_to_csv.py
│   │   ├── parse_residual_convergence.py
│   │   ├── extract2csv.py
│   │   ├── validate_drucker_prager.py
│   │   └── verify_loading_conditions.py
│   └── plotting/     # グラフ作成スクリプト
└── docs/             # ドキュメント

```

## 📊 利用可能なデータ

| ファイル | 内容 | 用途 |
|---------|------|------|
| `stress_strain_1elem_f.csv` | 応力-ひずみデータ（100ステップ） | S-Sカーブ作成 |
| `1elem_f_test_current_iterations.csv` | 反復回数データ | 累積反復プロット |
| `1elem_f_test_current_convergence_all.csv` | 残差収束データ | 2次収束確認 |

## 🔧 主要スクリプト

### データ抽出
```bash
# 応力-ひずみ抽出
python scripts/analysis/extract2csv.py output/RES_1elem_f.cml output/csv/stress_strain.csv

# 反復データ抽出
python scripts/analysis/parse_iterations_to_csv.py output/logs/logfile.log

# 収束データ抽出
python scripts/analysis/parse_residual_convergence.py output/logs/logfile.log
```

### 検証
```bash
# Drucker-Prager妥当性検証
python scripts/analysis/validate_drucker_prager.py output/csv/stress_strain.csv

# 載荷条件検証
python scripts/analysis/verify_loading_conditions.py
```

## 📈 Excelグラフ作成用CSVカラム

### 応力-ひずみ曲線
- X軸: `yystrain`（符号反転推奨）
- Y軸: `yystress`（符号反転推奨）

### 累積反復（Yamamoto図4タイプ）
- X軸: `Time_Normalized`
- Y軸: `Cumulative_Iterations`

### 残差収束（セミログ）
- X軸: `Iteration`
- Y軸: `Log10_Residual`

## ✅ クリーンアップ完了
- 古い解析結果を削除
- ネスト構造を簡素化
- CSVファイルを統合