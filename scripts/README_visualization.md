# 可視化ツール使用ガイド

## 概要
plstss240の有限要素解析結果を可視化するためのPythonツール群を整理しました。

## 主要ツール

### 1. Drucker-Prager基準解析ツール (推奨)
- **ファイル**: `drucker_prager_analysis.py`
- **機能**: RESファイルからDrucker-Prager降伏基準と累積塑性ひずみの関係を解析・可視化
- **出力**: 
  - グラフ: `drucker_prager_analysis_RES_cube250514.png`
  - データ: `drucker_prager_data_RES_cube250514.csv`

**使用方法**:
```bash
cd scripts
python drucker_prager_analysis.py
```

**プロット内容**:
1. Drucker-Prager基準値 vs 累積塑性ひずみ（メイン）
2. ステップ vs Drucker-Prager基準値
3. ステップ vs 累積塑性ひずみ
4. ステップごとの塑性ひずみ増分

### 2. 統合可視化ツール
- **ファイル**: `main_visualizer.py`
- **機能**: 従来の複数プロットツールの統合版
- **出力**: 複数のグラフファイル

### 3. データローダー
- **ファイル**: `data_loader.py`
- **機能**: 解析結果ファイルの読み込み関数群

## アーカイブされたツール

以下のツールは`archive/`ディレクトリに移動されました：

- `plot_*.py` - 個別プロットツール群
- `visualize_plastic_data.py` - 旧塑性データ可視化ツール
- `detailed_analysis.py` - 旧詳細解析ツール
- `test_paths.py` - パステストツール
- `check_stress_dp_signature.py` - チェックツール

## Drucker-Prager基準について

Drucker-Prager降伏基準は以下の式で定義されます：
```
DP = √J₂ + α·I₁
```

ここで：
- J₂: 偏差応力テンソルの第2不変量
- I₁: 応力テンソルの第1不変量（静水圧応力×3）
- α: Drucker-Prager係数（デフォルト: 0.3）

## 解析結果サマリー（cube250514）

- **要素**: Element 8（中央要素）
- **ステップ数**: 10ステップ
- **降伏開始**: ステップ4で塑性変形開始
- **最終DP値**: 249.06 MPa
- **最終累積塑性ひずみ**: 0.00213

## 使用上の注意

1. 日本語フォントの警告が表示されますが、グラフ生成には影響ありません
2. 要素番号は1-8の範囲で指定可能
3. 出力ファイルは`../output/`ディレクトリに保存されます

## ファイル構成

```
scripts/
├── drucker_prager_analysis.py  # メイン解析ツール
├── main_visualizer.py          # 統合可視化ツール
├── data_loader.py              # データ読み込み
├── archive/                    # 旧ツール群
└── README_visualization.md     # このファイル
``` 