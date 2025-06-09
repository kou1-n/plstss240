# 塑性変形可視化ガイド

plstss240での塑性変形解析結果を可視化するためのスクリプトとその使用方法を説明します。

## 概要

塑性変形解析では以下の項目を可視化できます：

1. **von Mises応力の時間変化**
2. **降伏応力との比較**
3. **累積塑性ひずみエネルギー**
4. **相当塑性ひずみ**
5. **応力-ひずみ関係**
6. **塑性変形の進展過程**

## 使用方法

### 1. 基本的な可視化

```bash
# WORKディレクトリで実行
python scripts/visualize_plastic_data.py
```

**出力ファイル：**
- `output/plastic_analysis.png` - 4つのグラフを含む総合図
- `output/plastic_data.csv` - 数値データ（Excel等で開ける）

### 2. 詳細解析

```bash
# より詳細な解析
python scripts/detailed_analysis.py
```

**出力ファイル：**
- `output/detailed_analysis.png` - 6つのグラフを含む詳細図
- `output/plastic_progression.csv` - 塑性変形進展データ

## 解析結果の見方

### 基本可視化の内容

1. **von Mises Stress vs Load Step**
   - 各要素のvon Mises応力の変化
   - 赤い点線：降伏応力（200 MPa）
   - 降伏開始ステップの特定

2. **Max/Avg von Mises Stress vs Load Step**
   - 全要素中の最大・平均von Mises応力
   - 構造全体の応力レベル把握

3. **Cumulative Plastic Energy vs Load Step**
   - 累積塑性ひずみエネルギーの変化
   - 塑性変形の蓄積量

4. **Equivalent Plastic Strain vs Load Step**
   - 相当塑性ひずみの変化
   - 塑性変形量の定量評価

### 詳細解析の内容

1. **von Mises Stress vs Displacement**
   - 要素5の荷重-変位関係
   - 降伏点の明確な特定

2. **Plastic Energy vs Displacement**
   - 塑性エネルギーの蓄積過程

3. **Equivalent Plastic Strain vs Displacement**
   - 塑性ひずみの増加率

4. **σzz vs εzz**
   - Z方向（主荷重方向）の応力-ひずみ関係
   - 材料の構成則確認

5. **Plastic Yield Detail (Steps 3-5)**
   - 降伏開始付近の詳細
   - 塑性化の瞬間を拡大表示

6. **Plastic Strain Distribution**
   - 最終ステップでの要素間塑性ひずみ分布
   - 最大塑性ひずみ要素の特定

## データファイルの活用

### CSV形式データ

生成されるCSVファイルは以下の列を含みます：

**plastic_data.csv:**
```
Step,Element,von_Mises_Stress,Total_Plastic_Energy,Equiv_Plastic_Strain
1,1,54.598000,0.008021,0.008956
1,2,52.880000,0.007168,0.008467
...
```

**plastic_progression.csv:**
```
Step,Plastic_Elements_Count,Max_von_Mises_Element,Max_von_Mises_Stress
1,8,1,54.598000
2,8,1,109.200000
...
```

### Excel等での二次解析

CSVファイルをExcelで開いて以下の解析が可能：

1. **ピボットテーブル**での要素別・ステップ別集計
2. **グラフ作成**でのカスタム可視化
3. **回帰分析**での材料パラメータ同定
4. **統計解析**での信頼性評価

## 解析結果の解釈

### cube250514.cmlの場合

- **材料：** 鋼材（E=200GPa, ν=0.3, σy=200MPa）
- **載荷：** Z方向圧縮（2.5mm/10ステップ）
- **要素：** 8要素六面体メッシュ

**典型的な結果：**

1. **ステップ1-3：** 弾性域（von Mises < 200 MPa）
2. **ステップ4：** 降伏開始（von Mises ≈ 200 MPa）
3. **ステップ5-10：** 塑性域（塑性ひずみ蓄積）

### 降伏判定基準

```
if von_Mises_stress >= 200.0:
    状態 = "塑性"
else:
    状態 = "弾性"
```

### 塑性ひずみの累積

相当塑性ひずみは以下の関係で推定：

```
ε_p_eq ≈ sqrt(2 * W_p / σ_y)
```

ここで：
- `W_p`：塑性ひずみエネルギー
- `σ_y`：降伏応力（200 MPa）

## カスタマイズ

### スクリプトの修正

1. **要素選択の変更**（`detailed_analysis.py`）：
```python
element_idx = 4  # 要素5を解析対象に
```

2. **降伏応力の変更**：
```python
yield_stress = 200.0  # MPa
```

3. **変位増分の変更**：
```python
displacement = step * 0.25e-3  # mm/step
```

### 材料パラメータの変更

異なる材料での解析時は、入力ファイル（`cube250514.cml`）の材料定数を変更し、
可視化スクリプト内の降伏応力も対応して変更してください。

## トラブルシューティング

### よくあるエラー

1. **ファイルが見つからない**
```
エラー: ファイルが見つかりません: output/RES_cube250514.cml
```
→ まずplstss240を実行して結果ファイルを生成

2. **matplotlibエラー**
```
ModuleNotFoundError: No module named 'matplotlib'
```
→ `pip install matplotlib numpy`でインストール

3. **データ形式エラー**
→ RES_*.cmlファイルが正常に生成されているか確認

### 対処方法

1. **解析の再実行**：
```bash
make clean
make
./plstss240 cube250514.cml
```

2. **出力ディレクトリの確認**：
```bash
ls output/
```

3. **ログファイルの確認**：
```bash
cat output/NOR_cube250514.txt
```

## 参考情報

- [plastic_monitoring.md](plastic_monitoring.md) - 塑性検出機能の詳細
- [README.md](../README.md) - plstss240の基本使用方法
- [CHANGELOG.md](../CHANGELOG.md) - 機能の変更履歴 