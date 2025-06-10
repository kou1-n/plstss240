# plstss240 可視化ツールファミリー

plstss240（Fortran有限要素解析プログラム）の結果を可視化するための統合ツールセットです。

## 🚀 クイックスタート

```bash
# 簡単起動（推奨）
python visualize.py

# 2軸プロットのみ
python visualize.py dual

# 包括的解析のみ
python visualize.py comprehensive

# 硬化曲線のみ
python visualize.py hardening
```

## 📁 ファイル構成

### メインツール
- **`plstss_visualizer.py`** - 統合可視化ツール（すべての機能を含む）
- **`visualize.py`** - 簡単起動スクリプト（エントリーポイント）
- **`hardening_simple.py`** - 硬化曲線専用ツール
- **`data_loader.py`** - データ読み込み共通ライブラリ

### アーカイブ
- **`archive/`** - 古いツールや開発用スクリプト

## 🎯 主な機能

### 1. 統合解析 (`plstss_visualizer.py`)
```bash
python plstss_visualizer.py --mode all --element 8
```

**出力:**
- 包括的解析図（4つのサブプロット）
- 2軸プロット（DP基準 vs 硬化応力）
- 理論硬化曲線
- CSVデータファイル

### 2. 2軸プロット解析
```bash
python plstss_visualizer.py --mode dual --element 8
```

**特徴:**
- Drucker-Prager降伏基準値（左軸、赤）
- 硬化後降伏応力（右軸、青）
- 初期降伏応力の参照線
- ステップ番号注釈

### 3. 包括的解析
```bash
python plstss_visualizer.py --mode comprehensive --element 8
```

**4つのサブプロット:**
1. DP基準 vs 塑性ひずみ
2. 硬化曲線
3. ステップ進行
4. 応力差分解析（弾性/塑性領域）

### 4. 硬化曲線解析
```bash
python plstss_visualizer.py --mode hardening
```

**出力:**
- 理論硬化曲線
- 初期降伏応力・漸近応力の参照線

## 🔧 コマンドオプション

### 基本オプション
```bash
--file, -f      RESファイル名 (デフォルト: RES_cube250514.cml)
--element, -e   要素番号 (デフォルト: 8)
--mode, -m      解析モード [all|dual|comprehensive|hardening]
--output, -o    出力ディレクトリ (デフォルト: ../output)
```

### 使用例
```bash
# 要素5の全解析
python visualize.py --element 5

# 要素10の2軸プロット
python visualize.py dual --element 10

# 異なるファイルの解析
python visualize.py --file RES_other.cml --element 1

# 出力先指定
python visualize.py --output ./results
```

## 📊 出力ファイル

### 画像ファイル
- `element{N}_dual_axis.png` - 2軸プロット
- `element{N}_comprehensive_analysis.png` - 包括的解析
- `theoretical_hardening_curve.png` - 理論硬化曲線
- `hardening_analysis_simple.png` - 硬化曲線比較

### データファイル
- `element{N}_analysis_RES_{filename}.csv` - 解析結果CSV

## 🧮 材料パラメータ

統合ツールは以下の材料パラメータを使用します：

```python
材料パラメータ = {
    'E': 200000.0,      # ヤング率 [MPa]
    'nu': 0.3,          # ポアソン比
    'sigma_y': 200.0,   # 初期降伏応力 [MPa]
    'hk': 200.0,        # 線形硬化係数 [MPa]
    'hpa': 400.0,       # 漸近応力 [MPa]
    'hpb': 10.0,        # 硬化指数
    'alpha_dp': 0.3     # Drucker-Prager係数
}
```

## 📈 解析結果の解釈

### 2軸プロット
- **赤線**: Drucker-Prager降伏基準値の推移
- **青線**: 硬化後降伏応力の推移
- **黒破線**: 初期降伏応力
- **注釈**: ステップ番号

### 包括的解析
1. **DP基準 vs 塑性ひずみ**: 降伏基準の進化
2. **硬化曲線**: 降伏応力の硬化
3. **ステップ進行**: 荷重ステップでの変化
4. **応力差分**: 弾性/塑性領域の可視化

### 数値出力
- **κ**: 等価塑性ひずみ
- **DP[MPa]**: Drucker-Prager基準値
- **σy[MPa]**: 硬化後降伏応力
- **DP-σy[MPa]**: 応力差分（負=弾性、正=塑性）

## 🔄 開発履歴

### Version 1.0 (現在)
- ✅ 統合可視化ツール (`plstss_visualizer.py`)
- ✅ 簡単起動スクリプト (`visualize.py`)
- ✅ 2軸プロット機能
- ✅ 包括的解析機能
- ✅ 硬化曲線解析
- ✅ CSV出力
- ✅ ファイル整理・アーカイブ

### アーカイブ済み
- 個別ツール（`drucker_prager_analysis.py`等）
- GUI版ツール
- 開発用スクリプト

## 🆘 トラブルシューティング

### よくある問題
1. **RESファイルが見つからない**
   ```bash
   # outputディレクトリに RES_cube250514.cml があることを確認
   ls ../output/RES_*.cml
   ```

2. **要素番号が存在しない**
   ```bash
   # 利用可能な要素番号を確認
   python visualize.py --element 1
   ```

3. **フォント警告**
   - 警告は無視して構いません（動作に影響なし）
   - 日本語フォントが利用できない環境での正常な動作

### デバッグモード
```bash
# 詳細なエラー情報を表示
python plstss_visualizer.py --mode dual --element 8 2>&1 | more
```

## 📞 サポート

問題や質問がある場合は、以下を確認してください：
1. RESファイルの存在確認
2. 要素番号の有効性
3. 出力ディレクトリの書き込み権限

---

**Author**: AI Assistant  
**Last Updated**: 2024年12月  
**Version**: 1.0 