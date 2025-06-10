#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Output Directory Organization Tool
outputディレクトリを分類・整理するツール
"""

import os
import shutil
from pathlib import Path
from datetime import datetime

def organize_output_directory():
    """outputディレクトリを整理する"""
    
    output_dir = Path("../output")
    
    if not output_dir.exists():
        print("outputディレクトリが見つかりません")
        return
    
    print("=== Output Directory Organization ===")
    print(f"整理対象: {output_dir.absolute()}")
    
    # 整理カテゴリの定義
    categories = {
        "fem_results": {
            "description": "FEM解析生ファイル",
            "patterns": ["RES_*.cml", "STS_*.txt", "DIS_*.txt", "NOR_*.txt", "ENE_*.txt", "TMP_*.txt"],
            "folder": "01_fem_raw_results"
        },
        "visualization": {
            "description": "可視化グラフ",
            "patterns": ["*.png"],
            "folder": "02_visualizations"
        },
        "data_analysis": {
            "description": "解析データ（CSV）",
            "patterns": ["*.csv"],
            "folder": "03_analysis_data"
        }
    }
    
    # サブカテゴリ（可視化内の詳細分類）
    viz_subcategories = {
        "hardening": {
            "description": "硬化関数分析",
            "patterns": ["hardening_*", "*hardening*"],
            "folder": "hardening_analysis"
        },
        "comprehensive": {
            "description": "包括分析",
            "patterns": ["*comprehensive*", "detailed_*", "deformation_increase_*"],
            "folder": "comprehensive_analysis"
        },
        "element_specific": {
            "description": "要素別分析",
            "patterns": ["element8_*", "drucker_prager_*"],
            "folder": "element_analysis"
        },
        "general": {
            "description": "一般分析",
            "patterns": ["plastic_*"],
            "folder": "general_analysis"
        }
    }
    
    # フォルダ作成
    create_directory_structure(output_dir, categories, viz_subcategories)
    
    # ファイル移動
    organize_files(output_dir, categories, viz_subcategories)
    
    # インデックス作成
    create_index_file(output_dir, categories, viz_subcategories)
    
    print("\n✅ 整理完了！")

def create_directory_structure(output_dir, categories, viz_subcategories):
    """ディレクトリ構造を作成"""
    
    print("\n📁 ディレクトリ構造作成中...")
    
    for cat_key, cat_info in categories.items():
        folder_path = output_dir / cat_info["folder"]
        folder_path.mkdir(exist_ok=True)
        print(f"  ✓ {cat_info['folder']} - {cat_info['description']}")
        
        # 可視化の場合はサブカテゴリも作成
        if cat_key == "visualization":
            for sub_key, sub_info in viz_subcategories.items():
                sub_folder_path = folder_path / sub_info["folder"]
                sub_folder_path.mkdir(exist_ok=True)
                print(f"    ├─ {sub_info['folder']} - {sub_info['description']}")

def organize_files(output_dir, categories, viz_subcategories):
    """ファイルを整理"""
    
    print("\n📋 ファイル移動中...")
    
    moved_files = []
    
    # 全ファイルリスト取得
    all_files = [f for f in output_dir.iterdir() if f.is_file()]
    
    for file_path in all_files:
        moved = False
        
        # メインカテゴリで分類
        for cat_key, cat_info in categories.items():
            if matches_patterns(file_path.name, cat_info["patterns"]):
                
                if cat_key == "visualization":
                    # 可視化ファイルはサブカテゴリで更に分類
                    for sub_key, sub_info in viz_subcategories.items():
                        if matches_patterns(file_path.name, sub_info["patterns"]):
                            dest_dir = output_dir / cat_info["folder"] / sub_info["folder"]
                            dest_path = dest_dir / file_path.name
                            shutil.move(str(file_path), str(dest_path))
                            print(f"  📊 {file_path.name} → {cat_info['folder']}/{sub_info['folder']}")
                            moved_files.append({
                                'file': file_path.name,
                                'category': cat_info['description'],
                                'subcategory': sub_info['description'],
                                'path': f"{cat_info['folder']}/{sub_info['folder']}"
                            })
                            moved = True
                            break
                    
                    # サブカテゴリにマッチしない場合は一般フォルダに
                    if not moved:
                        dest_dir = output_dir / cat_info["folder"] / "other"
                        dest_dir.mkdir(exist_ok=True)
                        dest_path = dest_dir / file_path.name
                        shutil.move(str(file_path), str(dest_path))
                        print(f"  📊 {file_path.name} → {cat_info['folder']}/other")
                        moved_files.append({
                            'file': file_path.name,
                            'category': cat_info['description'],
                            'subcategory': 'その他',
                            'path': f"{cat_info['folder']}/other"
                        })
                        moved = True
                else:
                    # 通常のカテゴリ移動
                    dest_dir = output_dir / cat_info["folder"]
                    dest_path = dest_dir / file_path.name
                    shutil.move(str(file_path), str(dest_path))
                    print(f"  📄 {file_path.name} → {cat_info['folder']}")
                    moved_files.append({
                        'file': file_path.name,
                        'category': cat_info['description'],
                        'subcategory': '',
                        'path': cat_info['folder']
                    })
                    moved = True
                break
        
        if not moved:
            print(f"  ❓ {file_path.name} - 分類不明（そのまま）")
    
    return moved_files

def matches_patterns(filename, patterns):
    """ファイル名がパターンにマッチするかチェック"""
    
    import fnmatch
    
    for pattern in patterns:
        if fnmatch.fnmatch(filename.lower(), pattern.lower()):
            return True
    return False

def create_index_file(output_dir, categories, viz_subcategories):
    """インデックスファイルを作成"""
    
    print("\n📑 インデックスファイル作成中...")
    
    index_content = f"""# Output Directory Index
Generated on: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}

## Directory Structure

### 📂 01_fem_raw_results - FEM解析生ファイル
plstss240の直接出力ファイル
- RES_*.cml - 解析結果メインファイル
- STS_*.txt - 応力結果
- DIS_*.txt - 変位結果  
- NOR_*.txt - 収束履歴
- ENE_*.txt - エネルギー履歴
- TMP_*.txt - 温度履歴

### 📊 02_visualizations - 可視化グラフ
#### hardening_analysis - 硬化関数分析
- 硬化曲線の解析結果
- 線形性評価
- 理論値との比較

#### comprehensive_analysis - 包括分析  
- 4分割総合分析
- 変形量増加分析
- 詳細解析結果

#### element_analysis - 要素別分析
- Element 8の個別分析
- Drucker-Prager基準分析
- デュアル軸プロット

#### general_analysis - 一般分析
- 塑性進展分析
- その他の解析結果

### 📈 03_analysis_data - 解析データ（CSV）
- element8_analysis_*.csv - 要素8解析データ
- drucker_prager_data_*.csv - DP基準データ
- hardening_*.csv - 硬化関数データ
- plastic_*.csv - 塑性進展データ

## Usage Tips

### 最新の結果を確認したい場合
1. `01_fem_raw_results/` - 最新のFEM解析結果
2. `02_visualizations/comprehensive_analysis/` - 総合的な可視化

### 硬化関数の詳細を調べたい場合
1. `02_visualizations/hardening_analysis/` - グラフ
2. `03_analysis_data/hardening_*.csv` - 数値データ

### 特定要素の解析結果
1. `02_visualizations/element_analysis/` - 要素8のグラフ
2. `03_analysis_data/element8_*.csv` - 要素8のデータ

---
Generated by plstss240 Output Organization Tool
"""
    
    index_path = output_dir / "README.md"
    with open(index_path, 'w', encoding='utf-8') as f:
        f.write(index_content)
    
    print(f"  ✓ インデックスファイル作成: {index_path.name}")

def show_organization_summary(output_dir):
    """整理結果のサマリー表示"""
    
    print("\n📋 整理結果サマリー:")
    
    total_files = 0
    for root, dirs, files in os.walk(output_dir):
        if files:
            rel_path = Path(root).relative_to(output_dir)
            print(f"  📁 {rel_path if str(rel_path) != '.' else 'ルート'}: {len(files)}ファイル")
            total_files += len(files)
    
    print(f"\n合計: {total_files}ファイル")

if __name__ == "__main__":
    organize_output_directory()
    show_organization_summary(Path("../output")) 