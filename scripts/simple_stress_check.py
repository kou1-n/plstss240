#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Simple Stress Check - von Mises問題詳細調査
"""

import numpy as np
import pandas as pd
from pathlib import Path

def main():
    """メイン関数"""
    
    print("=== 応力値詳細検証 ===")
    
    # STSファイルからデータ読み込み
    sts_file = Path("output/01_fem_raw_results/STS_cube250514.txt")
    
    if not sts_file.exists():
        print(f"エラー: {sts_file} が見つかりません")
        return
    
    print(f"解析ファイル: {sts_file}")
    
    # データ読み込み
    data = load_sts_data(sts_file)
    
    if data is not None:
        analyze_stress_data(data)

def load_sts_data(sts_file):
    """STSファイルからデータを読み込み"""
    
    try:
        # ヘッダーをスキップしてデータを読み込み
        with open(sts_file, 'r') as f:
            lines = f.readlines()
        
        # ヘッダーの確認
        print(f"ヘッダー1: {lines[0].strip()}")
        print(f"ヘッダー2: {lines[1].strip()}")
        
        # データ部分の抽出
        data_lines = [line.strip() for line in lines[2:] if line.strip()]
        
        if not data_lines:
            print("データが見つかりません")
            return None
        
        # データを解析
        steps = []
        strains = []  # xx, yy, zz
        stresses = [] # xx, yy, zz
        
        for line in data_lines:
            parts = line.split()
            if len(parts) >= 7:
                step = int(parts[0])
                strain_xx = float(parts[1])
                strain_yy = float(parts[2]) 
                strain_zz = float(parts[3])
                stress_xx = float(parts[4])
                stress_yy = float(parts[5])
                stress_zz = float(parts[6])
                
                steps.append(step)
                strains.append([strain_xx, strain_yy, strain_zz])
                stresses.append([stress_xx, stress_yy, stress_zz])
        
        return {
            'steps': np.array(steps),
            'strains': np.array(strains),
            'stresses': np.array(stresses)
        }
        
    except Exception as e:
        print(f"データ読み込みエラー: {e}")
        return None

def analyze_stress_data(data):
    """応力データの詳細分析"""
    
    steps = data['steps']
    strains = data['strains']  # [N, 3] - xx, yy, zz
    stresses = data['stresses'] # [N, 3] - xx, yy, zz
    
    print(f"\n📊 データ概要:")
    print(f"  ステップ数: {len(steps)}")
    print(f"  ステップ範囲: {steps[0]} → {steps[-1]}")
    
    print(f"\n📈 ひずみ範囲:")
    print(f"  εxx: {strains[:, 0].min():.6f} → {strains[:, 0].max():.6f}")
    print(f"  εyy: {strains[:, 1].min():.6f} → {strains[:, 1].max():.6f}")
    print(f"  εzz: {strains[:, 2].min():.6f} → {strains[:, 2].max():.6f}")
    
    print(f"\n📊 応力範囲:")
    print(f"  σxx: {stresses[:, 0].min():.2f} → {stresses[:, 0].max():.2f} MPa")
    print(f"  σyy: {stresses[:, 1].min():.2f} → {stresses[:, 1].max():.2f} MPa")
    print(f"  σzz: {stresses[:, 2].min():.2f} → {stresses[:, 2].max():.2f} MPa")
    
    # von Mises応力を計算（せん断応力成分なしの場合）
    print(f"\n🔢 von Mises応力計算:")
    
    von_mises_values = []
    
    for i in range(len(steps)):
        sxx, syy, szz = stresses[i]
        
        # せん断応力が0の場合のvon Mises応力
        von_mises = calculate_von_mises_principal(sxx, syy, szz)
        von_mises_values.append(von_mises)
        
        print(f"  Step {steps[i]:2d}: σxx={sxx:7.2f}, σyy={syy:7.2f}, σzz={szz:7.2f} → σvm={von_mises:7.2f} MPa")
        
        # 異常値チェック
        if von_mises < 0:
            print(f"    ⚠️  負のvon Mises応力!")
        
        if abs(szz) > 1000:
            print(f"    ⚠️  異常に大きなzz応力!")
    
    von_mises_values = np.array(von_mises_values)
    
    print(f"\n📊 von Mises応力統計:")
    print(f"  範囲: {von_mises_values.min():.2f} → {von_mises_values.max():.2f} MPa")
    print(f"  平均: {von_mises_values.mean():.2f} MPa")
    
    # 負の値があるかチェック
    negative_count = np.sum(von_mises_values < 0)
    if negative_count > 0:
        print(f"  ⚠️  負の値: {negative_count}個")
    else:
        print(f"  ✓ すべて正の値")
    
    # zz応力の挙動を詳細分析
    analyze_zz_stress_behavior(steps, stresses[:, 2])

def calculate_von_mises_principal(sxx, syy, szz):
    """主応力からvon Mises応力を計算（せん断応力=0の場合）"""
    
    # von Mises応力の計算式（せん断応力成分なし）
    # σvm = sqrt(0.5 * ((σxx-σyy)² + (σyy-σzz)² + (σzz-σxx)²))
    
    term1 = (sxx - syy)**2
    term2 = (syy - szz)**2
    term3 = (szz - sxx)**2
    
    von_mises = np.sqrt(0.5 * (term1 + term2 + term3))
    
    return von_mises

def analyze_zz_stress_behavior(steps, szz_stresses):
    """zz応力の挙動を詳細分析"""
    
    print(f"\n🔍 zz応力挙動分析:")
    
    # 変化率の計算
    dzz_values = np.diff(szz_stresses)
    
    print(f"  初期値: {szz_stresses[0]:.2f} MPa")
    print(f"  最終値: {szz_stresses[-1]:.2f} MPa")
    print(f"  総変化: {szz_stresses[-1] - szz_stresses[0]:.2f} MPa")
    
    # 変化パターンの確認
    increasing_count = np.sum(dzz_values > 0)
    decreasing_count = np.sum(dzz_values < 0)
    
    print(f"  増加ステップ: {increasing_count}")
    print(f"  減少ステップ: {decreasing_count}")
    
    # 異常な変化の検出
    large_changes = np.abs(dzz_values) > 50  # 50MPa以上の変化
    if np.any(large_changes):
        print(f"  ⚠️  大きな変化が {np.sum(large_changes)} ステップで発生")
        change_indices = np.where(large_changes)[0]
        for idx in change_indices:
            step_from = steps[idx]
            step_to = steps[idx+1]
            change = dzz_values[idx]
            print(f"    Step {step_from}→{step_to}: Δσzz = {change:.2f} MPa")
    else:
        print(f"  ✓ 変化は妥当な範囲内")

if __name__ == "__main__":
    main() 