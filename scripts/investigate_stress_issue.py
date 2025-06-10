#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Stress Issue Investigation Tool
応力値の異常調査ツール
"""

import numpy as np
import re
from pathlib import Path

def investigate_stress_issue():
    """応力問題の詳細調査"""
    
    print("=== Stress Issue Investigation ===")
    
    # RESファイルから要素8の応力を直接抽出
    res_file = Path("output/01_fem_raw_results/RES_cube250514.cml")
    
    if not res_file.exists():
        print(f"エラー: {res_file} が見つかりません")
        return
    
    print(f"解析ファイル: {res_file.name}")
    print()
    
    # RESファイルの解析
    element8_data = extract_element8_from_res(res_file)
    
    if element8_data:
        analyze_element8_stress(element8_data)
    
    # STSファイルの確認
    sts_file = Path("output/01_fem_raw_results/STS_cube250514.txt")
    if sts_file.exists():
        analyze_sts_output(sts_file)

def extract_element8_from_res(res_file):
    """RESファイルから要素8のデータを抽出"""
    
    print("📄 RESファイルから要素8データを抽出中...")
    
    with open(res_file, 'r') as f:
        content = f.read()
    
    # 要素8のストレスセクションを探す
    element_sections = []
    lines = content.split('\n')
    
    i = 0
    current_step = None
    while i < len(lines):
        line = lines[i].strip()
        
        # ステップ情報の検出
        if 'STEP=' in line:
            current_step = extract_step_number(line)
        
        # 要素8のストレスデータを探す
        if line.startswith('8') and current_step is not None:
            # 要素8の行から6つの応力成分を抽出
            stress_data = parse_stress_line(line)
            if stress_data:
                element_sections.append({
                    'step': current_step,
                    'stress': stress_data
                })
        
        i += 1
    
    print(f"  ✓ {len(element_sections)}ステップのデータを抽出")
    return element_sections

def extract_step_number(line):
    """ステップ番号を抽出"""
    match = re.search(r'STEP=\s*(\d+)', line)
    return int(match.group(1)) if match else None

def parse_stress_line(line):
    """応力行を解析"""
    parts = line.split()
    if len(parts) >= 7:  # 要素番号 + 6つの応力成分
        try:
            return [float(parts[i]) for i in range(1, 7)]  # σxx, σyy, σzz, τyz, τzx, τxy
        except ValueError:
            return None
    return None

def analyze_element8_stress(element8_data):
    """要素8の応力を分析"""
    
    print("\n🔍 要素8応力分析:")
    print("Step  σxx[MPa]   σyy[MPa]   σzz[MPa]   τyz[MPa]   τzx[MPa]   τxy[MPa]   σvm[MPa]")
    print("-" * 85)
    
    for data in element8_data:
        step = data['step']
        stress = data['stress']
        
        # von Mises応力を計算
        sxx, syy, szz, tyz, tzx, txy = stress
        von_mises = calculate_von_mises(sxx, syy, szz, tyz, tzx, txy)
        
        print(f"{step:4d}  {sxx:8.2f}   {syy:8.2f}   {szz:8.2f}   {tyz:8.2f}   {tzx:8.2f}   {txy:8.2f}   {von_mises:8.2f}")
        
        # 異常値の検出
        if von_mises < 0:
            print(f"  ⚠️  Step {step}: von Mises応力が負の値!")
        
        if abs(szz) > 1000:  # 異常に大きな値
            print(f"  ⚠️  Step {step}: zz方向応力が異常な値!")

def calculate_von_mises(sxx, syy, szz, tyz, tzx, txy):
    """von Mises応力を計算"""
    
    # von Mises応力の計算式
    # σvm = sqrt(0.5 * ((σxx-σyy)² + (σyy-σzz)² + (σzz-σxx)² + 6*(τxy² + τyz² + τzx²)))
    
    s1 = (sxx - syy)**2
    s2 = (syy - szz)**2  
    s3 = (szz - sxx)**2
    s4 = 6 * (txy**2 + tyz**2 + tzx**2)
    
    von_mises = np.sqrt(0.5 * (s1 + s2 + s3 + s4))
    
    return von_mises

def analyze_sts_output(sts_file):
    """STSファイルの出力を分析"""
    
    print(f"\n📊 STSファイル分析:")
    print(f"ファイル: {sts_file.name}")
    
    with open(sts_file, 'r') as f:
        lines = f.readlines()
    
    # ヘッダー行を確認
    if len(lines) >= 2:
        header1 = lines[0].strip()
        header2 = lines[1].strip()
        
        print(f"ヘッダー1: {header1}")
        print(f"ヘッダー2: {header2}")
        
        # 要素番号を確認
        if "1-1" in header2:
            print("⚠️  STSファイルは要素1の値を出力しています!")
            print("   .cmlファイルの変更が反映されていない可能性があります")
        elif "8-1" in header2:
            print("✓ STSファイルは要素8の値を出力しています")
        else:
            print("❓ 出力要素番号が不明です")

def check_cml_output_settings():
    """CMLファイルの出力設定を確認"""
    
    print(f"\n⚙️  CMLファイルの出力設定確認:")
    
    cml_files = ["cube250514.cml", "cube250514_dp.cml"]
    
    for cml_file in cml_files:
        if Path(cml_file).exists():
            print(f"\n📋 {cml_file}:")
            check_single_cml_output(cml_file)

def check_single_cml_output(cml_file):
    """単一CMLファイルの出力設定確認"""
    
    with open(cml_file, 'r') as f:
        lines = f.readlines()
    
    in_pstrs = False
    in_pstrn = False
    
    for i, line in enumerate(lines):
        line = line.strip()
        
        if line == "/PSTRS/":
            in_pstrs = True
            print(f"  応力出力設定 (行{i+1}):")
            continue
        elif line == "/PSTRN/":
            in_pstrn = True
            print(f"  ひずみ出力設定 (行{i+1}):")
            continue
        elif line.startswith("/") and (in_pstrs or in_pstrn):
            in_pstrs = False
            in_pstrn = False
            continue
        
        if (in_pstrs or in_pstrn) and line and not line.isdigit():
            if any(char.isdigit() for char in line):
                print(f"    {line}")

if __name__ == "__main__":
    investigate_stress_issue()
    check_cml_output_settings() 