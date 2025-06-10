#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Drucker-Prager降伏基準と累積塑性ひずみの統合可視化ツール
RESファイルを読み込み、Drucker-Prager基準値を計算して可視化
"""

import numpy as np
import matplotlib.pyplot as plt
import os
import sys
from pathlib import Path

# 日本語フォントの設定
plt.rcParams['font.family'] = ['DejaVu Sans', 'Arial Unicode MS', 'Hiragino Sans']
plt.rcParams['axes.unicode_minus'] = False


class DruckerPragerAnalyzer:
    def __init__(self, output_dir="../output"):
        self.output_dir = Path(output_dir)
        self.results = {}
        
    def read_res_file(self, filename):
        """RESファイルを読み込み、応力データを抽出"""
        res_path = self.output_dir / filename
        if not res_path.exists():
            print(f"Error: {res_path} が見つかりません")
            return None
            
        print(f"Reading {res_path}")
        
        stress_data = {}  # {step: {element: stress_tensor}}
        strain_data = {}  # {step: {element: plastic_strain}}
        
        with open(res_path, 'r') as f:
            lines = f.readlines()
            
        step_num = 0
        i = 0
        
        while i < len(lines):
            line = lines[i].strip()
            
            # /NODAL/セクションを検出（新しいステップ）
            if line == '/NODAL/':
                step_num += 1
                stress_data[step_num] = {}
                strain_data[step_num] = {}
                i += 1
                continue
            
            # /ELMTL/セクションで要素データを読み込み
            elif line == '/ELMTL/':
                i += 1
                if i < len(lines):
                    # ヘッダー行をスキップ
                    i += 2  # skip header lines
                    
                    # 要素データを読み込み
                    while i < len(lines) and lines[i].strip() and not lines[i].startswith('/'):
                        elem_line = lines[i].strip()
                        if elem_line and elem_line[0].isdigit():
                            parts = elem_line.split()
                            if len(parts) >= 7:
                                elem_num = int(parts[0])
                                
                                # 応力成分を読み込み（行1: SXX, SYY, SZZ, SXY, SYZ, SZX）
                                stress_tensor = {
                                    'SXX': float(parts[1]),
                                    'SYY': float(parts[2]), 
                                    'SZZ': float(parts[3]),
                                    'SXY': float(parts[4]),
                                    'SYZ': float(parts[5]),
                                    'SZX': float(parts[6])
                                }
                                stress_data[step_num][elem_num] = stress_tensor
                                
                                # 塑性ひずみ行（行2）
                                i += 1
                                if i < len(lines):
                                    strain_line = lines[i].strip()
                                    strain_parts = strain_line.split()
                                    if len(strain_parts) >= 6:
                                        # 有効塑性ひずみの計算
                                        ep_xx = float(strain_parts[0])
                                        ep_yy = float(strain_parts[1])
                                        ep_zz = float(strain_parts[2])
                                        ep_xy = float(strain_parts[3])
                                        ep_yz = float(strain_parts[4])
                                        ep_zx = float(strain_parts[5])
                                        
                                        plastic_strain = np.sqrt(2.0/3.0 * (ep_xx**2 + ep_yy**2 + ep_zz**2 + 
                                                                          2*(ep_xy**2 + ep_yz**2 + ep_zx**2)))
                                        strain_data[step_num][elem_num] = plastic_strain
                                
                                # von Mises応力行（行3）をスキップ
                                i += 1
                        
                        i += 1
                continue
            
            i += 1
                    
        return stress_data, strain_data
    
    def calculate_drucker_prager(self, stress_tensor, alpha=0.3):
        """
        Drucker-Prager降伏基準値を計算
        DP = sqrt(J2) + alpha * I1
        
        Parameters:
        - stress_tensor: 応力テンソル辞書 {SXX, SYY, SZZ, SXY, SYZ, SZX}
        - alpha: Drucker-Prager係数（デフォルト0.3）
        """
        if not stress_tensor:
            return 0.0
            
        # 応力成分の取得
        sxx = stress_tensor.get('SXX', 0.0)
        syy = stress_tensor.get('SYY', 0.0)
        szz = stress_tensor.get('SZZ', 0.0)
        sxy = stress_tensor.get('SXY', 0.0)
        syz = stress_tensor.get('SYZ', 0.0)
        szx = stress_tensor.get('SZX', 0.0)
        
        # 第1不変量 I1 = σxx + σyy + σzz
        I1 = sxx + syy + szz
        
        # 平均応力
        sm = I1 / 3.0
        
        # 偏差応力
        sxx_dev = sxx - sm
        syy_dev = syy - sm
        szz_dev = szz - sm
        
        # 第2偏差不変量 J2
        J2 = 0.5 * (sxx_dev**2 + syy_dev**2 + szz_dev**2) + sxy**2 + syz**2 + szx**2
        
        # Drucker-Prager基準値
        if J2 < 0:
            J2 = 0.0
        dp_value = np.sqrt(J2) + alpha * I1
        
        return dp_value
    
    def analyze_element(self, stress_data, strain_data, element_num=8):
        """指定要素について解析"""
        steps = []
        dp_values = []
        plastic_strains = []
        
        for step in sorted(stress_data.keys()):
            if (element_num in stress_data[step] and 
                element_num in strain_data[step]):
                
                stress_tensor = stress_data[step][element_num]
                plastic_strain = strain_data[step][element_num]
                
                dp_value = self.calculate_drucker_prager(stress_tensor)
                
                steps.append(step)
                dp_values.append(dp_value)
                plastic_strains.append(plastic_strain)
        
        return steps, dp_values, plastic_strains
    
    def plot_drucker_prager_analysis(self, filename, element_num=8):
        """Drucker-Prager基準と累積塑性ひずみの可視化"""
        # データ読み込み
        stress_data, strain_data = self.read_res_file(filename)
        if stress_data is None:
            return
            
        # 指定要素の解析
        steps, dp_values, plastic_strains = self.analyze_element(
            stress_data, strain_data, element_num)
        
        if not steps:
            print(f"Element {element_num} のデータが見つかりません")
            return
            
        # プロット作成
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(14, 10))
        fig.suptitle(f'Drucker-Prager基準解析 - Element {element_num}', 
                     fontsize=16, fontweight='bold')
        
        # 1. Drucker-Prager基準値 vs 累積塑性ひずみ（メインプロット）
        ax1.plot(plastic_strains, dp_values, 'bo-', linewidth=2, markersize=6)
        ax1.set_xlabel('累積塑性ひずみ')
        ax1.set_ylabel('Drucker-Prager基準値 [MPa]')
        ax1.set_title('Drucker-Prager基準 vs 累積塑性ひずみ')
        ax1.grid(True, alpha=0.3)
        
        # ステップ番号をアノテーション
        for i, (x, y, step) in enumerate(zip(plastic_strains, dp_values, steps)):
            if i % 2 == 0:  # 見やすくするため、2ステップおきに表示
                ax1.annotate(f'S{step}', (x, y), xytext=(5, 5), 
                           textcoords='offset points', fontsize=8)
        
        # 2. ステップ vs Drucker-Prager基準値
        ax2.plot(steps, dp_values, 'ro-', linewidth=2, markersize=6)
        ax2.set_xlabel('ステップ')
        ax2.set_ylabel('Drucker-Prager基準値 [MPa]')
        ax2.set_title('ステップ vs Drucker-Prager基準値')
        ax2.grid(True, alpha=0.3)
        
        # 3. ステップ vs 累積塑性ひずみ
        ax3.plot(steps, plastic_strains, 'go-', linewidth=2, markersize=6)
        ax3.set_xlabel('ステップ')
        ax3.set_ylabel('累積塑性ひずみ')
        ax3.set_title('ステップ vs 累積塑性ひずみ')
        ax3.grid(True, alpha=0.3)
        
        # 4. 塑性ひずみ増分
        strain_increments = [0]
        for i in range(1, len(plastic_strains)):
            increment = plastic_strains[i] - plastic_strains[i-1]
            strain_increments.append(increment)
        
        ax4.bar(steps, strain_increments, alpha=0.7, color='purple')
        ax4.set_xlabel('ステップ')
        ax4.set_ylabel('塑性ひずみ増分')
        ax4.set_title('ステップごとの塑性ひずみ増分')
        ax4.grid(True, alpha=0.3)
        
        plt.tight_layout()
        
        # 画像保存
        output_path = self.output_dir / f"drucker_prager_analysis_{filename.replace('.cml', '')}.png"
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        print(f"グラフを保存しました: {output_path}")
        
        # データ保存
        self.save_analysis_data(steps, dp_values, plastic_strains, strain_increments, 
                               filename, element_num)
        
        plt.show()
        
        return steps, dp_values, plastic_strains
    
    def save_analysis_data(self, steps, dp_values, plastic_strains, strain_increments,
                          filename, element_num):
        """解析データをCSVファイルに保存"""
        data_path = self.output_dir / f"drucker_prager_data_{filename.replace('.cml', '')}.csv"
        
        with open(data_path, 'w') as f:
            f.write(f"# Drucker-Prager基準解析データ\n")
            f.write(f"# ファイル: {filename}\n")
            f.write(f"# 要素: {element_num}\n")
            f.write("Step,Drucker_Prager_Value[MPa],Cumulative_Plastic_Strain,Strain_Increment\n")
            
            for i, (step, dp, strain, inc) in enumerate(zip(steps, dp_values, plastic_strains, strain_increments)):
                f.write(f"{step},{dp:.6f},{strain:.8f},{inc:.8f}\n")
        
        print(f"データを保存しました: {data_path}")


def main():
    """メイン実行関数"""
    analyzer = DruckerPragerAnalyzer()
    
    # 使用可能なRESファイルを確認
    res_files = list(analyzer.output_dir.glob("RES_*.cml"))
    
    if not res_files:
        print("RESファイルが見つかりません")
        return
    
    print("利用可能なRESファイル:")
    for i, file in enumerate(res_files):
        print(f"{i+1}. {file.name}")
    
    # デフォルトでRES_cube250514.cmlを解析
    default_file = "RES_cube250514.cml"
    if (analyzer.output_dir / default_file).exists():
        print(f"\n{default_file} を解析します...")
        analyzer.plot_drucker_prager_analysis(default_file, element_num=8)
    else:
        print(f"{default_file} が見つかりません")


if __name__ == "__main__":
    main() 