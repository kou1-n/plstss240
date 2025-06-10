#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
plstss240 統合可視化ツール
===============================

機能:
1. Drucker-Prager降伏基準解析
2. 硬化曲線解析  
3. 2軸プロット（DP基準 vs 硬化応力）
4. 応力-ひずみ関係の可視化
5. ステップ別解析

Author: AI Assistant
Version: 1.0
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')  # 非対話モード
import matplotlib.pyplot as plt
from pathlib import Path
import os
import sys
import argparse
from typing import Dict, List, Tuple, Optional
import warnings
warnings.filterwarnings('ignore')

# フォント設定
plt.rcParams['font.family'] = ['DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['figure.figsize'] = [12, 8]
plt.rcParams['font.size'] = 10


class PlstssVisualizer:
    """plstss240の統合可視化ツール"""
    
    def __init__(self, output_dir="../output"):
        self.output_dir = Path(output_dir)
        self.results = {}
        self.material_params = self._load_material_params()
        
    def _load_material_params(self):
        """材料パラメータをロード（デフォルト値）"""
        return {
            'E': 200000.0,      # Young's modulus [MPa]
            'nu': 0.3,          # Poisson's ratio
            'sigma_y': 200.0,   # Initial yield stress [MPa]
            'hk': 200.0,        # Linear hardening [MPa]
            'hpa': 400.0,       # Asymptotic stress [MPa]
            'hpb': 10.0,        # Hardening exponent
            'alpha_dp': 0.3     # Drucker-Prager parameter
        }
    
    def read_res_file(self, filename: str) -> Tuple[Dict, Dict]:
        """RESファイルを読み込み、応力・ひずみデータを抽出"""
        res_path = self.output_dir / filename
        if not res_path.exists():
            raise FileNotFoundError(f"RES file not found: {res_path}")
            
        print(f"Reading {res_path.name}...")
        
        stress_data = {}  # {step: {element: stress_tensor}}
        strain_data = {}  # {step: {element: plastic_strain}}
        
        with open(res_path, 'r') as f:
            lines = f.readlines()
            
        step_num = 0
        i = 0
        
        while i < len(lines):
            line = lines[i].strip()
            
            if line == '/NODAL/':
                step_num += 1
                stress_data[step_num] = {}
                strain_data[step_num] = {}
                i += 1
                continue
            
            elif line == '/ELMTL/':
                i += 1
                if i < len(lines):
                    i += 2  # skip header lines
                    
                    while i < len(lines) and lines[i].strip() and not lines[i].startswith('/'):
                        elem_line = lines[i].strip()
                        if elem_line and elem_line[0].isdigit():
                            parts = elem_line.split()
                            if len(parts) >= 7:
                                elem_num = int(parts[0])
                                
                                # 応力成分
                                stress_tensor = {
                                    'SXX': float(parts[1]), 'SYY': float(parts[2]), 
                                    'SZZ': float(parts[3]), 'SXY': float(parts[4]),
                                    'SYZ': float(parts[5]), 'SZX': float(parts[6])
                                }
                                stress_data[step_num][elem_num] = stress_tensor
                                
                                # 塑性ひずみ
                                i += 1
                                if i < len(lines):
                                    strain_line = lines[i].strip()
                                    strain_parts = strain_line.split()
                                    if len(strain_parts) >= 6:
                                        ep_xx = float(strain_parts[0])
                                        ep_yy = float(strain_parts[1])
                                        ep_zz = float(strain_parts[2])
                                        ep_xy = float(strain_parts[3])
                                        ep_yz = float(strain_parts[4])
                                        ep_zx = float(strain_parts[5])
                                        
                                        plastic_strain = np.sqrt(2.0/3.0 * (
                                            ep_xx**2 + ep_yy**2 + ep_zz**2 + 
                                            2*(ep_xy**2 + ep_yz**2 + ep_zx**2)))
                                        strain_data[step_num][elem_num] = plastic_strain
                                
                                i += 1  # skip von Mises stress line
                        
                        i += 1
                continue
            
            i += 1
                    
        return stress_data, strain_data
    
    def calculate_drucker_prager(self, stress_tensor: Dict, alpha: float = None) -> float:
        """Drucker-Prager降伏基準値を計算"""
        if alpha is None:
            alpha = self.material_params['alpha_dp']
            
        if not stress_tensor:
            return 0.0
            
        sxx = stress_tensor.get('SXX', 0.0)
        syy = stress_tensor.get('SYY', 0.0)
        szz = stress_tensor.get('SZZ', 0.0)
        sxy = stress_tensor.get('SXY', 0.0)
        syz = stress_tensor.get('SYZ', 0.0)
        szx = stress_tensor.get('SZX', 0.0)
        
        # 第1不変量と偏差応力
        I1 = sxx + syy + szz
        sm = I1 / 3.0
        sxx_dev = sxx - sm
        syy_dev = syy - sm
        szz_dev = szz - sm
        
        # 第2偏差不変量
        J2 = 0.5 * (sxx_dev**2 + syy_dev**2 + szz_dev**2) + sxy**2 + syz**2 + szx**2
        if J2 < 0:
            J2 = 0.0
        
        return np.sqrt(J2) + alpha * I1
    
    def calculate_hardened_yield_stress(self, kappa: float) -> float:
        """硬化後の降伏応力を計算"""
        params = self.material_params
        sigma_y = params['sigma_y']
        hk = params['hk']
        hpa = params['hpa']
        hpb = params['hpb']
        
        return sigma_y + hk * kappa + (hpa - sigma_y) * (1 - np.exp(-hpb * kappa))
    
    def analyze_element(self, stress_data: Dict, strain_data: Dict, element_num: int = 8):
        """指定要素の解析"""
        steps = []
        dp_values = []
        plastic_strains = []
        yield_stresses = []
        
        for step in sorted(stress_data.keys()):
            if (element_num in stress_data[step] and 
                element_num in strain_data[step]):
                
                stress_tensor = stress_data[step][element_num]
                plastic_strain = strain_data[step][element_num]
                
                dp_value = self.calculate_drucker_prager(stress_tensor)
                yield_stress = self.calculate_hardened_yield_stress(plastic_strain)
                
                steps.append(step)
                dp_values.append(dp_value)
                plastic_strains.append(plastic_strain)
                yield_stresses.append(yield_stress)
        
        return {
            'steps': steps,
            'dp_values': dp_values,
            'plastic_strains': plastic_strains,
            'yield_stresses': yield_stresses
        }
    
    def plot_drucker_prager_main(self, filename: str, element_num: int = 8):
        """メインのDrucker-Prager解析プロット"""
        stress_data, strain_data = self.read_res_file(filename)
        results = self.analyze_element(stress_data, strain_data, element_num)
        
        if not results['steps']:
            print(f"No data found for element {element_num}")
            return
        
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(16, 12))
        fig.suptitle(f'Element {element_num}: Drucker-Prager Analysis', 
                     fontsize=16, fontweight='bold')
        
        steps = results['steps']
        dp_values = results['dp_values']
        plastic_strains = results['plastic_strains']
        yield_stresses = results['yield_stresses']
        
        # 1. DP基準 vs 塑性ひずみ
        ax1.plot(plastic_strains, dp_values, 'bo-', linewidth=2, markersize=6)
        ax1.set_xlabel('Equivalent Plastic Strain κ')
        ax1.set_ylabel('Drucker-Prager Criterion [MPa]')
        ax1.set_title('DP Criterion vs Plastic Strain')
        ax1.grid(True, alpha=0.3)
        
        # 2. 硬化曲線
        ax2.plot(plastic_strains, yield_stresses, 'ro-', linewidth=2, markersize=6)
        ax2.axhline(y=self.material_params['sigma_y'], color='black', 
                   linestyle='--', alpha=0.6, label='Initial σy')
        ax2.set_xlabel('Equivalent Plastic Strain κ')
        ax2.set_ylabel('Hardened Yield Stress [MPa]')
        ax2.set_title('Hardening Curve')
        ax2.grid(True, alpha=0.3)
        ax2.legend()
        
        # 3. ステップ進行
        ax3.plot(steps, dp_values, 'go-', linewidth=2, markersize=6, label='DP Criterion')
        ax3.plot(steps, yield_stresses, 'mo-', linewidth=2, markersize=6, label='Yield Stress')
        ax3.set_xlabel('Load Step')
        ax3.set_ylabel('Stress [MPa]')
        ax3.set_title('Step Progression')
        ax3.grid(True, alpha=0.3)
        ax3.legend()
        
        # 4. 応力差
        stress_diff = np.array(dp_values) - np.array(yield_stresses)
        ax4.plot(plastic_strains, stress_diff, 'o-', color='purple', linewidth=2, markersize=6)
        ax4.axhline(y=0, color='black', linestyle='-', alpha=0.5)
        ax4.fill_between(plastic_strains, stress_diff, 0, 
                        where=stress_diff < 0, alpha=0.3, color='blue', label='Elastic')
        ax4.fill_between(plastic_strains, stress_diff, 0, 
                        where=stress_diff >= 0, alpha=0.3, color='red', label='Plastic')
        ax4.set_xlabel('Equivalent Plastic Strain κ')
        ax4.set_ylabel('Stress Difference (DP - σy) [MPa]')
        ax4.set_title('Stress Difference Analysis')
        ax4.grid(True, alpha=0.3)
        ax4.legend()
        
        plt.tight_layout()
        
        # 保存
        output_path = self.output_dir / f"element{element_num}_comprehensive_analysis.png"
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        print(f"Comprehensive analysis saved: {output_path.name}")
        plt.close()
        
        return results
    
    def plot_dual_axis(self, filename: str, element_num: int = 8):
        """2軸プロット（DP基準 vs 硬化応力）"""
        stress_data, strain_data = self.read_res_file(filename)
        results = self.analyze_element(stress_data, strain_data, element_num)
        
        if not results['steps']:
            print(f"No data found for element {element_num}")
            return
        
        steps = results['steps']
        dp_values = results['dp_values']
        plastic_strains = results['plastic_strains']
        yield_stresses = results['yield_stresses']
        
        # 2軸プロット
        fig, ax1 = plt.subplots(figsize=(14, 8))
        
        # 主軸：DP基準
        color1 = 'tab:red'
        ax1.set_xlabel('Equivalent Plastic Strain κ', fontsize=12)
        ax1.set_ylabel('Drucker-Prager Criterion [MPa]', color=color1, fontsize=12)
        line1 = ax1.plot(plastic_strains, dp_values, 'o-', color=color1, 
                        linewidth=3, markersize=8, label='DP Criterion', alpha=0.8)
        ax1.tick_params(axis='y', labelcolor=color1)
        ax1.grid(True, alpha=0.3)
        
        # 副軸：硬化応力
        ax2 = ax1.twinx()
        color2 = 'tab:blue'
        ax2.set_ylabel('Hardened Yield Stress [MPa]', color=color2, fontsize=12)
        line2 = ax2.plot(plastic_strains, yield_stresses, 's-', color=color2,
                        linewidth=3, markersize=8, label='Hardened Yield Stress', alpha=0.8)
        ax2.tick_params(axis='y', labelcolor=color2)
        
        # 初期降伏応力の参照線
        initial_yield = self.material_params['sigma_y']
        ax2.axhline(y=initial_yield, color='black', linestyle='--', alpha=0.6,
                   linewidth=2, label=f'Initial σy = {initial_yield} MPa')
        
        # ステップ注釈
        for i, (strain, dp, ys, step) in enumerate(zip(plastic_strains, dp_values, yield_stresses, steps)):
            if i % 2 == 0:
                ax1.annotate(f'S{step}', (strain, dp), textcoords="offset points",
                            xytext=(5, 10), ha='left', fontsize=9, color=color1, alpha=0.8)
        
        # 凡例
        lines = line1 + line2 + [plt.Line2D([0], [0], color='black', linestyle='--', alpha=0.6)]
        labels = ['DP Criterion', 'Hardened Yield Stress', f'Initial σy = {initial_yield} MPa']
        ax1.legend(lines, labels, loc='upper left', bbox_to_anchor=(0.02, 0.98), fontsize=11)
        
        plt.title(f'Element {element_num}: Dual-Axis Analysis\n'
                 f'Drucker-Prager Criterion & Hardened Yield Stress vs Plastic Strain', 
                 fontsize=14, fontweight='bold', pad=20)
        
        # 軸範囲調整
        ax1.set_xlim(0, max(plastic_strains) * 1.05)
        ax1.set_ylim(0, max(dp_values) * 1.05)
        ax2.set_ylim(initial_yield * 0.998, max(yield_stresses) * 1.002)
        
        plt.tight_layout()
        
        # 保存
        output_path = self.output_dir / f"element{element_num}_dual_axis.png"
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        print(f"Dual-axis plot saved: {output_path.name}")
        plt.close()
        
        # 数値データ表示
        self._print_analysis_table(steps, plastic_strains, dp_values, yield_stresses)
        
        return results
    
    def plot_hardening_curve(self, max_strain: float = 0.01):
        """硬化曲線の理論値プロット"""
        strains = np.linspace(0, max_strain, 1000)
        stresses = [self.calculate_hardened_yield_stress(k) for k in strains]
        
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.plot(strains, stresses, 'b-', linewidth=2.5, label='Theoretical Hardening Curve')
        ax.axhline(y=self.material_params['sigma_y'], color='red', 
                  linestyle='--', alpha=0.7, label='Initial Yield Stress')
        ax.axhline(y=self.material_params['hpa'], color='green', 
                  linestyle='--', alpha=0.7, label='Asymptotic Stress')
        
        ax.set_xlabel('Equivalent Plastic Strain κ')
        ax.set_ylabel('Hardened Yield Stress [MPa]')
        ax.set_title('Theoretical Hardening Curve\n'
                    f'σy={self.material_params["sigma_y"]}, hk={self.material_params["hk"]}, '
                    f'hpa={self.material_params["hpa"]}, hpb={self.material_params["hpb"]}')
        ax.grid(True, alpha=0.3)
        ax.legend()
        
        plt.tight_layout()
        
        output_path = self.output_dir / "theoretical_hardening_curve.png"
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        print(f"Theoretical hardening curve saved: {output_path.name}")
        plt.close()
    
    def _print_analysis_table(self, steps, plastic_strains, dp_values, yield_stresses):
        """解析結果テーブルを表示"""
        print(f"\n=== Analysis Results ===")
        print("Step  κ         DP[MPa]   σy[MPa]   DP-σy[MPa]")
        print("-" * 50)
        for step, strain, dp, ys in zip(steps, plastic_strains, dp_values, yield_stresses):
            stress_diff = dp - ys
            print(f"{step:2d}    {strain:.6f}  {dp:.1f}     {ys:.1f}     {stress_diff:+.1f}")
        
        stress_differences = [dp - ys for dp, ys in zip(dp_values, yield_stresses)]
        print(f"\nStatistics:")
        print(f"DP Criterion range: {min(dp_values):.1f} - {max(dp_values):.1f} MPa")
        print(f"Yield stress range: {min(yield_stresses):.1f} - {max(yield_stresses):.1f} MPa")
        print(f"Stress difference range: {min(stress_differences):+.1f} - {max(stress_differences):+.1f} MPa")
    
    def save_analysis_data(self, results: Dict, filename: str, element_num: int):
        """解析結果をCSVファイルに保存"""
        output_path = self.output_dir / f"element{element_num}_analysis_{filename.replace('.cml', '.csv')}"
        
        with open(output_path, 'w') as f:
            f.write("Step,PlasticStrain,DruckerPrager,YieldStress,StressDifference\n")
            for step, strain, dp, ys in zip(results['steps'], results['plastic_strains'], 
                                          results['dp_values'], results['yield_stresses']):
                stress_diff = dp - ys
                f.write(f"{step},{strain:.8f},{dp:.6f},{ys:.6f},{stress_diff:.6f}\n")
        
        print(f"Analysis data saved: {output_path.name}")


def main():
    """メイン関数"""
    parser = argparse.ArgumentParser(description='plstss240 Unified Visualization Tool')
    parser.add_argument('--file', '-f', type=str, default='RES_cube250514.cml',
                       help='RES file name (default: RES_cube250514.cml)')
    parser.add_argument('--element', '-e', type=int, default=8,
                       help='Element number to analyze (default: 8)')
    parser.add_argument('--mode', '-m', choices=['all', 'dual', 'comprehensive', 'hardening'],
                       default='all', help='Analysis mode (default: all)')
    parser.add_argument('--output', '-o', type=str, default='../output',
                       help='Output directory (default: ../output)')
    
    args = parser.parse_args()
    
    # 可視化ツール初期化
    visualizer = PlstssVisualizer(args.output)
    
    print(f"=== plstss240 Unified Visualization Tool ===")
    print(f"File: {args.file}")
    print(f"Element: {args.element}")
    print(f"Mode: {args.mode}")
    print(f"Output: {args.output}")
    print()
    
    try:
        if args.mode in ['all', 'comprehensive']:
            print("Creating comprehensive analysis...")
            results = visualizer.plot_drucker_prager_main(args.file, args.element)
            visualizer.save_analysis_data(results, args.file, args.element)
        
        if args.mode in ['all', 'dual']:
            print("Creating dual-axis analysis...")
            visualizer.plot_dual_axis(args.file, args.element)
        
        if args.mode in ['all', 'hardening']:
            print("Creating theoretical hardening curve...")
            visualizer.plot_hardening_curve()
        
        print(f"\nAll analyses completed successfully!")
        print(f"Results saved in: {visualizer.output_dir}")
        
    except FileNotFoundError as e:
        print(f"Error: {e}")
        return 1
    except Exception as e:
        print(f"Unexpected error: {e}")
        return 1
    
    return 0


if __name__ == "__main__":
    sys.exit(main()) 