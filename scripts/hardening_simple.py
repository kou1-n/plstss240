#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
シンプルな硬化曲線解析ツール
stress_vm.fの硬化モデルに基づいて降伏応力の硬化を可視化（非対話モード）
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')  # 非対話モード
import matplotlib.pyplot as plt
import os
import sys
from pathlib import Path

# 日本語フォントの設定
plt.rcParams['font.family'] = ['DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False


class SimpleHardeningAnalyzer:
    def __init__(self, output_dir="../output"):
        self.output_dir = Path(output_dir)
        
        # 材料パラメータ（cube250514.cmlから）
        self.E = 200000.0      # Young's modulus [MPa]
        self.nu = 0.3          # Poisson's ratio
        self.sigma_y = 200.0   # Initial yield stress [MPa]
        self.hk = 200.0        # Hardening parameter [MPa]
        self.hpa = 400.0       # Asymptotic hardening stress [MPa]
        self.hpb = 10.0        # Hardening exponent
        self.hpd = 0.0         # Kinematic hardening parameter
        
    def H_iso(self, kappa):
        """
        等方硬化関数 - stress_vm.fのH_iso関数に対応
        K(κ) = σy + hk*κ + (hpa - σy)*(1 - exp(-hpb*κ))
        """
        return (self.sigma_y + self.hk * kappa + 
                (self.hpa - self.sigma_y) * (1.0 - np.exp(-self.hpb * kappa)))
    
    def dH_iso_dk(self, kappa):
        """
        等方硬化関数の微分 dK/dκ
        """
        return (self.hk + (self.hpa - self.sigma_y) * self.hpb * np.exp(-self.hpb * kappa))
    
    def read_plastic_strain_history(self, filename="drucker_prager_data_RES_cube250514.csv"):
        """
        既存の解析結果から相当塑性ひずみの履歴を読み込み
        """
        data_path = self.output_dir / filename
        if not data_path.exists():
            print(f"Warning: {data_path} not found")
            return None, None
            
        steps = []
        plastic_strains = []
        
        with open(data_path, 'r') as f:
            lines = f.readlines()
            
        for line in lines[4:]:  # ヘッダーをスキップ
            if line.strip() and not line.startswith('#'):
                parts = line.strip().split(',')
                if len(parts) >= 3:
                    steps.append(int(parts[0]))
                    plastic_strains.append(float(parts[2]))
        
        return steps, plastic_strains
    
    def analyze_hardening_behavior(self):
        """
        stress_vm.fの硬化動作の詳細解析
        """
        print("=== stress_vm.f Hardening Model Analysis ===")
        print(f"Material Parameters:")
        print(f"  Initial yield stress σy = {self.sigma_y} MPa")
        print(f"  Linear hardening coeff hk = {self.hk} MPa")
        print(f"  Asymptotic stress hpa = {self.hpa} MPa")
        print(f"  Hardening exponent hpb = {self.hpb}")
        print(f"  Kinematic hardening hpd = {self.hpd} MPa")
        print()
        
        print("Hardening Function:")
        print("  Isotropic: H_iso(κ) = σy + hk*κ + (hpa - σy)*(1 - exp(-hpb*κ))")
        print("  Kinematic: H_kin(κ) = (1 - θ)*hk*κ  (currently inactive: hpd=0)")
        print()
        
        # サンプル計算
        kappa_samples = [0.0, 0.0005, 0.001, 0.0015, 0.002, 0.0025]
        print("Sample Calculations:")
        print("κ        H_iso(κ)  Δσy     dH/dκ")
        print("-" * 40)
        for kappa in kappa_samples:
            h_iso = self.H_iso(kappa)
            delta_sy = h_iso - self.sigma_y
            dh_dk = self.dH_iso_dk(kappa)
            print(f"{kappa:.4f}   {h_iso:.2f}    {delta_sy:.2f}   {dh_dk:.2f}")
        
        return kappa_samples
    
    def plot_hardening_curves(self):
        """
        硬化曲線の可視化（非対話モード）
        """
        # 理論曲線用の相当塑性ひずみ範囲
        kappa_max = 0.003
        kappa_theory = np.linspace(0, kappa_max, 1000)
        
        # 理論硬化曲線
        iso_hardening = self.H_iso(kappa_theory)
        iso_hardening_increment = iso_hardening - self.sigma_y
        hardening_rate = self.dH_iso_dk(kappa_theory)
        
        # 実データ読み込み
        steps, actual_plastic_strains = self.read_plastic_strain_history()
        
        if steps is not None:
            actual_yield_stress = [self.H_iso(kappa) for kappa in actual_plastic_strains]
            actual_iso_increment = [stress - self.sigma_y for stress in actual_yield_stress]
        
        # プロット作成
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(15, 12))
        fig.suptitle('Yield Stress Hardening Analysis (Element 8)', fontsize=16, fontweight='bold')
        
        # 1. 降伏応力 vs 相当塑性ひずみ
        ax1.plot(kappa_theory, iso_hardening, 'b-', linewidth=2, label='Theoretical curve')
        if steps is not None:
            ax1.plot(actual_plastic_strains, actual_yield_stress, 'ro-', markersize=6, 
                    linewidth=2, label='Analysis results')
        ax1.axhline(y=self.sigma_y, color='k', linestyle='--', alpha=0.7, 
                   label=f'Initial yield: {self.sigma_y} MPa')
        ax1.axhline(y=self.hpa, color='g', linestyle='--', alpha=0.7, 
                   label=f'Asymptotic: {self.hpa} MPa')
        ax1.set_xlabel('Equivalent plastic strain κ')
        ax1.set_ylabel('Yield stress σy [MPa]')
        ax1.set_title('Yield Stress Hardening')
        ax1.grid(True, alpha=0.3)
        ax1.legend()
        
        # 2. 等方硬化増分 vs 相当塑性ひずみ
        ax2.plot(kappa_theory, iso_hardening_increment, 'b-', linewidth=2, label='Theoretical')
        if steps is not None:
            ax2.plot(actual_plastic_strains, actual_iso_increment, 'ro-', markersize=6,
                    linewidth=2, label='Analysis results')
        ax2.set_xlabel('Equivalent plastic strain κ')
        ax2.set_ylabel('Isotropic hardening increment [MPa]')
        ax2.set_title('Isotropic Hardening Component')
        ax2.grid(True, alpha=0.3)
        ax2.legend()
        
        # 3. 硬化率 dσy/dκ vs 相当塑性ひずみ
        ax3.plot(kappa_theory, hardening_rate, 'g-', linewidth=2, label='Theoretical rate')
        ax3.axhline(y=self.hk, color='r', linestyle='--', alpha=0.7, 
                   label=f'Linear rate: {self.hk} MPa')
        ax3.set_xlabel('Equivalent plastic strain κ')
        ax3.set_ylabel('Hardening rate dσy/dκ [MPa]')
        ax3.set_title('Hardening Rate Evolution')
        ax3.grid(True, alpha=0.3)
        ax3.legend()
        
        # 4. 異なる硬化パラメータの比較
        hpb_values = [5.0, 10.0, 20.0]
        for hpb_val in hpb_values:
            temp_hardening = (self.sigma_y + self.hk * kappa_theory + 
                            (self.hpa - self.sigma_y) * (1.0 - np.exp(-hpb_val * kappa_theory)))
            ax4.plot(kappa_theory, temp_hardening, linewidth=2, 
                    label=f'hpb = {hpb_val}')
        
        ax4.axhline(y=self.sigma_y, color='k', linestyle='--', alpha=0.7)
        ax4.axhline(y=self.hpa, color='k', linestyle='--', alpha=0.7)
        ax4.set_xlabel('Equivalent plastic strain κ')
        ax4.set_ylabel('Yield stress [MPa]')
        ax4.set_title('Effect of hardening parameter hpb')
        ax4.grid(True, alpha=0.3)
        ax4.legend()
        
        plt.tight_layout()
        
        # 画像保存
        output_path = self.output_dir / "hardening_analysis_simple.png"
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        print(f"Hardening curve plot saved: {output_path}")
        
        plt.close()  # メモリ解放
        
        # 理論vs実データの詳細比較
        self.plot_theory_vs_actual(kappa_theory, iso_hardening, 
                                  actual_plastic_strains, actual_yield_stress)
    
    def plot_theory_vs_actual(self, kappa_theory, iso_hardening, 
                             actual_plastic_strains, actual_yield_stress):
        """
        理論曲線と実データの詳細比較
        """
        if actual_plastic_strains is None:
            print("No actual data available for comparison")
            return
        
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
        fig.suptitle('Theoretical vs Actual Hardening Comparison', fontsize=14, fontweight='bold')
        
        # 1. 硬化曲線の比較
        ax1.plot(kappa_theory, iso_hardening, 'b-', linewidth=2, label='Theoretical')
        ax1.plot(actual_plastic_strains, actual_yield_stress, 'ro-', markersize=8, 
                linewidth=2, label='FEM Analysis')
        ax1.axhline(y=self.sigma_y, color='k', linestyle='--', alpha=0.7, 
                   label=f'σy = {self.sigma_y} MPa')
        ax1.axhline(y=self.hpa, color='g', linestyle='--', alpha=0.7, 
                   label=f'σ∞ = {self.hpa} MPa')
        ax1.set_xlabel('Equivalent plastic strain κ')
        ax1.set_ylabel('Yield stress [MPa]')
        ax1.set_title('Hardening Curve Comparison')
        ax1.grid(True, alpha=0.3)
        ax1.legend()
        
        # 2. 硬化増分の比較
        theoretical_at_actual = [self.H_iso(kappa) for kappa in actual_plastic_strains]
        actual_increment = [stress - self.sigma_y for stress in actual_yield_stress]
        theory_increment = [stress - self.sigma_y for stress in theoretical_at_actual]
        
        x_pos = range(len(actual_plastic_strains))
        width = 0.35
        
        ax2.bar([x - width/2 for x in x_pos], theory_increment, width, 
               label='Theoretical', alpha=0.7, color='blue')
        ax2.bar([x + width/2 for x in x_pos], actual_increment, width, 
               label='FEM Analysis', alpha=0.7, color='red')
        
        ax2.set_xlabel('Analysis Step')
        ax2.set_ylabel('Hardening increment [MPa]')
        ax2.set_title('Hardening Increment by Step')
        ax2.set_xticks(x_pos)
        ax2.set_xticklabels([f'S{i+4}' for i in range(len(x_pos))])  # ステップ4から開始
        ax2.grid(True, alpha=0.3)
        ax2.legend()
        
        plt.tight_layout()
        
        # 保存
        output_path = self.output_dir / "hardening_comparison.png"
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        print(f"Hardening comparison plot saved: {output_path}")
        
        plt.close()
        
        # 数値比較表示
        print("\nNumerical Comparison (Theory vs FEM):")
        print("Step  κ_actual   σy_theory  σy_actual  Error[%]")
        print("-" * 50)
        for i, kappa in enumerate(actual_plastic_strains):
            theory_stress = self.H_iso(kappa)
            actual_stress = actual_yield_stress[i]
            error = abs(theory_stress - actual_stress) / actual_stress * 100
            print(f"{i+4:2d}    {kappa:.6f}  {theory_stress:.2f}     {actual_stress:.2f}     {error:.2f}")


def main():
    """メイン実行関数"""
    analyzer = SimpleHardeningAnalyzer()
    
    print("Starting yield stress hardening curve analysis...")
    print()
    
    # 硬化モデルの詳細解析
    analyzer.analyze_hardening_behavior()
    
    # 硬化曲線のプロット
    analyzer.plot_hardening_curves()
    
    print("\nAnalysis completed successfully!")


if __name__ == "__main__":
    main() 