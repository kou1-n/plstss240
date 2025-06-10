#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
降伏応力硬化曲線解析ツール
stress_vm.fの硬化モデルに基づいて降伏応力の硬化を可視化
"""

import numpy as np
import matplotlib.pyplot as plt
import os
import sys
from pathlib import Path

# 日本語フォントの設定
plt.rcParams['font.family'] = ['DejaVu Sans', 'Arial Unicode MS', 'Hiragino Sans']
plt.rcParams['axes.unicode_minus'] = False


class HardeningAnalyzer:
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
        
        # Elastic constants
        self.G = self.E / (2.0 * (1.0 + self.nu))  # Shear modulus
        self.K = self.E / (3.0 * (1.0 - 2.0*self.nu))  # Bulk modulus
        
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
    
    def H_kin(self, kappa, theta=1.0):
        """
        運動硬化関数 - stress_vm.fのH_kin関数に対応
        H_kin = (1 - θ) * hk * κ
        """
        return (1.0 - theta) * self.hk * kappa
    
    def current_yield_stress(self, kappa):
        """
        現在の降伏応力 = 初期降伏応力 + 等方硬化成分
        """
        return self.H_iso(kappa)
    
    def read_plastic_strain_history(self, filename="drucker_prager_data_RES_cube250514.csv"):
        """
        既存の解析結果から相当塑性ひずみの履歴を読み込み
        """
        data_path = self.output_dir / filename
        if not data_path.exists():
            print(f"Error: {data_path} が見つかりません")
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
    
    def plot_hardening_curves(self):
        """
        硬化曲線の可視化
        """
        # 理論曲線用の相当塑性ひずみ範囲
        kappa_max = 0.003  # 最大相当塑性ひずみ
        kappa_theory = np.linspace(0, kappa_max, 1000)
        
        # 理論硬化曲線
        iso_hardening = self.H_iso(kappa_theory)
        iso_hardening_increment = iso_hardening - self.sigma_y
        kin_hardening = self.H_kin(kappa_theory, theta=1.0)  # theta=1で運動硬化なし
        
        # 実際の解析結果読み込み
        steps, actual_plastic_strains = self.read_plastic_strain_history()
        
        if steps is not None:
            actual_yield_stress = [self.current_yield_stress(kappa) for kappa in actual_plastic_strains]
            actual_iso_increment = [stress - self.sigma_y for stress in actual_yield_stress]
        
        # プロット作成
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(15, 12))
        fig.suptitle('降伏応力硬化曲線解析 (Element 8)', fontsize=16, fontweight='bold')
        
        # 1. 降伏応力 vs 相当塑性ひずみ
        ax1.plot(kappa_theory, iso_hardening, 'b-', linewidth=2, label='理論曲線')
        if steps is not None:
            ax1.plot(actual_plastic_strains, actual_yield_stress, 'ro-', markersize=6, 
                    linewidth=2, label='解析結果')
        ax1.axhline(y=self.sigma_y, color='k', linestyle='--', alpha=0.7, label=f'初期降伏応力: {self.sigma_y} MPa')
        ax1.axhline(y=self.hpa, color='g', linestyle='--', alpha=0.7, label=f'漸近応力: {self.hpa} MPa')
        ax1.set_xlabel('相当塑性ひずみ κ')
        ax1.set_ylabel('降伏応力 σy [MPa]')
        ax1.set_title('降伏応力の硬化')
        ax1.grid(True, alpha=0.3)
        ax1.legend()
        
        # 2. 等方硬化増分 vs 相当塑性ひずみ
        ax2.plot(kappa_theory, iso_hardening_increment, 'b-', linewidth=2, label='理論曲線')
        if steps is not None:
            ax2.plot(actual_plastic_strains, actual_iso_increment, 'ro-', markersize=6,
                    linewidth=2, label='解析結果')
        ax2.set_xlabel('相当塑性ひずみ κ')
        ax2.set_ylabel('等方硬化増分 [MPa]')
        ax2.set_title('等方硬化成分')
        ax2.grid(True, alpha=0.3)
        ax2.legend()
        
        # 3. 硬化率 dσy/dκ vs 相当塑性ひずみ
        hardening_rate = self.dH_iso_dk(kappa_theory)
        ax3.plot(kappa_theory, hardening_rate, 'g-', linewidth=2, label='理論硬化率')
        ax3.axhline(y=self.hk, color='r', linestyle='--', alpha=0.7, label=f'線形硬化率: {self.hk} MPa')
        ax3.set_xlabel('相当塑性ひずみ κ')
        ax3.set_ylabel('硬化率 dσy/dκ [MPa]')
        ax3.set_title('硬化率の変化')
        ax3.grid(True, alpha=0.3)
        ax3.legend()
        
        # 4. 硬化パラメータの影響比較
        hpb_values = [5.0, 10.0, 20.0]  # 異なるhpb値
        for hpb_val in hpb_values:
            temp_hardening = (self.sigma_y + self.hk * kappa_theory + 
                            (self.hpa - self.sigma_y) * (1.0 - np.exp(-hpb_val * kappa_theory)))
            ax4.plot(kappa_theory, temp_hardening, linewidth=2, 
                    label=f'hpb = {hpb_val}')
        
        ax4.axhline(y=self.sigma_y, color='k', linestyle='--', alpha=0.7)
        ax4.axhline(y=self.hpa, color='k', linestyle='--', alpha=0.7)
        ax4.set_xlabel('相当塑性ひずみ κ')
        ax4.set_ylabel('降伏応力 [MPa]')
        ax4.set_title('硬化パラメータhpbの影響')
        ax4.grid(True, alpha=0.3)
        ax4.legend()
        
        plt.tight_layout()
        
        # 画像保存
        output_path = self.output_dir / "hardening_curve_analysis.png"
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        print(f"硬化曲線グラフを保存しました: {output_path}")
        
        # データ保存
        self.save_hardening_data(kappa_theory, iso_hardening, hardening_rate,
                                actual_plastic_strains, actual_yield_stress)
        
        plt.show()
    
    def save_hardening_data(self, kappa_theory, iso_hardening, hardening_rate,
                           actual_plastic_strains=None, actual_yield_stress=None):
        """
        硬化曲線データをCSVファイルに保存
        """
        # 理論曲線データ
        theory_path = self.output_dir / "hardening_theory_curve.csv"
        with open(theory_path, 'w') as f:
            f.write("# 硬化曲線理論データ\n")
            f.write("# 材料パラメータ:\n")
            f.write(f"# σy = {self.sigma_y} MPa\n")
            f.write(f"# hk = {self.hk} MPa\n")
            f.write(f"# hpa = {self.hpa} MPa\n")
            f.write(f"# hpb = {self.hpb}\n")
            f.write("Equiv_Plastic_Strain,Yield_Stress[MPa],Hardening_Increment[MPa],Hardening_Rate[MPa]\n")
            
            for i, kappa in enumerate(kappa_theory):
                f.write(f"{kappa:.8f},{iso_hardening[i]:.6f},{iso_hardening[i]-self.sigma_y:.6f},{hardening_rate[i]:.6f}\n")
        
        print(f"理論硬化曲線データを保存しました: {theory_path}")
        
        # 実際の解析結果データ
        if actual_plastic_strains is not None and actual_yield_stress is not None:
            actual_path = self.output_dir / "hardening_actual_data.csv"
            with open(actual_path, 'w') as f:
                f.write("# 実際の解析結果による硬化データ (Element 8)\n")
                f.write("Equiv_Plastic_Strain,Yield_Stress[MPa],Hardening_Increment[MPa]\n")
                
                for i, kappa in enumerate(actual_plastic_strains):
                    increment = actual_yield_stress[i] - self.sigma_y
                    f.write(f"{kappa:.8f},{actual_yield_stress[i]:.6f},{increment:.6f}\n")
            
            print(f"実際の硬化データを保存しました: {actual_path}")
    
    def compare_hardening_models(self):
        """
        異なる硬化モデルの比較
        """
        kappa = np.linspace(0, 0.003, 1000)
        
        # 線形硬化
        linear_hardening = self.sigma_y + self.hk * kappa
        
        # 指数硬化（Voce型）
        voce_hardening = self.H_iso(kappa)
        
        # 複合硬化（線形+指数）
        combined_hardening = self.sigma_y + self.hk * kappa + (self.hpa - self.sigma_y) * (1.0 - np.exp(-self.hpb * kappa))
        
        plt.figure(figsize=(12, 8))
        
        plt.subplot(2, 2, 1)
        plt.plot(kappa, linear_hardening, 'b-', linewidth=2, label='線形硬化')
        plt.plot(kappa, voce_hardening, 'r-', linewidth=2, label='指数硬化 (Voce)')
        plt.axhline(y=self.sigma_y, color='k', linestyle='--', alpha=0.7, label='初期降伏応力')
        plt.axhline(y=self.hpa, color='g', linestyle='--', alpha=0.7, label='漸近応力')
        plt.xlabel('相当塑性ひずみ κ')
        plt.ylabel('降伏応力 [MPa]')
        plt.title('硬化モデルの比較')
        plt.grid(True, alpha=0.3)
        plt.legend()
        
        # 硬化率の比較
        plt.subplot(2, 2, 2)
        linear_rate = np.full_like(kappa, self.hk)
        voce_rate = self.dH_iso_dk(kappa)
        
        plt.plot(kappa, linear_rate, 'b-', linewidth=2, label='線形硬化率')
        plt.plot(kappa, voce_rate, 'r-', linewidth=2, label='指数硬化率')
        plt.xlabel('相当塑性ひずみ κ')
        plt.ylabel('硬化率 dσy/dκ [MPa]')
        plt.title('硬化率の比較')
        plt.grid(True, alpha=0.3)
        plt.legend()
        
        # 硬化増分の比較
        plt.subplot(2, 2, 3)
        plt.plot(kappa, linear_hardening - self.sigma_y, 'b-', linewidth=2, label='線形硬化増分')
        plt.plot(kappa, voce_hardening - self.sigma_y, 'r-', linewidth=2, label='指数硬化増分')
        plt.xlabel('相当塑性ひずみ κ')
        plt.ylabel('硬化増分 [MPa]')
        plt.title('硬化増分の比較')
        plt.grid(True, alpha=0.3)
        plt.legend()
        
        # パラメータ依存性
        plt.subplot(2, 2, 4)
        hk_values = [100, 200, 300]
        for hk_val in hk_values:
            temp_hardening = self.sigma_y + hk_val * kappa + (self.hpa - self.sigma_y) * (1.0 - np.exp(-self.hpb * kappa))
            plt.plot(kappa, temp_hardening, linewidth=2, label=f'hk = {hk_val} MPa')
        
        plt.axhline(y=self.sigma_y, color='k', linestyle='--', alpha=0.7)
        plt.xlabel('相当塑性ひずみ κ')
        plt.ylabel('降伏応力 [MPa]')
        plt.title('線形硬化係数hkの影響')
        plt.grid(True, alpha=0.3)
        plt.legend()
        
        plt.tight_layout()
        
        # 保存
        output_path = self.output_dir / "hardening_model_comparison.png"
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        print(f"硬化モデル比較グラフを保存しました: {output_path}")
        
        plt.show()
    
    def analyze_stress_vm_hardening(self):
        """
        stress_vm.fの硬化処理の詳細解析
        """
        print("=== stress_vm.f 硬化モデル解析 ===")
        print(f"材料パラメータ:")
        print(f"  初期降伏応力 σy = {self.sigma_y} MPa")
        print(f"  線形硬化係数 hk = {self.hk} MPa")
        print(f"  漸近硬化応力 hpa = {self.hpa} MPa")
        print(f"  硬化指数 hpb = {self.hpb}")
        print(f"  運動硬化係数 hpd = {self.hpd} MPa")
        print()
        
        print("硬化関数:")
        print("  等方硬化: H_iso(κ) = σy + hk*κ + (hpa - σy)*(1 - exp(-hpb*κ))")
        print("  運動硬化: H_kin(κ) = (1 - θ)*hk*κ  (現在の設定ではhpd=0により無効)")
        print()
        
        # サンプル計算
        kappa_samples = [0.0, 0.001, 0.002, 0.003]
        print("サンプル計算:")
        print("κ\t\tH_iso(κ)\tΔσy\t\tdH/dκ")
        print("-" * 50)
        for kappa in kappa_samples:
            h_iso = self.H_iso(kappa)
            delta_sy = h_iso - self.sigma_y
            dh_dk = self.dH_iso_dk(kappa)
            print(f"{kappa:.3f}\t\t{h_iso:.2f}\t\t{delta_sy:.2f}\t\t{dh_dk:.2f}")


def main():
    """メイン実行関数"""
    analyzer = HardeningAnalyzer()
    
    print("降伏応力硬化曲線解析を開始します...")
    
    # 硬化モデルの詳細解析
    analyzer.analyze_stress_vm_hardening()
    
    # 主要な硬化曲線解析
    analyzer.plot_hardening_curves()
    
    # 硬化モデルの比較
    analyzer.compare_hardening_models()
    
    print("解析完了！")


if __name__ == "__main__":
    main() 