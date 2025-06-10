#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
変形量増加後の解析結果詳細分析
前回結果との比較と硬化関数の非線形性評価
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path

def analyze_increased_deformation():
    """変形量増加後の詳細分析"""
    
    print("=== 変形量増加後の解析結果分析 ===")
    
    # 新しい解析結果（変形量10倍）
    new_results = {
        'steps': [1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
        'plastic_strains': [0.002117, 0.004323, 0.006551, 0.008790, 0.011039, 
                           0.013294, 0.015553, 0.017817, 0.020083, 0.022352],
        'dp_criterion': [245.7, 290.5, 327.7, 361.6, 394.0, 
                        425.7, 457.0, 488.0, 518.7, 549.3],
        'yield_stresses': [204.6, 209.3, 214.0, 218.6, 223.1, 
                          227.6, 231.9, 236.2, 240.4, 244.5],
        'stress_differences': [41.1, 81.2, 113.7, 143.0, 170.9, 
                              198.2, 225.1, 251.8, 278.3, 304.7]
    }
    
    # 以前の解析結果（変形量小）
    old_results = {
        'plastic_strains': [0.00021868, 0.00043736, 0.00065604, 0.00086767, 0.00105206,
                           0.00125781, 0.00147247, 0.00168997, 0.00190826, 0.00212708],
        'yield_stresses': [200.480618, 200.960282, 201.438994, 201.901367, 202.303502,
                          202.751427, 203.217858, 203.689534, 204.161988, 204.634650]
    }
    
    # 材料パラメータ
    sigma_y = 200.0
    hk = 200.0
    hpa = 400.0
    hpb = 10.0
    
    print("\n🔍 解析範囲比較:")
    print(f"変形量調整前: κ = {min(old_results['plastic_strains']):.6f} → {max(old_results['plastic_strains']):.6f}")
    print(f"変形量調整後: κ = {min(new_results['plastic_strains']):.6f} → {max(new_results['plastic_strains']):.6f}")
    print(f"範囲拡大倍率: {max(new_results['plastic_strains'])/max(old_results['plastic_strains']):.1f}倍")
    
    print(f"\n📊 硬化降伏応力範囲:")
    print(f"調整前: {min(old_results['yield_stresses']):.1f} → {max(old_results['yield_stresses']):.1f} MPa")
    print(f"調整後: {min(new_results['yield_stresses']):.1f} → {max(new_results['yield_stresses']):.1f} MPa")
    print(f"応力上昇量: {max(old_results['yield_stresses'])-min(old_results['yield_stresses']):.1f} → {max(new_results['yield_stresses'])-min(new_results['yield_stresses']):.1f} MPa")
    
    # 硬化関数の理論値計算
    def H_iso(kappa):
        return sigma_y + hk*kappa + (hpa - sigma_y)*(1.0 - np.exp(-hpb*kappa))
    
    def dH_dk(kappa):
        """硬化率"""
        return hk + (hpa - sigma_y) * hpb * np.exp(-hpb * kappa)
    
    # 非線形性の評価
    print(f"\n🔄 非線形性評価:")
    
    # 新しい結果での硬化率変化
    new_rates = [dH_dk(k) for k in new_results['plastic_strains']]
    old_rates = [dH_dk(k) for k in old_results['plastic_strains']]
    
    print(f"硬化率変化（調整前）: {max(old_rates) - min(old_rates):.1f} MPa")
    print(f"硬化率変化（調整後）: {max(new_rates) - min(new_rates):.1f} MPa")
    print(f"非線形性向上: {(max(new_rates) - min(new_rates))/(max(old_rates) - min(old_rates)):.1f}倍")
    
    # 線形近似との誤差
    print(f"\n📏 線形近似との乖離:")
    for i, (kappa, actual_stress) in enumerate(zip(new_results['plastic_strains'], new_results['yield_stresses'])):
        linear_approx = sigma_y + hk * kappa
        full_theory = H_iso(kappa)
        error_vs_linear = abs(actual_stress - linear_approx) / actual_stress * 100
        error_vs_theory = abs(actual_stress - full_theory) / actual_stress * 100
        
        if i in [0, 4, 9]:  # 代表的なステップ
            print(f"Step {i+1:2d}: κ={kappa:.4f}, 実測={actual_stress:.1f}, 線形={linear_approx:.1f}, 理論={full_theory:.1f}")
            print(f"         線形誤差={error_vs_linear:.2f}%, 理論誤差={error_vs_theory:.3f}%")
    
    # 可視化
    create_comparison_plots(new_results, old_results, sigma_y, hk, hpa, hpb)
    
    print(f"\n✅ 変形量増加の効果:")
    print(f"1. プラスチックひずみ範囲が10倍拡大")
    print(f"2. 硬化関数の非線形性が明確に観察可能")
    print(f"3. 線形近似との乖離が拡大 ({abs(new_results['yield_stresses'][-1] - (sigma_y + hk*new_results['plastic_strains'][-1])):.1f} MPa)")
    print(f"4. 材料の真の挙動特性が可視化")

def create_comparison_plots(new_results, old_results, sigma_y, hk, hpa, hpb):
    """比較プロット作成"""
    
    def H_iso(kappa):
        return sigma_y + hk*kappa + (hpa - sigma_y)*(1.0 - np.exp(-hpb*kappa))
    
    def H_linear(kappa):
        return sigma_y + hk*kappa
    
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(16, 12))
    fig.suptitle('Deformation Increase Analysis: Before vs After', fontsize=16, fontweight='bold')
    
    # プロット1: 硬化曲線の比較
    kappa_theory = np.linspace(0, max(new_results['plastic_strains'])*1.1, 1000)
    theory_curve = [H_iso(k) for k in kappa_theory]
    linear_curve = [H_linear(k) for k in kappa_theory]
    
    ax1.plot(kappa_theory, theory_curve, 'b-', linewidth=2, label='Theory (Nonlinear)')
    ax1.plot(kappa_theory, linear_curve, 'r--', linewidth=2, label='Linear approximation')
    ax1.plot(old_results['plastic_strains'], old_results['yield_stresses'], 'go', 
             markersize=6, label='Before (small strain)')
    ax1.plot(new_results['plastic_strains'], new_results['yield_stresses'], 'ro', 
             markersize=8, label='After (large strain)')
    
    ax1.set_xlabel('Equivalent Plastic Strain κ')
    ax1.set_ylabel('Hardened Yield Stress [MPa]')
    ax1.set_title('Hardening Curve Comparison')
    ax1.grid(True, alpha=0.3)
    ax1.legend()
    
    # プロット2: 相対誤差比較
    old_errors = []
    new_errors = []
    
    for kappa, stress in zip(old_results['plastic_strains'], old_results['yield_stresses']):
        linear = H_linear(kappa)
        error = abs(stress - linear) / stress * 100
        old_errors.append(error)
    
    for kappa, stress in zip(new_results['plastic_strains'], new_results['yield_stresses']):
        linear = H_linear(kappa)
        error = abs(stress - linear) / stress * 100
        new_errors.append(error)
    
    ax2.plot(old_results['plastic_strains'], old_errors, 'g-s', 
             linewidth=2, markersize=6, label='Before (small strain)')
    ax2.plot(new_results['plastic_strains'], new_errors, 'r-o', 
             linewidth=2, markersize=6, label='After (large strain)')
    
    ax2.set_xlabel('Equivalent Plastic Strain κ')
    ax2.set_ylabel('Error vs Linear Approximation [%]')
    ax2.set_title('Nonlinearity Detection')
    ax2.grid(True, alpha=0.3)
    ax2.legend()
    
    # プロット3: Drucker-Prager進展
    ax3.plot(new_results['plastic_strains'], new_results['dp_criterion'], 'b-o', 
             linewidth=2, markersize=6, label='DP Criterion')
    ax3.plot(new_results['plastic_strains'], new_results['yield_stresses'], 'r-s', 
             linewidth=2, markersize=6, label='Hardened Yield Stress')
    
    ax3.set_xlabel('Equivalent Plastic Strain κ')
    ax3.set_ylabel('Stress [MPa]')
    ax3.set_title('DP Criterion vs Yield Stress')
    ax3.grid(True, alpha=0.3)
    ax3.legend()
    
    # プロット4: 応力差の進展
    ax4.plot(new_results['steps'], new_results['stress_differences'], 'mo-', 
             linewidth=2, markersize=6)
    ax4.set_xlabel('Load Step')
    ax4.set_ylabel('DP - Yield Stress [MPa]')
    ax4.set_title('Stress Difference Evolution')
    ax4.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    # 保存
    output_path = Path("../output") / "deformation_increase_analysis.png"
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"\n📊 比較分析グラフを保存: {output_path.name}")
    plt.close()

if __name__ == "__main__":
    analyze_increased_deformation() 