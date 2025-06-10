#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
硬化関数の線形性解析ツール
なぜ硬化降伏応力が線形に見えるのかを詳細分析
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')  # 非対話モード
import matplotlib.pyplot as plt
from pathlib import Path

def H_iso(kappa, sigma_y=200.0, hk=200.0, hpa=400.0, hpb=10.0):
    """等方硬化関数"""
    return sigma_y + hk*kappa + (hpa - sigma_y)*(1.0 - np.exp(-hpb*kappa))

def H_iso_linear_only(kappa, sigma_y=200.0, hk=200.0):
    """線形項のみの硬化関数"""
    return sigma_y + hk*kappa

def H_iso_nonlinear_only(kappa, sigma_y=200.0, hpa=400.0, hpb=10.0):
    """非線形項のみの硬化関数"""
    return sigma_y + (hpa - sigma_y)*(1.0 - np.exp(-hpb*kappa))

def analyze_hardening_linearity():
    """硬化関数の線形性を詳細分析"""
    
    # 実際の解析データ（element 8）
    actual_plastic_strains = [
        0.00021868, 0.00043736, 0.00065604, 0.00086767, 0.00105206,
        0.00125781, 0.00147247, 0.00168997, 0.00190826, 0.00212708
    ]
    
    actual_yield_stresses = [
        200.480618, 200.960282, 201.438994, 201.901367, 202.303502,
        202.751427, 203.217858, 203.689534, 204.161988, 204.634650
    ]
    
    # 材料パラメータ
    sigma_y = 200.0
    hk = 200.0
    hpa = 400.0
    hpb = 10.0
    
    print("=== 硬化関数の線形性解析 ===")
    print(f"材料パラメータ:")
    print(f"  σy = {sigma_y} MPa")
    print(f"  hk = {hk} MPa") 
    print(f"  hpa = {hpa} MPa")
    print(f"  hpb = {hpb}")
    print()
    
    # 実際のデータ範囲での分析
    kappa_min = min(actual_plastic_strains)
    kappa_max = max(actual_plastic_strains)
    
    print(f"実際のデータ範囲:")
    print(f"  κ範囲: {kappa_min:.6f} → {kappa_max:.6f}")
    print(f"  Δκ = {kappa_max - kappa_min:.6f}")
    print()
    
    # 各項の寄与を計算
    print("各項の寄与分析:")
    print("κ        線形項    非線形項   合計     線形近似   誤差[%]")
    print("-" * 60)
    
    for kappa in actual_plastic_strains:
        linear_term = hk * kappa
        nonlinear_term = (hpa - sigma_y) * (1.0 - np.exp(-hpb * kappa))
        total_hardening = sigma_y + linear_term + nonlinear_term
        linear_approx = sigma_y + hk * kappa
        error_percent = abs(total_hardening - linear_approx) / total_hardening * 100
        
        print(f"{kappa:.6f}  {linear_term:.3f}     {nonlinear_term:.3f}     {total_hardening:.3f}    {linear_approx:.3f}      {error_percent:.3f}")
    
    print()
    
    # 非線形項の寄与率
    kappa_range = np.array(actual_plastic_strains)
    linear_contributions = hk * kappa_range
    nonlinear_contributions = (hpa - sigma_y) * (1.0 - np.exp(-hpb * kappa_range))
    
    avg_linear = np.mean(linear_contributions)
    avg_nonlinear = np.mean(nonlinear_contributions)
    
    print(f"平均寄与:")
    print(f"  線形項: {avg_linear:.3f} MPa ({avg_linear/(avg_linear+avg_nonlinear)*100:.1f}%)")
    print(f"  非線形項: {avg_nonlinear:.3f} MPa ({avg_nonlinear/(avg_linear+avg_nonlinear)*100:.1f}%)")
    print()
    
    # 最大誤差の計算
    full_hardening = [H_iso(k, sigma_y, hk, hpa, hpb) for k in actual_plastic_strains]
    linear_hardening = [H_iso_linear_only(k, sigma_y, hk) for k in actual_plastic_strains]
    
    max_abs_error = max(abs(f - l) for f, l in zip(full_hardening, linear_hardening))
    max_rel_error = max(abs(f - l)/f for f, l in zip(full_hardening, linear_hardening)) * 100
    
    print(f"線形近似の誤差:")
    print(f"  最大絶対誤差: {max_abs_error:.4f} MPa")
    print(f"  最大相対誤差: {max_rel_error:.3f}%")
    print()
    
    # 可視化
    create_detailed_plots(actual_plastic_strains, actual_yield_stresses,
                         sigma_y, hk, hpa, hpb)

def create_detailed_plots(actual_strains, actual_stresses, sigma_y, hk, hpa, hpb):
    """詳細なプロットを作成"""
    
    # 1. 実際のデータ範囲での詳細分析
    kappa_detail = np.linspace(0, max(actual_strains)*1.2, 1000)
    
    # 2. より広い範囲での比較
    kappa_wide = np.linspace(0, 0.01, 1000)  # 10倍広い範囲
    
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(16, 12))
    fig.suptitle('Hardening Function Linearity Analysis', fontsize=16, fontweight='bold')
    
    # プロット1: 実際のデータ範囲での比較
    full_hardening_detail = [H_iso(k, sigma_y, hk, hpa, hpb) for k in kappa_detail]
    linear_hardening_detail = [H_iso_linear_only(k, sigma_y, hk) for k in kappa_detail]
    
    ax1.plot(kappa_detail, full_hardening_detail, 'b-', linewidth=2, label='Full (Linear + Nonlinear)')
    ax1.plot(kappa_detail, linear_hardening_detail, 'r--', linewidth=2, label='Linear only')
    ax1.plot(actual_strains, actual_stresses, 'ko', markersize=6, label='Actual FEM data')
    
    ax1.set_xlabel('Equivalent Plastic Strain κ')
    ax1.set_ylabel('Hardened Yield Stress [MPa]')
    ax1.set_title('Actual Data Range (Small κ)')
    ax1.grid(True, alpha=0.3)
    ax1.legend()
    
    # プロット2: 広い範囲での比較
    full_hardening_wide = [H_iso(k, sigma_y, hk, hpa, hpb) for k in kappa_wide]
    linear_hardening_wide = [H_iso_linear_only(k, sigma_y, hk) for k in kappa_wide]
    nonlinear_hardening_wide = [H_iso_nonlinear_only(k, sigma_y, hpa, hpb) for k in kappa_wide]
    
    ax2.plot(kappa_wide, full_hardening_wide, 'b-', linewidth=2, label='Full (Linear + Nonlinear)')
    ax2.plot(kappa_wide, linear_hardening_wide, 'r--', linewidth=2, label='Linear only')
    ax2.plot(kappa_wide, nonlinear_hardening_wide, 'g:', linewidth=2, label='Nonlinear only')
    ax2.axhline(y=hpa, color='k', linestyle='-', alpha=0.5, label=f'Asymptotic σy = {hpa} MPa')
    
    # 実際のデータ範囲をハイライト
    ax2.axvspan(0, max(actual_strains), alpha=0.2, color='yellow', label='Actual data range')
    
    ax2.set_xlabel('Equivalent Plastic Strain κ')
    ax2.set_ylabel('Hardened Yield Stress [MPa]')
    ax2.set_title('Wide Range Comparison (10× larger κ)')
    ax2.grid(True, alpha=0.3)
    ax2.legend()
    
    # プロット3: 寄与率の分析
    kappa_contrib = np.array(actual_strains)
    linear_contrib = hk * kappa_contrib
    nonlinear_contrib = (hpa - sigma_y) * (1.0 - np.exp(-hpb * kappa_contrib))
    
    ax3.plot(actual_strains, linear_contrib, 'r-o', linewidth=2, markersize=6, label='Linear contribution')
    ax3.plot(actual_strains, nonlinear_contrib, 'g-s', linewidth=2, markersize=6, label='Nonlinear contribution')
    
    ax3.set_xlabel('Equivalent Plastic Strain κ')
    ax3.set_ylabel('Stress Contribution [MPa]')
    ax3.set_title('Individual Term Contributions')
    ax3.grid(True, alpha=0.3)
    ax3.legend()
    
    # プロット4: 相対誤差の分析
    full_values = [H_iso(k, sigma_y, hk, hpa, hpb) for k in actual_strains]
    linear_values = [H_iso_linear_only(k, sigma_y, hk) for k in actual_strains]
    relative_errors = [abs(f - l)/f * 100 for f, l in zip(full_values, linear_values)]
    
    ax4.plot(actual_strains, relative_errors, 'mo-', linewidth=2, markersize=6)
    ax4.set_xlabel('Equivalent Plastic Strain κ')
    ax4.set_ylabel('Relative Error [%]')
    ax4.set_title('Linear Approximation Error')
    ax4.grid(True, alpha=0.3)
    
    # 平均誤差ライン
    avg_error = np.mean(relative_errors)
    ax4.axhline(y=avg_error, color='r', linestyle='--', alpha=0.7, 
               label=f'Average error: {avg_error:.3f}%')
    ax4.legend()
    
    plt.tight_layout()
    
    # 保存
    output_path = Path("../output") / "hardening_linearity_analysis.png"
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"詳細解析グラフを保存: {output_path.name}")
    plt.close()
    
    # 追加: 硬化率の分析
    plot_hardening_rate_analysis(actual_strains, sigma_y, hk, hpa, hpb)

def plot_hardening_rate_analysis(actual_strains, sigma_y, hk, hpa, hpb):
    """硬化率の詳細分析"""
    
    def dH_dk(kappa, hk, hpa, sigma_y, hpb):
        """硬化率の導関数"""
        return hk + (hpa - sigma_y) * hpb * np.exp(-hpb * kappa)
    
    kappa_detail = np.linspace(0, max(actual_strains)*1.2, 1000)
    kappa_wide = np.linspace(0, 0.01, 1000)
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    fig.suptitle('Hardening Rate Analysis', fontsize=16, fontweight='bold')
    
    # 実際のデータ範囲での硬化率
    hardening_rates_detail = [dH_dk(k, hk, hpa, sigma_y, hpb) for k in kappa_detail]
    
    ax1.plot(kappa_detail, hardening_rates_detail, 'b-', linewidth=2, label='Hardening rate dσy/dκ')
    ax1.axhline(y=hk, color='r', linestyle='--', alpha=0.7, label=f'Linear rate: {hk} MPa')
    
    # 実際のデータポイントでの硬化率
    actual_rates = [dH_dk(k, hk, hpa, sigma_y, hpb) for k in actual_strains]
    ax1.plot(actual_strains, actual_rates, 'ro', markersize=6, label='Actual data points')
    
    ax1.set_xlabel('Equivalent Plastic Strain κ')
    ax1.set_ylabel('Hardening Rate dσy/dκ [MPa]')
    ax1.set_title('Actual Data Range')
    ax1.grid(True, alpha=0.3)
    ax1.legend()
    
    # 広い範囲での硬化率
    hardening_rates_wide = [dH_dk(k, hk, hpa, sigma_y, hpb) for k in kappa_wide]
    
    ax2.plot(kappa_wide, hardening_rates_wide, 'b-', linewidth=2, label='Hardening rate dσy/dκ')
    ax2.axhline(y=hk, color='r', linestyle='--', alpha=0.7, label=f'Asymptotic rate: {hk} MPa')
    
    # 初期硬化率
    initial_rate = dH_dk(0, hk, hpa, sigma_y, hpb)
    ax2.axhline(y=initial_rate, color='g', linestyle=':', alpha=0.7, 
               label=f'Initial rate: {initial_rate:.0f} MPa')
    
    # 実際のデータ範囲をハイライト
    ax2.axvspan(0, max(actual_strains), alpha=0.2, color='yellow', label='Actual range')
    
    ax2.set_xlabel('Equivalent Plastic Strain κ')
    ax2.set_ylabel('Hardening Rate dσy/dκ [MPa]')
    ax2.set_title('Wide Range (showing nonlinearity)')
    ax2.grid(True, alpha=0.3)
    ax2.legend()
    
    plt.tight_layout()
    
    # 保存
    output_path = Path("../output") / "hardening_rate_analysis.png"
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"硬化率解析グラフを保存: {output_path.name}")
    plt.close()
    
    # 硬化率の数値分析
    print("硬化率の数値分析:")
    print("κ        dσy/dκ[MPa]  vs線形[%]")
    print("-" * 35)
    for i, kappa in enumerate(actual_strains):
        rate = dH_dk(kappa, hk, hpa, sigma_y, hpb)
        vs_linear = (rate - hk) / hk * 100
        print(f"{kappa:.6f}  {rate:.1f}        {vs_linear:+.2f}")
    
    print(f"\n初期硬化率: {initial_rate:.0f} MPa")
    print(f"漸近硬化率: {hk} MPa")
    print(f"実際の範囲での硬化率変化: {max(actual_rates) - min(actual_rates):.1f} MPa")

if __name__ == "__main__":
    print("硬化関数の線形性解析を開始...")
    analyze_hardening_linearity()
    print("\n解析完了！") 