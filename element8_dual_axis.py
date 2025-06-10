#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
要素8番のDrucker-Prager降伏基準値と硬化降伏応力の2軸プロット
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')  # 非対話モード
import matplotlib.pyplot as plt
from pathlib import Path

# フォント設定
plt.rcParams['font.family'] = ['DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

def create_dual_axis_plot():
    """
    2軸プロットを作成する関数
    """
    output_dir = Path("../output")
    
    # データの手動入力（既存のCSVファイルから取得済み）
    # Drucker-Prager基準値データ
    plastic_strains = [
        0.00021868, 0.00043736, 0.00065604, 0.00086767, 0.00105206,
        0.00125781, 0.00147247, 0.00168997, 0.00190826, 0.00212708
    ]
    
    dp_values = [
        52.099710, 104.202753, 156.300128, 203.760310, 216.151956,
        223.581788, 230.399459, 237.026191, 243.233027, 249.063495
    ]
    
    # 硬化降伏応力データ
    yield_stresses = [
        200.480618, 200.960282, 201.438994, 201.901367, 202.303502,
        202.751427, 203.217858, 203.689534, 204.161988, 204.634650
    ]
    
    steps = [4, 5, 6, 7, 8, 9, 10, 11, 12, 13]  # 対応するステップ番号
    
    # 図の作成
    fig, ax1 = plt.subplots(figsize=(14, 8))
    
    # 主軸：Drucker-Prager基準値
    color1 = 'tab:red'
    ax1.set_xlabel('Equivalent Plastic Strain κ', fontsize=12)
    ax1.set_ylabel('Drucker-Prager Criterion [MPa]', color=color1, fontsize=12)
    line1 = ax1.plot(plastic_strains, dp_values, 'o-', color=color1, 
                    linewidth=3, markersize=8, label='DP Criterion', alpha=0.8)
    ax1.tick_params(axis='y', labelcolor=color1, labelsize=11)
    ax1.grid(True, alpha=0.3)
    
    # 副軸：硬化後降伏応力
    ax2 = ax1.twinx()
    color2 = 'tab:blue'
    ax2.set_ylabel('Hardened Yield Stress σy [MPa]', color=color2, fontsize=12)
    line2 = ax2.plot(plastic_strains, yield_stresses, 's-', color=color2,
                    linewidth=3, markersize=8, label='Hardened Yield Stress', alpha=0.8)
    ax2.tick_params(axis='y', labelcolor=color2, labelsize=11)
    
    # 初期降伏応力の参照線
    initial_yield = 200.0
    ax2.axhline(y=initial_yield, color='black', linestyle='--', alpha=0.6,
               linewidth=2, label=f'Initial σy = {initial_yield} MPa')
    
    # ステップ番号の注釈
    for i, (strain, dp, ys, step) in enumerate(zip(plastic_strains, dp_values, yield_stresses, steps)):
        if i % 2 == 0:  # 隔ステップで表示して見やすくする
            ax1.annotate(f'S{step}', (strain, dp), textcoords="offset points",
                        xytext=(5, 10), ha='left', fontsize=9, color=color1, alpha=0.8)
            ax2.annotate(f'S{step}', (strain, ys), textcoords="offset points",
                        xytext=(-5, -15), ha='right', fontsize=9, color=color2, alpha=0.8)
    
    # 凡例の統合
    lines = line1 + line2 + [plt.Line2D([0], [0], color='black', linestyle='--', alpha=0.6)]
    labels = ['DP Criterion', 'Hardened Yield Stress', f'Initial σy = {initial_yield} MPa']
    ax1.legend(lines, labels, loc='upper left', bbox_to_anchor=(0.02, 0.98), fontsize=11)
    
    # タイトルと書式設定
    plt.title('Element 8: Drucker-Prager Criterion & Hardened Yield Stress\nvs Equivalent Plastic Strain', 
              fontsize=14, fontweight='bold', pad=20)
    
    # 軸の範囲調整
    ax1.set_xlim(0, max(plastic_strains) * 1.05)
    ax1.set_ylim(0, max(dp_values) * 1.05)
    ax2.set_ylim(initial_yield * 0.998, max(yield_stresses) * 1.002)
    
    plt.tight_layout()
    
    # 保存
    output_path = output_dir / "element8_dual_axis_plot.png"
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Dual-axis plot saved: {output_path}")
    
    plt.close()
    
    # 数値データの表示
    print("\n=== Element 8 Dual-Axis Analysis ===")
    print("Step  κ         DP[MPa]   σy[MPa]   DP-σy[MPa]")
    print("-" * 50)
    for i, (step, strain, dp, ys) in enumerate(zip(steps, plastic_strains, dp_values, yield_stresses)):
        stress_diff = dp - ys
        print(f"{step:2d}    {strain:.6f}  {dp:.1f}     {ys:.1f}     {stress_diff:+.1f}")
    
    # 統計情報
    stress_differences = [dp - ys for dp, ys in zip(dp_values, yield_stresses)]
    print(f"\nStatistics:")
    print(f"DP Criterion range: {min(dp_values):.1f} - {max(dp_values):.1f} MPa")
    print(f"Yield stress range: {min(yield_stresses):.1f} - {max(yield_stresses):.1f} MPa")
    print(f"Stress difference range: {min(stress_differences):+.1f} - {max(stress_differences):+.1f} MPa")
    
    # 差分グラフも作成
    create_difference_plot(plastic_strains, stress_differences, steps, output_dir)

def create_difference_plot(plastic_strains, stress_differences, steps, output_dir):
    """
    応力差分の詳細プロット
    """
    fig, ax = plt.subplots(figsize=(12, 6))
    
    # 差分プロット
    ax.plot(plastic_strains, stress_differences, 'o-', color='purple', 
           linewidth=3, markersize=8, alpha=0.8)
    ax.axhline(y=0, color='black', linestyle='-', alpha=0.5, linewidth=1)
    
    # ステップ番号の注釈
    for i, (strain, diff, step) in enumerate(zip(plastic_strains, stress_differences, steps)):
        ax.annotate(f'S{step}', (strain, diff), textcoords="offset points",
                   xytext=(0, 12), ha='center', fontsize=9, alpha=0.8)
    
    ax.set_xlabel('Equivalent Plastic Strain κ', fontsize=12)
    ax.set_ylabel('Stress Difference (DP - σy) [MPa]', fontsize=12)
    ax.set_title('Element 8: Difference between DP Criterion and Hardened Yield Stress', 
                fontsize=13, fontweight='bold')
    ax.grid(True, alpha=0.3)
    
    # 負の領域（弾性域）と正の領域（塑性域）を色分け
    ax.fill_between(plastic_strains, stress_differences, 0, 
                   where=np.array(stress_differences) < 0, alpha=0.3, color='blue', 
                   label='Elastic range (DP < σy)')
    ax.fill_between(plastic_strains, stress_differences, 0, 
                   where=np.array(stress_differences) >= 0, alpha=0.3, color='red',
                   label='Plastic range (DP ≥ σy)')
    
    ax.legend(fontsize=11)
    plt.tight_layout()
    
    # 保存
    output_path = output_dir / "element8_stress_difference_plot.png"
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Stress difference plot saved: {output_path}")
    
    plt.close()

if __name__ == "__main__":
    print("Creating dual-axis plot for Element 8...")
    create_dual_axis_plot()
    print("Dual-axis plot creation completed!") 