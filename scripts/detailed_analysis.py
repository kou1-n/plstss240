#!/usr/bin/env python3
"""
詳細解析スクリプト：特定要素の応力-ひずみ関係や塑性挙動の詳細可視化
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
import os
import sys

# visualize_plastic_data.pyの関数を使用
sys.path.append('.')
from visualize_plastic_data import read_cml_file, calculate_equivalent_plastic_strain

rcParams['font.family'] = 'DejaVu Sans'
plt.rcParams['axes.unicode_minus'] = False

def plot_stress_strain_curves(step_data, output_dir='output'):
    """応力-ひずみ曲線をプロット"""
    
    steps = sorted(step_data.keys())
    num_elements = len(step_data[steps[0]])
    
    # ステップごとの変位を計算（Z方向の変位増分から）
    displacements = []
    for step in steps:
        displacement = step * 0.25e-3  # 各ステップ0.25mmの変位
        displacements.append(displacement)
    
    # 応力-ひずみ関係の可視化
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    
    # 中央要素（要素5）に注目
    element_idx = 4  # 0-indexedで要素5
    
    # 1. von Mises応力 vs 変位
    ax = axes[0, 0]
    von_mises_history = [step_data[step][element_idx]['von_mises'] for step in steps]
    ax.plot(displacements, von_mises_history, 'bo-', linewidth=2, markersize=6)
    ax.axhline(y=200, color='red', linestyle='--', linewidth=2, label='Yield Stress')
    ax.set_xlabel('Displacement (mm)')
    ax.set_ylabel('von Mises Stress (MPa)')
    ax.set_title(f'Element {element_idx+1}: von Mises Stress vs Displacement')
    ax.grid(True, alpha=0.3)
    ax.legend()
    
    # 2. 塑性ひずみエネルギー vs 変位
    ax = axes[0, 1]
    plastic_energy = [step_data[step][element_idx]['total_plastic_energy'] for step in steps]
    ax.plot(displacements, plastic_energy, 'go-', linewidth=2, markersize=6)
    ax.set_xlabel('Displacement (mm)')
    ax.set_ylabel('Plastic Energy')
    ax.set_title(f'Element {element_idx+1}: Plastic Energy vs Displacement')
    ax.grid(True, alpha=0.3)
    
    # 3. 相当塑性ひずみ vs 変位
    ax = axes[0, 2]
    equiv_plastic_strain = [step_data[step][element_idx]['equiv_plastic_strain'] for step in steps]
    ax.plot(displacements, equiv_plastic_strain, 'mo-', linewidth=2, markersize=6)
    ax.set_xlabel('Displacement (mm)')
    ax.set_ylabel('Equivalent Plastic Strain')
    ax.set_title(f'Element {element_idx+1}: Equiv. Plastic Strain vs Displacement')
    ax.grid(True, alpha=0.3)
    
    # 4. Z方向応力 vs Z方向ひずみ
    ax = axes[1, 0]
    z_stress = [step_data[step][element_idx]['stress'][2] for step in steps]  # σzz
    z_strain = [step_data[step][element_idx]['strain'][2] for step in steps]  # εzz
    ax.plot(z_strain, z_stress, 'co-', linewidth=2, markersize=6)
    ax.set_xlabel('Z-Strain')
    ax.set_ylabel('Z-Stress (MPa)')
    ax.set_title(f'Element {element_idx+1}: σzz vs εzz')
    ax.grid(True, alpha=0.3)
    
    # 5. 塑性変形開始の詳細（ステップ3-5）
    ax = axes[1, 1]
    plastic_start_steps = [3, 4, 5]
    plastic_start_von_mises = [step_data[step][element_idx]['von_mises'] for step in plastic_start_steps]
    plastic_start_displacement = [step * 0.25e-3 for step in plastic_start_steps]
    
    ax.plot(plastic_start_displacement, plastic_start_von_mises, 'ro-', linewidth=3, markersize=8)
    ax.axhline(y=200, color='red', linestyle='--', linewidth=2, label='Yield Stress')
    ax.set_xlabel('Displacement (mm)')
    ax.set_ylabel('von Mises Stress (MPa)')
    ax.set_title('Plastic Yield Detail (Steps 3-5)')
    ax.grid(True, alpha=0.3)
    ax.legend()
    
    # 6. 全要素の最終ステップでの塑性ひずみ分布
    ax = axes[1, 2]
    final_step = max(steps)
    final_plastic_strains = [step_data[final_step][i]['equiv_plastic_strain'] for i in range(num_elements)]
    element_numbers = list(range(1, num_elements + 1))
    
    bars = ax.bar(element_numbers, final_plastic_strains, color='skyblue', edgecolor='navy')
    ax.set_xlabel('Element Number')
    ax.set_ylabel('Equivalent Plastic Strain')
    ax.set_title(f'Step {final_step}: Plastic Strain Distribution')
    ax.grid(True, alpha=0.3, axis='y')
    
    # 最大値の要素をハイライト
    max_idx = np.argmax(final_plastic_strains)
    bars[max_idx].set_color('red')
    
    plt.tight_layout()
    
    # 保存
    output_file = os.path.join(output_dir, 'detailed_analysis.png')
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"詳細解析グラフを保存しました: {output_file}")
    
    plt.show()

def analyze_plastic_progression(step_data, output_dir='output'):
    """塑性変形の進展解析"""
    
    steps = sorted(step_data.keys())
    num_elements = len(step_data[steps[0]])
    
    print("\n=== 塑性変形進展解析 ===")
    print(f"解析対象：{num_elements}要素、{len(steps)}ステップ")
    
    # 各ステップで塑性変形している要素数をカウント
    plastic_elements_count = []
    for step in steps:
        plastic_count = 0
        for elem in step_data[step]:
            if elem['total_plastic_energy'] > 1e-12:  # 塑性変形の閾値
                plastic_count += 1
        plastic_elements_count.append(plastic_count)
        print(f"Step {step}: 塑性変形要素数 = {plastic_count}/{num_elements}")
    
    # 塑性変形開始ステップを特定
    first_plastic_step = None
    for i, count in enumerate(plastic_elements_count):
        if count > 0:
            first_plastic_step = steps[i]
            break
    
    if first_plastic_step:
        print(f"\n塑性変形開始ステップ: {first_plastic_step}")
        print(f"Step {first_plastic_step}での要素別詳細:")
        for i, elem in enumerate(step_data[first_plastic_step]):
            if elem['total_plastic_energy'] > 1e-12:
                print(f"  Element {i+1}: von Mises = {elem['von_mises']:.2f} MPa (降伏)")
            else:
                print(f"  Element {i+1}: von Mises = {elem['von_mises']:.2f} MPa (弾性)")
    
    # 最大応力要素の特定
    max_stress_data = []
    for step in steps:
        max_stress = 0
        max_element = 0
        for i, elem in enumerate(step_data[step]):
            if elem['von_mises'] > max_stress:
                max_stress = elem['von_mises']
                max_element = i + 1
        max_stress_data.append((step, max_element, max_stress))
        print(f"Step {step}: 最大von Mises応力 = {max_stress:.2f} MPa (Element {max_element})")
    
    # CSVで詳細データを出力
    detail_csv = os.path.join(output_dir, 'plastic_progression.csv')
    with open(detail_csv, 'w') as f:
        f.write('Step,Plastic_Elements_Count,Max_von_Mises_Element,Max_von_Mises_Stress\n')
        for i, step in enumerate(steps):
            step_num, max_elem, max_stress = max_stress_data[i]
            f.write(f'{step},{plastic_elements_count[i]},{max_elem},{max_stress:.6f}\n')
    
    print(f"\n詳細データを保存しました: {detail_csv}")

def main():
    """メイン関数"""
    cml_file = 'output/RES_cube250514.cml'
    
    if not os.path.exists(cml_file):
        print(f"エラー: ファイルが見つかりません: {cml_file}")
        return
    
    print("詳細解析を開始...")
    
    # データ読み込み
    step_data = read_cml_file(cml_file)
    step_data = calculate_equivalent_plastic_strain(step_data)
    
    # 詳細可視化
    plot_stress_strain_curves(step_data)
    
    # 塑性進展解析
    analyze_plastic_progression(step_data)

if __name__ == '__main__':
    main() 