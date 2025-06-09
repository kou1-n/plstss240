#!/usr/bin/env python3
"""
可視化スクリプト：降伏応力とvon Mises応力、累積塑性ひずみエネルギーの時間変化
RES_*.cmlファイルから要素データを抽出してグラフ化
"""

import re
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
import os

# 日本語フォントの設定
rcParams['font.family'] = 'DejaVu Sans'
plt.rcParams['axes.unicode_minus'] = False

def read_cml_file(filename):
    """RES_*.cmlファイルを読み込んでデータを抽出"""
    with open(filename, 'r') as f:
        content = f.read()
    
    # ステップデータを格納する辞書
    step_data = {}
    
    # 各ステップのデータを抽出
    step_pattern = r'/NODAL/\s*(\d+)\s*\d+\s*\d+.*?/ELMTL/\s*\1\s*\d+\s*\d+\s*(\d+)(.*?)(?=/NODAL/|/ENDOF/)'
    
    for match in re.finditer(step_pattern, content, re.DOTALL):
        step_num = int(match.group(1))
        num_elements = int(match.group(2))
        element_data_text = match.group(3)
        
        # 要素データを解析
        elements = []
        lines = element_data_text.strip().split('\n')
        
        i = 0
        while i < len(lines):
            line = lines[i].strip()
            if line and line.split()[0].isdigit():  # 要素番号行
                element_num = int(line.split()[0])
                stress_line = line
                strain_line = lines[i+1].strip() if i+1 < len(lines) else ""
                other_line = lines[i+2].strip() if i+2 < len(lines) else ""
                
                # 応力成分を抽出（要素番号を除く）
                stress_values = [float(x) for x in stress_line.split()[1:7]]
                # ひずみ成分を抽出
                strain_values = [float(x) for x in strain_line.split()[0:6]]
                # その他の値を抽出（von Mises応力、塑性エネルギーなど）
                other_values = [float(x) for x in other_line.split()]
                
                von_mises_stress = other_values[0] if len(other_values) > 0 else 0.0
                plastic_energy_incr = other_values[1] if len(other_values) > 1 else 0.0
                plastic_strain_incr = other_values[2] if len(other_values) > 2 else 0.0
                plastic_energy_incr2 = other_values[3] if len(other_values) > 3 else 0.0
                total_plastic_energy = other_values[4] if len(other_values) > 4 else 0.0
                
                elements.append({
                    'element': element_num,
                    'stress': stress_values,  # [σxx, σyy, σzz, τxy, τxz, τyz]
                    'strain': strain_values,  # [εxx, εyy, εzz, γxy, γxz, γyz]
                    'von_mises': von_mises_stress,
                    'plastic_energy_incr': plastic_energy_incr,
                    'plastic_strain_incr': plastic_strain_incr,
                    'total_plastic_energy': total_plastic_energy
                })
                
                i += 3
            else:
                i += 1
        
        step_data[step_num] = elements
    
    return step_data

def calculate_equivalent_plastic_strain(step_data):
    """累積相当塑性ひずみを計算"""
    # 各要素について累積相当塑性ひずみを計算
    for step_num in sorted(step_data.keys()):
        for element in step_data[step_num]:
            # 簡易的な累積塑性ひずみ計算（塑性ひずみエネルギーから推定）
            # σy = 200 MPa とすると、εp_eq ≈ sqrt(2*Wp/σy)
            yield_stress = 200.0  # MPa
            if element['total_plastic_energy'] > 0:
                element['equiv_plastic_strain'] = np.sqrt(2 * element['total_plastic_energy'] / yield_stress)
            else:
                element['equiv_plastic_strain'] = 0.0
    
    return step_data

def plot_data(step_data, output_dir='output'):
    """データを可視化"""
    
    # 降伏応力（一定値）
    yield_stress = 200.0  # MPa
    
    # ステップ番号とデータを準備
    steps = sorted(step_data.keys())
    
    # 各要素のデータを抽出
    num_elements = len(step_data[steps[0]])
    
    # 図1：von Mises応力の時間変化
    plt.figure(figsize=(12, 8))
    
    plt.subplot(2, 2, 1)
    for elem_idx in range(min(4, num_elements)):  # 最初の4要素を表示
        von_mises_history = []
        for step in steps:
            von_mises_history.append(step_data[step][elem_idx]['von_mises'])
        plt.plot(steps, von_mises_history, 'o-', label=f'Element {elem_idx+1}')
    
    plt.axhline(y=yield_stress, color='red', linestyle='--', linewidth=2, label='Yield Stress (200 MPa)')
    plt.xlabel('Load Step')
    plt.ylabel('von Mises Stress (MPa)')
    plt.title('von Mises Stress vs Load Step')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    # 図2：最大von Mises応力の変化
    plt.subplot(2, 2, 2)
    max_von_mises = []
    avg_von_mises = []
    
    for step in steps:
        von_mises_values = [elem['von_mises'] for elem in step_data[step]]
        max_von_mises.append(max(von_mises_values))
        avg_von_mises.append(np.mean(von_mises_values))
    
    plt.plot(steps, max_von_mises, 'ro-', label='Maximum von Mises')
    plt.plot(steps, avg_von_mises, 'bo-', label='Average von Mises')
    plt.axhline(y=yield_stress, color='red', linestyle='--', linewidth=2, label='Yield Stress (200 MPa)')
    plt.xlabel('Load Step')
    plt.ylabel('von Mises Stress (MPa)')
    plt.title('Max/Avg von Mises Stress vs Load Step')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    # 図3：累積塑性ひずみエネルギーの変化
    plt.subplot(2, 2, 3)
    for elem_idx in range(min(4, num_elements)):
        plastic_energy_history = []
        for step in steps:
            plastic_energy_history.append(step_data[step][elem_idx]['total_plastic_energy'])
        plt.plot(steps, plastic_energy_history, 'o-', label=f'Element {elem_idx+1}')
    
    plt.xlabel('Load Step')
    plt.ylabel('Cumulative Plastic Energy')
    plt.title('Cumulative Plastic Energy vs Load Step')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    # 図4：相当塑性ひずみの変化
    plt.subplot(2, 2, 4)
    for elem_idx in range(min(4, num_elements)):
        equiv_strain_history = []
        for step in steps:
            equiv_strain_history.append(step_data[step][elem_idx]['equiv_plastic_strain'])
        plt.plot(steps, equiv_strain_history, 'o-', label=f'Element {elem_idx+1}')
    
    plt.xlabel('Load Step')
    plt.ylabel('Equivalent Plastic Strain')
    plt.title('Equivalent Plastic Strain vs Load Step')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    # 保存
    output_file = os.path.join(output_dir, 'plastic_analysis.png')
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"グラフを保存しました: {output_file}")
    
    plt.show()
    
    # 数値データもCSVで出力
    csv_file = os.path.join(output_dir, 'plastic_data.csv')
    with open(csv_file, 'w') as f:
        f.write('Step,Element,von_Mises_Stress,Total_Plastic_Energy,Equiv_Plastic_Strain\n')
        for step in steps:
            for elem_idx, element in enumerate(step_data[step]):
                f.write(f'{step},{elem_idx+1},{element["von_mises"]:.6f},'
                       f'{element["total_plastic_energy"]:.6f},{element["equiv_plastic_strain"]:.6f}\n')
    
    print(f"数値データを保存しました: {csv_file}")

def main():
    """メイン関数"""
    # 結果ファイルのパス
    cml_file = 'output/RES_cube250514.cml'
    
    if not os.path.exists(cml_file):
        print(f"エラー: ファイルが見つかりません: {cml_file}")
        print("まずplstss240を実行して結果ファイルを生成してください。")
        return
    
    print(f"結果ファイルを読み込み中: {cml_file}")
    
    # データを読み込み
    step_data = read_cml_file(cml_file)
    
    # 累積相当塑性ひずみを計算
    step_data = calculate_equivalent_plastic_strain(step_data)
    
    print(f"読み込み完了: {len(step_data)} ステップ")
    
    # 可視化
    plot_data(step_data)
    
    # ステップ4でのデータを詳細表示（塑性変形開始ステップ）
    if 4 in step_data:
        print("\n--- Step 4での詳細データ（塑性変形開始） ---")
        for i, elem in enumerate(step_data[4]):
            print(f"Element {i+1}: von Mises = {elem['von_mises']:.2f} MPa, "
                  f"Plastic Energy = {elem['total_plastic_energy']:.6f}, "
                  f"Equiv Plastic Strain = {elem['equiv_plastic_strain']:.6f}")

if __name__ == '__main__':
    main() 