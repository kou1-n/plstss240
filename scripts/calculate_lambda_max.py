#!/usr/bin/env python3
"""
Drucker-Prager y軸圧縮・比例載荷設定
目標ステップ k* で f=0 になるよう λ_max を決定
"""

import math
import shutil
from datetime import datetime

# ===============================
# 入力パラメータ
# ===============================
k_star = 10  # 目標降伏ステップ（ステップ10で降伏f=0、ステップ11で塑性域）

# CMLファイルからの抽出値
N = 100  # 総ステップ数 (/SOLUT/)

# /MATER/からのDPパラメータ (line 22: φ, ψ in radians)
phi_rad = 1.74532E-01  # 10 degrees
psi_rad = 8.72664E-02  # 5 degrees
phi_deg = math.degrees(phi_rad)
psi_deg = math.degrees(psi_rad)

# 降伏応力 σ_y (cohesion相当) - /MATER/ line 20 の10番目
sigma_y = 200.0  # MPa

# 座標系（mm単位と仮定）
# 節点座標から載荷面を特定
# 節点3,4,7,8が載荷節点、節点3,4が上面(y=1.0)、節点7,8が下面(y=0.0)
# 載荷面として下面（y=0）を選択（圧縮側）
loading_nodes = [7, 8]  # 節点7,8が下面

# 載荷面の座標（x,z平面に投影）
# 節点7: (-0.5, 0.5, 0.0), 節点8: (0.5, 0.5, 0.0)
# これだけでは面積計算不可。全4節点（5,6,7,8）で下面を形成
bottom_face_coords = [
    (-0.5, -0.5, 0.0),  # 節点6
    (0.5, -0.5, 0.0),   # 節点5
    (0.5, 0.5, 0.0),    # 節点8
    (-0.5, 0.5, 0.0)    # 節点7
]

# 面積計算（x-z平面への投影）
def calculate_polygon_area(coords):
    """x-z平面に投影した多角形の面積を計算"""
    n = len(coords)
    area = 0.0
    for i in range(n):
        j = (i + 1) % n
        area += coords[i][0] * coords[j][2]
        area -= coords[j][0] * coords[i][2]
    return abs(area) / 2.0

# 実際は正方形なので
A = 1.0  # mm^2（1x1の正方形）
print(f"載荷面積 A = {A} mm^2")

# /LOADC/からの基準節点力（節点3,4,7,8に各-930.923）
# 下面の節点7,8の合計
F_y0 = -930.923 * 2  # N (2節点分の合計、負値は圧縮)
print(f"基準節点力 F_y^0 = {F_y0} N (圧縮)")

# ===============================
# Drucker-Pragerパラメータ計算
# ===============================

# η, ξ の計算（Drucker-Prager）
sin_phi = math.sin(phi_rad)
cos_phi = math.cos(phi_rad)
eta = 6 * sin_phi / (math.sqrt(3) * (3 - sin_phi))
xi = 6 * cos_phi / (math.sqrt(3) * (3 - sin_phi))

# k の計算（cohesion = σ_y として）
k = xi * sigma_y

print(f"\n=== Drucker-Pragerパラメータ ===")
print(f"φ = {phi_deg:.1f}°")
print(f"ψ = {psi_deg:.1f}°")
print(f"η = {eta:.6f}")
print(f"ξ = {xi:.6f}")
print(f"k = {k:.3f} MPa")

# ===============================
# 圧縮降伏応力の計算
# ===============================
# σ_y^(c) = k / (1 - η/3)
sigma_y_c = k / (1 - eta/3)
print(f"\n圧縮降伏応力 σ_y^(c) = {sigma_y_c:.3f} MPa")

# ===============================
# λ_max の計算
# ===============================
# λ_max = (N / k_star) * (A * σ_y^(c)) / |F_y^0|
# 単位整合性: MPa = N/mm^2、A[mm^2]、F_y^0[N] → OK

lambda_max = (N / k_star) * (A * sigma_y_c) / abs(F_y0)
print(f"\n=== 荷重係数 ===")
print(f"λ_max = {lambda_max:.6f}")

# ===============================
# 検証表の作成
# ===============================
print(f"\n=== 検証表 (目標ステップ k* = {k_star}) ===")

# データ作成
steps = range(1, N+1)
data = []

for s in steps:
    # ステップsでの節点力
    F_y_s = (s/N) * lambda_max * F_y0

    # 平均応力
    sigma_s = abs(F_y_s) / A  # MPa

    # 降伏関数（単軸圧縮近似）
    # f = σ(1 - η/3) - k
    f_s = sigma_s * (1 - eta/3) - k

    data.append({
        'step': s,
        'F_y(s) [N]': F_y_s,
        'σ(s) [MPa]': sigma_s,
        'f(s) [MPa]': f_s
    })

# k*付近の表示
print(f"{'step':>6} {'F_y(s) [N]':>15} {'σ(s) [MPa]':>15} {'f(s) [MPa]':>15}")
print("-" * 60)
for d in data[k_star-3:k_star+3]:
    print(f"{d['step']:6d} {d['F_y(s) [N]']:15.3f} {d['σ(s) [MPa]']:15.6f} {d['f(s) [MPa]']:15.6e}")

# f≈0の確認
f_at_kstar = None
for d in data:
    if d['step'] == k_star:
        f_at_kstar = d['f(s) [MPa]']
        break

print(f"\nステップ {k_star} での f = {f_at_kstar:.3e} MPa")
if abs(f_at_kstar) < 1e-9:
    print("✓ |f| < 1e-9 を満たす")
else:
    tolerance = abs(f_at_kstar / k)
    print(f"  相対誤差: {tolerance:.3%}")

# ===============================
# グラフ用データをCSV出力（後処理用）
# ===============================
print("\n=== グラフ用データ出力 ===")
with open('dp_verification_graph_data.csv', 'w') as f:
    f.write("step,F_y_N,sigma_MPa,f_MPa\n")
    for d in data:
        f.write(f"{d['step']},{d['F_y(s) [N]']:.6f},{d['σ(s) [MPa]']:.6f},{d['f(s) [MPa]']:.6e}\n")
print("グラフデータ保存: dp_verification_graph_data.csv")

# ===============================
# 出力サマリー
# ===============================
print("\n" + "="*60)
print("計算結果サマリー")
print("="*60)
print(f"総ステップ数 N = {N}")
print(f"目標降伏ステップ k* = {k_star}")
print(f"Drucker-Prager: η = {eta:.6f}, k = {k:.3f} MPa")
print(f"圧縮降伏応力 σ_y^(c) = {sigma_y_c:.3f} MPa")
print(f"載荷面積 A = {A} mm^2")
print(f"基準節点力 F_y^0 = {F_y0} N")
print(f"荷重係数 λ_max = {lambda_max:.6f}")
print("="*60)

# ===============================
# CMLファイル更新
# ===============================
print("\n=== CMLファイル更新 ===")

input_file = 'input_files/1elem_f.cml'
backup_file = 'input_files/1elem_f.cml.bak'
output_file = 'input_files/1elem_f_kstar10.cml'

# バックアップ作成
shutil.copy(input_file, backup_file)
print(f"バックアップ作成: {backup_file}")

# ファイル読み込みと更新
with open(input_file, 'r') as f:
    lines = f.readlines()

# /LOADC/セクションのFactor of Total Loadsを更新
# 現在の構造では明示的な係数がないので、力の値を直接スケーリング
new_lines = []
loadc_section = False
for i, line in enumerate(lines):
    # /LOADC/セクションの節点力定義行を探す
    if '/LOADC/' in line:
        loadc_section = True
        new_lines.append(line)
        continue

    if loadc_section and line.strip() and line.strip()[0].isdigit():
        # 節点力定義行の処理 (固定フォーマット)
        # フォーマット: 節点番号(8文字) + 各方向力(各12文字×6)
        try:
            # 固定位置から値を読み取る
            node = line[0:8].strip()
            # 各力成分は12文字幅（符号込み）
            fx_str = line[8:20].strip()
            fy_str = line[20:32].strip()
            fz_str = line[32:44].strip()
            mx_str = line[44:56].strip()
            my_str = line[56:68].strip()
            mz_str = line[68:80].strip()

            fx = 0.0 if not fx_str else float(fx_str.replace('D', 'E'))
            fy = 0.0 if not fy_str else float(fy_str.replace('D', 'E'))
            fz = 0.0 if not fz_str else float(fz_str.replace('D', 'E'))
            mx = 0.0 if not mx_str else float(mx_str.replace('D', 'E'))
            my = 0.0 if not my_str else float(my_str.replace('D', 'E'))
            mz = 0.0 if not mz_str else float(mz_str.replace('D', 'E'))

            # y方向力にλ_maxを適用
            fy_new = fy * lambda_max

            # フォーマットを保持して出力（固定幅）
            new_line = f"{node:>8} {fx:+.5E}{fy_new:+.5E} {fz:+.5E} {mx:+.5E} {my:+.5E} {mz:+.5E}\n"
            new_lines.append(new_line)
        except (ValueError, IndexError) as e:
            # パース失敗時は元の行を保持
            new_lines.append(line)
    else:
        if loadc_section and '/SOLUT/' in line:
            loadc_section = False
        new_lines.append(line)

# 新ファイル書き込み
with open(output_file, 'w') as f:
    f.writelines(new_lines)

print(f"更新済みCMLファイル: {output_file}")
print(f"  λ_max = {lambda_max:.6f} を適用")

# 検証用CSV出力
with open('dp_verification_table.csv', 'w') as f:
    f.write("step,F_y_N,sigma_MPa,f_MPa\n")
    for d in data:
        f.write(f"{d['step']},{d['F_y(s) [N]']:.6f},{d['σ(s) [MPa]']:.6f},{d['f(s) [MPa]']:.6e}\n")
print(f"\n検証表をCSVファイルに保存: dp_verification_table.csv")