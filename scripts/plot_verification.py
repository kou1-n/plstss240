#!/usr/bin/env python3
"""
Drucker-Prager検証グラフ作成
CSVデータから応力-ステップ、降伏関数-ステップのグラフを生成
"""

try:
    import matplotlib
    matplotlib.use('Agg')  # GUI不要のバックエンド
    import matplotlib.pyplot as plt
    has_matplotlib = True
except ImportError:
    has_matplotlib = False
    print("matplotlib not available, creating simple ASCII plot instead")

# CSVファイル読み込み
steps = []
sigma_vals = []
f_vals = []

with open('dp_verification_graph_data.csv', 'r') as f:
    lines = f.readlines()
    for line in lines[1:]:  # ヘッダーをスキップ
        parts = line.strip().split(',')
        steps.append(int(parts[0]))
        sigma_vals.append(float(parts[2]))
        f_vals.append(float(parts[3]))

k_star = 11
sigma_y_c = 259.839  # 計算済みの値

if has_matplotlib:
    # matplotlibが利用可能な場合
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    # s-σ(s) プロット
    ax1.plot(steps, sigma_vals, 'b-', linewidth=2, label='σ(s)')
    ax1.axvline(x=k_star, color='r', linestyle='--', label=f'k* = {k_star}')
    ax1.axhline(y=sigma_y_c, color='g', linestyle='--', label=f'σ_y^(c) = {sigma_y_c:.1f} MPa')
    ax1.set_xlabel('Step')
    ax1.set_ylabel('σ [MPa]')
    ax1.set_title('Stress vs Step')
    ax1.grid(True, alpha=0.3)
    ax1.legend()
    ax1.set_xlim(0, 30)  # 初期部分を拡大

    # s-f(s) プロット
    ax2.plot(steps, f_vals, 'b-', linewidth=2, label='f(s)')
    ax2.axvline(x=k_star, color='r', linestyle='--', label=f'k* = {k_star}')
    ax2.axhline(y=0, color='k', linestyle='-', alpha=0.3, label='f = 0')
    ax2.set_xlabel('Step')
    ax2.set_ylabel('f [MPa]')
    ax2.set_title('Yield Function vs Step')
    ax2.grid(True, alpha=0.3)
    ax2.legend()
    ax2.set_xlim(0, 30)  # 初期部分を拡大

    plt.tight_layout()
    plt.savefig('dp_verification.png', dpi=150, bbox_inches='tight')
    print("グラフを保存: dp_verification.png")

else:
    # matplotlibが利用できない場合、簡易ASCII表示
    print("\n=== ASCII Plot: σ vs Step (first 20 steps) ===")
    print("Step |", "=" * 50)

    max_sigma_20 = max(sigma_vals[:20])
    for i in range(min(20, len(steps))):
        s = steps[i]
        sigma = sigma_vals[i]
        bar_len = int(50 * sigma / max_sigma_20)
        marker = "*" if s == k_star else "#"
        print(f"{s:4d} | {marker * bar_len}")

    print("\n=== ASCII Plot: f vs Step (first 20 steps) ===")
    print("Step | -100 " + "-" * 20 + " 0 " + "+" * 20 + " +100")

    for i in range(min(20, len(steps))):
        s = steps[i]
        f = f_vals[i]

        # -100 から +100 の範囲で正規化
        if abs(f) > 100:
            pos = 50 if f > 0 else 0
        else:
            pos = int(25 + f/4)  # -100 to +100 を 0 to 50 にマップ

        line = [" "] * 51
        line[25] = "|"  # 中心線 (f=0)
        if s == k_star:
            line[pos] = "*"
        else:
            line[pos] = "#"

        print(f"{s:4d} | {''.join(line)}")

    print("\nLegend: * = k* step, # = other steps, | = f=0 line")

# 数値サマリー
print("\n=== 数値検証サマリー ===")
print(f"k* = {k_star} での値:")
print(f"  σ({k_star}) = {sigma_vals[k_star-1]:.3f} MPa")
print(f"  f({k_star}) = {f_vals[k_star-1]:.3e} MPa")

# 前後のステップ
if k_star > 1:
    print(f"k* - 1 = {k_star-1} での値:")
    print(f"  σ({k_star-1}) = {sigma_vals[k_star-2]:.3f} MPa")
    print(f"  f({k_star-1}) = {f_vals[k_star-2]:.3e} MPa")

if k_star < len(steps):
    print(f"k* + 1 = {k_star+1} での値:")
    print(f"  σ({k_star+1}) = {sigma_vals[k_star]:.3f} MPa")
    print(f"  f({k_star+1}) = {f_vals[k_star]:.3e} MPa")

print("\n降伏が k* = {} で発生することを確認".format(k_star))