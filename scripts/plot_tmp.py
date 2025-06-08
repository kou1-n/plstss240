# plot_tmp.py
import matplotlib.pyplot as plt
import pandas as pd

def plot_temperature_vs_step(df_tmp, temp_col_idx, fig_title='Temperature vs. Step', material_info=""):
    """
    温度のステップ変化をプロットする。
    temp_col_idx は Step列を除いたDataFrameの0始まりのインデックス
    """
    if df_tmp.empty or 'Step' not in df_tmp.columns:
        # print("TMPデータが空か、Step列が存在しません。")
        return
    
    cols_without_step = [col for col in df_tmp.columns if col.lower() != 'step']
    
    if not (0 <= temp_col_idx < len(cols_without_step)):
        print("TMPプロットエラー: 指定された温度列のインデックスが無効です。")
        # ... (エラーメッセージの詳細は省略)
        return

    temp_col_name = cols_without_step[temp_col_idx]

    if not pd.api.types.is_numeric_dtype(df_tmp[temp_col_name]):
        print(f"TMPプロットエラー: 温度列 '{temp_col_name}' が数値データではありません。")
        return

    fig, ax = plt.subplots(figsize=(8, 6)) # fig, ax を取得
    ax.plot(df_tmp['Step'], df_tmp[temp_col_name], marker='o', linestyle='-')
    ax.set_xlabel('Step')
    ax.set_ylabel(f'Temperature ({temp_col_name})')
    ax.set_title(fig_title)
    ax.grid(True)

    if material_info:
        plt.figtext(0.5, 0.01, material_info, ha="center", va="bottom", fontsize=8,
                    bbox={"facecolor":"white", "alpha":0.7, "pad":3, "edgecolor":"none"})
    plt.subplots_adjust(bottom=0.15 if material_info else 0.1)
    plt.show()
