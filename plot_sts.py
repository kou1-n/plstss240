# plot_sts.py
import matplotlib.pyplot as plt
import pandas as pd

def plot_stress_strain_curve(df_sts, strain_col_idx, stress_col_idx, fig_title='Stress-Strain Curve', material_info=""):
    """
    応力-ひずみ曲線をプロットする。
    strain_col_idx, stress_col_idx は Step列を除いたDataFrameの0始まりのインデックス
    """
    if df_sts.empty:
        # print("STSデータが空です。")
        return
    
    cols_without_step = [col for col in df_sts.columns if col.lower() != 'step']

    if not (0 <= strain_col_idx < len(cols_without_step) and \
            0 <= stress_col_idx < len(cols_without_step)):
        print("STSプロットエラー: 指定されたひずみまたは応力の列インデックスが無効です。")
        # ... (エラーメッセージの詳細は省略)
        return

    strain_col_name = cols_without_step[strain_col_idx]
    stress_col_name = cols_without_step[stress_col_idx]

    if not (pd.api.types.is_numeric_dtype(df_sts[strain_col_name]) and \
            pd.api.types.is_numeric_dtype(df_sts[stress_col_name])):
        print(f"STSプロットエラー: ひずみ列 '{strain_col_name}' または応力列 '{stress_col_name}' が数値データではありません。")
        return

    fig, ax = plt.subplots(figsize=(8, 6)) # fig, ax を取得
    ax.plot(df_sts[strain_col_name], df_sts[stress_col_name], marker='o', linestyle='-')
    ax.set_xlabel(f'Strain ({strain_col_name})')
    ax.set_ylabel(f'Stress ({stress_col_name})')
    ax.set_title(fig_title)
    ax.grid(True)

    if material_info:
        plt.figtext(0.5, 0.01, material_info, ha="center", va="bottom", fontsize=8,
                    bbox={"facecolor":"white", "alpha":0.7, "pad":3, "edgecolor":"none"})
    plt.subplots_adjust(bottom=0.15 if material_info else 0.1)
    plt.show()

def plot_stress_vs_step(df_sts, stress_col_indices, fig_title='Stress vs. Step', material_info=""):
    """指定された応力成分のステップ変化をプロットする。"""
    if df_sts.empty or 'Step' not in df_sts.columns:
        # print("STSデータが空か、Step列が存在しません。")
        return

    cols_without_step = [col for col in df_sts.columns if col.lower() != 'step']

    fig, ax = plt.subplots(figsize=(10, 6)) # fig, ax を取得
    plotted = False
    for idx in stress_col_indices:
        if 0 <= idx < len(cols_without_step):
            col_name = cols_without_step[idx]
            if pd.api.types.is_numeric_dtype(df_sts[col_name]):
                ax.plot(df_sts['Step'], df_sts[col_name], marker='.', linestyle='-', label=col_name)
                plotted = True
            # else:
                # print(f"STSプロット警告: 列 '{col_name}' は数値データではないためスキップされました。")
        # else:
            # print(f"STSプロット警告: 無効なインデックス {idx} が指定されたためスキップされました。")
            
    if plotted:
        ax.set_xlabel('Step')
        ax.set_ylabel('Stress')
        ax.set_title(fig_title)
        ax.legend()
        ax.grid(True)

        if material_info:
            plt.figtext(0.5, 0.01, material_info, ha="center", va="bottom", fontsize=8,
                        bbox={"facecolor":"white", "alpha":0.7, "pad":3, "edgecolor":"none"})
        plt.subplots_adjust(bottom=0.15 if material_info else 0.1)
        plt.show()
    else:
        # print("プロット可能な応力列が見つからなかったか、すべて数値データではありませんでした。")
        # if cols_without_step:
            # print(f"利用可能なStep列を除いた列名: {cols_without_step}")
        plt.close(fig) # プロットしない場合は図を閉じる
