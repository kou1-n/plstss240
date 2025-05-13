# plot_ene.py
import matplotlib.pyplot as plt
import pandas as pd

def plot_energy_vs_step(df_ene, fig_title='Energy vs. Step', material_info=""):
    """
    各種エネルギーのステップ変化をプロットする。
    """
    if df_ene.empty or 'Step' not in df_ene.columns:
        print("ENEデータが空か、Step列が存在しません。")
        return

    fig, ax = plt.subplots(figsize=(10, 6)) # fig, ax を取得
    plotted = False
    
    energy_cols_candidates = ['TotalEne.', 'Elastic Ene.', 'Plastic Ene.']
    energy_cols_to_plot = [col for col in energy_cols_candidates if col in df_ene.columns]
    
    if not energy_cols_to_plot:
        # print("ENEプロット情報: 標準的なエネルギー列名が見つかりませんでした。Step以外の数値列をプロットします。")
        energy_cols_to_plot = [col for col in df_ene.columns if col.lower() != 'step' and pd.api.types.is_numeric_dtype(df_ene[col])]

    if not energy_cols_to_plot:
        print("ENEプロットエラー: プロット可能な数値データ列が見つかりませんでした。")
        # print(f"利用可能な列: {df_ene.columns.tolist()}")
        plt.close(fig) # プロットできない場合は図を閉じる
        return

    for col_name in energy_cols_to_plot:
        if col_name in df_ene.columns and pd.api.types.is_numeric_dtype(df_ene[col_name]):
            ax.plot(df_ene['Step'], df_ene[col_name], marker='.', linestyle='-', label=col_name)
            plotted = True
            
    if plotted:
        ax.set_xlabel('Step')
        ax.set_ylabel('Energy')
        ax.set_title(fig_title)
        ax.legend()
        ax.grid(True)

        if material_info:
            plt.figtext(0.5, 0.01, material_info, ha="center", va="bottom", fontsize=8,
                        bbox={"facecolor":"white", "alpha":0.7, "pad":3, "edgecolor":"none"})
        plt.subplots_adjust(bottom=0.15 if material_info else 0.1) # 情報表示スペースを確保
        plt.show()
    else:
        # print("ENEプロットエラー: プロット可能なエネルギー列が見つかりませんでした。")
        plt.close(fig) # プロットできない場合は図を閉じる
