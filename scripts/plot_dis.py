# plot_dis.py
import matplotlib.pyplot as plt
import pandas as pd

def plot_force_displacement(df_dis, disp_col_idx, force_col_idx, fig_title='Force-Displacement Curve', material_info=""): # material_info 引数を追加
    """
    荷重-変位曲線をプロットする。
    disp_col_idx, force_col_idx は Step列を除いたDataFrameの0始まりのインデックス
    """
    if df_dis.empty:
        # print("DISデータが空です。") # mainで表示するので抑制
        return

    cols_without_step = [col for col in df_dis.columns if col.lower() != 'step']
    
    if not (0 <= disp_col_idx < len(cols_without_step) and \
            0 <= force_col_idx < len(cols_without_step)):
        print("DISプロットエラー: 指定された変位または荷重の列インデックスが無効です。")
        # ... (エラーメッセージの詳細は省略)
        return

    disp_col_name = cols_without_step[disp_col_idx]
    force_col_name = cols_without_step[force_col_idx]

    if not (pd.api.types.is_numeric_dtype(df_dis[disp_col_name]) and \
            pd.api.types.is_numeric_dtype(df_dis[force_col_name])):
        print(f"DISプロットエラー: 変位列 '{disp_col_name}' または荷重列 '{force_col_name}' が数値データではありません。")
        return

    fig, ax = plt.subplots(figsize=(8, 6)) # fig オブジェクトを取得
    ax.plot(df_dis[disp_col_name], df_dis[force_col_name], marker='o', linestyle='-')
    ax.set_xlabel(f'Displacement ({disp_col_name})')
    ax.set_ylabel(f'Force ({force_col_name})')
    ax.set_title(fig_title)
    ax.grid(True)

    # 材料情報をグラフの下部に表示
    if material_info:
        plt.figtext(0.5, 0.01, material_info, ha="center", va="bottom", fontsize=8,
                    bbox={"facecolor":"white", "alpha":0.7, "pad":3, "edgecolor":"none"})
    
    plt.subplots_adjust(bottom=0.15 if material_info else 0.1) # 情報表示スペースを確保
    plt.show()

def plot_displacement_vs_step(df_dis, disp_col_indices, fig_title='Displacement vs. Step', material_info=""): # material_info 引数を追加
    """指定された変位成分のステップ変化をプロットする。"""
    if df_dis.empty or 'Step' not in df_dis.columns:
        # print("DISデータが空か、Step列が存在しません。")
        return

    cols_without_step = [col for col in df_dis.columns if col.lower() != 'step']

    fig, ax = plt.subplots(figsize=(10, 6)) # fig オブジェクトを取得
    plotted = False
    for idx in disp_col_indices:
        if 0 <= idx < len(cols_without_step):
            col_name = cols_without_step[idx]
            if pd.api.types.is_numeric_dtype(df_dis[col_name]):
                ax.plot(df_dis['Step'], df_dis[col_name], marker='.', linestyle='-', label=col_name)
                plotted = True
            # else:
                # print(f"DISプロット警告: 列 '{col_name}' は数値データではないためスキップされました。")
        # else:
            # print(f"DISプロット警告: 無効なインデックス {idx} が指定されたためスキップされました。")
            
    if plotted:
        ax.set_xlabel('Step')
        ax.set_ylabel('Displacement')
        ax.set_title(fig_title)
        ax.legend()
        ax.grid(True)

        if material_info:
            plt.figtext(0.5, 0.01, material_info, ha="center", va="bottom", fontsize=8,
                        bbox={"facecolor":"white", "alpha":0.7, "pad":3, "edgecolor":"none"})
        plt.subplots_adjust(bottom=0.15 if material_info else 0.1)
        plt.show()
    else:
        # print("プロット可能な変位列が見つからなかったか、すべて数値データではありませんでした。")
        # if cols_without_step:
            # print(f"利用可能なStep列を除いた列名: {cols_without_step}")
        pass # エラーはmain側で表示されるのでここでは抑制も可
