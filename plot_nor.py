# plot_nor.py
import matplotlib.pyplot as plt
import pandas as pd

def plot_convergence_history(df_nor, target_step=1, fig_title_prefix='Convergence History', material_info=""):
    """
    指定されたステップの収束履歴をプロットする。
    fig_title_prefix は main_visualizer.py から完全なタイトルとして渡されることを想定。
    """
    if df_nor.empty or not all(col in df_nor.columns for col in ['Step', 'Iter']):
        print("NORデータが空か、Step/Iter列が存在しません。")
        return

    # Iter 列と Step 列を数値型に変換 (エラーはNaNに)
    df_nor['Iter'] = pd.to_numeric(df_nor['Iter'], errors='coerce')
    df_nor['Step'] = pd.to_numeric(df_nor['Step'], errors='coerce')
    
    # NaNを含む行を除外してからStepを整数型に (プロット前にNaNは扱いにくいため)
    df_nor_copy = df_nor.dropna(subset=['Iter', 'Step']).copy()
    if df_nor_copy.empty:
        print(f"NORデータに有効なStep/Iter情報を持つ行がありません (ステップ {target_step} 確認前)。")
        return
    df_nor_copy['Step'] = df_nor_copy['Step'].astype(int)


    df_step_nor = df_nor_copy[df_nor_copy['Step'] == target_step]
    
    if df_step_nor.empty:
        print(f"NORデータにステップ {target_step} の有効なデータが存在しません。")
        return

    fig, ax = plt.subplots(figsize=(10, 6)) # fig, ax を取得
    plotted = False
    
    norm_cols_candidates = ['unorm', 'tnorm']
    norm_cols_to_plot = [col for col in norm_cols_candidates if col in df_step_nor.columns]

    if not norm_cols_to_plot:
        # print(f"NORプロット情報 (Step {target_step}): 標準的なノルム列名が見つかりませんでした。Step/Iter以外の数値列をプロットします。")
        potential_cols = [col for col in df_step_nor.columns if col.lower() not in ['step', 'iter']]
        for pc in potential_cols:
            if pd.api.types.is_numeric_dtype(df_step_nor[pc]):
                norm_cols_to_plot.append(pc)
    
    if not norm_cols_to_plot:
        print(f"NORプロットエラー: ステップ {target_step} でプロット可能な数値データ列が見つかりませんでした。")
        # print(f"利用可能な列 (Step {target_step}): {df_step_nor.columns.tolist()}")
        plt.close(fig)
        return

    for col_name in norm_cols_to_plot:
        if col_name in df_step_nor.columns and pd.api.types.is_numeric_dtype(df_step_nor[col_name]):
            ax.plot(df_step_nor['Iter'], df_step_nor[col_name], marker='.', linestyle='-', label=col_name)
            plotted = True

    if plotted:
        ax.set_xlabel('Iteration')
        ax.set_ylabel('Norm Value')
        ax.set_title(fig_title_prefix) # mainから渡された完全なタイトルを使用
        
        can_log_scale = True
        for col_name in norm_cols_to_plot:
            if col_name in df_step_nor.columns and pd.api.types.is_numeric_dtype(df_step_nor[col_name]):
                # NaNを除外して正の値のみか確認
                if not (df_step_nor[col_name].dropna() > 0).all() and not df_step_nor[col_name].dropna().empty:
                    can_log_scale = False
                    break
        
        # プロットされたデータが存在し、かつ最大値が0より大きい場合にlogスケールを試みる
        if can_log_scale and any(df_step_nor[col].max() > 0 for col in norm_cols_to_plot if col in df_step_nor and pd.api.types.is_numeric_dtype(df_step_nor[col]) and not df_step_nor[col].dropna().empty):
             ax.set_yscale('log')
        
        ax.legend()
        ax.grid(True)

        if material_info:
            plt.figtext(0.5, 0.01, material_info, ha="center", va="bottom", fontsize=8,
                        bbox={"facecolor":"white", "alpha":0.7, "pad":3, "edgecolor":"none"})
        plt.subplots_adjust(bottom=0.15 if material_info else 0.1)
        plt.show()
    else:
        # print(f"NORプロットエラー: ステップ {target_step} でプロット可能なノルム列が見つかりませんでした。")
        plt.close(fig)
