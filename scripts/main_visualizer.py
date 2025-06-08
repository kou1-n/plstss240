# main_visualizer.py
import data_loader # data_loader.py の関数をインポート
import plot_dis
import plot_ene
import plot_nor
import plot_sts
import plot_tmp
import pandas as pd

def get_user_choice(prompt, options):
    """ユーザーに選択肢を提示し、選択を取得する関数"""
    print(prompt)
    for i, option in enumerate(options):
        print(f"{i+1}. {option}")
    while True:
        try:
            choice = int(input("番号を選択してください: "))
            if 1 <= choice <= len(options):
                return choice
            else:
                print("無効な番号です。もう一度入力してください。")
        except ValueError:
            print("数値を入力してください。")

def get_column_info_from_user(df, prompt_message, num_indices=1, allow_skip=False):
    """
    ユーザーにデータフレームの列インデックスを選択させ、
    選択された (Step列を除いた0始まりの)インデックスと列名のリストのタプルを返す。
    [(idx, col_name), (idx, col_name), ...]
    """
    if df.empty:
        # print("データフレームが空です。列を選択できません。") # mainループで表示するのでここでは抑制
        return [] 

    print(prompt_message)
    cols_without_step_info = []
    if 'step' in df.columns.str.lower(): 
        cols_without_step_info = [(i, col) for i, col in enumerate(df.columns) if col.lower() != 'step']
        print("利用可能なデータ列 (Step列を除く、0始まりのインデックス):")
        cols_display_list = [info[1] for info in cols_without_step_info]
        for i, col_name in enumerate(cols_display_list):
            print(f"  {i}: {col_name}")
    else: 
        cols_without_step_info = [(i, col) for i, col in enumerate(df.columns)]
        cols_display_list = [info[1] for info in cols_without_step_info]
        print("利用可能なデータ列 (0始まりのインデックス、Step列が見つかりませんでした):")
        for i, col_name in enumerate(cols_display_list):
            print(f"  {i}: {col_name}")

    if not cols_display_list: 
        print("Step列以外の利用可能なデータ列がありません。")
        return []

    selected_columns_info = [] 
    for i in range(num_indices):
        while True:
            try:
                user_input_prompt = f"{i+1}番目の列のインデックスを入力してください"
                if allow_skip and i > 0:
                    user_input_prompt += " (スキップする場合はEnterキーのみ): "
                else:
                    user_input_prompt += ": "
                user_input = input(user_input_prompt)
                if allow_skip and i > 0 and not user_input: 
                    return selected_columns_info 
                idx_val = int(user_input)
                if 0 <= idx_val < len(cols_display_list):
                    selected_columns_info.append((idx_val, cols_display_list[idx_val]))
                    break 
                else:
                    print(f"無効なインデックスです。0 から {len(cols_display_list)-1} の間で入力してください。")
            except ValueError:
                if allow_skip and i > 0 and not user_input: 
                    return selected_columns_info
                print("数値を入力してください。")
    return selected_columns_info


def main():
    case_prefix_input = input("解析ケースのプレフィックスを入力してください (例: cube250514, デフォルトは cube250514): ")
    case_prefix = case_prefix_input if case_prefix_input else "cube250514"
    print(f"解析ケース: '{case_prefix}'")

    base_path = "./"
    cml_file_path = f"{base_path}{case_prefix}.cml" # CMLファイルパス

    # --- 材料情報の読み込み ---
    material_info_text = data_loader.parse_cml_material_section(cml_file_path)
    print(f"\n--- 材料情報 ({cml_file_path}) ---")
    print(material_info_text)
    print("---------------------------------------\n")


    file_paths = {
        "DIS": f"{base_path}DIS_{case_prefix}.txt",
        "ENE": f"{base_path}ENE_{case_prefix}.txt",
        "NOR": f"{base_path}NOR_{case_prefix}.txt",
        "STS": f"{base_path}STS_{case_prefix}.txt",
        "TMP": f"{base_path}TMP_{case_prefix}.txt",
    }

    print(f"\n--- {case_prefix} のデータファイルを読み込み中 ---")
    dataframes = {}
    for key, path in file_paths.items():
        # print(f"'{path}' をロード中...") # data_loader内でエラー表示するので抑制も可
        df = data_loader.parse_header_and_load(path)
        dataframes[key] = df
        if not df.empty:
            # print(f"  {key}データサンプル:\n{df.head()}")
            # print(f"  {key} Columns: {df.columns.tolist()}")
            pass
        else:
            print(f"  警告: {key}データは空、または読み込めませんでした。 ({path})")
        # print("-" * 20) # ログが冗長ならコメントアウト
    
    print("---------------------------------------\n")

    while True:
        plot_choice = get_user_choice(
            "どのデータを可視化しますか？",
            [
                "変位/荷重データ (DIS)",
                "エネルギーデータ (ENE)",
                "収束履歴データ (NOR)",
                "応力/ひずみデータ (STS)",
                "温度データ (TMP)",
                "終了"
            ]
        )

        if plot_choice == 1: # DIS
            df_dis = dataframes.get("DIS")
            if df_dis is None or df_dis.empty:
                print("DISデータが読み込まれていないか空です。")
                continue
            dis_plot_type = get_user_choice(
                "DISデータのプロットタイプを選択してください:",
                ["荷重-変位曲線", "変位 vs. ステップ", "戻る"]
            )
            if dis_plot_type == 1:
                selected_info = get_column_info_from_user(df_dis, "荷重-変位曲線のために、変位と荷重の列をそれぞれ選択:", 2)
                if len(selected_info) == 2:
                    disp_idx, disp_name = selected_info[0]
                    force_idx, force_name = selected_info[1]
                    title = f"{case_prefix}: Force ({force_name}) vs. Displacement ({disp_name})"
                    plot_dis.plot_force_displacement(df_dis, disp_idx, force_idx, title, material_info_text)
            elif dis_plot_type == 2:
                selected_info = get_column_info_from_user(df_dis, "変位 vs. ステップのために、プロットしたい変位の列を選択 (最大2つ):", 2, allow_skip=True)
                if selected_info: 
                    indices = [info[0] for info in selected_info]
                    names = [info[1] for info in selected_info]
                    title = f"{case_prefix}: Displacement ({', '.join(names)}) vs. Step"
                    plot_dis.plot_displacement_vs_step(df_dis, indices, title, material_info_text)
            elif dis_plot_type == 3:
                continue

        elif plot_choice == 2: # ENE
            df_ene = dataframes.get("ENE")
            if df_ene is None or df_ene.empty:
                print("ENEデータが読み込まれていないか空です。")
                continue
            title = f"{case_prefix}: Energy vs. Step"
            plot_ene.plot_energy_vs_step(df_ene, title, material_info_text)

        elif plot_choice == 3: # NOR
            df_nor = dataframes.get("NOR")
            if df_nor is None or df_nor.empty:
                print("NORデータが読み込まれていないか空です。")
                continue
            if 'Step' in df_nor.columns and not df_nor['Step'].dropna().empty:
                valid_steps = pd.to_numeric(df_nor['Step'], errors='coerce').dropna().unique()
                if valid_steps.size > 0:
                    min_step, max_step = int(valid_steps.min()), int(valid_steps.max())
                    while True:
                        try:
                            target_step_str = input(f"収束履歴を表示するステップ番号を入力してください ({min_step}-{max_step}): ")
                            target_step = int(target_step_str)
                            if target_step in valid_steps: 
                                title = f"{case_prefix}: Convergence History (Step {target_step})"
                                plot_nor.plot_convergence_history(df_nor, target_step, title, material_info_text) 
                                break
                            else:
                                print(f"無効なステップ番号、またはそのステップのデータがありません。利用可能なステップ: {sorted([int(s) for s in valid_steps])}")
                        except ValueError:
                            print("数値を入力してください。")
                else:
                    print("NORデータに有効なStep情報がありません。")
            else:
                print("NORデータにStep列がないか、有効なStep情報がありません。")

        elif plot_choice == 4: # STS
            df_sts = dataframes.get("STS")
            if df_sts is None or df_sts.empty:
                print("STSデータが読み込まれていないか空です。")
                continue
            sts_plot_type = get_user_choice(
                "STSデータのプロットタイプを選択してください:",
                ["応力-ひずみ曲線", "応力 vs. ステップ", "戻る"]
            )
            if sts_plot_type == 1:
                selected_info = get_column_info_from_user(df_sts, "応力-ひずみ曲線のために、ひずみと応力の列をそれぞれ選択:", 2)
                if len(selected_info) == 2:
                    strain_idx, strain_name = selected_info[0]
                    stress_idx, stress_name = selected_info[1]
                    title = f"{case_prefix}: Stress ({stress_name}) vs. Strain ({strain_name})"
                    plot_sts.plot_stress_strain_curve(df_sts, strain_idx, stress_idx, title, material_info_text)
            elif sts_plot_type == 2:
                selected_info = get_column_info_from_user(df_sts, "応力 vs. ステップのために、プロットしたい応力の列を選択 (最大2つ):", 2, allow_skip=True)
                if selected_info:
                    indices = [info[0] for info in selected_info]
                    names = [info[1] for info in selected_info]
                    title = f"{case_prefix}: Stress ({', '.join(names)}) vs. Step"
                    plot_sts.plot_stress_vs_step(df_sts, indices, title, material_info_text)
            elif sts_plot_type == 3:
                continue

        elif plot_choice == 5: # TMP
            df_tmp = dataframes.get("TMP")
            if df_tmp is None or df_tmp.empty:
                print("TMPデータが読み込まれていないか空です。")
                continue
            cols_without_step_tmp = [col for col in df_tmp.columns if col.lower() != 'step']
            if not cols_without_step_tmp:
                print("TMPファイルにはステップ番号以外のデータ列がありません。")
                continue
            selected_info = get_column_info_from_user(df_tmp, "温度 vs. ステップのために、プロットしたい温度の列を選択:", 1)
            if selected_info: 
                temp_idx, temp_name = selected_info[0]
                title = f"{case_prefix}: Temperature ({temp_name}) vs. Step"
                plot_tmp.plot_temperature_vs_step(df_tmp, temp_idx, title, material_info_text)
        
        elif plot_choice == 6: # 終了
            print("プログラムを終了します。")
            break
        
        print("-" * 30)

if __name__ == "__main__":
    main()
