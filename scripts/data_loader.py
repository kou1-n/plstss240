# data_loader.py
import pandas as pd
import re # 正規表現モジュールをインポート

def parse_header_and_load(filepath):
    """
    ファイルのコメント行からヘッダー情報を解析し、データをDataFrameとして読み込む。
    ヘッダーが不完全な場合（例: TMPファイル）にも対応を試みる。
    """
    extracted_column_names = []
    comment_lines_count = 0
    header_line_content = "" # 最後に解析された有効なヘッダー行候補

    try:
        with open(filepath, 'r', encoding='utf-8') as f: # encodingを指定
            for i, line in enumerate(f):
                stripped_line = line.strip()
                if stripped_line.startswith('%'):
                    comment_lines_count = i + 1
                    potential_header = stripped_line.lstrip('%').strip()
                    if potential_header: # 空でないコメント行をヘッダー候補として保持
                        header_line_content = potential_header
                elif header_line_content: # ヘッダー候補が見つかった後の最初の非コメント行で終了
                    break
                elif not stripped_line: # 空行はスキップ
                    comment_lines_count = i + 1 # コメント行としてカウントアップ
                    continue
                else: # コメントでもヘッダーでもない最初のデータ行（ヘッダーなしと判断）
                    break
        
        if header_line_content:
            extracted_column_names = [name for name in header_line_content.split(' ') if name]

        # --- データ部分をまずヘッダーなしで読み込み、実際の列数を確認 ---
        # encoding を指定
        df = pd.read_csv(filepath, delim_whitespace=True, skiprows=comment_lines_count, header=None, comment='%', escapechar='\\', encoding='utf-8') 
        
        if df.empty:
            # print(f"情報 ({filepath}): ファイルが空か、データ部分が読み込めませんでした。")
            return pd.DataFrame()

        actual_num_columns = df.shape[1]
        final_column_names = []

        if extracted_column_names: 
            if len(extracted_column_names) == actual_num_columns:
                final_column_names = extracted_column_names
            elif len(extracted_column_names) < actual_num_columns:
                # print(f"情報 ({filepath}): ヘッダー列名数({len(extracted_column_names)})が実データ列数({actual_num_columns})より少ないため、不足分を補完します。")
                final_column_names = extracted_column_names[:] 
                for i in range(len(extracted_column_names), actual_num_columns):
                    if len(extracted_column_names) == 1 and extracted_column_names[0].lower() == 'step' and i == 1:
                        final_column_names.append('Temperature_1') 
                    else:
                        final_column_names.append(f'data_col_{i - len(extracted_column_names) +1}')
            else: 
                # print(f"警告 ({filepath}): ヘッダー列名数({len(extracted_column_names)})が実データ列数({actual_num_columns})より多いため、実データ列数に合わせます。")
                final_column_names = extracted_column_names[:actual_num_columns]
        elif actual_num_columns > 0:
            # print(f"情報 ({filepath}): ヘッダーが抽出できませんでした。実際のデータ列数 ({actual_num_columns}個) に基づいて仮の列名を付与します。")
            if actual_num_columns == 1:
                 final_column_names = ['Step'] 
            elif actual_num_columns > 1:
                final_column_names = ['Step'] + [f'data_col_{i}' for i in range(1, actual_num_columns)]
            else: 
                 return pd.DataFrame()
        elif not extracted_column_names and actual_num_columns == 0 : 
             return pd.DataFrame()


        # 本読み込み (列名を指定して再度読み込む)
        if final_column_names:
             # encoding を指定
            df = pd.read_csv(filepath, delim_whitespace=True, skiprows=comment_lines_count, names=final_column_names, comment='%', escapechar='\\', encoding='utf-8')
        else: # これでも列名が決まらない場合
             # encoding を指定
            df = pd.read_csv(filepath, delim_whitespace=True, skiprows=comment_lines_count, header=None, comment='%', escapechar='\\', encoding='utf-8')
            if not df.empty:
                df.columns = [f'col_{j}' for j in range(df.shape[1])] 


        if not df.empty:
            if 'Step' in df.columns:
                df['Step'] = pd.to_numeric(df['Step'], errors='coerce')
                df.dropna(subset=['Step'], inplace=True) 
            elif df.shape[1] > 0 and 'col_0' in df.columns : 
                df.rename(columns={'col_0': 'Step'}, inplace=True)
                df['Step'] = pd.to_numeric(df['Step'], errors='coerce')
                df.dropna(subset=['Step'], inplace=True)
            
            for col in df.columns:
                if col.lower() != 'step': 
                    df[col] = pd.to_numeric(df[col], errors='ignore')
        
        return df

    except pd.errors.EmptyDataError:
        # print(f"ファイル {filepath} は空か、読み込めるデータがありません。")
        return pd.DataFrame()
    except Exception as e:
        # print(f"ファイル {filepath} の読み込み中に予期せぬエラーが発生しました: {e}")
        try:
            # encoding を指定
            df = pd.read_csv(filepath, delim_whitespace=True, comment='%', header=None, escapechar='\\', encoding='utf-8')
            if not df.empty:
                first_data_row = 0
                for k_idx in range(len(df)):
                    row_series = df.iloc[k_idx]
                    if any(pd.to_numeric(s, errors='coerce') is not pd.NA for s in row_series) or \
                       any(isinstance(s, (int, float)) for s in row_series):
                        first_data_row = k_idx
                        break
                else: 
                     return pd.DataFrame()
                # encoding を指定
                df = pd.read_csv(filepath, delim_whitespace=True, skiprows=first_data_row, header=None, comment='%', escapechar='\\', encoding='utf-8')
                if not df.empty:
                    if df.shape[1] == 1:
                        df.columns = ['Step']
                    elif df.shape[1] > 1:
                        df.columns = ['Step'] + [f'data_col_{j}' for j in range(1, df.shape[1])]
                    else:
                        return pd.DataFrame()

                    if 'Step' in df.columns:
                        df['Step'] = pd.to_numeric(df['Step'], errors='coerce')
                        df.dropna(subset=['Step'], inplace=True)
                    for col in df.columns:
                        if col.lower() != 'step':
                            df[col] = pd.to_numeric(df[col], errors='ignore')
                return df
            else:
                return pd.DataFrame()
        except Exception as e_fallback:
            # print(f"最終フォールバック読み込みも失敗しました ({filepath}): {e_fallback}")
            return pd.DataFrame()

def parse_cml_material_section(cml_filepath):
    """
    CMLファイルから最初の /MATER/ セクションの主要な物性値を読み取り、
    整形された文字列として返す。
    """
    material_params = {}
    try:
        with open(cml_filepath, 'r', encoding='utf-8') as f: # encodingを指定
            in_mater_section = False
            mater_lines_count = 0
            num_materials_to_read = 0
            materials_read = 0

            for line in f:
                stripped_line = line.strip()
                if stripped_line.startswith('/MATER/'):
                    in_mater_section = True
                    mater_lines_count = 0
                    continue

                if in_mater_section:
                    mater_lines_count += 1
                    if mater_lines_count == 1: # Number of material sets
                        try:
                            num_materials_to_read = int(stripped_line.split()[0])
                        except (ValueError, IndexError):
                            # print(f"警告: /MATER/ の材料セット数を読み取れませんでした。: {stripped_line}")
                            return "Material info: Could not parse num_materials."
                        continue
                    
                    if materials_read < num_materials_to_read:
                        if mater_lines_count == 2: # Mat_ID, Model_ID (最初の材料セットのみ対象)
                            try:
                                parts = stripped_line.split()
                                material_params['MatID'] = parts[0]
                                material_params['Model'] = parts[1]
                            except IndexError:
                                # print(f"警告: 材料ID/モデルを読み取れませんでした。: {stripped_line}")
                                pass # 続行して物性値の読み取りを試みる
                        
                        # CMLフォーマットに基づき、特定の行から値を抽出
                        # (Isotropic linear elastic-plastic material model 1 を想定)
                        elif mater_lines_count == 3: # E, nu, rho, alpha, K_therm
                            values = stripped_line.split()
                            if len(values) >= 2:
                                material_params['E'] = values[0]
                                material_params['nu'] = values[1]
                        elif mater_lines_count == 4: # C_p, f_i, q1, q2, sigma_Y
                            values = stripped_line.split()
                            if len(values) >= 5:
                                material_params['sigma_Y'] = values[4]
                        elif mater_lines_count == 5: # H, s_inf, delta, n, f_F
                            values = stripped_line.split()
                            if len(values) >= 3:
                                material_params['H'] = values[0]
                                material_params['s_inf'] = values[1] # sigma_Y_infinity
                                material_params['delta'] = values[2]
                            # 1材料セット分の読み込みが完了したら終了 (今回は最初の材料のみ)
                            materials_read += 1
                            if materials_read >= num_materials_to_read: # 通常は1
                                in_mater_section = False # 次の/MATER/があっても読まない
                                break 
                    else: # 必要な材料セット数を読み終えたらセクション終了
                        in_mater_section = False
                        break
                
                if stripped_line.startswith('/') and not stripped_line.startswith('/MATER/'):
                    # 他のセクションヘッダーが見つかったら、/MATER/セクションは終了とみなす
                    if in_mater_section:
                        in_mater_section = False
                        break
        
        if not material_params:
            return "Material info: /MATER/ section not found or empty."

        # 整形された文字列を作成
        info_parts = []
        if 'MatID' in material_params: info_parts.append(f"MatID:{material_params['MatID']}")
        if 'Model' in material_params: info_parts.append(f"Model:{material_params['Model']}")
        if 'E' in material_params: info_parts.append(f"E:{material_params['E']}")
        if 'nu' in material_params: info_parts.append(f"nu:{material_params['nu']}")
        if 'sigma_Y' in material_params: info_parts.append(f"sigma_Y:{material_params['sigma_Y']}")
        if 'H' in material_params: info_parts.append(f"H:{material_params['H']}")
        if 's_inf' in material_params: info_parts.append(f"s_inf:{material_params['s_inf']}")
        if 'delta' in material_params: info_parts.append(f"delta:{material_params['delta']}")
        
        return "Mat: " + ", ".join(info_parts) if info_parts else "Material info: No key params found."

    except FileNotFoundError:
        return f"Material info: CML file '{cml_filepath}' not found."
    except Exception as e:
        return f"Material info: Error parsing CML '{cml_filepath}': {e}"

