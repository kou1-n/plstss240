#!/usr/bin/env python3
"""
パス設定テスト用スクリプト
"""

import os

def test_paths():
    # パス設定（main_visualizer.pyと同じロジック）
    script_dir = os.path.dirname(os.path.abspath(__file__))  # scripts/フォルダ
    work_dir = os.path.dirname(script_dir)  # WORKディレクトリ
    output_dir = os.path.join(work_dir, "output")  # WORK/outputディレクトリ
    
    print("=== パス設定テスト ===")
    print(f"スクリプトディレクトリ: {script_dir}")
    print(f"WORKディレクトリ: {work_dir}")
    print(f"OUTPUTディレクトリ: {output_dir}")
    print()
    
    # CMLファイルのチェック
    case_prefix = "cube250514"
    cml_file_path = os.path.join(work_dir, f"{case_prefix}.cml")
    print(f"CMLファイルパス: {cml_file_path}")
    print(f"CMLファイル存在: {os.path.exists(cml_file_path)}")
    print()
    
    # 結果ファイルのチェック
    file_types = ["DIS", "ENE", "NOR", "STS", "TMP"]
    result_files = {
        file_type: os.path.join(output_dir, f"{file_type}_{case_prefix}.txt")
        for file_type in file_types
    }
    
    print("=== 結果ファイル存在確認 ===")
    for file_type, file_path in result_files.items():
        exists = os.path.exists(file_path)
        print(f"{file_type}: {file_path} -> {exists}")
    
    # RES_*.cmlファイルのチェック
    res_file = os.path.join(output_dir, f"RES_{case_prefix}.cml")
    print(f"RES: {res_file} -> {os.path.exists(res_file)}")

if __name__ == "__main__":
    test_paths() 