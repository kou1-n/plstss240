#!/usr/bin/env python3
"""
PLSTss実行スクリプト - 解析結果をターミナル出力とファイルに同時に保存
Usage: python run_plstss.py input_name [output_dir]
"""

import subprocess
import sys
import os
from datetime import datetime

def run_plstss_with_logging(input_name, output_dir="output/analysis_logs"):
    """
    PLSTssを実行し、ターミナル出力をファイルに保存

    Args:
        input_name: 入力ファイル名（拡張子なし）
        output_dir: 出力ディレクトリ（デフォルト: output/analysis_logs）
    """
    # 出力ディレクトリを作成
    os.makedirs(output_dir, exist_ok=True)

    # タイムスタンプを含む出力ファイル名を生成
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_file = os.path.join(output_dir, f"{input_name}_{timestamp}.txt")

    # PLSTss実行コマンド
    plstss_path = os.path.expanduser("~/bin/plstss2393_nagasaku")

    print(f"=== Running PLSTss Analysis ===")
    print(f"Input: {input_name}")
    print(f"Output will be saved to: {output_file}")
    print(f"{'='*40}\n")

    # ファイルとコンソールの両方に出力
    with open(output_file, 'w', encoding='utf-8') as f:
        # ヘッダー情報を書き込み
        f.write(f"PLSTss Analysis Log\n")
        f.write(f"{'='*60}\n")
        f.write(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"Input file: {input_name}\n")
        f.write(f"{'='*60}\n\n")

        # PLSTssプロセスを実行
        process = subprocess.Popen(
            plstss_path,
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            universal_newlines=True,
            bufsize=1
        )

        # 入力ファイル名を送信
        process.stdin.write(f"{input_name}\n")
        process.stdin.flush()

        # 出力を逐次読み込み、画面とファイルに出力
        for line in process.stdout:
            print(line, end='')  # コンソールに表示
            f.write(line)        # ファイルに書き込み
            f.flush()            # バッファをフラッシュ

        # プロセスの終了を待つ
        process.wait()

        # フッター情報を書き込み
        f.write(f"\n{'='*60}\n")
        f.write(f"Analysis completed at: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"Return code: {process.returncode}\n")

    print(f"\n{'='*40}")
    print(f"Analysis complete. Output saved to: {output_file}")

    # 解析が成功した場合、主要な結果ファイルの存在を確認
    result_files = [
        f"output/STS_{input_name}.txt",
        f"output/DIS_{input_name}.txt",
        f"output/NOR_{input_name}.txt",
        f"output/ENE_{input_name}.txt"
    ]

    print("\nGenerated result files:")
    for rfile in result_files:
        if os.path.exists(rfile):
            print(f"  ✓ {rfile}")
        else:
            print(f"  ✗ {rfile} (not found)")

    return output_file, process.returncode

def main():
    """メイン関数"""
    if len(sys.argv) < 2:
        print("Usage: python run_plstss.py input_name [output_dir]")
        print("Example: python run_plstss.py 1elem_f")
        print("         python run_plstss.py 1elem_f output/my_logs")
        sys.exit(1)

    input_name = sys.argv[1]
    output_dir = sys.argv[2] if len(sys.argv) > 2 else "output/analysis_logs"

    # 入力ファイルの存在確認
    input_files = [
        f"{input_name}.cml",
        f"input_files/{input_name}.cml"
    ]

    input_exists = False
    for ifile in input_files:
        if os.path.exists(ifile):
            input_exists = True
            print(f"Found input file: {ifile}")
            break

    if not input_exists:
        print(f"Error: Input file {input_name}.cml not found")
        print("Searched in: ./  and  ./input_files/")
        sys.exit(1)

    # PLSTss実行
    output_file, return_code = run_plstss_with_logging(input_name, output_dir)

    if return_code == 0:
        print("\n✓ Analysis completed successfully")
    else:
        print(f"\n⚠ Analysis finished with return code: {return_code}")

    return return_code

if __name__ == "__main__":
    sys.exit(main())