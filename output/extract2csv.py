#!/usr/bin/env python3
"""Extract von Mises stress, equivalent plastic strain, yystress, and yystrain from a RES_*.cml file.

このスクリプトは有限要素法解析の結果ファイル（RES_*.cml）から材料の応力・ひずみ情報を抽出します。
抽出される主な物理量：
- von Mises応力（相当応力）：多軸応力状態を単一の値で表現
- 相当塑性ひずみ：塑性変形の累積量
- yy方向の応力・ひずみ：y軸方向の垂直応力・ひずみ成分

Usage: python extract2csv.py input.cml output.csv

The script reads /ELMTL/ blocks and outputs rows:
    step, element_id, von, eps, yystress, yystrain
in CSV format.
"""
import sys
import csv


def parse_res(filename):
    """
    RES_*.cmlファイルを解析して応力・ひずみデータを抽出する関数
    
    Args:
        filename (str): 入力するcmlファイルのパス
    
    Returns:
        list: (step, element_id, von, eps, yystress, yystrain)のタプルのリスト
    """
    results = []
    # ファイルを開いて全行を読み込み、各行の末尾の改行文字を削除
    with open(filename) as f:
        lines = [line.rstrip() for line in f]
    
    # ファイルを1行ずつ処理するためのインデックス
    i = 0
    while i < len(lines):
        # /ELMTL/ブロックの開始を検出（各ステップの要素材料データブロック）
        if lines[i].strip() == '/ELMTL/':
            # /ELMTL/の次の行：ステップ番号やその他の情報を含む
            # 例: "1 0 0 0" のような形式（最初の数値がステップ番号）
            step_info = lines[i + 1].split()
            if not step_info:
                i += 1
                continue
            step = int(step_info[0])  # ステップ番号を取得
            
            # /ELMTL/から2行目：要素数の情報
            # 例: "1 0 0" のような形式（最初の数値が要素数）
            nel = int(lines[i + 2].split()[0])  # 要素数を取得
            i += 3  # 次のデータ行へ移動
            
            # 各要素のデータを処理（要素数分ループ）
            for _ in range(nel):
                # 1行目：応力データ行
                # 形式: "要素ID xxstress yystress zzstress ..."
                # 例: "1 0.5 1.2 0.3 ..." （要素1のxx,yy,zz方向の応力など）
                stress_tokens = lines[i].split()
                if not stress_tokens:
                    break
                el_id = int(stress_tokens[0])  # 要素IDを取得（1列目）
                yystress = float(stress_tokens[2])  # yy方向の応力を取得（3列目）
                i += 1
                
                # 2行目：ひずみデータ行
                # 形式: "xxstrain yystrain zzstrain ..."
                # 例: "0.001 0.002 0.0015 ..." （xx,yy,zz方向のひずみなど）
                strain_tokens = lines[i].split()
                yystrain = float(strain_tokens[1])  # yy方向のひずみを取得（2列目）
                i += 1
                
                # 3行目：von Mises応力と相当塑性ひずみなどの情報
                # 形式: "von_mises_stress equivalent_plastic_strain その他..."
                # 例: "150.5 0.003 ..." （von Mises応力、相当塑性ひずみ、その他の値）
                vals = lines[i].split()
                if len(vals) >= 3:
                    von = float(vals[0])  # von Mises応力（相当応力）を取得（1列目）
                    eps = float(vals[1])  # 相当塑性ひずみを取得（2列目）
                    # 抽出したデータをタプルとして結果リストに追加
                    results.append((step, el_id, von, eps, yystress, yystrain))
                i += 1
        else:
            # /ELMTL/ブロック以外の行はスキップ
            i += 1
    return results


def main():
    """
    メイン関数：コマンドライン引数を処理してCSVファイルを出力
    """
    # コマンドライン引数のチェック（入力ファイルと出力ファイルの2つが必要）
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} input.cml output.csv")
        sys.exit(1)
    
    # 入力ファイル名と出力ファイル名を取得
    src, dest = sys.argv[1], sys.argv[2]
    
    # cmlファイルを解析してデータを抽出
    data = parse_res(src)
    
    # 抽出したデータをCSVファイルに書き込み
    with open(dest, 'w', newline='') as f:
        writer = csv.writer(f)
        # ヘッダー行を書き込み
        writer.writerow(['step', 'element', 'von', 'eps', 'yystress', 'yystrain'])
        # データ行を書き込み（各タプルが1行になる）
        writer.writerows(data)


if __name__ == '__main__':
    main()
