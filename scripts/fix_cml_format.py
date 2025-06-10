#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
CMLファイルのフォーマット修正スクリプト
"""

def fix_cml_format():
    """cube250514.cmlのインデントを修正"""
    
    cml_file = "cube250514.cml"
    
    print(f"📝 {cml_file}のフォーマットを修正中...")
    
    with open(cml_file, 'r') as f:
        lines = f.readlines()
    
    # 修正が必要な行を特定（2-9行目のインデント）
    fixed_lines = []
    
    for i, line in enumerate(lines):
        line_num = i + 1
        
        # 73-80行目のインデントを修正（2-9番のノード）
        if 73 <= line_num <= 80:
            # 行の先頭が数字から始まっている場合、7つのスペースを追加
            if line.strip() and line[0].isdigit():
                # 既存のスペースを除去してから正しいインデントを追加
                content = line.strip()
                fixed_line = "       " + content + "\n"
                fixed_lines.append(fixed_line)
                print(f"  修正: 行{line_num}: '{line.strip()}' → '       {content}'")
            else:
                fixed_lines.append(line)
        else:
            fixed_lines.append(line)
    
    # ファイルを書き戻し
    with open(cml_file, 'w') as f:
        f.writelines(fixed_lines)
    
    print(f"  ✓ {cml_file}のフォーマット修正完了")

if __name__ == "__main__":
    fix_cml_format() 