#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
plstss240 可視化ツール簡単起動スクリプト
=====================================

使用方法:
  python visualize.py                    # すべての解析を実行
  python visualize.py dual               # 2軸プロットのみ
  python visualize.py comprehensive      # 包括的解析のみ
  python visualize.py hardening         # 硬化曲線のみ
  python visualize.py --element 5        # 要素5を解析

Author: AI Assistant
"""

import sys
import subprocess
from pathlib import Path

def print_usage():
    """使用方法を表示"""
    print("""
plstss240 可視化ツール
====================

使用方法:
  python visualize.py [mode] [options]

モード:
  (なし)          - すべての解析を実行 (デフォルト)
  dual           - 2軸プロット (DP基準 vs 硬化応力)
  comprehensive  - 包括的解析 (4つのサブプロット)
  hardening      - 理論硬化曲線

オプション:
  --element N    - 解析する要素番号 (デフォルト: 8)
  --file F       - RES ファイル名 (デフォルト: RES_cube250514.cml)

例:
  python visualize.py                    # すべての解析
  python visualize.py dual               # 2軸プロットのみ
  python visualize.py --element 5        # 要素5の全解析
  python visualize.py dual --element 10  # 要素10の2軸プロット

出力場所: ../output/
""")

def main():
    """メイン関数"""
    args = sys.argv[1:]
    
    # ヘルプ表示
    if not args or args[0] in ['-h', '--help', 'help']:
        print_usage()
        return 0
    
    # 統合ツールのパス
    visualizer_script = Path(__file__).parent / "plstss_visualizer.py"
    
    if not visualizer_script.exists():
        print("Error: plstss_visualizer.py が見つかりません")
        return 1
    
    # モード判定
    mode = 'all'  # デフォルト
    other_args = []
    
    if args and args[0] in ['dual', 'comprehensive', 'hardening']:
        mode = args[0]
        other_args = args[1:]
    else:
        other_args = args
    
    # コマンド構築
    cmd = ['python', str(visualizer_script), '--mode', mode] + other_args
    
    print(f"実行中: {' '.join(cmd[1:])}")
    print()
    
    # 実行
    try:
        result = subprocess.run(cmd, check=False)
        return result.returncode
    except KeyboardInterrupt:
        print("\n中断されました")
        return 1
    except Exception as e:
        print(f"Error: {e}")
        return 1

if __name__ == "__main__":
    sys.exit(main()) 