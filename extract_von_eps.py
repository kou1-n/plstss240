#!/usr/bin/env python3
"""Extract von Mises stress and equivalent plastic strain from a RES_*.cml file.

Usage: python extract_von_eps.py input.cml output.csv

The script reads /ELMTL/ blocks and outputs rows:
    step, element_id, von, eps
in CSV format.
"""
import sys
import csv


def parse_res(filename):
    results = []
    with open(filename) as f:
        lines = [line.rstrip() for line in f]
    i = 0
    while i < len(lines):
        if lines[i].strip() == '/ELMTL/':
            # Step information and number of elements
            step_info = lines[i + 1].split()
            if not step_info:
                i += 1
                continue
            step = int(step_info[0])
            nel = int(lines[i + 2].split()[0])
            i += 3
            for _ in range(nel):
                # Stress line (contains element id)
                tokens = lines[i].split()
                if not tokens:
                    break
                el_id = int(tokens[0])
                i += 1
                # Skip strain line
                i += 1
                # Line with von and eps
                vals = lines[i].split()
                if len(vals) >= 2:
                    von = float(vals[0])
                    eps = float(vals[1])
                    results.append((step, el_id, von, eps))
                i += 1
        else:
            i += 1
    return results


def main():
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} input.cml output.csv")
        sys.exit(1)
    src, dest = sys.argv[1], sys.argv[2]
    data = parse_res(src)
    with open(dest, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['step', 'element', 'von', 'eps'])
        writer.writerows(data)


if __name__ == '__main__':
    main()
