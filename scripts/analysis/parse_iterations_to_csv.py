#!/usr/bin/env python3
"""
Parse PLSTss analysis log file and create CSV file with iteration data.
For creating cumulative iteration plot similar to Yamamoto et al. (2020) Figure 4.
"""

import re
import csv
from pathlib import Path
import argparse


def parse_log_file(log_path):
    """
    Parse PLSTss log file to extract iteration counts per loading step.

    Returns:
        data: List of dictionaries containing step data
    """
    with open(log_path, 'r', encoding='utf-8', errors='ignore') as f:
        content = f.read()

    # Find all loading steps
    step_pattern = r'\*\*\*\*\* Loading Step :\s+(\d+)/\s*(\d+)'
    steps_data = re.findall(step_pattern, content)

    if not steps_data:
        raise ValueError("No loading steps found in log file")

    total_steps = int(steps_data[0][1])

    # Split content by loading steps
    step_sections = re.split(r'\*\*\*\*\* Loading Step :', content)[1:]

    data = []
    cumulative_iterations = 0

    for section in step_sections:
        # Extract step number
        step_match = re.match(r'\s*(\d+)/\s*(\d+)', section)
        if not step_match:
            continue

        step_num = int(step_match.group(1))

        # Count iterations in this step
        iter_count = len(re.findall(r'^\s*Iteration\s*:', section, re.MULTILINE))

        # Update cumulative
        cumulative_iterations += iter_count

        # Calculate normalized time (0 to 1)
        time_normalized = step_num / total_steps

        data.append({
            'Step': step_num,
            'Time_Normalized': time_normalized,
            'Iterations': iter_count,
            'Cumulative_Iterations': cumulative_iterations
        })

    return data, total_steps


def write_csv(data, output_path):
    """
    Write iteration data to CSV file.
    """
    with open(output_path, 'w', newline='', encoding='utf-8') as csvfile:
        fieldnames = ['Step', 'Time_Normalized', 'Iterations', 'Cumulative_Iterations']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

        # Write header
        writer.writeheader()

        # Write data
        for row in data:
            writer.writerow(row)

    print(f"CSV file saved to: {output_path}")


def main():
    parser = argparse.ArgumentParser(description='Parse PLSTss log file and create CSV with iteration data')
    parser.add_argument('log_file', type=str, help='Path to log file')
    parser.add_argument('-o', '--output', type=str, help='Output CSV file path')

    args = parser.parse_args()

    log_path = Path(args.log_file)
    if not log_path.exists():
        print(f"Error: Log file not found: {log_path}")
        return 1

    try:
        # Parse log file
        print(f"Parsing log file: {log_path}")
        data, total_steps = parse_log_file(log_path)

        # Calculate statistics
        total_iterations = data[-1]['Cumulative_Iterations'] if data else 0
        avg_iterations = total_iterations / len(data) if data else 0
        max_iterations = max(row['Iterations'] for row in data) if data else 0
        max_iter_step = next((row['Step'] for row in data if row['Iterations'] == max_iterations), None)

        print(f"Found {len(data)} loading steps (out of {total_steps} total)")
        print(f"Total iterations: {total_iterations}")
        print(f"Average iterations per step: {avg_iterations:.2f}")
        print(f"Maximum iterations: {max_iterations} (at step {max_iter_step})")

        # Determine output path
        if args.output:
            output_path = Path(args.output)
        else:
            output_path = log_path.parent / f"{log_path.stem}_iterations.csv"

        # Write CSV
        write_csv(data, output_path)

        # Also write a summary file
        summary_path = output_path.parent / f"{output_path.stem}_summary.txt"
        with open(summary_path, 'w', encoding='utf-8') as f:
            f.write(f"PLSTss Iteration Analysis Summary\n")
            f.write(f"================================\n\n")
            f.write(f"Log file: {log_path}\n")
            f.write(f"Total loading steps analyzed: {len(data)}/{total_steps}\n")
            f.write(f"Total iterations: {total_iterations}\n")
            f.write(f"Average iterations per step: {avg_iterations:.2f}\n")
            f.write(f"Maximum iterations: {max_iterations} (at step {max_iter_step})\n\n")

            # Find steps with multiple iterations
            multi_iter_steps = [row['Step'] for row in data if row['Iterations'] > 1]
            if multi_iter_steps:
                f.write(f"Steps with multiple iterations: {multi_iter_steps}\n")
            else:
                f.write(f"All steps converged in 1 iteration\n")

        print(f"Summary saved to: {summary_path}")

        return 0

    except Exception as e:
        print(f"Error: {e}")
        return 1


if __name__ == '__main__':
    import sys

    # If no arguments provided, use default log file
    if len(sys.argv) == 1:
        default_log = '//akita/nagasaku/WORK/output/logs/1elem_f_test_current.log'
        if Path(default_log).exists():
            sys.argv.append(default_log)
            print(f"Using default log file: {default_log}")

    sys.exit(main())