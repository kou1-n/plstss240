#!/usr/bin/env python3
"""
Parse PLSTss log file to extract residual norms for convergence analysis.
Creates CSV data for plotting quadratic convergence of Return Mapping (RM).
"""

import re
import csv
from pathlib import Path
import argparse


def parse_residual_norms(log_path, target_step=None):
    """
    Parse PLSTss log file to extract residual norms per iteration.

    Args:
        log_path: Path to log file
        target_step: Specific step to analyze (e.g., step 13 with 2 iterations)
                    If None, analyzes all steps with multiple iterations

    Returns:
        data: Dictionary with step numbers as keys, containing iteration and residual data
    """
    with open(log_path, 'r', encoding='utf-8', errors='ignore') as f:
        content = f.read()

    # Split content by loading steps
    step_pattern = r'\*\*\*\*\* Loading Step :\s+(\d+)/\s*(\d+)'
    step_splits = re.split(step_pattern, content)[1:]  # Skip before first match

    data = {}

    # Process in groups of 3 (step_num, total_steps, content)
    for i in range(0, len(step_splits), 3):
        if i + 2 >= len(step_splits):
            break

        step_num = int(step_splits[i])
        total_steps = int(step_splits[i + 1])
        step_content = step_splits[i + 2]

        # Find all iterations and residuals in this step
        iter_pattern = r'Iteration\s*:\s*(\d+)'
        residual_pattern = r'rnorm\s*:\s*([\d.E+-]+).*?fnorm\s*:\s*([\d.E+-]+).*?RESIDUAL\s*:\s*([\d.E+-]+)'

        iterations = re.findall(iter_pattern, step_content)
        residuals = re.findall(residual_pattern, step_content)

        if len(iterations) > 1:  # Only interested in steps with multiple iterations
            step_data = {
                'step': step_num,
                'iterations': [],
                'rnorm': [],
                'fnorm': [],
                'residual': []
            }

            for idx, (iter_num, (rnorm, fnorm, residual)) in enumerate(zip(iterations, residuals)):
                step_data['iterations'].append(int(iter_num))
                step_data['rnorm'].append(float(rnorm))
                step_data['fnorm'].append(float(fnorm))
                step_data['residual'].append(float(residual))

            data[step_num] = step_data

    if target_step and target_step in data:
        return {target_step: data[target_step]}

    return data


def write_convergence_csv(data, output_path):
    """
    Write residual convergence data to CSV file.
    """
    with open(output_path, 'w', newline='', encoding='utf-8') as csvfile:
        fieldnames = ['Step', 'Iteration', 'Rnorm', 'Fnorm', 'Residual', 'Log10_Residual']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

        writer.writeheader()

        for step_num, step_data in sorted(data.items()):
            for i in range(len(step_data['iterations'])):
                # Calculate log10 for semi-log plots
                import math
                log_residual = math.log10(step_data['residual'][i]) if step_data['residual'][i] > 0 else -20

                writer.writerow({
                    'Step': step_num,
                    'Iteration': step_data['iterations'][i],
                    'Rnorm': f"{step_data['rnorm'][i]:.6E}",
                    'Fnorm': f"{step_data['fnorm'][i]:.6E}",
                    'Residual': f"{step_data['residual'][i]:.6E}",
                    'Log10_Residual': f"{log_residual:.4f}"
                })

    print(f"Convergence data saved to: {output_path}")


def analyze_convergence_rate(data):
    """
    Analyze convergence rate (linear, quadratic, etc.)
    """
    analysis = []

    for step_num, step_data in sorted(data.items()):
        residuals = step_data['residual']

        if len(residuals) >= 2:
            # Calculate convergence rate between iterations
            rates = []
            for i in range(1, len(residuals)):
                if residuals[i-1] > 0 and residuals[i] > 0:
                    # Quadratic convergence: r_{k+1} ≈ C * r_k^2
                    # Linear convergence: r_{k+1} ≈ C * r_k
                    linear_rate = residuals[i] / residuals[i-1] if residuals[i-1] != 0 else 0
                    quadratic_estimate = residuals[i] / (residuals[i-1]**2) if residuals[i-1] != 0 else 0

                    rates.append({
                        'iteration': i + 1,
                        'linear_rate': linear_rate,
                        'quadratic_estimate': quadratic_estimate
                    })

            analysis.append({
                'step': step_num,
                'total_iterations': len(residuals),
                'initial_residual': residuals[0],
                'final_residual': residuals[-1],
                'reduction_factor': residuals[-1] / residuals[0] if residuals[0] != 0 else 0,
                'convergence_rates': rates
            })

    return analysis


def main():
    parser = argparse.ArgumentParser(description='Parse residual norms for convergence analysis')
    parser.add_argument('log_file', type=str, help='Path to log file')
    parser.add_argument('-o', '--output', type=str, help='Output CSV file path')
    parser.add_argument('-s', '--step', type=int, help='Analyze specific step number')

    args = parser.parse_args()

    log_path = Path(args.log_file)
    if not log_path.exists():
        print(f"Error: Log file not found: {log_path}")
        return 1

    try:
        print(f"Parsing log file: {log_path}")

        # Parse residual norms
        data = parse_residual_norms(log_path, args.step)

        if not data:
            print("No steps with multiple iterations found.")
            return 0

        # Determine output path
        if args.output:
            output_path = Path(args.output)
        else:
            suffix = f"_step{args.step}" if args.step else "_all"
            output_path = log_path.parent / f"{log_path.stem}_convergence{suffix}.csv"

        # Write CSV
        write_convergence_csv(data, output_path)

        # Analyze convergence rates
        print("\nConvergence Analysis:")
        print("=" * 50)

        analysis = analyze_convergence_rate(data)
        for step_analysis in analysis:
            print(f"\nStep {step_analysis['step']}:")
            print(f"  Iterations: {step_analysis['total_iterations']}")
            print(f"  Initial residual: {step_analysis['initial_residual']:.6E}")
            print(f"  Final residual: {step_analysis['final_residual']:.6E}")
            print(f"  Total reduction: {step_analysis['reduction_factor']:.6E}")

            for rate_info in step_analysis['convergence_rates']:
                print(f"  Iteration {rate_info['iteration']}: "
                      f"linear_rate={rate_info['linear_rate']:.6E}, "
                      f"quadratic_est={rate_info['quadratic_estimate']:.6E}")

        # Save analysis to text file
        analysis_path = output_path.parent / f"{output_path.stem}_analysis.txt"
        with open(analysis_path, 'w', encoding='utf-8') as f:
            f.write("Residual Norm Convergence Analysis\n")
            f.write("=" * 50 + "\n\n")

            for step_analysis in analysis:
                f.write(f"Loading Step {step_analysis['step']}:\n")
                f.write(f"  Total iterations: {step_analysis['total_iterations']}\n")
                f.write(f"  Initial residual: {step_analysis['initial_residual']:.6E}\n")
                f.write(f"  Final residual: {step_analysis['final_residual']:.6E}\n")
                f.write(f"  Reduction factor: {step_analysis['reduction_factor']:.6E}\n\n")

                f.write("  Convergence rates:\n")
                for rate_info in step_analysis['convergence_rates']:
                    f.write(f"    Iteration {rate_info['iteration']}: ")
                    f.write(f"linear={rate_info['linear_rate']:.6E}, ")
                    f.write(f"quadratic={rate_info['quadratic_estimate']:.6E}\n")
                f.write("\n")

        print(f"\nAnalysis saved to: {analysis_path}")

        return 0

    except Exception as e:
        print(f"Error: {e}")
        import traceback
        traceback.print_exc()
        return 1


if __name__ == '__main__':
    import sys

    # If no arguments provided, use default
    if len(sys.argv) == 1:
        default_log = '//akita/nagasaku/WORK/output/logs/1elem_f_test_current.log'
        if Path(default_log).exists():
            sys.argv.extend([default_log, '-s', '13'])  # Focus on step 13 which has 2 iterations
            print(f"Using default: {default_log}, analyzing step 13")

    sys.exit(main())