#!/usr/bin/env python3
"""
Parse Block Newton convergence data showing simultaneous convergence of:
1. Equilibrium residual (||Rf||/||F||)
2. Yield function residual (||Rg||)

This demonstrates the key characteristic of Block Newton method according to
Yamamoto et al. (2021): both residuals converge simultaneously without local iteration.
"""

import re
import sys
import pandas as pd
import numpy as np

def parse_bn_log(filename):
    """
    Parse PLSTss log file with Block Newton output to extract dual residuals.

    Returns:
        DataFrame with columns: Step, Iteration, Rf_norm, Rg_norm, Converged
    """
    data = []
    current_step = 0
    current_iter = 0

    with open(filename, 'r', encoding='utf-8', errors='ignore') as f:
        for line in f:
            # Match loading step
            step_match = re.search(r'Loading Step\s*:\s*(\d+)/\s*(\d+)', line)
            if step_match:
                current_step = int(step_match.group(1))
                current_iter = 0
                continue

            # Match iteration number
            iter_match = re.search(r'Iteration\s*:\s*(\d+)', line)
            if iter_match:
                current_iter = int(iter_match.group(1))
                continue

            # Match Block Newton dual residuals
            # Format: ||Rf||/||F|| : value, ||Rg|| : value [Block Newton]
            bn_match = re.search(r'\|\|Rf\|\|/\|\|F\|\|\s*:\s*([\d.E+-]+),\s*\|\|Rg\|\|\s*:\s*([\d.E+-]+)\s*\[Block Newton\]', line)
            if bn_match:
                rf_norm = float(bn_match.group(1))
                rg_norm = float(bn_match.group(2))
                data.append({
                    'Step': current_step,
                    'Iteration': current_iter,
                    'Rf_norm': rf_norm,
                    'Rg_norm': rg_norm,
                    'Converged': False
                })
                continue

            # Match convergence message
            if 'Block Newton: Both criteria satisfied' in line:
                if data:
                    data[-1]['Converged'] = True

    return pd.DataFrame(data)

def analyze_convergence(df):
    """
    Analyze Block Newton convergence characteristics.
    """
    print("\n" + "="*70)
    print("Block Newton Simultaneous Convergence Analysis")
    print("="*70)

    # Group by step to analyze iterations per step
    steps_summary = df.groupby('Step').agg({
        'Iteration': 'max',
        'Rf_norm': ['first', 'last', 'min'],
        'Rg_norm': ['first', 'last', 'max']
    })

    print(f"\nTotal load steps analyzed: {df['Step'].nunique()}")
    print(f"Total iterations: {len(df)}")
    print(f"Average iterations per step: {len(df) / df['Step'].nunique():.2f}")

    # Find plastic steps (where Rg > 0)
    plastic_df = df[df['Rg_norm'] > 1e-14]
    if not plastic_df.empty:
        print(f"\nPlastic load steps (Rg > 0):")
        plastic_steps = plastic_df['Step'].unique()
        print(f"  Steps: {plastic_steps[:10]}{'...' if len(plastic_steps) > 10 else ''}")
        print(f"  Total: {len(plastic_steps)} steps")

        # Analyze convergence rates for plastic steps
        for step in plastic_steps[:3]:  # Show first 3 plastic steps
            step_data = df[df['Step'] == step]
            if len(step_data) > 1:
                print(f"\n  Step {step} convergence:")
                for _, row in step_data.iterrows():
                    print(f"    Iter {row['Iteration']}: ||Rf||/||F|| = {row['Rf_norm']:.3e}, ||Rg|| = {row['Rg_norm']:.3e}")

                # Calculate convergence rate
                if len(step_data) >= 2:
                    rf_rates = []
                    rg_rates = []
                    for i in range(1, len(step_data)):
                        prev = step_data.iloc[i-1]
                        curr = step_data.iloc[i]
                        if prev['Rf_norm'] > 1e-16:
                            rf_rate = np.log10(curr['Rf_norm'] / prev['Rf_norm'])
                            rf_rates.append(rf_rate)
                        if prev['Rg_norm'] > 1e-16:
                            rg_rate = np.log10(curr['Rg_norm'] / prev['Rg_norm'])
                            rg_rates.append(rg_rate)

                    if rf_rates:
                        print(f"    Rf convergence rate: {np.mean(rf_rates):.2f} (log scale)")
                    if rg_rates:
                        print(f"    Rg convergence rate: {np.mean(rg_rates):.2f} (log scale)")
    else:
        print("\nNo plastic steps detected (all Rg values near zero)")

    # Check simultaneity: both residuals should decrease together
    print("\n" + "-"*70)
    print("Simultaneity Check (key characteristic of Block Newton):")
    print("-"*70)

    multi_iter_steps = df[df['Iteration'] > 1]['Step'].unique()
    if len(multi_iter_steps) > 0:
        simultaneous_count = 0
        for step in multi_iter_steps:
            step_data = df[df['Step'] == step].sort_values('Iteration')
            if len(step_data) >= 2:
                # Check if both residuals decrease
                rf_decreasing = all(step_data['Rf_norm'].diff().dropna() <= 0)
                rg_decreasing = all(step_data['Rg_norm'].diff().dropna() <= 0) or all(step_data['Rg_norm'] < 1e-14)
                if rf_decreasing and rg_decreasing:
                    simultaneous_count += 1

        print(f"Steps with multiple iterations: {len(multi_iter_steps)}")
        print(f"Steps with simultaneous convergence: {simultaneous_count}")
        print(f"Simultaneity ratio: {100*simultaneous_count/len(multi_iter_steps):.1f}%")
    else:
        print("All steps converged in single iteration (characteristic of well-conditioned Block Newton)")

    return df

def save_to_csv(df, output_file):
    """Save convergence data to CSV for plotting."""
    df.to_csv(output_file, index=False)
    print(f"\nConvergence data saved to: {output_file}")

    # Create summary for plotting
    summary = df.groupby('Step').agg({
        'Iteration': 'max',
        'Rf_norm': 'min',
        'Rg_norm': 'max'
    }).reset_index()
    summary.columns = ['Step', 'Iterations', 'Min_Rf', 'Max_Rg']

    summary_file = output_file.replace('.csv', '_summary.csv')
    summary.to_csv(summary_file, index=False)
    print(f"Summary data saved to: {summary_file}")

def main():
    if len(sys.argv) != 2:
        print("Usage: python parse_bn_convergence.py <log_file>")
        print("Example: python parse_bn_convergence.py output/1elem_f_bn_full.log")
        sys.exit(1)

    log_file = sys.argv[1]

    try:
        df = parse_bn_log(log_file)

        if df.empty:
            print("No Block Newton convergence data found in log file.")
            print("Make sure the log file is from a Block Newton run (MATYPE=5 or 6).")
            return

        analyzed_df = analyze_convergence(df)

        # Save to CSV
        output_file = log_file.replace('.log', '_bn_convergence.csv')
        save_to_csv(analyzed_df, output_file)

        print("\n" + "="*70)
        print("Analysis complete!")

    except FileNotFoundError:
        print(f"Error: File '{log_file}' not found.")
        sys.exit(1)
    except Exception as e:
        print(f"Error processing log file: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()