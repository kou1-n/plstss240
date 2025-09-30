#!/usr/bin/env python3
"""
Parse Block Newton convergence data showing simultaneous convergence.
Simple version without external dependencies.
"""

import re
import sys

def parse_bn_log(filename):
    """
    Parse PLSTss log file with Block Newton output to extract dual residuals.
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
            bn_match = re.search(r'\|\|Rf\|\|/\|\|F\|\|\s*:\s*([\d.E+-]+),\s*\|\|Rg\|\|\s*:\s*([\d.E+-]+)\s*\[Block Newton\]', line)
            if bn_match:
                rf_norm = float(bn_match.group(1))
                rg_norm = float(bn_match.group(2))
                data.append({
                    'Step': current_step,
                    'Iteration': current_iter,
                    'Rf_norm': rf_norm,
                    'Rg_norm': rg_norm
                })

    return data

def analyze_convergence(data):
    """
    Analyze Block Newton convergence characteristics.
    """
    print("\n" + "="*70)
    print("Block Newton Simultaneous Convergence Analysis")
    print("="*70)

    if not data:
        print("No convergence data found.")
        return

    # Count steps and iterations
    steps = set(d['Step'] for d in data)
    total_steps = len(steps)
    total_iterations = len(data)

    print(f"\nTotal load steps analyzed: {total_steps}")
    print(f"Total iterations: {total_iterations}")
    print(f"Average iterations per step: {total_iterations / total_steps:.2f}")

    # Find plastic steps (where Rg > 0)
    plastic_data = [d for d in data if d['Rg_norm'] > 1e-14]
    if plastic_data:
        plastic_steps = set(d['Step'] for d in plastic_data)
        print(f"\nPlastic load steps (||Rg|| > 1e-14):")
        print(f"  Steps: {sorted(plastic_steps)[:10]}{'...' if len(plastic_steps) > 10 else ''}")
        print(f"  Total: {len(plastic_steps)} steps")

        # Show convergence for first few plastic steps
        print("\n" + "-"*70)
        print("Detailed convergence for plastic steps:")
        print("-"*70)

        for step in sorted(plastic_steps)[:5]:  # Show first 5 plastic steps
            step_data = [d for d in data if d['Step'] == step]
            if len(step_data) > 0:
                print(f"\nStep {step}:")
                for d in step_data:
                    print(f"  Iter {d['Iteration']}: ||Rf||/||F|| = {d['Rf_norm']:.3e}, ||Rg|| = {d['Rg_norm']:.3e}")

                # Check for simultaneous decrease
                if len(step_data) >= 2:
                    rf_decreased = all(step_data[i]['Rf_norm'] <= step_data[i-1]['Rf_norm']
                                     for i in range(1, len(step_data)))
                    rg_decreased = all(step_data[i]['Rg_norm'] <= step_data[i-1]['Rg_norm']
                                     for i in range(1, len(step_data)))

                    if rf_decreased and rg_decreased:
                        print("  -> Both residuals decrease simultaneously (Block Newton characteristic)")
                    elif rf_decreased:
                        print("  -> Only equilibrium residual decreases")
                    elif rg_decreased:
                        print("  -> Only yield residual decreases")
    else:
        print("\nNo plastic steps detected (all ||Rg|| values near zero)")

    # Check simultaneity for multi-iteration steps
    print("\n" + "-"*70)
    print("Multi-iteration steps (key for observing simultaneous convergence):")
    print("-"*70)

    multi_iter_steps = []
    for step in steps:
        step_data = [d for d in data if d['Step'] == step]
        max_iter = max(d['Iteration'] for d in step_data)
        if max_iter > 1:
            multi_iter_steps.append(step)

    if multi_iter_steps:
        print(f"\nSteps requiring multiple iterations: {multi_iter_steps[:10]}{'...' if len(multi_iter_steps) > 10 else ''}")
        print(f"Total: {len(multi_iter_steps)} steps")

        # Check simultaneity
        simultaneous_count = 0
        for step in multi_iter_steps:
            step_data = sorted([d for d in data if d['Step'] == step],
                             key=lambda x: x['Iteration'])
            if len(step_data) >= 2:
                rf_decreases = all(step_data[i]['Rf_norm'] <= step_data[i-1]['Rf_norm']
                                  for i in range(1, len(step_data)))
                rg_decreases = all(step_data[i]['Rg_norm'] <= step_data[i-1]['Rg_norm']
                                  for i in range(1, len(step_data)))
                if rf_decreases and rg_decreases:
                    simultaneous_count += 1

        print(f"\nSteps with simultaneous convergence: {simultaneous_count}/{len(multi_iter_steps)}")
        print(f"Simultaneity ratio: {100*simultaneous_count/len(multi_iter_steps):.1f}%")
        print("\nNote: High simultaneity ratio confirms Block Newton implementation.")
    else:
        print("\nAll steps converged in single iteration.")
        print("This indicates either:")
        print("  1. Problem is in elastic regime")
        print("  2. Tolerance is too loose")
        print("  3. Block Newton's excellent convergence (no local iteration needed)")

def save_to_csv(data, output_file):
    """Save convergence data to CSV."""
    with open(output_file, 'w') as f:
        f.write("Step,Iteration,Rf_norm,Rg_norm\n")
        for d in data:
            f.write(f"{d['Step']},{d['Iteration']},{d['Rf_norm']:.6e},{d['Rg_norm']:.6e}\n")
    print(f"\nConvergence data saved to: {output_file}")

def main():
    if len(sys.argv) != 2:
        print("Usage: python parse_bn_convergence_simple.py <log_file>")
        print("Example: python parse_bn_convergence_simple.py output/1elem_f_bn_full.log")
        sys.exit(1)

    log_file = sys.argv[1]

    try:
        data = parse_bn_log(log_file)

        if not data:
            print("No Block Newton convergence data found in log file.")
            print("Make sure the log file is from a Block Newton run (MATYPE=5 or 6).")
            return

        analyze_convergence(data)

        # Save to CSV
        output_file = log_file.replace('.log', '_bn_convergence.csv')
        save_to_csv(data, output_file)

        print("\n" + "="*70)
        print("Analysis complete!")
        print("="*70)

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