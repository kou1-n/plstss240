#!/usr/bin/env python3
"""
Parse PLSTss analysis log file and create cumulative iteration plot.
Similar to Figure 4 in Yamamoto et al. (2020) - Return Mapping (RM) portion.
"""

import re
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import argparse


def parse_log_file(log_path):
    """
    Parse PLSTss log file to extract iteration counts per loading step.

    Returns:
        steps: List of step numbers
        iterations: List of iteration counts for each step
        cumulative: Cumulative sum of iterations
    """
    with open(log_path, 'r') as f:
        content = f.read()

    # Find all loading steps
    step_pattern = r'\*\*\*\*\* Loading Step :\s+(\d+)/\s*(\d+)'
    steps_data = re.findall(step_pattern, content)

    if not steps_data:
        raise ValueError("No loading steps found in log file")

    total_steps = int(steps_data[0][1])

    # Split content by loading steps
    step_sections = re.split(r'\*\*\*\*\* Loading Step :', content)[1:]

    steps = []
    iterations = []

    for section in step_sections:
        # Extract step number
        step_match = re.match(r'\s*(\d+)/\s*(\d+)', section)
        if not step_match:
            continue

        step_num = int(step_match.group(1))

        # Count iterations in this step
        iter_count = len(re.findall(r'^\s*Iteration\s*:', section, re.MULTILINE))

        steps.append(step_num)
        iterations.append(iter_count)

    # Calculate cumulative iterations
    cumulative = np.cumsum(iterations)

    return steps, iterations, cumulative, total_steps


def plot_cumulative_iterations(steps, iterations, cumulative, total_steps, output_path=None, title=None):
    """
    Create cumulative iteration plot similar to Yamamoto et al. Figure 4.
    """
    fig, ax = plt.subplots(figsize=(10, 6))

    # Normalize time (loading steps) to [0, 1]
    time_normalized = np.array(steps) / total_steps

    # Plot cumulative iterations
    ax.plot(time_normalized, cumulative, 'b-', linewidth=2, label='Return Mapping (RM)',
            marker='o', markersize=3, markevery=5)

    # Styling
    ax.set_xlabel('Time (normalized)', fontsize=12)
    ax.set_ylabel('Cumulative Iterations', fontsize=12)

    if title:
        ax.set_title(title, fontsize=14)
    else:
        ax.set_title('Cumulative Iterations - Return Mapping Method', fontsize=14)

    ax.grid(True, alpha=0.3)
    ax.legend(loc='upper left', fontsize=11)

    # Set x-axis limits
    ax.set_xlim(0, 1.0)
    ax.set_ylim(0, max(cumulative) * 1.1)

    # Add some statistics as text
    avg_iter = cumulative[-1] / len(steps)
    max_iter_step = steps[np.argmax(iterations)]
    max_iter = max(iterations)

    stats_text = f'Total iterations: {cumulative[-1]}\n'
    stats_text += f'Average per step: {avg_iter:.2f}\n'
    stats_text += f'Max iterations: {max_iter} (step {max_iter_step})'

    ax.text(0.98, 0.05, stats_text, transform=ax.transAxes,
            fontsize=10, verticalalignment='bottom', horizontalalignment='right',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    plt.tight_layout()

    if output_path:
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        print(f"Plot saved to: {output_path}")

    return fig, ax


def plot_iterations_per_step(steps, iterations, output_path=None):
    """
    Create bar plot of iterations per loading step.
    """
    fig, ax = plt.subplots(figsize=(12, 6))

    # Create bar plot
    bars = ax.bar(steps, iterations, width=0.8, color='steelblue', edgecolor='black', linewidth=0.5)

    # Highlight steps with multiple iterations
    for i, (step, iter_count) in enumerate(zip(steps, iterations)):
        if iter_count > 1:
            bars[i].set_color('orange')

    ax.set_xlabel('Loading Step', fontsize=12)
    ax.set_ylabel('Number of Iterations', fontsize=12)
    ax.set_title('Iterations per Loading Step', fontsize=14)

    # Set integer y-axis
    ax.yaxis.set_major_locator(plt.MaxNLocator(integer=True))

    # Add grid
    ax.grid(True, axis='y', alpha=0.3)

    # Add legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='steelblue', label='1 iteration'),
        Patch(facecolor='orange', label='Multiple iterations')
    ]
    ax.legend(handles=legend_elements, loc='upper right')

    # Rotate x-axis labels if too many steps
    if len(steps) > 50:
        plt.xticks(steps[::5], rotation=45, ha='right')

    plt.tight_layout()

    if output_path:
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        print(f"Plot saved to: {output_path}")

    return fig, ax


def main():
    parser = argparse.ArgumentParser(description='Plot cumulative iterations from PLSTss log file')
    parser.add_argument('log_file', type=str, help='Path to log file')
    parser.add_argument('-o', '--output', type=str, help='Output plot file path')
    parser.add_argument('--bar', action='store_true', help='Also create bar plot of iterations per step')
    parser.add_argument('--show', action='store_true', help='Show plots interactively')

    args = parser.parse_args()

    log_path = Path(args.log_file)
    if not log_path.exists():
        print(f"Error: Log file not found: {log_path}")
        return 1

    try:
        # Parse log file
        print(f"Parsing log file: {log_path}")
        steps, iterations, cumulative, total_steps = parse_log_file(log_path)

        print(f"Found {len(steps)} loading steps")
        print(f"Total iterations: {cumulative[-1]}")
        print(f"Average iterations per step: {cumulative[-1]/len(steps):.2f}")

        # Determine output paths
        if args.output:
            output_path = Path(args.output)
        else:
            output_path = log_path.parent / f"{log_path.stem}_cumulative.png"

        # Create cumulative plot
        fig1, ax1 = plot_cumulative_iterations(steps, iterations, cumulative, total_steps, output_path)

        # Create bar plot if requested
        if args.bar:
            bar_output = output_path.parent / f"{output_path.stem}_bar.png"
            fig2, ax2 = plot_iterations_per_step(steps, iterations, bar_output)

        # Show plots if requested
        if args.show:
            plt.show()

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
            sys.argv.extend([default_log, '--show'])
            print(f"Using default log file: {default_log}")

    sys.exit(main())