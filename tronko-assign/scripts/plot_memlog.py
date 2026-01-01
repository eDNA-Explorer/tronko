#!/usr/bin/env python3
"""
Plot memory usage from tronko-assign TSV memory logs.
Usage: python plot_memlog.py log1.tsv [log2.tsv ...] [-o output.png]
"""
import sys
import argparse

def load_memlog(path):
    """Load a memory log TSV file."""
    try:
        import pandas as pd
        return pd.read_csv(path, sep='\t', comment='#',
                           names=['wall_time', 'phase', 'rss_mb', 'vm_mb',
                                  'peak_rss_mb', 'cpu_user', 'cpu_sys', 'extra'],
                           skiprows=1)  # Skip header row
    except ImportError:
        print("Error: pandas is required. Install with: pip install pandas")
        sys.exit(1)

def plot_logs(log_files, output='memory_comparison.png'):
    """Plot RSS over time for multiple log files."""
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("Error: matplotlib is required. Install with: pip install matplotlib")
        sys.exit(1)

    fig, axes = plt.subplots(2, 1, figsize=(12, 8))

    for path in log_files:
        try:
            df = load_memlog(path)
            label = path.split('/')[-1].replace('.tsv', '')

            # RSS over time
            axes[0].plot(df['wall_time'], df['rss_mb'], label=label, marker='.', markersize=4)
        except Exception as e:
            print(f"Warning: Could not load {path}: {e}")
            continue

    axes[0].set_xlabel('Time (s)')
    axes[0].set_ylabel('RSS (MB)')
    axes[0].legend()
    axes[0].set_title('Memory Usage Over Time')
    axes[0].grid(True, alpha=0.3)

    # Phase comparison (bar chart of key phases)
    key_phases = ['STARTUP', 'REFERENCE_LOADED', 'BWA_INDEX', 'THREADS_ALLOCATED']
    phase_data = {}

    for path in log_files:
        try:
            df = load_memlog(path)
            label = path.split('/')[-1].replace('.tsv', '')
            for phase in key_phases:
                phase_rows = df[df['phase'] == phase]
                if not phase_rows.empty:
                    if phase not in phase_data:
                        phase_data[phase] = {}
                    phase_data[phase][label] = phase_rows['rss_mb'].iloc[0]
        except:
            continue

    if phase_data:
        import pandas as pd
        phase_df = pd.DataFrame(phase_data).T
        phase_df.plot(kind='bar', ax=axes[1], alpha=0.7)
        axes[1].set_ylabel('RSS (MB)')
        axes[1].set_title('Memory by Phase')
        axes[1].legend()
        axes[1].tick_params(axis='x', rotation=45)

    plt.tight_layout()
    plt.savefig(output, dpi=150)
    print(f"Saved plot to {output}")

def main():
    parser = argparse.ArgumentParser(description='Plot memory usage from tronko-assign TSV logs')
    parser.add_argument('logs', nargs='+', help='TSV log files to plot')
    parser.add_argument('-o', '--output', default='memory_comparison.png', help='Output image file')
    args = parser.parse_args()

    plot_logs(args.logs, args.output)

if __name__ == '__main__':
    main()
