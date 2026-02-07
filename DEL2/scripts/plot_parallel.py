"""
plot_parallel_io.py — Visualization for parallel I/O benchmark comparison.

Compares sequential vs parallel I/O approaches for SpMV matrix loading.

Dependencies: pandas, numpy, matplotlib

Usage:
    python plot_parallel_io.py

Expects CSVs in ./results/parallel_io/:
    strong_scaling_sequential.csv
    strong_scaling_parallel.csv
    io_comparison.csv
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------
# Directory layout
# ---------------------------------------------------------------------------
RESULTS_DIR = "../results/parallel_io"
PLOTS_DIR   = "../plots/parallel_io_reading"

# Create output directory
os.makedirs(PLOTS_DIR, exist_ok=True)

# ---------------------------------------------------------------------------
# Color scheme
# ---------------------------------------------------------------------------
IO_COLOURS = {
    "sequential": "#E63946",
    "parallel":   "#06A77D"
}
IO_MARKERS = {
    "sequential": "o",
    "parallel":   "s"
}

# ===========================================================================
# 1. DATA LOADING & PREPROCESSING
# ===========================================================================
def load_csv(path):
    """Read CSV and strip whitespace from string columns."""
    df = pd.read_csv(path)
    df.columns = [c.strip() for c in df.columns]
    for c in df.select_dtypes(include="object").columns:
        df[c] = df[c].str.strip()
    return df


def calculate_90th_percentile(df):
    """
    Compute 90th percentile across ranks for each (num_procs, io_type, run, matrix).
    This captures the tail latency from the slowest rank.
    """
    group_cols = ["num_procs", "io_type", "run", "matrix"]

    p90 = (
        df.groupby(group_cols)
        .agg(
            elapsed_time = ("elapsed_time",  lambda x: np.percentile(x, 90)),
            comm_time    = ("comm_time",     lambda x: np.percentile(x, 90)),
            io_time      = ("io_time",       lambda x: np.percentile(x, 90)),
        )
        .reset_index()
    )
    p90["comp_time"] = (p90["elapsed_time"] - p90["comm_time"]).clip(lower=0.0)
    return p90


def aggregate_over_runs(p90_df):
    """
    Collapse repeated runs → median of p90 values.
    """
    agg = (
        p90_df.groupby(["num_procs", "io_type", "matrix"])
        .agg(
            elapsed_time = ("elapsed_time", "median"),
            io_time      = ("io_time",      "median"),
        )
        .reset_index()
    )
    return agg


# ===========================================================================
# 2. I/O TIME COMPARISON PLOTS
# ===========================================================================
def plot_io_time_comparison(io_comparison_csv):
    """
    Plot absolute I/O times: sequential vs parallel across all matrices.
    One subplot per matrix.
    """
    df = pd.read_csv(io_comparison_csv)
    matrices = sorted(df["matrix"].unique())
    
    n_matrices = len(matrices)
    fig, axes = plt.subplots(1, n_matrices, figsize=(7 * n_matrices, 6), squeeze=False)
    
    for idx, matrix in enumerate(matrices):
        ax = axes[0, idx]
        sub = df[df["matrix"] == matrix]
        
        for io_type in ["sequential", "parallel"]:
            chunk = sub[sub["io_type"] == io_type].sort_values("num_procs")
            if chunk.empty:
                continue
                
            procs = chunk["num_procs"].values
            times = chunk["io_time_seconds"].values
            
            ax.plot(procs, times,
                    color=IO_COLOURS[io_type],
                    marker=IO_MARKERS[io_type],
                    linewidth=2.5, markersize=9,
                    label=io_type.capitalize())
        
        ax.set_xscale("log", base=2)
        ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f"{int(x)}"))
        ax.set_xlabel("Number of Processes", fontsize=12, fontweight="bold")
        if idx == 0:
            ax.set_ylabel("I/O Time (seconds)", fontsize=12, fontweight="bold")
        ax.set_title(f"{matrix}", fontsize=13, fontweight="bold")
        ax.legend(fontsize=10)
        ax.grid(True, alpha=0.3)
    
    fig.suptitle("I/O Time Comparison: Sequential vs Parallel", 
                 fontsize=15, fontweight="bold", y=1.02)
    plt.tight_layout()
    plt.savefig(os.path.join(PLOTS_DIR, "io_time_comparison.png"),
                dpi=300, bbox_inches="tight")
    plt.close()
    print("    saved  io_time_comparison.png")


def plot_io_speedup_per_matrix(io_comparison_csv):
    """
    Plot I/O time speedup: sequential vs parallel, one plot per matrix.
    Speedup = T_sequential / T_parallel
    """
    df = pd.read_csv(io_comparison_csv)
    matrices = sorted(df["matrix"].unique())

    for matrix in matrices:
        sub = df[df["matrix"] == matrix]
        
        # Pivot: rows = num_procs, columns = io_type
        pivot = sub.pivot(index="num_procs", columns="io_type", values="io_time_seconds")
        
        if "sequential" not in pivot.columns or "parallel" not in pivot.columns:
            continue
            
        procs = pivot.index.values
        speedup = pivot["sequential"] / pivot["parallel"]
        
        fig, ax = plt.subplots(figsize=(10, 7))
        
        ax.plot(procs, speedup, 
                marker="D", linewidth=2.5, markersize=10,
                color="#2E86AB", label="I/O Speedup")
        ax.axhline(1.0, color="gray", linestyle="--", linewidth=2, 
                   label="No speedup", alpha=0.7)
        
        ax.set_xscale("log", base=2)
        ax.set_yscale("log", base=2)
        ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f"{int(x)}"))
        ax.set_xlabel("Number of Processes", fontsize=13, fontweight="bold")
        ax.set_ylabel("I/O Speedup (x)", fontsize=13, fontweight="bold")
        ax.set_title(f"I/O Speedup: Parallel vs Sequential ({matrix})",
                     fontsize=14, fontweight="bold")
        ax.legend(fontsize=11)
        ax.grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig(os.path.join(PLOTS_DIR, f"io_speedup_{matrix}.png"),
                    dpi=300, bbox_inches="tight")
        plt.close()
        print(f"    saved  io_speedup_{matrix}.png")


# ===========================================================================
# 3. TOTAL TIME BREAKDOWN (I/O + SpMV)
# ===========================================================================
def plot_io_spmv_breakdown(strong_agg):
    """
    Stacked bar chart: I/O time + SpMV time.
    Shows total time AND where the benefit comes from.
    This is the KEY plot showing total execution time impact.
    """
    matrices = sorted(strong_agg["matrix"].unique())
    
    for matrix in matrices:
        sub = strong_agg[strong_agg["matrix"] == matrix]
        
        fig, ax = plt.subplots(figsize=(12, 7))
        
        procs_sorted = sorted(sub["num_procs"].unique())
        n_procs = len(procs_sorted)
        
        io_types = ["sequential", "parallel"]
        bar_w = 0.35
        
        for ti, io_type in enumerate(io_types):
            chunk = sub[sub["io_type"] == io_type].set_index("num_procs")
            chunk = chunk.reindex(procs_sorted)
            
            io_time_ms = chunk["io_time"].values * 1000
            spmv_time_ms = chunk["elapsed_time"].values * 1000  # SpMV only
            
            x = np.arange(n_procs) + (ti - 0.5) * bar_w
            
            # Apply hatching to BOTH layers for parallel I/O
            hatch_pattern = "///" if io_type == "parallel" else None
            
            # I/O layer (bottom)
            ax.bar(x, io_time_ms, width=bar_w,
                   color="#FFA500", alpha=0.85,
                   edgecolor="white", linewidth=0.7,
                   hatch=hatch_pattern,
                   label=f"I/O" if ti == 0 else None)
            
            # SpMV layer (top)
            ax.bar(x, spmv_time_ms, width=bar_w, bottom=io_time_ms,
                   color="#4169E1", alpha=0.85,
                   edgecolor="white", linewidth=0.7,
                   hatch=hatch_pattern,
                   label=f"SpMV" if ti == 0 else None)
        
        ax.set_xticks(range(n_procs))
        ax.set_xticklabels([str(int(p)) for p in procs_sorted])
        ax.set_xlabel("Number of Processes", fontsize=13, fontweight="bold")
        ax.set_ylabel("Time (ms)", fontsize=13, fontweight="bold")
        ax.set_title(f"Total Execution Time: I/O + SpMV ({matrix})",
                     fontsize=14, fontweight="bold")
        
        # Add legend showing both I/O types and components
        from matplotlib.patches import Patch
        legend_elements = [
            Patch(facecolor='#FFA500', label='I/O Time'),
            Patch(facecolor='#4169E1', label='SpMV Time'),
            Patch(facecolor='white', edgecolor='black', label='Sequential', hatch=''),
            Patch(facecolor='white', edgecolor='black', label='Parallel', hatch='///')
        ]
        ax.legend(handles=legend_elements, fontsize=10, loc='upper left')
        ax.grid(True, alpha=0.3, axis="y")
        
        plt.tight_layout()
        plt.savefig(os.path.join(PLOTS_DIR, f"total_time_breakdown_{matrix}.png"),
                    dpi=300, bbox_inches="tight")
        plt.close()
        print(f"    saved  total_time_breakdown_{matrix}.png")

def plot_io_fraction(strong_agg):
    """
    Show I/O time as percentage of total execution time.
    Demonstrates how parallel I/O reduces the I/O bottleneck.
    """
    matrices = sorted(strong_agg["matrix"].unique())
    
    for matrix in matrices:
        sub = strong_agg[strong_agg["matrix"] == matrix]
        
        fig, ax = plt.subplots(figsize=(10, 7))
        
        for io_type in ["sequential", "parallel"]:
            chunk = sub[sub["io_type"] == io_type].sort_values("num_procs")
            if chunk.empty:
                continue
            
            procs = chunk["num_procs"].values
            # Total time = I/O + SpMV
            total_time = chunk["io_time"] + chunk["elapsed_time"]
            io_fraction = (chunk["io_time"] / total_time) * 100
            
            ax.plot(procs, io_fraction,
                    color=IO_COLOURS[io_type],
                    marker=IO_MARKERS[io_type],
                    linewidth=2.5, markersize=9,
                    label=io_type.capitalize())
        
        ax.set_xscale("log", base=2)
        ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f"{int(x)}"))
        ax.set_xlabel("Number of Processes", fontsize=13, fontweight="bold")
        ax.set_ylabel("I/O Time / Total Time (%)", fontsize=13, fontweight="bold")
        ax.set_title(f"I/O Overhead as Percentage of Total Time ({matrix})",
                     fontsize=14, fontweight="bold")
        ax.legend(fontsize=11)
        ax.grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig(os.path.join(PLOTS_DIR, f"io_fraction_{matrix}.png"),
                    dpi=300, bbox_inches="tight")
        plt.close()
        print(f"    saved  io_fraction_{matrix}.png")


# ===========================================================================
# 4. MAIN
# ===========================================================================
def main():
    print("=" * 70)
    print("  Parallel I/O Benchmark — Focused Visualization")
    print("=" * 70)

    # Load CSVs
    print("\nLoading results...")
    strong_seq = load_csv(os.path.join(RESULTS_DIR, "strong_scaling_sequential.csv"))
    strong_par = load_csv(os.path.join(RESULTS_DIR, "strong_scaling_parallel.csv"))
    
    # Combine for unified processing
    strong_df = pd.concat([strong_seq, strong_par], ignore_index=True)

    # 90th percentile reduction
    print("Computing 90th percentile across ranks...")
    strong_90 = calculate_90th_percentile(strong_df)

    # Aggregate over runs
    strong_agg = aggregate_over_runs(strong_90)

    matrices = sorted(strong_agg["matrix"].unique())
    print(f"  Matrices found: {matrices}")
    print(f"  I/O types: {sorted(strong_agg['io_type'].unique())}")
    print(f"  Proc counts: {sorted(strong_agg['num_procs'].unique())}\n")

    # ================================================================
    # PLOT GENERATION
    # ================================================================
    print("-" * 70)
    print("  Generating Plots")
    print("-" * 70)
    
    # I/O Time Comparison
    io_csv = os.path.join(RESULTS_DIR, "io_comparison.csv")
    if os.path.exists(io_csv):
        print("\n[1/4] I/O Time Comparison...")
        plot_io_time_comparison(io_csv)
        
        print("[2/4] I/O Speedup (per matrix)...")
        plot_io_speedup_per_matrix(io_csv)
    else:
        print("  [WARN] io_comparison.csv not found, skipping I/O time plots")

    # Total Time Breakdown (I/O + SpMV)
    print("[3/4] Total Time Breakdown (I/O + SpMV)...")
    plot_io_spmv_breakdown(strong_agg)
    
    # I/O Fraction
    print("[4/4] I/O Overhead Fraction...")
    plot_io_fraction(strong_agg)

    # ================================================================
    print("\n" + "=" * 70)
    print("  All plots complete!")
    print("=" * 70)
    print(f"\nOutput directory: {PLOTS_DIR}/")
    print("\nGenerated plots:")
    for f in sorted(os.listdir(PLOTS_DIR)):
        if f.endswith('.png'):
            print(f"  - {f}")
    print()


if __name__ == "__main__":
    main()