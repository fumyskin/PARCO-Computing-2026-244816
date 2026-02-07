"""
plot_results.py — Full benchmark visualisation for parallel SpMV.

Dependencies (all standard / PyPI-mainstream):
    pandas, numpy, matplotlib   (NO scipy, NO seaborn)

Usage:
    python plot_results.py

Expects the four CSVs produced by run_benchmark.sh to live in  ./results/ :
    strong_scaling_all.csv
    weak_scaling_all.csv
    strong_scaling_hybrid.csv
    weak_scaling_hybrid.csv
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------
# Directory layout  (mirrors the existing structure)
# ---------------------------------------------------------------------------
RESULTS_DIR = "../results/1D_and_2D"
PLOTS_DIR   = "../plots/1D_and_2D"

SUBDIRS = {
    "data":                os.path.join(PLOTS_DIR, "data_reduction"),
    "strong_speedup":      os.path.join(PLOTS_DIR, "speedup_strong"),
    "strong_comp_comm":    os.path.join(PLOTS_DIR, "comp_vs_comm_strong"),
    "strong_comm_volume":  os.path.join(PLOTS_DIR, "comm_volume_strong"),
    "weak":                os.path.join(PLOTS_DIR, "weak_scaling"),
    "comparison":          os.path.join(PLOTS_DIR, "comparison"),
}

for d in SUBDIRS.values():
    os.makedirs(d, exist_ok=True)

# ---------------------------------------------------------------------------
# Palette constants  —  kept close to the original code's colours
# ---------------------------------------------------------------------------
DIST_COLOURS = {"1D": "#2E86AB", "2D": "#E63946", "2D_cyclic": "#06A77D"}
DIST_MARKERS = {"1D": "o",       "2D": "s",       "2D_cyclic": "D"}
COMP_COLOUR  = "steelblue"
COMM_COLOUR  = "coral"

DIST_ORDER   = ["1D", "2D", "2D_cyclic"]   # canonical draw-order


# ===========================================================================
# 1.  DATA LOADING  &  90th-PERCENTILE REDUCTION
# ===========================================================================
def load_csv(path):
    """Read one CSV, strip any surrounding whitespace from string cols."""
    df = pd.read_csv(path)
    df.columns = [c.strip() for c in df.columns]
    for c in df.select_dtypes(include="object").columns:
        df[c] = df[c].str.strip()
    return df


def calculate_90th_percentile(df):
    """
    Correct p90 reduction:
        group by  (num_procs, distribution, run, [matrix if present])
        → take the 90th percentile of elapsed_time / comm_time  ACROSS RANKS
          (this captures the slowest-rank tail)
        Then derive comp_time.

    The original code grouped by (num_procs, rank), which computed p90 over
    repeated runs *per rank* and therefore never captured inter-rank imbalance.
    """
    group_cols = ["num_procs", "distribution", "run"]
    if "matrix" in df.columns:
        group_cols.append("matrix")

    p90 = (
        df.groupby(group_cols)
        .agg(
            elapsed_time    =("elapsed_time",    lambda x: np.percentile(x, 90)),
            comm_time       =("comm_time",       lambda x: np.percentile(x, 90)),
            local_nz        =("local_nz",        "mean"),
            ghost_entries   =("ghost_entries",   "mean"),
            local_flops     =("local_flops",     "mean"),
        )
        .reset_index()
    )
    p90["comp_time"] = (p90["elapsed_time"] - p90["comm_time"]).clip(lower=0.0)
    return p90


def aggregate_over_runs(p90_df, extra_group=None):
    """
    Collapse the 10 repeated runs → median of the p90 values.
    Also keep min / mean / max of ghost_entries and local_nz across runs
    for communication-volume and load-balance tables.

    extra_group: list of additional columns beyond (num_procs, distribution).
    """
    base = ["num_procs", "distribution"]
    if extra_group:
        base = base + extra_group

    agg = (
        p90_df.groupby(base)
        .agg(
            elapsed_time          =("elapsed_time",   "median"),
            comm_time             =("comm_time",      "median"),
            comp_time             =("comp_time",      "median"),
            ghost_entries_min     =("ghost_entries",  "min"),
            ghost_entries_mean    =("ghost_entries",  "mean"),
            ghost_entries_max     =("ghost_entries",  "max"),
            local_nz_min          =("local_nz",       "min"),
            local_nz_mean         =("local_nz",       "mean"),
            local_nz_max          =("local_nz",       "max"),
        )
        .reset_index()
    )
    return agg


# ===========================================================================
# 2.  AMDAHL  &  GUSTAFSON  (grid-search only — no scipy)
# ===========================================================================
def amdahl_speedup(p, parallel_fraction):
    """S(p) = 1 / ( (1-pf) + pf/p )"""
    s = 1.0 - parallel_fraction
    return 1.0 / (s + parallel_fraction / np.asarray(p, dtype=float))


def estimate_parallel_fraction(procs, speedups, n_steps=2000):
    """
    Grid-search over parallel_fraction in [0.01, 0.9999] minimising
    sum-of-squared-errors between empirical speedup and Amdahl prediction.
    """
    best_pf   = 0.95
    best_err  = float("inf")
    for pf in np.linspace(0.01, 0.9999, n_steps):
        pred = amdahl_speedup(procs, pf)
        err  = float(np.sum((speedups - pred) ** 2))
        if err < best_err:
            best_err = err
            best_pf  = pf
    return best_pf


def gustafson_scaled_speedup(p, serial_fraction):
    """
    Gustafson scaled speedup:   S(p) = s + p*(1-s)   where s = serial fraction.
    Equivalently  S(p) = p - (p-1)*s.
    """
    return np.asarray(p, dtype=float) - (np.asarray(p, dtype=float) - 1.0) * serial_fraction


def estimate_serial_fraction_weak(procs, scaled_speedups, n_steps=2000):
    """
    Grid-search serial_fraction in [0.001, 0.5] minimising SSE against
    empirical *scaled* speedup  (= p * T(1)/T(p) ).
    """
    best_s   = 0.05
    best_err = float("inf")
    for s in np.linspace(0.001, 0.50, n_steps):
        pred = gustafson_scaled_speedup(procs, s)
        err  = float(np.sum((scaled_speedups - pred) ** 2))
        if err < best_err:
            best_err = err
            best_s   = s
    return best_s


# ===========================================================================
# 3.  HELPER: baseline T(1) lookup  (handles missing 1-proc 2D runs)
# ===========================================================================
def get_baseline_time(agg_df, distribution, matrix=None):
    """
    Return T(1) for the given distribution.
    If distribution has no 1-proc row (2D/2D_cyclic are skipped at 1 proc
    in the bash script), fall back to the 1D baseline.
    """
    mask = (agg_df["distribution"] == distribution) & (agg_df["num_procs"] == 1)
    if matrix is not None:
        mask = mask & (agg_df["matrix"] == matrix)

    row = agg_df.loc[mask]
    if row.empty:
        # fall back to 1D
        mask2 = (agg_df["distribution"] == "1D") & (agg_df["num_procs"] == 1)
        if matrix is not None:
            mask2 = mask2 & (agg_df["matrix"] == matrix)
        row = agg_df.loc[mask2]

    if row.empty:
        return None
    return float(row["elapsed_time"].iloc[0])


# ===========================================================================
# 4.  PLOTS  ──  strong scaling
# ===========================================================================
# ---------- PLOT 1: Strong Speedup per matrix (all distributions + Amdahl) --
def plot_strong_speedup(strong_agg, label="MPI"):
    """One figure per matrix.  Three coloured lines (1D / 2D / 2D_cyclic),
    ideal linear, and one Amdahl envelope fitted on the 1D curve."""
    matrices = sorted(strong_agg["matrix"].unique())

    for matrix in matrices:
        fig, ax = plt.subplots(figsize=(10, 7))
        sub = strong_agg[strong_agg["matrix"] == matrix]

        # We fit Amdahl on 1D (most data points); overlay curve applies to all.
        t1_1d = get_baseline_time(sub, "1D", matrix)

        for dist in DIST_ORDER:
            chunk = sub[sub["distribution"] == dist].sort_values("num_procs")
            if chunk.empty:
                continue

            t1 = get_baseline_time(sub, dist, matrix)
            if t1 is None:
                continue

            procs   = chunk["num_procs"].values.astype(float)
            speedup = t1 / chunk["elapsed_time"].values

            ax.plot(procs, speedup,
                    color=DIST_COLOURS[dist], marker=DIST_MARKERS[dist],
                    linewidth=2.5, markersize=9, label=dist, zorder=3)

        # ── Amdahl reference (fit on 1D) ──
        if t1_1d is not None:
            ref = sub[sub["distribution"] == "1D"].sort_values("num_procs")
            if not ref.empty:
                p_ref = ref["num_procs"].values.astype(float)
                s_ref = t1_1d / ref["elapsed_time"].values
                pf    = estimate_parallel_fraction(p_ref, s_ref)

                p_cont = np.logspace(0, np.log2(float(sub["num_procs"].max())),
                                     200, base=2)
                ax.plot(p_cont, amdahl_speedup(p_cont, pf),
                        color="#888888", linestyle="-.", linewidth=2,
                        label=f"Amdahl (pf={pf:.3f})", zorder=2)

        # ── Ideal linear ──
        p_max = float(sub["num_procs"].max())
        ax.plot([1, p_max], [1, p_max], "k--", linewidth=2, label="Ideal (linear)", zorder=1)

        ax.set_xscale("log", base=2)
        ax.set_yscale("log", base=2)
        ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f"{int(x)}"))
        ax.set_xlabel("Number of Processes", fontsize=13, fontweight="bold")
        ax.set_ylabel("Speedup (x)", fontsize=13, fontweight="bold")
        ax.set_title(f"Strong Scaling – Speedup  [{label}]  ({matrix})",
                     fontsize=14, fontweight="bold")
        ax.legend(fontsize=11)
        ax.grid(True, alpha=0.3)

        plt.tight_layout()
        plt.savefig(os.path.join(SUBDIRS["strong_speedup"],
                                 f"speedup_{label}_{matrix}.png"),
                    dpi=300, bbox_inches="tight")
        plt.close()
        print(f"    saved  speedup_{label}_{matrix}.png")


# ---------- PLOT 2: Strong Comp vs Comm  (grouped stacked bars) -----------
def plot_strong_comp_comm(strong_agg, label="MPI"):
    """One figure per matrix.  For every proc-count we draw one narrow
    stacked bar per distribution (comp bottom, comm on top)."""
    matrices = sorted(strong_agg["matrix"].unique())

    for matrix in matrices:
        fig, ax = plt.subplots(figsize=(11, 7))
        sub = strong_agg[strong_agg["matrix"] == matrix].copy()

        procs_sorted = sorted(sub["num_procs"].unique())
        n_procs      = len(procs_sorted)

        dists_present = [d for d in DIST_ORDER if d in sub["distribution"].values]
        n_dists       = len(dists_present)
        bar_w         = 0.7 / max(n_dists, 1)

        for di, dist in enumerate(dists_present):
            chunk = sub[sub["distribution"] == dist].set_index("num_procs")
            chunk = chunk.reindex(procs_sorted)   # NaN where no data

            comp_ms = (chunk["comp_time"].values * 1000)
            comm_ms = (chunk["comm_time"].values * 1000)

            x = np.arange(n_procs) + (di - n_dists / 2.0 + 0.5) * bar_w

            ax.bar(x, comp_ms, width=bar_w,
                   color=COMP_COLOUR, alpha=0.85, edgecolor="white", linewidth=0.7,
                   label="Computation" if di == 0 else None)
            ax.bar(x, comm_ms, width=bar_w, bottom=comp_ms,
                   color=COMM_COLOUR, alpha=0.85, edgecolor="white", linewidth=0.7,
                   label="Communication" if di == 0 else None,
                   hatch="///" if dist != "1D" else None)


        ax.set_xticks(range(n_procs))
        ax.set_xticklabels([str(int(p)) for p in procs_sorted])
        
        # Add distribution labels below each process count group
        y_min = ax.get_ylim()[0]
        label_y = y_min - (ax.get_ylim()[1] - y_min) * 0.05  # 5% below bottom
        
        for pi in range(n_procs):
            # Group center for this process count
            group_center = pi
            for di, dist in enumerate(dists_present):
                x_pos = group_center + (di - n_dists / 2.0 + 0.5) * bar_w
                ax.text(x_pos, label_y, dist, ha="center", va="top",
                        fontsize=7, color=DIST_COLOURS[dist], fontweight="bold",
                        rotation=0)
        
        ax.set_xlabel("Number of Processes", fontsize=13, fontweight="bold")
        ax.set_ylabel("Time (ms)", fontsize=13, fontweight="bold")
        ax.set_title(f"Computation vs Communication  [{label}]  ({matrix})",
                     fontsize=14, fontweight="bold")
        ax.legend(fontsize=11)
        ax.grid(True, alpha=0.3, axis="y")

        plt.tight_layout()
        plt.savefig(os.path.join(SUBDIRS["strong_comp_comm"],
                                 f"comp_vs_comm_{label}_{matrix}.png"),
                    dpi=300, bbox_inches="tight")
        plt.close()
        print(f"    saved  comp_vs_comm_{label}_{matrix}.png")


# ---------- PLOT 2b: MPI vs Hybrid side-by-side (comp/comm) ---------------
def plot_strong_comp_comm_mpi_vs_hybrid(strong_agg_mpi, strong_agg_hyb):
    """For each matrix: grouped bars comparing MPI and Hybrid per proc count,
    keeping all distributions visible via hatching."""
    matrices = sorted(strong_agg_mpi["matrix"].unique())

    for matrix in matrices:
        mpi = strong_agg_mpi[strong_agg_mpi["matrix"] == matrix]
        hyb = strong_agg_hyb[strong_agg_hyb["matrix"] == matrix]
        if mpi.empty or hyb.empty:
            continue

        # only proc counts present in both
        common_procs = sorted(
            set(mpi["num_procs"].values) & set(hyb["num_procs"].values)
        )
        if not common_procs:
            continue

        fig, ax = plt.subplots(figsize=(11, 7))
        n_procs = len(common_procs)

        # two mega-bars per proc-count: MPI on the left, Hybrid on the right
        mega_w = 0.35

        for pi, p in enumerate(common_procs):
            mpi_row = mpi[mpi["num_procs"] == p]
            hyb_row = hyb[hyb["num_procs"] == p]
            if mpi_row.empty or hyb_row.empty:
                continue

            # take first distribution row as representative (already aggregated)
            comp_mpi = float(mpi_row["comp_time"].iloc[0]) * 1000
            comm_mpi = float(mpi_row["comm_time"].iloc[0]) * 1000
            comp_hyb = float(hyb_row["comp_time"].iloc[0]) * 1000
            comm_hyb = float(hyb_row["comm_time"].iloc[0]) * 1000

            x_mpi = pi - mega_w / 2.0
            x_hyb = pi + mega_w / 2.0

            ax.bar(x_mpi, comp_mpi, width=mega_w, color="steelblue", alpha=0.85,
                   edgecolor="white", linewidth=0.7,
                   label="MPI – Computation" if pi == 0 else None)
            ax.bar(x_mpi, comm_mpi, width=mega_w, bottom=comp_mpi,
                   color="coral", alpha=0.85, edgecolor="white", linewidth=0.7,
                   label="MPI – Communication" if pi == 0 else None)

            ax.bar(x_hyb, comp_hyb, width=mega_w, color="seagreen", alpha=0.85,
                   edgecolor="white", linewidth=0.7,
                   label="Hybrid – Computation" if pi == 0 else None)
            ax.bar(x_hyb, comm_hyb, width=mega_w, bottom=comp_hyb,
                   color="goldenrod", alpha=0.85, edgecolor="white", linewidth=0.7,
                   label="Hybrid – Communication" if pi == 0 else None)

        ax.set_xticks(range(n_procs))
        ax.set_xticklabels([str(int(p)) for p in common_procs])
        ax.set_xlabel("Number of Processes", fontsize=13, fontweight="bold")
        ax.set_ylabel("Time (ms)", fontsize=13, fontweight="bold")
        ax.set_title(f"Comp vs Comm – MPI vs MPI+OpenMP  ({matrix})",
                     fontsize=14, fontweight="bold")
        ax.legend(fontsize=10)
        ax.grid(True, alpha=0.3, axis="y")

        plt.tight_layout()
        plt.savefig(os.path.join(SUBDIRS["strong_comp_comm"],
                                 f"comp_vs_comm_mpi_vs_hybrid_{matrix}.png"),
                    dpi=300, bbox_inches="tight")
        plt.close()
        print(f"    saved  comp_vs_comm_mpi_vs_hybrid_{matrix}.png")


# ---------- PLOT 3: Communication-Volume heatmap table --------------------
def plot_comm_volume_heatmap(strong_agg, label="MPI"):
    """
    One figure per matrix.  Rows = distributions, columns = proc counts.
    Cell value = mean ghost_entries (summed across ranks already in the
    original raw CSV; here we have the per-rank mean from aggregation, so
    we multiply back by num_procs to recover total volume).
    """
    matrices = sorted(strong_agg["matrix"].unique())

    for matrix in matrices:
        sub = strong_agg[strong_agg["matrix"] == matrix].copy()
        # total ghost volume = mean_per_rank * num_procs
        sub["ghost_volume"] = sub["ghost_entries_mean"] * sub["num_procs"]

        procs_sorted = sorted(sub["num_procs"].unique())
        dists_present = [d for d in DIST_ORDER if d in sub["distribution"].values]

        n_rows = len(dists_present)
        n_cols = len(procs_sorted)

        # build the 2-D grid
        grid = np.full((n_rows, n_cols), np.nan)
        for ri, dist in enumerate(dists_present):
            for ci, p in enumerate(procs_sorted):
                cell = sub[(sub["distribution"] == dist) & (sub["num_procs"] == p)]
                if not cell.empty:
                    grid[ri, ci] = float(cell["ghost_volume"].iloc[0])

        # --- figure ---
        fig, ax = plt.subplots(figsize=(max(8, n_cols * 1.3), max(3.5, n_rows * 1.4)))

        im = ax.imshow(grid, aspect="auto", cmap="YlOrRd", origin="upper")

        # grid lines
        ax.set_xticks(np.arange(n_cols))
        ax.set_xticklabels([str(int(p)) for p in procs_sorted], fontsize=10)
        ax.set_yticks(np.arange(n_rows))
        ax.set_yticklabels(dists_present, fontsize=11)
        ax.set_xlabel("Number of Processes", fontsize=13, fontweight="bold")
        ax.set_title(f"Communication Volume  [{label}]  ({matrix})",
                     fontsize=14, fontweight="bold")

        # cell annotations
        vmax = np.nanmax(grid) if not np.all(np.isnan(grid)) else 1.0
        for ri in range(n_rows):
            for ci in range(n_cols):
                v = grid[ri, ci]
                if np.isnan(v):
                    continue
                colour = "white" if v > vmax * 0.55 else "black"
                ax.text(ci, ri, f"{int(v):,}", ha="center", va="center",
                        fontsize=9, color=colour, fontweight="bold")

        fig.colorbar(im, ax=ax, shrink=0.85, label="Total ghost entries")
        plt.tight_layout()
        plt.savefig(os.path.join(SUBDIRS["strong_comm_volume"],
                                 f"comm_volume_{label}_{matrix}.png"),
                    dpi=300, bbox_inches="tight")
        plt.close()
        print(f"    saved  comm_volume_{label}_{matrix}.png")

        # ── also dump the CSV table (kept from original) ──
        rows = []
        for ri, dist in enumerate(dists_present):
            row_d = sub[sub["distribution"] == dist].sort_values("num_procs")
            for _, r in row_d.iterrows():
                rows.append({
                    "distribution":       dist,
                    "num_procs":          int(r["num_procs"]),
                    "ghost_entries_min":  int(r["ghost_entries_min"]),
                    "ghost_entries_mean": round(float(r["ghost_entries_mean"]), 2),
                    "ghost_entries_max":  int(r["ghost_entries_max"]),
                    "ghost_volume_total": int(r["ghost_volume"]),
                    "load_imbalance":     round(float(r["ghost_entries_max"] /
                                                      max(r["ghost_entries_mean"], 1)), 3),
                })
        pd.DataFrame(rows).to_csv(
            os.path.join(SUBDIRS["strong_comm_volume"],
                         f"comm_volume_{label}_{matrix}.csv"),
            index=False)


# ---------- PLOT 4: Matrix comparison  (speedup + efficiency) -------------
def plot_strong_comparison(strong_agg, label="MPI"):
    """Two-panel: left = speedup, right = efficiency.  One curve per matrix,
    using the 1D distribution as the representative curve."""
    matrices = sorted(strong_agg["matrix"].unique())
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))

    for matrix in matrices:
        sub   = strong_agg[(strong_agg["matrix"] == matrix) &
                           (strong_agg["distribution"] == "1D")].sort_values("num_procs")
        if sub.empty:
            continue
        t1      = float(sub["elapsed_time"].iloc[0])
        procs   = sub["num_procs"].values.astype(float)
        speedup = t1 / sub["elapsed_time"].values

        ax1.plot(procs, speedup, marker="o", linewidth=2.5, markersize=8, label=matrix)
        ax2.plot(procs, (speedup / procs) * 100,
                 marker="s", linewidth=2.5, markersize=8, label=matrix)

    p_max = float(strong_agg["num_procs"].max())
    ax1.plot([1, p_max], [1, p_max], "k--", linewidth=2, label="Ideal")
    ax2.axhline(100, color="k", linestyle="--", linewidth=2, label="Ideal")

    for ax in (ax1, ax2):
        ax.set_xscale("log", base=2)
        ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f"{int(x)}"))
        ax.set_xlabel("Number of Processes", fontsize=13, fontweight="bold")
        ax.legend(fontsize=11)
        ax.grid(True, alpha=0.3)

    ax1.set_yscale("log", base=2)
    ax1.set_ylabel("Speedup (x)", fontsize=13, fontweight="bold")
    ax1.set_title(f"Strong Scaling – Speedup Comparison  [{label}]",
                  fontsize=14, fontweight="bold")

    ax2.set_ylim(0, 110)
    ax2.set_ylabel("Efficiency (%)", fontsize=13, fontweight="bold")
    ax2.set_title(f"Strong Scaling – Efficiency Comparison  [{label}]",
                  fontsize=14, fontweight="bold")

    plt.tight_layout()
    plt.savefig(os.path.join(SUBDIRS["comparison"],
                             f"strong_scaling_comparison_{label}.png"),
                dpi=300, bbox_inches="tight")
    plt.close()
    print(f"    saved  strong_scaling_comparison_{label}.png")


# ---------- PLOT 5: MPI vs Hybrid speedup per matrix ----------------------
def plot_strong_mpi_vs_hybrid(strong_agg_mpi, strong_agg_hyb):
    """Per-matrix figure, 1D distribution, two speedup curves."""
    matrices = sorted(strong_agg_mpi["matrix"].unique())

    for matrix in matrices:
        mpi = (strong_agg_mpi[(strong_agg_mpi["matrix"] == matrix) &
                              (strong_agg_mpi["distribution"] == "1D")]
               .sort_values("num_procs"))
        hyb = (strong_agg_hyb[(strong_agg_hyb["matrix"] == matrix) &
                              (strong_agg_hyb["distribution"] == "1D")]
               .sort_values("num_procs"))
        if mpi.empty or hyb.empty:
            continue

        common = sorted(set(mpi["num_procs"].values) & set(hyb["num_procs"].values))
        if not common:
            continue
        mpi = mpi[mpi["num_procs"].isin(common)]
        hyb = hyb[hyb["num_procs"].isin(common)]

        procs = np.array(common, dtype=float)
        t1_m  = float(mpi["elapsed_time"].iloc[0])
        t1_h  = float(hyb["elapsed_time"].iloc[0])

        fig, ax = plt.subplots(figsize=(10, 7))
        ax.plot(procs, t1_m / mpi["elapsed_time"].values,
                marker="o", linewidth=2.5, markersize=8,
                color="#1f77b4", label="MPI")
        ax.plot(procs, t1_h / hyb["elapsed_time"].values,
                marker="s", linewidth=2.5, markersize=8,
                color="#ff7f0e", label="MPI+OpenMP")
        ax.plot([1, procs[-1]], [1, procs[-1]], "k--", linewidth=2, label="Ideal")

        ax.set_xscale("log", base=2)
        ax.set_yscale("log", base=2)
        ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f"{int(x)}"))
        ax.set_xlabel("Number of Processes", fontsize=13, fontweight="bold")
        ax.set_ylabel("Speedup (x)", fontsize=13, fontweight="bold")
        ax.set_title(f"Strong Scaling – MPI vs MPI+OpenMP  ({matrix})",
                     fontsize=14, fontweight="bold")
        ax.legend(fontsize=11)
        ax.grid(True, alpha=0.3)

        plt.tight_layout()
        plt.savefig(os.path.join(SUBDIRS["comparison"],
                                 f"strong_mpi_vs_hybrid_{matrix}.png"),
                    dpi=300, bbox_inches="tight")
        plt.close()
        print(f"    saved  strong_mpi_vs_hybrid_{matrix}.png")


# ===========================================================================
# 5.  PLOTS  ──  weak scaling
# ===========================================================================
# ---------- PLOT 6: Weak Speedup + Gustafson  (all distributions) ----------
def plot_weak_speedup(weak_agg, label="MPI"):
    """
    Scaled speedup  =  p * T(1) / T(p).   Ideal = p  (linear).
    Gustafson fitted on the 1D curve; one overlay curve for reference.
    """
    fig, ax = plt.subplots(figsize=(10, 7))

    # 1D baseline for all distributions (consistent reference)
    row1 = weak_agg[(weak_agg["distribution"] == "1D") &
                    (weak_agg["num_procs"] == weak_agg["num_procs"].min())]
    if row1.empty:
        print("  [WARN] no 1-proc weak baseline — skipping weak speedup plot")
        plt.close()
        return
    t1_ref = float(row1["elapsed_time"].iloc[0])
    p_ref  = int(row1["num_procs"].iloc[0])

    p_max = float(weak_agg["num_procs"].max())
    p_cont = np.logspace(0, np.log2(p_max), 200, base=2)

    # Gustafson fit on 1D
    ref_1d = weak_agg[weak_agg["distribution"] == "1D"].sort_values("num_procs")
    if not ref_1d.empty:
        ps_1d = ref_1d["num_procs"].values.astype(float)
        ss_1d = ps_1d * (t1_ref / ref_1d["elapsed_time"].values)
        s_est = estimate_serial_fraction_weak(ps_1d, ss_1d)

        ax.plot(p_cont, gustafson_scaled_speedup(p_cont, s_est),
                color="#E63946", linestyle="-.", linewidth=2,
                label=f"Gustafson (s={s_est:.3f})", zorder=2)

    # Ideal  S = p
    ax.plot([1, p_max], [1, p_max], "k--", linewidth=2, label="Ideal (linear)", zorder=1)

    # Empirical curves per distribution
    for dist in DIST_ORDER:
        chunk = weak_agg[weak_agg["distribution"] == dist].sort_values("num_procs")
        if chunk.empty:
            continue
        procs   = chunk["num_procs"].values.astype(float)
        s_speedup = procs * (t1_ref / chunk["elapsed_time"].values)

        ax.plot(procs, s_speedup,
                color=DIST_COLOURS[dist], marker=DIST_MARKERS[dist],
                linewidth=2.5, markersize=10, label=dist, zorder=3)

    ax.set_xscale("log", base=2)
    ax.set_yscale("log", base=2)
    ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f"{int(x)}"))
    ax.set_xlabel("Number of Processes", fontsize=13, fontweight="bold")
    ax.set_ylabel("Scaled Speedup (x)", fontsize=13, fontweight="bold")
    ax.set_title(f"Weak Scaling – Scaled Speedup + Gustafson  [{label}]",
                 fontsize=14, fontweight="bold")
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(os.path.join(SUBDIRS["weak"], f"weak_scaling_speedup_{label}.png"),
                dpi=300, bbox_inches="tight")
    plt.close()
    print(f"    saved  weak_scaling_speedup_{label}.png")


# ---------- PLOT 7: Weak Efficiency ----------------------------------------
def plot_weak_efficiency(weak_agg, label="MPI"):
    fig, ax = plt.subplots(figsize=(10, 7))

    row1 = weak_agg[(weak_agg["distribution"] == "1D") &
                    (weak_agg["num_procs"] == weak_agg["num_procs"].min())]
    if row1.empty:
        plt.close()
        return
    t1_ref = float(row1["elapsed_time"].iloc[0])

    # Gustafson efficiency curve (reuse fit)
    ref_1d = weak_agg[weak_agg["distribution"] == "1D"].sort_values("num_procs")
    if not ref_1d.empty:
        ps_1d = ref_1d["num_procs"].values.astype(float)
        ss_1d = ps_1d * (t1_ref / ref_1d["elapsed_time"].values)
        s_est = estimate_serial_fraction_weak(ps_1d, ss_1d)

        p_max  = float(weak_agg["num_procs"].max())
        p_cont = np.logspace(0, np.log2(p_max), 200, base=2)
        gust_eff = (gustafson_scaled_speedup(p_cont, s_est) / p_cont) * 100
        ax.plot(p_cont, gust_eff, color="#E63946", linestyle="-.", linewidth=2,
                label=f"Gustafson (s={s_est:.3f})", zorder=2)

    ax.axhline(100, color="k", linestyle="--", linewidth=2, label="Ideal", zorder=1)

    for dist in DIST_ORDER:
        chunk = weak_agg[weak_agg["distribution"] == dist].sort_values("num_procs")
        if chunk.empty:
            continue
        procs = chunk["num_procs"].values.astype(float)
        eff   = (t1_ref / chunk["elapsed_time"].values) * 100

        ax.plot(procs, eff,
                color=DIST_COLOURS[dist], marker=DIST_MARKERS[dist],
                linewidth=2.5, markersize=10, label=dist, zorder=3)

    ax.set_xscale("log", base=2)
    ax.set_ylim(0, 110)
    ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f"{int(x)}"))
    ax.set_xlabel("Number of Processes", fontsize=13, fontweight="bold")
    ax.set_ylabel("Efficiency (%)", fontsize=13, fontweight="bold")
    ax.set_title(f"Weak Scaling – Efficiency  [{label}]",
                 fontsize=14, fontweight="bold")
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(os.path.join(SUBDIRS["weak"], f"weak_scaling_efficiency_{label}.png"),
                dpi=300, bbox_inches="tight")
    plt.close()
    print(f"    saved  weak_scaling_efficiency_{label}.png")


# ---------- PLOT 8: Weak Time Breakdown  (stacked bars, faceted by dist) --
def plot_weak_time_breakdown(weak_agg, label="MPI"):
    """One sub-plot per distribution.  Stacked comp/comm bars + total line."""
    dists_present = [d for d in DIST_ORDER if d in weak_agg["distribution"].values]
    n_d = len(dists_present)
    if n_d == 0:
        return

    fig, axes = plt.subplots(1, n_d, figsize=(6 * n_d, 6), squeeze=False)

    for di, dist in enumerate(dists_present):
        ax   = axes[0, di]
        chunk = weak_agg[weak_agg["distribution"] == dist].sort_values("num_procs")
        if chunk.empty:
            continue

        procs   = chunk["num_procs"].values
        comp_ms = chunk["comp_time"].values * 1000
        comm_ms = chunk["comm_time"].values * 1000
        tot_ms  = chunk["elapsed_time"].values * 1000

        x     = np.arange(len(procs))
        width = 0.55

        ax.bar(x, comp_ms, width, color=COMP_COLOUR, alpha=0.85,
               edgecolor="white", linewidth=0.7,
               label="Computation")
        ax.bar(x, comm_ms, width, bottom=comp_ms,
               color=COMM_COLOUR, alpha=0.85,
               edgecolor="white", linewidth=0.7,
               label="Communication")
        ax.plot(x, tot_ms, color="#6A4C93", marker="o", linewidth=2.5,
                markersize=8, label="Total (p90)", zorder=4)

        # baseline line (1-proc total)
        row1 = weak_agg[(weak_agg["distribution"] == "1D") &
                        (weak_agg["num_procs"] == weak_agg["num_procs"].min())]
        if not row1.empty:
            ax.axhline(float(row1["elapsed_time"].iloc[0]) * 1000,
                       color="gray", linestyle="--", linewidth=1.8, alpha=0.7,
                       label="Baseline (1 proc)")

        ax.set_xticks(x)
        ax.set_xticklabels([str(int(p)) for p in procs])
        ax.set_xlabel("Number of Processes", fontsize=13, fontweight="bold")
        if di == 0:
            ax.set_ylabel("Time (ms)", fontsize=13, fontweight="bold")
        ax.set_title(f"{dist}", fontsize=13, fontweight="bold")
        ax.legend(fontsize=9, loc="upper left")
        ax.grid(True, alpha=0.3, axis="y")

    fig.suptitle(f"Weak Scaling – Execution Time Breakdown  [{label}]",
                 fontsize=15, fontweight="bold", y=1.02)
    plt.tight_layout()
    plt.savefig(os.path.join(SUBDIRS["weak"],
                             f"weak_scaling_time_breakdown_{label}.png"),
                dpi=300, bbox_inches="tight")
    plt.close()
    print(f"    saved  weak_scaling_time_breakdown_{label}.png")


# ---------- PLOT 10: Weak Comp vs Comm  (one panel, all dists) ------------
def plot_weak_comp_comm(weak_agg, label="MPI"):
    fig, ax = plt.subplots(figsize=(10, 7))

    procs_sorted  = sorted(weak_agg["num_procs"].unique())
    n_procs       = len(procs_sorted)
    dists_present = [d for d in DIST_ORDER if d in weak_agg["distribution"].values]
    n_dists       = len(dists_present)
    bar_w         = 0.7 / max(n_dists, 1)

    for di, dist in enumerate(dists_present):
        chunk = weak_agg[weak_agg["distribution"] == dist].set_index("num_procs")
        chunk = chunk.reindex(procs_sorted)

        comp_ms = chunk["comp_time"].values * 1000
        comm_ms = chunk["comm_time"].values * 1000

        x = np.arange(n_procs) + (di - n_dists / 2.0 + 0.5) * bar_w
        ax.bar(x, comp_ms, width=bar_w,
               color=COMP_COLOUR, alpha=0.85, edgecolor="white", linewidth=0.7,
               label="Computation" if di == 0 else None)
        ax.bar(x, comm_ms, width=bar_w, bottom=comp_ms,
               color=COMM_COLOUR, alpha=0.85, edgecolor="white", linewidth=0.7,
               label="Communication" if di == 0 else None,
               hatch="///" if dist != "1D" else None)

    ax.set_xticks(range(n_procs))
    ax.set_xticklabels([str(int(p)) for p in procs_sorted])
    
    # Add distribution labels below each process count group
    y_min = ax.get_ylim()[0]
    label_y = y_min - (ax.get_ylim()[1] - y_min) * 0.05  # 5% below bottom
    
    for pi in range(n_procs):
        # Group center for this process count
        group_center = pi
        for di, dist in enumerate(dists_present):
            x_pos = group_center + (di - n_dists / 2.0 + 0.5) * bar_w
            ax.text(x_pos, label_y, dist, ha="center", va="top",
                    fontsize=7, color=DIST_COLOURS[dist], fontweight="bold",
                    rotation=0)
    
    ax.set_xlabel("Number of Processes", fontsize=13, fontweight="bold")
    ax.set_ylabel("Time (ms)", fontsize=13, fontweight="bold")
    ax.set_title(f"Weak Scaling – Computation vs Communication  [{label}]",
                 fontsize=14, fontweight="bold")
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3, axis="y")

    plt.tight_layout()
    plt.savefig(os.path.join(SUBDIRS["weak"],
                             f"weak_scaling_comp_vs_comm_{label}.png"),
                dpi=300, bbox_inches="tight")
    plt.close()
    print(f"    saved  weak_scaling_comp_vs_comm_{label}.png")


# ---------- PLOT 11: Weak Load Balance  -----------------------------------
def plot_weak_load_balance(weak_agg, label="MPI"):
    fig, ax = plt.subplots(figsize=(10, 7))

    # use 1D as representative (simplest story)
    sub = weak_agg[weak_agg["distribution"] == "1D"].sort_values("num_procs")
    if sub.empty:
        plt.close()
        return

    procs = sub["num_procs"].values.astype(float)
    ax.plot(procs, sub["local_nz_min"],  marker="v", linewidth=2.5, markersize=10,
            color="#E63946", label="Min NNZ")
    ax.plot(procs, sub["local_nz_mean"], marker="o", linewidth=2.5, markersize=10,
            color="#06FFA5", label="Avg NNZ")
    ax.plot(procs, sub["local_nz_max"],  marker="^", linewidth=2.5, markersize=10,
            color="#F77F00", label="Max NNZ")

    ax.set_xscale("log", base=2)
    ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f"{int(x)}"))
    ax.set_xlabel("Number of Processes", fontsize=13, fontweight="bold")
    ax.set_ylabel("Non-zeros per Rank", fontsize=13, fontweight="bold")
    ax.set_title(f"Weak Scaling – Load Balance (NNZ per Rank)  [{label}]",
                 fontsize=14, fontweight="bold")
    ax.legend(fontsize=12)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(os.path.join(SUBDIRS["weak"],
                             f"weak_scaling_load_balance_{label}.png"),
                dpi=300, bbox_inches="tight")
    plt.close()
    print(f"    saved  weak_scaling_load_balance_{label}.png")


# ---------- PLOT 12: Weak Comm Volume CSV table  --------------------------
def export_weak_comm_volume(weak_agg, label="MPI"):
    rows = []
    for _, r in weak_agg.sort_values(["distribution", "num_procs"]).iterrows():
        rows.append({
            "distribution":       r["distribution"],
            "num_procs":          int(r["num_procs"]),
            "ghost_entries_min":  int(r["ghost_entries_min"]),
            "ghost_entries_mean": round(float(r["ghost_entries_mean"]), 2),
            "ghost_entries_max":  int(r["ghost_entries_max"]),
            "load_imbalance":     round(float(r["ghost_entries_max"] /
                                              max(r["ghost_entries_mean"], 1)), 3),
        })
    out = os.path.join(SUBDIRS["weak"], f"comm_volume_weak_{label}.csv")
    pd.DataFrame(rows).to_csv(out, index=False)
    print(f"    saved  comm_volume_weak_{label}.csv")


# ---------- PLOT 9: MPI vs Hybrid weak (speedup + efficiency) -------------
def plot_weak_mpi_vs_hybrid(weak_agg_mpi, weak_agg_hyb):
    """Two-panel: scaled speedup (left) and efficiency (right)."""
    # restrict to 1D for a clean comparison
    mpi = weak_agg_mpi[weak_agg_mpi["distribution"] == "1D"].sort_values("num_procs")
    hyb = weak_agg_hyb[weak_agg_hyb["distribution"] == "1D"].sort_values("num_procs")
    if mpi.empty or hyb.empty:
        return

    common = sorted(set(mpi["num_procs"].values) & set(hyb["num_procs"].values))
    if not common:
        return
    mpi = mpi[mpi["num_procs"].isin(common)]
    hyb = hyb[hyb["num_procs"].isin(common)]

    procs   = np.array(common, dtype=float)
    t1_m    = float(mpi["elapsed_time"].iloc[0])
    t1_h    = float(hyb["elapsed_time"].iloc[0])
    ss_mpi  = procs * (t1_m / mpi["elapsed_time"].values)
    ss_hyb  = procs * (t1_h / hyb["elapsed_time"].values)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))

    # --- scaled speedup ---
    ax1.plot(procs, ss_mpi, marker="o", linewidth=2.5, markersize=8,
             color="#1f77b4", label="MPI")
    ax1.plot(procs, ss_hyb, marker="s", linewidth=2.5, markersize=8,
             color="#ff7f0e", label="MPI+OpenMP")
    ax1.plot([1, procs[-1]], [1, procs[-1]], "k--", linewidth=2, label="Ideal")

    ax1.set_xscale("log", base=2)
    ax1.set_yscale("log", base=2)
    ax1.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f"{int(x)}"))
    ax1.set_xlabel("Number of Processes", fontsize=13, fontweight="bold")
    ax1.set_ylabel("Scaled Speedup (x)", fontsize=13, fontweight="bold")
    ax1.set_title("Weak Scaling – Scaled Speedup  MPI vs MPI+OpenMP",
                  fontsize=14, fontweight="bold")
    ax1.legend(fontsize=11)
    ax1.grid(True, alpha=0.3)

    # --- efficiency ---
    ax2.plot(procs, (t1_m / mpi["elapsed_time"].values) * 100,
             marker="o", linewidth=2.5, markersize=8, color="#1f77b4", label="MPI")
    ax2.plot(procs, (t1_h / hyb["elapsed_time"].values) * 100,
             marker="s", linewidth=2.5, markersize=8, color="#ff7f0e", label="MPI+OpenMP")
    ax2.axhline(100, color="k", linestyle="--", linewidth=2, label="Ideal")

    ax2.set_xscale("log", base=2)
    ax2.set_ylim(0, 110)
    ax2.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f"{int(x)}"))
    ax2.set_xlabel("Number of Processes", fontsize=13, fontweight="bold")
    ax2.set_ylabel("Efficiency (%)", fontsize=13, fontweight="bold")
    ax2.set_title("Weak Scaling – Efficiency  MPI vs MPI+OpenMP",
                  fontsize=14, fontweight="bold")
    ax2.legend(fontsize=11)
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(os.path.join(SUBDIRS["comparison"], "weak_mpi_vs_hybrid.png"),
                dpi=300, bbox_inches="tight")
    plt.close()
    print("    saved  weak_mpi_vs_hybrid.png")


# ===========================================================================
# 6.  MAIN
# ===========================================================================
def main():
    print("=" * 70)
    print("  SpMV Benchmark — Full Plot Suite")
    print("=" * 70)

    # ── load raw CSVs ──
    print("\nLoading results...")
    strong_df       = load_csv(os.path.join(RESULTS_DIR, "strong_scaling_all.csv"))
    weak_df         = load_csv(os.path.join(RESULTS_DIR, "weak_scaling_all.csv"))
    strong_hybrid_df = load_csv(os.path.join(RESULTS_DIR, "strong_scaling_hybrid.csv"))
    weak_hybrid_df   = load_csv(os.path.join(RESULTS_DIR, "weak_scaling_hybrid.csv"))

    # ── 90th percentile across ranks ──
    print("\nComputing 90th percentile across ranks …")
    strong_90   = calculate_90th_percentile(strong_df)
    weak_90     = calculate_90th_percentile(weak_df)
    strong_h90  = calculate_90th_percentile(strong_hybrid_df)
    weak_h90    = calculate_90th_percentile(weak_hybrid_df)

    # persist reduced CSVs
    strong_90 .to_csv(os.path.join(SUBDIRS["data"], "strong_scaling_90th.csv"),        index=False)
    weak_90   .to_csv(os.path.join(SUBDIRS["data"], "weak_scaling_90th.csv"),          index=False)
    strong_h90.to_csv(os.path.join(SUBDIRS["data"], "strong_scaling_hybrid_90th.csv"), index=False)
    weak_h90  .to_csv(os.path.join(SUBDIRS["data"], "weak_scaling_hybrid_90th.csv"),   index=False)
    print("  MPI + Hybrid 90th-percentile CSVs saved to plots/data_reduction/")

    # ── collapse runs → median  (keeps distribution!) ──
    strong_agg      = aggregate_over_runs(strong_90,  extra_group=["matrix"])
    weak_agg        = aggregate_over_runs(weak_90)
    strong_hyb_agg  = aggregate_over_runs(strong_h90, extra_group=["matrix"])
    weak_hyb_agg    = aggregate_over_runs(weak_h90)

    matrices = sorted(strong_agg["matrix"].unique())
    print(f"\n  Matrices found: {matrices}")
    print(f"  Distributions : {sorted(strong_agg['distribution'].unique())}")
    print(f"  Proc counts   : {sorted(strong_agg['num_procs'].unique())}\n")

    # ================================================================
    # STRONG SCALING
    # ================================================================
    print("-" * 70)
    print("  STRONG SCALING  –  MPI")
    print("-" * 70)
    plot_strong_speedup(strong_agg, label="MPI")
    plot_strong_comp_comm(strong_agg, label="MPI")
    plot_comm_volume_heatmap(strong_agg, label="MPI")
    plot_strong_comparison(strong_agg, label="MPI")

    print("\n" + "-" * 70)
    print("  STRONG SCALING  –  MPI+OpenMP")
    print("-" * 70)
    plot_strong_speedup(strong_hyb_agg, label="Hybrid")
    plot_strong_comp_comm(strong_hyb_agg, label="Hybrid")
    plot_comm_volume_heatmap(strong_hyb_agg, label="Hybrid")
    plot_strong_comparison(strong_hyb_agg, label="Hybrid")

    print("\n" + "-" * 70)
    print("  STRONG SCALING  –  MPI vs MPI+OpenMP")
    print("-" * 70)
    plot_strong_mpi_vs_hybrid(strong_agg, strong_hyb_agg)
    plot_strong_comp_comm_mpi_vs_hybrid(strong_agg, strong_hyb_agg)

    # ================================================================
    # WEAK SCALING
    # ================================================================
    print("\n" + "-" * 70)
    print("  WEAK SCALING  –  MPI")
    print("-" * 70)
    plot_weak_speedup(weak_agg, label="MPI")
    plot_weak_efficiency(weak_agg, label="MPI")
    plot_weak_time_breakdown(weak_agg, label="MPI")
    plot_weak_comp_comm(weak_agg, label="MPI")
    plot_weak_load_balance(weak_agg, label="MPI")
    export_weak_comm_volume(weak_agg, label="MPI")

    print("\n" + "-" * 70)
    print("  WEAK SCALING  –  MPI+OpenMP")
    print("-" * 70)
    plot_weak_speedup(weak_hyb_agg, label="Hybrid")
    plot_weak_efficiency(weak_hyb_agg, label="Hybrid")
    plot_weak_time_breakdown(weak_hyb_agg, label="Hybrid")
    plot_weak_comp_comm(weak_hyb_agg, label="Hybrid")
    plot_weak_load_balance(weak_hyb_agg, label="Hybrid")
    export_weak_comm_volume(weak_hyb_agg, label="Hybrid")

    print("\n" + "-" * 70)
    print("  WEAK SCALING  –  MPI vs MPI+OpenMP")
    print("-" * 70)
    plot_weak_mpi_vs_hybrid(weak_agg, weak_hyb_agg)

    # ================================================================
    print("\n" + "=" * 70)
    print("  All plots complete.  Output tree:")
    print("=" * 70)
    for root, dirs, files in os.walk(PLOTS_DIR):
        level  = root.replace(PLOTS_DIR, "").count(os.sep)
        indent = " " * 2 * level
        print(f"{indent}{os.path.basename(root)}/")
        sub_indent = " " * 2 * (level + 1)
        for f in sorted(files):
            print(f"{sub_indent}{f}")


if __name__ == "__main__":
    main()