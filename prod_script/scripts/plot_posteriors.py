#!/usr/bin/env python3
"""
plot_posteriors.py

Plots the posterior probability distribution from a production array
assignment run. One PNG per multi-library sample.

For each sample, reads the assign_kinnex.py output and plots a histogram
of Top_Posterior values coloured by classification tier:
  HIGH_CONF  — confident assignments (posterior ≥ POSTERIOR_HIGH_CONF)
  LOW_CONF   — low-confidence assignments (posterior ≥ POSTERIOR_LOW_CONF)
  UNASSIGNED — below threshold or insufficient observations

Vertical lines mark the HIGH_CONF and LOW_CONF posterior thresholds.
Threshold values are read from the assignment file's provenance header
so the plot always reflects the parameters actually used, not hardcoded
defaults.

Unlike visualize_posteriors.py (test_script), this script does not
require ground truth data — it is designed for production runs on real
data where ground truth is not available.

Usage:
    python plot_posteriors.py <assign_file> <out_png>

Arguments:
    assign_file  Path to results/assigned/{sample}.txt
    out_png      Output PNG path (e.g. results/qc/posteriors/{sample}_posteriors.png)

Author:  KM
Created: 2026-03
"""

import sys
import os
import math
import argparse

import pandas as pd
import matplotlib
matplotlib.use("Agg")  # non-interactive backend for cluster use
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import utils

# ---------------------------------------------------------------------------
# Defaults — used only if the assignment file header cannot be parsed
# ---------------------------------------------------------------------------
DEFAULT_POSTERIOR_HIGH_CONF = 0.840
DEFAULT_POSTERIOR_LOW_CONF  = 0.500

# Classification tier colours — chosen for colour-blind accessibility
COLOUR_HIGH_CONF  = "#2E86AB"   # blue
COLOUR_LOW_CONF   = "#F4A261"   # amber
COLOUR_UNASSIGNED = "#BBBBBB"   # grey


def load_thresholds(assign_file):
    """
    Read posterior threshold values from the assignment file's provenance
    header. Falls back to hardcoded defaults if the header is absent or
    the values cannot be parsed.

    Parameters
    ----------
    assign_file : str
        Path to assign_kinnex.py output file.

    Returns
    -------
    tuple (float, float)
        (posterior_high_conf, posterior_low_conf)
    """
    header = utils.parse_assignment_header(assign_file)
    try:
        high = float(header["POSTERIOR_HIGH_CONF"])
        low  = float(header["POSTERIOR_LOW_CONF"])
        return high, low
    except (KeyError, ValueError):
        print(
            f"[plot_posteriors] WARNING: could not read thresholds from "
            f"assignment file header — using defaults "
            f"({DEFAULT_POSTERIOR_HIGH_CONF} / {DEFAULT_POSTERIOR_LOW_CONF})",
            file=sys.stderr,
        )
        return DEFAULT_POSTERIOR_HIGH_CONF, DEFAULT_POSTERIOR_LOW_CONF


def plot_posteriors(assign_file, out_png):
    """
    Generate and save the posterior distribution histogram for one sample.

    Parameters
    ----------
    assign_file : str
        Path to the assignment .txt file from assign_kinnex.py.
    out_png : str
        Output PNG path.
    """
    # ── Load data ──────────────────────────────────────────────────────────
    df = utils.load_assignments_df(assign_file)

    if df.empty:
        print(
            f"[plot_posteriors] WARNING: no assignments found in {assign_file} "
            f"— skipping plot.",
            file=sys.stderr,
        )
        # Write a minimal placeholder so Snakemake output is satisfied
        fig, ax = plt.subplots(figsize=(8, 4))
        ax.text(0.5, 0.5, "No assignments found", ha="center", va="center",
                transform=ax.transAxes, fontsize=14, color="#999999")
        os.makedirs(os.path.dirname(out_png), exist_ok=True) if os.path.dirname(out_png) else None
        fig.savefig(out_png, dpi=150, bbox_inches="tight")
        plt.close(fig)
        return

    df["Top_Posterior"] = pd.to_numeric(df["Top_Posterior"], errors="coerce")
    df = df.dropna(subset=["Top_Posterior"])

    # ── Thresholds ─────────────────────────────────────────────────────────
    thresh_high, thresh_low = load_thresholds(assign_file)

    # ── Derive sample name from filename for the plot title ────────────────
    sample = os.path.splitext(os.path.basename(assign_file))[0]

    # ── Split by classification tier ───────────────────────────────────────
    hc         = df[df["Classification"] == "HIGH_CONF"]["Top_Posterior"]
    lc         = df[df["Classification"] == "LOW_CONF"]["Top_Posterior"]
    unassigned = df[df["Classification"] == "UNASSIGNED"]["Top_Posterior"]

    n_total = len(df)
    n_hc    = len(hc)
    n_lc    = len(lc)
    n_ua    = len(unassigned)

    # ── Bin edges: 50 bins across [0, 1] ──────────────────────────────────
    # Posteriors for UNASSIGNED ZMWs are the top posterior score even though
    # they didn't meet the classification thresholds. Including them in the
    # histogram at their actual posterior value (often near 0.5) is more
    # honest than excluding them.
    bins = [i / 50 for i in range(51)]  # 0.00, 0.02, 0.04, ... 1.00

    # ── Figure ─────────────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(10, 5))
    fig.patch.set_facecolor("#f8f9fb")
    ax.set_facecolor("#f8f9fb")

    # Stacked histogram — plot in order: UNASSIGNED (bottom), LOW_CONF, HIGH_CONF
    # bottom= accumulates so bars stack correctly
    counts_ua, _ = ax.hist(
        unassigned, bins=bins,
        color=COLOUR_UNASSIGNED, edgecolor="white", linewidth=0.4,
        label=f"UNASSIGNED  ({n_ua:,}  {100*n_ua/n_total:.1f}%)",
    )[0:2] if len(unassigned) > 0 else ([], bins)

    # Build stacked bottoms manually for correct stacking
    import numpy as np
    bottom_lc = np.zeros(len(bins) - 1)
    if len(unassigned) > 0:
        bottom_lc, _ = np.histogram(unassigned, bins=bins)

    bottom_hc = bottom_lc.copy()
    if len(lc) > 0:
        lc_counts, _ = np.histogram(lc, bins=bins)
        ax.bar(
            [b + 0.01 for b in bins[:-1]],  # slight offset = bin centre approx
            lc_counts,
            width=0.02,
            bottom=bottom_lc,
            color=COLOUR_LOW_CONF, edgecolor="white", linewidth=0.4,
            label=f"LOW_CONF  ({n_lc:,}  {100*n_lc/n_total:.1f}%)",
            align="edge",
        )
        bottom_hc = bottom_lc + lc_counts

    if len(hc) > 0:
        hc_counts, _ = np.histogram(hc, bins=bins)
        ax.bar(
            [b + 0.01 for b in bins[:-1]],
            hc_counts,
            width=0.02,
            bottom=bottom_hc,
            color=COLOUR_HIGH_CONF, edgecolor="white", linewidth=0.4,
            label=f"HIGH_CONF  ({n_hc:,}  {100*n_hc/n_total:.1f}%)",
            align="edge",
        )

    # Redraw UNASSIGNED on top of the bar() calls so the first hist() call
    # is replaced — easier to use ax.bar() for all three tiers consistently.
    # Clear and redo all three as bar() calls.
    ax.cla()
    ax.set_facecolor("#f8f9fb")

    ua_counts = np.zeros(len(bins) - 1)
    lc_counts = np.zeros(len(bins) - 1)
    hc_counts = np.zeros(len(bins) - 1)

    if len(unassigned) > 0:
        ua_counts, _ = np.histogram(unassigned, bins=bins)
    if len(lc) > 0:
        lc_counts, _ = np.histogram(lc, bins=bins)
    if len(hc) > 0:
        hc_counts, _ = np.histogram(hc, bins=bins)

    bar_width  = bins[1] - bins[0]
    bar_lefts  = np.array(bins[:-1])

    ax.bar(bar_lefts, ua_counts,
           width=bar_width, align="edge",
           color=COLOUR_UNASSIGNED, edgecolor="white", linewidth=0.3,
           label=f"UNASSIGNED  ({n_ua:,},  {100*n_ua/n_total:.1f}%)")
    ax.bar(bar_lefts, lc_counts,
           width=bar_width, align="edge",
           bottom=ua_counts,
           color=COLOUR_LOW_CONF, edgecolor="white", linewidth=0.3,
           label=f"LOW_CONF  ({n_lc:,},  {100*n_lc/n_total:.1f}%)")
    ax.bar(bar_lefts, hc_counts,
           width=bar_width, align="edge",
           bottom=ua_counts + lc_counts,
           color=COLOUR_HIGH_CONF, edgecolor="white", linewidth=0.3,
           label=f"HIGH_CONF  ({n_hc:,},  {100*n_hc/n_total:.1f}%)")

    # ── Threshold lines ────────────────────────────────────────────────────
    ax.axvline(thresh_high, color="#1a1a2e", linewidth=1.4, linestyle="--", zorder=5)
    ax.axvline(thresh_low,  color="#1a1a2e", linewidth=1.0, linestyle=":",  zorder=5)

    y_max = ax.get_ylim()[1]
    ax.text(thresh_high + 0.01, y_max * 0.97,
            f"HIGH_CONF\n{thresh_high:.3f}",
            fontsize=8, va="top", color="#1a1a2e",
            fontfamily="monospace")
    ax.text(thresh_low + 0.01, y_max * 0.97,
            f"LOW_CONF\n{thresh_low:.3f}",
            fontsize=8, va="top", color="#1a1a2e",
            fontfamily="monospace")

    # ── Labels and formatting ──────────────────────────────────────────────
    ax.set_xlabel("Top Posterior Probability", fontsize=11, labelpad=8)
    ax.set_ylabel("ZMW count", fontsize=11, labelpad=8)
    ax.set_title(
        f"{sample}  —  ZMW posterior distribution  (n={n_total:,})",
        fontsize=12, fontweight="bold", pad=12, color="#1a1a2e",
    )
    ax.set_xlim(0, 1)
    ax.legend(
        frameon=True, framealpha=0.9, edgecolor="#cccccc",
        fontsize=9, loc="upper left",
    )
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color("#cccccc")
    ax.spines["bottom"].set_color("#cccccc")
    ax.tick_params(colors="#555555")
    ax.yaxis.set_major_formatter(
        matplotlib.ticker.FuncFormatter(lambda x, _: f"{int(x):,}")
    )

    fig.tight_layout()

    out_dir = os.path.dirname(out_png)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    fig.savefig(out_png, dpi=150, bbox_inches="tight", facecolor=fig.get_facecolor())
    plt.close(fig)
    print(f"[plot_posteriors] Written: {out_png}", file=sys.stderr)


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Plot posterior distribution histogram for one production sample."
    )
    parser.add_argument("assign_file", help="Path to results/assigned/{sample}.txt")
    parser.add_argument("out_png",     help="Output PNG path")
    args = parser.parse_args()

    if not os.path.exists(args.assign_file):
        print(f"ERROR: assignment file not found: {args.assign_file}", file=sys.stderr)
        sys.exit(1)

    plot_posteriors(args.assign_file, args.out_png)


if __name__ == "__main__":
    main()
