#!/usr/bin/env python3
"""Build the PDAC figures from aggregate analysis tables."""

from __future__ import annotations

import csv
import argparse
from pathlib import Path

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.text import Text
import numpy as np
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Rectangle
from matplotlib.lines import Line2D
from matplotlib.transforms import Bbox


HERE = Path(__file__).resolve().parent
DEFAULT_DATA_DIR = HERE / "source_data"
DEFAULT_OUTPUT_DIR = HERE / "output"
OUT = DEFAULT_OUTPUT_DIR

COHORTS = ["Fudan", "CPTAC", "APGI", "Ruijin"]
COHORT_LABELS = ["Fudan", "CPTAC-PDAC", "APGI", "Ruijin"]
SETS = ["R_protocol", "R_strict"]
SET_LABELS = {"R_protocol": "R", "R_strict": "R_strict"}
RULES = ["single_A", "single_B", "equal_weight", "continuous"]
RULE_LABELS = {
    "single_A": "Single A",
    "single_B": "Single B",
    "equal_weight": "Equal-weight",
    "continuous": "Continuous",
    "double_sig": "Double-significant-first",
}
RULE_COLORS = {
    "single_A": "#4477AA",
    "single_B": "#CC6677",
    "equal_weight": "#228833",
    "continuous": "#AA3377",
    "double_sig": "#66CCEE",
}
SET_COLORS = {"R_protocol": "#4477AA", "R_strict": "#CC6677"}
COHORT_COLORS = dict(zip(COHORTS, ["#4477AA", "#CC6677", "#228833", "#AA3377"]))
PROTEIN_INFO = {
    "Fudan": "Mass spectrometry; linear intensity",
    "CPTAC": "TMT; log2 MD abundance",
    "APGI": "DIA/MaxLFQ",
    "Ruijin": "DIA; released protein abundance",
}

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
    "font.size": 10,
    "font.weight": "bold",
    "axes.titlesize": 12,
    "axes.titleweight": "bold",
    "axes.labelsize": 10.5,
    "axes.labelweight": "bold",
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "legend.fontsize": 10,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.linewidth": 0.7,
    "pdf.fonttype": 42,
    "svg.fonttype": "none",
    "savefig.dpi": 600,
})


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def save_figure(fig: plt.Figure, name: str) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    if name.startswith("FigS"):
        for artist in fig.findobj(match=Text):
            artist.set_fontweight("bold")
        for ext in ("svg", "pdf", "png"):
            fig.savefig(OUT / f"{name}.{ext}", dpi=300 if ext == "png" else None, bbox_inches="tight")
        plt.close(fig)
        return
    fig.canvas.draw()
    for artist in fig.findobj(match=Text):
        if artist.get_text().strip():
            artist.set_fontweight("bold")
            if artist.get_fontsize() < 9.5:
                artist.set_fontsize(10)
    crop = Bbox.from_extents(0.015, 0.55, 7.285, 4.11) if name == "Fig1" else "tight"
    fig.savefig(OUT / f"{name}.svg", bbox_inches=crop, pad_inches=0.06)
    fig.savefig(OUT / f"{name}.pdf", bbox_inches=crop, pad_inches=0.06)
    fig.savefig(OUT / f"{name}.png", dpi=300, bbox_inches=crop, pad_inches=0.06)
    plt.close(fig)


FIGURE_DATA_FILES = {
    "Fig1_cohort_design": "Fig1_cohort_design.tsv",
    "Fig2_support_patient_bootstrap": "Fig2_support_patient_bootstrap.tsv",
    "Fig3_fixed_budget_four_rules": "Fig3_fixed_budget_four_rules.tsv",
    "Fig4_decomposition_k500": "Fig4_decomposition_k500.tsv",
    "FigS1_sensitivity_grid": "FigS1_sensitivity_grid.tsv",
    "FigS2_decomposition_patient_bootstrap": "FigS2_decomposition_patient_bootstrap.tsv",
    "FigS3_historical_purity_summary": "FigS3_historical_purity_summary.tsv",
    "FigS4_shared_testable_contrast": "FigS4_shared_testable_contrast.tsv",
    "FigS5_BioTExt_signed_effect": "FigS5_BioTExt_signed_effect.tsv",
    "FigS6_HPA_fixed_k500": "FigS6_HPA_fixed_k500.tsv",
}


def load_sources(data_dir: Path) -> dict[str, list[dict[str, str]]]:
    """Load the bundled, aggregate figure tables without upstream project files."""
    missing = [filename for filename in FIGURE_DATA_FILES.values() if not (data_dir / filename).is_file()]
    if missing:
        raise FileNotFoundError(
            f"Missing figure source table(s) in {data_dir}: {', '.join(missing)}"
        )
    sources = {
        name: read_tsv(data_dir / filename)
        for name, filename in FIGURE_DATA_FILES.items()
    }
    empty = [FIGURE_DATA_FILES[name] for name, rows in sources.items() if not rows]
    if empty:
        raise ValueError(f"Figure source table(s) are empty: {', '.join(empty)}")
    return sources

def box(ax: plt.Axes, x: float, y: float, w: float, h: float, title: str, lines: list[str], fill: str, edge: str) -> None:
    patch = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.005,rounding_size=0.010", facecolor=fill, edgecolor=edge, linewidth=1.15)
    ax.add_patch(patch)
    ax.text(x + 0.012, y + h - 0.035, title, weight="bold", fontsize=12, va="top")
    ax.text(x + 0.012, y + h - 0.093, "\n".join(lines), fontsize=10, weight="bold", va="top", linespacing=1.13)


def make_fig1(data: list[dict[str, object]]) -> None:
    fig = plt.figure(figsize=(7.3, 4.15))
    ax = fig.add_axes([0, 0, 1, 1])
    meta = data[0]
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
    ax.text(0.025, 0.955, "A", fontsize=18, weight="bold", va="top")
    ax.text(0.105, 0.952, "RNA signature definition", fontsize=13, weight="bold", va="top")
    ax.text(0.525, 0.955, "B", fontsize=18, weight="bold", va="top")
    ax.text(0.605, 0.952, "Paired protein evaluation", fontsize=13, weight="bold", va="top")
    box(ax, 0.035, 0.59, 0.195, 0.21, "GSE28735", [f"{int(meta['gse28735_pairs'])} tumour–normal", "pairs"], "#E8F1F8", "#4477AA")
    box(ax, 0.035, 0.31, 0.195, 0.21, "GSE15471", [f"{int(meta['gse15471_pairs'])} tumour–normal", "pairs"], "#E8F1F8", "#4477AA")
    box(ax, 0.282, 0.38, 0.222, 0.35, "Fixed RNA lists", [f"R = {int(meta['R_total']):,} genes", f"R_strict = {int(meta['R_strict_total']):,} genes", "Defined from", "RNA alone"], "#E7F2E5", "#228833")
    for y, tip in ((0.695, 0.593), (0.415, 0.517)):
        ax.add_patch(FancyArrowPatch((0.238, y), (0.278, tip), arrowstyle="-|>", mutation_scale=12, linewidth=1.3, color="#666666"))
    ax.add_patch(Rectangle((0.564, 0.16), 0.425, 0.66, fill=False, edgecolor="#A8A8A8", linewidth=0.75))
    ax.add_patch(FancyArrowPatch((0.511, 0.555), (0.564, 0.555), arrowstyle="-|>", mutation_scale=12, linewidth=1.3, color="#666666"))
    platform_labels = {"Fudan": "Label-free MS", "CPTAC": "TMT", "APGI": "DIA / MaxLFQ", "Ruijin": "DIA"}
    for index, (cohort, label) in enumerate(zip(COHORTS, COHORT_LABELS)):
        row = next(item for item in data if item["cohort"] == cohort)
        x = 0.579 + 0.208 * (index % 2)
        y = 0.59 - 0.28 * (index // 2)
        title = f"{label}*" if cohort == "APGI" else label
        box(ax, x, y, 0.192, 0.21, title, [f"{int(row['pairs'])} paired cases", platform_labels[cohort]], "#FFF7E6", "#CC8844")
    ax.text(0.579, 0.26, "* Pairs labelled PDAC in the source study.", fontsize=10, color="#444444", va="top", weight="bold")
    save_figure(fig, "Fig1")


def make_fig2(data: list[dict[str, object]]) -> None:
    fig, ax = plt.subplots(figsize=(7.3, 4.9))
    x = np.arange(4); width = 0.27
    for j, sig in enumerate(SETS):
        rows = [next(r for r in data if r["cohort"] == c and r["set"] == sig) for c in COHORTS]
        xx = x + (j - 0.5) * 0.34
        y = np.array([float(r["rate"]) for r in rows]) * 100
        lo = np.array([float(r["lo95"]) for r in rows]) * 100
        hi = np.array([float(r["hi95"]) for r in rows]) * 100
        ax.bar(xx, y, width, color=SET_COLORS[sig], label=SET_LABELS[sig], edgecolor="white", linewidth=0.5)
        ax.errorbar(xx, y, yerr=[y - lo, hi - y], fmt="none", ecolor="#222222", elinewidth=0.8, capsize=2.5, zorder=4)
        for xi, yi, upper, row in zip(xx, y, hi, rows):
            ax.text(
                xi,
                upper + 2.3,
                f"{yi:.1f}%\nn={int(row['n_testable']):,}",
                ha="center",
                va="bottom",
                fontsize=9.5,
                weight="bold",
                linespacing=1.0,
                clip_on=False,
                bbox={"facecolor": "white", "edgecolor": "none", "pad": 0.6},
            )
    ax.set_xticks(x, [f"{label}\n({int(next(r for r in data if r['cohort']==c)['pairs'])} pairs)" for label, c in zip(COHORT_LABELS, COHORTS)])
    ax.set_ylabel("Protein support among testable genes (%)", labelpad=8)
    ax.set_ylim(0, 104); ax.set_yticks(np.arange(0, 101, 20))
    ax.grid(axis="y", color="#DDDDDD", linewidth=0.45, zorder=0)
    ax.set_axisbelow(True)
    handles, legend_labels = ax.get_legend_handles_labels()
    ax.legend(handles, legend_labels, frameon=False, loc="lower center", bbox_to_anchor=(0.5, 1.025), ncol=2, prop={"weight": "bold", "size": 10.5})
    fig.text(0.04, 0.032, "n = testable genes. Whiskers show 95% patient-pair bootstrap intervals.", ha="left", fontsize=10, weight="bold", color="#444444")
    fig.tight_layout(rect=[0.015, 0.105, 0.985, 0.93], pad=0.65)
    save_figure(fig, "Fig2")


def make_fig3(data: list[dict[str, object]]) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(7.3, 6.1), sharex=True, sharey=True)
    ks = [100, 250, 500, 1000, 2000]
    markers = {"single_A": "o", "single_B": "s", "equal_weight": "^", "continuous": "D"}
    for panel_index, (ax, cohort, label) in enumerate(zip(axes.flat, COHORTS, COHORT_LABELS)):
        subset = [r for r in data if r["cohort"] == cohort]
        for rule in RULES:
            vals = sorted((r for r in subset if r["rule"] == rule), key=lambda r: int(r["k"]))
            ax.plot([int(r["k"]) for r in vals], [float(r["supported_per100"]) for r in vals], marker=markers[rule], markersize=5.4, linewidth=1.8, color=RULE_COLORS[rule], markerfacecolor="none" if rule == "continuous" else RULE_COLORS[rule], label=RULE_LABELS[rule])
        ax.axvline(500, color="#555555", linestyle="--", linewidth=1.0)
        panel_letter = chr(ord("A") + panel_index)
        ax.annotate(panel_letter, (0, 1), xycoords="axes fraction", xytext=(-22, 17), textcoords="offset points", fontsize=18, fontweight="bold", ha="left", va="center", annotation_clip=False)
        ax.set_title(label, loc="left", weight="bold", fontsize=12, pad=12)
        # Candidate budgets are strictly positive, so the log axis is defined.
        if any(int(r["k"]) <= 0 for r in subset):
            raise ValueError(f"Log-scaled candidate budget must be positive for {cohort}")
        ax.set_xscale("log"); ax.set_xticks(ks, ["100", "250", "500", "1,000", "2,000"])
        ax.set_ylim(0, 100); ax.set_yticks(np.arange(0, 101, 20))
        ax.grid(axis="y", color="#DDDDDD", linewidth=0.45)
    axes[0, 0].set_ylabel("Supported per 100 nominations", labelpad=7)
    axes[1, 0].set_ylabel("Supported per 100 nominations", labelpad=7)
    axes[1, 0].set_xlabel("Candidates nominated", labelpad=6)
    axes[1, 1].set_xlabel("Candidates nominated", labelpad=6)
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, frameon=False, loc="upper center", bbox_to_anchor=(0.5, 0.99), ncol=4, prop={"weight": "bold", "size": 9.5})
    fig.text(0.04, 0.025, "Dashed line marks 500 nominations. All nominations count toward the budget.", ha="left", fontsize=10, weight="bold", color="#444444")
    fig.subplots_adjust(left=0.14, right=0.977, top=0.84, bottom=0.145, wspace=0.22, hspace=0.38)
    save_figure(fig, "Fig3")


def make_fig4(data: list[dict[str, object]]) -> None:
    fig, ax = plt.subplots(figsize=(7.3, 4.55))
    centers = np.arange(4) * 1.85
    width = 0.32
    coverage = [float(next(r for r in data if r["cohort"] == c)["coverage_component_per100"]) for c in COHORTS]
    residual = [float(next(r for r in data if r["cohort"] == c)["residual_per100"]) for c in COHORTS]
    totals = [float(next(r for r in data if r["cohort"] == c)["total_gain_per100"]) for c in COHORTS]
    label_background = {"facecolor": "white", "edgecolor": "none", "pad": 0.5}
    for center, cov, res, total in zip(centers, coverage, residual, totals):
        x_cov, x_res, x_total = center - 0.55, center, center + 0.55
        ax.bar(x_cov, cov, width, bottom=0, color="#66C2A5", edgecolor="white", linewidth=0.5, zorder=3)
        ax.bar(x_res, res, width, bottom=cov, color="#E6AB5B", edgecolor="white", linewidth=0.5, zorder=3)
        ax.bar(x_total, total, width, bottom=0, color="#4477AA", edgecolor="white", linewidth=0.5, zorder=3)
        ax.plot([x_cov + width / 2, x_res - width / 2], [cov, cov], color="#888888", linewidth=0.65, zorder=2)
        ax.plot([x_res + width / 2, x_total - width / 2], [total, total], color="#888888", linewidth=0.65, zorder=2)
        ax.text(x_cov, cov + 0.18, f"{cov:+.2f}", ha="center", va="bottom", fontsize=10, weight="bold", bbox=label_background)
        ax.text(x_res, cov + res + 0.18, f"{res:+.2f}", ha="center", va="bottom", fontsize=10, weight="bold", bbox=label_background)
        total_offset = 0.15 if total >= 0 else -0.15
        ax.text(x_total, total + total_offset, f"{total:+.2f}", ha="center", va="bottom" if total >= 0 else "top", fontsize=10, weight="bold", bbox=label_background)
    ax.set_xticks(centers, COHORT_LABELS)
    ax.set_ylabel("Continuous minus Single A\n(change per 100 nominations)")
    bounds = [0.0]
    for cov, total in zip(coverage, totals):
        bounds.extend([cov, total])
    lower = min(bounds) - 1.1
    upper = max(bounds) + 1.1
    ax.set_ylim(lower, upper)
    ax.yaxis.set_major_locator(plt.MaxNLocator(6))
    ax.axhline(0, color="#333333", linewidth=0.65)
    ax.grid(axis="y", color="#DDDDDD", linewidth=0.45); ax.set_axisbelow(True)
    ax.set_xlim(centers[0] - 0.8, centers[-1] + 0.8)
    fig.text(0.13, 0.97, "Coverage and residual contributions at 500 nominations", va="top", weight="bold", fontsize=12)
    fig.legend(handles=[Rectangle((0, 0), 1, 1, color=c) for c in ["#66C2A5", "#E6AB5B", "#4477AA"]], labels=["Coverage", "Residual", "Total"], frameon=False, loc="upper left", bbox_to_anchor=(0.122, 0.888), ncol=3, prop={"weight": "bold", "size": 10.5})
    fig.text(0.04, 0.035, "Coverage and residual contributions sum to the observed difference in support.", ha="left", fontsize=10, weight="bold", color="#444444")
    fig.subplots_adjust(left=0.13, right=0.973, top=0.756, bottom=0.19)
    save_figure(fig, "Fig4")


SUPPLEMENT_STYLE = {
    "axes.labelsize": 10,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "legend.fontsize": 9,
    "legend.title_fontsize": 9,
}


def add_panel_tag(ax: plt.Axes, label: str) -> None:
    ax.text(
        -0.13, 1.035, label, transform=ax.transAxes,
        ha="left", va="bottom", fontsize=18, fontweight="bold", clip_on=False,
    )


def make_fig_s1(data: list[dict[str, object]]) -> None:
    rows: list[tuple[str, str, str, str]] = []
    for scope in ["frozen_set", "top500"]:
        for cohort in COHORTS:
            identifiers = sorted({str(r["identifier"]) for r in data if r["scope"] == scope and r["cohort"] == cohort})
            if scope == "frozen_set":
                identifiers = [sig for sig in SETS if sig in identifiers]
            else:
                # Spotlight the two rules cited in the Results while retaining
                # the complete threshold axis; the complete grid remains in Table S1.
                rule_order = ["single_A", "continuous"]
                identifiers = [rule for rule in rule_order if rule in identifiers]
            for identifier in identifiers:
                scope_label = "Fixed list" if scope == "frozen_set" else "Top 500"
                row_label = f"{COHORT_LABELS[COHORTS.index(cohort)]} · {SET_LABELS.get(identifier, RULE_LABELS.get(identifier, identifier))}\n{scope_label}"
                rows.append((scope, cohort, identifier, row_label))
    columns = [(5, "none"), (10, "none"), (15, "none"), (20, "none"), (10, "none"), (10, "0.25"), (10, "0.5")]
    cell_lookup = {(str(r["scope"]), str(r["cohort"]), str(r["identifier"]), int(r["min_pairs"]), str(r["effect_floor"])): float(r["supported_per100_itt"]) for r in data}
    matrix = np.full((len(rows), len(columns)), np.nan, dtype=float)
    for i, (scope, cohort, identifier, _) in enumerate(rows):
        for j, (min_pairs, effect_floor) in enumerate(columns):
            key = (scope, cohort, identifier, min_pairs, effect_floor)
            if key in cell_lookup:
                matrix[i, j] = cell_lookup[key]
    if np.isnan(matrix).any():
        raise ValueError("Incomplete FigS1 sensitivity grid for a required min-pair/effect-floor cell")

    # The sensitivity grid is dense; give each row more vertical space and
    # keep the exported page wide enough for its numeric labels to read cleanly.
    height = max(7.2, 0.43 * len(rows) + 1.9)
    fig, ax = plt.subplots(figsize=(7.5, height))
    cmap = mpl.colors.LinearSegmentedColormap.from_list("mcp_support", ["#F7F7F7", "#DCEAF4", "#78A8C1", "#355F7E"])
    image = ax.imshow(matrix, aspect="auto", interpolation="nearest", cmap=cmap, vmin=0, vmax=100)
    ax.set_yticks(np.arange(len(rows)), [row[3] for row in rows], fontsize=9.2)
    ax.set_xticks(np.arange(7), ["5", "10*", "15", "20", "none", "0.25", "0.50"], fontsize=9.5)
    ax.set_xlim(-0.5, 6.5)
    ax.set_ylim(len(rows) - 0.5, -1.25)
    ax.text(1.5, -0.85, "Minimum pairs (no effect floor)", ha="center", va="center", fontsize=9.5, weight="bold")
    ax.text(5.0, -0.85, "Effect floor at 10 pairs", ha="center", va="center", fontsize=9.5, weight="bold")
    ax.set_title("Protein support across pair and effect thresholds", loc="center", weight="bold", pad=26, fontsize=12)
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            value = matrix[i, j]
            background = cmap(value / 100)
            text_color = "#17232B" if sum(channel * weight for channel, weight in zip(background[:3], (0.2126, 0.7152, 0.0722))) > 0.51 else "#FFFFFF"
            ax.text(j, i, f"{value:.1f}", ha="center", va="center", fontsize=9.5, fontweight="bold", color=text_color)
    if any(scope == "frozen_set" for scope, _, _, _ in rows) and any(scope == "top500" for scope, _, _, _ in rows):
        split = sum(1 for scope, _, _, _ in rows if scope == "frozen_set")
        ax.axhline(split - 0.5, color="#555555", linewidth=0.9)
        ax.add_patch(Rectangle((0.5, -0.5), 1, len(rows), fill=False, edgecolor="#222222", linewidth=0.8))
        ax.add_patch(Rectangle((3.5, -0.5), 1, len(rows), fill=False, edgecolor="#222222", linewidth=0.8))
    ax.axvline(3.5, color="#555555", linewidth=0.8)
    colorbar = fig.colorbar(image, ax=ax, fraction=0.025, pad=0.02)
    colorbar.set_label("Supported per 100 nominations", fontsize=9.5, fontweight="bold")
    colorbar.ax.tick_params(labelsize=9)
    fig.text(0.015, 0.070, "* Both outlined columns use the primary setting: at least 10 complete pairs and no effect floor.", ha="left", fontsize=8.6, color="#333333")
    fig.text(0.015, 0.040, "Effect floors use each cohort's analysed scale. See Table S1 for all rules.", ha="left", fontsize=8.6, color="#333333")
    fig.subplots_adjust(left=0.31, right=0.96, top=0.88, bottom=0.15)
    save_figure(fig, "FigS1")


def make_fig_s2(data: list[dict[str, object]], fig1_data: list[dict[str, object]]) -> None:
    baselines = ["single_A", "single_B", "equal_weight"]
    labels = ["Single A", "Single B", "Equal-weight"]
    components = [("coverage", "Coverage", "#287E6B"), ("residual", "Residual", "#A96916")]
    fig, axes = plt.subplots(2, 2, figsize=(7.5, 5.9), sharey=True)
    for panel, (ax, cohort, cohort_label) in enumerate(zip(axes.flat, COHORTS, COHORT_LABELS)):
        cohort_row = next(r for r in fig1_data if r["cohort"] == cohort)
        n_pairs = int(cohort_row["pairs"])
        for i, baseline in enumerate(baselines):
            for comp_idx, (component, _, color) in enumerate(components):
                row = next(r for r in data if r["cohort"] == cohort and r["baseline"] == baseline and r["component"] == component)
                x = i + (-0.14 if comp_idx == 0 else 0.14)
                mean = float(row["bootstrap_mean"])
                point = float(row["point"])
                lo = float(row["lo95"])
                hi = float(row["hi95"])
                ax.errorbar(x, mean, yerr=[[mean - lo], [hi - mean]], fmt="D", markerfacecolor="white", markeredgecolor=color, ecolor=color, capsize=2.5, markersize=3.7, linewidth=1.0, zorder=3)
                ax.scatter([x - 0.10], [point], s=16, marker="o", color=color, edgecolors="white", linewidths=0.35, zorder=4)
        ax.axhline(0, color="#555555", linewidth=0.65)
        ax.set_xticks(np.arange(3), labels)
        ax.set_title(f"{cohort_label} (n={n_pairs} pairs)", loc="center", weight="bold", fontsize=11)
        add_panel_tag(ax, "ABCD"[panel])
        ax.grid(axis="y", color="#DDDDDD", linewidth=0.45)
    axes[0, 0].set_ylabel("")
    axes[1, 0].set_ylabel("")
    fig.text(0.012, 0.49, "Change in support per 100 nominations", rotation=90, va="center", ha="left", fontsize=9.5, fontweight="bold")
    comp_handles = [Rectangle((0, 0), 1, 1, color=color, label=label) for _, label, color in components]
    point_handle = Line2D([], [], marker="o", linestyle="none", color="#444444", label="Full-sample point")
    boot_handle = Line2D([], [], marker="D", markerfacecolor="white", markeredgecolor="#444444", linestyle="none", label="Bootstrap mean (95% CI)")
    fig.legend(handles=comp_handles + [point_handle, boot_handle], frameon=False, loc="upper center", bbox_to_anchor=(0.5, 0.955), ncol=4, fontsize=8.6)
    fig.suptitle("Coverage and residual contributions at 500 nominations", y=0.995, fontsize=12, weight="bold")
    fig.text(0.02, 0.045, "Each contrast compares Continuous with the named RNA rule.", ha="left", fontsize=8.6, color="#333333")
    fig.text(0.02, 0.018, "Whiskers show 95% confidence intervals from resampling patient pairs.", ha="left", fontsize=8.6, color="#333333")
    fig.subplots_adjust(left=0.12, right=0.98, top=0.83, bottom=0.13, hspace=0.55, wspace=0.34)
    save_figure(fig, "FigS2")


def make_fig_s3(data: list[dict[str, object]]) -> None:
    order = ["R", "R_strict", "background"]
    rows = {r["group"]: r for r in data}
    x = np.arange(3)
    values = [float(rows[k]["median_abs_rho_purity"]) for k in order]
    ns = [int(rows[k]["n"]) for k in order]
    fig, ax = plt.subplots(figsize=(6.2, 4.6))
    colors = ["#4477AA", "#CC6677", "#BBBBBB"]
    ax.bar(x, values, 0.56, color=colors, edgecolor="white", linewidth=0.5)
    for xi, val, n in zip(x, values, ns):
        ax.text(xi, val + 0.004, f"{val:.3f}\nn={n:,}", ha="center", va="bottom", fontsize=10, fontweight="bold", bbox={"facecolor": "white", "edgecolor": "none", "pad": 0.6})
    ax.set_xticks(x, ["R (historical list)", "R_strict (historical list)", "Background"])
    ax.set_ylabel("Median |Spearman ρ| with tumour purity", fontsize=10)
    ax.set_ylim(0, 0.135); ax.set_yticks(np.arange(0, 0.131, 0.03))
    ax.grid(axis="y", color="#DDDDDD", linewidth=0.45); ax.set_axisbelow(True)
    ax.set_title("Within-tumour purity correlations in the original CPTAC release", loc="center", weight="bold", fontsize=11)
    fig.text(0.015, 0.046, "n = proteins in each group. Confidence intervals are unavailable for this summary.", ha="left", fontsize=8.6, color="#333333")
    fig.tight_layout(rect=[0, 0.12, 1, 0.93])
    save_figure(fig, "FigS3")


def make_fig_s4(data: list[dict[str, object]]) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(7.5, 6.6), sharex=True, sharey=True)
    rules = ["single_A", "single_B", "equal_weight", "continuous", "double_sig"]
    labels = ["Single A", "Single B", "Equal", "Continuous", "Double-sig"]
    x = np.arange(5); width = 0.34
    for panel, (ax, cohort, cohort_label) in enumerate(zip(axes.flat, COHORTS, COHORT_LABELS)):
        rows = {r["rule"]: r for r in data if r["cohort"] == cohort}
        itt = [float(rows[rule]["primary_ITT_per100"]) for rule in rules]
        shared = [float(rows[rule]["shared_per100"]) for rule in rules]
        ax.bar(x - width / 2, itt, width, color="#4477AA", label="Primary ITT")
        for xi, rule in zip(x, rules):
            row = rows[rule]
            mean = float(row["shared_bootstrap_mean"])
            lo = float(row["shared_lo95"])
            hi = float(row["shared_hi95"])
            point = float(row["shared_per100"])
            ax.errorbar(xi + 0.30, mean, yerr=[[mean - lo], [hi - mean]], fmt="D", markerfacecolor="white", markeredgecolor="#8F5A16", ecolor="#8F5A16", capsize=1.8, markersize=3.1, linewidth=0.8, zorder=4)
            ax.scatter([xi + 0.08], [point], s=16, marker="o", color="#E6AB5B", edgecolors="#8F5A16", linewidths=0.45, zorder=5)
        ax.set_title(cohort_label, loc="center", weight="bold", fontsize=11)
        add_panel_tag(ax, "ABCD"[panel])
        ax.set_ylim(0, 100); ax.set_yticks(np.arange(0, 101, 20))
        ax.set_xlim(-0.55, 4.65)
        ax.grid(axis="y", color="#DDDDDD", linewidth=0.45)
    fig.text(0.012, 0.52, "Supported per 100 nominations", rotation=90, va="center", ha="left", fontsize=9.5)
    axes[1, 0].set_xticks(x, labels, rotation=25, ha="right")
    axes[1, 1].set_xticks(x, labels, rotation=25, ha="right")
    legend_handles = [
        Rectangle((0, 0), 1, 1, facecolor="#4477AA", edgecolor="none"),
        Line2D([], [], marker="o", linestyle="none", color="#E6AB5B", markeredgecolor="#8F5A16", markersize=4),
        Line2D([], [], marker="D", linestyle="none", markerfacecolor="white", markeredgecolor="#8F5A16", color="#8F5A16", markersize=4),
    ]
    legend_labels = ["Full RNA universe: observed", "Shared-testable genes: observed", "Shared-testable genes: bootstrap mean (95% CI)"]
    fig.legend(legend_handles, legend_labels, frameon=False, loc="upper center", bbox_to_anchor=(0.5, 0.94), ncol=3, fontsize=8.2)
    fig.suptitle("Protein support across two eligible gene sets", y=0.99, fontsize=12, weight="bold")
    universe_sizes = sorted({int(r["shared_universe_n"]) for r in data})
    universe_text = ", ".join(f"{n:,}" for n in universe_sizes)
    fig.text(0.015, 0.065, f"Each rule nominates 500 candidates. The sensitivity analysis starts from {universe_text} genes testable in all four cohorts.", ha="left", fontsize=8.3, color="#333333")
    fig.text(0.015, 0.035, "Whiskers show 95% confidence intervals from resampling patient pairs.", ha="left", fontsize=8.3, color="#333333")
    fig.subplots_adjust(left=0.12, right=0.98, top=0.81, bottom=0.16, hspace=0.48, wspace=0.34)
    save_figure(fig, "FigS4")


def make_fig_s5(data: list[dict[str, object]]) -> None:
    fig = plt.figure(figsize=(7.5, 6.8))
    layout = fig.add_gridspec(2, 2, height_ratios=[2.1, 1.0], hspace=0.54, wspace=0.36)
    first_ax = fig.add_subplot(layout[0, 0])
    second_ax = fig.add_subplot(layout[0, 1], sharey=first_ax)
    axes = [first_ax, second_ax]
    contrast_ax = fig.add_subplot(layout[1, :])
    tissue_order = ["Bulk", "Epithelial", "Stroma"]
    tissue_labels_by_key = {"Bulk": "Bulk", "Epithelial": "Epithelial-enriched", "Stroma": "Stromal-enriched"}
    pair_counts = {}
    for tissue in tissue_order:
        counts = {int(r["n_pairs"]) for r in data if r["tissue"] == tissue and r["set"] in SETS}
        if len(counts) != 1:
            raise ValueError(f"BioTExt matched-pair count differs across signatures for {tissue}: {counts}")
        pair_counts[tissue] = counts.pop()
    tissue_labels = [f"{tissue_labels_by_key[t]}\n(n={pair_counts[t]} pairs)" for t in tissue_order]
    tissue_colors = {"Bulk": "#4477AA", "Epithelial": "#CC6677", "Stroma": "#228833"}
    for panel, (ax, sig, title) in enumerate(zip(axes, SETS, ["R", "R_strict"])):
        rows = {str(r["tissue"]): r for r in data if r["set"] == sig}
        if not all(tissue in rows for tissue in tissue_order + ["Epithelial_minus_Bulk"]):
            raise ValueError(f"Incomplete BioTExt tissue summary for {sig}")
        for i, tissue in enumerate(tissue_order):
            row = rows[tissue]
            mean = float(row["mean_signed_log2FC"])
            lo = float(row["lo95"])
            hi = float(row["hi95"])
            ax.errorbar(i, mean, yerr=[[mean - lo], [hi - mean]], fmt="o", color=tissue_colors[tissue], ecolor=tissue_colors[tissue], capsize=3, markersize=4.5, linewidth=1.1, zorder=3)
            if tissue == "Stroma":
                ax.annotate(
                    f"{int(row['gene_BH05_hits'])} genes at\nBH q<0.05",
                    (i, hi), xytext=(2, 6), textcoords="offset points",
                    ha="right", va="bottom", fontsize=9.2, bbox={"facecolor": "white", "edgecolor": "none", "pad": 0.6},
                )
        diff = next(r for r in data if r["set"] == sig and r["tissue"] == "Epithelial_minus_Bulk")
        ax.axhline(0, color="#555555", linewidth=0.65)
        ax.set_xticks(np.arange(3), tissue_labels)
        ax.set_title(f"{title} (shared genes={int(diff['n_genes']):,})", loc="center", weight="bold", fontsize=11)
        add_panel_tag(ax, "AB"[panel])
        ax.grid(axis="y", color="#DDDDDD", linewidth=0.45)
    axes[0].set_ylabel("Mean RNA-aligned protein log2 fold change")
    axes[0].set_ylim(-0.25, 1.55)
    fig.suptitle("RNA-aligned protein changes by BioTExt tissue preparation", y=0.985, fontsize=12, weight="bold")
    contrast_rows = [next(r for r in data if r["set"] == sig and r["tissue"] == "Epithelial_minus_Bulk") for sig in SETS]
    contrast_y = np.arange(len(SETS))
    for yi, sig, row in zip(contrast_y, SETS, contrast_rows):
        mean = float(row["mean_signed_log2FC"])
        lo = float(row["lo95"])
        hi = float(row["hi95"])
        contrast_ax.errorbar(mean, yi, xerr=[[mean - lo], [hi - mean]], fmt="o", color=SET_COLORS[sig], ecolor=SET_COLORS[sig], capsize=2.5, markersize=4.2, linewidth=1.0, zorder=3)
    contrast_ax.axvline(0, color="#555555", linewidth=0.7)
    contrast_ax.set_yticks(contrast_y, ["R", "R_strict"])
    contrast_ax.set_ylim(1.5, -0.5)
    contrast_ax.set_xlim(-0.25, 0.45)
    contrast_ax.set_xlabel("Mean paired epithelial − bulk difference in RNA-aligned protein log2 fold change\n(95% bootstrap CI)", fontsize=9.5)
    contrast_ax.set_title("Paired epithelial–bulk contrast (n=15 pairs per signature)", loc="center", weight="bold", fontsize=10.5)
    contrast_ax.grid(axis="x", color="#DDDDDD", linewidth=0.45)
    contrast_ax.set_axisbelow(True)
    bh_family = int(next(r for r in data if r["set"] == "R_protocol" and r["tissue"] == "Stroma")["BH_family_size"])
    r_stroma = next(r for r in data if r["set"] == "R_protocol" and r["tissue"] == "Stroma")
    strict_stroma = next(r for r in data if r["set"] == "R_strict" and r["tissue"] == "Stroma")
    r_genes = int(r_stroma["n_genes"])
    strict_genes = int(strict_stroma["n_genes"])
    r_hits, strict_hits = int(r_stroma["gene_BH05_hits"]), int(strict_stroma["gene_BH05_hits"])
    fig.text(0.015, 0.070, "Stromal estimates use tumour–adjacent pairs. The within-patient contrast in C compares epithelial-enriched with bulk.", ha="left", fontsize=8.2, color="#333333")
    fig.text(0.015, 0.035, f"Stromal gene tests at BH q<0.05: {r_hits} of {r_genes:,} (R) and {strict_hits} of {strict_genes:,} (R_strict). Adjustment covered {bh_family:,} tests.", ha="left", fontsize=8.2, color="#333333")
    fig.subplots_adjust(left=0.13, right=0.98, top=0.87, bottom=0.20, hspace=0.54, wspace=0.36)
    a_bbox = axes[0].get_position()
    c_bbox = contrast_ax.get_position()
    c_x = a_bbox.x0 - 0.13 * a_bbox.width
    c_y = c_bbox.y1 + 0.02 * c_bbox.height
    fig.text(c_x, c_y, "C", ha="left", va="bottom", fontsize=18,
             fontweight="bold")
    save_figure(fig, "FigS5")


def make_fig_s6(data: list[dict[str, object]]) -> None:
    rules = ["single_A", "single_B", "equal_weight", "continuous", "double_sig"]
    rows = {r["rule"]: r for r in data}
    values = [(float(rows[k]["supported_per100_low"]) + float(rows[k]["supported_per100_high"])) / 2 for k in rules]
    lo = [float(rows[k]["supported_per100_low"]) for k in rules]
    hi = [float(rows[k]["supported_per100_high"]) for k in rules]
    fig, ax = plt.subplots(figsize=(7.5, 4.7))
    x = np.arange(len(rules))
    ax.bar(x, values, 0.62, color=[RULE_COLORS[k] for k in rules], edgecolor="white", linewidth=0.5)
    ax.errorbar(x, values, yerr=[np.asarray(values) - np.asarray(lo), np.asarray(hi) - np.asarray(values)], fmt="none", ecolor="#222222", capsize=2.0, linewidth=0.8, zorder=4)
    for xi, k, val in zip(x, rules, values):
        row = rows[k]
        low, high = int(row["supported_low"]), int(row["supported_high"])
        count = f"{low}/500" if low == high else f"{low}–{high}/500"
        ax.text(xi, val + 1.4, count, ha="center", va="bottom", fontsize=10.5, fontweight="bold")
    display_labels = ["Single A", "Single B", "Equal-weight", "Continuous", "Double-significant-\nfirst"]
    ax.set_xticks(x, display_labels, ha="center")
    ax.set_ylabel("Supported per 100 nominations", fontsize=10)
    ax.set_ylim(0, 86); ax.set_yticks(np.arange(0, 81, 20))
    ax.grid(axis="y", color="#DDDDDD", linewidth=0.45); ax.set_axisbelow(True)
    ax.set_title("Protein support in the ProteinTalks/HPA summary", loc="center", weight="bold", fontsize=12)
    fig.text(0.015, 0.045, "Bar labels give counts out of 500. Whiskers show rounding bounds for published values.", ha="left", fontsize=8.6, color="#333333")
    fig.tight_layout(rect=[0.02, 0.10, 0.98, 0.93])
    save_figure(fig, "FigS6")


def make_contact_sheet(names: list[str]) -> None:
    rows = (len(names) + 1) // 2
    fig, axes = plt.subplots(rows, 2, figsize=(11.7, 3.5 * rows), squeeze=False)
    for ax in axes.flat:
        ax.axis("off")
    for ax, name in zip(axes.flat, names):
        img = plt.imread(OUT / f"{name}.png")
        ax.imshow(img)
        ax.set_title(name, loc="left", fontsize=10, weight="bold")
    fig.tight_layout(pad=1.2)
    fig.savefig(OUT / "contact_sheet.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


FIGURE_NAMES = ["Fig1", "Fig2", "Fig3", "Fig4", "FigS1", "FigS2", "FigS3", "FigS4", "FigS5", "FigS6"]


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Render the ten figures from bundled aggregate source tables."
    )
    parser.add_argument(
        "--data-dir",
        type=Path,
        default=DEFAULT_DATA_DIR,
        help="Directory containing the bundled figure TSV files (default: package source_data).",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help="Directory for PDF/PNG figures and the contact sheet (default: package output).",
    )
    parser.add_argument(
        "--figures",
        nargs="+",
        choices=FIGURE_NAMES,
        help="Optional subset to render; by default all ten figures are rendered.",
    )
    parser.add_argument(
        "--font-path",
        type=Path,
        help="Optional local Arial Bold font file to register for exact font embedding.",
    )
    args = parser.parse_args()

    if args.font_path:
        font_path = args.font_path.resolve()
        if not font_path.is_file():
            raise FileNotFoundError(font_path)
        font_manager.fontManager.addfont(str(font_path))
        font_name = font_manager.FontProperties(fname=str(font_path)).get_name()
        mpl.rcParams["font.family"] = [font_name]
        mpl.rcParams["font.sans-serif"] = [font_name]
        print(f"Registered bold font: {font_name} ({font_path.name})")

    global OUT
    OUT = args.output_dir.resolve()
    sources = load_sources(args.data_dir.resolve())
    plotters = {
        "Fig1": lambda: make_fig1(sources["Fig1_cohort_design"]),
        "Fig2": lambda: make_fig2(sources["Fig2_support_patient_bootstrap"]),
        "Fig3": lambda: make_fig3(sources["Fig3_fixed_budget_four_rules"]),
        "Fig4": lambda: make_fig4(sources["Fig4_decomposition_k500"]),
        "FigS1": lambda: make_fig_s1(sources["FigS1_sensitivity_grid"]),
        "FigS2": lambda: make_fig_s2(
            sources["FigS2_decomposition_patient_bootstrap"],
            sources["Fig1_cohort_design"],
        ),
        "FigS3": lambda: make_fig_s3(sources["FigS3_historical_purity_summary"]),
        "FigS4": lambda: make_fig_s4(sources["FigS4_shared_testable_contrast"]),
        "FigS5": lambda: make_fig_s5(sources["FigS5_BioTExt_signed_effect"]),
        "FigS6": lambda: make_fig_s6(sources["FigS6_HPA_fixed_k500"]),
    }
    selected = args.figures or FIGURE_NAMES
    for name in selected:
        style = SUPPLEMENT_STYLE if name.startswith("FigS") else {}
        with mpl.rc_context(style):
            plotters[name]()
    make_contact_sheet(selected)
    print(f"Rendered {len(selected)} figure(s) to {OUT}")
    print("Outputs: PDF and 300-dpi PNG per figure, plus contact_sheet.png")


if __name__ == "__main__":
    main()
