"""Module 3c - circadian time-course line plot.

Mean expression (+/- SD) per time point with a fitted two-harmonic curve
(24 h + 12 h OLS), JTK_CYCLE annotation (P, fold change, peak CT) and a
rhythmicity tag per gene. Supports overlaying multiple genes and multiple
datasets (files get colors, genes get linestyles-free distinct colors).
"""

import argparse
import os

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

try:
    from .timeparse import parse_ct_headers, build_manual_mapping
    from .rhythm import jtk_rhythmicity, fit_two_harmonic_ols, two_harmonic_predict
    from .plotlib import DEFAULT_COLORS, save_figure
except ImportError:  # allow running as a plain script
    from timeparse import parse_ct_headers, build_manual_mapping
    from rhythm import jtk_rhythmicity, fit_two_harmonic_ols, two_harmonic_predict
    from plotlib import DEFAULT_COLORS, save_figure


def fmt_p(p):
    return "n.s." if np.isnan(p) else (format(p, ".2e") if p < 0.001 else format(p, ".3f"))


def fmt_val(v, d=2):
    return "n.a." if np.isnan(v) else format(v, "." + str(d) + "f")


def run_lineplot(args):
    use_ct_mode = not args.zt_mode

    input_files = [f.strip() for f in args.input.split(",")]
    n_files = len(input_files)
    selected_genes = [g.strip() for g in args.genes.split(",")] if args.genes else []

    fig, ax = plt.subplots(figsize=(7, 6))
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    combined_annotations = []
    dfs = [pd.read_csv(f) for f in input_files]

    if not selected_genes:
        common_genes = set(dfs[0].iloc[:, 0])
        for df in dfs[1:]:
            common_genes.intersection_update(df.iloc[:, 0])
        selected_genes = sorted(common_genes)

    use_file_colors = n_files > 1
    global_time_points = set()

    for gene_idx, gene_name in enumerate(selected_genes):
        for file_idx, (df, file_path) in enumerate(zip(dfs, input_files)):
            gene_col = df.columns[0]
            row_data = df[df[gene_col] == gene_name]
            if row_data.empty:
                continue

            row = row_data.iloc[0]
            headers = df.columns[1:]

            # absolute times (no folding) so multi-day courses display in full
            if use_ct_mode:
                mapping = parse_ct_headers(headers, fold=False)
            else:
                mapping = build_manual_mapping(
                    len(headers), args.replicates, args.interval, args.zero_time, fold=False)

            if not mapping:
                continue
            unique_times = sorted(mapping.keys())
            for t in unique_times:
                global_time_points.add(t)

            raw_values_orig = row.iloc[1:].values.astype(float)

            mean_val = np.mean(raw_values_orig)
            std_val = np.std(raw_values_orig)

            if args.zscore and std_val > 1e-4:
                raw_values_plot = (raw_values_orig - mean_val) / std_val
            else:
                raw_values_plot = raw_values_orig

            # JTK test on all replicate values
            all_times, all_values = [], []
            for t in unique_times:
                for col_idx in mapping[t]:
                    all_times.append(float(t))
                    all_values.append(raw_values_orig[col_idx])
            all_times = np.array(all_times)
            all_values = np.array(all_values)

            p_jtk, peak_time, fold_change, amplitude, mesor = jtk_rhythmicity(all_times, all_values)
            fit_params = fit_two_harmonic_ols(all_times, all_values)

            is_rhythmic = (p_jtk < 0.05) and (fold_change >= args.fold_change)
            rhythmic_tag = "+" if is_rhythmic else "-"

            if use_file_colors:
                color = DEFAULT_COLORS[file_idx % len(DEFAULT_COLORS)]
                base_name = os.path.splitext(os.path.basename(file_path))[0]
                legend_label = gene_name + " (" + base_name + ")"
                file_tag = "[" + base_name + "] "
            else:
                color = DEFAULT_COLORS[gene_idx % len(DEFAULT_COLORS)]
                legend_label = gene_name
                file_tag = ""

            y_means_plot = [np.mean(raw_values_plot[mapping[t]]) for t in unique_times]
            y_stds_plot = [np.std(raw_values_plot[mapping[t]], ddof=1)
                           if len(raw_values_plot[mapping[t]]) > 1 else 0.0
                           for t in unique_times]
            t_array_plot = np.array(unique_times)
            y_array_plot = np.array(y_means_plot)
            y_err_plot = np.array(y_stds_plot)

            t_fine = np.linspace(min(unique_times), max(unique_times), 300)

            if not args.no_curve and fit_params is not None:
                y_fine = two_harmonic_predict(t_fine, fit_params)
                if y_fine is not None:
                    if args.zscore and std_val > 1e-4:
                        y_fine = (y_fine - mean_val) / std_val
                    ax.plot(t_fine, y_fine, color=color, linewidth=1.8,
                            alpha=0.45, linestyle="--", zorder=2)

            ax.errorbar(t_array_plot, y_array_plot, yerr=y_err_plot, marker="o", color=color,
                        linewidth=2.5, markersize=6, capsize=3, capthick=1.5,
                        elinewidth=1.2, label=legend_label, zorder=3)

            label_text = (file_tag + gene_name + " " + rhythmic_tag
                          + " | P=" + fmt_p(p_jtk) + ", FC=" + fmt_val(fold_change, 2)
                          + ", Peak=CT" + fmt_val(peak_time, 1))
            if args.show_all:
                label_text += ", Mesor=" + fmt_val(mesor, 2) + ", Amp=" + fmt_val(amplitude, 2)
            combined_annotations.append((label_text, color))

    n_annotations = len(combined_annotations)
    y_start = 0.96
    y_step = 0.035
    if n_annotations > 8:
        y_step = max(0.02, 0.8 / n_annotations)
    for i, (txt, col) in enumerate(combined_annotations):
        ax.text(0.98, y_start - (i * y_step), txt, transform=ax.transAxes, color=col,
                fontsize=8, fontweight="bold", ha="right", va="top", family="monospace")

    ax.set_xlabel("Circadian time (CT)", fontsize=12, fontweight="bold")
    if args.zscore:
        ax.set_ylabel("Z-score", fontsize=12, fontweight="bold")
    else:
        ax.set_ylabel("Relative Expression", fontsize=12, fontweight="bold")

    if global_time_points:
        x_min = int(np.floor(min(global_time_points)))
        x_max = int(np.ceil(max(global_time_points)))
        ax.set_xticks(np.arange(x_min, x_max + 1, 3))
    ax.grid(axis="y", linestyle="--", alpha=0.3)
    ax.legend(loc="upper left", frameon=False, fontsize=9)
    plt.tight_layout()

    save_figure(fig, args.output, png_dpi=300)
    plt.close()
    print("Success! Line plot saved to: " + args.output)


def build_parser(add_help=True):
    parser = argparse.ArgumentParser(
        prog="rgannotator line",
        add_help=add_help,
        description="Module 3c - circadian time-course line plot with fitted curves",
    )
    parser.add_argument("-i", "--input", required=True, help="Input CSV path(s), comma-separated")
    parser.add_argument("-o", "--output", required=True, help="Output image (SVG/PNG/PDF)")
    parser.add_argument("-g", "--genes", type=str, default=None,
                        help="Comma-separated gene names. Default: all common genes.")
    parser.add_argument("-z", "--zscore", action="store_true",
                        help="Plot as Z-scores (statistics still use original data)")
    parser.add_argument("-nc", "--no_curve", action="store_true",
                        help="Hide fitted curves (data points only)")
    parser.add_argument("-fc", "--fold_change", type=float, default=1.5,
                        help="Minimum fold change for the rhythmic tag (default: 1.5)")
    parser.add_argument("-all", "--show_all", action="store_true",
                        help="Show Mesor and Amplitude in annotations")
    parser.add_argument("-ZT", "--zt_mode", action="store_true",
                        help="Manual mode instead of CT-column parsing")
    parser.add_argument("-zt0", "--zero_time", type=float, default=0)
    parser.add_argument("-r", "--replicates", type=int, default=1)
    parser.add_argument("-int", "--interval", type=int, default=4)
    return parser


def main(argv=None):
    args = build_parser().parse_args(argv)
    run_lineplot(args)


if __name__ == "__main__":
    main()
