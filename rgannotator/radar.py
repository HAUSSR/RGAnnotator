"""Module 3b - 24 h circadian radar (polar) plot.

All sampled time points are folded onto a single 0-24 h day (CT labels wrap
modulo 24), replicates are averaged per time point and drawn as a closed
polygon around the clock face. Each polygon is annotated with the JTK_CYCLE
rhythmicity test (P, fold change, peak CT) and the fitted peak is starred.
Ideal for comparing several genes (or the same gene across datasets) on one
clock face.
"""

import argparse
import os
import re

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


def run_radar(args):
    use_ct_mode = not args.zt_mode

    input_files = [f.strip() for f in args.input.split(",")]
    n_files = len(input_files)
    selected_genes = [g.strip() for g in args.genes.split(",")] if args.genes else []

    fig, ax = plt.subplots(figsize=(8, 8), subplot_kw=dict(polar=True))
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)

    combined_annotations = []
    dfs = [pd.read_csv(f) for f in input_files]

    if not selected_genes:
        common_genes = set(dfs[0].iloc[:, 0])
        for df in dfs[1:]:
            common_genes.intersection_update(df.iloc[:, 0])
        selected_genes = sorted(common_genes)

    use_file_colors = n_files > 1
    global_folded_times = set()

    # first pass: collect folded time points to set the axis ticks
    parsed = []
    for gene_idx, gene_name in enumerate(selected_genes):
        for file_idx, (df, file_path) in enumerate(zip(dfs, input_files)):
            gene_col = df.columns[0]
            row_data = df[df[gene_col] == gene_name]
            if row_data.empty:
                continue
            row = row_data.iloc[0]
            headers = df.columns[1:]

            if use_ct_mode:
                fold_mapping = parse_ct_headers(headers, fold=True)
            else:
                fold_mapping = build_manual_mapping(
                    len(headers), args.replicates, args.interval, args.zero_time, fold=True)
            if not fold_mapping:
                continue
            parsed.append((gene_idx, gene_name, file_idx, file_path, row, headers, fold_mapping))
            global_folded_times.update(fold_mapping.keys())

    # axis ticks from the actual sampling grid (not hardcoded 3 h)
    folded_times = sorted(global_folded_times)
    if folded_times:
        angles = np.array(folded_times) / 24 * 2 * np.pi
        ax.set_xticks(angles)
        ax.set_xticklabels(["CT" + (str(int(t)) if t == int(t) else str(t)) for t in folded_times],
                           fontsize=11 if len(folded_times) <= 8 else 8, fontweight="bold")

    ax.yaxis.grid(True, color="darkgray", alpha=0.7, linestyle="--", linewidth=1.2)
    ax.xaxis.grid(True, color="darkgray", alpha=0.8, linestyle="-", linewidth=1.2)
    ax.spines["polar"].set_visible(True)
    ax.spines["polar"].set_color("black")
    ax.spines["polar"].set_linewidth(1.5)
    ax.set_yticklabels([])

    for gene_idx, gene_name, file_idx, file_path, row, headers, fold_mapping in parsed:
        unique_times_24 = sorted(fold_mapping.keys())

        # statistics run on original (unfolded) times
        all_times, all_values = [], []
        raw_values = row.iloc[1:].values.astype(float)
        plot_means_24 = []
        for t_24 in unique_times_24:
            indices = fold_mapping[t_24]
            for idx in indices:
                if use_ct_mode:
                    m = re.search(r"CT(\d+\.?\d*)", str(headers[idx]))
                    if m:
                        all_times.append(float(m.group(1)))
                else:
                    all_times.append(args.zero_time + (idx // args.replicates) * args.interval)
            vals = [raw_values[idx] for idx in indices]
            all_values.extend(vals)
            plot_means_24.append(np.mean(vals))

        all_times = np.array(all_times)
        all_values = np.array(all_values)
        plot_means_24 = np.array(plot_means_24)

        mean_val = np.mean(all_values)
        std_val = np.std(all_values)
        if args.zscore and std_val > 1e-4:
            plot_means_24 = (plot_means_24 - mean_val) / std_val

        p_jtk, peak_time, fold_change, amplitude, mesor = jtk_rhythmicity(all_times, all_values)
        fit_params = fit_two_harmonic_ols(all_times, all_values)

        is_rhythmic = (p_jtk < 0.05) and (fold_change >= args.fold_change)
        rhythmic_tag = "+" if is_rhythmic else "-"

        if use_file_colors:
            color = DEFAULT_COLORS[(file_idx + 1) % len(DEFAULT_COLORS)]
            base_name = os.path.splitext(os.path.basename(file_path))[0]
            legend_label = gene_name + " (" + base_name + ")"
            file_tag = "[" + base_name + "] "
        else:
            color = DEFAULT_COLORS[1:][gene_idx % (len(DEFAULT_COLORS) - 1)]
            legend_label = gene_name
            file_tag = ""

        # polar plots cannot show negative radii (they mirror to the opposite
        # angle), so shift the display space when needed
        display_vals = plot_means_24
        r_offset = 0.0
        if display_vals.size and display_vals.min() < 0:
            data_range = display_vals.max() - display_vals.min()
            r_offset = -display_vals.min() + 0.1 * (data_range if data_range > 0 else 1.0)
            display_vals = display_vals + r_offset

        theta_points = np.array(unique_times_24) / 24 * 2 * np.pi
        theta_closed = np.append(theta_points, theta_points[0])
        y_closed = np.append(display_vals, display_vals[0])

        ax.plot(theta_closed, y_closed, color=color, linewidth=2.5, label=legend_label, zorder=3)
        ax.fill(theta_closed, y_closed, color=color, alpha=0.2, zorder=2)
        ax.scatter(theta_points, display_vals, color=color, s=50, zorder=4,
                   edgecolor="white", linewidth=1.5)

        # fitted-peak star
        if not np.isnan(peak_time) and fit_params is not None:
            peak_y = two_harmonic_predict(peak_time, fit_params)
            if peak_y is not None:
                if args.zscore and std_val > 1e-4:
                    peak_y = (peak_y - mean_val) / std_val
                peak_y = peak_y + r_offset
                peak_theta = peak_time / 24 * 2 * np.pi
                ax.scatter([peak_theta], [peak_y], color=color, marker="*", s=250,
                           zorder=5, edgecolor="black", linewidth=0.8)

        label_text = (file_tag + gene_name + " " + rhythmic_tag
                      + " | P=" + fmt_p(p_jtk) + ", FC=" + fmt_val(fold_change, 2)
                      + ", Peak=CT" + fmt_val(peak_time, 1))
        if args.show_all:
            label_text += ", Mesor=" + fmt_val(mesor, 2) + ", Amp=" + fmt_val(amplitude, 2)
        combined_annotations.append((label_text, color))

    # annotation block
    n_annotations = len(combined_annotations)
    y_start = 1.05
    y_step = 0.04
    if n_annotations > 8:
        y_step = max(0.02, 0.8 / n_annotations)
    for i, (txt, col) in enumerate(combined_annotations):
        fig.text(0.99, y_start - (i * y_step), txt, color=col, fontsize=8,
                 fontweight="bold", ha="right", va="top", family="monospace")

    fig.text(0.5, 0.98, "Circadian Expression Radar (24h Aggregated)",
             ha="center", fontsize=14, fontweight="bold")
    ax.legend(loc="upper left", bbox_to_anchor=(-0.1, 1.1), frameon=False, fontsize=9)

    plt.tight_layout()
    save_figure(fig, args.output, png_dpi=300)
    plt.close()
    print("Success! Radar plot saved to: " + args.output)


def build_parser(add_help=True):
    parser = argparse.ArgumentParser(
        prog="rgannotator radar",
        add_help=add_help,
        description="Module 3b - 24h circadian radar (polar) plot",
    )
    parser.add_argument("-i", "--input", required=True, help="Input CSV path(s), comma-separated")
    parser.add_argument("-o", "--output", required=True, help="Output image (SVG/PNG/PDF)")
    parser.add_argument("-g", "--genes", type=str, default=None,
                        help="Comma-separated gene names. Default: all common genes.")
    parser.add_argument("-z", "--zscore", action="store_true",
                        help="Plot as Z-scores (recommended when comparing genes/datasets)")
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
    run_radar(args)


if __name__ == "__main__":
    main()
