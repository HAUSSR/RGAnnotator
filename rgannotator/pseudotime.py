"""Module 3a - circadian pseudotime expression plot (CT/ZT as pseudotime).

Rows are genes (optionally ordered by peak time or hierarchical clustering),
columns are circadian time. Profiles are period-folded onto 0-24 h, cubic-
interpolated, Gaussian-smoothed and row Z-scored. Optional features: mark
specific genes on the left, attach a category annotation strip on the right,
and duplicate the cycle (-nc 2) for wrap-around display.
"""

import argparse
import os
import warnings

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.patches import Rectangle
from scipy.ndimage import gaussian_filter1d
from scipy.interpolate import interp1d
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import pdist

warnings.filterwarnings("ignore")

try:
    from .timeparse import parse_ct_headers, build_manual_mapping, detect_ct_mode
    from .plotlib import get_cmap, save_figure
except ImportError:  # allow running as a plain script
    from timeparse import parse_ct_headers, build_manual_mapping, detect_ct_mode
    from plotlib import get_cmap, save_figure


# ========================
# 1. Data loading
# ========================
def load_data(filepath, use_ct, zt0, interval, replicates, gene_filter):
    df = pd.read_csv(filepath)
    gene_col = df.columns[0]
    headers = df.columns[1:]

    if gene_filter:
        df = df[df[gene_col].isin(gene_filter)].copy()
        if df.empty:
            return None, None, None

    # Time parsing + folding of multi-day sampling onto 0-24 h
    if use_ct:
        time_mapping = parse_ct_headers(headers, fold=True)
    else:
        time_mapping = build_manual_mapping(len(headers), replicates, interval, zt0, fold=True)

    if not time_mapping:
        return None, None, None

    unique_times = sorted(time_mapping.keys())
    gene_names = df[gene_col].values

    # average replicates
    expr_matrix = np.zeros((len(df), len(unique_times)))
    for j, t in enumerate(unique_times):
        indices = time_mapping[t]
        vals = df.iloc[:, [idx + 1 for idx in indices]].values.astype(float)
        expr_matrix[:, j] = np.mean(vals, axis=1)

    return expr_matrix, unique_times, gene_names


# ========================
# 2. Periodic interpolation + Gaussian smoothing
# ========================
def smooth_periodic(expr_matrix, times, n_bins=240, sigma=2.0):
    n_genes = expr_matrix.shape[0]
    times_arr = np.array(times)

    # replicate one period left/right so the boundary stays continuous
    ext_times = np.concatenate([times_arr - 24, times_arr, times_arr + 24])
    ext_matrix = np.tile(expr_matrix, (1, 3))

    target_times = np.linspace(0, 24, n_bins, endpoint=False)

    smoothed = np.zeros((n_genes, n_bins))
    for i in range(n_genes):
        f = interp1d(ext_times, ext_matrix[i, :], kind="cubic", fill_value="extrapolate")
        smoothed[i, :] = f(target_times)

    smoothed = gaussian_filter1d(smoothed, sigma=sigma, axis=1, mode="wrap")
    return smoothed, target_times


# ========================
# 3. Row Z-score
# ========================
def zscore_rows(matrix):
    mean = np.mean(matrix, axis=1, keepdims=True)
    std = np.std(matrix, axis=1, keepdims=True)
    std[std < 1e-8] = 1.0
    return (matrix - mean) / std


# ========================
# 4. Gene ordering
# ========================
def order_genes(z_matrix, times, method="peak"):
    if method == "cluster":
        if z_matrix.shape[0] > 1:
            Z = linkage(pdist(z_matrix, "euclidean"), method="average")
            return leaves_list(Z)
        return np.array([0])
    peak_idx = np.argmax(z_matrix, axis=1)
    return np.argsort(times[peak_idx])


# ========================
# 5. Category palette
# ========================
def make_category_colors(categories):
    palettes = []
    for cmap_name in ["tab20", "tab20b", "tab20c", "Set3", "Paired", "Accent", "Dark2"]:
        cmap = get_cmap(cmap_name)
        if hasattr(cmap, "colors"):
            palettes.extend(list(cmap.colors))
        else:
            palettes.extend([cmap(i / 20) for i in range(20)])

    uniq, seen = [], set()
    for c in palettes:
        key = tuple(np.round(c[:3], 4))
        if key not in seen:
            uniq.append(c)
            seen.add(key)

    if len(categories) > len(uniq):
        extra = len(categories) - len(uniq)
        uniq.extend([plt.cm.hsv(i / max(extra, 1)) for i in range(extra)])

    return {cat: uniq[i] for i, cat in enumerate(sorted(categories))}


# ========================
# 6. Single pseudotime panel
# ========================
def plot_panel(fig, gs_cell, heat_z, gene_names, times_smooth,
                       annotation_df=None, anno_col=None, mark_genes=None,
                       title="Pseudotime", cmap_name="RdBu_r", z_clip=2.5, n_cycle=2):
    if n_cycle == 2:
        heat_z_plot = np.concatenate([heat_z, heat_z], axis=1)
    else:
        heat_z_plot = heat_z.copy()

    heat_z_plot = np.clip(heat_z_plot, -z_clip, z_clip)
    has_anno = annotation_df is not None and anno_col is not None
    has_mark = mark_genes is not None and len(mark_genes) > 0

    width_ratios = []
    if has_mark:
        width_ratios.append(0.25)
    width_ratios.append(1.0)
    if has_anno:
        width_ratios.append(0.05)

    inner = gs_cell.subgridspec(
        1, len(width_ratios),
        width_ratios=width_ratios,
        wspace=0.02,
    )

    col_idx = 0
    ax_mark = None
    if has_mark:
        ax_mark = fig.add_subplot(inner[0, col_idx])
        col_idx += 1

    ax_hm = fig.add_subplot(inner[0, col_idx])
    col_idx += 1

    ax_ann = None
    if has_anno:
        ax_ann = fig.add_subplot(inner[0, col_idx])
        col_idx += 1

    cmap = get_cmap(cmap_name)
    im = ax_hm.imshow(
        heat_z_plot, aspect="auto", interpolation="nearest",
        cmap=cmap, vmin=-z_clip, vmax=z_clip, origin="upper",
    )
    ax_hm.set_title(title, fontsize=11, fontweight="bold", pad=8)
    ax_hm.set_yticks([])

    n_genes = heat_z_plot.shape[0]
    ax_hm.set_ylim(n_genes - 0.5, -0.5)

    n_cols = heat_z_plot.shape[1]
    max_ct = 48 if n_cycle == 2 else 24
    tick_times = range(0, max_ct + 1, 6)
    ax_hm.set_xticks([tt / max_ct * (n_cols - 1) for tt in tick_times])
    ax_hm.set_xticklabels(["CT" + str(tt) for tt in tick_times], fontsize=8)
    ax_hm.set_xlabel("Circadian Time", fontsize=9)
    ax_hm.tick_params(length=0)
    for spine in ax_hm.spines.values():
        spine.set_visible(False)

    # left-side markers for genes of interest
    if has_mark:
        ax_mark.set_ylim(n_genes - 0.5, -0.5)
        ax_mark.set_xlim(0, 1)
        ax_mark.axis("off")
        mark_set = set(mark_genes)
        for i, g in enumerate(gene_names):
            if g in mark_set:
                ax_mark.text(0.95, i, g, ha="right", va="center", fontsize=7,
                             color="black", fontweight="bold", fontstyle="italic")
                ax_mark.scatter(0.98, i, marker=">", color="red", s=25, clip_on=False, zorder=10)

    # right-side category strip
    cat_colors, cats = None, None
    if has_anno:
        gene_to_cat = dict(zip(annotation_df.iloc[:, 0].astype(str), annotation_df[anno_col]))
        categories = [gene_to_cat.get(str(g), "Unknown") for g in gene_names]
        cats = sorted(set(categories))
        cat_colors = make_category_colors(cats)
        cat_codes = [cats.index(c) for c in categories]

        ax_ann.imshow(
            np.array(cat_codes).reshape(-1, 1), aspect="auto",
            interpolation="nearest", cmap=ListedColormap([cat_colors[c] for c in cats]),
            origin="upper",
        )
        ax_ann.set_ylim(n_genes - 0.5, -0.5)
        ax_ann.set_xticks([])
        ax_ann.set_yticks([])
        ax_ann.set_title(anno_col, fontsize=7, pad=4, rotation=90)
        for spine in ax_ann.spines.values():
            spine.set_visible(False)

    return im, cat_colors, cats


# ========================
# 7. Main
# ========================
def run_pseudotime(args):
    use_ct = not args.zt_mode
    input_files = [f.strip() for f in args.input.split(",")]
    n_files = len(input_files)

    gene_filter = None
    if args.gene_list:
        with open(args.gene_list) as f:
            gene_filter = set(line.strip() for line in f if line.strip())
        print("Loaded " + str(len(gene_filter)) + " genes from filter file")

    mark_genes = None
    if args.mark_genes:
        mark_genes = [g.strip() for g in args.mark_genes.split(",")]
        print("Will mark genes: " + str(mark_genes))

    annotation_df, anno_col = None, None
    if args.annotation:
        annotation_df = pd.read_csv(args.annotation)
        anno_col = annotation_df.columns[1]
        print("Annotation: " + str(len(annotation_df)) + " entries, category: '" + anno_col + "'")

    has_anno = annotation_df is not None

    # layout
    if n_files == 1:
        nrows, ncols = 1, 1
    elif n_files == 2:
        nrows, ncols = 1, 2
    elif n_files <= 4:
        nrows, ncols = 1, 3
    else:
        nrows, ncols = (n_files + 2) // 3, 3

    panel_w, panel_h = 4.5, 8
    fig_w = panel_w * ncols + (3 if has_anno else 1.5)
    fig_h = panel_h * nrows + 1.5
    fig = plt.figure(figsize=(fig_w, fig_h))

    if has_anno:
        total_cols = ncols + 2
        outer = fig.add_gridspec(1, total_cols, width_ratios=[1] * ncols + [0.04, 0.22], wspace=0.25)
        gs_heatmaps = [outer[0, i] for i in range(ncols)]
        gs_cbar = outer[0, ncols]
        gs_legend = outer[0, ncols + 1]
    else:
        outer = fig.add_gridspec(nrows, ncols, wspace=0.35, hspace=0.4)
        gs_heatmaps = [outer[r, c] for r in range(nrows) for c in range(ncols)]
        gs_cbar = gs_legend = None

    all_cat_colors, all_cats, last_im = None, None, None

    for file_idx, filepath in enumerate(input_files):
        if file_idx >= len(gs_heatmaps):
            break
        gs_cell = gs_heatmaps[file_idx]
        print("\nProcessing: " + filepath)

        expr_matrix, times, gene_names = load_data(
            filepath, use_ct, args.zero_time, args.interval,
            args.replicates, gene_filter,
        )
        if expr_matrix is None:
            print("  Skipping: no valid data")
            continue

        print("  Genes: " + str(len(gene_names)) + ", Time points: " + str(times))

        smoothed, times_smooth = smooth_periodic(expr_matrix, times, args.n_bins, args.sigma)
        z_matrix = zscore_rows(smoothed)
        order = order_genes(z_matrix, times_smooth, args.order_by)

        base_name = os.path.splitext(os.path.basename(filepath))[0]
        im, cat_colors, cats = plot_panel(
            fig, gs_cell,
            z_matrix[order], gene_names[order], times_smooth,
            annotation_df=annotation_df, anno_col=anno_col,
            mark_genes=mark_genes,
            title=base_name, cmap_name=args.colormap,
            z_clip=args.z_clip, n_cycle=args.n_cycle,
        )
        last_im = im
        if cat_colors:
            all_cat_colors, all_cats = cat_colors, cats

    # hide unused panels
    for idx in range(n_files, nrows * ncols):
        if idx < len(gs_heatmaps):
            fig.add_subplot(gs_heatmaps[idx]).axis("off")

    # colorbar
    if gs_cbar is not None:
        cax = fig.add_subplot(gs_cbar)
        cb = fig.colorbar(last_im, cax=cax)
        cb.set_label("Row Z-score", fontsize=10, fontweight="bold")
        cb.ax.tick_params(labelsize=8)
    else:
        fig.colorbar(last_im, ax=fig.axes, fraction=0.02, pad=0.02, label="Row Z-score")

    # category legend
    if gs_legend is not None and all_cats and all_cat_colors:
        lax = fig.add_subplot(gs_legend)
        lax.axis("off")
        lax.set_title(anno_col, fontsize=10, fontweight="bold", pad=8)
        dy = 1.0 / max(len(all_cats), 1)
        for i, cat in enumerate(all_cats):
            y = 1.0 - (i + 1) * dy + 0.005
            lax.add_patch(Rectangle((0, y), 0.15, dy * 0.7, color=all_cat_colors[cat],
                                    transform=lax.transAxes, clip_on=False))
            lax.text(0.20, y + dy * 0.35, cat, transform=lax.transAxes, va="center", fontsize=7)

    fig.suptitle("Circadian Pseudotime Expression", fontsize=14, fontweight="bold", y=0.98)
    plt.tight_layout(rect=[0, 0, 1, 0.95])

    save_figure(fig, args.output, png_dpi=1200)
    plt.close()

    print("\nSaved to: " + args.output)
    print("   Order: " + args.order_by + ", Cycles: " + str(args.n_cycle) + ", Cmap: " + args.colormap)


def build_parser(add_help=True):
    parser = argparse.ArgumentParser(
        prog="rgannotator pseudotime",
        add_help=add_help,
        description="Module 3a - circadian pseudotime expression plot (CT as pseudotime)",
    )
    parser.add_argument("-i", "--input", required=True, help="Input CSV, comma-separated for multiple panels")
    parser.add_argument("-o", "--output", required=True, help="Output image (PNG/SVG/PDF)")
    parser.add_argument("-g", "--mark_genes", type=str, default=None,
                        help="Comma-separated gene names to mark on the left (e.g. -g MsLHY)")
    parser.add_argument("-gl", "--gene_list", type=str, default=None,
                        help="Gene list file to filter rows (one gene per line)")
    parser.add_argument("-anno", "--annotation", type=str, default=None,
                        help="Annotation CSV: col1=gene, col2=classification")
    parser.add_argument("-cmap", "--colormap", type=str, default="RdBu_r")
    parser.add_argument("-zclip", "--z_clip", type=float, default=2.5)
    parser.add_argument("-bins", "--n_bins", type=int, default=240)
    parser.add_argument("-sigma", "--sigma", type=float, default=2.0)
    parser.add_argument("-nc", "--n_cycle", type=int, default=1, choices=[1, 2])
    parser.add_argument("-order", "--order_by", choices=["peak", "cluster"], default="peak")
    parser.add_argument("-ZT", "--zt_mode", action="store_true",
                        help="Manual mode instead of CT-column parsing")
    parser.add_argument("-zt0", "--zero_time", type=float, default=0)
    parser.add_argument("-int", "--interval", type=int, default=4)
    parser.add_argument("-r", "--replicates", type=int, default=1)
    return parser


def main(argv=None):
    args = build_parser().parse_args(argv)
    run_pseudotime(args)


if __name__ == "__main__":
    main()
