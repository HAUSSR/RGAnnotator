"""Shared plotting helpers (color palettes, colormap lookup, output format)."""

import os

import matplotlib
matplotlib.use("Agg")
import matplotlib as mpl

mpl.rcParams["svg.fonttype"] = "none"

# Okabe-Ito colorblind-safe palette
DEFAULT_COLORS = [
    "#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442",
    "#0072B2", "#D55E00", "#CC79A7",
]


def get_cmap(name):
    """Colormap lookup compatible with matplotlib < 3.9 and >= 3.9."""
    try:
        return mpl.colormaps[name]
    except AttributeError:  # matplotlib < 3.5 fallback
        return mpl.cm.get_cmap(name)


def save_figure(fig, output_path, png_dpi=300):
    """Save a figure, choosing the format from the output file extension."""
    ext = os.path.splitext(output_path)[1].lower().lstrip(".")
    fmt = {"svg": "svg", "png": "png", "pdf": "pdf", "jpg": "png", "jpeg": "png", "tiff": "png"}.get(ext, "png")
    dpi = png_dpi if fmt in ("png", "jpg", "jpeg", "tiff") else None
    fig.savefig(output_path, format=fmt, bbox_inches="tight", dpi=dpi)
