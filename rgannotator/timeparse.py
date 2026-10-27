"""Shared time-design parsing for expression matrices.

An RGAnnotator expression table is a CSV whose first column holds gene IDs and
whose remaining columns hold expression values. Two column layouts are supported:

CT mode      headers encode circadian time directly, e.g. ``CT0, CT4, CT8, ...
             CT15_2`` (optional ``_replicate`` suffix). Multi-day sampling such
             as ``CT24 ... CT72`` is accepted; pass ``fold=True`` to wrap times
             onto a single 0-24 h day (used by the radar and heatmap modules),
             or ``fold=False`` to keep absolute times (used by the statistics).
Manual mode  headers carry no time information; the time series is described by
             ``zero_time`` (start ZT), ``interval`` (hours between points) and
             ``replicates`` (consecutive columns per time point).
"""

import re

CT_PATTERN = re.compile(r"CT(\d+\.?\d*)")


def parse_ct_headers(headers, fold=False):
    """Map CT-style headers to time values.

    Returns ``{time_value: [column_index, ...]}`` where indices are relative to
    the first expression column (0-based). Non-matching headers are skipped.
    With ``fold=True`` times are reduced modulo 24 so extra circadian days
    merge onto one 0-24 h cycle.
    """
    mapping = {}
    for i, h in enumerate(headers):
        match = CT_PATTERN.search(str(h))
        if match:
            t_val = float(match.group(1))
            if fold:
                t_val %= 24
            if t_val not in mapping:
                mapping[t_val] = []
            mapping[t_val].append(i)
    return mapping


def build_manual_mapping(n_expr_cols, replicates, interval, zero_time, fold=False):
    """Map consecutive columns to evenly spaced time points.

    Returns ``{time_value: [column_index, ...]}`` mirroring :func:`parse_ct_headers`.
    """
    mapping = {}
    num_tp = n_expr_cols // replicates
    for i in range(num_tp):
        t_val = zero_time + i * interval
        if fold:
            t_val %= 24
        if t_val not in mapping:
            mapping[t_val] = []
        mapping[t_val].extend(range(i * replicates, (i + 1) * replicates))
    return mapping


def detect_ct_mode(headers):
    """True if at least two headers look like CT time labels."""
    return sum(1 for h in headers if CT_PATTERN.search(str(h))) >= 2
