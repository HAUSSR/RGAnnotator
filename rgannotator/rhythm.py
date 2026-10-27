"""Module 2 - circadian rhythm detection for expression matrices.

Two statistical engines are integrated and share one command-line interface:

* ``jtk`` (default): JTK_CYCLE non-parametric test. Kendall's tau between the
  expression values and a 24 h reference cosine scanned over 24 phases
  (Bonferroni-corrected across phases), computed on ALL individual replicate
  values. Robust to outliers and non-normal noise.
* ``cosinor``: single-cosine (24 h period) least-squares fit with an OLS
  F-test on the first-harmonic design matrix, fitted on per-time-point means.

Both engines return peak time (CT, 0-24 h), mesor, amplitude and fold change
(max/min time-point mean), and the runner adds BH-FDR correction plus a dual
rhythmicity criterion (Adj_P < 0.05 AND Fold_Change >= threshold).
"""

import argparse

import numpy as np
import pandas as pd
from concurrent.futures import ProcessPoolExecutor
from functools import partial
from scipy.stats import kendalltau
from scipy.optimize import curve_fit
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests

try:
    from .timeparse import parse_ct_headers, build_manual_mapping, detect_ct_mode
except ImportError:  # allow running as a plain script
    from timeparse import parse_ct_headers, build_manual_mapping, detect_ct_mode


# ========================
# JTK_CYCLE engine (non-parametric, all replicate values)
# ========================
def jtk_rhythmicity(all_times, all_values):
    """JTK_CYCLE scan over 24 phases. Returns (p, peak_time, fold_change, amplitude, mesor)."""
    if np.std(all_values) < 1e-4:
        return 1.0, np.nan, 1.0, 0.0, np.nan

    n_phases = 24  # 1 h resolution
    phases = np.linspace(0, 24, n_phases, endpoint=False)

    best_p = 1.0
    best_phase = 0.0
    best_tau = 0.0

    for phase in phases:
        ref = np.cos(2 * np.pi * (all_times - phase) / 24)
        tau, p = kendalltau(all_values, ref)
        if not np.isnan(p) and p < best_p:
            best_p = p
            best_phase = phase
            best_tau = tau

    p_corrected = min(1.0, best_p * n_phases)  # Bonferroni across phases
    peak_time = best_phase % 24 if best_tau >= 0 else (best_phase + 12) % 24

    # Effect sizes from time-point means
    unique_times = sorted(set(all_times.tolist()))
    tp_means = [np.mean(all_values[all_times == t]) for t in unique_times]
    max_mean = max(tp_means)
    min_mean = min(tp_means)

    amplitude = max_mean - min_mean
    mesor = np.mean(all_values)

    if min_mean > 0:
        fold_change = max_mean / min_mean
    elif max_mean > 0:
        fold_change = float("inf")
    else:
        fold_change = 1.0

    return p_corrected, peak_time, fold_change, amplitude, mesor


# ========================
# Cosinor engine (single 24 h cosine fit)
# ========================
def cosinor_func(time, mesor, amplitude, acrophase):
    return mesor + amplitude * np.cos(2 * np.pi * time / 24 - acrophase)


def cosinor_rhythmicity(times, expression):
    """Fit y = M + A*cos(2*pi*t/24 - phi) on per-time-point means and test the
    first harmonic with OLS. Returns (p, peak_time, fold_change, amplitude, mesor)
    with the same layout as :func:`jtk_rhythmicity`."""
    if np.std(expression) < 1e-4:
        return 1.0, np.nan, 1.0, 0.0, np.nan

    p0 = [np.mean(expression), (np.max(expression) - np.min(expression)) / 2, 0]
    params, _ = curve_fit(cosinor_func, times, expression, p0=p0, maxfev=10000)
    mesor, amplitude, acrophase = params

    omega = 2 * np.pi / 24
    X = np.column_stack([np.cos(omega * times), np.sin(omega * times)])
    X = sm.add_constant(X)
    ols_res = sm.OLS(expression, X).fit()
    p_val = ols_res.f_pvalue

    peak_hour = (acrophase * 24) / (2 * np.pi) % 24

    tp_means = expression
    max_mean, min_mean = np.max(tp_means), np.min(tp_means)
    if min_mean > 0:
        fold_change = max_mean / min_mean
    elif max_mean > 0:
        fold_change = float("inf")
    else:
        fold_change = 1.0

    return p_val, peak_hour, fold_change, abs(amplitude), mesor


# ========================
# Two-harmonic OLS (used for fitted curves in plots)
# ========================
def fit_two_harmonic_ols(all_times, all_values):
    """OLS with 24 h + 12 h harmonics; used to draw smooth fitted curves."""
    try:
        if np.std(all_values) < 1e-4:
            return None
        omega1 = 2 * np.pi / 24
        omega2 = 4 * np.pi / 24
        X = np.column_stack([
            np.cos(omega1 * all_times), np.sin(omega1 * all_times),
            np.cos(omega2 * all_times), np.sin(omega2 * all_times),
        ])
        X = sm.add_constant(X)
        ols_res = sm.OLS(all_values, X).fit()
        return tuple(ols_res.params)
    except Exception:
        return None


def two_harmonic_predict(time, params):
    if params is None:
        return None
    M, A1, B1, A2, B2 = params
    omega1 = 2 * np.pi / 24
    omega2 = 4 * np.pi / 24
    return (M + A1 * np.cos(omega1 * time) + B1 * np.sin(omega1 * time)
            + A2 * np.cos(omega2 * time) + B2 * np.sin(omega2 * time))


# ========================
# Per-gene workers
# ========================
def _empty_result(gene_id):
    return {
        "Gene_ID": gene_id, "P_value": np.nan, "Peak_Time": np.nan,
        "Fold_Change": np.nan, "Amplitude": np.nan, "Mesor": np.nan,
        "Rhythmic": False,
    }


def process_gene(row_data, time_points, mapping, num_replicates,
                 is_ct_mode, method, min_fold_change):
    gene_id = row_data[0]
    raw_values = np.array(row_data[1:], dtype=float)
    res = _empty_result(gene_id)

    try:
        if is_ct_mode:
            # raw_values indices line up with the header mapping
            all_times, all_values = [], []
            for t in sorted(mapping.keys()):
                for idx in mapping[t]:
                    all_times.append(float(t))
                    all_values.append(raw_values[idx])
            all_times = np.array(all_times)
            all_values = np.array(all_values)
        else:
            num_tp = len(time_points)
            all_times = np.repeat(time_points, num_replicates)
            all_values = raw_values

        if method == "cosinor":
            # cosinor fits per-time-point means
            unique_times = sorted(set(all_times.tolist()))
            expression = np.array([np.mean(all_values[all_times == t]) for t in unique_times])
            calc_times = np.array(unique_times)
            p_val, peak_time, fold_change, amp, mesor = cosinor_rhythmicity(calc_times, expression)
        else:
            p_val, peak_time, fold_change, amp, mesor = jtk_rhythmicity(all_times, all_values)

        is_rhythmic = (p_val < 0.05) and (fold_change >= min_fold_change)

        return {
            "Gene_ID": gene_id,
            "P_value": p_val,
            "Peak_Time": round(peak_time, 2) if not np.isnan(peak_time) else np.nan,
            "Fold_Change": round(fold_change, 4) if np.isfinite(fold_change) else float("inf"),
            "Amplitude": round(amp, 4),
            "Mesor": round(mesor, 4),
            "Rhythmic": is_rhythmic,
        }
    except Exception:
        return res


# ========================
# Runner
# ========================
def run_rhythm(args):
    df = pd.read_csv(args.input)
    gene_rows = df.values.tolist()
    headers = df.columns[1:]

    # Time layout: force CT with -CT, manual with -zt-mode, otherwise auto-detect
    if args.ct_mode:
        is_ct = True
    elif args.zt_mode:
        is_ct = False
    else:
        is_ct = detect_ct_mode(headers)

    mapping = None
    time_points = None
    if is_ct:
        mapping = parse_ct_headers(headers, fold=False)
        if len(mapping) < 3:
            raise SystemExit("ERROR: CT mode selected but fewer than 3 CT columns found in " + args.input)
        unique_times = sorted(mapping.keys())
        reps = [len(mapping[t]) for t in unique_times]
        print("--- CT Mode ---")
        print("Detected time points: " + str(unique_times))
        print("Replicates per time point: " + str(reps))
    else:
        num_tp = (df.shape[1] - 1) // args.replicates
        time_points = np.arange(args.zero_time, args.zero_time + num_tp * args.interval, args.interval)
        print("--- Manual Mode ---")
        print("Calculated time points: " + str(time_points))
        print("Replicates per time point: " + str(args.replicates))

    method_name = "JTK_CYCLE (non-parametric, all replicates)" if args.method == "jtk" \
        else "Cosinor (single 24 h cosine fit)"
    print("--- Analysis Parameters ---")
    print("Method: " + method_name)
    print("P-value correction: " + ("phase-Bonferroni + BH-FDR" if args.method == "jtk" else "BH-FDR"))
    print("Effect size filter: FC >= " + str(args.fold_change))

    worker = partial(
        process_gene, time_points=time_points, mapping=mapping,
        num_replicates=args.replicates, is_ct_mode=is_ct,
        method=args.method, min_fold_change=args.fold_change,
    )

    if args.threads > 1:
        with ProcessPoolExecutor(max_workers=args.threads) as executor:
            results = list(executor.map(worker, gene_rows))
    else:
        results = [worker(r) for r in gene_rows]

    res_df = pd.DataFrame([r for r in results if r is not None])

    if "P_value" in res_df and res_df["P_value"].notna().any():
        valid = res_df["P_value"].notna()
        _, adj_p, _, _ = multipletests(res_df.loc[valid, "P_value"], method="fdr_bh")
        res_df.loc[valid, "Adj_P_value"] = adj_p

    if "Adj_P_value" in res_df.columns:
        res_df["Rhythmic"] = (res_df["Adj_P_value"] < 0.05) & (res_df["Fold_Change"] >= args.fold_change)

    total = len(res_df)
    rhythmic_count = int(res_df["Rhythmic"].sum()) if "Rhythmic" in res_df.columns else 0
    print("\n--- Results Summary ---")
    print("Total genes: " + str(total))
    if total:
        print("Rhythmic (FDR<0.05 + FC>=" + str(args.fold_change) + "): "
              + str(rhythmic_count) + " (" + format(rhythmic_count / total * 100, ".1f") + "%)")

    cols_order = ["Gene_ID", "P_value", "Adj_P_value", "Peak_Time",
                  "Mesor", "Amplitude", "Fold_Change", "Rhythmic"]
    res_df = res_df[[c for c in cols_order if c in res_df.columns]]

    res_df.to_csv(args.output, index=False)
    print("\nSuccess! Saved to " + args.output)
    return res_df


def build_parser(add_help=True):
    parser = argparse.ArgumentParser(
        prog="rgannotator rhythm",
        add_help=add_help,
        description="Module 2 - circadian rhythm detection (JTK_CYCLE or Cosinor)",
    )
    parser.add_argument("-i", "--input", required=True, help="Input expression CSV (genes x time columns)")
    parser.add_argument("-o", "--output", required=True, help="Output result CSV")
    parser.add_argument("-m", "--method", choices=["jtk", "cosinor"], default="jtk",
                        help="Rhythm detection method (default: jtk)")
    parser.add_argument("-CT", "--ct_mode", action="store_true",
                        help="Force CT-column parsing (e.g. CT0, CT15_2). Auto-detected by default.")
    parser.add_argument("-ZT", "--zt_mode", action="store_true",
                        help="Force manual mode: time points built from -zt0/-int/-r")
    parser.add_argument("-t", "--threads", type=int, default=1, help="Parallel processes (default: 1)")
    parser.add_argument("-zt", "--zero_time", type=float, default=0, help="Manual mode: start ZT (default: 0)")
    parser.add_argument("-r", "--replicates", type=int, default=1, help="Manual mode: replicates per point (default: 1)")
    parser.add_argument("-int", "--interval", type=int, default=4, help="Manual mode: hours between points (default: 4)")
    parser.add_argument("-fc", "--fold_change", type=float, default=1.5,
                        help="Minimum fold change (max/min time-point mean). Default: 1.5")
    return parser


def main(argv=None):
    args = build_parser().parse_args(argv)
    run_rhythm(args)


if __name__ == "__main__":
    main()
