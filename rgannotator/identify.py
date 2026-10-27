"""Module 1 - identify circadian clock genes in a genome proteome.

Reciprocal-best-hit homology search against a curated panel of 49 Arabidopsis
clock genes:

1. makeblastdb on the target proteome
2. forward BLASTP: clock panel -> target proteome (candidate collection)
3. reverse BLASTP: candidates -> full TAIR10 proteome (specificity check)
4. keep candidates whose best TAIR10 hit is itself a clock gene
5. rename hits ``<prefix><ClockName>-<rank>`` (e.g. ``MsCCA1-1``) and write
   a per-species table + renamed FASTA, plus an all-species summary.
"""

import argparse
import glob
import os
import shutil
import subprocess
import tarfile

import pandas as pd
from Bio import SeqIO

try:
    from . import __version__
except ImportError:  # allow running as a plain script
    __version__ = "1.0.0"

# Bundled data shipped with the repository
PKG_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(os.path.dirname(PKG_DIR), "data")
DEFAULT_CIRCADIAN_FA = os.path.join(DATA_DIR, "circadian.fa")
DEFAULT_DB_PREFIX = os.path.join(DATA_DIR, "TAIR10_DB", "TAIR10_DB_new")
DB_TARBALL = os.path.join(DATA_DIR, "TAIR10_DB_new.tar.gz")

PROTEIN_EXTS = (".pep.fa", ".pep.fasta", ".protein.fa", ".fa", ".fasta", ".faa")


def run_cmd(cmd):
    subprocess.run(cmd, check=True)


def check_blast_tools():
    missing = [t for t in ("makeblastdb", "blastp") if shutil.which(t) is None]
    if missing:
        raise SystemExit(
            "ERROR: BLAST+ executables not found on PATH: " + ", ".join(missing)
            + "\nInstall NCBI BLAST+ (e.g. conda install -c bioconda blast)."
        )


def ensure_tair_db(db_prefix):
    """Extract the bundled TAIR10 BLAST database on first use."""
    if os.path.exists(db_prefix + ".pin") or os.path.exists(db_prefix + ".phr"):
        return db_prefix
    if not os.path.exists(DB_TARBALL):
        raise SystemExit(
            "ERROR: TAIR10 BLAST database not found.\nExpected either database files at "
            + db_prefix + ".* or the archive " + DB_TARBALL + ".\n"
            "Alternatively rebuild it from the TAIR10 proteome with:\n"
            "  makeblastdb -in TAIR10_pep_20210128 -dbtype prot -out " + db_prefix
        )
    dest_dir = os.path.dirname(db_prefix)
    os.makedirs(dest_dir, exist_ok=True)
    print("Extracting bundled TAIR10 BLAST database -> " + dest_dir)
    with tarfile.open(DB_TARBALL, "r:gz") as tar:
        tar.extractall(dest_dir)
    return db_prefix


def parse_circadian_info(fasta_path):
    """Map TAIR accession -> clock gene symbol from '>' header 'Symbol | ATxGxxxxx.n | Category'."""
    mapping = {}
    for record in SeqIO.parse(fasta_path, "fasta"):
        if "|" in record.description:
            symbol = record.description.split("|")[0].strip()
            ath_id = record.description.split("|")[1].strip()
            mapping[ath_id] = symbol
        else:
            mapping[record.id] = record.id
    return mapping


def get_prefix(filename):
    basename = os.path.basename(filename)
    name = basename
    for ext in PROTEIN_EXTS:
        if name.endswith(ext):
            name = name[: -len(ext)]
            break
    parts = name.split("_")
    if len(parts) >= 2:
        prefix = parts[0][0].upper() + parts[1][0].lower()
    else:
        prefix = parts[0][:2].capitalize()
    return prefix


def species_stem(filename):
    basename = os.path.basename(filename)
    for ext in PROTEIN_EXTS:
        if basename.endswith(ext):
            return basename[: -len(ext)]
    return basename


def collect_pep_files(input_path):
    """Accept a single proteome file or a directory of proteomes."""
    if os.path.isdir(input_path):
        files = []
        for ext in PROTEIN_EXTS:
            files.extend(glob.glob(os.path.join(input_path, "*" + ext)))
        return sorted(set(files))
    return [input_path]


def process_single_file(input_file, db_prefix, circadian_fa, circadian_map, circadian_ids,
                        threads, evalue, output_dir):
    species_name = species_stem(input_file)
    prefix = get_prefix(input_file)
    print("")
    print("=" * 60)
    print("  Processing: " + species_name + "  (Prefix: " + prefix + ")")
    print("=" * 60)

    tmp_dir = os.path.join(output_dir, "tmp_blast_" + species_name)
    if os.path.exists(tmp_dir):
        shutil.rmtree(tmp_dir)
    os.makedirs(tmp_dir)

    try:
        # Step 1: build target database
        print("  [1/6] Building database: " + input_file)
        mc_db = os.path.join(tmp_dir, "mc_db")
        run_cmd(["makeblastdb", "-in", input_file, "-dbtype", "prot", "-out", mc_db])

        # Step 2: forward BLAST (collect candidates)
        print("  [2/6] Forward BLAST (finding candidates)")
        fwd_out = os.path.join(tmp_dir, "fwd_results.txt")
        run_cmd([
            "blastp", "-query", circadian_fa, "-db", mc_db,
            "-evalue", str(evalue), "-num_threads", str(threads),
            "-outfmt", "6 qseqid sseqid", "-out", fwd_out,
        ])

        fwd_hits = pd.read_csv(fwd_out, sep="\t", names=["Ath_Query", "Mc_Candidate"])
        candidate_ids = fwd_hits["Mc_Candidate"].unique()

        if len(candidate_ids) == 0:
            print("  X " + species_name + ": No candidates found, skipping.")
            return None

        # Step 3: extract candidate sequences
        print("  [3/6] Extracting " + str(len(candidate_ids)) + " candidate sequences")
        candidate_fa = os.path.join(tmp_dir, "candidates.fa")
        mc_dict = SeqIO.to_dict(SeqIO.parse(input_file, "fasta"))
        with open(candidate_fa, "w") as f:
            for cid in candidate_ids:
                if cid in mc_dict:
                    SeqIO.write(mc_dict[cid], f, "fasta")

        # Step 4: reverse BLAST against full TAIR10 proteome (validation)
        print("  [4/6] Reverse BLAST (validation)")
        rev_out = os.path.join(tmp_dir, "rev_results.txt")
        run_cmd([
            "blastp", "-query", candidate_fa, "-db", db_prefix,
            "-evalue", str(evalue), "-max_target_seqs", "500",
            "-num_threads", str(threads),
            "-outfmt", "6 qseqid sseqid evalue bitscore", "-out", rev_out,
        ])

        # Step 5: final filtering and renaming
        print("  [5/6] Final filtering and renaming")
        rev_df = pd.read_csv(rev_out, sep="\t", names=["Mc_Gene", "Ath_Hit", "Evalue", "Bitscore"])
        rev_df = rev_df.sort_values(by=["Mc_Gene", "Bitscore"], ascending=[True, False])
        best_hits_df = rev_df.drop_duplicates(subset=["Mc_Gene"], keep="first").copy()

        def validate_and_map(hit_id):
            if hit_id in circadian_ids:
                return hit_id
            base = hit_id.split(".")[0]
            for cid in circadian_ids:
                if base in cid:
                    return cid
            return None

        best_hits_df["Matched_Ath_ID"] = best_hits_df["Ath_Hit"].apply(validate_and_map)
        final_df = best_hits_df.dropna(subset=["Matched_Ath_ID"]).copy()
        final_df["Clock_Name"] = final_df["Matched_Ath_ID"].map(circadian_map)

        if len(final_df) == 0:
            print("  X " + species_name + ": No validated circadian homologs found.")
            return None

        final_df = final_df.sort_values(by=["Clock_Name", "Evalue", "Bitscore"], ascending=[True, True, False])
        final_df["Rank"] = final_df.groupby("Clock_Name").cumcount() + 1
        final_df["New_Name"] = final_df.apply(lambda x: prefix + x["Clock_Name"] + "-" + str(x["Rank"]), axis=1)
        final_df["Species"] = species_name

        # Step 6: output
        print("  [6/6] Writing output")
        csv_path = os.path.join(output_dir, species_name + "_clockfinder.csv")
        fa_path = os.path.join(output_dir, species_name + "_clockfinder.fa")

        cols = ["Species", "Mc_Gene", "Matched_Ath_ID", "Clock_Name", "New_Name", "Evalue", "Bitscore"]
        final_df[cols].to_csv(csv_path, index=False)

        with open(fa_path, "w") as f_fasta:
            for _, row in final_df.iterrows():
                gene_id = row["Mc_Gene"]
                if gene_id in mc_dict:
                    record = mc_dict[gene_id]
                    record.id = row["New_Name"]
                    record.description = ""
                    SeqIO.write(record, f_fasta, "fasta")

        print("  OK " + species_name + ": " + str(len(final_df)) + " circadian homologs identified")
        print("     CSV:   " + csv_path)
        print("     FASTA: " + fa_path)

        return final_df[cols]

    finally:
        if os.path.exists(tmp_dir):
            shutil.rmtree(tmp_dir)


def run_identify(args):
    check_blast_tools()

    pep_files = collect_pep_files(args.input)
    if len(pep_files) == 0:
        raise SystemExit("ERROR: no protein FASTA found under " + args.input
                         + " (expected extensions: " + ", ".join(PROTEIN_EXTS) + ")")

    db_prefix = ensure_tair_db(args.tair_db)
    circadian_fa = args.circadian_fa
    if not os.path.exists(circadian_fa):
        raise SystemExit("ERROR: circadian gene FASTA not found: " + circadian_fa)

    print("Found " + str(len(pep_files)) + " protein file(s) to process:")
    for f in pep_files:
        print("  " + species_stem(f).ljust(40) + " -> Prefix: " + get_prefix(f))

    os.makedirs(args.output_dir, exist_ok=True)

    print("\n--- Loading circadian gene info ---")
    circadian_map = parse_circadian_info(circadian_fa)
    circadian_ids = set(circadian_map.keys())
    print("Total " + str(len(circadian_ids)) + " circadian query sequences")
    print("Reverse BLAST database: " + db_prefix)

    all_results = []
    success_count = 0
    fail_count = 0

    for i, pep_file in enumerate(pep_files, 1):
        print("\n[" + str(i) + "/" + str(len(pep_files)) + "] Starting...")
        try:
            result = process_single_file(
                pep_file, db_prefix, circadian_fa, circadian_map, circadian_ids,
                args.threads, args.evalue, args.output_dir,
            )
        except subprocess.CalledProcessError as e:
            print("  X BLAST failed for " + pep_file + ": " + str(e))
            result = None
        if result is not None and len(result) > 0:
            all_results.append(result)
            success_count += 1
        else:
            fail_count += 1

    if all_results:
        print("\n" + "=" * 60)
        print("  Combining all species results")
        print("=" * 60)
        combined_df = pd.concat(all_results, ignore_index=True)
        combined_csv = os.path.join(args.output_dir, "all_species_clockfinder_summary.csv")
        combined_df.to_csv(combined_csv, index=False)
        print("OK Summary: " + combined_csv)
        print("   Total " + str(len(combined_df)) + " records from " + str(success_count) + " species")

    print("\n" + "=" * 60)
    print("  Done!")
    print("  Success: " + str(success_count) + " species")
    print("  Failed:  " + str(fail_count) + " species")
    print("  Output:  " + args.output_dir + "/")
    print("=" * 60)


def build_parser(add_help=True):
    parser = argparse.ArgumentParser(
        prog="rgannotator identify",
        add_help=add_help,
        description="Module 1 - identify circadian clock genes via reciprocal BLASTP",
    )
    parser.add_argument("-i", "--input", required=True,
                        help="Target proteome FASTA, or a directory containing proteome FASTAs")
    parser.add_argument("-o", "--output_dir", default="clockfinder_results",
                        help="Output directory (default: clockfinder_results)")
    parser.add_argument("-t", "--threads", type=int, default=8, help="BLAST threads (default: 8)")
    parser.add_argument("-e", "--evalue", type=float, default=1e-10, help="E-value threshold (default: 1e-10)")
    parser.add_argument("--circadian_fa", default=DEFAULT_CIRCADIAN_FA,
                        help="Curated clock gene FASTA (default: bundled data/circadian.fa)")
    parser.add_argument("--tair_db", default=DEFAULT_DB_PREFIX,
                        help="TAIR10 BLASTP database prefix (default: bundled data/TAIR10_DB/TAIR10_DB_new, "
                             "auto-extracted from data/TAIR10_DB_new.tar.gz on first use)")
    return parser


def main(argv=None):
    args = build_parser().parse_args(argv)
    run_identify(args)


if __name__ == "__main__":
    main()
