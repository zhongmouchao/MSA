#!/usr/bin/env python3
"""
aggregate_counts.py -- Step 4 of the viroGym surveillance pipeline.

Merges the per-chunk motif count tables produced by the alignment jobs into
one training table per benchmark, with STRICT cleaning -- this replaces the
old pipeline's integration, which let garbage date headers spawn columns
like '1953-unknown' and scattered counts across them.

Per group (from task_plan.json), per motif:
  1. find every chunk file  aligned/<group>/**/<chunk>_<tag>_counts_by_<bin>_wide.csv
  2. merge rows by (Motif Nucleotide, Motif Peptide), summing counts
  3. DATE CLEANING -- a column is kept only if:
       - yearly : ^(19|20)NN$           monthly: ^(19|20)NN-(01..12)$
       - year within the plausible range for that virus
         (SCV2 >= 2019, HIV >= 1980, IAV >= 1900, others >= 1950,
          everything <= current year + 1)
     every other column ('1953-unknown', '20213', 'NA', years outside the
     range, ...) is folded into ONE single 'unknown' column
  4. VARIANT FILTER (EvoPU-grade): drop every variant whose Motif Peptide
     contains '-' or 'X' at more than 10% of its positions
     (note: run_pipeline already rejects such SEQUENCES upstream via
     --motif-min-called 0.9; this is the variant-level backstop)
  5. write training_data/<benchmark_stem>_training_counts.csv for EVERY
     benchmark attached to the motif (shared references are counted once,
     written once per benchmark name)

Also writes training_data/_aggregation_summary.csv: per motif, chunks
merged, variants kept/dropped, counts kept/dropped, date range, unknowns.

Run from the pipeline root:

    python3 aggregate_counts.py
    python3 aggregate_counts.py --plan task_plan.json --out-dir training_data
"""

import argparse
import csv
import glob
import os
import re
import sys
import time

YEAR_RE = re.compile(r"^(19|20)\d{2}$")
MONTH_RE = re.compile(r"^(19|20)\d{2}-(0[1-9]|1[0-2])$")

# earliest plausible collection year per virus (1918 flu is real!)
MIN_YEAR = {"SCV2": 2019, "HIV": 1980, "IAV": 1900}
DEFAULT_MIN_YEAR = 1950
MAX_GAP_FRAC = 0.10          # EvoPU filter: reject variants with >10% '-'/'X'
KEY_COLS = ("Motif Nucleotide", "Motif Peptide")


def max_year():
    return time.localtime().tm_year + 1


def col_ok(col, binning, min_year):
    if binning == "monthly":
        m = MONTH_RE.fullmatch(col)
    else:
        m = YEAR_RE.fullmatch(col)
    if not m:
        return False
    year = int(col[:4])
    return min_year <= year <= max_year()


def gap_frac(peptide):
    if not peptide:
        return 1.0
    bad = peptide.count("-") + peptide.count("X")
    return bad / len(peptide)


def find_count_files(root, gid, tag, binning):
    word = "month" if binning == "monthly" else "year"
    pat = os.path.join(root, "aligned", gid, "**",
                       f"*_{tag}_counts_by_{word}_wide.csv")
    return sorted(p for p in glob.glob(pat, recursive=True)
                  if os.path.isfile(p))


def merge_motif(paths, binning, min_year):
    """Merge chunk CSVs -> (merged_dict, kept_cols, unknown_total).

    merged_dict: (nucleotide, peptide) -> {col: count}
    """
    merged = {}
    colset = set()
    unknown_total = 0
    for path in paths:
        with open(path, newline="", encoding="utf-8-sig") as fh:
            reader = csv.DictReader(fh)
            fieldnames = reader.fieldnames or []
            keys = [k for k in KEY_COLS if k in fieldnames]
            if not keys:
                continue
            date_cols = [c for c in fieldnames if c not in KEY_COLS]
            for row in reader:
                nuc = (row.get("Motif Nucleotide") or "").strip()
                pep = (row.get("Motif Peptide") or "").strip()
                if not pep:
                    continue
                rec = merged.setdefault((nuc, pep), {})
                for c in date_cols:
                    v = (row.get(c) or "").strip()
                    if not v:
                        continue
                    try:
                        n = int(float(v))
                    except ValueError:
                        continue
                    if n == 0:
                        continue
                    if col_ok(c, binning, min_year):
                        colset.add(c)
                        rec[c] = rec.get(c, 0) + n
                    else:
                        unknown_total += n
                        rec["unknown"] = rec.get("unknown", 0) + n
    return merged, sorted(colset), unknown_total


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--plan", default="task_plan.json")
    ap.add_argument("--out-dir", default="training_data")
    args = ap.parse_args()

    root = os.getcwd()
    if not os.path.isfile(args.plan):
        print(f"[FATAL] {args.plan} not found -- run plan_tasks.py first")
        return 1
    import json
    with open(args.plan, encoding="utf-8") as fh:
        plan = json.load(fh)

    os.makedirs(args.out_dir, exist_ok=True)
    summary_rows = []
    n_tables = 0

    for g in plan.get("groups", []):
        gid = g["group_id"]
        binning = g.get("binning", "yearly")
        min_year = MIN_YEAR.get(g.get("virus_code", ""), DEFAULT_MIN_YEAR)

        for motif in g.get("motifs", []):
            tag = motif["tag"]
            paths = find_count_files(root, gid, tag, binning)
            if not paths:
                summary_rows.append({
                    "group_id": gid, "motif_tag": tag,
                    "benchmarks": ";".join(motif["benchmarks"]),
                    "chunks_merged": 0, "variants_kept": 0,
                    "variants_dropped_gapX": 0, "counts_kept": 0,
                    "counts_dropped_gapX": 0, "counts_unknown": 0,
                    "date_range": "", "status": "no count files yet"})
                continue

            merged, cols, unknown_total = merge_motif(paths, binning, min_year)

            kept, dropped = {}, {}
            for key, rec in merged.items():
                pep = key[1]
                (dropped if gap_frac(pep) > MAX_GAP_FRAC else kept)[key] = rec

            kept_counts = sum(sum(r.values()) for r in kept.values())
            dropped_counts = sum(sum(r.values()) for r in dropped.values())
            date_range = f"{cols[0]}..{cols[-1]}" if cols else ""

            # one table per benchmark attached to this motif
            for bf in motif["benchmarks"]:
                stem = bf[:-4] if bf.lower().endswith(".csv") else bf
                out_path = os.path.join(args.out_dir,
                                        f"{stem}_training_counts.csv")
                out_cols = [c for c in KEY_COLS] + cols + ["unknown"]
                with open(out_path, "w", newline="", encoding="utf-8") as fh:
                    w = csv.writer(fh)
                    w.writerow(out_cols)
                    for (nuc, pep), rec in sorted(kept.items(),
                                                  key=lambda kv: -sum(kv[1].values())):
                        w.writerow([nuc, pep] +
                                   [rec.get(c, 0) for c in cols] +
                                   [rec.get("unknown", 0)])
                n_tables += 1

            summary_rows.append({
                "group_id": gid, "motif_tag": tag,
                "benchmarks": ";".join(motif["benchmarks"]),
                "chunks_merged": len(paths),
                "variants_kept": len(kept),
                "variants_dropped_gapX": len(dropped),
                "counts_kept": kept_counts,
                "counts_dropped_gapX": dropped_counts,
                "counts_unknown": unknown_total,
                "date_range": date_range, "status": "OK"})

            print(f"  {gid} [{tag}]: {len(paths)} chunk(s), "
                  f"{len(kept):,} variants kept "
                  f"({len(dropped)} dropped by 10% gap/X rule), "
                  f"{kept_counts:,} counts, {unknown_total:,} unknown-date, "
                  f"{len(motif['benchmarks'])} benchmark table(s)")

    summary_path = os.path.join(args.out_dir, "_aggregation_summary.csv")
    with open(summary_path, "w", newline="", encoding="utf-8") as fh:
        fieldnames = ["group_id", "motif_tag", "benchmarks", "chunks_merged",
                      "variants_kept", "variants_dropped_gapX", "counts_kept",
                      "counts_dropped_gapX", "counts_unknown", "date_range",
                      "status"]
        w = csv.DictWriter(fh, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(summary_rows)

    print()
    print("=" * 72)
    print("AGGREGATION SUMMARY")
    print("=" * 72)
    ok = sum(1 for r in summary_rows if r["status"] == "OK")
    todo = sum(1 for r in summary_rows if r["status"] != "OK")
    print(f"motifs aggregated : {ok}")
    print(f"motifs with no data yet : {todo}")
    print(f"training tables written : {n_tables} -> {os.path.abspath(args.out_dir)}/")
    print(f"summary : {summary_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
