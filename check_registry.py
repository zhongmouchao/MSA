#!/usr/bin/env python3
"""
check_registry.py -- Step 0 of the viroGym surveillance pipeline.

Cross-checks the DMS benchmark folder against the registry CSV before any
alignment work is planned. The registry CSV lives INSIDE the benchmark folder
together with the DMS data files. Run it from the pipeline root (the folder
that contains Virogym_Benchmark/):

    python3 check_registry.py
    python3 check_registry.py --benchmark-dir Virogym_Benchmark
    python3 check_registry.py --registry path/to/other_registry.csv

It reports four categories:

  [OK]        DMS file present AND registered          -> will be planned
  [NEW]       DMS file present but NOT in the registry -> skipped until you
              add a row to the registry CSV
  [MISSING]   registry row with no DMS file            -> skipped (dataset not
              provided yet; harmless)
  [PROBLEM]   registry row fails a content check       -> must be fixed before
              the planner will accept this benchmark

Content checks on each registry row:
  - required columns present and non-empty:
      benchmark_file, dms_id, virus_code, host, target_seq,
      target_type, recommended_reference
  - target_seq contains only amino-acid letters plus '-', '.', '*', 'X'
      (gaps are allowed here; the planner de-gaps before use)
  - seq_len (if given) equals the length of target_seq as stored,
      gap characters included                                          [warning]
  - recommended_reference is 'target_as_reference' or 'full_protein_plus_motif'
  - no duplicate benchmark_file rows in the registry

Exit code: 0 when there are no [PROBLEM] rows, 1 otherwise.
[NEW] and [MISSING] are informational and never block.
"""

import argparse
import csv
import os
import sys

REQUIRED_COLUMNS = [
    "benchmark_file",
    "dms_id",
    "virus_code",
    "host",
    "target_seq",
    "target_type",
    "recommended_reference",
]

ALLOWED_REFERENCE_MODES = {"target_as_reference", "full_protein_plus_motif"}
VALID_SEQ_CHARS = set("ACDEFGHIKLMNPQRSTVWYBXZUOJ-. *")  # generous on purpose


def degap(seq: str) -> str:
    return seq.replace("-", "").replace(".", "").replace(" ", "").strip()


def load_registry(path: str):
    with open(path, newline="", encoding="utf-8-sig") as fh:
        reader = csv.DictReader(fh)
        rows = list(reader)
        cols = reader.fieldnames or []
    return cols, rows


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--benchmark-dir", default="Virogym_Benchmark",
                    help="folder containing the raw DMS csv files and the "
                         "registry CSV (default: Virogym_Benchmark)")
    ap.add_argument("--registry", default=None,
                    help="registry CSV path (default: "
                         "<benchmark-dir>/virogym_79_benchmarks_combined_with_target_type.csv)")
    args = ap.parse_args()

    if args.registry is None:
        args.registry = os.path.join(
            args.benchmark_dir,
            "virogym_79_benchmarks_combined_with_target_type.csv")

    # ---- locate inputs ----------------------------------------------------
    if not os.path.isdir(args.benchmark_dir):
        print(f"[FATAL] benchmark folder not found: {args.benchmark_dir}")
        return 1
    if not os.path.isfile(args.registry):
        print(f"[FATAL] registry CSV not found: {args.registry}")
        return 1

    registry_name = os.path.basename(args.registry)
    dms_files = sorted(f for f in os.listdir(args.benchmark_dir)
                       if f.lower().endswith(".csv") and f != registry_name)
    cols, rows = load_registry(args.registry)

    # ---- column sanity ----------------------------------------------------
    missing_cols = [c for c in REQUIRED_COLUMNS if c not in cols]
    if missing_cols:
        print(f"[FATAL] registry is missing required column(s): {missing_cols}")
        print(f"        columns found: {cols}")
        return 1

    by_file = {}
    problems = []          # (benchmark_file, message)
    dup_seen = set()
    for r in rows:
        bf = (r.get("benchmark_file") or "").strip()
        if not bf:
            problems.append(("(blank benchmark_file)", "row has no benchmark_file"))
            continue
        if bf in dup_seen:
            problems.append((bf, "duplicate benchmark_file row in registry"))
        dup_seen.add(bf)
        by_file[bf] = r

    # ---- per-row content checks -------------------------------------------
    for bf, r in by_file.items():
        for col in REQUIRED_COLUMNS:
            if not (r.get(col) or "").strip():
                problems.append((bf, f"empty required field: {col}"))

        raw_seq = (r.get("target_seq") or "").strip()
        if raw_seq:
            bad = set(raw_seq.upper()) - VALID_SEQ_CHARS
            if bad:
                problems.append((bf, f"target_seq has unexpected character(s): {sorted(bad)}"))
            dg = degap(raw_seq.upper())
            seq_len_field = (r.get("seq_len") or "").strip()
            stored_len = len(raw_seq.replace(" ", ""))
            if seq_len_field.isdigit() and int(seq_len_field) not in (stored_len, len(dg)):
                problems.append((
                    bf,
                    f"seq_len={seq_len_field} but target_seq is {stored_len} aa as stored "
                    f"({len(dg)} aa de-gapped) -- update seq_len or the sequence"))
            if "-" in raw_seq or "." in raw_seq:
                print(f"  [note] {bf}: target_seq contains gaps -- "
                      f"planner will de-gap before use")

        mode = (r.get("recommended_reference") or "").strip()
        if mode and mode not in ALLOWED_REFERENCE_MODES:
            problems.append((bf,
                             f"recommended_reference must be one of "
                             f"{sorted(ALLOWED_REFERENCE_MODES)}, got '{mode}'"))

    # ---- folder <-> registry cross-check ----------------------------------
    registered = set(by_file)
    on_disk = set(dms_files)

    ok = sorted(on_disk & registered)
    new = sorted(on_disk - registered)
    missing = sorted(registered - on_disk)

    # ---- report -------------------------------------------------------------
    print()
    print("=" * 72)
    print("REGISTRY CHECK SUMMARY")
    print("=" * 72)
    print(f"benchmark folder : {os.path.abspath(args.benchmark_dir)}")
    print(f"registry         : {os.path.abspath(args.registry)}")
    print(f"DMS files on disk: {len(dms_files)}")
    print(f"registry rows    : {len(rows)}")
    print("-" * 72)
    print(f"[OK]      registered and present : {len(ok)}")
    print(f"[NEW]     on disk, not registered: {len(new)}")
    print(f"[MISSING] registered, no file    : {len(missing)}")
    print(f"[PROBLEM] rows failing checks    : {len(problems)}")
    print("-" * 72)

    if new:
        print("\n[NEW] files skipped until you add a registry row:")
        for f in new:
            print(f"  - {f}")
    if missing:
        print("\n[MISSING] registry rows with no DMS file (skipped):")
        for f in missing:
            print(f"  - {f}  (dms_id={by_file[f].get('dms_id','')})")
    if problems:
        print("\n[PROBLEM] rows that must be fixed:")
        for bf, msg in problems:
            print(f"  - {bf}: {msg}")

    if not problems:
        print("\nAll registered rows pass content checks.")
        if new:
            print("(Unregistered [NEW] files exist but never block the pipeline.)")
        return 0
    print("\nFix the [PROBLEM] rows above, then re-run this check.")
    return 1


if __name__ == "__main__":
    sys.exit(main())
