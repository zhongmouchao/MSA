#!/usr/bin/env python3
"""
scan_inputs.py -- Step 1 of the viroGym surveillance pipeline.

Recursively scans the raw FASTA folder and records, for every FASTA file,
the facts the task planner needs:

  - molecule type : nucleotide | protein
  - state         : raw | aligned   (aligned = sequences contain gap '-'
                    AND all sampled sequences have the same length)
  - median length : of up to --sample sequences (nt for nucleotide,
                    aa for protein)
  - record count  : exact number of FASTA records (headers counted)

Run from the pipeline root:

    python3 scan_inputs.py                      # scans ./RAW_FASTA_DATA
    python3 scan_inputs.py --fasta-root RAW_FASTA_DATA --sample 1000

Output:
  input_inventory.json   (default; consumed by plan_tasks.py)
  + a printed summary table, one line per FASTA file.

The scanner never modifies anything. Re-run it whenever you add or replace
FASTA files -- it always rebuilds the inventory from scratch.
"""

import argparse
import json
import os
import sys
import time

FASTA_EXTS = (".fasta", ".fa", ".fna", ".faa", ".fas")

NT_CHARS = set("ACGTUNacgtun-.?")
GAP_CHARS = set("-.")


def iter_fasta(path, limit=None):
    """Yield (header, sequence) tuples. Stops after `limit` records if given."""
    header = None
    chunks = []
    n = 0
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks)
                    n += 1
                    if limit and n >= limit:
                        return
                header = line[1:].strip()
                chunks = []
            else:
                chunks.append(line.strip())
        if header is not None:
            yield header, "".join(chunks)


def count_records(path):
    """Exact record count via a fast buffered scan for '>'."""
    n = 0
    with open(path, "rb") as fh:
        while True:
            buf = fh.read(8 * 1024 * 1024)
            if not buf:
                break
            n += buf.count(b"\n>")
    # first record has no leading newline
    with open(path, "rb") as fh:
        first = fh.read(1)
    if first == b">":
        n += 1
    return n


def classify_sample(seqs):
    """Return (molecule, state, median_length) from a list of sequences."""
    if not seqs:
        return "unknown", "unknown", 0

    lengths = sorted(len(s) for s in seqs)
    median = lengths[len(lengths) // 2]

    # molecule: fraction of characters that look nucleotide-like
    total = 0
    nt_like = 0
    for s in seqs:
        total += len(s)
        nt_like += sum(1 for c in s if c in NT_CHARS)
    molecule = "nucleotide" if (total and nt_like / total >= 0.95) else "protein"

    # state: gaps present AND uniform length -> aligned
    has_gaps = any(any(c in GAP_CHARS for c in s) for s in seqs)
    uniform = lengths[0] == lengths[-1]
    state = "aligned" if (has_gaps and uniform) else "raw"

    return molecule, state, median


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fasta-root", default="RAW_FASTA_DATA",
                    help="root folder with raw FASTA files, scanned "
                         "recursively (default: RAW_FASTA_DATA)")
    ap.add_argument("--sample", type=int, default=1000,
                    help="sequences sampled per file for type/length "
                         "detection (default: 1000)")
    ap.add_argument("--out", default="input_inventory.json",
                    help="inventory output path (default: %(default)s)")
    args = ap.parse_args()

    if not os.path.isdir(args.fasta_root):
        print(f"[FATAL] FASTA root folder not found: {args.fasta_root}")
        return 1

    fasta_files = []
    for dirpath, _dirnames, filenames in os.walk(args.fasta_root):
        for fn in sorted(filenames):
            if fn.lower().endswith(FASTA_EXTS):
                fasta_files.append(os.path.join(dirpath, fn))
    fasta_files.sort()

    if not fasta_files:
        print(f"[FATAL] no FASTA files ({', '.join(FASTA_EXTS)}) found under "
              f"{args.fasta_root}")
        return 1

    print(f"Scanning {len(fasta_files)} FASTA file(s) under "
          f"{os.path.abspath(args.fasta_root)} ...")

    entries = []
    for i, path in enumerate(fasta_files, 1):
        rel = os.path.relpath(path, args.fasta_root)
        t0 = time.time()

        seqs = []
        for _h, s in iter_fasta(path, limit=args.sample):
            seqs.append(s)
        molecule, state, median = classify_sample(seqs)
        n_records = count_records(path)
        size_bytes = os.path.getsize(path)

        entries.append({
            "path": os.path.abspath(path),
            "rel_path": rel,
            "folder": os.path.dirname(rel),
            "molecule": molecule,
            "state": state,
            "median_length": median,
            "n_records": n_records,
            "size_bytes": size_bytes,
        })

        print(f"  [{i}/{len(fasta_files)}] {rel}: {molecule}, {state}, "
              f"median {median:,}, {n_records:,} records, "
              f"{size_bytes / 1e6:.1f} MB  ({time.time() - t0:.1f}s)")

    inventory = {
        "fasta_root": os.path.abspath(args.fasta_root),
        "scanned_at": time.strftime("%Y-%m-%d %H:%M:%S"),
        "sample_size": args.sample,
        "files": entries,
    }
    with open(args.out, "w", encoding="utf-8") as fh:
        json.dump(inventory, fh, indent=2)

    print("-" * 72)
    print(f"total files   : {len(entries)}")
    print(f"total records : {sum(e['n_records'] for e in entries):,}")
    print(f"total size    : {sum(e['size_bytes'] for e in entries) / 1e9:.2f} GB")
    print(f"inventory written to {os.path.abspath(args.out)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
