#!/usr/bin/env python3
"""
Split a large FASTA into N parts of (nearly) equal RECORD count.

    python3 split_fasta.py input/big.fasta                 # 4 parts, next to the original
    python3 split_fasta.py input/big.fasta --parts 8
    python3 split_fasta.py input/big.fasta --out-dir input

Writes <base>_part_<i>_of_<N>.fasta files, copying original lines VERBATIM
(no re-wrapping, headers untouched). Splits only at record boundaries —
safe for any FASTA. Streams line by line; memory use is negligible.
"""
import argparse
import math
import os
import sys


def count_records(path):
    n = 0
    with open(path, "rb") as h:
        for line in h:
            if line.startswith(b">"):
                n += 1
    return n


def main():
    ap = argparse.ArgumentParser(description="Split a FASTA into N equal-record parts")
    ap.add_argument("fasta", help="input FASTA file")
    ap.add_argument("--parts", type=int, default=4, help="number of parts (default 4)")
    ap.add_argument("--out-dir", default=None,
                    help="output directory (default: same as the input file)")
    ap.add_argument("--force", action="store_true", help="overwrite existing part files")
    args = ap.parse_args()

    if not os.path.exists(args.fasta):
        sys.exit(f"not found: {args.fasta}")
    if args.parts < 2:
        sys.exit("--parts must be >= 2")

    out_dir = args.out_dir or os.path.dirname(os.path.abspath(args.fasta))
    os.makedirs(out_dir, exist_ok=True)
    base = os.path.splitext(os.path.basename(args.fasta))[0]
    ext = os.path.splitext(args.fasta)[1] or ".fasta"

    n = count_records(args.fasta)
    if n == 0:
        sys.exit(f"no FASTA records found in {args.fasta}")
    per = math.ceil(n / args.parts)
    print(f"{args.fasta}: {n:,} records -> {args.parts} parts of ~{per:,} records")

    paths = [os.path.join(out_dir, f"{base}_part_{i}_of_{args.parts}{ext}")
             for i in range(1, args.parts + 1)]
    if not args.force:
        existing = [p for p in paths if os.path.exists(p)]
        if existing:
            sys.exit(f"part file(s) already exist (use --force to overwrite):\n  "
                     + "\n  ".join(existing))

    handles = [open(p, "wb") for p in paths]
    counts = [0] * args.parts
    rec = 0
    current = None
    try:
        with open(args.fasta, "rb") as h:
            for line in h:
                if line.startswith(b">"):
                    rec += 1
                    current = min((rec - 1) // per, args.parts - 1)
                    counts[current] += 1
                if current is not None:
                    handles[current].write(line)
    finally:
        for h in handles:
            h.close()

    for p, c in zip(paths, counts):
        print(f"  wrote {c:,} records -> {p}")
    assert sum(counts) == n, "record count mismatch!"


if __name__ == "__main__":
    main()
