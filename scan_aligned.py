#!/usr/bin/env python3
"""
scan_aligned.py -- Step 1.5 of the viroGym surveillance pipeline.

Inventories the aligned/ folder so the pipeline never redoes finished work,
and registers alignment runs that were produced outside this pipeline.

Two modes
---------

1) Scan mode (default):

    python3 scan_aligned.py

   Reads every group folder under aligned/. A group folder is usable when it
   contains alignment_manifest.json (written by run_tasks.py or by
   --register). For each chunk subfolder the scanner records:
     - .msa_complete marker present?          -> alignment reusable
     - <stem>_<tag>_counts_by_*_wide.csv      -> counting already done
   Writes aligned_inventory.json.

2) Register mode (for your EXISTING runs, e.g. aligned/sars_spike):

    python3 scan_aligned.py --register aligned/sars_spike \
        --config path/to/the/reference.txt/you/used

   - parses the config, hashes its >protein_reference block
   - verifies every chunk subfolder (.msa_complete + aligned proteins fasta)
   - writes alignment_manifest.json inside that folder
   - tells you the canonical group name (e.g. SCV2_1273aa_a1b2c3d4) and how
     to rename/symlink the folder so the pipeline finds it
   - if the registry is available, also reports whether any motif in the
     config is IDENTICAL to a benchmark target_seq -- those benchmarks can
     reuse the already-computed counts CSVs directly.

Reuse rule (locked design): an alignment is reused only when the reference
hash matches exactly. No coordinate mapping, no exceptions.
"""

import argparse
import glob
import hashlib
import json
import os
import re
import sys
import time

REGISTRY_NAME = "virogym_79_benchmarks_combined_with_target_type.csv"
COUNTS_RE = re.compile(r"^(?P<stem>.+)_(?P<tag>.+)_counts_by_(month|year)_wide\.csv$")


def degap(seq: str) -> str:
    return seq.replace("-", "").replace(".", "").replace(" ", "").strip().upper()


def ref_hash(seq: str) -> str:
    return hashlib.sha1(seq.encode()).hexdigest()[:8]


def parse_config(path):
    """Parse a run_pipeline config: returns (protein_ref, {tag: motif_seq})."""
    blocks = {}
    name = None
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                name = line[1:].strip()
                blocks[name] = []
            elif name and line:
                blocks[name].append(line)
    seqs = {k: degap("".join(v)) for k, v in blocks.items()}
    protein = seqs.get("protein_reference")
    motifs = {k[len("motif_reference_"):]: v
              for k, v in seqs.items() if k.startswith("motif_reference_")}
    return protein, motifs


def scan_chunk_folder(path):
    """Inspect one chunk folder; returns a status dict."""
    stem = os.path.basename(path)
    msa_done = os.path.isfile(os.path.join(path, ".msa_complete"))
    has_aligned_proteins = bool(glob.glob(
        os.path.join(path, "*_aligned_proteins.fasta"))) or bool(
        glob.glob(os.path.join(path, "*_aligned_protein.fasta")))
    counts = {}
    for csv_path in glob.glob(os.path.join(path, "*_counts_by_*_wide.csv")):
        m = COUNTS_RE.match(os.path.basename(csv_path))
        if m:
            counts[m.group("tag")] = csv_path
    return {
        "msa_complete": msa_done,
        "aligned_proteins": has_aligned_proteins,
        "counts_by_motif": counts,
        "reusable": msa_done and has_aligned_proteins,
    }


def load_manifest(group_dir):
    p = os.path.join(group_dir, "alignment_manifest.json")
    if not os.path.isfile(p):
        return None
    try:
        with open(p, encoding="utf-8") as fh:
            return json.load(fh)
    except (json.JSONDecodeError, OSError):
        return None


def cmd_scan(aligned_root, out_path):
    groups = {}
    unregistered = []
    if not os.path.isdir(aligned_root):
        print(f"[note] no aligned/ folder found at {aligned_root} -- "
              f"nothing aligned yet, inventory will be empty")
    else:
        for name in sorted(os.listdir(aligned_root)):
            gdir = os.path.join(aligned_root, name)
            if not os.path.isdir(gdir):
                continue
            manifest = load_manifest(gdir)
            if manifest is None:
                unregistered.append(name)
                continue
            chunks = {}
            for sub in sorted(os.listdir(gdir)):
                cdir = os.path.join(gdir, sub)
                if os.path.isdir(cdir):
                    chunks[sub] = scan_chunk_folder(cdir)
            groups[manifest["group_id"]] = {
                "folder": os.path.abspath(gdir),
                "ref_hash": manifest["ref_hash"],
                "ref_len": manifest["ref_len"],
                "config": manifest.get("config"),
                "chunks": chunks,
            }

    inv = {
        "aligned_root": os.path.abspath(aligned_root),
        "scanned_at": time.strftime("%Y-%m-%d %H:%M:%S"),
        "groups": groups,
        "unregistered_folders": unregistered,
    }
    with open(out_path, "w", encoding="utf-8") as fh:
        json.dump(inv, fh, indent=2)

    print()
    print("=" * 72)
    print("ALIGNED-DATA INVENTORY")
    print("=" * 72)
    if not groups and not unregistered:
        print("no aligned data found -- everything will be aligned from scratch")
    for gid, g in groups.items():
        done = sum(1 for c in g["chunks"].values() if c["reusable"])
        counted = sum(1 for c in g["chunks"].values() if c["counts_by_motif"])
        print(f"  {gid}: {done} reusable chunk(s), "
              f"{counted} with counts already computed")
    for name in unregistered:
        print(f"  [unregistered] {name}  -> run with "
              f"--register aligned/{name} --config <the config you used>")
    print(f"\ninventory written to {os.path.abspath(out_path)}")
    return 0


def cmd_register(folder, config, aligned_root, benchmark_dir, registry_path):
    if not os.path.isdir(folder):
        print(f"[FATAL] folder not found: {folder}")
        return 1
    if not os.path.isfile(config):
        print(f"[FATAL] config not found: {config}")
        return 1

    protein, motifs = parse_config(config)
    if not protein:
        print(f"[FATAL] no >protein_reference block found in {config}")
        return 1

    h = ref_hash(protein)
    # virus code guess: look at registry for a matching target
    virus = "UNK"
    motif_matches = []
    if registry_path and os.path.isfile(registry_path):
        import csv as _csv
        with open(registry_path, newline="", encoding="utf-8-sig") as fh:
            rows = list(_csv.DictReader(fh))
        for r in rows:
            t = degap(r.get("target_seq") or "")
            if not t:
                continue
            if t == protein:
                virus = (r.get("virus_code") or "UNK").strip() or "UNK"
            for tag, mseq in motifs.items():
                if t == mseq:
                    motif_matches.append((tag, r["benchmark_file"]))

    group_id = f"{virus}_{len(protein)}aa_{h}"

    # verify chunks
    chunks = {}
    n_reusable = 0
    for sub in sorted(os.listdir(folder)):
        cdir = os.path.join(folder, sub)
        if os.path.isdir(cdir):
            st = scan_chunk_folder(cdir)
            chunks[sub] = st
            n_reusable += st["reusable"]

    manifest = {
        "group_id": group_id,
        "ref_hash": h,
        "ref_len": len(protein),
        "config": os.path.abspath(config),
        "motif_tags": sorted(motifs),
        "registered_at": time.strftime("%Y-%m-%d %H:%M:%S"),
        "source_folder_name": os.path.basename(os.path.normpath(folder)),
    }
    with open(os.path.join(folder, "alignment_manifest.json"), "w",
              encoding="utf-8") as fh:
        json.dump(manifest, fh, indent=2)

    print()
    print("=" * 72)
    print("REGISTER EXISTING ALIGNMENT")
    print("=" * 72)
    print(f"folder           : {os.path.abspath(folder)}")
    print(f"config           : {os.path.abspath(config)}")
    print(f"protein reference: {len(protein)} aa, hash {h}")
    print(f"motifs in config : {', '.join(sorted(motifs)) or '(none)'}")
    print(f"chunk folders    : {len(chunks)}  ({n_reusable} reusable: "
          f".msa_complete + aligned proteins)")
    print("-" * 72)
    print(f"canonical group name: {group_id}")
    canonical = os.path.join(aligned_root, group_id)
    if os.path.abspath(folder) != os.path.abspath(canonical):
        print(f"\nto make the pipeline see it, either rename:")
        print(f"  mv {folder} {canonical}")
        print(f"or symlink (keeps your current name working):")
        print(f"  ln -s {os.path.abspath(folder)} {canonical}")
    else:
        print("folder already has the canonical name -- nothing to move")

    bad = {s: c for s, c in chunks.items() if not c["reusable"]}
    if bad:
        print(f"\nwarning: {len(bad)} chunk folder(s) NOT reusable "
              f"(missing .msa_complete or aligned proteins):")
        for s in bad:
            print(f"  - {s}")

    if motif_matches:
        print("\nmotifs in this config that EXACTLY match benchmark targets")
        print("(their already-computed counts CSVs can be reused):")
        for tag, bf in sorted(motif_matches):
            print(f"  - motif '{tag}' == target of {bf}")
    else:
        print("\nno motif in this config exactly matches any benchmark target;")
        print("counting will be redone by the pipeline where needed.")
    print(f"\nmanifest written inside {folder}")
    print("re-run:  python3 scan_aligned.py   to include it in the inventory")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--aligned-root", default="aligned")
    ap.add_argument("--out", default="aligned_inventory.json")
    ap.add_argument("--register", metavar="FOLDER", default=None,
                    help="register an existing aligned folder instead of scanning")
    ap.add_argument("--config", default=None,
                    help="the run_pipeline config used for that folder "
                         "(required with --register)")
    ap.add_argument("--benchmark-dir", default="Virogym_Benchmark")
    ap.add_argument("--registry", default=None)
    args = ap.parse_args()

    if args.register:
        if not args.config:
            print("[FATAL] --register needs --config <the config you used>")
            return 1
        registry = args.registry or os.path.join(args.benchmark_dir, REGISTRY_NAME)
        return cmd_register(args.register, args.config,
                            args.aligned_root, args.benchmark_dir, registry)
    return cmd_scan(args.aligned_root, args.out)


if __name__ == "__main__":
    sys.exit(main())
