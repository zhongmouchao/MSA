#!/usr/bin/env python3
"""
plan_tasks.py -- Step 2 of the viroGym surveillance pipeline.

Joins the registry CSV (what each benchmark needs) with input_inventory.json
(what FASTA data exists) and writes task_plan.json plus one ready-to-use
pipeline config file per alignment group.

Run from the pipeline root, AFTER check_registry.py and scan_inputs.py:

    python3 plan_tasks.py
    python3 plan_tasks.py --benchmark-dir Virogym_Benchmark \
        --inventory input_inventory.json --out task_plan.json

Grouping rules
--------------
- Benchmarks are grouped by their DE-GAPPED protein reference sequence.
  One alignment per (group, FASTA file); every distinct motif in the group
  is then extracted from that same alignment and counted separately.
- Benchmarks sharing BOTH reference and motif (e.g. replicate studies of the
  same protein) are counted once; the aggregator writes the result out under
  each benchmark's own name.
- recommended_reference = 'target_as_reference'
      protein_reference = motif = de-gapped target_seq (DMS coordinates)
- recommended_reference = 'full_protein_plus_motif'
      protein_reference = full-length UniProt protein, motif = target_seq.
      The full protein must be supplied as references/<stem>_fullprotein.fasta
      or as a 'full_protein_seq' column in the registry, else the benchmark
      is BLOCKED with instructions.

FASTA routing rules
-------------------
- SCV2            -> folder SCV2/
- IAV, human host -> IAV/human/<segment>, then IAV/human/, then IAV/all/<segment>,
                     then IAV/all/   (most specific host wins; never mixes
                     a host folder with all/ for the same benchmark)
- IAV, avian host -> same with avian first
- IAV, other host -> IAV/all/
- any other virus -> folder named by virus_code
- When an IAV file is NOT inside a dedicated segment subfolder, a segment
  hint is read from its filename (H1/H3/H5/HA, N1/NA, PB1, PB2, PA, NP);
  a file whose hint contradicts the benchmark's segment is excluded.

Suitability rules (per FASTA file vs reference of length L aa)
--------------------------------------------------------------
- gene match      : 0.6 <= median_len/3 / L <= 1.5   (nucleotide)
                    0.6 <= median_len   / L <= 1.5   (protein)
- parent match    : same ratio against the UniProt full-length parent,
                    only when the target is a sub-region of that parent
                    (lets spike-length files serve RBD benchmarks, while
                    still rejecting them for e.g. Mpro)
- genome match    : median nucleotide length >= 8000 nt
- aligned inputs  : recorded as mode 'aligned' (the runner skips the MSA step)

Output: task_plan.json + planned_configs/<group_id>.txt for every group
that has at least one suitable FASTA file.
"""

import argparse
import csv
import hashlib
import json
import os
import re
import sys
import time

# ---- tunables ---------------------------------------------------------------
GENOME_MIN_NT = 8000        # median nt length at/above which a file counts as genome
GENE_RATIO_MIN = 0.6
GENE_RATIO_MAX = 1.5
PARENT_MIN_FACTOR = 1.2     # parent must be >= 1.2x target to count as a real parent

# IAV segment detection from the registry 'protein' field
SEGMENT_KEYWORDS = [
    ("PB2", ["PB2"]),
    ("PB1", ["PB1"]),
    ("PA",  ["PA polymerase", "PA subunit"]),
    ("HA",  ["hemagglutinin", "HA "]),
    ("NA",  ["neuraminidase", "NA "]),
    ("NP",  ["nucleoprotein", "NP "]),
]

# IAV segment hints from FASTA filenames
HINT_PATTERNS = [
    ("PB2", re.compile(r"PB2", re.I)),
    ("PB1", re.compile(r"PB1", re.I)),
    ("PA",  re.compile(r"(?:^|[_\-])PA(?:[_\-]|$)", re.I)),
    ("NP",  re.compile(r"(?:^|[_\-])NP(?:[_\-]|$)|NUCLEOP", re.I)),
    ("HA",  re.compile(r"(?:^|[_\-])H\d+|HEMAG|(?:^|[_\-])HA(?:[_\-]|$)", re.I)),
    ("NA",  re.compile(r"(?:^|[_\-])N\d+|NEURAM|(?:^|[_\-])NA(?:[_\-]|$)", re.I)),
]

REGISTRY_NAME = "virogym_79_benchmarks_combined_with_target_type.csv"


def degap(seq: str) -> str:
    return seq.replace("-", "").replace(".", "").replace(" ", "").strip().upper()


def ref_hash(seq: str) -> str:
    return hashlib.sha1(seq.encode()).hexdigest()[:8]


def detect_segment(protein_field: str):
    text = protein_field or ""
    low = text.lower()
    for seg, keys in SEGMENT_KEYWORDS:
        for k in keys:
            if k.lower() in low:
                return seg
    return None


def filename_segment_hint(fname: str):
    base = os.path.basename(fname)
    for seg, pat in HINT_PATTERNS:
        if pat.search(base):
            return seg
    return None


def folder_tiers(virus_code: str, host: str, segment):
    """Ordered list of folder prefixes; first non-empty tier wins."""
    v = (virus_code or "").strip()
    h = (host or "").strip().lower()
    if v == "SCV2":
        return ["SCV2"]
    if v == "IAV":
        tiers = []
        if h in ("human", "avian"):
            if segment:
                tiers += [f"IAV/{h}/{segment}", f"IAV/all/{segment}"]
            tiers += [f"IAV/{h}"]
        if segment:
            tiers += [f"IAV/all/{segment}"]
        tiers += ["IAV/all", "IAV"]
        # dedupe, preserve order
        seen = set()
        return [t for t in tiers if not (t in seen or seen.add(t))]
    return [v]


def files_in_folder(entries, prefix):
    out = []
    for e in entries:
        folder = e["folder"].replace("\\", "/")
        if folder == prefix or folder.startswith(prefix + "/"):
            out.append(e)
    return out


def suitability(entry, ref_len_aa, parent_len_aa):
    """Return (suitable: bool, mode: str, reason: str)."""
    mol = entry["molecule"]
    state = entry["state"]
    med = entry["median_length"]

    if mol == "protein":
        ratio = med / ref_len_aa if ref_len_aa else 0
        if GENE_RATIO_MIN <= ratio <= GENE_RATIO_MAX:
            return True, ("protein_aligned" if state == "aligned" else "protein_raw"), ""
        return False, "", f"protein median {med} aa vs ref {ref_len_aa} aa (ratio {ratio:.2f})"

    if mol == "nucleotide":
        if med >= GENOME_MIN_NT:
            return True, ("genome_aligned" if state == "aligned" else "genome"), ""
        aa = med / 3.0
        ratio = aa / ref_len_aa if ref_len_aa else 0
        if GENE_RATIO_MIN <= ratio <= GENE_RATIO_MAX:
            return True, ("cds_aligned" if state == "aligned" else "cds"), ""
        if parent_len_aa and parent_len_aa >= PARENT_MIN_FACTOR * ref_len_aa:
            pratio = aa / parent_len_aa
            if GENE_RATIO_MIN <= pratio <= GENE_RATIO_MAX:
                return True, ("subregion_aligned" if state == "aligned" else "subregion"), ""
        return False, "", (f"nt median {med} ({aa:.0f} aa) vs ref {ref_len_aa} aa "
                           f"(ratio {ratio:.2f})")
    return False, "", "unknown molecule type"


def read_fasta_single(path):
    """Read a one-record fasta; return sequence or None."""
    seq = []
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                if seq:
                    break
                continue
            seq.append(line.strip())
    return "".join(seq) or None


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--benchmark-dir", default="Virogym_Benchmark")
    ap.add_argument("--registry", default=None,
                    help="default: <benchmark-dir>/" + REGISTRY_NAME)
    ap.add_argument("--inventory", default="input_inventory.json")
    ap.add_argument("--out", default="task_plan.json")
    ap.add_argument("--config-dir", default="planned_configs")
    args = ap.parse_args()

    registry_path = args.registry or os.path.join(args.benchmark_dir, REGISTRY_NAME)
    for p, what in [(registry_path, "registry CSV"), (args.inventory, "inventory JSON")]:
        if not os.path.isfile(p):
            print(f"[FATAL] {what} not found: {p}")
            return 1

    with open(args.inventory, encoding="utf-8") as fh:
        inventory = json.load(fh)
    entries = inventory.get("files", [])
    if not entries:
        print("[FATAL] inventory contains no files -- run scan_inputs.py first")
        return 1

    with open(registry_path, newline="", encoding="utf-8-sig") as fh:
        rows = list(csv.DictReader(fh))

    dms_on_disk = {f for f in os.listdir(args.benchmark_dir)
                   if f.lower().endswith(".csv") and f != REGISTRY_NAME} \
        if os.path.isdir(args.benchmark_dir) else set()

    # ---- pass 1: build groups keyed by de-gapped protein reference --------
    groups = {}          # ref hash -> group dict
    skipped = []         # benchmarks that will not get any task
    blocked = []

    for r in rows:
        bf = (r.get("benchmark_file") or "").strip()
        if not bf:
            continue
        stem = bf[:-4] if bf.lower().endswith(".csv") else bf
        virus = (r.get("virus_code") or "").strip()
        host = (r.get("host") or "").strip()
        mode = (r.get("recommended_reference") or "").strip() or "target_as_reference"
        target = degap(r.get("target_seq") or "")
        if not target:
            skipped.append({"benchmark_file": bf, "reason": "empty target_seq"})
            continue

        parent_len = None
        try:
            parent_len = int(float((r.get("uniprot_full_len") or "").strip()))
        except ValueError:
            pass

        if mode == "full_protein_plus_motif":
            full_seq = degap(r.get("full_protein_seq") or "")
            if not full_seq:
                ref_fasta = os.path.join("references", f"{stem}_fullprotein.fasta")
                if os.path.isfile(ref_fasta):
                    full_seq = degap(read_fasta_single(ref_fasta) or "")
            if not full_seq:
                blocked.append({
                    "benchmark_file": bf,
                    "reason": ("needs the full-length parent protein as "
                               f"references/{stem}_fullprotein.fasta or a "
                               "'full_protein_seq' registry column")})
                continue
            protein_ref = full_seq
            motif = target
        else:
            protein_ref = target
            motif = target

        h = ref_hash(protein_ref)
        if h not in groups:
            groups[h] = {
                "group_id": f"{virus}_{len(protein_ref)}aa_{h}",
                "virus_code": virus,
                "host": host,
                "segment": detect_segment(r.get("protein") or ""),
                "binning": "monthly" if virus == "SCV2" else "yearly",
                "reference_mode": mode,
                "ref_len": len(protein_ref),
                "ref_hash": h,
                "protein_reference_seq": protein_ref,
                "parent_len_aa": parent_len,
                "motifs": {},          # motif hash -> {tag, seq, benchmarks}
            }
        g = groups[h]
        mh = ref_hash(motif)
        if mh not in g["motifs"]:
            same_as_ref = (mh == g["ref_hash"])
            g["motifs"][mh] = {
                "tag": "whole" if same_as_ref else stem,
                "len": len(motif),
                "seq": motif,
                "benchmarks": [],
            }
        g["motifs"][mh]["benchmarks"].append(bf)

    # ---- pass 2: attach suitable FASTA files to each group -----------------
    plan_groups = []
    for h, g in sorted(groups.items(), key=lambda kv: kv[1]["group_id"]):
        tiers = folder_tiers(g["virus_code"], g["host"], g["segment"])
        chosen = []
        tier_used = None
        for t in tiers:
            cand = files_in_folder(entries, t)
            if cand:
                chosen = cand
                tier_used = t
                break

        suitable = []
        rejected = []
        in_segment_subfolder = bool(g["segment"]) and tier_used and \
            tier_used.rstrip("/").endswith("/" + g["segment"])
        for e in chosen:
            # IAV segment guard for files outside a dedicated segment subfolder
            if g["virus_code"] == "IAV" and g["segment"] and not in_segment_subfolder:
                hint = filename_segment_hint(e["rel_path"])
                if hint and hint != g["segment"]:
                    rejected.append({"file": e["rel_path"],
                                     "reason": f"filename suggests segment {hint}, "
                                               f"benchmark needs {g['segment']}"})
                    continue
            ok, fmode, why = suitability(e, g["ref_len"], g["parent_len_aa"])
            if ok:
                suitable.append({"path": e["path"], "rel_path": e["rel_path"],
                                 "mode": fmode, "state": e["state"],
                                 "n_records": e["n_records"]})
            else:
                rejected.append({"file": e["rel_path"], "reason": why})

        benchmark_names = sorted({b for m in g["motifs"].values()
                                  for b in m["benchmarks"]})
        if not chosen:
            for bf in benchmark_names:
                skipped.append({"benchmark_file": bf,
                                "reason": f"no FASTA folder matched "
                                          f"(tried: {', '.join(tiers)})"})
            continue
        if not suitable:
            for bf in benchmark_names:
                skipped.append({"benchmark_file": bf,
                                "reason": f"folder '{tier_used}' has no suitable "
                                          f"file for a {g['ref_len']} aa reference"})
            continue

        # write the run_pipeline config for this group
        os.makedirs(args.config_dir, exist_ok=True)
        cfg_path = os.path.join(args.config_dir, g["group_id"] + ".txt")
        with open(cfg_path, "w", encoding="utf-8") as fh:
            fh.write(">protein_reference\n")
            fh.write(g["protein_reference_seq"] + "\n")
            for m in g["motifs"].values():
                fh.write(f">motif_reference_{m['tag']}\n")
                fh.write(m["seq"] + "\n")

        plan_groups.append({
            "group_id": g["group_id"],
            "virus_code": g["virus_code"],
            "host": g["host"],
            "binning": g["binning"],
            "reference_mode": g["reference_mode"],
            "ref_len": g["ref_len"],
            "config": os.path.abspath(cfg_path),
            "motifs": [{"tag": m["tag"], "len": m["len"],
                        "benchmarks": m["benchmarks"]}
                       for m in g["motifs"].values()],
            "benchmarks": benchmark_names,
            "fasta_tier": tier_used,
            "fasta_files": suitable,
            "rejected_files": rejected,
        })

    plan = {
        "created_at": time.strftime("%Y-%m-%d %H:%M:%S"),
        "registry": os.path.abspath(registry_path),
        "inventory": os.path.abspath(args.inventory),
        "groups": plan_groups,
        "skipped": skipped,
        "blocked": blocked,
        "dms_missing_on_disk": sorted(
            bf for g in plan_groups for bf in g["benchmarks"]
            if bf not in dms_on_disk),
    }
    with open(args.out, "w", encoding="utf-8") as fh:
        json.dump(plan, fh, indent=2)

    # ---- report --------------------------------------------------------------
    n_jobs = sum(len(g["fasta_files"]) for g in plan_groups)
    n_benchmarks_planned = len({b for g in plan_groups for b in g["benchmarks"]})
    print()
    print("=" * 72)
    print("TASK PLAN SUMMARY")
    print("=" * 72)
    print(f"registry rows          : {len(rows)}")
    print(f"reference groups       : {len(groups)} unique de-gapped references")
    print(f"groups with work       : {len(plan_groups)}")
    print(f"benchmarks planned     : {n_benchmarks_planned}")
    print(f"alignment jobs         : {n_jobs}  (group x FASTA file)")
    print(f"benchmarks skipped     : {len(skipped)}")
    print(f"benchmarks blocked     : {len(blocked)}")
    print("-" * 72)
    for g in plan_groups:
        print(f"  {g['group_id']}: {len(g['benchmarks'])} benchmark(s), "
              f"{len(g['motifs'])} motif(s), {len(g['fasta_files'])} file(s) "
              f"[{g['fasta_tier']}], binning={g['binning']}")
    if skipped:
        print("\nskipped benchmarks:")
        for s in skipped:
            print(f"  - {s['benchmark_file']}: {s['reason']}")
    if blocked:
        print("\nblocked benchmarks:")
        for s in blocked:
            print(f"  - {s['benchmark_file']}: {s['reason']}")
    if plan["dms_missing_on_disk"]:
        print("\nwarning: planned benchmarks whose DMS file is not in "
              f"{args.benchmark_dir}:")
        for bf in plan["dms_missing_on_disk"]:
            print(f"  - {bf}")
    print(f"\nplan written to {os.path.abspath(args.out)}")
    print(f"configs in      {os.path.abspath(args.config_dir)}/")
    return 0


if __name__ == "__main__":
    sys.exit(main())
