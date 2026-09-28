#!/usr/bin/env python3
"""
MSA Pipeline — streamlined, cluster-ready sequence alignment + motif extraction + counting
============================================================================================

Usage (terminal / cluster):

    python3 run_pipeline.py                          # process every dataset in ./input
    python3 run_pipeline.py --monthly --start-year 2000
    python3 run_pipeline.py --dataset my_data.fasta  # process just this one dataset
                                                     # (integrated tables are rebuilt
                                                     #  cumulatively over ./output)

Dataset layout in the input directory
-------------------------------------
Each dataset is a PAIR of files sharing the same base name:

    input/
        my_data.fasta        # sequences (nucleotide or protein, raw or aligned)
        my_data.txt          # config: reference sequence(s), FASTA-style

Config file format (FASTA-style; protein entry optional, motif entries repeatable):

    >protein_reference
    MFVFLVLLPLVSSQCVNLTTRTQ...
    >motif_reference
    DPSKPSKRSFIEDLLFNKVTLADAGFIKQYGDCLGDIAARDLICAQKF
    >motif_reference_fusion
    GLFGAIAGFIENGWEGMIDGWYG
    >motif_reference_stem
    GVTAA+PTITDQESLY+PKVRDQAG

Motif rules:
    * '+' joins parts of a CONCATENATED motif: parts are located IN ORDER on
      the protein reference (or alignment consensus) and the extracted columns
      are the concatenation of the per-part columns (non-contiguous OK).
    * Motif placement fallback chain: protein reference -> first data sequence
      containing all parts in order (handles subtype-specific motifs) ->
      column consensus.
    * Multiple motif_reference entries = extract & count each motif separately
      from the SAME alignment. The header suffix after "motif_reference"
      (e.g. "_fusion") becomes the motif tag used in output file names.
    * With a single unnamed motif, output files keep the classic names
      (<name>_motif_*.fasta, <name>_counts_by_year_wide.csv,
      integrated_counts_by_year_wide.csv); named/multiple motifs get the tag
      in the file names (<name>_<tag>_counts_by_year_wide.csv,
      integrated_<tag>_counts_by_year_wide.csv).

Intent is inferred from the config:

    protein only  -> ORF (if genomes) -> MSA -> counting (no motif extraction)
    motif only    -> motif extraction -> counting   (input MUST be pre-aligned)
    both          -> ORF (if genomes) -> MSA -> motif extraction -> counting

Auto-detection (per dataset, from the first N sequences):
    * nucleotide vs protein           (alphabet)
    * raw vs already-aligned          (gaps + uniform length)
    * whole-genome vs gene/CDS        (median length vs protein reference length)

Headers must contain a collection date: YYYY, YYYY-MM, or YYYY-MM-DD
(e.g. GISAID style  >hCoV-19/USA/XX/2023|2023-05-12). Records without a
parseable date are counted in the "unknown" column.

Dependencies: biopython, pandas, HMMER (hmmbuild + hmmalign on PATH).
"""

import argparse
import logging
import os
import re
import subprocess
import sys
import tempfile
from collections import Counter

import numpy as np
import pandas as pd
from Bio import AlignIO, SeqIO
from Bio.Align import MultipleSeqAlignment
from Bio.Data import CodonTable
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# =========================
# CONSTANTS
# =========================
FASTA_EXTS = (".fasta", ".fa", ".fna", ".faa", ".fas")
SAMPLE_N = 1000                 # records sampled for auto-detection
GENOME_LEN_FACTOR = 1.5         # median nt len / 3 > factor * ref aa len  -> genome
REF_LEN_MIN_FACTOR = 20         # HMM input filter: drop proteins longer than
                                # max(ref + 100, ref * 20). With k-mer-verified ORF
                                # selection, long ORFs are legitimate genes that
                                # merely CONTAIN a short reference (e.g. a 306-aa
                                # nsp3 fragment inside 4,406-aa ORF1a); over-length
                                # material becomes insert columns and is trimmed.
MAX_CONSENSUS_MISMATCH_FRAC = 0.25  # motif placement on consensus must be better than this

_table = CodonTable.unambiguous_dna_by_id[1]
_FWD = _table.forward_table
_STOPS = set(_table.stop_codons)
_VALID = set("ACGT")

# Date in header: YYYY, YYYY-MM, YYYY-MM-DD.
# GISAID uses several header conventions that put the collection date in
# different places:
#   hCoV-19 style: >hCoV-19/USA/XX/2025|EPI_ISL_20125422|2025-07-20  (date LAST)
#   EpiFlu style:  >2009-08-28|A/Singapore/ON1900/2009|A_/_H3N2|EPI_ISL_122791
#                  (date FIRST, isolate IDs like "ON1900" inside the strain name)
# So position alone cannot be trusted. Priority: (1) a hyphenated full date
# YYYY-MM[-DD] anywhere, checking the first and last pipe fields before the
# middle ones; (2) the strain-name year suffix "/YYYY"; (3) a standalone bare
# year. Isolate IDs (ON1900), accessions (EPI_ISL_2099) and pipeline
# annotations ([cds_aa=1,273]) must never be mistaken for the date.
FULL_DATE = re.compile(r"(?<![\d/])((?:19|20)\d{2})-(\d{1,2})(?:-(\d{1,2}))?(?!\d)")
STRAIN_YEAR = re.compile(r"/((?:19|20)\d{2})(?=\D|$)")
BARE_YEAR = re.compile(r"(?<![\w/])((?:19|20)\d{2})(?!\d)")


# =========================
# SMALL UTILITIES
# =========================
# numpy-vectorized translation: byte -> base index (0-3, 4 = invalid),
# codon id = a*25 + b*5 + c -> amino acid lookup table
_B2I = np.full(256, 4, dtype=np.uint8)
_B2I[ord("A")] = 0
_B2I[ord("C")] = 1
_B2I[ord("G")] = 2
_B2I[ord("T")] = 3
_CODON_LUT = np.full(125, ord("X"), dtype=np.uint8)
for _a, _b1 in (("A", 0), ("C", 1), ("G", 2), ("T", 3)):
    for _b, _b2 in (("A", 0), ("C", 1), ("G", 2), ("T", 3)):
        for _c, _b3 in (("A", 0), ("C", 1), ("G", 2), ("T", 3)):
            _codon = _a + _b + _c
            _aa = "*" if _codon in _STOPS else _FWD.get(_codon, "X")
            _CODON_LUT[_b1 * 25 + _b2 * 5 + _b3] = ord(_aa)
del _a, _b1, _b, _b2, _c, _b3, _codon, _aa


def translate_ambig_to_X(nt: str) -> str:
    nt = nt.upper()
    nt = nt[:len(nt) - (len(nt) % 3)]
    if not nt:
        return ""
    arr = np.frombuffer(nt.encode("ascii"), dtype=np.uint8)
    idx = _B2I[arr].reshape(-1, 3)
    cid = idx[:, 0] * 25 + idx[:, 1] * 5 + idx[:, 2]
    return _CODON_LUT[cid].tobytes().decode("ascii")


def choose_frame_by_min_stops(nt: str, test_nt: int = 1200):
    nt = nt.upper()
    test = nt[:min(len(nt), test_nt)]
    best = None  # (stop_count, -longest_orf, frame)
    for frame in (0, 1, 2):
        s = test[frame:]
        s = s[:len(s) - (len(s) % 3)]
        if len(s) < 3:
            continue
        prot = translate_ambig_to_X(s)
        stops = prot.count("*")
        longest_orf = max((len(x) for x in prot.split("*")), default=0)
        cand = (stops, -longest_orf, frame)
        if best is None or cand < best:
            best = cand
    return best[2] if best else None


def longest_orf_with_coords(protein: str):
    best_len = 0
    best = None
    pos = 0
    for seg in protein.split("*"):
        if len(seg) > best_len:
            best_len = len(seg)
            best = (seg, pos, pos + len(seg))
        pos += len(seg) + 1
    return best


def hamming(a: str, b: str) -> int:
    return sum(c1 != c2 for c1, c2 in zip(a, b))


def _norm_date(y, m, d):
    """Normalize a regex match to (year, month) or None. GISAID writes unknown
    month/day as 00, so 0 -> None; out-of-range month/day -> reject the match."""
    year = int(y)
    month = int(m) if m is not None else None
    day = int(d) if d is not None else None
    if month == 0:
        month, day = None, None
    if day == 0:
        day = None
    if month is not None and not (1 <= month <= 12):
        return None
    if day is not None and not (1 <= day <= 31):
        return None
    return year, month


def parse_date(description: str):
    """Return (year:int|None, month:int|None).

    Priority (see the FULL_DATE comment block above):
      1. first valid hyphenated full date, scanning the first pipe field,
         then the last, then the middle fields (covers both EpiFlu
         'date-first' and hCoV-19 'date-last' conventions);
      2. the LAST strain-name year suffix '/YYYY' (the year is the final
         component of a GISAID strain name, so later matches beat isolate
         IDs like 'ON1900');
      3. the first standalone bare year (not touching letters, '_', '/')."""
    fields = [f.strip() for f in description.split("|")]
    order = [fields[0]] + [fields[-1]] + fields[1:-1] if len(fields) > 1 else fields
    for field in order:
        for m in FULL_DATE.finditer(field):
            r = _norm_date(*m.groups())
            if r:
                return r
    matches = STRAIN_YEAR.findall(description)
    if matches:
        return int(matches[-1]), None
    m = BARE_YEAR.search(description)
    if m:
        return int(m.group(1)), None
    return None, None


def write_fasta(records, path):
    with open(path, "w") as h:
        SeqIO.write(records, h, "fasta")


# =========================
# CONFIG PARSING
# =========================
MOTIF_HEADER = re.compile(r"(?i)^motif_reference[\s_\-:.]*(\S*)$")
AA_ALPHABET = set("ACDEFGHIKLMNPQRSTVWYX*")


class Motif:
    """One motif reference: a tag (for output naming) and its parts.
    Multiple parts ('AAA+BBB' in the config) mean a concatenated motif:
    the parts are located IN ORDER on the protein reference / consensus and
    the extracted columns are the concatenation of the per-part columns."""

    def __init__(self, tag, parts):
        self.tag = tag
        self.parts = parts
        self.seq = "".join(parts)

    def __len__(self):
        return len(self.seq)


def parse_config(config_path):
    """Return (protein_ref|None, motifs:list[Motif]). Raises on problems.

    Config format (FASTA-style):
        >protein_reference          (optional; needed for MSA / ORF picking)
        <protein sequence>
        >motif_reference            (optional; repeatable)
        <motif sequence, parts joined by '+' for a concatenated motif>
        >motif_reference_stem       (extra motifs get their own tag -> own outputs)
        <motif sequence>
    An empty entry means that reference is intentionally not provided."""
    protein_ref = None
    motifs = []
    header, buf = None, []

    def flush():
        nonlocal protein_ref
        if header is None:
            return
        seq = "".join(buf).strip().upper().replace(" ", "").replace("\t", "")
        if not seq:
            print(f"NOTE: config entry '{header}' is empty -> ignored")
            return
        if "protein" in header.lower():
            if protein_ref is not None:
                raise ValueError("Config has more than one protein reference.")
            if any(c not in AA_ALPHABET for c in seq):
                raise ValueError("protein reference contains non-amino-acid characters.")
            protein_ref = seq
            return
        m = MOTIF_HEADER.match(header)
        if m:
            parts = [p for p in seq.split("+") if p]
            for p in parts:
                if any(c not in AA_ALPHABET for c in p):
                    raise ValueError(f"motif reference '{header}' contains "
                                     f"non-amino-acid characters.")
            if not parts:
                raise ValueError(f"motif reference '{header}' has no usable parts.")
            motifs.append([m.group(1), parts])   # tag finalized below
            return
        raise ValueError(
            f"Config header '{header}' must contain 'protein' or start with "
            f"'motif_reference'.")

    with open(config_path) as h:
        for line in h:
            line = line.rstrip("\n")
            if line.startswith(">"):
                flush()
                header, buf = line[1:].strip(), []
            elif line.strip():
                buf.append(line.strip())
        flush()

    # finalize tags: sanitize for filenames; defaults depend on motif count
    seen = set()
    out = []
    for i, (tag, parts) in enumerate(motifs, 1):
        tag = re.sub(r"[^A-Za-z0-9_\-]", "_", tag)
        if not tag:
            tag = "motif" if len(motifs) == 1 else f"motif{i}"
        if tag in seen:
            raise ValueError(f"Duplicate motif tag '{tag}' in config.")
        seen.add(tag)
        out.append(Motif(tag, parts))

    if protein_ref is None and not out:
        raise ValueError("Config has no '>protein_reference' or '>motif_reference' entry.")
    return protein_ref, out


# =========================
# AUTO-DETECTION
# =========================
def detect_input(fasta_path):
    """Sample the file; return dict with molecule/alignment/length info."""
    n = 0
    lengths = []
    has_gap = False
    has_protein_letter = False
    nt_chars = 0
    total_chars = 0

    with open(fasta_path, "r", encoding="latin-1") as h:
        for rec in SeqIO.parse(h, "fasta"):
            s = str(rec.seq).upper()
            lengths.append(len(s))
            if "-" in s:
                has_gap = True
            core = s.replace("-", "")
            if any(c in "EFILPQXZ*" for c in core):
                has_protein_letter = True
            nt_chars += sum(1 for c in core if c in "ACGTUN")
            total_chars += len(core)
            n += 1
            if n >= SAMPLE_N:
                break

    if n == 0:
        raise ValueError("Input FASTA contains no sequences.")

    is_protein = has_protein_letter or (total_chars > 0 and nt_chars / total_chars < 0.98)
    aligned = has_gap and len(set(lengths)) == 1
    lengths_sorted = sorted(lengths)
    median_len = lengths_sorted[len(lengths_sorted) // 2]

    return {
        "molecule": "protein" if is_protein else "nucleotide",
        "aligned": aligned,
        "median_len": median_len,
        "sampled": n,
    }


# =========================
# STAGE 1: ORF EXTRACTION (genomes) / DIRECT TRANSLATION (CDS)
# =========================
ORF_KMER_K = 10            # k-mer size for ORF-reference similarity scoring
ORF_MIN_KMER_SCORE = 0.3   # minimum fraction of ORF k-mers found in the reference


def _kmers(seq: str, k: int) -> set:
    return {seq[i:i + k] for i in range(len(seq) - k + 1) if "X" not in seq[i:i + k]}


def extract_orfs(fasta_path, ref_seq, min_orf_aa, logger, stats=None):
    """Whole genomes -> the ORF that matches the reference protein.

    Collects all ORFs (>= min_orf_aa) across all 3 forward frames, then picks
    the one with the highest k-mer coverage of the reference (fraction of the
    reference's 10-mers present in the ORF). Similarity, not length, is the
    criterion: random non-gene ORFs of comparable length are common in
    genomes, so length-only selection picks the wrong ORF for short proteins.
    Coverage is normalized by the REFERENCE, so a short fragment reference
    still matches the full-length gene ORF that contains it, and fused ORFs
    score high too. ORFs scoring < ORF_MIN_KMER_SCORE are rejected.
    Returns (nt_records, aa_records)."""
    ref_aa_len = len(ref_seq)
    ref_kmers = _kmers(ref_seq, ORF_KMER_K)
    nt_out, aa_out = [], []
    n_in = n_no_orf = 0
    with open(fasta_path, "r", encoding="latin-1") as handle:
        for rec in SeqIO.parse(handle, "fasta"):
            n_in += 1
            nt_full = str(rec.seq).upper()
            best = None  # (score, orf_len, orf_aa, nt_start, nt_end, frame)
            for frame in (0, 1, 2):
                coding = nt_full[frame:]
                coding = coding[:len(coding) - (len(coding) % 3)]
                if len(coding) < 3:
                    continue
                protein = translate_ambig_to_X(coding)
                pos = 0
                for seg in protein.split("*"):
                    if len(seg) >= min_orf_aa:
                        km = _kmers(seg, ORF_KMER_K)
                        score = (len(km & ref_kmers) / len(ref_kmers)) if ref_kmers else 0.0
                        if score >= ORF_MIN_KMER_SCORE and \
                                (best is None or score > best[0]):
                            best = (score, len(seg), seg,
                                    frame + pos * 3, frame + (pos + len(seg)) * 3,
                                    frame)
                    pos += len(seg) + 1
            if best is None:
                n_no_orf += 1
                continue
            score, orf_len, orf_aa, nt_start, nt_end, best_frame = best
            nt_orf = nt_full[nt_start:nt_end]
            # note: thousands separators in annotations (e.g. "orf_aa=4,406")
            # ensure they can never be mistaken for a year by date parsers
            desc = (f"{rec.description} [orf_aa={orf_len:,}] [frame +{best_frame + 1}] "
                    f"[kmer_score={score:.2f}]")
            nt_out.append(SeqRecord(Seq(nt_orf), id=rec.id, description=desc))
            aa_out.append(SeqRecord(Seq(orf_aa), id=rec.id, description=desc))
    logger.info(f"ORF extraction: {n_in:,} in, {len(nt_out):,} kept "
                f"(best k-mer match to ref {ref_aa_len} aa), {n_no_orf:,} with no "
                f"matching ORF >= {min_orf_aa} aa")
    if stats is not None:
        stats["orf_in"] = n_in
        stats["orf_kept"] = len(nt_out)
        stats["orf_dropped_no_orf"] = n_no_orf
    return nt_out, aa_out


def translate_cds(fasta_path, ref_aa_len, min_aa, logger, stats=None):
    """Gene/CDS input -> best-frame direct translation. Returns (nt_records, aa_records)."""
    nt_out, aa_out = [], []
    n_in = n_stop = n_short = 0
    min_len = min_aa
    with open(fasta_path, "r", encoding="latin-1") as handle:
        for rec in SeqIO.parse(handle, "fasta"):
            n_in += 1
            nt = str(rec.seq).upper()
            frame = choose_frame_by_min_stops(nt)
            if frame is None:
                n_short += 1
                continue
            coding = nt[frame:]
            coding = coding[:len(coding) - (len(coding) % 3)]
            if len(coding) < 3:
                n_short += 1
                continue
            protein = translate_ambig_to_X(coding)
            if "*" in protein:  # truncate at first internal stop
                protein = protein.split("*")[0]
                coding = coding[:len(protein) * 3]
                n_stop += 1
            if len(protein) < min_len:
                n_short += 1
                continue
            desc = f"{rec.description} [cds_aa={len(protein):,}] [frame +{frame + 1}]"
            nt_out.append(SeqRecord(Seq(coding), id=rec.id, description=desc))
            aa_out.append(SeqRecord(Seq(protein), id=rec.id, description=desc))
    logger.info(f"CDS translation: {n_in:,} in, {len(nt_out):,} kept, "
                f"{n_stop:,} truncated at stop, {n_short:,} too short")
    if stats is not None:
        stats["cds_in"] = n_in
        stats["cds_kept"] = len(nt_out)
        stats["cds_dropped_short"] = n_short
    return nt_out, aa_out


# =========================
# STAGE 2: REF-ONLY HMM ALIGNMENT + CODON-AWARE BACK-TRANSLATION
# =========================
def hmm_align_proteins(aa_records, ref_seq, logger, stats=None):
    """Ref-only HMM alignment. Returns protein MSA including the REF row (first)."""
    threshold = max(len(ref_seq) + 100, len(ref_seq) * REF_LEN_MIN_FACTOR)
    with tempfile.TemporaryDirectory() as tmpdir:
        ref_sto = os.path.join(tmpdir, "ref.sto")
        with open(ref_sto, "w") as f:
            f.write("# STOCKHOLM 1.0\nREF_SEQ %s\n//\n" % ref_seq)
        hmm = os.path.join(tmpdir, "ref.hmm")
        r = subprocess.run(["hmmbuild", hmm, ref_sto], capture_output=True, text=True)
        if r.returncode != 0:
            raise RuntimeError(f"hmmbuild failed: {r.stderr.strip()[:500]}")

        query_fa = os.path.join(tmpdir, "queries.fasta")
        kept, dropped = [], 0
        for r in aa_records:
            if len(r.seq) > threshold:
                dropped += 1
                continue
            kept.append(SeqRecord(Seq(str(r.seq).replace("X", "x")),
                                  id=r.id, description=r.description))
        if dropped:
            logger.info(f"HMM input: dropped {dropped:,} proteins longer than {threshold} aa")
        write_fasta([SeqRecord(Seq(ref_seq), id="REF_SEQ")] + kept, query_fa)

        final_sto = os.path.join(tmpdir, "aligned.sto")
        # note: hmmalign is single-threaded; it has no --cpu option
        cmd = ["hmmalign", "-o", final_sto, hmm, query_fa]
        r = subprocess.run(cmd, capture_output=True, text=True)
        if r.returncode != 0:
            # negative returncode = killed by a signal (e.g. -9 = SIGKILL,
            # typically the SLURM/cgroup OOM killer -> raise --mem)
            raise RuntimeError(f"hmmalign failed (exit code {r.returncode}): "
                               f"{r.stderr.strip()[:500]}")
        aln = AlignIO.read(final_sto, "stockholm")

    # restore original IDs by position (REF first)
    ids = ["REF_SEQ"] + [r.id for r in kept]
    descs = {r.id: r.description for r in kept}
    fixed = MultipleSeqAlignment([
        SeqRecord(Seq(str(rec.seq).upper()), id=ids[i],
                  description="reference" if i == 0 else descs.get(ids[i], ""))
        for i, rec in enumerate(aln)
    ])
    logger.info(f"HMM alignment: {len(fixed):,} rows (incl. reference), "
                f"width {fixed.get_alignment_length():,}")
    if stats is not None:
        stats["hmm_dropped_long"] = dropped
        stats["aligned"] = len(fixed) - 1
    return fixed


def back_translate(prot_aln, nt_records, logger, keep_cols=None):
    """Codon-aware back-translation; pads missing codons with '---'.
    REF row excluded. If keep_cols is given (reference match columns), only
    those columns are emitted — insert-column codons are consumed but dropped,
    so the output has exactly len(keep_cols) codons per sequence."""
    nt = {r.id: str(r.seq).upper() for r in nt_records}
    aln_len = prot_aln.get_alignment_length()
    keep_list = sorted(keep_cols) if keep_cols is not None else list(range(aln_len))
    keep_np = np.asarray(keep_list, dtype=np.int64)
    out, missing_nt, padded = [], 0, 0
    for rec in prot_aln:
        if rec.id == "REF_SEQ":
            continue
        coding = nt.get(rec.id)
        if coding is None:
            missing_nt += 1
            continue
        coding = coding[:len(coding) - (len(coding) % 3)]
        codons = np.frombuffer(coding.encode("ascii"), dtype=np.uint8).reshape(-1, 3)
        arr = np.frombuffer(str(rec.seq).encode("ascii"), dtype=np.uint8)
        is_res = arr != 45                       # residue vs '-' per column
        cum = np.cumsum(is_res) - 1              # codon index of each residue column
        res_kept = is_res[keep_np]               # kept columns carrying a residue
        codon_idx = cum[keep_np[res_kept]]
        ok = codon_idx < len(codons)
        padded += int((~ok).sum())
        result = np.full((len(keep_list), 3), 45, dtype=np.uint8)  # default '---'
        rows = np.nonzero(res_kept)[0][ok]
        result[rows] = codons[codon_idx[ok]]
        out.append(SeqRecord(Seq(result.tobytes().decode("ascii")),
                             id=rec.id, description=rec.description))
    logger.info(f"Back-translation: {len(out):,} written, {missing_nt:,} missing NT, "
                f"{padded:,} padded codons, width {len(keep_list) * 3:,} nt")
    return MultipleSeqAlignment(out)


def trim_to_reference(prot_aln, logger):
    """Drop insert columns: keep only columns where the reference has a residue.
    Returns (trimmed_alignment, keep_cols). The reference row becomes gap-free."""
    ref = next(r for r in prot_aln if r.id == "REF_SEQ")
    keep = np.fromiter((i for i, c in enumerate(str(ref.seq)) if c != "-"),
                       dtype=np.int64)
    dropped = prot_aln.get_alignment_length() - len(keep)
    trimmed = MultipleSeqAlignment([
        SeqRecord(Seq(np.frombuffer(str(r.seq).encode("ascii"), dtype=np.uint8)[keep]
                      .tobytes().decode("ascii")),
                  id=r.id, description=r.description)
        for r in prot_aln
    ])
    logger.info(f"Trimmed to reference columns: width {len(keep):,} "
                f"({dropped:,} insert columns removed)")
    return trimmed, keep.tolist()


# =========================
# STAGE 3: MOTIF EXTRACTION
# =========================
def _parts_in_seq(parts, seq):
    """True if every motif part occurs in `seq`, in order (concatenated motif)."""
    pos = 0
    for p in parts:
        i = seq.find(p, pos)
        if i == -1:
            return False
        pos = i + len(p)
    return True


def find_motif_columns_in_data(prot_aln, parts, logger):
    """Data-driven fallback (for subtype-specific motifs that are absent from
    both the protein reference and the column consensus): use the first aligned
    row whose UNGAPPED sequence contains all parts in order as the coordinate
    reference, and map the parts' residues to alignment columns via that row."""
    for r in prot_aln:
        if r.id == "REF_SEQ":
            continue
        aln_s = str(r.seq)
        raw = aln_s.replace("-", "")
        pos, ranges, ok = 0, [], True
        for p in parts:
            i = raw.find(p, pos)
            if i == -1:
                ok = False
                break
            ranges.append((i, i + len(p)))
            pos = i + len(p)
        if not ok:
            continue
        raw_to_aln, raw_pos = {}, 0
        for aln_pos, aa in enumerate(aln_s):
            if aa != "-":
                raw_to_aln[raw_pos] = aln_pos
                raw_pos += 1
        cols = []
        for s, e in ranges:
            cols.extend(raw_to_aln[i] for i in range(s, e))
        logger.info(f"Motif located on data sequence '{r.id}' "
                    f"(parts found in order): columns {cols[0] + 1}-{cols[-1] + 1} "
                    f"({len(parts)} part(s), {len(cols)} columns)")
        return cols
    raise ValueError("No sequence in the alignment contains all motif parts in "
                     "order. Check the motif sequence(s).")


def find_motif_columns_in_ref(prot_aln, parts):
    """Motif columns via the reference row (internal alignments). Parts are
    located IN ORDER on the ungapped reference; the column list is the
    concatenation of the per-part ranges (non-contiguous for concat motifs)."""
    ref = next(r for r in prot_aln if r.id == "REF_SEQ")
    ref_aln = str(ref.seq)
    raw = ref_aln.replace("-", "")
    raw_to_aln, raw_pos = {}, 0
    for aln_pos, aa in enumerate(ref_aln):
        if aa != "-":
            raw_to_aln[raw_pos] = aln_pos
            raw_pos += 1
    cols = []
    start_search = 0
    for part in parts:
        idx = raw.find(part, start_search)
        if idx == -1:
            raise ValueError(f"Motif part '{part[:20]}...' not found in order in the "
                             f"protein reference. Check the motif sequence.")
        cols.extend(raw_to_aln[i] for i in range(idx, idx + len(part)))
        start_search = idx + len(part)
    return cols


def _place_part_on_consensus(consensus, part, start_search, label, logger):
    """Hamming-place one motif part on the consensus at/after start_search."""
    L = len(part)
    best = min(((hamming(consensus[i:i + L], part), i)
                for i in range(start_search, len(consensus) - L + 1)), default=None)
    if best is None:
        raise ValueError(f"Alignment shorter than motif part ({label}).")
    mis, start = best
    if mis > int(MAX_CONSENSUS_MISMATCH_FRAC * L):
        raise ValueError(f"Motif part ({label}) could not be placed on the alignment "
                         f"consensus ({mis}/{L} mismatches). Is this the right "
                         f"alignment/motif?")
    logger.info(f"Motif part ({label}) placed on consensus at columns "
                f"{start + 1}-{start + L} ({mis} mismatches)")
    return start


def find_motif_columns_by_consensus(prot_aln, parts, logger):
    """Motif columns via the column consensus (external aligned input, no
    reference row). Parts are placed IN ORDER, each after the previous one."""
    aln_len = prot_aln.get_alignment_length()
    consensus = []
    for c in range(aln_len):
        col = Counter(str(r.seq)[c].upper() for r in prot_aln)
        col.pop("-", None)
        consensus.append(col.most_common(1)[0][0] if col else "-")
    consensus = "".join(consensus)
    cols = []
    start_search = 0
    for k, part in enumerate(parts, 1):
        start = _place_part_on_consensus(consensus, part, start_search,
                                         f"part {k}/{len(parts)}", logger)
        cols.extend(range(start, start + len(part)))
        start_search = start + len(part)
    return cols


def extract_motif_aligned(prot_aln, nuc_aln, motif_cols, motif, logger, stats=None,
                          min_called_frac=0.9):
    """Per-residue column extraction. A fragment is kept when its called
    (non-gap, non-X) positions cover >= min_called_frac of the motif length;
    use min_called_frac=1.0 for the old zero-tolerance behaviour. NT output
    is de-gapped; Ns are allowed only where the AA is X (same positions)."""
    L = len(motif)
    nuc_by_id = {r.id: str(r.seq) for r in nuc_aln} if nuc_aln is not None else {}
    prot_out, nuc_out = [], []
    n_ambig = n_missing_nt = n_short_nt = 0

    for rec in prot_aln:
        if rec.id == "REF_SEQ":
            continue
        aa_aln = str(rec.seq).upper()
        aa_frag = "".join(aa_aln[c] for c in motif_cols)
        bad = sum(1 for ch in aa_frag if ch in "-X")
        if L - bad < min_called_frac * L:
            n_ambig += 1
            continue
        nt_frag = ""
        if nuc_aln is not None:
            nt_aln = nuc_by_id.get(rec.id)
            if nt_aln is None:
                n_missing_nt += 1
                continue
            max_c = motif_cols[-1]
            if len(nt_aln) < max_c * 3 + 3:
                n_short_nt += 1
                continue
            nt_frag = "".join(nt_aln[c * 3:c * 3 + 3] for c in motif_cols).replace("-", "").upper()
            # every called (non-gap, non-X) position must contribute a full codon
            if len(nt_frag) < 3 * (L - bad):
                n_short_nt += 1
                continue
        prot_out.append(SeqRecord(Seq(aa_frag), id=rec.id, description=rec.description))
        nuc_out.append(SeqRecord(Seq(nt_frag), id=rec.id, description=rec.description))

    logger.info(f"Motif extraction: kept {len(prot_out):,}; dropped: "
                f"below {min_called_frac:.0%} called {n_ambig:,}, "
                f"missing NT {n_missing_nt:,}, short NT {n_short_nt:,}")
    if stats is not None:
        stats["motif_kept"] = len(prot_out)
        stats["motif_dropped_ambig"] = n_ambig
        stats["motif_dropped_nt"] = n_missing_nt + n_short_nt
    return prot_out, nuc_out


def extract_motif_from_aligned_nt(nuc_aln_records, parts, logger, stats=None,
                                  min_called_frac=0.9):
    """Motif-only config + pre-aligned NT input: consensus placement (parts
    located in order, each after the previous), codon columns, per-sequence
    translation. A fragment is kept when its called (non-gap, non-X, non-stop)
    positions cover >= min_called_frac of the motif length.
    Returns (prot_records, nuc_records, motif_cols)."""
    aln_len = len(nuc_aln_records[0].seq)
    if aln_len % 3 != 0:
        raise ValueError(f"Aligned NT length {aln_len} is not a multiple of 3 — "
                         f"is this a codon-aware alignment?")
    aa_cols = aln_len // 3

    # consensus codon per column -> consensus protein
    consensus = []
    for c in range(aa_cols):
        codons = Counter(str(r.seq)[c * 3:c * 3 + 3].upper() for r in nuc_aln_records)
        codons.pop("---", None)
        aa = "X"
        for codon, _ in codons.most_common():
            if all(b in _VALID for b in codon):
                aa = "*" if codon in _STOPS else _FWD.get(codon, "X")
                break
        consensus.append(aa)
    consensus = "".join(consensus)

    motif_cols = []
    start_search = 0
    for k, part in enumerate(parts, 1):
        start = _place_part_on_consensus(consensus, part, start_search,
                                         f"part {k}/{len(parts)}", logger)
        motif_cols.extend(range(start, start + len(part)))
        start_search = start + len(part)
    L = len(motif_cols)

    prot_out, nuc_out = [], []
    n_ambig = 0
    for r in nuc_aln_records:
        s = str(r.seq).upper()
        codons = [s[c * 3:c * 3 + 3] for c in motif_cols]
        aa_chars = []
        bad = 0
        for cd in codons:
            if "-" in cd:
                aa_chars.append("-")
                bad += 1
                continue
            a = "X" if any(b not in _VALID for b in cd) \
                else ("*" if cd in _STOPS else _FWD.get(cd, "X"))
            aa_chars.append(a)
            if a in "X*":
                bad += 1
        if L - bad < min_called_frac * L:
            n_ambig += 1
            continue
        aa = "".join(aa_chars)
        nt_frag = "".join(codons).replace("-", "")
        prot_out.append(SeqRecord(Seq(aa), id=r.id, description=r.description))
        nuc_out.append(SeqRecord(Seq(nt_frag), id=r.id, description=r.description))
    logger.info(f"Motif extraction: kept {len(prot_out):,}; dropped: "
                f"below {min_called_frac:.0%} called {n_ambig:,}")
    if stats is not None:
        stats["motif_kept"] = len(prot_out)
        stats["motif_dropped_ambig"] = n_ambig
    return prot_out, nuc_out, motif_cols


# =========================
# STAGE 4: COUNTING (always runs)
# =========================
def count_variants(nt_records, aa_records, monthly, start_year, out_csv, logger):
    """Deduplicate (nt, aa) variants and pivot to wide year or year-month counts."""
    aa_by_id = {r.id: str(r.seq).upper() for r in aa_records}
    counts = Counter()
    observed = set()           # (year, month|None)
    n_missing_date = 0

    for r in nt_records:
        nt = str(r.seq).upper().replace("-", "")
        aa = aa_by_id.get(r.id, translate_ambig_to_X(nt) if nt else "")
        year, month = parse_date(r.description)
        if year is None:
            n_missing_date += 1
        else:
            observed.add((year, month))
        key_date = (year, month) if monthly else (year, None)
        counts[(nt, aa, key_date)] += 1

    if not counts:
        raise ValueError("Nothing to count (no records reached the counting stage).")

    years = sorted({y for (y, m) in observed if y is not None})
    if start_year is None:
        start_year = min(years) if years else 2000
    max_year = max(years) if years else start_year

    # build ordered date-column labels
    if monthly:
        columns = []
        for y in range(start_year, max_year + 1):
            for m in range(1, 13):
                columns.append((f"{y}-{m:02d}", (y, m)))
            columns.append((f"{y}-unknown", (y, None)))   # year known, month unknown
        columns.append(("unknown", (None, None)))
    else:
        columns = [(f"before {start_year}", ("__before__", None))]
        columns += [(str(y), (y, None)) for y in range(start_year, max_year + 1)]
        columns.append(("unknown", (None, None)))

    def col_of(key_date):
        year, month = key_date
        if year is None:
            return (None, None)
        if year < start_year and not monthly:
            return ("__before__", None)
        return (year, month if monthly else None)

    rows = {}
    for (nt, aa, key_date), c in counts.items():
        row = rows.setdefault((nt, aa), Counter())
        row[col_of(key_date)] += c

    table = []
    for (nt, aa), row in rows.items():
        entry = {"Motif Nucleotide": nt, "Motif Peptide": aa}
        total = 0
        for label, col_key in columns:
            v = row.get(col_key, 0)
            entry[label] = v
            total += v
        entry["Count"] = total
        table.append(entry)

    df = pd.DataFrame(table).sort_values("Count", ascending=False).reset_index(drop=True)
    df.to_csv(out_csv, index=False)
    logger.info(f"Counting: {sum(counts.values()):,} records -> {len(df):,} unique variants; "
                f"{n_missing_date:,} records without a parseable date; wrote {out_csv}")
    return {"counted_records": sum(counts.values()),
            "unique_variants": len(df),
            "unknown_date": n_missing_date}


# =========================
# DATASET DRIVER
# =========================
def process_dataset(fasta_path, config_path, out_dir, args, logger):
    name = os.path.splitext(os.path.basename(fasta_path))[0]
    stats = {"dataset": name}
    with open(fasta_path, "rb") as h:
        stats["input_sequences"] = sum(1 for line in h if line.startswith(b">"))

    protein_ref, motifs = parse_config(config_path)
    logger.info(f"Config: protein_reference={'yes (%d aa)' % len(protein_ref) if protein_ref else 'no'}, "
                f"motifs={len(motifs)}"
                + (": " + "; ".join(f"{m.tag} ({len(m)} aa, {len(m.parts)} part(s))"
                                    for m in motifs) if motifs else ""))
    stats["protein_ref"] = "yes" if protein_ref else "no"
    stats["motif_ref"] = "yes" if motifs else "no"
    stats["motifs"] = ";".join(m.tag for m in motifs)

    det = detect_input(fasta_path)
    logger.info(f"Detected: {det['molecule']}, {'aligned' if det['aligned'] else 'raw'}, "
                f"median length {det['median_len']:,} (sampled {det['sampled']:,})")
    stats["detected"] = f"{det['molecule']}/{'aligned' if det['aligned'] else 'raw'}"

    want_motif = bool(motifs)

    # --- consistency double-checks -------------------------------------
    if det["aligned"] and protein_ref:
        logger.info("Input is already aligned -> MSA skipped (protein reference not needed for alignment)")
    if not det["aligned"] and not protein_ref:
        raise ValueError("Input is raw (unaligned) but the config has no protein_reference — "
                         "cannot build the MSA. Add the protein reference or provide aligned input.")
    if want_motif and protein_ref:
        for m in motifs:
            if not _parts_in_seq(m.parts, protein_ref):
                logger.warning(f"Motif '{m.tag}' is not found in order in the protein "
                               f"reference; falling back to consensus-based placement.")

    nt_records = aa_records = None        # unaligned paired sequences
    prot_aln = nuc_aln = None             # alignments
    motif_cols_map = {}                   # motif tag -> alignment columns

    # --- stage: ORF / translation (raw NT only) -------------------------
    if not det["aligned"] and det["molecule"] == "nucleotide":
        # minimum ORF/CDS length: absolute override, else 50% of reference length
        min_aa = (args.min_orf_aa if args.min_orf_aa is not None
                  else max(1, int(args.min_orf_frac * len(protein_ref))))
        logger.info(f"Minimum ORF/CDS length: {min_aa} aa")
        stats["min_aa"] = min_aa
        is_genome = det["median_len"] / 3 > GENOME_LEN_FACTOR * len(protein_ref)
        if is_genome:
            logger.info("Genome-length input -> ORF extraction")
            nt_records, aa_records = extract_orfs(fasta_path, protein_ref,
                                                  min_aa, logger, stats)
            write_fasta(nt_records, os.path.join(out_dir, f"{name}_orf_nucleotides.fasta"))
            write_fasta(aa_records, os.path.join(out_dir, f"{name}_orf_proteins.fasta"))
        else:
            logger.info("Gene/CDS-length input -> direct best-frame translation")
            nt_records, aa_records = translate_cds(fasta_path, len(protein_ref),
                                                   min_aa, logger, stats)
    elif not det["aligned"]:
        with open(fasta_path, "r", encoding="latin-1") as h:
            aa_records = [SeqRecord(Seq(str(r.seq).upper()), id=r.id, description=r.description)
                          for r in SeqIO.parse(h, "fasta")]
        stats["protein_in"] = len(aa_records)

    # --- stage: MSA ------------------------------------------------------
    if not det["aligned"]:
        prot_aln_full = hmm_align_proteins(aa_records, protein_ref, logger, stats)
        prot_aln, keep_cols = trim_to_reference(prot_aln_full, logger)
        write_fasta(list(prot_aln), os.path.join(out_dir, f"{name}_aligned_proteins.fasta"))
        if det["molecule"] == "nucleotide":
            # back-translate from the FULL alignment (insert codons must be
            # consumed), emitting only reference columns
            nuc_aln = back_translate(prot_aln_full, nt_records, logger,
                                     keep_cols=keep_cols)
            write_fasta(list(nuc_aln), os.path.join(out_dir, f"{name}_aligned_nucleotides.fasta"))
        if want_motif:
            for m in motifs:
                if _parts_in_seq(m.parts, protein_ref):
                    cols = find_motif_columns_in_ref(prot_aln, m.parts)
                    logger.info(f"Motif '{m.tag}' located on reference: columns "
                                f"{cols[0] + 1}-{cols[-1] + 1} "
                                f"({len(m.parts)} part(s), {len(cols)} columns)")
                else:
                    try:
                        cols = find_motif_columns_in_data(prot_aln, m.parts, logger)
                        logger.info(f"Motif '{m.tag}' placed via a data sequence")
                    except ValueError:
                        cols = find_motif_columns_by_consensus(prot_aln, m.parts, logger)
                motif_cols_map[m.tag] = cols
    else:
        with open(fasta_path, "r", encoding="latin-1") as h:
            records = list(SeqIO.parse(h, "fasta"))
        if det["molecule"] == "protein":
            prot_aln = MultipleSeqAlignment(
                [SeqRecord(Seq(str(r.seq).upper()), id=r.id, description=r.description)
                 for r in records])
        else:
            nuc_aln = records  # keep as list of SeqRecord

    # --- stage: motif extraction (per motif) ------------------------------
    tag = "month" if args.monthly else "year"
    legacy_names = len(motifs) == 1 and motifs[0].tag == "motif"
    count_jobs = []        # (count_nt, count_aa, out_csv)
    motif_detail = []

    if want_motif:
        for m in motifs:
            mstats = {}
            if det["aligned"] and det["molecule"] == "nucleotide":
                prot_motif, nuc_motif, _ = extract_motif_from_aligned_nt(
                    nuc_aln, m.parts, logger, mstats,
                    min_called_frac=args.motif_min_called)
            else:
                cols = motif_cols_map.get(m.tag)
                if cols is None:
                    try:
                        cols = find_motif_columns_in_data(prot_aln, m.parts, logger)
                    except ValueError:
                        cols = find_motif_columns_by_consensus(prot_aln, m.parts, logger)
                prot_motif, nuc_motif = extract_motif_aligned(
                    prot_aln, nuc_aln, cols, m.seq, logger, mstats,
                    min_called_frac=args.motif_min_called)
            prefix = f"{name}_motif" if legacy_names else f"{name}_{m.tag}_motif"
            write_fasta(prot_motif, os.path.join(out_dir, f"{prefix}_proteins.fasta"))
            if nuc_motif and any(len(r.seq) for r in nuc_motif):
                write_fasta(nuc_motif, os.path.join(out_dir, f"{prefix}_nucleotides.fasta"))
            csv_path = (os.path.join(out_dir, f"{name}_counts_by_{tag}_wide.csv")
                        if legacy_names else
                        os.path.join(out_dir, f"{name}_{m.tag}_counts_by_{tag}_wide.csv"))
            count_jobs.append((nuc_motif, prot_motif, csv_path, m.tag))
            for k in ("motif_kept", "motif_dropped_ambig", "motif_dropped_nt"):
                stats[k] = stats.get(k, 0) + mstats.get(k, 0)
            motif_detail.append(f"{m.tag}: kept {mstats.get('motif_kept', 0):,}")
    else:
        # counting on full sequences (ungapped)
        if nuc_aln is not None:
            count_nt = [SeqRecord(Seq(str(r.seq).replace("-", "")), id=r.id,
                                  description=r.description) for r in nuc_aln]
            if prot_aln is not None:
                count_aa = [r for r in prot_aln if r.id != "REF_SEQ"]
            else:
                count_aa = []  # AA derived from NT during counting
        elif prot_aln is not None:
            count_aa = [r for r in prot_aln if r.id != "REF_SEQ"]
            count_nt = [SeqRecord(Seq(""), id=r.id, description=r.description) for r in count_aa]
        else:
            raise ValueError("Nothing to count: no sequences survived the pipeline.")
        count_jobs.append((count_nt, count_aa,
                           os.path.join(out_dir, f"{name}_counts_by_{tag}_wide.csv"), None))

    # --- stage: counting (always, once per motif) --------------------------
    total_records = total_variants = total_unknown = 0
    csv_paths = []
    for cnt_nt, cnt_aa, csv_path, mtag in count_jobs:
        count_stats = count_variants(cnt_nt, cnt_aa, args.monthly, args.start_year,
                                     csv_path, logger)
        total_records += count_stats["counted_records"]
        total_variants += count_stats["unique_variants"]
        total_unknown += count_stats["unknown_date"]
        csv_paths.append(csv_path)
        if mtag is not None:
            motif_detail.append(f"{mtag}: {count_stats['unique_variants']:,} variants")
    stats["counted_records"] = total_records
    stats["unique_variants"] = total_variants
    stats["unknown_date"] = total_unknown
    stats["motif_detail"] = "; ".join(motif_detail)
    stats["counts_csv"] = ";".join(csv_paths)
    return stats


# =========================
# RUN-LEVEL OUTPUTS
# =========================
SUMMARY_COLUMNS = [
    "dataset", "status", "detected", "protein_ref", "motif_ref", "motifs",
    "input_sequences", "orf_in", "orf_kept", "orf_dropped_no_orf",
    "cds_in", "cds_kept", "cds_dropped_short", "protein_in",
    "hmm_dropped_long", "aligned",
    "motif_kept", "motif_dropped_ambig", "motif_dropped_nt",
    "counted_records", "unique_variants", "unknown_date", "motif_detail", "error",
]


def write_run_summary(all_stats, out_path):
    """Per-dataset summary table + a TOTAL row for the numeric columns."""
    df = pd.DataFrame(all_stats)
    for c in SUMMARY_COLUMNS:
        if c not in df.columns:
            df[c] = ""
    df = df[SUMMARY_COLUMNS]
    total = {"dataset": "TOTAL", "status": ""}
    for c in SUMMARY_COLUMNS:
        if c in ("dataset", "status", "detected", "protein_ref", "motif_ref",
                 "motifs", "motif_detail", "error"):
            continue
        total[c] = pd.to_numeric(df[c], errors="coerce").sum()
    df = pd.concat([df, pd.DataFrame([total])], ignore_index=True)
    df.to_csv(out_path, index=False)
    return df


def integrate_counts(output_dir, csv_paths, out_path, tag="year"):
    """Merge per-dataset counts CSVs into one integrated wide CSV."""
    frames = []
    for p in csv_paths:
        if os.path.exists(p):
            # keep_default_na=False: protein-only inputs leave the nucleotide
            # column empty; "" must not become NaN or groupby drops those rows
            frames.append(pd.read_csv(p, keep_default_na=False))
    if not frames:
        return None

    # unify year/month columns across datasets
    if tag == "year":
        starts = []
        for df in frames:
            for c in df.columns:
                if c.startswith("before "):
                    starts.append(int(c.split()[1]))
        start = min(starts) if starts else None
        for i, df in enumerate(frames):
            remap = {}
            for c in df.columns:
                if c.startswith("before "):
                    remap[c] = f"before {start}"
                elif c.isdigit() and start is not None and int(c) < start:
                    remap[c] = f"before {start}"
            if remap:
                long = df.melt(id_vars=["Motif Nucleotide", "Motif Peptide"],
                               var_name="col", value_name="n")
                long["col"] = long["col"].map(lambda c: remap.get(c, c))
                df = (long.groupby(["Motif Nucleotide", "Motif Peptide", "col"],
                                   as_index=False)["n"].sum()
                          .pivot_table(index=["Motif Nucleotide", "Motif Peptide"],
                                       columns="col", values="n", fill_value=0)
                          .reset_index())
                df.columns.name = None
            frames[i] = df

    merged = pd.concat(frames, ignore_index=True)
    merged = merged.groupby(["Motif Nucleotide", "Motif Peptide"],
                            as_index=False).sum(numeric_only=True)

    # order columns: identifiers, date columns sorted, unknown, Count
    def col_key(c):
        if c.startswith("before "):
            return (0, c)
        if c == "unknown":
            return (2, "")
        if c == "Count":
            return (3, "")
        return (1, c)

    date_cols = sorted([c for c in merged.columns
                        if c not in ("Motif Nucleotide", "Motif Peptide", "Count")],
                       key=col_key)
    merged = merged[["Motif Nucleotide", "Motif Peptide"] + date_cols]
    merged[date_cols] = merged[date_cols].fillna(0).astype(int)
    merged["Count"] = merged[date_cols].sum(axis=1)
    merged = merged.sort_values("Count", ascending=False).reset_index(drop=True)

    merged.to_csv(out_path, index=False)
    return out_path, len(merged), int(merged["Count"].sum())


# =========================
# MAIN
# =========================
def main():
    here = os.path.dirname(os.path.abspath(__file__))
    ap = argparse.ArgumentParser(description="Batch MSA + motif extraction + counting pipeline")
    ap.add_argument("--input-dir", default=os.path.join(here, "input"))
    ap.add_argument("--output-dir", default=os.path.join(here, "output"))
    ap.add_argument("--monthly", action="store_true",
                    help="split counts by year AND month (YYYY-MM columns)")
    ap.add_argument("--start-year", type=int, default=None,
                    help="first year column (default: earliest year observed)")
    ap.add_argument("--min-orf-frac", type=float, default=0.5,
                    help="minimum ORF/CDS length as a fraction of the reference "
                         "protein length (default 0.5)")
    ap.add_argument("--min-orf-aa", type=int, default=None,
                    help="absolute minimum ORF/CDS length in aa (overrides --min-orf-frac)")
    ap.add_argument("--motif-min-called", type=float, default=0.9,
                    help="minimum fraction of motif positions that must be called "
                         "(non-gap, non-X) to keep a sequence (default 0.9; "
                         "use 1.0 for zero tolerance)")
    ap.add_argument("--dataset", default=None,
                    help="process only this one dataset (FASTA file name or base "
                         "name in the input dir) instead of every FASTA found")
    args = ap.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    fastas = sorted(f for f in os.listdir(args.input_dir)
                    if f.lower().endswith(FASTA_EXTS))
    if args.dataset:
        base = os.path.splitext(args.dataset)[0]
        match = [f for f in fastas
                 if f == args.dataset or os.path.splitext(f)[0] == base]
        if not match:
            print(f"Dataset '{args.dataset}' not found in {args.input_dir}.\n"
                  f"Available: {', '.join(fastas)}")
            sys.exit(1)
        fastas = match
    if not fastas:
        print(f"No FASTA files found in {args.input_dir}")
        sys.exit(1)

    print(f"Found {len(fastas)} dataset(s) in {args.input_dir}\n")
    ok, failed, all_stats = [], [], []
    for fa in fastas:
        name = os.path.splitext(fa)[0]
        config = os.path.join(args.input_dir, name + ".txt")
        fasta_path = os.path.join(args.input_dir, fa)
        out_dir = os.path.join(args.output_dir, name)
        os.makedirs(out_dir, exist_ok=True)

        logger = logging.getLogger(name)
        logger.setLevel(logging.INFO)
        logger.handlers.clear()
        fmt = logging.Formatter("%(asctime)s [%(levelname)s] %(message)s", "%H:%M:%S")
        fh = logging.FileHandler(os.path.join(out_dir, f"{name}.log"), mode="w")
        fh.setFormatter(fmt)
        sh = logging.StreamHandler(sys.stdout)
        sh.setFormatter(fmt)
        logger.addHandler(fh)
        logger.addHandler(sh)

        print(f"{'=' * 60}\nDATASET: {name}\n{'=' * 60}")
        try:
            if not os.path.exists(config):
                # fall back to a single shared config in the input dir
                shared = next((os.path.join(args.input_dir, c) for c in
                               ("reference.txt", "config.txt")
                               if os.path.exists(os.path.join(args.input_dir, c))), None)
                if shared is None:
                    raise FileNotFoundError(
                        f"missing config file: {name}.txt (or a shared reference.txt)")
                config = shared
                logger.info(f"Using shared config: {os.path.basename(shared)}")
            stats = process_dataset(fasta_path, config, out_dir, args, logger)
            stats["status"] = "OK"
            stats["error"] = ""
            all_stats.append(stats)
            ok.append(name)
            print(f"--> {name}: DONE\n")
        except Exception as e:
            logger.error(f"FAILED: {e}")
            failed.append((name, str(e)))
            all_stats.append({"dataset": name, "status": "FAILED", "error": str(e)})
            print(f"--> {name}: FAILED ({e})\n")

    # ---- run-level outputs: summary table + integrated counts -----------
    tag = "month" if args.monthly else "year"
    summary_path = os.path.join(args.output_dir, "run_summary.csv")
    summary_df = write_run_summary(all_stats, summary_path)

    print(f"{'=' * 60}\nRUN SUMMARY\n{'=' * 60}")
    print(summary_df.to_string(index=False))

    # Integration scans EVERY dataset dir in the output folder (not just the
    # ones processed in this run), so single-dataset runs (--dataset, e.g. one
    # SLURM job per FASTA) rebuild the integrated tables cumulatively.
    # Group per-dataset count CSVs by motif infix: "" (legacy/single motif)
    # or "_<motif tag>" -> one integrated CSV per motif.
    groups = {}
    suffix = f"_counts_by_{tag}_wide.csv"
    for n in sorted(os.listdir(args.output_dir)):
        d = os.path.join(args.output_dir, n)
        if not os.path.isdir(d):
            continue
        for f in sorted(os.listdir(d)):
            if f == n + suffix or (f.startswith(n + "_") and f.endswith(suffix)):
                infix = f[len(n):len(f) - len(suffix)]   # "" or "_fusion"
                groups.setdefault(infix, []).append(os.path.join(d, f))
    for infix, paths in sorted(groups.items()):
        out_path = os.path.join(args.output_dir,
                                f"integrated{infix}_counts_by_{tag}_wide.csv")
        result = integrate_counts(args.output_dir, paths, out_path, tag)
        if result:
            out_path, n_variants, total = result
            label = infix.lstrip("_") or "motif"
            print(f"\nIntegrated counts [{label}]: {n_variants:,} unique variants, "
                  f"{total:,} total records -> {out_path}")

    print(f"\nRun summary written to {summary_path}")
    print(f"SUMMARY: {len(ok)} succeeded, {len(failed)} failed")
    for n, e in failed:
        print(f"  FAILED {n}: {e}")
    print(f"Results in {args.output_dir}")
    sys.exit(1 if failed and not ok else 0)


if __name__ == "__main__":
    main()
