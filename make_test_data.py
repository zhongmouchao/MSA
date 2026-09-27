#!/usr/bin/env python3
"""Generate small synthetic datasets to test run_pipeline.py end-to-end.
Writes into ./input. Safe to delete the generated files afterwards."""

import os
import random

random.seed(42)
HERE = os.path.dirname(os.path.abspath(__file__))
INPUT = os.path.join(HERE, "input")
os.makedirs(INPUT, exist_ok=True)

REF = "MFVFLVLLPLVSSQCVNLTTRTQLPPAYTNSFTRGVYYPDKVFRSSVLHSTQDLFLPFFSNVTWFHAIHVSGTNGTKRFDNPVLPFNDGVYFASTEKSNIIRGWIFGTTLDSKTQSLLIVNNATNVVIKVCEFQFCNDPFLGVYYHKNNKSWMESEFRVYSSANNCTFEYVSQPFLMDLEGKQGNFKNLREFVFKNIDGYFKIYSKHTPINLVRDLPQGFSALEPLVDLPIGINITRFQTLLALHRSYLTPGDSSSGWTAGAAAYYVGYLQPRTFLLKYNENGTITDAVDCALDPLSETKCTLKSFTVEKGIYQTSNFRVQPTESIVRFPNITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSKVGGNYNYLYRLFRKSNLKPFERDISTEIYQAGSTPCNGVEGFNCYFPLQSYGFQPTNGVGYQPYRVVVLSFELLHAPATVCGPKKSTNLVKNKCVNFNFNGLTGTGVLTESNKKFLPFQQFGRDIADTTDAVRDPQTLEILDITPCSFGGVSVITPGTNTSNQVAVLYQDVNCTEVPVAIHADQLTPTWRVYSTGSNVFQTRAGCLIGAEHVNNSYECDIPIGAGICASYQTQTNSPRRARSVASQSIIAYTMSLGAENSVAYSNNSIAIPTNFTISVTTEILPVSMTKTSVDCTMYICGDSTECSNLLLQYGSFCTQLNRALTGIAVEQDKNTQEVFAQVKQIYKTPPIKDFGGFNFSQILPDPSKPSKRSFIEDLLFNKVTLADAGFIKQYGDCLGDIAARDLICAQKFNGLTVLPPLLTDEMIAQYTSALLAGTITSGWTFGAGAALQIPFAMQMAYRFNGIGVTQNVLYENQKLIANQFNSAIGKIQDSLSSTASALGKLQDVVNQNAQALNTLVKQLSSNFGAISSVLNDILSRLDKVEAEVQIDRLITGRLQSLQTYVTQQLIRAAEIRASANLAATKMSECVLGQSKRVDFCGKGYHLMSFPQSAPHGVVFLHVTYVPAQEKNFTTAPAICHDGKAHFPREGVFVSNGTHWFVTQRNFYEPQIITTDNTFVSGNCDVVIGIVNNTVYDPLQPELDSFKEELDKYFKNHTSPDVDLGDISGINASVVNIQKEIDRLNEVAKNLNESLIDLQELGKYEQYIKWPWYIWLGFIAGLIAIVMVTIMLCCMTSCCSCLKGCCSCGSCCKFDEDDSEPVLKGVKLHYT"
MOTIF = "DPSKPSKRSFIEDLLFNKVTLADAGFIKQYGDCLGDIAARDLICAQKF"

CODON = {  # one common codon per amino acid
    "A": "GCT", "R": "CGT", "N": "AAT", "D": "GAT", "C": "TGT", "Q": "CAA",
    "E": "GAA", "G": "GGT", "H": "CAT", "I": "ATT", "L": "CTT", "K": "AAA",
    "M": "ATG", "F": "TTT", "P": "CCT", "S": "TCT", "T": "ACT", "W": "TGG",
    "Y": "TAT", "V": "GTT",
}
AAS = sorted(CODON)

DATES = ["2020-03-15", "2021-06-02", "2022-11-20", "2023-01-30", "2023-05",
         "2024", "2024-07-14", "2025-02-09", None]


def encode(protein):
    return "".join(CODON[a] for a in protein)


def mutate(protein, rate=0.01):
    prot = list(protein)
    for i in range(len(prot)):
        if random.random() < rate:
            prot[i] = random.choice(AAS)
    return "".join(prot)


def header(i):
    d = random.choice(DATES)
    return f"testseq/{i:04d}|EPI_ISL_{100000 + i}" + (f"|{d}" if d else "")


def write(name, records):
    with open(os.path.join(INPUT, name), "w") as h:
        for hdr, seq in records:
            h.write(f">{hdr}\n{seq}\n")


def write_config(name, protein=True, motif=True):
    with open(os.path.join(INPUT, name), "w") as h:
        if protein:
            h.write(f">protein_reference\n{REF}\n")
        if motif:
            h.write(f">motif_reference\n{MOTIF}\n")


# --- 1. example_demo: 60 raw CDS nucleotide sequences, both references ----
recs = []
for i in range(60):
    prot = mutate(REF, 0.005)
    nt = encode(prot)
    if i % 17 == 16:  # a few sequences with Ns
        p = random.randrange(0, len(nt) - 3)
        nt = nt[:p] + "N" + nt[p + 1:]
    recs.append((header(i), nt))
write("example_demo.fasta", recs)
write_config("example_demo.txt", protein=True, motif=True)

# --- 2. test_genome: 20 genome-like sequences (CDS + random flanks) -------
recs = []
for i in range(20):
    prot = mutate(REF, 0.005)
    flank5 = "".join(random.choices("ACGT", k=2500))
    flank3 = "".join(random.choices("ACGT", k=2500))
    recs.append((header(100 + i), flank5 + encode(prot) + flank3))
write("test_genome.fasta", recs)
write_config("test_genome.txt", protein=True, motif=True)

# --- 3. test_proteinonly: 20 raw CDS, protein reference only --------------
recs = []
for i in range(20):
    prot = mutate(REF, 0.005)
    recs.append((header(200 + i), encode(prot)))
write("test_proteinonly.fasta", recs)
write_config("test_proteinonly.txt", protein=True, motif=False)

print("Wrote example_demo, test_genome, test_proteinonly to", INPUT)
print("(test_aligned is created from example_demo's aligned output after the first run)")
