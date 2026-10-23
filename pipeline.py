#!/usr/bin/env python3
"""pipeline.py - sequence analysis pipeline, runs end to end with one command:

    python pipeline.py

Steps (each reads the previous step's file and keeps its own output in OUT_DIR):
    1. clean    data/example.fasta  ->  01_clean.fasta
    2. stats    01_clean.fasta      ->  02_stats.csv
    3. orfs     01_clean.fasta      ->  03_orfs.csv         (both strands, 3 frames each)
    4. proteins 03_orfs.csv         ->  04_top_proteins.fasta
    5. gc       01_clean.fasta      ->  05_gc_windows.csv
    6. report   02 + 03 + params    ->  report.md, parameters.json

Standard library only, nothing to install.
"""

# ============================ PARAMETERS ============================
INPUT_FASTA = "data/example.fasta"   # any DNA FASTA file
OUT_DIR = "work"                     # all intermediate + final files go here
MIN_SEQ_LEN = 50                     # step 1: drop sequences shorter than this (bp)
MIN_ORF_NT = 90                      # step 3: minimum ORF length incl. stop codon (nt)
TOP_N_PROTEINS = 5                   # step 4: how many longest proteins to export
GC_WINDOW = 100                      # step 5: sliding window size (bp)
GC_STEP = 50                         # step 5: window step (bp)
# ====================================================================

import argparse
import csv
import datetime
import json
from collections import Counter
from pathlib import Path

# Standard genetic code (NCBI table 1); codon index order is T, C, A, G
_AA = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
_B = "TCAG"
CODON = {a + b + c: _AA[16 * i + 4 * j + k]
         for i, a in enumerate(_B) for j, b in enumerate(_B) for k, c in enumerate(_B)}
_COMP = str.maketrans("ACGTN", "TGCAN")


# ------------------------------------------------------------ helpers
def read_fasta(path):
    records, header, chunks = [], None, []
    for raw in open(path):
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                records.append((header, "".join(chunks)))
            header, chunks = line[1:].strip(), []
        elif header is not None:
            chunks.append(line)
    if header is not None:
        records.append((header, "".join(chunks)))
    return records


def write_fasta(path, records, width=70):
    with open(path, "w") as f:
        for header, seq in records:
            f.write(f">{header}\n")
            for i in range(0, len(seq), width):
                f.write(seq[i:i + width] + "\n")


def seq_id(header):
    return header.split()[0]


def revcomp(seq):
    return seq.translate(_COMP)[::-1]


def translate(seq):
    return "".join(CODON.get(seq[i:i + 3], "X") for i in range(0, len(seq) - 2, 3))


def gc_percent(seq):
    c = Counter(seq)
    acgt = c["A"] + c["C"] + c["G"] + c["T"]
    return 100.0 * (c["G"] + c["C"]) / acgt if acgt else 0.0


def read_csv(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f))


# -------------------------------------------------------------- steps
def step1_clean(in_path, out_path):
    """Upper-case, U->T, replace non-ACGT symbols by N, drop short sequences."""
    kept, dropped = [], 0
    for header, seq in read_fasta(in_path):
        seq = seq.upper().replace("U", "T")
        seq = "".join(ch if ch in "ACGT" else "N" for ch in seq if ch.isalpha() or ch == "-")
        if len(seq) < MIN_SEQ_LEN:
            dropped += 1
            continue
        kept.append((header, seq))
    if not kept:
        raise SystemExit("No sequences left after cleaning - check INPUT_FASTA / MIN_SEQ_LEN")
    write_fasta(out_path, kept)
    print(f"  kept {len(kept)} sequence(s), dropped {dropped} shorter than {MIN_SEQ_LEN} bp")


def step2_stats(clean_path, out_path):
    with open(out_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["id", "length", "gc_percent", "A", "C", "G", "T", "N"])
        for header, seq in read_fasta(clean_path):
            c = Counter(seq)
            w.writerow([seq_id(header), len(seq), f"{gc_percent(seq):.2f}",
                        c["A"], c["C"], c["G"], c["T"], c["N"]])


def find_orfs(seq):
    """ATG...stop ORFs in all 6 frames. Coordinates are 1-based on the forward strand."""
    L, found = len(seq), []
    for strand, s in (("+", seq), ("-", revcomp(seq))):
        for frame in range(3):
            start = None
            for i in range(frame, L - 2, 3):
                cod = s[i:i + 3]
                if start is None:
                    if cod == "ATG":
                        start = i
                elif CODON.get(cod) == "*":
                    end = i + 3
                    if end - start >= MIN_ORF_NT:
                        a, b = (start + 1, end) if strand == "+" else (L - end + 1, L - start)
                        found.append((strand, frame + 1, a, b, end - start,
                                      translate(s[start:end - 3])))
                    start = None
    return found


def step3_orfs(clean_path, out_path):
    total = 0
    with open(out_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["id", "strand", "frame", "start", "end", "length_nt", "length_aa", "protein"])
        for header, seq in read_fasta(clean_path):
            orfs = sorted(find_orfs(seq), key=lambda o: -o[4])
            for strand, frame, a, b, n, prot in orfs:
                w.writerow([seq_id(header), strand, frame, a, b, n, len(prot), prot])
            total += len(orfs)
    print(f"  found {total} ORF(s) >= {MIN_ORF_NT} nt")


def step4_proteins(orfs_path, out_path):
    rows = sorted(read_csv(orfs_path), key=lambda r: -int(r["length_aa"]))[:TOP_N_PROTEINS]
    records = [(f"{r['id']}|ORF{i}|{r['strand']}{r['frame']}|{r['start']}-{r['end']}|{r['length_aa']}aa",
                r["protein"]) for i, r in enumerate(rows, start=1)]
    write_fasta(out_path, records)
    print(f"  wrote top {len(records)} protein(s)")


def step5_gc(clean_path, out_path):
    with open(out_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["id", "window_start", "window_end", "gc_percent"])
        for header, seq in read_fasta(clean_path):
            for i in range(0, max(len(seq) - GC_WINDOW, 0) + 1, GC_STEP):
                win = seq[i:i + GC_WINDOW]
                w.writerow([seq_id(header), i + 1, i + len(win), f"{gc_percent(win):.2f}"])


def step6_report(out_dir, stats_path, orfs_path):
    params = {
        "date_run": datetime.datetime.now().strftime("%Y-%m-%d %H:%M"),
        "INPUT_FASTA": INPUT_FASTA, "OUT_DIR": str(out_dir), "MIN_SEQ_LEN": MIN_SEQ_LEN,
        "MIN_ORF_NT": MIN_ORF_NT, "TOP_N_PROTEINS": TOP_N_PROTEINS,
        "GC_WINDOW": GC_WINDOW, "GC_STEP": GC_STEP,
    }
    (out_dir / "parameters.json").write_text(json.dumps(params, indent=2) + "\n")

    orfs_by_id = {}
    for r in read_csv(orfs_path):
        orfs_by_id.setdefault(r["id"], []).append(r)

    lines = ["# Pipeline report", "", "## Parameters", "", "| Parameter | Value |", "|---|---|"]
    lines += [f"| {k} | {v} |" for k, v in params.items()]
    lines += ["", "## Per-sequence summary", "",
              "| id | length (bp) | GC % | ORFs | longest ORF |", "|---|---|---|---|---|"]
    for r in read_csv(stats_path):
        orfs = orfs_by_id.get(r["id"], [])
        if orfs:
            top = orfs[0]
            longest = f"{top['length_aa']} aa ({top['strand']} strand, {top['start']}-{top['end']})"
        else:
            longest = "none"
        lines.append(f"| {r['id']} | {r['length']} | {r['gc_percent']} | {len(orfs)} | {longest} |")
    lines += ["", "## Intermediate files", "",
              "- `01_clean.fasta` cleaned input", "- `02_stats.csv` length, GC, composition",
              "- `03_orfs.csv` all ORFs with coordinates and protein sequence",
              "- `04_top_proteins.fasta` longest proteins", "- `05_gc_windows.csv` sliding-window GC",
              "- `parameters.json` parameters used for this run"]
    (out_dir / "report.md").write_text("\n".join(lines) + "\n")


# --------------------------------------------------------------- main
def main():
    global INPUT_FASTA
    ap = argparse.ArgumentParser(description="Sequence analysis pipeline (parameters: top of this file)")
    ap.add_argument("--input", default=INPUT_FASTA, help=f"input FASTA (default {INPUT_FASTA})")
    ap.add_argument("--out-dir", default=OUT_DIR, help=f"output folder (default {OUT_DIR})")
    args = ap.parse_args()
    INPUT_FASTA = args.input

    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    p = {n: out / n for n in ["01_clean.fasta", "02_stats.csv", "03_orfs.csv",
                              "04_top_proteins.fasta", "05_gc_windows.csv"]}

    print("[1/6] clean");    step1_clean(args.input, p["01_clean.fasta"])
    print("[2/6] stats");    step2_stats(p["01_clean.fasta"], p["02_stats.csv"])
    print("[3/6] orfs");     step3_orfs(p["01_clean.fasta"], p["03_orfs.csv"])
    print("[4/6] proteins"); step4_proteins(p["03_orfs.csv"], p["04_top_proteins.fasta"])
    print("[5/6] gc");       step5_gc(p["01_clean.fasta"], p["05_gc_windows.csv"])
    print("[6/6] report");   step6_report(out, p["02_stats.csv"], p["03_orfs.csv"])
    print(f"Done. Results in {out}/ (start with report.md)")


if __name__ == "__main__":
    main()
