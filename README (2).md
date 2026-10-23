# Task 03 — Scripted analysis pipeline

A multi-step sequence analysis that runs end to end from **one command**:

```bash
python pipeline.py
```

Python 3.7+ only, standard library, nothing to install. It works from a clean checkout: the input file is in `data/`.

## What the pipeline does

Each step reads the previous step's file, and every intermediate file is kept in `work/`:

| Step | Input | Output | What it does |
|---|---|---|---|
| 1. clean | `data/example.fasta` | `01_clean.fasta` | upper-case, U→T, odd symbols → N, drops sequences shorter than `MIN_SEQ_LEN` |
| 2. stats | `01_clean.fasta` | `02_stats.csv` | length, GC %, base counts |
| 3. orfs | `01_clean.fasta` | `03_orfs.csv` | ORFs (ATG→stop) in all 6 frames, with coordinates and protein sequence |
| 4. proteins | `03_orfs.csv` | `04_top_proteins.fasta` | the longest proteins |
| 5. gc | `01_clean.fasta` | `05_gc_windows.csv` | sliding-window GC % |
| 6. report | `02` + `03` | `report.md`, `parameters.json` | summary table and the parameters used |

## Parameters

All parameters are in one block at the **top of `pipeline.py`** (`INPUT_FASTA`, `OUT_DIR`, `MIN_SEQ_LEN`, `MIN_ORF_NT`, `TOP_N_PROTEINS`, `GC_WINDOW`, `GC_STEP`). The values used for a run are saved in `work/parameters.json`. Input and output folder can also be overridden:

```bash
python pipeline.py --input my_sequences.fasta --out-dir my_results
```

## Using your own data

Replace `data/example.fasta` with any DNA FASTA file (for example a public gene or genome record from NCBI) and run `python pipeline.py` again.

## Decisions I made

- **No dependencies**, so the pipeline runs anywhere Python is installed.
- An ORF here means `ATG` to the first in-frame stop codon, at least `MIN_ORF_NT` long (default 90 nt = 30 codons), on either strand. Overlapping ORFs in different frames are all reported, and coordinates are always on the forward strand.
- Ambiguous bases (`N`) are kept in length but ignored in GC %.
- `data/example.fasta` is a **synthetic demo sequence** with two planted ORFs (one on each strand), used to check that the ORF finder returns the exact expected coordinates. Replace it with real data for real results.

## Check that it is reproducible

```bash
rm -rf work
python pipeline.py
```

Re-running gives identical files (only the date in `parameters.json` and `report.md` changes).
