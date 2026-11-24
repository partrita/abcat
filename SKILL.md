---
name: abcat
description: >-
  Analyze antibody heavy and light chain amino acid sequences with abcat.
  Extracts Variable domain (VH, VL) numbering and CDR1/CDR2/CDR3 loops via ANARCII,
  and classifies Constant domain (CH1, CH2, CH3, CL) isotype, subclass, IMGT allele,
  isoallotype, and allotype markers.
---

# abcat — Antibody Sequence Analyzer

Use this skill when the user wants to analyze antibody amino acid sequences for:
- Variable domain (`VH`, `VL`, `VK`) numbering (`imgt`, `kabat`, `martin`, `chothia`, `aho`)
- `CDR1`, `CDR2`, `CDR3`, and `FR1`–`FR4` extraction
- Constant domain (`CH1`, `Hinge`, `CH2`, `CH3`, `CL`) isotype and subclass classification
- IMGT allele, isoallotype, and allotype marker identification (such as `G1m17`, `G1m1`, `G1m3`, `G2m23`, `G3m*`, `nG4m(a)/(b)`, `Km1/Km3`)

## Requirements

- Python 3.14 or newer
- `abcat` installed (`uv pip install abcat` or `uv run abcat` inside the project)

## CLI Commands

### 1. Analyze a Single Sequence or FASTA File

```bash
abcat analyze --sequence "<AMINO_ACID_SEQUENCE>" --scheme imgt --verbose
```

Options:
- `--sequence`, `-s`: Single amino acid sequence. Wrap multiline sequences in quotes.
- `--input`, `-i`: Path to an input FASTA file.
- `--output`, `-o`: Path to save JSON or CSV output. Without `--output`, prints a summary table.
- `--format`, `-f`: Output format (`json` or `csv`, default `json`).
- `--scheme`: Numbering scheme (`imgt`, `kabat`, `martin`, `chothia`, `aho`, default `imgt`).
- `--verbose`: Show allotypes and `CDR1`/`CDR2`/`CDR3` sequences in the table.
- `--sequence-id`: Filter table output to one sequence ID.

### 2. Batch Process a FASTA File

```bash
abcat batch --input input.fasta --output results.csv --format csv --scheme imgt
```

### 3. Analyze Variable Domain Only

```bash
abcat vdomain --sequence "<AMINO_ACID_SEQUENCE>" --scheme kabat
```

Also accepts `--input` / `-i` for a FASTA file (uses the first record).

### 4. Analyze Constant Domain Only

```bash
abcat cdomain --sequence "<AMINO_ACID_SEQUENCE>"
```

## Python API

```python
from abcat import analyze_chain, analyze_vdomain, analyze_cdomain

result = analyze_chain(seq, scheme="imgt")
v_res = analyze_vdomain(seq, scheme="kabat")
c_res = analyze_cdomain(seq)
```

## Limitations

- Input must be amino acid sequences, not nucleotide sequences.
- Constant domain allele and allotype reference data targets human heavy and light chain constant genes.
