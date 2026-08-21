# abcat

[English](./README.en.md) | [한국어](./README.md)

`abcat` is a Python CLI and library designed for unified sequence analysis of antibody Heavy and Light chains across both Variable Domains (VH, VL) and Constant Domains (CH1, CH2, CH3, CL), leveraging [ANARCII](https://github.com/oxpig/ANARCII) as a core dependency.

For Variable Domains, `abcat` uses ANARCII to perform numbering (IMGT, Kabat, Chothia, Martin, AHo) and to delineate CDR1, CDR2, and CDR3 loops. For Constant Domains, it performs sequence alignments and fingerprint matching to determine Isotype, Subclass, IMGT Allele, Isoallotype, and Allotype.

## Key Features

1. **Variable Domain (VH / VL) Analysis (Powered by ANARCII):**
   - Precise Variable Domain boundary detection and IMGT/Kabat/Chothia/Martin/AHo numbering.
   - Automatic extraction of CDR1, CDR2, CDR3 and Framework (FR1, FR2, FR3, FR4) sequences and lengths.
   - Heavy chain (VH) and Light chain (VK, VL) classification with confidence scores.

2. **Constant Domain Analysis (Isotype, Subclass, IMGT Allele, Isoallotype, Allotype):**
   - Automatic identification of Constant Regions (CH1, Hinge, CH2, CH3 / CL).
   - Isotype determination: Heavy chain (`IgG`, `IgA`, `IgM`, `IgE`, `IgD`), Light chain (`Kappa`, `Lambda`).
   - Subclass classification: `IgG1`, `IgG2`, `IgG3`, `IgG4`, `IgA1`, `IgA2`, `IGKC`, `IGLC1`–`IGLC7`.
   - IMGT Allele classification based on EU numbering key polymorphic positions:
     - `IGHG1*01`–`IGHG1*13`, `IGHG2*01`–`IGHG2*06`, `IGHG3*01`–`IGHG3*19`, `IGHG4*01`–`IGHG4*04`, `IGKC*01`–`IGKC*04`, `IGHE*01`–`IGHE*04`.
   - Allotype and Isoallotype determination via IMGT polymorphic residue fingerprints:
     - Heavy chain: `G1m1`, `G1m2`, `G1m3`, `G1m17`, `G1m27`, `G1m28`, `G2m..`, `G2m23`, `G3m5`–`G3m28`, `nG1m1`, `nG1m17`, `nG3m5`, `nG3m11`, `nG3m21`, `nG4m(a)`, `nG4m(b)`, `IGHE*01`–`IGHE*04`, etc.
     - Light chain: `Km1`, `Km1,2`, `Km3`, etc.

3. **CLI and Python API Support:**
   - Single sequence analysis and FASTA batch file processing.
   - Comprehensive reporting in CSV, JSON, and rich console table formats.

## Tech Stack

- Python: 3.14+
- Dependency Management: `uv`
- V-Domain Analysis: `anarcii` (https://github.com/oxpig/ANARCII)
- C-Domain Alignment & Sequence Comparison: `biopython`
- Data Models: `pydantic`
- CLI Framework: `typer` / `rich`
- Testing & Quality: `pytest`, `ruff`, `mypy`, `ty`

## Architecture & Workflow

```text
[Input Amino Acid Sequence (Full Heavy/Light Chain or V+C)]
                     │
                     ▼
       ┌──────────────────────────┐
       │     ANARCII Wrapper      │  (V-domain Numbering & Classification)
       └─────────────┬────────────┘
                     │
         ┌───────────┴───────────┐
         ▼                       ▼
┌──────────────────┐    ┌───────────────────────────────────┐
│ V-Domain Engine  │    │      Constant Domain Engine       │
│                  │    │                                   │
│ - VH / VL Check  │    │ - Constant Region Extraction      │
│ - IMGT Numbering │    │ - C-Gene Profile Alignment        │
│ - CDR1/2/3 Bounds│    │ - Isotype & Subclass Matching     │
│ - FR1/2/3/4      │    │ - IMGT Allele Classification      │
│                  │    │ - Allotype/Isoallotype Fingerprint│
└────────┬─────────┘    └─────────────────┬─────────────────┘
         │                                │
         └────────────────┬───────────────┘
                          ▼
            ┌──────────────────────────┐
            │ Unified Result (JSON/CSV)│
            └──────────────────────────┘
```

## IMGT Constant Domain Alleles (Reference)

`abcat` identifies IMGT alleles through key amino acid residue combinations (fingerprints) based on EU numbering:

### 1. IGHG1 (IgG1) Alleles

| IMGT Allele | Subclass | Allotype Association | Key EU Positions & Residues                          |
| :---------- | :------- | :------------------- | :--------------------------------------------------- |
| `IGHG1*01`  | IgG1     | G1m17,1              | 214K, 309L, 356D, 358L, 384N                         |
| `IGHG1*02`  | IgG1     | G1m17,1              | 214K, 309L, 356D, 358L, 384N (silent variant of *01) |
| `IGHG1*03`  | IgG1     | G1m3, nG1m1          | 199I, 214R, 309L, 356E, 358M, 384N                   |
| `IGHG1*04`  | IgG1     | G1m17,1,27           | 214K, 309L, 356D, 358L, 384N, 422I                   |
| `IGHG1*05`  | IgG1     | G1m17,1,28           | 214K, 309L, 356D, 358L, 384N, 435R, 436Y             |
| `IGHG1*07`  | IgG1     | G1m17,1,2            | 214K, 309L, 356D, 358L, 384N, 431G                   |
| `IGHG1*08`  | IgG1     | G1m3,1               | 199I, 214R, 309L, 356D, 358L, 384N                   |
| `IGHG1*11`  | IgG1     | G1m17,1 (309V)       | 214K, 309V, 356D, 358L, 384N                         |
| `IGHG1*13`  | IgG1     | G1m17,1 (296F)       | 214K, 296F, 309L, 356D, 358L, 384N                   |

### 2. IGHG2, IGHG3, IGHG4 & Light Chain Alleles

| Subclass | IMGT Allele       | Allotype / Marker | Key Residues (EU numbering)  |
| :------- | :---------------- | :---------------- | :--------------------------- |
| **IgG2** | `IGHG2*01`        | G2m.. (G2m-)      | 192S, 193L, 282V, 309V       |
|          | `IGHG2*02`        | G2m23             | 282M                         |
|          | `IGHG2*04`        | G2m(ny)           | 192N, 193F, 282V             |
|          | `IGHG2*06`        | -                 | 282V                         |
| **IgG3** | `IGHG3*01`        | G3m5,26           | 291P, 384S, 435R, 436F       |
|          | `IGHG3*04`        | G3m21             | 291L, 384S, 435H, 436Y       |
|          | `IGHG3*11`        | G3m5,13,14        | 397V, 419Q, 435R, 436F       |
|          | `IGHG3*12`        | G3m15,16          | 292W, 378M, 384N, 435H, 436Y |
|          | `IGHG3*13`        | G3m6,24           | 384S, 419E, 435H, 436Y       |
|          | `IGHG3*14`        | G3m10,27,28       | 384S, 422I, 435R, 436Y       |
|          | `IGHG3*17`        | G3m11             | 384S, 435H, 436Y             |
|          | `IGHG3*18`        | nG3m11            | 384N, 435H, 436Y             |
|          | `IGHG3*19`        | G3m16,5           | 292W, 384S, 435R, 436F       |
| **IgG4** | `IGHG4*01`        | nG4m(a)           | 309L                         |
|          | `IGHG4*02`        | nG4m(b)           | 309V                         |
|          | `IGHG4*03`        | -                 | 309L, 409K                   |
|          | `IGHG4*04`        | -                 | 309L, 445P                   |
| **IGKC** | `IGKC*01` / `*04` | Km3               | 153A, 191V                   |
|          | `IGKC*02`         | Km1,2             | 153A, 191L                   |
|          | `IGKC*03`         | Km1               | 153V, 191L                   |

## Key Allotype Fingerprints (Reference)

| Subclass / Chain | Marker / Allotype         | IMGT / EU Position              | Key Polymorphisms                                   |
| :--------------- | :------------------------ | :------------------------------ | :-------------------------------------------------- |
| IgG1 (CH1)       | `G1m17` vs `G1m3`         | CH1 IMGT 103, 120 (EU 199, 214) | K214 = G1m17, I199/R214 = G1m3, R214 = nG1m17       |
| IgG1 (CH3)       | `G1m1` vs `nG1m1`         | CH3 IMGT 12, 14 (EU 356, 358)   | D356/L358 = G1m1, E356/M358 = nG1m1                 |
| IgG1 (CH3)       | `G1m2`                    | CH3 IMGT 110 (EU 431)           | G431 = G1m2                                         |
| IgG1 (CH3)       | `G1m27`                   | CH3 IMGT 101 (EU 422)           | I422 = G1m27                                        |
| IgG1 (CH3)       | `G1m28`                   | CH3 IMGT 115, 116 (EU 435, 436) | R435/Y436 = G1m28                                   |
| IgG2 (CH2)       | `G2m23` vs `G2m..`        | CH2 EU 282                      | M282 = G2m23, V282 = G2m..                          |
| IgG3 (CH3)       | `G3m5` vs `nG3m5`         | CH3 EU 435, 436                 | R435/F436 = G3m5, H435/Y436 = nG3m5                 |
| IgG3 (CH3)       | `G3m26`                   | CH3 EU 436                      | R436 = G3m26                                        |
| IgG4 (CH2)       | `nG4m(a)` vs `nG4m(b)`    | CH2 EU 309                      | L309 = nG4m(a), V309 = nG4m(b)                      |
| IgE (CH1, CH2)   | `IGHE*01`–`IGHE*04`       | CH1 IMGT 41, CH2 IMGT 41        | C141/W246 (*01), W141/W246 (*02), C141/L246 (*03)   |
| Kappa (CL)       | `Km1` vs `Km1,2` vs `Km3` | CL IMGT 45, 83 (EU 153, 191)    | V153/L191 = Km1, A153/L191 = Km1,2, A153/V191 = Km3 |

> **Note on Engineered Antibodies (e.g., Trastuzumab):**  
> Trastuzumab's Heavy Chain Fc differs from the natural G1m1,17 form by harboring an engineered `E356-M358` motif in the CH3 region to remove the G1m1 epitope. Consequently, the analyzer identifies only the `G1m17` allotype (`G1m17 only`), and classifies the `E356/M358` segment as the `nG1m1` isoallotype.

## Installation

### 1. From PyPI

```bash
uv pip install abcat
```

### 2. Development Setup (using `uv`)

```bash
# Clone repository
git clone https://github.com/user/abcat.git
cd abcat

# Sync virtual environment including development dependencies
uv sync --extra dev

# Run tests and CLI
uv run pytest
uv run abcat analyze --sequence EVQLVES...
```

## CLI & Python API Usage Example

### 1. CLI Examples

```bash
# Basic analysis (IMGT scheme)
abcat analyze --sequence EVQLVESGGGLVQPGGSLRLSCAASGFTFSDHYMDWVRQAPGKGLEWVGRIRSKANSYATAYAASVKGRFTISRDDSKNTLYLQMNSLRAEDTAVYYCARFDAYWGQGTLVTVSSASTKGPSVFPLAPSSKSTSGGTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSSVVTVPSSSLGTQTYICNVNHKPSNTKVDKKVEPKSCDKTHTCPPCPAPELLGGPSVFLFPPKPKDTLMISRTPEVTCVVVDVSHEDPEVKFNWYVDGVEVHNAKTKPREEQYNSTYRVVSVLTVLHQDWLNGKEYKCKVSNKALPAPIEKTISKAKGQPREPQVYTLPPSRDELTKNQVSLTCLVKGFYPSDIAVEWESNGQPENNYKTTPPVLDSDGSFFLYSKLTVDKSRWQQGNVFSCSVMHEALHNHYTQKSLSLSPGK

# Sequence containing line breaks (wrap with quotes)
abcat analyze --sequence "EVQLLESGGGLVQPGGSLRLSCAASGIDLSTYAMGWVRQAPGKGLEWVGLIHRSGRTYYA
TWAKGRFTISKDSSKNTLYLQMNSLRAEDTAVYYCTRSYPDYSATASIWGQGTTVTVSSA
STKGPSVFPLAPSSKSTSGGTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSG
LYSLSSVVTVPSSSLGTQTYICNVNHKPSNTKVDKKVEPKSCDKTHTCPPCPAPELLGGP
SVFLFPPKPKDTLMISRTPEVTCVVVDVSHEDPEVKFNWYVDGVEVHNAKTKPREEQYNS
TYRVVSVLTVLHQDWLNGKEYKCKVSNKALPAPIEKTISKAKGQPREPQVYTLPPSREEM
TKNQVSLTCLVKGFYPSDIAVEWESNGQPENNYKTTPPVLDSDGSFFLYSKLTVDKSRWQ
QGNVFSCSVMHEALHNHYTQKSLSLSPGK"

# Specify numbering scheme (supports imgt, kabat, martin, chothia, aho)
abcat analyze --sequence EVQLVES... --scheme kabat
abcat vdomain --sequence EVQLVES... --scheme martin
abcat vdomain --sequence EVQLVES... --scheme chothia
abcat vdomain --sequence EVQLVES... --scheme aho

# Batch analysis from FASTA file to CSV
abcat batch --input examples/full_antibodies.fasta --output results.csv --format csv --scheme imgt
```

### 2. Python API Example

```python
from abcat import analyze_chain, analyze_vdomain

seq = "EVQLVESGGGLVQPGGSLRLSCAASGFTFSDHYMDWVRQAPGKGLEWVGRIRSKANSYATAYAASVKGRFTISRDDSKNTLYLQMNSLRAEDTAVYYCARFDAYWGQGTLVTVSSASTKGPSVFPLAPSSKSTSGGTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSSVVTVPSSSLGTQTYICNVNHKPSNTKVDKKVEPKSCDKTHTCPPCPAPELLGGPSVFLFPPKPKDTLMISRTPEVTCVVVDVSHEDPEVKFNWYVDGVEVHNAKTKPREEQYNSTYRVVSVLTVLHQDWLNGKEYKCKVSNKALPAPIEKTISKAKGQPREPQVYTLPPSRDELTKNQVSLTCLVKGFYPSDIAVEWESNGQPENNYKTTPPVLDSDGSFFLYSKLTVDKSRWQQGNVFSCSVMHEALHNHYTQKSLSLSPGK"

# Analyze with Martin scheme
result_martin = analyze_chain(seq, scheme="martin")
print("Scheme:", result_martin.v_analysis.scheme)
print("CDR1:", result_martin.v_analysis.cdrs["CDR1"].sequence)
print("CDR2:", result_martin.v_analysis.cdrs["CDR2"].sequence)
print("CDR3:", result_martin.v_analysis.cdrs["CDR3"].sequence)

# Analyze only V-Domain with Kabat scheme
v_kabat = analyze_vdomain(seq, scheme="kabat")
print("Kabat CDR3:", v_kabat.cdrs["CDR3"].sequence)
```

## PyPI Deployment (GitHub Actions)

When a GitHub Release is created, the workflow builds wheels and source distributions using `uv build` and automatically publishes them to PyPI via PyPI Trusted Publisher (OIDC).

## License

[MIT](./LICENSE)

## References

- GM Allotypes Reference:
  Currently testable (serologically) GM allotypes and amino acid substitutions.
  *J Immunol.* 2025 Dec 1;214(12):3181–3187. doi: [10.1093/jimmun/vkaf190](https://doi.org/10.1093/jimmun/vkaf190).
- Allelic Diversity & Isoallotypes Reference:
  Warrender AK, Kelton W. Beyond Allotypes: The Influence of Allelic Diversity in Antibody Constant Domains. *Front Immunol.* 2020 Aug 18;11:2016. doi: [10.3389/fimmu.2020.02016](https://doi.org/10.3389/fimmu.2020.02016). PMID: 32973808; PMCID: PMC7461860.
- Lefranc MP, Lefranc G. Using IMGT unique numbering for IG allotypes and Fc-engineered variants of effector properties and half-life of therapeutic antibodies. Immunol Rev. 2024 Nov;328(1):473-506. doi: 10.1111/imr.13399. Epub 2024 Oct 4. PMID: 39367563; PMCID: PMC11659927.
