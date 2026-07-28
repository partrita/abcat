# abmap

`abmap`는 [ANARCII](https://github.com/oxpig/ANARCII)를 핵심 의존성으로 활용하여 항체 Heavy chain 및 Light chain 서열의 Variable Domain(VH, VL) 및 Constant Domain(CH1, CH2, CH3, CL)을 통합 분석하는 Python CLI/라이브러리 도구입니다.

Variable Domain에 대해서는 ANARCII를 통한 번호 부여(IMGT, Kabat, Chothia 등) 및 CDR1, CDR2, CDR3 루프 분리를 수행하고, Constant Domain에 대해서는 정밀 서열 비교를 통해 Isotype, Subclass, Isoallotype, Allotype을 판정합니다.

## Key Features

1. Variable Domain (VH / VL) 분석 (ANARCII 기반):
   - ANARCII 통합으로 정확한 Variable Domain 경계 인식 및 IMGT/Kabat/Chothia/AHo 번호 부여.
   - CDR1, CDR2, CDR3 및 Framework (FR1, FR2, FR3, FR4) 구간 서열 및 길이 자동 추출.
   - Heavy chain(VH), Light chain(VK, VL) 자동 구분 및 confidence score 제공.

2. Constant Domain 분석 (Isotype, Subclass, Isoallotype, Allotype):
   - Constant Region (CH1, Hinge, CH2, CH3 / CL) 자동 도출.
   - Isotype 판정: Heavy chain (`IgG`, `IgA`, `IgM`, `IgE`, `IgD`), Light chain (`Kappa`, `Lambda`).
   - Subclass 판정: `IgG1`, `IgG2`, `IgG3`, `IgG4`, `IgA1`, `IgA2`, `IGKC`, `IGLC1`~`IGLC7`.
   - Allotype 및 Isoallotype 판정: IMGT 다형성 위치(Polymorphic positions) 분석을 통한 알로타입 마커 추출.
     - Heavy chain: `G1m17`, `G1m3`, `G1m1`, `G1m2`, `G2m23`, `G3m5`, `nG1m1`, `nG4m(a)`, `nG4m(b)` isoallotypes 등.
     - Light chain: `Km1`, `Km1,2`, `Km2`, `Km3` 등.

3. CLI 및 Python API 지원:
   - 단일 서열 및 FASTA 배치 파일 분석 지원.
   - CSV, JSON, 콘솔 테이블 형태의 풍부한 리포팅.

## Tech Stack

- Python: 3.14+
- Dependency Management: `uv`
- V-Domain Analysis: `anarcii` (https://github.com/oxpig/ANARCII)
- C-Domain Alignment & Sequence Comparison: `biopython`
- Data Models: `pydantic`
- CLI Framework: `typer` / `rich`
- Testing & Quality: `pytest`, `ruff`, `mypy`

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
│ - FR1/2/3/4      │    │ - Allotype/Isoallotype Fingerprint│
└────────┬─────────┘    └─────────────────┬─────────────────┘
         │                                │
         └────────────────┬───────────────┘
                          ▼
            ┌──────────────────────────┐
            │ Unified Result (JSON/CSV)│
            └──────────────────────────┘
```

## Key Allotype Fingerprints (Reference)

| Subclass / Chain | Marker / Allotype         | IMGT Position / Sequence Motif  | Key Polymorphisms                             |
| :--------------- | :------------------------ | :------------------------------ | :-------------------------------------------- |
| IgG1 (CH1)       | `G1m17` vs `G1m3`         | CH1 IMGT 120 (EU 214)           | K120 = G1m17, R120 = G1m3                     |
| IgG1 (CH3)       | `G1m1` vs `nG1m1`         | CH3 IMGT 356, 358 (EU 356, 358) | D356/L358 = G1m1, E356/M358 = nG1m1           |
| IgG1 (CH3)       | `G1m2`                    | CH3 IMGT 431 (EU 431)           | G431 = G1m2+                                  |
| IgG2 (CH2)       | `G2m23`                   | CH2 IMGT 282 (EU 282)           | V282 = G2m23+                                 |
| IgG4 (CH2)       | `nG4m(a)` vs `nG4m(b)`    | CH2 IMGT 115 (EU 309)           | L309 = nG4m(a), V309 = nG4m(b)                |
| Kappa (CL)       | `Km1` vs `Km1,2` vs `Km3` | CL IMGT 45, 83 (EU 153, 191)    | V45/L83 = Km1, A45/L83 = Km1,2, A45/V83 = Km3 |

> Note on Engineered Antibodies (e.g., Trastuzumab):
> Trastuzumab의 Heavy Chain Fc는 자연형 G1m1,17 형태가 아닌, CH3 영역이 `E356-M358`로 엔지니어링되어 G1m1 에피토프가 제거된 비정상(allotypic / engineered) 형태입니다. 따라서 판정 시스템은 `G1m17` 알로타입만 감지하며(`G1m17 only`), CH3의 `E356/M358` 서열은 `nG1m1` isoallotype으로 분류됩니다.

## CLI & Python API Usage Example

### 1. CLI Examples

```bash
# 기본 분석 (IMGT scheme)
uv run abmap analyze --sequence EVQLVESGGGLVQPGGSLRLSCAASGFTFSDHYMDWVRQAPGKGLEWVGRIRSKANSYATAYAASVKGRFTISRDDSKNTLYLQMNSLRAEDTAVYYCARFDAYWGQGTLVTVSSASTKGPSVFPLAPSSKSTSGGTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSSVVTVPSSSLGTQTYICNVNHKPSNTKVDKKVEPKSCDKTHTCPPCPAPELLGGPSVFLFPPKPKDTLMISRTPEVTCVVVDVSHEDPEVKFNWYVDGVEVHNAKTKPREEQYNSTYRVVSVLTVLHQDWLNGKEYKCKVSNKALPAPIEKTISKAKGQPREPQVYTLPPSRDELTKNQVSLTCLVKGFYPSDIAVEWESNGQPENNYKTTPPVLDSDGSFFLYSKLTVDKSRWQQGNVFSCSVMHEALHNHYTQKSLSLSPGK

# 줄바꿈이 포함된 서열 분석 (큰따옴표 "..." 로 감싸서 입력)
uv run abmap analyze --sequence "EVQLLESGGGLVQPGGSLRLSCAASGIDLSTYAMGWVRQAPGKGLEWVGLIHRSGRTYYA
TWAKGRFTISKDSSKNTLYLQMNSLRAEDTAVYYCTRSYPDYSATASIWGQGTTVTVSSA
STKGPSVFPLAPSSKSTSGGTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSG
LYSLSSVVTVPSSSLGTQTYICNVNHKPSNTKVDKKVEPKSCDKTHTCPPCPAPELLGGP
SVFLFPPKPKDTLMISRTPEVTCVVVDVSHEDPEVKFNWYVDGVEVHNAKTKPREEQYNS
TYRVVSVLTVLHQDWLNGKEYKCKVSNKALPAPIEKTISKAKGQPREPQVYTLPPSREEM
TKNQVSLTCLVKGFYPSDIAVEWESNGQPENNYKTTPPVLDSDGSFFLYSKLTVDKSRWQ
QGNVFSCSVMHEALHNHYTQKSLSLSPGK"

# 넘버링 체계(Scheme) 변경 예제 (imgt, kabat, martin, chothia, aho 지원)
uv run abmap analyze --sequence EVQLVES... --scheme kabat
uv run abmap vdomain --sequence EVQLVES... --scheme martin
uv run abmap vdomain --sequence EVQLVES... --scheme chothia
uv run abmap vdomain --sequence EVQLVES... --scheme aho

# FASTA 배치 파일 분석 및 CSV 저장
uv run abmap batch --input examples/full_antibodies.fasta --output results.csv --format csv --scheme imgt
```

### 2. Python API Example

```python
from abmap import analyze_chain, analyze_vdomain

seq = "EVQLVESGGGLVQPGGSLRLSCAASGFTFSDHYMDWVRQAPGKGLEWVGRIRSKANSYATAYAASVKGRFTISRDDSKNTLYLQMNSLRAEDTAVYYCARFDAYWGQGTLVTVSSASTKGPSVFPLAPSSKSTSGGTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSSVVTVPSSSLGTQTYICNVNHKPSNTKVDKKVEPKSCDKTHTCPPCPAPELLGGPSVFLFPPKPKDTLMISRTPEVTCVVVDVSHEDPEVKFNWYVDGVEVHNAKTKPREEQYNSTYRVVSVLTVLHQDWLNGKEYKCKVSNKALPAPIEKTISKAKGQPREPQVYTLPPSRDELTKNQVSLTCLVKGFYPSDIAVEWESNGQPENNYKTTPPVLDSDGSFFLYSKLTVDKSRWQQGNVFSCSVMHEALHNHYTQKSLSLSPGK"

# Martin scheme으로 분석
result_martin = analyze_chain(seq, scheme="martin")
print("Scheme:", result_martin.v_analysis.scheme)
print("CDR1:", result_martin.v_analysis.cdrs["CDR1"].sequence)
print("CDR2:", result_martin.v_analysis.cdrs["CDR2"].sequence)
print("CDR3:", result_martin.v_analysis.cdrs["CDR3"].sequence)

# Kabat scheme으로 V-Domain만 분석
v_kabat = analyze_vdomain(seq, scheme="kabat")
print("Kabat CDR3:", v_kabat.cdrs["CDR3"].sequence)
```

## PyPI Deployment (GitHub Actions)

본 프로젝트는 GitHub Release 발급 시 `uv build`를 실행하고, PyPI Trusted Publisher (OIDC)를 통해 최신 버전을 PyPI에 자동 게시하는 GitHub Action 워크플로우를 제공합니다.

## License

[MIT](./LICENSE)

## References

- GM Allotypes Reference:
  Currently testable (serologically) GM allotypes and amino acid substitutions.
  *J Immunol.* 2025 Dec 1;214(12):3181–3187. doi: [10.1093/jimmun/vkaf190](https://doi.org/10.1093/jimmun/vkaf190).
- Allelic Diversity & Isoallotypes Reference:
  Warrender AK, Kelton W. Beyond Allotypes: The Influence of Allelic Diversity in Antibody Constant Domains. *Front Immunol.* 2020 Aug 18;11:2016. doi: [10.3389/fimmu.2020.02016](https://doi.org/10.3389/fimmu.2020.02016). PMID: 32973808; PMCID: PMC7461860.