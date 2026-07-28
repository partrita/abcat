import csv
import io
from pathlib import Path

from abmap.schemas import BatchAnalysisResult, FullChainAnalysis


def clean_sequence(seq: str) -> str:
    """Cleans an input amino acid sequence by removing whitespace, line breaks, tabs, non-alphabet characters, and FASTA header if present."""
    if not seq:
        return ""
    import re

    lines = seq.strip().splitlines()
    filtered_lines = []
    for line in lines:
        line_str = line.strip()
        if line_str.startswith(">"):
            continue
        filtered_lines.append(line_str)

    raw = "".join(filtered_lines)
    return "".join(re.findall(r"[A-Za-z]+", raw)).upper()


def parse_fasta(fasta_content: str) -> list[tuple[str, str]]:
    """Parses a FASTA string into a list of (header, sequence) tuples."""
    records = []
    current_header = None
    current_seq = []

    for line in fasta_content.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if current_header:
                records.append((current_header, clean_sequence("".join(current_seq))))
            current_header = line[1:].strip()
            current_seq = []
        else:
            current_seq.append(line)

    if current_header:
        records.append((current_header, clean_sequence("".join(current_seq))))

    return records


def read_fasta_file(filepath: str | Path) -> list[tuple[str, str]]:
    """Reads a FASTA file from disk and parses records."""
    path = Path(filepath)
    if not path.exists():
        raise FileNotFoundError(f"FASTA file not found: {filepath}")
    with open(path, "r", encoding="utf-8") as f:
        return parse_fasta(f.read())


def export_json(result: FullChainAnalysis | BatchAnalysisResult) -> str:
    """Exports Pydantic result model to formatted JSON string."""
    return result.model_dump_json(indent=2)


def export_csv(result: FullChainAnalysis | BatchAnalysisResult) -> str:
    """Exports analysis results to CSV format."""
    output = io.StringIO()
    writer = csv.writer(output)

    headers = [
        "Sequence_ID",
        "Chain_Type",
        "Isotype",
        "Subclass",
        "Allotypes",
        "Isoallotypes",
        "CDR1_Seq",
        "CDR2_Seq",
        "CDR3_Seq",
        "V_Domain_Seq",
        "C_Region_Seq",
    ]
    writer.writerow(headers)

    analyses: list[FullChainAnalysis] = []
    if isinstance(result, BatchAnalysisResult):
        analyses = result.results
    else:
        analyses = [result]

    for item in analyses:
        seq_id = item.sequence_id
        chain_type = item.v_analysis.chain_type.value if item.v_analysis else "Unknown"
        isotype = item.c_analysis.isotype if item.c_analysis else "Unknown"
        subclass = item.c_analysis.subclass if item.c_analysis else "Unknown"

        present_allotypes = [
            a.allotype
            for a in (item.c_analysis.allotypes if item.c_analysis else [])
            if a.status == "present"
        ]
        present_isoallotypes = [
            a.allotype
            for a in (item.c_analysis.isoallotypes if item.c_analysis else [])
            if a.status == "present"
        ]

        allotype_str = ";".join(present_allotypes) if present_allotypes else "None"
        isoallotype_str = (
            ";".join(present_isoallotypes) if present_isoallotypes else "None"
        )

        cdr1 = (
            item.v_analysis.cdrs.get("CDR1").sequence
            if (item.v_analysis and "CDR1" in item.v_analysis.cdrs)
            else ""
        )
        cdr2 = (
            item.v_analysis.cdrs.get("CDR2").sequence
            if (item.v_analysis and "CDR2" in item.v_analysis.cdrs)
            else ""
        )
        cdr3 = (
            item.v_analysis.cdrs.get("CDR3").sequence
            if (item.v_analysis and "CDR3" in item.v_analysis.cdrs)
            else ""
        )
        v_seq = item.v_analysis.v_domain_sequence if item.v_analysis else ""
        c_seq = item.c_analysis.c_region_sequence if item.c_analysis else ""

        row = [
            seq_id,
            chain_type,
            isotype,
            subclass,
            allotype_str,
            isoallotype_str,
            cdr1,
            cdr2,
            cdr3,
            v_seq,
            c_seq,
        ]
        writer.writerow(row)

    return output.getvalue()
