import json
import logging
from pathlib import Path
from typing import Any

from Bio.Align import PairwiseAligner

from abmap.schemas import AllotypeMarkerCall, CAnalysisResult, ChainType
from abmap.utils import clean_sequence

logger = logging.getLogger(__name__)

# Data directory path
DATA_DIR = Path(__file__).parent.parent.parent / "data"

_C_GENE_DB: dict[str, Any] | None = None
_ALLOTYPE_DB: dict[str, Any] | None = None


def _load_databases() -> tuple[dict[str, Any], dict[str, Any]]:
    global _C_GENE_DB, _ALLOTYPE_DB
    if _C_GENE_DB is None or _ALLOTYPE_DB is None:
        c_gene_path = DATA_DIR / "human_c_genes.json"
        allotype_path = DATA_DIR / "allotype_markers.json"

        if c_gene_path.exists():
            with open(c_gene_path, "r", encoding="utf-8") as f:
                _C_GENE_DB = json.load(f)
        else:
            _C_GENE_DB = {"heavy": {}, "light": {}}

        if allotype_path.exists():
            with open(allotype_path, "r", encoding="utf-8") as f:
                _ALLOTYPE_DB = json.load(f)
        else:
            _ALLOTYPE_DB = {"heavy": {}, "light": {}}

    return _C_GENE_DB, _ALLOTYPE_DB


def analyze_cdomain(
    sequence: str,
    v_domain_len: int = 0,
    chain_hint: ChainType | None = None,
) -> CAnalysisResult:
    """
    Analyzes the Constant Domain of an antibody sequence.
    Determines Isotype, Subclass, Isoallotype, and Allotype markers.
    """
    cleaned = clean_sequence(sequence)

    # Extract constant domain sequence if full sequence is provided
    if v_domain_len > 0 and len(cleaned) > v_domain_len:
        c_seq = cleaned[v_domain_len:]
    else:
        c_seq = cleaned

    c_gene_db, allotype_db = _load_databases()
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -5
    aligner.extend_gap_score = -1

    best_match_key = None
    best_c_info = None
    best_score = -1e9
    best_alignment = None
    chain_category = "heavy"

    # Filter targets based on chain hint if provided
    search_groups = []
    if chain_hint == ChainType.HEAVY:
        search_groups = [("heavy", c_gene_db.get("heavy", {}))]
    elif chain_hint in (ChainType.KAPPA, ChainType.LAMBDA):
        search_groups = [("light", c_gene_db.get("light", {}))]
    else:
        search_groups = [
            ("heavy", c_gene_db.get("heavy", {})),
            ("light", c_gene_db.get("light", {})),
        ]

    for cat_name, genes in search_groups:
        for gene_name, gene_info in genes.items():
            ref_seq = gene_info["sequence"]
            alignments = aligner.align(c_seq, ref_seq)
            if alignments:
                top = alignments[0]
                if top.score > best_score:
                    best_score = top.score
                    best_match_key = gene_name
                    best_c_info = gene_info
                    best_alignment = top
                    chain_category = cat_name

    if not best_match_key or not best_c_info or not best_alignment:
        return CAnalysisResult(
            c_region_sequence=c_seq,
            isotype="Unknown",
            subclass="Unknown",
            matched_c_gene="None",
            alignment_identity=0.0,
            allotypes=[],
            isoallotypes=[],
        )

    ref_len = len(best_c_info["sequence"])
    identity = round(min(1.0, max(0.0, best_score / (2.0 * ref_len))), 4)

    isotype = best_c_info["isotype"]
    subclass = best_c_info["subclass"]

    # Call Allotypes and Isoallotypes
    allotype_calls, isoallotype_calls = _call_allotypes(
        query_c_seq=c_seq,
        matched_subclass=subclass,
        chain_category=chain_category,
        allotype_db=allotype_db,
        alignment=best_alignment,
    )

    return CAnalysisResult(
        c_region_sequence=c_seq,
        isotype=isotype,
        subclass=subclass,
        matched_c_gene=best_match_key,
        alignment_identity=identity,
        allotypes=allotype_calls,
        isoallotypes=isoallotype_calls,
    )


def _get_query_index_for_ref_offset(
    alignment: Any, ref_idx: int, default_len: int
) -> int | None:
    """
    Maps 0-indexed reference C-gene position (ref_idx) to 0-indexed position in query_c_seq using PairwiseAligner alignment.
    alignment.aligned[0] is query blocks, alignment.aligned[1] is target (ref) blocks.
    """
    if hasattr(alignment, "aligned") and len(alignment.aligned) >= 2:
        query_blocks, ref_blocks = alignment.aligned[0], alignment.aligned[1]
        for (q_start, q_end), (r_start, r_end) in zip(query_blocks, ref_blocks):
            if r_start <= ref_idx < r_end:
                return q_start + (ref_idx - r_start)
        return None
    if 0 <= ref_idx < default_len:
        return ref_idx
    return None


def _call_allotypes(
    query_c_seq: str,
    matched_subclass: str,
    chain_category: str,
    allotype_db: dict[str, Any],
    alignment: Any,
) -> tuple[list[AllotypeMarkerCall], list[AllotypeMarkerCall]]:
    """Scans key polymorphic position markers to assign Allotypes and Isoallotypes."""
    allotypes: list[AllotypeMarkerCall] = []
    isoallotypes: list[AllotypeMarkerCall] = []

    subclass_markers = allotype_db.get(chain_category, {}).get(matched_subclass, [])
    if not subclass_markers:
        return allotypes, isoallotypes

    # Align query C-region to reference
    for rule in subclass_markers:
        allotype_name = rule["allotype"]
        opposing = rule.get("opposing_allotype")
        domain = rule.get("domain", "C")
        category = rule.get("type", "allotype")
        markers = rule.get("markers", [])

        matches = {}
        all_matched = True

        for m in markers:
            ref_offset = (
                m["subclass_offset"] - 1
            )  # 0-indexed offset in reference C-gene
            expected_aa = m["amino_acid"]
            eu_pos = m["eu_position"]

            query_idx = _get_query_index_for_ref_offset(
                alignment, ref_offset, len(query_c_seq)
            )

            if query_idx is not None and 0 <= query_idx < len(query_c_seq):
                actual_aa = query_c_seq[query_idx]
                matches[f"EU_{eu_pos}"] = actual_aa
                if actual_aa != expected_aa:
                    all_matched = False
            else:
                all_matched = False

        status = "present" if all_matched else "absent"

        call = AllotypeMarkerCall(
            allotype=allotype_name,
            opposing_allotype=opposing,
            category=category,
            domain=domain,
            status=status,
            matched_residues=matches,
        )

        if category == "isoallotype":
            isoallotypes.append(call)
        else:
            allotypes.append(call)

    return allotypes, isoallotypes
