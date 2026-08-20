import json
import logging
from pathlib import Path
from typing import Any

from Bio.Align import PairwiseAligner

from abcat.schemas import AllotypeMarkerCall, CAnalysisResult, ChainType
from abcat.utils import clean_sequence

logger = logging.getLogger(__name__)

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


def _make_aligner() -> PairwiseAligner:
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -5
    aligner.extend_gap_score = -1
    return aligner


def _isotype_candidates(
    chain_category: str,
    c_gene_db: dict[str, Any],
) -> dict[str, dict[str, tuple[str, dict[str, Any]]]]:
    """Group constant-gene references by isotype before subclass matching."""
    grouped: dict[str, dict[str, tuple[str, dict[str, Any]]]] = {}
    for gene_name, gene_info in c_gene_db.get(chain_category, {}).items():
        isotype = gene_info.get("isotype", "Unknown")
        grouped.setdefault(isotype, {})[gene_name] = (gene_name, gene_info)
    return grouped


def _classify_isotype(
    c_seq: str,
    chain_category: str,
    c_gene_db: dict[str, Any],
    aligner: PairwiseAligner,
) -> tuple[str, float]:
    """Determine isotype first, using the N-terminal constant-region signature."""
    grouped = _isotype_candidates(chain_category, c_gene_db)
    signature_len = min(70, len(c_seq))
    signature = c_seq[:signature_len]

    best_isotype = "Unknown"
    best_score = float("-inf")
    for isotype, genes in grouped.items():
        iso_score = float("-inf")
        for _, gene_info in genes.values():
            ref_signature = gene_info["sequence"][:signature_len]
            alignments = aligner.align(signature, ref_signature)
            if alignments:
                iso_score = max(iso_score, alignments[0].score)
        if iso_score > best_score:
            best_score = iso_score
            best_isotype = isotype

    return best_isotype, best_score


def _select_subclass(
    c_seq: str,
    chain_category: str,
    isotype: str,
    c_gene_db: dict[str, Any],
    aligner: PairwiseAligner,
) -> tuple[str | None, dict[str, Any] | None, float, Any | None]:
    """Match the constant sequence only against genes belonging to the chosen isotype."""
    best_gene: str | None = None
    best_info: dict[str, Any] | None = None
    best_score = float("-inf")
    best_alignment: Any | None = None

    for gene_name, gene_info in c_gene_db.get(chain_category, {}).items():
        if gene_info.get("isotype") != isotype:
            continue
        alignments = aligner.align(c_seq, gene_info["sequence"])
        if alignments and alignments[0].score > best_score:
            best_score = alignments[0].score
            best_gene = gene_name
            best_info = gene_info
            best_alignment = alignments[0]

    return best_gene, best_info, best_score, best_alignment


def _infer_eu_start(reference_c_seq: str, markers: list[dict[str, Any]]) -> int | None:
    """Infer the EU numbering origin for the reference constant sequence."""
    del reference_c_seq  # Reserved for future explicit EU-numbered reference tables.
    starts = {
        int(m["subclass_offset"]) - int(m["eu_position"]) + 1
        for rule in markers
        for m in rule.get("markers", [])
        if "subclass_offset" in m and "eu_position" in m
    }
    if len(starts) == 1:
        return starts.pop()
    if not starts:
        logger.warning("No EU origin metadata available for allotype markers")
    else:
        logger.warning("Inconsistent EU origins in allotype markers: %s", sorted(starts))
    return None


def _get_query_index_for_eu_position(
    alignment: Any,
    eu_position: int,
    reference_eu_start: int | None,
) -> int | None:
    """Map a canonical EU position through reference and query pairwise alignment."""
    if reference_eu_start is None:
        return None
    ref_idx = eu_position - reference_eu_start
    if ref_idx < 0:
        return None
    return _get_query_index_for_ref_offset(alignment, ref_idx, 0)


def _get_query_index_for_ref_offset(
    alignment: Any, ref_idx: int, default_len: int
) -> int | None:
    """Map a 0-indexed reference position to query position using aligned blocks."""
    if hasattr(alignment, "aligned") and len(alignment.aligned) >= 2:
        query_blocks, ref_blocks = alignment.aligned[0], alignment.aligned[1]
        for (q_start, q_end), (r_start, r_end) in zip(query_blocks, ref_blocks):
            if r_start <= ref_idx < r_end:
                return q_start + (ref_idx - r_start)
        return None
    if 0 <= ref_idx < default_len:
        return ref_idx
    return None


def analyze_cdomain(
    sequence: str,
    v_domain_len: int = 0,
    chain_hint: ChainType | None = None,
) -> CAnalysisResult:
    """Analyze the constant domain in order: isotype -> subclass -> EU allotype."""
    cleaned = clean_sequence(sequence)
    c_seq = cleaned[v_domain_len:] if v_domain_len > 0 and len(cleaned) > v_domain_len else cleaned

    c_gene_db, allotype_db = _load_databases()
    aligner = _make_aligner()

    if chain_hint == ChainType.HEAVY:
        search_groups = ["heavy"]
    elif chain_hint in (ChainType.KAPPA, ChainType.LAMBDA):
        search_groups = ["light"]
    else:
        search_groups = ["heavy", "light"]

    # Stage 1: determine isotype without letting allotype/subclass differences influence the class call.
    best_chain: str | None = None
    best_isotype = "Unknown"
    best_iso_score = float("-inf")
    for category in search_groups:
        isotype, score = _classify_isotype(c_seq, category, c_gene_db, aligner)
        if score > best_iso_score:
            best_iso_score = score
            best_isotype = isotype
            best_chain = category

    if best_chain is None or best_isotype == "Unknown":
        return CAnalysisResult(
            c_region_sequence=c_seq,
            isotype="Unknown",
            subclass="Unknown",
            matched_c_gene="None",
            alignment_identity=0.0,
            allotypes=[],
            isoallotypes=[],
        )

    # Stage 2: subclass is matched only inside the selected isotype.
    best_match_key, best_c_info, best_score, best_alignment = _select_subclass(
        c_seq, best_chain, best_isotype, c_gene_db, aligner
    )
    if not best_match_key or not best_c_info or best_alignment is None:
        return CAnalysisResult(
            c_region_sequence=c_seq,
            isotype=best_isotype,
            subclass="Unknown",
            matched_c_gene="None",
            alignment_identity=0.0,
            allotypes=[],
            isoallotypes=[],
        )

    ref_len = len(best_c_info["sequence"])
    identity = round(min(1.0, max(0.0, best_score / (2.0 * ref_len))), 4)
    subclass = best_c_info["subclass"]

    # Stage 3: determine allotypes from canonical EU-numbered positions.
    allotype_calls, isoallotype_calls = _call_allotypes(
        query_c_seq=c_seq,
        matched_subclass=subclass,
        chain_category=best_chain,
        allotype_db=allotype_db,
        alignment=best_alignment,
        reference_c_seq=best_c_info["sequence"],
    )

    present_allotypes = [a.allotype for a in allotype_calls if a.status == "present"]
    summary_str = format_allotype_summary(present_allotypes)

    return CAnalysisResult(
        c_region_sequence=c_seq,
        isotype=best_isotype,
        subclass=subclass,
        matched_c_gene=best_match_key,
        alignment_identity=identity,
        allotypes=allotype_calls,
        isoallotypes=isoallotype_calls,
        allotype_summary=summary_str,
    )


def format_allotype_summary(present_allotypes: list[str]) -> str:
    """Formats a list of present allotype names into IMGT notation (e.g. G1m17,1)."""
    if not present_allotypes:
        return "None"

    import re

    groups: dict[str, list[str]] = {}
    for allo in present_allotypes:
        match = re.match(r"^([A-Za-z]+[0-9]*m|\w+[\*\_]?)(.*)$", allo)
        if match:
            prefix, suffix = match.groups()
            groups.setdefault(prefix, []).append(suffix)
        else:
            groups.setdefault("other", []).append(allo)

    formatted_parts = []
    for prefix, suffixes in groups.items():
        if prefix == "other":
            formatted_parts.extend(suffixes)
        elif len(suffixes) == 1:
            formatted_parts.append(f"{prefix}{suffixes[0]}")
        else:
            formatted_parts.append(f"{prefix}{','.join(suffixes)}")

    return ", ".join(formatted_parts)


def _call_allotypes(
    query_c_seq: str,
    matched_subclass: str,
    chain_category: str,
    allotype_db: dict[str, Any],
    alignment: Any,
    reference_c_seq: str,
) -> tuple[list[AllotypeMarkerCall], list[AllotypeMarkerCall]]:
    """Call allotype/isoallotype markers using canonical EU positions."""
    allotypes: list[AllotypeMarkerCall] = []
    isoallotypes: list[AllotypeMarkerCall] = []

    subclass_markers = allotype_db.get(chain_category, {}).get(matched_subclass, [])
    if not subclass_markers:
        return allotypes, isoallotypes

    eu_start = _infer_eu_start(reference_c_seq, subclass_markers)
    if eu_start is None:
        return allotypes, isoallotypes

    for rule in subclass_markers:
        allotype_name = rule["allotype"]
        opposing = rule.get("opposing_allotype")
        domain = rule.get("domain", "C")
        category = rule.get("type", "allotype")
        markers = rule.get("markers", [])

        matches: dict[str, str] = {}
        all_matched = True
        any_observed = False

        for marker in markers:
            eu_pos = int(marker["eu_position"])
            expected_aa = marker["amino_acid"]
            query_idx = _get_query_index_for_eu_position(best_alignment := alignment, eu_pos, eu_start)
            if query_idx is None or not (0 <= query_idx < len(query_c_seq)):
                all_matched = False
                continue

            actual_aa = query_c_seq[query_idx]
            matches[f"EU_{eu_pos}"] = actual_aa
            any_observed = True
            if actual_aa != expected_aa:
                all_matched = False

        if all_matched and any_observed:
            status = "present"
        elif any_observed:
            status = "absent"
        else:
            status = "inconclusive"

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
