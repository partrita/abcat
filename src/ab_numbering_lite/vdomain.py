import logging
from typing import Dict, Tuple, Any
from Bio.Align import PairwiseAligner

from ab_numbering_lite.schemas import (
    VAnalysisResult,
    ChainType,
    CDRRegion,
    FrameworkRegion,
)

logger = logging.getLogger(__name__)

SCHEME_BOUNDARIES = {
    "imgt": {
        "FR1": (1, 26),
        "CDR1": (27, 38),
        "FR2": (39, 55),
        "CDR2": (56, 65),
        "FR3": (66, 104),
        "CDR3": (105, 117),
        "FR4": (118, 128),
    },
    "kabat": {
        "FR1": (1, 30),
        "CDR1": (31, 35),
        "FR2": (36, 49),
        "CDR2": (50, 65),
        "FR3": (66, 94),
        "CDR3": (95, 102),
        "FR4": (103, 113),
    },
    "chothia": {
        "FR1": (1, 25),
        "CDR1": (26, 32),
        "FR2": (33, 51),
        "CDR2": (52, 56),
        "FR3": (57, 94),
        "CDR3": (95, 102),
        "FR4": (103, 113),
    },
    "martin": {
        "FR1": (1, 25),
        "CDR1": (26, 35),
        "FR2": (36, 49),
        "CDR2": (50, 58),
        "FR3": (59, 94),
        "CDR3": (95, 102),
        "FR4": (103, 113),
    },
    "aho": {
        "FR1": (1, 26),
        "CDR1": (27, 40),
        "FR2": (41, 57),
        "CDR2": (58, 68),
        "FR3": (69, 106),
        "CDR3": (107, 138),
        "FR4": (139, 150),
    },
}

HUMAN_CONSENSUS_V = {
    ChainType.HEAVY: "EVQLVESGGGLVQPGGSLRLSCAASGFTFSDHYMDWVRQAPGKGLEWVGRIRSKANSYATAYAASVKGRFTISRDDSKNTLYLQMNSLRAEDTAVYYCARFDAYWGQGTLVTVSS",
    ChainType.KAPPA: "DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK",
    ChainType.LAMBDA: "QSALTQPASVSGSPGQSITISCTGTSSDVGGYNYVSWYQQHPGKAPKLMIYEVSNRPSGVSNRFSGSKSGNTASLTISGLQAEDEADYYCSSYTSSSTWVFGGGTKLTVL",
}

VALID_SCHEMES = ["imgt", "kabat", "martin", "chothia", "aho"]


def analyze_vdomain(sequence: str, scheme: str = "imgt") -> VAnalysisResult:
    """
    Analyzes the variable domain of an antibody sequence using ANARCII (with fallback alignment).
    Extracts numbering (imgt, kabat, martin, chothia, aho), chain type, CDR1/2/3, and FR1/2/3/4 regions.
    """
    cleaned_seq = sequence.strip().upper().replace(" ", "").replace("\n", "")
    scheme_lower = scheme.lower()
    if scheme_lower not in VALID_SCHEMES:
        logger.warning(
            f"Unknown scheme '{scheme}'. Supported schemes: {VALID_SCHEMES}. Defaulting to 'imgt'."
        )
        scheme_lower = "imgt"

    try:
        import anarcii  # type: ignore

        # Use ANARCII Class API if present
        if hasattr(anarcii, "Anarcii"):
            model = anarcii.Anarcii(seq_type="antibody")
            results = model.number([cleaned_seq])
            if scheme_lower != "imgt" and hasattr(model, "to_scheme"):
                try:
                    converted = model.to_scheme(scheme_lower)
                    if converted:
                        results = converted
                except Exception as scheme_err:
                    logger.warning(
                        f"Failed to convert to scheme '{scheme_lower}': {scheme_err}"
                    )

            if results:
                return _parse_anarcii_results(cleaned_seq, results, scheme_lower)

        # Function interface fallback
        results = None
        if hasattr(anarcii, "run_anarcii"):
            results = anarcii.run_anarcii([("query", cleaned_seq)], scheme=scheme_lower)
        elif hasattr(anarcii, "anarcii"):
            results = anarcii.anarcii([("query", cleaned_seq)], scheme=scheme_lower)

        if results:
            return _parse_anarcii_results(cleaned_seq, results, scheme_lower)
    except Exception as e:
        logger.debug(
            f"ANARCII API call failed or not installed ({e}), using built-in alignment fallback."
        )

    return _fallback_vdomain_analysis(cleaned_seq, scheme_lower)


def _parse_anarcii_results(full_seq: str, results: Any, scheme: str) -> VAnalysisResult:
    """Parses standard ANARCII python output tuple/list structure."""
    try:
        domain_numbering = None
        chain_letter = None

        if (
            isinstance(results, tuple)
            and len(results) == 3
            and isinstance(results[1], list)
        ):
            numbered_seqs, alignment_details, _ = results
            if numbered_seqs and numbered_seqs[0] and numbered_seqs[0][0]:
                domain_numbering = numbered_seqs[0][0][0]
                if (
                    alignment_details
                    and alignment_details[0]
                    and alignment_details[0][0]
                ):
                    chain_letter = alignment_details[0][0].get("chain_type", "H")
        elif isinstance(results, list) and len(results) > 0:
            first_item = results[0]
            if isinstance(first_item, tuple) and len(first_item) == 3:
                domain_numbering = first_item[0]
            elif isinstance(first_item, list):
                domain_numbering = first_item

        if not domain_numbering:
            return _fallback_vdomain_analysis(full_seq, scheme)

        numbering_dict: Dict[str, str] = {}
        v_seq_chars = []

        if isinstance(domain_numbering, (list, tuple)):
            for item in domain_numbering:
                if isinstance(item, (list, tuple)) and len(item) == 2:
                    pos_info, aa = item
                    if isinstance(pos_info, (list, tuple)) and len(pos_info) == 2:
                        pos_num, ins_code = pos_info
                        if aa != "-":
                            pos_str = f"{pos_num}{str(ins_code).strip()}"
                            numbering_dict[pos_str] = aa
                            v_seq_chars.append(aa)

        if not numbering_dict:
            return _fallback_vdomain_analysis(full_seq, scheme)

        v_domain_seq = "".join(v_seq_chars)

        if not chain_letter:
            import re

            if "118" in numbering_dict and numbering_dict["118"] in ["W", "w"]:
                chain_type = ChainType.HEAVY
            elif re.search(r"W[GQAE]G[A-Z]TLTVSS|WG[A-Z]G", v_domain_seq):
                chain_type = ChainType.HEAVY
            else:
                chain_type = ChainType.KAPPA
        else:
            chain_map = {
                "H": ChainType.HEAVY,
                "K": ChainType.KAPPA,
                "L": ChainType.LAMBDA,
            }
            chain_type = chain_map.get(chain_letter, ChainType.HEAVY)

        cdrs, frameworks = _extract_regions_from_numbering(
            numbering_dict, scheme=scheme
        )

        return VAnalysisResult(
            chain_type=chain_type,
            scheme=scheme,
            v_domain_sequence=v_domain_seq,
            numbering=numbering_dict,
            cdrs=cdrs,
            frameworks=frameworks,
            confidence=0.99,
        )
    except Exception as err:
        logger.debug(f"Error parsing ANARCII results ({err}), using fallback")
        return _fallback_vdomain_analysis(full_seq, scheme)


def _fallback_vdomain_analysis(sequence: str, scheme: str) -> VAnalysisResult:
    """Robust built-in landmark-anchored V-domain alignment fallback."""
    import re

    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 3
    aligner.mismatch_score = -1
    aligner.open_gap_score = -5
    aligner.extend_gap_score = -0.5
    aligner.end_insertion_score = 0
    aligner.end_deletion_score = 0

    best_chain = ChainType.UNKNOWN
    best_score = -1e9

    for chain, ref_seq in HUMAN_CONSENSUS_V.items():
        alignments = aligner.align(sequence, ref_seq)
        if alignments and alignments[0].score > best_score:
            best_score = alignments[0].score
            best_chain = chain

    if best_score < 50:
        return VAnalysisResult(
            chain_type=ChainType.UNKNOWN,
            scheme=scheme,
            v_domain_sequence="",
            confidence=0.0,
            error_message="Sequence does not match antibody V-domain profile",
        )

    # Determine FR4 J-motif end
    if best_chain == ChainType.HEAVY:
        j_match = re.search(r"W[GQAE]G[A-Z]TLTVSS|WG[A-Z]G", sequence)
    else:
        j_match = re.search(r"F[GQAE]G[A-Z]TKVEIK|FG[A-Z]G", sequence)

    v_end_idx = j_match.end() if j_match else min(len(sequence), 128)
    v_domain_seq = sequence[:v_end_idx]

    # Find key anchor residues: Cys23, Trp41, Cys104, Trp118
    cys_indices = [m.start() for m in re.finditer(r"C", v_domain_seq)]
    c23_idx = cys_indices[0] if len(cys_indices) > 0 else 22
    c104_idx = (
        cys_indices[1]
        if len(cys_indices) > 1
        else max(len(v_domain_seq) - 15, c23_idx + 60)
    )

    w_indices = [m.start() for m in re.finditer(r"W", v_domain_seq)]
    w41_idx = next((i for i in w_indices if c23_idx < i < c104_idx), c23_idx + 15)
    w118_idx = j_match.start() if j_match else len(v_domain_seq) - 11

    numbering_dict: Dict[str, str] = {}

    for i, aa in enumerate(v_domain_seq):
        if i <= c23_idx:
            pos = int(23 - (c23_idx - i))
        elif c23_idx < i < w41_idx:
            total_gap = max(1, w41_idx - c23_idx)
            pos = int(23 + round((i - c23_idx) / total_gap * 18))
        elif i == w41_idx:
            pos = 41
        elif w41_idx < i < c104_idx:
            total_gap = max(1, c104_idx - w41_idx)
            pos = int(41 + round((i - w41_idx) / total_gap * 63))
        elif i == c104_idx:
            pos = 104
        elif c104_idx < i < w118_idx:
            pos = int(105 + (i - c104_idx - 1))
        else:
            pos = int(118 + (i - w118_idx))

        numbering_dict[str(pos)] = aa

    cdrs, frameworks = _extract_regions_from_numbering(numbering_dict, scheme=scheme)

    return VAnalysisResult(
        chain_type=best_chain,
        scheme=scheme,
        v_domain_sequence=v_domain_seq,
        numbering=numbering_dict,
        cdrs=cdrs,
        frameworks=frameworks,
        confidence=round(min(1.0, max(0.5, best_score / 300.0)), 2),
    )


def _extract_regions_from_numbering(
    numbering: Dict[str, str],
    scheme: str = "imgt",
) -> Tuple[Dict[str, CDRRegion], Dict[str, FrameworkRegion]]:
    """Separates numbered dictionary into CDR1/2/3 and FR1/2/3/4 regions based on scheme boundaries."""
    cdrs: Dict[str, CDRRegion] = {}
    frameworks: Dict[str, FrameworkRegion] = {}

    scheme_lower = scheme.lower()
    boundaries = SCHEME_BOUNDARIES.get(scheme_lower, SCHEME_BOUNDARIES["imgt"])

    for region_name, (start_num, end_num) in boundaries.items():
        seq_chars = []
        for pos_str, aa in numbering.items():
            num_part = int("".join(filter(str.isdigit, pos_str)))
            if start_num <= num_part <= end_num:
                seq_chars.append(aa)

        region_seq = "".join(seq_chars)

        if region_name.startswith("CDR"):
            cdrs[region_name] = CDRRegion(
                name=region_name,
                sequence=region_seq,
                start_pos=str(start_num),
                end_pos=str(end_num),
                length=len(region_seq),
            )
        else:
            frameworks[region_name] = FrameworkRegion(
                name=region_name,
                sequence=region_seq,
            )

    return cdrs, frameworks
