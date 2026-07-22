from ab_numbering_lite.schemas import FullChainAnalysis, BatchAnalysisResult
from ab_numbering_lite.vdomain import analyze_vdomain
from ab_numbering_lite.cdomain import analyze_cdomain


def analyze_chain(
    sequence: str,
    sequence_id: str = "query",
    scheme: str = "imgt",
) -> FullChainAnalysis:
    """
    Performs full integrated analysis of an antibody sequence:
    1. V-domain numbering & CDR1/2/3 extraction via ANARCII (with alignment fallback).
    2. Constant region extraction & Isotype / Subclass / Isoallotype / Allotype determination.
    """
    v_res = analyze_vdomain(sequence, scheme=scheme)
    v_domain_len = len(v_res.v_domain_sequence)
    c_res = analyze_cdomain(
        sequence, v_domain_len=v_domain_len, chain_hint=v_res.chain_type
    )

    return FullChainAnalysis(
        sequence_id=sequence_id,
        full_sequence=sequence,
        v_analysis=v_res,
        c_analysis=c_res,
        success=True,
    )


__all__ = [
    "analyze_chain",
    "analyze_vdomain",
    "analyze_cdomain",
    "FullChainAnalysis",
    "BatchAnalysisResult",
]
