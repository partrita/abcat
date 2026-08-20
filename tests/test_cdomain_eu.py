from copy import deepcopy

from abcat.cdomain import (
    _call_allotypes,
    _load_databases,
    _make_aligner,
    analyze_cdomain,
)
from abcat.schemas import ChainType


def test_isotype_is_selected_before_subclass():
    c_gene_db, _ = _load_databases()
    igg1 = c_gene_db["heavy"]["IGHG1"]["sequence"]

    result = analyze_cdomain(igg1, chain_hint=ChainType.HEAVY)

    assert result.isotype == "IgG"
    assert result.subclass == "IgG1"
    assert result.matched_c_gene == "IGHG1"


def test_allotype_matching_uses_eu_position_not_legacy_offset():
    c_gene_db, allotype_db = _load_databases()
    reference = c_gene_db["heavy"]["IGHG1"]["sequence"]
    markers = deepcopy(allotype_db["heavy"]["IgG1"])

    # EU position is canonical; legacy subclass offsets should not be required.
    for rule in markers:
        for marker in rule.get("markers", []):
            marker.pop("subclass_offset", None)

    alignment = _make_aligner().align(reference, reference)[0]
    allotypes, isoallotypes = _call_allotypes(
        query_c_seq=reference,
        matched_subclass="IgG1",
        chain_category="heavy",
        allotype_db={"heavy": {"IgG1": markers}},
        alignment=alignment,
        reference_c_seq=reference,
    )

    present = {a.allotype for a in allotypes if a.status == "present"}
    assert "G1m17" in present
    g1m17 = next(a for a in allotypes if a.allotype == "G1m17")
    assert g1m17.matched_residues["EU_214"] == "K"
    assert not any(a.allotype == "G1m17" for a in isoallotypes)
