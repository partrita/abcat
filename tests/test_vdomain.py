from ab_numbering_lite.vdomain import analyze_vdomain
from ab_numbering_lite.schemas import ChainType

TRASTUZUMAB_VH = (
    "EVQLVESGGGLVQPGGSLRLSCAASGFTFTDYTMDWVRQAPGKGLEWVADVNPNSGGSIYNQRFKGRFTLS"
    "VDRSKNTLYLQMNSLRAEDTAVYYCARNLGPSFYFDYWGQGTLVTVSS"
)

TRASTUZUMAB_VL = (
    "DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTDFTLT"
    "ISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK"
)


def test_vdomain_trastuzumab_vh():
    result = analyze_vdomain(TRASTUZUMAB_VH, scheme="imgt")
    assert result.chain_type == ChainType.HEAVY
    assert "CDR1" in result.cdrs
    assert "CDR2" in result.cdrs
    assert "CDR3" in result.cdrs
    assert result.cdrs["CDR3"].sequence != ""
    assert len(result.v_domain_sequence) > 100


def test_vdomain_trastuzumab_vl():
    result = analyze_vdomain(TRASTUZUMAB_VL, scheme="imgt")
    assert result.chain_type == ChainType.KAPPA
    assert "CDR1" in result.cdrs
    assert "CDR2" in result.cdrs
    assert "CDR3" in result.cdrs


def test_vdomain_non_antibody():
    result = analyze_vdomain("MKTLLILAVIMIFLLQGQAKSEVN", scheme="imgt")
    assert result.chain_type == ChainType.UNKNOWN or result.confidence < 0.5


def test_vdomain_alternate_schemes():
    for scheme in ["kabat", "martin", "chothia", "aho"]:
        res = analyze_vdomain(TRASTUZUMAB_VH, scheme=scheme)
        assert res.scheme == scheme
        assert res.chain_type == ChainType.HEAVY
