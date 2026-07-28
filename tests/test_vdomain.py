from abmap.schemas import ChainType
from abmap.vdomain import analyze_vdomain

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


def test_vdomain_multiline_sequence():
    multiline_seq = (
        "EVQLVESGGGLVQPGGSLRLSCAASGFTFT\n"
        "DYTMDWVRQAPGKGLEWVADVNPNSGGSIY\r\n"
        "\tNQRFKGRFTLSVDRSKNTLYLQMNSLRAED\n"
        "TAVYYCARNLGPSFYFDYWGQGTLVTVSS"
    )
    single_res = analyze_vdomain(TRASTUZUMAB_VH, scheme="imgt")
    multi_res = analyze_vdomain(multiline_seq, scheme="imgt")
    assert multi_res.chain_type == single_res.chain_type
    assert multi_res.v_domain_sequence == single_res.v_domain_sequence
    assert multi_res.cdrs["CDR3"].sequence == single_res.cdrs["CDR3"].sequence
