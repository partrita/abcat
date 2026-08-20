from abcat.cdomain import analyze_cdomain
from abcat.schemas import ChainType

TRASTUZUMAB_HEAVY_FULL = (
    "EVQLVESGGGLVQPGGSLRLSCAASGFTFTDYTMDWVRQAPGKGLEWVADVNPNSGGSIYNQRFKGRFTLSVDRSKNTLYLQ"
    "MNSLRAEDTAVYYCARNLGPSFYFDYWGQGTLVTVSSASTKGPSVFPLAPSSKSTSGGTAALGCLVKDYFPEPVTVSWNSGAL"
    "TSGVHTFPAVLQSSGLYSLSSVVTVPSSSLGTQTYICNVNHKPSNTKVDKKVEPKSCDKTHTCPPCPAPELLGGPSVFLFPPK"
    "PKDTLMISRTPEVTCVVVDVSHEDPEVKFNWYVDGVEVHNAKTKPREEQYNSTYRVVSVLTVLHQDWLNGKEYKCKVSNKALP"
    "APIEKTISKAKGQPREPQVYTLPPSREEMTKNQVSLTCLVKGFYPSDIAVEWESNGQPENNYKTTPPVLDSDGSFFLYSKLTV"
    "DKSRWQQGNVFSCSVMHEALHNHYTQKSLSLSPGK"
)

TRASTUZUMAB_LIGHT_FULL = (
    "DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTDFTLTISSLQPEDF"
    "ATYYCQQHYTTPPTFGQGTKVEIKRTVAAPSVFIFPPSDEQLKSGTASVVCLLNNFYPREAKVQWKVDNALQSGNSQESVTEQ"
    "DSKDSTYSLSSTLTLSKADYEKHKVYACEVTHQGLSSPVTKSFNRGEC"
)


# Wild-type natural G1m1,17 IgG1 heavy chain sequence (contains D356/L358 and K214)
WILDTYPE_G1M1_17_HEAVY = (
    "EVQLVESGGGLVQPGGSLRLSCAASGFTFTDYTMDWVRQAPGKGLEWVADVNPNSGGSIYNQRFKGRFTLSVDRSKNTLYLQ"
    "MNSLRAEDTAVYYCARNLGPSFYFDYWGQGTLVTVSSASTKGPSVFPLAPSSKSTSGGTAALGCLVKDYFPEPVTVSWNSGAL"
    "TSGVHTFPAVLQSSGLYSLSSVVTVPSSSLGTQTYICNVNHKPSNTKVDKKVEPKSCDKTHTCPPCPAPELLGGPSVFLFPPK"
    "PKDTLMISRTPEVTCVVVDVSHEDPEVKFNWYVDGVEVHNAKTKPREEQYNSTYRVVSVLTVLHQDWLNGKEYKCKVSNKALP"
    "APIEKTISKAKGQPREPQVYTLPPSRDELTKNQVSLTCLVKGFYPSDIAVEWESNGQPENNYKTTPPVLDSDGSFFLYSKLTV"
    "DKSRWQQGNVFSCSVMHEALHNHYTQKSLSLSPGK"
)


def test_cdomain_trastuzumab_heavy():
    res = analyze_cdomain(
        TRASTUZUMAB_HEAVY_FULL, v_domain_len=120, chain_hint=ChainType.HEAVY
    )
    assert res.isotype == "IgG"
    assert res.subclass == "IgG1"
    assert res.matched_c_gene == "IGHG1"
    assert res.alignment_identity > 0.9

    allotype_names = [a.allotype for a in res.allotypes if a.status == "present"]
    isoallotype_names = [i.allotype for i in res.isoallotypes if i.status == "present"]

    assert "G1m17" in allotype_names
    assert "G1m1" not in allotype_names
    assert "nG1m1" in isoallotype_names


def test_cdomain_wildtype_g1m1_17_heavy():
    """Tests natural wild-type IgG1 heavy chain possessing both G1m1 (D356/L358) and G1m17 (K214)."""
    # Replace E356/M358 (REEMTKN) in Trastuzumab with natural D356/L358 (RDELTKN)
    wt_seq = TRASTUZUMAB_HEAVY_FULL.replace("REEMTKN", "RDELTKN")
    res = analyze_cdomain(wt_seq, v_domain_len=120, chain_hint=ChainType.HEAVY)
    assert res.isotype == "IgG"
    assert res.subclass == "IgG1"

    allotype_names = [a.allotype for a in res.allotypes if a.status == "present"]
    assert "G1m17" in allotype_names
    assert "G1m1" in allotype_names


def test_cdomain_trastuzumab_light():
    res = analyze_cdomain(
        TRASTUZUMAB_LIGHT_FULL, v_domain_len=108, chain_hint=ChainType.KAPPA
    )
    assert res.isotype == "Kappa"
    assert res.subclass == "IGKC"
    assert res.matched_c_gene == "IGKC"
    assert res.alignment_identity > 0.9


def test_cdomain_partial_marker_match_rejected():
    """Tests that partial matching of allotype markers does not call the allotype."""
    # Trastuzumab heavy has E356 and M358 (nG1m1)
    # If sequence has D356 and M358 (partial match for G1m1), G1m1 must NOT be present
    partial_seq = TRASTUZUMAB_HEAVY_FULL.replace("REEMTKN", "RDEMTKN")
    res = analyze_cdomain(partial_seq, v_domain_len=120, chain_hint=ChainType.HEAVY)
    allotype_names = [a.allotype for a in res.allotypes if a.status == "present"]
    assert "G1m1" not in allotype_names
    assert "nG1m1" not in allotype_names


def test_cdomain_g1m27_and_g1m28():
    """Tests IgG1 G1m27 (I422) and G1m28 (R435, Y436) allotypes."""
    # Modify Trastuzumab heavy sequence to introduce I422 (G1m27) and R435/Y436 (G1m28)
    # Trastuzumab has V422 (GNVFSCS) -> I422 (GNIFSCS) and H435/Y436 (HYTQK) -> R435/Y436 (RYTQK)
    seq = TRASTUZUMAB_HEAVY_FULL.replace("GNVFSCS", "GNIFSCS").replace("HYTQK", "RYTQK")
    res = analyze_cdomain(seq, v_domain_len=120, chain_hint=ChainType.HEAVY)
    allotype_names = [a.allotype for a in res.allotypes if a.status == "present"]
    assert "G1m27" in allotype_names
    assert "G1m28" in allotype_names


def test_cdomain_ige_allotypes():
    """Tests IgE constant region allotype calling for IGHE alleles."""
    from abcat.cdomain import _load_databases

    c_gene_db, _ = _load_databases()
    ige_seq = c_gene_db["heavy"]["IGHE"]["sequence"]

    res = analyze_cdomain(ige_seq, chain_hint=ChainType.HEAVY)
    assert res.isotype == "IgE"
    assert res.subclass == "IgE"

    allotype_names = [a.allotype for a in res.allotypes if a.status == "present"]
    assert "IGHE*02" in allotype_names or "IGHE*04" in allotype_names

    # Test IGHE*03 variant by replacing W149 with L149 (INITWLED -> INITLLED) and W43 with C43 (PEPVMVTW -> PEPVMVTC)
    ige_03_seq = ige_seq.replace("PEPVMVTW", "PEPVMVTC").replace("INITWLED", "INITLLED")
    res_03 = analyze_cdomain(ige_03_seq, chain_hint=ChainType.HEAVY)
    allotype_names_03 = [a.allotype for a in res_03.allotypes if a.status == "present"]
    assert "IGHE*03" in allotype_names_03


def test_cdomain_multiple_simultaneous_allotypes():
    """Tests simultaneous calling of multiple allotypes (e.g., G1m17 and G1m1 -> G1m17,1)."""
    # Trastuzumab heavy chain has G1m17 (K214) and nG1m1 (E356/M358)
    # Restore natural G1m1 (D356/L358): REEMTKN -> RDELTKN
    seq = TRASTUZUMAB_HEAVY_FULL.replace("REEMTKN", "RDELTKN")
    res = analyze_cdomain(seq, v_domain_len=120, chain_hint=ChainType.HEAVY)

    present_allotypes = [a.allotype for a in res.allotypes if a.status == "present"]
    assert "G1m17" in present_allotypes
    assert "G1m1" in present_allotypes
    assert res.allotype_summary == "G1m17,1"


def test_example_fasta_igg1_human():
    """Tests that IGG1_human sequence from examples fasta yields Isotype: IgG, Subclass: 1 (IgG1), and Allotype: G1m17,1."""
    from pathlib import Path

    from abcat import analyze_chain
    from abcat.utils import read_fasta_file

    fasta_path = Path(__file__).parent.parent / "examples" / "full_antibodies.fasta"
    records = dict(read_fasta_file(fasta_path))

    # Match IGG1_human sequence (case-insensitive header match)
    matching_key = next((k for k in records if "igg1_human" in k.lower()), None)
    assert matching_key is not None, (
        "IGG1_human sequence not found in examples/full_antibodies.fasta"
    )

    seq = records[matching_key]
    res = analyze_chain(seq, sequence_id=matching_key)

    assert res.c_analysis is not None
    assert res.c_analysis.isotype == "IgG"
    assert res.c_analysis.subclass == "IgG1" or res.c_analysis.subclass == "1"
    assert "1" in res.c_analysis.subclass
    assert res.c_analysis.allotype_summary == "G1m17,1"

    present_allotypes = [
        a.allotype for a in res.c_analysis.allotypes if a.status == "present"
    ]
    assert "G1m17" in present_allotypes
    assert "G1m1" in present_allotypes
