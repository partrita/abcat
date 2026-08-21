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

    c_gene_db, _, _ = _load_databases()
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


def test_cdomain_imgt_alleles_ighg1():
    """Tests IMGT allele calling for IGHG1 variants (*01, *03, *04, *07, *08, *11, *13)."""
    # IGHG1*01: 214K, 309L, 356D, 358L, 384N
    res_01 = analyze_cdomain(
        WILDTYPE_G1M1_17_HEAVY, v_domain_len=119, chain_hint=ChainType.HEAVY
    )
    assert res_01.imgt_allele is not None
    assert "IGHG1*01" in res_01.imgt_alleles

    # IGHG1*03: 199I, 214R, 309L, 356E, 358M, 384N
    # In WILDTYPE_G1M1_17_HEAVY: NTKVDKKVEPK -> NIKVDKRVEPK (199I, 214R), RDELTKN -> REEMTKN (356E, 358M)
    seq_03 = WILDTYPE_G1M1_17_HEAVY.replace("NTKVDKKVEPK", "NIKVDKRVEPK").replace(
        "RDELTKN", "REEMTKN"
    )
    res_03 = analyze_cdomain(seq_03, v_domain_len=119, chain_hint=ChainType.HEAVY)
    assert res_03.imgt_allele == "IGHG1*03"
    assert "IGHG1*03" in res_03.imgt_alleles

    # IGHG1*04: 214K, 309L, 356D, 358L, 384N, 422I
    # In WILDTYPE_G1M1_17_HEAVY: GNVFSCS -> GNIFSCS (422I)
    seq_04 = WILDTYPE_G1M1_17_HEAVY.replace("GNVFSCS", "GNIFSCS")
    res_04 = analyze_cdomain(seq_04, v_domain_len=119, chain_hint=ChainType.HEAVY)
    assert res_04.imgt_allele == "IGHG1*04"

    # IGHG1*07: 214K, 309L, 356D, 358L, 384N, 431G
    # In WILDTYPE_G1M1_17_HEAVY: HEALHNH -> HEGLHNH (431G)
    seq_07 = WILDTYPE_G1M1_17_HEAVY.replace("HEALHNH", "HEGLHNH")
    res_07 = analyze_cdomain(seq_07, v_domain_len=119, chain_hint=ChainType.HEAVY)
    assert res_07.imgt_allele == "IGHG1*07"

    # IGHG1*08: 199I, 214R, 309L, 356D, 358L, 384N
    # In WILDTYPE_G1M1_17_HEAVY: NTKVDKKVEPK -> NIKVDKRVEPK (199I, 214R), keeps 356D, 358L
    seq_08 = WILDTYPE_G1M1_17_HEAVY.replace("NTKVDKKVEPK", "NIKVDKRVEPK")
    res_08 = analyze_cdomain(seq_08, v_domain_len=119, chain_hint=ChainType.HEAVY)
    assert res_08.imgt_allele == "IGHG1*08"

    # IGHG1*11: 214K, 309V, 356D, 358L, 384N
    # In WILDTYPE_G1M1_17_HEAVY: LTVLHQ -> LTVVHQ (309V)
    seq_11 = WILDTYPE_G1M1_17_HEAVY.replace("LTVLHQ", "LTVVHQ")
    res_11 = analyze_cdomain(seq_11, v_domain_len=119, chain_hint=ChainType.HEAVY)
    assert res_11.imgt_allele == "IGHG1*11"

    # IGHG1*13: 214K, 296F, 309L, 356D, 358L, 384N
    # In WILDTYPE_G1M1_17_HEAVY: EQYNSTY -> EQFNSTY (296F)
    seq_13 = WILDTYPE_G1M1_17_HEAVY.replace("EQYNSTY", "EQFNSTY")
    res_13 = analyze_cdomain(seq_13, v_domain_len=119, chain_hint=ChainType.HEAVY)
    assert res_13.imgt_allele == "IGHG1*13"


def test_cdomain_imgt_alleles_ighg2_ighg3_ighg4_igkc():
    """Tests IMGT allele calling for IGHG2, IGHG3, IGHG4, and IGKC."""
    from abcat.cdomain import _load_databases

    c_gene_db, _, _ = _load_databases()

    # IGHG2
    g2_seq = c_gene_db["heavy"]["IGHG2"]["sequence"]
    res_g2_04 = analyze_cdomain(g2_seq, chain_hint=ChainType.HEAVY)
    assert "IGHG2*04" in res_g2_04.imgt_alleles

    # IGHG2*01: replace 192N/193F (SSNFGT) with 192S/193L (SSSLGT)
    g2_01_seq = g2_seq.replace("SSNFGT", "SSSLGT")
    res_g2_01 = analyze_cdomain(g2_01_seq, chain_hint=ChainType.HEAVY)
    assert "IGHG2*01" in res_g2_01.imgt_alleles

    # IGHG2*02 (282M): replace DGVEVH with DGMEVH
    g2_02_seq = g2_seq.replace("DGVEVH", "DGMEVH")
    res_g2_02 = analyze_cdomain(g2_02_seq, chain_hint=ChainType.HEAVY)
    assert "IGHG2*02" in res_g2_02.imgt_alleles

    # IGHG3
    g3_seq = c_gene_db["heavy"]["IGHG3"]["sequence"]
    res_g3 = analyze_cdomain(g3_seq, chain_hint=ChainType.HEAVY)
    assert "IGHG3*01" in res_g3.imgt_alleles

    # IGHG3*04 (291L, 435H, 436Y): replace TKPREEQYN with TKLREEQYN, and HNRFT with HNHYT
    g3_04_seq = g3_seq.replace("TKPREEQYN", "TKLREEQYN").replace("HNRFT", "HNHYT")
    res_g3_04 = analyze_cdomain(g3_04_seq, chain_hint=ChainType.HEAVY)
    assert "IGHG3*04" in res_g3_04.imgt_alleles

    # IGHG4
    g4_seq = c_gene_db["heavy"]["IGHG4"]["sequence"]
    res_g4_01 = analyze_cdomain(g4_seq, chain_hint=ChainType.HEAVY)
    assert "IGHG4*01" in res_g4_01.imgt_alleles

    # IGHG4*02 (309V): replace LTVLHQ with LTVVHQ
    g4_02_seq = g4_seq.replace("LTVLHQ", "LTVVHQ")
    res_g4_02 = analyze_cdomain(g4_02_seq, chain_hint=ChainType.HEAVY)
    assert "IGHG4*02" in res_g4_02.imgt_alleles

    # IGKC (Km3: 153A, 191V -> IGKC*01)
    kc_seq = c_gene_db["light"]["IGKC"]["sequence"]
    res_kc_01 = analyze_cdomain(kc_seq, chain_hint=ChainType.KAPPA)
    assert "IGKC*01" in res_kc_01.imgt_alleles
