from abmap.cdomain import analyze_cdomain
from abmap.schemas import ChainType

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
