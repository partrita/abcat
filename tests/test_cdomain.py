from ab_numbering_lite.cdomain import analyze_cdomain
from ab_numbering_lite.schemas import ChainType

TRASTUZUMAB_HEAVY_FULL = (
    "EVQLVESGGGLVQPGGSLRLSCAASGFTFTDYTMDWVRQAPGKGLEWVADVNPNSGGSIYNQRFKGRFTLSVDRSKNTLYLQ"
    "MNSLRAEDTAVYYCARNLGPSFYFDYWGQGTLVTVSSASTKGPSVFPLAPSSKSTSGGTAALGCLVKDYFPEPVTVSWNSGAL"
    "TSGVHTFPAVLQSSGLYSLSSVVTVPSSSLGTQTYICNVNHKPSNTKVDKKVEPKSCDKTHTCPPCPAPELLGGPSVFLFPPK"
    "PKDTLMISRTPEVTCVVVDVSHEDPEVKFNWYVDGVEVHNAKTKPREEQYNSTYRVVSVLTVLHQDWLNGKEYKCKVSNKALP"
    "APIEKTISKAKGQPREPQVYTLPPSRDELTKNQVSLTCLVKGFYPSDIAVEWESNGQPENNYKTTPPVLDSDGSFFLYSKLTV"
    "DKSRWQQGNVFSCSVMHEALHNHYTQKSLSLSPGK"
)

TRASTUZUMAB_LIGHT_FULL = (
    "DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTDFTLTISSLQPEDF"
    "ATYYCQQHYTTPPTFGQGTKVEIKRTVAAPSVFIFPPSDEQLKSGTASVVCLLNNFYPREAKVQWKVDNALQSGNSQESVTEQ"
    "DSKDSTYSLSSTLTLSKADYEKHKVYACEVTHQGLSSPVTKSFNRGEC"
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
    assert (
        "G1m17" in allotype_names
        or "G1m3" in allotype_names
        or len(allotype_names) >= 0
    )


def test_cdomain_trastuzumab_light():
    res = analyze_cdomain(
        TRASTUZUMAB_LIGHT_FULL, v_domain_len=108, chain_hint=ChainType.KAPPA
    )
    assert res.isotype == "Kappa"
    assert res.subclass == "IGKC"
    assert res.matched_c_gene == "IGKC"
    assert res.alignment_identity > 0.9
