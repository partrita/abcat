from typer.testing import CliRunner

from ab_numbering_lite.cli import app

runner = CliRunner()

TRASTUZUMAB_VH = (
    "EVQLVESGGGLVQPGGSLRLSCAASGFTFTDYTMDWVRQAPGKGLEWVADVNPNSGGSIYNQRFKGRFTLS"
    "VDRSKNTLYLQMNSLRAEDTAVYYCARNLGPSFYFDYWGQGTLVTVSS"
)


def test_cli_analyze():
    result = runner.invoke(app, ["analyze", "--sequence", TRASTUZUMAB_VH])
    assert result.exit_code == 0
    assert "Heavy" in result.stdout or "Seq ID" in result.stdout


def test_cli_vdomain():
    result = runner.invoke(app, ["vdomain", "--sequence", TRASTUZUMAB_VH])
    assert result.exit_code == 0
    assert "Chain Type" in result.stdout


def test_cli_cdomain():
    result = runner.invoke(
        app,
        [
            "cdomain",
            "--sequence",
            "ASTKGPSVFPLAPSSKSTSGGTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSSVVTVPSSSLGTQTYICNVNHKPSNTKVDKKVEPKSCDK",
        ],
    )
    assert result.exit_code == 0
    assert "Isotype" in result.stdout


def test_export_csv():
    from ab_numbering_lite import analyze_chain
    from ab_numbering_lite.utils import export_csv

    res = analyze_chain(TRASTUZUMAB_VH, sequence_id="test1")
    csv_out = export_csv(res)
    assert "Sequence_ID,Chain_Type,Isotype" in csv_out
    assert "test1,Heavy" in csv_out
