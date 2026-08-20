from typer.testing import CliRunner

from abcat.cli import app

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
    from abcat import analyze_chain
    from abcat.utils import export_csv

    res = analyze_chain(TRASTUZUMAB_VH, sequence_id="test1")
    csv_out = export_csv(res)
    assert "Sequence_ID,Chain_Type,Isotype" in csv_out
    assert "test1,Heavy" in csv_out


def test_cli_analyze_fasta(tmp_path):
    import json
    from pathlib import Path

    fasta_path = Path(__file__).parent.parent / "examples" / "full_antibodies.fasta"
    out_file = tmp_path / "result.json"
    result = runner.invoke(
        app,
        [
            "analyze",
            "--input",
            str(fasta_path),
            "--output",
            str(out_file),
            "--format",
            "json",
        ],
    )
    assert result.exit_code == 0
    assert out_file.exists()

    data = json.loads(out_file.read_text(encoding="utf-8"))
    results_map = {r["sequence_id"]: r for r in data["results"]}
    assert "IGG1_human(P0DOX5)" in results_map

    igg1_res = results_map["IGG1_human(P0DOX5)"]["c_analysis"]
    assert igg1_res["isotype"] == "IgG"
    assert igg1_res["subclass"] == "IgG1"
    assert igg1_res["allotype_summary"] == "G1m17,1"
