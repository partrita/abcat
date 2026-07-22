from pathlib import Path
from typing import Optional
import typer
from rich.console import Console
from rich.table import Table

from ab_numbering_lite import analyze_chain, analyze_vdomain, analyze_cdomain
from ab_numbering_lite.schemas import BatchAnalysisResult
from ab_numbering_lite.utils import read_fasta_file, export_json, export_csv

app = typer.Typer(
    name="abnl",
    help="ab-numbering-lite: ANARCII V-domain CDR analyzer and Constant domain Allotype classifier",
    add_completion=False,
)
console = Console()


@app.command("analyze")
def analyze_command(
    sequence: Optional[str] = typer.Option(
        None, "--sequence", "-s", help="Single amino acid sequence"
    ),
    input_file: Optional[Path] = typer.Option(
        None, "--input", "-i", help="Path to input FASTA file"
    ),
    output: Optional[Path] = typer.Option(
        None, "--output", "-o", help="Output file path"
    ),
    output_format: str = typer.Option(
        "json", "--format", "-f", help="Output format: json or csv"
    ),
    scheme: str = typer.Option(
        "imgt",
        "--scheme",
        help="Numbering scheme: imgt, kabat, martin, chothia, aho (default: imgt)",
    ),
):
    """Analyzes antibody chain sequence(s) for V-domain CDRs and Constant region Allotypes."""
    records = []

    if sequence:
        records.append(("query", sequence))
    elif input_file and input_file.exists():
        records = read_fasta_file(input_file)
    else:
        console.print(
            "[red]Error:[/red] Please provide either --sequence or a valid --input FASTA file."
        )
        raise typer.Exit(code=1)

    results = []
    for seq_id, seq_str in records:
        res = analyze_chain(seq_str, sequence_id=seq_id, scheme=scheme)
        results.append(res)

    batch_res = BatchAnalysisResult(
        total_sequences=len(results),
        successful_analyses=len(results),
        results=results,
    )

    formatted_output = (
        export_json(batch_res)
        if output_format.lower() == "json"
        else export_csv(batch_res)
    )

    if output:
        output.write_text(formatted_output, encoding="utf-8")
        console.print(f"[green]Saved analysis results to {output}[/green]")
    else:
        _render_rich_table(results)


@app.command("batch")
def batch_command(
    input_file: Path = typer.Option(
        ..., "--input", "-i", help="Path to input FASTA file"
    ),
    output: Path = typer.Option(..., "--output", "-o", help="Output file path"),
    output_format: str = typer.Option(
        "csv", "--format", "-f", help="Output format: csv or json"
    ),
    scheme: str = typer.Option(
        "imgt",
        "--scheme",
        help="Numbering scheme: imgt, kabat, martin, chothia, aho (default: imgt)",
    ),
):
    """Batch processes a FASTA file of antibody sequences."""
    if not input_file.exists():
        console.print(f"[red]Error:[/red] Input file not found: {input_file}")
        raise typer.Exit(code=1)

    records = read_fasta_file(input_file)
    results = [
        analyze_chain(seq_str, sequence_id=seq_id, scheme=scheme)
        for seq_id, seq_str in records
    ]

    batch_res = BatchAnalysisResult(
        total_sequences=len(results),
        successful_analyses=len(results),
        results=results,
    )

    formatted_output = (
        export_json(batch_res)
        if output_format.lower() == "json"
        else export_csv(batch_res)
    )
    output.write_text(formatted_output, encoding="utf-8")
    console.print(
        f"[bold green]Batch processing completed.[/bold green] Wrote {len(results)} records to [cyan]{output}[/cyan]"
    )


@app.command("vdomain")
def vdomain_command(
    sequence: Optional[str] = typer.Option(None, "--sequence", "-s"),
    input_file: Optional[Path] = typer.Option(None, "--input", "-i"),
    scheme: str = typer.Option(
        "imgt",
        "--scheme",
        help="Numbering scheme: imgt, kabat, martin, chothia, aho (default: imgt)",
    ),
):
    """Runs ANARCII V-domain numbering and CDR extraction only."""
    seq_str = sequence
    if not seq_str and input_file and input_file.exists():
        records = read_fasta_file(input_file)
        seq_str = records[0][1]

    if not seq_str:
        console.print("[red]Error:[/red] Provide sequence or FASTA file.")
        raise typer.Exit(code=1)

    res = analyze_vdomain(seq_str, scheme=scheme)
    console.print(f"[bold cyan]Chain Type:[/bold cyan] {res.chain_type.value}")
    console.print(f"[bold cyan]Scheme:[/bold cyan] {res.scheme.upper()}")
    console.print(
        f"[bold cyan]V-Domain Length:[/bold cyan] {len(res.v_domain_sequence)}"
    )

    for cdr_name, cdr in res.cdrs.items():
        console.print(
            f"  [yellow]{cdr_name}:[/yellow] {cdr.sequence} (Len: {cdr.length})"
        )


@app.command("cdomain")
def cdomain_command(
    sequence: str = typer.Option(
        ..., "--sequence", "-s", help="Constant region or full chain sequence"
    ),
):
    """Runs Constant Region Isotype, Subclass, and Allotype analysis only."""
    res = analyze_cdomain(sequence)
    console.print(f"[bold green]Isotype:[/bold green] {res.isotype}")
    console.print(f"[bold green]Subclass:[/bold green] {res.subclass}")
    console.print(f"[bold green]Matched Gene:[/bold green] {res.matched_c_gene}")

    if res.allotypes:
        console.print("[bold yellow]Allotypes:[/bold yellow]")
        for a in res.allotypes:
            status_color = "green" if a.status == "present" else "dim"
            console.print(
                f"  [{status_color}]• {a.allotype}[/{status_color}] ({a.domain}): {a.status}"
            )


def _render_rich_table(results):
    table = Table(title="Antibody Sequence Analysis Summary")
    table.add_column("Seq ID", style="cyan", no_wrap=True)
    table.add_column("Chain", style="magenta")
    table.add_column("Scheme", style="yellow")
    table.add_column("Isotype", style="green")
    table.add_column("Subclass", style="bold green")
    table.add_column("Allotypes", style="yellow")
    table.add_column("CDR1", style="blue")
    table.add_column("CDR2", style="blue")
    table.add_column("CDR3", style="bold blue")

    for r in results:
        chain = r.v_analysis.chain_type.value if r.v_analysis else "-"
        scheme_str = r.v_analysis.scheme.upper() if r.v_analysis else "IMGT"
        iso = r.c_analysis.isotype if r.c_analysis else "-"
        sub = r.c_analysis.subclass if r.c_analysis else "-"
        allos = (
            ",".join(
                [
                    a.allotype
                    for a in (r.c_analysis.allotypes if r.c_analysis else [])
                    if a.status == "present"
                ]
            )
            or "None"
        )
        cdr1 = (
            r.v_analysis.cdrs.get("CDR1").sequence
            if (r.v_analysis and "CDR1" in r.v_analysis.cdrs)
            else "-"
        )
        cdr2 = (
            r.v_analysis.cdrs.get("CDR2").sequence
            if (r.v_analysis and "CDR2" in r.v_analysis.cdrs)
            else "-"
        )
        cdr3 = (
            r.v_analysis.cdrs.get("CDR3").sequence
            if (r.v_analysis and "CDR3" in r.v_analysis.cdrs)
            else "-"
        )

        table.add_row(
            r.sequence_id, chain, scheme_str, iso, sub, allos, cdr1, cdr2, cdr3
        )

    console.print(table)


if __name__ == "__main__":
    app()
