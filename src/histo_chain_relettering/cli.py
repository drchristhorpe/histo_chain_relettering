"""Command line interface for histo_chain_relettering."""

from __future__ import annotations

from pathlib import Path

import click
from rich.console import Console
from rich.table import Table

from histo_chain_relettering.core import ChainReletterer, RelettererError, save_structure, structure_format

console = Console()


def _parse_mapping(entries: tuple[str, ...]) -> dict[str, str]:
    chain_types: dict[str, str] = {}
    for entry in entries:
        if "=" not in entry:
            raise click.BadParameter(f"expected CHAIN=CHAIN_TYPE, got {entry!r}", param_hint="'--map'")
        chain_id, chain_type = entry.split("=", 1)
        chain_id = chain_id.strip()
        chain_type = chain_type.strip()
        if not chain_id or not chain_type:
            raise click.BadParameter(f"expected CHAIN=CHAIN_TYPE, got {entry!r}", param_hint="'--map'")
        if chain_id in chain_types:
            raise click.BadParameter(f"chain {chain_id!r} mapped more than once", param_hint="'--map'")
        chain_types[chain_id] = chain_type
    return chain_types


@click.command()
@click.argument("filename", type=click.Path(exists=True, dir_okay=False, path_type=Path))
@click.option(
    "--map",
    "-m",
    "mapping",
    multiple=True,
    required=True,
    metavar="CHAIN=CHAIN_TYPE",
    help="Source chain id to chain type, e.g. -m P=peptide. "
    "Repeatable. CHAIN_TYPE must be a key of chain_letters.json. "
    "Chains not mapped are left unchanged.",
)
@click.option(
    "--format",
    "-f",
    "output_format",
    type=click.Choice(["pdb", "cif"]),
    default=None,
    help="Output structure format (default: same as the input file).",
)
@click.option(
    "--output",
    "-o",
    type=click.Path(dir_okay=False, path_type=Path),
    default=None,
    help="Structure output path (defaults to '<stem>_relettered.<format>').",
)
def main(filename: Path, mapping: tuple[str, ...], output_format: str | None, output: Path | None) -> None:
    """Reletter FILENAME's chains to the standardized single-letter scheme.

    FILENAME is a .cif/.mmcif or .pdb/.ent structure file.
    """
    chain_types = _parse_mapping(mapping)

    try:
        reletterer = ChainReletterer(filename)
        result = reletterer.reletter(chain_types)
    except RelettererError as exc:
        raise click.ClickException(str(exc)) from exc

    fmt = output_format or structure_format(filename)
    output_path = output or filename.with_name(f"{filename.stem}_relettered.{fmt}")
    save_structure(result.structure, output_path, fmt=fmt)

    table = Table(title=f"Chain relettering: {filename}")
    table.add_column("source chain")
    table.add_column("chain type")
    table.add_column("target letter")
    for source_id, chain_type in result.chain_types.items():
        table.add_row(source_id, chain_type, result.mapping[source_id])
    console.print(table)

    if result.unmapped_chains:
        console.print(f"[yellow]Left unmapped:[/yellow] {', '.join(result.unmapped_chains)}")

    console.print(f"[bold green]Wrote relettered structure to[/bold green] {output_path}")


if __name__ == "__main__":
    main()
