from pathlib import Path

from click.testing import CliRunner

from histo_chain_relettering.cli import main

FIXTURE = Path(__file__).parent / "fixtures" / "8gvi_1_aligned.cif"

FULL_MAP_ARGS = [
    "--map", "H=MHC alpha",
    "--map", "L=Beta-2 microglobulin",
    "--map", "P=Peptide",
    "--map", "A=TCR alpha",
    "--map", "B=TCR beta",
]


def test_cli_writes_default_output_same_format_as_input(tmp_path):
    runner = CliRunner()
    result = runner.invoke(main, [str(FIXTURE), *FULL_MAP_ARGS])
    assert result.exit_code == 0, result.output

    expected = FIXTURE.with_name(f"{FIXTURE.stem}_relettered.cif")
    try:
        assert expected.exists()
        assert "Wrote relettered structure to" in result.output
    finally:
        expected.unlink(missing_ok=True)


def test_cli_format_and_output_override(tmp_path):
    runner = CliRunner()
    out_path = tmp_path / "out.pdb"
    result = runner.invoke(
        main, [str(FIXTURE), *FULL_MAP_ARGS, "--format", "pdb", "--output", str(out_path)]
    )
    assert result.exit_code == 0, result.output
    assert out_path.exists()


def test_cli_partial_mapping_reports_unmapped(tmp_path):
    runner = CliRunner()
    out_path = tmp_path / "out.cif"
    result = runner.invoke(main, [str(FIXTURE), "--map", "P=Peptide", "--output", str(out_path)])
    assert result.exit_code == 0, result.output
    assert "Left unmapped" in result.output
    assert out_path.exists()


def test_cli_malformed_map_entry_errors():
    runner = CliRunner()
    result = runner.invoke(main, [str(FIXTURE), "--map", "P-Peptide"])
    assert result.exit_code != 0


def test_cli_unknown_role_errors():
    runner = CliRunner()
    result = runner.invoke(main, [str(FIXTURE), "--map", "P=Nonsense"])
    assert result.exit_code != 0
    assert "Unknown chain role" in result.output


def test_cli_missing_file_errors():
    runner = CliRunner()
    result = runner.invoke(main, ["does_not_exist.cif", "--map", "P=Peptide"])
    assert result.exit_code != 0


def test_cli_requires_at_least_one_map():
    runner = CliRunner()
    result = runner.invoke(main, [str(FIXTURE)])
    assert result.exit_code != 0
