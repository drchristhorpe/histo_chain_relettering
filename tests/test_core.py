from pathlib import Path

import pytest

from histo_chain_relettering import (
    ChainReletterer,
    RelettererError,
    load_chain_letters,
    load_structure,
    reletter_chains,
    structure_format,
)

FIXTURE = Path(__file__).parent / "fixtures" / "8gvi_1_aligned.cif"

FULL_ROLES = {
    "H": "MHC alpha",
    "L": "Beta-2 microglobulin",
    "P": "Peptide",
    "A": "TCR alpha",
    "B": "TCR beta",
}

# {original chain id: residue count}, confirmed by inspection of the fixture.
ORIGINAL_RESIDUE_COUNTS = {"A": 203, "B": 241, "H": 274, "L": 99, "P": 8}


def test_structure_format_by_extension():
    assert structure_format("x.cif") == "cif"
    assert structure_format("x.mmcif") == "cif"
    assert structure_format("x.pdb") == "pdb"
    assert structure_format("x.ent") == "pdb"
    with pytest.raises(RelettererError, match="Unrecognised structure file extension"):
        structure_format("x.txt")


def test_load_structure_missing_file_raises():
    with pytest.raises(RelettererError, match="No such file"):
        load_structure("does_not_exist.cif")


def test_load_chain_letters_bundled_scheme():
    scheme = load_chain_letters()
    assert scheme == {
        "A": ["MHC alpha"],
        "B": ["MHC beta", "Beta-2 microglobulin"],
        "C": ["Peptide"],
        "D": ["TCR alpha"],
        "E": ["TCR beta"],
    }


def test_resolve_role_exact_and_case_insensitive():
    reletterer = ChainReletterer(FIXTURE)
    assert reletterer.resolve_role("Peptide") == "C"
    assert reletterer.resolve_role("peptide") == "C"
    assert reletterer.resolve_role("  TCR alpha  ") == "D"


def test_resolve_unknown_role_raises_with_valid_list():
    reletterer = ChainReletterer(FIXTURE)
    with pytest.raises(RelettererError, match="Unknown chain role"):
        reletterer.resolve_role("Nonsense role")


def test_full_mapping_resolves_and_relabels():
    reletterer = ChainReletterer(FIXTURE)
    result = reletterer.reletter(FULL_ROLES)

    assert result.mapping == {"H": "A", "L": "B", "P": "C", "A": "D", "B": "E"}
    assert result.unmapped_chains == []

    model = result.structure[0]
    assert {chain.id for chain in model} == {"A", "B", "C", "D", "E"}
    assert [chain.id for chain in model] == ["A", "B", "C", "D", "E"]

    # New letter -> original source id, to check residue counts moved correctly.
    new_to_original = {target: source for source, target in result.mapping.items()}
    for new_id, original_id in new_to_original.items():
        assert len(list(model[new_id])) == ORIGINAL_RESIDUE_COUNTS[original_id]


def test_partial_mapping_leaves_others_untouched():
    reletterer = ChainReletterer(FIXTURE)
    result = reletterer.reletter({"P": "Peptide"})

    assert result.mapping == {"P": "C"}
    assert result.unmapped_chains == ["A", "B", "H", "L"]

    model = result.structure[0]
    assert {chain.id for chain in model} == {"A", "B", "H", "L", "C"}
    assert [chain.id for chain in model] == ["A", "B", "C", "H", "L"]
    assert len(list(model["C"])) == ORIGINAL_RESIDUE_COUNTS["P"]
    assert len(list(model["A"])) == ORIGINAL_RESIDUE_COUNTS["A"]


def test_missing_source_chain_raises():
    reletterer = ChainReletterer(FIXTURE)
    with pytest.raises(RelettererError, match="No such chain 'Z'"):
        reletterer.reletter({"Z": "Peptide"})


def test_two_sources_same_target_letter_raises():
    reletterer = ChainReletterer(FIXTURE)
    with pytest.raises(RelettererError, match="resolve to the same target letter"):
        reletterer.reletter({"H": "MHC beta", "L": "Beta-2 microglobulin"})


def test_target_collides_with_unmapped_chain_raises():
    reletterer = ChainReletterer(FIXTURE)
    # H -> "MHC alpha" -> target "A", but source chain "A" (TCR alpha) is
    # deliberately left unmapped -> would collide.
    with pytest.raises(RelettererError, match="would collide with chain"):
        reletterer.reletter({"H": "MHC alpha"})


def test_reletter_chains_convenience_function():
    result = reletter_chains(FIXTURE, FULL_ROLES)
    assert result.mapping == {"H": "A", "L": "B", "P": "C", "A": "D", "B": "E"}


def test_save_structure_pdb_and_cif_roundtrip(tmp_path):
    from histo_chain_relettering import save_structure

    result = reletter_chains(FIXTURE, FULL_ROLES)

    cif_path = save_structure(result.structure, tmp_path / "out.cif", fmt="cif")
    reloaded_cif = load_structure(cif_path)
    assert {chain.id for chain in reloaded_cif[0]} == {"A", "B", "C", "D", "E"}
    assert [chain.id for chain in reloaded_cif[0]] == ["A", "B", "C", "D", "E"]

    pdb_path = save_structure(result.structure, tmp_path / "out.pdb", fmt="pdb")
    reloaded_pdb = load_structure(pdb_path)
    assert {chain.id for chain in reloaded_pdb[0]} == {"A", "B", "C", "D", "E"}
