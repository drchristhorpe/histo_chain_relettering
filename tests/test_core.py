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
    "H": "class_i_alpha",
    "L": "beta2m",
    "P": "peptide",
    "A": "tcr_alpha",
    "B": "tcr_beta",
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


def test_the_bundled_scheme_is_keyed_by_chain_type():
    """Keyed by chain type, not by letter, and the value is a list because a
    complex can hold more than one chain of a type."""
    scheme = load_chain_letters()

    assert scheme["class_i_alpha"] == ["A"]
    assert scheme["beta2m"] == ["B"]
    assert scheme["peptide"] == ["C"]
    assert scheme["tcr_alpha"] == ["D"]
    assert scheme["tcr_beta"] == ["E"]
    # Types that can appear twice carry a letter for each.
    assert len(scheme["cd8a"]) == 2


def test_letters_for_returns_the_whole_list():
    reletterer = ChainReletterer(FIXTURE)

    assert reletterer.letters_for("peptide") == ["C"]
    assert reletterer.letters_for("cd8a") == ["F", "G"]


def test_an_unknown_chain_type_raises():
    """Chain types are slugs from a fixed vocabulary, so this is a caller bug
    rather than something to normalise away."""
    reletterer = ChainReletterer(FIXTURE)

    with pytest.raises(RelettererError, match="Unknown chain type"):
        reletterer.letters_for("Nonsense role")


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
    result = reletterer.reletter({"P": "peptide"})

    assert result.mapping == {"P": "C"}
    assert result.unmapped_chains == ["A", "B", "H", "L"]

    model = result.structure[0]
    assert {chain.id for chain in model} == {"A", "B", "H", "L", "C"}
    assert [chain.id for chain in model] == ["A", "B", "C", "H", "L"]
    assert len(list(model["C"])) == ORIGINAL_RESIDUE_COUNTS["P"]
    assert len(list(model["A"])) == ORIGINAL_RESIDUE_COUNTS["A"]


def test_missing_source_chain_raises():
    reletterer = ChainReletterer(FIXTURE)
    with pytest.raises(RelettererError, match="No such chain"):
        reletterer.reletter({"Z": "peptide"})


def test_two_chain_types_sharing_a_letter_raises():
    """Two chains of the *same* type is normal — they take successive letters.
    Two different types that both want `D` is a real conflict."""
    reletterer = ChainReletterer(FIXTURE)

    with pytest.raises(RelettererError, match="resolve to the same target letter"):
        reletterer.reletter({"H": "tcr_alpha", "L": "ab_heavy"})


def test_repeated_chain_types_take_successive_letters():
    """The reason the scheme stores a list at all."""
    reletterer = ChainReletterer(FIXTURE)

    result = reletterer.reletter({"H": "cd8a", "L": "cd8a", "P": "peptide"})

    assert sorted(result.mapping.values()) == ["C", "F", "G"]
    assert result.mapping["P"] == "C"


def test_running_out_of_letters_raises():
    reletterer = ChainReletterer(FIXTURE)

    with pytest.raises(RelettererError, match="only 1 letter"):
        reletterer.reletter({"H": "peptide", "L": "peptide"})


def test_records_are_the_published_shape():
    reletterer = ChainReletterer(FIXTURE)

    records = reletterer.reletter(FULL_ROLES).records

    assert records[0] == {
        "chain_type": "class_i_alpha",
        "consistent_chain_letter": "A",
        "pdb_chain_letter": "H",
    }
    assert [r["consistent_chain_letter"] for r in records] == ["A", "B", "C", "D", "E"]


def test_allocation_does_not_depend_on_mapping_order():
    """Letters are allocated in *file* order, so the same file always gives the
    same answer however the caller happened to build its dict."""
    forward = ChainReletterer(FIXTURE).reletter({"H": "cd8a", "L": "cd8a"})
    reverse = ChainReletterer(FIXTURE).reletter({"L": "cd8a", "H": "cd8a"})

    assert forward.mapping == reverse.mapping


def test_target_collides_with_unmapped_chain_raises():
    reletterer = ChainReletterer(FIXTURE)
    # H -> "class_i_alpha" -> target "A", but source chain "A" (TCR alpha) is
    # deliberately left unmapped -> would collide.
    with pytest.raises(RelettererError, match="would collide with chain"):
        reletterer.reletter({"H": "class_i_alpha"})


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
