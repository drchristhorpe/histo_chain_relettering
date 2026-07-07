"""Structure loading and two-step chain relettering."""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

from Bio.PDB.mmcifio import MMCIFIO
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.PDBIO import PDBIO
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Structure import Structure

_CIF_SUFFIXES = {".cif", ".mmcif"}
_PDB_SUFFIXES = {".pdb", ".ent"}

_WRITERS = {"pdb": PDBIO, "cif": MMCIFIO}


class RelettererError(ValueError):
    """Raised for problems loading structures, resolving roles, or reletter collisions."""


def structure_format(path: str | Path) -> str:
    """'cif' or 'pdb', chosen from ``path``'s extension (case-insensitive)."""
    suffix = Path(path).suffix.lower()
    if suffix in _CIF_SUFFIXES:
        return "cif"
    if suffix in _PDB_SUFFIXES:
        return "pdb"
    raise RelettererError(
        f"Unrecognised structure file extension {suffix!r} for {path}; "
        "expected one of .cif, .mmcif, .pdb, .ent"
    )


def load_structure(path: str | Path, structure_id: str | None = None) -> Structure:
    """Parse a PDB or mmCIF file into a Bio.PDB Structure.

    Format is chosen from the file extension (case-insensitive):
    ``.cif``/``.mmcif`` -> mmCIF, ``.pdb``/``.ent`` -> legacy PDB.
    """
    path = Path(path)
    if not path.is_file():
        raise RelettererError(f"No such file: {path}")

    fmt = structure_format(path)
    sid = structure_id or path.stem
    parser = MMCIFParser(QUIET=True) if fmt == "cif" else PDBParser(QUIET=True)

    structure = parser.get_structure(sid, str(path))
    if len(structure) == 0:
        raise RelettererError(f"No models found in {path}")
    return structure


def save_structure(structure: Structure, path: str | Path, fmt: str = "pdb") -> Path:
    """Write a structure to a PDB or mmCIF file."""
    path = Path(path)
    if fmt not in _WRITERS:
        raise RelettererError(f"fmt must be one of {sorted(_WRITERS)}, got {fmt!r}")
    io = _WRITERS[fmt]()
    io.set_structure(structure)
    io.save(str(path))
    return path


def default_chain_letters_path() -> Path:
    return Path(__file__).parent / "data" / "chain_letters.json"


def load_chain_letters() -> dict[str, list[str]]:
    """The bundled standardized ``{target_letter: [role names]}`` scheme."""
    return json.loads(default_chain_letters_path().read_text())


def _role_to_letter_map(chain_letters: dict[str, list[str]]) -> dict[str, str]:
    """Reverse lookup: role name (lowercased, trimmed) -> target letter."""
    return {
        role.strip().lower(): letter for letter, roles in chain_letters.items() for role in roles
    }


@dataclass(frozen=True)
class RelabelResult:
    """The outcome of relettering a structure's chains.

    ``structure`` is the parsed Bio.PDB Structure, relettered in place.
    """

    path: Path
    structure: Structure
    roles: dict[str, str]  # {source_id: role_name}, as given
    mapping: dict[str, str]  # {source_id: target_letter}, resolved
    unmapped_chains: list[str]  # chain ids present but not in `roles`

    def to_dict(self) -> dict[str, object]:
        """A JSON-serializable view of the result (omits ``structure``)."""
        return {
            "path": str(self.path),
            "roles": dict(self.roles),
            "mapping": dict(self.mapping),
            "unmapped_chains": list(self.unmapped_chains),
        }


class ChainReletterer:
    """Loads a structure and the bundled standardized chain-letter scheme
    once; ``.reletter()`` relabels any number of chains per an explicit
    ``{source_id: role_name}`` mapping.

    >>> reletterer = ChainReletterer("8gvi_1_aligned.cif")
    >>> reletterer.resolve_role("Peptide")
    'C'
    >>> result = reletterer.reletter({
    ...     "H": "MHC alpha", "L": "Beta-2 microglobulin", "P": "Peptide",
    ...     "A": "TCR alpha", "B": "TCR beta",
    ... })
    >>> result.mapping
    {'H': 'A', 'L': 'B', 'P': 'C', 'A': 'D', 'B': 'E'}
    """

    def __init__(self, path: str | Path):
        self.path = Path(path)
        self.structure = load_structure(self.path)
        self.chain_letters = load_chain_letters()
        self._role_to_letter = _role_to_letter_map(self.chain_letters)

    def resolve_role(self, role: str) -> str:
        """The standardized target letter for a biological role name
        (case-insensitive, whitespace-trimmed match against
        ``chain_letters.json``)."""
        try:
            return self._role_to_letter[role.strip().lower()]
        except KeyError:
            valid = sorted(self._role_to_letter, key=str.lower)
            raise RelettererError(
                f"Unknown chain role {role!r}; valid roles are: {', '.join(valid)}"
            ) from None

    def reletter(self, roles: dict[str, str]) -> RelabelResult:
        """Reletter chains per ``roles`` (source chain id -> role name).

        Chains not present in ``roles`` are left untouched. Raises
        ``RelettererError`` if a source id doesn't exist in the
        structure, a role name doesn't resolve, two source chains
        resolve to the same target letter, or a target letter collides
        with a chain deliberately left unmapped.

        On return, every chain in the model (mapped or not) is reordered
        by its (possibly new) id, so writing the structure lists chains
        in id order rather than original/rename order.
        """
        model = self.structure[0]
        available = {chain.id for chain in model}

        mapping: dict[str, str] = {}
        for source_id, role in roles.items():
            if source_id not in available:
                raise RelettererError(
                    f"No such chain {source_id!r} in structure; available chains: "
                    f"{', '.join(sorted(available))}"
                )
            mapping[source_id] = self.resolve_role(role)

        sources_by_target: dict[str, list[str]] = {}
        for source_id, target in mapping.items():
            sources_by_target.setdefault(target, []).append(source_id)
        clashes = {t: ids for t, ids in sources_by_target.items() if len(ids) > 1}
        if clashes:
            detail = "; ".join(f"{t} <- {', '.join(sorted(ids))}" for t, ids in clashes.items())
            raise RelettererError(f"Multiple source chains resolve to the same target letter: {detail}")

        unmapped = sorted(available - set(mapping))
        unmapped_conflicts = set(mapping.values()) & set(unmapped)
        if unmapped_conflicts:
            raise RelettererError(
                "Target letter(s) "
                f"{', '.join(sorted(unmapped_conflicts))} would collide with chain(s) "
                f"left unmapped: {', '.join(sorted(unmapped_conflicts))}"
            )

        # Step 1: proxy rename (source id repeated 4x) so every mapped
        # chain vacates its original id before any final target id is
        # assigned. See docs/PLAN.md §5 for why this is required.
        for source_id in mapping:
            model[source_id].id = source_id * 4

        # Step 2: proxy -> resolved target letter.
        for source_id, target in mapping.items():
            model[source_id * 4].id = target

        # Reorder for writing: chain id order, not original/rename order.
        # Only `child_list` (plain list, drives PDBIO/MMCIFIO iteration)
        # needs sorting; `child_dict` lookups are unaffected by list order.
        model.child_list.sort(key=lambda chain: chain.id)

        return RelabelResult(
            path=self.path,
            structure=self.structure,
            roles=dict(roles),
            mapping=mapping,
            unmapped_chains=unmapped,
        )


def reletter_chains(path: str | Path, roles: dict[str, str]) -> RelabelResult:
    """Convenience wrapper: reletter a structure file's chains without
    holding onto a ``ChainReletterer`` instance."""
    return ChainReletterer(path).reletter(roles)
