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
    """Raised for problems loading structures, resolving chain types, or collisions."""


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


def load_chain_letters(path: str | Path | None = None) -> dict[str, list[str]]:
    """The standardized ``{chain_type: [letters]}`` scheme.

    Keyed by *chain type* — ``class_i_alpha``, ``beta2m``, ``tcr_alpha`` — the
    vocabulary the histo pipelines already speak, rather than prose role names.
    The value is a list because a complex can hold more than one chain of a
    type: the second CD8 alpha in an assembly takes the second letter, so the
    order of the list is load-bearing.
    """
    return json.loads(Path(path or default_chain_letters_path()).read_text())


@dataclass(frozen=True)
class RelabelResult:
    """The outcome of relettering a structure's chains.

    ``structure`` is the parsed Bio.PDB Structure, relettered in place.
    """

    path: Path
    structure: Structure
    chain_types: dict[str, str]  # {source_id: chain_type}, as given
    mapping: dict[str, str]  # {source_id: target_letter}, resolved
    unmapped_chains: list[str]  # chain ids present but not in `chain_types`

    @property
    def records(self) -> list[dict[str, str]]:
        """The mapping as one record per chain, sorted by target letter.

        This is the shape the histo pipelines publish, and the shape the
        previous pipeline's ``chain_relettering.json`` used.
        """
        return sorted(
            (
                {
                    "chain_type": self.chain_types[source_id],
                    "consistent_chain_letter": target,
                    "pdb_chain_letter": source_id,
                }
                for source_id, target in self.mapping.items()
            ),
            key=lambda record: record["consistent_chain_letter"],
        )

    def to_dict(self) -> dict[str, object]:
        """A JSON-serializable view of the result (omits ``structure``)."""
        return {
            "path": str(self.path),
            "chain_types": dict(self.chain_types),
            "mapping": dict(self.mapping),
            "unmapped_chains": list(self.unmapped_chains),
            "records": self.records,
        }


class ChainReletterer:
    """Loads a structure and the standardized chain-letter scheme once;
    ``.reletter()`` relabels any number of chains per an explicit
    ``{source_id: chain_type}`` mapping.

    >>> reletterer = ChainReletterer("8gvi_1_aligned.cif")
    >>> reletterer.letters_for("peptide")
    ['C']
    >>> result = reletterer.reletter({
    ...     "H": "class_i_alpha", "L": "beta2m", "P": "peptide",
    ...     "A": "tcr_alpha", "B": "tcr_beta",
    ... })
    >>> result.mapping
    {'H': 'A', 'L': 'B', 'P': 'C', 'A': 'D', 'B': 'E'}
    """

    def __init__(self, path: str | Path, chain_letters_path: str | Path | None = None):
        self.path = Path(path)
        self.structure = load_structure(self.path)
        self.chain_letters = load_chain_letters(chain_letters_path)

    def letters_for(self, chain_type: str) -> list[str]:
        """The standardized letters a chain type may take, in occurrence order."""
        try:
            return self.chain_letters[chain_type]
        except KeyError:
            raise RelettererError(
                f"Unknown chain type {chain_type!r}; it has no entry in the "
                "chain-letter scheme, so its assembly cannot be relettered"
            ) from None

    def reletter(self, chain_types: dict[str, str]) -> RelabelResult:
        """Reletter chains per ``chain_types`` (source chain id -> chain type).

        A chain type may appear more than once in an assembly — two CD8 alphas,
        four copies of a trimer — so letters are allocated by *occurrence*: the
        nth chain of a type takes the nth letter of its list. Occurrence order
        is the order the chains appear in the file, which is what makes the
        result reproducible.

        Chains not present in ``chain_types`` are left untouched. Raises
        ``RelettererError`` if a source id is not in the structure, a chain type
        has no entry, a type runs out of letters, two chains resolve to the same
        letter, or a target collides with a chain deliberately left unmapped.
        """
        model = self.structure[0]
        available = {chain.id for chain in model}

        unknown = sorted(set(chain_types) - available)
        if unknown:
            raise RelettererError(
                f"No such chain(s) {', '.join(unknown)} in structure; available "
                f"chains: {', '.join(sorted(available))}"
            )

        # Allocate in file order, so the same file always gives the same answer.
        seen: dict[str, int] = {}
        mapping: dict[str, str] = {}
        for chain in model:
            source_id = chain.id
            if source_id not in chain_types:
                continue
            chain_type = chain_types[source_id]
            letters = self.letters_for(chain_type)
            index = seen.get(chain_type, 0)
            if index >= len(letters):
                raise RelettererError(
                    f"Chain type {chain_type!r} appears {index + 1} times but has "
                    f"only {len(letters)} letter(s) ({', '.join(letters)}); the "
                    "scheme needs another letter for it"
                )
            seen[chain_type] = index + 1
            mapping[source_id] = letters[index]

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
                f"Target letter(s) {', '.join(sorted(unmapped_conflicts))} would "
                f"collide with chain(s) left unmapped: {', '.join(sorted(unmapped_conflicts))}"
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
            chain_types=dict(chain_types),
            mapping=mapping,
            unmapped_chains=unmapped,
        )


def reletter_chains(
    path: str | Path,
    chain_types: dict[str, str],
    chain_letters_path: str | Path | None = None,
) -> RelabelResult:
    """Convenience wrapper: reletter a structure file's chains without
    holding onto a ``ChainReletterer`` instance."""
    return ChainReletterer(path, chain_letters_path).reletter(chain_types)
