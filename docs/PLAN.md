# histo_chain_relettering — Design & Implementation Plan

## 1. Purpose

`histo_chain_relettering` reletters a structure file's (PDB/mmCIF) chain
IDs to the family's standardized single-letter scheme, defined in
`chain_letters.json`:

```json
{
    "A": ["MHC alpha"],
    "B": ["MHC beta", "Beta-2 microglobulin"],
    "C": ["Peptide"],
    "D": ["TCR alpha"],
    "E": ["TCR beta"]
}
```

The caller states, per source chain, which biological role it plays
(e.g. `H` is `"MHC alpha"`), and the tool resolves that role to its
standardized letter and relabels the chain. Chains not mentioned are
left untouched. Ships as:

1. A Python library (`import histo_chain_relettering`)
2. A CLI tool (`histo-chain-relettering`), built with Click
3. A Claude Code / Claude Desktop skill wrapping the CLI

## 2. Tooling

- Python **3.14**, managed with **uv**
- **Biopython** (`Bio.PDB`) for structure parsing/writing, following the
  `load_structure`/`save_structure` pattern established in
  `histo_altloc_cleaner`/`histo_aligner`
- **Click** for the CLI, **Rich** for console table output
- `pytest` for tests, run against a real committed fixture

## 3. Mapping input: role name, not raw letter

Confirmed with the user: `--map` takes a source chain id and a
**biological role name** (`--map H="MHC alpha"`), not a raw target
letter directly. The role string must match — case-insensitively,
whitespace-trimmed — one of the strings listed under some letter in
`chain_letters.json`; the tool resolves it to that letter. This is
self-documenting (the CLI invocation reads as "chain H is the MHC
alpha") and validates every mapping against the family's controlled
vocabulary, rather than requiring the caller to already know/compute
the target letter by hand. An unknown role raises `RelettererError`
listing every valid role name.

Chains present in the structure but absent from `--map` are left
completely unchanged — the tool doesn't require full chain coverage
(e.g. water/heteroatom chains don't need mapping).

## 4. Reference data (`histo_chain_relettering/data/chain_letters.json`)

Bundled inside the package (`src/histo_chain_relettering/data/`, not the
repo root) per `CONVENTIONS.md` §3, so it ships with the wheel. Loaded
via plain `pathlib`, mirroring `histo_aligner.reference_data`'s pattern
(no `importlib.resources`):

```python
def default_chain_letters_path() -> Path:
    return Path(__file__).parent / "data" / "chain_letters.json"
```

Not overridable via a CLI flag — it's the family-wide standard, not a
per-run configuration knob.

## 5. Two-step relettering to avoid id collisions

A single-pass rename (`chain.id = target_letter` for every mapped
chain, in whatever order) is unsafe whenever a target letter is already
in use by another chain in the structure that hasn't been renamed yet.
This is the common case, not an edge case: 8gvi's own chains include
`A` (TCR alpha, mapped to target `D`) *and* `H` (MHC alpha, mapped to
target `A`) — renaming `H` to `A` directly would collide with the
still-present original chain `A`.

Checked against Biopython's `Bio.PDB.Entity.id` setter: assigning an id
that collides with an existing sibling does **not** raise — it emits a
`BiopythonWarning` and silently overwrites the parent's `child_dict`
entry, leaving the structure in an inconsistent state (the old chain
object becomes unreachable by `model[id]` even though it's still in
`child_list`). This must be avoided entirely, not caught after the fact.

The confirmed fix — a two-pass rename, exactly as specified by the
user:

1. **Proxy pass**: every mapped chain's id is renamed to its own id
   repeated 4 times (`P` -> `PPPP`, `H` -> `HHHH`, ...). A 4-character
   proxy can never collide with any single-character source or target
   id in this scheme, so every proxy assignment is guaranteed
   collision-free regardless of order.
2. **Final pass**: every proxy is renamed to its resolved target letter
   (`PPPP` -> `C`, `HHHH` -> `A`, ...). By this point every source chain
   that would have collided with a target letter has already vacated
   its original id in pass 1, so this pass is also collision-free.

Both passes run entirely in memory before anything is written — the
intermediate 4-letter ids are never serialized to disk (which matters
for classic PDB format's 1-character chain id column limit).

## 6. Validation

Beyond role resolution (§3), `reletter()` checks two collision cases
before mutating anything:

- **Two source chains resolve to the same target letter** (e.g. two
  chains both mapped to `"MHC beta"`/`"Beta-2 microglobulin"`, which
  share letter `B`) — raises `RelettererError` naming both source
  chains and the shared target.
- **A target letter collides with an chain deliberately left
  unmapped** (e.g. mapping only `H` -> `"MHC alpha"` (target `A`) while
  chain `A` itself is left out of `--map`) — raises `RelettererError`,
  since the two-step process only protects mapped chains from each
  other, not from chains the caller chose not to touch.

## 7. Output format

`--format {pdb,cif}`, **default: same as the input file's format** (not
a fixed default) — the user's stated requirement is "return the same
type of file it was given, with an optional flag" for the other format.
Independent of the two-step in-memory rename, which never touches
serialization.

## 8. Library API (`histo_chain_relettering/core.py`)

```python
class RelettererError(ValueError):
    """Raised for problems loading structures, resolving roles, or reletter collisions."""

def load_structure(path: str | Path, structure_id: str | None = None) -> Structure: ...
def save_structure(structure: Structure, path: str | Path, fmt: str = "pdb") -> Path: ...
def structure_format(path: str | Path) -> str:
    """'cif' or 'pdb', from the file extension."""

def default_chain_letters_path() -> Path: ...
def load_chain_letters() -> dict[str, list[str]]: ...

@dataclass(frozen=True)
class RelabelResult:
    path: Path
    structure: Structure
    roles: dict[str, str]            # {source_id: role_name}, as given
    mapping: dict[str, str]          # {source_id: target_letter}, resolved
    unmapped_chains: list[str]       # chain ids present but not in `roles`

    def to_dict(self) -> dict: ...   # JSON-serializable, omits `structure`

class ChainReletterer:
    """Parses a structure and the bundled chain_letters.json once;
    `.reletter()` performs the two-step rename.

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
    def __init__(self, path: str | Path) -> None: ...
    def resolve_role(self, role: str) -> str: ...
    def reletter(self, roles: dict[str, str]) -> RelabelResult: ...

def reletter_chains(path: str | Path, roles: dict[str, str]) -> RelabelResult:
    """One-shot convenience wrapper."""
```

## 9. CLI (`histo_chain_relettering/cli.py`)

```
histo-chain-relettering FILENAME --map CHAIN=ROLE [--map CHAIN=ROLE ...]
                         [-f/--format {pdb,cif}] [-o/--output PATH]
```

| Flag | Short | Default | Notes |
|---|---|---|---|
| `FILENAME` | (positional) | required | `.cif`/`.mmcif` or `.pdb`/`.ent` |
| `--map` | `-m` | required, repeatable | `CHAIN=ROLE`, e.g. `-m H="MHC alpha"` |
| `--format` | `-f` | same as input | `pdb` or `cif` |
| `--output` | `-o` | `<stem>_relettered.<format>` | structure output path |

CLI prints (Rich) a table of source chain id / role / resolved target
letter, a line noting any chains left unmapped, then confirms the
output path.

## 10. Claude skill

`skills/histo-chain-relettering/SKILL.md` — wraps the CLI: relettering
a structure's chains to the standardized scheme given a role for each
chain to touch.

## 11. Package layout

```
histo_chain_relettering/
  pyproject.toml
  README.md
  CLAUDE.md
  CHANGELOG.md
  docs/PLAN.md
  .gitignore
  src/histo_chain_relettering/
    __init__.py
    core.py               # load/save_structure, ChainReletterer, reletter_chains()
    cli.py                # Click CLI, entry point: histo-chain-relettering
    py.typed
    data/
      chain_letters.json  # bundled standardized role -> letter scheme
  tests/
    fixtures/
      8gvi_1_aligned.cif  # copy of histo_com's fixture (chains A/B/H/L/P)
    test_core.py
    test_cli.py
  skills/histo-chain-relettering/SKILL.md
  tmp/
    .gitkeep
```

## 12. Test fixture

`8gvi_1_aligned.cif`, copied from `histo_com/tests/fixtures/` (same
source, per-repo copy convention — siblings aren't assumed checked out
next to each other). Confirmed by inspection: chains `A` (203 res, TCR
alpha), `B` (241 res, TCR beta), `H` (274 res, MHC alpha), `L` (99 res,
beta-2 microglobulin), `P` (8 res, peptide) — a genuine pMHC-TCR
complex that exercises the full standardized scheme (`A`-`E`) in one
fixture, and specifically exercises the collision case in §5 (`A` is
both a source id and a target letter).

## 13. Testing plan

Real fixture only, no synthetic structures:

- `structure_format()`/`load_structure()`: cif/pdb dispatch by
  extension, missing file, unknown extension.
- `load_chain_letters()` returns the bundled 5-entry scheme.
- `resolve_role()`: exact match, case-insensitive match, unknown role
  raises `RelettererError` listing valid roles.
- Full 8gvi mapping (`H`/`L`/`P`/`A`/`B` -> all 5 roles) resolves to
  `{"H": "A", "L": "B", "P": "C", "A": "D", "B": "E"}`; after
  `.reletter()`, the structure has exactly chains `{A, B, C, D, E}` with
  residue counts matching the *original* chain each now represents
  (e.g. new chain `A` has 274 residues — MHC's original count from old
  chain `H`).
  `unmapped_chains` is empty.
- Partial mapping leaves the unmapped chain's id and residues
  untouched.
- Collision cases raise `RelettererError`: two source chains mapped to
  roles sharing one target letter; a target letter colliding with a
  chain deliberately left unmapped.
- `--format` defaults to the input's own format; explicit `cif`/`pdb`
  both round-trip and parse back with the expected relettered chain ids.
- CLI (`CliRunner`): expected output file written, malformed `--map`
  entry (no `=`) errors, unknown role errors with non-zero exit, missing
  file errors.

## 14. Workflow

1. Write this plan.
2. Move `data/chain_letters.json` into
   `src/histo_chain_relettering/data/chain_letters.json`; scaffold
   `pyproject.toml`, `src/` skeleton.
3. Implement `core.py`, `cli.py`.
4. Copy the 8gvi fixture, write and run tests.
5. Write `README.md`, `CLAUDE.md`, `CHANGELOG.md`,
   `skills/histo-chain-relettering/SKILL.md`, `.gitignore`.
6. Manually exercise `uv run histo-chain-relettering ...` against the
   8gvi fixture with the full role mapping from §12, inspect the output
   structure in `tmp/`.
7. Present for approval, then commit.
