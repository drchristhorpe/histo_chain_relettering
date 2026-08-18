# histo-chain-relettering

Reletter a 3D biological structure file's (PDB/mmCIF) chains to the
family's standardized single-letter scheme, defined in
[`chain_letters.json`](src/histo_chain_relettering/data/chain_letters.json):

| Chain type | Letter(s) |
|---|---|
| `class_i_alpha`, `truncated_class_i_alpha`, `hybrid_class_i_alpha`, `mr1`, `cd1d`, … | `A` |
| `beta2m` | `B` |
| `peptide` | `C` |
| `tcr_alpha`, `ab_heavy`, `nanobody` | `D` |
| `tcr_beta`, `ab_light` | `E` |
| `cd8a`, `cd8b` | `F`, `G` |
| `nkg2a`, `nkg2d`, `ly49a`, `cd94`, `kir` | `K`, `L` |

56 chain types in all. You state which **chain type** each source chain
is; the tool resolves that to a letter and relabels the chain in two
steps to avoid id collisions (see "Notes and limitations").

The value is a *list* because an assembly can hold more than one chain of
a type — two CD8 alphas, four copies of a trimer. The **nth chain of a
type takes the nth letter**, allocated in the order the chains appear in
the file so the same file always gives the same answer.

Chain types are the vocabulary the histo pipelines already use, so a
caller that has a `chain_types.json` can pass it through unchanged.

Built on [Biopython](https://biopython.org/), it ships as:

- a Python library — `import histo_chain_relettering`
- a CLI tool — `histo-chain-relettering`
- a [Claude Code / Claude Desktop skill](skills/histo-chain-relettering/SKILL.md)

Requires Python 3.14+.

## Install

```bash
uv sync                 # dev environment, from a checkout
uv tool install .       # install the `histo-chain-relettering` CLI globally
# or
pip install .
```

## CLI usage

```
histo-chain-relettering FILENAME --map CHAIN=CHAIN_TYPE [--map CHAIN=CHAIN_TYPE ...] [--format pdb|cif] [--output PATH]
```

- `FILENAME` — a `.cif`/`.mmcif` or `.pdb`/`.ent` structure file.
- `--map`, `-m` (required, repeatable) — `CHAIN=CHAIN_TYPE`, e.g. `-m P=peptide`.
  `ROLE` must match one of the role strings in `chain_letters.json`
  (case-insensitive). Chains not mentioned are left unchanged.
- `--format`, `-f` — output structure format, `pdb` or `cif` (default:
  same as the input file).
- `--output`, `-o` — output path (default `<stem>_relettered.<format>`).

### Example

```bash
$ histo-chain-relettering 8gvi_1_aligned.cif \
    --map H="MHC alpha" \
    --map L="Beta-2 microglobulin" \
    --map P=Peptide \
    --map A="TCR alpha" \
    --map B="TCR beta"
             Chain relettering: 8gvi_1_aligned.cif
┏━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━┓
┃ source chain ┃ role                  ┃ target letter ┃
┡━━━━━━━━━━━━━╇━━━━━━━━━━━━━━━━━━━━━━━╇━━━━━━━━━━━━━━━┩
│ H            │ MHC alpha             │ A             │
│ L            │ Beta-2 microglobulin  │ B             │
│ P            │ Peptide               │ C             │
│ A            │ TCR alpha             │ D             │
│ B            │ TCR beta              │ E             │
└──────────────┴───────────────────────┴───────────────┘
Wrote relettered structure to 8gvi_1_aligned_relettered.cif
```

Note in this example that source chain `H` (MHC alpha) and source chain
`A` (TCR alpha) both need relettering, and their target letters (`A`
and `D` respectively) overlap with *other* source chain ids already in
the structure — this is exactly the case the two-step process (below)
exists for.

## Library usage

```python
from histo_chain_relettering import ChainReletterer, reletter_chains, save_structure

reletterer = ChainReletterer("8gvi_1_aligned.cif")
reletterer.resolve_role("Peptide")            # -> "C"

result = reletterer.reletter({
    "H": "MHC alpha",
    "L": "Beta-2 microglobulin",
    "P": "Peptide",
    "A": "TCR alpha",
    "B": "TCR beta",
})
result.mapping             # -> {"H": "A", "L": "B", "P": "C", "A": "D", "B": "E"}
result.unmapped_chains     # -> []
save_structure(result.structure, "out.cif", fmt="cif")
```

A one-shot convenience function is also available:

```python
from histo_chain_relettering import reletter_chains

result = reletter_chains("8gvi_1_aligned.cif", {"P": "Peptide"})
```

## Notes and limitations

- Relettering happens in **two steps** to avoid id collisions: each
  mapped chain is first renamed to its own id repeated 4 times (`P` ->
  `PPPP`), then every proxy is renamed to its resolved target letter
  (`PPPP` -> `C`). This is necessary because a chain's target letter is
  often already in use by *another* chain in the structure that hasn't
  been relettered yet (e.g. 8gvi's chain `A` is TCR alpha, but the
  target letter for MHC alpha is also `A`). Both steps happen entirely
  in memory before anything is written.
- Chains not mentioned in the mapping are left completely unchanged.
- Two source chains resolving to the same target letter, or a target
  letter colliding with a chain deliberately left unmapped, both raise
  an error rather than silently overwriting a chain.
- Only the **first model** in a file is used.
- The standardized role/letter scheme (`chain_letters.json`) is bundled
  with the package and is not overridable via the CLI — it's the
  family-wide standard, not a per-run configuration knob.

## Development

```bash
uv sync
uv run pytest
```

The test fixture, `tests/fixtures/8gvi_1_aligned.cif`, is copied from
`histo_com`'s own fixtures (same underlying source; a genuine pMHC-TCR
complex with chains `A`/`B`/`H`/`L`/`P`).

See [docs/PLAN.md](docs/PLAN.md) for the design rationale and
[CHANGELOG.md](CHANGELOG.md) for release history.
