# CLAUDE.md

Guidance for Claude Code when working in this repository. Follows the
shared conventions in [../CONVENTIONS.md](../CONVENTIONS.md).

## What this is

`histo_chain_relettering` reletters a structure file's (PDB/mmCIF)
chain IDs to the family's standardized single-letter biological-role
scheme, bundled at `src/histo_chain_relettering/data/chain_letters.json`.
It's a Python library, a Click-based CLI (`histo-chain-relettering`)
with Rich console output, and a Claude skill wrapping that CLI
(`skills/histo-chain-relettering/`). See [README.md](README.md) for
user-facing usage and [docs/PLAN.md](docs/PLAN.md) for the design
rationale.

## Environment

- Python 3.14, managed with `uv`. Use `uv sync`, `uv run <cmd>`, `uv run pytest`.
- Don't invoke a bare `python`/`pip` — always go through `uv run` /
  `uv add` so the lockfile stays authoritative.

## Layout

```
src/histo_chain_relettering/
  core.py   # load/save_structure, ChainReletterer, reletter_chains()
  cli.py    # Click CLI (entry point: histo-chain-relettering)
  data/
    chain_letters.json   # bundled standardized role -> letter scheme
tests/
  fixtures/8gvi_1_aligned.cif   # real pMHC-TCR complex (chains A/B/H/L/P)
  test_core.py
  test_cli.py
skills/histo-chain-relettering/SKILL.md
```

## Key invariants — don't break these

- **Relettering is always two-step, never a direct single-pass
  rename.** Checked directly against Biopython's `Bio.PDB.Entity.id`
  setter: assigning a chain id that's already used by a sibling chain
  does **not** raise — it emits a `BiopythonWarning` and silently
  overwrites the parent's `child_dict` entry, leaving the old chain
  object unreachable via `model[id]` even though it's still iterated by
  `child_list`. This must be prevented up front, not caught afterwards.
  Step 1 renames every mapped chain to its own id repeated 4 times
  (`P` -> `PPPP`); step 2 renames every proxy to its resolved target
  letter. A 4-character proxy can never collide with any single-
  character source or target id in this scheme, so ordering within
  each pass never matters. See `docs/PLAN.md` §5.
- Both rename passes happen **entirely in memory before any file is
  written** — the 4-character proxy ids are never serialized. This
  matters because classic PDB format only has a 1-character chain id
  column; only the final single-letter ids ever reach disk.
- Chain roles are resolved against the bundled
  `data/chain_letters.json` **only** — there's no CLI/library override.
  It's the family-wide standard, not a per-run config knob.
- After relettering, `.reletter()` sorts the model's `child_list` (the
  plain list `Bio.PDB.Entity.__iter__`/`PDBIO`/`MMCIFIO` all iterate for
  write order) by chain id — **every** chain, mapped or not, ends up in
  id order in the output file, not original/rename order. Only
  `child_list` is touched; `child_dict` (used for `model[id]` lookups)
  is a dict and unaffected by list order.
- Role name matching is case-insensitive and whitespace-trimmed, but
  otherwise exact — no fuzzy/abbreviation matching (e.g. `"MHC"` or
  `"B2M"` do **not** resolve; the caller must use the full role string
  from `chain_letters.json`, e.g. `"MHC alpha"`/`"Beta-2
  microglobulin"`). This was a deliberate scope decision, not an
  oversight — don't add fuzzy matching without checking first.
- Chains not mentioned in the `roles` mapping passed to `.reletter()`
  are left **completely untouched** — the tool never requires full
  chain coverage.
- `.reletter()` validates two collision cases *before* mutating
  anything: (1) two source chains resolving to the same target letter,
  and (2) a target letter colliding with a chain deliberately left
  unmapped. Both raise `RelettererError` rather than silently
  overwriting a chain (see the Biopython behavior above — silent
  corruption is exactly what these checks exist to prevent).
- Only the first model (`structure[0]`) is used everywhere — sufficient
  for X-ray/cryo-EM, a documented limitation for NMR ensembles.
- Output format defaults to the **input file's own format** (not a
  fixed default like `cif`) — `--format`/`fmt=` overrides explicitly.

## Testing

- `uv run pytest` — the fixture is already committed in
  `tests/fixtures/`; tests don't hit the network.
- `8gvi_1_aligned.cif` (copied from `histo_com/tests/fixtures/`) is a
  real pMHC-TCR complex with chains `A` (203 res, TCR alpha), `B` (241
  res, TCR beta), `H` (274 res, MHC alpha), `L` (99 res, beta-2
  microglobulin), `P` (8 res, peptide) — confirmed by inspection. It
  exercises the full standardized scheme in one structure and
  specifically exercises the collision case above: `A` is both a source
  chain id (TCR alpha) and the target letter for a different role (MHC
  alpha).

## Scope

The CLI intentionally exposes exactly three options: `--map`
(required, repeatable), `--format`, and `--output` (plus the positional
`FILENAME`). Don't add further options (e.g. a `chain_letters.json`
override flag, fuzzy/abbreviation role matching, batch processing of
multiple structures in one invocation) without checking with the user
first — these were deliberate constraints from the initial design
conversation, not oversights.
