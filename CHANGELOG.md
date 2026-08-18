# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

## 0.2.0

### Changed

- **`chain_letters.json` is keyed by chain type, not by letter.** Was
  `{"A": ["MHC alpha"]}` over 5 letters and prose role names; now
  `{"class_i_alpha": ["A"]}` over 56 chain types, the vocabulary the histo
  pipelines already speak. A caller with a `chain_types.json` can pass it
  through unchanged instead of translating at the boundary.
- `ChainReletterer.resolve_role(role)` -> `letters_for(chain_type)`, returning
  the whole list rather than one letter.
- `reletter()` takes `{chain_id: chain_type}`; `RelabelResult.roles` is now
  `.chain_types`.
- CLI `--map CHAIN=ROLE` -> `--map CHAIN=CHAIN_TYPE`.

### Added

- **Multi-copy allocation.** The nth chain of a type takes the nth letter, in
  file order. An assembly with two CD8 alphas gets `F` and `G`; previously a
  second chain of a type was a collision error.
- `RelabelResult.records` — the mapping as
  `{chain_type, consistent_chain_letter, pdb_chain_letter}` per chain, the
  shape the histo pipelines publish.
- `load_chain_letters(path)` and `ChainReletterer(path, chain_letters_path)`
  accept an external scheme, so a pipeline can supply its own.

## [Unreleased]

## [0.1.0] - 2026-07-07

### Added

- `histo_chain_relettering` Python library: `ChainReletterer` class that
  parses a PDB/mmCIF structure and the bundled standardized
  `chain_letters.json` role->letter scheme once, resolves an explicit
  `{source_chain_id: role_name}` mapping to target letters
  (`resolve_role()`, case-insensitive/whitespace-trimmed exact match),
  and relabels chains via `.reletter()`.
- Relettering is two-step to avoid id collisions: each mapped chain is
  first renamed to its own id repeated 4 times (`P` -> `PPPP`), then
  every proxy is renamed to its resolved target letter (`PPPP` -> `C`) —
  required because Biopython's `Bio.PDB.Entity.id` setter silently
  overwrites a colliding sibling id rather than raising.
- Validation before any mutation: two source chains resolving to the
  same target letter, or a target letter colliding with a chain
  deliberately left unmapped, both raise `RelettererError`.
- Chains not present in the mapping are left completely unchanged.
- Written chain order always follows chain id (`A`, `B`, `C`, ...), not
  the original file's/rename order — `.reletter()` sorts the in-memory
  model's chain list before returning, covering mapped and unmapped
  chains alike.
- `--format {pdb,cif}` output format, defaulting to the **input file's
  own format** (not a fixed default).
- `histo-chain-relettering` CLI (Click-based, Rich console table output)
  with `--map` (required, repeatable `CHAIN=ROLE`), `--format`, and
  `--output` options.
- Claude Code / Claude Desktop skill
  (`skills/histo-chain-relettering/`) wrapping the CLI.
- Test suite (pytest) against a committed real pMHC-TCR fixture
  (`8gvi_1_aligned.cif`, copied from `histo_com`'s fixtures) that
  exercises the full standardized letter scheme and the id-collision
  case directly.
- `README.md`, `CLAUDE.md`, and design plan (`docs/PLAN.md`).
- Bundled reference data moved from the repo-root `data/` directory to
  `src/histo_chain_relettering/data/chain_letters.json` so it ships
  with the wheel, per `CONVENTIONS.md` §3.
