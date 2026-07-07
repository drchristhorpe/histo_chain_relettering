---
name: histo-chain-relettering
description: Reletter a structure file's (PDB/mmCIF) chains to the standardized single-letter biological-role scheme (A=MHC alpha, B=MHC beta/Beta-2 microglobulin, C=Peptide, D=TCR alpha, E=TCR beta), given which role each chain plays. Use when asked to reletter, relabel, rename, or standardize a structure's chain IDs, or to convert a structure's chains to the standard pMHC/TCR letter scheme.
---

# histo-chain-relettering

`histo-chain-relettering` is a CLI tool (installed from the
`histo_chain_relettering` package) that relabels a PDB/mmCIF structure's
chain IDs to the family's standardized scheme, using Biopython. Invoke
it with the Bash tool.

## When to use this skill

The user provides (or references) a `.cif`/`.mmcif` or `.pdb`/`.ent`
structure file and asks to reletter/relabel/standardize its chain IDs —
typically to the pMHC/TCR scheme where MHC alpha is `A`, MHC beta/beta-
2 microglobulin is `B`, peptide is `C`, TCR alpha is `D`, and TCR beta
is `E`.

## Checking availability

```bash
histo-chain-relettering --help
```

If this fails with "command not found", install it first:

```bash
uv tool install histo_chain_relettering   # or: pip install histo_chain_relettering
```

(If working from a checkout of the `histo_chain_relettering` source
repo instead of an installed package, use `uv run
histo-chain-relettering ...` there instead.)

## Usage

```
histo-chain-relettering FILENAME --map CHAIN=ROLE [--map CHAIN=ROLE ...] [--format pdb|cif] [--output PATH]
```

- `FILENAME`: the structure file to reletter.
- `--map`/`-m` (required, repeatable): `CHAIN=ROLE`, e.g. `-m
  H="MHC alpha"`. `ROLE` must be one of the exact role strings in the
  standardized scheme: `MHC alpha`, `MHC beta`, `Beta-2 microglobulin`,
  `Peptide`, `TCR alpha`, `TCR beta` (case-insensitive). If you don't
  already know which chain in the file plays which role, inspect the
  structure first (e.g. via `histo-com` or by reading the file) —
  this tool doesn't infer roles automatically.
- `--format`/`-f` (optional): output structure format, `pdb` or `cif`.
  Defaults to the same format as the input file.
- `--output`/`-o` (optional): output path. Defaults to
  `<stem>_relettered.<format>`.

Chains not mentioned in `--map` are left completely unchanged — you
don't need to map every chain in the file, only the ones you want
relettered.

Example: relettering a pMHC-TCR complex where `H`=MHC alpha, `L`=beta-2
microglobulin, `P`=peptide, `A`=TCR alpha, `B`=TCR beta:

```bash
histo-chain-relettering structure.cif \
  --map H="MHC alpha" \
  --map L="Beta-2 microglobulin" \
  --map P=Peptide \
  --map A="TCR alpha" \
  --map B="TCR beta"
```

This produces `structure_relettered.cif` with chains renamed to
`A`/`B`/`C`/`D`/`E` respectively (MHC alpha -> `A`, beta-2 microglobulin
-> `B`, peptide -> `C`, TCR alpha -> `D`, TCR beta -> `E`).

## Output

One relettered structure file, plus a Rich console table (source chain
id / role / resolved target letter) and a note of any chains left
unmapped.

Report the result back to the user in whatever form they asked for
(the output file path, the resolved mapping, etc.) — this skill only
tells you how to obtain the relettered structure. An unknown role name
or a chain-id collision (e.g. two source chains resolving to the same
target letter) errors out with a clear message rather than silently
producing a wrong structure.
