# abd-prep

Prepare and validate input files for **Amazon Bio Discovery (ABD)** -- the
AWS antibody design service -- before you upload them:

| Input | What it is | Where ABD uses it | Command group |
|---|---|---|---|
| **PDB** | 3D structure of a protein (atom coordinates) | target antigen for *de novo* design, antibody framework / seed structure | `abd-prep pdb` |
| **FASTA** | amino-acid sequence as plain text | seed antibody for *directed evolution* and *property prediction* | `abd-prep fasta` |
| **CSV** | table of measurements (one row per antibody variant) | fine-tuning and results analysis with your own wet-lab data | `abd-prep csv` |

ABD rejects files that look right in a viewer but break its parser:
alternate atom locations (`PDB contains alternate locations`), stray
`TER` records, unlabelled heavy/light chains, Windows line ends, a `*`
at the end of a sequence, a UTF-8 byte-order mark written by Excel, a
`;`-separated "CSV". `abd-prep` fixes the *formatting* problems, reports
everything else as a finding with a stable code, and **never changes the
science**: no residue is renumbered or renamed, no sequence is truncated
or edited, no value is converted.

Python 3.8+, **standard library only**. Nothing to install: clone and
run `python3 -m abd_prep`.

```bash
# simple mode -- every file a first ABD run needs, in one step
python3 -m abd_prep quickstart
#   -> ./abd-quickstart/target.pdb            antigen (HER2 domain IV, 1N8Z chain C 489-607)
#      ./abd-quickstart/antibody_fv.pdb       seed antibody structure (trastuzumab Fv, chains L + H)
#      ./abd-quickstart/antibody_seed.fasta   the same seed as sequence (>H, >L)
#      ./abd-quickstart/assay_example.csv     assay-table template (placeholder row, not data)
#      ./abd-quickstart/README.txt            which file goes to which ABD recipe
#      ./abd-quickstart/quickstart.report.json

# the structure step on its own: download RCSB 1N8Z, cut HER2 domain IV + trastuzumab Fv, clean, validate
python3 -m abd_prep pdb split --preset her2-trastuzumab --out ./out
#   -> ./out/her2-trastuzumab_target.pdb   (antigen, chain C 489-607)
#      ./out/her2-trastuzumab_fv.pdb       (antibody Fv, chains L + H)
#      ./out/her2-trastuzumab.report.json  (what was done, SHA-256 of every file)
```

## Contents

- [Install](#install)
- [Quickstart: everything for a first run](#quickstart-everything-for-a-first-run)
- [Walkthrough: HER2 / trastuzumab](#walkthrough-her2--trastuzumab)
- [Which file does each ABD recipe need?](#which-file-does-each-abd-recipe-need)
- [Where do these files come from?](#where-do-these-files-come-from)
- [Command reference](#command-reference)
  - [`quickstart` -- simple mode](#quickstart--simple-mode)
  - [Common options and exit codes](#common-options-and-exit-codes)
  - [`pdb` -- structures](#pdb--structures)
  - [`fasta` -- seed antibody sequences](#fasta--seed-antibody-sequences)
  - [`csv` -- assay and results tables](#csv--assay-and-results-tables)
- [report.json](#reportjson)
- [Finding codes](#finding-codes)
- [Limitations](#limitations)
- [Conventions not yet confirmed against the ABD console](#conventions-not-yet-confirmed-against-the-abd-console)
- [Documentation](#documentation)
- [Development](#development)
- [Data sources and licenses](#data-sources-and-licenses)

## Install

No installation is required:

```bash
git clone https://github.com/hmkim/abd-prep.git
cd abd-prep
python3 -m abd_prep --help
```

To get the `abd-prep` command on your `PATH`:

```bash
pip install .                 # or: pipx install git+https://github.com/hmkim/abd-prep.git
abd-prep --version
```

The only network access is `pdb fetch`, `pdb split --preset/--pdb` and
`quickstart`, which download from RCSB (`files.rcsb.org`,
`data.rcsb.org`) over HTTPS into `~/.cache/abd-prep/` (override with
`--cache DIR` or `$ABD_PREP_CACHE`). `HTTPS_PROXY` is honoured. The
`fasta` and `csv` groups never touch the network.

## Quickstart: everything for a first run

`abd-prep quickstart` is the simple mode: one command, one folder, every
file a first ABD experiment needs, each already validated and each with
its purpose written down. It adds no rule of its own -- it runs `pdb
split`, `pdb inspect --fasta` + `fasta validate --against` and `csv
validate --against` for you and keeps their findings under the prefixes
`target/`, `fv/`, `seed/` and `csv/`.

```bash
python3 -m abd_prep quickstart                      # -> ./abd-quickstart/
python3 -m abd_prep quickstart --out ./my-test      # another folder
python3 -m abd_prep quickstart --preset her2-trastuzumab-6lbx-6b9z
python3 -m abd_prep quickstart --no-csv --strict    # skip the CSV template; warnings fail
```

```
[quickstart] her2-trastuzumab  ->  abd-quickstart
  WARNING PDB_RESOLUTION        1N8Z                       resolution 2.52 A is coarser than the 2.5 A the ABD guide recommends
  INFO    PDB_HETATM_REMOVED    target/file                removed HETATM records: HOH x79, NAG x28, SO4 x5
  INFO    PDB_GAP_NUMBERING     target/C:580-591           residue numbering jumps 580 -> 591 (expected)
  INFO    PDB_GAP_CA_DISTANCE   target/C:580-591           CA-CA distance 6.7 A between consecutive residues (> 4.2 A): chain break -- expected for this preset
  INFO    PDB_HETATM_REMOVED    fv/file                    removed HETATM records: HOH x79, NAG x28, SO4 x5
  INFO    FASTA_PDB_IDENTITY    seed/H:H<-antibody_fv.pdb  PDB chain H vs record 'antibody_fv.pdb': 120/120 aligned residues identical (PDB residue 1 = FASTA position 1)
  INFO    FASTA_PDB_IDENTITY    seed/L:L<-antibody_fv.pdb  PDB chain L vs record 'antibody_fv.pdb': 107/107 aligned residues identical (PDB residue 1 = FASTA position 1)
  INFO    CSV_CROSSREF_MATCHED  csv/against                1/1 row(s) match antibody_seed.fasta (0 by id, 1 by sequence)
  summary: 0 error(s), 1 warning(s), 7 info
  report: abd-quickstart/quickstart.report.json

Files in abd-quickstart:
  target.pdb           antigen structure (chain C 109 res)
                       -> target antigen -- de novo design (required); directed evolution with antigen (optional)
  antibody_fv.pdb      seed antibody structure (chain L 107 res, chain H 120 res)
                       -> seed antibody as a structure -- directed evolution, property prediction
  antibody_seed.fasta  seed antibody sequence (>H 120 res, >L 107 res)
                       -> seed antibody as a sequence -- directed evolution, property prediction (same antibody as antibody_fv.pdb)
  assay_example.csv    assay-table template (1 placeholder row, not a measurement)
                       -> assay / results table template -- fine-tuning, results analysis; replace the placeholder row with your data
  README.txt           what each file is, where it came from, what to do next
Next: read README.txt, then create an experiment in the ABD console and upload the file(s) its recipe asks for.
```

| File | What it is | ABD recipe that takes it |
|---|---|---|
| `target.pdb` | antigen structure, chain ID kept from the source entry (`C` for 1N8Z) | *de novo* design (required); *directed evolution* with antigen (optional) |
| `antibody_fv.pdb` | seed antibody structure, chains `L` + `H` | *directed evolution*, *property prediction* (seed as a structure) |
| `antibody_seed.fasta` | the same seed as sequence, records `>H` then `>L`, cross-checked against `antibody_fv.pdb` | *directed evolution*, *property prediction* (seed as a sequence) |
| `assay_example.csv` | **template**: the standard columns and one row -- the seed itself with a placeholder value, labelled as such in `notes` | fine-tuning / results analysis, after you replace the row with your measurements |
| `README.txt` | file list with purposes, RCSB source and revision date, next steps | keep next to the files |
| `quickstart.report.json` | the usual report: sources, SHA-256 of every file, every finding | keep next to the files |

Running it twice into the same folder stops with exit code 2 ("already
exist ... pass `--force`") so an edited seed file is never overwritten by
accident. The command needs a preset with a target and a heavy chain
(all shipped presets qualify); for your own structure use `pdb split
--input` and the `fasta` / `csv` commands directly.

## Walkthrough: HER2 / trastuzumab

This is the manual PyMOL procedure from the ABD walkthrough blog post
(download 1N8Z, keep the HER2 domain IV and the trastuzumab Fv, remove
alternate locations and waters, save two files), done in one command and
with a reproducibility record.

**1. Cut, clean and validate the structures**

```bash
python3 -m abd_prep pdb split --preset her2-trastuzumab --out ./out
```

```
[split] her2-trastuzumab  ->  target: ./out/her2-trastuzumab_target.pdb, fv: ./out/her2-trastuzumab_fv.pdb
  WARNING PDB_RESOLUTION       1N8Z              resolution 2.52 A is coarser than the 2.5 A the ABD guide recommends
  INFO    PDB_HETATM_REMOVED   target/file       removed HETATM records: HOH x79, NAG x28, SO4 x5
  INFO    PDB_GAP_NUMBERING    target/C:580-591  residue numbering jumps 580 -> 591 (expected)
  INFO    PDB_GAP_CA_DISTANCE  target/C:580-591  CA-CA distance 6.7 A between consecutive residues (> 4.2 A): chain break -- expected for this preset
  INFO    PDB_HETATM_REMOVED   fv/file           removed HETATM records: HOH x79, NAG x28, SO4 x5
  summary: 0 error(s), 1 warning(s), 4 info
  report: ./out/her2-trastuzumab.report.json
```

The preset encodes what the blog post does by hand: HER2 extracellular
domain IV is chain C residues 489-607 (mature-protein numbering), the
antibody variable fragment (Fv -- the part of an antibody that binds the
antigen) is light chain A 1-107 and heavy chain B 1-120, written as
chains `L` and `H`. The one warning is real: 1N8Z is a 2.52 A structure
and the ABD guide recommends 2.5 A or better.

Without a preset, spell out the selection yourself -- for any RCSB entry
or a local file:

```bash
python3 -m abd_prep pdb split --pdb 1N8Z --target C:489-607 --heavy B:1-120 --light A:1-107 \
    --rename-heavy H --rename-light L --out ./out
python3 -m abd_prep pdb split --input ./my_complex.pdb --target B:530-601 --out ./out
```

**2. Upload**

- `her2-trastuzumab_target.pdb` is the **target antigen** file (every
  recipe that takes an antigen).
- `her2-trastuzumab_fv.pdb` is the **antibody framework / seed
  structure** (directed evolution and property prediction accept a seed
  as PDB or FASTA).

`./out/her2-trastuzumab.report.json` records the RCSB revision date, the
exact command, every finding and the SHA-256 of each file, so the
experiment can be traced back to the exact input.

**3. Seed sequence as FASTA**

If you would rather upload the seed as a sequence, extract it from the
structure and confirm the two agree:

```bash
python3 -m abd_prep pdb inspect ./out/her2-trastuzumab_fv.pdb --fasta > seed.fasta
python3 -m abd_prep fasta validate seed.fasta --framework igg --against ./out/her2-trastuzumab_fv.pdb
```

Or clean a sequence you already have (from a sequencing report, a paper,
a database):

```bash
python3 -m abd_prep fasta clean trastuzumab_from_paper.fasta seed.fasta --framework igg
```

**4. Your own assay table**

```bash
python3 -m abd_prep csv validate results.csv --map id=Clone,sequence=AA_seq,value=KD,unit=Unit
python3 -m abd_prep csv clean results.csv results_clean.csv --map id=Clone,sequence=AA_seq,value=KD,unit=Unit
python3 -m abd_prep csv validate results_clean.csv --against seed.fasta     # do the rows refer to the seeds?
```

## Which file does each ABD recipe need?

ABD groups its workflows into *recipes*. The deciding question is whether
the recipe starts from an **existing antibody** (a *seed*) or not.

| Recipe | Target antigen | Seed antibody | Notes |
|---|---|---|---|
| **Antibody generation -- de novo** | required (PDB) | not used | pick a framework: Nanobody, IgG or scFv |
| **Antibody generation -- directed evolution** | optional (PDB) | required (FASTA or PDB) | with an antigen: binding-aware optimisation; without: optimise the seed on its own |
| **Property prediction** | optional | required (FASTA or PDB) | affinity, stability, developability predictions |

The *framework* decides how many chains the seed has, and `abd-prep fasta`
checks the file against it (`--framework`):

| Framework | What it is | Chains in the seed |
|---|---|---|
| **IgG** | the standard Y-shaped antibody; the variable region is two domains | one **heavy** (VH) + one **light** (VL) record |
| **Nanobody (VHH)** | a single-domain antibody derived from camelid heavy-chain-only antibodies | one record (heavy only) |
| **scFv** | VH and VL fused into one chain through a flexible linker (typically `(GGGGS)n`) | one record containing both domains and the linker |

Every upload must be smaller than 100 MB and encoded as UTF-8; the tool
checks both for every file.

## Where do these files come from?

**PDB** files come from the [RCSB Protein Data Bank](https://www.rcsb.org/)
(public crystal / cryo-EM structures) or from your own structure
determination or modelling. Public entries contain far more than ABD
wants: the whole antibody Fab or a second copy of the complex, waters,
ions, glycans and buffer molecules (`HETATM` records), and *alternate
locations* -- two sets of coordinates for one atom where the
crystallographer saw it in two positions. The `pdb` group selects the
chains and residue ranges you need and resolves the rest.

**FASTA** files are the antibody as text: one header line starting with
`>` and the amino-acid sequence in one-letter code. They come from a
sequencing report of a hit antibody (hybridoma or phage display, usually
delivered as FASTA or Excel by the discovery platform or CRO), from a
paper, patent or database (IMGT, SAbDab, UniProt, INN lists), or are
typed into a text editor. Typical problems: the header does not say which
record is the heavy and which the light chain; a signal peptide or the
constant region is still attached; alignment gaps, numbers or lower-case
letters are embedded; a translation stop `*` is left at the end; `X`, `B`
or `Z` stand in for unknown residues. The `fasta` group fixes the
formatting and points out the rest.

**CSV** files are the laboratory results: one row per clone or variant,
columns for the measured values (affinity KD, melting temperature Tm,
expression titre ...). They are exported from the instrument software
(SPR, BLI, ELISA readers), arrive from a CRO as Excel, or are assembled
by hand. Typical problems: every lab names its columns differently
(`KD`, `Kd (nM)`, `affinity`), units are mixed, missing values are
written a dozen ways, Excel saves UTF-8 with a byte-order mark and -- in
many locales -- with `;` instead of `,`. The `csv` group maps your
columns onto a standard set, normalises encoding and delimiter, and
reports what cannot be fixed mechanically (missing values, non-numeric
cells, mixed units).

## Command reference

Every group has `clean` (write a corrected copy), `validate` (check
without writing) and `inspect` (summarise). `pdb` adds `fetch` and
`split`. `quickstart` sits outside the groups because it spans all three
(see [above](#quickstart-everything-for-a-first-run)). Run
`abd-prep GROUP ACTION --help` for the full option list.

### `quickstart` -- simple mode

```bash
abd-prep quickstart [--preset NAME] [--out DIR] [--force] [--no-csv] [--cache DIR] [--refresh]
```

`--preset` defaults to `her2-trastuzumab`; `--out` to `./abd-quickstart`.
Writes `target.pdb`, `antibody_fv.pdb`, `antibody_seed.fasta`,
`assay_example.csv` (unless `--no-csv`), `README.txt` and
`quickstart.report.json`; refuses to overwrite existing files without
`--force`. Takes the common options below except `--encoding` (it reads
no user text file).

### Common options and exit codes

| Option | Meaning |
|---|---|
| `--json` | print the report as JSON instead of the text table |
| `--strict` | treat warnings as failures (exit code 1) -- for CI and scripts |
| `--encoding ENC` | input text encoding for FASTA/CSV (default `utf-8`; e.g. `cp949`, `latin-1`). `clean` converts to UTF-8 |
| `--report PATH` | write `report.json` to PATH |
| `--no-report` | do not write `report.json` (`clean` writes `<output stem>.report.json` by default; `validate` writes one only with `--report`) |
| `--report-paths basename|full` | record file paths as basename (default; reports travel in mail and tickets) or absolute |

Exit codes: **0** pass (warnings and info allowed), **1** at least one
validation error (or a warning under `--strict`), **2** input or option
error (unreadable file, unknown chain, bad range, `.xlsx` given where
CSV is expected). `clean` writes its output even when there are
validation errors, so you can look at the result; the exit code still
says 1. The one exception is an input that cannot be decoded at all.

### `pdb` -- structures

```bash
abd-prep pdb fetch    1N8Z 6LBX [--out DIR] [--cache DIR] [--refresh] [--no-meta]
abd-prep pdb split    --preset NAME | --pdb ID | --input FILE  [--target SPEC] [--heavy SPEC] [--light SPEC]
                      [--rename-target ID] [--rename-heavy ID] [--rename-light ID] [--name NAME] [--out DIR]
                      [--keep-het] [--std-residues] [--strip-h] [--renumber]
abd-prep pdb clean    IN.pdb OUT.pdb [--keep SPEC ...] [--rename A=L,B=H] [--keep-het] [--std-residues] [--strip-h] [--renumber]
abd-prep pdb validate FILE... [--preset NAME] [--expect-chains N] [--allow-het]
abd-prep pdb inspect  FILE... [--fasta]
```

**Selecting chains and residues.** `--keep`, `--target`, `--heavy` and
`--light` take `CHAIN`, `CHAIN:LO-HI`, `CHAIN:LO-HI,LO-HI` or a single
`CHAIN:RES`. Numbers are the residue numbers *as written in the file*
(PDB columns 23-27), including insertion codes, so Kabat-numbered
antibodies (`82A`, `82B`, `100A` ...) select correctly:

- a bound written **without** an insertion code includes every inserted
  residue of that number: `H:1-82` keeps 82, 82A, 82B and 82C;
- a bound written **with** one is exact: `H:1-82A` keeps 82 and 82A only.

Negative residue numbers (`-3`) are accepted.

**What `clean` does, in order**

1. Keeps the first `MODEL` only (NMR ensembles, multi-model files) and warns.
2. Applies the chain / residue selection.
3. Resolves **alternate locations**: for each atom the conformer with the
   highest occupancy is kept and column 17 is blanked (ties go to the
   first one in the file). When two *different residues* are modelled at
   one position (e.g. `ACSO` / `BCYS` at 6B9Z B175), the residue name with
   the higher occupancy is chosen first, and a conformer that the HETATM
   filter would remove never displaces one that would stay. This is a
   different rule from PyMOL's `remove not alt ''+A`, which keeps
   conformer `A` regardless of occupancy; results can differ in the
   affected residues.
4. Removes **HETATM** records (water, ions, glycans, ligands) unless
   `--keep-het`. A HETATM residue that is part of the protein backbone
   (has N, CA and C and sits inside the chain -- selenomethionine `MSE`,
   oxidised cysteine `CSO`, pyroglutamate `PCA` ...) would leave a hole in
   the chain, so its removal is a **warning**. With `--std-residues` such
   residues are converted to their standard parent instead:

   | Modified | Written as | Atoms |
   |---|---|---|
   | `MSE` | `MET` | `SE` renamed `SD`, element column updated |
   | `CSO`, `CSD`, `OCS`, `CME` | `CYS` | oxygen / adduct atoms dropped |
   | `SEP`, `TPO`, `PTR` | `SER`, `THR`, `TYR` | phosphate dropped |
   | `HYP` | `PRO` | `OD1` dropped |
   | `MLY`, `KCX` | `LYS` | methyl / carboxyl atoms dropped |
   | `PCA` | `GLN` | `OE` renamed `OE1`; side chain stays incomplete |

   Any atom not in the standard residue's atom set is dropped and listed
   in the report. Modified residues without a rule are removed with a
   warning (`--keep-het` retains them).
5. Removes hydrogen and deuterium atoms with `--strip-h` (element column,
   falling back to the atom name when the column is blank).
6. Renames chains (`--rename A=L,B=H`). Renaming is applied line by line
   at write time, so swapping two IDs (`A=B,B=A`) is safe. Two source
   chains mapped to one output ID is an error (`PDB_RENAME_COLLISION`).
7. Regenerates **`TER` and `END`**: one `TER` per chain boundary naming
   the last residue, original `TER` records are discarded (a leftover
   `TER` from a removed chain is a common reason strict parsers reject
   edited files). Atom serial numbers are kept as in the source unless
   `--renumber`; `validate` warns about duplicates.

Everything else in `ATOM`/`HETATM` records is copied byte for byte.
`HEADER`, `REMARK`, `SEQRES`, `SSBOND`, `LINK`, `CRYST1` and all other
record types are **not** written (see [Limitations](#limitations)).

**`split`** is `fetch` + `clean` + `validate` producing the two files ABD
wants: `<name>_target.pdb` (one chain) and `<name>_fv.pdb` (heavy +
light), plus `<name>.report.json` with the RCSB entry metadata (revision
date, resolution, method). Target and antibody may come from different
entries (`--pdb` for one entry; a preset can combine two).

**Presets**

| Preset | Entry | Target | Fv | Notes |
|---|---|---|---|---|
| `her2-trastuzumab` (default) | 1N8Z (HER2 ECD / Herceptin Fab, 2.52 A) | `C:489-607` -> chain `C`, 109 residues | `A:1-107` -> `L`, `B:1-120` -> `H` | Matches the blog walkthrough and the reference trastuzumab sequence. Chain C is numbered as the mature protein (UniProt P04626 numbering minus 22). Loop 581-590 is not modelled in 1N8Z; `validate` reports it as an expected gap. |
| `her2-trastuzumab-6lbx-6b9z` | 6LBX (target) + 6B9Z (Fv) | `B:530-601` -> chain `B`, 72 residues | `A:1-107` -> `L`, `B:1-120` -> `H` | The 6B9Z antibody differs from trastuzumab at several framework positions; treat it as a trastuzumab-derived variant. Target and antibody come from different crystals, so their coordinate frames are unrelated. |

Presets also carry expectations (residue counts, known gaps) that
`validate --preset NAME` applies to any file whose chain IDs match a
preset part. The target chain keeps its original ID in both presets
(`--rename-target` changes it) until the ABD chain-label requirement is
confirmed -- see [Conventions not yet confirmed](#conventions-not-yet-confirmed-against-the-abd-console).

**`validate`** checks what ABD's parser trips over: alternate locations,
`HETATM` records, `TER`/`END` placement, chain count, residue counts
plausible for a VL (105-115) and VH (110-125), chain breaks (consecutive
residues with a CA-CA distance above 4.2 A -- a peptide bond gives about
3.8 A), residue-number jumps, non-standard residue names, duplicate
serials, non-ASCII bytes. **`inspect`** prints chains, ranges and
sequences; `--fasta` writes them as FASTA with the ABD heavy/light
headers so the output can go straight into `fasta validate`.

### `fasta` -- seed antibody sequences

```bash
abd-prep fasta clean    IN.fasta OUT.fasta --framework igg|nanobody|scfv [--heavy-id ID] [--light-id ID]
                        [--replace-ambiguous X=G,...] [--no-heuristics] [--wrap N]
abd-prep fasta validate FILE... [--framework ...] [--heavy-id ID] [--light-id ID] [--replace-ambiguous ...]
                        [--no-heuristics] [--against FV.pdb] [--chain-map H=B,L=A]
abd-prep fasta inspect  FILE... [--framework ...] [--heavy-id ID] [--light-id ID]
```

**Which record is heavy, which is light?** The tool is deliberately
conservative: uploading the two chains swapped is the most expensive
mistake, so it never guesses from sequence content. A role is assigned
from, in priority order:

1. `--heavy-id ID` / `--light-id ID` -- the record ID (first token of the
   header line);
2. the framework when it implies a single record: `nanobody` -> heavy,
   `scfv` -> one record holding both domains;
3. header tokens, split on anything that is not a letter or digit (so
   RCSB headers such as `1N8Z_2|Chains B|... heavy chain` work):
   heavy -- `H`, `HC`, `VH`, `VHH`, `heavy`, `heavychain`, `IGH`, `nanobody`,
   `Nb`, `sdAb`, or an IMGT gene name starting `IGHV`/`IGHG`/...; light --
   `L`, `LC`, `VL`, `VK`, `light`, `lightchain`, `kappa`, `lambda`, `IGK`,
   `IGL`, or `IGKV`/`IGLV`/...; scFv -- `scFv`, `singlechain`.

A header with both heavy and light tokens is ambiguous (unless it also
says `scFv`), and a record whose role cannot be told is an **error** for
`--framework igg` and a warning otherwise.

**What `clean` changes -- formatting only**

- CR/LF line ends -> LF; UTF-8 without BOM; `--encoding` converts other encodings.
- Headers are rewritten to the ABD convention with the original header
  preserved as the description: `>H heavy trastuzumab VH ...`,
  `>L VL_1 kappa ...`, `>scFv ...`.
- Blanks, digits (from numbered alignments), alignment gaps `-` and `.`
  are removed from the sequence; lower-case letters are upper-cased.
- A terminal `*` (translation stop) is removed. A `*` **inside** the
  sequence is an **error**, not removed: it means the translation hit an
  internal stop codon, and the reading frame should be checked instead.
- Ambiguous letters (`X`, `B`, `Z`, `J`, `U`, `O`) are never replaced
  automatically; `--replace-ambiguous X=G,B=N` applies an explicit
  substitution and records it. (`validate` applies it in memory only; a
  file is written by `clean` alone.)

The sequence is otherwise written exactly as read. Nothing is truncated
to the variable domain and no residue is changed.

**What `validate` checks**

- Alphabet: only the 20 standard amino acids; a record made solely of
  A/C/G/T/U/N (20 residues or more) is reported as a nucleotide sequence.
- Framework consistency: `igg` needs exactly one heavy and one light
  record, `nanobody` and `scfv` exactly one record.
- Heuristics (warnings only, `--no-heuristics` turns them off): length
  outside the typical range for the role -- VL 95-130, VH 105-140, VHH
  100-140, scFv 220-280 (deliberately wider than the structure-side VL
  105-115 / VH 110-125, because a seed sequence may legitimately carry a
  few extra residues); a leading `M` followed by a hydrophobic stretch
  (signal peptide not removed?); two records with the same sequence. For
  scFv records the `(G3-4S)n` linker is located and reported as info.
- `--against FV.pdb` (**the same molecule?**): compares each record with
  the corresponding PDB chain, position by position, using the PDB
  residue numbers. Residues missing from the structure (unmodelled loops)
  are *info*; the seed being longer than the structure at either end
  (signal peptide, constant region) is *info*; the seed being shorter is
  a *warning*; a residue that differs is an **error**. The alignment
  offset is found by scanning every start position, so a full-length
  sequence still lines up with an Fv structure and a real point mutation
  is still reported. Below 50 % identity the record and chain are
  declared non-corresponding in one line. PDB chains labelled `H`/`L`
  map to the roles automatically; otherwise give `--chain-map H=B,L=A`.

### `csv` -- assay and results tables

```bash
abd-prep csv clean    IN.csv OUT.csv [--map STD=HEADER,...] [--delimiter CHAR] [--drop-missing] [--drop-unmapped]
abd-prep csv validate FILE... [--map STD=HEADER,...] [--delimiter CHAR] [--against seed.fasta|fv.pdb]
abd-prep csv inspect  FILE... [--map STD=HEADER,...] [--delimiter CHAR]
```

**Standard columns.** Lab tables all name their columns differently, so
`--map` projects your headers onto a standard set; a header that already
equals a standard name is mapped automatically.

| Standard column | Content | Requirement |
|---|---|---|
| `id` | clone / variant identifier | **required** (rows must be identifiable) |
| `sequence` | full antibody sequence -- or | one of the two forms is **recommended**; without any, `--against` matches by `id` only |
| `heavy` + `light` | the two chains in separate columns | |
| `value` | the measured number (KD, Tm, titre ...) | optional |
| `unit` | unit of `value` (nM, pM, degC ...) | optional; mixing is reported, never converted |
| `target`, `assay`, `notes` | antigen name, assay type (SPR, BLI, ELISA ...), free text | optional |

```bash
abd-prep csv clean results.csv results_clean.csv --map id=Clone,sequence=AA_seq,value=KD,unit=Unit
abd-prep csv clean results.csv results_clean.csv --map id=Clone --map "value=KD, nM"   # header containing a comma
```

Unmapped columns pass through unchanged after the standard ones
(`--drop-unmapped` removes them). A `--map` that names a header not in
the file is an option error (exit 2) listing the real headers.

**What `clean` changes -- encoding, delimiter and header only**

- Encoding: UTF-8 without BOM. Other encodings are converted when named
  with `--encoding` (e.g. `cp949`, `latin-1`); without it a non-UTF-8 file
  is an error that names the byte position.
- Delimiter: `,` in the output. The input delimiter (`,` `;` `|` tab) is
  detected with `csv.Sniffer`; if detection fails and the header line
  contains exactly one candidate it is used, otherwise `,` is assumed with
  a warning. `--delimiter ';'` (or `tab`) overrides.
- Header: standard columns renamed to their standard names, in standard
  order; quoting normalised to the minimum.
- Rows whose cell count differs from the header (an unquoted delimiter
  inside a value is the usual cause) are reported as errors with their
  line numbers and are **not** written.

Cell values are written exactly as read. No number is reformatted, no
unit converted, no sequence case-folded.

**What `validate` checks**

Header present, no empty or duplicate names; `id` present; a sequence
form present (warning); missing values per standard column (`''`, `NA`,
`N/A`, `NaN`, `None`, `null`, `-`, `--`, `ND`, `N.D.`, case-insensitive --
`--drop-missing` removes those rows); `value` cells that are not plain
numbers (`<0.1`, `>1000`, text) and cells with a comma inside the number
(`1,2` -- European decimal mark or thousands separator?); mixed units;
duplicate ids; sequence columns checked with the FASTA alphabet rules.
`--against seed.fasta` or `--against fv.pdb` cross-references the rows:
a row matches when its `id` equals a record / chain ID **or** its
sequence (or `heavy`+`light`) equals a record / chain sequence; unmatched
rows are listed.

## report.json

Every `clean` (and `split`) writes a report next to its output, and
`validate` writes one with `--report`. One schema for every group:

```json
{
  "tool":    {"name": "abd-prep", "version": "0.1.0", "python": "3.9.25"},
  "command": "abd-prep pdb clean 1N8Z_mini.pdb fv.pdb --keep A --keep B --rename A=L,B=H",
  "generated_at": "2026-10-26T05:53:28Z",
  "inputs":  [{"kind": "pdb", "path": "1N8Z_mini.pdb", "sha256": "9a9d571b...", "bytes": 10555}],
  "outputs": [{"kind": "pdb", "path": "fv.pdb", "sha256": "419c6ec4...", "bytes": 6723}],
  "findings": [
    {"level": "info", "code": "PDB_ALTLOC_RESOLVED", "where": "file",
     "message": "2 atoms in 2 residues had alternate locations; kept the highest-occupancy conformer, dropped 2 atoms"},
    {"level": "info", "code": "PDB_HETATM_REMOVED", "where": "file", "message": "removed HETATM records: HOH x2"}
  ],
  "details": {"atoms_in": 123, "atoms_out": 80, "hetatm_removed": {"HOH": 2}, "rename": {"A": "L", "B": "H"},
              "chains_out": {"L": {"residues": 5, "atoms": 40, "first": "1", "last": "5"}, "H": {"...": "..."}},
              "ter_records_written": 2}
}
```

`details` is group specific: `split` adds `sources[]` (PDB ID, RCSB
revision date, resolution, method) and `parts{target, fv}`; `fasta`
adds `records[]` and, with `--against`, per-chain `offset`, `identical`,
`mismatches` and `pdb_gaps`; `csv` adds the column mapping, unmapped
columns, missing-value counts, units, `rows_read` / `rows_written` and
the `--against` match counts. Paths are basenames unless
`--report-paths full`.

## Finding codes

Codes are stable identifiers meant for scripts and for searching this
README. **error** fails the run (exit 1); **warning** is reported, and
fails only under `--strict`; **info** records what was done.

**File level (all groups)**

| Code | Level | Meaning |
|---|---|---|
| `FILE_TOO_LARGE` | error | 100 MB or more; ABD rejects the upload |
| `FILE_NOT_UTF8` | error | not valid UTF-8 and no `--encoding` given |
| `FILE_ENCODING` | error | does not decode with the `--encoding` given |
| `FILE_ENCODING_INPUT` | warning | (`validate`) decodes as the given encoding but ABD needs UTF-8; run `clean --encoding` |
| `FILE_ENCODING_CONVERTED` | info | (`clean`) input decoded as the given encoding, output written as UTF-8 |
| `FILE_ENCODING_ASCII` | info | `--encoding` given but the file is plain ASCII (already UTF-8) |
| `FILE_BOM` | info | UTF-8 byte-order mark present (Excel); `clean` removes it |
| `FILE_EXTENSION` | warning | unusual extension for the group; the file is still parsed (`.xlsx` is an input error with a "save as CSV UTF-8" hint) |

**`pdb`**

| Code | Level | Meaning |
|---|---|---|
| `PDB_NO_ATOMS` | error | no `ATOM`/`HETATM` records |
| `PDB_ALTLOC_PRESENT` | error | column 17 not blank -- the "contains alternate locations" rejection |
| `PDB_END_MISSING` | error | no `END` record |
| `PDB_CHAIN_COUNT` | error | chain count differs from `--expect-chains` |
| `PDB_TER_MISSING` / `PDB_TER_MISPLACED` / `PDB_TER_MISMATCH` / `PDB_TER_EXTRA` | error | `TER` absent at a chain boundary / left over from a removed chain / naming a different residue than the atom before it / several at one boundary |
| `PDB_CHAIN_NOT_FOUND` | error | a selected chain is not in the file |
| `PDB_SELECTION_EMPTY` | error | the range selected no atoms (or only `HETATM` records that were removed -- `--keep-het`) |
| `PDB_RENAME_COLLISION` | error | two source chains map to one output chain ID |
| `PDB_HETATM_PRESENT` | warning / info | `HETATM` records remain (info with `--allow-het`) |
| `PDB_GAP_CA_DISTANCE` | warning / info | consecutive residues with CA-CA > 4.2 A -- a chain break (info when a preset lists it as expected) |
| `PDB_MODIFIED_RESIDUE_REMOVED` | warning | a modified backbone residue inside the selection was removed with the `HETATM` records -> gap; use `--std-residues` |
| `PDB_MODIFIED_RESIDUE_NO_RULE` | warning | as above but no conversion rule exists; `--keep-het` retains it |
| `PDB_NONSTANDARD_RESNAME` | warning | `ATOM` record with a non-standard residue name |
| `PDB_RESIDUE_COUNT` | warning | chain H/L (or a preset expectation) outside the plausible Fv size |
| `PDB_RESIDUE_DUPLICATE` | warning | same residue number twice in a chain with different names |
| `PDB_MULTI_MODEL` | warning | more than one `MODEL`; `clean` keeps the first |
| `PDB_CHAIN_BLANK` | warning | blank chain ID |
| `PDB_SERIAL_DUPLICATE` | warning | duplicate atom/`TER` serial numbers (`--renumber`) |
| `PDB_NON_ASCII` | warning | `ATOM`/`HETATM` record with non-ASCII bytes (written unchanged) |
| `PDB_RESOLUTION` | warning | (`split`) entry resolution coarser than the 2.5 A the ABD guide recommends |
| `PDB_EXTENSION` | warning | file name does not end in `.pdb`/`.ent` |
| `PDB_ALTLOC_RESOLVED` | info | alternate locations resolved: atoms affected, atoms dropped |
| `PDB_RESIDUE_HETEROGENEITY` | info | two residue types modelled at one position; which one was kept |
| `PDB_MODIFIED_CONFORMER_REMOVED` | info | a modified-residue conformer was removed while the standard conformer at that position stays (no gap) |
| `PDB_MODIFIED_RESIDUE_CONVERTED` | info | `--std-residues` conversion applied (atoms renamed / dropped) |
| `PDB_HETATM_REMOVED` | info | `HETATM` records removed, by residue name |
| `PDB_HYDROGENS_REMOVED` / `PDB_HYDROGENS_PRESENT` | info | hydrogens removed by `--strip-h` / present in the file |
| `PDB_GAP_NUMBERING` | info | residue numbers skip or decrease (Kabat codes, unmodelled loops, chimeras) |
| `PDB_PRESET_NO_MATCH` | info | chain IDs match no part of the preset; generic rules applied |
| `PDB_ENTRY_META_UNAVAILABLE` | info | RCSB entry metadata could not be fetched; files are unaffected |

**`fasta`**

| Code | Level | Meaning |
|---|---|---|
| `FASTA_NO_RECORDS` | error | no line starts with `>` |
| `FASTA_NO_HEADER` | error | sequence text before the first header |
| `FASTA_EMPTY_HEADER` | error | header line without a record ID |
| `FASTA_EMPTY_SEQUENCE` | error | header without sequence |
| `FASTA_TOO_LONG` | error | more than 10 000 residues -- not an antibody construct |
| `FASTA_INTERNAL_STOP` | error | `*` inside the sequence: internal stop codon, check the reading frame |
| `FASTA_BAD_ALPHABET` | error | letters outside the 20 standard amino acids (positions listed); `--replace-ambiguous` to substitute explicitly |
| `FASTA_NUCLEIC` | error | A/C/G/T/U/N only: a DNA/RNA sequence, not the translated protein |
| `FASTA_ROLE_UNKNOWN` | error / warning | heavy or light cannot be told from the header (error under `--framework igg`) |
| `FASTA_ROLE_AMBIGUOUS` | error | header carries both heavy and light tokens |
| `FASTA_ROLE_CONFLICT` | error | a light-chain token in a nanobody (single heavy domain) file |
| `FASTA_FRAMEWORK_MISMATCH` | error | record count / roles do not fit `--framework` |
| `FASTA_PDB_MISMATCH` | error | (`--against`) residues that differ between record and PDB chain |
| `FASTA_PDB_NO_CORRESPONDENCE` | error | (`--against`) record and chain do not correspond (identity below 50 %) |
| `FASTA_PDB_NOT_COMPARED` | error | (`--against`) no record carries a role, nothing was compared |
| `FASTA_DUPLICATE_ID` | warning | the same record ID appears more than once |
| `FASTA_LENGTH_RANGE` | warning | length outside the typical range for the role (truncated? constant region or signal peptide included?) |
| `FASTA_SIGNAL_PEPTIDE` | warning | starts with `M` + hydrophobic stretch; mature variable domains usually start with E/Q/D |
| `FASTA_DUPLICATE_SEQUENCE` | warning | two records with the same sequence |
| `FASTA_PDB_UNCOVERED` | warning | (`--against`) PDB residues outside the record: the seed is shorter than the structure |
| `FASTA_PDB_ROLE_MISSING` | warning | (`--against`) a PDB chain is mapped to a role no record has |
| `FASTA_CRLF` / `FASTA_LOWERCASE` / `FASTA_WHITESPACE_REMOVED` / `FASTA_GAPS_REMOVED` / `FASTA_TRAILING_STOP` | info | formatting fixed by `clean` |
| `FASTA_AMBIGUOUS_REPLACED` | info | an explicit `--replace-ambiguous` substitution was applied |
| `FASTA_ROLE_OVERRIDDEN` | info | header suggests one role but `--heavy-id/--light-id` assigns another |
| `FASTA_FRAMEWORK_UNSPECIFIED` | info | `--framework` not given; count/role rules skipped |
| `FASTA_SCFV_LINKER` | info | linker located (or no common pattern found; sequence unchanged) |
| `FASTA_HEURISTICS_OFF` | info | `--no-heuristics` given |
| `FASTA_PDB_IDENTITY` / `FASTA_PDB_GAP` / `FASTA_EXTRA_RESIDUES` | info | (`--against`) identity summary / residues unmodelled in the structure / residues the record has beyond the structure at either end |

**`csv`**

| Code | Level | Meaning |
|---|---|---|
| `CSV_NO_HEADER` | error | empty file or no header row |
| `CSV_HEADER_EMPTY_NAME` / `CSV_HEADER_DUPLICATE` | error | header column without a name / duplicate header names |
| `CSV_ROW_LENGTH` | error | rows whose cell count differs from the header (not written by `clean`) |
| `CSV_ID_MISSING` | error | no `id` column after mapping |
| `CSV_SEQUENCE_ALPHABET` / `CSV_SEQUENCE_NUCLEIC` | error | sequence cells with non-standard letters / that look like DNA |
| `CSV_SEQUENCE_MISSING` | warning | no `sequence` (or `heavy`+`light`) column; `--against` will match by id only |
| `CSV_DELIMITER_GUESSED` | warning | delimiter not detected, `,` assumed |
| `CSV_NO_DATA_ROWS` | warning | header only |
| `CSV_MISSING_VALUES` | warning | missing values per standard column (`--drop-missing`) |
| `CSV_VALUE_NOT_NUMERIC` | warning | `value` cells that are not plain numbers (`<0.1`, text); kept as written |
| `CSV_DECIMAL_COMMA` | warning | a comma inside a number: decimal mark or thousands separator |
| `CSV_UNIT_MIXED` | warning | several units in `unit`; not converted |
| `CSV_ID_DUPLICATE` | warning | repeated `id` values |
| `CSV_CROSSREF_UNMATCHED` / `CSV_CROSSREF_EMPTY` | warning | (`--against`) rows matching no record / the reference has nothing to match against |
| `CSV_DELIMITER` / `CSV_BLANK_ROWS` / `CSV_ROWS_DROPPED` / `CSV_SEQUENCE_LOWERCASE` / `CSV_CROSSREF_MATCHED` | info | detected delimiter / blank rows skipped / rows removed by `--drop-missing` / lower-case sequence letters (left as-is) / match summary |

## Limitations

- **Legacy PDB format only**, fixed columns. mmCIF is not read; entries
  that RCSB serves only as mmCIF (very large assemblies, more than 62
  chains or 99 999 atoms) are reported as such by `fetch`.
- Only `ATOM`, `HETATM`, `TER` and `END` are written. `HEADER`,
  `REMARK`, `SEQRES`, `SSBOND`, `LINK`, `CRYST1`, `HELIX`/`SHEET` and
  every other record type are dropped. The report keeps the provenance
  (entry ID, revision date, SHA-256) instead.
- Alternate locations are resolved **per atom by occupancy**, not by
  whole-residue conformer, and differently from PyMOL's "keep `A`".
- Fv boundaries are not computed. Presets use literature values for
  trastuzumab; for any other antibody pass the ranges yourself. The
  FASTA group warns about likely signal peptides and constant regions
  but never truncates; antibody numbering tools (ANARCI) are not included.
- The FASTA heuristics (length ranges, signal peptide) can misfire on
  unusual antibodies (very long CDR-H3, humanised variants); they are
  warnings and can be switched off with `--no-heuristics`.
- `--against` aligns by a single offset per chain; it is designed for
  "same molecule, maybe missing residues", not for aligning different
  antibodies (it then reports non-correspondence in one line).
- CSV: `.xlsx` is not read (save as *CSV UTF-8* from Excel). Values are
  validated for *format* only -- no unit conversion, no range or
  plausibility checks. Multi-line cells are handled, but a header line is
  required.
- `--strict` turns every warning into a failure, including the honest
  `PDB_RESOLUTION` warning on 1N8Z (2.52 A) and the heuristics above.
- Not a substitute for the ABD validation itself: a clean exit means the
  file has none of the known formatting problems, not that ABD will
  accept the biology.

## Conventions not yet confirmed against the ABD console

The public ABD user guide states the requirements only in general terms
("proper chain identification (H/L/T)", "verify sequence format and
headers"). Three conventions are therefore implemented as the most likely
reading and isolated in one place each, so a correction is a one-file
change:

| Convention | Current default | Where to change |
|---|---|---|
| Chain IDs ABD expects in a PDB framework / target file | Fv written as `L`/`H`; target keeps its original ID (`--rename-target`) | `abd_prep/pdb/presets.py` (`Part.out_chain`) |
| FASTA header for heavy / light / scFv records | `>H ...`, `>L ...`, `>scFv ...` with the original header kept as description | `abd_prep/fasta/roles.py` (`ABD_HEADER_ID`, `fasta_header`) |
| CSV column schema for fine-tuning / results upload | the standard columns above, tool-defined | `abd_prep/csvdata/schema.py` |

If you have confirmed any of these against the console, please open an
issue.

## Documentation

This README is the complete English reference. Three longer documents
are available in Korean under [`docs/`](docs/):

| Document | Audience | Covers |
|---|---|---|
| [`docs/easy-tutorial.ko.md`](docs/easy-tutorial.ko.md) | first-time users, little or no terminal experience | 20-minute path: check Python, download, run `quickstart`, which file goes to which ABD recipe, swapping in your own PDB / FASTA / CSV, reading messages, the six places people get stuck |
| [`docs/user-guide.ko.md`](docs/user-guide.ko.md) | researchers preparing ABD uploads, no programming background assumed | what each file format contains and where it comes from, a 5-minute HER2 / trastuzumab walkthrough, every command with its options and messages explained in plain language, reading `report.json`, CI usage, FAQ, glossary, cheat sheet |
| [`docs/developer-guide.ko.md`](docs/developer-guide.ko.md) | contributors | design principles, architecture and data flow, module-by-module internals, the `report.json` schema, how to add presets / finding codes / actions / format groups, test strategy, release checklist, security notes, decision history |

## Development

Tests use `unittest` and run without the network or any third-party
package (pytest collects them too):

```bash
python3 -m unittest discover -s tests -t .
python3 tests/make_seq_fixtures.py --check     # FASTA/CSV fixtures are byte-identical to their generator
```

The regression tests against the full 1N8Z / 6LBX / 6B9Z entries run
only when those files are present in `~/.cache/abd-prep/`
(`python3 -m abd_prep pdb fetch 1N8Z 6LBX 6B9Z`); otherwise they are
skipped.

**Test fixtures** (`tests/fixtures/`) contain public sequences and
synthetic data only -- never measured results:

| Fixture | Source | Purpose |
|---|---|---|
| `1N8Z_mini.pdb` | RCSB 1N8Z chains A 1-5, B 1-5, C 489-493 + 2 waters; altlocs injected at two atoms (one 0.60/0.40, one 0.50/0.50 tie) | altloc selection and tie rule, chain selection, `TER`/`END` regeneration |
| `6B9Z_cso_mini.pdb` | RCSB 6B9Z chain B 173-177 (B175 is `ACSO` 0.75 / `BCYS` 0.25); B177 rewritten as a synthetic `MSE`; wrapped in two `MODEL` blocks | modified residues, `--std-residues`, residue heterogeneity, first-model-only |
| `kabat_mini.pdb` | RCSB 1N8Z chain B 80-86 renumbered Kabat-style (82, 82A, 82B, 82C) | insertion-code ranges, numbering-gap rules |
| `seed_igg.fasta` | trastuzumab VH (1N8Z B 1-120) and VL (1N8Z A 1-107) with CRLF, lower case, gaps, trailing `*` | normal `igg` path |
| `seed_nanobody.fasta` | **synthetic** VHH: 1N8Z VH with camelid hallmark substitutions V37F/G44E/L45R/W47G (labelled as such in the header) | `nanobody` path |
| `seed_scfv.fasta` | 1N8Z VH-(G4S)3-VL, synthetic linker | `scfv` path, linker detection |
| `seed_warn_fulllength.fasta` | **synthetic** 19-residue leader + 1N8Z VH + synthetic 68-residue C-terminal extension; 1N8Z VL | signal-peptide and length warnings |
| `seed_bad_alphabet.fasta`, `seed_bad_nucleic.fasta`, `seed_bad_count.fasta` | 1N8Z sequences with injected `X`/`B`/`Z` and an internal `*`; a synthetic nucleotide pattern; a third record | error paths |
| `assay_ok.csv` | 1N8Z-derived variant sequences with **made-up** numbers | clean table |
| `assay_messy_utf8bom.csv` | same, as Excel writes it: UTF-8 BOM, `;` delimiter, missing values, `<0.1`, `1,2`, mixed units, duplicate id | messy-table path |
| `assay_messy_cp949.csv` | same with Korean column names in CP949 (EUC-KR) and one row with an unquoted comma | `--encoding` conversion, `CSV_ROW_LENGTH` |

The FASTA/CSV fixtures are generated by `tests/make_seq_fixtures.py`
(so CRLF, BOM and CP949 bytes are never edited by hand) and the PDB
fixtures by `tests/make_fixtures.py` from the cached RCSB files; each
fixture names its source in its header or `REMARK` lines.

Layout:

```
abd_prep/
  cli.py              argument parsing, one function per action, exit codes
  quickstart.py       the simple mode: pdb split + seed FASTA + example CSV + README.txt in one folder
  common/             finding/result model, report.json, SHA-256 / size / encoding checks, residue tables
  pdb/                pdbfile (fixed-column records, TER/END), ranges, clean, validate, fetch, presets, split
  fasta/              fastafile (parser/normaliser), roles (header tokens, ABD headers), validate, crosscheck
  csvdata/            schema (standard columns, missing tokens), csvfile (delimiter, header, rows), validate, crosscheck
tests/                unittest suite + fixtures + fixture generators
```

## Data sources and licenses

- Code: [MIT](LICENSE).
- Structures 1N8Z, 6LBX and 6B9Z are from the RCSB Protein Data Bank and
  are redistributed under
  [CC0 1.0](https://www.rcsb.org/pages/usage-policy). Please cite the
  primary publications of the entries you use.
- Test fixtures contain public trastuzumab sequences derived from 1N8Z,
  clearly labelled synthetic constructs, and invented numbers. No
  experimental data is included.
- This is an independent, community tool. It is not an Amazon Bio
  Discovery component and is not endorsed by AWS; the ABD input
  requirements it checks are taken from the public ABD user guide.
