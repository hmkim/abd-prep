"""``abd-prep quickstart``: every file a first ABD run needs, in one step.

The simple mode. One command fills a folder with a target antigen
structure, a seed antibody structure, the same seed as a FASTA sequence,
an example assay table and a README that says which ABD recipe takes
which file::

    <out>/
      README.txt              which file goes where, where it came from, next steps
      target.pdb              antigen structure           de novo design; directed evolution with antigen
      antibody_fv.pdb         seed antibody structure     directed evolution, property prediction
      antibody_seed.fasta     the same seed as sequence   directed evolution, property prediction
      assay_example.csv       assay-table template        fine-tuning / results (placeholder row)
      quickstart.report.json  sources, SHA-256 of every file, every finding

Nothing here is a new rule. The module composes the existing paths --
``pdb split`` (fetch -> clean -> validate), ``pdb inspect --fasta`` +
``fasta validate --against`` (F14-F18) and ``csv validate --against``
(F19-F22) -- and keeps their findings under the prefixes ``target/``,
``fv/``, ``seed/`` and ``csv/``, so the report reads like the four
commands run by hand.

The example CSV is a *template*: its only row is the seed antibody itself
with a placeholder number, and the ``notes`` cell says so. It exists so
the expected shape can be seen and the CSV path exercised end to end;
nothing in it is a measurement (plan §1: the tool never invents data that
could be mistaken for results).
"""
import csv
import os
from dataclasses import dataclass, field
from datetime import datetime, timezone
from typing import Dict, List, Optional, Sequence, Tuple

from . import __version__
from .common.report import Report
from .common.result import InputError, ValidationResult
from .csvdata import crosscheck as csv_crosscheck
from .csvdata import validate as csv_validate
from .fasta import crosscheck as fasta_crosscheck
from .fasta import validate as fasta_validate
from .fasta.fastafile import Record, write_fasta
from .fasta.roles import ROLE_HEAVY, ROLE_LIGHT
from .pdb import presets as pdb_presets
from .pdb.fetch import Opener, default_opener
from .pdb.split import SplitOptions, SplitOutcome, run_split
from .pdb.validate import inspect_file

OUT_DIR_DEFAULT = "abd-quickstart"
REPORT_NAME = "quickstart.report.json"

#: file key -> file name inside the output directory (fixed, descriptive;
#: provenance lives in README.txt and the report, not in the file name)
FILES: Dict[str, str] = {
    "target": "target.pdb",
    "fv": "antibody_fv.pdb",
    "seed": "antibody_seed.fasta",
    "csv": "assay_example.csv",
    "readme": "README.txt",
}
#: file key -> report.json ``kind``
KINDS: Dict[str, str] = {"target": "pdb", "fv": "pdb", "seed": "fasta", "csv": "csv", "readme": "text"}
#: file key -> which ABD recipe takes it (README.txt and the console summary)
USES: Dict[str, str] = {
    "target": "target antigen -- de novo design (required); directed evolution with antigen (optional)",
    "fv": "seed antibody as a structure -- directed evolution, property prediction",
    "seed": "seed antibody as a sequence -- directed evolution, property prediction (same antibody as antibody_fv.pdb)",
    "csv": "assay / results table template -- fine-tuning, results analysis; replace the placeholder row with your data",
    "readme": "what each file is, where it came from, what to do next",
}

#: placeholder cells of the example CSV row (never a measurement)
EXAMPLE_VALUE = "1.0"
EXAMPLE_UNIT = "nM"
EXAMPLE_ASSAY = "SPR"
EXAMPLE_NOTES = ("EXAMPLE ROW -- the value is a placeholder, not a measurement; "
                 "replace this row with your own data before uploading")


@dataclass
class QuickstartOptions:
    preset: pdb_presets.Preset
    out_dir: str = OUT_DIR_DEFAULT
    force: bool = False
    with_csv: bool = True
    cache: Optional[str] = None
    refresh: bool = False
    path_mode: str = "basename"
    report_path: Optional[str] = None      # None -> <out>/quickstart.report.json; "" -> none


@dataclass
class QuickstartOutcome:
    report: Report
    files: Dict[str, str] = field(default_factory=dict)      # file key -> path written
    report_path: Optional[str] = None
    framework: Optional[str] = None
    summary: Dict[str, str] = field(default_factory=dict)    # file key -> one-line content summary

    @property
    def result(self) -> ValidationResult:
        return self.report.result


# -- planning -------------------------------------------------------------------------

def framework_for(preset: pdb_presets.Preset) -> Optional[str]:
    """``igg`` for heavy + light, ``nanobody`` for heavy only, ``None`` without a heavy chain."""
    heavy, light = preset.part("heavy"), preset.part("light")
    if heavy and light:
        return "igg"
    if heavy:
        return "nanobody"
    return None


def check_preset(preset: pdb_presets.Preset) -> str:
    if preset.part("target") is None:
        raise InputError(f"preset {preset.name!r} has no target part; quickstart needs a target antigen "
                         "and an antibody (use 'pdb split' for partial presets)")
    fw = framework_for(preset)
    if fw is None:
        raise InputError(f"preset {preset.name!r} has no heavy chain; quickstart needs a target antigen "
                         "and an antibody (use 'pdb split' for partial presets)")
    return fw


def planned_files(opts: QuickstartOptions) -> Dict[str, str]:
    """file key -> path, in the order they are written."""
    keys = ["target", "fv", "seed"] + (["csv"] if opts.with_csv else []) + ["readme"]
    return {k: os.path.join(opts.out_dir, FILES[k]) for k in keys}


def check_overwrite(paths: Sequence[str], force: bool) -> None:
    existing = [os.path.basename(p) for p in paths if os.path.exists(p)]
    if existing and not force:
        raise InputError(f"{len(existing)} file(s) already exist in {os.path.dirname(paths[0]) or '.'}: "
                         f"{', '.join(existing)}; pass --force to overwrite them, or choose another --out")


def _prefix(result: ValidationResult, group: str) -> None:
    for f in result.findings:
        f.where = f"{group}/{f.where}" if f.where else group


# -- steps ------------------------------------------------------------------------------

def seed_records(fv_path: str, preset: pdb_presets.Preset) -> List[Record]:
    """Heavy then light, read back from the written Fv file (so FASTA == PDB by construction)."""
    info = inspect_file(fv_path)
    chains: Dict[str, Dict[str, object]] = info["chains"]  # type: ignore[assignment]
    src = os.path.basename(fv_path)
    records: List[Record] = []
    for role, part_role in ((ROLE_HEAVY, "heavy"), (ROLE_LIGHT, "light")):
        part = preset.part(part_role)
        if part is None:
            continue
        ch = part.written_chain
        if ch not in chains:
            raise InputError(f"chain {ch} ({part_role}) is missing from {src}; the preset selection "
                             f"{part.selection!r} produced no residues")
        s = chains[ch]
        records.append(Record(id=src, description=f"chain {ch} {s['first']}-{s['last']} ({s['residues']} res)",
                              sequence=str(s["sequence"]), role=role))
    return records


def write_seed(fv_path: str, seed_path: str, preset: pdb_presets.Preset, framework: str,
               result: ValidationResult) -> Tuple[List[Record], Dict[str, object]]:
    """``pdb inspect --fasta`` + ``fasta validate --framework ... --against fv`` in one step."""
    records = seed_records(fv_path, preset)
    write_fasta(records, seed_path)
    res = ValidationResult()
    fasta_validate.validate_records(records, fasta_validate.FastaOptions(framework=framework), res)
    against = fasta_crosscheck.crosscheck(records, fv_path, {}, res)
    _prefix(res, "seed")
    result.extend(res)
    return records, {"framework": framework, "records": [r.to_dict() for r in records], "against": against}


def example_rows(records: Sequence[Record], preset: pdb_presets.Preset) -> Tuple[List[str], List[List[str]]]:
    """Header and the single placeholder row of the example CSV."""
    target = preset.part("target")
    target_name = f"{target.pdb_id} {target.selection}" if target else ""
    by_role = {r.role: r for r in records}
    if ROLE_LIGHT in by_role:
        header = ["id", "heavy", "light", "value", "unit", "target", "assay", "notes"]
        seq_cells = [by_role[ROLE_HEAVY].sequence, by_role[ROLE_LIGHT].sequence]
    else:
        header = ["id", "sequence", "value", "unit", "target", "assay", "notes"]
        seq_cells = [by_role[ROLE_HEAVY].sequence]
    row = [f"{preset.name}-seed"] + seq_cells + [EXAMPLE_VALUE, EXAMPLE_UNIT, target_name, EXAMPLE_ASSAY, EXAMPLE_NOTES]
    return header, [row]


def write_example_csv(csv_path: str, seed_path: str, records: Sequence[Record], preset: pdb_presets.Preset,
                      result: ValidationResult) -> Dict[str, object]:
    """Write the template and run ``csv validate --against seed`` on it."""
    header, rows = example_rows(records, preset)
    with open(csv_path, "w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh, lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)
    res = ValidationResult()
    hook = csv_crosscheck.CrossRef(csv_crosscheck.load_reference(seed_path))
    details = csv_validate.process_csv(csv_path, "utf-8", csv_validate.CsvOptions(), res, row_hook=hook)
    details["against"] = hook.report(res)
    details["placeholder_rows"] = len(rows)
    _prefix(res, "csv")
    result.extend(res)
    return details


def readme_text(opts: QuickstartOptions, split: SplitOutcome, files: Dict[str, str], framework: str,
                command: str) -> str:
    preset = opts.preset
    now = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
    lines: List[str] = [
        "abd-prep quickstart -- files for a first Amazon Bio Discovery (ABD) run",
        "=" * 72,
        f"generated {now} by abd-prep {__version__}",
        f"preset:  {preset.name}",
        f"         {preset.description}",
        "",
        "SOURCES",
        "-------",
    ]
    for src in split.report.details.get("sources", []):
        entry = src.get("entry") or {}
        bits = [f"RCSB {src.get('pdb_id')}"]
        if entry.get("title"):
            bits.append(str(entry["title"]))
        if entry.get("resolution_angstrom") is not None:
            bits.append(f"{entry['resolution_angstrom']} A")
        if entry.get("revision_date"):
            bits.append(f"revised {entry['revision_date']}")
        lines.append("  " + " -- ".join(bits))
    lines += ["  Structures come from the RCSB Protein Data Bank (CC0). Cite the primary publication.", ""]

    lines += ["FILES", "-----"]
    width = max(len(FILES[k]) for k in files)
    for key in files:
        lines.append(f"  {FILES[key]:{width}s}  {USES[key]}")
    lines += [
        "",
        "WHICH FILE FOR WHICH ABD RECIPE",
        "-------------------------------",
        f"  Antibody generation, de novo:              {FILES['target']} as the target antigen "
        "(choose the framework in the console)",
        f"  Antibody generation, directed evolution:   {FILES['seed']} (or {FILES['fv']}) as the seed antibody; "
        f"add {FILES['target']} for antigen-aware optimisation",
        f"  Property prediction:                       {FILES['seed']} (or {FILES['fv']}) as the seed antibody",
    ]
    if "csv" in files:
        lines.append(f"  Fine-tuning / results analysis:            {FILES['csv']} shows the expected table shape; "
                     "replace its placeholder row with your own measurements")
    target_part = preset.part("target")
    target_chain = target_part.written_chain if target_part else "<ID>"
    lines += [
        "",
        "BEFORE YOU UPLOAD",
        "-----------------",
        f"  * The antibody framework is '{framework}' "
        + ("(one heavy + one light chain, written as chains H and L / records >H and >L)."
           if framework == "igg" else "(one heavy-chain domain, written as chain H / record >H)."),
        f"  * The target keeps its chain ID from the source entry ({target_chain}); if the console asks for a "
        f"specific ID such as T, rename it: abd-prep pdb clean {FILES['target']} target_T.pdb --rename {target_chain}=T",
    ]
    if "csv" in files:
        lines.append(f"  * {FILES['csv']} contains no measurement. Its one row is the seed itself with a placeholder "
                     "value so you can see the columns; the 'notes' cell says so.")
    lines += [
        f"  * {REPORT_NAME} records every step, the RCSB revision used and the SHA-256 of each file. "
        "Keep it next to the files you upload.",
        "",
        "REPRODUCE",
        "---------",
        f"  {command or 'abd-prep quickstart --preset ' + preset.name}",
        "",
        "Three ABD conventions are not confirmed against the console yet (chain IDs, FASTA headers, CSV "
        "columns); the project README lists them and the one-line fixes. abd-prep is an independent tool, "
        "not an ABD component.",
        "",
    ]
    return "\n".join(lines)


# -- driver -----------------------------------------------------------------------------

def run_quickstart(opts: QuickstartOptions, opener: Opener = default_opener, command: str = "") -> QuickstartOutcome:
    framework = check_preset(opts.preset)
    files = planned_files(opts)
    check_overwrite(list(files.values()), opts.force)
    os.makedirs(opts.out_dir, exist_ok=True)

    report = Report(command=command, path_mode=opts.path_mode)
    outcome = QuickstartOutcome(report=report, framework=framework)

    # 1. structures: fetch -> clean -> validate (pdb split)
    sopts = SplitOptions.from_preset(opts.preset, opts.out_dir, cache=opts.cache, refresh=opts.refresh,
                                     path_mode=opts.path_mode, report_path="",
                                     filenames={"target": FILES["target"], "fv": FILES["fv"]})
    split = run_split(sopts, opener=opener, command=command)
    report.inputs.extend(split.report.inputs)
    report.result.extend(split.result)
    split_details = dict(split.report.details)
    split_details.pop("preset", None)
    report.details = {"preset": opts.preset.to_dict(), "framework": framework, "split": split_details}
    for group in ("target", "fv"):
        outcome.files[group] = split.written[group]
        report.add_output(KINDS[group], split.written[group])

    # 2. seed sequence: inspect --fasta + fasta validate --against
    records, seed_details = write_seed(files["fv"], files["seed"], opts.preset, framework, report.result)
    report.details["seed"] = seed_details
    outcome.files["seed"] = files["seed"]
    report.add_output(KINDS["seed"], files["seed"])

    # 3. example table: write + csv validate --against seed
    if opts.with_csv:
        report.details["csv"] = write_example_csv(files["csv"], files["seed"], records, opts.preset, report.result)
        outcome.files["csv"] = files["csv"]
        report.add_output(KINDS["csv"], files["csv"])

    # 4. README.txt
    with open(files["readme"], "w", encoding="utf-8", newline="\n") as fh:
        fh.write(readme_text(opts, split, files, framework, command))
    outcome.files["readme"] = files["readme"]
    report.add_output(KINDS["readme"], files["readme"])

    report.details["files"] = {k: {"path": FILES[k], "use": USES[k]} for k in outcome.files}
    outcome.summary = _summaries(report.details, outcome.files)

    if opts.report_path == "":
        outcome.report_path = None
    else:
        outcome.report_path = opts.report_path or os.path.join(opts.out_dir, REPORT_NAME)
        report.write(outcome.report_path)
    return outcome


def _summaries(details: Dict[str, object], files: Dict[str, str]) -> Dict[str, str]:
    """One line of content per file for the console summary."""
    out: Dict[str, str] = {}
    parts = details.get("split", {}).get("parts", {})  # type: ignore[union-attr]
    for group, label in (("target", "antigen structure"), ("fv", "seed antibody structure")):
        chains = parts.get(group, {}).get("clean", {}).get("chains_out", {})
        desc = ", ".join(f"chain {ch} {d['residues']} res" for ch, d in chains.items())
        out[group] = f"{label} ({desc})" if desc else label
    recs = details.get("seed", {}).get("records", [])  # type: ignore[union-attr]
    out["seed"] = "seed antibody sequence (" + ", ".join(f">{r['role']} {r['length']} res" for r in recs) + ")"
    if "csv" in files:
        n = details.get("csv", {}).get("placeholder_rows", 0)  # type: ignore[union-attr]
        out["csv"] = f"assay-table template ({n} placeholder row, not a measurement)"
    out["readme"] = "what each file is, where it came from, what to do next"
    return out


def format_files(outcome: QuickstartOutcome, out_dir: str) -> str:
    """Console block printed after the findings table."""
    lines = [f"Files in {out_dir}:"]
    width = max(len(FILES[k]) for k in outcome.files)
    for key in outcome.files:
        lines.append(f"  {FILES[key]:{width}s}  {outcome.summary.get(key, '')}")
        if key != "readme":
            lines.append(f"  {'':{width}s}  -> {USES[key]}")
    lines.append(f"Next: read {FILES['readme']}, then create an experiment in the ABD console and upload the "
                 "file(s) its recipe asks for.")
    return "\n".join(lines)
