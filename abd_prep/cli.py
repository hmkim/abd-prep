"""Command line: ``abd-prep <group> <action> ...`` (plan §3.1).

Groups are input formats (``pdb``, ``fasta``, ``csv``); actions are
``clean`` / ``validate`` / ``inspect`` plus the PDB-only ``fetch`` and
``split`` (M2). ``abd-prep quickstart`` is the one cross-format command:
the simple mode that prepares every file a first ABD run needs in one
step (see :mod:`abd_prep.quickstart`). Exit codes: 0 pass (warnings
allowed), 1 validation error (or warnings with ``--strict``), 2
input/usage error.

Common options accepted by every ``clean``/``validate`` (and ``quickstart``):

``--json``                 print the report as JSON instead of text
``--strict``               warnings also fail (exit 1)
``--encoding CODEC``       input text encoding (FASTA/CSV; ignored for PDB)
``--report PATH``          where to write report.json
``--no-report``            do not write report.json
``--report-paths MODE``    ``basename`` (default) or ``full`` paths in the report
"""
import argparse
import os
import sys
from typing import Callable, Dict, List, Optional, Sequence, Tuple

from . import __version__, quickstart
from .common import textfile
from .common.report import PATH_MODES, Report, render_command
from .common.result import EXIT_INPUT, InputError, ValidationResult
from .csvdata import crosscheck as csv_crosscheck
from .csvdata import schema as csv_schema
from .csvdata import validate as csv_validate
from .csvdata.csvfile import parse_delimiter
from .fasta import crosscheck as fasta_crosscheck
from .fasta import validate as fasta_validate
from .fasta.fastafile import read_fasta, write_fasta
from .fasta.roles import FRAMEWORKS
from .pdb import clean as pdb_clean
from .pdb import fetch as pdb_fetch
from .pdb import presets as pdb_presets
from .pdb import split as pdb_split
from .pdb import validate as pdb_validate
from .pdb.ranges import ChainSelection, RangeSyntaxError, parse_rename

PDB_EXTENSIONS = (".pdb", ".ent", ".pdb1")


# -- argument parsing -----------------------------------------------------------

def _common_options(parser: argparse.ArgumentParser, with_report: bool = True, with_encoding: bool = True) -> None:
    g = parser.add_argument_group("common options")
    g.add_argument("--json", action="store_true", help="print the report as JSON instead of a text table")
    g.add_argument("--strict", action="store_true", help="treat warnings as failures (exit code 1)")
    if with_encoding:
        g.add_argument("--encoding", default="utf-8",
                       help="input text encoding for FASTA/CSV (default utf-8; PDB is ASCII and ignores this)")
    if with_report:
        g.add_argument("--report", metavar="PATH", help="write report.json to PATH")
        g.add_argument("--no-report", action="store_true", help="do not write report.json")
        g.add_argument("--report-paths", choices=PATH_MODES, default="basename",
                       help="record file paths in the report as basename (default) or full")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="abd-prep",
        description="Prepare and validate Amazon Bio Discovery inputs (PDB structures, FASTA seed "
                    "antibodies, assay CSVs). Standard library only. New here? Run: abd-prep quickstart")
    parser.add_argument("--version", action="version", version=f"abd-prep {__version__}")
    groups = parser.add_subparsers(dest="group", metavar="GROUP")
    groups.required = True
    preset_names = sorted(pdb_presets.PRESETS)

    # ======================================================================
    # quickstart (simple mode: every file for a first ABD run, one command)
    # ======================================================================
    p = groups.add_parser("quickstart", help="simple mode: prepare every file a first ABD run needs in one step "
                                             "(target PDB, antibody Fv PDB, seed FASTA, example CSV, README)")
    p.add_argument("--preset", choices=preset_names, metavar="NAME", default=pdb_presets.DEFAULT_PRESET,
                   help=f"named target + antibody selection (default {pdb_presets.DEFAULT_PRESET}); one of: "
                        + ", ".join(preset_names))
    p.add_argument("--out", metavar="DIR", default=quickstart.OUT_DIR_DEFAULT,
                   help=f"output directory (default ./{quickstart.OUT_DIR_DEFAULT})")
    p.add_argument("--force", action="store_true", help="overwrite files that already exist in --out")
    p.add_argument("--no-csv", action="store_true", help="do not write the example assay table")
    p.add_argument("--cache", metavar="DIR", help="RCSB download cache directory")
    p.add_argument("--refresh", action="store_true", help="re-download cached entries")
    _common_options(p, with_encoding=False)
    p.set_defaults(func=cmd_quickstart)

    pdb = groups.add_parser("pdb", help="protein structures (target antigen, antibody Fv)")
    actions = pdb.add_subparsers(dest="action", metavar="ACTION")
    actions.required = True

    # pdb fetch --------------------------------------------------------------
    p = actions.add_parser("fetch", help="download PDB-format file(s) from RCSB into the local cache")
    p.add_argument("ids", nargs="+", metavar="PDB_ID", help="4-character RCSB entry ID(s), e.g. 1N8Z")
    p.add_argument("--out", metavar="DIR", help="also copy the file(s) into DIR")
    p.add_argument("--cache", metavar="DIR",
                   help=f"cache directory (default {pdb_fetch.DEFAULT_CACHE_DIR}, or ${pdb_fetch.CACHE_ENV})")
    p.add_argument("--refresh", action="store_true", help="re-download even if the file is cached")
    p.add_argument("--no-meta", action="store_true", help="skip the entry-metadata request")
    p.add_argument("--json", action="store_true", help="print the result as JSON")
    p.set_defaults(func=cmd_pdb_fetch)

    # pdb split --------------------------------------------------------------
    p = actions.add_parser("split", help="cut target and antibody Fv files from RCSB entries or a local file, "
                                         "clean them, regenerate TER/END and validate (fetch + clean + validate)")
    src = p.add_mutually_exclusive_group()
    src.add_argument("--preset", choices=preset_names, metavar="NAME",
                     help="named selection; one of: " + ", ".join(preset_names))
    src.add_argument("--pdb", metavar="PDB_ID", help="RCSB entry to cut from (used with --target/--heavy/--light)")
    src.add_argument("--input", metavar="FILE", help="local PDB file to cut from instead of an RCSB entry")
    p.add_argument("--target", metavar="SPEC", help="target chain and range, e.g. C:489-607")
    p.add_argument("--heavy", metavar="SPEC", help="heavy-chain range, e.g. B:1-120")
    p.add_argument("--light", metavar="SPEC", help="light-chain range, e.g. A:1-107")
    p.add_argument("--rename-target", metavar="ID", help="chain ID to write for the target (default: keep)")
    p.add_argument("--rename-heavy", metavar="ID", help="chain ID to write for the heavy chain (default: keep)")
    p.add_argument("--rename-light", metavar="ID", help="chain ID to write for the light chain (default: keep)")
    p.add_argument("--name", metavar="NAME",
                   help="output file stem: NAME_target.pdb / NAME_fv.pdb / NAME.report.json "
                        "(default: preset name, or the entry ID / input file name)")
    p.add_argument("--out", metavar="DIR", default=".", help="output directory (default: current directory)")
    p.add_argument("--cache", metavar="DIR", help="RCSB download cache directory")
    p.add_argument("--refresh", action="store_true", help="re-download cached entries")
    p.add_argument("--keep-het", action="store_true", help="keep HETATM records (default: remove)")
    p.add_argument("--std-residues", action="store_true",
                   help="convert modified backbone residues (CSO, MSE, ...) to their standard parent")
    p.add_argument("--strip-h", action="store_true", help="remove hydrogen (and deuterium) atoms")
    p.add_argument("--renumber", action="store_true", help="rewrite atom and TER serial numbers 1..N")
    _common_options(p)
    p.set_defaults(func=cmd_pdb_split)

    # pdb clean --------------------------------------------------------------
    p = actions.add_parser("clean", help="remove altlocs/HETATM, select chains and residues, rename chains, "
                                         "regenerate TER/END (prototype clean_pdb.py compatible)")
    p.add_argument("input", help="input PDB file")
    p.add_argument("output", help="output PDB file")
    p.add_argument("--keep", action="append", default=[], metavar="SPEC",
                   help="chain or chain:ranges to keep, e.g. B or A:1-107 or C:489-607,610-612 or H:1-113A "
                        "(repeatable; default keeps every chain)")
    p.add_argument("--rename", default="", metavar="OLD=NEW,...", help="rename chains, e.g. A=L,B=H")
    p.add_argument("--keep-het", action="store_true", help="keep HETATM records (default: remove)")
    p.add_argument("--std-residues", action="store_true",
                   help="convert modified backbone residues (CSO, MSE, SEP, TPO, PTR, PCA ...) to their "
                        "standard parent instead of removing them with the HETATM records")
    p.add_argument("--strip-h", action="store_true", help="remove hydrogen (and deuterium) atoms")
    p.add_argument("--renumber", action="store_true", help="rewrite atom and TER serial numbers 1..N")
    _common_options(p)
    p.set_defaults(func=cmd_pdb_clean)

    # pdb validate -----------------------------------------------------------
    p = actions.add_parser("validate", help="check ABD readiness: altlocs, HETATM, chains, TER/END, gaps")
    p.add_argument("files", nargs="+", help="PDB file(s) to check")
    p.add_argument("--allow-het", action="store_true", help="HETATM records are expected (info, not warning)")
    p.add_argument("--expect-chains", type=int, metavar="N", help="fail unless the file has exactly N chains")
    p.add_argument("--preset", choices=preset_names, metavar="NAME",
                   help="apply a preset's expectations (residue counts, known gaps); each file is matched to "
                        "the preset's target or Fv part by its chain IDs")
    _common_options(p)
    p.set_defaults(func=cmd_pdb_validate)

    # pdb inspect ------------------------------------------------------------
    p = actions.add_parser("inspect", help="summarize chains, ranges and sequences "
                                           "(prototype inspect_pdb.py compatible)")
    p.add_argument("files", nargs="+", help="PDB file(s)")
    p.add_argument("--fasta", action="store_true", help="print the chain sequences as FASTA instead")
    p.add_argument("--json", action="store_true", help="print the summary as JSON")
    p.set_defaults(func=cmd_pdb_inspect)

    # ======================================================================
    # fasta group (seed antibody sequences, M2.5)
    # ======================================================================
    fasta = groups.add_parser("fasta", help="seed antibody sequences (directed evolution, property prediction)")
    factions = fasta.add_subparsers(dest="action", metavar="ACTION")
    factions.required = True

    def _fasta_role_options(parser: argparse.ArgumentParser, framework_required: bool) -> None:
        parser.add_argument("--framework", choices=FRAMEWORKS, required=framework_required,
                            help="antibody format: igg (one heavy + one light record), nanobody (one VHH record), "
                                 "scfv (one record, linker kept)" + ("" if framework_required else
                                 "; when omitted the record-count/role rules are skipped"))
        parser.add_argument("--heavy-id", metavar="ID", help="record ID (first header token) of the heavy chain")
        parser.add_argument("--light-id", metavar="ID", help="record ID (first header token) of the light chain")
        parser.add_argument("--replace-ambiguous", metavar="X=G,...",
                            help="replace ambiguous letters explicitly, e.g. X=G,B=N (never done automatically)")
        parser.add_argument("--no-heuristics", action="store_true",
                            help="skip the length / signal-peptide / duplicate-sequence warnings")

    p = factions.add_parser("clean", help="normalize a seed FASTA (line ends, case, gaps, headers) for upload")
    p.add_argument("input", help="input FASTA file")
    p.add_argument("output", help="output FASTA file (UTF-8, LF, ABD headers)")
    _fasta_role_options(p, framework_required=True)
    p.add_argument("--wrap", type=int, default=0, metavar="N", help="wrap sequence lines at N characters (default: one line)")
    _common_options(p)
    p.set_defaults(func=cmd_fasta_clean)

    p = factions.add_parser("validate", help="check a seed FASTA: alphabet, chain roles, framework, heuristics")
    p.add_argument("files", nargs="+", help="FASTA file(s) to check")
    _fasta_role_options(p, framework_required=False)
    p.add_argument("--against", metavar="PDB", help="compare the records with the chains of this PDB file (F18)")
    p.add_argument("--chain-map", metavar="H=B,L=A",
                   help="which PDB chain is heavy/light when the PDB chains are not labelled H/L")
    _common_options(p)
    p.set_defaults(func=cmd_fasta_validate)

    p = factions.add_parser("inspect", help="summarize records: id, inferred role, length, alphabet")
    p.add_argument("files", nargs="+", help="FASTA file(s)")
    p.add_argument("--framework", choices=FRAMEWORKS, help="framework used for role assignment")
    p.add_argument("--heavy-id", metavar="ID", help="record ID of the heavy chain")
    p.add_argument("--light-id", metavar="ID", help="record ID of the light chain")
    p.add_argument("--encoding", default="utf-8", help="input text encoding (default utf-8)")
    p.add_argument("--json", action="store_true", help="print the summary as JSON")
    p.set_defaults(func=cmd_fasta_inspect)

    # ======================================================================
    # csv group (assay / results tables, M2.5)
    # ======================================================================
    csvg = groups.add_parser("csv", help="assay / results tables (fine-tuning data, wet-lab measurements)")
    cactions = csvg.add_subparsers(dest="action", metavar="ACTION")
    cactions.required = True

    def _csv_options(parser: argparse.ArgumentParser) -> None:
        parser.add_argument("--map", action="append", default=[], metavar="STD=HEADER,...",
                            help="project source headers onto standard columns, e.g. "
                                 "id=Clone,sequence=AA_seq,value=KD,unit=Unit (repeatable). Standard columns: "
                                 + ", ".join(csv_schema.STANDARD_COLUMNS))
        parser.add_argument("--delimiter", metavar="CHAR",
                            help="input delimiter (',' ';' '|' or tab); detected when omitted")

    p = cactions.add_parser("clean", help="normalize a CSV: UTF-8 without BOM, ',' delimiter, standard column names")
    p.add_argument("input", help="input CSV file")
    p.add_argument("output", help="output CSV file")
    _csv_options(p)
    p.add_argument("--drop-missing", action="store_true", help="drop rows with a missing value in a standard column")
    p.add_argument("--drop-unmapped", action="store_true", help="drop columns that are not mapped to a standard column")
    _common_options(p)
    p.set_defaults(func=cmd_csv_clean)

    p = cactions.add_parser("validate", help="check a CSV: encoding, delimiter, header, standard columns, values")
    p.add_argument("files", nargs="+", help="CSV file(s) to check")
    _csv_options(p)
    p.add_argument("--against", metavar="FILE", help="cross-reference ids/sequences with a seed FASTA or Fv PDB (F22)")
    _common_options(p)
    p.set_defaults(func=cmd_csv_validate)

    p = cactions.add_parser("inspect", help="summarize columns, fill rates and sample values")
    p.add_argument("files", nargs="+", help="CSV file(s)")
    _csv_options(p)
    p.add_argument("--encoding", default="utf-8", help="input text encoding (default utf-8)")
    p.add_argument("--json", action="store_true", help="print the summary as JSON")
    p.set_defaults(func=cmd_csv_inspect)

    return parser


# -- helpers -----------------------------------------------------------------------

def _require_file(path: str) -> None:
    if not os.path.isfile(path):
        raise InputError(f"file not found: {path}")


def _extension_warning(path: str, result: ValidationResult) -> None:
    if not path.lower().endswith(PDB_EXTENSIONS):
        result.warning("PDB_EXTENSION", f"{os.path.basename(path)!r} does not have a .pdb extension; "
                       "the file is still parsed as PDB", "file")


def _default_report_path(output: str) -> str:
    root, _ = os.path.splitext(output)
    return root + ".report.json"


def _emit(report: Report, args: argparse.Namespace, title: str) -> None:
    if getattr(args, "json", False):
        print(report.to_json())
    else:
        print(report.result.format_table(title))
    if getattr(args, "no_report", False):
        return
    path = getattr(args, "report", None)
    if path:
        report.write(path)
        if not args.json:
            print(f"  report: {path}")


# -- commands -----------------------------------------------------------------------

def cmd_quickstart(args: argparse.Namespace, argv: Sequence[str]) -> int:
    preset = pdb_presets.get_preset(args.preset)
    opts = quickstart.QuickstartOptions(
        preset=preset, out_dir=args.out, force=args.force, with_csv=not args.no_csv,
        cache=args.cache, refresh=args.refresh, path_mode=args.report_paths,
        report_path="" if args.no_report else args.report)
    outcome = quickstart.run_quickstart(opts, command=render_command(argv, args.report_paths))
    if args.json:
        print(outcome.report.to_json())
    else:
        print(outcome.result.format_table(f"[quickstart] {preset.name}  ->  {args.out}"))
        if outcome.report_path:
            print(f"  report: {outcome.report_path}")
        print()
        print(quickstart.format_files(outcome, args.out))
    return outcome.result.exit_code(args.strict)


def cmd_pdb_clean(args: argparse.Namespace, argv: Sequence[str]) -> int:
    _require_file(args.input)
    try:
        opts = pdb_clean.CleanOptions(
            selection=ChainSelection.parse(args.keep),
            rename=parse_rename(args.rename),
            keep_het=args.keep_het,
            std_residues=args.std_residues,
            strip_h=args.strip_h,
            renumber=args.renumber,
        )
    except RangeSyntaxError as exc:
        raise InputError(str(exc))
    if os.path.abspath(args.input) == os.path.abspath(args.output):
        raise InputError("input and output are the same file")

    report = Report(command=render_command(argv, args.report_paths), path_mode=args.report_paths)
    report.add_input("pdb", args.input)
    _extension_warning(args.input, report.result)
    outcome = pdb_clean.clean_file(args.input, args.output, opts)
    report.result.extend(outcome.result)
    report.details = outcome.details
    report.add_output("pdb", args.output)
    if args.report is None and not args.no_report:
        args.report = _default_report_path(args.output)
    chains = ", ".join(f"{ch}({d['residues']} res)" for ch, d in outcome.details["chains_out"].items())
    _emit(report, args, f"[clean] {args.input} -> {args.output}  chains: {chains or 'none'}")
    return report.result.exit_code(args.strict)


def _preset_options_for(preset: pdb_presets.Preset, path: str) -> Tuple[Optional[str], pdb_validate.ValidateOptions]:
    """Match a file to the preset's target or Fv part by chain IDs; return (part, options)."""
    info = pdb_validate.inspect_file(path)
    chains = set(info["chains"])
    for group, parts in pdb_presets.output_groups(preset.parts).items():
        out_chains = {p.written_chain for p in parts}
        if chains == out_chains:
            return group, pdb_validate.ValidateOptions(
                expect_chains=len(parts),
                expect_residues={ch: (n, n) for ch, n in preset.expect_residues.items() if ch in out_chains},
                expected_gaps={g for g in preset.expected_gaps if g[0] in out_chains},
            )
    return None, pdb_validate.ValidateOptions()


def cmd_pdb_validate(args: argparse.Namespace, argv: Sequence[str]) -> int:
    for f in args.files:
        _require_file(f)
    preset = pdb_presets.get_preset(args.preset) if args.preset else None
    report = Report(command=render_command(argv, args.report_paths), path_mode=args.report_paths)
    per_file: Dict[str, object] = {}
    worst = 0
    for path in args.files:
        report.add_input("pdb", path)
        if preset is not None:
            part, vopts = _preset_options_for(preset, path)
            vopts.allow_het = args.allow_het
            if args.expect_chains is not None:
                vopts.expect_chains = args.expect_chains
        else:
            part = None
            vopts = pdb_validate.ValidateOptions(allow_het=args.allow_het, expect_chains=args.expect_chains)
        res, details = pdb_validate.validate_file(path, vopts)
        if preset is not None:
            if part is None:
                res.info("PDB_PRESET_NO_MATCH", f"chain IDs do not match any part of preset {preset.name!r}; "
                         "generic rules applied", "file")
            else:
                details["preset_part"] = part
        _extension_warning(path, res)
        per_file[os.path.basename(path) if args.report_paths == "basename" else path] = details
        # prefix the location with the file name when several files are checked
        if len(args.files) > 1:
            for f in res.findings:
                f.where = f"{os.path.basename(path)}:{f.where}" if f.where else os.path.basename(path)
        report.result.extend(res)
        worst = max(worst, res.exit_code(args.strict))
    report.details = {"files": per_file, "preset": preset.name if preset else None}
    _emit(report, args, f"[validate] {' '.join(args.files)}")
    return worst


def cmd_pdb_fetch(args: argparse.Namespace, argv: Sequence[str]) -> int:
    import json
    results = []
    for raw_id in args.ids:
        fr = pdb_fetch.fetch(raw_id, cache=args.cache, refresh=args.refresh, with_meta=not args.no_meta)
        dest = pdb_fetch.copy_to(fr, args.out) if args.out else fr.path
        entry = fr.to_dict("full")
        entry["path"] = dest
        results.append(entry)
        if not args.json:
            origin = "cache" if fr.from_cache else "RCSB"
            print(f"{fr.pdb_id}: {dest}  ({fr.bytes} bytes, {origin})")
            if fr.meta:
                m = fr.meta
                res = f"{m.resolution:.2f} A" if m.resolution is not None else "n/a"
                print(f"  {m.title}")
                print(f"  method {m.experimental_method or 'n/a'}, resolution {res}, "
                      f"released {m.initial_release_date or 'n/a'}, revised {m.revision_date or 'n/a'}")
            elif fr.meta_error:
                print(f"  entry metadata unavailable: {fr.meta_error}")
    if args.json:
        print(json.dumps(results if len(results) > 1 else results[0], indent=2))
    return 0


def _split_options(args: argparse.Namespace) -> pdb_split.SplitOptions:
    common = dict(cache=args.cache, refresh=args.refresh, keep_het=args.keep_het,
                  std_residues=args.std_residues, strip_h=args.strip_h, renumber=args.renumber,
                  path_mode=args.report_paths,
                  report_path="" if args.no_report else args.report)
    if args.preset:
        if any([args.target, args.heavy, args.light, args.rename_target, args.rename_heavy, args.rename_light]):
            raise InputError("--preset cannot be combined with --target/--heavy/--light/--rename-*; "
                             "spell the selection out instead of using the preset")
        preset = pdb_presets.get_preset(args.preset)
        opts = pdb_split.SplitOptions.from_preset(preset, args.out, **common)
        if args.name:
            opts.name = args.name
        return opts

    if args.input:
        _require_file(args.input)
        source = pdb_split.LOCAL_SOURCE
    elif args.pdb:
        source = pdb_fetch.normalize_pdb_id(args.pdb)
    else:
        raise InputError("give --preset, or --pdb PDB_ID / --input FILE together with --target/--heavy/--light")

    parts: List[pdb_presets.Part] = []
    for role, spec, rename in (("target", args.target, args.rename_target),
                               ("heavy", args.heavy, args.rename_heavy),
                               ("light", args.light, args.rename_light)):
        if spec:
            parts.append(pdb_presets.Part(role, source, spec, rename))
        elif rename:
            raise InputError(f"--rename-{role} given without --{role}")
    parts_t = tuple(parts)
    pdb_presets.validate_parts(parts_t)
    name = args.name or pdb_split.derive_name(parts_t, args.input)
    return pdb_split.SplitOptions(name=name, parts=parts_t, out_dir=args.out, input_path=args.input, **common)


def cmd_pdb_split(args: argparse.Namespace, argv: Sequence[str]) -> int:
    opts = _split_options(args)
    outcome = pdb_split.run_split(opts, command=render_command(argv, args.report_paths))
    report = outcome.report
    if args.json:
        print(report.to_json())
    else:
        files = ", ".join(f"{g}: {p}" for g, p in outcome.written.items())
        print(report.result.format_table(f"[split] {opts.name}  ->  {files}"))
        if outcome.report_path:
            print(f"  report: {outcome.report_path}")
    return report.result.exit_code(args.strict)


def cmd_pdb_inspect(args: argparse.Namespace, argv: Sequence[str]) -> int:
    import json
    infos = []
    for path in args.files:
        _require_file(path)
        infos.append(pdb_validate.inspect_file(path))
    if args.json:
        print(json.dumps(infos if len(infos) > 1 else infos[0], indent=2))
    elif args.fasta:
        sys.stdout.write("".join(pdb_validate.to_fasta(i) for i in infos))
    else:
        print("\n".join(pdb_validate.format_inspect(i) for i in infos))
    return 0


# -- fasta commands -------------------------------------------------------------------

def _fasta_options(args: argparse.Namespace) -> fasta_validate.FastaOptions:
    return fasta_validate.FastaOptions(
        framework=getattr(args, "framework", None),
        heavy_id=getattr(args, "heavy_id", None),
        light_id=getattr(args, "light_id", None),
        replace_ambiguous=fasta_validate.parse_replacements(getattr(args, "replace_ambiguous", None)),
        heuristics=not getattr(args, "no_heuristics", False),
    )


def _load_fasta(path: str, encoding: str, action: str, opts: fasta_validate.FastaOptions, result: ValidationResult):
    """precheck + parse + F14-F17; returns the FastaFile or None when undecodable."""
    textfile.check_extension(path, textfile.FASTA_EXTENSIONS, "FASTA", result)
    codec = textfile.precheck(path, encoding, result, action)
    if codec is None:
        return None
    fasta = read_fasta(path, result, encoding=codec)
    fasta_validate.validate_records(fasta.records, opts, result)
    return fasta


def cmd_fasta_clean(args: argparse.Namespace, argv: Sequence[str]) -> int:
    _require_file(args.input)
    if os.path.abspath(args.input) == os.path.abspath(args.output):
        raise InputError("input and output are the same file")
    opts = _fasta_options(args)
    report = Report(command=render_command(argv, args.report_paths), path_mode=args.report_paths)
    report.add_input("fasta", args.input)
    fasta = _load_fasta(args.input, args.encoding, "clean", opts, report.result)
    written = False
    if fasta is not None:
        write_fasta(fasta.records, args.output, wrap=args.wrap)
        report.add_output("fasta", args.output)
        written = True
        report.details = {"framework": opts.framework, "records": [r.to_dict() for r in fasta.records]}
    if args.report is None and not args.no_report:
        args.report = _default_report_path(args.output)
    summary = ", ".join(f"{r.role or '?'}({r.length} res)" for r in fasta.records) if fasta else "not written"
    _emit(report, args, f"[fasta clean] {args.input} -> {args.output if written else '(no output)'}  records: {summary}")
    return report.result.exit_code(args.strict)


def cmd_fasta_validate(args: argparse.Namespace, argv: Sequence[str]) -> int:
    for f in args.files:
        _require_file(f)
    if args.against and len(args.files) != 1:
        raise InputError("--against compares one FASTA file at a time; give a single input file")
    if args.chain_map and not args.against:
        raise InputError("--chain-map only makes sense together with --against")
    if args.against:
        _require_file(args.against)
    chain_map = fasta_crosscheck.parse_chain_map(args.chain_map)
    opts = _fasta_options(args)
    report = Report(command=render_command(argv, args.report_paths), path_mode=args.report_paths)
    per_file: Dict[str, object] = {}
    worst = 0
    for path in args.files:
        report.add_input("fasta", path)
        res = ValidationResult()
        fasta = _load_fasta(path, args.encoding, "validate", opts, res)
        details: Dict[str, object] = {"framework": opts.framework}
        if fasta is not None:
            details["records"] = [r.to_dict() for r in fasta.records]
            if args.against:
                report.add_input("pdb", args.against)
                details["against"] = fasta_crosscheck.crosscheck(fasta.records, args.against, chain_map, res)
        per_file[os.path.basename(path) if args.report_paths == "basename" else path] = details
        if len(args.files) > 1:
            for f in res.findings:
                f.where = f"{os.path.basename(path)}:{f.where}" if f.where else os.path.basename(path)
        report.result.extend(res)
        worst = max(worst, res.exit_code(args.strict))
    report.details = {"files": per_file}
    _emit(report, args, f"[fasta validate] {' '.join(args.files)}")
    return worst


def cmd_fasta_inspect(args: argparse.Namespace, argv: Sequence[str]) -> int:
    import json
    infos = []
    for path in args.files:
        _require_file(path)
        res = ValidationResult()
        textfile.check_extension(path, textfile.FASTA_EXTENSIONS, "FASTA", res)
        codec = textfile.precheck(path, args.encoding, res, "validate")
        if codec is None:
            raise InputError("; ".join(f.message for f in res.errors))
        fasta = read_fasta(path, res, encoding=codec)
        fasta_validate.validate_records(fasta.records, _fasta_options(args), res)
        info = fasta_validate.inspect_records(fasta.records, path)
        info["findings"] = res.to_list()
        infos.append(info)
    if args.json:
        print(json.dumps(infos if len(infos) > 1 else infos[0], indent=2, ensure_ascii=False))
    else:
        for info in infos:
            print(fasta_validate.format_inspect(info))
    return 0


# -- csv commands ---------------------------------------------------------------------

def _csv_options(args: argparse.Namespace) -> csv_validate.CsvOptions:
    return csv_validate.CsvOptions(
        mapping=csv_schema.parse_map(getattr(args, "map", None)),
        delimiter=parse_delimiter(getattr(args, "delimiter", None)),
        drop_missing=getattr(args, "drop_missing", False),
        drop_unmapped=getattr(args, "drop_unmapped", False),
    )


def cmd_csv_clean(args: argparse.Namespace, argv: Sequence[str]) -> int:
    _require_file(args.input)
    if os.path.abspath(args.input) == os.path.abspath(args.output):
        raise InputError("input and output are the same file")
    opts = _csv_options(args)
    report = Report(command=render_command(argv, args.report_paths), path_mode=args.report_paths)
    report.add_input("csv", args.input)
    textfile.check_extension(args.input, textfile.CSV_EXTENSIONS, "CSV", report.result)
    codec = textfile.precheck(args.input, args.encoding, report.result, "clean")
    written = False
    if codec is not None:
        report.details = csv_validate.process_csv(args.input, codec, opts, report.result, out_path=args.output)
        written = bool(report.details.get("output_written"))
        if written:
            report.add_output("csv", args.output)
    if args.report is None and not args.no_report:
        args.report = _default_report_path(args.output)
    rows = report.details.get("rows_written", 0) if written else 0
    _emit(report, args, f"[csv clean] {args.input} -> {args.output if written else '(no output)'}  rows: {rows}")
    return report.result.exit_code(args.strict)


def cmd_csv_validate(args: argparse.Namespace, argv: Sequence[str]) -> int:
    for f in args.files:
        _require_file(f)
    if args.against and len(args.files) != 1:
        raise InputError("--against cross-references one CSV file at a time; give a single input file")
    opts = _csv_options(args)
    report = Report(command=render_command(argv, args.report_paths), path_mode=args.report_paths)
    per_file: Dict[str, object] = {}
    worst = 0
    for path in args.files:
        report.add_input("csv", path)
        res = ValidationResult()
        textfile.check_extension(path, textfile.CSV_EXTENSIONS, "CSV", res)
        codec = textfile.precheck(path, args.encoding, res, "validate")
        details: Dict[str, object] = {}
        if codec is not None:
            hook = None
            if args.against:
                ref = csv_crosscheck.load_reference(args.against)
                report.add_input(ref.kind, args.against)
                hook = csv_crosscheck.CrossRef(ref)
            details = csv_validate.process_csv(path, codec, opts, res, row_hook=hook)
            if hook is not None:
                details["against"] = hook.report(res)
        per_file[os.path.basename(path) if args.report_paths == "basename" else path] = details
        if len(args.files) > 1:
            for f in res.findings:
                f.where = f"{os.path.basename(path)}:{f.where}" if f.where else os.path.basename(path)
        report.result.extend(res)
        worst = max(worst, res.exit_code(args.strict))
    report.details = {"files": per_file}
    _emit(report, args, f"[csv validate] {' '.join(args.files)}")
    return worst


def cmd_csv_inspect(args: argparse.Namespace, argv: Sequence[str]) -> int:
    import json
    opts = _csv_options(args)
    infos = []
    for path in args.files:
        _require_file(path)
        res = ValidationResult()
        textfile.check_extension(path, textfile.CSV_EXTENSIONS, "CSV", res)
        codec = textfile.precheck(path, args.encoding, res, "validate")
        if codec is None:
            raise InputError("; ".join(f.message for f in res.errors))
        info = csv_validate.inspect_csv(path, codec, opts.delimiter, res, opts.mapping)
        info["findings"] = res.to_list()
        infos.append(info)
    if args.json:
        print(json.dumps(infos if len(infos) > 1 else infos[0], indent=2, ensure_ascii=False))
    else:
        for info in infos:
            print(csv_validate.format_inspect(info))
    return 0


# -- entry point --------------------------------------------------------------------

def main(argv: Optional[Sequence[str]] = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    parser = build_parser()
    try:
        args = parser.parse_args(argv)
    except SystemExit as exc:          # argparse already printed the message
        return int(exc.code) if isinstance(exc.code, int) else EXIT_INPUT
    func: Callable[[argparse.Namespace, Sequence[str]], int] = args.func
    try:
        return func(args, argv)
    except InputError as exc:
        print(f"abd-prep: error: {exc}", file=sys.stderr)
        return EXIT_INPUT
    except UnicodeError as exc:
        print(f"abd-prep: error: text encoding problem: {exc}", file=sys.stderr)
        return EXIT_INPUT
    except OSError as exc:
        print(f"abd-prep: error: {exc}", file=sys.stderr)
        return EXIT_INPUT


if __name__ == "__main__":
    sys.exit(main())
