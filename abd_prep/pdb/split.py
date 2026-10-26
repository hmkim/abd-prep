"""``abd-prep pdb split``: fetch -> select -> clean -> write -> validate (plan §3.1).

One run produces up to two files in ``--out``:

* ``<name>_target.pdb``  the antigen part
* ``<name>_fv.pdb``      heavy + light chains (same source entry)

plus ``<name>.report.json`` with the §3.5 schema. ``<name>`` is the preset
name, ``--name``, or derived from the entry/file -- never a fixed string
like ``trastuzumab_fv.pdb`` (plan R5: the 6B9Z antibody may not be the
reference trastuzumab sequence, so a fixed name could lie).
"""
import os
from dataclasses import dataclass, field
from typing import Dict, FrozenSet, List, Optional, Sequence, Tuple

from ..common.report import Report
from ..common.result import InputError, ValidationResult
from .clean import CleanOptions, clean_file
from .fetch import RECOMMENDED_RESOLUTION, FetchResult, Opener, default_opener, fetch, normalize_pdb_id
from .presets import Part, Preset, output_filename, output_groups, validate_parts
from .ranges import ChainSelection
from .validate import ValidateOptions, validate_file

LOCAL_SOURCE = "input"


@dataclass
class SplitOptions:
    name: str
    parts: Tuple[Part, ...]
    out_dir: str = "."
    preset: Optional[Preset] = None
    input_path: Optional[str] = None          # source for parts whose pdb_id == "input"
    expect_residues: Dict[str, int] = field(default_factory=dict)
    expected_gaps: FrozenSet[Tuple[str, int, int]] = frozenset()
    cache: Optional[str] = None
    refresh: bool = False
    keep_het: bool = False
    std_residues: bool = False
    strip_h: bool = False
    renumber: bool = False
    path_mode: str = "basename"
    report_path: Optional[str] = None         # None -> <out>/<name>.report.json; "" -> none
    #: group ("target"/"fv") -> file name inside ``out_dir``; default ``<name>_<group>.pdb``
    filenames: Dict[str, str] = field(default_factory=dict)

    @classmethod
    def from_preset(cls, preset: Preset, out_dir: str, **kw) -> "SplitOptions":
        return cls(name=preset.name, parts=preset.parts, out_dir=out_dir, preset=preset,
                   expect_residues=dict(preset.expect_residues), expected_gaps=preset.expected_gaps, **kw)


@dataclass
class SplitOutcome:
    report: Report
    written: Dict[str, str]                      # group ("target"/"fv") -> path
    sources: Dict[str, str]                      # pdb_id -> local path used
    report_path: Optional[str] = None

    @property
    def result(self) -> ValidationResult:
        return self.report.result


def derive_name(parts: Sequence[Part], input_path: Optional[str]) -> str:
    """Default ``<name>`` when neither a preset nor ``--name`` is given."""
    ids = []
    for p in parts:
        if p.pdb_id not in ids:
            ids.append(p.pdb_id)
    if ids == [LOCAL_SOURCE] and input_path:
        return os.path.splitext(os.path.basename(input_path))[0]
    return "-".join(i.lower() for i in ids if i != LOCAL_SOURCE) or "split"


def _prefix(result: ValidationResult, group: str) -> None:
    for f in result.findings:
        f.where = f"{group}/{f.where}" if f.where else group


def resolve_sources(opts: SplitOptions, report: Report, opener: Opener) -> Dict[str, str]:
    """Fetch (or locate) every distinct source entry; record inputs + metadata."""
    paths: Dict[str, str] = {}
    sources_detail: List[Dict[str, object]] = []
    for pid in dict.fromkeys(p.pdb_id for p in opts.parts):
        if pid == LOCAL_SOURCE:
            if not opts.input_path:
                raise InputError("a part refers to the local input but --input was not given")
            if not os.path.isfile(opts.input_path):
                raise InputError(f"file not found: {opts.input_path}")
            paths[pid] = opts.input_path
            sources_detail.append({"pdb_id": None, "path": os.path.basename(opts.input_path)
                                   if opts.path_mode == "basename" else os.path.abspath(opts.input_path)})
        else:
            fr: FetchResult = fetch(pid, opts.cache, opts.refresh, opener)
            paths[pid] = fr.path
            sources_detail.append(fr.to_dict(opts.path_mode))
            if fr.meta_error:
                report.result.info("PDB_ENTRY_META_UNAVAILABLE",
                                   f"RCSB entry metadata not recorded: {fr.meta_error}", pid)
            elif fr.meta and fr.meta.resolution is not None and fr.meta.resolution > RECOMMENDED_RESOLUTION:
                report.result.warning("PDB_RESOLUTION",
                                      f"resolution {fr.meta.resolution:.2f} A is coarser than the "
                                      f"{RECOMMENDED_RESOLUTION} A the ABD guide recommends", pid)
        report.add_input("pdb", paths[pid])
    report.details["sources"] = sources_detail
    return paths


def run_split(opts: SplitOptions, opener: Opener = default_opener, command: str = "") -> SplitOutcome:
    validate_parts(opts.parts)
    for p in opts.parts:
        if p.pdb_id != LOCAL_SOURCE:
            normalize_pdb_id(p.pdb_id)
    if not opts.name or os.sep in opts.name or opts.name in (".", ".."):
        raise InputError(f"bad output name {opts.name!r}")
    for group, fname in opts.filenames.items():
        if not fname or os.sep in fname or fname in (".", ".."):
            raise InputError(f"bad output file name {fname!r} for {group}")
    os.makedirs(opts.out_dir, exist_ok=True)

    report = Report(command=command, path_mode=opts.path_mode)
    report.details = {
        "name": opts.name,
        "preset": opts.preset.to_dict() if opts.preset else None,
        "parts": {},
    }
    sources = resolve_sources(opts, report, opener)
    written: Dict[str, str] = {}

    for group, parts in output_groups(opts.parts).items():
        src_id = parts[0].pdb_id
        selection = ChainSelection.parse([p.selection for p in parts])
        rename = {p.chain: p.out_chain for p in parts if p.out_chain and p.out_chain != p.chain}
        copts = CleanOptions(selection=selection, rename=rename, keep_het=opts.keep_het,
                             std_residues=opts.std_residues, strip_h=opts.strip_h, renumber=opts.renumber)
        out_path = os.path.join(opts.out_dir, opts.filenames.get(group) or output_filename(opts.name, group))
        if os.path.abspath(out_path) == os.path.abspath(sources[src_id]):
            raise InputError(f"output {out_path} would overwrite the source file")
        outcome = clean_file(sources[src_id], out_path, copts)
        _prefix(outcome.result, group)
        report.result.extend(outcome.result)

        out_chains = [p.written_chain for p in parts]
        vopts = ValidateOptions(
            expect_chains=len(parts),
            expect_residues={ch: (n, n) for ch, n in opts.expect_residues.items() if ch in out_chains},
            expected_gaps={g for g in opts.expected_gaps if g[0] in out_chains},
        )
        vres, vdetails = validate_file(out_path, vopts)
        _prefix(vres, group)
        report.result.extend(vres)

        report.details["parts"][group] = {
            "source": None if src_id == LOCAL_SOURCE else src_id,
            "roles": {p.role: {"selection": p.selection, "out_chain": p.written_chain} for p in parts},
            "file": os.path.basename(out_path) if opts.path_mode == "basename" else os.path.abspath(out_path),
            "clean": outcome.details,
            "validation": vdetails,
        }
        report.add_output("pdb", out_path)
        written[group] = out_path

    report_path: Optional[str]
    if opts.report_path == "":
        report_path = None
    else:
        report_path = opts.report_path or os.path.join(opts.out_dir, f"{opts.name}.report.json")
        report.write(report_path)
    return SplitOutcome(report, written, sources, report_path)
