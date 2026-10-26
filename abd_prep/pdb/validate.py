"""Validation rules (plan §3.3.1) and the ``inspect`` summary (F7, F8).

Each rule emits one finding code so tests, README and ``report.json``
can refer to the same name:

======================= ======== =================================================
code                    level    rule
======================= ======== =================================================
PDB_NO_ATOMS            error    no ATOM/HETATM records
PDB_ALTLOC_PRESENT      error    column 17 not blank (ABD: "contains alternate locations")
PDB_END_MISSING         error    no END record
PDB_CHAIN_COUNT         error    chain count != ``expect_chains``
PDB_TER_MISSING         error    chain boundary without a TER
PDB_TER_MISPLACED       error    TER not at a chain boundary (left over from a removed chain)
PDB_TER_MISMATCH        error    TER names a different residue than the atom before it
PDB_TER_EXTRA           error    several TERs at one boundary
PDB_HETATM_PRESENT      warning  HETATM records remain (info with ``allow_het``)
PDB_GAP_CA_DISTANCE     warning  consecutive residues with CA-CA > 4.2 A (info if expected)
PDB_NONSTANDARD_RESNAME warning  ATOM record with a non-standard residue name
PDB_RESIDUE_COUNT       warning  chain H/L (or an expectation) outside the plausible Fv size
PDB_MULTI_MODEL         warning  more than one MODEL block
PDB_CHAIN_BLANK         warning  blank chain ID
PDB_SERIAL_DUPLICATE    warning  duplicate atom/TER serial numbers
PDB_RESIDUE_DUPLICATE   warning  same (resSeq, iCode) twice in a chain with different names
PDB_NON_ASCII           warning  ATOM/HETATM record with non-ASCII bytes
PDB_GAP_NUMBERING       info     residue numbers skip or decrease (Kabat codes, disordered loops ...)
PDB_HYDROGENS_PRESENT   info     hydrogen atoms present
PDB_EXTENSION           warning  file name does not end in .pdb/.ent (CLI)
======================= ======== =================================================
"""
import math
import os
from collections import Counter, OrderedDict
from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple

from ..common import residues as R
from ..common.result import ValidationResult
from ..fasta.roles import fasta_header
from .pdbfile import Atom, PdbFile, iter_residues, read_pdb

#: Consecutive-residue CA-CA distance above which the chain is treated as
#: broken (a peptide bond gives ~3.8 A; 4.2 A leaves room for strain).
CA_GAP_ANGSTROM = 4.2

#: Plausible Fv residue counts for chains labelled H/L (plan §3.3.1).
FV_RESIDUE_RANGE: Dict[str, Tuple[int, int]] = {"L": (105, 115), "H": (110, 125)}


@dataclass
class ValidateOptions:
    allow_het: bool = False
    expect_chains: Optional[int] = None
    #: chain -> (min, max) residue count; overrides FV_RESIDUE_RANGE (presets, M2)
    expect_residues: Dict[str, Tuple[int, int]] = field(default_factory=dict)
    #: expected gaps as ``{(chain, resseq_before, resseq_after)}`` (presets, M2)
    expected_gaps: Set[Tuple[str, int, int]] = field(default_factory=set)


@dataclass
class ChainSummary:
    chain: str
    first: str
    last: str
    residues: int
    atoms: int
    sequence: str
    nonstandard: List[str]

    def to_dict(self) -> Dict[str, object]:
        return {"chain": self.chain, "first": self.first, "last": self.last, "residues": self.residues,
                "atoms": self.atoms, "sequence": self.sequence, "nonstandard": self.nonstandard}


def summarize_chains(atoms: Iterable[Atom]) -> "OrderedDict[str, ChainSummary]":
    """Per-chain range, counts and one-letter sequence over ATOM records."""
    out: "OrderedDict[str, ChainSummary]" = OrderedDict()
    seqs: Dict[str, List[str]] = {}
    for key, resname, res_atoms in iter_residues(a for a in atoms if not a.is_hetatm):
        ch = key[0]
        label = f"{key[1]}{key[2].strip()}"
        s = out.get(ch)
        if s is None:
            s = out[ch] = ChainSummary(ch, label, label, 0, 0, "", [])
            seqs[ch] = []
        s.last = label
        s.residues += 1
        s.atoms += len(res_atoms)
        seqs[ch].append(R.one_letter(resname))
        if resname not in R.STANDARD_AA3 and resname not in s.nonstandard:
            s.nonstandard.append(resname)
    for ch, s in out.items():
        s.sequence = "".join(seqs[ch])
    return out


def _ca(res_atoms: Sequence[Atom]) -> Optional[Tuple[float, float, float]]:
    for a in res_atoms:
        if a.name == "CA":
            return a.xyz
    return None


def _dist(p: Tuple[float, float, float], q: Tuple[float, float, float]) -> float:
    return math.sqrt(sum((x - y) ** 2 for x, y in zip(p, q)))


def validate_atoms(pdb: PdbFile, opts: Optional[ValidateOptions] = None) -> Tuple[ValidationResult, Dict[str, object]]:
    """Apply every §3.3.1 rule to a parsed file."""
    opts = opts or ValidateOptions()
    res = ValidationResult()
    atoms = pdb.first_model_atoms()
    details: Dict[str, object] = {"models": pdb.n_models, "ter_records": len(pdb.ters), "has_end": pdb.has_end}

    if not atoms:
        res.error("PDB_NO_ATOMS", "no ATOM/HETATM records", "file")
        return res, details

    if pdb.n_models > 1:
        res.warning("PDB_MULTI_MODEL", f"{pdb.n_models} MODEL blocks; ABD expects a single model "
                    "(clean keeps only the first)", "file")

    # altloc -----------------------------------------------------------------
    alt_atoms = [a for a in atoms if a.altloc != " "]
    if alt_atoms:
        residues = sorted({f"{a.chain}:{a.resseq}{a.icode.strip()}" for a in alt_atoms})
        shown = ", ".join(residues[:5]) + (" ..." if len(residues) > 5 else "")
        res.error("PDB_ALTLOC_PRESENT",
                  f"{len(alt_atoms)} atoms in {len(residues)} residues carry alternate locations "
                  f"(ABD rejects: 'PDB contains alternate locations'); residues: {shown}", "file")

    # HETATM -----------------------------------------------------------------
    het = Counter(a.resname for a in atoms if a.is_hetatm)
    details["hetatm"] = dict(sorted(het.items()))
    if het:
        summary = ", ".join(f"{k} x{v}" for k, v in sorted(het.items()))
        if opts.allow_het:
            res.info("PDB_HETATM_PRESENT", f"HETATM records present (allowed): {summary}", "file")
        else:
            res.warning("PDB_HETATM_PRESENT", f"HETATM records present: {summary}; ABD inputs normally contain "
                        "protein atoms only (clean removes them by default)", "file")

    # chains -----------------------------------------------------------------
    summaries = summarize_chains(atoms)
    details["chains"] = OrderedDict((ch, s.to_dict()) for ch, s in summaries.items())
    chain_ids = list(summaries)
    if " " in {a.chain for a in atoms}:
        res.warning("PDB_CHAIN_BLANK", "some atoms have a blank chain ID; ABD needs labelled chains", "file")
    if opts.expect_chains is not None and len(chain_ids) != opts.expect_chains:
        res.error("PDB_CHAIN_COUNT", f"{len(chain_ids)} protein chain(s) ({''.join(chain_ids)}), "
                  f"expected {opts.expect_chains}", "file")

    for ch, s in summaries.items():
        if s.nonstandard:
            res.warning("PDB_NONSTANDARD_RESNAME",
                        f"non-standard residue name(s) in ATOM records: {', '.join(s.nonstandard)}", ch)
        rng = opts.expect_residues.get(ch) or FV_RESIDUE_RANGE.get(ch)
        if rng and not (rng[0] <= s.residues <= rng[1]):
            res.warning("PDB_RESIDUE_COUNT",
                        f"chain {ch} has {s.residues} residues, expected {rng[0]}-{rng[1]}", ch)

    # gaps, duplicates -----------------------------------------------------
    prev: Dict[str, Tuple[Tuple[str, int, str], Optional[Tuple[float, float, float]]]] = {}
    seen_keys: Dict[str, Dict[Tuple[int, str], str]] = {}
    for key, resname, res_atoms in iter_residues(a for a in atoms if not a.is_hetatm):
        ch, num, ic = key
        bucket = seen_keys.setdefault(ch, {})
        if (num, ic) in bucket and bucket[(num, ic)] != resname:
            res.warning("PDB_RESIDUE_DUPLICATE",
                        f"residue {num}{ic.strip()} appears twice with different names "
                        f"({bucket[(num, ic)]}, {resname})", ch)
        bucket.setdefault((num, ic), resname)
        ca = _ca(res_atoms)
        if ch in prev:
            (pch, pnum, pic), pca = prev[ch]
            where = f"{ch}:{pnum}{pic.strip()}-{num}{ic.strip()}"
            expected = (ch, pnum, num) in opts.expected_gaps
            if num - pnum > 1:
                res.info("PDB_GAP_NUMBERING",
                         f"residue numbering jumps {pnum} -> {num}" + (" (expected)" if expected else ""), where)
            elif num < pnum:
                res.info("PDB_GAP_NUMBERING",
                         f"residue numbering decreases {pnum} -> {num} (chimeric construct or renumbered segment?)",
                         where)
            if ca is not None and pca is not None:
                d = _dist(pca, ca)
                if d > CA_GAP_ANGSTROM:
                    level = "info" if expected else "warning"
                    res.add(level, "PDB_GAP_CA_DISTANCE",
                            f"CA-CA distance {d:.1f} A between consecutive residues (> {CA_GAP_ANGSTROM} A): "
                            f"chain break" + (" -- expected for this preset" if expected else ""), where)
        prev[ch] = (key, ca)

    # TER consistency (P1 regression guard) ----------------------------------
    # after_index counts atoms across every model, so the check is only
    # meaningful for single-model files (clean output always is).
    if pdb.n_models <= 1:
        _check_ters(pdb, atoms, res)

    # serials, END, hydrogens -------------------------------------------------
    serials = Counter(a.serial for a in atoms if a.serial is not None)
    serials.update(t.serial for t in pdb.ters if t.serial is not None)
    dup = sorted(s for s, n in serials.items() if n > 1)
    if dup:
        shown = ", ".join(str(s) for s in dup[:5]) + (" ..." if len(dup) > 5 else "")
        res.warning("PDB_SERIAL_DUPLICATE", f"{len(dup)} duplicate serial number(s): {shown} "
                    "(clean --renumber rewrites them)", "file")
    if not pdb.has_end:
        res.error("PDB_END_MISSING", "file does not end with an END record", "file")
    n_h = sum(1 for a in atoms if a.is_hydrogen)
    if n_h:
        res.info("PDB_HYDROGENS_PRESENT", f"{n_h} hydrogen atoms present (clean --strip-h removes them)", "file")
    n_bad = sum(1 for a in atoms if not a.is_ascii)
    if n_bad:
        res.warning("PDB_NON_ASCII", f"{n_bad} ATOM/HETATM record(s) contain non-ASCII bytes; PDB is an ASCII "
                    "format and parsers may reject them", "file")

    details["atoms"] = len(atoms)
    return res, details


def _check_ters(pdb: PdbFile, atoms: List[Atom], res: ValidationResult) -> None:
    """TER must sit at every chain boundary, once, naming the preceding residue."""
    boundaries: List[Tuple[int, Atom]] = []
    for i in range(1, len(atoms)):
        if atoms[i].chain != atoms[i - 1].chain:
            boundaries.append((i, atoms[i - 1]))
    boundaries.append((len(atoms), atoms[-1]))
    expected_at = {pos: last for pos, last in boundaries}
    ters_at: Counter = Counter(t.after_index for t in pdb.ters)
    for pos, last in boundaries:
        if pos not in ters_at:
            res.error("PDB_TER_MISSING", f"no TER after residue {last.resseq}{last.icode.strip()}", last.chain)
    for t in pdb.ters:
        if t.after_index not in expected_at:
            res.error("PDB_TER_MISPLACED",
                      f"TER ({t.resname} {t.chain} {t.resseq}) is not at a chain boundary -- left over from a "
                      f"removed chain or cut residue", f"line {t.lineno}")
            continue
        last = expected_at[t.after_index]
        if (t.resname, t.chain, t.resseq, t.icode) != (last.resname, last.chain, last.resseq, last.icode):
            res.error("PDB_TER_MISMATCH",
                      f"TER says {t.resname} {t.chain} {t.resseq}{t.icode.strip()} but the preceding atom is "
                      f"{last.resname} {last.chain} {last.resseq}{last.icode.strip()}", f"line {t.lineno}")
    for pos, n in ters_at.items():
        if n > 1:
            where = expected_at[pos].chain if pos in expected_at else f"after atom {pos}"
            res.error("PDB_TER_EXTRA", f"{n} TER records at one boundary", where)


def validate_file(path: str, opts: Optional[ValidateOptions] = None) -> Tuple[ValidationResult, Dict[str, object]]:
    return validate_atoms(read_pdb(path), opts)


# -- inspect (prototype inspect_pdb.py) --------------------------------------

def inspect_atoms(pdb: PdbFile) -> Dict[str, object]:
    atoms = pdb.first_model_atoms()
    summaries = summarize_chains(atoms)
    het = Counter(a.resname for a in atoms if a.is_hetatm)
    return {
        "path": pdb.path,
        "models": pdb.n_models,
        "altloc_atoms": sum(1 for a in atoms if a.altloc != " "),
        "hetatm": dict(sorted(het.items())),
        "ter_records": len(pdb.ters),
        "has_end": pdb.has_end,
        "chains": OrderedDict((ch, s.to_dict()) for ch, s in summaries.items()),
    }


def inspect_file(path: str) -> Dict[str, object]:
    return inspect_atoms(read_pdb(path))


def format_inspect(info: Dict[str, object]) -> str:
    """Text layout close to the prototype's ``inspect_pdb.py`` output."""
    het = info["hetatm"]
    het_text = ", ".join(f"{k} x{v}" for k, v in het.items()) if het else "none"
    lines = [f"== {info['path']}  altloc atoms: {info['altloc_atoms']}  HETATM: {het_text}  "
             f"MODEL: {info['models']}  TER: {info['ter_records']}  END: {'yes' if info['has_end'] else 'no'}"]
    for ch, s in info["chains"].items():
        extra = f"  non-standard: {','.join(s['nonstandard'])}" if s["nonstandard"] else ""
        lines.append(f"  chain {ch}: {s['first']}-{s['last']} ({s['residues']} res, {s['atoms']} atoms){extra}")
        lines.append(f"     {s['sequence']}")
    return "\n".join(lines)


def to_fasta(info: Dict[str, object], width: int = 0) -> str:
    """FASTA text for every chain (F8); headers follow :mod:`abd_prep.fasta.roles`."""
    out: List[str] = []
    src = os.path.basename(info.get("path") or "")
    for ch, s in info["chains"].items():
        desc = f"{src} chain {ch} {s['first']}-{s['last']} ({s['residues']} res)".strip()
        out.append(fasta_header(ch, desc))
        seq = s["sequence"]
        if width and width > 0:
            out.extend(seq[i:i + width] for i in range(0, len(seq), width))
        else:
            out.append(seq)
    return "\n".join(out) + ("\n" if out else "")
