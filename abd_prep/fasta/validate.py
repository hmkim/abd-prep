"""FASTA validation rules F15-F17 (plan §3.3.2) and the ``inspect`` summary.

``validate_records`` is the single entry point used by ``fasta clean``,
``fasta validate`` and ``fasta inspect``: it assigns roles (F14), checks
the alphabet (F15), the framework consistency (F16) and, unless disabled,
the length / signal-peptide / duplicate heuristics (F17). Sequences are
never modified here except for the explicit ``--replace-ambiguous``
substitutions the user asked for.

Finding codes:

``FASTA_BAD_ALPHABET``          error    letters outside the 20 standard amino acids
``FASTA_NUCLEIC``               error    record looks like DNA/RNA (A/C/G/T/U/N only)
``FASTA_AMBIGUOUS_REPLACED``    info     ``--replace-ambiguous`` applied
``FASTA_FRAMEWORK_UNSPECIFIED`` info     --framework not given; F16 skipped
``FASTA_FRAMEWORK_MISMATCH``    error    record count / roles do not fit the framework
``FASTA_LENGTH_RANGE``          warning  outside the plausible range for the role
``FASTA_SIGNAL_PEPTIDE``        warning  starts with M + hydrophobic stretch
``FASTA_DUPLICATE_SEQUENCE``    warning  two records have the same sequence
``FASTA_SCFV_LINKER``           info     common (G4S)n-type linker found / not found
``FASTA_HEURISTICS_OFF``        info     --no-heuristics given
"""
import re
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

from ..common import residues as R
from ..common.result import InputError, ValidationResult
from .fastafile import Record
from .roles import ROLE_HEAVY, ROLE_LIGHT, ROLE_SCFV, assign_roles, check_framework

#: Plausible residue counts per role (plan §3.3.2). Wider than the PDB
#: rules (§3.3.1: VL 105-115, VH 110-125) on purpose: a seed sequence is
#: whatever the user has, not an Fv cut from a crystal structure, so
#: both ends vary more.
LENGTH_RANGE: Dict[str, Tuple[int, int]] = {
    ROLE_LIGHT: (95, 130),
    ROLE_HEAVY: (105, 140),
    "VHH": (100, 140),
    ROLE_SCFV: (220, 280),
}

#: Minimum length before the "looks like DNA/RNA" rule applies: short
#: peptides made of A/C/G/T letters are real proteins.
NUCLEIC_MIN_LENGTH = 20

#: scFv linker: two or more (G3-4 S) units, e.g. (G4S)3 = GGGGSGGGGSGGGGS.
LINKER_RE = re.compile(r"(?:G{3,4}S){2,}")

#: Signal-peptide heuristic: sequence starts with M and, within the first
#: 30 residues, has a window of 8 with at least 7 hydrophobic residues.
HYDROPHOBIC = frozenset("AILMFVWC")
SIGNAL_SCAN = 30
SIGNAL_WINDOW = 8
SIGNAL_MIN_HYDROPHOBIC = 7


@dataclass
class FastaOptions:
    framework: Optional[str] = None             # igg | nanobody | scfv | None
    heavy_id: Optional[str] = None
    light_id: Optional[str] = None
    replace_ambiguous: Dict[str, str] = field(default_factory=dict)   # e.g. {"X": "G"}
    heuristics: bool = True


def parse_replacements(text: Optional[str]) -> Dict[str, str]:
    """``"X=G,B=N"`` -> ``{"X": "G", "B": "N"}`` (F15 explicit substitutions)."""
    out: Dict[str, str] = {}
    if not text:
        return out
    for item in text.split(","):
        item = item.strip()
        if not item:
            continue
        if "=" not in item:
            raise InputError(f"--replace-ambiguous: expected FROM=TO, got {item!r}")
        src, dst = (p.strip().upper() for p in item.split("=", 1))
        if len(src) != 1 or len(dst) != 1:
            raise InputError(f"--replace-ambiguous: single letters only, got {item!r}")
        if not src.isalpha():
            raise InputError(f"--replace-ambiguous: {src!r} is not a letter; an internal '*' is an error to fix at the "
                             "source, not a character to replace")
        if src in R.STANDARD_AA1:
            raise InputError(f"--replace-ambiguous: {src} is a standard amino acid and is never replaced")
        if dst not in R.STANDARD_AA1:
            raise InputError(f"--replace-ambiguous: replacement {dst!r} is not one of the 20 standard amino acids")
        out[src] = dst
    return out


def looks_like_signal_peptide(seq: str) -> bool:
    if not seq.startswith("M") or len(seq) < SIGNAL_WINDOW + 1:
        return False
    head = seq[1:SIGNAL_SCAN]
    for i in range(0, len(head) - SIGNAL_WINDOW + 1):
        window = head[i:i + SIGNAL_WINDOW]
        if sum(1 for ch in window if ch in HYDROPHOBIC) >= SIGNAL_MIN_HYDROPHOBIC:
            return True
    return False


def _positions(items: List[Tuple[int, str]], limit: int = 8) -> str:
    shown = ", ".join(f"{ch}@{pos}" for pos, ch in items[:limit])
    return shown + (f" ... ({len(items)} total)" if len(items) > limit else "")


def check_alphabet(rec: Record, result: ValidationResult, replace: Optional[Dict[str, str]] = None) -> bool:
    """F15 for one record; returns True when the sequence is acceptable."""
    where = rec.id or f"line {rec.lineno}"
    seq = rec.sequence
    if replace:
        counts: Dict[str, int] = {}
        for src, dst in replace.items():
            n = seq.count(src)
            if n:
                counts[src] = n
                seq = seq.replace(src, dst)
        if counts:
            rec.sequence = seq
            desc = ", ".join(f"{src}->{replace[src]} x{n}" for src, n in counts.items())
            result.info("FASTA_AMBIGUOUS_REPLACED", f"explicit replacement applied: {desc}", where)
    bad = [(i + 1, ch) for i, ch in enumerate(seq) if ch not in R.STANDARD_AA1 and ch != "*"]
    _, nucleic = R.check_protein_alphabet(seq.replace("*", ""))
    if nucleic and len(seq) >= NUCLEIC_MIN_LENGTH:
        result.error("FASTA_NUCLEIC",
                     f"{len(seq)} residues made only of A/C/G/T/U/N: this looks like a DNA/RNA sequence; "
                     "ABD needs the translated amino-acid sequence", where)
        return False
    if bad:
        result.error("FASTA_BAD_ALPHABET",
                     f"non-standard letter(s) {_positions(bad)}: only the 20 standard amino acids are "
                     "accepted (X/B/Z/J/U/O are ambiguous). Replace them explicitly with "
                     "--replace-ambiguous X=G or fix the source sequence", where)
        return False
    return True


def length_key(rec: Record, framework: Optional[str]) -> Optional[str]:
    if rec.role == ROLE_HEAVY and framework == "nanobody":
        return "VHH"
    return rec.role if rec.role in LENGTH_RANGE else None


def check_heuristics(records: Sequence[Record], framework: Optional[str], result: ValidationResult) -> None:
    """F17: length range, signal peptide, duplicate sequences (warnings only)."""
    seen: Dict[str, str] = {}
    for rec in records:
        where = rec.id or f"line {rec.lineno}"
        seq = rec.sequence
        if not seq:
            continue
        key = length_key(rec, framework)
        if key is not None:
            lo, hi = LENGTH_RANGE[key]
            label = "VHH" if key == "VHH" else {"H": "VH", "L": "VL", "HL": "scFv"}[key]
            if len(seq) < lo:
                result.warning("FASTA_LENGTH_RANGE",
                               f"{len(seq)} residues is shorter than a typical {label} ({lo}-{hi}); truncated?", where)
            elif len(seq) > hi:
                result.warning("FASTA_LENGTH_RANGE",
                               f"{len(seq)} residues is longer than a typical {label} ({lo}-{hi}); constant region "
                               "or signal peptide included? ABD seeds are variable domains", where)
        if looks_like_signal_peptide(seq):
            result.warning("FASTA_SIGNAL_PEPTIDE",
                           "starts with M followed by a hydrophobic stretch: signal peptide not removed? Mature "
                           "variable domains usually start with E/Q/D", where)
        if seq in seen:
            result.warning("FASTA_DUPLICATE_SEQUENCE", f"same sequence as record {seen[seq]!r}", where)
        else:
            seen[seq] = rec.id


def check_scfv_linker(records: Sequence[Record], framework: Optional[str], result: ValidationResult) -> None:
    for rec in records:
        if framework != "scfv" and rec.role != ROLE_SCFV:
            continue
        where = rec.id or f"line {rec.lineno}"
        m = LINKER_RE.search(rec.sequence)
        if m:
            result.info("FASTA_SCFV_LINKER",
                        f"G/S linker {m.group(0)} at positions {m.start() + 1}-{m.end()} (kept as-is)", where)
        else:
            result.info("FASTA_SCFV_LINKER", "no common (G3-4S)n linker pattern found; other linkers are fine, "
                        "the sequence is left unchanged", where)


def validate_records(records: Sequence[Record], opts: Optional[FastaOptions], result: ValidationResult) -> None:
    """Run F14-F17 over parsed records, recording into ``result``."""
    opts = opts or FastaOptions()
    assign_roles(records, result, opts.framework, opts.heavy_id, opts.light_id)
    for rec in records:
        check_alphabet(rec, result, opts.replace_ambiguous)
    check_framework(records, opts.framework, result)
    check_scfv_linker(records, opts.framework, result)
    if opts.heuristics:
        check_heuristics(records, opts.framework, result)
    else:
        result.info("FASTA_HEURISTICS_OFF", "length/signal-peptide/duplicate heuristics disabled", "file")


# -- inspect -----------------------------------------------------------------------

def inspect_records(records: Sequence[Record], path: str = "") -> Dict[str, object]:
    rows = []
    for rec in records:
        bad, nucleic = R.check_protein_alphabet(rec.sequence.replace("*", ""))
        rows.append({
            "id": rec.id, "description": rec.description, "role": rec.role, "role_source": rec.role_source,
            "length": rec.length, "nonstandard": len(bad), "looks_nucleic": nucleic and rec.length >= NUCLEIC_MIN_LENGTH,
            "internal_stop": "*" in rec.sequence, "signal_peptide": looks_like_signal_peptide(rec.sequence),
            "sequence": rec.sequence,
        })
    return {"path": path, "records": rows}


def format_inspect(info: Dict[str, object]) -> str:
    rows: List[Dict[str, object]] = info["records"]  # type: ignore[assignment]
    lines = [f"== {info['path']}  records: {len(rows)}"]
    for r in rows:
        flags = []
        if r["nonstandard"]:
            flags.append(f"non-standard letters: {r['nonstandard']}")
        if r["looks_nucleic"]:
            flags.append("looks like DNA/RNA")
        if r["internal_stop"]:
            flags.append("internal '*'")
        if r["signal_peptide"]:
            flags.append("signal peptide?")
        role = r["role"] or "?"
        src = f" ({r['role_source']})" if r["role_source"] else ""
        extra = f"  [{'; '.join(flags)}]" if flags else ""
        lines.append(f"  {r['id'] or '(no id)'}: role {role}{src}, {r['length']} res{extra}")
        lines.append(f"     {r['sequence']}")
    return "\n".join(lines)
