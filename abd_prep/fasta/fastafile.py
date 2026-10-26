"""FASTA records: parsing, F13 normalization and writing.

A record is a ``>`` header line followed by one or more sequence lines.
Normalization (plan F13) touches *formatting only*; it never changes,
cuts or re-orders residues:

* CR/LF line endings and a UTF-8 BOM are tolerated (reported as info)
* white space and digits inside sequence lines are removed (GenPept and
  Word paste-ins carry position numbers)
* alignment gaps ``-`` and ``.`` are removed
* lower-case letters become upper-case
* a ``*`` at the very end (translation stop) is removed; a ``*`` anywhere
  else is an **error** -- it marks an internal stop codon, i.e. a wrong
  translation, and silently deleting it would let a broken sequence pass

Finding codes:

``FASTA_NO_RECORDS``          error    no ``>`` header in the file
``FASTA_NO_HEADER``           error    sequence text before the first header
``FASTA_EMPTY_HEADER``        error    ``>`` with no record ID
``FASTA_EMPTY_SEQUENCE``      error    header without residues
``FASTA_TOO_LONG``            error    more than ``MAX_SEQUENCE_LENGTH`` residues
``FASTA_INTERNAL_STOP``       error    ``*`` inside the sequence
``FASTA_DUPLICATE_ID``        warning  two records share an ID
``FASTA_CRLF``                info     Windows line endings
``FASTA_LOWERCASE``           info     lower-case letters upper-cased
``FASTA_GAPS_REMOVED``        info     ``-``/``.`` removed
``FASTA_WHITESPACE_REMOVED``  info     blanks/digits removed inside sequence lines
``FASTA_TRAILING_STOP``       info     terminal ``*`` removed
"""
import os
import re
from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Optional, Sequence, TextIO, Tuple, Union

from ..common.result import ValidationResult
from .roles import fasta_header

#: Longer than any antibody construct; plan §2 non-functional requirements.
MAX_SEQUENCE_LENGTH = 10_000

_WS_DIGITS_RE = re.compile(r"[\s\d]+")
_GAP_RE = re.compile(r"[-.]+")


@dataclass
class Record:
    """One FASTA record after normalization."""

    id: str                       # first token after '>'
    description: str              # rest of the header line
    sequence: str                 # normalized residues (may still contain '*' when internal)
    lineno: int = 0               # line of the header
    role: Optional[str] = None    # 'H', 'L' or 'HL' (scFv), assigned by roles.assign_roles
    role_source: str = ""         # 'option', 'header', 'framework'

    @property
    def header(self) -> str:
        return f"{self.id} {self.description}".strip()

    @property
    def length(self) -> int:
        return len(self.sequence)

    def to_dict(self) -> Dict[str, object]:
        return {"id": self.id, "description": self.description, "length": self.length,
                "role": self.role, "role_source": self.role_source, "lineno": self.lineno}


@dataclass
class FastaFile:
    records: List[Record] = field(default_factory=list)
    path: str = ""
    crlf_lines: int = 0

    def by_id(self) -> Dict[str, Record]:
        out: Dict[str, Record] = {}
        for r in self.records:
            out.setdefault(r.id, r)
        return out


def read_fasta(source: Union[str, TextIO, Iterable[str]], result: ValidationResult, path: str = "",
               encoding: str = "utf-8-sig") -> FastaFile:
    """Parse and normalize a FASTA file, recording F13 findings in ``result``."""
    if isinstance(source, str):
        with open(source, "r", encoding=encoding, newline="") as fh:
            return _read_lines(fh, result, path or source)
    return _read_lines(source, result, path)


def _read_lines(lines: Iterable[str], result: ValidationResult, path: str) -> FastaFile:
    fasta = FastaFile(path=path)
    header: Optional[Tuple[str, str, int]] = None   # (id, description, lineno)
    chunks: List[str] = []
    orphan_lines: List[int] = []
    for lineno, raw in enumerate(lines, start=1):
        if raw.endswith("\r\n") or raw.endswith("\r"):
            fasta.crlf_lines += 1
        line = raw.rstrip("\r\n")
        if lineno == 1 and line.startswith("\ufeff"):
            line = line[1:]
        if line.startswith(">"):
            if header is not None:
                fasta.records.append(_finish(header, chunks, result))
            parts = line[1:].strip().split(None, 1)       # ID = first token (space OR tab separated)
            ident = parts[0] if parts else ""
            desc = parts[1].strip() if len(parts) > 1 else ""
            header = (ident, desc, lineno)
            chunks = []
        elif line.strip():
            if header is None:
                orphan_lines.append(lineno)
            else:
                chunks.append(line)
    if header is not None:
        fasta.records.append(_finish(header, chunks, result))

    if orphan_lines:
        result.error("FASTA_NO_HEADER",
                     f"sequence text before the first '>' header (line {orphan_lines[0]}"
                     f"{' and ' + str(len(orphan_lines) - 1) + ' more' if len(orphan_lines) > 1 else ''}); "
                     "add a header line such as '>heavy'", "file")
    if not fasta.records:
        result.error("FASTA_NO_RECORDS", "no FASTA records found (no line starts with '>')", "file")
    if fasta.crlf_lines:
        result.info("FASTA_CRLF", f"{fasta.crlf_lines} line(s) end with CR/LF (Windows); written as LF", "file")

    seen: Dict[str, int] = {}
    for r in fasta.records:
        seen[r.id] = seen.get(r.id, 0) + 1
    for ident, n in seen.items():
        if n > 1 and ident:
            result.warning("FASTA_DUPLICATE_ID", f"record ID {ident!r} appears {n} times; --heavy-id/--light-id "
                           "cannot tell them apart", ident)
    return fasta


def _finish(header: Tuple[str, str, int], chunks: Sequence[str], result: ValidationResult) -> Record:
    ident, desc, lineno = header
    where = ident or f"line {lineno}"
    if not ident:
        result.error("FASTA_EMPTY_HEADER", f"header on line {lineno} has no record ID", where)
    raw = "".join(chunks)
    seq, n_ws = _WS_DIGITS_RE.subn("", raw)
    if n_ws:
        result.info("FASTA_WHITESPACE_REMOVED", f"{n_ws} run(s) of blanks/digits removed inside the sequence", where)
    seq, n_gap = _GAP_RE.subn("", seq)
    if n_gap:
        result.info("FASTA_GAPS_REMOVED", f"alignment gap character(s) '-'/'.' removed ({n_gap} run(s))", where)
    n_lower = sum(1 for ch in seq if ch.islower())
    if n_lower:
        result.info("FASTA_LOWERCASE", f"{n_lower} lower-case letter(s) upper-cased", where)
        seq = seq.upper()
    stripped = seq.rstrip("*")
    if len(stripped) != len(seq):
        n_stop = len(seq) - len(stripped)
        seq = stripped
        result.info("FASTA_TRAILING_STOP",
                    "terminal '*' (translation stop) removed" + (f" ({n_stop} characters)" if n_stop > 1 else ""), where)
    stops = [i + 1 for i, ch in enumerate(seq) if ch == "*"]
    if stops:
        shown = ", ".join(str(p) for p in stops[:5]) + (" ..." if len(stops) > 5 else "")
        result.error("FASTA_INTERNAL_STOP",
                     f"'*' inside the sequence at position(s) {shown}: an internal stop codon means the "
                     "translation is wrong; check the DNA reading frame instead of deleting the '*'", where)
    if not seq:
        result.error("FASTA_EMPTY_SEQUENCE", "header without sequence", where)
    elif len(seq) > MAX_SEQUENCE_LENGTH:
        result.error("FASTA_TOO_LONG", f"{len(seq)} residues; longer than any antibody construct "
                     f"(limit {MAX_SEQUENCE_LENGTH})", where)
    return Record(ident, desc, seq, lineno)


def format_record(record: Record, wrap: int = 0, rewrite_header: bool = True) -> str:
    """One record as text (header + sequence, LF terminated).

    With ``rewrite_header`` the header follows the ABD convention from
    :mod:`abd_prep.fasta.roles` (``>H``/``>L``) and the original header is
    kept after it as the description, so provenance is not lost.
    """
    if rewrite_header and record.role:
        head = fasta_header(record.role, record.header)
    else:
        head = ">" + record.header
    lines = [head]
    seq = record.sequence
    if wrap and wrap > 0:
        lines.extend(seq[i:i + wrap] for i in range(0, len(seq), wrap))
    else:
        lines.append(seq)
    return "\n".join(lines) + "\n"


def write_fasta(records: Sequence[Record], out: Union[str, TextIO], wrap: int = 0,
                rewrite_header: bool = True) -> None:
    """Write records as UTF-8 text with LF line endings."""
    text = "".join(format_record(r, wrap, rewrite_header) for r in records)
    if isinstance(out, str):
        os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok=True)
        with open(out, "w", encoding="utf-8", newline="\n") as fh:
            fh.write(text)
    else:
        out.write(text)
