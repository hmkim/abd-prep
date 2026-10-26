"""CSV schema: standard columns, ``--map`` projection, missing tokens, numbers.

**This is the only module to edit when the real ABD CSV schema is
confirmed (plan R8).** The README says so too. Until then the standard
column set below is the tool's own definition, chosen so a lab result
table can be projected onto it with ``--map``:

``id``        clone / variant identifier (REQUIRED -- rows must be identifiable)
``sequence``  full antibody sequence, or
``heavy`` + ``light``  the two chains in separate columns (one of the two
              forms is RECOMMENDED; without any, F22 can still match by id)
``value``     measured number (KD, Tm, titre ...)
``unit``      unit of ``value`` (nM, pM, degC ...) -- mixing is reported, never converted
``target``    antigen name
``assay``     assay type (SPR, BLI, ELISA ...)
``notes``     free text

Values are never changed by this tool: no unit conversion, no number
reformatting, no case folding of sequences (plan §1 non-goals).
"""
import re
from typing import Dict, FrozenSet, Iterable, List, Optional, Tuple

from ..common.result import InputError

#: Standard column names in output order.
STANDARD_COLUMNS: Tuple[str, ...] = ("id", "sequence", "heavy", "light", "value", "unit", "target", "assay", "notes")
REQUIRED_COLUMNS: Tuple[str, ...] = ("id",)
#: Any one of these column groups satisfies the "has a sequence" recommendation.
SEQUENCE_GROUPS: Tuple[Tuple[str, ...], ...] = (("sequence",), ("heavy", "light"))
#: Columns whose values are protein sequences (F15 alphabet check applies).
SEQUENCE_COLUMNS: Tuple[str, ...] = ("sequence", "heavy", "light")
#: Columns whose values should parse as numbers.
NUMERIC_COLUMNS: Tuple[str, ...] = ("value",)

#: Tokens that mean "no value" (compared case-insensitively after stripping
#: blanks). Referenced by the validate rules and the README; defined once.
MISSING_TOKENS: FrozenSet[str] = frozenset({"", "na", "n/a", "nan", "none", "null", "-", "--", "nd", "n.d."})

_NUMBER_RE = re.compile(r"^[+-]?(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?$")
_DECIMAL_COMMA_RE = re.compile(r"^[+-]?\d+,\d+(?:[eE][+-]?\d+)?$")
_CENSORED_RE = re.compile(r"^[<>]=?\s*[+-]?(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?$")

KIND_NUMBER = "number"
KIND_MISSING = "missing"
KIND_DECIMAL_COMMA = "decimal_comma"
KIND_CENSORED = "censored"
KIND_TEXT = "text"


def is_missing(value: str) -> bool:
    return value.strip().lower() in MISSING_TOKENS


def classify_number(value: str) -> Tuple[str, Optional[float]]:
    """Classify a cell meant to hold a number.

    Returns ``(kind, number)`` where kind is one of ``number`` (parsed),
    ``missing`` (a missing token), ``decimal_comma`` (``1,2`` -- European
    decimal mark), ``censored`` (``<0.1``, ``>1000`` -- assay limit) or
    ``text``. Only ``number`` carries a float; the cell text itself is
    never rewritten.
    """
    text = value.strip()
    if is_missing(text):
        return KIND_MISSING, None
    if _NUMBER_RE.match(text):
        try:
            return KIND_NUMBER, float(text)
        except ValueError:  # pragma: no cover - regex already guarantees this
            return KIND_TEXT, None
    if _DECIMAL_COMMA_RE.match(text):
        return KIND_DECIMAL_COMMA, None
    if _CENSORED_RE.match(text):
        return KIND_CENSORED, None
    return KIND_TEXT, None


def parse_map(specs: Optional[Iterable[str]]) -> Dict[str, str]:
    """``["id=Clone,sequence=AA_seq", "value=KD"]`` -> ``{"id": "Clone", ...}``.

    Keys are standard column names; values are source header names
    (case-sensitive, matched after stripping blanks). Repeated ``--map``
    options are joined, so a header that itself contains a comma can be
    given in its own ``--map``.
    """
    out: Dict[str, str] = {}
    used: Dict[str, str] = {}
    for spec in specs or []:
        items = [spec] if spec.count("=") == 1 else spec.split(",")
        for item in items:
            item = item.strip()
            if not item:
                continue
            if "=" not in item:
                raise InputError(f"--map: expected STANDARD=SOURCE_HEADER, got {item!r}")
            std, src = (p.strip() for p in item.split("=", 1))
            std = std.lower()
            if std not in STANDARD_COLUMNS:
                raise InputError(f"--map: {std!r} is not a standard column; choose from {', '.join(STANDARD_COLUMNS)}")
            if not src:
                raise InputError(f"--map: empty source header for {std}")
            if std in out:
                raise InputError(f"--map: standard column {std!r} mapped twice ({out[std]!r} and {src!r})")
            if src in used:
                raise InputError(f"--map: source header {src!r} mapped to both {used[src]!r} and {std!r}")
            out[std] = src
            used[src] = std
    return out


def resolve_columns(header: List[str], mapping: Dict[str, str]) -> Tuple[Dict[str, int], List[int]]:
    """Return ``(standard_column -> index, unmapped indices)`` for a header.

    Explicit mappings win; a header that equals a standard name
    (case-insensitively) is mapped automatically unless that standard
    column was mapped elsewhere. A mapping that names a header the file
    does not have is an input error (exit 2) -- it is an option mistake.
    """
    stripped = [h.strip() for h in header]
    resolved: Dict[str, int] = {}
    taken: set = set()
    for std, src in mapping.items():
        if src not in stripped:
            raise InputError(f"--map {std}={src}: the file has no column {src!r}; header is: "
                             + ", ".join(repr(h) for h in stripped))
        idx = stripped.index(src)
        resolved[std] = idx
        taken.add(idx)
    lower = [h.lower() for h in stripped]
    for std in STANDARD_COLUMNS:
        if std in resolved:
            continue
        if std in lower:
            idx = lower.index(std)
            if idx not in taken:
                resolved[std] = idx
                taken.add(idx)
    unmapped = [i for i in range(len(header)) if i not in taken]
    return resolved, unmapped


def has_sequence_columns(resolved: Dict[str, int]) -> bool:
    return any(all(c in resolved for c in group) for group in SEQUENCE_GROUPS)
