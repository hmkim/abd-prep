"""Chain / residue range selection (F2; fixes prototype problem P2).

Grammar of one ``--keep`` / ``--target`` / ``--heavy`` / ``--light`` value::

    CHAIN                      whole chain            A
    CHAIN:LO-HI                one range              B:1-120
    CHAIN:LO-HI,LO-HI          several ranges         C:489-607,610-612
    CHAIN:RES                  single residue         H:100
    CHAIN:LO-HIi               with insertion code    H:1-113A

Residues are keyed by ``(resSeq, iCode)`` exactly as they appear in PDB
columns 23-26 / 27, so Kabat-numbered antibodies (82A/82B/82C, 100A...)
select correctly.

Insertion-code semantics at the bounds (documented in README):

* a bound written **without** an insertion code includes every inserted
  residue of that number -- ``1-82`` keeps 82, 82A, 82B, 82C;
* a bound written **with** one is exact -- ``1-82A`` keeps 82 and 82A only.

PDB residue numbers may be negative (``-3``); the grammar handles it.
"""
import re
from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Optional, Tuple

ResidueKey = Tuple[int, str]  # (resSeq, iCode) with iCode " " when absent

#: Sorts after every real insertion code (printable ASCII), used for the
#: "include all insertions of this residue number" upper bound.
_ICODE_MAX = "\x7f"

_RES_RE = re.compile(r"^(-?\d+)([A-Za-z]?)$")
_RANGE_RE = re.compile(r"^(-?\d+[A-Za-z]?)-(-?\d+[A-Za-z]?)$")


class RangeSyntaxError(ValueError):
    """Raised for an unparsable range specification."""


def parse_residue(text: str) -> Tuple[int, str]:
    """``"113A"`` -> ``(113, "A")``; ``"5"`` -> ``(5, "")`` (no code given)."""
    m = _RES_RE.match(text)
    if not m:
        raise RangeSyntaxError(f"bad residue {text!r} (expected e.g. 100 or 100A)")
    return int(m.group(1)), m.group(2).upper()


@dataclass(frozen=True)
class ResidueRange:
    """Inclusive range over ``(resSeq, iCode)`` keys."""

    lo: ResidueKey
    hi: ResidueKey

    @classmethod
    def parse(cls, text: str) -> "ResidueRange":
        text = text.strip()
        m = _RANGE_RE.match(text)
        if m:
            lo_n, lo_i = parse_residue(m.group(1))
            hi_n, hi_i = parse_residue(m.group(2))
        else:
            lo_n, lo_i = parse_residue(text)   # single residue
            hi_n, hi_i = lo_n, lo_i
        lo = (lo_n, lo_i or " ")
        hi = (hi_n, hi_i or _ICODE_MAX)         # no code given -> all insertions
        if lo > hi:
            raise RangeSyntaxError(f"range {text!r}: start is after end")
        return cls(lo, hi)

    def contains(self, key: ResidueKey) -> bool:
        return self.lo <= key <= self.hi

    def __str__(self) -> str:
        lo = f"{self.lo[0]}{self.lo[1].strip()}"
        hi_code = "" if self.hi[1] == _ICODE_MAX else self.hi[1].strip()
        hi = f"{self.hi[0]}{hi_code}"
        return lo if (self.lo[0] == self.hi[0] and self.lo[1] == " " and self.hi[1] == _ICODE_MAX) else f"{lo}-{hi}"


@dataclass
class ChainSelection:
    """Which chains/residues to keep. Empty ``chains`` means keep everything.

    ``chains[chain_id]`` is a list of ranges; an empty list keeps the whole
    chain. Chain IDs are compared as written (case-sensitive, PDB col 22).
    """

    chains: Dict[str, List[ResidueRange]] = field(default_factory=dict)

    @classmethod
    def parse(cls, specs: Iterable[str]) -> "ChainSelection":
        sel = cls()
        for spec in specs:
            sel.add(spec)
        return sel

    def add(self, spec: str) -> None:
        spec = spec.strip()
        if not spec:
            raise RangeSyntaxError("empty selection")
        if ":" in spec:
            chain, rest = spec.split(":", 1)
            chain = chain.strip()
            if len(chain) != 1:
                raise RangeSyntaxError(f"chain id must be one character in {spec!r} (PDB format limit)")
            ranges = [ResidueRange.parse(part) for part in rest.split(",") if part.strip()]
            if not ranges:
                raise RangeSyntaxError(f"no residue range after ':' in {spec!r}")
            self.chains.setdefault(chain, []).extend(ranges)
        else:
            if len(spec) != 1:
                raise RangeSyntaxError(f"chain id must be one character: {spec!r}")
            self.chains.setdefault(spec, [])   # whole chain (empty range list)

    @property
    def is_empty(self) -> bool:
        return not self.chains

    def wants_chain(self, chain: str) -> bool:
        return self.is_empty or chain in self.chains

    def wants(self, chain: str, resseq: int, icode: str = " ") -> bool:
        if self.is_empty:
            return True
        ranges = self.chains.get(chain)
        if ranges is None:
            return False
        if not ranges:
            return True
        key = (resseq, icode or " ")
        return any(r.contains(key) for r in ranges)

    def describe(self) -> Dict[str, str]:
        return {ch: ",".join(str(r) for r in rs) if rs else "*" for ch, rs in self.chains.items()}


def parse_rename(text: Optional[str]) -> Dict[str, str]:
    """``"A=L,B=H"`` -> ``{"A": "L", "B": "H"}`` (prototype ``--rename`` syntax)."""
    out: Dict[str, str] = {}
    if not text:
        return out
    for part in text.split(","):
        part = part.strip()
        if not part:
            continue
        if "=" not in part:
            raise RangeSyntaxError(f"bad rename {part!r} (expected OLD=NEW)")
        old, new = (p.strip() for p in part.split("=", 1))
        if len(old) != 1 or len(new) != 1:
            raise RangeSyntaxError(f"chain ids must be one character in rename {part!r}")
        if old in out:
            raise RangeSyntaxError(f"chain {old!r} renamed twice")
        out[old] = new
    return out
