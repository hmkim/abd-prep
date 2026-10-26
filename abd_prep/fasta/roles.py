"""Chain roles and the ABD FASTA header convention (F14; R7).

Two things live here and nowhere else:

1. **The ABD header convention.** The exact header ABD expects for seed
   antibodies is not documented beyond "proper chain identification
   (H/L)". Until it is confirmed on the console (plan R7) the convention
   is a bare role letter: ``>H`` heavy, ``>L`` light, and ``>scFv`` for a
   single-chain record. **Change only** :data:`ABD_HEADER_ID` when R7 is
   settled -- ``pdb inspect --fasta``, ``fasta clean`` and the tests all
   read it.

2. **Role inference from header tokens.** Deliberately conservative:
   a header is split on non-alphanumeric characters and each token is
   looked up in :data:`HEAVY_TOKENS` / :data:`LIGHT_TOKENS` (plus the
   ``IGHV*`` / ``IGKV*`` / ``IGLV*`` gene-name prefixes). A header with no
   known token is **not guessed** -- the user names the records with
   ``--heavy-id`` / ``--light-id``. Uploading heavy and light swapped is
   the most expensive mistake this tool can prevent.

Roles are ``"H"``, ``"L"`` and, for an scFv record that carries both
variable domains, ``"HL"``. There is no target/antigen role: ABD takes
the antigen as a PDB structure, not a FASTA record (plan §8.1/§8.2).
"""
import re
from typing import TYPE_CHECKING, Dict, FrozenSet, Iterable, List, Optional, Sequence, Tuple

from ..common.result import InputError, ValidationResult

if TYPE_CHECKING:  # pragma: no cover - import for type checkers only
    from .fastafile import Record

#: Role -> header ID written after ``>``. Keys are the roles the tool
#: knows; the PDB group maps chain IDs ``H``/``L`` straight onto them.
ABD_HEADER_ID: Dict[str, str] = {
    "H": "H",
    "L": "L",
    "HL": "scFv",
}

ROLE_HEAVY = "H"
ROLE_LIGHT = "L"
ROLE_SCFV = "HL"

FRAMEWORKS = ("igg", "nanobody", "scfv")

#: Header tokens (lower-case) that mark a heavy chain / VH / VHH.
HEAVY_TOKENS: FrozenSet[str] = frozenset({
    "h", "hc", "vh", "vhh", "heavy", "heavychain", "hchain", "igh", "nanobody", "nb", "sdab",
})
#: Header tokens (lower-case) that mark a light chain / VL.
LIGHT_TOKENS: FrozenSet[str] = frozenset({
    "l", "lc", "vl", "vk", "light", "lightchain", "lchain", "kappa", "lambda", "igk", "igl",
})
#: Gene-name prefixes (IMGT nomenclature), e.g. IGHV3-23, IGKV1-39, IGLV2-14.
HEAVY_PREFIXES: Tuple[str, ...] = ("ighv", "ighg", "ighm", "igha", "ighe", "ighd", "ighj")
LIGHT_PREFIXES: Tuple[str, ...] = ("igkv", "iglv", "igkc", "iglc", "igkj", "iglj")
#: Header tokens that mark a single-chain Fv (both domains in one record).
SCFV_TOKENS: FrozenSet[str] = frozenset({"scfv", "scfvs", "singlechain"})

_TOKEN_RE = re.compile(r"[^A-Za-z0-9]+")


def fasta_header(role_or_chain: str, description: str = "") -> str:
    """Build a FASTA header line (without the newline).

    ``role_or_chain`` is ``H``/``L``/``HL`` (mapped through
    :data:`ABD_HEADER_ID`) or any other chain ID, which is written as-is.
    ``description`` is appended after a space; standard FASTA parsers
    treat only the first token as the record ID, so it is safe to include.
    """
    ident = ABD_HEADER_ID.get(role_or_chain, role_or_chain)
    return f">{ident} {description}".rstrip()


def header_tokens(header: str) -> List[str]:
    return [t.lower() for t in _TOKEN_RE.split(header) if t]


def infer_role(header: str) -> Tuple[Optional[str], List[str]]:
    """Infer ``H``/``L``/``HL`` from a header; return ``(role, matched_tokens)``.

    An scFv token wins outright (an scFv header legitimately names both
    VH and VL). Otherwise ``role`` is ``None`` when no known token is
    present, and ``"?"`` when both heavy and light tokens are present.
    """
    hits_h: List[str] = []
    hits_l: List[str] = []
    hits_s: List[str] = []
    for tok in header_tokens(header):
        if tok in SCFV_TOKENS:
            hits_s.append(tok)
        elif tok in HEAVY_TOKENS or tok.startswith(HEAVY_PREFIXES):
            hits_h.append(tok)
        elif tok in LIGHT_TOKENS or tok.startswith(LIGHT_PREFIXES):
            hits_l.append(tok)
    if hits_s:
        return ROLE_SCFV, hits_s
    if hits_h and hits_l:
        return "?", hits_h + hits_l
    if hits_h:
        return ROLE_HEAVY, hits_h
    if hits_l:
        return ROLE_LIGHT, hits_l
    return None, []


def _find(records: Sequence["Record"], ident: str, option: str) -> "Record":
    hits = [r for r in records if r.id == ident]
    if not hits:
        ids = ", ".join(repr(r.id) for r in records[:10]) or "(none)"
        raise InputError(f"{option} {ident!r}: no record has that ID; records are: {ids}")
    if len(hits) > 1:
        raise InputError(f"{option} {ident!r}: {len(hits)} records share that ID; make the IDs unique first")
    return hits[0]


def assign_roles(records: Sequence["Record"], result: ValidationResult, framework: Optional[str] = None,
                 heavy_id: Optional[str] = None, light_id: Optional[str] = None) -> None:
    """Set ``record.role`` on every record (F14).

    Precedence: explicit ``--heavy-id``/``--light-id`` > single-record
    frameworks (``nanobody`` -> H, ``scfv`` -> HL) > header tokens. A
    record that ends up without a role is an error under ``--framework
    igg`` and a warning otherwise (the framework rules are skipped then).
    A light-chain token on a nanobody record is a conflict error.
    """
    if framework is not None and framework not in FRAMEWORKS:
        raise InputError(f"unknown framework {framework!r}; choose one of {', '.join(FRAMEWORKS)}")
    if heavy_id is not None and light_id is not None and heavy_id == light_id:
        raise InputError("--heavy-id and --light-id name the same record")

    for r in records:
        r.role, r.role_source = None, ""
    if heavy_id is not None:
        rec = _find(records, heavy_id, "--heavy-id")
        rec.role, rec.role_source = ROLE_HEAVY, "option"
    if light_id is not None:
        rec = _find(records, light_id, "--light-id")
        rec.role, rec.role_source = ROLE_LIGHT, "option"

    for r in records:
        inferred, tokens = infer_role(r.header)
        where = r.id or f"line {r.lineno}"
        if r.role is not None:
            if inferred not in (None, "?", r.role):
                result.info("FASTA_ROLE_OVERRIDDEN",
                            f"header suggests {inferred} ({', '.join(tokens)}) but the option assigns {r.role}", where)
            continue
        if framework == "scfv":
            r.role, r.role_source = ROLE_SCFV, "framework"
            continue
        if framework == "nanobody":
            if inferred in (ROLE_LIGHT, ROLE_SCFV):
                what = "a light chain" if inferred == ROLE_LIGHT else "an scFv"
                result.error("FASTA_ROLE_CONFLICT",
                             f"header names {what} ({', '.join(tokens)}) but a nanobody (VHH) is a single heavy-chain "
                             "domain; check the framework or the record", where)
            r.role, r.role_source = ROLE_HEAVY, "framework"
            continue
        if inferred == "?":
            result.error("FASTA_ROLE_AMBIGUOUS",
                         f"header has both heavy and light tokens ({', '.join(tokens)}); use --heavy-id/--light-id", where)
        elif inferred is None:
            if framework == "igg":
                result.error("FASTA_ROLE_UNKNOWN",
                             "cannot tell heavy from light by the header; name the records with "
                             f"--heavy-id {r.id or 'ID'} or --light-id {r.id or 'ID'} (known header tokens: "
                             "heavy/VH/VHH/IGHV*, light/VL/kappa/lambda/IGKV*/IGLV*)", where)
            else:
                result.warning("FASTA_ROLE_UNKNOWN",
                               "chain role (H/L) not recognized from the header; pass --heavy-id/--light-id "
                               "or use recognizable headers before uploading", where)
        else:
            r.role, r.role_source = inferred, "header"


def check_framework(records: Sequence["Record"], framework: Optional[str], result: ValidationResult) -> None:
    """F16 record-count / role consistency against ``--framework``."""
    if framework is None:
        result.info("FASTA_FRAMEWORK_UNSPECIFIED",
                    "--framework not given; record-count/role consistency rules skipped", "file")
        return
    n = len(records)
    roles = [r.role for r in records]
    if framework == "igg":
        n_h, n_l = roles.count(ROLE_HEAVY), roles.count(ROLE_LIGHT)
        n_s = roles.count(ROLE_SCFV)
        if not (n == 2 and n_h == 1 and n_l == 1):
            other = f", {n_s} scFv" if n_s else ""
            result.error("FASTA_FRAMEWORK_MISMATCH",
                         f"igg expects exactly one heavy and one light record; found {n} record(s) "
                         f"({n_h} heavy, {n_l} light{other}, {n - n_h - n_l - n_s} without a role)", "file")
    elif framework in ("nanobody", "scfv"):
        if n != 1:
            result.error("FASTA_FRAMEWORK_MISMATCH",
                         f"{framework} expects exactly one record; found {n}", "file")


def roles_summary(records: Iterable["Record"]) -> Dict[str, int]:
    out: Dict[str, int] = {}
    for r in records:
        key = r.role or "unassigned"
        out[key] = out.get(key, 0) + 1
    return out
