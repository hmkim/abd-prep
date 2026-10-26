"""Findings and validation results shared by every ``validate``/``clean``.

Severity levels and exit-code contract (plan §3.1):

* ``info``    -- recorded, never affects the exit code
* ``warning`` -- exit 0 unless ``--strict`` is given (then 1)
* ``error``   -- exit 1

Exit code 2 (input/usage error) is raised by the CLI layer as
:class:`InputError`, not stored as a finding.
"""
import json
from dataclasses import asdict, dataclass, field
from typing import Dict, Iterable, List, Optional

LEVELS = ("info", "warning", "error")
EXIT_OK = 0
EXIT_VALIDATION = 1
EXIT_INPUT = 2


class InputError(Exception):
    """Bad input file or options; the CLI maps this to exit code 2."""


@dataclass
class Finding:
    """One validation observation.

    ``code`` is a stable upper-case identifier (one per rule in plan §3.3,
    e.g. ``PDB_ALTLOC_PRESENT``) so tests, README and report.json agree.
    ``where`` names the location (chain, residue, record id, row...).
    """

    level: str
    code: str
    message: str
    where: str = ""

    def __post_init__(self) -> None:
        if self.level not in LEVELS:
            raise ValueError(f"unknown finding level {self.level!r}")


@dataclass
class ValidationResult:
    """Ordered collection of findings plus a few counters."""

    findings: List[Finding] = field(default_factory=list)

    # -- recording -----------------------------------------------------
    def add(self, level: str, code: str, message: str, where: str = "") -> Finding:
        f = Finding(level, code, message, where)
        self.findings.append(f)
        return f

    def info(self, code: str, message: str, where: str = "") -> Finding:
        return self.add("info", code, message, where)

    def warning(self, code: str, message: str, where: str = "") -> Finding:
        return self.add("warning", code, message, where)

    def error(self, code: str, message: str, where: str = "") -> Finding:
        return self.add("error", code, message, where)

    def extend(self, other: "ValidationResult") -> None:
        self.findings.extend(other.findings)

    # -- querying ------------------------------------------------------
    def by_level(self, level: str) -> List[Finding]:
        return [f for f in self.findings if f.level == level]

    @property
    def errors(self) -> List[Finding]:
        return self.by_level("error")

    @property
    def warnings(self) -> List[Finding]:
        return self.by_level("warning")

    @property
    def has_errors(self) -> bool:
        return any(f.level == "error" for f in self.findings)

    @property
    def has_warnings(self) -> bool:
        return any(f.level == "warning" for f in self.findings)

    def codes(self) -> List[str]:
        return [f.code for f in self.findings]

    def exit_code(self, strict: bool = False) -> int:
        if self.has_errors or (strict and self.has_warnings):
            return EXIT_VALIDATION
        return EXIT_OK

    # -- rendering -----------------------------------------------------
    def to_list(self) -> List[Dict[str, str]]:
        return [asdict(f) for f in self.findings]

    def to_json(self, indent: Optional[int] = 2) -> str:
        return json.dumps(self.to_list(), indent=indent, ensure_ascii=False)

    def format_table(self, title: str = "") -> str:
        """Plain-text table: LEVEL  CODE  WHERE  MESSAGE, plus a summary line."""
        lines: List[str] = []
        if title:
            lines.append(title)
        if not self.findings:
            lines.append("  (no findings)")
        else:
            w_code = max(len(f.code) for f in self.findings)
            w_where = max(len(f.where) for f in self.findings)
            for f in self.findings:
                lines.append(f"  {f.level.upper():7s} {f.code:{w_code}s}  {f.where:{w_where}s}  {f.message}")
        n_err, n_warn, n_info = len(self.errors), len(self.warnings), len(self.by_level("info"))
        lines.append(f"  summary: {n_err} error(s), {n_warn} warning(s), {n_info} info")
        return "\n".join(lines)


def merge(results: Iterable[ValidationResult]) -> ValidationResult:
    out = ValidationResult()
    for r in results:
        out.extend(r)
    return out
