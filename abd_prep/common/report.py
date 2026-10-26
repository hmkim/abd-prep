"""``report.json`` -- one schema for every group and action (plan §3.5).

::

    {
      "tool":      {"name": "abd-prep", "version": "...", "python": "3.x.y"},
      "command":   "abd-prep pdb clean in.pdb out.pdb --keep A:1-107",
      "generated_at": "2026-10-23T06:30:00Z",
      "inputs":    [{"kind": "pdb", "path": "in.pdb",  "sha256": "...", "bytes": 123}],
      "outputs":   [{"kind": "pdb", "path": "out.pdb", "sha256": "...", "bytes": 120}],
      "findings":  [{"level": "warning", "code": "...", "where": "...", "message": "..."}],
      "details":   {...group specific...}
    }

Paths are recorded as basenames by default (``--report-paths basename``)
because reports travel in mail and tickets and absolute paths leak the
user name and home directory. ``--report-paths full`` keeps them whole.
"""
import json
import os
import platform
import shlex
from dataclasses import dataclass, field
from datetime import datetime, timezone
from typing import Any, Dict, List, Optional, Sequence

from .. import __version__
from .hashing import file_size, sha256_file
from .result import ValidationResult

TOOL_NAME = "abd-prep"
PATH_MODES = ("basename", "full")


def _render_path(path: str, mode: str) -> str:
    if mode == "full":
        return os.path.abspath(path)
    return os.path.basename(path)


def file_entry(kind: str, path: str, path_mode: str = "basename") -> Dict[str, Any]:
    """Describe one input/output file (kind, path, sha256, bytes)."""
    return {
        "kind": kind,
        "path": _render_path(path, path_mode),
        "sha256": sha256_file(path),
        "bytes": file_size(path),
    }


def render_command(argv: Sequence[str], path_mode: str = "basename") -> str:
    """Re-join argv for the report; path-like arguments follow ``path_mode``."""
    parts: List[str] = [TOOL_NAME]
    for a in argv:
        if path_mode == "basename" and (os.sep in a or a.startswith("~")) and not a.startswith("-"):
            a = os.path.basename(a)
        parts.append(shlex.quote(a))
    return " ".join(parts)


@dataclass
class Report:
    command: str = ""
    inputs: List[Dict[str, Any]] = field(default_factory=list)
    outputs: List[Dict[str, Any]] = field(default_factory=list)
    result: ValidationResult = field(default_factory=ValidationResult)
    details: Dict[str, Any] = field(default_factory=dict)
    path_mode: str = "basename"

    def add_input(self, kind: str, path: str) -> None:
        self.inputs.append(file_entry(kind, path, self.path_mode))

    def add_output(self, kind: str, path: str) -> None:
        self.outputs.append(file_entry(kind, path, self.path_mode))

    def to_dict(self) -> Dict[str, Any]:
        return {
            "tool": {"name": TOOL_NAME, "version": __version__, "python": platform.python_version()},
            "command": self.command,
            "generated_at": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
            "inputs": self.inputs,
            "outputs": self.outputs,
            "findings": self.result.to_list(),
            "details": self.details,
        }

    def to_json(self, indent: Optional[int] = 2) -> str:
        return json.dumps(self.to_dict(), indent=indent, ensure_ascii=False)

    def write(self, path: str) -> None:
        os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
        with open(path, "w", encoding="utf-8") as fh:
            fh.write(self.to_json())
            fh.write("\n")
