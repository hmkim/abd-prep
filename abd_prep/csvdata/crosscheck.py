"""F22: do the CSV rows refer to the sequences in a seed FASTA / Fv PDB?

Links in-silico predictions to wet-lab results: each row matches the
reference when its ``id`` equals a FASTA record ID (or a PDB chain ID),
or when its ``sequence`` -- or ``heavy`` and ``light`` both -- equals a
reference sequence. Rows that match neither way are listed.

Finding codes:

``CSV_CROSSREF_MATCHED``    info     k/n rows matched (by id / by sequence)
``CSV_CROSSREF_UNMATCHED``  warning  rows that matched nothing (first few listed)
``CSV_CROSSREF_EMPTY``      warning  the reference has no usable IDs/sequences
"""
import os
from dataclasses import dataclass, field
from typing import Dict, List, Set

from ..common.result import InputError, ValidationResult
from ..fasta.fastafile import read_fasta
from ..pdb.validate import inspect_file as pdb_inspect_file

LIST_LIMIT = 10
PDB_EXTENSIONS = (".pdb", ".ent", ".pdb1")


@dataclass
class Reference:
    path: str
    ids: Set[str] = field(default_factory=set)
    sequences: Set[str] = field(default_factory=set)
    kind: str = "fasta"


def load_reference(path: str, encoding: str = "utf-8-sig") -> Reference:
    """IDs and sequences from a FASTA (records) or a PDB (chains)."""
    if not os.path.isfile(path):
        raise InputError(f"--against: file not found: {path}")
    ref = Reference(path=path)
    if path.lower().endswith(PDB_EXTENSIONS):
        ref.kind = "pdb"
        info = pdb_inspect_file(path)
        for ch, s in info["chains"].items():  # type: ignore[union-attr]
            ref.ids.add(ch)
            if s["sequence"]:
                ref.sequences.add(s["sequence"].upper())
        return ref
    scratch = ValidationResult()
    try:
        fasta = read_fasta(path, scratch, encoding=encoding)
    except UnicodeDecodeError as exc:
        raise InputError(f"--against {os.path.basename(path)}: not valid UTF-8 ({exc.reason} near byte {exc.start}); "
                         "run 'fasta clean' on the reference first")
    for rec in fasta.records:
        if rec.id:
            ref.ids.add(rec.id)          # a cleaned seed file has IDs H/L, so a CSV keyed "H"/"L" matches it
        if rec.sequence:
            ref.sequences.add(rec.sequence.upper())
    return ref


class CrossRef:
    """Row hook for :func:`abd_prep.csvdata.validate.process_csv`."""

    def __init__(self, ref: Reference) -> None:
        self.ref = ref
        self.rows = 0
        self.by_id = 0
        self.by_sequence = 0
        self.unmatched: List[str] = []

    def __call__(self, lineno: int, values: Dict[str, str]) -> None:
        self.rows += 1
        ident = values.get("id", "").strip()
        matched = False
        if ident and ident in self.ref.ids:
            self.by_id += 1
            matched = True
        seq = values.get("sequence", "").strip().upper()
        heavy = values.get("heavy", "").strip().upper()
        light = values.get("light", "").strip().upper()
        if (seq and seq in self.ref.sequences) or (heavy and light and heavy in self.ref.sequences
                                                   and light in self.ref.sequences):
            self.by_sequence += 1
            matched = True
        if not matched:
            self.unmatched.append(f"line {lineno}" + (f" id {ident!r}" if ident else ""))

    def report(self, result: ValidationResult) -> Dict[str, object]:
        name = os.path.basename(self.ref.path)
        if not self.ref.ids and not self.ref.sequences:
            result.warning("CSV_CROSSREF_EMPTY", f"{name} has no record IDs or sequences to match against", "against")
        matched = self.rows - len(self.unmatched)
        result.info("CSV_CROSSREF_MATCHED",
                    f"{matched}/{self.rows} row(s) match {name} ({self.by_id} by id, {self.by_sequence} by sequence)",
                    "against")
        if self.unmatched:
            shown = "; ".join(self.unmatched[:LIST_LIMIT]) + (" ..." if len(self.unmatched) > LIST_LIMIT else "")
            result.warning("CSV_CROSSREF_UNMATCHED",
                           f"{len(self.unmatched)} row(s) match neither an ID nor a sequence in {name}: {shown}",
                           "against")
        return {"reference": name, "kind": self.ref.kind, "rows": self.rows, "matched": matched,
                "by_id": self.by_id, "by_sequence": self.by_sequence, "unmatched": len(self.unmatched)}
