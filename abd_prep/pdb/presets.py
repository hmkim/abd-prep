"""Presets: named target + antibody selections (F6, plan §3.4).

A preset is a list of *parts* -- ``target``, ``heavy``, ``light`` -- each
naming the RCSB entry, the chain and residue range to cut, and the chain
ID to write. Target and antibody may come from different entries (the
6LBX + 6B9Z preset); heavy and light must share one entry because they
are written into one file and must sit in the same coordinate frame.

Expectations (``expect_residues``, ``expected_gaps``) were measured on the
RCSB files current in 2026-10 and feed ``validate`` so a preset run is
checked against known-good numbers. Residue numbering: 1N8Z chain C uses
mature-protein numbering (UniProt P04626 minus the 22-residue signal
peptide); the Fab chains A/B are sequential with no insertion codes.
"""
from dataclasses import dataclass, field
from typing import Dict, FrozenSet, List, Optional, Tuple

from ..common.result import InputError
from .ranges import ChainSelection, RangeSyntaxError

ROLES = ("target", "heavy", "light")
FV_ROLES = ("heavy", "light")


@dataclass(frozen=True)
class Part:
    """One chain cut from one entry, written under ``out_chain``."""

    role: str                 # target | heavy | light
    pdb_id: str               # RCSB ID, or "input" for a local file
    selection: str            # e.g. "C:489-607" (one chain, any number of ranges)
    out_chain: Optional[str] = None   # None keeps the source chain ID

    @property
    def chain(self) -> str:
        return self.selection.split(":", 1)[0].strip()

    @property
    def written_chain(self) -> str:
        return self.out_chain or self.chain

    def check(self) -> None:
        if self.role not in ROLES:
            raise InputError(f"unknown part role {self.role!r}")
        try:
            sel = ChainSelection.parse([self.selection])
        except RangeSyntaxError as exc:
            raise InputError(f"{self.role}: {exc}")
        if len(sel.chains) != 1:
            raise InputError(f"{self.role}: selection {self.selection!r} must name exactly one chain")
        if self.out_chain is not None and len(self.out_chain) != 1:
            raise InputError(f"{self.role}: output chain ID must be one character, got {self.out_chain!r}")


@dataclass(frozen=True)
class Preset:
    name: str
    description: str
    parts: Tuple[Part, ...]
    #: output chain -> exact residue count measured on the current RCSB file
    expect_residues: Dict[str, int] = field(default_factory=dict)
    #: (output chain, residue before, residue after) of known unmodelled loops
    expected_gaps: FrozenSet[Tuple[str, int, int]] = frozenset()
    notes: Tuple[str, ...] = ()

    def part(self, role: str) -> Optional[Part]:
        for p in self.parts:
            if p.role == role:
                return p
        return None

    @property
    def pdb_ids(self) -> List[str]:
        seen: Dict[str, None] = {}
        for p in self.parts:
            seen.setdefault(p.pdb_id)
        return list(seen)

    def to_dict(self) -> Dict[str, object]:
        return {
            "name": self.name,
            "description": self.description,
            "parts": [{"role": p.role, "pdb_id": p.pdb_id, "selection": p.selection,
                       "out_chain": p.written_chain} for p in self.parts],
            "expect_residues": dict(self.expect_residues),
            "expected_gaps": sorted([list(g) for g in self.expected_gaps]),
            "notes": list(self.notes),
        }


PRESETS: Dict[str, Preset] = {
    "her2-trastuzumab": Preset(
        name="her2-trastuzumab",
        description="HER2 extracellular domain IV (target) + trastuzumab Fv, both from 1N8Z "
                    "(HER2 ECD / Herceptin Fab complex, 2.52 A). Matches the blog walkthrough.",
        parts=(
            Part("target", "1N8Z", "C:489-607"),          # mature numbering; chain ID kept (R1 open)
            Part("light", "1N8Z", "A:1-107", "L"),
            Part("heavy", "1N8Z", "B:1-120", "H"),
        ),
        expect_residues={"C": 109, "L": 107, "H": 120},
        expected_gaps=frozenset({("C", 580, 591)}),        # residues 581-590 unmodelled in 1N8Z
        notes=(
            "1N8Z chain C is numbered as the mature protein (UniProt P04626 numbering minus 22).",
            "Fab chains A (light, 1-214) and B (heavy, 1-220) are numbered sequentially, no insertion codes; "
            "Fv ends 107 (VL) / 120 (VH) are literature values for trastuzumab.",
            "The 10-residue loop 581-590 is not modelled in 1N8Z; validate reports it as an expected gap.",
            "Target chain keeps its original ID C until the ABD chain-label requirement is confirmed; "
            "pass --rename-target to change it.",
        ),
    ),
    "her2-trastuzumab-6lbx-6b9z": Preset(
        name="her2-trastuzumab-6lbx-6b9z",
        description="HER2 domain IV fragment from 6LBX (chain B 530-601, repebody partner removed) + "
                    "trastuzumab-derived Fv from 6B9Z. The pair actually produced with the prototype scripts.",
        parts=(
            Part("target", "6LBX", "B:530-601"),
            Part("light", "6B9Z", "A:1-107", "L"),
            Part("heavy", "6B9Z", "B:1-120", "H"),
        ),
        expect_residues={"B": 72, "L": 107, "H": 120},
        expected_gaps=frozenset(),
        notes=(
            "6B9Z antibody differs from the 1N8Z trastuzumab at several framework positions "
            "(e.g. VL SPILLSASV vs SPSSLSASV); treat it as a trastuzumab-derived variant, not the "
            "reference sequence (plan R5).",
            "6B9Z residue B175 is CSO/CYS heterogeneity outside the Fv range; the Fv cut is unaffected.",
            "Target and antibody come from different crystals, so their coordinate frames are unrelated.",
        ),
    ),
}

DEFAULT_PRESET = "her2-trastuzumab"


def get_preset(name: str) -> Preset:
    try:
        return PRESETS[name]
    except KeyError:
        raise InputError(f"unknown preset {name!r}; available: {', '.join(sorted(PRESETS))}")


def output_groups(parts: Tuple[Part, ...]) -> Dict[str, List[Part]]:
    """Group parts into output files: ``target`` -> one file, ``fv`` -> heavy+light."""
    groups: Dict[str, List[Part]] = {}
    for p in parts:
        key = "target" if p.role == "target" else "fv"
        groups.setdefault(key, []).append(p)
    return groups


def output_filename(name: str, group: str) -> str:
    return f"{name}_{group}.pdb"


def validate_parts(parts: Tuple[Part, ...]) -> None:
    """Shared sanity rules for presets and ``split`` arguments."""
    if not parts:
        raise InputError("nothing to split: give --target and/or --heavy/--light (or --preset)")
    roles = [p.role for p in parts]
    for r in roles:
        if roles.count(r) > 1:
            raise InputError(f"part {r!r} given twice")
    for p in parts:
        p.check()
    fv = [p for p in parts if p.role in FV_ROLES]
    if len({p.pdb_id for p in fv}) > 1:
        raise InputError("heavy and light chains must come from the same entry (they share one file "
                         "and one coordinate frame)")
    written = [p.written_chain for p in fv]
    if len(written) != len(set(written)):
        raise InputError(f"heavy and light would both be written as chain {written[0]!r}")
    for p in fv:
        for q in fv:
            if p is not q and p.pdb_id == q.pdb_id and p.chain == q.chain:
                raise InputError(f"heavy and light select the same source chain {p.chain!r}")


for _preset in PRESETS.values():
    validate_parts(_preset.parts)
