"""Residue and amino-acid tables shared by the pdb, fasta and csv groups.

Everything here is data: the 20 standard amino acids, their canonical
heavy-atom names (used when converting modified residues back to the
standard form), the table of modified residues the ``--std-residues``
option knows how to convert, and the protein alphabet check used by
FASTA and CSV sequence columns (F15).
"""
from typing import Dict, FrozenSet, Iterable, List, NamedTuple, Optional, Tuple

#: Three-letter -> one-letter code for the 20 standard amino acids.
AA3_TO_1: Dict[str, str] = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V",
}

AA1_TO_3: Dict[str, str] = {one: three for three, one in AA3_TO_1.items()}

#: The 20 standard residue names as they appear in PDB column 18-20.
STANDARD_AA3: FrozenSet[str] = frozenset(AA3_TO_1)

#: The 20 standard one-letter codes (FASTA / CSV sequence alphabet).
STANDARD_AA1: FrozenSet[str] = frozenset(AA3_TO_1.values())

#: One-letter codes that are legal IUPAC symbols but ambiguous or
#: non-standard for ABD (rejected by F15 unless explicitly replaced).
AMBIGUOUS_AA1: FrozenSet[str] = frozenset("XBZJUO")

#: Nucleotide letters; a record made only of these is probably DNA/RNA.
NUCLEIC_LETTERS: FrozenSet[str] = frozenset("ACGTUN")

_BACKBONE = ("N", "CA", "C", "O", "OXT")

#: Canonical heavy-atom names per standard residue (PDB v3 nomenclature).
#: Used to strip atoms that have no counterpart after a modified-residue
#: conversion (e.g. the extra OD of CSO once it becomes CYS).
STANDARD_ATOMS: Dict[str, FrozenSet[str]] = {
    "ALA": frozenset(_BACKBONE + ("CB",)),
    "ARG": frozenset(_BACKBONE + ("CB", "CG", "CD", "NE", "CZ", "NH1", "NH2")),
    "ASN": frozenset(_BACKBONE + ("CB", "CG", "OD1", "ND2")),
    "ASP": frozenset(_BACKBONE + ("CB", "CG", "OD1", "OD2")),
    "CYS": frozenset(_BACKBONE + ("CB", "SG")),
    "GLN": frozenset(_BACKBONE + ("CB", "CG", "CD", "OE1", "NE2")),
    "GLU": frozenset(_BACKBONE + ("CB", "CG", "CD", "OE1", "OE2")),
    "GLY": frozenset(_BACKBONE),
    "HIS": frozenset(_BACKBONE + ("CB", "CG", "ND1", "CD2", "CE1", "NE2")),
    "ILE": frozenset(_BACKBONE + ("CB", "CG1", "CG2", "CD1")),
    "LEU": frozenset(_BACKBONE + ("CB", "CG", "CD1", "CD2")),
    "LYS": frozenset(_BACKBONE + ("CB", "CG", "CD", "CE", "NZ")),
    "MET": frozenset(_BACKBONE + ("CB", "CG", "SD", "CE")),
    "PHE": frozenset(_BACKBONE + ("CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ")),
    "PRO": frozenset(_BACKBONE + ("CB", "CG", "CD")),
    "SER": frozenset(_BACKBONE + ("CB", "OG")),
    "THR": frozenset(_BACKBONE + ("CB", "OG1", "CG2")),
    "TRP": frozenset(_BACKBONE + ("CB", "CG", "CD1", "CD2", "NE1", "CE2", "CE3", "CZ2", "CZ3", "CH2")),
    "TYR": frozenset(_BACKBONE + ("CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ", "OH")),
    "VAL": frozenset(_BACKBONE + ("CB", "CG1", "CG2")),
}


class ModifiedResidue(NamedTuple):
    """How to turn one modified residue back into its standard parent."""

    parent: str                       # standard three-letter code
    rename: Dict[str, str]            # atom name -> new atom name
    note: str                         # shown as an info finding


#: Modified residues ``--std-residues`` can convert (F4 / P3).
#: After the renames listed here, any atom whose name is not in
#: ``STANDARD_ATOMS[parent]`` is dropped and reported.
MODIFIED_RESIDUES: Dict[str, ModifiedResidue] = {
    "MSE": ModifiedResidue("MET", {"SE": "SD"}, "selenomethionine -> methionine (SE renamed SD)"),
    "CSO": ModifiedResidue("CYS", {}, "S-hydroxycysteine -> cysteine (OD dropped)"),
    "CSD": ModifiedResidue("CYS", {}, "3-sulfinoalanine -> cysteine (OD1/OD2 dropped)"),
    "OCS": ModifiedResidue("CYS", {}, "cysteine sulfonic acid -> cysteine (OD1-3 dropped)"),
    "CME": ModifiedResidue("CYS", {}, "S,S-(2-hydroxyethyl)thiocysteine -> cysteine (adduct dropped)"),
    "SEP": ModifiedResidue("SER", {}, "phosphoserine -> serine (phosphate dropped)"),
    "TPO": ModifiedResidue("THR", {}, "phosphothreonine -> threonine (phosphate dropped)"),
    "PTR": ModifiedResidue("TYR", {}, "phosphotyrosine -> tyrosine (phosphate dropped)"),
    "HYP": ModifiedResidue("PRO", {}, "hydroxyproline -> proline (OD1 dropped)"),
    "MLY": ModifiedResidue("LYS", {}, "N-dimethyl-lysine -> lysine (methyls dropped)"),
    "KCX": ModifiedResidue("LYS", {}, "lysine NZ-carboxylic acid -> lysine (carboxyl dropped)"),
    "PCA": ModifiedResidue("GLN", {"OE": "OE1"},
                           "pyroglutamic acid -> glutamine (OE renamed OE1; side chain incomplete, NE2 absent)"),
}

#: Elements the standard residues contain; used to rewrite PDB columns
#: 77-78 after an atom rename (first letter of the canonical atom name).
PROTEIN_ELEMENTS: FrozenSet[str] = frozenset("CNOS")


def one_letter(resname: str, unknown: str = "X") -> str:
    """Return the one-letter code for a three-letter residue name."""
    return AA3_TO_1.get(resname.upper(), unknown)


def check_protein_alphabet(seq: str) -> Tuple[List[Tuple[int, str]], bool]:
    """Validate a one-letter protein sequence (F15).

    Returns ``(bad_positions, looks_nucleic)`` where ``bad_positions`` is a
    list of ``(1-based index, character)`` for every character outside
    the 20 standard letters, and ``looks_nucleic`` is True when the whole
    sequence is made of nucleotide letters only (A/C/G/T/U/N) -- a strong
    hint that a DNA/RNA sequence was supplied instead of a protein.
    """
    bad = [(i + 1, ch) for i, ch in enumerate(seq) if ch not in STANDARD_AA1]
    letters = set(seq)
    looks_nucleic = bool(letters) and letters <= NUCLEIC_LETTERS
    return bad, looks_nucleic


def element_for_atom(atom_name: str) -> Optional[str]:
    """Best-effort element symbol for a protein heavy atom name.

    Standard residues contain only C/N/O/S, so the first alphabetic
    character of the canonical name is the element. Returns None when
    the name does not start with one of those letters.
    """
    for ch in atom_name:
        if ch.isalpha():
            return ch if ch in PROTEIN_ELEMENTS else None
    return None


#: Two-letter element symbols that start with ``H`` and can appear as the
#: whole atom name of a HETATM (metals): never hydrogens.
_HEAVY_H_NAMES: FrozenSet[str] = frozenset({"HG", "HO", "HF"})


def is_hydrogen_name(atom_name: str, resname: str = "") -> bool:
    """Name-based hydrogen test used when PDB columns 77-78 are blank.

    Leading digits are skipped (``1HB``, ``2HG1``); the first letter must
    be ``H``. In standard amino acids every such name is a hydrogen
    (``HE`` in ARG, ``HG`` in SER/CYS, ``HH`` in TYR ...). For a
    non-standard residue an exact two-letter name ``HG``/``HO``/``HF`` is a
    metal (mercury, holmium, hafnium) and is kept.
    """
    name = atom_name.strip()
    i = 0
    while i < len(name) and name[i].isdigit():
        i += 1
    core = name[i:]
    if not core or core[0] != "H":
        return False
    if resname.upper() not in STANDARD_AA3 and core in _HEAVY_H_NAMES:
        return False
    return True


def standard_atoms(resname: str) -> Optional[FrozenSet[str]]:
    """Canonical heavy-atom set for a standard residue, or None."""
    return STANDARD_ATOMS.get(resname.upper())


def sequence_from_residues(resnames: Iterable[str]) -> str:
    """Join three-letter names into a one-letter string (unknown -> X)."""
    return "".join(one_letter(r) for r in resnames)
