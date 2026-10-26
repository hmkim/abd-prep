"""Legacy PDB fixed-column records: parse, edit in place, write (F10, F11).

The whole tool works on the raw 80-column line. An :class:`Atom` keeps
its original text and edits only the columns it has to (altloc, chain,
record name, residue/atom name, element), so untouched coordinates and
B-factors are written back byte-for-byte -- this is what makes the
regression comparison against the prototype's outputs exact.

Column reference (1-based, PDB format v3.3)::

     1- 6  record name      ATOM  / HETATM
     7-11  serial
    13-16  atom name
    17     altLoc
    18-20  resName
    22     chainID
    23-26  resSeq
    27     iCode
    31-38 / 39-46 / 47-54   x / y / z
    55-60  occupancy
    61-66  tempFactor
    77-78  element
    79-80  charge

Files are read and written as ``latin-1``: every byte maps to exactly one
character, so a stray non-ASCII byte in a record survives the round trip
unchanged instead of crashing the writer (clean reports it as
``PDB_NON_ASCII``). Well-formed PDB files are pure ASCII, for which
latin-1 is identical.
"""
from collections import Counter
from dataclasses import dataclass, field
from typing import Dict, Iterable, Iterator, List, Optional, Sequence, TextIO, Tuple, Union

from ..common.residues import is_hydrogen_name
from ..common.result import InputError

LINE_WIDTH = 80
COORD_RECORDS = ("ATOM", "HETATM")
FILE_ENCODING = "latin-1"


def _pad(line: str) -> str:
    return line if len(line) >= LINE_WIDTH else line.ljust(LINE_WIDTH)


def format_atom_name(name: str) -> str:
    """Place a bare atom name into the 4-character field (columns 13-16).

    Names of 1-3 characters start in column 14 (``" CA "``, ``" SD "``);
    4-character names fill the field (``"HD11"``). Element symbols of two
    letters (``FE``) would start in column 13, but standard residues
    contain none, so the 3-or-fewer rule is enough here.
    """
    name = name.strip()
    if len(name) >= 4:
        return name[:4]
    return (" " + name).ljust(4)


@dataclass
class Atom:
    """One ATOM/HETATM record, edited through the raw line."""

    line: str                     # without trailing newline
    model: int = 0                # 0 = outside any MODEL block
    lineno: int = 0               # 1-based position in the source file

    # -- read-only views ------------------------------------------------
    @property
    def record(self) -> str:
        return self.line[0:6].strip()

    @property
    def is_hetatm(self) -> bool:
        return self.line.startswith("HETATM")

    @property
    def serial(self) -> Optional[int]:
        return _int_or_none(self.line[6:11])

    @property
    def name(self) -> str:
        return self.line[12:16].strip()

    @property
    def name_field(self) -> str:
        return self.line[12:16]

    @property
    def altloc(self) -> str:
        return self.line[16] if len(self.line) > 16 else " "

    @property
    def resname(self) -> str:
        return self.line[17:20].strip()

    @property
    def chain(self) -> str:
        return self.line[21] if len(self.line) > 21 else " "

    @property
    def resseq(self) -> int:
        try:
            return int(self.line[22:26])
        except ValueError:
            raise InputError(f"line {self.lineno}: residue number {self.line[22:26]!r} is not an integer "
                             "(hybrid-36 / very large structures are outside the PDB-format scope)")

    @property
    def icode(self) -> str:
        c = self.line[26] if len(self.line) > 26 else " "
        return c if c.strip() else " "

    @property
    def residue_key(self) -> Tuple[str, int, str]:
        return (self.chain, self.resseq, self.icode)

    @property
    def xyz(self) -> Tuple[float, float, float]:
        try:
            return float(self.line[30:38]), float(self.line[38:46]), float(self.line[46:54])
        except ValueError:
            raise InputError(f"line {self.lineno}: coordinates are not numeric")

    @property
    def occupancy(self) -> Optional[float]:
        try:
            return float(self.line[54:60])
        except (ValueError, IndexError):
            return None

    @property
    def element(self) -> str:
        return self.line[76:78].strip() if len(self.line) >= 78 else ""

    @property
    def is_hydrogen(self) -> bool:
        """Hydrogen/deuterium by the element column, falling back to the name.

        Shared by ``clean --strip-h`` and ``validate`` so both agree (F12).
        """
        el = self.element
        if el:
            return el in ("H", "D")
        return is_hydrogen_name(self.name, self.resname)

    @property
    def is_ascii(self) -> bool:
        return all(ord(ch) < 128 for ch in self.line)

    # -- edits (in place) -----------------------------------------------
    def _splice(self, start: int, end: int, text: str) -> None:
        line = _pad(self.line) if end > len(self.line) else self.line
        self.line = line[:start] + text + line[end:]

    def set_altloc(self, code: str = " ") -> None:
        self._splice(16, 17, code[:1] or " ")

    def set_chain(self, chain: str) -> None:
        self._splice(21, 22, chain[:1])

    def set_record(self, record: str) -> None:
        self._splice(0, 6, record.ljust(6)[:6])

    def set_resname(self, resname: str) -> None:
        self._splice(17, 20, resname.rjust(3)[:3])

    def set_name(self, name: str) -> None:
        self._splice(12, 16, format_atom_name(name))

    def set_element(self, symbol: str) -> None:
        self._splice(76, 78, symbol.rjust(2)[:2])

    def set_serial(self, serial: int) -> None:
        self._splice(6, 11, _serial_text(serial))

    def copy(self) -> "Atom":
        return Atom(self.line, self.model, self.lineno)


def _int_or_none(text: str) -> Optional[int]:
    try:
        return int(text)
    except ValueError:
        return None


def _serial_text(serial: int) -> str:
    return f"{serial:5d}" if serial <= 99999 else "*****"


@dataclass
class TerRecord:
    """A TER line and where it sat: after ``after_index`` atoms."""

    line: str
    after_index: int              # number of ATOM/HETATM records seen before this TER
    lineno: int = 0

    @property
    def serial(self) -> Optional[int]:
        return _int_or_none(self.line[6:11])

    @property
    def resname(self) -> str:
        return self.line[17:20].strip()

    @property
    def chain(self) -> str:
        return self.line[21] if len(self.line) > 21 else " "

    @property
    def resseq(self) -> Optional[int]:
        return _int_or_none(self.line[22:26])

    @property
    def icode(self) -> str:
        c = self.line[26] if len(self.line) > 26 else " "
        return c if c.strip() else " "


@dataclass
class PdbFile:
    """Everything the tool needs from one PDB file."""

    atoms: List[Atom] = field(default_factory=list)   # all models, file order
    ters: List[TerRecord] = field(default_factory=list)
    n_models: int = 0
    has_end: bool = False
    record_counts: Counter = field(default_factory=Counter)
    path: str = ""

    def first_model_atoms(self) -> List[Atom]:
        """Atoms of the first MODEL (or every atom when there is no MODEL)."""
        if self.n_models <= 1:
            return list(self.atoms)
        return [a for a in self.atoms if a.model <= 1]

    def chains(self, atoms: Optional[Sequence[Atom]] = None) -> List[str]:
        seen: Dict[str, None] = {}
        for a in (atoms if atoms is not None else self.atoms):
            seen.setdefault(a.chain)
        return list(seen)


def read_pdb(source: Union[str, TextIO, Iterable[str]], path: str = "") -> PdbFile:
    """Parse a PDB file from a path, an open text file or an iterable of lines."""
    if isinstance(source, str):
        with open(source, "r", encoding=FILE_ENCODING) as fh:
            return _read_lines(fh, path or source)
    return _read_lines(source, path)


def _read_lines(lines: Iterable[str], path: str) -> PdbFile:
    pdb = PdbFile(path=path)
    model = 0
    for lineno, raw in enumerate(lines, start=1):
        line = raw.rstrip("\r\n")
        rec = line[0:6].strip()
        if not rec:
            continue
        pdb.record_counts[rec] += 1
        if rec in COORD_RECORDS:
            pdb.atoms.append(Atom(line, model, lineno))
        elif rec == "MODEL":
            model += 1
            pdb.n_models = model
        elif rec == "ENDMDL":
            continue
        elif rec == "TER":
            pdb.ters.append(TerRecord(line, len(pdb.atoms), lineno))
        elif rec == "END":
            pdb.has_end = True
    return pdb


def format_ter(serial: int, resname: str, chain: str, resseq: int, icode: str = " ") -> str:
    """Build a TER record for the residue that ends a chain."""
    return _pad(f"TER   {_serial_text(serial)}      {resname.rjust(3)[:3]} {chain[:1]}{resseq:4d}{(icode or ' ')[:1]}")


END_LINE = _pad("END")


def write_pdb(atoms: Sequence[Atom], out: Union[str, TextIO], renumber: bool = False) -> List[str]:
    """Write atoms with freshly generated TER records and a final END (F10).

    A TER is emitted every time the (output) chain ID changes and after the
    last atom. Its serial is ``last atom serial + 1`` -- the PDB convention
    -- or the next sequential number when ``renumber`` is set, in which
    case atom serials are rewritten 1..N as well.

    Returns the TER lines written (useful for tests and reports). Use
    :func:`ter_serial_collisions` to find TER serials that duplicate an
    atom serial, which happens when the source had no TER at a boundary.
    """
    if isinstance(out, str):
        with open(out, "w", encoding=FILE_ENCODING, newline="\n") as fh:
            return write_pdb(atoms, fh, renumber)

    ters: List[str] = []
    counter = 0
    prev: Optional[Atom] = None
    last_serial = 0
    for atom in atoms:
        if prev is not None and atom.chain != prev.chain:
            counter += 1
            ter_serial = counter if renumber else last_serial + 1
            ter = format_ter(ter_serial, prev.resname, prev.chain, prev.resseq, prev.icode)
            ters.append(ter)
            out.write(ter + "\n")
        if renumber:
            counter += 1
            atom = atom.copy()
            atom.set_serial(counter)
        out.write(atom.line + "\n")
        last_serial = atom.serial if atom.serial is not None else last_serial + 1
        prev = atom
    if prev is not None:
        counter += 1
        ter_serial = counter if renumber else last_serial + 1
        ter = format_ter(ter_serial, prev.resname, prev.chain, prev.resseq, prev.icode)
        ters.append(ter)
        out.write(ter + "\n")
    out.write(END_LINE + "\n")
    return ters


def ter_serial_collisions(atoms: Sequence[Atom], ters: Sequence[str]) -> List[int]:
    """TER serials (from :func:`write_pdb`) that duplicate an atom serial.

    With ``renumber=False`` a TER takes ``last atom serial + 1``; if the
    source file had no TER at that boundary, the next chain's first atom
    already owns that number. ``--renumber`` resolves it.
    """
    atom_serials = {a.serial for a in atoms if a.serial is not None}
    hits: List[int] = []
    for t in ters:
        s = _int_or_none(t[6:11])
        if s is not None and s in atom_serials:
            hits.append(s)
    return hits


def iter_residues(atoms: Iterable[Atom]) -> Iterator[Tuple[Tuple[str, int, str], str, List[Atom]]]:
    """Group consecutive atoms into residues: ``(key, resname, atoms)``."""
    current_key = None
    current_name = ""
    bucket: List[Atom] = []
    for atom in atoms:
        key = atom.residue_key
        if key != current_key:
            if bucket:
                yield current_key, current_name, bucket  # type: ignore[misc]
            current_key, current_name, bucket = key, atom.resname, []
        bucket.append(atom)
    if bucket:
        yield current_key, current_name, bucket  # type: ignore[misc]
