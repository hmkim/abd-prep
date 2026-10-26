"""Shared helpers: fixture paths and a minimal ATOM-line builder."""
import os
import tempfile
from typing import List, Optional

FIXTURES = os.path.join(os.path.dirname(os.path.abspath(__file__)), "fixtures")
CACHE = os.path.expanduser("~/.cache/abd-prep")


def fixture(name: str) -> str:
    return os.path.join(FIXTURES, name)


def atom_line(serial: int, name: str, resname: str, chain: str, resseq: int, x: float = 0.0, y: float = 0.0,
              z: float = 0.0, altloc: str = " ", icode: str = " ", occ: float = 1.0, element: Optional[str] = None,
              hetatm: bool = False) -> str:
    """Build one 80-column ATOM/HETATM record."""
    rec = "HETATM" if hetatm else "ATOM  "
    field = name if len(name) >= 4 else (" " + name).ljust(4)
    el = (element if element is not None else name.strip()[0]).rjust(2)
    return (f"{rec}{serial:5d} {field}{altloc}{resname:>3s} {chain}{resseq:4d}{icode}   "
            f"{x:8.3f}{y:8.3f}{z:8.3f}{occ:6.2f}{0.0:6.2f}          {el}  ")


def residue(serial: int, resname: str, chain: str, resseq: int, x: float = 0.0, icode: str = " ",
            names=("N", "CA", "C", "O"), hetatm: bool = False, altloc: str = " ", occ: float = 1.0) -> List[str]:
    """A residue with backbone atoms spaced so CA-CA to the next residue is 3.8 A when x steps by 3.8."""
    return [atom_line(serial + i, n, resname, chain, resseq, x + i * 0.1, 0.0, 0.0, altloc=altloc, icode=icode,
                      occ=occ, hetatm=hetatm) for i, n in enumerate(names)]


def ter_line(serial: int, resname: str, chain: str, resseq: int, icode: str = " ") -> str:
    return f"TER   {serial:5d}      {resname:>3s} {chain}{resseq:4d}{icode}".ljust(80)


def write_temp(lines: List[str], suffix: str = ".pdb") -> str:
    fd, path = tempfile.mkstemp(suffix=suffix)
    with os.fdopen(fd, "w", newline="\n") as fh:
        fh.write("\n".join(lines) + "\n")
    return path
