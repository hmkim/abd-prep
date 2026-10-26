"""Tests pinned to the M1 code-review findings (plan §5.1, M2 review)."""
import io
import os
import tempfile
import unittest

from abd_prep.common.hashing import check_utf8
from abd_prep.pdb.clean import CleanOptions, clean_atoms, clean_file
from abd_prep.pdb.pdbfile import read_pdb, ter_serial_collisions, write_pdb
from abd_prep.pdb.ranges import ChainSelection
from abd_prep.pdb.validate import validate_atoms
from tests.helpers import atom_line, fixture, residue, ter_line


def parse(lines):
    return read_pdb(io.StringIO("\n".join(lines) + "\n"))


class Utf8OffsetTest(unittest.TestCase):
    def test_offset_points_at_the_bad_byte(self):
        fd, path = tempfile.mkstemp()
        with os.fdopen(fd, "wb") as fh:
            fh.write(b"abc\xff\xfedef")
        ok, detail = check_utf8(path)
        self.assertFalse(ok)
        self.assertIn("near byte 3", detail)

    def test_unknown_encoding(self):
        ok, detail = check_utf8(fixture("1N8Z_mini.pdb"), encoding="no-such-codec")
        self.assertFalse(ok)
        self.assertIn("unknown encoding", detail)


class NonAsciiRoundTripTest(unittest.TestCase):
    def test_byte_preserved_and_warned(self):
        lines = residue(1, "ALA", "A", 1)
        lines[0] = lines[0][:72] + "\xe9" + lines[0][73:]
        fd, src = tempfile.mkstemp(suffix=".pdb")
        with os.fdopen(fd, "wb") as fh:
            fh.write(("\n".join(lines) + "\nEND\n").encode("latin-1"))
        out = os.path.join(tempfile.mkdtemp(), "o.pdb")
        outcome = clean_file(src, out, CleanOptions())
        self.assertIn("PDB_NON_ASCII", outcome.result.codes())
        with open(out, "rb") as fh:
            data = fh.read()
        self.assertIn(b"\xe9", data)                      # byte survives, no crash
        res, _ = validate_atoms(read_pdb(out))
        self.assertIn("PDB_NON_ASCII", res.codes())


class HydrogenParityTest(unittest.TestCase):
    def test_validate_and_clean_agree_on_name_only_hydrogens(self):
        lines = residue(1, "ALA", "A", 1) + [atom_line(5, "HB1", "ALA", "A", 1, element="")]
        pdb = parse(lines + [ter_line(6, "ALA", "A", 1), "END"])
        res, _ = validate_atoms(pdb)
        h = [f for f in res.findings if f.code == "PDB_HYDROGENS_PRESENT"]
        self.assertEqual(len(h), 1)
        self.assertIn("1 hydrogen", h[0].message)
        out = clean_atoms(pdb, CleanOptions(strip_h=True))
        self.assertEqual(out.details["hydrogens_removed"], 1)


class TerSerialCollisionTest(unittest.TestCase):
    def test_collision_detected_and_renumber_fixes_it(self):
        lines = residue(1, "ALA", "A", 1) + residue(5, "GLY", "B", 1)      # no TER between chains
        pdb = parse(lines + ["END"])
        buf = io.StringIO()
        ters = write_pdb(pdb.atoms, buf)
        self.assertEqual(ter_serial_collisions(pdb.atoms, ters), [5])
        buf = io.StringIO()
        write_pdb(pdb.atoms, buf, renumber=True)
        written = read_pdb(io.StringIO(buf.getvalue()))
        res, _ = validate_atoms(written)
        self.assertNotIn("PDB_SERIAL_DUPLICATE", res.codes())

    def test_clean_file_warns(self):
        lines = residue(1, "ALA", "A", 1) + residue(5, "GLY", "B", 1) + ["END"]
        fd, src = tempfile.mkstemp(suffix=".pdb")
        with os.fdopen(fd, "w") as fh:
            fh.write("\n".join(lines) + "\n")
        outcome = clean_file(src, os.path.join(tempfile.mkdtemp(), "o.pdb"), CleanOptions())
        self.assertIn("PDB_SERIAL_DUPLICATE", outcome.result.codes())
        outcome = clean_file(src, os.path.join(tempfile.mkdtemp(), "o.pdb"), CleanOptions(renumber=True))
        self.assertNotIn("PDB_SERIAL_DUPLICATE", outcome.result.codes())


class NumberingDecreaseTest(unittest.TestCase):
    def test_decrease_is_info(self):
        lines = residue(1, "ALA", "A", 100) + residue(5, "GLY", "A", 95, x=3.8) + [ter_line(9, "GLY", "A", 95), "END"]
        res, _ = validate_atoms(parse(lines))
        gaps = [f for f in res.findings if f.code == "PDB_GAP_NUMBERING"]
        self.assertEqual(len(gaps), 1)
        self.assertEqual(gaps[0].level, "info")
        self.assertIn("decreases 100 -> 95", gaps[0].message)


class SelectionEmptyMessageTest(unittest.TestCase):
    def test_hetatm_only_selection_mentions_keep_het(self):
        pdb = read_pdb(fixture("1N8Z_mini.pdb"))                    # waters are HOH A 215-216
        out = clean_atoms(pdb, CleanOptions(selection=ChainSelection.parse(["A:215-216"])))
        msg = next(f.message for f in out.result.findings if f.code == "PDB_SELECTION_EMPTY")
        self.assertIn("--keep-het", msg)
        out = clean_atoms(pdb, CleanOptions(selection=ChainSelection.parse(["A:900-950"])))
        msg = next(f.message for f in out.result.findings if f.code == "PDB_SELECTION_EMPTY")
        self.assertNotIn("--keep-het", msg)


if __name__ == "__main__":
    unittest.main()
