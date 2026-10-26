import io
import os
import unittest

from abd_prep.pdb.pdbfile import Atom, format_atom_name, format_ter, read_pdb, write_pdb
from tests.helpers import atom_line, fixture, residue, ter_line

LINE = "ATOM      1  N   ASP L   1      22.146  36.117  40.236  1.00 14.85           N  "


class AtomTest(unittest.TestCase):
    def test_column_views(self):
        a = Atom(LINE)
        self.assertEqual(a.record, "ATOM")
        self.assertFalse(a.is_hetatm)
        self.assertEqual(a.serial, 1)
        self.assertEqual(a.name, "N")
        self.assertEqual(a.altloc, " ")
        self.assertEqual(a.resname, "ASP")
        self.assertEqual(a.chain, "L")
        self.assertEqual(a.resseq, 1)
        self.assertEqual(a.icode, " ")
        self.assertEqual(a.residue_key, ("L", 1, " "))
        self.assertAlmostEqual(a.xyz[0], 22.146)
        self.assertAlmostEqual(a.occupancy, 1.0)
        self.assertEqual(a.element, "N")

    def test_edits_touch_only_their_columns(self):
        a = Atom(LINE)
        a.set_chain("H")
        a.set_altloc(" ")
        self.assertEqual(a.line[21], "H")
        self.assertEqual(a.line[:21] + "L" + a.line[22:], LINE)
        a.set_record("HETATM")
        self.assertTrue(a.line.startswith("HETATM"))
        a.set_resname("MSE")
        self.assertEqual(a.resname, "MSE")
        a.set_name("SD")
        self.assertEqual(a.line[12:16], " SD ")
        a.set_element("S")
        self.assertEqual(a.line[76:78], " S")
        a.set_serial(99999)
        self.assertEqual(a.serial, 99999)
        self.assertEqual(len(a.line), len(LINE))

    def test_short_line_is_padded_on_edit(self):
        a = Atom(LINE[:54])            # no occupancy/element columns
        self.assertIsNone(a.occupancy)
        self.assertEqual(a.element, "")
        a.set_element("C")
        self.assertEqual(len(a.line), 80)
        self.assertEqual(a.element, "C")

    def test_format_atom_name(self):
        self.assertEqual(format_atom_name("CA"), " CA ")
        self.assertEqual(format_atom_name("N"), " N  ")
        self.assertEqual(format_atom_name("HD11"), "HD11")
        self.assertEqual(format_atom_name(" SD "), " SD ")


class ReaderTest(unittest.TestCase):
    def test_fixture_bookkeeping(self):
        pdb = read_pdb(fixture("1N8Z_mini.pdb"))
        self.assertEqual(pdb.n_models, 0)
        self.assertEqual(len(pdb.ters), 3)
        self.assertTrue(pdb.has_end)
        self.assertEqual(pdb.chains(), ["A", "B", "C"])
        self.assertEqual(pdb.record_counts["REMARK"], 3)
        self.assertEqual(len(pdb.first_model_atoms()), len(pdb.atoms))

    def test_models(self):
        pdb = read_pdb(fixture("6B9Z_cso_mini.pdb"))
        self.assertEqual(pdb.n_models, 2)
        first = pdb.first_model_atoms()
        self.assertEqual(len(first) * 2, len(pdb.atoms))
        self.assertTrue(all(a.model == 1 for a in first))

    def test_ter_position(self):
        lines = residue(1, "ALA", "A", 1) + [ter_line(5, "ALA", "A", 1)] + residue(6, "GLY", "B", 1) + ["END"]
        pdb = read_pdb(io.StringIO("\n".join(lines) + "\n"))
        self.assertEqual([t.after_index for t in pdb.ters], [4])
        self.assertEqual(pdb.ters[0].chain, "A")
        self.assertEqual(pdb.ters[0].resseq, 1)

    def test_bad_residue_number_is_input_error(self):
        from abd_prep.common.result import InputError
        bad = LINE[:22] + "A0B1" + LINE[26:]
        pdb = read_pdb(io.StringIO(bad + "\n"))
        with self.assertRaises(InputError):
            _ = pdb.atoms[0].resseq


class WriterTest(unittest.TestCase):
    def test_format_ter_columns(self):
        t = format_ter(1646, "CYS", "A", 214)
        self.assertEqual(len(t), 80)
        self.assertEqual(t[:27], "TER    1646      CYS A 214 ")
        self.assertEqual(format_ter(5, "LYS", "H", 100, "A")[22:27], " 100A")

    def test_ter_and_end_regenerated(self):
        atoms = [Atom(l) for l in residue(10, "ALA", "L", 1) + residue(14, "GLY", "L", 2) + residue(100, "SER", "H", 1)]
        buf = io.StringIO()
        ters = write_pdb(atoms, buf)
        out = buf.getvalue().splitlines()
        self.assertEqual(len(ters), 2)
        self.assertEqual(ters[0][:27], "TER      18      GLY L   2 ")   # last L serial 17 + 1
        self.assertEqual(ters[1][:27], "TER     104      SER H   1 ")
        self.assertEqual(out[-1].rstrip(), "END")
        self.assertEqual(sum(1 for l in out if l.startswith("TER")), 2)
        self.assertEqual(sum(1 for l in out if l.startswith("ATOM")), 12)
        # atom text untouched
        self.assertEqual(out[0], atoms[0].line)

    def test_renumber(self):
        atoms = [Atom(l) for l in residue(500, "ALA", "A", 7) + residue(900, "GLY", "B", 8)]
        buf = io.StringIO()
        ters = write_pdb(atoms, buf, renumber=True)
        serials = [int(l[6:11]) for l in buf.getvalue().splitlines() if l[:3] in ("ATO", "TER")]
        self.assertEqual(serials, list(range(1, 11)))
        self.assertEqual(ters[0][6:11], "    5")
        # the caller's atoms are not mutated
        self.assertEqual(atoms[0].serial, 500)

    def test_write_to_path(self):
        import tempfile
        atoms = [Atom(l) for l in residue(1, "ALA", "A", 1)]
        path = os.path.join(tempfile.mkdtemp(), "x.pdb")
        write_pdb(atoms, path)
        with open(path) as fh:
            text = fh.read()
        self.assertTrue(text.endswith("END".ljust(80) + "\n"))
        self.assertIn("TER       5      ALA A   1", text)


if __name__ == "__main__":
    unittest.main()
