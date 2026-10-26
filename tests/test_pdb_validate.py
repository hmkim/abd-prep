import io
import unittest

from abd_prep.pdb.clean import CleanOptions, clean_atoms
from abd_prep.pdb.pdbfile import Atom, read_pdb, write_pdb
from abd_prep.pdb.ranges import ChainSelection, parse_rename
from abd_prep.pdb.validate import (ValidateOptions, format_inspect, inspect_atoms, inspect_file, summarize_chains,
                                   to_fasta, validate_atoms, validate_file)
from tests.helpers import atom_line, fixture, residue, ter_line


def parse(lines):
    return read_pdb(io.StringIO("\n".join(lines) + "\n"))


def clean_text(fixture_name, keep=(), rename=""):
    """Clean a fixture through the library and re-read the written text."""
    pdb = read_pdb(fixture(fixture_name))
    out = clean_atoms(pdb, CleanOptions(selection=ChainSelection.parse(list(keep)), rename=parse_rename(rename)))
    buf = io.StringIO()
    write_pdb(out.atoms, buf)
    return read_pdb(io.StringIO(buf.getvalue()))


def good_two_chain():
    """Two chains, proper TERs, END; CA-CA 3.8 A along x."""
    lines = []
    serial = 1
    for i in range(3):
        lines += residue(serial, "ALA", "L", i + 1, x=3.8 * i)
        serial += 4
    lines.append(ter_line(serial, "ALA", "L", 3))
    serial += 1
    for i in range(3):
        lines += residue(serial, "GLY", "H", i + 1, x=3.8 * i)
        serial += 4
    lines.append(ter_line(serial, "GLY", "H", 3))
    lines.append("END".ljust(80))
    return lines


class RuleTest(unittest.TestCase):
    def test_clean_output_passes(self):
        pdb = clean_text("1N8Z_mini.pdb", keep=["A", "B"], rename="A=L,B=H")
        res, details = validate_atoms(pdb)
        self.assertFalse(res.has_errors, res.format_table())
        # L/H have only 5 residues in the fixture -> residue count warnings only
        self.assertEqual(set(res.codes()), {"PDB_RESIDUE_COUNT"})
        self.assertEqual(list(details["chains"]), ["L", "H"])

    def test_altloc_is_error(self):
        res, _ = validate_atoms(read_pdb(fixture("1N8Z_mini.pdb")))
        self.assertIn("PDB_ALTLOC_PRESENT", res.codes())
        self.assertTrue(res.has_errors)

    def test_hetatm_warning_and_allow(self):
        pdb = read_pdb(fixture("1N8Z_mini.pdb"))
        res, _ = validate_atoms(pdb)
        f = next(x for x in res.findings if x.code == "PDB_HETATM_PRESENT")
        self.assertEqual(f.level, "warning")
        res, _ = validate_atoms(pdb, ValidateOptions(allow_het=True))
        f = next(x for x in res.findings if x.code == "PDB_HETATM_PRESENT")
        self.assertEqual(f.level, "info")

    def test_expect_chains(self):
        pdb = parse(good_two_chain())
        self.assertNotIn("PDB_CHAIN_COUNT", validate_atoms(pdb, ValidateOptions(expect_chains=2))[0].codes())
        self.assertIn("PDB_CHAIN_COUNT", validate_atoms(pdb, ValidateOptions(expect_chains=1))[0].codes())

    def test_prototype_style_ter_leftovers(self):
        """TERs copied from the source: wrong residue, extra TERs of removed chains (P1)."""
        lines = good_two_chain()
        # replace the two correct TERs with prototype-style leftovers
        lines = [l for l in lines if not l.startswith("TER")]
        idx = next(i for i, l in enumerate(lines) if l[21] == "H")
        lines.insert(idx, ter_line(900, "CYS", "A", 214))                 # mismatch at L boundary
        lines.insert(len(lines) - 1, ter_line(901, "LYS", "B", 221))      # mismatch at end
        lines.insert(len(lines) - 1, ter_line(902, "GLY", "E", 81))       # extra
        res, _ = validate_atoms(parse(lines))
        codes = res.codes()
        self.assertEqual(codes.count("PDB_TER_MISMATCH"), 3)
        self.assertIn("PDB_TER_EXTRA", codes)

    def test_ter_missing_and_misplaced(self):
        lines = [l for l in good_two_chain() if not l.startswith("TER")]
        res, _ = validate_atoms(parse(lines))
        self.assertEqual(res.codes().count("PDB_TER_MISSING"), 2)
        lines = good_two_chain()
        lines.insert(4, ter_line(50, "ALA", "L", 1))     # inside chain L after residue 1
        res, _ = validate_atoms(parse(lines))
        self.assertIn("PDB_TER_MISPLACED", res.codes())

    def test_numbering_gap_info_and_ca_distance_warning(self):
        lines = residue(1, "ALA", "A", 1, x=0) + residue(5, "ALA", "A", 2, x=3.8) + residue(9, "ALA", "A", 4, x=7.6) \
            + residue(13, "ALA", "A", 5, x=20.0) + [ter_line(17, "ALA", "A", 5), "END"]
        res, _ = validate_atoms(parse(lines))
        gaps = [(f.code, f.where, f.level) for f in res.findings if f.code.startswith("PDB_GAP")]
        self.assertIn(("PDB_GAP_NUMBERING", "A:2-4", "info"), gaps)
        self.assertNotIn("PDB_GAP_CA_DISTANCE", [g[0] for g in gaps if g[1] == "A:2-4"])   # 3.8 A: numbering only
        self.assertIn(("PDB_GAP_CA_DISTANCE", "A:4-5", "warning"), gaps)
        # expected gap -> info
        res, _ = validate_atoms(parse(lines), ValidateOptions(expected_gaps={("A", 4, 5)}))
        f = next(x for x in res.findings if x.code == "PDB_GAP_CA_DISTANCE")
        self.assertEqual(f.level, "info")

    def test_kabat_insertions_are_not_gaps(self):
        res, _ = validate_atoms(read_pdb(fixture("kabat_mini.pdb")))
        self.assertFalse(any(c.startswith("PDB_GAP") for c in res.codes()), res.format_table())

    def test_end_missing(self):
        lines = good_two_chain()[:-1]
        self.assertIn("PDB_END_MISSING", validate_atoms(parse(lines))[0].codes())

    def test_multi_model_warning(self):
        res, details = validate_atoms(read_pdb(fixture("6B9Z_cso_mini.pdb")))
        self.assertIn("PDB_MULTI_MODEL", res.codes())
        self.assertEqual(details["models"], 2)

    def test_nonstandard_resname_in_atom_records(self):
        lines = residue(1, "ALA", "A", 1) + residue(5, "UNK", "A", 2, x=3.8) + [ter_line(9, "UNK", "A", 2), "END"]
        self.assertIn("PDB_NONSTANDARD_RESNAME", validate_atoms(parse(lines))[0].codes())

    def test_residue_count_for_hl(self):
        pdb = parse(good_two_chain())
        res, _ = validate_atoms(pdb)
        self.assertEqual(res.codes().count("PDB_RESIDUE_COUNT"), 2)
        res, _ = validate_atoms(pdb, ValidateOptions(expect_residues={"L": (1, 5), "H": (1, 5)}))
        self.assertNotIn("PDB_RESIDUE_COUNT", res.codes())

    def test_duplicate_serial_and_blank_chain(self):
        lines = residue(1, "ALA", " ", 1) + residue(1, "GLY", " ", 2, x=3.8) + [ter_line(9, "GLY", " ", 2), "END"]
        codes = validate_atoms(parse(lines))[0].codes()
        self.assertIn("PDB_SERIAL_DUPLICATE", codes)
        self.assertIn("PDB_CHAIN_BLANK", codes)

    def test_hydrogens_info(self):
        lines = residue(1, "ALA", "A", 1) + [atom_line(5, "H", "ALA", "A", 1, element="H"), ter_line(6, "ALA", "A", 1), "END"]
        self.assertIn("PDB_HYDROGENS_PRESENT", validate_atoms(parse(lines))[0].codes())

    def test_residue_duplicate(self):
        lines = residue(1, "ALA", "A", 1) + residue(5, "GLY", "A", 2, x=3.8) + residue(9, "SER", "A", 1, x=7.6) \
            + [ter_line(13, "SER", "A", 1), "END"]
        self.assertIn("PDB_RESIDUE_DUPLICATE", validate_atoms(parse(lines))[0].codes())

    def test_empty_file(self):
        res, _ = validate_atoms(parse(["END"]))
        self.assertEqual(res.codes(), ["PDB_NO_ATOMS"])


class InspectTest(unittest.TestCase):
    def test_summary(self):
        info = inspect_file(fixture("1N8Z_mini.pdb"))
        self.assertEqual(list(info["chains"]), ["A", "B", "C"])
        self.assertEqual(info["chains"]["A"]["sequence"], "DIQMT")
        self.assertEqual(info["chains"]["C"]["first"], "489")
        self.assertEqual(info["altloc_atoms"], 4)
        self.assertEqual(info["hetatm"], {"HOH": 2})
        text = format_inspect(info)
        self.assertIn("chain A: 1-5 (5 res", text)
        self.assertIn("DIQMT", text)

    def test_fasta_headers_follow_roles_table(self):
        pdb = clean_text("1N8Z_mini.pdb", keep=["A", "B"], rename="A=L,B=H")
        fasta = to_fasta(inspect_atoms(pdb))
        lines = fasta.splitlines()
        self.assertTrue(lines[0].startswith(">L "))
        self.assertEqual(lines[1], "DIQMT")
        self.assertTrue(lines[2].startswith(">H "))
        self.assertEqual(lines[3], "EVQLV")

    def test_fasta_wrap(self):
        pdb = read_pdb(fixture("1N8Z_mini.pdb"))
        lines = to_fasta(inspect_atoms(pdb), width=2).splitlines()
        self.assertEqual(lines[1:4], ["DI", "QM", "T"])


if __name__ == "__main__":
    unittest.main()
