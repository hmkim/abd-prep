import io
import os
import tempfile
import unittest

from abd_prep.pdb.clean import CleanOptions, clean_atoms, clean_file, polymer_like_hetero_residues, select_altlocs
from abd_prep.pdb.pdbfile import read_pdb
from abd_prep.pdb.ranges import ChainSelection, parse_rename
from tests.helpers import atom_line, fixture, residue, ter_line


def parse(lines):
    return read_pdb(io.StringIO("\n".join(lines) + "\n"))


def opts(keep=(), rename="", **kw):
    return CleanOptions(selection=ChainSelection.parse(list(keep)), rename=parse_rename(rename), **kw)


def residues_of(atoms):
    seen = []
    for a in atoms:
        k = (a.chain, a.resseq, a.icode, a.resname)
        if k not in seen:
            seen.append(k)
    return seen


class AltlocTest(unittest.TestCase):
    def test_highest_occupancy_wins(self):
        pdb = read_pdb(fixture("1N8Z_mini.pdb"))
        out = clean_atoms(pdb, opts(keep=["A"]))
        cb = [a for a in out.atoms if a.resseq == 2 and a.name == "CB"]
        self.assertEqual(len(cb), 1)
        self.assertAlmostEqual(cb[0].occupancy, 0.60)
        self.assertEqual(cb[0].altloc, " ")
        self.assertIn("PDB_ALTLOC_RESOLVED", out.result.codes())
        self.assertFalse(any(a.altloc != " " for a in out.atoms))

    def test_tie_keeps_first_seen(self):
        pdb = read_pdb(fixture("1N8Z_mini.pdb"))
        out = clean_atoms(pdb, opts(keep=["B"]))
        cb = [a for a in out.atoms if a.resseq == 3 and a.name == "CB"]
        self.assertEqual(len(cb), 1)
        # first altloc (A) has the original x; the injected B was shifted by +0.5
        src = [a for a in pdb.atoms if a.chain == "B" and a.resseq == 3 and a.name == "CB"]
        self.assertAlmostEqual(cb[0].xyz[0], src[0].xyz[0])

    def test_select_altlocs_key_matches_prototype(self):
        lines = [atom_line(1, "CB", "ALA", "A", 1, altloc="A", occ=0.3), atom_line(2, "CB", "ALA", "A", 1, altloc="B", occ=0.7)]
        best = select_altlocs([a for a in parse(lines).atoms])
        self.assertEqual(list(best.values()), [(0.7, "B")])

    def test_residue_heterogeneity_cso_cys(self):
        pdb = read_pdb(fixture("6B9Z_cso_mini.pdb"))
        # default: CSO (HETATM) is removed, so the CYS conformer (occ 0.25) must survive intact
        out = clean_atoms(pdb, opts())
        r175 = [a for a in out.atoms if a.resseq == 175]
        self.assertEqual({a.resname for a in r175}, {"CYS"})
        self.assertEqual(sorted(a.name for a in r175), ["C", "CA", "CB", "N", "O", "SG"])
        self.assertIn("PDB_MODIFIED_CONFORMER_REMOVED", out.result.codes())
        # the gap warning in this fixture belongs to the synthetic MSE at B177, never to B175
        removed = [f.where for f in out.result.findings if f.code == "PDB_MODIFIED_RESIDUE_REMOVED"]
        self.assertEqual(removed, ["B:177"])
        # with --std-residues the higher-occupancy CSO conformer wins and becomes CYS
        out2 = clean_atoms(pdb, opts(std_residues=True))
        r175 = [a for a in out2.atoms if a.resseq == 175]
        self.assertEqual({a.resname for a in r175}, {"CYS"})
        self.assertTrue(all(abs(a.occupancy - 0.75) < 1e-6 for a in r175))
        self.assertEqual(sorted(a.name for a in r175), ["C", "CA", "CB", "N", "O", "SG"])
        self.assertIn("PDB_RESIDUE_HETEROGENEITY", out2.result.codes())
        self.assertTrue(all(a.record == "ATOM" for a in r175))


class HetatmTest(unittest.TestCase):
    def test_hetatm_removed_by_default_and_counted(self):
        pdb = read_pdb(fixture("1N8Z_mini.pdb"))
        out = clean_atoms(pdb, opts())
        self.assertFalse(any(a.is_hetatm for a in out.atoms))
        self.assertEqual(out.details["hetatm_removed"], {"HOH": 2})

    def test_keep_het(self):
        pdb = read_pdb(fixture("1N8Z_mini.pdb"))
        out = clean_atoms(pdb, opts(keep_het=True))
        self.assertEqual(sum(1 for a in out.atoms if a.is_hetatm), 2)

    def test_polymer_like_detection(self):
        pdb = read_pdb(fixture("6B9Z_cso_mini.pdb"))
        poly = polymer_like_hetero_residues(pdb.first_model_atoms())
        self.assertEqual(sorted(poly.values()), ["CSO", "MSE"])
        water = parse(residue(1, "ALA", "A", 1) + [atom_line(9, "O", "HOH", "A", 300, hetatm=True)])
        self.assertEqual(polymer_like_hetero_residues(water.atoms), {})

    def test_modified_residue_gap_warning(self):
        # a lone MSE inside the selection, no standard conformer -> warning
        lines = residue(1, "ALA", "A", 1) + residue(5, "MSE", "A", 2, x=3.8, names=("N", "CA", "C", "O", "CB", "SE"), hetatm=True) \
            + residue(11, "GLY", "A", 3, x=7.6)
        out = clean_atoms(parse(lines), opts(keep=["A:1-3"]))
        self.assertIn("PDB_MODIFIED_RESIDUE_REMOVED", out.result.codes())
        self.assertEqual([r[1] for r in residues_of(out.atoms)], [1, 3])
        # outside the selection: silent
        out = clean_atoms(parse(lines), opts(keep=["A:1"]))
        self.assertNotIn("PDB_MODIFIED_RESIDUE_REMOVED", out.result.codes())

    def test_modified_residue_without_rule(self):
        lines = residue(1, "ALA", "A", 1) + residue(5, "XYZ", "A", 2, x=3.8, names=("N", "CA", "C", "O", "CB", "XX"), hetatm=True)
        out = clean_atoms(parse(lines), opts())
        self.assertIn("PDB_MODIFIED_RESIDUE_NO_RULE", out.result.codes())

    def test_std_residues_conversion_mse(self):
        pdb = read_pdb(fixture("6B9Z_cso_mini.pdb"))
        out = clean_atoms(pdb, opts(std_residues=True))
        r177 = [a for a in out.atoms if a.resseq == 177]
        self.assertEqual({a.resname for a in r177}, {"MET"})
        self.assertIn("SD", [a.name for a in r177])
        self.assertNotIn("SE", [a.name for a in r177])
        sd = next(a for a in r177 if a.name == "SD")
        self.assertEqual(sd.element, "S")
        self.assertEqual(sd.line[12:16], " SD ")
        self.assertTrue(all(a.record == "ATOM" for a in r177))
        codes = out.result.codes()
        self.assertEqual(codes.count("PDB_MODIFIED_RESIDUE_CONVERTED"), 2)   # CSO and MSE
        self.assertEqual([d["from"] for d in out.details["modified_residues_converted"]], ["CSO", "MSE"])

    def test_std_residues_drops_only_extra_atoms(self):
        pdb = read_pdb(fixture("6B9Z_cso_mini.pdb"))
        out = clean_atoms(pdb, opts(std_residues=True))
        r175 = [a.name for a in out.atoms if a.resseq == 175]
        self.assertNotIn("OD", r175)
        self.assertEqual(len(r175), 6)


class SelectionTest(unittest.TestCase):
    def test_insertion_codes(self):
        pdb = read_pdb(fixture("kabat_mini.pdb"))
        keys = [(r[1], r[2]) for r in residues_of(clean_atoms(pdb, opts(keep=["B:80-82"])).atoms)]
        self.assertEqual(keys, [(80, " "), (81, " "), (82, " "), (82, "A"), (82, "B"), (82, "C")])
        keys = [(r[1], r[2]) for r in residues_of(clean_atoms(pdb, opts(keep=["B:80-82A"])).atoms)]
        self.assertEqual(keys, [(80, " "), (81, " "), (82, " "), (82, "A")])
        keys = [(r[1], r[2]) for r in residues_of(clean_atoms(pdb, opts(keep=["B:82B-83"])).atoms)]
        self.assertEqual(keys, [(82, "B"), (82, "C"), (83, " ")])

    def test_whole_chain_and_order_preserved(self):
        pdb = read_pdb(fixture("1N8Z_mini.pdb"))
        out = clean_atoms(pdb, opts(keep=["C", "A"]))
        self.assertEqual([r[0] for r in residues_of(out.atoms)][:1], ["A"])    # file order, not option order
        self.assertEqual({r[0] for r in residues_of(out.atoms)}, {"A", "C"})

    def test_chain_not_found_is_error(self):
        pdb = read_pdb(fixture("1N8Z_mini.pdb"))
        out = clean_atoms(pdb, opts(keep=["Z"]))
        self.assertIn("PDB_CHAIN_NOT_FOUND", out.result.codes())
        self.assertTrue(out.result.has_errors)

    def test_empty_range_is_error(self):
        pdb = read_pdb(fixture("1N8Z_mini.pdb"))
        out = clean_atoms(pdb, opts(keep=["A:900-950"]))
        self.assertIn("PDB_SELECTION_EMPTY", out.result.codes())

    def test_rename_and_collision(self):
        pdb = read_pdb(fixture("1N8Z_mini.pdb"))
        out = clean_atoms(pdb, opts(keep=["A", "B"], rename="A=L,B=H"))
        self.assertEqual(list(out.details["chains_out"]), ["L", "H"])
        self.assertNotIn("PDB_RENAME_COLLISION", out.result.codes())
        out = clean_atoms(pdb, opts(keep=["A", "B"], rename="A=B"))   # A becomes B while B stays B
        self.assertIn("PDB_RENAME_COLLISION", out.result.codes())

    def test_swap_rename(self):
        pdb = read_pdb(fixture("1N8Z_mini.pdb"))
        out = clean_atoms(pdb, opts(keep=["A", "B"], rename="A=B,B=A"))
        self.assertEqual(list(out.details["chains_out"]), ["B", "A"])
        self.assertNotIn("PDB_RENAME_COLLISION", out.result.codes())


class ModelAndHydrogenTest(unittest.TestCase):
    def test_first_model_only(self):
        pdb = read_pdb(fixture("6B9Z_cso_mini.pdb"))
        out = clean_atoms(pdb, opts())
        self.assertIn("PDB_MULTI_MODEL", out.result.codes())
        self.assertEqual(out.details["models_in_file"], 2)
        # every output residue appears once
        keys = [(r[1], r[2]) for r in residues_of(out.atoms)]
        self.assertEqual(len(keys), len(set(keys)))

    def test_strip_h_by_element_and_by_name(self):
        lines = residue(1, "ALA", "A", 1) + [
            atom_line(5, "H", "ALA", "A", 1, element="H"),
            atom_line(6, "HB1", "ALA", "A", 1, element=""),   # blank element column -> name rule
            atom_line(7, "HG", "SER", "A", 2, element=""),    # hydrogen gamma in a standard residue
        ]
        out = clean_atoms(parse(lines), opts(strip_h=True))
        self.assertEqual(len(out.atoms), 4)
        self.assertEqual(out.details["hydrogens_removed"], 3)
        out = clean_atoms(parse(lines), opts())
        self.assertEqual(len(out.atoms), 7)


class CleanFileTest(unittest.TestCase):
    def test_end_to_end_regenerates_ter(self):
        pdb_path = fixture("1N8Z_mini.pdb")
        out_path = os.path.join(tempfile.mkdtemp(), "out.pdb")
        outcome = clean_file(pdb_path, out_path, opts(keep=["A:1-3", "B:1-3"], rename="A=L,B=H"))
        self.assertEqual(len(outcome.ters), 2)
        with open(out_path) as fh:
            lines = fh.read().splitlines()
        ters = [l for l in lines if l.startswith("TER")]
        self.assertEqual([t[17:26] for t in ters], ["GLN L   3", "GLN H   3"])
        self.assertEqual(lines[-1].rstrip(), "END")
        self.assertFalse(any(l.startswith(("REMARK", "HETATM")) for l in lines))


if __name__ == "__main__":
    unittest.main()
