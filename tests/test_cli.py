import contextlib
import io
import json
import os
import tempfile
import unittest

from abd_prep import __version__
from abd_prep.cli import main
from tests.helpers import fixture, residue, ter_line, write_temp


def run(*argv):
    out, err = io.StringIO(), io.StringIO()
    with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
        code = main(list(argv))
    return code, out.getvalue(), err.getvalue()


class CleanCommandTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.out = os.path.join(self.tmp, "fv.pdb")

    def test_prototype_arguments(self):
        code, out, err = run("pdb", "clean", fixture("1N8Z_mini.pdb"), self.out,
                             "--keep", "A:1-3", "--keep", "B:1-3", "--rename", "A=L,B=H")
        self.assertEqual(code, 0, err)
        self.assertIn("chains: L(3 res), H(3 res)", out)
        self.assertTrue(os.path.exists(self.out))
        report_path = os.path.join(self.tmp, "fv.report.json")
        self.assertTrue(os.path.exists(report_path))
        with open(report_path) as fh:
            rep = json.load(fh)
        self.assertEqual(set(rep), {"tool", "command", "generated_at", "inputs", "outputs", "findings", "details"})
        self.assertEqual(rep["tool"], {"name": "abd-prep", "version": __version__, "python": rep["tool"]["python"]})
        self.assertEqual(rep["inputs"][0]["path"], "1N8Z_mini.pdb")     # basename by default
        self.assertEqual(rep["outputs"][0]["path"], "fv.pdb")
        self.assertEqual(len(rep["inputs"][0]["sha256"]), 64)
        self.assertNotIn(self.tmp, rep["command"])
        self.assertIn("abd-prep pdb clean", rep["command"])
        self.assertEqual(rep["details"]["rename"], {"A": "L", "B": "H"})

    def test_full_paths_and_custom_report(self):
        rp = os.path.join(self.tmp, "r.json")
        code, _, _ = run("pdb", "clean", fixture("1N8Z_mini.pdb"), self.out, "--report", rp, "--report-paths", "full")
        self.assertEqual(code, 0)
        with open(rp) as fh:
            rep = json.load(fh)
        self.assertEqual(rep["outputs"][0]["path"], os.path.abspath(self.out))
        self.assertFalse(os.path.exists(os.path.join(self.tmp, "fv.report.json")))

    def test_no_report(self):
        code, _, _ = run("pdb", "clean", fixture("1N8Z_mini.pdb"), self.out, "--no-report")
        self.assertEqual(code, 0)
        self.assertEqual([f for f in os.listdir(self.tmp) if f.endswith(".json")], [])

    def test_json_output(self):
        code, out, _ = run("pdb", "clean", fixture("1N8Z_mini.pdb"), self.out, "--json", "--no-report")
        self.assertEqual(code, 0)
        rep = json.loads(out)
        self.assertIn("PDB_ALTLOC_RESOLVED", [f["code"] for f in rep["findings"]])

    def test_exit_codes(self):
        # validation error (chain missing) -> 1, file still written
        code, _, _ = run("pdb", "clean", fixture("1N8Z_mini.pdb"), self.out, "--keep", "Z", "--no-report")
        self.assertEqual(code, 1)
        # warning only -> 0; with --strict -> 1
        code, _, _ = run("pdb", "clean", fixture("6B9Z_cso_mini.pdb"), self.out, "--no-report")
        self.assertEqual(code, 0)
        code, _, _ = run("pdb", "clean", fixture("6B9Z_cso_mini.pdb"), self.out, "--no-report", "--strict")
        self.assertEqual(code, 1)
        # input errors -> 2
        code, _, err = run("pdb", "clean", os.path.join(self.tmp, "missing.pdb"), self.out)
        self.assertEqual(code, 2)
        self.assertIn("file not found", err)
        code, _, err = run("pdb", "clean", fixture("1N8Z_mini.pdb"), self.out, "--keep", "A:5-1")
        self.assertEqual(code, 2)
        code, _, err = run("pdb", "clean", fixture("1N8Z_mini.pdb"), fixture("1N8Z_mini.pdb"))
        self.assertEqual(code, 2)
        self.assertIn("same file", err)

    def test_std_residues_flag(self):
        code, out, _ = run("pdb", "clean", fixture("6B9Z_cso_mini.pdb"), self.out, "--std-residues", "--no-report")
        self.assertEqual(code, 0)
        self.assertIn("PDB_MODIFIED_RESIDUE_CONVERTED", out)
        code, out, _ = run("pdb", "validate", self.out)
        self.assertEqual(code, 0, out)

    def test_extension_warning(self):
        src = fixture("1N8Z_mini.pdb")
        odd = os.path.join(self.tmp, "structure.txt")
        with open(src) as a, open(odd, "w") as b:
            b.write(a.read())
        code, out, _ = run("pdb", "clean", odd, self.out, "--no-report")
        self.assertEqual(code, 0)
        self.assertIn("PDB_EXTENSION", out)


class ValidateCommandTest(unittest.TestCase):
    def test_broken_ter_fails(self):
        lines = residue(1, "ALA", "L", 1) + [ter_line(5, "CYS", "A", 214), "END"]
        path = write_temp(lines)
        code, out, _ = run("pdb", "validate", path)
        self.assertEqual(code, 1)
        self.assertIn("PDB_TER_MISMATCH", out)

    def test_multiple_files_and_expect_chains(self):
        good = write_temp(residue(1, "ALA", "A", 1) + [ter_line(5, "ALA", "A", 1), "END"])
        code, out, _ = run("pdb", "validate", good, good, "--expect-chains", "1")
        self.assertEqual(code, 0, out)
        code, out, _ = run("pdb", "validate", good, "--expect-chains", "2")
        self.assertEqual(code, 1)
        self.assertIn("PDB_CHAIN_COUNT", out)

    def test_json_and_strict(self):
        code, out, _ = run("pdb", "validate", fixture("kabat_mini.pdb"), "--json")
        self.assertEqual(code, 0)
        rep = json.loads(out)
        self.assertEqual(rep["inputs"][0]["kind"], "pdb")
        # fixture has a chain named B with 7 residues: no H/L count rule -> no warnings -> strict still 0
        code, _, _ = run("pdb", "validate", fixture("kabat_mini.pdb"), "--strict")
        self.assertEqual(code, 0)
        # 1N8Z_mini has altlocs (error) and HETATM (warning)
        code, _, _ = run("pdb", "validate", fixture("1N8Z_mini.pdb"))
        self.assertEqual(code, 1)


class InspectCommandTest(unittest.TestCase):
    def test_text_fasta_json(self):
        code, out, _ = run("pdb", "inspect", fixture("1N8Z_mini.pdb"))
        self.assertEqual(code, 0)
        self.assertIn("chain C: 489-493", out)
        code, out, _ = run("pdb", "inspect", fixture("1N8Z_mini.pdb"), "--fasta")
        self.assertEqual(out.count(">"), 3)
        self.assertIn("\nDIQMT\n", out)
        code, out, _ = run("pdb", "inspect", fixture("1N8Z_mini.pdb"), fixture("kabat_mini.pdb"), "--json")
        self.assertEqual(len(json.loads(out)), 2)


class SplitCommandTest(unittest.TestCase):
    def setUp(self):
        self.out = tempfile.mkdtemp()

    def test_local_input_explicit_parts(self):
        code, out, err = run("pdb", "split", "--input", fixture("1N8Z_mini.pdb"), "--target", "C:489-493",
                             "--heavy", "B:1-5", "--light", "A:1-5", "--rename-heavy", "H", "--rename-light", "L",
                             "--name", "mini", "--out", self.out)
        self.assertEqual(code, 0, err + out)
        self.assertEqual(sorted(os.listdir(self.out)), ["mini.report.json", "mini_fv.pdb", "mini_target.pdb"])
        self.assertIn("[split] mini", out)
        with open(os.path.join(self.out, "mini.report.json")) as fh:
            rep = json.load(fh)
        self.assertTrue(rep["command"].startswith("abd-prep pdb split"))
        self.assertNotIn(self.out, rep["command"])

    def test_default_name_from_input_and_json(self):
        code, out, _ = run("pdb", "split", "--input", fixture("kabat_mini.pdb"), "--heavy", "B:80-83",
                           "--rename-heavy", "H", "--out", self.out, "--json", "--no-report")
        self.assertEqual(code, 0, out)
        rep = json.loads(out)
        self.assertEqual([o["path"] for o in rep["outputs"]], ["kabat_mini_fv.pdb"])
        self.assertEqual(os.listdir(self.out), ["kabat_mini_fv.pdb"])

    def test_usage_errors(self):
        code, _, err = run("pdb", "split", "--preset", "her2-trastuzumab", "--target", "C:1-5", "--out", self.out)
        self.assertEqual(code, 2)
        self.assertIn("--preset cannot be combined", err)
        code, _, err = run("pdb", "split", "--pdb", "1N8Z", "--out", self.out)
        self.assertEqual(code, 2)
        self.assertIn("nothing to split", err)
        code, _, err = run("pdb", "split", "--input", fixture("1N8Z_mini.pdb"), "--rename-heavy", "H", "--out", self.out)
        self.assertEqual(code, 2)
        self.assertIn("without --heavy", err)
        code, _, err = run("pdb", "split", "--pdb", "bad-id", "--target", "A", "--out", self.out)
        self.assertEqual(code, 2)
        self.assertIn("not a PDB ID", err)
        code, _, err = run("pdb", "split", "--preset", "nope", "--out", self.out)
        self.assertEqual(code, 2)                               # argparse choices


class FetchCommandTest(unittest.TestCase):
    def test_invalid_id_is_exit_2_without_network(self):
        code, _, err = run("pdb", "fetch", "../etc")
        self.assertEqual(code, 2)
        self.assertIn("not a PDB ID", err)


class ValidatePresetTest(unittest.TestCase):
    def test_preset_match_and_no_match(self):
        out = tempfile.mkdtemp()
        # build a two-chain H/L file from the fixture: matches the preset's fv part
        code, _, _ = run("pdb", "clean", fixture("1N8Z_mini.pdb"), os.path.join(out, "fv.pdb"),
                         "--keep", "A:1-5", "--keep", "B:1-5", "--rename", "A=L,B=H", "--no-report")
        self.assertEqual(code, 0)
        code, text, _ = run("pdb", "validate", os.path.join(out, "fv.pdb"), "--preset", "her2-trastuzumab", "--json")
        self.assertEqual(code, 0, text)
        rep = json.loads(text)
        self.assertEqual(rep["details"]["files"]["fv.pdb"]["preset_part"], "fv")
        codes = [f["code"] for f in rep["findings"]]
        self.assertIn("PDB_RESIDUE_COUNT", codes)             # 5 residues vs the preset's 107/120
        self.assertNotIn("PDB_PRESET_NO_MATCH", codes)
        # chain C alone is the target part; chain B alone matches nothing
        code, text, _ = run("pdb", "validate", fixture("kabat_mini.pdb"), "--preset", "her2-trastuzumab")
        self.assertEqual(code, 0)
        self.assertIn("PDB_PRESET_NO_MATCH", text)


class ParserTest(unittest.TestCase):
    def test_version_and_usage(self):
        code, out, _ = run("--version")
        self.assertEqual(code, 0)
        self.assertIn(__version__, out)
        code, _, err = run()
        self.assertEqual(code, 2)
        code, _, err = run("pdb")
        self.assertEqual(code, 2)
        code, _, err = run("fasta", "clean", "x", "y")      # group not implemented yet (M2.5)
        self.assertEqual(code, 2)


if __name__ == "__main__":
    unittest.main()
