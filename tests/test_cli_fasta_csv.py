"""CLI: ``abd-prep fasta ...`` and ``abd-prep csv ...`` (plan §6 #5-#8) + fixture provenance."""
import contextlib
import io
import json
import os
import subprocess
import sys
import tempfile
import unittest

from abd_prep.cli import main
from tests.helpers import CACHE, fixture
from tests.make_seq_fixtures import VH, VL, build, FIXTURES


def run(*argv):
    out, err = io.StringIO(), io.StringIO()
    with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
        code = main(list(argv))
    return code, out.getvalue(), err.getvalue()


class FastaCleanCliTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.out = os.path.join(self.tmp, "seed_clean.fasta")

    def test_dod5_seed_igg(self):
        code, out, err = run("fasta", "clean", fixture("seed_igg.fasta"), self.out, "--framework", "igg")
        self.assertEqual(code, 0, err)
        with open(self.out) as fh:
            lines = fh.read().splitlines()
        self.assertEqual(lines[0].split()[0], ">H")
        self.assertEqual(lines[1], VH)
        self.assertEqual(lines[2].split()[0], ">L")
        self.assertEqual(lines[3], VL)
        self.assertIn("records: H(120 res), L(107 res)", out)
        with open(os.path.join(self.tmp, "seed_clean.report.json")) as fh:
            rep = json.load(fh)
        self.assertEqual(rep["inputs"][0]["kind"], "fasta")
        self.assertEqual(rep["outputs"][0]["path"], "seed_clean.fasta")
        self.assertEqual([r["role"] for r in rep["details"]["records"]], ["H", "L"])
        self.assertIn("FASTA_CRLF", [f["code"] for f in rep["findings"]])

    def test_dod5_fulllength_warnings_and_strict(self):
        args = ("fasta", "clean", fixture("seed_warn_fulllength.fasta"), self.out, "--framework", "igg",
                "--heavy-id", "heavy_full", "--light-id", "light_full", "--no-report")
        code, out, _ = run(*args)
        self.assertEqual(code, 0)
        self.assertIn("2 warning(s)", out)
        self.assertEqual(run(*args, "--strict")[0], 1)

    def test_framework_required_for_clean(self):
        code, _, err = run("fasta", "clean", fixture("seed_igg.fasta"), self.out)
        self.assertEqual(code, 2)

    def test_errors_still_write_output_and_exit_1(self):
        code, _, _ = run("fasta", "clean", fixture("seed_bad_alphabet.fasta"), self.out, "--framework", "igg",
                         "--no-report")
        self.assertEqual(code, 1)
        self.assertTrue(os.path.exists(self.out))          # like pdb clean: written, but reported as failing

    def test_replace_ambiguous_and_wrap(self):
        code, out, _ = run("fasta", "clean", fixture("seed_bad_alphabet.fasta"), self.out, "--framework", "igg",
                           "--replace-ambiguous", "X=G,B=N,Z=Q", "--wrap", "60", "--no-report")
        self.assertEqual(code, 1)                           # internal '*' remains an error
        self.assertIn("FASTA_AMBIGUOUS_REPLACED", out)
        with open(self.out) as fh:
            self.assertTrue(all(len(ln) <= 60 for ln in fh.read().splitlines() if not ln.startswith(">")))

    def test_same_file_and_missing_file(self):
        self.assertEqual(run("fasta", "clean", fixture("seed_igg.fasta"), fixture("seed_igg.fasta"), "--framework", "igg")[0], 2)
        self.assertEqual(run("fasta", "clean", "/nonexistent.fasta", self.out, "--framework", "igg")[0], 2)

    def test_not_utf8_input_exit_1_no_output(self):
        bad = os.path.join(self.tmp, "cp949.fasta")
        with open(bad, "wb") as fh:
            fh.write(">중쇄\nEVQL\n".encode("cp949"))
        code, out, _ = run("fasta", "clean", bad, self.out, "--framework", "nanobody", "--no-report")
        self.assertEqual(code, 1)
        self.assertIn("FILE_NOT_UTF8", out)
        self.assertFalse(os.path.exists(self.out))
        code, out, _ = run("fasta", "clean", bad, self.out, "--framework", "nanobody", "--no-report",
                           "--encoding", "cp949")
        self.assertEqual(code, 0, out)
        self.assertIn("FILE_ENCODING_CONVERTED", out)


class FastaValidateCliTest(unittest.TestCase):
    def test_dod6_error_fixtures(self):
        expected = {"seed_bad_alphabet.fasta": "FASTA_BAD_ALPHABET", "seed_bad_nucleic.fasta": "FASTA_NUCLEIC",
                    "seed_bad_count.fasta": "FASTA_FRAMEWORK_MISMATCH"}
        for name, code_name in expected.items():
            code, out, _ = run("fasta", "validate", fixture(name), "--framework", "igg", "--no-report")
            self.assertEqual(code, 1, name)
            self.assertIn(code_name, out, name)
        code, out, _ = run("fasta", "validate", fixture("seed_bad_alphabet.fasta"), "--framework", "igg", "--no-report")
        self.assertIn("FASTA_INTERNAL_STOP", out)

    def test_multiple_files_prefix_where_and_json(self):
        code, out, _ = run("fasta", "validate", fixture("seed_nanobody.fasta"), fixture("seed_scfv.fasta"),
                           "--json", "--no-report")
        self.assertEqual(code, 0)
        rep = json.loads(out)
        self.assertEqual(len(rep["inputs"]), 2)
        self.assertTrue(all(f["where"].split(":")[0].endswith(".fasta") for f in rep["findings"]))
        self.assertEqual(set(rep["details"]["files"]), {"seed_nanobody.fasta", "seed_scfv.fasta"})

    def test_extension_warning(self):
        tmp = tempfile.mkdtemp()
        path = os.path.join(tmp, "seed.txt")
        with open(path, "w") as fh:
            fh.write(">VHH\n" + VH + "\n")
        code, out, _ = run("fasta", "validate", path, "--framework", "nanobody", "--no-report")
        self.assertEqual(code, 0)
        self.assertIn("FILE_EXTENSION", out)
        self.assertEqual(run("fasta", "validate", path, "--framework", "nanobody", "--no-report", "--strict")[0], 1)

    def test_against_option_errors(self):
        self.assertEqual(run("fasta", "validate", fixture("seed_igg.fasta"), fixture("seed_scfv.fasta"),
                             "--against", fixture("1N8Z_mini.pdb"))[0], 2)
        self.assertEqual(run("fasta", "validate", fixture("seed_igg.fasta"), "--chain-map", "H=B")[0], 2)
        self.assertEqual(run("fasta", "validate", fixture("seed_igg.fasta"), "--against", "/nonexistent.pdb")[0], 2)

    def test_against_mini_fixture_needs_chain_map(self):
        # 1N8Z_mini has chains A/B/C with the first residues of each 1N8Z chain
        code, _, err = run("fasta", "validate", fixture("seed_igg.fasta"), "--against", fixture("1N8Z_mini.pdb"),
                           "--no-report")
        self.assertEqual(code, 2)
        self.assertIn("--chain-map", err)
        code, out, _ = run("fasta", "validate", fixture("seed_igg.fasta"), "--against", fixture("1N8Z_mini.pdb"),
                           "--chain-map", "H=B,L=A", "--framework", "igg", "--json", "--no-report")
        rep = json.loads(out)
        codes = [f["code"] for f in rep["findings"]]
        self.assertNotIn("FASTA_PDB_MISMATCH", codes)
        self.assertEqual(codes.count("FASTA_PDB_IDENTITY"), 2)
        self.assertIn("against", rep["details"]["files"]["seed_igg.fasta"])
        # the --against path is absolute here; basename mode must not leak its directory anywhere in the report
        self.assertEqual(rep["details"]["files"]["seed_igg.fasta"]["against"]["pdb"], "1N8Z_mini.pdb")
        self.assertNotIn(os.path.dirname(fixture("1N8Z_mini.pdb")), out)

    def test_inspect_text_and_json(self):
        code, out, _ = run("fasta", "inspect", fixture("seed_bad_count.fasta"), "--framework", "igg")
        self.assertEqual(code, 0)
        self.assertIn("records: 3", out)
        self.assertIn("heavy: role H (header), 120 res", out)
        code, out, _ = run("fasta", "inspect", fixture("seed_scfv.fasta"), "--json")
        info = json.loads(out)
        self.assertEqual(info["records"][0]["length"], 242)


@unittest.skipUnless(os.path.isfile(os.path.join(CACHE, "1N8Z.pdb")), "1N8Z not in the local cache")
class FastaAgainstPresetCliTest(unittest.TestCase):
    def test_dod7_cli(self):
        tmp = tempfile.mkdtemp()
        self.assertEqual(run("pdb", "split", "--preset", "her2-trastuzumab", "--out", tmp, "--no-report")[0], 0)
        clean = os.path.join(tmp, "seed_clean.fasta")
        self.assertEqual(run("fasta", "clean", fixture("seed_igg.fasta"), clean, "--framework", "igg", "--no-report")[0], 0)
        code, out, _ = run("fasta", "validate", clean, "--against", os.path.join(tmp, "her2-trastuzumab_fv.pdb"),
                           "--no-report")
        self.assertEqual(code, 0, out)
        self.assertIn("120/120", out)
        self.assertIn("107/107", out)
        self.assertEqual(run("fasta", "validate", clean, "--against", os.path.join(CACHE, "1N8Z.pdb"), "--no-report")[0], 2)
        self.assertEqual(run("fasta", "validate", clean, "--against", os.path.join(CACHE, "1N8Z.pdb"),
                             "--chain-map", "H=B,L=A", "--no-report")[0], 0)


class CsvCliTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.out = os.path.join(self.tmp, "assay_clean.csv")

    def test_dod8_messy_utf8bom(self):
        code, out, _ = run("csv", "clean", fixture("assay_messy_utf8bom.csv"), self.out,
                           "--map", "id=Clone,sequence=AA_seq,value=KD,unit=Unit")
        self.assertEqual(code, 0, out)
        for code_name in ("FILE_BOM", "CSV_DELIMITER", "CSV_MISSING_VALUES", "CSV_VALUE_NOT_NUMERIC",
                          "CSV_DECIMAL_COMMA", "CSV_UNIT_MIXED", "CSV_ID_DUPLICATE"):
            self.assertIn(code_name, out, code_name)
        with open(self.out, "rb") as fh:
            raw = fh.read()
        self.assertTrue(raw.startswith(b"id,sequence,value,unit,Note\n"))
        with open(os.path.join(self.tmp, "assay_clean.report.json")) as fh:
            rep = json.load(fh)
        self.assertEqual(rep["inputs"][0]["kind"], "csv")
        self.assertEqual(rep["details"]["rows_written"], 8)

    def test_dod8_cp949(self):
        args = ("csv", "clean", fixture("assay_messy_cp949.csv"), self.out, "--map", "id=클론,sequence=서열,value=KD",
                "--no-report")
        code, out, _ = run(*args, "--encoding", "cp949")
        self.assertEqual(code, 1)
        self.assertEqual(out.count("ERROR"), 1)
        self.assertIn("CSV_ROW_LENGTH", out)
        with open(self.out, encoding="utf-8") as fh:
            self.assertEqual(fh.readline().rstrip("\n"), "id,sequence,value,단위,비고")
        os.remove(self.out)
        code, out, _ = run(*args)
        self.assertEqual(code, 1)
        self.assertIn("FILE_NOT_UTF8", out)
        self.assertFalse(os.path.exists(self.out))

    def test_dod8_assay_ok_validate(self):
        code, out, _ = run("csv", "validate", fixture("assay_ok.csv"), "--map",
                           "id=id,sequence=sequence,value=value,unit=unit", "--no-report")
        self.assertEqual(code, 0)
        self.assertIn("0 error(s), 0 warning(s)", out)

    def test_validate_against_fasta_and_pdb(self):
        path = os.path.join(self.tmp, "x.csv")
        with open(path, "w") as fh:
            fh.write(f"id,sequence,value\nheavy,{VH},1\nAbX,{VH[:50]},2\n")
        code, out, _ = run("csv", "validate", path, "--against", fixture("seed_igg.fasta"), "--json", "--no-report")
        self.assertEqual(code, 0)
        rep = json.loads(out)
        against = rep["details"]["files"]["x.csv"]["against"]
        self.assertEqual((against["matched"], against["unmatched"], against["kind"]), (1, 1, "fasta"))
        self.assertEqual(rep["inputs"][1]["kind"], "fasta")
        code, out, _ = run("csv", "validate", path, "--against", fixture("1N8Z_mini.pdb"), "--json", "--no-report")
        self.assertEqual(json.loads(out)["details"]["files"]["x.csv"]["against"]["kind"], "pdb")
        self.assertEqual(run("csv", "validate", path, fixture("assay_ok.csv"), "--against", fixture("seed_igg.fasta"))[0], 2)

    def test_option_errors_exit_2(self):
        self.assertEqual(run("csv", "clean", fixture("assay_ok.csv"), self.out, "--map", "clone=id")[0], 2)
        self.assertEqual(run("csv", "clean", fixture("assay_ok.csv"), self.out, "--map", "id=NoSuchColumn")[0], 2)
        self.assertEqual(run("csv", "clean", fixture("assay_ok.csv"), self.out, "--delimiter", "ab")[0], 2)
        self.assertEqual(run("csv", "validate", fixture("assay_ok.csv"), "--encoding", "klingon")[0], 2)
        xlsx = os.path.join(self.tmp, "assay.xlsx")
        with open(xlsx, "wb") as fh:
            fh.write(b"PK\x03\x04")
        code, _, err = run("csv", "validate", xlsx)
        self.assertEqual(code, 2)
        self.assertIn("CSV UTF-8", err)

    def test_inspect(self):
        code, out, _ = run("csv", "inspect", fixture("assay_messy_utf8bom.csv"), "--map", "id=Clone")
        self.assertEqual(code, 0)
        self.assertIn("delimiter: ;  rows: 8", out)
        self.assertIn("Clone -> id", out)
        code, out, _ = run("csv", "inspect", fixture("assay_messy_cp949.csv"), "--encoding", "cp949", "--json")
        info = json.loads(out)
        self.assertEqual(info["bad_length_rows"], 1)
        self.assertEqual(info["columns"][0]["name"], "클론")

    def test_tab_delimiter_and_tsv_extension(self):
        path = os.path.join(self.tmp, "a.tsv")
        with open(path, "w") as fh:
            fh.write("id\tvalue\nAb1\t1\n")
        code, out, _ = run("csv", "validate", path, "--delimiter", "tab", "--no-report")
        self.assertEqual(code, 0, out)
        self.assertNotIn("FILE_EXTENSION", out)
        self.assertIn("delimiter tab", out)


class FixtureProvenanceTest(unittest.TestCase):
    """The committed FASTA/CSV fixtures are exactly what tests/make_seq_fixtures.py generates."""

    def test_fixtures_match_generator(self):
        for name, data in build().items():
            with open(os.path.join(FIXTURES, name), "rb") as fh:
                self.assertEqual(fh.read(), data, name)

    def test_generator_check_mode(self):
        proc = subprocess.run([sys.executable, os.path.join(os.path.dirname(FIXTURES), "make_seq_fixtures.py"), "--check"],
                              capture_output=True, text=True)
        self.assertEqual(proc.returncode, 0, proc.stdout + proc.stderr)

    def test_fixture_byte_features(self):
        with open(fixture("seed_igg.fasta"), "rb") as fh:
            self.assertIn(b"\r\n", fh.read())
        with open(fixture("assay_messy_utf8bom.csv"), "rb") as fh:
            self.assertTrue(fh.read().startswith(b"\xef\xbb\xbf"))
        with open(fixture("assay_messy_cp949.csv"), "rb") as fh:
            raw = fh.read()
        self.assertFalse(raw.startswith(b"\xef\xbb\xbf"))
        with self.assertRaises(UnicodeDecodeError):
            raw.decode("utf-8")


if __name__ == "__main__":
    unittest.main()
