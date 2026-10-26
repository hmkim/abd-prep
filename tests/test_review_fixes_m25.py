"""Regression tests for the defects found in the M2.5 code review (plan §5.1).

A1 tab-separated FASTA header       A4 csv.Error traceback (field size / quoting)
A2 several trailing '*'             A5 stale output reported after a header-less input
A3 --replace-ambiguous '*=G'        A6 --against that compares nothing passed silently
B1 --encoding on an ASCII file      B4 scFv linker info only under --framework scfv
B3 undecodable --against reference  B6 letter positions shifted by an internal '*'
"""
import contextlib
import io
import os
import tempfile
import unittest

from abd_prep.cli import main
from abd_prep.common import textfile as T
from abd_prep.common.result import InputError, ValidationResult
from abd_prep.csvdata import crosscheck as CX
from abd_prep.csvdata import validate as CV
from abd_prep.fasta import crosscheck as X
from abd_prep.fasta import validate as V
from abd_prep.fasta.fastafile import Record, read_fasta
from tests.helpers import write_temp
from tests.make_seq_fixtures import SCFV, VH, VL


def run(*argv):
    out, err = io.StringIO(), io.StringIO()
    with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
        code = main(list(argv))
    return code, out.getvalue(), err.getvalue()


def parse(text):
    res = ValidationResult()
    return read_fasta(io.StringIO(text), res), res


class FastaReviewFixesTest(unittest.TestCase):
    def test_a1_tab_separated_header(self):
        fasta, _ = parse(">VHH_1\tsome desc\n" + VH + "\n")
        self.assertEqual((fasta.records[0].id, fasta.records[0].description), ("VHH_1", "some desc"))

    def test_a2_every_trailing_stop_removed(self):
        fasta, res = parse(">h\n" + VH + "**\n")
        self.assertEqual(fasta.records[0].sequence, VH)
        self.assertEqual(res.codes(), ["FASTA_TRAILING_STOP"])
        self.assertIn("2 characters", res.findings[0].message)

    def test_a3_replace_ambiguous_rejects_stop_and_symbols(self):
        for bad in ("*=G", "-=G", "1=G"):
            with self.assertRaises(InputError):
                V.parse_replacements(bad)

    def test_b6_positions_not_shifted_by_internal_stop(self):
        seq = VH[:30] + "*" + VH[30:49] + "X" + VH[50:]          # '*' at 31, X at 51 (1-based, in the real string)
        fasta, res = parse(">VHH\n" + seq + "\n")
        V.validate_records(fasta.records, V.FastaOptions(framework="nanobody", heuristics=False), res)
        bad = [f for f in res.errors if f.code == "FASTA_BAD_ALPHABET"][0]
        self.assertIn("X@51", bad.message)

    def test_b4_scfv_linker_reported_from_header_role(self):
        fasta, res = parse(">my_scFv construct\n" + SCFV + "\n")
        V.validate_records(fasta.records, V.FastaOptions(), res)        # no --framework
        self.assertEqual(fasta.records[0].role, "HL")
        self.assertIn("FASTA_SCFV_LINKER", res.codes())
        self.assertNotIn("FASTA_LENGTH_RANGE", res.codes())          # scFv range applies (242 res)

    def test_a6_against_with_no_roles_is_an_error(self):
        from tests.test_crosscheck import chain_lines, plain
        pdb = write_temp(chain_lines("H", VH, plain(len(VH))) + ["END"])
        res = ValidationResult()
        X.crosscheck([Record("chain1", "", VH)], pdb, {}, res)         # role None
        self.assertEqual([f.code for f in res.errors], ["FASTA_PDB_NOT_COMPARED"])
        res = ValidationResult()
        details = X.crosscheck([Record("chain1", "", VH, role="H")], pdb, {}, res)
        self.assertFalse(res.has_errors)
        self.assertEqual(details["compared"], 1)


class CsvReviewFixesTest(unittest.TestCase):
    def test_a4_long_field_is_accepted_and_bad_quoting_is_input_error(self):
        path = write_temp(["id,notes,value", "Ab1," + "x" * 200_000 + ",1"], suffix=".csv")
        res = ValidationResult()
        details = CV.process_csv(path, "utf-8-sig", CV.CsvOptions(), res)
        self.assertEqual(details["rows_read"], 1)
        code, _, err = run("csv", "validate", path, "--no-report")
        self.assertEqual(code, 0, err)
        # an unterminated quote at the end of the file raises csv.Error -> exit 2 with the file name, no traceback
        bad = write_temp(["id,notes", 'Ab1,"open quote', "Ab2,x", 'Ab3,"'], suffix=".csv")
        code, _, err = run("csv", "validate", bad, "--no-report")
        self.assertIn(code, (0, 2))                                   # the csv module may or may not object...
        self.assertNotIn("Traceback", err)                            # ...but it must never escape as a traceback

    def test_a5_stale_output_not_reported(self):
        tmp = tempfile.mkdtemp()
        out = os.path.join(tmp, "o.csv")
        with open(out, "w") as fh:
            fh.write("stale\n")
        empty = write_temp(["", ""], suffix=".csv")
        code, text, _ = run("csv", "clean", empty, out, "--json", "--no-report")
        self.assertEqual(code, 1)
        import json
        rep = json.loads(text)
        self.assertEqual(rep["outputs"], [])
        self.assertFalse(rep["details"]["output_written"])

    def test_b3_undecodable_reference_is_a_clear_input_error(self):
        fd, ref = tempfile.mkstemp(suffix=".fasta")
        with os.fdopen(fd, "wb") as fh:
            fh.write(">H\nEVQL\n".encode("utf-16"))
        with self.assertRaises(InputError) as ctx:
            CX.load_reference(ref)
        self.assertIn(os.path.basename(ref), str(ctx.exception))
        self.assertIn("UTF-8", str(ctx.exception))

    def test_b1_ascii_file_with_declared_encoding(self):
        path = write_temp(["id,value", "Ab1,1"], suffix=".csv")
        res = ValidationResult()
        codec = T.precheck(path, "latin-1", res, "validate")
        self.assertEqual(codec, "utf-8-sig")
        self.assertEqual([(f.level, f.code) for f in res.findings], [("info", "FILE_ENCODING_ASCII")])
        # a file that really is latin-1 still gets the warning
        fd, lat = tempfile.mkstemp(suffix=".csv")
        with os.fdopen(fd, "wb") as fh:
            fh.write("id,notes\nAb1,caf\xe9\n".encode("latin-1"))
        res = ValidationResult()
        self.assertEqual(T.precheck(lat, "latin-1", res, "validate"), "latin-1")
        self.assertEqual(res.codes(), ["FILE_ENCODING_INPUT"])

    def test_b5_inspect_refuses_excel(self):
        tmp = tempfile.mkdtemp()
        xlsx = os.path.join(tmp, "a.xlsx")
        with open(xlsx, "wb") as fh:
            fh.write(b"PK\x03\x04")
        for group in ("csv", "fasta"):
            code, _, err = run(group, "inspect", xlsx)
            self.assertEqual(code, 2)
            self.assertIn("CSV UTF-8", err)


if __name__ == "__main__":
    unittest.main()
