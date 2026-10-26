"""fasta/fastafile.py + fasta/validate.py: F13 normalization, F15-F17 rules, writer."""
import io
import os
import tempfile
import unittest

from abd_prep.common.result import InputError, ValidationResult
from abd_prep.fasta import validate as V
from abd_prep.fasta.fastafile import MAX_SEQUENCE_LENGTH, Record, read_fasta, write_fasta
from tests.helpers import fixture
from tests.make_seq_fixtures import LEADER, SCFV, VH, VHH, VL


def parse(text: str):
    res = ValidationResult()
    fasta = read_fasta(io.StringIO(text), res, path="mem")
    return fasta, res


def validate(text: str, **kw):
    fasta, res = parse(text)
    V.validate_records(fasta.records, V.FastaOptions(**kw), res)
    return fasta, res


class ParseNormalizeTest(unittest.TestCase):
    def test_plain_two_records(self):
        fasta, res = parse(">heavy desc here\n" + VH + "\n>light\n" + VL + "\n")
        self.assertEqual(res.codes(), [])
        self.assertEqual([r.id for r in fasta.records], ["heavy", "light"])
        self.assertEqual(fasta.records[0].description, "desc here")
        self.assertEqual(fasta.records[0].sequence, VH)
        self.assertEqual(fasta.records[0].lineno, 1)
        self.assertEqual(fasta.records[1].lineno, 3)

    def test_crlf_bom_whitespace_digits_gaps_lowercase(self):
        text = ("\ufeff>h\r\n" + "   1 " + VH[:10].lower() + " 11 " + VH[10:20] + "\r\n" + VH[20:30] + "--" + VH[30:] + "\r\n")
        fasta, res = parse(text)
        self.assertEqual(fasta.records[0].sequence, VH)
        codes = res.codes()
        for code in ("FASTA_CRLF", "FASTA_WHITESPACE_REMOVED", "FASTA_GAPS_REMOVED", "FASTA_LOWERCASE"):
            self.assertIn(code, codes)
        self.assertTrue(all(f.level == "info" for f in res.findings))

    def test_trailing_stop_removed_internal_stop_error(self):
        fasta, res = parse(">h\n" + VH + "*\n")
        self.assertEqual(fasta.records[0].sequence, VH)
        self.assertEqual(res.codes(), ["FASTA_TRAILING_STOP"])
        fasta, res = parse(">h\n" + VH[:60] + "*" + VH[60:] + "\n")
        self.assertIn("FASTA_INTERNAL_STOP", res.codes())
        self.assertTrue(res.has_errors)
        self.assertIn("61", res.errors[0].message)
        self.assertIn("*", fasta.records[0].sequence)       # never silently deleted

    def test_no_header_no_records_empty(self):
        _, res = parse(VH + "\n")
        self.assertIn("FASTA_NO_HEADER", res.codes())
        self.assertIn("FASTA_NO_RECORDS", res.codes())
        _, res = parse("")
        self.assertEqual(res.codes(), ["FASTA_NO_RECORDS"])
        fasta, res = parse(">h\n>l\n" + VL + "\n")
        self.assertIn("FASTA_EMPTY_SEQUENCE", res.codes())
        self.assertEqual(res.errors[0].where, "h")

    def test_empty_header_and_duplicate_id(self):
        _, res = parse(">\n" + VH + "\n>x\n" + VL + "\n>x\n" + VL + "\n")
        self.assertIn("FASTA_EMPTY_HEADER", res.codes())
        dup = [f for f in res.findings if f.code == "FASTA_DUPLICATE_ID"]
        self.assertEqual(len(dup), 1)
        self.assertEqual(dup[0].level, "warning")

    def test_too_long(self):
        _, res = parse(">h\n" + "A" * (MAX_SEQUENCE_LENGTH + 1) + "\n")
        self.assertIn("FASTA_TOO_LONG", res.codes())

    def test_header_id_is_first_token_only(self):
        fasta, _ = parse(">1N8Z_2|Chains B|Herceptin Fab heavy chain|Homo sapiens\n" + VH + "\n")
        self.assertEqual(fasta.records[0].id, "1N8Z_2|Chains")
        self.assertEqual(fasta.records[0].description, "B|Herceptin Fab heavy chain|Homo sapiens")


class AlphabetTest(unittest.TestCase):
    def test_standard_only_ok(self):
        _, res = validate(">heavy\n" + VH + "\n>light\n" + VL + "\n", framework="igg")
        self.assertEqual(res.codes(), [])

    def test_ambiguous_letters_error_with_positions(self):
        bad = VH[:9] + "X" + VH[10:49] + "B" + VH[50:]
        _, res = validate(">heavy\n" + bad + "\n", framework="nanobody")
        errs = [f for f in res.errors if f.code == "FASTA_BAD_ALPHABET"]
        self.assertEqual(len(errs), 1)
        self.assertIn("X@10", errs[0].message)
        self.assertIn("B@50", errs[0].message)

    def test_explicit_replacement_is_info(self):
        bad = VH[:9] + "X" + VH[10:]
        fasta, res = validate(">heavy\n" + bad + "\n", framework="nanobody", replace_ambiguous={"X": "G"})
        self.assertEqual(fasta.records[0].sequence, VH)
        self.assertEqual(res.codes(), ["FASTA_AMBIGUOUS_REPLACED"])

    def test_parse_replacements(self):
        self.assertEqual(V.parse_replacements("X=G, b=n"), {"X": "G", "B": "N"})
        self.assertEqual(V.parse_replacements(None), {})
        for bad in ("X", "XX=G", "X=Z", "A=G"):
            with self.assertRaises(InputError):
                V.parse_replacements(bad)

    def test_nucleic_detected_only_when_long(self):
        _, res = validate(">heavy\n" + "ACGT" * 20 + "\n", framework="nanobody", heuristics=False)
        self.assertEqual([f.code for f in res.errors], ["FASTA_NUCLEIC"])
        _, res = validate(">heavy\nACGTACG\n", framework="nanobody", heuristics=False)
        self.assertNotIn("FASTA_NUCLEIC", res.codes())         # short peptide of A/C/G/T is a protein

    def test_internal_stop_not_double_reported_as_alphabet(self):
        _, res = validate(">heavy\n" + VH[:60] + "*" + VH[60:] + "\n", framework="nanobody")
        self.assertIn("FASTA_INTERNAL_STOP", res.codes())
        self.assertNotIn("FASTA_BAD_ALPHABET", res.codes())


class FrameworkTest(unittest.TestCase):
    def test_igg_ok(self):
        fasta, res = validate(">heavy\n" + VH + "\n>light\n" + VL + "\n", framework="igg")
        self.assertEqual([r.role for r in fasta.records], ["H", "L"])
        self.assertEqual(res.codes(), [])

    def test_igg_wrong_count(self):
        _, res = validate(">heavy\n" + VH + "\n", framework="igg")
        self.assertEqual(res.codes(), ["FASTA_FRAMEWORK_MISMATCH"])
        _, res = validate(">heavy\n" + VH + "\n>light\n" + VL + "\n>VHH_x\n" + VHH + "\n", framework="igg")
        self.assertEqual([f.code for f in res.errors], ["FASTA_FRAMEWORK_MISMATCH"])

    def test_igg_unknown_role_is_error_with_hint(self):
        _, res = validate(">chain1\n" + VH + "\n>chain2\n" + VL + "\n", framework="igg")
        codes = [f.code for f in res.errors]
        self.assertEqual(codes.count("FASTA_ROLE_UNKNOWN"), 2)
        self.assertIn("--heavy-id chain1", res.errors[0].message)

    def test_no_framework_unknown_role_is_warning(self):
        _, res = validate(">chain1\n" + VH + "\n")
        self.assertEqual([(f.level, f.code) for f in res.findings if f.code != "FASTA_LENGTH_RANGE"],
                         [("warning", "FASTA_ROLE_UNKNOWN"), ("info", "FASTA_FRAMEWORK_UNSPECIFIED")])

    def test_nanobody_single_record_role_h(self):
        fasta, res = validate(">VHH_x\n" + VHH + "\n", framework="nanobody")
        self.assertEqual(fasta.records[0].role, "H")
        self.assertEqual(fasta.records[0].role_source, "framework")
        self.assertEqual(res.codes(), [])
        _, res = validate(">light_x\n" + VL + "\n", framework="nanobody")
        self.assertIn("FASTA_ROLE_CONFLICT", res.codes())

    def test_scfv_single_record_with_linker(self):
        fasta, res = validate(">scFv_1\n" + SCFV + "\n", framework="scfv")
        self.assertEqual(fasta.records[0].role, "HL")
        self.assertEqual(res.codes(), ["FASTA_SCFV_LINKER"])
        self.assertIn("121-135", res.findings[0].message)
        _, res = validate(">scFv_1\n" + VH + VL + "\n", framework="scfv")
        self.assertIn("no common", [f.message for f in res.findings if f.code == "FASTA_SCFV_LINKER"][0])

    def test_scfv_two_records_mismatch(self):
        _, res = validate(">a\n" + SCFV + "\n>b\n" + SCFV + "\n", framework="scfv")
        self.assertIn("FASTA_FRAMEWORK_MISMATCH", res.codes())


class HeuristicsTest(unittest.TestCase):
    def test_length_ranges(self):
        _, res = validate(">heavy\n" + VH[:80] + "\n>light\n" + VL + "\n", framework="igg")
        self.assertEqual([f.code for f in res.warnings], ["FASTA_LENGTH_RANGE"])
        self.assertIn("shorter", res.warnings[0].message)
        _, res = validate(">heavy\n" + VH + VL + "\n>light\n" + VL + "\n", framework="igg")
        self.assertIn("longer", res.warnings[0].message)

    def test_vhh_range_used_for_nanobody(self):
        self.assertEqual(V.LENGTH_RANGE["VHH"], (100, 140))
        _, res = validate(">x\n" + VHH[:99] + "\n", framework="nanobody")
        self.assertIn("VHH", res.warnings[0].message)

    def test_scfv_range(self):
        _, res = validate(">x\n" + SCFV + "\n", framework="scfv")
        self.assertEqual(res.warnings, [])
        _, res = validate(">x\n" + SCFV + SCFV + "\n", framework="scfv")
        self.assertEqual([f.code for f in res.warnings], ["FASTA_LENGTH_RANGE"])

    def test_signal_peptide(self):
        self.assertTrue(V.looks_like_signal_peptide(LEADER + VH))
        self.assertFalse(V.looks_like_signal_peptide(VH))
        self.assertFalse(V.looks_like_signal_peptide("M" + VL))          # M but no hydrophobic core
        _, res = validate(">heavy\n" + LEADER + VH + "\n>light\n" + VL + "\n", framework="igg")
        self.assertEqual([f.code for f in res.warnings], ["FASTA_SIGNAL_PEPTIDE"])       # 139 res: still in range
        _, res = validate(">heavy\n" + LEADER + VH + VL[:30] + "\n>light\n" + VL + "\n", framework="igg")
        self.assertEqual(sorted(f.code for f in res.warnings), ["FASTA_LENGTH_RANGE", "FASTA_SIGNAL_PEPTIDE"])

    def test_duplicate_sequence_warning(self):
        _, res = validate(">heavy\n" + VH + "\n>light\n" + VH + "\n", framework="igg")
        self.assertIn("FASTA_DUPLICATE_SEQUENCE", res.codes())

    def test_heuristics_off(self):
        _, res = validate(">heavy\n" + VH[:80] + "\n>light\n" + VL + "\n", framework="igg", heuristics=False)
        self.assertEqual(res.codes(), ["FASTA_HEURISTICS_OFF"])


class WriteTest(unittest.TestCase):
    def test_headers_rewritten_and_wrap(self):
        fasta, _ = validate(">VL_1 kappa some text\n" + VL + "\n>heavy\n" + VH + "\n", framework="igg")
        buf = io.StringIO()
        write_fasta(fasta.records, buf, wrap=60)
        lines = buf.getvalue().splitlines()
        self.assertEqual(lines[0], ">L VL_1 kappa some text")
        self.assertEqual(lines[1], VL[:60])
        self.assertEqual(lines[2], VL[60:])
        self.assertEqual(lines[3], ">H heavy")
        self.assertTrue(all(len(ln) <= 60 for ln in lines[1:]))

    def test_round_trip_file(self):
        fasta, _ = validate(">heavy\n" + VH + "\n>light\n" + VL + "\n", framework="igg")
        path = os.path.join(tempfile.mkdtemp(), "out.fasta")
        write_fasta(fasta.records, path)
        with open(path, "rb") as fh:
            data = fh.read()
        self.assertNotIn(b"\r", data)
        self.assertTrue(data.startswith(b">H heavy\n"))
        again, res = parse(data.decode("utf-8"))
        self.assertEqual([r.sequence for r in again.records], [VH, VL])
        self.assertEqual(res.codes(), [])

    def test_scfv_header_constant(self):
        rec = Record("scFv_1", "", SCFV, role="HL")
        buf = io.StringIO()
        write_fasta([rec], buf)
        self.assertTrue(buf.getvalue().startswith(">scFv scFv_1\n"))


class FixtureTest(unittest.TestCase):
    """The committed fixtures behave as plan §6 #5/#6 describe."""

    def load(self, name, **kw):
        res = ValidationResult()
        fasta = read_fasta(fixture(name), res)
        V.validate_records(fasta.records, V.FastaOptions(**kw), res)
        return fasta, res

    def test_seed_igg(self):
        fasta, res = self.load("seed_igg.fasta", framework="igg")
        self.assertEqual([r.role for r in fasta.records], ["H", "L"])
        self.assertEqual([r.sequence for r in fasta.records], [VH, VL])
        self.assertFalse(res.has_errors)
        self.assertFalse(res.has_warnings)
        for code in ("FASTA_CRLF", "FASTA_LOWERCASE", "FASTA_GAPS_REMOVED", "FASTA_TRAILING_STOP"):
            self.assertIn(code, res.codes())

    def test_seed_warn_fulllength(self):
        _, res = self.load("seed_warn_fulllength.fasta", framework="igg", heavy_id="heavy_full", light_id="light_full")
        self.assertFalse(res.has_errors)
        self.assertEqual(sorted(f.code for f in res.warnings), ["FASTA_LENGTH_RANGE", "FASTA_SIGNAL_PEPTIDE"])
        self.assertEqual(res.exit_code(strict=False), 0)
        self.assertEqual(res.exit_code(strict=True), 1)

    def test_bad_alphabet(self):
        _, res = self.load("seed_bad_alphabet.fasta", framework="igg")
        codes = [f.code for f in res.errors]
        self.assertEqual(codes.count("FASTA_BAD_ALPHABET"), 2)
        self.assertIn("FASTA_INTERNAL_STOP", codes)

    def test_bad_nucleic(self):
        _, res = self.load("seed_bad_nucleic.fasta", framework="igg")
        self.assertEqual([f.code for f in res.errors], ["FASTA_NUCLEIC"])

    def test_bad_count(self):
        _, res = self.load("seed_bad_count.fasta", framework="igg")
        self.assertEqual([f.code for f in res.errors], ["FASTA_FRAMEWORK_MISMATCH"])

    def test_nanobody_and_scfv_clean(self):
        fasta, res = self.load("seed_nanobody.fasta", framework="nanobody")
        self.assertEqual(fasta.records[0].sequence, VHH)
        self.assertEqual(res.codes(), [])
        fasta, res = self.load("seed_scfv.fasta", framework="scfv")
        self.assertEqual(fasta.records[0].sequence, SCFV)
        self.assertEqual(res.codes(), ["FASTA_SCFV_LINKER"])


if __name__ == "__main__":
    unittest.main()
