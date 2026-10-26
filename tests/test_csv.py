"""csvdata: schema.py, csvfile.py, validate.py (F19-F21) and the CSV fixtures (plan §6 #8)."""
import csv
import os
import tempfile
import unittest

from abd_prep.common.result import InputError, ValidationResult
from abd_prep.csvdata import csvfile as CF
from abd_prep.csvdata import schema as S
from abd_prep.csvdata import validate as CV
from tests.helpers import fixture, write_temp
from tests.make_seq_fixtures import VH, VL


def run(path, out=None, codec="utf-8-sig", **kw):
    res = ValidationResult()
    details = CV.process_csv(path, codec, CV.CsvOptions(**kw), res, out_path=out)
    return res, details


def tmp_csv(lines):
    return write_temp(lines, suffix=".csv")


class SchemaTest(unittest.TestCase):
    def test_parse_map(self):
        self.assertEqual(S.parse_map(["id=Clone,sequence=AA_seq", "value=KD"]),
                         {"id": "Clone", "sequence": "AA_seq", "value": "KD"})
        self.assertEqual(S.parse_map(["notes=Comment, free text"]), {"notes": "Comment, free text"})   # single pair keeps commas
        self.assertEqual(S.parse_map(None), {})
        for bad in (["id"], ["clone=id"], ["id=A,id=B"], ["id=A,value=A"], ["id="]):
            with self.assertRaises(InputError):
                S.parse_map(bad)

    def test_resolve_columns_explicit_auto_and_unmapped(self):
        header = ["Clone", "sequence", "KD", "Unit ", "Note"]
        resolved, unmapped = S.resolve_columns(header, {"id": "Clone", "value": "KD"})
        self.assertEqual(resolved, {"id": 0, "sequence": 1, "value": 2, "unit": 3})   # 'sequence'/'Unit' auto-mapped
        self.assertEqual(unmapped, [4])
        with self.assertRaises(InputError) as ctx:
            S.resolve_columns(header, {"id": "clone"})                                   # case-sensitive source names
        self.assertIn("'Clone'", str(ctx.exception))

    def test_auto_map_does_not_steal_explicitly_mapped_column(self):
        resolved, _ = S.resolve_columns(["id", "value"], {"notes": "id"})
        self.assertEqual(resolved, {"notes": 0, "value": 1})

    def test_missing_tokens(self):
        for tok in ("", "  ", "NA", "n/a", "NaN", "None", "null", "-", "ND"):
            self.assertTrue(S.is_missing(tok), tok)
        for tok in ("0", "0.0", "none-such", "NAC"):
            self.assertFalse(S.is_missing(tok), tok)

    def test_classify_number(self):
        self.assertEqual(S.classify_number(" 2.3 "), ("number", 2.3))
        self.assertEqual(S.classify_number("1e-9")[0], "number")
        self.assertEqual(S.classify_number("-.5")[0], "number")
        self.assertEqual(S.classify_number("1,2"), ("decimal_comma", None))
        self.assertEqual(S.classify_number("<0.1"), ("censored", None))
        self.assertEqual(S.classify_number(">= 1000"), ("censored", None))
        self.assertEqual(S.classify_number("NA"), ("missing", None))
        self.assertEqual(S.classify_number("weak"), ("text", None))
        self.assertEqual(S.classify_number("1.2.3"), ("text", None))


class DelimiterTest(unittest.TestCase):
    def test_parse_delimiter(self):
        self.assertEqual(CF.parse_delimiter("tab"), "\t")
        self.assertEqual(CF.parse_delimiter("\\t"), "\t")
        self.assertEqual(CF.parse_delimiter(";"), ";")
        self.assertIsNone(CF.parse_delimiter(None))
        with self.assertRaises(InputError):
            CF.parse_delimiter(";;")

    def test_sniff(self):
        self.assertEqual(CF.sniff_delimiter("a;b;c\n1;2;3\n4;5;6\n"), ";")
        self.assertEqual(CF.sniff_delimiter("a\tb\n1\t2\n"), "\t")
        self.assertEqual(CF.sniff_delimiter("id,value\nAb1,1\n"), ",")
        self.assertIsNone(CF.sniff_delimiter(""))
        # Sniffer gives up on a single header-only line; the header fallback decides when one candidate is present
        self.assertEqual(CF.sniff_delimiter("id;value\n"), ";")
        self.assertIsNone(CF.sniff_delimiter("id\n"))

    def test_guessed_delimiter_warning_and_explicit_info(self):
        path = tmp_csv(["id", "Ab1", "Ab2"])
        res, _ = run(path)
        self.assertIn("CSV_DELIMITER_GUESSED", res.codes())
        path = tmp_csv(["id;value", "Ab1;1", "Ab2;2"])
        res, details = run(path, delimiter=";")
        self.assertEqual([f.code for f in res.findings if f.code.startswith("CSV_DELIM")], ["CSV_DELIMITER"])
        self.assertEqual(details["delimiter"], ";")

    def test_wrong_explicit_delimiter_surfaces_as_missing_id(self):
        path = tmp_csv(["id;value", "Ab1;1"])
        res, _ = run(path, delimiter=",")
        self.assertIn("CSV_ID_MISSING", res.codes())


class HeaderRowTest(unittest.TestCase):
    def test_empty_file_and_blank_header(self):
        res, _ = run(tmp_csv([""]))
        self.assertEqual(res.codes(), ["CSV_NO_HEADER"])
        res, _ = run(tmp_csv(["", "", "id,value", "Ab1,1"]))              # leading blank lines skipped
        self.assertNotIn("CSV_NO_HEADER", res.codes())

    def test_header_problems(self):
        res, _ = run(tmp_csv(["id,,value,value", "a,b,c,d"]))
        self.assertIn("CSV_HEADER_EMPTY_NAME", res.codes())
        self.assertIn("CSV_HEADER_DUPLICATE", res.codes())

    def test_row_length_error_lists_line_and_excludes_row(self):
        out = os.path.join(tempfile.mkdtemp(), "o.csv")
        res, details = run(tmp_csv(["id,value", "Ab1,1", "Ab2,2,extra", "", "Ab3,3"]), out=out)
        err = [f for f in res.errors if f.code == "CSV_ROW_LENGTH"][0]
        self.assertIn("line 3 (3 cells)", err.message)
        self.assertEqual(details["rows_read"], 2)
        self.assertEqual(details["rows_written"], 2)
        self.assertIn("CSV_BLANK_ROWS", res.codes())
        with open(out, newline="") as fh:
            self.assertEqual(list(csv.reader(fh)), [["id", "value"], ["Ab1", "1"], ["Ab3", "3"]])

    def test_quoted_multiline_field_keeps_line_numbers(self):
        path = tmp_csv(['id,notes,value', 'Ab1,"two', 'lines",1', 'Ab2,x,oops', 'Ab3,y,3,extra'])
        res, _ = run(path)
        err = [f for f in res.errors if f.code == "CSV_ROW_LENGTH"][0]
        self.assertIn("line 5", err.message)
        self.assertIn("CSV_VALUE_NOT_NUMERIC", res.codes())

    def test_header_only(self):
        res, _ = run(tmp_csv(["id,value"]))
        self.assertIn("CSV_NO_DATA_ROWS", res.codes())


class ColumnRulesTest(unittest.TestCase):
    def test_id_required_sequence_recommended(self):
        res, _ = run(tmp_csv(["Clone,KD", "Ab1,1"]))
        self.assertEqual([f.code for f in res.errors], ["CSV_ID_MISSING"])
        self.assertIn("CSV_SEQUENCE_MISSING", [f.code for f in res.warnings])
        res, _ = run(tmp_csv(["Clone,KD", "Ab1,1"]), mapping={"id": "Clone", "value": "KD"})
        self.assertFalse(res.has_errors)
        res, _ = run(tmp_csv(["id,heavy,light", "Ab1,EVQL,DIQM"]))
        self.assertNotIn("CSV_SEQUENCE_MISSING", res.codes())

    def test_value_rules(self):
        # with a ',' delimiter a decimal comma makes the row too long -> row-length error, not a value warning
        res, _ = run(tmp_csv(["id,value", "a,1.0", "c,1,2"]))
        self.assertIn("CSV_ROW_LENGTH", res.codes())
        path = tmp_csv(["id;value", "a;1.0", "b;<0.1", "c;1,2", "d;NA", "e;", "f;strong"])
        res, details = run(path, delimiter=";")
        codes = res.codes()
        self.assertIn("CSV_VALUE_NOT_NUMERIC", codes)
        self.assertIn("CSV_DECIMAL_COMMA", codes)
        self.assertIn("CSV_MISSING_VALUES", codes)
        self.assertEqual(details["missing"], {"value": 2})
        nn = [f for f in res.findings if f.code == "CSV_VALUE_NOT_NUMERIC"][0]
        self.assertIn("'<0.1'", nn.message)
        self.assertIn("'strong'", nn.message)
        self.assertTrue(all(f.level == "warning" for f in res.findings if f.code.startswith("CSV_")
                            and f.code not in ("CSV_DELIMITER",)))

    def test_units_and_duplicates(self):
        res, details = run(tmp_csv(["id,value,unit", "a,1,nM", "b,2,pM", "a,3,nM"]))
        self.assertIn("CSV_UNIT_MIXED", res.codes())
        self.assertIn("CSV_ID_DUPLICATE", res.codes())
        self.assertEqual(details["units"], {"nM": 2, "pM": 1})
        res, _ = run(tmp_csv(["id,value,unit", "a,1,nM", "b,2,nM"]))
        self.assertEqual(res.codes(), ["CSV_SEQUENCE_MISSING"])            # no sequence column: recommendation only

    def test_sequence_column_rules(self):
        res, _ = run(tmp_csv(["id,sequence", "a," + VH, "b," + VH[:30].lower(), "c,EVQLXZ", "d," + "ACGT" * 10]))
        self.assertEqual(sorted(f.code for f in res.errors), ["CSV_SEQUENCE_ALPHABET", "CSV_SEQUENCE_NUCLEIC"])
        self.assertIn("CSV_SEQUENCE_LOWERCASE", res.codes())
        bad = [f for f in res.errors if f.code == "CSV_SEQUENCE_ALPHABET"][0]
        self.assertIn("'XZ'", bad.message)

    def test_drop_missing_and_unmapped(self):
        out = os.path.join(tempfile.mkdtemp(), "o.csv")
        path = tmp_csv(["Clone,KD,Note", "a,1,x", "b,NA,y", "c,3,z"])
        res, details = run(path, out=out, mapping={"id": "Clone", "value": "KD"}, drop_missing=True, drop_unmapped=True)
        self.assertEqual(details["dropped"], 1)
        self.assertIn("CSV_ROWS_DROPPED", res.codes())
        with open(out, newline="") as fh:
            self.assertEqual(list(csv.reader(fh)), [["id", "value"], ["a", "1"], ["c", "3"]])


class CleanOutputTest(unittest.TestCase):
    def test_standard_columns_first_then_unmapped_values_unchanged(self):
        out = os.path.join(tempfile.mkdtemp(), "o.csv")
        path = tmp_csv(["Note;KD;Clone;Unit", "hello, world;1,2;Ab1;nM", "x;<0.1;Ab2;pM"])
        run(path, out=out, mapping={"id": "Clone", "value": "KD", "unit": "Unit"})
        with open(out, "rb") as fh:
            raw = fh.read()
        self.assertFalse(raw.startswith(b"\xef\xbb\xbf"))
        self.assertNotIn(b"\r", raw)
        with open(out, newline="") as fh:
            rows = list(csv.reader(fh))
        self.assertEqual(rows[0], ["id", "value", "unit", "Note"])
        self.assertEqual(rows[1], ["Ab1", "1,2", "nM", "hello, world"])      # values untouched, comma quoted
        self.assertIn(b'"hello, world"', raw)
        self.assertIn(b'"1,2"', raw)

    def test_no_output_when_encoding_fails_is_callers_job(self):
        # process_csv itself is only called after precheck succeeded; it never deletes an existing file
        out = os.path.join(tempfile.mkdtemp(), "o.csv")
        run(tmp_csv(["id", "a"]), out=out)
        self.assertTrue(os.path.exists(out))


class EncodingTest(unittest.TestCase):
    def test_cp949_runtime_generated(self):
        fd, path = tempfile.mkstemp(suffix=".csv")
        with os.fdopen(fd, "wb") as fh:
            fh.write("클론,서열,KD\nAb1,EVQL,1.0\n".encode("cp949"))
        out = os.path.join(tempfile.mkdtemp(), "o.csv")
        res, details = run(path, out=out, codec="cp949", mapping={"id": "클론", "sequence": "서열", "value": "KD"})
        self.assertFalse(res.has_errors)
        with open(out, encoding="utf-8", newline="") as fh:
            self.assertEqual(list(csv.reader(fh)), [["id", "sequence", "value"], ["Ab1", "EVQL", "1.0"]])

    def test_utf8_bom_header_is_clean(self):
        fd, path = tempfile.mkstemp(suffix=".csv")
        with os.fdopen(fd, "wb") as fh:
            fh.write(b"\xef\xbb\xbfid,value\nAb1,1\n")
        res, details = run(path)
        self.assertEqual(res.codes(), ["CSV_SEQUENCE_MISSING"])
        self.assertEqual(list(details["columns"].items()), [("id", "id"), ("value", "value")])


class InspectTest(unittest.TestCase):
    def test_inspect_fill_rates_and_mapping(self):
        res = ValidationResult()
        info = CV.inspect_csv(fixture("assay_messy_utf8bom.csv"), "utf-8-sig", None, res, {"id": "Clone", "value": "KD"})
        self.assertEqual(info["rows"], 8)
        cols = {c["name"]: c for c in info["columns"]}
        self.assertEqual(cols["Clone"]["standard"], "id")
        self.assertEqual(cols["KD"]["filled"], 6)
        self.assertEqual(cols["AA_seq"]["filled"], 7)
        text = CV.format_inspect(info)
        self.assertIn("Clone -> id: 8/8", text)


class FixtureTest(unittest.TestCase):
    """Plan §6 #8 on the committed fixtures."""

    def test_assay_ok_no_findings(self):
        res, details = run(fixture("assay_ok.csv"), mapping={"id": "id", "sequence": "sequence", "value": "value",
                                                              "unit": "unit"})
        self.assertEqual(res.codes(), [])
        self.assertEqual(details["rows_read"], 20)

    def test_messy_utf8bom(self):
        out = os.path.join(tempfile.mkdtemp(), "o.csv")
        res, details = run(fixture("assay_messy_utf8bom.csv"), out=out,
                           mapping={"id": "Clone", "sequence": "AA_seq", "value": "KD", "unit": "Unit"})
        self.assertFalse(res.has_errors)
        self.assertEqual(details["delimiter"], ";")
        for code in ("CSV_MISSING_VALUES", "CSV_VALUE_NOT_NUMERIC", "CSV_DECIMAL_COMMA", "CSV_UNIT_MIXED",
                     "CSV_ID_DUPLICATE", "CSV_SEQUENCE_LOWERCASE"):
            self.assertIn(code, res.codes(), code)
        with open(out, newline="", encoding="utf-8") as fh:
            rows = list(csv.reader(fh))
        self.assertEqual(rows[0], ["id", "sequence", "value", "unit", "Note"])
        self.assertEqual(len(rows), 9)

    def test_messy_cp949(self):
        res, details = run(fixture("assay_messy_cp949.csv"), codec="cp949",
                           mapping={"id": "클론", "sequence": "서열", "value": "KD"})
        self.assertEqual([f.code for f in res.errors], ["CSV_ROW_LENGTH"])
        self.assertIn("line 4 (6 cells)", res.errors[0].message)
        self.assertEqual(details["rows_read"], 3)


if __name__ == "__main__":
    unittest.main()
