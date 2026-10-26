"""common/textfile.py: size, encoding, BOM and extension pre-checks."""
import os
import tempfile
import unittest

from abd_prep.common import textfile as T
from abd_prep.common.result import InputError, ValidationResult


def write_bytes(data: bytes, suffix: str = ".fasta") -> str:
    fd, path = tempfile.mkstemp(suffix=suffix)
    with os.fdopen(fd, "wb") as fh:
        fh.write(data)
    return path


class EncodingTest(unittest.TestCase):
    def test_utf8_plain_ok(self):
        res = ValidationResult()
        codec = T.precheck(write_bytes(b">h\nEVQL\n"), "utf-8", res)
        self.assertEqual(codec, "utf-8-sig")
        self.assertEqual(res.codes(), [])

    def test_utf8_bom_is_info_only(self):
        res = ValidationResult()
        T.precheck(write_bytes(b"\xef\xbb\xbfid,value\n"), None, res, "clean")
        self.assertEqual(res.codes(), ["FILE_BOM"])
        self.assertIn("removed", res.findings[0].message)
        res2 = ValidationResult()
        T.precheck(write_bytes(b"\xef\xbb\xbfid,value\n"), None, res2, "validate")
        self.assertIn("clean removes it", res2.findings[0].message)

    def test_not_utf8_without_encoding_is_error_and_stops(self):
        res = ValidationResult()
        codec = T.precheck(write_bytes("클론,값\n".encode("cp949"), ".csv"), "utf-8", res)
        self.assertIsNone(codec)
        self.assertEqual(res.codes(), ["FILE_NOT_UTF8"])
        self.assertIn("--encoding", res.findings[0].message)

    def test_declared_encoding_validate_warns_clean_informs(self):
        path = write_bytes("클론,값\n".encode("cp949"), ".csv")
        res = ValidationResult()
        self.assertEqual(T.precheck(path, "cp949", res, "validate"), "cp949")
        self.assertEqual(res.codes(), ["FILE_ENCODING_INPUT"])
        self.assertEqual(res.findings[0].level, "warning")
        res = ValidationResult()
        self.assertEqual(T.precheck(path, "cp949", res, "clean"), "cp949")
        self.assertEqual(res.codes(), ["FILE_ENCODING_CONVERTED"])
        self.assertEqual(res.findings[0].level, "info")

    def test_wrong_declared_encoding_is_error(self):
        res = ValidationResult()
        # 0xff 0xfe is not valid in cp949 (lead byte 0xff is never used)
        self.assertIsNone(T.precheck(write_bytes(b"id\n\xff\xfe\n", ".csv"), "cp949", res))
        self.assertEqual(res.codes(), ["FILE_ENCODING"])

    def test_unknown_encoding_is_input_error(self):
        with self.assertRaises(InputError):
            T.normalize_encoding("klingon-8")

    def test_utf8_aliases(self):
        for name in ("utf-8", "UTF8", "utf_8", "utf-8-sig"):
            self.assertEqual(T.normalize_encoding(name), "utf-8-sig")
        self.assertEqual(T.normalize_encoding("latin-1"), "latin-1")


class SizeTest(unittest.TestCase):
    def test_too_large(self):
        path = write_bytes(b"x")
        with open(path, "wb") as fh:
            fh.truncate(T.ABD_MAX_BYTES)       # sparse file, no real disk use
        res = ValidationResult()
        T.precheck(path, None, res)
        self.assertIn("FILE_TOO_LARGE", res.codes())
        os.remove(path)


class ExtensionTest(unittest.TestCase):
    def test_usual_extensions_silent(self):
        for ext in T.FASTA_EXTENSIONS:
            res = ValidationResult()
            T.check_extension("seed" + ext, T.FASTA_EXTENSIONS, "FASTA", res)
            self.assertEqual(res.codes(), [])
        for ext in T.CSV_EXTENSIONS:
            res = ValidationResult()
            T.check_extension("a" + ext, T.CSV_EXTENSIONS, "CSV", res)
            self.assertEqual(res.codes(), [])

    def test_unusual_extension_warns_but_continues(self):
        res = ValidationResult()
        T.check_extension("seed.seq", T.FASTA_EXTENSIONS, "FASTA", res)
        self.assertEqual([(f.level, f.code) for f in res.findings], [("warning", "FILE_EXTENSION")])

    def test_excel_is_input_error_with_save_as_hint(self):
        for ext in (".xlsx", ".xls", ".XLSM"):
            with self.assertRaises(InputError) as ctx:
                T.check_extension("assay" + ext, T.CSV_EXTENSIONS, "CSV", ValidationResult())
            self.assertIn("CSV UTF-8", str(ctx.exception))


if __name__ == "__main__":
    unittest.main()
