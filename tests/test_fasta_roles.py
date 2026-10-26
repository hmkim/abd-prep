"""fasta/roles.py: header-token inference, explicit IDs, framework roles, header constant."""
import unittest

from abd_prep.common.result import InputError, ValidationResult
from abd_prep.fasta import roles as R
from abd_prep.fasta.fastafile import Record
from abd_prep.pdb.validate import to_fasta


def rec(ident, desc="", seq="EVQL"):
    return Record(ident, desc, seq)


class InferTest(unittest.TestCase):
    def test_heavy_tokens(self):
        for header in ("heavy", "VH_1", "Ab1 heavy chain", "H", "hc", "IGHV3-23*01 germline", "VHH_nb",
                       "1N8Z_2|Chains B|Herceptin Fab heavy chain|Homo sapiens", "my-nanobody"):
            self.assertEqual(R.infer_role(header)[0], "H", header)

    def test_light_tokens(self):
        for header in ("light", "VL_1 kappa", "L", "lc", "IGKV1-39", "IGLV2-14", "lambda chain", "vk",
                       "1N8Z_1|Chains A|Herceptin Fab light chain|Homo sapiens"):
            self.assertEqual(R.infer_role(header)[0], "L", header)

    def test_unknown_not_guessed(self):
        for header in ("chain1", "Ab001", "seq_A", "trastuzumab", "H1", "heavyweight", "1N8Z_1|Chains A"):
            self.assertIsNone(R.infer_role(header)[0], header)

    def test_ambiguous(self):
        role, tokens = R.infer_role("VH_VL_fusion heavy and light")
        self.assertEqual(role, "?")
        self.assertIn("vh", tokens)
        self.assertIn("vl", tokens)

    def test_scfv_token_wins_over_vh_vl(self):
        self.assertEqual(R.infer_role("scFv_1 VH-(G4S)3-VL construct")[0], "HL")
        self.assertEqual(R.infer_role("my_scfv")[0], "HL")

    def test_tokenizer_splits_on_non_alphanumerics(self):
        self.assertEqual(R.header_tokens("A|B-c_d e.f"), ["a", "b", "c", "d", "e", "f"])


class AssignTest(unittest.TestCase):
    def assign(self, records, **kw):
        res = ValidationResult()
        R.assign_roles(records, res, **kw)
        return res

    def test_header_inference_igg(self):
        recs = [rec("heavy"), rec("VL_1", "kappa")]
        res = self.assign(recs, framework="igg")
        self.assertEqual([(r.role, r.role_source) for r in recs], [("H", "header"), ("L", "header")])
        self.assertEqual(res.codes(), [])

    def test_explicit_ids_win_and_note_override(self):
        recs = [rec("a", "light"), rec("b")]
        res = self.assign(recs, framework="igg", heavy_id="a", light_id="b")
        self.assertEqual([(r.role, r.role_source) for r in recs], [("H", "option"), ("L", "option")])
        self.assertEqual(res.codes(), ["FASTA_ROLE_OVERRIDDEN"])

    def test_explicit_id_not_found_or_duplicate(self):
        with self.assertRaises(InputError):
            self.assign([rec("a"), rec("b")], framework="igg", heavy_id="zzz")
        with self.assertRaises(InputError):
            self.assign([rec("a"), rec("a")], framework="igg", heavy_id="a")
        with self.assertRaises(InputError):
            self.assign([rec("a"), rec("b")], framework="igg", heavy_id="a", light_id="a")

    def test_unknown_role_error_under_igg_warning_otherwise(self):
        recs = [rec("x"), rec("y")]
        res = self.assign(recs, framework="igg")
        self.assertEqual([f.level for f in res.findings], ["error", "error"])
        self.assertEqual([r.role for r in recs], [None, None])
        res = self.assign([rec("x")])
        self.assertEqual([f.level for f in res.findings], ["warning"])

    def test_ambiguous_is_error(self):
        res = self.assign([rec("VH_VL")], framework="igg")
        self.assertEqual([f.code for f in res.errors], ["FASTA_ROLE_AMBIGUOUS"])

    def test_scfv_and_nanobody_roles_from_framework(self):
        recs = [rec("anything")]
        self.assign(recs, framework="scfv")
        self.assertEqual((recs[0].role, recs[0].role_source), ("HL", "framework"))
        recs = [rec("anything")]
        self.assign(recs, framework="nanobody")
        self.assertEqual((recs[0].role, recs[0].role_source), ("H", "framework"))

    def test_nanobody_light_token_conflict(self):
        recs = [rec("VL_1")]
        res = self.assign(recs, framework="nanobody")
        self.assertEqual([f.code for f in res.errors], ["FASTA_ROLE_CONFLICT"])
        self.assertEqual(recs[0].role, "H")              # still written as H so the file stays usable

    def test_unknown_framework(self):
        with self.assertRaises(InputError):
            self.assign([rec("a")], framework="fab")


class HeaderConstantTest(unittest.TestCase):
    def test_header_ids(self):
        self.assertEqual(R.fasta_header("H"), ">H")
        self.assertEqual(R.fasta_header("L", "desc"), ">L desc")
        self.assertEqual(R.fasta_header("HL", "x"), ">scFv x")
        self.assertEqual(R.fasta_header("C", "target"), ">C target")     # unknown chain written as-is

    def test_pdb_inspect_fasta_uses_same_constant(self):
        info = {"path": "fv.pdb", "chains": {"H": {"first": "1", "last": "3", "residues": 3, "sequence": "EVQ"},
                                             "L": {"first": "1", "last": "3", "residues": 3, "sequence": "DIQ"}}}
        text = to_fasta(info)
        self.assertTrue(text.startswith(">" + R.ABD_HEADER_ID["H"] + " "))
        self.assertIn("\n>" + R.ABD_HEADER_ID["L"] + " ", text)


if __name__ == "__main__":
    unittest.main()
