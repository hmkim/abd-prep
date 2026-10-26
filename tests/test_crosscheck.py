"""fasta/crosscheck.py (F18) and csvdata/crosscheck.py (F22)."""
import os
import tempfile
import unittest

from abd_prep.common.residues import AA1_TO_3
from abd_prep.common.result import InputError, ValidationResult
from abd_prep.csvdata import crosscheck as CX
from abd_prep.csvdata import validate as CV
from abd_prep.fasta import crosscheck as X
from abd_prep.fasta.fastafile import Record, read_fasta
from tests.helpers import CACHE, fixture, residue, write_temp
from tests.make_seq_fixtures import VH, VL


def chain_lines(chain, letters, numbering, serial_start=1):
    """ATOM lines for ``letters`` placed at the given (resseq, icode) numbering."""
    lines = []
    serial = serial_start
    for i, (aa, (resseq, icode)) in enumerate(zip(letters, numbering)):
        lines.extend(residue(serial, AA1_TO_3[aa], chain, resseq, x=i * 3.8, icode=icode))
        serial += 4
    return lines


def plain(n, start=1):
    return [(start + i, " ") for i in range(n)]


def rec(ident, seq, role):
    return Record(ident, "", seq, role=role)


class ChainResiduesTest(unittest.TestCase):
    def test_gap_aware_index_and_insertion_codes(self):
        numbering = [(1, " "), (2, " "), (5, " "), (6, " "), (6, "A"), (6, "B"), (7, " ")]
        path = write_temp(chain_lines("H", VH[:7], numbering) + ["END"])
        chains = X.pdb_chain_residues(path)
        cr = chains["H"]
        self.assertEqual(cr.letters, VH[:7])
        self.assertEqual(cr.index, [0, 1, 4, 5, 6, 7, 8])
        self.assertEqual(cr.labels, ["1", "2", "5", "6", "6A", "6B", "7"])
        self.assertEqual(cr.gaps, [("2", "5", 2, 4)])
        self.assertEqual(cr.span, 9)

    def test_hetatm_polymer_residue_included_ligand_skipped(self):
        lines = chain_lines("A", "EV", plain(2))
        lines += residue(9, "MSE", "A", 3, x=7.6, hetatm=True)        # selenomethionine in the chain
        lines += chain_lines("A", "L", [(4, " ")], serial_start=13)
        lines += residue(17, "HOH", "A", 101, x=50, names=("O",), hetatm=True)
        path = write_temp(lines + ["END"])
        cr = X.pdb_chain_residues(path)["A"]
        self.assertEqual(cr.letters, "EVML")
        self.assertEqual(cr.gaps, [])


class AlignTest(unittest.TestCase):
    def chain(self, letters, numbering, chain="H"):
        path = write_temp(chain_lines(chain, letters, numbering) + ["END"])
        return X.pdb_chain_residues(path)[chain]

    def test_identical(self):
        al = X.align(self.chain(VH[:30], plain(30)), VH[:30])
        self.assertEqual((al.offset, al.overlap, al.matches, al.mismatches), (0, 30, 30, []))
        self.assertEqual((al.extra_n, al.extra_c, al.uncovered), (0, 0, []))

    def test_pdb_gap_is_reported_not_mismatched(self):
        numbering = plain(10) + plain(10, start=14)                   # 11-13 missing
        al = X.align(self.chain(VH[:10] + VH[13:23], numbering), VH[:23])
        self.assertEqual(al.mismatches, [])
        self.assertEqual(al.gaps_in_range, [("10", "14", 3, 11)])

    def test_point_mutation(self):
        mutated = VH[:4] + ("A" if VH[4] != "A" else "G") + VH[5:30]
        al = X.align(self.chain(VH[:30], plain(30)), mutated)
        self.assertEqual(len(al.mismatches), 1)
        label, pdb_aa, fasta_pos, fasta_aa = al.mismatches[0]
        self.assertEqual((label, pdb_aa, fasta_pos), ("5", VH[4], 5))
        self.assertEqual(al.matches, 29)

    def test_fasta_with_leader_and_tail(self):
        al = X.align(self.chain(VH[:40], plain(40)), "MGWSLILL" + VH[:40] + "ASTK")
        self.assertEqual((al.offset, al.extra_n, al.extra_c, al.matches), (8, 8, 4, 40))

    def test_fasta_shorter_than_pdb(self):
        al = X.align(self.chain(VH[:40], plain(40)), VH[:25])
        self.assertEqual(al.matches, 25)
        self.assertEqual(al.uncovered, [str(i) for i in range(26, 41)])

    def test_kabat_numbering_starts_at_offset(self):
        numbering = plain(5, start=82)                                 # 82..86, no insertion codes
        al = X.align(self.chain(VH[:5], numbering), "ZZZ" + VH[:5] + "ZZ")
        self.assertEqual((al.offset, al.matches), (3, 5))


class CrosscheckTest(unittest.TestCase):
    def setUp(self):
        gap_numbering = plain(50) + plain(57, start=64)                # H: 1-50, 64-120 (51-63 missing)
        lines = chain_lines("H", VH[:50] + VH[63:], gap_numbering)
        lines += chain_lines("L", VL, plain(len(VL)), serial_start=5000)
        lines.append("END")
        self.fv_path = write_temp(lines)
        self.raw_path = write_temp(chain_lines("B", VH, plain(len(VH)))
                                   + chain_lines("A", VL, plain(len(VL)), serial_start=5000) + ["END"])

    def run_check(self, records, path, chain_map=None):
        res = ValidationResult()
        details = X.crosscheck(records, path, X.parse_chain_map(chain_map), res)
        return res, details

    def test_identical_with_gap(self):
        res, details = self.run_check([rec("h", VH, "H"), rec("l", VL, "L")], self.fv_path)
        self.assertFalse(res.has_errors)
        self.assertFalse(res.has_warnings)
        codes = res.codes()
        self.assertEqual(codes.count("FASTA_PDB_IDENTITY"), 2)
        self.assertEqual(codes.count("FASTA_PDB_GAP"), 1)
        gap = [f for f in res.findings if f.code == "FASTA_PDB_GAP"][0]
        self.assertIn("13 missing", gap.message)
        self.assertIn("51-63", gap.message)
        self.assertEqual(details["chain_map"], {"H": "H", "L": "L"})

    def test_mismatch_is_error_with_positions(self):
        mutated = VL[:20] + ("A" if VL[20] != "A" else "G") + VL[21:]
        res, _ = self.run_check([rec("h", VH, "H"), rec("l", mutated, "L")], self.fv_path)
        errs = [f for f in res.errors]
        self.assertEqual([f.code for f in errs], ["FASTA_PDB_MISMATCH"])
        self.assertIn("PDB 21", errs[0].message)
        self.assertEqual(errs[0].where, "L:L<-l")

    def test_unrelated_sequence_no_correspondence(self):
        res, _ = self.run_check([rec("h", "ACDEFGHIKLMNPQRSTVWY" * 6, "H"), rec("l", VL, "L")], self.fv_path)
        self.assertEqual([f.code for f in res.errors], ["FASTA_PDB_NO_CORRESPONDENCE"])

    def test_missing_role_is_warning(self):
        res, _ = self.run_check([rec("h", VH, "H")], self.fv_path)
        self.assertEqual([f.code for f in res.warnings], ["FASTA_PDB_ROLE_MISSING"])

    def test_scfv_record_compared_to_both_chains(self):
        res, details = self.run_check([rec("s", VH + "GGGGSGGGGSGGGGS" + VL, "HL")], self.fv_path)
        self.assertFalse(res.has_errors)
        self.assertEqual(res.codes().count("FASTA_PDB_IDENTITY"), 2)
        self.assertEqual(details["chains"]["L:L<-s"]["extra_n"], 120 + 15)

    def test_raw_chain_ids_need_chain_map(self):
        with self.assertRaises(InputError) as ctx:
            self.run_check([rec("h", VH, "H"), rec("l", VL, "L")], self.raw_path)
        self.assertIn("--chain-map", str(ctx.exception))
        res, details = self.run_check([rec("h", VH, "H"), rec("l", VL, "L")], self.raw_path, "H=B,L=A")
        self.assertFalse(res.has_errors)
        self.assertEqual(details["chain_map"], {"H": "B", "L": "A"})
        with self.assertRaises(InputError):
            self.run_check([rec("h", VH, "H")], self.raw_path, "H=Z")

    def test_parse_chain_map(self):
        self.assertEqual(X.parse_chain_map("h=B, L=A"), {"H": "B", "L": "A"})
        self.assertEqual(X.parse_chain_map(None), {})
        for bad in ("H", "T=A", "H=AB", "H=A,H=B"):
            with self.assertRaises(InputError):
                X.parse_chain_map(bad)


@unittest.skipUnless(os.path.isfile(os.path.join(CACHE, "1N8Z.pdb")), "1N8Z not in the local cache")
class RealDataTest(unittest.TestCase):
    """Plan §6 #7: seed_igg.fasta vs the her2-trastuzumab preset Fv (both 1N8Z-derived)."""

    def test_seed_igg_against_preset_fv_and_raw_entry(self):
        from abd_prep.pdb import presets as P
        from abd_prep.pdb import split as S
        out = tempfile.mkdtemp()
        opts = S.SplitOptions.from_preset(P.get_preset("her2-trastuzumab"), out, report_path="")
        outcome = S.run_split(opts)
        fv_pdb = outcome.written["fv"]
        res = ValidationResult()
        fasta = read_fasta(fixture("seed_igg.fasta"), res)
        from abd_prep.fasta.roles import assign_roles
        assign_roles(fasta.records, res, "igg")
        X.crosscheck(fasta.records, fv_pdb, {}, res)
        self.assertEqual([f.code for f in res.errors], [])
        self.assertEqual([f.code for f in res.warnings], [])
        self.assertEqual(res.codes().count("FASTA_PDB_IDENTITY"), 2)
        self.assertNotIn("FASTA_PDB_GAP", res.codes())          # 1N8Z Fv chains A/B have no missing residues
        with self.assertRaises(InputError):
            X.crosscheck(fasta.records, os.path.join(CACHE, "1N8Z.pdb"), {}, ValidationResult())
        res2 = ValidationResult()
        X.crosscheck(fasta.records, os.path.join(CACHE, "1N8Z.pdb"), {"H": "B", "L": "A"}, res2)
        self.assertEqual([f.code for f in res2.errors], [])
        self.assertEqual([f.code for f in res2.warnings], ["FASTA_PDB_UNCOVERED", "FASTA_PDB_UNCOVERED"])


class CsvCrossRefTest(unittest.TestCase):
    def setUp(self):
        self.ref = CX.load_reference(fixture("seed_igg.fasta"))

    def test_reference_from_fasta(self):
        self.assertEqual(self.ref.kind, "fasta")
        self.assertEqual(self.ref.ids, {"heavy", "VL_1"})
        self.assertEqual(self.ref.sequences, {VH, VL})

    def test_reference_from_pdb(self):
        path = write_temp(chain_lines("H", VH[:10], plain(10)) + ["END"])
        os.rename(path, path + ".pdb")
        ref = CX.load_reference(path + ".pdb")
        self.assertEqual((ref.kind, ref.ids, ref.sequences), ("pdb", {"H"}, {VH[:10]}))

    def test_missing_reference_file(self):
        with self.assertRaises(InputError):
            CX.load_reference("/nonexistent/ref.fasta")

    def test_match_by_id_sequence_and_pair(self):
        hook = CX.CrossRef(self.ref)
        hook(2, {"id": "heavy", "sequence": "AAAA"})                     # id only
        hook(3, {"id": "x1", "sequence": VH.lower()})                    # sequence (case-insensitive)
        hook(4, {"id": "x2", "heavy": VH, "light": VL})                  # heavy+light pair
        hook(5, {"id": "x3", "heavy": VH, "light": "AAAA"})              # half a pair -> unmatched
        hook(6, {"id": "x4"})                                            # nothing
        res = ValidationResult()
        details = hook.report(res)
        self.assertEqual((details["matched"], details["by_id"], details["by_sequence"]), (3, 1, 2))
        self.assertEqual([f.code for f in res.warnings], ["CSV_CROSSREF_UNMATCHED"])
        self.assertIn("'x3'", res.warnings[0].message)
        self.assertIn("'x4'", res.warnings[0].message)

    def test_through_process_csv(self):
        path = write_temp(["id,sequence,value", f"heavy,{VH},1.0", f"VL_1,{VL},2.0", f"Ab9,{VH[:50]},3"], suffix=".csv")
        hook = CX.CrossRef(self.ref)
        res = ValidationResult()
        CV.process_csv(path, "utf-8-sig", CV.CsvOptions(), res, row_hook=hook)
        details = hook.report(res)
        self.assertEqual((details["rows"], details["matched"]), (3, 2))


if __name__ == "__main__":
    unittest.main()
