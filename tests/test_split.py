"""Presets (F6) and ``pdb split`` (plan §3.1, DoD #1-#3)."""
import json
import os
import tempfile
import unittest

from abd_prep.common.result import InputError
from abd_prep.pdb import fetch as F
from abd_prep.pdb import presets as P
from abd_prep.pdb import split as S
from abd_prep.pdb.pdbfile import read_pdb
from abd_prep.pdb.validate import validate_atoms
from tests.helpers import CACHE, fixture
from tests.test_fetch import fake_opener

REFERENCE_DIR = os.environ.get(
    "ABD_PREP_REFERENCE_DIR",
    os.path.normpath(os.path.join(os.path.dirname(__file__), "..", "..", "..", "data")))

#: preset -> (group -> prototype reference file)
REFERENCES = {
    "her2-trastuzumab": {"target": "alt_HER2_domainIV_from1N8Z.pdb", "fv": "alt_trastuzumab_Fv_from1N8Z.pdb"},
    "her2-trastuzumab-6lbx-6b9z": {"target": "HER2_domainIV_target.pdb", "fv": "trastuzumab_Fv_HLT.pdb"},
}


def coord_records(path):
    with open(path) as fh:
        return [l.rstrip("\n") for l in fh if l.startswith(("ATOM", "HETATM"))]


class PresetDefinitionTest(unittest.TestCase):
    def test_presets_are_well_formed(self):
        self.assertIn(P.DEFAULT_PRESET, P.PRESETS)
        for name, preset in P.PRESETS.items():
            with self.subTest(preset=name):
                self.assertEqual(preset.name, name)
                P.validate_parts(preset.parts)
                roles = {p.role for p in preset.parts}
                self.assertEqual(roles, {"target", "heavy", "light"})
                groups = P.output_groups(preset.parts)
                self.assertEqual(set(groups), {"target", "fv"})
                fv_chains = {p.written_chain for p in groups["fv"]}
                self.assertEqual(fv_chains, {"H", "L"})
                for ch in preset.expect_residues:
                    self.assertIn(ch, {p.written_chain for p in preset.parts})
                for g in preset.expected_gaps:
                    self.assertIn(g[0], {p.written_chain for p in preset.parts})
                d = preset.to_dict()
                json.dumps(d)                                      # serialisable for report.json
                self.assertEqual(len(d["parts"]), 3)

    def test_default_preset_is_1n8z(self):
        preset = P.get_preset("her2-trastuzumab")
        self.assertEqual(preset.pdb_ids, ["1N8Z"])
        self.assertEqual(preset.part("target").selection, "C:489-607")
        self.assertIsNone(preset.part("target").out_chain)        # chain label decision (R1) still open
        self.assertEqual(preset.part("heavy").out_chain, "H")
        self.assertEqual(preset.part("light").out_chain, "L")
        self.assertIn(("C", 580, 591), preset.expected_gaps)

    def test_unknown_preset(self):
        with self.assertRaises(InputError):
            P.get_preset("nope")

    def test_part_rules(self):
        with self.assertRaises(InputError):
            P.validate_parts(())
        with self.assertRaises(InputError):               # two chains in one part
            P.validate_parts((P.Part("target", "1N8Z", "A:1-5,B:1-5"),))
        with self.assertRaises(InputError):               # heavy/light from different entries
            P.validate_parts((P.Part("heavy", "1N8Z", "B:1-120", "H"), P.Part("light", "6B9Z", "A:1-107", "L")))
        with self.assertRaises(InputError):               # written chain collision
            P.validate_parts((P.Part("heavy", "1N8Z", "B:1-120"), P.Part("light", "1N8Z", "A:1-107", "B")))
        with self.assertRaises(InputError):               # same source chain twice
            P.validate_parts((P.Part("heavy", "1N8Z", "B:1-120", "H"), P.Part("light", "1N8Z", "B:1-107", "L")))
        with self.assertRaises(InputError):               # role twice
            P.validate_parts((P.Part("target", "1N8Z", "C:1-5"), P.Part("target", "1N8Z", "A:1-5")))
        with self.assertRaises(InputError):               # bad range
            P.validate_parts((P.Part("target", "1N8Z", "C:9-1"),))
        self.assertEqual(P.output_filename("x", "fv"), "x_fv.pdb")


class SplitLocalInputTest(unittest.TestCase):
    """The mini fixture holds A 1-5, B 1-5, C 489-493 (+ 2 waters), so expectations are tiny."""

    def setUp(self):
        self.out = tempfile.mkdtemp()
        self.parts = (P.Part("target", "input", "C:489-493"),
                      P.Part("light", "input", "A:1-5", "L"),
                      P.Part("heavy", "input", "B:1-5", "H"))

    def test_files_report_and_chains(self):
        opts = S.SplitOptions(name="mini", parts=self.parts, out_dir=self.out, input_path=fixture("1N8Z_mini.pdb"),
                              expect_residues={"C": 5, "L": 5, "H": 5})
        outcome = S.run_split(opts, command="abd-prep pdb split ...")
        self.assertEqual(sorted(outcome.written), ["fv", "target"])
        self.assertEqual(os.path.basename(outcome.written["target"]), "mini_target.pdb")
        self.assertEqual(os.path.basename(outcome.written["fv"]), "mini_fv.pdb")
        self.assertEqual(outcome.report_path, os.path.join(self.out, "mini.report.json"))
        self.assertFalse(outcome.result.has_errors, outcome.result.format_table())
        fv = read_pdb(outcome.written["fv"])
        self.assertEqual(fv.chains(), ["L", "H"])
        self.assertEqual(len(fv.ters), 2)
        self.assertTrue(fv.has_end)
        self.assertFalse(any(a.altloc != " " for a in fv.atoms))       # fixture altlocs resolved
        target = read_pdb(outcome.written["target"])
        self.assertEqual(target.chains(), ["C"])
        # findings carry the part prefix
        wheres = {f.where for f in outcome.result.findings}
        self.assertTrue(any(w.startswith("fv/") for w in wheres), wheres)
        with open(outcome.report_path) as fh:
            rep = json.load(fh)
        self.assertEqual([i["kind"] for i in rep["inputs"]], ["pdb"])
        self.assertEqual([o["path"] for o in rep["outputs"]], ["mini_target.pdb", "mini_fv.pdb"])
        self.assertEqual(rep["details"]["name"], "mini")
        self.assertIsNone(rep["details"]["preset"])
        self.assertEqual(rep["details"]["parts"]["fv"]["roles"]["heavy"]["out_chain"], "H")
        self.assertEqual(rep["details"]["parts"]["target"]["clean"]["chains_out"]["C"]["residues"], 5)
        self.assertEqual(rep["details"]["sources"][0]["pdb_id"], None)
        self.assertNotIn(self.out, json.dumps(rep))                     # basename mode leaks no paths

    def test_residue_expectation_mismatch_is_reported(self):
        opts = S.SplitOptions(name="mini", parts=self.parts, out_dir=self.out, input_path=fixture("1N8Z_mini.pdb"),
                              expect_residues={"C": 109}, report_path="")
        outcome = S.run_split(opts)
        self.assertIsNone(outcome.report_path)
        codes = [f.code for f in outcome.result.findings if f.where.startswith("target/")]
        self.assertIn("PDB_RESIDUE_COUNT", codes)

    def test_target_only_and_heavy_only(self):
        opts = S.SplitOptions(name="t", parts=(self.parts[0],), out_dir=self.out,
                              input_path=fixture("1N8Z_mini.pdb"), report_path="")
        self.assertEqual(list(S.run_split(opts).written), ["target"])
        opts = S.SplitOptions(name="nb", parts=(self.parts[2],), out_dir=self.out,
                              input_path=fixture("1N8Z_mini.pdb"), report_path="")
        outcome = S.run_split(opts)
        self.assertEqual(list(outcome.written), ["fv"])
        self.assertEqual(read_pdb(outcome.written["fv"]).chains(), ["H"])

    def test_derive_name(self):
        self.assertEqual(S.derive_name(self.parts, "/x/y/my_structure.pdb"), "my_structure")
        self.assertEqual(S.derive_name((P.Part("target", "1N8Z", "C:1-5"),), None), "1n8z")
        self.assertEqual(S.derive_name((P.Part("target", "6LBX", "B:1-5"), P.Part("heavy", "6B9Z", "B:1-5")), None),
                         "6lbx-6b9z")

    def test_missing_input_and_bad_name(self):
        opts = S.SplitOptions(name="x", parts=self.parts, out_dir=self.out, input_path=None)
        with self.assertRaises(InputError):
            S.run_split(opts)
        opts = S.SplitOptions(name="../x", parts=self.parts, out_dir=self.out, input_path=fixture("1N8Z_mini.pdb"))
        with self.assertRaises(InputError):
            S.run_split(opts)


class SplitMockedFetchTest(unittest.TestCase):
    def test_preset_run_with_injected_download(self):
        cache = tempfile.mkdtemp()
        out = tempfile.mkdtemp()
        with open(fixture("1N8Z_mini.pdb"), "rb") as fh:
            body = fh.read()
        entry = {"rcsb_accession_info": {"revision_date": "2024-10-30T00:00:00.000+00:00"},
                 "rcsb_entry_info": {"resolution_combined": [3.1]}, "struct": {"title": "mini"}}
        op = fake_opener({F.pdb_url("1N8Z"): (200, body), F.entry_url("1N8Z"): (200, json.dumps(entry).encode())})
        preset = P.get_preset("her2-trastuzumab")
        opts = S.SplitOptions.from_preset(preset, out, cache=cache)
        outcome = S.run_split(opts, opener=op)
        codes = outcome.result.codes()
        self.assertIn("PDB_RESOLUTION", codes)                   # 3.1 A > 2.5 A recommendation
        self.assertIn("PDB_RESIDUE_COUNT", codes)                # mini fixture != 109/107/120
        self.assertEqual(sorted(outcome.written), ["fv", "target"])
        self.assertTrue(os.path.basename(outcome.written["target"]).startswith("her2-trastuzumab_"))
        with open(outcome.report_path) as fh:
            rep = json.load(fh)
        self.assertEqual(rep["details"]["preset"]["name"], "her2-trastuzumab")
        self.assertEqual(rep["details"]["sources"][0]["entry"]["revision_date"], "2024-10-30")
        self.assertEqual(rep["details"]["sources"][0]["pdb_id"], "1N8Z")


def _cache_and_references_available():
    ok = all(os.path.isfile(os.path.join(CACHE, f"{pid}.pdb")) for pid in ("1N8Z", "6LBX", "6B9Z"))
    return ok and all(os.path.isfile(os.path.join(REFERENCE_DIR, f))
                      for refs in REFERENCES.values() for f in refs.values())


@unittest.skipUnless(_cache_and_references_available(), "RCSB originals / prototype reference outputs not available")
class PresetRegressionTest(unittest.TestCase):
    """DoD #1-#3 offline: every preset reproduces the prototype's ATOM records and validates clean."""

    def test_presets_reproduce_prototype_outputs(self):
        out = tempfile.mkdtemp()
        for name, refs in REFERENCES.items():
            with self.subTest(preset=name):
                preset = P.get_preset(name)
                opts = S.SplitOptions.from_preset(preset, out, cache=CACHE)
                outcome = S.run_split(opts, opener=fake_opener({}))   # cache hit only: any download = failure
                self.assertFalse(outcome.result.has_errors, outcome.result.format_table())
                self.assertFalse([f for f in outcome.result.warnings if f.code != "PDB_RESOLUTION"],
                                 outcome.result.format_table())
                for group, ref in refs.items():
                    self.assertEqual(coord_records(os.path.join(REFERENCE_DIR, ref)),
                                     coord_records(outcome.written[group]), f"{name}/{group}")
                    res, _ = validate_atoms(read_pdb(outcome.written[group]))
                    self.assertFalse(res.has_errors, res.format_table())
                self.assertTrue(os.path.isfile(os.path.join(out, f"{name}.report.json")))
                if name == "her2-trastuzumab":
                    # the 1N8Z target gap must be reported as *expected* (info), never as a warning
                    gap = [f for f in outcome.result.findings if f.code == "PDB_GAP_CA_DISTANCE"]
                    self.assertEqual([f.level for f in gap], ["info"])
                    self.assertIn("C:580-591", gap[0].where)


if __name__ == "__main__":
    unittest.main()
