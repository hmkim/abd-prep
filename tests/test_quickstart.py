"""``abd-prep quickstart`` -- the one-step simple mode (every file for a first ABD run)."""
import csv
import json
import os
import shutil
import tempfile
import unittest

from abd_prep import quickstart as Q
from abd_prep.common.result import InputError
from abd_prep.fasta.fastafile import read_fasta
from abd_prep.common.result import ValidationResult
from abd_prep.pdb import fetch as F
from abd_prep.pdb import presets as P
from abd_prep.pdb import split as S
from abd_prep.pdb.pdbfile import read_pdb
from abd_prep.pdb.validate import inspect_file
from tests.helpers import CACHE, fixture
from tests.make_seq_fixtures import VH, VL
from tests.test_cli import run
from tests.test_fetch import fake_opener

ENTRY_MINI = {"rcsb_accession_info": {"revision_date": "2024-10-30T00:00:00.000+00:00"},
              "rcsb_entry_info": {"resolution_combined": [3.1]}, "struct": {"title": "mini"}}

EXPECTED_FILES = ["README.txt", "antibody_fv.pdb", "antibody_seed.fasta", "assay_example.csv",
                  "quickstart.report.json", "target.pdb"]


def mini_opener():
    with open(fixture("1N8Z_mini.pdb"), "rb") as fh:
        body = fh.read()
    return fake_opener({F.pdb_url("1N8Z"): (200, body), F.entry_url("1N8Z"): (200, json.dumps(ENTRY_MINI).encode())})


def seed_cache_with_mini(cache):
    """Pre-populate a cache dir so the CLI path needs no network."""
    os.makedirs(cache, exist_ok=True)
    shutil.copyfile(fixture("1N8Z_mini.pdb"), os.path.join(cache, "1N8Z.pdb"))
    meta = F.EntryMeta(pdb_id="1N8Z", title="mini", revision_date="2024-10-30", resolution=3.1)
    with open(os.path.join(cache, "1N8Z.entry.json"), "w", encoding="utf-8") as fh:
        json.dump(meta.to_dict(), fh)


class QuickstartMockedTest(unittest.TestCase):
    """The mini fixture holds A 1-5, B 1-5, C 489-493, so sizes are tiny but every step runs."""

    def setUp(self):
        self.cache = tempfile.mkdtemp()
        self.out = os.path.join(tempfile.mkdtemp(), "kit")
        self.preset = P.get_preset("her2-trastuzumab")

    def _opts(self, **kw):
        return Q.QuickstartOptions(preset=self.preset, out_dir=self.out, cache=self.cache, **kw)

    def test_every_file_is_written_and_consistent(self):
        outcome = Q.run_quickstart(self._opts(), opener=mini_opener(), command="abd-prep quickstart --out kit")
        self.assertEqual(sorted(os.listdir(self.out)), EXPECTED_FILES)
        self.assertEqual(list(outcome.files), ["target", "fv", "seed", "csv", "readme"])
        self.assertEqual(outcome.framework, "igg")
        self.assertFalse(outcome.result.has_errors, outcome.result.format_table())

        # structures: same content the split path produces, under the fixed names
        fv = read_pdb(outcome.files["fv"])
        self.assertEqual(fv.chains(), ["L", "H"])
        self.assertEqual(len(fv.ters), 2)
        self.assertEqual(read_pdb(outcome.files["target"]).chains(), ["C"])

        # seed FASTA: heavy first, sequences identical to the Fv chains, ABD headers
        fasta = read_fasta(outcome.files["seed"], ValidationResult())
        self.assertEqual([r.id for r in fasta.records], ["H", "L"])
        chains = inspect_file(outcome.files["fv"])["chains"]
        self.assertEqual(fasta.records[0].sequence, chains["H"]["sequence"])
        self.assertEqual(fasta.records[1].sequence, chains["L"]["sequence"])
        self.assertIn("antibody_fv.pdb chain H 1-5", fasta.records[0].description)

        # example CSV: header, one placeholder row that matches the seed, labelled as such
        with open(outcome.files["csv"], encoding="utf-8", newline="") as fh:
            rows = list(csv.reader(fh))
        self.assertEqual(rows[0], ["id", "heavy", "light", "value", "unit", "target", "assay", "notes"])
        self.assertEqual(len(rows), 2)
        self.assertEqual(rows[1][0], "her2-trastuzumab-seed")
        self.assertEqual(rows[1][1], fasta.records[0].sequence)
        self.assertEqual(rows[1][2], fasta.records[1].sequence)
        self.assertEqual(rows[1][5], "1N8Z C:489-607")
        self.assertIn("placeholder", rows[1][7])

        # README names every file and the placeholder nature of the CSV
        with open(outcome.files["readme"], encoding="utf-8") as fh:
            readme = fh.read()
        for name in EXPECTED_FILES:
            self.assertIn(name, readme)
        self.assertIn("contains no measurement", readme)
        self.assertIn("--rename C=T", readme)
        self.assertIn("revised 2024-10-30", readme)

        # findings keep their group prefixes and the cross-checks ran
        wheres = [f.where for f in outcome.result.findings]
        self.assertTrue(any(w.startswith("target/") for w in wheres), wheres)
        self.assertTrue(any(w.startswith("fv/") for w in wheres), wheres)
        self.assertTrue(any(w.startswith("seed/") for w in wheres), wheres)
        self.assertTrue(any(w.startswith("csv/") for w in wheres), wheres)
        codes = outcome.result.codes()
        self.assertEqual(codes.count("FASTA_PDB_IDENTITY"), 2)
        self.assertIn("CSV_CROSSREF_MATCHED", codes)
        self.assertNotIn("CSV_CROSSREF_UNMATCHED", codes)
        self.assertIn("PDB_RESOLUTION", codes)                  # 3.1 A entry metadata

        # report: inputs = the source entry, outputs = the five files, basename mode leaks no path
        self.assertEqual(outcome.report_path, os.path.join(self.out, "quickstart.report.json"))
        with open(outcome.report_path) as fh:
            rep = json.load(fh)
        self.assertEqual([i["path"] for i in rep["inputs"]], ["1N8Z.pdb"])
        self.assertEqual([(o["kind"], o["path"]) for o in rep["outputs"]],
                         [("pdb", "target.pdb"), ("pdb", "antibody_fv.pdb"), ("fasta", "antibody_seed.fasta"),
                          ("csv", "assay_example.csv"), ("text", "README.txt")])
        self.assertEqual(sorted(rep["details"]), ["csv", "files", "framework", "preset", "seed", "split"])
        self.assertEqual(rep["details"]["preset"]["name"], "her2-trastuzumab")
        self.assertEqual(rep["details"]["split"]["sources"][0]["pdb_id"], "1N8Z")
        self.assertNotIn("preset", rep["details"]["split"])
        self.assertEqual(rep["details"]["csv"]["against"]["matched"], 1)
        self.assertEqual(rep["details"]["seed"]["against"]["chain_map"], {"H": "H", "L": "L"})
        self.assertNotIn(self.out, json.dumps(rep))

    def test_refuses_to_overwrite_unless_forced(self):
        Q.run_quickstart(self._opts(), opener=mini_opener())
        with self.assertRaises(InputError) as cm:
            Q.run_quickstart(self._opts(), opener=mini_opener())
        self.assertIn("--force", str(cm.exception))
        self.assertIn("target.pdb", str(cm.exception))
        outcome = Q.run_quickstart(self._opts(force=True), opener=mini_opener())
        self.assertEqual(sorted(os.listdir(self.out)), EXPECTED_FILES)
        self.assertFalse(outcome.result.has_errors)

    def test_no_csv_and_no_report(self):
        outcome = Q.run_quickstart(self._opts(with_csv=False, report_path=""), opener=mini_opener())
        self.assertEqual(sorted(os.listdir(self.out)), ["README.txt", "antibody_fv.pdb", "antibody_seed.fasta", "target.pdb"])
        self.assertIsNone(outcome.report_path)
        self.assertNotIn("csv", outcome.files)
        self.assertNotIn("csv", outcome.report.details)
        self.assertFalse(any(f.where.startswith("csv/") for f in outcome.result.findings))
        with open(outcome.files["readme"], encoding="utf-8") as fh:
            self.assertNotIn("assay_example.csv", fh.read())

    def test_custom_report_path_and_full_paths(self):
        rp = os.path.join(tempfile.mkdtemp(), "r.json")
        outcome = Q.run_quickstart(self._opts(report_path=rp, path_mode="full"), opener=mini_opener())
        self.assertEqual(outcome.report_path, rp)
        with open(rp) as fh:
            rep = json.load(fh)
        self.assertTrue(all(os.path.isabs(o["path"]) for o in rep["outputs"]))

    def test_preset_without_antibody_or_target_is_rejected(self):
        target_only = P.Preset("t", "target only", (P.Part("target", "1N8Z", "C:489-493"),))
        with self.assertRaises(InputError) as cm:
            Q.run_quickstart(Q.QuickstartOptions(preset=target_only, out_dir=self.out, cache=self.cache),
                             opener=mini_opener())
        self.assertIn("no heavy chain", str(cm.exception))
        fv_only = P.Preset("f", "fv only", (P.Part("heavy", "1N8Z", "B:1-5", "H"), P.Part("light", "1N8Z", "A:1-5", "L")))
        with self.assertRaises(InputError) as cm:
            Q.run_quickstart(Q.QuickstartOptions(preset=fv_only, out_dir=self.out, cache=self.cache),
                             opener=mini_opener())
        self.assertIn("no target part", str(cm.exception))
        self.assertFalse(os.path.exists(self.out))                  # nothing written before the check

    def test_nanobody_style_preset_writes_one_record_and_sequence_column(self):
        vhh = P.Preset("nb", "heavy only", (P.Part("target", "1N8Z", "C:489-493"), P.Part("heavy", "1N8Z", "B:1-5", "H")))
        outcome = Q.run_quickstart(Q.QuickstartOptions(preset=vhh, out_dir=self.out, cache=self.cache), opener=mini_opener())
        self.assertEqual(outcome.framework, "nanobody")
        fasta = read_fasta(outcome.files["seed"], ValidationResult())
        self.assertEqual([r.id for r in fasta.records], ["H"])
        with open(outcome.files["csv"], encoding="utf-8", newline="") as fh:
            rows = list(csv.reader(fh))
        self.assertEqual(rows[0][:2], ["id", "sequence"])
        self.assertEqual(rows[1][1], fasta.records[0].sequence)
        self.assertFalse(outcome.result.has_errors, outcome.result.format_table())

    def test_split_filename_override_is_validated(self):
        opts = S.SplitOptions.from_preset(self.preset, self.out, cache=self.cache, filenames={"fv": "a/b.pdb"})
        with self.assertRaises(InputError):
            S.run_split(opts, opener=mini_opener())


class QuickstartCliTest(unittest.TestCase):
    def setUp(self):
        self.cache = tempfile.mkdtemp()
        seed_cache_with_mini(self.cache)
        self.out = os.path.join(tempfile.mkdtemp(), "kit")

    def test_text_output_lists_files_and_next_step(self):
        code, out, err = run("quickstart", "--out", self.out, "--cache", self.cache)
        self.assertEqual(code, 0, err)
        self.assertIn("[quickstart] her2-trastuzumab  ->  " + self.out, out)
        self.assertIn("summary: 0 error(s)", out)
        self.assertIn("report: " + os.path.join(self.out, "quickstart.report.json"), out)
        self.assertIn("Files in " + self.out, out)
        for name in ("target.pdb", "antibody_fv.pdb", "antibody_seed.fasta", "assay_example.csv", "README.txt"):
            self.assertIn(name, out)
        self.assertIn("Next: read README.txt", out)
        self.assertEqual(out.count("what each file is, where it came from"), 1)     # README line printed once
        self.assertEqual(sorted(os.listdir(self.out)), EXPECTED_FILES)

    def test_json_output_is_the_report(self):
        code, out, err = run("quickstart", "--out", self.out, "--cache", self.cache, "--json", "--no-csv")
        self.assertEqual(code, 0, err)
        rep = json.loads(out)
        self.assertEqual(rep["command"], f"abd-prep quickstart --out kit --cache {os.path.basename(self.cache)} --json --no-csv")
        self.assertEqual([o["path"] for o in rep["outputs"]], ["target.pdb", "antibody_fv.pdb", "antibody_seed.fasta", "README.txt"])

    def test_second_run_without_force_is_exit_2(self):
        self.assertEqual(run("quickstart", "--out", self.out, "--cache", self.cache)[0], 0)
        code, out, err = run("quickstart", "--out", self.out, "--cache", self.cache)
        self.assertEqual(code, 2)
        self.assertIn("--force", err)
        code, out, err = run("quickstart", "--out", self.out, "--cache", self.cache, "--force")
        self.assertEqual(code, 0, err)

    def test_strict_fails_on_the_resolution_warning(self):
        code, out, err = run("quickstart", "--out", self.out, "--cache", self.cache, "--strict")
        self.assertEqual(code, 1)                                   # 3.1 A metadata -> PDB_RESOLUTION warning

    def test_unknown_preset_is_a_usage_error(self):
        code, out, err = run("quickstart", "--preset", "nope", "--out", self.out)
        self.assertEqual(code, 2)


@unittest.skipUnless(os.path.isfile(os.path.join(CACHE, "1N8Z.pdb")), "RCSB original 1N8Z not cached")
class QuickstartRegressionTest(unittest.TestCase):
    """Offline against the real 1N8Z: the seed is trastuzumab and every check passes."""

    def test_default_preset_from_cache(self):
        out = os.path.join(tempfile.mkdtemp(), "kit")
        opts = Q.QuickstartOptions(preset=P.get_preset(P.DEFAULT_PRESET), out_dir=out, cache=CACHE)
        outcome = Q.run_quickstart(opts, opener=fake_opener({}))      # cache hit only: any download = failure
        self.assertFalse(outcome.result.has_errors, outcome.result.format_table())
        self.assertEqual([f.code for f in outcome.result.warnings], ["PDB_RESOLUTION"])
        fasta = read_fasta(outcome.files["seed"], ValidationResult())
        self.assertEqual([(r.id, r.sequence) for r in fasta.records], [("H", VH), ("L", VL)])
        self.assertEqual(outcome.report.details["csv"]["against"]["by_sequence"], 1)
        self.assertEqual(outcome.report.details["split"]["parts"]["target"]["clean"]["chains_out"]["C"]["residues"], 109)


if __name__ == "__main__":
    unittest.main()
