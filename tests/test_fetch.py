import json
import os
import tempfile
import unittest

from abd_prep.common.result import InputError
from abd_prep.pdb import fetch as F
from tests.helpers import fixture


def fake_opener(table):
    """``{url: (status, body)}`` -> opener; unknown URLs answer 404."""
    calls = []

    def opener(url):
        calls.append(url)
        return table.get(url, (404, b""))
    opener.calls = calls
    return opener


ENTRY_1N8Z = {
    "rcsb_accession_info": {"initial_release_date": "2003-02-18T00:00:00.000+00:00",
                            "revision_date": "2024-10-30T00:00:00.000+00:00"},
    "rcsb_entry_info": {"resolution_combined": [2.52], "experimental_method": "X-ray",
                        "deposited_model_count": 1},
    "struct": {"title": "Crystal structure of extracellular domain of human HER2 complexed with Herceptin Fab"},
    "exptl": [{"method": "X-RAY DIFFRACTION"}],
    "audit_author": [{"name": "should not be kept"}],
}


class PdbIdTest(unittest.TestCase):
    def test_valid_ids_upper_cased(self):
        self.assertEqual(F.normalize_pdb_id("1n8z"), "1N8Z")
        self.assertEqual(F.normalize_pdb_id(" 6B9Z "), "6B9Z")

    def test_invalid_ids_rejected(self):
        for bad in ("", "0ABC", "1N8", "1N8ZZ", "../x", "1N8Z/", "1N8Z.pdb", "A1N8", "1-8Z"):
            with self.subTest(bad=bad):
                with self.assertRaises(InputError):
                    F.normalize_pdb_id(bad)

    def test_urls_are_fixed_hosts(self):
        self.assertEqual(F.pdb_url("1n8z"), "https://files.rcsb.org/download/1N8Z.pdb")
        self.assertEqual(F.entry_url("1n8z"), "https://data.rcsb.org/rest/v1/core/entry/1N8Z")


class EntryMetaTest(unittest.TestCase):
    def test_summarize_keeps_only_known_fields(self):
        meta = F.summarize_entry(ENTRY_1N8Z, "1N8Z")
        self.assertEqual(meta.revision_date, "2024-10-30")
        self.assertEqual(meta.initial_release_date, "2003-02-18")
        self.assertEqual(meta.resolution, 2.52)
        self.assertEqual(meta.experimental_method, "X-ray")
        self.assertEqual(meta.model_count, 1)
        self.assertIn("Herceptin", meta.title)
        d = meta.to_dict()
        self.assertNotIn("audit_author", json.dumps(d))
        self.assertEqual(F.EntryMeta.from_dict(d), meta)

    def test_summarize_tolerates_missing_fields(self):
        meta = F.summarize_entry({}, "1ABC")
        self.assertEqual(meta.pdb_id, "1ABC")
        self.assertIsNone(meta.resolution)
        self.assertEqual(meta.revision_date, "")
        meta = F.summarize_entry({"exptl": [{"method": "SOLUTION NMR"}], "rcsb_entry_info": {}}, "1ABC")
        self.assertEqual(meta.experimental_method, "SOLUTION NMR")


class FetchPdbTest(unittest.TestCase):
    def setUp(self):
        self.cache = tempfile.mkdtemp()
        with open(fixture("1N8Z_mini.pdb"), "rb") as fh:
            self.body = fh.read()
        self.table = {
            F.pdb_url("1N8Z"): (200, self.body),
            F.entry_url("1N8Z"): (200, json.dumps(ENTRY_1N8Z).encode()),
            F.entry_url("4V9D"): (200, b"{}"),            # entry exists but no .pdb -> mmCIF only
        }

    def test_download_then_cache(self):
        op = fake_opener(self.table)
        r = F.fetch("1n8z", cache=self.cache, opener=op)
        self.assertFalse(r.from_cache)
        self.assertEqual(r.path, os.path.join(self.cache, "1N8Z.pdb"))
        self.assertEqual(r.bytes, len(self.body))
        self.assertEqual(r.meta.revision_date, "2024-10-30")
        self.assertTrue(os.path.isfile(os.path.join(self.cache, "1N8Z.entry.json")))
        self.assertFalse(os.path.exists(r.path + ".part"))
        n_calls = len(op.calls)
        r2 = F.fetch("1N8Z", cache=self.cache, opener=op)
        self.assertTrue(r2.from_cache)
        self.assertEqual(len(op.calls), n_calls)                 # no network on a cache hit
        r3 = F.fetch("1N8Z", cache=self.cache, opener=op, refresh=True)
        self.assertFalse(r3.from_cache)
        self.assertGreater(len(op.calls), n_calls)

    def test_metadata_failure_is_not_fatal(self):
        op = fake_opener({F.pdb_url("1N8Z"): (200, self.body)})
        r = F.fetch("1N8Z", cache=self.cache, opener=op)
        self.assertIsNone(r.meta)
        self.assertIn("404", r.meta_error)

    def test_mmcif_only_entry(self):
        with self.assertRaises(F.FetchError) as cm:
            F.fetch_pdb("4V9D", cache=self.cache, opener=fake_opener(self.table))
        self.assertIn("mmCIF only", str(cm.exception))
        self.assertEqual(os.listdir(self.cache), [])

    def test_unknown_entry(self):
        with self.assertRaises(F.FetchError) as cm:
            F.fetch_pdb("9ZZ9", cache=self.cache, opener=fake_opener(self.table))
        self.assertIn("does not exist", str(cm.exception))

    def test_server_error_and_bad_body(self):
        with self.assertRaises(F.FetchError) as cm:
            F.fetch_pdb("1N8Z", cache=self.cache, opener=fake_opener({F.pdb_url("1N8Z"): (503, b"")}))
        self.assertIn("HTTP 503", str(cm.exception))
        with self.assertRaises(F.FetchError) as cm:
            F.fetch_pdb("1N8Z", cache=self.cache, opener=fake_opener({F.pdb_url("1N8Z"): (200, b"<html>oops")}))
        self.assertIn("does not look like a PDB", str(cm.exception))
        self.assertEqual(os.listdir(self.cache), [])          # nothing half-written

    def test_cache_dir_env_and_override(self):
        old = os.environ.get(F.CACHE_ENV)
        try:
            os.environ[F.CACHE_ENV] = os.path.join(self.cache, "env")
            self.assertEqual(F.cache_dir(), os.path.join(self.cache, "env"))
            self.assertEqual(F.cache_dir(os.path.join(self.cache, "arg")), os.path.join(self.cache, "arg"))
        finally:
            if old is None:
                del os.environ[F.CACHE_ENV]
            else:
                os.environ[F.CACHE_ENV] = old

    def test_copy_to(self):
        r = F.fetch("1N8Z", cache=self.cache, opener=fake_opener(self.table), with_meta=False)
        out = os.path.join(self.cache, "out")
        dest = F.copy_to(r, out)
        self.assertEqual(dest, os.path.join(out, "1N8Z.pdb"))
        with open(dest, "rb") as fh:
            self.assertEqual(fh.read(), self.body)


if __name__ == "__main__":
    unittest.main()
