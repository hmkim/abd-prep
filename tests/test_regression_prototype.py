"""Regression against the prototype's four reference outputs (plan §6 #2).

Needs the original RCSB files in ``~/.cache/abd-prep/`` and the reference
outputs produced by ``data/clean_pdb.py`` (set ``ABD_PREP_REFERENCE_DIR``;
defaults to the workspace ``data/`` directory). Skipped otherwise, so the
suite stays network-free.
"""
import os
import tempfile
import unittest

from abd_prep.pdb.clean import CleanOptions, clean_file
from abd_prep.pdb.pdbfile import read_pdb
from abd_prep.pdb.ranges import ChainSelection, parse_rename
from abd_prep.pdb.validate import validate_atoms
from tests.helpers import CACHE

REFERENCE_DIR = os.environ.get(
    "ABD_PREP_REFERENCE_DIR",
    os.path.normpath(os.path.join(os.path.dirname(__file__), "..", "..", "..", "data")))

CASES = [
    ("HER2_domainIV_target.pdb", "6LBX.pdb", ["B:530-601"], ""),
    ("trastuzumab_Fv_HLT.pdb", "6B9Z.pdb", ["A:1-107", "B:1-120"], "A=L,B=H"),
    ("alt_HER2_domainIV_from1N8Z.pdb", "1N8Z.pdb", ["C:489-607"], ""),
    ("alt_trastuzumab_Fv_from1N8Z.pdb", "1N8Z.pdb", ["A:1-107", "B:1-120"], "A=L,B=H"),
]


def _available():
    return all(os.path.isfile(os.path.join(CACHE, src)) for _, src, _, _ in CASES) and \
        all(os.path.isfile(os.path.join(REFERENCE_DIR, ref)) for ref, _, _, _ in CASES)


def coord_records(path):
    with open(path) as fh:
        return [l.rstrip("\n") for l in fh if l.startswith(("ATOM", "HETATM"))]


@unittest.skipUnless(_available(), "RCSB originals / prototype reference outputs not available")
class PrototypeRegressionTest(unittest.TestCase):
    def test_atom_records_identical_and_ter_regenerated(self):
        tmp = tempfile.mkdtemp()
        for ref_name, src, keep, rename in CASES:
            with self.subTest(ref=ref_name):
                out = os.path.join(tmp, ref_name)
                outcome = clean_file(os.path.join(CACHE, src), out,
                                     CleanOptions(selection=ChainSelection.parse(keep), rename=parse_rename(rename)))
                self.assertFalse(outcome.result.has_errors, outcome.result.format_table())
                self.assertEqual(coord_records(os.path.join(REFERENCE_DIR, ref_name)), coord_records(out))
                # one TER per output chain, each naming the last residue written
                chains = list(outcome.details["chains_out"])
                self.assertEqual(len(outcome.ters), len(chains))
                for ter, ch in zip(outcome.ters, chains):
                    self.assertEqual(ter[21], ch)
                    self.assertEqual(ter[22:26].strip(), outcome.details["chains_out"][ch]["last"])
                res, _ = validate_atoms(read_pdb(out))
                self.assertFalse(res.has_errors, res.format_table())

    def test_prototype_outputs_fail_ter_validation(self):
        """The P1 defect in the reference files is what validate must catch."""
        for ref_name, _, _, _ in CASES:
            res, _ = validate_atoms(read_pdb(os.path.join(REFERENCE_DIR, ref_name)))
            self.assertTrue(res.has_errors, ref_name)
            self.assertTrue({"PDB_TER_MISMATCH", "PDB_TER_EXTRA", "PDB_TER_MISPLACED"} & set(res.codes()), ref_name)


if __name__ == "__main__":
    unittest.main()
