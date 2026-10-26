import unittest

from abd_prep.pdb.ranges import ChainSelection, RangeSyntaxError, ResidueRange, parse_rename


class ResidueRangeTest(unittest.TestCase):
    def test_plain_range(self):
        r = ResidueRange.parse("1-120")
        self.assertTrue(r.contains((1, " ")))
        self.assertTrue(r.contains((120, " ")))
        self.assertFalse(r.contains((121, " ")))
        self.assertFalse(r.contains((0, " ")))

    def test_upper_bound_without_code_includes_insertions(self):
        r = ResidueRange.parse("1-82")
        for code in (" ", "A", "B", "C"):
            self.assertTrue(r.contains((82, code)), code)
        self.assertFalse(r.contains((83, " ")))

    def test_upper_bound_with_code_is_exact(self):
        r = ResidueRange.parse("1-82A")
        self.assertTrue(r.contains((82, " ")))
        self.assertTrue(r.contains((82, "A")))
        self.assertFalse(r.contains((82, "B")))

    def test_lower_bound_code(self):
        r = ResidueRange.parse("82B-83")
        self.assertFalse(r.contains((82, " ")))
        self.assertFalse(r.contains((82, "A")))
        self.assertTrue(r.contains((82, "B")))
        self.assertTrue(r.contains((82, "C")))
        self.assertTrue(r.contains((83, " ")))

    def test_single_residue(self):
        r = ResidueRange.parse("100")
        self.assertTrue(r.contains((100, " ")))
        self.assertTrue(r.contains((100, "A")))   # all insertions of 100
        self.assertFalse(r.contains((101, " ")))
        self.assertEqual(str(r), "100")

    def test_negative_numbers(self):
        r = ResidueRange.parse("-3-10")
        self.assertTrue(r.contains((-3, " ")))
        self.assertTrue(r.contains((0, " ")))
        self.assertTrue(r.contains((10, " ")))

    def test_reversed_range_rejected(self):
        with self.assertRaises(RangeSyntaxError):
            ResidueRange.parse("120-1")

    def test_garbage_rejected(self):
        for bad in ("1-", "-", "a-b", "1--", "1-2-3", ""):
            with self.assertRaises(RangeSyntaxError, msg=bad):
                ResidueRange.parse(bad)

    def test_str_roundtrip(self):
        self.assertEqual(str(ResidueRange.parse("489-607")), "489-607")
        self.assertEqual(str(ResidueRange.parse("1-113A")), "1-113A")


class ChainSelectionTest(unittest.TestCase):
    def test_empty_selection_keeps_everything(self):
        sel = ChainSelection.parse([])
        self.assertTrue(sel.is_empty)
        self.assertTrue(sel.wants("Z", 999))

    def test_whole_chain(self):
        sel = ChainSelection.parse(["A"])
        self.assertTrue(sel.wants("A", 1))
        self.assertTrue(sel.wants("A", 9999, "B"))
        self.assertFalse(sel.wants("B", 1))

    def test_prototype_syntax(self):
        sel = ChainSelection.parse(["A:1-107", "B:1-120"])
        self.assertTrue(sel.wants("A", 107))
        self.assertFalse(sel.wants("A", 108))
        self.assertTrue(sel.wants("B", 120))
        self.assertFalse(sel.wants("B", 121))
        self.assertFalse(sel.wants("C", 500))

    def test_multiple_ranges_one_chain(self):
        sel = ChainSelection.parse(["C:489-607,610-612"])
        self.assertTrue(sel.wants("C", 607))
        self.assertFalse(sel.wants("C", 608))
        self.assertTrue(sel.wants("C", 611))
        self.assertEqual(sel.describe(), {"C": "489-607,610-612"})

    def test_repeated_chain_accumulates(self):
        sel = ChainSelection.parse(["C:489-607", "C:610-612"])
        self.assertTrue(sel.wants("C", 611))

    def test_describe_whole_chain(self):
        self.assertEqual(ChainSelection.parse(["A"]).describe(), {"A": "*"})

    def test_errors(self):
        for bad in (["AB"], ["AB:1-5"], ["A:"], [""], ["A:1-5-7"]):
            with self.assertRaises(RangeSyntaxError, msg=bad):
                ChainSelection.parse(bad)


class RenameTest(unittest.TestCase):
    def test_parse(self):
        self.assertEqual(parse_rename("A=L,B=H"), {"A": "L", "B": "H"})
        self.assertEqual(parse_rename(""), {})
        self.assertEqual(parse_rename(None), {})
        self.assertEqual(parse_rename(" A = L , "), {"A": "L"})

    def test_swap_is_allowed(self):
        self.assertEqual(parse_rename("A=B,B=A"), {"A": "B", "B": "A"})

    def test_errors(self):
        for bad in ("A", "AB=C", "A=LL", "A=L,A=H"):
            with self.assertRaises(RangeSyntaxError, msg=bad):
                parse_rename(bad)


if __name__ == "__main__":
    unittest.main()
