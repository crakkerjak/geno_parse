#!/usr/bin/python3
import unittest
import geno_parse
from collections import defaultdict


class TestRanges(unittest.TestCase):
    def test_read_gaps(self):
        actual = defaultdict(list)
        geno_parse.read_gaps('test_gap.bed', actual)
        expected = defaultdict(list)
        expected = {
            'chr1': [(25000, 26000), (26000, 27000), (29900, 29905),
                    (30500, 30505), (32000, 33000), (33000, 34000)],
            'chr2': [(27000, 28000), (28000, 29000), (29950, 29960)]
        }
        self.assertEqual(actual, expected)

    # merge neighboring ranges
    def test_merge_ranges_1(self):
        actual = defaultdict(list)
        actual['chr1'] = [(0,99),(100,200)]
        expected = defaultdict(list)
        expected['chr1'] = [(0,200)]
        geno_parse.merge_ranges(actual)
        self.assertEqual(actual, expected)

    # do not merge disjoint
    def test_merge_ranges_2(self):
        actual = defaultdict(list)
        actual['chr1'] = [(0,98),(100,200)]
        expected = defaultdict(list)
        expected['chr1'] = [(0,98),(100,200)]
        geno_parse.merge_ranges(actual)
        self.assertEqual(actual, expected)

    # merge 3
    def test_merge_ranges_3(self):
        actual = defaultdict(list)
        actual['chr1'] = [(0,99),(100,200),(200,300)]
        expected = defaultdict(list)
        expected['chr1'] = [(0,300)]
        geno_parse.merge_ranges(actual)
        self.assertEqual(actual, expected)

    # merge 5
    def test_merge_ranges_3(self):
        actual = defaultdict(list)
        actual['chr1'] = [(0,99),(100,200),(200,300),(300,400),(401,500)]
        expected = defaultdict(list)
        expected['chr1'] = [(0,500)]
        geno_parse.merge_ranges(actual)
        self.assertEqual(actual, expected)

    # do not merge from different chromosome names
    def test_merge_ranges_4(self):
        actual = defaultdict(list)
        actual['chr1'] = [(0,99),(200,300)]
        actual['chr2'] = [(100,200)]
        expected = defaultdict(list)
        expected['chr1'] = [(0,99),(200,300)]
        expected['chr2'] = [(100,200)]
        geno_parse.merge_ranges(actual)
        self.assertEqual(actual, expected)

    # merge madness!!!
    def test_merge_ranges_5(self):
        actual = defaultdict(list)
        actual['chr1'] = [
                (0, 100), (101, 200), (23000, 24125), (24125, 27778),
                (45597, 45882), (91794, 162650), (162650,165776),
                (165777, 487528), (579465, 579475), (635348, 661995)]
        expected = defaultdict(list)
        expected['chr1'] = [
                (0, 200), (23000, 27778), (45597, 45882), (91794, 487528),
                (579465, 579475), (635348, 661995)]
        geno_parse.merge_ranges(actual)
        self.assertEqual(actual, expected)

    # line below range
    def test_filter_ranges_1(self):
        remaining = ['chr1\t0\t99']
        ranges_to_remove = defaultdict(list)
        ranges_to_remove['chr1'] = [(100,200)]

        removed = geno_parse.filter_ranges(remaining, ranges_to_remove, 0)

        expected_remaining = ['chr1\t0\t99']
        expected_removed = []

        self.assertEqual(remaining, expected_remaining)
        self.assertEqual(removed, expected_removed)

    # intersection length of only 1, low
    def test_filter_ranges_2(self):
        remaining = ['chr1\t0\t100']
        ranges_to_remove = defaultdict(list)
        ranges_to_remove['chr1'] = [(100,200)]

        removed = geno_parse.filter_ranges(remaining, ranges_to_remove, 0)

        expected_remaining = []
        expected_removed = ['chr1\t0\t100']

        self.assertEqual(remaining, expected_remaining)
        self.assertEqual(removed, expected_removed)

    # line above range
    def test_filter_ranges_3(self):
        remaining = ['chr1\t201\t300']
        ranges_to_remove = defaultdict(list)
        ranges_to_remove['chr1'] = [(100,200)]

        removed = geno_parse.filter_ranges(remaining, ranges_to_remove, 0)

        expected_remaining = ['chr1\t201\t300']
        expected_removed = []

        self.assertEqual(remaining, expected_remaining)
        self.assertEqual(removed, expected_removed)

    # intersection length of only 1, high
    def test_filter_ranges_4(self):
        remaining = ['chr1\t200\t300']
        ranges_to_remove = defaultdict(list)
        ranges_to_remove['chr1'] = [(100,200)]

        removed = geno_parse.filter_ranges(remaining, ranges_to_remove, 0)

        expected_remaining = []
        expected_removed = ['chr1\t200\t300']

        self.assertEqual(remaining, expected_remaining)
        self.assertEqual(removed, expected_removed)

    # line fully inside range
    def test_filter_ranges_5(self):
        remaining = ['chr1\t150\t150']
        ranges_to_remove = defaultdict(list)
        ranges_to_remove['chr1'] = [(100,200)]

        removed = geno_parse.filter_ranges(remaining, ranges_to_remove, 0)

        expected_remaining = []
        expected_removed = ['chr1\t150\t150']

        self.assertEqual(remaining, expected_remaining)
        self.assertEqual(removed, expected_removed)

    # max_overlap 1 with overlap of 1
    def test_filter_ranges_6(self):
        remaining = ['chr1\t200\t300']
        ranges_to_remove = defaultdict(list)
        ranges_to_remove['chr1'] = [(100,200)]

        removed = geno_parse.filter_ranges(remaining, ranges_to_remove, 1)

        expected_remaining = ['chr1\t200\t300']
        expected_removed = []

        self.assertEqual(remaining, expected_remaining)
        self.assertEqual(removed, expected_removed)

    # max_overlap 1 with overlap of 2
    def test_filter_ranges_7(self):
        remaining = ['chr1\t199\t300']
        ranges_to_remove = defaultdict(list)
        ranges_to_remove['chr1'] = [(100,200)]

        removed = geno_parse.filter_ranges(remaining, ranges_to_remove, 1)

        expected_remaining = []
        expected_removed = ['chr1\t199\t300']

        self.assertEqual(remaining, expected_remaining)
        self.assertEqual(removed, expected_removed)

    # multiline all remain
    def test_filter_ranges_8(self):
        remaining = ['chr1\t0\t99', 'chr1\t201\t300', 'chr1\t350\t400']
        ranges_to_remove = defaultdict(list)
        ranges_to_remove['chr1'] = [(100,200)]

        removed = geno_parse.filter_ranges(remaining, ranges_to_remove, 0)

        expected_remaining = ['chr1\t0\t99', 'chr1\t201\t300', 'chr1\t350\t400']
        expected_removed = []

        self.assertEqual(remaining, expected_remaining)
        self.assertEqual(removed, expected_removed)

    # multiline all remove
    def test_filter_ranges_9(self):
        remaining = ['chr1\t0\t100', 'chr1\t200\t300', 'chr1\t150\t400']
        ranges_to_remove = defaultdict(list)
        ranges_to_remove['chr1'] = [(100,200)]

        removed = geno_parse.filter_ranges(remaining, ranges_to_remove, 0)

        expected_remaining = []
        expected_removed = ['chr1\t0\t100', 'chr1\t200\t300', 'chr1\t150\t400']

        self.assertEqual(remaining, expected_remaining)
        self.assertEqual(removed, expected_removed)


    # multiline remove middle
    def test_filter_ranges_10(self):
        remaining = ['chr1\t0\t99', 'chr1\t200\t300', 'chr1\t201\t400']
        ranges_to_remove = defaultdict(list)
        ranges_to_remove['chr1'] = [(100,200)]

        removed = geno_parse.filter_ranges(remaining, ranges_to_remove, 0)

        expected_remaining = ['chr1\t0\t99', 'chr1\t201\t400']
        expected_removed = ['chr1\t200\t300']

        self.assertEqual(remaining, expected_remaining)
        self.assertEqual(removed, expected_removed)

    # multiline remove first
    def test_filter_ranges_11(self):
        remaining = ['chr1\t200\t300', 'chr1\t0\t99', 'chr1\t201\t400']
        ranges_to_remove = defaultdict(list)
        ranges_to_remove['chr1'] = [(100,200)]

        removed = geno_parse.filter_ranges(remaining, ranges_to_remove, 0)

        expected_remaining = ['chr1\t0\t99', 'chr1\t201\t400']
        expected_removed = ['chr1\t200\t300']

        self.assertEqual(remaining, expected_remaining)
        self.assertEqual(removed, expected_removed)

    # multiline remove last
    def test_filter_ranges_10(self):
        remaining = ['chr1\t201\t400', 'chr1\t0\t99', 'chr1\t200\t300']
        ranges_to_remove = defaultdict(list)
        ranges_to_remove['chr1'] = [(100,200)]

        removed = geno_parse.filter_ranges(remaining, ranges_to_remove, 0)

        expected_remaining = ['chr1\t201\t400', 'chr1\t0\t99']
        expected_removed = ['chr1\t200\t300']

        self.assertEqual(remaining, expected_remaining)
        self.assertEqual(removed, expected_removed)

    # different chromosome test
    def test_filter_ranges_11(self):
        remaining = ['chr2\t0\t150']
        ranges_to_remove = defaultdict(list)
        ranges_to_remove['chr1'] = [(100,200)]

        removed = geno_parse.filter_ranges(remaining, ranges_to_remove, 0)

        expected_remaining = ['chr2\t0\t150']
        expected_removed = []

        self.assertEqual(remaining, expected_remaining)
        self.assertEqual(removed, expected_removed)

    # different chromosome test 2
    def test_filter_ranges_12(self):
        remaining = ['chr2\t0\t150']
        ranges_to_remove = defaultdict(list)
        ranges_to_remove['chr1'] = [(100,200)]
        ranges_to_remove['chr2'] = [(200,300)]

        removed = geno_parse.filter_ranges(remaining, ranges_to_remove, 0)

        expected_remaining = ['chr2\t0\t150']
        expected_removed = []

        self.assertEqual(remaining, expected_remaining)
        self.assertEqual(removed, expected_removed)

    # different chromosome test 3
    def test_filter_ranges_13(self):
        remaining = ['chr2\t0\t150', 'chr1\t0\t150']
        ranges_to_remove = defaultdict(list)
        ranges_to_remove['chr1'] = [(100,200)]
        ranges_to_remove['chr2'] = [(200,300)]

        removed = geno_parse.filter_ranges(remaining, ranges_to_remove, 0)

        expected_remaining = ['chr2\t0\t150']
        expected_removed = ['chr1\t0\t150']

        self.assertEqual(remaining, expected_remaining)
        self.assertEqual(removed, expected_removed)

    # intersection in different index in ranges_to_remove
    def test_filter_ranges_14(self):
        remaining = ['chr1\t50\t150']
        ranges_to_remove = defaultdict(list)
        ranges_to_remove['chr1'] = [(0,20),(100,200),(300,250)]

        removed = geno_parse.filter_ranges(remaining, ranges_to_remove, 0)

        expected_remaining = []
        expected_removed = ['chr1\t50\t150']

        self.assertEqual(remaining, expected_remaining)
        self.assertEqual(removed, expected_removed)


if __name__ == '__main__':
    unittest.main()
