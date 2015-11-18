#!/usr/bin/env python
import unittest2

from six import string_types
import pysvtools
import pysvtools.models


class TestModels(unittest2.TestCase):
    def test_event_successinstance(self):
        chrA = "chr1"
        chrB = "chr1"
        chrApos = 1
        chrBpos = 5
        event = pysvtools.models.Event(chrA, chrApos, chrB, chrBpos)

        self.assertEqual(event.chrA, chrA)
        self.assertEqual(event.chrB, chrB)
        self.assertEqual(event.chrApos, chrApos)
        self.assertEqual(event.chrBpos, chrBpos)

    @unittest2.expectedFailure
    def test_event_failinginstance(self):
        chrA = "chr1"
        chrB = "chr1"
        chrApos = "b"
        chrBpos = 5
        event = pysvtools.models.Event(chrA, chrApos, chrB, chrBpos)

    def test_exclusionregion_success(self):
        chrA = "chr1"
        chrApos = 1
        chrBpos = 5
        exregion = pysvtools.models.ExclusionRegion(chrA, chrApos, chrBpos)

        self.assertEqual(exregion.chromosome, chrA)
        self.assertEqual(exregion.start, chrApos)
        self.assertEqual(exregion.end, chrBpos)



    def test_exclusionregion_overlap(self):
        exregionA = pysvtools.models.ExclusionRegion("chrA", 1, 10)

        self.assertFalse(exregionA.overlaps("chrB", 1))
        self.assertFalse(exregionA.overlaps("chrB", 5))
        self.assertFalse(exregionA.overlaps("chrB", 10))

        self.assertFalse(exregionA.overlaps("chrA", 0))
        self.assertFalse(exregionA.overlaps("chrA", 50))

        self.assertTrue(exregionA.overlaps("chrA", 1))
        self.assertTrue(exregionA.overlaps("chrA", 5))
        self.assertTrue(exregionA.overlaps("chrA", 10))

