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

    def test_event_equal(self):
        chrA = "chr1"
        chrB = "chr1"
        chrApos = 1
        chrBpos = 5
        event = pysvtools.models.Event(chrA, chrApos, chrB, chrBpos, sv_type="DEL")
        self.assertEqual(event, event)

    def test_event_equal_a_in_b(self):
        chrA = "chr1"
        chrB = "chr1"
        chrApos = 1
        chrBpos = 5
        eventA = pysvtools.models.Event(chrA, chrApos + 3, chrB, chrBpos + 3, sv_type="DEL")
        eventB = pysvtools.models.Event(chrA, chrApos, chrB, chrBpos, sv_type="DEL")
        self.assertEqual(eventA, eventB)

    def test_event_equal_a_surround_b(self):
        chrA = "chr1"
        chrB = "chr1"
        chrApos = 2
        chrBpos = 5
        eventA = pysvtools.models.Event(chrA, chrApos - 1, chrB, chrBpos + 1, sv_type="DEL")
        eventB = pysvtools.models.Event(chrA, chrApos, chrB, chrBpos, sv_type="DEL")
        self.assertEqual(eventA, eventB)

    def test_event_equal_b_surround_a(self):
        chrA = "chr1"
        chrB = "chr1"
        chrApos = 2
        chrBpos = 5
        eventA = pysvtools.models.Event(chrA, chrApos, chrB, chrBpos, sv_type="DEL")
        eventB = pysvtools.models.Event(chrA, chrApos - 1, chrB, chrBpos + 1, sv_type="DEL")
        self.assertEqual(eventA, eventB)

    def test_event_equal_no_match_by_flanking(self):
        chrA = "chr1"
        chrB = "chr1"
        chrApos = 1
        chrBpos = 5
        eventA = pysvtools.models.Event(chrA, chrApos, chrB, chrBpos, sv_type="DEL")
        eventB = pysvtools.models.Event(chrA, chrApos + 5, chrB, chrBpos + 5, sv_type="DEL")
        self.assertEqual(eventA, eventB)

    def test_event_size(self):
        chrA = "chr1"
        chrB = "chr1"
        chrApos = 1
        chrBpos = 5
        eventA = pysvtools.models.Event(chrA, chrApos, chrB, chrBpos, sv_type="DEL")
        self.assertTrue(eventA.size >= 0)
        self.assertEqual(eventA.size, 4)

    def test_event_notequal(self):
        chrA = "chr1"
        chrB = "chr5"
        chrApos = 1
        chrBpos = 5
        eventA = pysvtools.models.Event(chrA, chrApos, chrA, chrBpos, sv_type="DEL")
        eventB = pysvtools.models.Event(chrA, chrApos, chrB, chrBpos, sv_type="DEL")
        self.assertNotEqual(eventA, eventB)



    def test_event_vcfalt(self):
        chrA = "chr1"
        chrB = "chr1"
        chrApos = 1
        chrBpos = 5
        event = pysvtools.models.Event(chrA, chrApos, chrB, chrBpos, sv_type="DEL")
        self.assertIsInstance(event.vcf_alt, string_types)
        self.assertEqual(event.vcf_alt, ".")

    def test_event_vcfalt_itx(self):
        chrA = "chr1"
        chrB = "chr1"
        chrApos = 1
        chrBpos = 5
        event = pysvtools.models.Event(chrA, chrApos, chrB, chrBpos, sv_type="ITX")
        self.assertIsInstance(event.vcf_alt, string_types)
        self.assertEqual(event.vcf_alt, "N[chr1:5[")

    def test_event_hexdigest(self):
        chrA = "chr1"
        chrB = "chr1"
        chrApos = 1
        chrBpos = 5
        event = pysvtools.models.Event(chrA, chrApos, chrB, chrBpos)
        self.assertIsInstance(event.hexdigest, string_types)

    def test_event_reprstring(self):
        chrA = "chr1"
        chrB = "chr1"
        chrApos = 1
        chrBpos = 5
        event = pysvtools.models.Event(chrA, chrApos, chrB, chrBpos)
        self.assertIsInstance(event.__repr__(), string_types)
        self.assertIn(chrA, event.__repr__())
        self.assertIn(chrApos.__str__(), event.__repr__())
        self.assertIn(chrB, event.__repr__())
        self.assertIn(chrBpos.__str__(), event.__repr__())

    def test_event_reprstring_chrB(self):
        chrA = "chr1"
        chrB = "chr2"
        chrApos = 1
        chrBpos = 5
        event = pysvtools.models.Event(chrA, chrApos, chrB, chrBpos)
        self.assertIsInstance(event.__repr__(), string_types)
        self.assertIn(chrA, event.__repr__())
        self.assertIn(chrApos.__str__(), event.__repr__())
        self.assertIn(chrB, event.__repr__())
        self.assertIn(chrBpos.__str__(), event.__repr__())

    def test_event_string(self):
        chrA = "chr1"
        chrB = "chr1"
        chrApos = 1
        chrBpos = 5
        event = pysvtools.models.Event(chrA, chrApos, chrB, chrBpos)
        self.assertIsInstance(event.__str__(), string_types)
        self.assertIn(chrA, event.bedRow)
        self.assertIn(chrApos.__str__(), event.__repr__())
        self.assertIn(chrB, event.__repr__())
        self.assertIn(chrBpos.__str__(), event.__repr__())

    def test_event_bedrow(self):
        chrA = "chr1"
        chrB = "chr1"
        chrApos = 1
        chrBpos = 5
        event = pysvtools.models.Event(chrA, chrApos, chrB, chrBpos)
        self.assertIsInstance(event.bedRow, string_types)
        self.assertIn(chrA, event.bedRow)
        self.assertIn(chrApos.__str__(), event.bedRow)
        self.assertIn(chrB, event.bedRow)
        self.assertIn(chrBpos.__str__(), event.bedRow)

    def test_event_bedrow_chrB(self):
        chrA = "chr1"
        chrB = "chr2"
        chrApos = 1
        chrBpos = 5
        event = pysvtools.models.Event(chrA, chrApos, chrB, chrBpos)
        self.assertIsInstance(event.bedRow, string_types)
        self.assertIn(chrA, event.bedRow)
        self.assertIn(chrApos.__str__(), event.bedRow)
        self.assertIn(chrB, event.bedRow)
        self.assertIn(chrBpos.__str__(), event.bedRow)

    def test_exclusionregion_success(self):
        chrA = "chr1"
        chrApos = 1
        chrBpos = 5
        exregion = pysvtools.models.ExclusionRegion(chrA, chrApos, chrBpos)

        self.assertEqual(exregion.chromosome, chrA)
        self.assertEqual(exregion.start, chrApos)
        self.assertEqual(exregion.end, chrBpos)

    def test_exclusionregion_reprstring(self):
        exregionA = pysvtools.models.ExclusionRegion("chrA", 1, 10)
        self.assertIsInstance(exregionA.__repr__(), string_types)

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
