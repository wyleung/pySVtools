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
        assert event.chrA == chrA
        assert event.chrB == chrB
        assert event.chrApos == chrApos
        assert event.chrBpos == chrBpos

    @unittest2.expectedFailure
    def test_event_failinginstance(self):
        chrA = "chr1"
        chrB = "chr1"
        chrApos = "b"
        chrBpos = 5
        event = pysvtools.models.Event(chrA, chrApos, chrB, chrBpos)
