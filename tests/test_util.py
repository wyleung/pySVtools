#!/usr/bin/env python

from six import string_types
import pysvtools.utils


class TestUtils(object):
    def test_vcfheader(self):
        assert isinstance(pysvtools.utils.vcfHeader(), string_types)
