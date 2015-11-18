#!/usr/bin/env python
import types
import pysvtools.utils


class TestUtils(object):
    def test_vcfheader(self):
        assert isinstance(pysvtools.utils.vcfHeader(), types.StringType)
