#!/usr/bin/env python

from six import string_types
import pysvtools


class TestPackage(object):

    def test_version_tupple(self):
        assert isinstance(pysvtools.__version_info__, tuple)

    def test_version_string(self):
        assert isinstance(pysvtools.__version__, string_types)
