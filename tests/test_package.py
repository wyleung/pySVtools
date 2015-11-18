#!/usr/bin/env python
import types
import pysvtools


class TestPackage(object):
    def test_versionname_string(self):
        assert isinstance(pysvtools.version("pysvtools"), types.StringType)

    def test_version_tupple(self):
        assert isinstance(pysvtools.__version_info__, types.TupleType)

    def test_version_string(self):
        assert isinstance(pysvtools.__version__, types.StringType)
