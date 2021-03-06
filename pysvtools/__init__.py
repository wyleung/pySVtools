"""
pySVtools: VCF toolkit for SV analysis


Copyright (c) 2013-2014 Leiden University Medical Center <sasc@lumc.nl>
Copyright (c) 2013-2015 Wai Yi Leung <w.y.leung@lumc.nl>

Licensed under the MIT license, see the LICENSE file.
"""

# On the event of a new release, we update the __version_info__ package
# global and set RELEASE to True.
# Before a release, a development version is denoted by a __version_info__
# ending with a 'dev' item and RELEASE is set to False.
#
# We follow a versioning scheme compatible with setuptools [1] where the
# __version_info__ variable always contains the version of the upcomming
# release (and not that of the previous release), post-fixed with a 'dev'
# item. Only in a release commit, this 'dev' item is removed (and added
# again in the next commit).
#
# [1] http://peak.telecommunity.com/DevCenter/setuptools#specifying-your-project-s-version

RELEASE = False

__version_info__ = ('0', '1', '3')

__version__ = '.'.join(__version_info__)
__author__ = 'LUMC, Wai Yi Leung'
__contact__ = 'w.y.leung@lumc.nl'
__homepage__ = 'https://github.com/wyleung/pySVtools'

usage = __doc__.split("\n\n\n")