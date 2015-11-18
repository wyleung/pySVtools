import sys
from setuptools import setup

if sys.version_info < (2, 6):
    raise Exception('pySVtools requires Python 2.6 or higher.')

# Todo: How does this play with pip freeze requirement files?
requires = ['pyvcf',
            'natsort']

# Python 2.6 does not include the argparse module.
try:
    import argparse
except ImportError:
    requires.append('argparse')

import pysvtools as distmeta

setup(
    name='pysvtools',
    version=distmeta.__version__,
    description='VCF toolkit for SV analysis',
    long_description=distmeta.__doc__,
    author=distmeta.__author__,
    author_email=distmeta.__contact__,
    url=distmeta.__homepage__,
    license='MIT License',
    platforms=['any'],
    packages=['pysvtools', 'pysvtools.models'],
    install_requires=requires,
    test_suite='nose2.collector.collector',
    entry_points = {
        'console_scripts': [
            'mergevcf = pysvtools.merge:main'
        ]
    },
    classifiers = [
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'Intended Audience :: Developers',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Topic :: Scientific/Engineering',
    ],
    keywords='bioinformatics'
)
