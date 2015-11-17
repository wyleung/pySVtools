# pySVTools

A collection of usefull operations on VCF files containing structural variants calls.

[![PyPI](https://img.shields.io/pypi/v/pySVtools.svg)](https://pypi.python.org/pypi/pySVtools)  [![PyPI](https://img.shields.io/pypi/wheel/pySVtools.svg)](https://pypi.python.org/pypi/pySVtools)  [![GitHub issues](https://img.shields.io/github/issues/wyleung/pySVtools.svg)](https://github.com/wyleung/pySVtools/issues)  [![GitHub stars](https://img.shields.io/github/stars/wyleung/pySVtools.svg)](https://github.com/wyleung/pySVtools/stargazers)  [![Coverage Status](https://coveralls.io/repos/wyleung/pySVtools/badge.svg?branch=master&service=github)](https://coveralls.io/github/wyleung/pySVtools?branch=master)  [![Build Status](https://travis-ci.org/wyleung/pySVtools.svg?branch=master)](https://travis-ci.org/wyleung/pySVtools)

# Installation

    # Stable version
    pip install -U pySVtools

    # Bleeding edge version:
    pip install -U git+https://github.com/wyleung/pySVtools.git#egg=pysvtools

# Dependencies

Installation of the dependencies is done if the installation is done using `easy_install` or `pip`. However, when used directly from the source, you should install the following external libraries:

( or use `pip install -r requirements.txt`)


 - [pyvcf](https://github.com/jamescasbon/PyVCF)
 


# Intersecting IN/DEL events

For `DEL` and `INS` events, you can intersect 2 or more `VCF`-files using the following command:

    mergevcf -f 100 -i sample1.vcf sample2.vcf \
                    -o intersected.tsv -b intersected.bed -v intersected.vcf

The resulting `tsv` file is a matrix listing the:

 - Intersected hits, with both breakpoints (5' and 3'), coverage (DP) and size.
 - Hit `location` in each sample, and size of the event

# Intersecting TRA/ CTX events

In order to merge translocations, one should set the flanking margin to a higher number.
Recommended setting is to try out with `-t -f 2000` first, this will give some confident calls.
One can allow more flanking by increasing the `-f` value. F.e.g.: `-t -f 5000` to allow 5kb difference in the centerpoint.




# Help

```bash
usage: merge.py [-h] [-c EXCLUSION_REGIONS] [-f FLANKING] [-t]
                [-i VCF [VCF ...]] [-o OUTPUT] [-b BEDOUTPUT] [-v VCFOUTPUT]
                [-r REGIONS_OUT]

optional arguments:
  -h, --help            show this help message and exit
  -c EXCLUSION_REGIONS, --exclusion_regions EXCLUSION_REGIONS
                        Exclusion regions file in BED format
  -f FLANKING, --flanking FLANKING
                        Centerpoint flanking [100]
  -t, --translocation_only
                        Do translocations only
  -i VCF [VCF ...], --vcf VCF [VCF ...]
                        The VCF(s) to compare, can be supplied multiple times
  -o OUTPUT, --output OUTPUT
                        Output summary to [sample.tsv]
  -b BEDOUTPUT, --bedoutput BEDOUTPUT
                        Output bed file to [sample.bed]
  -v VCFOUTPUT, --vcfoutput VCFOUTPUT
                        Output summary to [sample.vcf]
  -r REGIONS_OUT, --regions_out REGIONS_OUT
                        Output all regions to [regions_out.bed]
```


# Features

 1. Intersecting SV events, using multiple VCF. Usefull for finding recuring event accros multiple 'samples'
 1. Filtering events using known regions using ``bed`` files describing f.e.g. GC-rich regions and/or known CNV regions
 1. Summarizing SV events in a nice LaTeX table

# Contact

Please file an issue report on Github if there are any questions using this library/tool

