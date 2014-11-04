# pySVTools

A collection of usefull operations on VCF files containing structural variants calls.

# Installation

    # Stable version
    pip install pySVtools

    # Bleeding edge version:
    pip install git+https://github.com/wyleung/pySVtools.git#egg=pysvtools
    

# Merging / intersecting

For `DEL` and `INS` events, you can intersect 2 or more `VCF` -files using the following command:

    mergeVCF -f 100 -i sample1.vcf sample2.vcf -o intersected.tsv -b intersected.bed > intersected.vcf


# Help

```bash
usage: merge.py [-h] [-c CENTROMERS] [-f FLANKING] [-t] [-i VCF [VCF ...]]
                [-o OUTPUT] [-b BEDOUTPUT]

optional arguments:
  -h, --help            show this help message and exit
  -c CENTROMERS, --centromers CENTROMERS
                        Centromers definitions file in BED format
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
```


# Features

 1. Summarizing SV events
 1. Intersecting SV events, using multiple VCF. Usefull for finding recuring event accros multiple 'samples'
 1. Filtering events using known regions using ``bed`` files describing f.e.g. GC-rich regions and/or known CNV regions

# Contact

Please file an issue report on Github if there are any questions using this library/tool

