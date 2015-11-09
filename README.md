# pySVTools

A collection of usefull operations on VCF files containing structural variants calls.

# Installation

    # Stable version
    pip install pySVtools

    # Bleeding edge version:
    pip install git+https://github.com/wyleung/pySVtools.git#egg=pysvtools

# Dependencies

Installation of the dependencies is done if the installation is done using `easy_install` or `pip`. However, when used directly from the source, you should install the following external libraries:

( or use `pip install -r requirements.txt`)


 - [pyvcf](https://github.com/jamescasbon/PyVCF)
 


# Intersecting IN/DEL events

For `DEL` and `INS` events, you can intersect 2 or more `VCF`-files using the following command:

    mergevcf -f 100 -i sample1.vcf sample2.vcf \
                    -o intersected.tsv -b intersected.bed > intersected.vcf

The resulting `tsv` file is a matrix listing the:

 - Intersected hits, with both breakpoints (5' and 3'), coverage (DP) and size.
 - Hit `location` in each sample, and size of the event



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

