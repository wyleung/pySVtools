#!/usr/bin/env python2

from __future__ import print_function

__desc__ = """
    SV events merging
    Takes VCF files as input.
"""
__author__ = "Wai Yi Leung <w.y.leung@lumc.nl>"

import argparse
import collections
import datetime
import itertools
import logging
import os
import pprint
import sys
import re

try:
    import vcf
except:
    print("No PyVCF installation was found, please install with:\n\tpip install pyvcf")
    sys.exit(1)

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

from genomic_regions import ExclusionRegion, buildExclusion
from models import *
# read all samples in memory

def extractTXmate(alt):
    try:
        return re.findall(r"([\d\w\_]+)\:([\d]+)", alt, re.I | re.M)[0]
    except:
        raise IndexError


def extractTXmateINFOFIELD(breakpoints):
    if breakpoints == []:
        return ("0", 1)
    if type(breakpoints) == type(list()):
        try:
            breakpoints = breakpoints[1]
        except:
            print(breakpoints)

    breakpoints = breakpoints.replace('"', '')
    try:
        return re.findall(r"([\d\w\_]+)\:([\d]+)", breakpoints, re.I | re.M)[0]
    except:
        raise IndexError


def getDP(vcf_record):
    if 'DP' in vcf_record.INFO.keys():
        if type(vcf_record.INFO['DP']) == type(list):
            return vcf_record.INFO['DP'][0]
        return vcf_record.INFO['DP']
    elif len(vcf_record.samples):
        return getattr(vcf_record.samples[0].data, 'DP', 0)
    return 0


def firstFromList(arr):
    if type(arr) == type([]):
        return arr[0]
    return arr


def loadEventFromVCF(s, vcf_reader, edb, centerpointFlanking, transonly):
    svDB = collections.OrderedDict()
    skipped_events = 0
    for rec in vcf_reader:
        SVTYPE = rec.INFO['SVTYPE'][0]

        if transonly and SVTYPE != 'CTX':
            continue

        skip = True in list(filter((lambda y: y == True), list(map((lambda x: x.overlaps(rec.CHROM, rec.POS)), edb))))
        if skip:
            skipped_events += 1
            continue

        SVLEN = 0
        if "SVLEN" in rec.INFO.keys():
            if type(rec.INFO['SVLEN']) == type([]):
                SVLEN = rec.INFO['SVLEN'][0]
            elif type(rec.INFO['SVLEN']) == type(0):
                SVLEN = rec.INFO['SVLEN']

        if SVTYPE == 'bITX':
            # intrachromosomal events
            try:
                end = rec.INFO['SVEND'][0]
            except:
                end = rec.INFO['BREAKPOINTS'][0].replace('"', '').split('-')[1]
            finally:
                t = Event(rec.CHROM, rec.POS, rec.CHROM,
                          end, sv_type='ITX',
                          cp_flank=centerpointFlanking,
                          dp=getDP(rec))
                svDB[t.virtualChr] = svDB.get(t.virtualChr, [])
                svDB[t.virtualChr].append(t)
        elif SVTYPE == 'CTX':
            # interchromosomal events
            # check chromosome B:
            try:
                chrB, chrBpos = extractTXmate(str(rec.ALT.pop()))
            except:
                # this is only valid for the old YAMSVC output
                chrB, chrBpos = extractTXmateINFOFIELD(rec.INFO.get('BREAKPOINTS', []))

                t = Event(rec.CHROM, rec.POS, chrB, chrBpos,
                          sv_type='CTX',
                          cp_flank=centerpointFlanking,
                          dp=getDP(rec))
            else:
                t = Event(rec.CHROM, rec.POS, chrB, chrBpos,
                          sv_type='CTX',
                          cp_flank=centerpointFlanking,
                          dp=getDP(rec))
            finally:
                try:
                    svDB[t.virtualChr] = svDB.get(t.virtualChr, [])
                except:
                    pass
                else:
                    svDB[t.virtualChr].append(t)
        elif SVTYPE == 'DEL':
            try:
                if "SVEND" in rec.INFO.keys():
                    end = rec.INFO['SVEND'][0]
                elif "END" in rec.INFO.keys():
                    end = rec.INFO['END'][0]
                elif "SVLEN" in rec.INFO.keys():
                    end = rec.POS + abs(SVLEN)
                t = Event(rec.CHROM, rec.POS, rec.CHROM, end,
                          sv_type=rec.INFO['SVTYPE'][0],
                          cp_flank=centerpointFlanking,
                          dp=getDP(rec))
            except:
                print("Unexpected error:", sys.exc_info()[0])
                raise
            else:
                svDB[t.virtualChr] = svDB.get(t.virtualChr, [])
                svDB[t.virtualChr].append(t)
        else:
            # all other events not covered in this analysis, we only check the overlap
            try:
                if "SVEND" in rec.INFO.keys():
                    if type(rec.INFO['SVEND']) == type([]):
                        end = rec.INFO['SVEND'][0]
                    else:
                        end = rec.INFO['SVEND']
                elif "END" in rec.INFO.keys():
                    end = rec.INFO['END'][0]
                elif "SVLEN" in rec.INFO.keys():
                    end = rec.POS + abs(SVLEN)
                t = Event(rec.CHROM, rec.POS, rec.CHROM, end,
                          sv_type=rec.INFO['SVTYPE'][0],
                          cp_flank=centerpointFlanking,
                          dp=getDP(rec))
            except:
                print("Unexpected error:", sys.exc_info()[0])
                raise
            else:
                svDB[t.virtualChr] = svDB.get(t.virtualChr, [])
                svDB[t.virtualChr].append(t)
    logger.info("Skipped {} events overlapping excluded regions.".format(skipped_events))
    return svDB


def vcfHeader():
    # print the VCF header
    TS_NOW = datetime.datetime.now()
    VCF_HEADER = """##fileformat=VCFv4.1
##fileDate={filedate}
##source=yamsvp-transmerge
##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length of variation">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">""".format(filedate=TS_NOW.strftime("%Y%m%d"))
    return VCF_HEADER + "\n" + "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	default"


# main functions
def startMerge(vcf_files, exclusion_regions, output_file, centerpointFlanking, bedoutput, transonly=False,
               regions_out="regions_out.bed", vcf_output="output.vcf"):
    regions_out_file = open(regions_out, "w")
    vcf_output_file = open(vcf_output, "w")

    samplelist = vcf_files

    sampleDB = collections.OrderedDict()
    svDB = collections.OrderedDict()

    commonhits = collections.OrderedDict()
    edb = []
    if type(exclusion_regions) != type([]):
        exclusion_regions = []
    for exclusion_region in exclusion_regions:
        edb += buildExclusion(exclusion_region)

    for s in samplelist:
        logger.info('Reading SV-events from sample: {} '.format(s))
        sampleDB[s] = vcf.Reader(open(s, 'r'))
        svDB[s] = loadEventFromVCF(s, sampleDB[s], edb, centerpointFlanking, transonly)
        n_events = sum([len(calls) for chromlist, calls in svDB[s].items()])
        logger.info('Loaded SV-events from sample: {} '.format(n_events))

    pairs_to_check = itertools.combinations(svDB.keys(), 2)
    # collect all chromosomes seen:
    chromosomes_to_check = []
    for sample in svDB.keys():
        chromosomes_to_check += svDB[sample].keys()

    chromosomes_to_check = list(set(chromosomes_to_check))
    chromosomes_to_check.sort()

    for (s1, s2) in pairs_to_check:
        _s1 = os.path.basename(s1)
        _s2 = os.path.basename(s2)
        logger.debug('Pairwise compare: {} x {}'.format(_s1, _s2))
        for _chromosome in chromosomes_to_check:
            s1_calls_in_chromosome = svDB[s1].get(_chromosome, [])
            s2_calls_in_chromosome = svDB[s2].get(_chromosome, [])
            s1_n_calls = len(s1_calls_in_chromosome)
            s2_n_calls = len(s2_calls_in_chromosome)

            _match = 0

            for i, t1 in enumerate(s1_calls_in_chromosome):
                for j, t2 in enumerate(s2_calls_in_chromosome):
                    if not (t1 and t2):
                        continue
                    if t1 == t2:
                        # determine the object with the most DP and size
                        if t1.size >= t2.size:
                            _m = t1
                        else:
                            _m = t2

                        # get the hashes from all hits, check weither one of them was already evaluated and thus in the table

                        if t1.matched_in or t2.matched_in:
                            m = t1.matched_in or t2.matched_in
                        else:
                            m = _m.hexdigest

                        t1.matched_in = m
                        t2.matched_in = m

                        # TODO: write getter method for the Event instead now by accessing the internal class variable
                        virtualchrom = _m.virtualChr

                        # split out per chromosome storage
                        commonhits[virtualchrom] = commonhits.get(virtualchrom, collections.OrderedDict())
                        commonhits[virtualchrom][m] = commonhits[virtualchrom].get(m, collections.OrderedDict())
                        commonhits[virtualchrom][m][s1] = t1
                        commonhits[virtualchrom][m][s2] = t2
                        _match += 1

                        # when found, continue? It is the best approach to break this long list comparison?
                        break
            logger.debug("Common hits so far in {}: {} / {} vs {}".format(_chromosome, _match, s1_n_calls, s2_n_calls))


    # vcf header
    print(vcfHeader(), file=vcf_output_file)

    # tsv file
    samplecols = "\t".join(map(lambda x: "{}\tsize".format(os.path.basename(x)), samplelist))
    header_line = "\t".join(['ChrA', 'ChrApos', 'ChrB', 'ChrBpos', 'SVTYPE', 'DP', 'Size', samplecols])

    structural_events = open(output_file, 'w')
    structural_events.write("{}\n".format(header_line))

    bed_structural_events = open(bedoutput, 'w')

    # the keys contain the virtual chromosomes, we need sorted output!
    commonhits_keys = commonhits.keys()
    commonhits_keys.sort()

    all_locations = []

    for virtualChr in commonhits_keys:
        for s, items in commonhits[virtualChr].items():
            if len(items):
                # check which samples has the same
                locations_found = []
                for s in samplelist:
                    if s in items.keys():
                        locations_found.append("{}\t{}".format(items[s], items[s].size))
                        # track all locations found for later intersecting or complementing the set of found/not-found
                        all_locations.append(items[s])
                    else:
                        locations_found.append("\t")
                # get the key with the highest DP
                sorted_by_dp = sorted(items.items(), key=lambda t: t[1].dp, reverse=True)
                fKey = sorted_by_dp[0][0]
                t = items[fKey]

                # write the new vcf record on the commandline
                # TODO: write the DP for each of the callers/sample

                INFOFIELDS = "IMPRECISE;SVTYPE={};END={}".format(
                    t.sv_type,
                    t.chrBpos
                )
                FORMATFIELDS = ":".join(map(str, [
                    '1/.',
                    t.dp])
                )
                print("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\tGT:DP\t{}".format(
                    t.chrA,
                    t.chrApos,
                    '.',
                    '.',
                    t.vcf_alt,
                    '100',
                    'PASS',
                    INFOFIELDS,
                    FORMATFIELDS
                ), file=vcf_output_file
                )

                structural_events.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                    t.chrA,
                    t.chrApos,
                    t.chrB,
                    t.chrBpos,
                    t.sv_type,
                    t.dp,
                    t.size,
                    "\t".join(locations_found)))
                if t.sv_type in ["INS", 'DEL', 'ITX']:
                    bed_structural_events.write("{chrom}\t{start}\t{end}\t{annot}\n".format(
                        chrom=t.chrA,
                        start=t.chrApos,
                        end=t.chrBpos,
                        annot="SVTYPE={svtype};DP={dp};SIZE={size}".format(
                            svtype=t.sv_type,
                            dp=t.dp,
                            size=t.size
                        )
                    ))
                else:
                    # for CTX events, we write 2 tracks, one for the left breakpoint and one for th right breakpoint
                    # both sites are flanked with 50 bases to generate a valid bed file and used in IGV to mark the region
                    bed_structural_events.write("{chrom}\t{start}\t{end}\t{annot}\n".format(
                        chrom=t.chrA,
                        start=t.chrApos - 50,
                        end=t.chrApos + 50,
                        annot="SVTYPE={svtype};DP={dp};SIZE={size};MATE={chrb}:{chrbpos}".format(
                            svtype=t.sv_type,
                            dp=t.dp,
                            size=t.size,
                            chrb=t.chrB,
                            chrbpos=t.chrBpos
                        )
                    ))
                    bed_structural_events.write("{chrom}\t{start}\t{end}\t{annot}\n".format(
                        chrom=t.chrB,
                        start=t.chrBpos - 50,
                        end=t.chrBpos + 50,
                        annot="SVTYPE={svtype};DP={dp};SIZE={size};MATE={chrb}:{chrbpos}".format(
                            svtype=t.sv_type,
                            dp=t.dp,
                            size=t.size,
                            chrb=t.chrA,
                            chrbpos=t.chrApos
                        )
                    ))

    for loc in all_locations:
        print(loc.bedRow, file=regions_out_file)

    structural_events.close()
    regions_out_file.close()

class Intersection(object):
    def __init__(self):
        pass



def main():
    parser = argparse.ArgumentParser()

    parser.add_argument('-c', '--exclusion_regions', action='append',
                        help='Exclusion regions file in BED format')

    parser.add_argument('-f', '--flanking', type=int,
                        help='Centerpoint flanking [100]', default=100)

    parser.add_argument('-t', '--translocation_only', action='store_true',
                        help='Do translocations only', required=False, default=False)

    parser.add_argument('-i', '--vcf', nargs='+',
                        help='The VCF(s) to compare, can be supplied multiple times')
    parser.add_argument('-o', '--output',
                        help='Output summary to [sample.tsv]', default='sample.tsv')
    parser.add_argument('-b', '--bedoutput',
                        help='Output bed file to [sample.bed]', default='sample.bed')
    parser.add_argument('-v', '--vcfoutput',
                        help='Output summary to [sample.vcf]', default='sample.vcf')
    parser.add_argument('-r', '--regions_out',
                        help='Output all regions to [regions_out.bed]', default='regions_out.bed')
    args = parser.parse_args()

    if args.vcf == None or len(args.vcf) < 2:
        logger.error("Please supply at least 2 VCF files to merge")
        sys.exit(1)

    startMerge(args.vcf, args.exclusion_regions,
               args.output, args.flanking, args.bedoutput, args.translocation_only, args.regions_out, args.vcfoutput)


if __name__ == "__main__":
    main()
