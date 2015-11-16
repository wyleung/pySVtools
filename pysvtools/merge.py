#!/usr/bin/env python2

from __future__ import print_function

import pprint

__desc__ = """
    Merging procedure for Structural Variation events.
    Follows the idea of centerpoint matching to allow flexible match vs. reciprocal overlap.
"""
__author__ = "Wai Yi Leung <w.y.leung@lumc.nl>"

import argparse
import collections
import itertools
from natsort import natsorted
import logging
import os
import sys

try:
    import vcf
except:
    print("No PyVCF installation was found, please install with:\n\tpip install pyvcf")
    sys.exit(1)

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

from pysvtools.models import Event, ExclusionRegion
from pysvtools.utils import extractTXmate, extractDPFromRecord, getSVType, getSVLEN, formatBedTrack, \
    formatVCFRecord, vcfHeader, build_exclusion


class VCFEventLoader(object):
    """
        Load VCF File and transform VCF record into an `Event`
    """

    def __init__(self):
        pass


class SVMerger(object):
    def __init__(self):
        pass


class ReportExport(object):
    def __init__(self):
        pass


# read all samples in memory
def loadEventFromVCF(s, vcf_reader, edb, centerpointFlanking, transonly, svmethod=""):
    """
        Loading VCF records and transform them to `Event`
    """
    svDB = collections.OrderedDict()
    skipped_events = 0

    for record in vcf_reader:
        SVTYPE = getSVType(record)

        if transonly and SVTYPE not in ['CTX', 'TRA']:
            continue

        skip = True in list(
            filter((lambda y: y == True), list(map((lambda x: x.overlaps(record.CHROM, record.POS)), edb))))
        if skip:
            skipped_events += 1
            continue

        SVLEN = getSVLEN(record)

        if SVTYPE == 'bITX':
            # intrachromosomal events
            try:
                end = record.INFO['SVEND'][0]
            except:
                end = record.INFO['BREAKPOINTS'][0].replace('"', '').split('-')[1]
            finally:
                t = Event(record.CHROM, record.POS, record.CHROM,
                          end, sv_type="TRA",
                          cp_flank=centerpointFlanking,
                          dp=extractDPFromRecord(record),
                          svmethod=svmethod)
                svDB[t.virtualChr] = svDB.get(t.virtualChr, [])
                svDB[t.virtualChr].append(t)
        elif SVTYPE in ['CTX', 'TRA']:
            # interchromosomal events
            # check chromosome B:
            chrB, chrBpos = extractTXmate(record)

            t = Event(record.CHROM, record.POS, chrB, chrBpos,
                      sv_type='TRA',
                      cp_flank=centerpointFlanking,
                      dp=extractDPFromRecord(record),
                      svmethod=svmethod)
            try:
                svDB[t.virtualChr] = svDB.get(t.virtualChr, [])
            except:
                pass
            else:
                svDB[t.virtualChr].append(t)

        elif SVTYPE == 'DEL':
            try:
                if "SVEND" in record.INFO.keys():
                    if type(record.INFO['SVEND']) == type([]):
                        end = record.INFO['SVEND'][0]
                    else:
                        end = record.INFO['SVEND']
                elif "END" in record.INFO.keys():
                    if type(record.INFO['END']) == type([]):
                        end = record.INFO['END'][0]
                    else:
                        end = record.INFO['END']
                elif "SVLEN" in record.INFO.keys():
                    end = record.POS + abs(SVLEN)
                t = Event(record.CHROM, record.POS, record.CHROM, end,
                          sv_type=SVTYPE,
                          cp_flank=centerpointFlanking,
                          dp=extractDPFromRecord(record),
                          svmethod=svmethod)
            except:
                print("Unexpected error:", sys.exc_info()[0])
                raise
            else:
                svDB[t.virtualChr] = svDB.get(t.virtualChr, [])
                svDB[t.virtualChr].append(t)
        else:
            # all other events not covered in this analysis, we only check the overlap
            try:
                if "SVEND" in record.INFO.keys():
                    if type(record.INFO['SVEND']) == type([]):
                        end = record.INFO['SVEND'][0]
                    else:
                        end = record.INFO['SVEND']
                elif "END" in record.INFO.keys():
                    if type(record.INFO['END']) == type([]):
                        end = record.INFO['END'][0]
                    else:
                        end = record.INFO['END']
                elif "SVLEN" in record.INFO.keys():
                    end = record.POS + abs(SVLEN)
                t = Event(record.CHROM, record.POS, record.CHROM, end,
                          sv_type=SVTYPE,
                          cp_flank=centerpointFlanking,
                          dp=extractDPFromRecord(record),
                          svmethod=svmethod)
            except:
                print("Unexpected error:", sys.exc_info()[0])
                raise
            else:
                svDB[t.virtualChr] = svDB.get(t.virtualChr, [])
                svDB[t.virtualChr].append(t)
    logger.info("Skipped {} events overlapping excluded regions.".format(skipped_events))
    return svDB


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
        edb += build_exclusion(exclusion_region)

    for s in samplelist:
        logger.info('Reading SV-events from sample: {} '.format(s))
        sampleDB[s] = vcf.Reader(open(s, 'r'))

        # extract SV caller from header
        sv_caller = sampleDB[s].metadata.get('source', [os.path.basename(s).strip('.vcf')]).pop(0).split(' ').pop(0)

        svDB[s] = loadEventFromVCF(s, sampleDB[s], edb, centerpointFlanking, transonly, sv_caller)
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
    samplecols = "\t".join(map(lambda x: "{}\tsize".format(os.path.basename(x).strip(".vcf")), samplelist))
    header_line = "\t".join(['ChrA', 'ChrApos', 'ChrB', 'ChrBpos', 'SVTYPE', 'DP', 'Size', samplecols])

    tsv_report_output = open(output_file, 'w')
    tsv_report_output.write("{}\n".format(header_line))

    bed_structural_events = open(bedoutput, 'w')
    all_locations = []

    for virtualChr in natsorted(commonhits.keys()):
        for s, items in sorted(commonhits[virtualChr].items(), key=lambda hit: hit[1].items()[0][1].chrApos):
            if len(items):
                # check which samples has the same
                locations_found = []
                for sample in samplelist:
                    if sample in items.keys():
                        locations_found.append("{}\t{}".format(items[sample], items[sample].size))
                        # track all locations found for later intersecting or complementing the set of found/not-found
                        all_locations.append(items[sample])
                    else:
                        locations_found.append("\t")
                # get the key with the highest DP
                sorted_by_dp = sorted(items.items(), key=lambda hit: hit[1].dp, reverse=True)
                fKey = sorted_by_dp[0][0]
                t = items[fKey]

                print(formatVCFRecord(t), file=vcf_output_file)

                tsv_report_output.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                    t.chrA,
                    t.chrApos,
                    t.chrB,
                    t.chrBpos,
                    t.sv_type,
                    t.dp,
                    t.size,
                    "\t".join(locations_found)))
                bed_structural_events.write(formatBedTrack(t))

    for loc in all_locations:
        print(loc.bedRow, file=regions_out_file)

    tsv_report_output.close()
    regions_out_file.close()


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
