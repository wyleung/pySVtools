#!/usr/bin/env python2

from __future__ import print_function

from merge import loadEventFromVCF

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

def matchEvents(sampleName, s2_calls_in_chromosome, allhits, virtualChromosome):
    s1_calls_in_chromosome = allhits.get(virtualChromosome, [])
    s1_n_calls = len(s1_calls_in_chromosome)
    s2_n_calls = len(s2_calls_in_chromosome)
    _match = 0

    for i, eventA in enumerate(s2_calls_in_chromosome):
        # for all events in the new "list"
        # find a match in the overall list, try to match or add

        matched = False
        for j, matchKey in enumerate(s1_calls_in_chromosome):
            _eventB = [x for x in s1_calls_in_chromosome[matchKey].values() if x.hexdigest == matchKey]
            if len(_eventB):
                eventB = _eventB[0]
            else:
                print([x.hexdigest for x in s1_calls_in_chromosome[matchKey].values()])
                print([x for x in s1_calls_in_chromosome[matchKey].values()])

                print(matchKey)
                print(_eventB)

            if not (eventA and eventB):
                continue
            if eventA == eventB:
                matchedInObject = matchKey
                # split out per chromosome storage
                allhits[virtualChromosome] = allhits.get(virtualChromosome, collections.OrderedDict())
                allhits[virtualChromosome][matchedInObject] = allhits[virtualChromosome].get(matchedInObject,
                                                                                                   collections.OrderedDict())
                allhits[virtualChromosome][matchedInObject][sampleName] = eventA
                _match += 1

                # when found, continue? It is the best approach to break this long list comparison?
                matched = True
                break

        if not matched:
            # add eventA to 'allhits[virtualChromosome]'
            if eventA.virtualChr == "chr17chr17":
                print("{} - {}".format(eventA.hexdigest, eventA))
            matchedInObject = eventA.hexdigest
            eventA.matched_in = matchedInObject
            virtualChromosome = eventA.virtualChr
            allhits[virtualChromosome] = allhits.get(virtualChromosome, collections.OrderedDict())
            allhits[virtualChromosome][matchedInObject] = allhits[virtualChromosome].get(matchedInObject,
                                                                                         collections.OrderedDict())
            allhits[virtualChromosome][matchedInObject][sampleName] = eventA
            # logger.debug("Adding new match in {} now {} items".format(
            #     virtualChromosome,
            #     len(allhits[virtualChromosome])
            # ))

    logger.debug("Common hits so far in {}: {} / {} vs {}".format(virtualChromosome, _match, s1_n_calls, s2_n_calls))
    return allhits[virtualChromosome]


def report(sampleNames, vcf_output_file, tsv_output_file, bed_output_file, commonhits, regions_out_file):
    # vcf header
    print(vcfHeader(sampleList=sampleNames), file=vcf_output_file)

    # tsv file
    samplecols = "\t".join(map(lambda x: "{}\tsize".format(x), sampleNames))
    header_line = "\t".join(['ChrA', 'ChrApos', 'ChrB', 'ChrBpos', 'SVTYPE', 'DP', 'Size', samplecols])

    tsv_report_output = open(tsv_output_file, 'w')
    tsv_report_output.write("{}\n".format(header_line))

    bed_structural_events = open(bed_output_file, 'w')
    all_locations = []

    for virtualChr in natsorted(commonhits.keys()):
        for sampleName, hits in sorted(commonhits[virtualChr].items(), key=lambda hit: hit[1].items()[0][1].chrApos):
            if len(hits):
                # check which samples has the same
                locations_found = []
                for sample in sampleNames:
                    if sample in hits.keys():
                        locations_found.append("{}\t{}".format(hits[sample], hits[sample].size))
                        # track all locations found for later intersecting or complementing the set of found/not-found
                        all_locations.append(hits[sample])
                    else:
                        locations_found.append("\t")
                # get the key with the highest DP
                sorted_by_dp = sorted(hits.items(), key=lambda hit: hit[1].dp, reverse=True)
                fKey = sorted_by_dp[0][0]
                outputRecord = hits[fKey]

                logger.debug(outputRecord)
                logger.debug(hits)

                print(formatVCFRecord(outputRecord, hits, sampleNames), file=vcf_output_file)

                tsv_report_output.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                    outputRecord.chrA,
                    outputRecord.chrApos,
                    outputRecord.chrB,
                    outputRecord.chrBpos,
                    outputRecord.sv_type,
                    outputRecord.dp,
                    outputRecord.size,
                    "\t".join(locations_found)))
                bed_structural_events.write(formatBedTrack(outputRecord))

    for loc in all_locations:
        print(loc.bedRow, file=regions_out_file)

    tsv_report_output.close()
    regions_out_file.close()

def startMerge(vcf_files, exclusion_regions, output_file, centerpointFlanking, bedoutput, transonly=False,
               regions_out="regions_out.bed", vcf_output="output.vcf"):
    regions_out_file = open(regions_out, "w")
    vcf_output_file = open(vcf_output, "w")

    # sampleList stores paths
    sampleList = vcf_files
    # sampleNames stores the real samplenames (now based on filenames only) Only single sample/vcf support now.
    # TODO: extract sample name from the samplename columns in the vcf
    sampleNames = [
        os.path.basename(x).replace(".vcf", "").replace(".realign", "").replace(".baserecal", "").replace(".dedup", "")
        for x in vcf_files]

    # sampleDB stores the file read handlers to the vcf files
    sampleDB = collections.OrderedDict()
    # svDB stores the actual content / parsed contents of the sv calls
    svDB = collections.OrderedDict()

    allhits = collections.OrderedDict()
    chromosomes_to_check = []


    edb = []
    if type(exclusion_regions) != type([]):
        exclusion_regions = []
    for exclusion_region in exclusion_regions:
        edb += build_exclusion(exclusion_region)

    for i, samplePath in enumerate(sampleList):
        sampleName = sampleNames[i]
        logger.info('Reading SV-events from sample: {} '.format(samplePath))
        sampleDB[sampleName] = vcf.Reader(open(samplePath, 'r'))

        # extract SV caller from header
        sv_caller = sampleDB[sampleName].metadata.get('source', [os.path.basename(samplePath).strip('.vcf')]).pop(
            0).split(' ').pop(0)

        svDB[sampleName] = loadEventFromVCF(sampleName, sampleDB[sampleName], edb, centerpointFlanking, transonly,
                                            sv_caller)
        n_events = sum([len(calls) for chromlist, calls in svDB[sampleName].items()])
        logger.info('Loaded SV-events from sample: {} '.format(n_events))

        # these are virtual chromosomes already
        chromosomes_to_check += svDB[sampleName].keys()
        chromosomes_to_check = list(set(chromosomes_to_check))
        chromosomes_to_check.sort()

        for _chromosome in chromosomes_to_check:
            s2_calls_in_chromosome = svDB[sampleName].get(_chromosome, [])
            allhits[_chromosome] = matchEvents(sampleName, s2_calls_in_chromosome, allhits, _chromosome)

        del svDB[sampleName]

    report(sampleNames, vcf_output_file, output_file, bedoutput, allhits, regions_out_file)

def main():
    parser = argparse.ArgumentParser()

    parser.add_argument('-c', '--exclusion_regions', action='append',
                        help='Exclusion regions file in BED format')

    parser.add_argument('-f', '--flanking', type=int,
                        help='Centerpoint flanking [100]', default=100)

    parser.add_argument('-s', '--sizeflanking', type=int,
                        help='Size Deviation [50]', default=50)

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
