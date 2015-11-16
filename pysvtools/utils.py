#!/usr/bin/env python2

import datetime
import re

import vcf.model
from pysvtools import __version__
from pysvtools.models.exclusionregion import ExclusionRegion


def extractTXmate(record):
    """
    :param record: pyVCF.vcf.model._Record
    :return: [ chrB, chrBpos ]: Returning mate chromosome if possible, otherwise [None, None]
    """
    chrB = None
    chrBpos = None

    # first check in record.alt
    try:
        alt = record.ALT.pop()
    except:
        alt = ""

    if type(alt) == vcf.model._SV:
        # try to extract from `record.INFO` on the keys: [CHR2, END]
        chrB = record.INFO.get('CHR2', None)
        chrBpos = record.sv_end
    elif type(alt) == vcf.model._Breakend:
        chrB = alt.chr
        chrBpos = alt.pos
    return [chrB, chrBpos]


def vcfHeader():
    ts_now = datetime.datetime.now()
    vcf_header = """##fileformat=VCFv4.1
##fileDate={filedate}
##source=pysvtools-{version}
##ALT=<ID=DEL,Description="Deletion">
##ALT=<ID=DUP,Description="Duplication">
##ALT=<ID=INS,Description="Insertion">
##ALT=<ID=INV,Description="Inversion">
##ALT=<ID=TRA,Description="Translocation">
##FILTER=<ID=LowQual,Description="PE support below 3 or mapping quality below 20.">
##FORMAT=<ID=DR,Number=1,Type=Integer,Description="# high-quality reference pairs">
##FORMAT=<ID=DV,Number=1,Type=Integer,Description="# high-quality variant pairs">
##FORMAT=<ID=FT,Number=1,Type=String,Description="Per-sample genotype filter">
##FORMAT=<ID=GL,Number=G,Type=Float,Description="Log10-scaled genotype likelihoods for RR,RA,AA genotypes">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification">
##FORMAT=<ID=RC,Number=1,Type=Integer,Description="Normalized high-quality read count for the SV">
##FORMAT=<ID=RR,Number=1,Type=Integer,Description="# high-quality reference junction reads">
##FORMAT=<ID=RV,Number=1,Type=Integer,Description="# high-quality variant junction reads">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
##INFO=<ID=CHR2,Number=1,Type=String,Description="Chromosome for END coordinate in case of a translocation">
##INFO=<ID=CIEND,Number=2,Type=Integer,Description="PE confidence interval around END">
##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="PE confidence interval around POS">
##INFO=<ID=CONSENSUS,Number=1,Type=String,Description="Split-read consensus sequence">
##INFO=<ID=CT,Number=1,Type=String,Description="Paired-end signature induced connection type">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the structural variant">
##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variation">
##INFO=<ID=MAPQ,Number=1,Type=Integer,Description="Median mapping quality of paired-ends">
##INFO=<ID=PE,Number=1,Type=Integer,Description="Paired-end support of the structural variant">
##INFO=<ID=PRECISE,Number=0,Type=Flag,Description="Precise structural variation">
##INFO=<ID=SR,Number=1,Type=Integer,Description="Split-read support">
##INFO=<ID=SRQ,Number=1,Type=Float,Description="Split-read consensus alignment quality">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length of the SV">
##INFO=<ID=SVMETHOD,Number=1,Type=String,Description="Type of approach used to detect SV">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">""".format(filedate=ts_now.strftime("%Y%m%d"),
                                                                          version=__version__)
    return vcf_header + "\n" + "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	default"


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


def extractDPFromRecord(record):
    if 'DP' in record.INFO.keys():
        if type(record.INFO['DP']) == type(list):
            return record.INFO['DP'][0]
        return record.INFO['DP']
    elif len(record.samples):
        return getattr(record.samples[0].data, 'DP', 0)
    return 0


def firstFromList(arr):
    if type(arr) == type([]):
        return arr[0]
    return arr


def getSVType(record):
    if type(record.INFO['SVTYPE']) == type([]):
        SVTYPE = record.INFO['SVTYPE'][0]
    else:
        SVTYPE = record.INFO['SVTYPE']
    return SVTYPE


def getSVLEN(record):
    SVLEN = 0
    if "SVLEN" in record.INFO.keys():
        if type(record.INFO['SVLEN']) == type([]):
            SVLEN = record.INFO['SVLEN'][0]
        elif type(record.INFO['SVLEN']) == type(0):
            SVLEN = record.INFO['SVLEN']
    return SVLEN


def formatBedTrack(mergedHit):
    formatted_bed = ""

    if mergedHit.sv_type in ["INS", 'DEL', 'ITX']:
        formatted_bed = "{chrom}\t{start}\t{end}\t{annot}\n".format(
            chrom=mergedHit.chrA,
            start=mergedHit.chrApos,
            end=mergedHit.chrBpos,
            annot="SVTYPE={svtype};DP={dp};SIZE={size}".format(
                svtype=mergedHit.sv_type,
                dp=mergedHit.dp,
                size=mergedHit.size
            )
        )
    else:
        # for CTX events, we write 2 tracks, one for the left breakpoint and one for th right breakpoint
        # both sites are flanked with 50 bases to generate a valid bed file and used in IGV to mark the region
        formatted_bed = "{chrom}\t{start}\t{end}\t{annot}\n".format(
            chrom=mergedHit.chrA,
            start=(mergedHit.chrApos > 50) and (mergedHit.chrApos - 50) or mergedHit.chrApos,
            end=mergedHit.chrApos + 50,
            annot="SVTYPE={svtype};DP={dp};SIZE={size};MATE={chrb}:{chrbpos}".format(
                svtype=mergedHit.sv_type,
                dp=mergedHit.dp,
                size=mergedHit.size,
                chrb=mergedHit.chrB,
                chrbpos=mergedHit.chrBpos
            )
        )
        formatted_bed += "{chrom}\t{start}\t{end}\t{annot}\n".format(
            chrom=mergedHit.chrB,
            start=mergedHit.chrBpos - 50,
            end=mergedHit.chrBpos + 50,
            annot="SVTYPE={svtype};DP={dp};SIZE={size};MATE={chrb}:{chrbpos}".format(
                svtype=mergedHit.sv_type,
                dp=mergedHit.dp,
                size=mergedHit.size,
                chrb=mergedHit.chrA,
                chrbpos=mergedHit.chrApos
            )
        )
    return formatted_bed


def formatVCFRecord(mergedHit):
    # TODO: write the DP for each of the callers/sample
    INFOFIELDS = "IMPRECISE;SVTYPE={};CHR2={};END={};SVMETHOD={svmethod}".format(
        mergedHit.sv_type,
        mergedHit.chrB,
        mergedHit.chrBpos,
        svmethod=mergedHit.svmethod
    )
    FORMATFIELDS = ":".join(map(str, [
        '1/.',
        mergedHit.dp]))

    formattedVCFRecord = "{chrA}\t{pos}\t{id}\t{ref}\t<{alt}>\t{qual}\t{filter}\t{info}\tGT:DP\t{format}".format(
        chrA=mergedHit.chrA,
        pos=mergedHit.chrApos,
        id='.',
        ref='N',
        alt=mergedHit.sv_type,
        qual='100',
        filter='PASS',
        info=INFOFIELDS,
        format=FORMATFIELDS
    )
    return formattedVCFRecord


def build_exclusion(bed_exclude=None):
    exclusiondb= []
    with open(bed_exclude, 'r') as fd:
        for r in fd:
            row = r.strip().split("\t")
            cols = dict(zip(['chromosome', 'start', 'end', 'band', 'color'], row))
            exclusion = ExclusionRegion(**cols)
            exclusiondb.append(exclusion)
    return exclusiondb
