import re


class Intersection(object):
    def __init__(self):
        pass


def extractTXmate(record):
    """
    :param record: pyVCF.vcf.model._Record
    :return: [ chrB, chrBpos ]: Returning mate chromosome if possible, otherwise [None, None]
    """
    chrB = None
    chrBpos = None
    # first check in record.alt
    alt = record.alt.pop()

    try:
        res = re.findall(r"([\d\w\_]+)\:([\d]+)", alt, re.I | re.M)[0]
    except:
        # try to extract from `record.INFO` on the keys: [CHR2, END]
        chrB = record.INFO.get('CHR2', None)
        chrBpos = record.INFO.get('END', None)
    else:
        chrB = res[0]
        chrBpos = res[1]

    return [chrB, chrBpos]


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
            chrom=t.chrA,
            start=t.chrApos,
            end=t.chrBpos,
            annot="SVTYPE={svtype};DP={dp};SIZE={size}".format(
                svtype=t.sv_type,
                dp=t.dp,
                size=t.size
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
        formatted_bed += "\n{chrom}\t{start}\t{end}\t{annot}\n".format(
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
    INFOFIELDS = "IMPRECISE;SVTYPE={};END={}".format(
        mergedHit.sv_type,
        mergedHit.chrBpos
    )
    FORMATFIELDS = ":".join(map(str, [
        '1/.',
        mergedHit.dp]))
    formattedVCFRecord = "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\tGT:DP\t{}".format(
        mergedHit.chrA,
        mergedHit.chrApos,
        '.',
        '.',
        mergedHit.vcf_alt,
        '100',
        'PASS',
        INFOFIELDS,
        FORMATFIELDS
    )
    return formattedVCFRecord
