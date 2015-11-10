import re


class Intersection(object):
    def __init__(self):
        pass


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
