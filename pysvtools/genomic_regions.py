#!/usr/bin/env python2

import vcf

class ExclusionRegion(object):
    def __init__( self, chromosome, start, end, *args, **kwargs ):
        self.chromosome = chromosome
        self.start = start
        self.end = end

    def overlaps( self, qChr, qPos ):
        # simple check on chromosome
        if qChr != self.chromosome:
            return False
        if qPos < self.start:
            return False
        if qPos > self.end:
            return False

        if qPos > self.start and qPos < self.end:
            return True

        # to be safe for any other combinations .
        return False

    def __eq__(self, other):
        # a little bit different, abuse the eq function to see whether a region
        # overlaps
        pass

    def __repr__( self ):
        return "<ExclusionRegion {}:{}-{}>".format(self.chromosome, self.start, self.end)


def buildExclusion( bed_exclude=None ):
    exclusionDB=[]
    with open( bed_exclude, 'r') as fd:
        for r in fd:
            row = r.strip().split("\t")
            cols = dict(zip(['chromosome', 'start', 'end', 'band', 'color' ], row))
            exclusion = ExclusionRegion( **cols )
            exclusionDB.append(exclusion)
    return exclusionDB


