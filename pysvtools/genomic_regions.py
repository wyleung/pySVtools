#!/usr/bin/env python2

class ExclusionRegion(object):
    def __init__( self, chromosome, start, end, *args, **kwargs ):
        self.chromosome = chromosome
        self.start = int(start)
        self.end = int(end)

    def overlaps( self, qChr, qPos ):
        # simple check on chromosome
        return (qChr == self.chromosome and qPos >= self.start and qPos <= self.end)

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


