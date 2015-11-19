#!/usr/bin/env python2


class ExclusionRegion(object):
    def __init__(self, chromosome, start, end, *args, **kwargs):
        self.chromosome = chromosome
        self.start = int(start)
        self.end = int(end)

    def overlaps(self, qChr, qPos):
        # simple check on chromosome
        return qChr == self.chromosome and self.start <= qPos <= self.end

    def __repr__(self):
        return "<ExclusionRegion {chromosome}:{start}-{end}>".format(
            chromosome=self.chromosome,
            start=self.start,
            end=self.end)
