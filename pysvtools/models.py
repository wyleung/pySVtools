#!/usr/bin/env python2

from __future__ import print_function

__desc__= """"""
__author__ = "Wai Yi Leung <w.y.leung@lumc.nl>"


class Event(object):
    _centerpointFlanking = 100

    def __init__( self, chrA, chrApos, chrB, chrBpos, sv_type=None, cp_flank=None, dp=0 ):
        (self.chrA, self.chrApos), (self.chrB, self.chrBpos) = sorted( [(chrA, chrApos), (chrB, chrBpos)] )

        self.chrApos = int(self.chrApos)
        self.chrBpos = int(self.chrBpos)

        self.seen = False
        self.matched_in = None
        
        # number of reads supporting this breakpoint
        self.support = 0

        _vChr = [self.chrA, self.chrB]
        _vChr.sort()
        self.virtualChr = "".join( _vChr )

        self.centerpoint = self.get_centerpoint
        self.sv_type = sv_type
        self.centerpointFlanking = cp_flank or self._centerpointFlanking
        self.dp = dp
        self._hash = None
        

    @property
    def get_centerpoint(self):
        cnt = int(self.size/2)
        positions = [self.chrApos, self.chrBpos]
        positions.sort()
        centerpoint = positions[0] + cnt
        return centerpoint

    @property
    def size( self ):
        positions = [self.chrApos, self.chrBpos]
        positions.sort()
        return abs( positions[1] - positions[0] )

    def __eq__( self, other ):
        # if the current sv event is already visited, then skip it.
#       This predicate is False under the condition where we have more than 2 samples
#       This logic is now kept in comments, to show that this is wrong assumption
#        if self.seen or other.seen:
#            return False

#       #FIXME: We also need to check on the SVTYPE of each
#
#

#        if not ((self.chrA == other.chrA) and (self.chrB == other.chrB)):
#            return False
        # check whether we are working on the same chromosome set.
        if not ( self.virtualChr == other.virtualChr ):
            return False
    
        if abs(self.centerpoint - other.centerpoint) > self.centerpointFlanking:
            return False
    
        centerpointA = self.centerpoint
        centerpointB = other.centerpoint
        
        lftA = centerpointA - self.centerpointFlanking
        rgtA = centerpointA + self.centerpointFlanking

        lftB = centerpointB - self.centerpointFlanking
        rgtB = centerpointB + self.centerpointFlanking

        # rules defined: centerpoints should overlap each other by at least 1 bp
        # look for hypothesis first:
        # A l--------------r
        # B            l--------------r
        if lftB < rgtA and rgtB > rgtA:
            # we have a overlap
#            print "A end to B start"
            return True
        # A            l--------------r
        # B l--------------r
        if lftB < lftA and rgtB > lftA:
#            print "B end to A start"
            return True
            
        # conditions where A fully overlaps B or the otherway around
        # A l--------------------r
        # B     l----------r
        # or
        # A     l----------r
        # B l--------------------r
        if (lftB > lftA and rgtB < rgtA) or (lftB < lftA and rgtB > rgtA):
#            print "overlapping things"
            return True

        # say false for other not (yet) tested hypotheses
        return False

    def naiveEQ( self, other ):
        pass

    @property
    def vcf_alt( self ):
        if self.chrA == self.chrB and self.sv_type != 'ITX':
            return '.'
        # return for the CTX events
        return 'N[{}:{}['.format(self.chrB, self.chrBpos)

    @property
    def hexdigest( self ):
        if not self._hash:
            self._hash = hashlib.sha1( str(self) ).hexdigest()
        return self._hash

    def __repr__( self ):
        return "<Event {0}:{1}-{2}:{3}>".format( 
            self.chrA,
            self.chrApos,
            self.chrB,
            self.chrBpos
        )

    def __str__( self ):
        if self.chrA == self.chrB:
            fmtstring = "{0}:{1}-{3}"
        else:
            fmtstring = "{0}:{1}-{2}:{3}"
        return fmtstring.format(
            self.chrA,
            self.chrApos,
            self.chrB,
            self.chrBpos
        )

