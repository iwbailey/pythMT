#!/usr/bin/env python

# Output info about moment tensors to check our moment tensor class is working
 
import sys
from MomentTensor import readpsmecaSm  

Ndata = 0
# read in from stdin assuming psmeca format
while 1:
    thisline = sys.stdin.readline()

    if thisline != '':

        # Read line as SymMT object
        MT, extra = readpsmecaSm( thisline )
        Ndata += 1

        print "\nData point Number %i" % Ndata

        print "Longitude: \t%.4f^o" % MT.c[0]
        print "Latitude: \t%.4f^o" % MT.c[1]
        print "Depth: \t\t%.2f km" % MT.c[2]
        print "Scalar Moment: \t%.4e dyne cm" % MT.getM0()
        print "Mw: \t\t%.2f" % MT.getMw()
        print "Trace: \t\t%.4e dyne cm" % ((MT.M(1,1) + MT.M(2,2) + MT.M(3,3)))
        print "f_CLVD: \t%.2f " % MT.getFclvd()

        # get faulting style
        [ptb, vals] = MT.getEig()
        if( abs( ptb[2,0] ) > abs( ptb[2,1] ) and abs( ptb[2,0] ) > abs( ptb[2,2] ) ):
            mostVert='P'
        elif( abs( ptb[2,1] ) > abs( ptb[2,2] ) ):
            mostVert = 'B'
        else:
            mostVert = 'T'
        print "Most Vertical Axis: %c" % mostVert

        print "in:  %s" % thisline.rstrip('\n')
        (x,y,z,mrr,mtt,mff, mrt, mrf, mtf, exp) = MT.getPsmecaSm()
        print "out: %8.3f %8.3f %8.2f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %2i " % \
                                         (x,y,z,mrr,mtt,mff, mrt, mrf, mtf,exp)
        print "Extra info not processed: %s" % extra


    else: break
            
###########################################################
