#!/usr/bin/env python

# Output info about moment tensors to check our moment tensor class is working
 
import os
import sys

# get file to test
srcpath = os.path.join( os.getcwd(), '..' )
sys.path.append( srcpath )
from iofuncs import readPsmecaList

ifname = "../../sample_data/socal_cmts.psmeca"

mtlist, alltxt = readPsmecaList( ifname ) 

nmt = len(mtlist)
print "read ",nmt, "tensors"


# loop through each 
for i in range(0,nmt):
    MT = mtlist[i]
    print "\n"
    print "Line Number: \t%i" % (i+1) 
    print "Longitude: \t%.4f^o" % MT.c[0]
    print "Latitude: \t%.4f^o" % MT.c[1]
    print "Depth: \t\t%.2f km" % MT.c[2]
    print "Scalar Moment: \t%.4e dyne cm" % MT.M0()
    print "Mw: \t\t%.2f" % MT.Mw()
    print "f_CLVD: \t%.2f" % MT.f_clvd()
    print "Trace/3: \t%.4e dyne cm" % ((MT.M(1,1) + MT.M(2,2) + MT.M(3,3))/3)
   


