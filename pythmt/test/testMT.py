#!/usr/bin/env python

# test that the moment tensor class is working by running some examples

import os
import sys
import numpy as np

# get file to test
srcpath = os.path.join( os.getcwd(), '..' )
sys.path.append( srcpath )
import momenttensor as mt

# # print documentation
# print mt.SymMT.__doc__

# define a standard double couple
mij = np.array( [1, 0, -1, 0, 0, 0] )
print "input:", mij
mij = mt.SymMT( mij )

print "smt\n", mij.smt()
print "Moment:", mij.M0()
print "Mw", mij.Mw()
print "eig:", mij.pbt()

# test default constructor
print "\ntest default constructor..."
mij = mt.SymMT()

print "smt\n", mij.smt()
print "Moment:", mij.M0()
print "Mw", mij.Mw()
print "eig:", mij.pbt()
