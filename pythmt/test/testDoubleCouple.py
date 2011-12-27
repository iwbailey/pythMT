#!/usr/bin/env python

# check that the double couple is working fine
 
import os
import sys
import numpy as np
from math import pi

# add next directory up to path
srcpath = os.path.join( os.getcwd(), '..' )
sys.path.append( srcpath )

# import module to test
import doublecouple 
from doublecouple import DoubleCouple as DC

# degrees to radians conversion
D2R = pi/180.0;

np.set_printoptions(precision=3, suppress=True)

print "Testing the default constructor"
dc0 = DC()
print "pbt = "
print dc0.pbt
print ""

#strike, dip, rake = 0, 45, -90 
strike, dip, rake = 10, 70, 25 
print( "Strike: %.1f^o, Dip: %.1f^o, Rake: %.1f^o" %  (strike, dip, rake) )
dc = DC( strike=D2R*strike, dip=D2R*dip, rake=D2R*rake )
print "pbt:"
print dc.pbt

# test the quaternion representation
(q1, q2, q3, q4)  = dc.quats()
print "q:" ,q1, q2, q3, q4

# test the moment tensor output
print "Mij:"
print dc.M()

print""

sys.exit()


###########################################################
