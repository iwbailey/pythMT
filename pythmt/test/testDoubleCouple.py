#!/usr/bin/env python

# check that the double couple is working fine
 
import sys
import numpy as np
from math import pi
from pythmt.MomentTensorDC import DoubleCouple  


# degrees to radians conversion
D2R = pi/180.0;

np.set_printoptions(precision=3, suppress=True)

print "Testing the default constructor"
dc0 = DoubleCouple()
print "pb = "
print dc0.pbt
print ""

print "Testing the conversion from strike, dip and rake"
stk = [0,   0, 45 ]
dip = [90, 45, 90 ]
rak = [0, -90,  0 ]

n = len(stk)

for i in range(0,n):
    print( "Strike: %.1f^o, Dip: %.1f^o, Rake: %.1f^o" % 
           (stk[i], dip[i], rak[i]) )
    dc = DoubleCouple( strike=D2R*stk[i], dip=D2R*dip[i], rake=D2R*rak[i] )
    print "pbt = "
    print dc.pbt
    print "moment tensor = "
    print dc.MTmat()

###########################################################
