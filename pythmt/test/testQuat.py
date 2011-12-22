#!/usr/bin/env python

import sys
import numpy as np
from numpy import pi, cos, sin
from quaternion import quatn as quat

# from cgkit.cgtypes import quat as cgquat
D2R = pi/180.0

q1 = [ 1/np.sqrt(2), 0, 1/np.sqrt(2), 0 ]
q2 = [ np.cos(0.5*45*D2R) , 0, 0, np.sin(0.5*45*D2R) ]

q1 = quat( q1 )
q2 = quat( q2 )

print "q1 = ", q1
print "q2 = ", q2
print "rotation by q2 then q1, i.e., q1*q2 = ", q1*q2
print "q2*q1.conjugate() = ", q2*q1.conjugate()

# case where reference dc is 1,1,1
dc1 = np.array( [[ 1, 0, 0 ],
                 [ 0, 1, 0 ],
                 [ 0, 0, 1 ]] )
#print quat( A=dc1 )
print ""

# dc1 = np.array( [[ 1, 0, 0 ],
#                  [ 0, 0, 1 ],
#                  [ 0, 1, 0 ]] , dtype=float)
dc1 = np.array( [[ 0, 0, 1 ],
                 [ 1, 0, 0 ],
                 [ 0, 1, 0 ]] )
q = quat( A=dc1 )
print "q=", q
q=q.normalize()
print q
(a, u) = q.toAngleAxis()
print "a =", a, "u=",u
print q.toMat3()

# print quat( A=dc1 )



# print cgquat( q1 ).toMat3()
# print quat( q2 ).toMat3()
sys.exit()
