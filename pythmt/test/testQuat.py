#!/usr/bin/env python

import os 
import sys
import numpy as np
from numpy import pi, cos, sin, arccos as acos

# get file to test
srcpath = os.path.join( os.getcwd(), '..' )
sys.path.append( srcpath )
from quaternions import quatn as quat

D2R = pi/180.0
R2D = 180/pi

# reference 
A0 = np.array( [[ 1, 0, 0 ],
                [ 0, 1, 0 ],
                [ 0, 0, 1 ]], dtype=float )

print "reference rotation matrix"
print A0
q0 = quat( A=A0 )
print "q = ", q0
print ""

# make x1->-x3, x3->x1
A1 = dc1 = np.array( [[ 0,  0, 1 ],
                      [ 0,  1, 0 ],
                      [ -1, 0, 0 ]] , dtype=float)

print "posn 1"
print A1
q1 = quat( A=A1 )
print "q1 = ", q1
(a, u) = q1.toAngleAxis()
print "angle", a*R2D, "ACW around axis",u
print ""

# make x1->-x3, x2->-x1, x3->x2
A2 = dc1 = np.array( [[  0, -1, 0 ],
                      [  0,  0, 1 ],
                      [ -1,  0, 0 ]] , dtype=float)

print "posn 2"
print A2
q2 = quat( A=A2 )
print "q2 = ", q2
(a, u) = q2.toAngleAxis()
print "angle", a*R2D, "ACW around axis",u
print ""

#  Combined rotation
print "q2 then q1,  q1*q2"
q3 = q1*q2
print "q1*q2 = ", q3
(a, u) = q3.toAngleAxis()
print "angle", a*R2D, "ACW around axis",u
print ""

#  Combined rotation
print "q1 then q2,  q2*q1"
q3 = q2*q1
print "q2*q1 = ", q3
(a, u) = q3.toAngleAxis()
print "angle", a*R2D, "ACW around axis",u
print ""

# rotation between
print "from q1 to q2"

# equivalent to rotating back to the reference then doing q2
q3 = q2*q1.conjugate()
print "q2*(q1^-1) = ", q3
(a, u) = q3.toAngleAxis()
print "angle", a*R2D, "ACW around axis",u
print ""

print "from q2 to q1"

q3 = q1*q2.conjugate()
print "(q2^-1)*q1 = ", q3
(a, u) = q3.toAngleAxis()
print "angle", a*R2D, "ACW around axis",u
print ""

sys.exit()

q1 = [ 1/np.sqrt(2), 0, 1/np.sqrt(2), 0 ]
q2 = [ np.cos(0.5*45*D2R) , 0, 0, np.sin(0.5*45*D2R) ]

q1 = quat( q1 )
q2 = quat( q2 )

print "q1 = ", q1
print "q2 = ", q2
print "rotation by q2 then q1, i.e., q1*q2 = ", q1*q2
print "q2*q1.conjugate() = ", q2*q1.conjugate()

qdiff = q2*q1.conjugate()
print qdiff.conjugate() * q2
print qdiff * q1

sys.exit()
# case where reference dc is 1,1,1
dc1 = np.array( [[ 1, 0, 0 ],
                 [ 0, 1, 0 ],
                 [ 0, 0, 1 ]] )
q=quat( A=dc1 )
print "from ",dc1," to", q
(a, u) = q.toAngleAxis()
print "angle", a*R2D, "ACW around axis",u
print ""

# 90 degree rotations

# around x axis
dc1 = np.array( [[ 1, 0, 0 ],
                 [ 0, 0, 1 ],
                 [ 0, -1, 0 ]] , dtype=float)

# around y axis
# dc1 = np.array( [[ 0, 0, 1 ],
#                  [ 0, 1, 0 ],
#                  [ 1, 0, 0 ]] , dtype=float)

# around z axis 
# dc1 = np.array( [[ 0, -1, 0 ],
#                  [ 1, 0, 0 ],
#                  [ 0, 0, 1 ]] , dtype=float)

# 120 deg rotations
# dc1 = np.array( [[ 0, 0, 1 ],
#                  [ 1, 0, 0 ],
#                  [ 0, 1, 0 ]] )
q = quat( A=dc1 )
print "from ",dc1," to", q
(a, u) = q.toAngleAxis()
print "angle", a*R2D, "ACW around axis",u
print ""

sys.exit()
