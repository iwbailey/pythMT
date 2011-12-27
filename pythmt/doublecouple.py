#!/usr/bin/env python
# doublecouple.py --- 
# 
# Filename: doublecouple.py
#
# Description: Python class for double-couple moment tensor, derived
# class of SymMT (see MomentTensor.py)
#
# Author: IW Bailey
# Maintainer: IW Bailey
# Created: Mon Dec 19 11:30:37 2011 (-0800)
# Version: 1
# Last-Updated: Tue Dec 27 14:27:17 2011 (-0800)
#           By: Iain William Bailey
#     Update #: 151


# Change Log:
#
#
#
#
#

# standard libraries
import numpy as NP
from numpy import sin, cos
from math import sqrt
# import sys

# other libs, must be in path
from momenttensor import EigMT, voigt2mat, mateig, mag2m0
from quaternions import quatn as quat

# FUNCTIONS ##################################################

def sdr2mtvect( strike, dip, rake ):
    """
    Convert strike, dip, rake (radians) into mt...
    mt = [ mrr, mtt, mpp, mtp, mrp, mrt ] where r=UP, t=S, p=E  
    """
    is2 = 1/sqrt(2.0)

    mrr =  is2*( sin(2*dip) * sin(rake) );

    mtt = -is2*( sin(dip)  * cos(rake) * sin(2*strike) +     
                 sin(2*dip) * sin(rake) * sin(strike)**2 );

    mpp =  is2*( sin(dip)  * cos(rake) * sin(2*strike) -
                 sin(2*dip) * sin(rake) * cos(strike)**2 );

    mtp = -is2*( sin(dip)  * cos(rake) * cos(2*strike) + 
                 0.5*sin(2*dip) * sin(rake) * sin(2*strike) );

    mrp =  is2*( cos(dip)  * cos(rake) * sin(strike)  -
                 cos(2*dip) * sin(rake) * cos(strike));

    mrt = -is2*( cos(dip)  * cos(rake) * cos(strike)  +    
                 cos(2*dip) * sin(rake) * sin(strike) );

    return NP.array([ mrr, mtt, mpp, mtp, mrp, mrt ] )

#------------------------------------------------------------
def sdr2smt( strike, dip, rake ):
    """
    Convert strike, dip, rake (radians) into unit norm mt matrix...
    mt = [ mrr, mrt, mrp ;
           mrt, mtt, mtp ;
           mrp, mtp, mpp ],  where r=UP, t=S, p=E  
    """
    # get in vect form
    mtvect = sdr2mtvect( strike, dip, rake )

    # convert to matrix form
    return voigt2mat( mtvect )


#------------------------------------------------------------
def sdr2ptax( strike, dip, rake ):
    """
    Convert strike, dip, rake (radians) into p and t-axis vectors...
    p = [ pr, pt, pp ]
    t = [ tr, tt, tp ] where r=UP, t=S, p=E  
    """
    # get unit norm matrix
    mt = sdr2smt( strike, dip, rake )
    
    # get eigen solution
    (ptb, vals) = mateig( mt )

    paxis = ptb[:,0]
    taxis = ptb[:,2]
    return (paxis, taxis)

#------------------------------------------------------------
def dcquat( dc0, dc1 ):
    """
    Get the quaternions describing rotation from one double couple to
    another.  dc0 is the one we're rotating from.
    """

    q0 = quat( A=dc0.pbt )  # make quaternions based on pbt as a rotation matrix
    q1 = quat( A=dc1.pbt )

    # quaternion for rotation from dc0 to dc1 is the same as inverse of dc0 then dc1
    qdiff = q1*q0.conjugate()   # equivelent to q1 then q0
    
    # return the four different rotations
    return qdiff, qdiff*quat([0,1,0,0]), qdiff*quat([0,0,1,0]), qdiff*quat([0,0,0,1]) 

##############################################################

# Python class representing a double-couple, inherited from moment tensor
class DoubleCouple( EigMT ):
    """
    A double couple defined by the orientation of its P and T axes.
    The default is located at [0,0,0] with p-axis in vertical(r)
    direction, t axis in east (phi) direction
    """

    # -------------------------------------------------- 
    # Default Constructor,
    def __init__( self 
                  , p = NP.array([1,0,0], dtype=float ) # vector for p-axis 
                  , t = NP.array([0,0,1], dtype=float ) # vector for t-axis
                  , strike = None, dip = None, rake = None # fault params in radians
                  , m0 = 1
                  , c = None # centroid location
                  , h = None  # hypocenter location
                  , mag = None
                  ):
        """
        Constructor for double-couple class
        """

        # set orientation of double-couple from strike, dip and rake
        if strike != None and dip != None and rake != None:
            # if the strike, dip and rake are set, use them
            (p, t) = sdr2ptax( strike, dip, rake)

        # set mag of double couple as scalar moment
        if mag != None: self.m0 = mag2m0( mag )
        else: self.m0 = m0

        # add the double couple constraint to the moment tensor rep
        pbtvals = NP.array( [ -m0, 0.0, m0 ] )

        # initialize the moment tensor part
        EigMT.__init__( self, p=p, t=t, pbtvals=pbtvals, c=c, h=h )

    # -------------------------------------------------- 
    def quats( self ):
        """
        Get the 4 quaternions describing rotation from the reference
        orientation where p = (1,0,0), t=(0,0,1)
        """
        # get direct mapping from rotation matrix
        q0 = quat( A=self.pbt )

        # get other four options by quaternion multiplication
        q1, q2, q3 = q0*quat([0,1,0,0]), q0*quat([0,0,1,0]), q0*quat([0,0,0,1]) 

        # force the scalar part to be positive
        for q in ( q0, q1, q2, q3 ):
            if q.w < 0: q*= -1

        return q0, q1, q2, q3

##############################################################
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License as
# published by the Free Software Foundation; either version 3, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; see the file COPYING.  If not, write to
# the Free Software Foundation, Inc., 51 Franklin Street, Fifth
# Floor, Boston, MA 02110-1301, USA.
#
# doublecouple.py ends here


