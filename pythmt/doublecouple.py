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
# Last-Updated: Wed Dec 21 17:00:51 2011 (-0800)
#           By: Iain William Bailey
#     Update #: 102


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
from momenttensor import voigt2mat, mateig
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
    q0 = quat( dc0.pbt )  # make quaternions based on pbt as a rotation matrix
    q1 = quat( dc1.pbt )

    # get the difference
    qdiff = q1*q0.conjugate()   # equivelent to q1 then q0
    
    # return the four different rotations
    return qdiff, qdiff*quat([0,1,0,0]), qdiff*quat([0,0,1,0]), qdiff*quat([0,0,0,1]) 

##############################################################

# Python class representing a double-couple
class DoubleCouple:
    """
    A double couple defined by the orientation of its P and T axes.
    The default is located at [0,0,0] with p-axis in vertical(r)
    direction, t axis in east (phi) direction
    """

    # -------------------------------------------------- 
    # Default Constructor,
    def __init__( self 
                  , p = NP.array([1,0,0]) # vector for p-axis 
                  , t = NP.array([0,0,1]) # vector for t-axis
                  , strike = None, dip = None, rake = None # fault params in radians
                  , m0 = 1
                  , c = None # centroid location
                  , h = None  # hypocenter location
                  , mag = None
                  ):
        """
        Default constructor
        """

        self.EPS = 1e-8 # accuracy for checks

        # set orientation of double-couple
        if strike != None and dip != None and rake != None:
            # if the strike, dip and rake are set, use them
            (p, t) = sdr2ptax( strike, dip, rake)
        else:
            # check p and t are orthogonal
            if( NP.abs( NP.dot( p,t ) ) > self.EPS ):
                print >> sys.stderr, "P and T-axes are not orthogonal"
                sys.exit()

        b = NP.cross( t, p ) # make such that p x b = t 
        self.pbt = NP.c_[ p, b, t ]

        # set mag of double couple as scalar moment
        if mag != None: self.m0 = mag2m0( mag )
        else: self.m0 = m0

        # set position of double-couple
        if c == None and h == None :
            # default position 
            h = NP.array( [0.0, 0.0, 0.0] )
        if c == None: c = h
        if h == None: h = c
        self.c = c
        self.h = h

    # -------------------------------------------------- 
    def MTmat( self ):
        """
        Get the moment tensor in matrix form
        """
        U = NP.zeros( (3,3) )
        U[0,0] = -1*self.m0 # size of the p axis
        U[2,2] = self.m0 # size of t axis

        # p, t, b axes orientations are a rotation matrix
        return NP.dot( NP.dot( self.pbt, U ), self.pbt.transpose() )
 
    # -------------------------------------------------- 
    
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


