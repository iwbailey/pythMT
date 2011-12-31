#!/usr/bin/env python
# summed_mt.py --- 
# 
# Filename: summed_mt.py
#
# Description: Python class for summed moment tensor, derived class of
# SymMT (see momenttensor.py)
#
# Author: IW Bailey
# Maintainer: IW Bailey
# Created: Fri Nov  5 10:06:42 2010 (-0700)
# Version: 1
# Last-Updated: Fri Dec 30 16:41:09 2011 (-0800)
#           By: Iain Bailey
#     Update #: 181

# Change Log:
#
#
#
#

# standard libraries
import numpy as NP
from math import sqrt

# other libs, must be in path
from momenttensor import SymMT
#import sys

# Python class representing a summed moment tensor, derived class of SymMT
class MTsum(SymMT):
    """
    Summed moment tensor object, basically contains a symmetric moment
    tensor with extra functions for adding up
    """

    # --------------------------------------------------
    # Default Constructor: create empty tensor
    def __init__( self  ):
        SymMT.__init__(self)
        self.cumNorm = 0.0
        self.count = 0.0 # number of tensors

    # --------------------------------------------------
    def add( self, newmt ):       
        """
        Add a tensor to the sum
        """
        summt = self.mhat*self.Norm + newmt.mhat*newmt.Norm

        # update the location. 
        self.c = ( self.cumNorm*self.c + 
                   newmt.Norm*newmt.c )/(self.cumNorm + newmt.Norm)

        self.h = ( self.cumNorm*self.h + 
                   newmt.Norm*newmt.h )/(self.cumNorm + newmt.Norm)

        # update the norm
        self.Norm = sqrt( NP.sum( summt*summt ) + NP.sum( summt[2:6]*summt[2:6] ) )
        self.cumNorm += newmt.Norm
         
        # update the MT
        self.mhat = summt / self.Norm
        
        # update number
        self.count += 1;


# #
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
# summed_mt.py ends here


