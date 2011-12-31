# iofuncs.py --- 
# 
# Filename: iofuncs.py
# Description: 
#
# Functions for related input and output in this package
#
# Author: Iain William Bailey
# Created: Wed Dec 21 10:00:03 2011 (-0800)
# Version: 1
# Last-Updated: Fri Dec 30 16:54:11 2011 (-0800)
#           By: Iain Bailey
#     Update #: 288

# Change Log:
# 
# 
# 
# ***********License stuff at bottom ************

# Code:

import sys
import numpy as NP
from math import sqrt, pi

from momenttensor import SymMT, EigMT
from doublecouple import DoubleCouple, sdr2smt

#--------------------------------------------------
def psmeca2SymMT( string ):
    """
    Convert from a psmeca -Sm format string to a SymMT object
    """ 

    # convert first 12 columns into array
    mtline = NP.fromstring( string, count=12, sep =' ', dtype=float )
    
    # get the location part
    c = mtline[0:3] # assume lon/lat/depth are centroid
    h = NP.array( [mtline[10], mtline[11], mtline[2]] ) # assume second lon/lat are hypocenter

    # get moment tensor part, voigt notation 
    # [mrr mtt mpp mrt mrp mtp] -> [mrr mtt mpp mtp mrp mrt]
    mhat = NP.r_[mtline[3:6], mtline[8], mtline[7], mtline[6]] 
    exp = mtline[9]
    
    # remove the norm, take care with off diagonals
    norm = NP.sqrt( NP.sum( mhat**2 ) + NP.sum( mhat[3:]**2 ) )
    mhat /= norm

    # get the true norm
    norm *= 10**exp

    # make the moment tensor object
    return SymMT( mhat=mhat, Norm=norm, c=c, h=h )
    
#--------------------------------------------------
def psmeca2EigMT( string ):
    """
    Convert from a psmeca -Sm format string to a EigMT object
    """ 
    
    # convert first 12 columns into array
    mtline = NP.fromstring( string, count=12, sep =' ', dtype=float )
    
    # get the location part
    c = mtline[0:3] # assume lon/lat/depth are centroid
    h = NP.array( [mtline[10], mtline[11], mtline[2]] ) # assume second lon/lat are hypocenter

    # get moment tensor part, voigt notation
    exp = mtline[9]
    # [mrr mtt mpp mrt mrp mtp] -> [mrr mtt mpp mtp mrp mrt]
    m = ( NP.r_[mtline[3:6], mtline[8], mtline[7], mtline[6]] )* 10**exp

    # make the moment tensor object
    return EigMT( m=m, c=c, h=h )


#--------------------------------------------------
def read_psmecalist( istream , isEig=False ):
    """
    From an input file or stdin, read a list of moment tensors in psmeca form

    Expected format
    lon/lat/z/mrr/mtt/mpp/mrt/mrp/mtp/exp/lon0/lat0/str/anything else

    if isEig is false, return list of SymMT objects, otherwise get a
    list of EigMT objects
    """

    mtlist=[] # this will be the output list

    # read everything
    alltxt = NP.genfromtxt( istream, delimiter='\n' , dtype=str)
    try: 
        istream.close()
    except:
        tmp=1

    # loop through all tensors
    n = len(alltxt)

    # check for desired output type
    if isEig:
        for i in range(0,n):
            mtlist.append( psmeca2EigMT( alltxt[i] ) )
    else:
        for i in range(0,n):
            mtlist.append( psmeca2SymMT( alltxt[i] ) )

    
    return mtlist, alltxt

#--------------------------------------------------
def sdr2dc( string ):
    """
    Convert from a lon/lat/depth/strike/dip/rake/mag string to a double couple type
    """
    # convert first 7 columns into array
    try:
        data = NP.fromstring( string, count=7, sep =' ', dtype=float )
    except IndexError:
        print >> sys.stderr, "Error: Require 7 columns: lon/lat/depth/strike/dip/rake/mag"
        return None

    # get variables
    hypo = data[0:3 ] # lon , lat, z
    sdr = data[3:6]*pi/180 # in radians
    mag = data[6]

    return DoubleCouple( c=hypo, h=hypo,
                         strike = sdr[0], dip = sdr[1], rake=sdr[2],
                         mag=mag )

#--------------------------------------------------
def read_sdrlist( istream ):
    """
    From an input file or stdin, read a list of focal mechanisms in
    strike/dip/rake form

    lon/lat/depth/strike/dip/rake/mag
    """

    mtlist = []

    # read everything
    alltxt = NP.genfromtxt( istream, delimiter='\n' , dtype=str)
    try: 
        istream.close()
    except:
        tmp=1

    # loop through all events
    n = len(alltxt)
    for i in range(0,n):
        mtlist.append( sdr2dc( alltxt[i] ) )

    return mtlist, alltxt
    
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
# 
# iofuncs.py ends here
