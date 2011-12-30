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
# Last-Updated: Fri Dec 30 12:14:57 2011 (-0800)
#           By: Iain William Bailey
#     Update #: 231

# Change Log:
# 
# 
# 
# ***********License stuff at bottom ************

# Code:

import sys
import numpy as NP
from math import sqrt, pi

from momenttensor import SymMT
from doublecouple import DoubleCouple, sdr2smt

# #--------------------------------------------------
# def readpsmecaSm( thisline , lcount=1):
#     """
#     Take a string as an input argument.  The string should be a psmeca
#     entry.  Return a moment tensor object

#     expected input file: 
#     x, y, z, mrr, mtt, mff, mrt, mrf, mtf, exp
#     r = up, t is south, f is east

#     the last three parts (lon', lat', id) and anything after are
#     output as a string
#     """

#     # Get matrix of floats without spaces
#     # this will crash or get things wrong if you have tabs ('\t') with no spaces
#     tmp = NP.array( thisline.split() )

#     # check the number of columns
#     nc = 9 # minimum number of required columns
#     if len( tmp ) < nc :
#         print >> sys.stderr, ( 'Line %i of input. Expecting >= %i cols, got %i' %
#                                ( lcount, nc, len(tmp) ) )
#         raise IOError(69,'Line %i of input. Expecting >= %i cols, got %i' %
#                       ( lcount, nc, len(tmp) ) )
#         return

#     # get coordinates, assume they represent the centroid
#     centroid = NP.array( [ float(tmp[0]), float(tmp[1]), float(tmp[2]) ] )

#     # get mechanism values
#     mrr, mtt, mff, mrt, mrf, mtf = float(tmp[3]), float(tmp[4]), float(tmp[5]), \
#         float(tmp[6]), float(tmp[7]), float(tmp[8])
#     expon = float(tmp[9]) # exponent

#     # get whatever is left as a string. reinsert the spaces 
#     endstr = "";
#     for i in range(10,len(tmp)): endstr += ' ' + tmp[i] 
#     endstr = endstr.rstrip('\n')
    
#     # get tensor norm without exponent
#     norm = sqrt( mrr**2 + mtt**2 + mff**2 + 2*( mrt**2 + mrf**2 + mtf**2 ) )

#     # Get the normalised tensor
#     m = NP.array( [mrr, mtt, mff, mtf, mrf, mrt] )/norm 
    
#     # Make a moment tensor object
#     MT = SymMT( m, norm*(10**expon), centroid )

#     return MT, endstr
def psmecaSmstring2SymMT( string ):
    
    # convert first 12 columns into array
    mtline = NP.fromstring( string, count=12, sep =' ', dtype=float )
    
    # get the location part
    c = mtline[0:3] # assume lon/lat/depth are centroid
    h = NP.array( [mtline[10], mtline[11], mtline[2]] ) # assume second lon/lat are hypocenter

    # get moment tensor part 
    mhat = mtline[3:9]
    exp = mtline[9]

    # remove the norm, take care with off diagonals, but save it
    norm = NP.sqrt( NP.sum( mhat**2 ) + NP.sum( mhat[3:]**2 ) )
    mhat /= norm
    norm *= 10**exp
        
    # make the moment tensor object
    return SymMT( mhat=mhat, Norm=norm, c=c, h=h )
    


#--------------------------------------------------
def readPsmecaList( istream ):
    """
    From an input file or stdin, read a list of moment tensors in psmeca form

    Expected format
    lon/lat/z/mrr/mtt/mpp/mrt/mrp/mtp/exp/lon0/lat0/str/anything else
    """

    mtlist=[] # this will be the output list

    # read everything
    alltxt = NP.genfromtxt( istream, delimiter='\n' , dtype=str)
    try: 
        istream.close()

    # loop through all tensors
    n = len(alltxt)
    for i in range(0,n):
        mtlist.append( psmecaSmstring2SymMT( alltxt[i] ) )

    
    return mtlist, alltxt

#--------------------------------------------------
def read_sdr( istream ):
    """
    From an input file or stdin, read a list of focal mechanisms in
    strike/dip/rake form

    lon/lat/depth/strike/dip/rake/mag
    """

    lcount = 0
    mtlist = []

    # get data in numpy array
    try:
        data = NP.loadtxt( istream, usecols=(0, 1, 2, 3, 4, 5, 6) ) #, autostrip=True )
    except IndexError:
        print >> sys.stderr, "Error: Require 7 columns: lon/lat/depth/strike/dip/rake/mag"
        return

    # Force 2d array and count the number of lines
    data = NP.array( data, ndmin = 2 )
    nline = NP.size( data, 0 )

    # loop through each line
    for i in range(0,nline):

        # get variables
        hypo = data[i, 0:3 ] # lon , lat, z
        sdr = data[i, 3:6]*pi/180 # in radians
        mag = data[i, 6]

        #TODO: error checking 

        # Make a double couple object
        MT = DoubleCouple( c=hypo, h=hypo,
                           strike = sdr[0], dip = sdr[1], rake=sdr[2],
                           mag=mag )

        mtlist.append( MT )
        
    return mtlist
    
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
