#!/usr/bin/env python

# genRandMT.py --- 
# 
# Filename: genRandMT.py
# Description: 
#   Generate a population of uniformly distributed random moment tensors
#   output in psmeca format
#
# Author: IW Bailey
# Maintainer: IW Bailey
# Created: Tue Apr 26 11:29:44 2011 (-0700)
# Version: 1
# Last-Updated: Mon Dec 19 11:17:47 2011 (-0800)
#           By: Iain Bailey
#     Update #: 118
# 
# 

# Change Log:
# 
# 
# 
# 

import sys
import numpy as NP
import argparse as AP

from pythmt.randomFuncs import randLonLat, randSphereSurf, randPBT
from pythmt import SymMT

#######################################################################    


### input arguments
parser = AP.ArgumentParser( "Generate uniform random distribution of moment tensors ")

# define region for centroid locations
parser.add_argument("-x", "--xmin", dest="x0",type = float, default=0.0 ,
                  help = "Minimum x  [default = 0.0]." )
parser.add_argument("-X", "--xmax", dest="x1",type = float, default=360.0 ,
                  help = "Maximum x  [default = 360.0]." )
parser.add_argument("-y", "--ymin", dest="y0",type = float, default=-90.0 ,
                  help = "Minimum y  [default = -90.0]." )
parser.add_argument("-Y", "--ymax", dest="y1",type = float, default=90.0 ,
                  help = "Maximum y  [default = 90.0]." )
parser.add_argument("-z", "--zmin", dest="z0",type = float, default=0.0 ,
                  help = "Minimum depth (km) [default = 0.0]." )
parser.add_argument("-Z", "--zmax", dest="z1",type = float, default=700.0 ,
                  help = "Maximum depth (km) [default = 1000.0]." )

# number of moment tensors
parser.add_argument( "-n", "--number", dest="nRand", type = int, default=100 ,
                     help = "Number of random tensors [100].")

# b-value
parser.add_argument("-b", "--bval", dest="b",type = float, default=0.0 ,
                    help = "b-value in G-R distribution [0.0]." )

# flag for only deviatoric mechanisms
parser.add_argument("--dev", action="store_true", dest="isDev", default=False,
                    help="Only deviatoric mechanisms")

# flag for only DC mechanisms
parser.add_argument("--dc", action="store_true", dest="isDC", default=False,
                    help="Only deviatoric, double-couple mechanisms")

args = parser.parse_args()

# nRand = number to generate 
# x0/x1/y0/y1/z0/z1: min max of region in lon/lat/depth
# b = b value
# isDev, isDC : flag for making deviatoric and or double couple


# set the a seed 
NP.random.seed(12345)

### get the coordinates 
(lon,lat) = randLonLat(args.nRand, args.x0, args.x1, args.y0, args.y1)
z = args.z0 + (args.z1-args.z0) * NP.random.rand( args.nRand, 1 )

### get the moments
# TODO
m0 = NP.ones( args.nRand )

### get the orientations
if ( args.isDC ):
    ( p, b, t) = randPBT(args.nRand)

    print NP.sum( p*b, 1 ), NP.sum( t*b, 1 ), NP.sum( p*t, 1 )
    print "TODO: still need to add DC option" 
    sys.exit()
    # TODO quaternion to mt
else:
    mtvecs = randSphereSurf(args.nRand,6)
    mtvecs[:,3:] /= NP.sqrt(2.0)
    if ( args.isDev ):
        # remove isotropic and re-normalise
        for i in range( 0,args.nRand ):
            mtvecs[i,0:3] -= NP.sum( mtvecs[i,0:3] )/3

            mtvecs[i,:] /= NP.sqrt( NP.sum( mtvecs[i,:]**2 ) + NP.sum( mtvecs[i,3:]**2 ) )



### print the outputs
for i in range(0,args.nRand):

    # make moment tensor
    mt = SymMT( mtvecs[i,:], m0[i], [lon[i,0], lat[i,0], z[i,0]] )
            
    # print output 
    print '%-10.4f %-10.4f %-8.2f %9.6f %9.6f %9.6f %9.6f %9.6f %9.6f %3i ' \
        % mt.getPsmecaSm()


#######################################################################    
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
######################################################################    
