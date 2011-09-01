#!/usr/bin/env python
# binAndSum.py --- 
# 
# Filename: binAndSum.py
# Description: 
#
# Read in moment tensors in psmeca format from stdin.  Allocate to
# bins on a regular grid based on specified parameters.  Compute the
# summed tensor for each bin, option to normalise each tensor before
# summation.  Output to stdout in user-specified format.
#
# type 'binAndSum.py -h' for help
#
# Author: IW Bailey
# Maintainer: IW Bailey
# Created: Fri Mar 11 15:31:41 2011 (-0800)
# Version: 1
# Last-Updated: Fri Aug  5 16:09:04 2011 (-0700)
#           By: Iain Bailey
#     Update #: 32
# 
# Change Log:
#
# Fri May  6 2011: Fixed a bug that was making things at z=z0 go into the last bin
# 
# 

# std python libraries used 
import numpy as NP
from math import sqrt, pi, atan2, asin, log10
import sys
from optparse import OptionParser

# personal libraries used, these need to be in the same directory or python path
import genParser as ARG
from ioFunctions import readpsmecaSm  
import SummedMomentTensor as SMT

# define constants
R2D=180/pi # Rad to deg conversion
progname='binAndSum'

# get the user input
# x0, x1, y0, y1, z0, z1 are dimensions of region
# nx, ny, nz are number of bins in each dimension
# stype = 0 then compute summed moment tensor
# stype = 1 then compute  summed source mechanism tensor
# otype = output type
parser = ARG.binAndSum()
(opt, args)=parser.parse_args()

# grid dimensions 
nx, ny, nz = opt.nx, opt.ny, opt.nz

# Get the dimensions of a single bin
dx = abs(opt.x1 - opt.x0)/nx
dy = abs(opt.y1 - opt.y0)/ny
dz = abs(opt.z1 - opt.z0)/nz

if( opt.verbose ):
    sys.stderr.write('%s: \tlimits(x) = %6.2f, %6.2f\n' % (progname,opt.x0, opt.x1))
    sys.stderr.write('\tlimits(y) = %6.2f, %6.2f\n' % (opt.y0, opt.y1))
    sys.stderr.write('\tlimits(z) = %6.2f, %6.2f\n' % (opt.z0, opt.z1))
    sys.stderr.write('\tdx = %6.4f, dy = %6.4f, dz = %6.4f km\n' % (dx, dy, dz))

# generate the empty moment tensors in each bin
bins = [] 
for i in range(0,nx*ny*nz):
    bins.append( SMT.MTsum() )

# read in data and sum
Ndata , Nremoved = 0 , 0

# read in from stdin assuming psmeca format, allocate to bins
while 1:
    thisline = sys.stdin.readline()

    if thisline != '':

        # Read psmeca input into symmetric moment tensor object
        mt, extra = readpsmecaSm( thisline )
        Ndata += 1

        # check if in region
        if( mt.c[0] >= opt.x0 and mt.c[0] < opt.x1
            and mt.c[1] >= opt.y0 and mt.c[1] < opt.y1
            and mt.c[2] >= opt.z0 and mt.c[2] < opt.z1 ):

            # find bin indices, the max condition prevents indices
            # being -1 when exactly x0, y0, z0
            i = max( int( round ( nx*( mt.c[0] - opt.x0 ) / (opt.x1-opt.x0) - 0.5 ) ), 0)
            j = max( int( round ( ny*( mt.c[1] - opt.y0 ) / (opt.y1-opt.y0) - 0.5 ) ), 0)
            k = max( int( round ( nz*( mt.c[2] - opt.z0 ) / (opt.z1-opt.z0) - 0.5 ) ), 0)
            
            # Add to the existing summed tensor in that bin
            if( opt.stype == 1):  
                # normalized summation
                mt.Norm = 1
                bins[ k*nx*ny + j*nx + i ].add( mt )
            else: 
                # kostrov summation
                bins[ k*nx*ny + j*nx + i ].add( mt )
        else:
            Nremoved += 1

    else: break
            
if( opt.verbose ):
    sys.stderr.write('%6i events read in.\n' % Ndata )
    sys.stderr.write('%6i events were not in region.\n' % Nremoved )

# write bin info to stdout
for k in range(0,nz):
    z = opt.z0 + k*dz + 0.5*dz # bin coords

    for j in range(0,ny):
        y = opt.y0 + j*dy + 0.5*dy
        for i in range(0, nx):
            x = opt.x0 + i*dx + 0.5*dx

            if( bins[k*nx*ny +j*nx + i].count > 0 ):

                # get the centroid locations
                c = bins[k*nx*ny +j*nx + i].getCentroid();
                Cx, Cy, Cz = c[0], c[1], c[2]
 
                # print centroid location
                if opt.otype == 4:
                    # cross section
                    sys.stdout.write('%8.2f %8.2f %8.2f ' % (Cx,Cz,Cy))
                else:
                    # map view
                    sys.stdout.write('%8.2f %8.2f %8.2f ' % (Cx, Cy, Cz))

                # print the tensor
                if opt.otype == 1: # mxx, myy...
                    for ti in range(0,6):
                        sys.stdout.write('%10.6f ' % bins[k*nx*ny +j*nx + i].smt(ti) )
                    sys.stdout.write('%10.6e ' % bins[k*nx*ny +j*nx + i].getM0() )

                elif opt.otype == 2: # T tx ty ...
                    vecs,vals = bins[k*nx*ny +j*nx + i].getEig()
                    for ti in range(0,3):
                        sys.stdout.write('%14.6e ' % vals[2-ti] )
                        sys.stdout.write('%10.6f %10.6f %10.6f ' % \
                                        (vecs[0,2-ti],vecs[1,2-ti],vecs[2,2-ti]) )
                    sys.stdout.write('%10.6e ' % bins[k*nx*ny +j*nx + i].getM0() )

                elif opt.otype == 3:# T taz tpl ...
                    vecs,vals = bins[k*nx*ny +j*nx + i].getEig()
                    exp = int( log10(0.5*(abs(vals[0])+abs(vals[2]))) + 0.5 )
                    fct = 10**exp # factor to remove
                    for ti in range(0,3):
                        sys.stdout.write('%9.5f ' % (vals[2-ti]/fct) )
                        if( vecs[2,2-ti] > 0 ): 
                            vecs[:,2-ti] *= -1
                        az = R2D*atan2(vecs[0,2-ti],vecs[1,2-ti])
                        pl = -R2D*asin(vecs[2,2-ti])
                        sys.stdout.write('%6.2f %6.2f' % (az,pl) )
                    sys.stdout.write('%3i ' % exp )
 
                elif opt.otype == 4: # cross section type psmeca
                    sys.stdout.write('%9.6f %9.6f %9.6f %9.6f %9.6f %9.6f %3i ' % \
                                         ( bins[k*nx*ny +j*nx + i].getPsmecaSide() ) )

                else: # default is option 0, psmeca
                    (X, Y, depth, mrr, mtt, mff, mrt, mrf, mtf, exp ) = \
                        bins[k*nx*ny +j*nx + i].getPsmecaSm()
                    sys.stdout.write('%9.6f %9.6f %9.6f %9.6f %9.6f %9.6f %3i ' % \
                                         ( mrr, mtt, mff, mrt, mrf, mtf, exp ) )
                   
                 # print bin coordinates
                if opt.otype == 4:
                    # cross section
                    sys.stdout.write('%8.2f %8.2f %8.2f ' % (x,z,y))
                else:
                    # map view
                    sys.stdout.write('%8.2f %8.2f %8.2f ' % (x,y,z))

                # print the number of events
                sys.stdout.write('%4i\n' % bins[k*nx*ny +j*nx + i].count)

######################################################################    
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
# binAndSum.py ends here

