#!/usr/bin/env python
# mtSum.py --- 
# 
# Filename: mtSum.py
# Description: 
#
# Read in moment tensors, various formats, output the summed tensor
#
# type 'mtSum.py -h' for help
#
# Author: IW Bailey
# Maintainer: IW Bailey
# Created: Fri Mar 11 15:31:41 2011 (-0800)
# Version: 1
# Last-Updated: Fri Aug  5 16:11:00 2011 (-0700)
#           By: Iain Bailey
#     Update #: 92
# 
# Change Log:
#
# Fri May  6 2011: Fixed a bug that was making things at z=z0 go into the last bin
# 
# 

# std python libraries used 
import numpy as NP
import sys
import argparse as AP

# personal libraries used, these need to be in the same directory or python path
import ioFunctions as IO
import SummedMomentTensor as SMT

# define constants
progname='mtSum'

#------------------------------
# command line arguments
parser = AP.ArgumentParser( description='Sum a set of moment tensors' )

parser.add_argument('ifile', nargs='?', type=AP.FileType('r'), default=sys.stdin, 
                    help = "Input file [stdin]")

parser.add_argument('ofile', nargs='?', type=AP.FileType('w'), default=sys.stdout,
                    help = "Output file [stdout]")

parser.add_argument("--ifmt", dest ="ifmt", type=int, default=0,
                    help = "Input data format; 0=psmeca/1=strike,dip,rake [0]" )

parser.add_argument("--stype", dest="stype", type=int, default=0,
                    help = "Summation type; 0=Kostrov, 1=Sum Normalised [0]" )

parser.add_argument("-v" , "--verbose", action="store_true", dest="isVb", default=False,
                    help="Verbose output")

args = parser.parse_args()

#------------------------------

# read the moment tensors
if( args.ifmt == 1 ):
    if( args.isVb ): sys.stderr.write('Reading Strike/dip/rake.\n' )
    mtlist = IO.readSDRList( args.ifile )
else:
    mtlist = IO.readPsmecaList( args.ifile )



ndata = len(mtlist)
if( args.isVb ): sys.stderr.write('%i events read in.\n' % ndata )

# generate the empty moment tensor summation
mtSum = SMT.MTsum()

# normalise correction
if( args.stype == 1):  
    for i in range( 0, ndata ):
        mtlist[i].Norm = 1

# add the tensors
for i in range( 0, ndata ):
    mtSum.add( mtlist[i] )

# output the result

# location            
for i in range(0,3): args.ofile.write('%8.2f ' % (mtSum.getCentroid(i)))

# get psmeca format 
(X, Y, depth, mrr, mtt, mff, mrt, mrf, mtf, exp ) = mtSum.getPsmecaSm()

args.ofile.write('%9.6f %9.6f %9.6f %9.6f %9.6f %9.6f %3i ' % 
                 ( mrr, mtt, mff, mrt, mrf, mtf, exp ) )

# print the number of events
args.ofile.write('%4i\n' % mtSum.count)

# # write bin info to stdout
# for k in range(0,nz):
#     z = opt.z0 + k*dz + 0.5*dz # bin coords

#     for j in range(0,ny):
#         y = opt.y0 + j*dy + 0.5*dy
#         for i in range(0, nx):
#             x = opt.x0 + i*dx + 0.5*dx

#             if( bins[k*nx*ny +j*nx + i].count > 0 ):

#                 # get the centroid locations
#                 c = bins[k*nx*ny +j*nx + i].getCentroid();
#                 Cx, Cy, Cz = c[0], c[1], c[2]
 
#                 # print centroid location
#                 if opt.otype == 4:
#                     # cross section
#                     sys.stdout.write('%8.2f %8.2f %8.2f ' % (Cx,Cz,Cy))
#                 else:
#                     # map view
#                     sys.stdout.write('%8.2f %8.2f %8.2f ' % (Cx, Cy, Cz))

#                 # print the tensor
#                 if opt.otype == 1: # mxx, myy...
#                     for ti in range(0,6):
#                         sys.stdout.write('%10.6f ' % bins[k*nx*ny +j*nx + i].smt(ti) )
#                     sys.stdout.write('%10.6e ' % bins[k*nx*ny +j*nx + i].getM0() )

#                 elif opt.otype == 2: # T tx ty ...
#                     vecs,vals = bins[k*nx*ny +j*nx + i].getEig()
#                     for ti in range(0,3):
#                         sys.stdout.write('%14.6e ' % vals[2-ti] )
#                         sys.stdout.write('%10.6f %10.6f %10.6f ' % \
#                                         (vecs[0,2-ti],vecs[1,2-ti],vecs[2,2-ti]) )
#                     sys.stdout.write('%10.6e ' % bins[k*nx*ny +j*nx + i].getM0() )

#                 elif opt.otype == 3:# T taz tpl ...
#                     vecs,vals = bins[k*nx*ny +j*nx + i].getEig()
#                     exp = int( log10(0.5*(abs(vals[0])+abs(vals[2]))) + 0.5 )
#                     fct = 10**exp # factor to remove
#                     for ti in range(0,3):
#                         sys.stdout.write('%9.5f ' % (vals[2-ti]/fct) )
#                         if( vecs[2,2-ti] > 0 ): 
#                             vecs[:,2-ti] *= -1
#                         az = R2D*atan2(vecs[0,2-ti],vecs[1,2-ti])
#                         pl = -R2D*asin(vecs[2,2-ti])
#                         sys.stdout.write('%6.2f %6.2f' % (az,pl) )
#                     sys.stdout.write('%3i ' % exp )
 
#                 elif opt.otype == 4: # cross section type psmeca
#                     sys.stdout.write('%9.6f %9.6f %9.6f %9.6f %9.6f %9.6f %3i ' % \
#                                          ( bins[k*nx*ny +j*nx + i].getPsmecaSide() ) )

#                 else: # default is option 0, psmeca
#                     (X, Y, depth, mrr, mtt, mff, mrt, mrf, mtf, exp ) = \
#                         bins[k*nx*ny +j*nx + i].getPsmecaSm()
#                     sys.stdout.write('%9.6f %9.6f %9.6f %9.6f %9.6f %9.6f %3i ' % \
#                                          ( mrr, mtt, mff, mrt, mrf, mtf, exp ) )
                   
#                  # print bin coordinates
#                 if opt.otype == 4:
#                     # cross section
#                     sys.stdout.write('%8.2f %8.2f %8.2f ' % (x,z,y))
#                 else:
#                     # map view
#                     sys.stdout.write('%8.2f %8.2f %8.2f ' % (x,y,z))

#                 # print the number of events
#                 sys.stdout.write('%4i\n' % bins[k*nx*ny +j*nx + i].count)

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
# mtSum.py ends here

