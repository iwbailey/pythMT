#!/usr/bin/env python
# mtBinSum.py --- 
# 
# Filename: mtBinSum.py
# Description: 
#
# Read in moment tensors, various formats, output the summed tensor
#
# type 'mtBinSum.py -h' for help
#
# Author: IW Bailey
# Maintainer: IW Bailey
# Created: Fri Mar 11 15:31:41 2011 (-0800)
# Version: 1
# Last-Updated: Mon Dec 19 10:57:03 2011 (-0800)
#           By: Iain Bailey
#     Update #: 208
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

# personal libraries used
import pythmt.ioFunctions as IO
import pythmt.EqkBin as EB
import pythmt.SummedMomentTensor as SMT

# define constants
progname='mtBinSum'

#------------------------------
# command line arguments
parser = AP.ArgumentParser( description='Sum a set of moment tensors' )

parser.add_argument('ifile', nargs='?', type=AP.FileType('r'), default=sys.stdin, 
                    help = "Input file [stdin]")

parser.add_argument('ofile', nargs='?', type=AP.FileType('w'), default=sys.stdout,
                    help = "Output file [stdout]")

parser.add_argument('--m01', nargs=2, type=float, default =[-99.9, 99.9],
                    help = "Magnitude range [-99 to 99]" )
parser.add_argument('--nm', dest = "nm", type=int, default=1,
                    help = "Number of magnitude bins [1]" )

parser.add_argument('--x01', nargs=2, type=float, default =[0.0, 360.0],
                    help = "Longitude or x range [0 to 360]" )
parser.add_argument('--nx', dest = "nx", type=int, default=1,
                    help = "Number of x bins [1]" )

parser.add_argument('--y01', nargs=2, type=float, default =[-90.0, 90.0],
                    help = "Latitude or y range [-90 to 90]" )
parser.add_argument('--ny', dest = "ny", type=int, default=1,
                    help = "Number of y bins [1]" )

parser.add_argument('--z01', nargs=2, type=float, default =[0.0, 700.0],
                    help = "Depth or z range [0.0 to 700]" )
parser.add_argument('--nz', dest = "nz", type=int, default=1,
                    help = "Number of z bins [1]" )

parser.add_argument("--ifmt", dest ="ifmt", type=int, default=0,
                    help = "Input data format; 0=psmeca / 1=strike,dip,rake [0]" )

parser.add_argument("--stype", dest="stype", type=int, default=0,
                    help = "Summation type; 0=Kostrov, 1=Sum Normalised [0]" )

parser.add_argument("-v" , "--verbose", action="store_true", dest="isVb", default=False,
                    help="Verbose output")

args = parser.parse_args()

#------------------------------

# read the moment tensors
if( args.ifmt == 1 ):
    if( args.isVb ): sys.stderr.write('Reading Strike/dip/rake.\n' )
    mtlist = IO.readSDRfile( args.ifile )
else:
    mtlist = IO.readPsmecaList( args.ifile )

ndata = len(mtlist)
if( args.isVb ): sys.stderr.write('%i events read in.\n' % ndata )

# normalise correction
if( args.stype == 1):  
    for i in range( 0, ndata ): mtlist[i].Norm = 1

# bin in magnitude
binx = 0.5*(args.x01[0]+args.x01[1])
biny = 0.5*(args.y01[0]+args.y01[1])
binz = 0.5*(args.z01[0]+args.z01[1])
( binlist1, ninside) = EB.binMT_mag( mtlist, args.m01[0], args.m01[1], args.nm, 
                                    NP.array([binx, biny, binz]) )

if( args.isVb ): sys.stderr.write('%i events within Mag range.\n' % ninside )

# check if also spatial binning
nxyz = args.nx * args.ny *args.nz 

if nxyz == 0:
    binlist = binlist1
else:
    # bin in xyz
    binlist =  []
    nmbins = len(binlist1) # actual number of magnitude bins containing data
    ninside = 0

    # loop through mag bins and further subdivide
    for i in range(0, nmbins):
        if( binlist1[i] == None ): continue # skip empty bins
        (tmp_binlist, tmp_nin) = EB.binMT_xyz( binlist1[i].eqklist, 
                                               args.x01[0], args.x01[1], args.nx, 
                                               args.y01[0], args.y01[1], args.ny, 
                                               args.z01[0], args.z01[1], args.nz)
        binlist.extend( tmp_binlist )
        ninside+=tmp_nin

if( args.isVb ): sys.stderr.write('%i events within XYZ range.\n' % ninside )

# compute summed tensors
nbins = len( binlist )
if( args.isVb ): sys.stderr.write('%i bins in total considered.\n' % nbins )

# loop through all bins
mtsumlist = []
for i in range(0,nbins):
    # skip if empty
    if( binlist[i] == None ):
        mtsumlist.append(None)
        continue

    # make summed tensor and add to list
    mtSum = SMT.MTsum()
    for mt in binlist[i].eqklist[:]: mtSum.add( mt )

    mtsumlist.append(mtSum)

# output
for i in range(0,nbins):

    # skip if no data in bin
    if( binlist[i] == None ): continue

    # get psmeca format for 1st 10 columns
    (X, Y, depth, mrr, mtt, mff, mrt, mrf, mtf, exp ) = mtsumlist[i].getPsmecaSm()
    args.ofile.write('%-10.6f %10.6f %6.2f %9.6f %9.6f %9.6f %9.6f %9.6f %9.6f %3i ' % 
                     ( X, Y, depth, mrr, mtt, mff, mrt, mrf, mtf, exp ) )

    # write bin info, keeping psmeca form
    args.ofile.write('%10.4f %10.4f Bin_%i %8.2f %6.2f %6i\n' % 
                     ( binlist[i].xyz[0], binlist[i].xyz[1], i, binlist[i].xyz[2], 
                       binlist[i].mag, binlist[i].n ) )


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
# mtBinSum.py ends here

