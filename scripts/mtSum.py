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
# Last-Updated: Sun Sep  4 16:46:42 2011 (-0700)
#           By: Iain William Bailey
#     Update #: 98
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

# generate the empty moment tensor summation
mtSum = SMT.MTsum()

# normalise correction
if( args.stype == 1):  
    for i in range( 0, ndata ): mtlist[i].Norm = 1

# add the tensors
for i in range( 0, ndata ):
    mtSum.add( mtlist[i] )

# output the result

# get psmeca format for 1st 10 columns
(X, Y, depth, mrr, mtt, mff, mrt, mrf, mtf, exp ) = mtSum.getPsmecaSm()
args.ofile.write('%-10.6f %10.6f %6.2f %9.6f %9.6f %9.6f %9.6f %9.6f %9.6f %3i ' % 
                 ( X, Y, depth, mrr, mtt, mff, mrt, mrf, mtf, exp ) )


# print the number of events
args.ofile.write('%4i\n' % mtSum.count)


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

