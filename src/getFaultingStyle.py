#!/usr/bin/env python
# getFaultingStyle.py --- 
# 
# Filename: getFaultingStyle.py
# Description: 
#
#  Read in moment tensors in a psmeca input format, compute the
#  faulting style (reverse, normal, strike-slip, oblique) based on
#  which of the eigen vectors is most vertical.  Output the original
#  data with two extra columns definigng the faulting style
#
# Author: IW Bailey
# Maintainer: IW Bailey
# Created: Fri Mar 11 15:31:41 2011 (-0800)
# Version: 1
# Last-Updated: Sat Sep  3 10:58:28 2011 (-0700)
#           By: Iain William Bailey
#     Update #: 54
#  
# Change Log:
# 
# 
# 
# 

# Standard libraries used
import sys
from numpy import sin, pi, argmax, abs
import argparse as AP

# personal libraries used, these need to be in the same directory or python path
import ioFunctions as IO


#------------------------------
# command line arguments
parser = AP.ArgumentParser( description='Get the faulting style for a population of focal mechanisms' )

parser.add_argument('ifile', nargs='?', type=AP.FileType('r'), default=sys.stdin, 
                    help = "Input file [stdin]")

parser.add_argument('ofile', nargs='?', type=AP.FileType('w'), default=sys.stdout,
                    help = "Output file [stdout]")

parser.add_argument("--ifmt", dest ="ifmt", type=int, default=0,
                    help = "Input data format; 0=psmeca/1=strike,dip,rake [0]" )

parser.add_argument("-p", "--pmin", dest="pmin", type=float, default=None ,
                    help = "Minimum plunge for an axis to be defined as vertical [take min of 3 axes]." )

parser.add_argument("-v" , "--verbose", action="store_true", dest="isVb", default=False,
                    help="Verbose output [False]")

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

# convert min plunge to min abs val of z component
if( args.pmin != None ):
    minZ = sin(abs(args.pmin)*pi/180)
else:
    # no oblique category if pmin is not defined
    minZ = 0.0


# loop through each tensor
for i in range(0,ndata):

    # get eigen solution
    # ptb rows are r; theta; phi, cols are p, b, t
    [pbt, vals] = mtlist[i].getEig()
  
    # find which axis is most vertical
    maxzi = argmax( abs( pbt[0,:] ) )

    if abs(pbt[2,maxzi]) < minZ:
        # none of the axes are vertical enough
        sys.stdout.write("1 %-12s " % ("oblique") )
    elif maxzi == 0:
        # p (min) axis is most vertical
        sys.stdout.write("1 %-12s " % ("normal") )
    elif maxzi == 1:
        # b (intermediate) axis is most vertical
        sys.stdout.write("2 %-12s " % ("strike-slip") )
    elif maxzi == 2:
        # t (max) axis is most vertical
        sys.stdout.write("3 %-12s " % ("reverse") )
    else:
        print >>sys.stderr(), "ERROR: don't understand eig vals"

    print "%9.6f %9.6f %9.6f" % (pbt[0,0],pbt[0,1],pbt[0,2])
sys.exit()    
            
###########################################################
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
###########################################################
# getFaultingStyle.py ends here
