#!/usr/bin/env python
# getMTeig.py --- 
# 
# Filename: getMTeig.py
# Description: 
#
#  Read in moment tensors in a psmeca input format, extract user
#  specified quantities related to the eigen values/ vectors for each
#
# Author: IW Bailey
# Maintainer: IW Bailey
# Created: Thu May  5 12:03:02 2011 (-0700)
# Version: 1
# Last-Updated: Mon Dec 19 11:22:59 2011 (-0800)
#           By: Iain Bailey
#     Update #: 37
#  
# Change Log:
# 
# 
# 
# 

# Standard libraries used
import sys
import argparse as AP

from optparse import OptionParser
from math import log10, sqrt
import numpy as NP

# Personal libraries used
from ioFunctions import readpsmecaSm  

# constants 
sqrt2 = sqrt(2.0)

# get command line options
parser = AP.ArgumentParser( "Generate uniform random distribution of moment tensors ")

# TODO
parser.add_option("--pos",action="store_true", dest="isPos", default=False,
                  help="Get position")
parser.add_option("-p",action="store_true", dest="isP", default=False,
                  help="Get Output only the p axis up/S/E order")
parser.add_option("-t",action="store_true", dest="isT", default=False,
                  help="Get Output only the t axis up/S/E order")
parser.add_option("-b",action="store_true", dest="isB", default=False,
                  help="Get Output only the b axis up/S/E order")
parser.add_option("--m0",action="store_true", dest="isM0", default=False,
                  help="Get scalar moment")
parser.add_option("--mw",action="store_true", dest="isMw", default=False,
                  help="Get moment magnitude")
parser.add_option("--check",action="store_true", dest="isCheck", default=False,
                  help="Random used for debugging")
(opt, args)=parser.parse_args()


# read in from stdin assuming psmeca format
Ndata = 0
while 1:
    thisline = sys.stdin.readline()

    if thisline != '':

        # Read line as SymMT object
        MT, extra = readpsmecaSm( thisline )
        Ndata += 1

        # get vals and vecs in E-N-Up coord system
        (tval, t, bval, b, pval, p) = MT.getPTB()

        # debugging option
        if opt.isCheck:
            # if( sqrt2*MT.mhat[2] > 1 ): 
            #     print NP.sum( MT.mhat[0:3]*MT.mhat[0:3] ) + 2*NP.sum( MT.mhat[3:6]*MT.mhat[3:6] )
            print MT.mhat[0:3]
#                print thisline.rstrip('\n') 
            continue

        # optional output the position
        if opt.isPos:
            sys.stdout.write("%-8.4f %8.4f %6.2f " %  (MT.c[0],MT.c[1],MT.c[2]) )

        # Output the t-axis
        if opt.isT or ( not opt.isP and not opt.isT and not opt.isB ):
            tr, tt, tf = t[2], -t[1], t[0] 
            sys.stdout.write("%-10.6g %10.6f %10.6f %10.6f " %  (tval, tr, tt, tf) )
        # Output the b-axis
        if opt.isB or ( not opt.isP and not opt.isT and not opt.isB ):
            br, bt, bf = b[2], -b[1], b[0] 
            sys.stdout.write("%-10.6g %10.6f %10.6f %10.6f " %  (bval, br, bt, bf) )
        # Output the p-axis
        if opt.isP or ( not opt.isP and not opt.isT and not opt.isB ):
            pr, pt, pf = p[2], -p[1], p[0] 
            sys.stdout.write("%-10.6g %10.6f %10.6f %10.6f " %  (pval, pr, pt, pf) )
            
        if opt.isM0:
            sys.stdout.write("%-8.6e " %  MT.getM0() )
        if opt.isMw:
            sys.stdout.write("%-6.3f " %  MT.getMw() )
        sys.stdout.write("\n")

    else: break
            
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
# getMTeig.py ends here
