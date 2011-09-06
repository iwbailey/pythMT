#!/usr/bin/env python
# mtDecomp.py --- 
# 
# Filename: mtDecomp.py
# Description: 
#
#  Read in moment tensors in a psmeca input format, decompose into DC,
#  iso, etc. according to options.
#  Output in the same format
#
# Author: IW Bailey
# Maintainer: IW Bailey
# Created: Fri Jun  3 16:23:39 2011 (-0700)
# Version: 1
# Last-Updated: Sun Sep  4 18:22:44 2011 (-0700)
#           By: Iain William Bailey
#     Update #: 97
#  
# Change Log:
# 
# 
# 
# 

# Standard libraries used
import sys
from optparse import OptionParser
from math import log10, sqrt
import numpy as NP

# Personal libraries used
from ioFunctions import readpsmecaSm  

# constants 
sqrt2 = sqrt(2.0)

# get command line options
parser = OptionParser()
parser.add_option("--dev",action="store_true", dest="isDev", default=False,
                  help="Get deviatoric moment tensor")
parser.add_option("--maxdc",action="store_true", dest="isMaxDC", default=False,
                  help="Get  max DC, by taking away min CLVD")
parser.add_option("--mindc",action="store_true", dest="isMaxDC", default=False,
                  help="Get minimum DC, by taking away max DC")
parser.add_option("--maxclvd",action="store_true", dest="isMaxCLVD", default=False,
                  help="Get solution with maximum CLVD")
parser.add_option("--check",action="store_true", dest="isCheck", default=False,
                  help="Random used for debugging")
(opt, args)=parser.parse_args()


# read in from stdin assuming psmeca format
Ndata = 0
lcount = 0

while 1:
    thisline = sys.stdin.readline()

    lcount += 1
    if thisline != '':
        
        # Read line as SymMT object
        try:
            MT, extra = readpsmecaSm( thisline, lcount )
            Ndata += 1
        except:
            sys.exit()
            continue


        # debugging option
        if opt.isCheck:
            # check decomposition is working 
            (MTdc, MTclvd, MTiso) = MT.decompose( 2 )
            
            print "ALL:", MT.getMvec()
            print "sum:", MTdc.getMvec() + MTclvd.getMvec() + MTiso.getMvec()
            print "DC:",MTdc.getMvec(), MTdc.getFclvd()
            print "CLVD:",MTclvd.getMvec(), MTclvd.getFclvd()
            print "ISO:", MTiso.getMvec()

            sys.exit()

        # do the decomposition
        if opt.isMaxDC:
            # take the 1/2 the size of the b axis away from the P and T-axes
            (MTdc, MTclvd, MTiso) = MT.decompose( 2 )
            (X, Y, depth, mrr, mtt, mff, mrt, mrf, mtf, exp ) = MTdc.getPsmecaSm()
            sys.stdout.write('%10.6f %10.6f %10.6f %9.6f %9.6f %9.6f %9.6f %9.6f %9.6f %3i\n' % 
                             ( X, Y, depth, mrr, mtt, mff, mrt, mrf, mtf, exp ) )
        elif opt.isMinDC:
           # take the 2* the size of the b axis away from the P and T-axes
            (MTdc, MTclvd, MTiso) = MT.decompose( 1 )
            (X, Y, depth, mrr, mtt, mff, mrt, mrf, mtf, exp ) = MTdc.getPsmecaSm()
            sys.stdout.write('%10.6f %10.6f %10.6f %9.6f %9.6f %9.6f %9.6f %9.6f %9.6f %3i\n' % 
                             ( X, Y, depth, mrr, mtt, mff, mrt, mrf, mtf, exp ) )
        elif opt.isMaxCLVD:
           # take the 2* the size of the b axis away from the P and T-axes
            (MTdc, MTclvd, MTiso) = MT.decompose( 1 )
            (X, Y, depth, mrr, mtt, mff, mrt, mrf, mtf, exp ) = MTclvd.getPsmecaSm()
            sys.stdout.write('%10.6f %10.6f %10.6f %9.6f %9.6f %9.6f %9.6f %9.6f %9.6f %3i\n' % 
                             ( X, Y, depth, mrr, mtt, mff, mrt, mrf, mtf, exp ) )

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
# mtDecomp.py ends here
