#!/usr/bin/env python
# getInfo.py --- 
# 
# Filename: getInfo.py
# Description: 
#
#  Read in moment tensors in a psmeca input format, extract user
#  specified quantities for each
#
# Author: IW Bailey
# Maintainer: IW Bailey
# Created: Fri Mar 11 15:31:41 2011 (-0800)
# Version: 1
# Last-Updated: Fri May  6 14:07:53 2011 (-0700)
#           By: Iain Bailey
#     Update #: 157
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
from MomentTensor import readpsmecaSm  

# constants 
sqrt2 = sqrt(2.0)
invs2= 1/sqrt2

# get command line options
parser = OptionParser()
parser.add_option("--all",action="store_true", dest="all", default=False,
                  help="Print all with other fields added on end")
parser.add_option("--pos",action="store_true", dest="pos", default=False,
                  help="Get position")
parser.add_option("--smt",action="store_true", dest="smt", default=False,
                  help="Get source mechanism tensor in rr/tt/ff/rt/rf/tf order")
parser.add_option("--m0",action="store_true", dest="m0", default=False,
                  help="Get scalar moment")
parser.add_option("--mw",action="store_true", dest="mw", default=False,
                  help="Get moment magnitude")
parser.add_option("--lognm",action="store_true", dest="lognm", default=False,
                  help="Get log_10 of the norm of the tensor")
parser.add_option("--fclvd",action="store_true", dest="fclvd", default=False,
                  help="Get f_clvd")
parser.add_option("--Gamma",action="store_true", dest="gamma", default=False,
                  help="Get Gamma (CLVD measure)")
parser.add_option("--mrr",action="store_true", dest="mrr", default=False,
                  help="Get M_rr/||M||")
parser.add_option("--check",action="store_true", dest="check", default=False,
                  help="Random used for debugging")
(opt, args)=parser.parse_args()

if (not opt.all and 
    not opt.pos and 
    not opt.smt and
    not opt.m0 and 
    not opt.mw and 
    not opt.lognm and 
    not opt.fclvd and
    not opt.mrr and 
    not opt.gamma and
    not opt.check ): 
    sys.stderr.write( "Nothing reqested.  Exiting...\n" )
    sys.exit()

# read in from stdin assuming psmeca format
Ndata = 0
while 1:
    thisline = sys.stdin.readline()

    if thisline != '':

        # Read line as SymMT object
        MT, extra = readpsmecaSm( thisline )
        Ndata += 1

        # output the original line with faulting style tacked on the end
        if opt.all:
            sys.stdout.write("%-s " % thisline.rstrip('\n') )
        elif opt.pos:
            sys.stdout.write("%-8.4f %8.4f %6.2f " %  (MT.c[0],MT.c[1],MT.c[2]) )

        if opt.smt:
            # get tensor components, r = up, t = south, f = east
            mrr, mtt, mff = MT.mhat[2], MT.mhat[1], MT.mhat[0]
            mrt, mrf, mtf = -MT.mhat[3], MT.mhat[4], -MT.mhat[5]
            sys.stdout.write("%-8.6f %8.6f %8.6f %8.6f %8.6f %8.6f" %  
                             ( mrr, mtt, mff, mrt, mrf, mtf ))
        if opt.m0:
            sys.stdout.write("%-8.6e " %  MT.getM0() )
        if opt.mw:
            sys.stdout.write("%-6.3f " %  MT.getMw() )
        if opt.lognm:
            sys.stdout.write("%-5.3f " %  log10(MT.Norm) )
        if opt.fclvd:
            sys.stdout.write("%-5.3f " %  MT.getFclvd() )
        if opt.gamma:
            sys.stdout.write("%-5.3f " %  MT.getGamma() )
        if( opt.mrr ):
            sys.stdout.write("%5.3f " %  ( sqrt2*MT.smt(3,3) ) )
        if opt.check:
            # if( sqrt2*MT.mhat[2] > 1 ): 
            #     print NP.sum( MT.mhat[0:3]*MT.mhat[0:3] ) + 2*NP.sum( MT.mhat[3:6]*MT.mhat[3:6] )
            print MT.mhat[0:3]
#                print thisline.rstrip('\n') 
            continue
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
# getInfo.py ends here
