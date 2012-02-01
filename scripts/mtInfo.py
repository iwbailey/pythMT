#!/usr/bin/env python
# mtInfo.py --- 
# 
# Filename: mtInfo.py
# Description: 
#
#  Read in moment tensors in a psmeca input format, extract user
#  specified quantities for each
#
# Author: IW Bailey
# Maintainer: IW Bailey
# Created: Fri Mar 11 15:31:41 2011 (-0800)
# Version: 2
# Last-Updated: Wed Feb  1 14:24:43 2012 (-0800)
#           By: Iain Bailey
#     Update #: 286
#  
# Change Log:
# Thu Dec 29 15:38:54 2011 : Changed name of file, optparse -> argparse
# Sun Sep  4 17:14:40 2011 (-0700): Fixed for r-theta-phi changes
# 
# 

# Standard libraries used
import sys
import argparse
from math import log10, sqrt
import numpy as NP

# Personal libraries used
from pythmt.iofuncs import read_psmecalist

# constants 
sqrt2 = sqrt(2.0)
invs2= 1/sqrt2

#--------------------------------------------------
# get command line options
parser = argparse.ArgumentParser( 
    description='Extract quantities from a psmeca style input file' )

parser.add_argument('ifile', nargs='?', type=argparse.FileType('r'), 
                    default=sys.stdin, 
                    help = "Input file, psmeca -Sm format [stdin]")

parser.add_argument('ofile', nargs='?', type=argparse.FileType('w'), 
                    default=sys.stdout, 
                    help = "Output file. [stdout]")

parser.add_argument("--paste",action="store_true", dest="isPaste", default=False,
                  help="Paste the new values on the end")
parser.add_argument("--pos",action="store_true", dest="pos", default=False,
                  help="Get position")
parser.add_argument("--smt",action="store_true", dest="smt", default=False,
                  help="Get source mechanism tensor in rr/tt/ff/rt/rf/tf order")
parser.add_argument("--m0",action="store_true", dest="m0", default=False,
                  help="Get scalar moment")
parser.add_argument("--mw",action="store_true", dest="mw", default=False,
                  help="Get moment magnitude")
parser.add_argument("--lognm",action="store_true", dest="lognm", default=False,
                  help="Get log_10 of the norm of the tensor")
parser.add_argument("--fclvd",action="store_true", dest="fclvd", default=False,
                  help="Get f_clvd")
parser.add_argument("--Gamma",action="store_true", dest="gamma", default=False,
                  help="Get Gamma (CLVD measure)")
parser.add_argument("--frr",action="store_true", dest="frr", default=False,
                  help="Get M_rr/max(P,T)")

parser.add_argument("-v" , "--verbose", action="store_true", dest="isVb", default=False,
                    help="Verbose output")

args=parser.parse_args()

#--------------------------------------------------

# check some info was provided
if args.ifile.isatty():
    print >> sys.stderr, "Error: No ifile and nothing in stdin. Exiting..."
    sys.exit()

# read the input into moment tensors
(mtlist, alltxt) = read_psmecalist( args.ifile )

# loop through all
n = len(mtlist)
if( args.isVb ): print "Read %i moment tensors" % n

sth=False # flag so we write a new line

for i in range(0,n):
    
    # output the original line 
    if args.isPaste:
        sth = True
        args.ofile.write( "%s " % alltxt[i] )
    if args.pos: 
        # output the position
        sth = True
        args.ofile.write("%-8.4f %8.4f %6.2f " %  (mtlist[i].c[0],mtlist[i].c[1],mtlist[i].c[2]) )

    if args.smt:
        # output the source mech tensor
        sth = True
        # get tensor components, r = up, t = south, f = east
        mrr, mtt, mff = mtlist[i].mhat[0], mtlist[i].mhat[1], mtlist[i].mhat[2]
        mrt, mrf, mtf = mtlist[i].mhat[5], mtlist[i].mhat[4], mtlist[i].mhat[3]
        args.ofile.write("%-9.6f %9.6f %9.6f %9.6f %9.6f %9.6f " %  
                             ( mrr, mtt, mff, mrt, mrf, mtf ))
    if args.m0:
        sth = True
        args.ofile.write("%-8.6e " %  mtlist[i].M0() )
    if args.mw:
        sth = True
        args.ofile.write("%6.3f " %  mtlist[i].Mw() )
    if args.lognm:
        sth = True
        args.ofile.write("%6.3f " %  log10(mtlist[i].Norm) )
    if args.fclvd:
        sth = True
        args.ofile.write("%6.3f " %  mtlist[i].f_clvd() )
    if args.gamma:
        sth = True
        args.ofile.write("%6.3f " %  mtlist[i].Gamma() )
    if( args.frr ):
        sth = True
        (tval, t, bval, b, pval, p) = mtlist[i].pbt()
        mrr =  mtlist[i].Norm*mtlist[i].smt(0,0)
        args.ofile.write("%6.3f " %  ( mrr/NP.max([tval,-pval]) ) )

    if sth: args.ofile.write("\n")
    else: 
        print "Nothing requested"
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
# mtInfo.py ends here
