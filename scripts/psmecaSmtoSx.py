#!/usr/bin/env python
# psmecaSmtoSx.py --- 
# 
# Filename: psmecaSmtoSx.py
# Description: 
#
#  Read in moment tensors in a psmeca input Sm (moment tensor) format,
#  convert to Sx (eigen value eigen vector format)
#
# Author: IW Bailey
# Maintainer: IW Bailey
# Created: Fri Mar 11 15:31:41 2011 (-0800)
# Version: 1
# Last-Updated: Wed Feb  1 15:21:41 2012 (-0800)
#           By: Iain Bailey
#     Update #: 74
#  
# Change Log:
#
# Wed Feb  1 2012 - changed to updated functions, changed file name
# 
# 
# 

# Standard libraries used
import sys
import argparse

# Personal libraries used
from pythmt.iofuncs import psmeca2SymMT

# get command line options
parser = argparse.ArgumentParser(
    description='From stdin read psmeca -Sm format, output psmeca -Sx format')

(args)=parser.parse_args()

# read in from stdin assuming psmeca format
Ndata = 0
while 1:
    thisline = sys.stdin.readline()

    if thisline != '':
        MT = psmeca2SymMT( thisline )

        (x,y,z,t,taz,tp,b,baz,bp,p,paz,pp,exp) = MT.getPsmecaSx()
        sys.stdout.write("%8.4f %8.4f %6.2f %10.6f %8.3f %8.3f %10.6f %8.3f %8.3f %10.6f %8.3f %8.3f %3i\n" % 
                         (x,y,z,t,taz,tp,b,baz,bp,p,paz,pp,exp))
        

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
# psmecaSmtoSx.py.py ends here
