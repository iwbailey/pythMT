#!/usr/bin/env python
# psmecaSmtoEig.py.py --- 
# 
# Filename: psmecaSmtoSx.py.py
# Description: 
#
#  Read in moment tensors in a psmeca input Sm (moment tensor) format,
#  convert to Sx (eigen value eigen vector format)
#
# Author: IW Bailey
# Maintainer: IW Bailey
# Created: Fri Mar 11 15:31:41 2011 (-0800)
# Version: 1
# Last-Updated: Fri Aug  5 16:11:13 2011 (-0700)
#           By: Iain Bailey
#     Update #: 57
#  
# Change Log:
# 
# 
# 
# 

# Standard libraries used
import sys
from optparse import OptionParser

# Personal libraries used
from ioFunctions import readpsmecaSm  

# get command line options
parser = OptionParser()
parser.add_option("--Sx",action="store_true", dest="sx", default=False,
                  help="Output psmeca Sx format")
(opt, args)=parser.parse_args()

if (not opt.sx): 
    opt.sx = True


# read in from stdin assuming psmeca format
Ndata = 0
while 1:
    thisline = sys.stdin.readline()

    if thisline != '':

        # Read line as SymMT object
        MT, extra = readpsmecaSm( thisline )
        Ndata += 1

        # output the original line with faulting style tacked on the end
        if opt.sx:
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
