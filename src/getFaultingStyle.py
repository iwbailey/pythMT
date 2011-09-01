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
# Last-Updated: Fri Aug  5 16:09:16 2011 (-0700)
#           By: Iain Bailey
#     Update #: 2
#  
# Change Log:
# 
# 
# 
# 

# Standard libraries used
import sys
from math import sin, pi
from optparse import OptionParser

# Personal libraries used
from ioFunctions import readpsmecaSm  

# Get the minimum plunge
# we regard an axis as vertical if plunge greater than this 
parser = OptionParser()
parser.add_option("-p", "--pmin", dest="pmin", type = "float", default=60.0 ,
                  help = "Minimum plunge for an axis to be defined as vertical [60.0^o]." )
(opt, args)=parser.parse_args()

# convert min plunge to min abs val of z component
minZ = sin(abs(opt.pmin)*pi/180)

# read in from stdin assuming psmeca format
Ndata = 0
while 1:
    thisline = sys.stdin.readline()

    if thisline != '':

        # Read line as SymMT object
        MT, extra = readpsmecaSm( thisline )
        Ndata += 1

        # get faulting style
        [ptb, vals] = MT.getEig()

        if( abs(ptb[2,0])<minZ and  abs(ptb[2,1])<minZ and abs(ptb[2,2])<minZ ):
            # oblique if no axis is within maxplunge of vertical
            ftype='O'
            ftypei=0
        elif( abs( ptb[2,0] ) > abs( ptb[2,1] ) and
              abs( ptb[2,0] ) > abs( ptb[2,2] ) ):
            # normal if the P axis is most vertical
            ftype='N'
            ftypei=1
        elif( abs( ptb[2,1] ) > abs( ptb[2,2] ) ):
            # strike slip if the B axis is most vertical 
            ftype = 'S'
            ftypei = 2
        else:
            # Reverse if the T axis is most vertical
            ftype = 'R'
            ftypei = 3

        # output the original line with faulting style tacked on the end
        print "%s %c %i" % ( thisline.rstrip('\n'), ftype, ftypei )

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
# getFaultingStyle.py ends here
