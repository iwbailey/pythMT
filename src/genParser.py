#!/usr/bin/env python
# genParser.py --- 
# 
# Filename: genParser.py
# Description: 
#
#  Define the terminal inputs for MT summation programs
#
# Author: IW Bailey
# Maintainer: IW Bailey
# Created: Fri Mar 11 15:31:41 2011 (-0800)
# Version: 1
# Last-Updated: Tue Apr 26 13:01:15 2011 (-0700)
#           By: Iain Bailey
#     Update #: 30
# 
# 

# Change Log:
# 
# Tue Apr 26 2011: added parser for genRandMT.py
# 
# 
# 
import sys
from optparse import OptionParser

#----------------------------------------------------------------------
def binAndSum():
    """
    Function generates the input parser for the binAndSum program
    """
    parser = OptionParser()

    # define region
    parser.add_option("-x", "--xmin", dest="x0",type = "float", default=0.0 ,
                      help = "Minimum x  [default = 0.0]." )
    parser.add_option("-X", "--xmax", dest="x1",type = "float", default=360.0 ,
                      help = "Maximum x  [default = 360.0]." )
    parser.add_option("-y", "--ymin", dest="y0",type = "float", default=-90.0 ,
                      help = "Minimum y  [default = -90.0]." )
    parser.add_option("-Y", "--ymax", dest="y1",type = "float", default=90.0 ,
                      help = "Maximum y  [default = 90.0]." )
    parser.add_option("-z", "--zmin", dest="z0",type = "float", default=0.0 ,
                      help = "Minimum depth (km) [default = 0.0]." )
    parser.add_option("-Z", "--zmax", dest="z1",type = "float", default=700.0 ,
                      help = "Maximum depth (km) [default = 1000.0]." )

    # number of bins
    parser.add_option( "--nx", dest="nx", type = "int", default=10 ,
                       help = "Number of bins in x direction [default = 10].")
    parser.add_option( "--ny", dest="ny", type = "int", default=10 ,
                       help = "Number of bins in y direction [default = 10].")
    parser.add_option( "--nz", dest="nz", type = "int", default=10 ,
                       help = "Number of bins in down direction [default = 10].")

    # type of summed tensor
    parser.add_option( "-t", "--sumtype", dest="stype",type = "int", default=0,
                       help = """Output option, (0=summed moment tensor (default) | \
                          1= summed source mechanism tensor)""" )
    # output type
    parser.add_option("-o", "--otype", dest="otype", type="int", default=0,
                      help="""Output format : \
                           0=[cx/cy/cz/mrr/mtt/mff/mrt/mrf/mtf/exp/x/y/z/n], \
                           1=[cx/cy/cz/mxx/myy/mzz/myz/mxz/mxy/M0/x/y/z/n], \
                           2=[cx/cy/cz/T/tx/ty/tz/B/bx/by/bz/P/px/py/pz/M0/x/y/z/n], \
                           3=[cx/cy/cz/T/ta/tp/B/ba/bp/P/pa/pp/exp/x/y/z/n], \
                           4=[cx/cz/cy/myy/mzz/mxx/-myz/mxy/-mxz/exp/x/y/z/n(for gmt x-section)] \
                           where [cx,cy,cz] is the weighted avg location, [x,y,z] is the bin center location; \
                           r/z=up,t=S,f/x=E,y=N """)

    # surpress output
    parser.add_option("-q", "--quiet",
                      action="store_false", dest="verbose", default=True,
                      help="don't print status messages to stderr")

    return parser

#----------------------------------------------------------------------
def genRandMT():
    """
    Function generates the input parser for the genRand program
    """
    parser = OptionParser()

    # define region
    parser.add_option("-x", "--xmin", dest="x0",type = "float", default=0.0 ,
                      help = "Minimum x  [default = 0.0]." )
    parser.add_option("-X", "--xmax", dest="x1",type = "float", default=360.0 ,
                      help = "Maximum x  [default = 360.0]." )
    parser.add_option("-y", "--ymin", dest="y0",type = "float", default=-90.0 ,
                      help = "Minimum y  [default = -90.0]." )
    parser.add_option("-Y", "--ymax", dest="y1",type = "float", default=90.0 ,
                      help = "Maximum y  [default = 90.0]." )
    parser.add_option("-z", "--zmin", dest="z0",type = "float", default=0.0 ,
                      help = "Minimum depth (km) [default = 0.0]." )
    parser.add_option("-Z", "--zmax", dest="z1",type = "float", default=700.0 ,
                      help = "Maximum depth (km) [default = 1000.0]." )

    # number of moment tensors
    parser.add_option( "-n", "--number", dest="nRand", type = "int", default=100 ,
                       help = "Number of random tensors [100].")

    # b-value
    parser.add_option("-b", "--bval", dest="b",type = "float", default=0.0 ,
                      help = "b-value in G-R distribution [0.0]." )

    # flag for only deviatoric mechanisms
    parser.add_option("--dev", action="store_true", dest="isDev", default=False,
                      help="Only deviatoric mechanisms")

    # flag for only DC mechanisms
    parser.add_option("--dc", action="store_true", dest="isDC", default=False,
                      help="Only deviatoric, double-couple mechanisms")

    (opt, args)=parser.parse_args()

    return (opt, args)

#######################################################################    
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
# genParser.py ends here
