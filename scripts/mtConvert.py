#!/usr/bin/env python
# 
# Filename: mtConvert.py
# Description: Convert between different moment tensor representations
# Author: Iain William Bailey
# Created: Tue Dec 27 10:24:02 2011 (-0800)
# Version: 1
# Last-Updated: Fri Dec 30 14:48:13 2011 (-0800)
#           By: Iain Bailey
#     Update #: 172
# Compatibility: 2.7.2+
# 
# 

# Change Log:
# 
# 
# 

# Code:

# Standard libraries used
import sys
import argparse 
from optparse import OptionParser
from math import log10, sqrt, pi
import numpy as NP

# Personal libraries used
import pythmt.iofuncs as iofuncs
from pythmt.momenttensor import mat2voigt

# constants 
R2D = 180/pi


# FUNCTIONS ############################################################

# input options all take the input file as an argument
#----------------------------------------------------------
def in_mt( ifile ):
    """
    Read in from moment tensor format, return EigMT's
    """
    print "MT input not yet added."
    return None

#----------------------------------------------------------
def in_pbt( ifile ):
    """
    Read in from PBT format, return EigMT's
    """
    print "P-B-T input not yet added."
    return None

#----------------------------------------------------------
def in_sdr( ifile ):
    """
    Read in from strike, dip, rake format
    """
    mtlist = iofuncs.read_sdrlist( ifile )
    return mtlist

#----------------------------------------------------------
def in_quat( ifile ):
    """
    Read in from Quaternion format
    """
    print "Quaternion input not yet added."
    return None

#----------------------------------------------------------

# output options take a list of moment tensors or double couples
#----------------------------------------------------------
def out_mt( mtlist, ofile ):
    """
    Write to moment tensor format
    """
    for mt in mtlist: 
        mt.c.tofile( ofile, sep=" ", format='%-9.4f')
        ofile.write(' ')
        mat2voigt(mt.M()).tofile( ofile, sep=" ", format='%8.4g') 
        ofile.write('\n')

    return

#----------------------------------------------------------
def out_pbt( mtlist, ofile ):
    """
    Write to P-B-T format
    """
    for mt in mtlist: 
        mt.c.tofile( ofile, sep=" ", format='%-9.4f')
        for j in (0,1,2):
            ofile.write(' %8.4g' % mt.pbtvals[j] )
            mt.pbt[:,j].tofile( ofile, sep=" ", format='%8.5f')
        ofile.write('\n');
    return

#----------------------------------------------------------
def out_sdr( dclist, ofile ):
    """
    Write in strike-dip-rake format
    """
    # TODO: add this
    print "Strike-dip rake output not yet added."
    return

#----------------------------------------------------------
def out_quat( dclist, ofile ):
    """
    Write quaternion output, only applies for double-couples
    """
    # TODO: add if for converting from moment tensor to double couple

    for dc in dclist: 
        dc.h.tofile( ofile, sep=" ", format='%-9.4f')
        for q in  dc.quats(): 
            ofile.write( " %-8.5f %8.5f %8.5f %8.5f" % (q.w, q.x, q.y, q.z) )
        ofile.write(" %8.4g\n" % dc.m0 )

    return
#----------------------------------------------------------
def out_rotax( dclist, ofile ):
    """
    Write rotation axis output, only applies for double-couples
    """
    # TODO: add if for converting from moment tensor to double couple

    for dc in dclist: 
        dc.h.tofile( ofile, sep=" ", format='%-9.4f')
        for q in  dc.quats(): 
            (a, x) = q.toAngleAxis()
            ofile.write( " %8.3f %8.5f %8.5f %8.5f" % (R2D*a, x[0], x[1], x[2]) )
        ofile.write(" %8.4g\n" % dc.m0 )

    return


#------------------------------
# command line arguments
parser = argparse.ArgumentParser( description='Convert between different moment tensor representations' )

parser.add_argument('ifile', nargs='?', type=argparse.FileType('r'), default=sys.stdin, 
                      help = "Input file [stdin]")

parser.add_argument('ofile', nargs='?', type=argparse.FileType('w'), default=sys.stdout,
                    help = "Output file [stdout]")

optstring = "0=moment tensor| 1=Pp,Bb,Tt| 2=strike,dip,rake(degs)| 3=quaternion| 4=rotn axis"

parser.add_argument( '-i', '--itype', type=int, default=0, 
                     help = 'Representation for input: '+optstring )

parser.add_argument( '-o', '--otype', type=int, default=0, 
                     help = "Representation for output: "+optstring ) 

args = parser.parse_args()
#------------------------------

in_options = {0: in_mt, 1: in_pbt, 2: in_sdr, 3: in_quat}  
out_options = {0: out_mt, 1: out_pbt, 2: out_sdr, 3: out_quat, 4: out_rotax}  

# get the input
try: 
    mtlist = in_options[args.itype](args.ifile)
except KeyError: 
    print "Input type (itype=%i) not recognised" % args.itype
    sys.exit()

if mtlist == None: sys.exit()

# output
try: 
    out_options[args.otype](mtlist, args.ofile)
except KeyError: 
    print "Output type (otype=%i) not recognised" % args.otype
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


# 
# mtConvert.py ends here
