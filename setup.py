#!/usr/bin/env python
# setup.py --- 
# 
# Filename: setup.py
# Description: 
# Author: Iain Bailey
# Maintainer: 
# Created: Fri Dec 16 15:10:30 2011 (-0800)
# Version: 
# Last-Updated: Wed Feb  1 14:44:38 2012 (-0800)
#           By: Iain Bailey
#     Update #: 56
# URL: 
# Keywords: 
# Compatibility: 
# 
# 

# Change Log:
# 
# 
# 
# 
# 
# 

# Code:

from distutils.core import setup, Extension
from distutils.sysconfig import get_python_lib, get_config_vars
import sys, shutil, os, os.path, time, glob, stat

######################################################################

BOOST_BASE = None
BOOST_LIB = "boost_python"
BOOST_DLL = "boost_python.dll"

# Name of the entire package (you can't use just this variable to rename
# the package, all the directories and source files have to be updated
# accordingly)
PACKAGE_NAME = "pythmt"

# The following variables are used to customize the compile process

INC_DIRS = []
LIB_DIRS = []
LIBS     = []
CC_ARGS  = []
LINK_ARGS = []
MACROS    = []
data_files = []
scripts = []

# list of scripts
scripts = ["scripts/mtRandom.py"
           , "scripts/mtFaultStyle.py"
           , "scripts/mtInfo.py"
           , "scripts/mtBinSum.py"
           , "scripts/mtDecomp.py"
           , "scripts/mtSum.py" 
           , "scripts/mtConvert.py"  # converts between different representations 
           , "scripts/psmecaSmtoSx.py"
            ]

setup( name = PACKAGE_NAME,
       version = "2.0",
       description = "Python Moment Tensors",
       author = "Iain Bailey",
       author_email = "iainbailey@gmail.com",
       url = "http://earth.usc.edu/~iwbailey",
       license = "GPL",

       # Process all pure Python modules in the source directory
       packages = [PACKAGE_NAME ],

       # add scripts
       scripts = scripts  )


print ("... finished setup")

# ***** BEGIN LICENSE BLOCK *****
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
# ***** END LICENSE BLOCK *****
# 
# setup.py ends here
