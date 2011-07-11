#!/usr/bin/env python
# MomentTensor.py --- 
# 
# Filename: MomentTensor.py
# Description: Python class for symmetric moment tensor
# Author: IW Bailey
# Maintainer: IW Bailey
# Created: Fri Nov  5 10:06:42 2010 (-0700)
# Version: 1
# Last-Updated: Sun Jun  5 13:02:31 2011 (-0700)
#           By: Iain Bailey
#     Update #: 356

# Commentary:
#
#
#
#

# Change Log:
#
#
#
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
#

# Code:

import numpy as NP
from  math import sqrt, log10, atan2, pi
import sys

##################################################
def readpsmecaSm( thisline , lcount=1):
    """
    Take a string as an input argument.  The string should be a psmeca
    entry.  Return a moment tensor object

    expected input file: 
    x, y, z, mrr, mtt, mff, mrt, mrf, mtf, exp
    r = up, t is south, f is east

    the last three parts (lon', lat', id) and anything after are
    output as a string
    """

    # Get matrix of floats without spaces
    # this will crash or get things wrong if you have tabs ('\t') with no spaces
    tmp = NP.array( thisline.split(' ') )
    tmp = tmp[ NP.where( tmp != "" ) ]

    # check the number of columns
    nc = 9 # minimum number of required columns
    if len( tmp ) < nc :
        print >> sys.stderr, ( 'Line %i of input. Expecting >= %i cols, got %i' %
                               ( lcount, nc, len(tmp) ) )
        raise IOError(69,'Line %i of input. Expecting >= %i cols, got %i' %
                      ( lcount, nc, len(tmp) ) )
        return


    # get coordinates, assume they represent the centroid
    centroid = NP.array( [ float(tmp[0]), float(tmp[1]), float(tmp[2]) ] )

    # get mechanism values
    mrr, mtt, mff, mrt, mrf, mtf = float(tmp[3]), float(tmp[4]), float(tmp[5]), \
        float(tmp[6]), float(tmp[7]), float(tmp[8])
    expon = float(tmp[9]) # exponent

    # get whatever is left as a string. reinsert the spaces 
    endstr = "";
    for i in range(10,len(tmp)): endstr += ' ' + tmp[i] 
    endstr = endstr.rstrip('\n')

    # get tensor norm without exponent
    norm = sqrt( mrr**2 + mtt**2 + mff**2 + 2*( mrt**2 + mrf**2 + mtf**2 ) )

    # Get the normalised tensor
    # transfer to xx, yy, zz, yz, xz, xy 
    # where x = east, y = north, z = up
    m = NP.array( [mff, mtt, mrr, -mrt, mrf, -mtf] )/norm 
    
    # Make a moment tensor object
    MT = SymMT( m, norm*(10**expon), centroid )

    return MT, endstr

##################################################
def azimplunge( x ):
    """ 
    Get the azimuth and plunge in degrees of a vector where x=[x,y,z]
    where x=E, y=N, z=up
    """
    if x[2]>0 : x*=-1;  # make downward pointing
    xh = sqrt(x[0]**2 + x[1]**2) # get horizontal length
    plunge = 180.0*atan2( -x[2], xh )/pi  # plunge should be positive
    azim = 180*atan2(x[0], x[1])/pi

    return (azim, plunge)

##################################################
def matrix2mt(Mmat, c, h):
    """
    Convert a numpy array into a SymMT object
    """
    # convert to vector with voigt notation
    m = NP.array([Mmat[0,0], Mmat[1,1], Mmat[2,2], Mmat[1,2], Mmat[0,2], Mmat[0,1]] )

    # get the norm
    norm = sqrt( NP.sum(m*m) + NP.sum(m[3:6]*m[3:6]) )

    # make the tensor object
    MT = SymMT( m/norm, norm, c, h )

    return MT

##################################################
class SymMT:
    """
    Define a class for a symmetric point source moment tensor.

    Instantiate with...

    mt = SymMT( m, M0, centroid , hypo, cov )

    where

    m = normalized symmetric tensor (source mech tensor) using voigt
      notation in the form [ mxx, myy, mzz, myz, mxz, mxy ] and x =
      East, y = North, z=up

    M0 = scalar moment [A.U], such that M_ij = sqrt(2) M0 m_ij

    centroid = centroid location [lon(deg), lat(deg), depth(km)] or
      [ x(km), y(km), z(km) ]

    hypocenter = hypocenter location [lon(deg), lat(deg), depth(km)] or
      [ x(km), y(km), z(km) ]


    cov = the covariance matrix of the full moment tensor sqrt(2)*M0*m
    ... this part will probably change... haven't used it yet

    If arguments are not given they will be set to default values
    """

    # --------------------------------------------------
    # Constructor
    def __init__( self
                  , mhat = NP.zeros( 6 )  # normalised tensor/source mech tensor
                  , Norm = 0.0  # norm of tensor, such that M = Norm*smt
                  , c = NP.array( [0.0, 0.0, 0.0] )  # centroid location
                  , h = NP.array( [0.0, 0.0, 0.0] )  # hypocenter location
                  , cov = NP.zeros( 36 ) # covariance matrix
                  ):

        self.mhat = mhat # mxx, myy, mzz, myz, mxz, mxy
        self.Norm = Norm # euclidean norm of moment tensor
        self.c = c # centroid lon lat depth
        self.h = h # hypocenter lon lat depth
        self.cov = cov # covariance matrix for M

    # --------------------------------------------------
    def smt( self, i, j ):
        """
        Get ij component of source mechanism tensor, where i=1,2,3
        """
        if( i == j ): return self.mhat[ i-1 ]
        else: return self.mhat[ 6 - i - j + 2 ]

    # --------------------------------------------------
    def M( self, i, j ):
        """ 
        Get ij component of moment tensor; i,j =1,2,3
        """
        return self.Norm * self.smt( i, j )

    # --------------------------------------------------
    def getMvec( self ):
        """
        Get the whole moment tensor in 6 component form
        [mxx, myy, mzz, myz, mxz, mxy]
        """
        return self.Norm * self.mhat

    # --------------------------------------------------
    def getSMTmat( self ):
        """
        Get the full 9 compnents of the source mechanism tensor
        [ mxx mxy mxz ]
        [ myx myy myz ]
        [ mzx mzy mzz ]
        """
        return NP.array( [
                [ self.mhat[0], self.mhat[5], self.mhat[4] ],
                [ self.mhat[5], self.mhat[1], self.mhat[3] ],
                [ self.mhat[4], self.mhat[3], self.mhat[2] ] ])

    # --------------------------------------------------
    def getM0( self ):
        """
        Get scalar moment
        """
        return self.Norm / sqrt(2.0)

    # --------------------------------------------------
    def getMw( self ):
        """
        Get moment magnitude
        MW = (2/3)*(log M0 - 16.1)
        """
#        return 2*(log10( self.getM0() )- 16.1)/3
        return 2*(log10( self.getM0() ))/3 - 10.7
    # --------------------------------------------------
    def getEig( self ):
        """
        Return the eigenvectors and eigenvalues
        P, B, T as columns of matrix
        eigenvalues as a horizontal vector with same order
        (1st column will be the smallest eigenvalue)
        """
        M = self.getSMTmat()
        
        # get min, int, max eigen vals, vects in that order
        (vals, vecs) = NP.linalg.eigh( M , UPLO='U')

        # sort the eigs from low to high
        ind = NP.argsort(vals)
        vals = vals[ind]
        vecs = vecs[:, ind]

        # prescribe orthogonality, b = t x p
        vecs[:,1] = NP.cross( vecs[:,2] , vecs[:,0] )

        # return vectors
        return ( vecs , vals*self.Norm )

    # --------------------------------------------------
    def getPTB( self ):
        """
        Return the Eigen solution in manageable parts

        OUT
        tval = eigen value for t axis
        t = numpy array for t axis vector (East-North-Up coord system)
        bval, b, pval, p = etc.
        """
        # get the eigenvalues and vectors
        (V,D) = self.getEig()
        p = V[:,0] 
        b = V[:,1]
        t = V[:,2]

        pval = D[0]
        bval = D[1]
        tval = D[2]

        return tval, t, bval, b, pval, p

    # --------------------------------------------------
    def getMhatIso( self ):
        """
        Get the isotropic component of the normalised tensor
        """
        return NP.sum( self.mhat[0:3] ) / 3

    # --------------------------------------------------
    def getIso( self ):
        """
        Get the isotropic component of the tensor
        """
        return self.Norm * self.getMhatIso() 

    # --------------------------------------------------
    def sepDevIso( self ):
        """
        Separate the deviatoric and the components of the tensor
        """
        iso = self.getIso()
        norm = sqrt( 3*(iso**2) ) 
        MTiso = SymMT( NP.array([ iso, iso, iso, 0,0,0])/norm , 
                       norm, self.c, self.h )

        mdev = self.getMvec() - NP.array([ iso, iso, iso, 0,0,0])
        norm = NP.sqrt( sum(mdev**2) +  sum(mdev[3:]**2) )
        MTdev = SymMT( mdev/norm , norm, self.c, self.h )

        return MTdev, MTiso


    # --------------------------------------------------
    def decompose( self , idecomp=1 ):
        """
        Decompose into 3 tensors that have the same eigen vectors
        idecomp = 1 CLVD B axis same as MT baxis
        idecomp = 2 CLVD P or T axis same as MT b-axis
        idecomp = 3 2 DCs 
        """

        (MTdev, MTiso) = self.sepDevIso()

        # get the eigen vectors
        (V,D) = MTdev.getEig()

        # clvd
        if idecomp == 1:
            # b-axis is clvd b-axis
            D2 = NP.array([ D[1], D[1], -2*D[1] ])
            D2 = D2[ NP.argsort(D2) ]
        elif idecomp == 2: 
            # b-axis is clvd p or t-axis
            D2 = NP.array([ -0.5*D[1], D[1], -0.5*D[1] ])
        elif idecomp == 3:
            # b-axis is alternative DC p or t-axis
            if( D[1] > 0 ): D2 = NP.array([ -D[1], D[1], 0.0 ])
            else: D2 = NP.array([ 0.0, D[1], -D[1] ])

        MT2  = matrix2mt( NP.dot( NP.dot(V, NP.diag( D2) ) , V.transpose() ),
                          self.c, self.h)

        # main dc
        Ddc = D - D2

        MTdc = matrix2mt( NP.dot( NP.dot(V, NP.diag( Ddc) ) , V.transpose() ),
                          self.c, self.h)

        return (MTdc, MT2, MTiso)
        
        

    # --------------------------------------------------
    def getRclvd( self ):
        """
        Return the rCLVD value
        """
        (V,D) = self.getEig()
        D -= self.getMhatIso() # get deviatoric
        
        eig2 = D[1]

        return 0.5*math.sqrt(6.0)*eig2 / sqrt( NP.sum( D**2 ) )

    # --------------------------------------------------
    def getFclvd( self ):
        """
        return th fCLVD value
        """
        (V,D) = self.getEig()

        D -= NP.sum(D)/3 # get deviatoric

        return -1.0 * D[1] / max( [ abs(D[0]), abs(D[2]) ] )

    # --------------------------------------------------
    def getGamma( self ):
        """ 
        Return Kagan's gamma value (esentially the determinant) which
        measures the CLVD size.
        See ...

        Frohlich (1995), "Characteristics of well determined... " PEPI 
        Kagan (2009), "On the geometric complexity..." PEPI

        This was incorrectly defined by a fact of -0.5 in Bailey et al (2009) GJI
        """
        M = self.getSMTmat()
        
        M -= self.getMhatIso()*NP.eye(3) # get deviatoric

        norm = NP.sqrt( NP.sum( M*M ) )
        det = NP.linalg.det(M)

        # note the 3 sqrt(6) differs from Frohlich (3/2) sqrt(3) b/c
        # of using the norm rather than the scalar moment
        gamma = 3*sqrt(6)*det/(norm*norm*norm);
 
        return gamma


    # --------------------------------------------------
    def getPsmecaSm( self ):
        """
        Return the gmt part for option -Sm
        longitude, latitude, depth in km, 
        mrr, mtt, mff, mrt, mrf, mtf in 10*exponent dynes-cm, exponent
        """ 
        # get location
        X, Y, depth = self.c[0], self.c[1], self.c[2]

        # get exponent
        exp = NP.round( log10( self.Norm ) )

        # get remaining factor
        fac =  self.Norm / 10**exp

        # get tensor components, r = up, t = south, f = east
        mrr, mtt, mff = self.mhat[2], self.mhat[1], self.mhat[0]
        mrt, mrf, mtf = -self.mhat[3], self.mhat[4], -self.mhat[5]

        return X, Y, depth, \
            fac*mrr, fac*mtt, fac*mff, fac*mrt, fac*mrf, fac*mtf, exp

    # --------------------------------------------------
    def getPsmecaSide( self ):
        """ 
        Return the cross section in gmt format
        t component points towards us (t->r), r component points up the page (r->-t)
        """
        X, Y, depth, mrr, mtt, mff, mrt, mrf, mtf, exp = psmecaSm( self )

        # convert so that where t is expected we give down (-r)
        # where r is expected we give N (-t)
        return X, Y, depth, mtt, mrr, mff, mrt, -mtf, -mrf, exp

    # --------------------------------------------------
    def getPsmecaSx( self ):
        """
        Return the gmt part for option -Sx
        longitude, latitude, depth in km,  
        value (in 10*exponent dynes-cm), azimuth, plunge of T, N, P axis, 
        exponent         
        """
        # get location
        X, Y, depth = self.c[0], self.c[1], self.c[2]

        # get the eigenvalues and vectors
        (V,D) = self.getEig()
        
        # convert vectors to azim and plunge
        (paz, ppl) = azimplunge( V[:,0] )
        (baz, bpl) = azimplunge( V[:,1] )
        (taz, tpl) = azimplunge( V[:,2] )

        # get exponent
        exp = NP.round( log10( self.Norm ) )

        # remove from eig values
        D = D / 10**exp

        return X, Y, depth, \
            D[2], taz, tpl, D[1], baz, bpl, D[0], paz, ppl, exp

    # --------------------------------------------------
    def rotate( self, R ):
        """
        Rotate the moment tensor using a rotation matrix
        """

        # get tensor as matrix
        M = self.SMTmat()

        # perform rotation
        M = dot( dot( R.transpose(), M ), R )

        # replace elements
        self.mhat = [ M[0,0], M[1,1], M[2,2], M[1,2], M[0,2], M[0,1] ]

        return

#     # Rotate the moment tensor around the x-axis so y moves towards z
#     def rotateX( self, theta ):
#         # convert theta to degrees
#         theta *= D2R

#         # make rotation matrix
#         R = NP.array([  [1,     0,          0           ],
#                         [0,     cos(theta), -sin(theta) ],
#                         [0,     sin(theta), cos(theta)  ]] )

#         # do the rotation
#         self.rotate( R )

#         return

#     # Rotate the moment tensor around the y-axis so z moves towards y
#     def rotateY( self, theta ):
#         # convert theta to degrees
#         theta *= D2R

#         # make rotation matrix
#         R = NP.array([  [cos(theta),    0, sin(theta)  ],
#                         [0,             1, 0           ],
#                         [-sin(theta),   0, cos(theta)  ]] )

#         # do the rotation
#         self.rotate( R )

#         return

#     # Rotate the moment tensor around the z-axis so x moves towards y
#     def rotateZ( self, theta ):
#         # convert theta to degrees
#         theta *= D2R

#         # make rotation matrix
#         R = NP.array([  [cos(theta),-sin(theta),0 ],
#                         [sin(theta),cos(theta), 0 ],
#                         [0,         0,          1 ]] )

#         # do the rotation
#         self.rotate( R )

#         return
# #
# MomentTensor.py ends here
