#!/usr/bin/env python
# momenttensor.py --- 
# 
# Filename: momenttensor.py
# Description: 
# 
# Python class for symmetric moment tensor and associated functions
#
# Author: IW Bailey
# Maintainer: IW Bailey
# Created: Fri Nov  5 10:06:42 2010 (-0700)
# Version: 1
# Last-Updated: Thu Dec 29 18:05:40 2011 (-0800)
#           By: Iain Bailey
#     Update #: 603

# Commentary:
#
#
#
#

# Change Log:
#
# Wed Dec 21 2011 : Removed strike and dip stuff, since added double
# couple as separate entity
#
# Tue Nov 15 2011 : mrt had a - where it should have had a +
#
#
#
#

# Code:

import numpy as NP
from numpy import sin, cos
from  math import sqrt, log10, atan2, pi
import sys

# Constants
EPS = 1e-8 # effective zero for accuracy checks



# FUNCTIONS ##################################################################

def azimplunge( x ):
    """ 
    Get the azimuth and plunge in degrees of a vector where x=[r,t,p]
    where r=up, t=South, p=East
    """
    if x[0]>0 : x*=-1;  # make downward pointing

    xh = sqrt(x[1]**2 + x[2]**2) # get horizontal length
    plunge = 180.0*atan2( -x[0], xh )/pi  # plunge should be positive
    azim = 180*atan2(x[2], -x[1])/pi # angle cw from n

    return (azim, plunge)

#---------------------------------------------------------------------
def voigt2mat( v ):
    """
    convert from vector with voigt notation to a matrix
    """
    return NP.array( [
            [ v[0], v[5], v[4] ],
            [ v[5], v[1], v[3] ],
            [ v[4], v[3], v[2] ] ])

#---------------------------------------------------------------------
def mat2voigt( M ):
    """
    convert from symmteric matrix as Numpy array to vector with voigt
    ordering
    """
    return NP.array( [M[0,0], M[1,1], M[2,2], M[1,2], M[0,2], M[0,1]] )

#---------------------------------------------------------------------
def voigtnorm( v ):
    """
    Get the norm of a 6 component matrix in voigt form
    """
    
    # diagonals once, off diagonals twice 
    return sqrt( NP.sum(v*v) + NP.sum(v[3:6]*v[3:6]) )

#---------------------------------------------------------------------
def mateig( mat ):
    """
    Get the eigen solution of a 3x3 moment tensor matrix 
    given as the ptb matrix and eigen values
    """
  
    # get min, int, max eigen vals, vects in that order
    (vals, vecs) = NP.linalg.eigh( mat , UPLO='U')

    # sort the eigs from low to high
    ind = NP.argsort(vals)
    vals = vals[ind]
    vecs = vecs[:, ind]

    # prescribe orthogonality, b = t x p
    vecs[:,1] = NP.cross( vecs[:,2] , vecs[:,0] )

    return ( vecs , vals )

#---------------------------------------------------------------------
def mag2m0( mag ):
    """
    Use Hanks and Kanamori reln to compute moment
    """
    return 10**(1.5*(mag+10.7))

#---------------------------------------------------------------------
def m02mag( m0 ):
    """
    Use Hanks and Kanamori reln to compute magnitude
    """
    return (log10(m0)/1.5) - 10.7 



###############################################################################

class SymMT:
    """
    Define a class for a symmetric point source moment tensor.

    Instantiate with...

    mt = SymMT( m, M0, centroid , hypo, cov )

    where

    m = normalized symmetric tensor (source mech tensor) using voigt
      notation in the form [ mrr, mtt, mpp, mtp, mrp, mrt ] and r =
      Up, t=South and p=East

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
                  , Norm = None  # norm of tensor, such that M = Norm*smt
                  , MT = None  # 3x3 numpy array option 
                  , c = NP.array( [0.0, 0.0, 0.0] )  # centroid location
                  , h = NP.array( [0.0, 0.0, 0.0] )  # hypocenter location
                  , cov = NP.zeros( 36 ) # covariance matrix
                  , mag = None
                  , rigidity = 30 # rigidity of source region (GPa)
                  ):

        if ( MT != None ):
            # initiate from the numpy matrix
            __frommatrix( MT )
        else:
            # use the vector form
            if Norm == None:
                # case where we normalise
                self.Norm = voigtnorm( mhat )
                self.mhat = mhat/self.Norm # mxx, myy, mzz, myz, mxz, mxy
            else:
                # force normalisation
                self.mhat = mhat/voigtnorm( mhat )
                self.Norm = Norm

        # magnitude not necessarily linked to moment
        if mag != None:
            self.mag = mag
        else:
            self.mag = m02mag( sqrt(2.0)*self.Norm )

        self.c = c # centroid lon lat depth
        self.h = h # hypocenter lon lat depth

        self.cov = cov # covariance matrix for M

        self.rigidity = rigidity # ridgidity in GPa
        
    #--------------------------------------------------------------------
    def __frommatrix( Mmat, c, h):
        """
        Convert a numpy array into a SymMT object
        """
        
        # convert to vector with voigt notation
        m = mat2voigt(Mmat)

        # get the norm
        norm = voigtnorm(m)

        # make the tensor object
        self.mhat = m/norm
        self.Norm = norm

        return 

    # --------------------------------------------------
    def smt( self, i=None, j=None ):
        """
        Get full Source mech tensor or ij component of source
        mechanism tensor, where i=1,2,3
        """
        if( i==None and j==None):
            # no indices defined, return a numpy matrix
            return voigt2mat( self.mhat )
        else:
            # else return a component
            if( i == j ): return self.mhat[ i-1 ]
            else: return self.mhat[ 6 - i - j + 2 ]

    # --------------------------------------------------
    def M( self, i=None, j=None ):
        """ 
        Get full moment tensor or ij component of moment tensor; i,j =1,2,3
        """
        if( i==None and j==None):
            # no indices defined, return a numpy matrix
            return self.Norm*voigt2mat( self.mhat )
        else:
            return self.Norm * self.smt( i, j )

    # --------------------------------------------------
    def M0( self ):
        """
        Get scalar moment as defined by Silver and Jordan
        """
        return self.Norm / sqrt(2.0)

    # --------------------------------------------------
    def Mw( self ):
        """
        Get moment magnitude
        """
        return m02mag( self.M0() )

    # --------------------------------------------------
    def eig( self ):
        """
        Return the eigenvectors [P, B, T] as columns of matrix,
        rows will be r, theta, phi components.
        eigenvalues as a horizontal vector with same order
        (1st column will be the smallest eigenvalue)
        """
        (vecs, vals) = mateig( voigt2mat(self.mhat) )

        # return vectors and values
        return ( vecs , vals*self.Norm )

    # --------------------------------------------------
    def pbt( self ):
        """
        Return the Eigen solution in manageable parts

        OUT
        pval = eigen value for p axis
        p = numpy array for p axis vector (up, South, East coord system)
        bval, b, tval, t = etc.
        """

        # get the eigenvalues and vectors
        (V,D) = self.eig()

        # extract 3 component r, theta, phi vectors
        p, b, t = V[:,0], V[:,1], V[:,2] # phi

        # get the values
        pval, bval, tval = D[0], D[1], D[2]

        return pval, p, bval, b, tval, t

    # --------------------------------------------------
    def iso_smt( self ):
        """
        Get the isotropic component of the normalised tensor
        """
        return NP.sum( self.mhat[0:3] ) / 3

    # --------------------------------------------------
    def iso( self ):
        """
        Get the isotropic component of the tensor
        """
        return self.Norm * self.getMhatIso() 

    # --------------------------------------------------
    def sepDevIso( self ):
        """
        Separate the deviatoric and the components of the tensor
        """
        EPS=1e-5;
        iso = self.getIso()
        norm = sqrt( 3*(iso**2) ) 
        if( norm < EPS ):
            MTiso = SymMT( NP.array([ 0,0,0, 0,0,0]) , 
                           0.0, self.c, self.h )
        else:
            MTiso = SymMT( NP.array([ iso, iso, iso, 0,0,0])/norm , 
                           norm, self.c, self.h )

        mdev = self.getMvec() - NP.array([ iso, iso, iso, 0,0,0])
        norm = NP.sqrt( sum(mdev**2) +  sum(mdev[3:]**2) )
        if( norm < EPS ):
            MTiso = SymMT( NP.array([ 0,0,0, 0,0,0]) , 
                           0.0, self.c, self.h )
        else:
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
    def f_clvd( self ):
        """
        return th fCLVD value
        """
        (V,D) = self.eig()

        D -= NP.sum(D)/3 # get deviatoric

        return -1.0 * D[1] / max( [ abs(D[0]), abs(D[2]) ] )

    # --------------------------------------------------
    def Gamma( self ):
        """ 
        Return Kagan's gamma value (esentially the determinant) which
        measures the CLVD size.
        See ...

        Frohlich (1995), "Characteristics of well determined... " PEPI 
        Kagan (2009), "On the geometric complexity..." PEPI

        This was incorrectly defined by a factor of -0.5 in Bailey et al (2009) GJI
        """
        M = self.smt()
        
        M -= self.iso_smt()*NP.eye(3) # get deviatoric

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
        mrr, mtt, mff = self.mhat[0], self.mhat[1], self.mhat[2]
        mtf, mrf, mrt = self.mhat[3], self.mhat[4], self.mhat[5]

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

###############################################################################

class EigMT:
    """
    class that stores moment tensor in its eigen decomposition format
    """
    def __init__(self
                 , m = [-1, 0, 1, 0, 0, 0]  # moment tensor in vector form
                 , p = None # 3x1 numpy array for p-axis vector
                 , t = None # 3x1 numpy array for t-axis vector
                 , pbtvals = None # 1x3 numpy array 
                 , c = None # centroid location
                 , h = None  # hypocenter location
                 ):
        """
        Constructor for moment tensor stored in its eigen decomposition
        """ 

        # set orientation and size part
        if ( p != None and t !=None, pbtvals != None ):
            # Set directly from arguments

            # check p and t are orthogonal
            if( NP.abs( NP.dot( p,t ) ) > EPS ):
                print >> sys.stderr, "P and T-axes are not orthogonal"
                sys.exit()

            # make the orientation part
            b = NP.cross( t, p ) # make such that p x b = t 
            self.pbt = NP.c_[ p, b, t ] # numpy 3x3 array

            # size part
            self.pbtvals = pbtvals  # nupy 3x1 array

        else:
            # convert from matrix vector representation
            Mij = voigt2mat( m )
            (self.pbt, self.pbtvals) = mateig( Mij )

        # set position
        if c == None and h == None :
            # default position 
            h = NP.array( [0.0, 0.0, 0.0] )
        if c == None: c = h
        if h == None: h = c
        self.c = c
        self.h = h

     # --------------------------------------------------
    def p( self, i=None ):
        """
        Get the vector for the p-axis or a component
        """
        if( i==None ):
            # no indices defined, return a numpy array
            return self.pbt[:,0]
        else:
            return self.pbt[i,0]

    # --------------------------------------------------
    def b( self, i=None ):
        """
        Get the vector for the b-axis or a component
        """
        if( i==None ):
            # no indices defined, return a numpy array
            return self.pbt[:,1]
        else:
            return self.pbt[i,1]

    # --------------------------------------------------
    def t( self, i=None ):
        """
        Get the vector for the t-axis or a component
        """
        if( i==None ):
            # no indices defined, return a numpy array
            return self.pbt[:,2]
        else:
            return self.pbt[i,2]

    # --------------------------------------------------
    def M( self, i=None, j=None ):
        """ 
        Get full moment tensor or ij component of moment tensor; i,j =0,1,2
        """
        M = self.pbt * NP.diag( self.pbtvals ) * self.pbt.transpose()

        if( i==None and j==None):
            # no indices defined, return a numpy matrix
            return M
        else:
            return M[ i, j ]

    # --------------------------------------------------
    def smt( self, i=None, j=None ):
        """ 
        Get source mechanism tensor or ij component of moment tensor; i,j =0,1,2
        """
        M = self.pbt * NP.diag( self.pbtvals ) * self.pbt.transpose()

        if( i==None and j==None):
            # no indices defined, return a numpy matrix
            return M
        else:
            return M[ i, j ]


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
# momenttensor.py ends here
