
from math import pi
from numpy import sin, arcsin as asin
from sys import stderr
import numpy as NP
import numpy.random as NPR


#--------------------------------------------------
def randLonLat( n=1, lon1=-180.0, lon2=180.0, lat1=-90.0, lat2=90.0 ):
    """
    Get n random lon/lat points in degrees
    """
    # deg/radian conversion
    D2R=pi/180.0
    R2D=1/D2R

    # get lons
    lons = lon1 + (lon2-lon1)*NPR.rand( n, 1 )

    # correct for area and get rand lat
    z1, z2 = sin( D2R*lat1 ) , sin( R2D*lat2 )
    lats = R2D*asin( z1 + (z2-z1)*NPR.rand( n,1 ) )
    
    return (lons, lats)

#--------------------------------------------------
def randSphere( n=1 , dim=3):
    """
    Get n random vector points within a unit sphere
    """
    
    # get starting vector
    vec = 2*NPR.rand( n, dim )- 1

    norm2 = NP.sum( vec**2, 1 ) 
    
    # get within unit sphere    
    for i in range(0,n):
        while norm2[i] > 1.0:
            vec[i,:] = 2*NPR.rand( 1, dim )- 1
            norm2[i] = NP.sum( vec[i,:]**2 )
            
        if( i  % 10000 == 0 ): stderr.write( "finished %i\n" %  i )

    return (vec, norm2)

#--------------------------------------------------
def randSphereSurf( n=1 , dim=3):
    """
    Get n random vector points on surface of a unit sphere
    """
    (vec, norm2) = randSphere( n , dim)

    # normalise
    norm = NP.sqrt(norm2)
    for i in range(0,dim): vec[:,i] /= norm.transpose()

    return vec

#--------------------------------------------------
def randPBT( n=1 ):
    """
    Get n random vector points on surface of a unit sphere
    """

    # get p
    p = randSphereSurf( n, 3 ) 

    # use dot product
    txy = randSphereSurf( n, 2 ) 
    t = NP.c_[ txy[:,0], txy[:,1], - (p[:,0]*txy[:,0] + p[:,1]*txy[:,1])/p[:,2] ] 

    t /= NP.tile( NP.sqrt( NP.sum( t**2 , 1) ), (3,1) ).transpose()

    # cross product
    b = NP.c_[ t[:,1]*p[:,2]-t[:,2]*p[:,1] , 
               -(t[:,0]*p[:,2] - t[:,2]*p[:,0]) , 
               t[:,0]*p[:,1] - t[:,1]*p[:,0] ]

    return (p, b, t)
