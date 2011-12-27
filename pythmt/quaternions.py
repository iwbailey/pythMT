
import numpy as np
import sys
from cgkit.cgtypes import quat as cgquat, mat3 
from math import sqrt

class quatn( cgquat ):
    """
    Quaternion class that basically links with the cgkit 
    """
    def __init__( self , q=[1,0,0,0], A=None ):

        if ( A != None ):
            # case where we have provided a numpy array as a rotation matrix

            # convert to mat3
            A = mat3( [ A[0,0], A[0,1], A[0,2], 
                        A[1,0], A[1,1], A[1,2], 
                        A[2,0], A[2,1], A[2,2] ] )

            # initialise the quaternion
            cgquat.__init__( self, A )

        else:
            # default is the four elements
            cgquat.__init__( self, q )

    #----------------------------------------------------------------------
    def __frommatrix( self, A ):
        """
        Convert from rotation matrix to quaternion
        """
        t = np.sum( np.diag(A) ) + 1.0 
    
        if( t > sys.float_info.epsilon ):
             s = 0.5/sqrt(t)
             w = 0.25 / s
             x1 = ( A[2,1] - A[1,2] ) * s
             x2 = ( A[0,2] - A[2,0] ) * s
             x3 = ( A[1,0] - A[0,1] ) * s
        else:
            if ( A[0,0] > A[1,1] and A[0,0] > A[2,2] ):
                s = 2.0 *sqrt( 1.0 + A[0,0] - A[1,1] - A[2,2] )
                x1 = 0.25 * s
                x2 = (A[0,1] + A[1,0]) / s
                x3 = (A[0,2] + A[2,0]) / s
                w = (A[2,1] - A[1,2]) / s
            elif (A[1,1] > A[2,2]):
                s = 2.0 *  sqrt( 1.0 + A[1,1] - A[0,0] - A[2,2] )
                x1 = (A[0,1] + A[1,0] ) / s
                x2 = 0.25 * s
                x3 = (A[1,2] + A[2,1] ) / s
                w = (A[0,2] - A[2,0] ) / s
            else:
                s = 2.0 *  sqrt( 1.0 + A[2,2] - A[0,0] - A[1,1] )
                x1 = (A[0,2] + A[2,0]) / s
                x2 = (A[1,2] + A[2,1]) / s
                x3 = 0.25 * s
                w = (A[1,0] - A[0,1]) / s

        # assign
        cgquat.__init__( self, w, x1, x2, x3 )
    #----------------------------------------------------------------------

