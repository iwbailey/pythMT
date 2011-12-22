
from cgkit.cgtypes import quat as cgquat, mat3 
import sys

class quatn( cgquat ):
    """
    Quaternion class that basically links with the cgkit 
    """
    def __init__( self , q=[1,0,0,0], A=None ):

        # case where we have provided a numpy array as a rotation matrix
        if ( A != None ):
            A = 1.0*A # make sure float
            # convert to mat3
            A = mat3( [ A[0,0], A[0,1], A[0,2], 
                        A[1,0], A[1,1], A[1,2], 
                        A[2,0], A[2,1], A[2,2] ] )
            print "A", A
            cgquat.__init__( self, A )
        else:
            # default is the four elements
            cgquat.__init__( self, q )
