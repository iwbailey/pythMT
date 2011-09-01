import sys
import numpy as NP
from math import sqrt
from MomentTensor import SymMT

#--------------------------------------------------
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

#--------------------------------------------------
def readPsmecaList( istream ):
    """
    From an input file or stdin, read a list of moment tensors in psmeca form
    """
    lcount = 0
    mtlist = []
    endstrlist = []

    for line in istream:
        (mt, endstr) = readpsmecaSm( line, lcount )
        mtlist.append( mt )
        lcount += 1

    return mtlist

#--------------------------------------------------
def readSDRList( istream ):
    """
    From an input file or stdin, read a list of moment tensors in psmeca form
    """
    for line in istream:
        (mt, enstr) = readSDR( line, lcount )
        mtlist.append( mt )
        lcount += 1

    return 1
