import sys
import numpy as NP
from math import sqrt, pi
from MomentTensor import SymMT, sdr2mt

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
    tmp = NP.array( thisline.split() )

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
    m = NP.array( [mrr, mtt, mff, mtf, mrf, mrt] )/norm 
    
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
def readSDRfile( istream ):
    """
    From an input file or stdin, read a list of focal mechanisms in
    strike/dip/rake form

    lon/lat/depth/strike/dip/rake/mag
    """
    lcount = 0
    mtlist = []

    # get data in numpy array
    data = NP.genfromtxt( istream, usecols=(0, 1, 2, 3, 4, 5, 6) ) #, autostrip=True )

    nline = NP.size( data, 0 )
    for i in range(0,nline):

        # get variables
        hypo = data[i, 0:3 ] # lon , lat, z
        sdr = data[i, 3:6]*pi/180 # in radians
        mag = data[i, 6]

        #TODO: error checking 

        # Make a moment tensor object
        MT = SymMT( c = hypo, h = hypo,
                    strike = sdr[0], dip = sdr[1], rake=sdr[2],
                    mag = mag )

        mtlist.append( MT )
        
    return mtlist
#--------------------------------------------------