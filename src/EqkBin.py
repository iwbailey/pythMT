# Functions that take a list of moment tensors and but them in bins
import numpy as np

#--------------------------------------------------
def binMT_mag( mtlist, m1=0.0, m2=10.0, nm=10, binlocn=np.array([0.0, 0.0, 0.0]) ):
    """
    Bin a list of moment tensors into equally spaced magnitude bins
    """

    # get bin locations in magnitude
    dm = (m2 - m1)/float(nm)
    bincenters = np.linspace( m1+0.5*dm, m2-0.5*dm, nm )

    # get the magnitudes from the list of moment tensors
    ndat = len(mtlist)
    mags = np.zeros( ndat )
    for i in range(0,ndat):
        mags[i] = mtlist[i].mag

    # convert to bin indices
    magi = np.floor( nm*( mags - m1 )/(m2 - m1) )

    # make bins
    binlist = []
    binnlist = []
    ntot = 0
    for i in range(0,nm):
        
        # find the indices of the mtlist
        idx = np.where( magi == i )
        nbini = len( idx[0] )
        binnlist.append( nbini )
        ntot += nbini
        
        # append to bins
        if nbini == 0: 
            binlist.append( None )
        else:
            thismtlist = []
            for j in idx[0][:]: thismtlist.append( mtlist[j] )
            thisbin = eqkBin( xyz=binlocn, mag=bincenters[i], eqklist=thismtlist ) 
            binlist.append( thisbin )

    return (binlist, ntot)

#--------------------------------------------------
class eqkBin:
    """
    A bin of earthquake data
    """
    def __init__( self
                  , xyz = np.array( [0.0, 0.0, 0.0] ) # 3-D bin center location
                  , mag = 0.0 # bin magnitude
                  , eqklist = [] # list of earthquake objects
                  ):
        self.xyz = xyz
        self.mag = mag
        self.eqklist = eqklist
        self.n = len(eqklist) # number of earthquakes



        
