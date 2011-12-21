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
def binMT_xyz( mtlist # list of moment tensors
               , x1=0.0, x2=360.0, nx=10
               , y1=-90, y2=90.0, ny=10
               , z1=0.0, z2=700, nz=10
               , binmag = 0.0 ):
    """
    Bin a list of moment tensors into an evenly spaced 3-d grid
    """

    # get bin coordinates
    dx = (x2 - x1)/float(nx)
    binx = np.linspace( x1+0.5*dx, x2-0.5*dx, nx )
    dy = (y2 - y1)/float(ny)
    biny = np.linspace( y1+0.5*dy, y2-0.5*dy, ny )
    dz = (z2 - z1)/float(nz)
    binz = np.linspace( z1+0.5*dz, z2-0.5*dz, nz )

    # get the coordinates from the list of moment tensors
    ndat = len(mtlist)
    pos = np.zeros( [ndat, 3] )
    for i in range(0,ndat):
        pos[i,:] = mtlist[i].c

    # convert to bin indices
    xi = np.floor( nx*( pos[:,0] - x1 )/(x2 - x1) )
    yi = np.floor( ny*( pos[:,1] - y1 )/(y2 - y1) )
    zi = np.floor( nz*( pos[:,2] - z1 )/(z2 - z1) )
    xyzi = xi*(ny*nz) + yi*nz + zi

    # make bins
    binlist = []
    binnlist = []
    ntot = 0
    for i in range(0,nx):
        for j in range(0,ny):
            for k in range(0,nz):
                # find the indices of the mtlist
                idx = np.where( np.logical_and( np.logical_and(  xi == i , 
                                                                 yi == j  ), 
                                                zi == k ) )
                nbini = len( idx[0] )
                binnlist.append( nbini )
                ntot += nbini
        
                # append to bins
                if nbini == 0: 
                    binlist.append( None )
                else:
                    thismtlist = []
                    for m in idx[0][:]: thismtlist.append( mtlist[m] )
                    thisbin = eqkBin( xyz=np.array([ binx[i], biny[j], binz[k] ] ), 
                                      mag=binmag, eqklist=thismtlist ) 
                    binlist.append( thisbin )

    return (binlist, ntot)

##--------------------------------------------------
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



        
