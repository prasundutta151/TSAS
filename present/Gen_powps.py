import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import sys
import pyfits
import os.path
from scipy.ndimage.filters import *
import numpy.ma as ma
from numpy.random import *
import scipy.ndimage.interpolation as intp

def FITS_READ(infile):

    file   = pyfits.open(infile)
    data   = file[0].data
    header = file[0].header
    dxy    = header['CDELT2']
    BMaj   = header['BMAJ']
    BMin   = header['BMIN']
    NGrid  = header['naxis1']
    ddim   = len(data.shape)
    data   = data.reshape(NGrid, NGrid) 

    return data, NGrid, dxy, BMaj, BMin

def FITS_WRITE(outfile, data, dxy, BMaj, BMin):

    hdu = pyfits.PrimaryHDU(data)
    hdulist = pyfits.HDUList([hdu])
    hdulist[0].header['CDELT1'] = dxy
    hdulist[0].header['CDELT2'] = dxy
    hdulist[0].header['BMAJ']   = BMaj
    hdulist[0].header['BMIN']   = BMin
    hdulist.writeto(outfile)
    
DTR = np.pi/180.
EXP_NARG=7

if(len(sys.argv)!=EXP_NARG):
    print
    print "USAGE: python/run", sys.argv[0], "<OUTFILE> <sigma> <alpha> <UMIN> <UMAX> <NBIN>"
    print "OUTFILE : Output file with power spectrum [U (klam), P(u)]"
    print "sigma   : Standard deviation of the field"
    print "alpha   : Slope of the power spectrum"
    print "Umin    : Value of Umin (Klambda)"
    print "Umax    : Value of Umax (Klambda)"
    print "NBin    : Number of bins"
    
    sys.exit()

OUTFILE  = sys.argv[1]
sigma    = np.float(sys.argv[2])
alpha    = np.float(sys.argv[3])
UMin     = np.float(sys.argv[4])
UMax     = np.float(sys.argv[5])
NBin     = np.int(sys.argv[6])

beta   = np.exp(np.log(UMax/UMin)/(1.*NBin))
xdummy = np.linspace(1,NBin, NBin)

UU     = UMin*beta**(1.*xdummy)

AA     = sigma**2*(alpha+2)/2./np.pi/(UMax**(alpha+2) - UMin**(alpha+2))
PU     = AA*UU**alpha

np.savetxt(OUTFILE, np.vstack((UU.T, PU.T)).T)
