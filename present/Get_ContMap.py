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
EXP_NARG=8
decC = 2.201447303442E+01

if(len(sys.argv)!=EXP_NARG):
    print
    print "USAGE: python/run", sys.argv[0], "<CONTFILE> <DELIFILE> <REGFILE> <OUTPFILE> <SMSIG> <NAVCHAN>"
    print "CONTFILE: Continium FITS file"
    print "DELIFILE: FITS File with Del I(\theta)"
    print "REGIFILE: Definition of region to exclude (only ellipse)"
    print "OUTPFILE: Output FITS file"
    print "SMSIG   : Smoothing length sclaes (units of pixel)"
    print "NAVCHAN : No of channels taken for the continium"
    print "SCALEW  : Scaling of window"
    sys.exit()

CONTFILE = sys.argv[1]
DELIFILE = sys.argv[2]
REGIFILE = sys.argv[3]
OUTPFILE = sys.argv[4]
SSIG     = np.float(sys.argv[5])
NAVCHAN  = np.int(sys.argv[6])
SCALEW   = np.float(sys.argv[7])

if(os.path.exists(OUTPFILE)):
    print "\nWARNING: File %s Exists, not created again\n" %OUTPFILE
    sys.exit(0)
if(os.path.exists(REGIFILE)==False):
    print "\nWARNING: File %s does not Exists\n" %REGIFILE
    sys.exit(0)
if(os.path.exists(DELIFILE)==False):
    print "\nWARNING: File %s does not Exists\n" %DELIFILE
    sys.exit(0)
if(os.path.exists(CONTFILE)==False):
    print "\nWARNING: File %s does not Exists\n" %CONTFILE
    sys.exit(0)

dataC, NGridC, dxyC, BMajC, BMinC = FITS_READ(CONTFILE)
dataI, NGridI, dxyI, BMajI, BMinI = FITS_READ(DELIFILE)

ellparm = np.loadtxt(REGIFILE)
#ellparm = np.array(ellparm[8:-1].split(',')).astype(int)
ex0 = ellparm.T[0]
ey0 = ellparm.T[1]
exL = ellparm.T[2]/1.25
eyL = ellparm.T[3]/1.25
eth = (np.pi/2.-ellparm.T[4])*DTR


xy = np.linspace(-ex0, -ex0+NGridC, NGridC).astype(int)
xx,yy = np.meshgrid(xy,xy)
xp = (xx*np.cos(eth) - yy*np.sin(eth))/exL
yp = (xx*np.sin(eth) + yy*np.cos(eth))/eyL

rr  = np.hypot(xx,yy)**2
rel = np.hypot(xp,yp)**2

xx = xx + ex0
yy = yy + ey0

Omean = np.mean(dataC[np.where(rel>1.)])
Ostdv = np.std(dataC[np.where(rel>1.)])
dataC = dataC - Omean

# TESTING ONLY
'''
dataO = dataC
dxyO = dxyC
BMajO = BMajC
BMinO = BMinC
FITS_WRITE("II.fits", dataO, dxyO, BMajO, BMinO)
'''
# TESTING DONE

BeamArea     = np.pi*BMajC*BMinC/4./np.log(2.)
PixelPerBeamC = BeamArea/dxyC**2
dataC = dataC/NAVCHAN

dataC = gaussian_filter(input=dataC, sigma=SSIG,  mode='reflect', truncate=7.) 
I0val = np.sum(dataC[np.where(rel<1.)])/PixelPerBeamC

dataC = dataC/I0val
stdC  = np.std(dataC)*I0val

I0val = I0val*SCALEW

print "BeamArea (Cont): %.3e sq deg" %BeamArea
print "NPix/Beam(Cont): %.3e" %PixelPerBeamC
print "I0val      (Jy): %.3e" %I0val
print "stdC       (Jy): %.3e" %stdC

BeamArea     = np.pi*BMajI*BMinI/4./np.log(2.)
PixelPerBeam = BeamArea/dxyI**2
#dataI = dataI/PixelPerBeam

Norm  = np.sum(dataC)/PixelPerBeam/np.pi**2./4.
dataI = dataI/Norm
stdI  = np.std(dataI)

print "BeamArea (DelI): %.3e sq deg" %BeamArea
print "NPix/Beam(DelI): %.3e" %PixelPerBeam
print "stdI       (Jy): %.3e" %stdI

if (dxyC*NGridC < dxyI*NGridI):

    print "Choosing the C image"
    dxyO   = dxyC
    NGridO = NGridC
    BMajO  = BMajC
    BMinO  = BMinC
    zratio = dxyC/dxyI
    dataIo = intp.zoom(dataI, zratio, order=3, mode='constant', cval=0.0, prefilter=True)/zratio**2
    xmin   = dataIo.shape[0]/2 - dataC.shape[0]/2
    xmax   = dataIo.shape[0]/2 + dataC.shape[0]/2
    ymin   = dataIo.shape[1]/2 - dataC.shape[1]/2
    ymax   = dataIo.shape[1]/2 + dataC.shape[1]/2
    dataI  = dataIo[xmin:xmax, ymin:ymax]

else:
    print "Choosing the I image"
    dxyO   = dxyI
    NGridO = NGridI
    BMajO  = BMajI
    BMinO  = BMinI
    zratio = dxyI/dxyC
    dataCo = intp.zoom(dataC, zratio, order=3, mode='constant', cval=0.0, prefilter=True)/zratio**2
    xmin   = dataCo.shape[0]/2 - dataI.shape[0]/2
    xmax   = dataCo.shape[0]/2 + dataI.shape[0]/2
    ymin   = dataCo.shape[1]/2 - dataI.shape[1]/2
    ymax   = dataCo.shape[1]/2 + dataI.shape[1]/2
    dataC  = dataCo[xmin:xmax, ymin:ymax]

print np.std(dataI), np.sum(dataI)
print np.std(dataC), np.sum(dataC)

dataO = dataC*(I0val + dataI)
FITS_WRITE("contmod.fits", dataO, dxyO, BMajO, BMinO)

dataO = I0val*dataC
FITS_WRITE("II.fits", dataO.T, dxyO, BMajO, BMinO)

dataO = dataC 
FITS_WRITE("Win.fits", dataO, dxyO, BMajO, BMinO)

dataO = dataI*dataC
FITS_WRITE("delI.fits", dataO, dxyO, BMajO, BMinO)

