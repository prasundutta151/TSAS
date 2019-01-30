import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import sys
import pyfits
import os.path
import scipy.interpolate as intrp

EXP_NARG=9

if(len(sys.argv)!=EXP_NARG):
    print
    print "USAGE: python/run", sys.argv[0], "<OBSFILE> <MODFILE> <UMIN> <UMAX> <TAU0> <SIGMA> <ALPHA> <CHIFILE>"
    print "OBSFILE : Observed Line Channel Power SPectra"
    print "MODFILE : Model Line Channel Mower Spectra"
    print "Umin    : Value of Umin (Klambda)"
    print "Umax    : Value of Umax (Klambda)"
    print "Tau0    : Value of Tau0"
    print "sigma   : Value of sigma"
    print "alpha   : Value of alpha"
    print "CHIFILE : File to write ChiSq/DF values"
    
    sys.exit()

OBSFILE  = sys.argv[1]
MODFILE  = sys.argv[2]
UMin     = np.float(sys.argv[3])
UMax     = np.float(sys.argv[4])
tau0     = np.float(sys.argv[5])
sigma    = np.float(sys.argv[6])
alpha    = np.float(sys.argv[7])
CHIFILE  = sys.argv[8]

data = np.loadtxt(OBSFILE, usecols=(0,1,3))
UO   = data.T[0]
PO   = data.T[1]
EO   = data.T[2]

data = np.loadtxt(MODFILE)
UM   = data.T[0]
PM   = data.T[1]

if(UMin < np.min(UM)):
    print "Model does not have required UMin"
    sys.exit(1)

if(UMax > np.max(UM)):
    print "Model does not have required UMax"
    sys.exit(1)

PMIF = intrp.interp1d(UM,PM)

index = np.where(UO>UMin)
UO = UO[index]
PO = PO[index]

index = np.where(UO<UMax)
UO = UO[index]
PO = PO[index]

index = np.where(PO>0.)
UO = UO[index]
PO = PO[index]

PMI = PMIF(UO)

DOF = np.size(UO) - 3

MATRIX = (PMI-PO)/PO
CHIVAL = np.sum(MATRIX**2)/DOF

data = np.vstack((tau0, sigma, alpha, CHIVAL))
f_handle = file(CHIFILE, 'a')
np.savetxt(f_handle, data.T,  fmt = "%.2f")
