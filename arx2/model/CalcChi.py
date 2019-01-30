import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import sys
import pyfits
import os.path
import scipy.interpolate as intrp
import re 

EXP_NARG=4

if(len(sys.argv)!=EXP_NARG):
    print
    print "USAGE: python/run", sys.argv[0], "<ALPHA>  <SIGMA> <TAU0>"
    print "alpha   : Value of alpha"
    print "sigma   : Value of sigma"
    print "Tau0    : Value of Tau0"
    sys.exit()

alpha    = np.float(sys.argv[1])
sigma    = np.float(sys.argv[2])
tau0     = np.float(sys.argv[3])


# Reading inputs from the inputfile
inputfile="tau.input"
fp = open(inputfile, mode='r')

tempr = fp.readline()
tempr = fp.readline()
OBSFILE = fp.readline().split(':')[0].split('\t')[0]
tempr = fp.readline()
WORKDIR = fp.readline().split(':')[0].split('\t')[0]
CHIFILE = fp.readline().split(':')[0].split('\t')[0]
tempr = fp.readline().split(':')[0]
tempr = re.findall(r"(\d+\.\d+)", tempr)
UMin = np.float(tempr[0])
UMax = np.float(tempr[1])
NBin = np.int(fp.readline().split(':')[0].split('\t')[1])
tempr = fp.readline()
tempr = fp.readline()
tempr = fp.readline()
NReal = np.int(fp.readline().split(':')[0].split('\t')[0])

real = 0
filename = WORKDIR + '/REAL/LIN_' + "%d" %(real) + ".POW"
data = np.loadtxt(filename, usecols=(0,1))

NBin = np.size(data.T[0])
UM = np.zeros(NBin)
PM = np.zeros(NBin)
EM = np.zeros(NBin)

for real in range(NReal):

    filename = WORKDIR + '/REAL/LIN_' + "%d" %(real) + ".POW"
    
    data = np.loadtxt(filename, usecols=(0,1))

    UM = UM + data.T[0]
    PM = PM + data.T[1]
    EM = EM + data.T[1]**2

UM = UM/(1.*NReal)
PM = PM/(1.*NReal)
EM = EM/(1.*NReal)
EM = np.sqrt(EM - PM**2)

data = np.loadtxt(OBSFILE, usecols=(0,1,3))
UO   = data.T[0]
PO   = data.T[1]
EO   = data.T[2]

'''
if(UMin < np.min(UM)):
    print "Model does not have required UMin"
    sys.exit(1)

if(UMax > np.max(UM)):
    print "Model does not have required UMax"
    sys.exit(1)
'''
PMIF = intrp.interp1d(UM,PM)
EMIF = intrp.interp1d(UM,EM)

index = np.where(UO>UMin)
UO = UO[index]
PO = PO[index]
EO = EO[index]

index = np.where(UO<UMax)
UO = UO[index]
PO = PO[index]
EO = EO[index]

index = np.where(PO>0.)
UO = UO[index]
PO = PO[index]

PMI = PMIF(UO)
EMI = EMIF(UO)
DOF = np.size(UO) - 3

MATRIX = (PMI-PO)/EMI
CHIVAL = np.sum(MATRIX**2)/DOF

'''
MATRIX1 = (PMI-PO)/EO
CHIVAL1 = np.sum(MATRIX**2)/DOF
'''

data = np.vstack((tau0, sigma, alpha, CHIVAL))
f_handle = file(CHIFILE, 'a')
np.savetxt(f_handle, data.T,  fmt = "%.2f")
f_handle.close()

PSFILE = WORKDIR + '/REAL/LIN_' + "%.2f_%.2f_%.2f" %(tau0, sigma, alpha) + ".POW"
data = np.vstack((UM, PM, EM))
np.savetxt(PSFILE, data.T)
