import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import sys
import pyfits
import os.path
import scipy.interpolate as intrp
import re 
EXP_NARG=2

if(len(sys.argv)!=EXP_NARG):
    print
    print "USAGE: python/run", sys.argv[0], "<ChiFile>"
    print "ChiFile: File with Chisq data"
    sys.exit()

ChiFile = sys.argv[1]

data = np.loadtxt(ChiFile)

tau0 = data.T[0]
sigm = data.T[1]
alph = data.T[2]
Chiv = data.T[3]

min_chi_index = np.where(Chiv == np.min(Chiv))

tau0_min = tau0[min_chi_index]

index_2d = np.where(tau0==tau0_min)
Chiv2D  = Chiv[index_2d]
sigm2D  = sigm[index_2d]
alph2D  = alph[index_2d]

NG = np.size(Chiv2D)
NG = np.int(np.sqrt(NG))

Chiv2D = np.reshape(Chiv2D, (NG, NG))
sigm2D = sigm2D[:NG]
alph2D = np.reshape(alph2D, (NG, NG)).T[0][:NG]

dx = alph2D[1] - alph2D[0]
dy = sigm2D[1] - sigm2D[0]

plt.clf()
plt.figure(1)
plt.xlim(alph2D[0]+dx/2., alph2D[NG-1]+dx/2.)
plt.ylim(sigm2D[0]+dy/2., sigm2D[NG-1]+dy/2.)
plt.pcolor(alph2D+dx/2, sigm2D+dy/2, np.log(Chiv2D))
plt.colorbar()

sig0_min = sig0[min_chi_index]
index_2d = np.where(sig0==sig0_min)
Chiv2D  = Chiv[index_2d]
sigm2D  = sigm[index_2d]
alph2D  = alph[index_2d]

NG = np.size(Chiv2D)
NG = np.int(np.sqrt(NG))

Chiv2D = np.reshape(Chiv2D, (NG, NG))
sigm2D = sigm2D[:NG]
alph2D = np.reshape(alph2D, (NG, NG)).T[0][:NG]

dx = alph2D[1] - alph2D[0]
dy = sigm2D[1] - sigm2D[0]

plt.clf()
plt.xlim(alph2D[0]+dx/2., alph2D[NG-1]+dx/2.)
plt.ylim(sigm2D[0]+dy/2., sigm2D[NG-1]+dy/2.)
plt.pcolor(alph2D+dx/2, sigm2D+dy/2, np.log(Chiv2D))
plt.colorbar()
