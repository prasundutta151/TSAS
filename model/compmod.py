import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

data = np.loadtxt("CRAB_CONT350_64.POW")
U1 = data.T[0]
P1 = data.T[1]

data = np.loadtxt("contmod.pow")
U2 = data.T[0]
P2 = data.T[1]

plt.clf()
plt.xlim(0.1,40.)
plt.ylim(0.01,1.e6)
plt.xlabel(r"U [k $\lambda$]", fontsize=20)
plt.ylabel(r"P (U) [arbitrary]", fontsize=20)
plt.loglog(U2,P2, "-", color='black', linewidth=3, label='Model')
plt.loglog(U1,P1, "o", color='blue',  label='Observed')

plt.legend()
