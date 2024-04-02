import numpy as np
# cgs
e = 4.80320425e-10
c = 2.99792458e10
me = 9.1093897e-28
mp = 1.6726231e-24
mec2 = me * c**2
mpc2 = mp * c**2
h = 6.62e-27
hbar = h/2/np.pi
lamCe = h/me/c
lambarCe = hbar/me/c
alphaF = 1/137.036
sigma_T = 6.65245854533e-25
G = 6.67259e-8
alphaG = G*mp**2/hbar/c
kB = 1.380658e-16
sigmaSB = 5.67051e-5
aSB = 7.5646e-15

Bc = me**2 * c**3 / e / hbar
Bcp = mp**2 * c**3 / e / hbar

eV2erg = 1.60218e-12
erg2eV = 1/eV2erg

AU = 1.496e13
pc = 3.086e18
ly = 9.463e17
msun = 1.989e33  # g
Nsun = msun/mp
Rsun = 6.96e10
Lsun = 3.9e33
Tsun = 5780

mEarth = 5.976e27
rEarth = 6.378e8
