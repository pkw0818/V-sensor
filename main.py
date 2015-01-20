# -*- coding: utf-8 -*-
"""
Created on Sun Jan 18 22:42:59 2015

@author: Philip
"""

from __future__ import division
import pylab as plt
import numpy as np
from sub import parameter, geometry, matrix, eigen

material = 'CdSe'
#material = 'ZnSe_CdS'
homospace = True
buffer = True # for heterostructure

"""
Define Parameter
"""
hbar = 6.626e-34/2/np.pi
e = 1.6e-19
eo = 8.85e-12;
m = 9.1e-31;
delta = 1e-15;              # use for avoid divergence in V(q)

if material == 'CdSe':
    param = parameter.CdSe()
    convgn = 'slowconv'
elif material == 'ZnSe_CdS':
    param = parameter.ZnSe_CdS()
    convgn = 'fastconv'
    
Eg = param[3]
k = param[5]

"""
Define Geometry 
"""
length = 2e-9
Vm_sweep = 0.2
Qw = length
L = Qw + 2e-9 
msize = 0.5e-10
n = np.linspace(-L, L, 2*L/msize + 1).size

print "Extra space is 2 nm and mesh size is %.1f A" % (msize * 1e10)

"""
Construct 1D Potential Energy
"""
cb1d, vb1d, er1d, me1d, mh1d = geometry.geo1d(material, L, Qw, msize)

Vm = Vm_sweep*L/Qw
Vm1d = np.linspace(0,Vm,n)
if Vm_sweep == 0:
    Vm1d = np.zeros(n)
Vm1d = Vm1d - (Vm/2)
CbEg = cb1d + Eg/2
VbEg = vb1d + Eg/2
CbEgVm = cb1d + Eg/2 + Vm1d
VbEgVm = vb1d + Eg/2 - Vm1d


"""
Construct Hamiltonian & Solve initial Hamiltonian
"""
# construct 2D + cylindrical matrix
KE, Apoisson = matrix.kinetic(material, me1d, mh1d, er1d, n, msize, hbar, m, e, eo, k)
PE = matrix.potential(n, CbEg, VbEg)

Hamilton = KE + PE # Vm is not applied yet. 
ev, ef = np.linalg.eigh(Hamilton)  # Eigen solver   
eE, hE, ewf_addr, hwf_addr = eigen.Energy(ev, cb1d, vb1d)     # Get eigen value 
psi_e, psi_h, psi_esq, psi_hsq = eigen.wfnormal(ef, n, msize) # normalization

plt.figure()
plt.plot(psi_esq)
plt.plot(psi_hsq)

distance = matrix.invdist(n)   # to get 1/|r1-r2|

"""
self-consistent iteration
"""
for i in range(0):
    Ve, Vh = matrix.Coulomb(material, n, distance, psi_esq, psi_hsq, e, eo, er1d)
    
    if convgn == 'slowconv':   # increase field 10% each gradually
        if i < 10:
            slowV = Vm1d*(i+1)/10            
            CbEgVm = CbEg + slowV
            VbEgVm = VbEg - slowV
    CbVm = np.diag(CbEgVm)  #Conduction band matrix
    VbVm = np.diag(VbEgVm)  # Valenceband matrix
    
    PE = matrix.potential(n, CbVm, VbVm)
    Vcoul = matrix.potential(n, -Vh, -Ve) 
    Hamilton = KE + PE + Vcoul
    ev, ef = np.linalg.eigh(Hamilton)  # Eigen solver   
    
    Ve1d = np.diag(Ve)
    Vh1d = np.diag(Vh)
    cbcorrect = CbEgVm - Vh1d
    vbcorrect = VbEgVm - Ve1d
    
    eE, hE, ewf_addr, hwf_addr = eigen.Energy(ev, cbcorrect, vbcorrect)     # Get eigen value 
    psi_e, psi_h, psi_esq, psi_hsq = eigen.wfnormal(ef, n, msize, ewf_addr, hwf_addr) # normalization

    
    plt.figure()
    plt.plot(psi_esq)
    plt.plot(psi_hsq)



