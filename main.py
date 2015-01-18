# -*- coding: utf-8 -*-
"""
Created on Sun Jan 18 22:42:59 2015

@author: Philip
"""

from __future__ import division
import pylab as plt
import numpy as np
from sub import parameter, geometry, matrix

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
elif material == 'ZnSe_CdS':
    param = parameter.ZnSe_CdS()
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
KE, Apoisson = matrix.kinetic(material, me1d, mh1d, er1d, n, msize, hbar, m, e, eo, k)
PE = matrix.potential(n, CbEg, VbEg)

Hamilton = KE + PE # Vm is not applied yet. 
ev, ef = np.linalg.eigh(Hamilton)   
   

Cb_temp = abs(ev - min(cb1d))
Vb_temp = abs(ev + min(vb1d))
Cb_E_address = np.where(Cb_temp == min(abs(Cb_temp)))[0][0]
Vb_E_address = np.where(Vb_temp == min(abs(Vb_temp)))[0][0]
Cb_energy = ev[Cb_E_address]
Vb_energy = ev[Vb_E_address]


psi_e = ef[0:n, n]  
psi_h = ef[n:2*n, n-1]
psi_e_sq = psi_e*np.conjugate(psi_e)   
psi_h_sq = psi_h*np.conjugate(psi_h)
norm_e = sum(psi_e_sq) *msize
norm_h = sum(psi_h_sq) *msize
psi_e_sq_norm = psi_e_sq /msize
psi_h_sq_norm = psi_h_sq /msize
psi_e_norm = psi_e /np.sqrt(norm_e)   #normalized wf
psi_h_norm = psi_h /np.sqrt(norm_h)
