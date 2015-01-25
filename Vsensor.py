# -*- coding: utf-8 -*-
"""
Created on Tue Jul 30 22:02:50 2013

@author: KyoungWon
""" 
from __future__ import division
import pylab as plt
import numpy as np
import main

length_sweep = np.linspace(2e-9,6e-9, 5) # [nm]
Vm_sweep = np.linspace(-0.2, 0.2, 21)  # [V]
dE = np.zeros(shape=(Vm_sweep.size, length_sweep.size, ), dtype=complex)
taum = np.zeros(shape=(Vm_sweep.size, length_sweep.size), dtype=complex)
overlap = np.zeros(shape=(Vm_sweep.size, length_sweep.size), dtype=complex)
kAmat = np.zeros(shape=(Vm_sweep.size, length_sweep.size), dtype=complex)

#output_hwf2=np.zeros(shape=(501,length_sweep.size), dtype=complex)
#output_dx2_int=np.zeros(shape=(801,sweep.size), dtype=complex)
#material = 'CdSe'
material = 'ZnSe_CdS'

for j in range(length_sweep.size):
    for i in range(Vm_sweep.size):        
        delta_E, tau, overlap_integral, kA = main.QCSE(material, length_sweep[j], Vm_sweep[i])
        dE[i,j] = delta_E
        taum[i,j] = tau   # radiative lifetime
        overlap[i,j] = overlap_integral  
        kAmat[i,j] = kA   # Auger rate
        #output[i,3]=overlap_integral2;
        #aa=abs(psi_h2_norm)
        #output_hwf2[:,i]=aa#[:,0]
        #output_dx2_int[:,i]=dx2_int
        print("   ")
        print("length is %.2f nm and Vm is %d mV") % (length_sweep[j]*1e9, Vm_sweep[i]*1e3)
        print("   ")

#plt.semilogy(sweep,output[:,1])
dE = np.real(dE)
kAmat = np.real(kAmat)
overlap = np.real(overlap)
taum = np.real(taum)