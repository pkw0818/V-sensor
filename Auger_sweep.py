
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 30 22:02:50 2013

@author: KyoungWon
""" 
from __future__ import division
#ip = get_ipython()
#ip.magic("reset -f")
from pylab import * 

length_sweep=np.linspace(4e-9,4e-9,1) # [nm]
Vm_sweep=np.linspace(-0.2, 0.2, 3)  # [V]
dE=zeros(shape=(Vm_sweep.size, length_sweep.size, ), dtype=complex)
kAmat=zeros(shape=(Vm_sweep.size, length_sweep.size), dtype=complex)
overlap=zeros(shape=(Vm_sweep.size, length_sweep.size), dtype=complex)

#output_hwf2=np.zeros(shape=(501,length_sweep.size), dtype=complex)
#output_dx2_int=np.zeros(shape=(801,sweep.size), dtype=complex)

#import CdSe
#import sc_CdSe
import sc_ZnSe_CdS
#import ZnSe_CdS
#import TwocolorQD2
for j in range(length_sweep.size):
    for i in range(Vm_sweep.size):
        #delta_E, kA, overlap_integral =sc_CdSe.sc_CdSe(length_sweep[j], Vm_sweep[i]);
        delta_E, kA, overlap_integral =sc_ZnSe_CdS.sc_ZnSe_CdS(length_sweep,Vm_sweep[i]);
        #delta_E, delta_E2, overlap_integral, overlap_integral2, psi_h2_norm =TwocolorQD2.TwocolorQD2(sweep[i]);
        dE[i,j]=delta_E
        kAmat[i,j]=kA
        overlap[i,j] = overlap_integral;
        #output[i,3]=overlap_integral2;
        #aa=abs(psi_h2_norm)
        #output_hwf2[:,i]=aa#[:,0]
        #output_dx2_int[:,i]=dx2_int
        print j,i

#plt.semilogy(sweep,output[:,1])
dE = real(dE)
kAmat = real(kAmat)
overlap = real(overlap)