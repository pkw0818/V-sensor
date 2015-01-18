# -*- coding: utf-8 -*-
"""
Created on Sat Mar 22 12:12:49 2014

@author: Philip
"""

# -*- coding: utf-8 -*-
"""
Created on Wed Jan 29 05:49:52 2014

@author: KyoungWon
"""


# -*- coding: utf-8 -*-
"""
Created on Tue Jul 30 22:02:50 2013

@author: KyoungWon
""" 
ip = get_ipython()
ip.magic("reset -f")
from pylab import * 

  
length_sweep=np.linspace(4.5e-9,10e-9,56) # [nm], half of NR length
#output=zeros(shape=(Vm_sweep.size,3), dtype=complex)
#output_hwf2=np.zeros(shape=(501,length_sweep.size), dtype=complex)
#output_dx2_int=np.zeros(shape=(801,sweep.size), dtype=complex)
energy_data=zeros(length_sweep.size, dtype=complex)
kA_data=zeros(length_sweep.size, dtype=complex)
overlap_data=zeros(length_sweep.size, dtype=complex)
tau_data=zeros(length_sweep.size, dtype=complex)

#import CdSe
#import sc_CdSe
import sc_ZnSe_CdS
#import ZnSe_CdS
#import TwocolorQD2

for i in range(length_sweep.size):
    Vm1=-0.08/(length_sweep[i]*2)   # electric field when Vm=-80mV
    Vm2=0.08/(length_sweep[i]*2)
    Vm_sweep=0
        #Vm_sweep=np.linspace(Vm1,Vm2,17)  #17 is number of Vm sweep from -80mV to 80mV
    #delta_E, kA, overlap_integral, tau =sc_CdSe.sc_CdSe(length_sweep[i],Vm_sweep)        
    delta_E, kA, overlap_integral, tau =sc_ZnSe_CdS.sc_ZnSe_CdS(length_sweep[i],Vm_sweep)        
    #delta_E, kA, overlap_integral =CdSe.CdSe(length_sweep[i],Vm_sweep[j])
    energy_data[i]=delta_E
    kA_data[i]=kA
    overlap_data[i]=overlap_integral
    tau_data[i]=tau
    #delta_E, kA, overlap_integral =ZnSe_CdS.ZnSe_CdS(length_sweep,Vm_sweep[i]);
    #delta_E, delta_E2, overlap_integral, overlap_integral2, psi_h2_norm =TwocolorQD2.TwocolorQD2(sweep[i]);
    #output[i,3]=overlap_integral2;
    #aa=abs(psi_h2_norm)
    #output_hwf2[:,i]=aa#[:,0]
    #output_dx2_int[:,i]=dx2_int
    print i

#plt.semilogy(sweep,output[:,1])
energy_data=real(energy_data)
kA_data=real(kA_data)
overlap_data=real(overlap_data)
tau_data=real(tau_data)
    