# -*- coding: utf-8 -*-
"""
Created on Sat Sep 07 15:08:57 2013

@author: KyoungWon
"""
ev_num=479
if ev_num >= size(ev)/2:
    psi=ef[0:n, ev_num]  
else:
    psi=ef[n:2*n, ev_num]
        
psi_sq=psi*conjugate(psi)   
norm_e=sum(psi_sq)*space
psi_norm=psi_sq/space
psi_sq_norm=psi_sq/space
psi_norm=psi/sqrt(norm_e)  
plot(psi_norm)