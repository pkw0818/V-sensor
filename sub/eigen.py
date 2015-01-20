# -*- coding: utf-8 -*-
"""
Created on Mon Jan 19 11:57:03 2015

@author: Philip
"""
from __future__ import division
import numpy as np

def wfnormal(ef, n, msize, *args):
    
    psi_e = ef[0:n, n]  
    psi_h = ef[n:2*n, n-1]

    if args is not None:
        psi_e = ef[0:n, args[0]]  
        psi_h = ef[n:2*n, args[0]]
        
    psi_e_sq = psi_e*np.conjugate(psi_e)   
    psi_h_sq = psi_h*np.conjugate(psi_h)
    norm_e = sum(psi_e_sq) *msize
    norm_h = sum(psi_h_sq) *msize
    psi_e_sq_norm = psi_e_sq /msize
    psi_h_sq_norm = psi_h_sq /msize
    psi_e_norm = psi_e /np.sqrt(norm_e)   #normalized wf
    psi_h_norm = psi_h /np.sqrt(norm_h)
    
    return psi_e_norm, psi_h_norm, psi_e_sq_norm, psi_h_sq_norm
    
def Energy(ev, cb, vb):
    Cb_temp = abs(ev - min(cb))
    Vb_temp = abs(ev + min(vb))
    Cb_E_address = np.where(Cb_temp == min(abs(Cb_temp)))[0][0]
    Vb_E_address = np.where(Vb_temp == min(abs(Vb_temp)))[0][0]
    Cb_E = ev[Cb_E_address]
    Vb_E = ev[Vb_E_address]
    
    return Cb_E, Vb_E, Cb_E_address, Vb_E_address