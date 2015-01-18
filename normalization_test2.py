# -*- coding: utf-8 -*-
"""
Created on Mon Aug 19 23:21:19 2013

@author: KyoungWon
"""
phi_kf2=phi_kf_n*conjugate(phi_kf_n)
phi_kf2_rdr=zeros((m,theta.size), dtype=complex)
for i in range(m):
    phi_kf2_rdr[i,:]=phi_kf2[i,:]*r[i]
norm_phi=sum(phi_kf2_rdr)*dr*2*d_theta*z[-1]