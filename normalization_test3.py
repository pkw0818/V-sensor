# -*- coding: utf-8 -*-
"""
Created on Wed Aug 21 15:55:38 2013

@author: KyoungWon
"""
aa=phi_kf_expanding*conjugate(phi_kf_expanding)
phi_kf2_rdr=zeros((r_limit,n-2*z_limit,theta.size), dtype=complex)

for i in range(r_limit):
    phi_kf2_rdr[i,:,:]=aa[i,:,:]*r[i]

norm_phi=sum(phi_kf2_rdr)*dr*2*d_theta*dz



bb=cc_hwf_expanding*conjugate(cc_hwf_expanding)
cc_hwf_expanding_rdr=zeros((r_limit,n-2*z_limit,theta.size), dtype=complex)

for i in range(r_limit):
    cc_hwf_expanding_rdr[i,:,:]=bb[i,:,:]*r[i+1]

norm_cc=sum(cc_hwf_expanding_rdr)*dr*2*d_theta*dz