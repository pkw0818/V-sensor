# -*- coding: utf-8 -*-
"""
Created on Mon Aug 19 16:48:15 2013

@author: KyoungWon
"""
eef2=eef_n*conjugate(eef_n)
eef2_rdr=zeros(((m-1)*n))
for i in range(m-1):
    eef2_rdr[i*n:(i+1)*n]=eef2[i*n:(i+1)*n]*r[i+1]
norm_e=sum(eef2_rdr)*dr*dz*2*pi