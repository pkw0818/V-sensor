# -*- coding: utf-8 -*-
"""
Created on Mon Jan 19 00:13:52 2015

@author: Philip
"""

from __future__ import division
import numpy as np
import parameter

def kinetic(material, me, mh, er, n, msize, hbar, m, e, eo, k):    
    
    Ae = np.mat(np.zeros((n,n)))
    Ah = np.mat(np.zeros((n,n)))
    Apoisson = np.mat(np.zeros((n,n)))   
    off_diag = np.mat(np.zeros((n,n)))
    
    Ae[0,0] = -(1/me[0] + 1/me[1]);
    Ae[1,0] = 1/me[0];
    Ae[n-1,n-1] = -(1/me[n-2] + 1/me[n-1]);
    #Ae[n-1,n-2]=Ae[n-1,n-2]*2  % Dirichlet-Neumann boundary condition for 
    Ae[n-2,n-1] = 1/me[n-1];
    
    Ah[0,0] = -(1/mh[0] + 1/mh[1]);
    Ah[1,0] = 1/mh[0];
    Ah[n-1,n-1] = -(1/mh[n-2] + 1/mh[n-1]);
    #Ah[n-1,n-2]=Ah[n-1,n-2]*2  % Dirichlet-Neumann boundary condition for 
    Ah[n-2,n-1] = 1/mh[n-1];
    
    Apoisson[0,0] = -(er[0]+er[1]);
    Apoisson[1,0] = er[0];
    Apoisson[n-1,n-1] = -(er[n-2] + er[n-1]);
    #Apoisson[n-1,n-2]=Apoisson[n-1,n-2]*2  % Dirichlet-Neumann boundary condition for 
    Apoisson[n-2,n-1] = er[n-1];
    
    off_diag[0,1] = 1;
    off_diag[1,0] = -1;
    off_diag[n-1,n-2] = -1;
    off_diag[n-2,n-1] = 1;
    
    for i in range(1,n-1):
        Ae[i,i] = -(1/me[i-1] + 1/me[i+1])
        Ah[i,i] = -(1/mh[i-1] + 1/mh[i+1])
        Apoisson[i,i] = -(er[i-1] + er[i+1])
        Ae[i+1,i] = 1/me[i-1]
        Ah[i+1,i] = 1/mh[i-1]
        Apoisson[i+1,i] = er[i-1]
        Ae[i-1,i] = 1/me[i+1]
        Ah[i-1,i] = 1/mh[i+1]
        Apoisson[i-1,i] = er[i+1]
        off_diag[i,i+1] = 1
        off_diag[i+1,i] = -1
    
    KE = np.zeros(shape=(2*n,2*n), dtype=complex)
    KE[0:n, 0:n] = -Ae/ m/ e* (hbar**2)/ 2/ (msize**2) 
    KE[n:2*n, n:2*n] = Ah/ m/ e* (hbar**2)/ 2/ (msize**2) 
    KE[0:n, n:2*n] = off_diag *(-1j) *hbar *k / 2 /msize
    KE[n:2*n, 0:n] = KE[0:n, n:2*n]
    
    Apoisson = Apoisson* eo/ (msize**2);
    
    return KE, Apoisson
    
def potential(n, Cb, Vb):
    PE = np.zeros(shape=(2*n,2*n), dtype=complex)
    PE[0:n, 0:n] = np.diag(Cb)
    PE[n:2*n, n:2*n] = -np.diag(Vb)
    
    return PE

def invdist(n):
    distance_lower = np.zeros((n,n))
    for i in range(n):
        temp = np.eye(n, k = i)*i  # this k is not k dot p's k
        distance_lower = distance_lower + temp
    distance_upper = distance_lower.T
    distance = distance_lower + distance_upper
    distance = distance + np.eye(n, k = 0)
    distance = 1/distance - np.eye(n, k = 0)
    
    return distance
    
def Coulomb (material, n, distance, psi_esq, psi_hsq, e, eo, er1d):
    if material == 'CdSe':
        param = parameter.CdSe()
        er_n = param[2]
    elif material =='ZnSe_CdS':
        param = parameter.ZnSe_CdS()
        er_1 = param[2]
        er_2 = param[10]
        er_n = (er_1 + er_2)/2
            
    Velectron = np.zeros(n)
    Vhole = np.zeros(n)
    e_over_r = np.zeros(n)
    h_over_r = np.zeros(n)
    for j in range(n):
        e_over_r = np.multiply(distance[j,:], np.real(psi_esq))
        h_over_r = np.multiply(distance[j,:], np.real(psi_hsq))
        Velectron[j] = sum(e_over_r) / er1d[j]
        Vhole[j] = sum(h_over_r) / er1d[j]
    Ve = Velectron*e/4/np.pi/eo/er_n
    Vh = Vhole*e/4/np.pi/eo/er_n
    #Ve = np.diag(Velectron)
    #Vh = np.diag(Vhole)
    
    return Ve, Vh