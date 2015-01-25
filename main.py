# -*- coding: utf-8 -*-
"""
Created on Sun Jan 18 22:42:59 2015

@author: Philip
"""

from __future__ import division
import pylab as plt
import numpy as np
#from scipy.sparse import linalg, coo_matrix

from sub import parameter, geometry, matrix, eigen

def QCSE(material, length, Vm_sweep):
    #material = 'CdSe'
    #material = 'ZnSe_CdS'
    #homospace = True
    #buffer = True # for heterostructure
    
    """ Define Parameter """
    
    hbar = 6.626e-34/2/np.pi
    e = 1.6e-19
    eo = 8.85e-12;
    m = 9.1e-31;
    
    if material == 'CdSe':
        param = parameter.CdSe()
        convgn = 'slowconv'
    elif material == 'ZnSe_CdS':
        param = parameter.ZnSe_CdS()
        convgn = 'fastconv'
        
    mh1 = param[1]
    er1 = param[2]
    Eg = param[3]
    Ep = param[4]
    k = param[5]
    hole_barrier = param[7]
    
    """ Define Geometry  """
    
    #length = 2e-9
    #Vm_sweep = 0.2
    Qw = length + 0.05e-9
    L = Qw + 2e-9 
    msize = 0.5e-10
    X = np.linspace(-L , L, 2*L/msize +1)
    n = X.size
    
    print "Extra space is 2 nm and mesh size is %.1f A" % (msize * 1e10)
    
    """ Construct 1D Potential Energy """
    
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
    
    
    """ Construct Hamiltonian & Solve initial Hamiltonian """
    
    # construct 2D + cylindrical matrix
    KE, Apoisson = matrix.kinetic(material, me1d, mh1d, er1d, n, msize, hbar, m, e, eo, k)
    PE = matrix.potential(n, CbEg, VbEg)
    
    Hamilton = KE + PE # Vm is not applied yet. 
    #sHamilton = coo_matrix(Hamilton)  # sparse matrix
    ev, ef = np.linalg.eigh(Hamilton)  # Eigen solver   
    
    eE, hE, ewf_addr, hwf_addr = eigen.Energy(n, ev, CbEg, VbEg)     # Get eigen value 
    psi_e, psi_h, psi_esq, psi_hsq = eigen.wfnormal(ef, n, msize, ewf_addr, hwf_addr) # normalization
    
    print('1st iteration.  Electron Energy = %.3f') % eE
    
    #plt.figure()
    #plt.plot(psi_esq)
    #plt.plot(psi_hsq)
    
    distance = matrix.invdist(n)   # to get 1/|r1-r2|
    
    """ self-consistent iteration """
    
    bfeE = 0
    for i in range(20):
        Ve, Vh = matrix.Coulomb(material, n, distance, psi_esq, psi_hsq, e, eo, er1d)
        
        if convgn == 'slowconv':   # increase field 10% each gradually
            if i < 10:
                slowV = Vm1d*(i+1)/10            
                CbEgVm = CbEg + slowV
                VbEgVm = VbEg - slowV
        
        PE = matrix.potential(n, CbEgVm, VbEgVm)
        Vcoul = matrix.potential(n, -Vh, -Ve) 
        Hamilton = KE + PE + Vcoul
        ev, ef = np.linalg.eigh(Hamilton)  # Eigen solver   
        
        cbcorrect = CbEgVm - Vh
        vbcorrect = VbEgVm - Ve
        
        eE, hE, ewf_addr, hwf_addr = eigen.Energy(n, ev, cbcorrect, vbcorrect)  # Get eigen value 
        psi_e, psi_h, psi_esq, psi_hsq = \
        eigen.wfnormal(ef, n, msize, ewf_addr, hwf_addr) # normalization
    
        #plt.figure()
        #plt.plot(psi_esq)
        #plt.plot(psi_hsq)
        print('%d st iteration.  Electron Energy = %.3f') % (i+2,  eE)
        
        if convgn == 'slowconv' and i >10:
            if abs(float(eE - bfeE)) < 0.001:
                print('Congratulation. Converged!')
                break
        elif convgn == 'fastconv':
            if abs(float(eE - bfeE)) < 0.001:
                print('Congratulation. Converged!')
                break
        
        bfeE = eE
    
    """ Get dE, overlap integral and radiative lifetime"""
    
    cc_e = np.conjugate(psi_e)
    cc_h = np.conjugate(psi_h)
    
    delta_E = eE - hE  # dE
    e_h_multiple = cc_e *psi_h
    spatial_sum = sum(e_h_multiple) *msize
    overlap_integral = abs(spatial_sum)**2
    tau = 2*np.pi*eo*m*(3e8**3)*(hbar**2)/np.sqrt(2.5)/(e**2)/(Eg+delta_E)/Ep/e/e/overlap_integral
    
    """ Calculation Auger rate """
    eq_x1 = cc_e * cc_h
    Eex = Eg- hole_barrier + (eE-Eg/2) + 2*abs(hE+Eg/2)  # excited energy of hole
    kf = np.sqrt(2*m*mh1*Eex*e)/ hbar  # [1/m]
    phi_F = 1/np.sqrt(2*L)*(np.e)**(-1j*kf*X)       #[1/sqrt(m)]
    
    Vdistance = distance / msize     #[1/m]
    
    " First integral over x2"
    dx2_int = np.zeros(n, dtype=complex)
    temp = np.diag(phi_F*cc_h)                       #[1/m]
    for x1 in range(n):
        temp2 = Vdistance[x1,:]                #[1/m]
        dx2_int[x1] = sum(np.dot(temp,temp2))      #[1m^2]   
    
    " Second integral over x1"
    Mif = sum(eq_x1*dx2_int)*msize**(2)*(e**2)*np.sqrt(2)/4/np.pi/eo/er1        #[1/m * C^2 / C^2 * J*m= J]
    dos_Ef = np.sqrt(2)*(m*mh1)**(1.5)/np.pi/np.pi/hbar/hbar/hbar*np.sqrt(Eex*e)*2*Qw*e  #3D DOS * Qw * 2
    kA = 2*np.pi/hbar*abs(Mif)**(2)*dos_Ef              #[s/J * (J)^2 *  * J^-1]
    Aug_life = 1/kA
    #Eex=kf**(2)*hbar**(2)/2/m/mh1/e              #[eV]
    #dos_Ef=1/pi/hbar*sqrt(m*mh1/2/Eex/e)*2*Qw           #[1/J)]
    #dos_Ef=m*mh1/pi/hbar/hbar*e   
    return delta_E, tau, overlap_integral, kA