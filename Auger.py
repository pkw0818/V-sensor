# -*- coding: utf-8 -*-
"""
Spyder Editor

This temporary script file is located here:
C:\Users\KyoungWon\.spyder2\.temp.py
"""

ip = get_ipython()
ip.magic("reset -f")
from pylab import *   #include numpy, plotting function
from scipy.interpolate import interp1d
from scipy.special import expi, exp1

# Define Constant 
hbar=6.626e-34/2/math.pi
e=1.6e-19
eo=8.85e-12;
m=9.1e-31;
me1=0.04;   # CdSe me=0.13 mh=0.45       # CdS me:0.21 mh:0.8
mh1=1.0;    # CdTe me=0.1  mh=0.4
#me2=0.04*9.1e-31;   
#mh2=1*9.1e-31;
er1=6.0;           # CdSe: 9.3(WZ) 10.2(ZB)        # CdTe: 10.2(ZB)
er2=12.0;           # CdS: 8.9 or 7 or 6.3
Eg=1.0;
Cbo=0.3;           # CdSe qX: 4.6     Eg: 2.1
Vbo=0.3;             # CdTe qX: 4.1     Eg: 1.8
k=math.sqrt(21./2);           # K dot P matrix element
delta=1e-11;              # use for avoid divergence in V(q)

#Define Geometry
L=30e-9;                # Sim. length L >> QW (confined area)
Qw = 2e-9;
space=0.5e-10;
X = np.linspace(-L , L, 2*L/space+1)
n= X.size
region=[-L,-Qw,Qw,L]
indx=np.zeros(len(region))
for var in range(len(region)):
    indx[var]=(region[var]-region[0])/space

CB_value=[Cbo,0,Cbo]
VB_value=[Vbo,0,Vbo]
er_value=[er1,er1,er1]
me_value=[me1,me1,me1]
mh_value=[mh1,mh1,mh1]

import piecewise2      #make a piecewise function with x as region, y as value
Cb=piecewise2.piecewise2(indx,CB_value)
Vb=piecewise2.piecewise2(indx,VB_value)
er=piecewise2.piecewise2(indx,er_value)
me=piecewise2.piecewise2(indx,me_value)
mh=piecewise2.piecewise2(indx,mh_value)

Field=0
PE=np.linspace(0,Field,n)
if Field==0:
    PE=np.zeros(n)

PE=PE-(Field/2)
Cb=Cb+PE+Eg/2
Vb=Vb-PE+Eg/2

#Define Hamiltonian cell
Ae=mat(zeros((n,n)))
Ah=mat(zeros((n,n)))
Apoisson=mat(zeros((n,n)))
Pe=diag(Cb)
Ph=diag(Vb)
off_diag=mat(zeros((n,n)))

Ae[0,0]=-(1/me[0]+1/me[1]);
Ae[1,0]=1/me[0];
Ae[n-1,n-1]=-(1/me[n-2]+1/me[n-1]);
#Ae[n-1,n-2]=Ae[n-1,n-2]*2  % Dirichlet-Neumann boundary condition for 
Ae[n-2,n-1]=1/me[n-1];

Ah[0,0]=-(1/mh[0]+1/mh[1]);
Ah[1,0]=1/mh[0];
Ah[n-1,n-1]=-(1/mh[n-2]+1/mh[n-1]);
#Ah[n-1,n-2]=Ah[n-1,n-2]*2  % Dirichlet-Neumann boundary condition for 
Ah[n-2,n-1]=1/mh[n-1];

Apoisson[0,0]=-(eo*er[0]+eo*er[1]);
Apoisson[1,0]=eo*er[0];
Apoisson[n-1,n-1]=-(eo*er[n-2]+eo*er[n-1]);
#Apoisson[n-1,n-2]=Apoisson[n-1,n-2]*2  % Dirichlet-Neumann boundary condition for 
Apoisson[n-2,n-1]=eo*er[n-1];

off_diag[0,1]=1;
off_diag[1,0]=-1;
off_diag[n-1,n-2]=-1;
off_diag[n-2,n-1]=1;

for i in range(1,n-1):
    Ae[i,i]=-(1/me[i-1]+1/me[i+1])
    Ah[i,i]=-(1/mh[i-1]+1/mh[i+1])
    Apoisson[i,i]=-(eo*er[i-1]+eo*er[i+1])
    Ae[i+1,i]=1/me[i-1]
    Ah[i+1,i]=1/mh[i-1]
    Apoisson[i+1,i]=eo*er[i-1]
    Ae[i-1,i]=1/me[i+1]
    Ah[i-1,i]=1/mh[i+1]
    Apoisson[i-1,i]=eo*er[i+1]
    off_diag[i,i+1]=1
    off_diag[i+1,i]=-1

#kp=matrix(zeros((2*n,2*n)),dtype=complex)
kp=np.zeros(shape=(2*n,2*n), dtype=complex)
kp[0:n, 0:n]=-Ae/m/e*(hbar**2)/2/(space**2) + Pe
kp[n:2*n, n:2*n]= Ah/m/e*(hbar**2)/2/(space**2) - Ph
kp[0:n, n:2*n]=off_diag*(-1j)*hbar*k/2/space
kp[n:2*n, 0:n]=kp[0:n, n:2*n]
Apoisson=Apoisson/(space**2);

ev, ef = linalg.eigh(kp)      #eigen value of Hermitian operator
Cb_temp=abs(ev-min(Cb))
Vb_temp=abs(ev+min(Vb))
Cb_E_address=np.where(Cb_temp==min(abs(Cb_temp)))
Vb_E_address=np.where(Vb_temp==min(abs(Vb_temp)))
Cb_energy=ev[Cb_E_address]
Vb_energy=ev[Vb_E_address]
Ee_address=Cb_E_address[0]
Eh_address=Vb_E_address[0]

psi_e=ef[0:n, Ee_address]  
psi_h=ef[n:2*n, Eh_address]
psi_e_sq=psi_e*conjugate(psi_e)   
psi_h_sq=psi_h*conjugate(psi_h)
norm_e=sum(psi_e_sq)*space
norm_h=sum(psi_h_sq)*space
psi_e_sq_norm=psi_e_sq/space
psi_h_sq_norm=psi_h_sq/space
psi_e_norm=psi_e/sqrt(norm_e)   #normalized wf
psi_h_norm=psi_h/sqrt(norm_h)

cc_psi_e=conjugate(psi_e_norm)
cc_psi_h=conjugate(psi_h_norm)
eq_x1=cc_psi_e*cc_psi_h

import fourier
eq_x1_FT, f0= fourier.fourier (eq_x1,space)

eq_x1_FT=eq_x1_FT/sqrt(2*pi)

Vq=zeros(f0.size)
for q in range(f0.size):
    xq=1j*f0[q]*delta
    Vq[q]=-expi(-xq)*math.e**(xq)-expi(xq)*math.e**(-xq)
Vq_prime=Vq*eq_x1_FT
Vq_prime=Vq_prime[:, 0]

psi_h_kf, f0=fourier.fourier(cc_psi_h,space)
psi_h_kf=psi_h_kf/sqrt(2*pi)
kf=sqrt(2*m*e*(Cb_energy-Eg/2+2*abs(Vb_energy+Eg/2)+Eg))/hbar
kf_address=np.where(abs(f0-kf)==min(abs(f0-kf)))
psi_h_kf_value=psi_h_kf[kf_address]

Vq_prime_integral_range_upper=np.where(abs(f0-1/Qw)==min(abs(f0-1/Qw)))
Vq_prime_integral_range_lower=np.where(abs(f0+1/Qw)==min(abs(f0+1/Qw)))
dq=f0[1]-f0[0]
Vq_prime_integral=sum(Vq_prime[Vq_prime_integral_range_lower[0]:Vq_prime_integral_range_upper[0]])*dq
Mif=(e**2)*sqrt(2)/8/sqrt(2*L)/(pi**2)/eo/er1*psi_h_kf_value*Vq_prime_integral
Eex=f0[kf_address]**(2)*hbar**(2)/2/m/mh1/e-Vbo
dos_Ef=1/pi/hbar*sqrt(m*mh1/2/Eex/e)
kA=2*pi/hbar*abs(Mif)**(2)*dos_Ef
Aug_life=1/kA
e_h_multiple=conjugate(psi_e_norm)*psi_h_norm
spatial_sum=sum(e_h_multiple)*space
overlap_integral=abs(spatial_sum)**2
delta_E=Cb_energy-Cb[(n-1)/2]-Vb_energy-Vb[(n-1)/2]

#import nextpow2    
#Nfft=nextpow2.nextpow2(n)
#Fy = np.fft.fft(eq_x1,Nfft,0)*space
#y0 = fftshift(Fy)
#f0 = np.linspace(-Nfft/2,Nfft/2,Nfft)/space/Nfft*2*pi
#plot(f0,abs(y0))
