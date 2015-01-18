# -*- coding: utf-8 -*-
"""
Created on Sat Jan 18 18:32:28 2014

@author: KyoungWon
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
me1=0.13;   # CdSe me=0.13 mh=0.45       # CdS me:0.21 mh:0.8
mh1=0.45;    # CdTe me=0.1  mh=0.4
#me2=0.16;   
#mh2=0.57;
er1=9.3;           # CdSe: 9.3(WZ) 10.2(ZB)        # CdTe: 10.2(ZB)
#er2=8.6;           # CdS: 8.9 or 7 or 6.3
Eg=1.75;             # ZnSe: me0.14 mh0.53 Eg2.72   CBO=0.8   qX:4.09 er=9.2
                     # CdS: me.18 mh0.6 Eg2.45       VBO=0.52    qX=4.79 er=8.6
Cbo=0;           # CdSe qX: 4.95     Eg: 1.75 er=9.5
Vbo=0;             # CdTe qX: 4.1     Eg: 1.8
k=math.sqrt(20./2);           # K dot P matrix element
delta=1e-15;              # use for avoid divergence in V(q)
hole_barrier=1     # 1eV

#Define Geometry
L=10e-9;                # Sim. length L >> QW (confined area)
Qw = 4e-9
space=0.5e-10;
X = np.linspace(-L , L, 2*L/space+1)
n= X.size
region=[-L,-Qw,Qw,L]
indx=np.zeros(len(region))
for var in range(len(region)):
    indx[var]=(region[var]-region[0])/space

CB_value=[4.95,0,4.95]
VB_value=[4.95,0,4.95]
er_value=[80,er1,80]
er_value_homo=[er1,er1,er1]  # for image charge
me_value=[me1,me1,me1]
mh_value=[mh1,mh1,mh1]

import piecewise2      #make a piecewise function with x as region, y as value
Cb=piecewise2.piecewise2(indx,CB_value)
Vb=piecewise2.piecewise2(indx,VB_value)
er=piecewise2.piecewise2(indx,er_value)
er_homo=piecewise2.piecewise2(indx,er_value_homo)
me=piecewise2.piecewise2(indx,me_value)
mh=piecewise2.piecewise2(indx,mh_value)
#Field=0
Field=8e-2/Qw*L
PE=np.linspace(0,Field,n)
if Field==0:
    PE=np.zeros(n)

PE=PE-(Field/2)
Cb1=Cb+Eg/2
Vb1=Vb+Eg/2

Cb2=Cb+PE+Eg/2
Vb2=Vb-PE+Eg/2


#Define Hamiltonian cell
Ae=mat(zeros((n,n)))
Ah=mat(zeros((n,n)))
Apoisson=mat(zeros((n,n)))
Pe1=diag(Cb1)
Ph1=diag(Vb1)
Pe2=diag(Cb2)
Ph2=diag(Vb2)
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
kp[0:n, 0:n]=-Ae/m/e*(hbar**2)/2/(space**2) + Pe1
kp[n:2*n, n:2*n]= Ah/m/e*(hbar**2)/2/(space**2) - Ph1
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
Ee_address=int(Ee_address[0])
Eh_address=int(Eh_address[0])

psi_e=ef[0:n, 401]  
psi_h=ef[n:2*n, 400]
psi_e_sq=psi_e*conjugate(psi_e)   
psi_h_sq=psi_h*conjugate(psi_h)
norm_e=sum(psi_e_sq)*space
norm_h=sum(psi_h_sq)*space
psi_e_sq_norm=psi_e_sq/space
psi_h_sq_norm=psi_h_sq/space
psi_e_norm=psi_e/sqrt(norm_e)   #normalized wf
psi_h_norm=psi_h/sqrt(norm_h)


## constructing distance matrix : |r-X| for calculating coulomb potential
distance_lower=zeros((n,n))
for i in range(n):
    temp=eye(n,k=i)*i
    distance_lower=distance_lower+temp
distance_upper=distance_lower.T
distance=distance_lower+distance_upper
distance=distance+eye(n,k=0)
distance=1/distance-eye(n,k=0)


## self-consistency
old_Cb_E=0
for i in range(10):
    
    Velectron=zeros((n))
    Vhole=zeros((n))
    e_over_r=zeros((n))
    h_over_r=zeros((n))
    for j in range(n):
        e_over_r=multiply(distance[j,:],psi_e_sq_norm)
        h_over_r=multiply(distance[j,:],psi_h_sq_norm)
        Velectron[j]=sum(e_over_r)/er[j]
        Vhole[j]=sum(h_over_r)/er[j]
    Velectron=Velectron*e/4/pi/eo
    Vhole=Vhole*e/4/pi/eo

    
    
    Velectron_mat=diag(Velectron)
    Vhole_mat=diag(Vhole)
    kkp=np.zeros(shape=(2*n,2*n), dtype=complex)
    kkp[0:n, 0:n]=-Ae/m/e*(hbar**2)/2/(space**2) + Pe2 - Vhole_mat
    kkp[n:2*n, n:2*n]= Ah/m/e*(hbar**2)/2/(space**2) - Ph2 + Velectron_mat
    kkp[0:n, n:2*n]=off_diag*(-1j)*hbar*k/2/space
    kkp[n:2*n, 0:n]=kkp[0:n, n:2*n]
    
    #kkp[0:n, 0:n]=kp[0:n, 0:n]-Vhole_mat
    #kkp[n:2*n, n:2*n]= kp[n:2*n, n:2*n]+Velectron_mat
    
    ev, ef = linalg.eigh(kkp)      #eigen value of Hermitian operator
    
    Cb_temp=min(Cb)
    Vb_temp=min(Vb)
    Cb_temp2=ev-Cb_temp
    Vb_temp2=ev*-1-Vb_temp
    Cb_temp3=100
    Vb_temp3=100
    for k1 in range (n,n*2):
        e_abs=abs(ef[0:n, k1])
        e_derivative=zeros(n-1)
        for kk in range(n-1):
            if e_abs[kk+1]-e_abs[kk] > 0:
                e_derivative[kk]=1
            else:
                e_derivative[kk]=-1
             
        ll=0
        for l in range(n-2):
            if e_derivative[l]*e_derivative[l+1] < 0:
                ll=ll+1
        if ll==1:
            Cb_E_address=k1
    
    for k2 in range (0,n):
        h_abs=abs(ef[n:n*2, k2])
        h_derivative=zeros(n-1)
        for kk in range(n-1):
            if h_abs[kk+1]-h_abs[kk] > 0:
                h_derivative[kk]=1
            else:
                h_derivative[kk]=-1
        ll=0
        for l in range(n-2):
            if h_derivative[l]*h_derivative[l+1] < 0:
                ll=ll+1
        if ll==1:
            Vb_E_address=k2
      
    
        #if Cb_temp2[k] > -0.2 and Cb_temp2[k] < Cb_temp3:
        #    Cb_temp3=Cb_temp2[k]
        #    Cb_E_address=k
        #if Vb_temp2[k] > -0.2 and Vb_temp2[k] < Vb_temp3:
        #    Vb_temp3=Vb_temp2[k]
        #    Vb_E_address=k
    if i == 0:
        Vb_E_address=Vb_E_address[0]
        Vb_E_address=int(Vb_E_address)
    #Cb_E_address=int(Cb_E_address)    
    #Cb_temp=abs(ev-min(Cb))
    #Vb_temp=abs(ev+min(Vb))
    #Cb_E_address=np.where(Cb_temp==min(abs(Cb_temp)))
    #Vb_E_address=np.where(Vb_temp==min(abs(Vb_temp)))
    Cb_energy=ev[Cb_E_address]
    Vb_energy=ev[Vb_E_address]
    #Ee_address=Cb_E_address[0]
    #Eh_address=Vb_E_address[0]
    #Ee_address=float(Ee_address[0])
    #Eh_address=float(Eh_address[0])
    
    
    psi_e=ef[0:n, Cb_E_address]  
    psi_h=ef[n:2*n, Vb_E_address]
    


    
    psi_e_sq=psi_e*conjugate(psi_e)   
    psi_h_sq=psi_h*conjugate(psi_h)
    norm_e=sum(psi_e_sq)*space
    norm_h=sum(psi_h_sq)*space
    psi_e_sq_norm=psi_e_sq/space
    psi_h_sq_norm=psi_h_sq/space
    psi_e_norm=psi_e/sqrt(norm_e)   #normalized wf
    psi_h_norm=psi_h/sqrt(norm_h)
    print Cb_energy
    plot(psi_e_sq_norm)

    if abs(float(Cb_energy-old_Cb_E)) < 0.01:
        break
    old_Cb_E=Cb_energy
        
Velectron=zeros((n))
Vhole=zeros((n))
for j in range(n):
    e_over_r=multiply(distance[j,:],psi_e_sq_norm)
    h_over_r=multiply(distance[j,:],psi_h_sq_norm)
    Velectron[j]=sum(e_over_r)/er[j]
    Vhole[j]=sum(h_over_r)/er[j]
Velectron=Velectron*e/4/pi/eo
Vhole=Vhole*e/4/pi/eo
figure(2)
plot(psi_e_sq_norm)
plot(psi_h_sq_norm)
figure(3)
plot(Velectron)
plot(Vhole)
figure(4)
plot(Cb-Vhole)
plot(Vb-Velectron)




cc_psi_e=conjugate(psi_e_norm)
cc_psi_h=conjugate(psi_h_norm)
eq_x1=cc_psi_e*cc_psi_h

Eex=Eg-hole_barrier+(Cb_energy-Eg/2)+2*abs(Vb_energy+Eg/2)
kf=sqrt(2*m*mh1*Eex*e)/hbar  # [1/m]
phi_F=1/sqrt(2*Qw)*(math.e)**(-1j*kf*X)                      #[1/sqrt(m)]

Ve_lower=zeros(shape=(n,n), dtype=complex)
for x1 in range(n):
    for x2 in range(x1):
        #dx2_int=cc_psi_h[x2]*phi_f[x2]*1/(abs(x1-x2)+delta)
        Ve_lower[x1,x2]=1/(abs(X[x1]-X[x2]))
Ve_upper=np.triu(Ve_lower.T,1)
Ve=Ve_lower+Ve_upper              #[1/m]

#Ve=ones(shape=(n,n), dtype=complex)
dx2_int=zeros(n, dtype=complex)
temp=diag(phi_F*cc_psi_h)                       #[1/m]
for x1 in range(n):
    temp2=Ve[x1,:]                #[1/m]
    dx2_int[x1]=sum(temp*temp2)      #[1m^2]   
#dx2_int=dx2_int                    #[1/m]
Mif=sum(diag(eq_x1*dx2_int))*space**(2)*(e**2)*sqrt(2)/4/pi/eo/er1        #[1/m * C^2 / C^2 * J*m= J]
aaa=temp*temp2
        
#Mif=(e**2)*sqrt(2)/8/sqrt(2*L)/(pi**2)/eo/er1*psi_h_kf_value*Vq_prime_integral
#Eex=kf**(2)*hbar**(2)/2/m/mh1/e              #[eV]
#dos_Ef=1/pi/hbar*sqrt(m*mh1/2/Eex/e)*2*Qw           #[1/J)]
dos_Ef=sqrt(2)*(m*mh1)**(1.5)/pi/pi/hbar/hbar/hbar*sqrt(Eex*e)*2*Qw*e  #3D DOS * Qw * 2
#dos_Ef=m*mh1/pi/hbar/hbar*e   
kA=2*pi/hbar*abs(Mif)**(2)*dos_Ef              #[s/J * (J)^2 *  * J^-1]
Aug_life=1/kA
e_h_multiple=conjugate(psi_e_norm)*psi_h_norm
spatial_sum=sum(e_h_multiple)*space
overlap_integral=abs(spatial_sum)**2
delta_E=Cb_energy-Cb[(n-1)/2]-Vb_energy-Vb[(n-1)/2]