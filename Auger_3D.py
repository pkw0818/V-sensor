# -*- coding: utf-8 -*-
"""
Created on Mon Aug 05 16:47:28 2013

@author: KyoungWon
"""

ip = get_ipython()
ip.magic("reset -f")
from pylab import *   #include numpy, plotting function
from mpl_toolkits.mplot3d import Axes3D  # 3D plotting
from matplotlib import cm                # color map, I use 'jet'
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
from scipy.sparse.linalg import eigs
from scipy.sparse import csc_matrix, eye
from scipy.special import jn
# Define Constant 
V=0     # to define electric field
hbar=6.626e-34/2/math.pi
plank=6.626e-34
e=1.6e-19
eo=8.85e-12;
me=9.1e-31;
me1=0.13*me;   # CdSe me=0.13 mh=0.45       # CdS me:0.21 mh:0.8
mh1=0.45*me;    # CdTe me=0.1  mh=0.4
#me2=0.04*9.1e-31;   
#mh2=1*9.1e-31;
er1=9.3*eo;           # CdSe: 9.3(WZ) 10.2(ZB)        # CdTe: 10.2(ZB)
er2=12.0*eo;           # CdS: 8.9 or 7 or 6.3
Eg=1.0;
Cb1=4;           # CdSe qX: 4.6     Eg: 2.1
Vb1=4;             # CdTe qX: 4.1     Eg: 1.8
k_num=math.sqrt(21./2);           # K dot P matrix element
delta=1e-12;        

#Define dimensiton and mesh
r=4e-9
dr=2e-10
r=np.linspace(0, r, r/dr+1)
z=12e-9
dz=2e-10
z=np.linspace(0,z, z/dz+1)
r_bound=2e-9  #r_max boundary of NR,    too reduce computation time
z_bound=2e-9  #Z_min and max boundary of NR
r_limit=int(r_bound/dr)
z_limit=int(z_bound/dz)
m=r.size
n=z.size
drdz=dr/dz
dzdr=dz/dr
theta=np.linspace(0,180,17)

Cb=np.zeros(shape=(m,n))
Vb=np.zeros(shape=(m,n))
er_org=np.ones(shape=(m,n))*er1
me_org=np.ones(shape=(m,n))*me1
mh_org=np.ones(shape=(m,n))*mh1
er=np.ones(shape=(m+2,n+2))*er1
me=np.ones(shape=(m+2,n+2))*me1
mh=np.ones(shape=(m+2,n+2))*mh1


geo=np.ones(shape=(m,n))*2
for i in range(m):
    for j in range(n):
        if r[i]**2 + (z[j]-4e-9)**2 <= r_bound**2:
            geo[i,j]=1
        elif r[i]**2 + (z[j]-8e-9)**2 <= r_bound**2:
            geo[i,j]=1
geo[0:r_limit+1, z_limit:n-z_limit]=1

for i in range(m):
    for j in range(n):
        if geo[i,j]==1:   # CdSe
            Cb[i,j]=0
            Vb[i,j]=0
            er_org[i,j]=er1
            me_org[i,j]=me1
            mh_org[i,j]=mh1
        elif geo[i,j]==2:   #Vacuum
            Cb[i,j]=Cb1
            Vb[i,j]=Vb1
            er_org[i,j]=eo
            me_org[i,j]=me1
            mh_org[i,j]=mh1

Cb=Cb+Eg/2
Vb=Vb+Eg/2

er_array=zeros((m-1)*n)
for i in range(m-1):
    er_array[i*n:(i+1)*n]=er_org[i+1,:]
            
bndr=zeros(shape=(m,n))
for i in range(m-1):
    for j in range(n-2):
        if er_org[i,j] > er_org[i+1,j]:
            bndr[i+1,j]=1
        elif er_org[i,j] > er_org[i,j+1]:
            bndr[i,j+1]=1
        elif er_org[i,j] < er_org[i,j+1]:
            bndr[i,j]=1
bndr=bndr*((er1+er1)/2-er1)
er_org=er_org+bndr
er[1:m+1,1:n+1]=er_org;
er[0,1:n+1]=er_org[0,:];
me[1:m+1,1:n+1]=me_org;
me[0,1:n+1]=me_org[0,:];
mh[1:m+1,1:n+1]=mh_org;
mh[0,1:n+1]=mh_org[0,:];

Vx=np.linspace(0,V,n)
PE_m=np.ones(m)
PE=PE_m.reshape(m,1)*Vx.reshape(1,n)

# Define Hamiltonian

Ae=np.zeros((n,n,m))
Ah=np.zeros((n,n,m))
r_mid=np.zeros((m+2))
r_mid[0]=0
r_mid[1:m+1]=r+dr/2
r_mid[m+1]=r_mid[m]+dr


for k in range(m): # m X m (radial), Each element is nXn (axial) submatrix
    a=k+1
    for j in range(1,n-1):
        b=j+1
        Ae[j,j,k]=-(1/dr**2/me[a-1,b]+1/dr**2/me[a+1,b]+1/dz**2/me[a,b-1]+1/dz**2/me[a,b+1])
        Ae[j,j-1,k]=1/dz**2/me[a,b-1]
        Ae[j,j+1,k]=1/dz**2/me[a,b+1]
        Ah[j,j,k]=-(1/dr**2/mh[a-1,b]+1/dr**2/mh[a+1,b]+1/dz**2/mh[a,b-1]+1/dz**2/mh[a,b+1])
        Ah[j,j-1,k]=1/dz**2/mh[a,b-1]
        Ah[j,j+1,k]=1/dz**2/mh[a,b+1]
    Ae[0,0,k]=-(1/dr**2/me[a-1,1]+1/dr**2/me[a+1,1]+1/dz**2/me[a,0]+1/dz**2/me[a,2]);  
    Ae[0,1,k]=1/me[a,2]/dz**2
    Ae[n-1,n-1,k]=-(1/dr**2/me[a-1,n+1]+1/dr**2/me[a+1,n+1]+1/dz**2/me[a,n]+1/dz**2/me[a,n+1]);
    Ae[n-1,n-2,k]=1/me[a,n]/dz**2
    Ah[0,0,k]=-(1/dr**2/mh[a-1,1]+1/dr**2/mh[a+1,1]+1/dz**2/mh[a,0]+1/dz**2/mh[a,2]);  
    Ah[0,1,k]=1/mh[a,2]/dz**2
    Ah[n-1,n-1,k]=-(1/dr**2/mh[a-1,n+1]+1/dr**2/mh[a+1,n+1]+1/dz**2/mh[a,n]+1/dz**2/mh[a,n+1]);
    Ah[n-1,n-2,k]=1/mh[a,n]/dz**2

Be_down=np.zeros((n,n,m))
Bh_down=np.zeros((n,n,m))
Be_up=np.zeros((n,n,m))
Bh_up=np.zeros((n,n,m))

for k in range(1,m):
    for j in range(n):
        a=k+1
        b=j+1
        Be_down[j,j,k]=1/me[a-1,b]/dr**2
        Bh_down[j,j,k]=1/mh[a-1,b]/dr**2

for k in range(m):
    for j in range(n):
        a=k+1
        b=j+1
        Be_up[j,j,k]=1/me[a+1,b]/dr**2
        Bh_up[j,j,k]=1/mh[a+1,b]/dr**2

CeM=np.zeros((m*n,m*n))
ChM=np.zeros((m*n,m*n))
CbMa=np.zeros((m*n,m*n))
VbMa=np.zeros((m*n,m*n))
Pfield=np.zeros((m*n,m*n))
for i in range(m):
    for j in range(m):
        if i==j:
            CeM[i*n:(i+1)*n,i*n:(i+1)*n]=Ae[:,:,i]
            ChM[i*n:(i+1)*n,i*n:(i+1)*n]=Ah[:,:,i]
            CbMa[i*n:(i+1)*n,i*n:(i+1)*n]=diag(Cb[i,:])
            VbMa[i*n:(i+1)*n,i*n:(i+1)*n]=diag(Vb[i,:])
            Pfield[i*n:(i+1)*n,i*n:(i+1)*n]=diag(PE[i,:])
        elif i==j+1:
            CeM[i*n:(i+1)*n,j*n:(j+1)*n]=Be_down[:,:,i]
            ChM[i*n:(i+1)*n,j*n:(j+1)*n]=Bh_down[:,:,i]
        elif i==j-1:
            CeM[i*n:(i+1)*n,j*n:(j+1)*n]=Be_up[:,:,i]
            ChM[i*n:(i+1)*n,j*n:(j+1)*n]=Bh_up[:,:,i]
            
CeM[n:2*n,n:2*n]=CeM[n:2*n,n:2*n]+CeM[n:2*n,:n]*4/3
CeM[n:2*n,2*n:3*n]=CeM[n:2*n,2*n:3*n]-CeM[n:2*n,:n]/3

ChM[n:2*n,n:2*n]=ChM[n:2*n,n:2*n]+ChM[n:2*n,:n]*4/3
ChM[n:2*n,2*n:3*n]=ChM[n:2*n,2*n:3*n]-ChM[n:2*n,:n]/3
            
Ke2=CeM[n:m*n,n:m*n]
Kh2=ChM[n:m*n,n:m*n]
CbM=CbMa[n:m*n,n:m*n]
VbM=VbMa[n:m*n,n:m*n]
PfieldM=Pfield[n:m*n,n:m*n]

He=-Ke2/e*hbar**2/2 + CbM + PfieldM
Hh=Kh2/e*hbar**2/2  - VbM + PfieldM

#He=Ke+PfieldM+CbM
#Hh=-Kh+PfieldM-VbM

off_1=eye(n*(m-1), k=1).todense()
off_2=eye(n*(m-1), k=-1).todense()
off_3=eye(n*(m-1), k=n).todense()
off_4=eye(n*(m-1), k=-n).todense()
off=off_1/dz-off_2/dz+off_3/dr-off_4/dr
#off=off.todense()
#off[n:2*n,n:2*n]=off[n:2*n,n:2*n]+off[n:2*n,:n]*4/3
#off[n:2*n,2*n:3*n]=off[n:2*n,2*n:3*n]-off[n:2*n,:n]/3
off=off*(-1j)*hbar*k_num/2

kp=np.zeros((2*n*(m-1),2*n*(m-1)), dtype=complex)
kp[:n*(m-1), :n*(m-1)]=He
kp[n*(m-1):2*n*(m-1),n*(m-1):2*n*(m-1)]=Hh
kp[:n*(m-1),n*(m-1):2*n*(m-1)]=off
kp[n*(m-1):2*n*(m-1),:n*(m-1)]=off

sp_kp=csc_matrix(kp)
#ev, ef = linalg.eig(kp) 

#spkp=coo_matrix(kp)
#sHe=coo_matrix(He)
#sHh=coo_matrix(Hh)
del He, Hh, Ke2, Kh2, CbM, VbM, PfieldM, CeM, ChM, CbMa, VbMa, Ae, Ah, Be_down, Bh_down, Be_up, Bh_up
del Cb, Vb, er, me, mh, er_org, me_org, mh_org, geo, er_array, off
#ev, ef = linalg.eig(kp) 
#evh, efh  = linalg.eig(Hh) 
eev, eef= eigs(sp_kp, 1, sigma=Eg/2, which='LM')
hev, hef= eigs(sp_kp, 1, sigma=-Eg/2, which='LM')
#hev, hef = eigs(sHh,1, which='SR')
eef=eef[:(m-1)*n,0]
hef=hef[(m-1)*n:(m-1)*n*2,0]

eef2=eef*conjugate(eef)
hef2=eef*conjugate(hef)
eef2_rdr=zeros(((m-1)*n), dtype=complex)
hef2_rdr=zeros(((m-1)*n), dtype=complex)
for i in range(m-1):
    eef2_rdr[i*n:(i+1)*n]=eef2[i*n:(i+1)*n]*r[i+1]
    hef2_rdr[i*n:(i+1)*n]=hef2[i*n:(i+1)*n]*r[i+1]
norm_e=sum(eef2_rdr)*dr*dz*2*pi
norm_h=sum(hef2_rdr)*dr*dz*2*pi
eef2_n=eef2/norm_e
hef2_n=hef2/norm_h
eef_n=eef/sqrt(norm_e)          #normalized wavefunction
hef_n=hef/sqrt(norm_h)
cc_eef=conjugate(eef_n)         # complex conjugated wavefunction
cc_hef=conjugate(hef_n)
ewf=eef_n.reshape((m-1),n)
hwf=hef_n.reshape((m-1),n)
cc_ewf=cc_eef.reshape((m-1),n)
cc_hwf=cc_hef.reshape((m-1),n)  #cc_hwf[0,:] is r at r[1]


cc_ewf_expanding=zeros((r_limit,n-2*z_limit,theta.size), dtype=complex)
cc_hwf_expanding=zeros((r_limit,n-2*z_limit,theta.size), dtype=complex)  # just for calculation later
for i in range(theta.size):                         #expanding along the theta
    cc_hwf_expanding[:,:,i]=cc_hwf[:r_limit,z_limit:n-z_limit]       #ignore r=0
    cc_ewf_expanding[:,:,i]=cc_ewf[:r_limit,z_limit:n-z_limit]
 

distance1=zeros((r_limit+1,r_limit+1,theta.size))  #z1=z2, 
for i in range(r_limit+1):     #r1
    for j in range (i+1):          #r2
        for k in range(theta.size):         #theta=theta1=theta2
            distance1[i,j,k]=r[i]**2+r[j]**2-2*r[i]*r[j]*cos(radians(theta[k]))

for i in range(theta.size):
    temp=distance1[:,:,i]
    temp=temp+triu(temp.T,1)
    distance1[:,:,i]=temp

Vcoul=zeros((r_limit+1,r_limit+1,theta.size,n-2*z_limit))
for r1 in range(r_limit+1):    #since r[0]=0   #r1 at r=0 should be taken
    for r2 in range(r_limit+1):     #r[0] should be ignored
        for delta_theta in range(theta.size):
            for delta_z in range(n-2*z_limit):
                Vcoul[r1,r2,delta_theta,delta_z]=r[r2]/sqrt(distance1[r1,r2,delta_theta]+z[delta_z]**2)
            
for i in range(r_limit+1):
    Vcoul[i,i,0,0]=0     #r1=r2, z1=z2, theta1=theta2, then Vcoul goes zero
for i in range(theta.size):
    Vcoul[0,0,i,0]=0     #r1=r2=0, z1=z2,  then Vcoul goes zero


kf=sqrt(2*mh1*e*(real(eev)-Eg/2+2*abs(hev+Eg/2)+Eg))/hbar  # [1/m]    
kf=kf[0]            
phi_kf=zeros((m,theta.size), dtype=complex)
#plane wave in cylindrical coordinate, expanding in fourier series of bessel function
phi_kf_fs=zeros((m,theta.size,theta.size), dtype=complex)
for i in range(theta.size):        # nth orde
    for j in range(m):             # r[j]
        for k in range(theta.size):     #theta[k]
            if i==0:
                phi_kf_fs[j,k,i]=jn(i,r[j]*kf)
            else:
                phi_kf_fs[j,k,i]=1j**i*jn(i,r[j]*kf)*2*cos(radians(i*theta[k]))
    phi_kf=phi_kf+phi_kf_fs[:,:,i]  # not normalized yet

#phi_kf normalization
d_theta=radians(2*pi/radians(theta.size-1))
phi_kf2=phi_kf*conjugate(phi_kf)
phi_kf2_rdr=zeros((m,theta.size), dtype=complex)
for i in range(m):
    phi_kf2_rdr[i,:]=phi_kf2[i,:]*r[i]
norm_phi=sum(phi_kf2_rdr)*dr*2*d_theta*z[-1]
phi_kf_n=phi_kf/sqrt(norm_phi)


            
# expanding along the Z        
phi_kf_expanding=zeros((r_limit,n-2*z_limit,theta.size), dtype=complex)
for i in range(n-2*z_limit):
    phi_kf_expanding[:,i,:]=phi_kf_n[1:r_limit+1,:]
   
import dx2_int

dx2_val=zeros((r_limit,n-2*z_limit,theta.size))
for i in range(r_limit):
    for j in range(n-2*z_limit):
        for k in range(theta.size):
            dx2_val[i,j,k]=dx2_int.dx2_int(i,j,k,m,n,theta,Vcoul,r_limit,z_limit,dr,dz,phi_kf_expanding,cc_hwf_expanding)

dx2_val_rdr=zeros((r_limit, n-2*z_limit,theta.size), dtype=complex)
for i in range(r_limit):
    dx2_val_rdr[i,:,:]=dx2_val[i,:,:]*r[i+1]

Mif=sum(cc_ewf_expanding*cc_hwf_expanding*dx2_val_rdr)*2*dr*dz*d_theta*(e**2)*sqrt(2)/4/pi/er1

Eex=kf**(2)*hbar**(2)/2/mh1/e              #[eV]
vol=pi*r[-1]**2*z[-1]
dos_Ef=8*pi*sqrt(2)/(plank**3)*sqrt(mh1**3)*sqrt(Eex*e)*vol         #[1/J)]
kA=2*pi/hbar*abs(Mif)**(2)*dos_Ef              #[s/J * (J)^2 *  * J^-1]
Aug_life=1/kA
#e_h_multiple=conjugate(psi_e_norm)*psi_h_norm
#spatial_sum=sum(e_h_multiple)*space
#overlap_integral=abs(spatial_sum)**2
delta_E=real(eev)-Eg/2+2*abs(hev+Eg/2)

