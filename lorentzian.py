# -*- coding: utf-8 -*-
"""
Created on Sun Sep 01 21:29:54 2013

@author: KyoungWon
"""

ip = get_ipython()
ip.magic("reset -f")
from pylab import *   #include numpy, plotting function
data = genfromtxt('C:\Users\KyoungWon\Google Drive\Python\Auger\CdSe_length.csv', delimiter = ',') 

mean=500
new_mean=np.linspace(400,600,201)
FWHM=30
X=np.linspace(0 , 1000, 10001)


dx=X[-1]/X.size
lorentzian=FWHM/2/pi/((X-mean)**2+(FWHM/2)**2)

#define dFWHM over dE relation  from experimental data
dE=np.linspace(-100,100,201)
FWHM_ratio=0.35/50*dE+1.02
dFWHM=FWHM*FWHM_ratio

#plot(dE, dFWHM)
#plot(X,lorentzian)
left_int=zeros((dE.size))
right_int=zeros((dE.size))
ratio_int=zeros((dE.size))
new_spectrum=zeros((X.size,dE.size))
for i in range(dE.size):
    new_spectrum[:,i]=dFWHM[i]/2/pi/((X-(mean+dE[i]))**2+(dFWHM[i]/2)**2)
    left_int[i]=sum(new_spectrum[:5000,i])*dx
    right_int[i]=sum(new_spectrum[5001:,i])*dx
    ratio_int[i]=right_int[i]/left_int[i]
    
#ROI_sp_change=ratio_int[50:150]

#plt.semilogy(right_int)

plt.plot(X[3000:8000],new_spectrum[3000:8000,100], lw=4, color='red')
plt.plot(X[3000:8000],new_spectrum[3000:8000,150], lw=4, color='blue')
XX=X[3000:8000]
YY=new_spectrum[3000:8000,150]
plt.fill_between(XX,YY,where=XX<500, color='green', alpha=0.5)
plt.fill_between(XX,YY,where=XX>=500, color='blue', alpha=0.5)
#plt.fill_between(X,new_spectrum[:,130],where=X>500, color='blue', alpha=0.5)   
