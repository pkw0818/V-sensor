# -*- coding: utf-8 -*-
"""
Created on Sat Jul 27 10:03:07 2013

@author: KyoungWon
"""
ip = get_ipython()
ip.magic("reset -f")
from pylab import *   #include numpy, plotting function
from scipy.interpolate import interp1d

t=np.linspace(0,1000,1001)/1000
x = 0.7*np.sin(2*pi*50*t) + np.sin(2*pi*120*t)
y=x+4*random_sample(size=(1001))
#plot(y[0:50])
NFFT=1024
Yy=np.fft.fft(y,1024)/1000
f=500*np.linspace(0,1,513)
Yy2=fftshift(Yy)
plot(2*abs(Yy[0:513]))

#aaa=np.fft.fft(eq_x1,2048)*space
#bbb=np.fft.fftfreq(2048)