# -*- coding: utf-8 -*-
"""
Created on Sun Jul 28 16:28:14 2013

@author: KyoungWon
"""

def fourier(input, dx):
    ip = get_ipython()
    from pylab import *   #include numpy, plotting function    
    import nextpow2    
    n=input.size
    Nfft=nextpow2.nextpow2(n)
    Fy = np.fft.fft(input,Nfft,0)*dx
    y0 = fftshift(Fy)
    f0 = np.linspace(-Nfft/2,Nfft/2,Nfft)/dx/Nfft*2*pi
    plot(f0,abs(y0))
    return y0, f0
    