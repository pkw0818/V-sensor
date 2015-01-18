# -*- coding: utf-8 -*-
"""
Created on Sun Sep 01 22:45:02 2013

@author: KyoungWon
"""


ip = get_ipython()
ip.magic("reset -f")
from pylab import *   #include numpy, plotting function
data = genfromtxt('C:\Users\KyoungWon\Google Drive\Python\Auger\CdSe_length.csv', delimiter = ',') 

X=np.linspace(0,41,41)
log_data=log10(data)
m,b = polyfit(X,log_data,1)
Y=m*X+b
eY=10**(Y)
plt.semilogy(eY)
plt.semilogy(data)

