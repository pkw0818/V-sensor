# -*- coding: utf-8 -*-
"""
Created on Wed Jan 29 04:29:38 2014

@author: KyoungWon
"""

# -*- coding: utf-8 -*-
"""
Created on Sun Sep 08 20:12:27 2013

@author: KyoungWon
"""

ip = get_ipython()
ip.magic("reset -f")
from pylab import *   #include numpy, plotting function
data = genfromtxt('C:/Users/KyoungWon/SkyDrive/V sensor/2d_data/CdSe_ratio_lamda.csv', delimiter = ',') 


length_size=size(data,0)
QY_size=size(data,1)

length=linspace(4,12,9)
QY=linspace(-80,80,17)

X , Y = meshgrid(length,QY)     #X, Y should be 2D

c = plt.contourf(X,Y,data, linspace(100,119,20))
plt.rcParams['font.size'] = 20
b = plt.colorbar(c, orientation='vertical')
lx = plt.xlabel("Length [nm]", fontsize=20)
#ly = plt.ylabel(r'$\frac{\Delta \tau}{\tau}$'   "[%]", fontsize=40)
#ax = plt.axis([0,10,0,8])
plt.show()
