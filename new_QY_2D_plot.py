# -*- coding: utf-8 -*-
"""
Created on Sun Sep 08 20:12:27 2013

@author: KyoungWon
"""

ip = get_ipython()
ip.magic("reset -f")
from pylab import *   #include numpy, plotting function
data = genfromtxt('C:\Users\KyoungWon\Google Drive\Python\Auger\New_QY_Zn_ki_blue.csv', delimiter = ',') 

length_size=size(data,0)
QY_size=size(data,1)

length=linspace(4,12,42)
QY=linspace(0.05,0.95,19)

X , Y = meshgrid(length[:-1],QY)     #X, Y should be 2D

c = plt.contourf(X,Y,data.T, linspace(0,2.1,22))
plt.rcParams['font.size'] = 20
b = plt.colorbar(c, orientation='vertical')
lx = plt.xlabel("Length [nm]", fontsize=20)
ly = plt.ylabel(r'$original X QY$', fontsize=20)
#ax = plt.axis([0,10,0,8])
plt.show()

