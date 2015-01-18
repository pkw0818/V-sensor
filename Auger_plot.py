# -*- coding: utf-8 -*-
"""
Created on Wed Aug 21 16:18:53 2013

@author: KyoungWon
"""


    
    

X , Y = meshgrid(r[:-1],z)     #X, Y should be 2D
fig1=plt.figure()
ax = fig1.gca(projection='3d')

surf = ax.plot_surface(X.T, Y.T, ewf, rstride=1, cstride=1, cmap=cm.jet,
        linewidth=0, antialiased=False)
fig1.colorbar(surf, shrink=0.5, aspect=5)   #shrink: size of colorbar

fig2=plt.figure()
ax = fig2.gca(projection='3d')
surf = ax.plot_surface(X.T, Y.T, hwf, rstride=1, cstride=1, cmap=cm.jet,
        linewidth=0, antialiased=False)
fig2.colorbar(surf, shrink=0.5, aspect=5)

X , Y = meshgrid(r,theta)     #X, Y should be 2D
fig3=plt.figure()
ax = fig3.gca(projection='3d')
surf = ax.plot_surface(X, Y, abs(phi_kf_n.T), rstride=1, cstride=1, cmap=cm.jet,
        linewidth=0, antialiased=False)
fig3.colorbar(surf, shrink=0.5, aspect=5)   #shrink: size of colorbar    