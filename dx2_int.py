# -*- coding: utf-8 -*-
"""
Created on Tue Aug 13 16:53:44 2013

@author: KyoungWon
"""
from pylab import * 
def dx2_int(r1,z1,theta1,m,n,theta,Vcoul,r_limit,z_limit,dr,dz,phi_kf_expanding,cc_hwf_expanding):     #i,j,k represent index of r1,z1,theta1    
    temp=zeros((r_limit,n-2*z_limit,theta.size))
    d_theta=radians(2*pi/radians(theta.size-1))

    for r2 in range(r_limit):    #since r[0]=0
        for z2 in range(n-2*z_limit):
            delta_z=abs(z1-z2)
            for theta2 in range(theta.size):
                delta_theta=abs(theta1-theta2)
                temp[r2,z2,theta2]=Vcoul[r1+1,r2+1,delta_theta,delta_z]

    dx2_int_value=sum(phi_kf_expanding*cc_hwf_expanding*temp)*2*dr*dz*d_theta               

    return dx2_int_value


