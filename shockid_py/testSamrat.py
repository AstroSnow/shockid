# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

from pipreadmods import pipread
import matplotlib.pyplot as plt
import numpy as np
from shockid import shockid


fname='/home/bs428/Downloads/thermal_tearing2d/adiabatic'

with open('/media/ben/datadrive1/Samrat_Data/non_adiabatic/t58/v1_t58.dat') as file:
    next(file)
    v1 = [[float(digit) for digit in line.split()] for line in file]
    
with open('/media/ben/datadrive1/Samrat_Data/non_adiabatic/t58/v2_t58.dat') as file:
    next(file)
    v2 = [[float(digit) for digit in line.split()] for line in file]
    
with open('/media/ben/datadrive1/Samrat_Data/non_adiabatic/t58/rho_t58.dat') as file:
    next(file)
    rho = [[float(digit) for digit in line.split()] for line in file]
    
v1=np.array(v1)
v2=np.array(v2)

margin=2
egx=np.size(v1[0,:])
egy=np.size(v1[:,0])
dx=1
dy=1

xg=np.linspace(0,egx,egx)
yg=np.linspace(0,egy,egy)

divv=(v1[margin+1:egy-margin+1,margin:egx-margin]-\
      v1[margin-1:egy-margin-1,margin:egx-margin])/(2.0*dy) \
    +(v2[margin:egy-margin,margin+1:egx-margin+1]-\
      v2[margin:egy-margin,margin-1:egx-margin-1])/(2.0*dx) 

#plt.contourf(divv)
#plt.contourf(divv[:,950:1150],levels=np.linspace(-0.055,-0.00001,100))
plt.contourf(np.log10(-divv[:,950:1095].T),levels=101)

[col,row]=shockid(xg,yg,0.0,rho,v2,v1,0.0,0.0,0.0,0.0,0.0)

plt.plot(col,row-950,color='r',marker='+',linestyle='',markersize=0.5)