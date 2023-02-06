# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

from pipreadmods import pipread
import matplotlib.pyplot as plt
import numpy as np


fname='/home/bs428/Downloads/thermal_tearing2d/adiabatic'

with open('/home/bs428/Downloads/thermal_tearing2d/adiabatic/v1_t23.5.dat') as file:
    next(file)
    v1 = [[float(digit) for digit in line.split()] for line in file]
    
with open('/home/bs428/Downloads/thermal_tearing2d/adiabatic/v2_t23.5.dat') as file:
    next(file)
    v2 = [[float(digit) for digit in line.split()] for line in file]
    
v1=np.array(v1)
v2=np.array(v2)

margin=2
egx=np.size(v1[0,:])
egy=np.size(v1[:,0])
dx=1
dy=1

divv=(v2[margin+1:egx-margin+1,margin:egy-margin]-\
      v2[margin-1:egx-margin-1,margin:egy-margin])/(2.0*dy) \
    +(v1[margin:egx-margin,margin+1:egy-margin+1]-\
      v1[margin:egx-margin,margin-1:egy-margin-1])/(2.0*dx) 

#plt.contourf(divv)
#plt.contourf(divv[:,950:1150],levels=np.linspace(-0.055,-0.00001,100))
plt.contourf(divv[:,950:1100],nlevels=101)