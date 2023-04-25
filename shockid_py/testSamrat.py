# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

from pipreadmods import pipread
import matplotlib.pyplot as plt
import numpy as np
from shockid import shockid


#fname='/home/bs428/Downloads/thermal_tearing2d/adiabatic'
fname='non_adiabatic/t47/'

with open(''.join((fname,'v1_t47.dat'))) as file:
    next(file)
    v1 = [[float(digit) for digit in line.split()] for line in file]
    
with open(''.join((fname,'v2_t47.dat'))) as file:
    next(file)
    v2 = [[float(digit) for digit in line.split()] for line in file]
    
with open(''.join((fname,'rho_t47.dat'))) as file:
    next(file)
    rho = [[float(digit) for digit in line.split()] for line in file]
    
with open(''.join((fname,'b1_t47.dat'))) as file:
    next(file)
    b1 = [[float(digit) for digit in line.split()] for line in file]
    
with open(''.join((fname,'b2_t47.dat'))) as file:
    next(file)
    b2 = [[float(digit) for digit in line.split()] for line in file]
    
with open(''.join((fname,'p_t47.dat'))) as file:
    next(file)
    pr = [[float(digit) for digit in line.split()] for line in file]

v1=np.asarray(v1)
v2=np.asarray(v2)
b1=np.asarray(b1)
b2=np.asarray(b2)
rho=np.asarray(rho)
pr=np.asarray(pr)

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


fig, ax = plt.subplots(figsize=(9, 6))
#plt.contourf(divv)
#plt.contourf(divv[:,950:1150],levels=np.linspace(-0.055,-0.00001,100))
#ax.contourf(np.log10(-divv.T),levels=101)
ax.contourf(np.log10(rho).T,levels=101,cmap='Greys')
plt.ylim([950, 1100])

shocks=shockid(xg,yg,0.0,rho.T,v1.T,v2.T,0.0,b1.T,b2.T,0.0,pr.T,smthfac=0,nproc=6)

plt.plot(shocks['slow'][:,1],shocks['slow'][:,0],color='r',marker='+',linestyle='',markersize=2.5)
plt.plot(shocks['fast'][:,1],shocks['fast'][:,0],color='b',marker='+',linestyle='',markersize=2.5)
plt.plot(shocks['int4'][:,1],shocks['int4'][:,0],color='g',marker='+',linestyle='',markersize=2.5)

savename='samrat_test_47.png'
plt.savefig(savename)