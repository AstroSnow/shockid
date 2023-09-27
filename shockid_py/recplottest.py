#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  6 17:26:37 2023

@author: bs428
"""
from pipreadmods import pipread
import matplotlib.pyplot as plt
import numpy as np
from shockid import shockid

#fname='../../Reconnection/MHDtest/Data/'
#fname='../../Reconnection/recdata/mhd_test_5/'
fname='../../Reconnection/recdata/mhd_rad/'
#fname='../../shockstab/Data/'

ds=pipread(fname,70)
#stopplt
nproc=30
    
v1=ds['vx_p']
v2=ds['vy_p']

margin=2
egx=np.size(v1[0,:])
egy=np.size(v1[:,0])
dx=ds['xgrid'][1]-ds['xgrid'][0]
dy=ds['ygrid'][1]-ds['ygrid'][0]

divv=(v1[margin+1:egy-margin+1,margin:egx-margin]-\
      v1[margin-1:egy-margin-1,margin:egx-margin])/(2.0*dy) \
    +(v2[margin:egy-margin,margin+1:egx-margin+1]-\
      v2[margin:egy-margin,margin-1:egx-margin-1])/(2.0*dx) 


T=ds['pr_p']/ds['ro_p']

#plt.contourf(ds['ro_p'],levels=101) ;stop
#plt.contourf(divv[:,950:1150],levels=np.linspace(-0.055,-0.00001,100))
#plt.contourf(ds['pr_p'][400:650,:],levels=101)
#plt.contourf(ds['pr_p']/ds['ro_p'],levels=101)
#plt.contourf(T[350:650,1000:1800],levels=101)

#stop
#[col,row]=shockid(ds['xgrid'],ds['ygrid'],0.0,ds['ro_p'],ds['vx_p'],ds['vy_p'],0.0,ds['bx'],ds['by'],ds['bz'],ds['pr_p'])
#plt.plot(row,col,color='r',linestyle='',marker='.',markersize=2.8)

#subset the data here
xs=800
xe=2000
ys=0000
ye=2000

xg=ds['xgrid'][xs:xe]
yg=ds['ygrid'][ys:ye]
ro=ds['ro_p'][ys:ye,xs:xe]
pr=ds['pr_p'][ys:ye,xs:xe]
vx=ds['vx_p'][ys:ye,xs:xe]
vy=ds['vy_p'][ys:ye,xs:xe]
bx=ds['bx'][ys:ye,xs:xe]
by=ds['by'][ys:ye,xs:xe]
bz=ds['bz'][ys:ye,xs:xe]


#plt.contourf(pr/ro,levels=101,cmap='Greys')
#stop
shocks=shockid(xg,yg,0.0,ro,vx,vy,0.0,bx,by,bz,pr,smthfac=3,nproc=nproc,convl=0.0001,avecyl=5)
fig, ax = plt.subplots(figsize=(9, 6))
plt.contourf(np.log10(ro),levels=101,cmap='Greys')
#stop
plt.plot(shocks['slow'][:,1],shocks['slow'][:,0],color='r',linestyle='',marker='.',markersize=2.8)
plt.plot(shocks['fast'][:,1],shocks['fast'][:,0],color='b',linestyle='',marker='.',markersize=2.8)
plt.plot(shocks['int1'][:,1],shocks['int1'][:,0],color='y',linestyle='',marker='.',markersize=2.8)
plt.plot(shocks['int2'][:,1],shocks['int2'][:,0],color='m',linestyle='',marker='.',markersize=2.8)
plt.plot(shocks['int3'][:,1],shocks['int3'][:,0],color='c',linestyle='',marker='.',markersize=2.8)
plt.plot(shocks['int4'][:,1],shocks['int4'][:,0],color='g',linestyle='',marker='.',markersize=2.8)

#savename='MHDtest5_30_full.png'
#plt.savefig(savename)
