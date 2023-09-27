#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 25 10:32:03 2023

@author: ben
"""
from pipreadmods import pipread
import matplotlib.pyplot as plt
import numpy as np
from shockid import shockid
import h5py
import gc

#fname='/home/bs428/Downloads/thermal_tearing2d/adiabatic/t47/'
#fname='/media/ben/datadrive1/Samrat_Data/non_adiabatic/t47/'
#fname='non_adiabatic/t47/'

#fname='/media/ben/datadrive1/Reconnection/OZ/testdata.0070.h5'
fname='/media/snow/datadrive/OZdata/MHD_OZ_hr_h5/testdata.0060.h5'

with h5py.File(fname, "r") as f:
    # List all groups
    #print("Keys: %s" % f.keys())

    # Get the data
    #data={}
#    for param in f:
    ref_data = np.array(f['ro_p'])
    ro=np.squeeze(ref_data)
    del ref_data
    gc.collect()
    ref_data = np.array(f['vx_p'])
    vx=np.squeeze(ref_data)
    del ref_data
    gc.collect()
    ref_data = np.array(f['bx'])
    bx=np.squeeze(ref_data)
    del ref_data
    gc.collect()
    ref_data = np.array(f['by'])
    by=np.squeeze(ref_data)
    del ref_data
    gc.collect()
    ref_data = np.array(f['vy_p'])
    vy=np.squeeze(ref_data)
    del ref_data
    gc.collect()
    ref_data = np.array(f['pr_p'])
    pr=np.squeeze(ref_data)
    del ref_data
    gc.collect()

margin=2
egx=np.size(ro[0,:])
egy=np.size(ro[:,0])
dx=1
dy=1

xg=np.linspace(0,egx,egx)
yg=np.linspace(0,egy,egy)

#divv=(v1[margin+1:egy-margin+1,margin:egx-margin]-\
#      v1[margin-1:egy-margin-1,margin:egx-margin])/(2.0*dy) \
#    +(v2[margin:egy-margin,margin+1:egx-margin+1]-\
#      v2[margin:egy-margin,margin-1:egx-margin-1])/(2.0*dx) 


trxs=4500
trxe=8200#egx#8000
trys=7700#6000
trye=9500#egy#8000


fig, ax = plt.subplots(figsize=(9, 6),dpi=300)
#plt.contourf(divv)
#plt.contourf(divv[:,950:1150],levels=np.linspace(-0.055,-0.00001,100))
#ax.contourf(np.log10(-divv.T),levels=101)
ax.contourf(ro[trys:trye,trxs:trxe],levels=101,cmap='Greys')
#plt.ylim([950, 1100])

#stop

nproc=30
shocks=shockid(xg[trxs:trxe],yg[trys:trye],0.0,ro[trys:trye,trxs:trxe],
			   vx[trys:trye,trxs:trxe],vy[trys:trye,trxs:trxe],0.0,
			   bx[trys:trye,trxs:trxe],by[trys:trye,trxs:trxe],0.0,
			   pr[trys:trye,trxs:trxe],smthfac=2,nproc=nproc,convl=0.0001,avecyl=5)

sl=ax.plot(shocks['slow'][:,1],shocks['slow'][:,0],color='r',linestyle='',marker='.',markersize=2.8)
fa=ax.plot(shocks['fast'][:,1],shocks['fast'][:,0],color='b',linestyle='',marker='.',markersize=2.8)
i1=ax.plot(shocks['int1'][:,1],shocks['int1'][:,0],color='y',linestyle='',marker='.',markersize=2.8)
i2=ax.plot(shocks['int2'][:,1],shocks['int2'][:,0],color='m',linestyle='',marker='.',markersize=2.8)
i3=ax.plot(shocks['int3'][:,1],shocks['int3'][:,0],color='c',linestyle='',marker='.',markersize=2.8)
i4=ax.plot(shocks['int4'][:,1],shocks['int4'][:,0],color='g',linestyle='',marker='.',markersize=2.8)
#if nproc == 1:
#plt.plot(shocks['slow'][:,1],shocks['slow'][:,0],color='r',marker='+',linestyle='',markersize=2.5)
#plt.plot(shocks['fast'][:,1],shocks['fast'][:,0],color='b',marker='+',linestyle='',markersize=2.5)
#plt.plot(shocks['int4'][:,1],shocks['int4'][:,0],color='g',marker='+',linestyle='',markersize=2.5)
"""if nproc > 1:
	for j in range (0,nproc):
		plt.plot(shocks['slow'][j][:,1],shocks['slow'][j][:,0],color='r',marker='+',linestyle='',markersize=2.5)
		plt.plot(shocks['fast'][j][:,1],shocks['fast'][j][:,0],color='b',marker='+',linestyle='',markersize=2.5)
		#plt.plot(shocks['int4'][j][:,1],shocks['int4'][j][:,0],color='g',marker='+',linestyle='',markersize=2.5)
"""
savename=''.join((fname,'shocks.png'))
plt.savefig(savename)

#Create a h5 file of the shock data
hf = h5py.File(''.join((fname,'shocks.h5')), 'w')
for f in shocks:
	hf.create_dataset(f, data=shocks[f])
hf.close()