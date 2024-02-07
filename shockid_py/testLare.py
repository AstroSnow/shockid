# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import matplotlib.pyplot as plt
import numpy as np
from shockid import shockid
import h5py

#Load the data in


margin=2
egx=np.size(v1[0,:])
egy=np.size(v1[:,0])
dx=1
dy=1

xg=np.linspace(0,egx,egx)
yg=np.linspace(0,egy,egy)

#divv=(v1[margin+1:egy-margin+1,margin:egx-margin]-\
#      v1[margin-1:egy-margin-1,margin:egx-margin])/(2.0*dy) \
#    +(v2[margin:egy-margin,margin+1:egx-margin+1]-\
#      v2[margin:egy-margin,margin-1:egx-margin-1])/(2.0*dx) 


fig, ax = plt.subplots(figsize=(9, 6))
ax.contourf(np.log10(rho).T,levels=101,cmap='Greys')
#plt.ylim([950, 1100])

#
nproc=1
shocks=shockid(xg,yg,0.0,rho,v1,v2,0.0,b1,b2,0.0,pr,smthfac=2,nproc=nproc,avecyl=5)

#if nproc == 1:
plt.plot(shocks['slow'][:,1],shocks['slow'][:,0],color='r',marker='+',linestyle='',markersize=2.5)
plt.plot(shocks['fast'][:,1],shocks['fast'][:,0],color='b',marker='+',linestyle='',markersize=2.5)
plt.plot(shocks['int4'][:,1],shocks['int4'][:,0],color='g',marker='+',linestyle='',markersize=2.5)
"""if nproc > 1:
	for j in range (0,nproc):
		plt.plot(shocks['slow'][j][:,1],shocks['slow'][j][:,0],color='r',marker='+',linestyle='',markersize=2.5)
		plt.plot(shocks['fast'][j][:,1],shocks['fast'][j][:,0],color='b',marker='+',linestyle='',markersize=2.5)
		#plt.plot(shocks['int4'][j][:,1],shocks['int4'][j][:,0],color='g',marker='+',linestyle='',markersize=2.5)
"""
savename='testLare_plot.png'
plt.savefig(savename)

#Create a h5 file of the shock data
hf = h5py.File(''.join((fname,'shocks.h5')), 'w')
for f in shocks:
	hf.create_dataset(f, data=shocks[f])
hf.close()