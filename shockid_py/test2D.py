# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

from pipreadmods import pipread
import matplotlib.pyplot as plt
import numpy as np
from shockid import shockid

fname='../MHD_2D_rec/'
fname='../../shockstab/Data/'

ds=pipread(fname,10)
    
v1=ds['vx_p']
v2=ds['vy_p']

margin=2
egx=np.size(v1[0,:])
egy=np.size(v1[:,0])
dx=ds['xgrid'][1]-ds['xgrid'][0]
dy=ds['ygrid'][1]-ds['ygrid'][0]

divv=(v2[margin+1:egy-margin+1,margin:egx-margin]-\
      v2[margin-1:egy-margin-1,margin:egx-margin])/(2.0*dy) \
    +(v1[margin:egy-margin,margin+1:egx-margin+1]-\
      v1[margin:egy-margin,margin-1:egx-margin-1])/(2.0*dx) 

#plt.contourf(divv)
#plt.contourf(divv[:,950:1150],levels=np.linspace(-0.055,-0.00001,100))
plt.contourf(divv,levels=101,vmax=-0.01)

#[col,row]=shockid(ds['xgrid'],ds['ygrid'],0.0,ds['ro_p'],ds['vx_p'],ds['vy_p'],0.0,ds['bx'],ds['by'],ds['bz'],ds['pr_p'])
#plt.plot(row,col,color='r',linestyle='',marker='.',markersize=2.8)

shocks=shockid(ds['xgrid'],ds['ygrid'],0.0,ds['ro_p'],ds['vx_p'],ds['vy_p'],0.0,ds['bx'],ds['by'],ds['bz'],ds['pr_p'],smthfac=2)
plt.plot(shocks['slow'][:,1],shocks['slow'][:,0],color='r',linestyle='',marker='.',markersize=2.8)
