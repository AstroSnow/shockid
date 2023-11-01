#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 31 07:20:28 2023

@author: ben
"""
import numpy as np
import h5py as h5
import matplotlib.pyplot as plt
from shockid import shockid

gamma=5.0/3.0

nx=100
ny=100
nz=100

dx=1.0/nx
dy=1.0/ny
dz=1.0/nz

xg=np.linspace(0,1,nx)
yg=np.linspace(0,1,ny)
zg=np.linspace(0,1,nz)

ro=np.zeros((nz,ny,nx))
vx=np.zeros((nz,ny,nx))
vy=np.zeros((nz,ny,nx))
vz=np.zeros((nz,ny,nx))
bx=np.zeros((nz,ny,nx))
by=np.zeros((nz,ny,nx))
bz=np.zeros((nz,ny,nx))
pr=np.zeros((nz,ny,nx))

M=2
r=M**2*(gamma+1.0)/(2.0+M**2*(gamma-1))

ro[:,:,:]=1
pr[:,:,:]=1.0/gamma
vx[:,:,:]=M

for i in range(0,nx):
	for j in range(0,ny):
		for k in range(0,nz):
			if xg[i] > 0.5:
				ro[k,j,i]=r
				vx[k,j,i]=M/r
				pr[k,j,i]=(1.0+gamma*M*M*(1.0-1.0/r))


margin=2
divv=np.zeros((nz,ny,nx))
divv[margin:nz-margin,margin:ny-margin,margin:nx-margin]=(vy[margin:nz-margin,margin+1:ny-margin+1,margin:nx-margin]-\
vy[margin:nz-margin,margin-1:ny-margin-1,margin:nx-margin])/(2.0*dy) \
   +(vx[margin:nz-margin,margin:ny-margin,margin+1:nx-margin+1]-\
vx[margin:nz-margin,margin:ny-margin,margin-1:nx-margin-1])/(2.0*dz) \
   +(vz[margin+1:nz-margin+1,margin:ny-margin,margin:nx-margin]-\
vz[margin-1:nz-margin-1,margin:ny-margin,margin:nx-margin])/(2.0*dx)

gradrox=np.gradient(ro,axis=2)
gradroy=np.gradient(ro,axis=1)
gradroz=np.gradient(ro,axis=0)

xs=0
xe=-1
ys=0
ye=-1
zs=0
ze=-1

shocks=shockid(xg[xs:xe],yg[ys:ye],zg[zs:ze],ro[zs:ze,ys:ye,xs:xe],
			   vx[zs:ze,ys:ye,xs:xe],vy[zs:ze,ys:ye,xs:xe],vz[zs:ze,ys:ye,xs:xe],
			   bx[zs:ze,ys:ye,xs:xe],by[zs:ze,ys:ye,xs:xe],bz[zs:ze,ys:ye,xs:xe],
			   pr[zs:ze,ys:ye,xs:xe],ndim=3,smthfac=0,nproc=1,convl=0.000,avecyl=5)