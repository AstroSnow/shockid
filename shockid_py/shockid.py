# -*- coding: utf-8 -*-
"""
PYTHON code for shockID
To Do:
1) everything
"""

def shockid(gridx,gridy,gridz,rog,vxg,vyg,vzg,bxg,byg,bzg,prg,ndim=2):
	import numpy as np
	#Parameters for shock limits
	convl=0.01 #Convergence threshold
	avecyl=4 #Cylinder to average over
	smthfac=1 #smoothing factor (1=no smoothing)
	shocktol=0.05 #tolerence for the shock transitions
	bulkspeed=0 # Use the bulk sound and alfven speeds or individual
	#species='plasma' # plasma or neutral
	species='neutral' # plasma or neutral
	data_subset=1
	
	#Stratified mean density profiles
	strat=1
	
	#Define simulation grid
	print('Removing ghost cells')
	margin=2

	#Define the end of the grid	
	egx=np.size(gridx)-1
	egy=np.size(gridy)-1
	if (ndim == 3):
		egz=np.size(gridz)-1

	#Define the grid
	x=gridx[margin:egx-margin]
	y=gridy[margin:egy-margin]
	if (ndim == 3):
		z=gridz[margin:egz-margin]

	#Smooth the data (if needed)
	ds=smoothdata(rog,vxg,vyg,vzg,bxg,byg,bzg,prg,ndim,species,margin,smthfac)

	#Divergence of velocity field
	print('Calculating velocity divergence')
	divv=divergence(ds,ndim=2)
	
	#Calculate the density gradients
	print('Calculating density gradients')
	gradrox=np.gradient(ds['ro'],axis=1)
	gradroy=np.gradient(ds['ro'],axis=2)
	gradroz=0.0
	if (ndim == 3):
		gradroz=np.gradient(ds['ro'],axis=3)
	gradmag=np.sqrt(gradrox**2+gradroy**2+gradroz**2)

	#identify candidate cells based on divv
	print('Identify candidate cells')
	if (ndim ==2):
		[col,row] = np.argwhere(divv < -convl)
	if (ndim ==3):
		[col,row,zrow] = np.argwhere(divv < -convl)
		
		
	#Filter the data 
	if data_subset == 1:
		print('Subsetting data')
		xminf=0.1#x(95);0.1
		xmaxf=0.9#x(105);0.9
		yminf=0.1
		ymaxf=0.9
		zminf=z[95]#0.45
		zmaxf=z[105]#0.55

		col2=[]
		row2=[]
		zrow2=[]

		if ndim == 3:
			xt=np.where((x[col] >= xminf) & (x[col] <= xmaxf) \
				&  (y[row] >= yminf) & (y[row] <= ymaxf) \
				&  (z[zrow] >= zminf) & (z[zrow] <= zmaxf))
			col2=col[xt]
			row2=row[xt]
			zrow2=zrow[xt]
			
		if ndim == 2:
			xt=np.where((x[col] >= xminf) & (x[col] <= xmaxf) \
				&  (y[row] >= yminf) & (y[row] <= ymaxf))
			col2=col[xt]
			row2=row[xt]
			
		col=col2
		row=row2
		if ndim == 3:
			zrow=zrow2
			
	#Remove non-local maximun gradient
	[col,row,zrow]=removeNonMax()
	
	#Allocate arrays for the shock locations
	#SHOULD BE A DICTIONARY
	fastx=[]
	fasty=[]
	fastz=[]
	fastm=[]

	slowx=[]
	slowy=[]
	slowz=[]
	slowm=[]

	int1x=[]
	int1y=[]
	int1z=[]

	int2x=[]
	int2y=[]
	int2z=[]

	int3x=[]
	int3y=[]
	int3z=[]

	int4x=[]
	int4y=[]
	int4z=[]

	print('Begin loop for n=',np.size(col))
	
#####################################################################################################
#Data smoothing routine
def smoothdata(ro,vx,vy,vz,bx,by,bz,pr,ndim,species,margin,smthfac):
	print('Smoothing doesnt work')
	print('Fit it yourself or do without')
	#Appy some smoothing	
	#Gaussian smoothing?
	if (ndim == 2):
#		ro=smooth(rog(margin:egx-margin,margin:egy-margin),smthfac)
#		vx=smooth(vxg(margin:egx-margin,margin:egy-margin),smthfac)
#		vy=smooth(vyg(margin:egx-margin,margin:egy-margin),smthfac)
#		vz=smooth(vzg(margin:egx-margin,margin:egy-margin),smthfac)
#		pr=smooth(prg(margin:egx-margin,margin:egy-margin),smthfac)
		ds={'ro':ro,'vx':vx,'vy':vy,'vz':vz,'pr':pr}
		if species == 'plasma':
#			bx=smooth(bxg(margin:egx-margin,margin:egy-margin),smthfac)
#			by=smooth(byg(margin:egx-margin,margin:egy-margin),smthfac)
#			bz=smooth(bzg(margin:egx-margin,margin:egy-margin),smthfac)
			ds[bx]=bx
			ds[by]=by
			ds[bz]=bz
	if (ndim == 3):
#		ro=smooth(rog(margin:egx-margin,margin:egy-margin,margin:egz-margin),smthfac)
#		vx=smooth(vxg(margin:egx-margin,margin:egy-margin,margin:egz-margin),smthfac)
#		vy=smooth(vyg(margin:egx-margin,margin:egy-margin,margin:egz-margin),smthfac)
#		vz=smooth(vzg(margin:egx-margin,margin:egy-margin,margin:egz-margin),smthfac)
#		pr=smooth(prg(margin:egx-margin,margin:egy-margin,margin:egz-margin),smthfac)
		ds={'ro':ro,'vx':vx,'vy':vy,'vz':vz,'pr':pr}
		if species == 'plasma':
#			bx=smooth(bxg(margin:egx-margin,margin:egy-margin,margin:egz-margin),smthfac)
#			by=smooth(byg(margin:egx-margin,margin:egy-margin,margin:egz-margin),smthfac)
#			bz=smooth(bzg(margin:egx-margin,margin:egy-margin,margin:egz-margin),smthfac)
			ds[bx]=bx
			ds[by]=by
			ds[bz]=bz
	return(ds)
################################################################
#Calculate the divergence of the velocity field
def divergence(ds,ndim,margin,egx,egy,egz,dx,dy,dz):
	if (ndim == 2):
		divv=(ds['vx'][margin+1:egx-margin+1,margin:egy-margin]-\
		ds['vx'][margin-1:egx-margin-1,margin:egy-margin])/(2.0*dx) \
	    +(ds['vy'][margin:egx-margin,margin+1:egy-margin+1]-\
		ds['vy'][margin:egx-margin,margin-1:egy-margin-1])/(2.0*dy) 
#	if (ndim == 3):
#	    divv=(ds['vx'](margin+1:egx-margin+1,margin:egy-margin,margin:egz-margin)-\
#		ds['vx'](margin-1:egx-margin-1,margin:egy-margin,margin:egz-margin))/(2.0*dx) \
#	    +(ds['vy'](margin:egx-margin,margin+1:egy-margin+1,margin:egz-margin)-\
#		ds['vy'](margin:egx-margin,margin-1:egy-margin-1,margin:egz-margin))/(2.0*dy) \
#	    +(ds['vz'](margin:egx-margin,margin:egy-margin,margin+1:egz-margin+1)-\
#		ds['vz'](margin:egx-margin,margin:egy-margin,margin-1:egz-margin-1))/(2.0*dz)
	return(divv)
################################################################
def removeNonMax(col,row,zrow):
	import numpy as np
	#Remove cells that are not a local maximum
	print('Removing non-local maximum density gradients from candidate cells n=',np.size(col))

	col2=[]
	row2=[]
	zrow2=[]

	for i in range(0,np.size(col)):
		#Step ii: find shock normal based on density gradient
		if ndim == 2:
		    normx=gradx(col(i),row(i))/gradmag(col(i),row(i))
		    normy=grady(col(i),row(i))/gradmag(col(i),row(i))
		
		if ndim == 3:
		    normx=gradx(col(i),row(i),zrow(i))/gradmag(col(i),row(i),zrow(i))
		    normy=grady(col(i),row(i),zrow(i))/gradmag(col(i),row(i),zrow(i))
		    normz=gradz(col(i),row(i),zrow(i))/gradmag(col(i),row(i),zrow(i))
				

		tempx=indgen(2*2+1)-2+col(i)
		tempy=indgen(2*2+1)-2+row(i)
		tempx2=indgen(2*2+1)-2
		tempy2=indgen(2*2+1)-2

		if ndim == 3:
			tempz=indgen(2*2+1)-2+zrow(i)
			tempz2=indgen(2*2+1)-2
		
		#use periodic BC to fix negative values
		for ii in range (0,2*2+1):
			if ndim == 2:
				if tempx(ii) < 0:
					tempx[ii]=np.size(divv[:,0])+tempx[ii]
				if tempy(ii) < 0:
					tempy[ii]=np.size(divv[0,:])+tempy[ii]
#		    if ndim == 3:
#		        if tempx(ii) lt 0 then tempx(ii)=N_elements(divv(*,0,0))+tempx(ii)
#		        if tempy(ii) lt 0 then tempy(ii)=N_elements(divv(0,*,0))+tempy(ii)
#		        if tempz(ii) lt 0 then tempz(ii)=N_elements(divv(0,0,*))+tempz(ii)
		    
		#use periodic BC to fix outside grid values
		for ii in range(0,2*2+1):
			if ndim == 2:
				if tempx[ii] > np.size(divv[:,0])-1:
					tempx[ii]=tempx[ii]-np.size(divv[:,0])
				if tempy[ii] > np.size(divv[0,:])-1:
					tempy[ii]=tempy[ii]-np.size(divv[0,:])
#			if ndim == 3:
#		        if tempx(ii) gt n_elements(divv(*,0,0))-1 then tempx(ii)=tempx(ii)-N_elements(divv(*,0,0))
#		        if tempy(ii) gt n_elements(divv(0,*,0))-1 then tempy(ii)=tempy(ii)-N_elements(divv(0,*,0))
#		        if tempz(ii) gt n_elements(divv(0,0,*))-1 then tempz(ii)=tempz(ii)-N_elements(divv(0,0,*))
		    
		#Interpolate the values normal to the shock NOTE: THES ARE NOT PARALLEL TO SHOCK YET!
		if ndim == 2:
		    ronorm=shocknormvals(ro,tempx,tempy,tempx2,tempy2,2,normx,normy)
		
#		if ndim == 3:
#		    tempxyz=dblarr(6,2*2+1)
#		    tempxyz(0,*)=tempx
#		    tempxyz(1,*)=tempy
#		    tempxyz(2,*)=tempz
#		    tempxyz(3,*)=tempx2
#		    tempxyz(4,*)=tempy2
#		    tempxyz(5,*)=tempz2
#		    ronorm=shocknormvals3D(ro,tempxyz,2,normx,normy,normz)
#		    if strat == 1:
#		        mrotemp=interpol(mro,z,tempz)
		    

		a=max(abs(deriv(ronorm)),b)
#		if strat eq 1 then a=max(abs(deriv(ronorm-mrotemp)),b)

		if b == 2:
			row2=[row2,row(i)]
			col2=[col2,col(i)]
#		if ndim eq 3 then zrow2=[zrow2,zrow(i)]
		
	row=row2
	col=col2
	if ndim==3:
		zrow=zrow2
