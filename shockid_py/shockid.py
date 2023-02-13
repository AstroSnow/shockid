# -*- coding: utf-8 -*-
"""
PYTHON code for shockID
To Do:
1) everything
"""
import numpy as np

def shockid(gridx,gridy,gridz,rog,vxg,vyg,vzg,bxg,byg,bzg,prg,ndim=2):
#	import numpy as np
	#Parameters for shock limits
	convl=0.01 #Convergence threshold
	avecyl=4 #Cylinder to average over
	smthfac=2 #smoothing factor (1=no smoothing)
	shocktol=0.05 #tolerence for the shock transitions
	bulkspeed=0 # Use the bulk sound and alfven speeds or individual
	#species='plasma' # plasma or neutral
	species='neutral' # plasma or neutral
	data_subset=0
	
	#Stratified mean density profiles
	strat=0
	
	#Define simulation grid
	print('Removing ghost cells')
	margin=2

	#Define the end of the grid	
	egx=np.size(gridx)-1
	egy=np.size(gridy)-1
	egz=0.0
	if (ndim == 3):
		egz=np.size(gridz)-1

	#Define the grid
	x=gridx[margin:egx-margin]
	y=gridy[margin:egy-margin]
	if (ndim == 3):
		z=gridz[margin:egz-margin]

	dx=x[1]-x[0]
	dy=y[1]-y[0]
	dz=0.0#z[1]-z[0]
	#Smooth the data (if needed)
	ds=smoothdata(rog,vxg,vyg,vzg,bxg,byg,bzg,prg,ndim,species,margin,smthfac)

	#Divergence of velocity field
	print('Calculating velocity divergence')
	divv=divergence(ds,2,margin,egx,egy,egz,dx,dy,dz)
	
	#Calculate the density gradients
	print('Calculating density gradients')
	gradrox=np.gradient(ds['ro'],axis=1)
	gradroy=np.gradient(ds['ro'],axis=0)
	gradroz=0.0
	if (ndim == 3):
		gradroz=np.gradient(ds['ro'],axis=3)
	gradmag=np.sqrt(gradrox**2+gradroy**2+gradroz**2)

	#identify candidate cells based on divv
	print('Identify candidate cells')
	if (ndim ==2):
		temp=np.argwhere(divv < -convl)
		row=temp[:,1]
		col=temp[:,0]
		zrow=0
		#[col,row] = np.argwhere(divv < -convl)
	if (ndim ==3):
		[col,row,zrow] = np.argwhere(divv < -convl)
		
	

	"""#Filter the data 
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
	"""		
	#Remove non-local maximun gradient
	[col,row]=removeNonMax(col,row,zrow,ds['ro'],gradrox,gradroy,gradroz,gradmag,divv,ndim,2)
	
	return(col,row)	
	"""
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
"""	
#####################################################################################################
#Data smoothing routine
def smoothdata(ro,vx,vy,vz,bx,by,bz,pr,ndim,species,margin,smthfac):
	from scipy.ndimage import gaussian_filter
	print('Smoothing isnt overly tested and is simulation specific')
	#Appy some smoothing	
	#Gaussian smoothing?
	if (ndim == 2):
		ro=gaussian_filter(ro,smthfac)
		vx=gaussian_filter(vx,smthfac)
		vy=gaussian_filter(vy,smthfac)
		pr=gaussian_filter(pr,smthfac)	
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
			ds[bx]=gaussian_filter(bx,smthfac)
			ds[by]=gaussian_filter(by,smthfac)
			ds[bz]=gaussian_filter(bz,smthfac)
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
		divv=(ds['vy'][margin+1:egy-margin+1,margin:egx-margin]-\
		ds['vy'][margin-1:egy-margin-1,margin:egx-margin])/(2.0*dy) \
	    +(ds['vx'][margin:egy-margin,margin+1:egx-margin+1]-\
		ds['vx'][margin:egy-margin,margin-1:egx-margin-1])/(2.0*dx) 
#	if (ndim == 3):
#	    divv=(ds['vx'](margin+1:egx-margin+1,margin:egy-margin,margin:egz-margin)-\
#		ds['vx'](margin-1:egx-margin-1,margin:egy-margin,margin:egz-margin))/(2.0*dx) \
#	    +(ds['vy'](margin:egx-margin,margin+1:egy-margin+1,margin:egz-margin)-\
#		ds['vy'](margin:egx-margin,margin-1:egy-margin-1,margin:egz-margin))/(2.0*dy) \
#	    +(ds['vz'](margin:egx-margin,margin:egy-margin,margin+1:egz-margin+1)-\
#		ds['vz'](margin:egx-margin,margin:egy-margin,margin-1:egz-margin-1))/(2.0*dz)
	return(divv)
################################################################
def removeNonMax(col,row,zrow,ro,gradx,grady,gradz,gradmag,divv,ndim,avecyl):
	import numpy as np
	#Remove cells that are not a local maximum
	print('Removing non-local maximum density gradients from candidate cells n=',np.size(col))

	col2=[]
	row2=[]
	zrow2=[]

	for i in range(0,np.size(col)):
		#Step ii: find shock normal based on density gradient
		if ndim == 2:
		    normx=gradx[col[i],row[i]]/gradmag[col[i],row[i]]
		    normy=grady[col[i],row[i]]/gradmag[col[i],row[i]]
		
		"""if ndim == 3:
		    normx=gradx(col(i),row(i),zrow(i))/gradmag(col(i),row(i),zrow(i))
		    normy=grady(col(i),row(i),zrow(i))/gradmag(col(i),row(i),zrow(i))
		    normz=gradz(col(i),row(i),zrow(i))/gradmag(col(i),row(i),zrow(i))
		"""		

		tempx=np.linspace(col[i]-2,col[i]+2,5)
		tempy=np.linspace(row[i]-2,row[i]+2,5)
		tempx2=np.linspace(-2,2,5)
		tempy2=np.linspace(-2,2,5)

#		if ndim == 3:
#			tempz=indgen(2*2+1)-2+zrow(i)
#			tempz2=indgen(2*2+1)-2
		
		#use periodic BC to fix negative values
		for ii in range(0,5):
			if ndim == 2:
				if tempx[ii] < 0:
					tempx[ii]=np.size(divv[:,0])+tempx[ii]
				if tempy[ii] < 0:
					tempy[ii]=np.size(divv[0,:])+tempy[ii]
#		    if ndim == 3:
#		        if tempx(ii) lt 0 then tempx(ii)=N_elements(divv(*,0,0))+tempx(ii)
#		        if tempy(ii) lt 0 then tempy(ii)=N_elements(divv(0,*,0))+tempy(ii)
#		        if tempz(ii) lt 0 then tempz(ii)=N_elements(divv(0,0,*))+tempz(ii)
		    
		#use periodic BC to fix outside grid values
		for ii in range(0,5):
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
		    ronorm=shocknormvals(ro,tempx,tempy,tempx2,tempy2,2,normx,normy,avecyl)
		
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
		    

		#a=max(abs(deriv(ronorm)),b)
		b=np.argmax(np.abs(np.gradient(ronorm)))
#		if strat eq 1 then a=max(abs(deriv(ronorm-mrotemp)),b)
		#print(np.gradient(ronorm))
		if b == 2:
			row2.append(row[i])#=[row2,row[i]]
			col2.append(col[i])#=[col2,col[i]]
#		if ndim eq 3 then zrow2=[zrow2,zrow(i)]
		
	row=row2
	col=col2
	if ndim==3:
		zrow=zrow2
		
	return(col,row)

def shocknormvals(var,tempx,tempy,tempx2,tempy2,ndim,normx,normy,avecyl):
    #from scipy.interpolate import RegularGridInterpolator
    import scipy.interpolate as spint
    RGI = spint.RegularGridInterpolator
	#Calculate the values normal to the shock
    tempa=np.zeros([2*avecyl+1,2*avecyl+1])
    
    #populate the array (density)
    for ii in range(0,2*avecyl+1):
        for jj in range(0,2*avecyl+1):
            #print(ii,jj,tempx[ii].astype(int),tempy[jj].astype(int))
            tempa[ii,jj]=var[tempx[ii].astype(int),tempy[jj].astype(int)]

    #find the normal line
    normlinex=normx*np.linspace(-avecyl,avecyl,avecyl*2+1)#(indgen(2*avecyl+1)-avecyl)
    normliney=normy*np.linspace(-avecyl,avecyl,avecyl*2+1)#*(indgen(2*avecyl+1)-avecyl)   
    pnts=[np.linspace(-avecyl,avecyl,avecyl*2+1),np.linspace(-avecyl,avecyl,avecyl*2+1)]
								
    #Interpolate the values normal to the shock
    rgi=RGI(points=pnts,values=tempa)
    #var2=rgi([normlinex,normliney])
    var2=np.zeros(2*avecyl+1)
    for ii in range(0,2*avecyl+1):
        var2[ii]=rgi([normliney[ii],normlinex[ii]])
#	var2=interp2d(tempa,tempx2,tempy2,normlinex,normliney)
    #print(normliney)
    return(var2)
