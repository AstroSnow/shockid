# -*- coding: utf-8 -*-
"""
PYTHON code for shockID
To Do:
1) everything
"""
import numpy as np
import time

def shockid(gridx,gridy,gridz,rog,vxg,vyg,vzg,bxg,byg,bzg,prg,ndim=2,smthfac=0,nproc=1,convl=0.0001,avecyl=2):
#	import numpy as np
	#Parameters for shock limits
#	convl=0.0001 #Convergence threshold
#	avecyl=2 #Cylinder to average over
#	smthfac=0 #smoothing factor (1=no smoothing)
	shocktol=0.05 #tolerence for the shock transitions
	bulkspeed=0 # Use the bulk sound and alfven speeds or individual
	#species='plasma' # plasma or neutral
	species='plasma' # plasma or neutral
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
	#print('Removing non-local maximum candidates')
	time1=time.perf_counter()
	if nproc == 1:
		[col,row]=removeNonMax(col,row,zrow,ds['ro'],gradrox,gradroy,gradroz,gradmag,divv,ndim,2)
	if nproc > 1:
		[col,row]=removeNonMaxPar(col,row,zrow,ds['ro'],gradrox,gradroy,gradroz,gradmag,divv,ndim,2,nproc)
	time2=time.perf_counter()
	
	print('Removing non-local max took ',time2-time1,' on ',nproc,' cores')
	#return(col,row)	
	
	#print(np.size(col))
	#stop
	#Allocate arrays for the shock locations
	print('Begin loop for n=',np.size(col))
	time1=time.perf_counter()
	if nproc == 1:
		#Series loop. Can be removed at some point
		shocks={}
		fast=[0,0]
		slow=[0,0]
		int1=[0,0]
		int2=[0,0]
		int3=[0,0]
		int4=[0,0]
	
		for i in range(0,np.size(col)):
			#Calcuate data along the LOS
			normarr=getNormVals(col[i],row[i],zrow,ds,gradrox,gradroy,gradroz,gradmag,divv,ndim,avecyl,gcalc=False)
			#get indecies of pre and post shock states
			[ipre,ipos]=prepostIndex(normarr['ro'],avecyl)
			#vsa=getShockFrame(normarr['ro'][ipos],normarr['ro'][ipre],normarr['vperp'][ipos],normarr['vperp'][ipre])
			vsa=getShockFrame(normarr['ro'][ipos],normarr['ro'][ipre],normarr['vperp'][ipos],normarr['vperp'][ipre],
						normarr['vpar'][ipos],normarr['vpar'][ipre],normarr['bpar'][ipos],normarr['bpar'][ipre],normarr['bperp'][ipos],normarr['bperp'][ipre])
			speeds=getWaveSpeeds(normarr['ro'],normarr['pr'],normarr['bx'],normarr['by'],normarr['bz'],normarr['bperp'], normarr['ang'])
			#Put velocity in shock frame
			vfpos=normarr['vperp'][ipos]+vsa
			vfpre=normarr['vperp'][ipre]+vsa
			
			vslowpre =np.abs(vfpre/np.sqrt(speeds['vslow2'][ipre]))
			valfpre  =np.abs(vfpre/np.sqrt(speeds['vap2'][ipre]))
			vfastpre =np.abs(vfpre/np.sqrt(speeds['vfast2'][ipre]))
			prestate=getState(vslowpre,valfpre,vfastpre)
			
			vslowpos =np.abs(vfpos/np.sqrt(speeds['vslow2'][ipos]))
			valfpos  =np.abs(vfpos/np.sqrt(speeds['vap2'][ipos]))
			vfastpos =np.abs(vfpos/np.sqrt(speeds['vfast2'][ipos]))
			posstate=getState(vslowpos,valfpos,vfastpos)
		
			#print(col[i],row[i])
			#Get the transitions
			if (prestate == 1) and (posstate==2):
				#Fast shocks
				fast=np.vstack((fast,[col[i],row[i]]))
				#print('fast shock')
			if (prestate == 3) and (posstate==4):
				#Fast shocks
				slow=np.vstack((slow,[col[i],row[i]]))
			if (prestate == 1) and (posstate==3):
				int1=np.vstack((int1,[col[i],row[i]]))
			if (prestate == 1) and (posstate==4):
				int2=np.vstack((int2,[col[i],row[i]]))
			if (prestate == 2) and (posstate==3):
				int3=np.vstack((int3,[col[i],row[i]]))
			if (prestate == 2) and (posstate==4):
				int4=np.vstack((int4,[col[i],row[i]]))
							
		shocks['slow']=slow
		shocks['fast']=fast
		shocks['int1']=int1
		shocks['int2']=int2
		shocks['int3']=int3
		shocks['int4']=int4
		
	if nproc > 1:
		#Parallel loop
		import multiprocessing as mp
		pool=mp.get_context('fork').Pool(nproc)
		print('using ',nproc,' cores')
		pitrange=np.zeros((nproc,2))
		st=int(0)
		for i in range(0,nproc-1):
			pitrange[i,0]=int(st)
			pitrange[i,1]=int(st+np.floor(np.size(col)/nproc))
			st=pitrange[i,1]+1
		pitrange[nproc-1,0]=int(st)
		pitrange[nproc-1,1]=int(np.size(col)-1)
		print(pitrange[:,0],pitrange[:,1])
		print('Each processor doing ',int(pitrange[1,0]-pitrange[0,0]),' elements')
		#sol=shockClassLoop(pitrange[0,0],pitrange[0,1],col,row,zrow,ds,gradrox,gradroy,gradroz,gradmag,divv,ndim,avecyl)
		sol=pool.starmap(shockClassLoop,[(pitrange[j,0],pitrange[j,1],col,row,zrow,ds,gradrox,gradroy,gradroz,gradmag,divv,ndim,avecyl) for j in range(0,nproc)])
		pool.close()
		pool.join()
		#NEED TO REJOIN EVERYTHING
		shocks={}
		shocks['slow']=np.vstack([sol[j]['slow'] for j in range(0,nproc)])
		shocks['fast']=np.vstack([sol[j]['fast'] for j in range(0,nproc)])
		shocks['int1']=np.vstack([sol[j]['int1'] for j in range(0,nproc)])
		shocks['int2']=np.vstack([sol[j]['int2'] for j in range(0,nproc)])
		shocks['int3']=np.vstack([sol[j]['int3'] for j in range(0,nproc)])
		shocks['int4']=np.vstack([sol[j]['int4'] for j in range(0,nproc)])

		shocks['slowmach']=np.vstack([sol[j]['slowmach'] for j in range(0,nproc)])
		shocks['fastmach']=np.vstack([sol[j]['fastmach'] for j in range(0,nproc)])
		shocks['int1mach']=np.vstack([sol[j]['int1mach'] for j in range(0,nproc)])
		shocks['int2mach']=np.vstack([sol[j]['int2mach'] for j in range(0,nproc)])
		shocks['int3mach']=np.vstack([sol[j]['int3mach'] for j in range(0,nproc)])
		shocks['int4mach']=np.vstack([sol[j]['int4mach'] for j in range(0,nproc)])
		
	time2=time.perf_counter()
	
	print('Finding shocks took ',time2-time1,' on ',nproc,' cores')
	return(shocks)
###############################################################################
def shockClassLoop(istart,iend,col,row,zrow,ds,gradrox,gradroy,gradroz,gradmag,divv,ndim,avecyl):
	#Loop for identifying shocks for use in parallel
	shocks={}
	fast=[0,0]
	slow=[0,0]
	int1=[0,0]
	int2=[0,0]
	int3=[0,0]
	int4=[0,0]
	
	#arrays for storing mach numbers
	slowmach=[0,0]
	fastmach=[0,0]
	int1mach=[0,0]
	int2mach=[0,0]
	int3mach=[0,0]
	int4mach=[0,0]
	for i in range(int(istart),int(iend)+1):
		#print(i)
		#Calcuate data along the LOS
		normarr=getNormVals(col[i],row[i],zrow,ds,gradrox,gradroy,gradroz,gradmag,divv,ndim,avecyl,gcalc=False)
		#get indecies of pre and post shock states
		[ipre,ipos]=prepostIndex(normarr['ro'],avecyl)
#		vsa=getShockFrame(normarr['ro'][ipos],normarr['ro'][ipre],normarr['vperp'][ipos],normarr['vperp'][ipre])
		vsa=getShockFrame(normarr['ro'][ipos],normarr['ro'][ipre],normarr['vperp'][ipos],normarr['vperp'][ipre],
					normarr['vpar'][ipos],normarr['vpar'][ipre],normarr['bpar'][ipos],normarr['bpar'][ipre],normarr['bperp'][ipos],normarr['bperp'][ipre])
		speeds=getWaveSpeeds(normarr['ro'],normarr['pr'],normarr['bx'],normarr['by'],normarr['bz'],normarr['bperp'], normarr['ang'])
		#Put velocity in shock frame
		vfpos=normarr['vperp'][ipos]+vsa
		vfpre=normarr['vperp'][ipre]+vsa
		
		vslowpre =np.abs(vfpre/np.sqrt(speeds['vslow2'][ipre]))
		valfpre  =np.abs(vfpre/np.sqrt(speeds['vap2'][ipre]))
		vfastpre =np.abs(vfpre/np.sqrt(speeds['vfast2'][ipre]))
		prestate=getState(vslowpre,valfpre,vfastpre)
		
		vslowpos =np.abs(vfpos/np.sqrt(speeds['vslow2'][ipos]))
		valfpos  =np.abs(vfpos/np.sqrt(speeds['vap2'][ipos]))
		vfastpos =np.abs(vfpos/np.sqrt(speeds['vfast2'][ipos]))
		posstate=getState(vslowpos,valfpos,vfastpos)
	
		#print(col[i],row[i])
		#Get the transitions
		if (prestate == 1) and (posstate==2):
			#Fast shocks
			fast=np.vstack((fast,[col[i],row[i]]))
			fastmach=np.vstack((fastmach,[valfpre,valfpos]))
			#print('fast shock')
		if (prestate == 3) and (posstate==4):
			#Fast shocks
			slow=np.vstack((slow,[col[i],row[i]]))
			slowmach=np.vstack((slowmach,[valfpre,valfpos]))
		if (prestate == 1) and (posstate==3):
			int1=np.vstack((int1,[col[i],row[i]]))
			int1mach=np.vstack((int1mach,[valfpre,valfpos]))
		if (prestate == 1) and (posstate==4):
			int2=np.vstack((int2,[col[i],row[i]]))
			int2mach=np.vstack((int2mach,[valfpre,valfpos]))
		if (prestate == 2) and (posstate==3):
			int3=np.vstack((int3,[col[i],row[i]]))
			int3mach=np.vstack((int3mach,[valfpre,valfpos]))
		if (prestate == 2) and (posstate==4):
			int4=np.vstack((int4,[col[i],row[i]]))
			int4mach=np.vstack((int4mach,[valfpre,valfpos]))
			
	shocks['slow']=slow
	shocks['fast']=fast
	shocks['int1']=int1
	shocks['int2']=int2
	shocks['int3']=int3
	shocks['int4']=int4
	
	shocks['slowmach']=slowmach
	shocks['fastmach']=fastmach
	shocks['int1mach']=int1mach
	shocks['int2mach']=int2mach
	shocks['int3mach']=int3mach
	shocks['int4mach']=int4mach
	return(shocks)

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
			ds['bx']=gaussian_filter(bx,smthfac)
			ds['by']=gaussian_filter(by,smthfac)
			ds['bz']=gaussian_filter(bz,smthfac)
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
def getNormVals(col,row,zrow,var,gradx,grady,gradz,gradmag,divv,ndim,avecyl,gcalc=False):
	#Step ii: find shock normal based on density gradient
	if ndim == 2:
	    normx=gradx[col,row]/gradmag[col,row]
	    normy=grady[col,row]/gradmag[col,row]
	"""if ndim == 3:
	    normx=gradx(col(i),row(i),zrow(i))/gradmag(col(i),row(i),zrow(i))
	    normy=grady(col(i),row(i),zrow(i))/gradmag(col(i),row(i),zrow(i))
	    normz=gradz(col(i),row(i),zrow(i))/gradmag(col(i),row(i),zrow(i))
	"""		

	tempx=np.linspace(col-avecyl,col+avecyl,2*avecyl+1)
	tempy=np.linspace(row-avecyl,row+avecyl,2*avecyl+1)
	tempx2=np.linspace(-avecyl,avecyl,2*avecyl+1)
	tempy2=np.linspace(-avecyl,avecyl,2*avecyl+1)

#		if ndim == 3:
#			tempz=indgen(2*2+1)-2+zrow(i)
#			tempz2=indgen(2*2+1)-2
	
	#use periodic BC to fix negative values
	for ii in range(0,2*avecyl+1):
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
	for ii in range(0,2*avecyl+1):
		if ndim == 2:
			if tempx[ii] > np.size(divv[:,0])-1:
				tempx[ii]=tempx[ii]-np.size(divv[:,0])
			if tempy[ii] > np.size(divv[0,:])-1:
				tempy[ii]=tempy[ii]-np.size(divv[0,:])
				
	#Interpolate the values normal to the shock NOTE: THES ARE NOT PARALLEL TO SHOCK YET!
	normvals={}
	if ndim == 2:
		if gcalc == True:	
			ronorm=shocknormvals(var,tempx,tempy,tempx2,tempy2,2,normx,normy,avecyl)   
			return(ronorm)
		if gcalc == False:
			ronorm=shocknormvals(var['ro'],tempx,tempy,tempx2,tempy2,2,normx,normy,avecyl)            
			normvals['ro']=ronorm
			vxnorm=shocknormvals(var['vx'],tempx,tempy,tempx2,tempy2,2,normx,normy,avecyl)
			normvals['vx']=vxnorm
			vynorm=shocknormvals(var['vy'],tempx,tempy,tempx2,tempy2,2,normx,normy,avecyl)
			normvals['vy']=vynorm
			prnorm=shocknormvals(var['pr'],tempx,tempy,tempx2,tempy2,2,normx,normy,avecyl)
			normvals['pr']=prnorm
			bxnorm=shocknormvals(var['bx'],tempx,tempy,tempx2,tempy2,2,normx,normy,avecyl)
			normvals['bx']=bxnorm
			bynorm=shocknormvals(var['by'],tempx,tempy,tempx2,tempy2,2,normx,normy,avecyl)
			normvals['by']=bynorm
			bznorm=var['bz']
			if np.size(var['bz']) >1:
				bznorm=shocknormvals(var['bz'],tempx,tempy,tempx2,tempy2,2,normx,normy,avecyl)
			normvals['bz']=bznorm
			normvals['normx']=normx
			normvals['normy']=normy
			normvals['perpx']=(-normy/normx)/np.sqrt(normy**2/normx**2+1.0)
			normvals['perpy']=1.0/np.sqrt(normy**2/normx**2+1.0)
			normvals['bperp']=normvals['normx']*normvals['bx']+normvals['normy']*normvals['by']    
			normvals['bpar']=normvals['perpx']*normvals['bx']+normvals['perpy']*normvals['by']
			normvals['vperp']=normvals['normx']*normvals['vx']+normvals['normy']*normvals['vy']    
			normvals['vpar']=normvals['perpx']*normvals['vx']+normvals['perpy']*normvals['vy']
			normvals['ang']=np.arctan(normvals['bpar']/normvals['bperp']) 
			return(normvals)
	
################################################################
def removeNonMax(col,row,zrow,ro,gradx,grady,gradz,gradmag,divv,ndim,avecyl):
	import numpy as np
	#Remove cells that are not a local maximum
	print('Removing non-local maximum density gradients from candidate cells n=',np.size(col))

	col2=[]
	row2=[]
	zrow2=[]

	for i in range(0,np.size(col)):
		ronorm=getNormVals(col[i],row[i],zrow,ro,gradx,grady,gradz,gradmag,divv,ndim,avecyl,gcalc=True)

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
################################################################
def removeNonMaxPar(col,row,zrow,ro,gradx,grady,gradz,gradmag,divv,ndim,avecyl,nproc):
	import numpy as np
	import multiprocessing as mp
	#Remove cells that are not a local maximum in parallel
	print('Removing non-local maximum density gradients from candidate cells n=',np.size(col))

#	col2=[]
#	row2=[]
	zrow2=[]

#	solarr=np.zeros((np.size(col),nproc))	
	
	#mp.set_start_method('fork')
	pool=mp.get_context('fork').Pool(nproc)
#	pool=mp.Pool(nproc)
	print('using ',nproc,' cores')
#	itrange=list(range(0,np.size(col)))
	pitrange=np.zeros((nproc,2))
	st=int(0)
	for i in range(0,nproc-1):
		pitrange[i,0]=int(st)
		pitrange[i,1]=int(st+np.floor(np.size(col)/nproc))
		st=pitrange[i,1]+1
	pitrange[nproc-1,0]=int(st)
	pitrange[nproc-1,1]=int(np.size(col)-1)
	print('Each processor doing ',pitrange[1,0]-pitrange[0,0],' elements')

	sol=pool.starmap(removeNonMaxParLoop,[(pitrange[j,0],pitrange[j,1],col,row,zrow,ro,gradx,grady,gradz,gradmag,divv,ndim,avecyl) for j in range(0,nproc)])
	pool.close()
	pool.join()

	sol=np.asarray(sol)
	sol=np.sum(sol,axis=0)
	#print(np.size(sol))
	row=row[np.argwhere(sol == 1)]
	col=col[np.argwhere(sol == 1)]
	#row=row2
	#col=col2
	if ndim==3:
		zrow=zrow2
	#print(np.size(row),np.size(col))	
	#print(row.reshape(np.size(row)))
	#stop
	return(col.reshape(np.size(row)),row.reshape(np.size(row)))

###############################################################################
def removeNonMaxParLoop(istart,iend,col,row,zrow,ro,gradx,grady,gradz,gradmag,divv,ndim,avecyl):
	sol=np.zeros(np.size(col))
#	row2=[]
#	col2=[]
#	print(istart,iend)
	for i in range(int(istart),int(iend)+1):
#		print(i)
		ronorm=getNormVals(col[i],row[i],zrow,ro,gradx,grady,gradz,gradmag,divv,ndim,avecyl,True)
		b=np.argmax(np.abs(np.gradient(ronorm)))
#		if strat eq 1 then a=max(abs(deriv(ronorm-mrotemp)),b)
#print(np.gradient(ronorm))
		if b == 2:
			sol[i]=1
	return(sol)

###############################################################################
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

###############################################################################
def prepostIndex(ro,avecyl):
	if ro[0] < ro[2*avecyl]: 
		ipre=np.argmin(ro[0:avecyl-1])
		ipos=np.argmax(ro[avecyl+1:avecyl*2])
		ipos=ipos+avecyl+1
	else:
		ipos=np.argmax(ro[0:avecyl-1])
		ipre=np.argmin(ro[avecyl+1:avecyl*2])
		ipre=ipre+avecyl+1
		
	return(ipre,ipos)

###############################################################################
def getShockFrame(ropos,ropre,vperppos,vperppre,vparpos,vparpre,bparpos,bparpre,bperppos,bperppre):
	#Just assuming mass conservation. There are better ways to do this.
	vsa=(ropos*vperppos-ropre*vperppre)/(ropre-ropos)
	#Electric field
#	vsa=(vperppos*bparpos-vperppre*bparpre-vparpos*bperppos+vparpre*bperppre)/(bparpos-bparpre)
	return(vsa)

###############################################################################
def getWaveSpeeds(ro,pr,bx,by,bz,bperp,ang):
	cs2=5.0/3.0*pr/ro
	va2=(bx**2+by**2)/ro
	if np.size(bz) > 1:
		va2=(bx**2+by**2+bz**2)/ro
	vap2=bperp**2/ro
	vslow2=0.5*(cs2+va2-np.sqrt((cs2+va2)**2 - 4.0*va2*cs2*(np.cos(ang))**2))
	vfast2=0.5*(cs2+va2+np.sqrt((cs2+va2)**2 - 4.0*va2*cs2*(np.cos(ang))**2))
	speeds={'cs2':cs2,'va2':va2,'vap2':vap2,'vslow2':vslow2,'vfast2':vfast2}
	return(speeds)
	
###############################################################################
def getState(vslow,valp,vfast):
	state=4
	if vfast >= 1:
		state=1
	if valp >=1 and vfast < 1:
		state=2
	if vslow >=1 and valp < 1:
		state=3
	return(state)