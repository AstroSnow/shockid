;Test of the 3d shockid function for Henrik

RESOLVE_ROUTINE, 'shockid3d_fun2';, /IS_FUNCTION

;Location of your data files
;fname='../PIPstab3D/Data/'
;fname='../../datadrive/PIPstab3D/PIP3D/'
;fname='../../datadrive/simdata/MHD_OZ_3D'
;fname='../Convection3D_288x288x168/sim_clean/c3d_0055900_T_32bits.h5'
fname='../Convection3D_576x576x328/sim_ambi/c3d_0056210_T_32bits.h5'
;Time step to read in
;t=10

;Read in the data
;rdmpi,ds,datapath=fname,time_step=t

file_id=H5F_OPEN(fname)
ro=H5D_read(H5D_open(file_id,'rho'))
pr=H5D_read(H5D_open(file_id,'pe'))

vx=H5D_read(H5D_open(file_id,'vx'))
vy=H5D_read(H5D_open(file_id,'vy'))
vz=H5D_read(H5D_open(file_id,'vz'))

mu0=4.0*!pi*1.e-7

bx=H5D_read(H5D_open(file_id,'bx'))/sqrt(mu0)
by=H5D_read(H5D_open(file_id,'by'))/sqrt(mu0)
bz=H5D_read(H5D_open(file_id,'bz'))/sqrt(mu0)

x=indgen(n_elements(ro(*,0,0)))
x=1.0*x/max(x)
y=indgen(n_elements(ro(0,*,0)))
y=1.0*y/max(y)
z=indgen(n_elements(ro(0,0,*)))
z=1.0*z/max(z)

betal=pr/((bx^2+by^2+bz^2)/2.0)
;stop

;Some info about the numerical grid
margin=1 ; number of ghost cells either side
egx=N_elements(x)-1
egy=N_elements(y)-1
egz=N_elements(z)-1
smthfac=1

;Run the shock detection routine 'shockid3d_fun2'
;Outputs: shocks - structure of all the detected shocks
;Inputs: ds.x - 1D array of your x data
;        ds.y - 1D array of your y data
;        ds.z - 1D array of your z data
;        ds.ro_p - 3D array of your density data
;        ds.vx_p - 3D array of your x-velocity data
;        ds.vy_p - 3D array of your y-velocity data
;        ds.vz_p - 3D array of your z-velocity data
;        ds.bx - 3D array of your x-magnetic field data
;        ds.by - 3D array of your y-magnetic field data
;        ds.bz - 3D array of your z-magnetic field data
;        ds.pr_p - 3D array of your gas pressure data
;        ndim=3 - this is should be left as is
shockid3d_fun2,shocks,x,y,z,ro,vx,vy,vz,bx,by,bz,pr,ndim=3

; Print out the number of shocks
print,'Slow',n_elements(shocks.slowx)
print,'Fast',n_elements(shocks.fastx)

;Save the shocks
save,shocks,filename='testelena.dat'

;Calculate divergence of velocity
;divv=(smooth(ds.vx_p(margin+1:egx-margin+1,margin:egy-margin,margin:egz-margin),smthfac)-$
;smooth(ds.vx_p(margin-1:egx-margin-1,margin:egy-margin,margin:egz-margin),smthfac))/(ds.x(10)-ds.x(8)) $
;    +(smooth(ds.vy_p(margin:egx-margin,margin+1:egy-margin+1,margin:egz-margin),smthfac)-$
;smooth(ds.vy_p(margin:egx-margin,margin-1:egy-margin-1,margin:egz-margin),smthfac))/(ds.y(10)-ds.y(8))$
;    +(smooth(ds.vz_p(margin:egx-margin,margin:egy-margin,margin+1:egz-margin+1),smthfac)-$
;smooth(ds.vz_p(margin:egx-margin,margin:egy-margin,margin-1:egz-margin-1),smthfac))/(ds.z(10)-ds.z(8))

;Some basic plotting routines

;;Plot a 2D slice of the data at z=100 grid cell
;fx=shocks.fastx(where((shocks.fastz le 110) and (shocks.fastz ge 90)))
;fy=shocks.fasty(where((shocks.fastz le 110) and (shocks.fastz ge 90)))
;sx=shocks.slowx(where((shocks.slowz le 110) and (shocks.slowz ge 90)))
;sy=shocks.slowy(where((shocks.slowz le 110) and (shocks.slowz ge 90)))

;Plot a 2D slice of the data at z=100 grid cell
fx=shocks.fasty(where((shocks.fastx le 110) and (shocks.fastx ge 90)))
fy=shocks.fastz(where((shocks.fastx le 110) and (shocks.fastx ge 90)))
sx=shocks.slowy(where((shocks.slowx le 110) and (shocks.slowx ge 90)))
sy=shocks.slowz(where((shocks.slowx le 110) and (shocks.slowx ge 90)))
i4x=shocks.int4y(where((shocks.int4x le 110) and (shocks.int4x ge 90)))
i4y=shocks.int4z(where((shocks.int4x le 110) and (shocks.int4x ge 90)))

;Use the density as a background
c=image(reform(ro(100,*,*)),rgb_tab=0)
p=plot(fx+margin,fy+margin,'b.',/overplot,sym_thick=2) ;The 6 is there because of the ghost cells in my code that the shockid routine assumes. We can change this quite easily.
p=plot(sx+margin,sy+margin,'r.',/overplot,sym_thick=2)
p=plot(i4x+margin,i4y+margin,'g.',/overplot,sym_thick=2)

c=contour(reform((betal(100,*,*))),color='r',c_value=1,/overplot) 
;Use the divergence of v as a background
;c=contour(divv(*,*,100),/fill,n_levels=101,rgb_tab=0)
;p=plot(fx+margin,fy,'b.',/overplot,sym_thick=2) ;Margin not needed here
;p=plot(sx+margin,sy,'r.',/overplot,sym_thick=2)


;ps=plot3d(shocks.slowx,shocks.slowy,shocks.slowz,'r.')
;ps.sym_size=4 
;pf=plot3d(shocks.fastx,shocks.fasty,shocks.fastz,'b.',/overplot)
;pf.sym_size=4 
;pi4=plot3d(shocks.int4x,shocks.int4y,shocks.int4z,'g.',/overplot)
;pi4.sym_size=4 
stop

;Save the buffered plot
p.save,'test3d2.pdf'

END
