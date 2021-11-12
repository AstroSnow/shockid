;Test of the 3d shockid function for Henrik

RESOLVE_ROUTINE, 'shockid3d_fun2';, /IS_FUNCTION

;Location of your data files
;fname='../PIPstab3D/Data/'
;fname='../../datadrive/PIPstab3D/PIP3D/'
fname='../../datadrive/simdata/MHD_OZ_3D'

;Time step to read in
t=2

;Read in the data
rdmpi,ds,datapath=fname,time_step=t

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
shockid3d_fun2,shocks,ds.x,ds.y,ds.z,ds.ro_p,ds.vx_p,ds.vy_p,ds.vz_p,ds.bx,ds.by,ds.bz,ds.pr_p,ndim=3

; Print out the sizes of the elements of a
help,shocks

;Save the shocks
save,shocks,filename='test3d2.dat'

;Some basic plotting routines

;Plot a 2D slice of the data at z=100 grid cell
fx=shocks.fastx(where(shocks.fastz eq 100))
fy=shocks.fasty(where(shocks.fastz eq 100))
sx=shocks.slowx(where(shocks.slowz eq 100))
sy=shocks.slowy(where(shocks.slowz eq 100))
c=contour(ds.ro_p(*,*,100),/fill,n_levels=101,rgb_tab=0)
p=plot(fx+6,fy+6,'b.',/overplot) ;The 6 is there because of the ghost cells in my code that the shockid routine assumes. We can change this quite easily.
p=plot(sx+6,sy+6,'r.',/overplot)

stop

;3D plot of the shocks - VERY SLOW TO RUN
p=plot3d(shocks.slowx(0)*[1,1],shocks.slowy(0)*[1,1],shocks.slowz(0)*[1,1],'ro',/buffer)
for j=1,n_elements(shocks.slowx)-1 do begin
    p=plot3d(shocks.slowx(j)*[1,1],shocks.slowy(j)*[1,1],shocks.slowz(j)*[1,1],'ro',/overplot)
end
for j=0,n_elements(shocks.fastx)-1 do begin
    p=plot3d(shocks.fastx(j)*[1,1],shocks.fasty(j)*[1,1],shocks.fastz(j)*[1,1],'bo',/overplot)
end

;Save the buffered plot
p.save,'test3d2.pdf'

END
