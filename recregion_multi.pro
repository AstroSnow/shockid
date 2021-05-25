;Plot the current density with shock ID's

;use t=70,MHD
resolve_routine,'shockid_fun'

w=window(dim=[900,900],/buffer)

;fname='Sims/MHD_test' & tread=70
fname='/home/snow/datadrive/simdata/MHD_OZ_1024/Data' 
tread=[0,58,59,60,61,62,63,64,65,66]

for ti=1,9 do begin

shockid_fun,shocks,fname=fname,tread=tread(ti)

rdmpi,ds,datapath=fname,time_step=tread(ti),current=1,var=['bx','by','bz']
i=image(ds.j_z,ds.x,ds.y,xr=[0.3,0.7],yr=[0.3,0.7],axis_sty=2,dim=[970,600],xtitle='x', ytitle='y',layout=[3,3,ti],/current,margin=[0,0,0,0],/buffer)
;v=streamline(ds.bx,ds.by,ds.x+5*dx,ds.y+5*dx,streamline_stepsize=0.01,x_stream=101,y_stream=101,arrow_size=0.1,/overplot,streamline_nstep=200)

xshift=5
yshift=5

symth=3

ps=plot(1.0*shocks.slowx/(1.0*n_elements(shocks.x)-xshift),1.0*shocks.slowy/(1.0*n_elements(shocks.y)-yshift),/overplot,'.r',sym_thick=symth)

pf=plot(1.0*shocks.fastx/(1.0*n_elements(shocks.x)-xshift),1.0*shocks.fasty/(1.0*n_elements(shocks.y)-yshift),/overplot,'.b',sym_thick=symth)

pi1=plot(1.0*shocks.int1x/(1.0*n_elements(shocks.x)-xshift),1.0*shocks.int1y/(1.0*n_elements(shocks.y)-yshift),/overplot,'.y',sym_thick=symth)

pi2=plot(1.0*shocks.int2x/(1.0*n_elements(shocks.x)-xshift),1.0*shocks.int2y/(1.0*n_elements(shocks.y)-yshift),/overplot,'.m',sym_thick=symth)

pi3=plot(1.0*shocks.int3x/(1.0*n_elements(shocks.x)-xshift),1.0*shocks.int3y/(1.0*n_elements(shocks.y)-yshift),/overplot,'.c',sym_thick=symth)

pi4=plot(1.0*shocks.int4x/(1.0*n_elements(shocks.x)-xshift),1.0*shocks.int4y/(1.0*n_elements(shocks.y)-yshift),/overplot,'.g',sym_thick=symth)

endfor

i.save,'recregion_multi.pdf'

END
