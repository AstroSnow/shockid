;Plot the current density with shock ID's

;use t=70,MHD
resolve_routine,'shockid_fun'

;fname='Sims/MHD_test' & tread=70
fname='/home/snow/datadrive/simdata/MHD_OZ_1024/Data' & tread=70

shockid_fun,shocks,fname=fname,tread=tread

rdmpi,ds,datapath=fname,time_step=tread,current=1,var=['bx','by','bz']
i=image(ds.j_z,ds.x,ds.y,xr=[0.7,1.0],yr=[0.2,0.5],axis_sty=2,dim=[970,600],xtitle='x',ytitle='y')
;v=streamline(ds.bx,ds.by,ds.x+5*dx,ds.y+5*dx,streamline_stepsize=0.01,x_stream=101,y_stream=101,arrow_size=0.1,/overplot,streamline_nstep=200)

xshift=5
yshift=5

ps=plot(1.0*shocks.slowx/(1.0*n_elements(shocks.x)-xshift),1.0*shocks.slowy/(1.0*n_elements(shocks.y)-yshift),/overplot,'.r',sym_thick=6)

ps=plot(1.0*shocks.fastx/(1.0*n_elements(shocks.x)-xshift),1.0*shocks.fasty/(1.0*n_elements(shocks.y)-yshift),/overplot,'.b',sym_thick=6)

pi1=plot(1.0*shocks.int1x/(1.0*n_elements(shocks.x)-xshift),1.0*shocks.int1y/(1.0*n_elements(shocks.y)-yshift),/overplot,'.y',sym_thick=6)

pi2=plot(1.0*shocks.int2x/(1.0*n_elements(shocks.x)-xshift),1.0*shocks.int2y/(1.0*n_elements(shocks.y)-yshift),/overplot,'.m',sym_thick=6)

pi3=plot(1.0*shocks.int3x/(1.0*n_elements(shocks.x)-xshift),1.0*shocks.int3y/(1.0*n_elements(shocks.y)-yshift),/overplot,'.c',sym_thick=6)

pi4=plot(1.0*shocks.int4x/(1.0*n_elements(shocks.x)-xshift),1.0*shocks.int4y/(1.0*n_elements(shocks.y)-yshift),/overplot,'.g',sym_thick=6)
END
