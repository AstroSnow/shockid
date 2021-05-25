;plasmoid_ex

restore,'MHD_shockdata_t_72.sav'
i=image(ds.j_z,ds.x,ds.y,xr=[0.8,0.91],yr=[0.3,0.4],axis_sty=2,dim=[970,600],xtitle='x'$
, ytitle='y',font_size=14)

ps=plot(1.0*shocks.slowx/(1.0*n_elements(shocks.x)),1.0*shocks.slowy/(1.0*n_elements(shocks.y)),/overplot,'.r',sym_thick=10)

ps=plot(1.0*shocks.fastx/(1.0*n_elements(shocks.x)),1.0*shocks.fasty/(1.0*n_elements(shocks.y)),/overplot,'.b',sym_thick=10)

pi1=plot(1.0*shocks.int1x/(1.0*n_elements(shocks.x)),1.0*shocks.int1y/(1.0*n_elements(shocks.y)),/overplot,'.y',sym_thick=10)

pi2=plot(1.0*shocks.int2x/(1.0*n_elements(shocks.x)),1.0*shocks.int2y/(1.0*n_elements(shocks.y)),/overplot,'.m',sym_thick=10)

pi3=plot(1.0*shocks.int3x/(1.0*n_elements(shocks.x)),1.0*shocks.int3y/(1.0*n_elements(shocks.y)),/overplot,'.c',sym_thick=10)

pi4=plot(1.0*shocks.int4x/(1.0*n_elements(shocks.x)),1.0*shocks.int4y/(1.0*n_elements(shocks.y)),/overplot,'.g',sym_thick=10)


END
