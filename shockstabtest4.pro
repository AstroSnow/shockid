;Identify and clasify shocks based on the SHOCKFIND algorythm (https://arxiv.org/abs/1608.02050)

;Compile a few necessary function
RESOLVE_ROUTINE, 'exist', /IS_FUNCTION
RESOLVE_ROUTINE, 'where_vector', /IS_FUNCTION
RESOLVE_ROUTINE, 'interp2d', /IS_FUNCTION
RESOLVE_ROUTINE, 'datatype', /IS_FUNCTION
RESOLVE_ROUTINE, 'shocknormvals', /IS_FUNCTION
RESOLVE_ROUTINE, 'shockstabtestfunc', /IS_FUNCTION

w=window(dim=[970,600])

;Basic properties (not yet passed to the function
;Parameters for shock limits
convl=0.1 ;Convergence threshold
avecyl=5 ;Cylinder to average over
smthfac=5 ;smoothing factor (1=no smoothing)

;Input directory
;fname='../MHD_test'
;fname='Sims/MHD_test' & tread=100
;fname='../Shock_stab/PIP/Sims/Data_4000x200_MHD_test' & tread=10
fname='../Shock_stab/PIP/Sims/MHD_sf_o4_hr' & tread=30
i=0
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;Read in the data
print,'Reading data'
rdmpi,ds,datapath=fname,time_step=tread

;Define simulation grid
print,'Removing ghost cells'
margin=6

egx=N_elements(ds.x)-1
egy=N_elements(ds.y)-1

x=ds.x(margin:egx-margin)
y=ds.y(margin:egy-margin)
ro=smooth(ds.ro_p(margin:egx-margin,margin:egy-margin),smthfac)

;use rgb_table=0 (black),62(red),49(blue)
c=image(ro,x,y,rgb_table=0,$
position=[0.1,0.1,0.9,0.9],/current,xr=[-4,4],yr=[0,10],axis_style=2,ytitle='y',aspect_ratio=0.3, font_size=16,xtitle='x')

aslow=shockstabtestfunc(fname,tread)

;c=image(ro,x,y,rgb_tab=0,xtitle='x',ytitle='y',dim=[970,600], position=[0.1,0.1,0.7,0.9],axis_style=2,aspect_ratio=1,xr=[-4,4])
ps=plot(aslow(*,0),aslow(*,1),/overplot,'.r',name='slow',sym_thick=2)


end
