;Identify and clasify shocks based on the SHOCKFIND algorythm (https://arxiv.org/abs/1608.02050)

;Compile a few necessary function
RESOLVE_ROUTINE, 'exist', /IS_FUNCTION
RESOLVE_ROUTINE, 'where_vector', /IS_FUNCTION
RESOLVE_ROUTINE, 'interp2d', /IS_FUNCTION
RESOLVE_ROUTINE, 'datatype', /IS_FUNCTION
RESOLVE_ROUTINE, 'shocknormvals', /IS_FUNCTION

;Input directory
;fname='../MHD_test'
fname='Sims/MHD_test' & tread=100
;fname='../Shock_stab/PIP/Sims/Data_4000x200_MHD_test' & tread=50


;Read in the data
print,'Reading data'
rdmpi,ds,datapath=fname,time_step=tread

;Parameters for shock limits
convl=400.0 ;Convergence threshold
avecyl=5 ;Cylinder to average over
smthfac=1 ;smoothing factor (1=no smoothing)

;Define simulation grid
print,'Removing ghost cells'
margin=6

egx=N_elements(ds.x)-1
egy=N_elements(ds.y)-1

x=ds.x(margin:egx-margin)
y=ds.y(margin:egy-margin)
ro=smooth(ds.ro_p(margin:egx-margin,margin:egy-margin),smthfac)
vx=smooth(ds.vx_p(margin:egx-margin,margin:egy-margin),smthfac)
vy=smooth(ds.vy_p(margin:egx-margin,margin:egy-margin),smthfac)
vz=smooth(ds.vz_p(margin:egx-margin,margin:egy-margin),smthfac)
bx=smooth(ds.bx(margin:egx-margin,margin:egy-margin),smthfac)
by=smooth(ds.by(margin:egx-margin,margin:egy-margin),smthfac)
bz=smooth(ds.bz(margin:egx-margin,margin:egy-margin),smthfac)
pr=smooth(ds.pr_p(margin:egx-margin,margin:egy-margin),smthfac)

;2D divergence of velocity
print,'Calculating divergence'
divv=dblarr(n_elements(x),n_elements(y))
;for i=0,n_elements(x)-1 do begin
;for j=0,n_elements(y)-1 do begin
;    divv(i,j)=(ds.vx_p(margin+i+1,margin+j)-ds.vx_p(margin+i-1,margin+j))/(ds.x(margin+i+1)-ds.x(margin+i-1)) $
;+(ds.vy_p(margin+i,margin+j+1)-ds.vy_p(margin+i,margin+j-1))/(ds.y(margin+j+1)-ds.y(margin+j-1))
;endfor
;endfor
;In vector form (so quicker)
divv=(smooth(ds.vx_p(margin+1:egx-margin+1,margin:egy-margin),smthfac)-smooth(ds.vx_p(margin-1:egx-margin-1,margin:egy-margin),smthfac))/(ds.x(10)-ds.x(8)) $
+(smooth(ds.vy_p(margin:egx-margin,margin+1:egy-margin+1),smthfac)-smooth(ds.vy_p(margin:egx-margin,margin-1:egy-margin-1),smthfac))/(ds.y(10)-ds.y(8))


;Maximum gradient of density
print,'Calculating density gradients'
gradx=dblarr(n_elements(x),n_elements(y))
grady=dblarr(n_elements(x),n_elements(y))
gradz=dblarr(n_elements(x),n_elements(y))
;for i=0,n_elements(x)-1 do begin
;for j=0,n_elements(y)-1 do begin
;    gradx(i,j)=(ds.ro_p(margin+i+1,margin+j)-ds.ro_p(margin+i-1,margin+j))/(ds.x(margin+i+1)-ds.x(margin+i-1))
;    grady(i,j)=(ds.ro_p(margin+i,margin+j+1)-ds.ro_p(margin+i,margin+j-1))/(ds.y(margin+j+1)-ds.y(margin+j-1))
;    gradz(i,j)=0.0
;endfor
;endfor
;In vector form (so quicker)
gradx=(smooth(ds.ro_p(margin+1:egx-margin+1,margin:egy-margin),smthfac)-smooth(ds.ro_p(margin-1:egx-margin-1,margin:egy-margin),smthfac))/(ds.x(10)-ds.x(8))
grady=(smooth(ds.ro_p(margin:egx-margin,margin+1:egy-margin+1),smthfac)-smooth(ds.ro_p(margin:egx-margin,margin-1:egy-margin-1),smthfac))/(ds.y(10)-ds.y(8))
gradz=0.0*dblarr(n_elements(x),n_elements(y))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;Step i: identify candidate cells based on 
print,'Identify candidate cells'
gradmag=sqrt(gradx^2+grady^2+gradz^2)

;Find candidate cells
index=where(divv LT -convl)
s = SIZE(divv)
ncol = s(1)
col = index MOD ncol
row = index / ncol
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

fastx=[]
fasty=[]
fastm=[]

slowx=[]
slowy=[]
slowm=[]

int1x=[]
int1y=[]

int2x=[]
int2y=[]

int3x=[]
int3y=[]

int4x=[]
int4y=[]

;create arrays to store shock normal data
arrronorm=dblarr(N_elements(col),avecyl*2+1)
arrvxnorm=dblarr(N_elements(col),avecyl*2+1)
arrvynorm=dblarr(N_elements(col),avecyl*2+1)
arrbxnorm=dblarr(N_elements(col),avecyl*2+1)
arrbynorm=dblarr(N_elements(col),avecyl*2+1)
arrprnorm=dblarr(N_elements(col),avecyl*2+1)
arrvpar=dblarr(N_elements(col),avecyl*2+1)
arrvperp=dblarr(N_elements(col),avecyl*2+1)
arrvslow2=dblarr(N_elements(col),avecyl*2+1)
arrvalf2=dblarr(N_elements(col),avecyl*2+1)
arrvfast2=dblarr(N_elements(col),avecyl*2+1)
arrtype=dblarr(N_elements(col))
arrvs=dblarr(N_elements(col))
arrnx=dblarr(N_elements(col))
arrny=dblarr(N_elements(col))
arrpx=dblarr(N_elements(col))
arrpy=dblarr(N_elements(col))

print,'Begin loop'
;loop over candidates
for i=0,n_elements(col)-1 do begin

    print,i,(n_elements(col)-1)
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ;Step ii: find shock normal based on density gradient
    normx=gradx(col(i),row(i))/gradmag(col(i),row(i))
    normy=grady(col(i),row(i))/gradmag(col(i),row(i))

    ;Find perpendicular direction too
    perpx=(-normy/normx)/sqrt(normy^2/normx^2+1.0)
    perpy=1.0/sqrt(normy^2/normx^2+1.0)

    ;print,normx,normy,gradx(col(i),row(i)),grady(col(i),row(i))

    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ;Step iii) 

    ;First create a (2*avecyl)^2 box of values
;    tempa=dblarr(2*avecyl+1,2*avecyl+1)
    tempx=indgen(2*avecyl+1)-avecyl+col(i)
    tempy=indgen(2*avecyl+1)-avecyl+row(i)
    tempx2=indgen(2*avecyl+1)-avecyl
    tempy2=indgen(2*avecyl+1)-avecyl

    ;use periodic BC to fix negative values
    for ii=0,2*avecyl do begin
        if tempx(ii) lt 0 then tempx(ii)=N_elements(divv(*,0))+tempx(ii)
        if tempy(ii) lt 0 then tempy(ii)=N_elements(divv(0,*))+tempy(ii)
    endfor

    ;use periodic BC to fix outside grid values
    for ii=0,2*avecyl do begin
        if tempx(ii) gt n_elements(divv(*,0))-1 then tempx(ii)=tempx(ii)-N_elements(divv(*,0))
        if tempy(ii) gt n_elements(divv(0,*))-1 then tempy(ii)=tempy(ii)-N_elements(divv(0,*))
    endfor


    ;Interpolate the values normal to the shock NOTE: THES ARE NOT PARALLEL TO SHOCK YET!
    ronorm=shocknormvals(ro,tempx,tempy,tempx2,tempy2,avecyl,normx,normy)
    vxnorm=shocknormvals(vx,tempx,tempy,tempx2,tempy2,avecyl,normx,normy)
    vynorm=shocknormvals(vy,tempx,tempy,tempx2,tempy2,avecyl,normx,normy)
    bxnorm=shocknormvals(bx,tempx,tempy,tempx2,tempy2,avecyl,normx,normy)
    bynorm=shocknormvals(by,tempx,tempy,tempx2,tempy2,avecyl,normx,normy)
    prnorm=shocknormvals(pr,tempx,tempy,tempx2,tempy2,avecyl,normx,normy)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Check for maximum density gradient along the normal
    ;rgnorm=dblarr(2*avecyl+1)
    ;for ii=1,2*avecyl-2 do begin
    ;    rgnorm(ii)=abs(ronorm(ii+1)-ronorm(ii-1))
    ;endfor

    ;if (max(rgnorm)-abs(rgnorm(avecyl))) lt 1.0e-6 then begin

    ;a=max(abs(deriv(ronorm)),b)
    ;if b eq avecyl then begin

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;I think these are the wrong way round
    ;Calculate the parallel and perpendicular velocity and magnetic field
;    vpar=normx*vxnorm+normy*vynorm
    ;vperp=sqrt(vxnorm^2+vynorm^2-vpar^2)
;    vperp=perpx*vxnorm+perpy*vynorm
;    bpar=normx*bxnorm+normy*bynorm
    ;bperp=sqrt(bxnorm^2+bynorm^2-bpar^2)
;    bperp=perpx*bxnorm+perpy*bynorm
;    ang=atan(bpar/bperp) ;DOUBLE CHECK THIS
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

    ;Calculate the parallel and perpendicular velocity and magnetic field
    vperp=normx*vxnorm+normy*vynorm
    ;vperp=sqrt(vxnorm^2+vynorm^2-vpar^2)
    vpar=perpx*vxnorm+perpy*vynorm
    bperp=normx*bxnorm+normy*bynorm
    ;bperp=sqrt(bxnorm^2+bynorm^2-bpar^2)
    bpar=perpx*bxnorm+perpy*bynorm
    ang=atan(bpar/bperp) ;DOUBLE CHECK THIS

    ;define pre and post shock indicies
    ipre=0
    ipos=2*avecyl

    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ;Estimate shock speed 
;    vs=(vpar(ipre)-vpar(ipos))/(1.0-ronorm(ipre)/ronorm(ipos))
    vs=(vperp(ipre)-vperp(ipos))/(1.0-ronorm(ipre)/ronorm(ipos))
    ;Calculate Mach numbers
    cs2=5.0/3.0*prnorm/ronorm
    va2=(bxnorm^2+bynorm^2)/ronorm
    vap2=bperp^2/ronorm
    vslow2=0.5*(cs2+va2-sqrt((cs2+va2)^2 - 4.0*va2*cs2*(cos(ang))^2))
    vfast2=0.5*(cs2+va2+sqrt((cs2+va2)^2 - 4.0*va2*cs2*(cos(ang))^2))

;    p=plot((vpar-vs)^2)
;    p=plot(cs2,/overplot,'r')
;    p=plot(va2,/overplot,'b')
;    p=plot(vslow2,/overplot,'g')
;    p=plot(vfast2,/overplot,'m')

    vf=vperp
    if vs LT 0.0 then vf=vperp-vs
    if vs GT 0.0 then vf=vperp+vs
    
    ;if ronorm(ipre) LT ronorm(ipos) then vf=vperp-vs
    ;if ronorm(ipos) LT ronorm(ipre) then vf=vperp+vs

    ;if min(abs(vperp)) LT vs then vf=vperp+vs
    ;if min(abs(vperp)) GT vs then vf=vperp-vs

    ;alocate data to perminant array
    arrronorm(i,*)=ronorm
    arrvxnorm(i,*)=vxnorm
    arrvynorm(i,*)=vynorm
    arrbxnorm(i,*)=bxnorm
    arrbynorm(i,*)=bynorm
    arrprnorm(i,*)=prnorm
    arrvpar(i,*)=vpar
    arrvperp(i,*)=vperp
    arrvslow2(i,*)=vslow2
    arrvfast2(i,*)=vfast2
    arrvalf2(i,*)=vap2
    arrvs(i)=vs
    arrnx(i)=normx
    arrny(i)=normy
    arrpx(i)=perpx
    arrpy(i)=perpy

    ;Test if it is a shock
;    if (max(vf/sqrt(vslow2)) GT 1) and (min(vf/sqrt(vslow2)) LT 1) then begin
;        slowx=[slowx,col(i)]
;        slowy=[slowy,row(i)]
;    endif

;    if (max(vf/sqrt(vfast2)) GT 1) and (min(vf/sqrt(vfast2)) LT 1) then begin
;        fastx=[fastx,col(i)]
;        fasty=[fasty,row(i)]
;    endif 

    ;Test if it is a shock
;    if (vf(ipre)/sqrt(vslow2(ipre)) GT 1) and (vf(ipos)/sqrt(vslow2(ipos)) LT 1) then begin
;        slowx=[slowx,col(i)]
;        slowy=[slowy,row(i)]
;    endif

;    if (vf(ipre)/sqrt(vfast2(ipre)) GT 1) and (vf(ipos)/sqrt(vfast2(ipos)) LT 1) then begin
;        fastx=[fastx,col(i)]
;        fasty=[fasty,row(i)]
;    endif 
    
;    if (max(vf(ipre:avecyl)/sqrt(vslow2(ipre:avecyl))) GT 1) and (min(vf(avecyl+1:ipos)/sqrt(vslow2(avecyl+1:ipos))) LT 1) then begin
;        slowx=[slowx,col(i)]
;        slowy=[slowy,row(i)]
;    endif

;    if (max(vf(ipre:avecyl)/sqrt(vfast2(ipre:avecyl))) GT 1) and (min(vf(avecyl+1:ipos)/sqrt(vfast2(avecyl+1:ipos))) LT 1) then begin
;        fastx=[fastx,col(i)]
;        fasty=[fasty,row(i)]
;    endif 


;Clasification based on transitions
    vslowpre =mean(abs(vf(ipre:avecyl-3)/sqrt(vslow2(ipre:avecyl-3))))
    valfpre  =mean(abs(vf(ipre:avecyl-3)/sqrt(vap2(ipre:avecyl-3))))
    vfastpre =mean(abs(vf(ipre:avecyl-3)/sqrt(vfast2(ipre:avecyl-3))))

    tpre=4
    if vfastpre gt 1 then tpre=1 else $
    if (valfpre gt 1) and (vfastpre lt 1) then tpre=2 else $
    if (vslowpre gt 1) and (valfpre lt 1) then tpre=3 

    vslowpos =mean(abs(vf(avecyl+4:ipos)/sqrt(vslow2(avecyl+4:ipos))))
    valfpos  =mean(abs(vf(avecyl+4:ipos)/sqrt(vap2(avecyl+4:ipos))))
    vfastpos =mean(abs(vf(avecyl+4:ipos)/sqrt(vfast2(avecyl+4:ipos))))

    tpos=4
    if vfastpos gt 1 then tpos=1 else $
    if (valfpos gt 1) and (vfastpos lt 1) then tpos=2 else $
    if (vslowpos gt 1) and (valfpos lt 1) then tpos=3 

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Check for maximum density gradient along the normal
    a=max(abs(deriv(ronorm)),b)
    if b eq avecyl then begin



    if (tpre eq 1) and (tpos eq 2) then begin
        fastx=[fastx,col(i)]
        fasty=[fasty,row(i)]
        arrtype(i)=6
    endif

    if (tpre eq 2) and (tpos eq 1) then begin
        fastx=[fastx,col(i)]
        fasty=[fasty,row(i)]
        arrtype(i)=6
    endif

    if (tpre eq 3) and (tpos eq 4) then begin
        slowx=[slowx,col(i)]
        slowy=[slowy,row(i)]
        arrtype(i)=1
    endif 

    if (tpre eq 4) and (tpos eq 3) then begin
        slowx=[slowx,col(i)]
        slowy=[slowy,row(i)]
        arrtype(i)=1
    endif  

    if (tpre eq 1) and (tpos eq 3) then begin
        int1x=[int1x,col(i)]
        int1y=[int1y,row(i)]
        arrtype(i)=4
    endif   

    if (tpre eq 3) and (tpos eq 1) then begin
        int1x=[int1x,col(i)]
        int1y=[int1y,row(i)]
        arrtype(i)=4
    endif   

    if (tpre eq 1) and (tpos eq 4) then begin
        int2x=[int2x,col(i)]
        int2y=[int2y,row(i)]
        arrtype(i)=3
    endif  

    if (tpre eq 4) and (tpos eq 1) then begin
        int2x=[int2x,col(i)]
        int2y=[int2y,row(i)]
        arrtype(i)=3
    endif   

    if (tpre eq 2) and (tpos eq 3) then begin
        int3x=[int3x,col(i)]
        int3y=[int3y,row(i)]
        arrtype(i)=5
    endif   

    if (tpre eq 3) and (tpos eq 2) then begin
        int3x=[int3x,col(i)]
        int3y=[int3y,row(i)]
        arrtype(i)=5
    endif   

    if (tpre eq 2) and (tpos eq 4) then begin
        int4x=[int4x,col(i)]
        int4y=[int4y,row(i)]
        arrtype(i)=2
    endif 

    if (tpre eq 4) and (tpos eq 2) then begin
        int4x=[int4x,col(i)]
        int4y=[int4y,row(i)]
        arrtype(i)=2
    endif    

    endif

endfor

;c=image(divv,x,y,rgb_tab=0,xtitle='x',ytitle='y',dim=[970,600], position=[0.1,0.1,0.7,0.9],axis_style=2,xr=[0,1],yr=[0,1],/buffer)
;c=contour(divv,x,y,/fill,n_levels=101,rgb_tab=0,xtitle='x',ytitle='y',dim=[970,600], position=[0.1,0.1,0.7,0.9],xr=[0,1],yr=[0,1])
diva=divv                     
;diva[where(divv gt 0)]=1.0/0.0
diva[where(divv gt -convl)]=1.0/0.0
c=image(alog10(-diva),x,y,rgb_tab=0,xtitle='x',ytitle='y',dim=[970,600], position=[0.1,0.1,0.7,0.9],axis_style=2,xr=[0,1],yr=[0,1])

if N_elements(slowx) gt 0 then $
ps=plot(1.0*slowx/(n_elements(x)-1),1.0*slowy/(n_elements(y)-1),/overplot,'.r',name='slow',thick=0.1)
if N_elements(fastx) gt 0 then $
pf=plot(1.0*fastx/(n_elements(x)-1),1.0*fasty/(n_elements(y)-1),/overplot,'.b',name='fast',thick=0.1)
if N_elements(int1x) gt 0 then $
pi1=plot(1.0*int1x/(n_elements(x)-1),1.0*int1y/(n_elements(y)-1),/overplot,'.y',name='Int (1-3)',thick=0.1)
if N_elements(int2x) gt 0 then $
pi2=plot(1.0*int2x/(n_elements(x)-1),1.0*int2y/(n_elements(y)-1),/overplot,'.m',name='Int (1-4)',thick=0.1)
if N_elements(int3x) gt 0 then $
pi3=plot(1.0*int3x/(n_elements(x)-1),1.0*int3y/(n_elements(y)-1),/overplot,'.c',name='Int (2-3)',thick=0.1)
if N_elements(int4x) gt 0 then $
pi4=plot(1.0*int4x/(n_elements(x)-1),1.0*int4y/(n_elements(y)-1),/overplot,'.g',name='Int (2-4)',thick=0.1)

;to get the labels right:
ps=plot([-1,-1],[-1,-1],/overplot,'r',name='slow (3-4)',thick=2)
pf=plot([-1,-1],[-1,-1],/overplot,'b',name='fast (1-2)',thick=2)
pi1=plot([-1,-1],[-1,-1],/overplot,'y',name='Int (1-3)',thick=2)
pi2=plot([-1,-1],[-1,-1],/overplot,'m',name='Int (1-4)',thick=2)
pi3=plot([-1,-1],[-1,-1],/overplot,'c',name='Int (2-3)',thick=2)
pi4=plot([-1,-1],[-1,-1],/overplot,'g',name='Int (2-4)',thick=2)

l=legend(target=[ps,pf,pi1,pi2,pi3,pi4])
l.position=[0.9,0.9]

;Annotate the number of shocks
ts=text(0.75,0.6,'Slow = '+strtrim(n_elements(slowx),1),/overplot)
ts=text(0.75,0.55,'Fast = '+strtrim(n_elements(fastx),1),/overplot)
ts=text(0.75,0.5,'Int (1-3) = '+strtrim(n_elements(int1x),1),/overplot)
ts=text(0.75,0.45,'Int (1-4) = '+strtrim(n_elements(int2x),1),/overplot)
ts=text(0.75,0.4,'Int (2-3) = '+strtrim(n_elements(int3x),1),/overplot)
ts=text(0.75,0.35,'Int (2-4) = '+strtrim(n_elements(int4x),1),/overplot)

c.save,'shockid3_r.pdf'

save,arrronorm,arrvxnorm,arrvynorm,arrbxnorm,arrbynorm,$
arrprnorm,arrvpar,arrvperp,arrvslow2,arrvalf2,arrvfast2,$
arrvs,arrtype,col,row,filename='shockid3_data.dat'
;Order the perminant arrays
arrronorms=arrronorm(sort(arrtype),*)
arrvxnorms=arrvxnorm(sort(arrtype),*)
arrvynorms=arrvynorm(sort(arrtype),*)
arrbxnorms=arrbxnorm(sort(arrtype),*)
arrbynorms=arrbynorm(sort(arrtype),*)
arrprnorms=arrprnorm(sort(arrtype),*)
arrvpars=arrvpar(sort(arrtype),*)
arrvperps=arrvperp(sort(arrtype),*)
arrvslow2s=arrvslow2(sort(arrtype),*)
arrvfast2s=arrvfast2(sort(arrtype),*)
arrvalf2s=arrvalf2(sort(arrtype),*)
arrvss=arrvs(sort(arrtype),*)
arrtypes=arrtype(sort(arrtype))


;p=plot(1.0*col/(n_elements(x)-1.0),1.0*row/(n_elements(y)-1.0),/overplot,'.')


END
