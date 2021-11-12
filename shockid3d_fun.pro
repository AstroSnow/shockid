;Identify and clasify shocks based on the SHOCKFIND algorythm (https://arxiv.org/abs/1608.02050)
pro shockid3d_fun,shocks,fname=fname,tread=tread,ndim=ndim
  Compile_Opt DEFINT32  

;Compile a few necessary function
RESOLVE_ROUTINE, 'exist', /IS_FUNCTION
RESOLVE_ROUTINE, 'where_vector', /IS_FUNCTION
RESOLVE_ROUTINE, 'interp2d', /IS_FUNCTION
RESOLVE_ROUTINE, 'datatype', /IS_FUNCTION
RESOLVE_ROUTINE, 'shocknormvals', /IS_FUNCTION
RESOLVE_ROUTINE, 'shocknormvals3d', /IS_FUNCTION
RESOLVE_ROUTINE, 'shockplane3d', /IS_FUNCTION

;Input directory
;fname='../MHD_test'
;fname='Sims/MHD_test' & tread=70
;fname='Sims/PIP_test' & tread=70
;fname='Sims/PIP_test_ac100' & tread=100
;fname='../Shock_stab/PIP/Sims/Data_4000x200_MHD_test' & tread=50
;fname='../Shock_stab/PIP/Sims/MHD_sf_o4_hr' & tread=30
;fname='/media/snow/datadrive/OZdata/MHD_OZ_16384/Data' & tread=100

;Read in the data
print,'Reading data'
rdmpi,ds,datapath=fname,time_step=tread;,var=['ro_p','vx_p','vy_p','vz_p']
;Parameters for shock limits
convl=0.01 ;Convergence threshold
avecyl=5 ;Cylinder to average over
smthfac=0 ;smoothing factor (1=no smoothing)
shocktol=0.05 ;tolerence for the shock transitions
bulkspeed=0 ; Use the bulk sound and alfven speeds or individual
species='plasma' ; plasma or neutral
;species='neutral' ; plasma or neutral

;Data handeling flags
data_subset=0 ;flag to subset the data
data_save=0


;Define simulation grid
print,'Removing ghost cells'
margin=6

egx=N_elements(ds.x)-1
egy=N_elements(ds.y)-1
if (ndim eq 3) then egz=N_elements(ds.z)-1

x=ds.x(margin:egx-margin)
y=ds.y(margin:egy-margin)
if (ndim eq 3) then z=ds.z(margin:egz-margin)

;MHD variables
if (ds.fl_pip eq 0) and (ds.fl_mhd eq 1) then begin
    print,'Using MHD quantities'
    if (ndim eq 3) then begin
        ro=smooth(ds.ro_p(margin:egx-margin,margin:egy-margin,margin:egz-margin),smthfac)
        vx=smooth(ds.vx_p(margin:egx-margin,margin:egy-margin,margin:egz-margin),smthfac)
        vy=smooth(ds.vy_p(margin:egx-margin,margin:egy-margin,margin:egz-margin),smthfac)
        vz=smooth(ds.vz_p(margin:egx-margin,margin:egy-margin,margin:egz-margin),smthfac)
        bx=smooth(ds.bx(margin:egx-margin,margin:egy-margin,margin:egz-margin),smthfac)
        by=smooth(ds.by(margin:egx-margin,margin:egy-margin,margin:egz-margin),smthfac)
        bz=smooth(ds.bz(margin:egx-margin,margin:egy-margin,margin:egz-margin),smthfac)
        pr=smooth(ds.pr_p(margin:egx-margin,margin:egy-margin,margin:egz-margin),smthfac)
    endif
    if (ndim eq 2) then begin
        ro=smooth(ds.ro_p(margin:egx-margin,margin:egy-margin),smthfac)
        vx=smooth(ds.vx_p(margin:egx-margin,margin:egy-margin),smthfac)
        vy=smooth(ds.vy_p(margin:egx-margin,margin:egy-margin),smthfac)
        vz=smooth(ds.vz_p(margin:egx-margin,margin:egy-margin),smthfac)
        bx=smooth(ds.bx(margin:egx-margin,margin:egy-margin),smthfac)
        by=smooth(ds.by(margin:egx-margin,margin:egy-margin),smthfac)
        bz=smooth(ds.bz(margin:egx-margin,margin:egy-margin),smthfac)
        pr=smooth(ds.pr_p(margin:egx-margin,margin:egy-margin),smthfac)
    endif
endif
;PIP plasma only variables
if (ds.fl_pip eq 1) and (bulkspeed eq 0) and (species eq 'plasma') then begin
    print,'Using plasma only quantities'
    if (ndim eq 3) then begin
        ro=smooth(ds.ro_p(margin:egx-margin,margin:egy-margin,margin:egz-margin),smthfac)
        vx=smooth(ds.vx_p(margin:egx-margin,margin:egy-margin,margin:egz-margin),smthfac)
        vy=smooth(ds.vy_p(margin:egx-margin,margin:egy-margin,margin:egz-margin),smthfac)
        vz=smooth(ds.vz_p(margin:egx-margin,margin:egy-margin,margin:egz-margin),smthfac)
        bx=smooth(ds.bx(margin:egx-margin,margin:egy-margin,margin:egz-margin),smthfac)
        by=smooth(ds.by(margin:egx-margin,margin:egy-margin,margin:egz-margin),smthfac)
        bz=smooth(ds.bz(margin:egx-margin,margin:egy-margin,margin:egz-margin),smthfac)
        pr=smooth(ds.pr_p(margin:egx-margin,margin:egy-margin,margin:egz-margin),smthfac)
    endif
    if (ndim eq 2) then begin
        ro=smooth(ds.ro_p(margin:egx-margin,margin:egy-margin),smthfac)
        vx=smooth(ds.vx_p(margin:egx-margin,margin:egy-margin),smthfac)
        vy=smooth(ds.vy_p(margin:egx-margin,margin:egy-margin),smthfac)
        vz=smooth(ds.vz_p(margin:egx-margin,margin:egy-margin),smthfac)
        bx=smooth(ds.bx(margin:egx-margin,margin:egy-margin),smthfac)
        by=smooth(ds.by(margin:egx-margin,margin:egy-margin),smthfac)
        bz=smooth(ds.bz(margin:egx-margin,margin:egy-margin),smthfac)
        pr=smooth(ds.pr_p(margin:egx-margin,margin:egy-margin),smthfac)
    endif
endif
;PIP plasma bulk variables
if (ds.fl_pip eq 1) and (bulkspeed eq 1) and (species eq 'plasma') then begin
    print,'Using plasma bulk quantities'
    ro=smooth(ds.ro_p(margin:egx-margin,margin:egy-margin),smthfac)+$
        smooth(ds.ro_n(margin:egx-margin,margin:egy-margin),smthfac)
    vx=smooth(ds.vx_p(margin:egx-margin,margin:egy-margin),smthfac)
    vy=smooth(ds.vy_p(margin:egx-margin,margin:egy-margin),smthfac)
    vz=smooth(ds.vz_p(margin:egx-margin,margin:egy-margin),smthfac)
    bx=smooth(ds.bx(margin:egx-margin,margin:egy-margin),smthfac)
    by=smooth(ds.by(margin:egx-margin,margin:egy-margin),smthfac)
    bz=smooth(ds.bz(margin:egx-margin,margin:egy-margin),smthfac)
    pr=smooth(ds.pr_p(margin:egx-margin,margin:egy-margin),smthfac)+$
        smooth(ds.pr_n(margin:egx-margin,margin:egy-margin),smthfac)
endif
;PIP plasma only variables
if (ds.fl_pip eq 1) and (bulkspeed eq 0) and (species eq 'neutral') then begin
    print,'Using neutral only quantities'
    ro=smooth(ds.ro_n(margin:egx-margin,margin:egy-margin),smthfac)
    vx=smooth(ds.vx_n(margin:egx-margin,margin:egy-margin),smthfac)
    vy=smooth(ds.vy_n(margin:egx-margin,margin:egy-margin),smthfac)
    vz=smooth(ds.vz_n(margin:egx-margin,margin:egy-margin),smthfac)
    pr=smooth(ds.pr_n(margin:egx-margin,margin:egy-margin),smthfac)
endif

;2D divergence of velocity
print,'Calculating divergence'
if (ndim eq 2) then begin
    divv=dblarr(n_elements(x),n_elements(y))
    if (species ne 'neutral') or (ds.fl_pip eq 0) then begin
    divv=(smooth(ds.vx_p(margin+1:egx-margin+1,margin:egy-margin),smthfac)-smooth(ds.vx_p(margin-1:egx-margin-1,margin:egy-margin),smthfac))/(ds.x(10)-ds.x(8)) $
    +(smooth(ds.vy_p(margin:egx-margin,margin+1:egy-margin+1),smthfac)-smooth(ds.vy_p(margin:egx-margin,margin-1:egy-margin-1),smthfac))/(ds.y(10)-ds.y(8))
    endif
    if (species eq 'neutral') and (ds.fl_pip eq 1) then begin
    divv=(smooth(ds.vx_n(margin+1:egx-margin+1,margin:egy-margin),smthfac)-smooth(ds.vx_n(margin-1:egx-margin-1,margin:egy-margin),smthfac))/(ds.x(10)-ds.x(8)) $
    +(smooth(ds.vy_n(margin:egx-margin,margin+1:egy-margin+1),smthfac)-smooth(ds.vy_n(margin:egx-margin,margin-1:egy-margin-1),smthfac))/(ds.y(10)-ds.y(8))
    endif
endif

if (ndim eq 3) then begin
    divv=dblarr(n_elements(x),n_elements(y),n_elements(z))
    if (species ne 'neutral') or (ds.fl_pip eq 0) then begin
    divv=(smooth(ds.vx_p(margin+1:egx-margin+1,margin:egy-margin,margin:egz-margin),smthfac)-$
smooth(ds.vx_p(margin-1:egx-margin-1,margin:egy-margin,margin:egz-margin),smthfac))/(ds.x(10)-ds.x(8)) $
    +(smooth(ds.vy_p(margin:egx-margin,margin+1:egy-margin+1,margin:egz-margin),smthfac)-$
smooth(ds.vy_p(margin:egx-margin,margin-1:egy-margin-1,margin:egz-margin),smthfac))/(ds.y(10)-ds.y(8))$
    +(smooth(ds.vz_p(margin:egx-margin,margin:egy-margin,margin+1:egz-margin+1),smthfac)-$
smooth(ds.vz_p(margin:egx-margin,margin:egy-margin,margin-1:egz-margin-1),smthfac))/(ds.z(10)-ds.z(8))
    endif

endif

;Maximum gradient of density
print,'Calculating density gradients'
gradx=dblarr(n_elements(x),n_elements(y))
grady=dblarr(n_elements(x),n_elements(y))
gradz=dblarr(n_elements(x),n_elements(y))
if (species ne 'neutral') or (ds.fl_pip eq 0) then begin
if ndim eq 2 then begin
gradx=(smooth(ds.ro_p(margin+1:egx-margin+1,margin:egy-margin),smthfac)-smooth(ds.ro_p(margin-1:egx-margin-1,margin:egy-margin),smthfac))/(ds.x(10)-ds.x(8))
grady=(smooth(ds.ro_p(margin:egx-margin,margin+1:egy-margin+1),smthfac)-smooth(ds.ro_p(margin:egx-margin,margin-1:egy-margin-1),smthfac))/(ds.y(10)-ds.y(8))
gradz=0.0*dblarr(n_elements(x),n_elements(y))
endif
if ndim eq 3 then begin
gradx=(smooth(ds.ro_p(margin+1:egx-margin+1,margin:egy-margin,margin:egz-margin),smthfac)-$
smooth(ds.ro_p(margin-1:egx-margin-1,margin:egy-margin,margin:egz-margin),smthfac))/(ds.x(10)-ds.x(8))
grady=(smooth(ds.ro_p(margin:egx-margin,margin+1:egy-margin+1,margin:egz-margin),smthfac)$
-smooth(ds.ro_p(margin:egx-margin,margin-1:egy-margin-1,margin:egz-margin),smthfac))/(ds.y(10)-ds.y(8))
gradz=(smooth(ds.ro_p(margin:egx-margin,margin:egy-margin,margin+1:egz-margin+1),smthfac)$
-smooth(ds.ro_p(margin:egx-margin,margin:egy-margin,margin-1:egz-margin-1),smthfac))/(ds.z(10)-ds.z(8))
endif
endif
if (species eq 'neutral') and (ds.fl_pip eq 1) then begin
gradx=(smooth(ds.ro_n(margin+1:egx-margin+1,margin:egy-margin),smthfac)-smooth(ds.ro_n(margin-1:egx-margin-1,margin:egy-margin),smthfac))/(ds.x(10)-ds.x(8))
grady=(smooth(ds.ro_n(margin:egx-margin,margin+1:egy-margin+1),smthfac)-smooth(ds.ro_n(margin:egx-margin,margin-1:egy-margin-1),smthfac))/(ds.y(10)-ds.y(8))
gradz=0.0*dblarr(n_elements(x),n_elements(y))
endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;Step i: identify candidate cells based on 
print,'Identify candidate cells'
gradmag=sqrt(gradx^2+grady^2+gradz^2)

;Find candidate cells
if ndim eq 2 then begin
index=where(divv LT -convl)
s = SIZE(divv)
ncol = s(1)
col = index MOD ncol
row = index / ncol
endif
if ndim eq 3 then begin
index = WHERE(divv LT -convl)
s = SIZE(divv)
ncol = s[1]
nrow = s[2]
col = index mod ncol
row = (index / ncol) mod nrow
zrow = index / (nrow*ncol)
;p=plot3D(col,row,zrow,'.')
;stop
endif
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;Remove cells that are not a local maximum
print,'Removing non-local maximum density gradients from candidate cells n=',n_elements(col)

col2=[]
row2=[]
zrow2=[]

for i=0,n_elements(col)-1 do begin
;for i=0,1000 do begin

    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ;Step ii: find shock normal based on density gradient
    if ndim eq 2 then begin
        normx=gradx(col(i),row(i))/gradmag(col(i),row(i))
        normy=grady(col(i),row(i))/gradmag(col(i),row(i))
    endif
    if ndim eq 3 then begin
        normx=gradx(col(i),row(i),zrow(i))/gradmag(col(i),row(i),zrow(i))
        normy=grady(col(i),row(i),zrow(i))/gradmag(col(i),row(i),zrow(i))
        normz=gradz(col(i),row(i),zrow(i))/gradmag(col(i),row(i),zrow(i))
    endif

    ;Find perpendicular direction too (not needed here)
;    perpx=(-normy/normx)/sqrt(normy^2/normx^2+1.0)
;    perpy=1.0/sqrt(normy^2/normx^2+1.0)

    tempx=indgen(2*avecyl+1)-avecyl+col(i)
    tempy=indgen(2*avecyl+1)-avecyl+row(i)
    tempx2=indgen(2*avecyl+1)-avecyl
    tempy2=indgen(2*avecyl+1)-avecyl

    if ndim eq 3 then begin
        tempz=indgen(2*avecyl+1)-avecyl+zrow(i)
        tempz2=indgen(2*avecyl+1)-avecyl
    endif

    ;use periodic BC to fix negative values
    for ii=0,2*avecyl do begin
        if ndim eq 2 then begin
            if tempx(ii) lt 0 then tempx(ii)=N_elements(divv(*,0))+tempx(ii)
            if tempy(ii) lt 0 then tempy(ii)=N_elements(divv(0,*))+tempy(ii)
        endif
        if ndim eq 3 then begin
            if tempx(ii) lt 0 then tempx(ii)=N_elements(divv(*,0,0))+tempx(ii)
            if tempy(ii) lt 0 then tempy(ii)=N_elements(divv(0,*,0))+tempy(ii)
            if tempz(ii) lt 0 then tempz(ii)=N_elements(divv(0,0,*))+tempz(ii)
        endif
    endfor

    ;use periodic BC to fix outside grid values
    for ii=0,2*avecyl do begin
        if ndim eq 2 then begin
            if tempx(ii) gt n_elements(divv(*,0))-1 then tempx(ii)=tempx(ii)-N_elements(divv(*,0))
            if tempy(ii) gt n_elements(divv(0,*))-1 then tempy(ii)=tempy(ii)-N_elements(divv(0,*))
        endif
        if ndim eq 3 then begin
            if tempx(ii) gt n_elements(divv(*,0,0))-1 then tempx(ii)=tempx(ii)-N_elements(divv(*,0,0))
            if tempy(ii) gt n_elements(divv(0,*,0))-1 then tempy(ii)=tempy(ii)-N_elements(divv(0,*,0))
            if tempz(ii) gt n_elements(divv(0,0,*))-1 then tempz(ii)=tempz(ii)-N_elements(divv(0,0,*))
        endif
    endfor


    ;Interpolate the values normal to the shock NOTE: THES ARE NOT PARALLEL TO SHOCK YET!
    if ndim eq 2 then $
    ronorm=shocknormvals(ro,tempx,tempy,tempx2,tempy2,avecyl,normx,normy)

    if ndim eq 3 then begin
    tempxyz=dblarr(6,2*avecyl+1)
    tempxyz(0,*)=tempx
    tempxyz(1,*)=tempy
    tempxyz(2,*)=tempz
    tempxyz(3,*)=tempx2
    tempxyz(4,*)=tempy2
    tempxyz(5,*)=tempz2
    ronorm=shocknormvals3D(ro,tempxyz,avecyl,normx,normy,normz)
    endif

    a=max(abs(deriv(ronorm)),b)
    if b eq avecyl then begin
        row2=[row2,row(i)]
        col2=[col2,col(i)]
	if ndim eq 3 then zrow2=[zrow2,zrow(i)]
    endif
endfor
row=row2
col=col2
if ndim eq 3 then zrow=zrow2

;p=plot3D(col,row,zrow,'.',xr=[0,n_elements(x)],yr=[0,n_elements(y)],zr=[0,n_elements(z)])
;print,'ONLY DONE UP TO HERE'
;stop

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;Filter the candidate cells for testing
if data_subset eq 1 then begin
print,'Subsetting data'
xminf=0.28
xmaxf=0.35
yminf=0.22
ymaxf=0.3
zminf=0.1
zmaxf=0.2

col2=[]
row2=[]
zrow2=[]
for i=0,n_elements(row)-1 do begin
    if ndim eq 2 then begin
    if (x(col(i)) LT xmaxf) and (x(col(i)) GT xminf) and (y(row(i)) LT ymaxf) and (y(row(i)) GT yminf) then begin
            row2=[row2,row(i)]
            col2=[col2,col(i)]
    endif
    endif
    if ndim eq 3 then begin
    if (x(col(i)) LT xmaxf) and (x(col(i)) GT xminf) and $
	(y(row(i)) LT ymaxf) and (y(row(i)) GT yminf) and $
	(z(zrow(i)) LT zmaxf) and (z(zrow(i)) GT zminf) then begin
            row2=[row2,row(i)]
            col2=[col2,col(i)]
	    zrow2=[zrow2,zrow(i)]
    endif
    endif
endfor

col=col2
row=row2
if ndim eq 3 then zrow=zrow2

endif
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;Arrays for the shock locations
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

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
if data_save eq 1 then begin
;create arrays to store shock normal data
arrronorm=dblarr(N_elements(col),avecyl*2+1)
;arrvxnorm=dblarr(N_elements(col),avecyl*2+1)
;arrvynorm=dblarr(N_elements(col),avecyl*2+1)
;arrbxnorm=dblarr(N_elements(col),avecyl*2+1)
;arrbynorm=dblarr(N_elements(col),avecyl*2+1)
arrprnorm=dblarr(N_elements(col),avecyl*2+1)
arrvpar=dblarr(N_elements(col),avecyl*2+1)
arrvperp=dblarr(N_elements(col),avecyl*2+1)
arrbpar=dblarr(N_elements(col),avecyl*2+1)
arrbperp=dblarr(N_elements(col),avecyl*2+1)
arrvslow2=dblarr(N_elements(col),avecyl*2+1)
arrvalf2=dblarr(N_elements(col),avecyl*2+1)
arrvfast2=dblarr(N_elements(col),avecyl*2+1)
arrmf=dblarr(N_elements(col),avecyl*2+1)
arrma=dblarr(N_elements(col),avecyl*2+1)
arrms=dblarr(N_elements(col),avecyl*2+1)
arrtype=dblarr(N_elements(col))
arrvs=dblarr(N_elements(col))
;arrnx=dblarr(N_elements(col))
;arrny=dblarr(N_elements(col))
;arrpx=dblarr(N_elements(col))
;arrpy=dblarr(N_elements(col))
;arrrescode=dblarr(N_elements(col))
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

print,'Begin loop for n=',n_elements(col)
;loop over candidates
for i=0,n_elements(col)-1 do begin
;print,'REMEMBER TO CHANGE START TO 0'
;for i=0,225 do begin

    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ;Step ii: find shock normal based on density gradient
    if ndim eq 2 then begin
        normx=gradx(col(i),row(i))/gradmag(col(i),row(i))
        normy=grady(col(i),row(i))/gradmag(col(i),row(i))
    endif
    if ndim eq 3 then begin
        normx=gradx(col(i),row(i),zrow(i))/gradmag(col(i),row(i),zrow(i))
        normy=grady(col(i),row(i),zrow(i))/gradmag(col(i),row(i),zrow(i))
        normz=gradz(col(i),row(i),zrow(i))/gradmag(col(i),row(i),zrow(i))
    endif

    ;Find perpendicular direction too
    perpx=(-normy/normx)/sqrt(normy^2/normx^2+1.0)
    perpy=1.0/sqrt(normy^2/normx^2+1.0)
    if ndim eq 2 then begin
        perpx=(-normy/normx)/sqrt(normy^2/normx^2+1.0)
        perpy=1.0/sqrt(normy^2/normx^2+1.0)
    endif
    if ndim eq 3 then begin
    ;Get the plane perpendicular to the shock 
    pplane=shockplane3d(normx,normy,normz)
    ;Get the density in this plane
    ppro=dblarr(3,3)
    for j=0,2 do begin
    for k=0,2 do begin
        ppro(j,k)=interpolate(ro,pplane(0,j,k)+col(i),pplane(1,j,k)+row(i),pplane(2,j,k)+zrow(i))
    end
    end
    ;get the gradient of thedensity along the plane
    ppgrox=ppro(2,1)-ppro(0,1)
    ppgroy=ppro(1,2)-ppro(1,0)
    if abs(ppgrox) lt 1.0e-8 then begin
        ppx=0.0
        ppy=1.0
    endif else if abs(ppgroy) lt 1.0e-8 then begin
        ppx=1.0
        ppy=0.0
    endif else begin
        ppx=ppgrox/sqrt(ppgrox^2+ppgroy^2)
        ppy=ppgroy/sqrt(ppgrox^2+ppgroy^2)
    endelse
    ;get the perp vectors
    perpx=normy*ppx-normz*ppy
    perpy=-normx*ppx
    perpz=normx*ppy
    ;Make sure its a unit vector
    perpt=sqrt(perpx^2+perpy^2+perpz^2)
    perpx=perpx/perpt
    perpy=perpy/perpt
    perpz=perpz/perpt    
;print,perpx,perpy,perpz,perpt
	;perpx=(-normy-normz)/normx/sqrt(2.0+((-normy-normz)/normx)^2)
	;perpy=1.0/sqrt(2.0+((-normy-normz)/normx)^2)
	;perpz=1.0/sqrt(2.0+((-normy-normz)/normx)^2)
    endif

;test plot of the plane
;p=plot3d(pplane(0,0,0)*[1,1],pplane(1,0,0)*[1,1],pplane(2,0,0)*[1,1],'bo')
;for j=1,2 do begin
;for k=1,2 do begin
;    p=plot3d(pplane(0,j,k)*[1,1],pplane(1,j,k)*[1,1],pplane(2,j,k)*[1,1],/overplot,'bo')
;end
;end

;stop
    ;print,normx,normy,gradx(col(i),row(i)),grady(col(i),row(i))

    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ;Step iii) 

    ;First create a (2*avecyl)^2 box of values
;    tempa=dblarr(2*avecyl+1,2*avecyl+1)
    tempx=indgen(2*avecyl+1)-avecyl+col(i)
    tempy=indgen(2*avecyl+1)-avecyl+row(i)
    tempx2=indgen(2*avecyl+1)-avecyl
    tempy2=indgen(2*avecyl+1)-avecyl

    if ndim eq 3 then begin
        tempz=indgen(2*avecyl+1)-avecyl+zrow(i)
        tempz2=indgen(2*avecyl+1)-avecyl
    endif

    ;use periodic BC to fix negative values
    for ii=0,2*avecyl do begin
        if ndim eq 2 then begin
            if tempx(ii) lt 0 then tempx(ii)=N_elements(divv(*,0))+tempx(ii)
            if tempy(ii) lt 0 then tempy(ii)=N_elements(divv(0,*))+tempy(ii)
        endif
        if ndim eq 3 then begin
            if tempx(ii) lt 0 then tempx(ii)=N_elements(divv(*,0,0))+tempx(ii)
            if tempy(ii) lt 0 then tempy(ii)=N_elements(divv(0,*,0))+tempy(ii)
            if tempz(ii) lt 0 then tempz(ii)=N_elements(divv(0,0,*))+tempz(ii)
        endif
    endfor

    ;use periodic BC to fix outside grid values
    for ii=0,2*avecyl do begin
        if ndim eq 2 then begin
            if tempx(ii) gt n_elements(divv(*,0))-1 then tempx(ii)=tempx(ii)-N_elements(divv(*,0))
            if tempy(ii) gt n_elements(divv(0,*))-1 then tempy(ii)=tempy(ii)-N_elements(divv(0,*))
        endif
        if ndim eq 3 then begin
            if tempx(ii) gt n_elements(divv(*,0,0))-1 then tempx(ii)=tempx(ii)-N_elements(divv(*,0,0))
            if tempy(ii) gt n_elements(divv(0,*,0))-1 then tempy(ii)=tempy(ii)-N_elements(divv(0,*,0))
            if tempz(ii) gt n_elements(divv(0,0,*))-1 then tempz(ii)=tempz(ii)-N_elements(divv(0,0,*))
        endif
    endfor

    ;Interpolate the values normal to the shock NOTE: THES ARE NOT PARALLEL TO SHOCK YET!
    if ndim eq 2 then begin
	ronorm=shocknormvals(ro,tempx,tempy,tempx2,tempy2,avecyl,normx,normy)
	vxnorm=shocknormvals(vx,tempx,tempy,tempx2,tempy2,avecyl,normx,normy)
	vynorm=shocknormvals(vy,tempx,tempy,tempx2,tempy2,avecyl,normx,normy)
	prnorm=shocknormvals(pr,tempx,tempy,tempx2,tempy2,avecyl,normx,normy)

	if (species ne 'neutral') or (ds.fl_pip eq 0) then begin
	bxnorm=shocknormvals(bx,tempx,tempy,tempx2,tempy2,avecyl,normx,normy)
	bynorm=shocknormvals(by,tempx,tempy,tempx2,tempy2,avecyl,normx,normy)
	endif
    endif
    if ndim eq 3 then begin
    tempxyz=dblarr(6,2*avecyl+1)
    tempxyz(0,*)=tempx
    tempxyz(1,*)=tempy
    tempxyz(2,*)=tempz
    tempxyz(3,*)=tempx2
    tempxyz(4,*)=tempy2
    tempxyz(5,*)=tempz2
    ronorm=shocknormvals3d(ro,tempxyz,avecyl,normx,normy,normz)
    vxnorm=shocknormvals3d(vx,tempxyz,avecyl,normx,normy,normz)
    vynorm=shocknormvals3d(vy,tempxyz,avecyl,normx,normy,normz)
    vznorm=shocknormvals3d(vz,tempxyz,avecyl,normx,normy,normz)
    prnorm=shocknormvals3d(pr,tempxyz,avecyl,normx,normy,normz)
	if (species ne 'neutral') or (ds.fl_pip eq 0) then begin
        bxnorm=shocknormvals3d(bx,tempxyz,avecyl,normx,normy,normz)
        bynorm=shocknormvals3d(by,tempxyz,avecyl,normx,normy,normz)
        bznorm=shocknormvals3d(bz,tempxyz,avecyl,normx,normy,normz)
	endif
    endif


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ;Calculate the parallel and perpendicular velocity and magnetic field
    if ndim eq 2 then vperp=normx*vxnorm+normy*vynorm
    if ndim eq 3 then vperp=normx*vxnorm+normy*vynorm+normy*vznorm
    if ndim eq 2 then vpar=perpx*vxnorm+perpy*vynorm
    if ndim eq 3 then vpar=perpx*vxnorm+perpy*vynorm+perpz*vznorm


    if (species ne 'neutral') or (ds.fl_pip eq 0) then begin
	if ndim eq 2 then begin
	    bperp=normx*bxnorm+normy*bynorm
	    ;bperp=sqrt(bxnorm^2+bynorm^2-bpar^2)
	    bpar=perpx*bxnorm+perpy*bynorm
	    ang=atan(bpar/bperp) ;DOUBLE CHECK THIS
	endif
	if ndim eq 3 then begin
	    bperp=normx*bxnorm+normy*bynorm+normz*bznorm
	    bpar=perpx*bxnorm+perpy*bynorm+perpz*bznorm
;	    bmag2=bxnorm^2+bynorm^2+bznorm^2
;	    vmag2=vxnorm^2+vynorm^2+vznorm^2
;	    sinang=sqrt(((bmag2-bperp^2)*vperp - vmag2+vperp^2)/(bmag2+vperp^2-vmag2-bperp^2))
;	    sinang(where (sinang gt  1.0))=1.0
;	    sinang(where (sinang lt -1.0))=-1.0;
;	    ang2=asin(sinang)
;            vpar=sqrt(vmag2-vperp^2)*sin(ang2)
;            bpar=sqrt(bmag2-bperp^2)*sin(ang2)
	    ang=atan(bpar/bperp)
        for bi=0,n_elements(bperp)-1 do begin
            if abs(bperp(bi)) lt 1.0e-8 then ang(bi)=1.5707964
        end
;print,ang
	endif
    endif

    ;define pre and post shock indicies
    if ronorm(0) LT ronorm(2*avecyl) then begin 
;            ipre=0
;            ipos=2*avecyl
            rominpre=min(ronorm(0:avecyl-1),ipre)
            rominpos=max(ronorm(avecyl+1:avecyl*2),ipos)
            ipos=ipos+avecyl+1
    endif else begin
;            ipos=0
;            ipre=2*avecyl
            rominpos=max(ronorm(0:avecyl-1),ipos)
            rominpre=min(ronorm(avecyl+1:avecyl*2),ipre)
            ipre=ipre+avecyl+1
    endelse

    ;calculate conserved quantities
    con_mas=(ronorm(ipos)*vperp(ipos)-ronorm(ipre)*vperp(ipre))/(ronorm(ipre)-ronorm(ipos))


    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ;Estimate shock speed 
    vs=con_mas


    ;Calculate Mach numbers
    cs2=5.0/3.0*prnorm/ronorm
    if (species ne 'neutral') or (ds.fl_pip eq 0) then begin
    if ndim eq 2 then va2=(bxnorm^2+bynorm^2)/ronorm
    if ndim eq 3 then va2=(bxnorm^2+bynorm^2+bznorm^2)/ronorm
    vap2=bperp^2/ronorm
    vslow2=0.5*(cs2+va2-sqrt((cs2+va2)^2 - 4.0*va2*cs2*(cos(ang))^2))
    vfast2=0.5*(cs2+va2+sqrt((cs2+va2)^2 - 4.0*va2*cs2*(cos(ang))^2))
    endif

;print,vap2,vslow2,vfast2

    vf=vperp

    vf=vperp+vs 

    if data_save eq 1 then begin
    ;alocate data to perminant array
    arrronorm(i,*)=ronorm
    arrprnorm(i,*)=prnorm
    arrvpar(i,*)=vpar
    arrvperp(i,*)=vperp
    arrbpar(i,*)=bpar
    arrbperp(i,*)=bperp
    arrvslow2(i,*)=vslow2
    arrvfast2(i,*)=vfast2
    arrvalf2(i,*)=vap2
    arrvs(i)=vs
    arrms(i,*)=abs(vf/sqrt(vslow2))
    arrma(i,*)=abs(vf/sqrt(vap2))
    arrmf(i,*)=abs(vf/sqrt(vfast2))
    endif

;stop
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;Plasma transisitons
if (species ne 'neutral') or (ds.fl_pip eq 0) then begin
    vslowpre =abs(vf(ipre)/sqrt(vslow2(ipre)))
    valfpre  =abs(vf(ipre)/sqrt(vap2(ipre)))
    vfastpre =abs(vf(ipre)/sqrt(vfast2(ipre)))

    tpre=4
;    if vfastpre gt 1.0 then tpre=1 else $
;    if (valfpre gt 1.0) and (vfastpre lt 1.0) then tpre=2 else $
;    if (vslowpre gt 1.0) and (valfpre lt 1.0) then tpre=3 

    diffmin=1000.0
    if vfastpre gt (1.0-shocktol) then begin 
        tpre=1 
        if (vfastpre-1.0) lt shocktol then diffmin=abs(vfastpre-1.0)
    endif

    if (valfpre gt (1.0-shocktol)) and (vfastpre lt (1.0+shocktol)) then begin 
        if tpre eq 4 then begin
            tpre=2 
            if (valfpre-1.0) lt shocktol then diffmin=abs(valfpre-1.0)
        endif else begin
            alerr=abs(valfpre-1.0)
            if (alerr lt shocktol) and (alerr lt diffmin) then diffmin=alerr & tpre=2
        endelse
    endif

    if (vslowpre gt (1.0-shocktol)) and (valfpre lt (1.0+shocktol)) then begin 
        if tpre eq 4 then begin
            tpre=3 
            if (vslowpre-1.0) lt shocktol then diffmin=abs(vslowpre-1.0)
        endif else begin
            alerr=abs(vslowpre-1.0)
            if (alerr lt shocktol) and (alerr lt diffmin) then diffmin=alerr & tpre=3
        endelse
    endif

;    vslowpos =mean(abs(vf(avecyl+4:ipos)/sqrt(vslow2(avecyl+4:ipos))))
;    valfpos  =mean(abs(vf(avecyl+4:ipos)/sqrt(vap2(avecyl+4:ipos))))
;    vfastpos =mean(abs(vf(avecyl+4:ipos)/sqrt(vfast2(avecyl+4:ipos))))
    vslowpos =abs(vf(ipos)/sqrt(vslow2(ipos)))
    valfpos  =abs(vf(ipos)/sqrt(vap2(ipos)))
    vfastpos =abs(vf(ipos)/sqrt(vfast2(ipos)))

    tpos=4
;    if vfastpos gt 1.0 then tpos=1 else $
;    if (valfpos gt 1.0) and (vfastpos lt 1.0) then tpos=2 else $
;    if (vslowpos gt 1.0) and (valfpos lt 1.0) then tpos=3 

    diffmin=1000.0
    if vfastpos gt (1.0-shocktol) then begin 
        tpos=1 
        if (vfastpos-1.0) lt shocktol then diffmin=abs(vfastpos-1.0)
    endif

    if (valfpos gt (1.0-shocktol)) and (vfastpos lt (1.0+shocktol)) then begin 
        if tpos eq 4 then begin
            tpos=2 
            if (valfpos-1.0) lt shocktol then diffmin=abs(valfpos-1.0)
        endif else begin
            alerr=abs(valfpos-1.0)
            if (alerr lt shocktol) and (alerr lt diffmin) then diffmin=alerr & tpos=2
        endelse
    endif

    if (vslowpos gt (1.0-shocktol)) and (valfpos lt (1.0+shocktol)) then begin 
        if tpos eq 4 then begin
            tpos=3 
            if (vslowpos-1.0) lt shocktol then diffmin=abs(vslowpos-1.0)
        endif else begin
            alerr=abs(vslowpos-1.0)
            if (alerr lt shocktol) and (alerr lt diffmin) then diffmin=alerr & tpos=3
        endelse
    endif
endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;Neutral transitions
if (species eq 'neutral') and (ds.fl_pip eq 1) then begin

    sonicm=abs(vf/sqrt(cs2))

    tpos=2
    if sonicm(ipos) gt 1 then tpos = 1

    tpre=2
    if sonicm(ipre) gt 1 then tpre = 1

endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Check for maximum density gradient along the normal
    rescode='Not local maximum density gradient'
    a=max(abs(deriv(ronorm)),b)
    if b eq avecyl then begin

    rescode='No transition satisfied'

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;Compressible transitions

;Plasma transisitons
if (species ne 'neutral') or (ds.fl_pip eq 0) then begin
    if (tpre eq 1) and (tpos eq 2) then begin
        fastx=[fastx,col(i)]
        fasty=[fasty,row(i)]
	if ndim eq 3 then fastz=[fastz,zrow(i)]
        if data_save eq 1 then arrtype(i)=6
        rescode='fast'
    endif

    if (tpre eq 3) and (tpos eq 4) then begin
        slowx=[slowx,col(i)]
        slowy=[slowy,row(i)]
	if ndim eq 3 then slowz=[slowz,zrow(i)]
        if data_save eq 1 then arrtype(i)=1
        rescode='slow'
    endif 

    if (tpre eq 1) and (tpos eq 3) and (bpar(ipos)*bpar(ipre) LT 0) then begin
        int1x=[int1x,col(i)]
        int1y=[int1y,row(i)]
	if ndim eq 3 then int1z=[int1z,zrow(i)]
        if data_save eq 1 then arrtype(i)=4
        rescode='int'
    endif  

    if (tpre eq 1) and (tpos eq 4) and (bpar(ipos)*bpar(ipre) LT 0) then begin
        int2x=[int2x,col(i)]
        int2y=[int2y,row(i)]
	if ndim eq 3 then int2z=[int2z,zrow(i)]
        if data_save eq 1 then arrtype(i)=3
        rescode='int'
    endif  
 

    if (tpre eq 2) and (tpos eq 3) and (bpar(ipos)*bpar(ipre) LT 0) then begin
        int3x=[int3x,col(i)]
        int3y=[int3y,row(i)]
	if ndim eq 3 then int3z=[int3z,zrow(i)]
        if data_save eq 1 then arrtype(i)=5
        rescode='int'
    endif     

    if (tpre eq 2) and (tpos eq 4) and (bpar(ipos)*bpar(ipre) LT 0) then begin
        int4x=[int4x,col(i)]
        int4y=[int4y,row(i)]
	if ndim eq 3 then int4z=[int4z,zrow(i)]
        if data_save eq 1 then arrtype(i)=2
        rescode='int'
    endif 
endif

;Neutral transitions
if (species eq 'neutral') and (ds.fl_pip eq 1) then begin
    if (tpre eq 1) and (tpos eq 2) then begin
        slowx=[slowx,col(i)]
        slowy=[slowy,row(i)]
        if data_save eq 1 then arrtype(i)=1
        rescode='sonic'
    endif
endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;Expansive transitions (non-shock?)

;    if (tpre eq 2) and (tpos eq 1) then begin
;        fastx=[fastx,col(i)]
;        fasty=[fasty,row(i)]
;        arrtype(i)=6
;    endif

;    if (tpre eq 4) and (tpos eq 3) then begin
;        slowx=[slowx,col(i)]
;        slowy=[slowy,row(i)]
;        arrtype(i)=1
;    endif    

;    if (tpre eq 3) and (tpos eq 1) then begin
;        int1x=[int1x,col(i)]
;        int1y=[int1y,row(i)]
;        arrtype(i)=4
;    endif  

;    if (tpre eq 4) and (tpos eq 1) then begin
;        int2x=[int2x,col(i)]
;        int2y=[int2y,row(i)]
;        arrtype(i)=3
;    endif  

;    if (tpre eq 3) and (tpos eq 2) then begin
;        int3x=[int3x,col(i)]
;        int3y=[int3y,row(i)]
;        arrtype(i)=5
;    endif 

;    if (tpre eq 4) and (tpos eq 2) then begin
;        int4x=[int4x,col(i)]
;        int4y=[int4y,row(i)]
;        arrtype(i)=2
;    endif    

    endif

;    print,i,(n_elements(col)-1),rescode
endfor

if n_elements(slowx) eq 0 then slowx=[0]
if n_elements(slowy) eq 0 then slowy=[0]
if n_elements(slowz) eq 0 then slowz=[0]

if n_elements(fastx) eq 0 then fastx=[0]
if n_elements(fasty) eq 0 then fasty=[0]
if n_elements(fastz) eq 0 then fastz=[0]

if n_elements(int1x) eq 0 then int1x=[0]
if n_elements(int1y) eq 0 then int1y=[0]
if n_elements(int1z) eq 0 then int1z=[0]

if n_elements(int2x) eq 0 then int2x=[0]
if n_elements(int2y) eq 0 then int2y=[0]
if n_elements(int2z) eq 0 then int2z=[0]

if n_elements(int3x) eq 0 then int3x=[0]
if n_elements(int3y) eq 0 then int3y=[0]
if n_elements(int3z) eq 0 then int3z=[0]

if n_elements(int4x) eq 0 then int4x=[0]
if n_elements(int4y) eq 0 then int4y=[0]
if n_elements(int4z) eq 0 then int4z=[0]

if ndim eq 2 then $
shocks=create_struct(['x','y','slowx','slowy','fastx','fasty',$
'int1x','int1y','int2x','int2y','int3x','int3y','int4x','int4y'],$
x,y,slowx,slowy,fastx,fasty,int1x,int1y,int2x,int2y,int3x,int3y,int4x,int4y)

if ndim eq 3 then $
shocks=create_struct(['x','y','z','slowx','slowy','slowz','fastx','fasty','fastz',$
'int1x','int1y','int1z','int2x','int2y','int2z','int3x','int3y','int3z','int4x','int4y','int4z'],$
x,y,z,slowx,slowy,slowz,fastx,fasty,fastz,int1x,int1y,int1z,int2x,int2y,int2z,$
int3x,int3y,int3z,int4x,int4y,int4z)


;p=plot3d(slowx,slowy,slowz,'.r',xr=[0,n_elements(x)],yr=[0,n_elements(y)],zr=[0,n_elements(z)])


if data_save eq 1 then begin
save,arrronorm,arrvxnorm,arrvynorm,arrbxnorm,arrbynorm,$
arrprnorm,arrvpar,arrvperp,arrvslow2,arrvalf2,arrvfast2,$
arrvs,arrtype,col,row,filename='shockid_v7_data.dat'
;Order the perminant arrays
arrronorms=arrronorm(sort(arrtype),*)
arrvxnorms=arrvxnorm(sort(arrtype),*)
arrvynorms=arrvynorm(sort(arrtype),*)
arrbxnorms=arrbxnorm(sort(arrtype),*)
arrbynorms=arrbynorm(sort(arrtype),*)
arrprnorms=arrprnorm(sort(arrtype),*)
arrvpars=arrvpar(sort(arrtype),*)
arrvperps=arrvperp(sort(arrtype),*)
arrbpars=arrbpar(sort(arrtype),*)
arrbperps=arrbperp(sort(arrtype),*)
arrvslow2s=arrvslow2(sort(arrtype),*)
arrvfast2s=arrvfast2(sort(arrtype),*)
arrvalf2s=arrvalf2(sort(arrtype),*)
arrvss=arrvs(sort(arrtype),*)
arrtypes=arrtype(sort(arrtype))
endif

;p=plot(1.0*col/(n_elements(x)-1.0),1.0*row/(n_elements(y)-1.0),/overplot,'.')

close,/all

END
