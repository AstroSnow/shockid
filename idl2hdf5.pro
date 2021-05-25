; IDL routine to save HDF5 files 

;Input directory
fname='../MHD_test'

;Define save directory MUST ALREADY EXIST
savdir=fname+'/HDF5/'

;Variables to save
outvar=['ro_p','vx_p','vy_p','vz_p','bx','by','bz','divv','grad','v','b','pr_p','grad3']

;Read in the data
rdmpi,ds,datapath=fname,time_step=100

;Define simulation grid
margin=6

egx=N_elements(ds.x)-1
egy=N_elements(ds.y)-1

x=ds.x(margin:egx-margin)
y=ds.y(margin:egy-margin)
ro=ds.ro_p(margin:egx-margin,margin:egy-margin)
vx=ds.vx_p(margin:egx-margin,margin:egy-margin)
vy=ds.vy_p(margin:egx-margin,margin:egy-margin)
vz=ds.vz_p(margin:egx-margin,margin:egy-margin)
bx=ds.bx(margin:egx-margin,margin:egy-margin)
by=ds.by(margin:egx-margin,margin:egy-margin)
bz=ds.bz(margin:egx-margin,margin:egy-margin)
pr_p=ds.pr_p(margin:egx-margin,margin:egy-margin)

v=dblarr(3,n_elements(x),n_elements(y))
v(0,*,*)=vx
v(1,*,*)=vy
v(2,*,*)=vz

b=dblarr(3,n_elements(x),n_elements(y))
b(0,*,*)=bx
b(1,*,*)=by
b(2,*,*)=bz

;2D divergence of velocity
divv=dblarr(n_elements(x),n_elements(y))
for i=0,n_elements(x)-1 do begin
for j=0,n_elements(y)-1 do begin
    divv(i,j)=(ds.vx_p(margin+i+1,margin+j)-ds.vx_p(margin+i-1,margin+j))/(ds.x(margin+i+1)-ds.x(margin+i-1)) $
+(ds.vy_p(margin+i,margin+j+1)-ds.vy_p(margin+i,margin+j+1))/(ds.y(margin+j+1)-ds.y(margin+j-1))
endfor
endfor

;Maximum gradient of density
grad=dblarr(n_elements(x),n_elements(y))
for i=0,n_elements(x)-1 do begin
for j=0,n_elements(y)-1 do begin
    grad(i,j)=sqrt(((ds.ro_p(margin+i+1,margin+j)-ds.ro_p(margin+i-1,margin+j))/(ds.x(margin+i+1)-ds.x(margin+i-1)))^2 $
+((ds.ro_p(margin+i,margin+j+1)-ds.ro_p(margin+i,margin+j+1))/(ds.y(margin+j+1)-ds.y(margin+j-1)))^2)
endfor
endfor

;Maximum gradient of density
grad3=dblarr(3,n_elements(x),n_elements(y))
for i=0,n_elements(x)-1 do begin
for j=0,n_elements(y)-1 do begin
    grad3(0,i,j)=(ds.ro_p(margin+i+1,margin+j)-ds.ro_p(margin+i-1,margin+j))/(ds.x(margin+i+1)-ds.x(margin+i-1))
    grad3(1,i,j)=(ds.ro_p(margin+i,margin+j+1)-ds.ro_p(margin+i,margin+j+1))/(ds.y(margin+j+1)-ds.y(margin+j-1))
    grad3(2,i,j)=0.0
endfor
endfor


for i=0,n_elements(outvar)-1 do begin

    if outvar(i) eq 'ro_p' then data = ro
    if outvar(i) eq 'vx_p' then data = vx
    if outvar(i) eq 'vy_p' then data = vy
    if outvar(i) eq 'vz_p' then data = vz
    if outvar(i) eq 'bx' then data = bx
    if outvar(i) eq 'by' then data = by
    if outvar(i) eq 'bz' then data = bz
    if outvar(i) eq 'divv' then data = divv
    if outvar(i) eq 'grad' then data = grad
    if outvar(i) eq 'v' then data = v
    if outvar(i) eq 'b' then data = b
    if outvar(i) eq 'pr_p' then data = pr_p
    if outvar(i) eq 'grad3' then data = grad3

    ;; get data type and space, needed to create the dataset

    file = savdir+outvar(i)+'.h5'
    fid = H5F_CREATE(file)

    datatype_id = H5T_IDL_CREATE(data)
    dataspace_id = H5S_CREATE_SIMPLE(size(data,/DIMENSIONS))

    ;; create dataset in the output file

    dataset_id = H5D_CREATE(fid,$

    outvar(i),datatype_id,dataspace_id)

    ;; write data to dataset

    H5D_WRITE,dataset_id,data

    ;; close all open identifiers

    H5D_CLOSE,dataset_id

    H5S_CLOSE,dataspace_id

    H5T_CLOSE,datatype_id

    H5F_CLOSE,fid

endfor

END
