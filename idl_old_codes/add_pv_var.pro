pro add_pv_var,pv,vars,flag_mhd,var_names,mpi_x,mpi_y,mpi_z,margin,ix_m,jx_m,kx_m,$
                      time_step,datapath
  nv=n_elements(vars)
  gm=pv.info.gm
  if flag_mhd eq 1 then begin
     vnames=["ro_p","en_p","mx_p","my_p","mz_p","bx","by","bz"]
  endif else begin
     vnames=["ro_n","en_n","mx_n","my_n","mz_n"] 
  endelse

;make array for needed components
vartarr=[]
for v=0,nv-1 do begin
	if (vars(v) eq 'bx') or (vars(v) eq 'by') or (vars(v) eq 'bz') or (vars(v) eq 'ro_p') or (vars(v) eq 'ro_n') then vartarr=[vartarr,vars(v)]
	if (vars(v) eq 'vx_p') then vartarr=[vartarr,'mx_p','ro_p']
	if (vars(v) eq 'vy_p') then vartarr=[vartarr,'my_p','ro_p']
	if (vars(v) eq 'vz_p') then vartarr=[vartarr,'mz_p','ro_p']
	if (vars(v) eq 'pr_p') then vartarr=[vartarr,'en_p','mx_p','my_p','mz_p','ro_p','bx','by','bz']

	if (vars(v) eq 'vx_n') then vartarr=[vartarr,'mx_n','ro_n']
	if (vars(v) eq 'vy_n') then vartarr=[vartarr,'my_n','ro_n']
	if (vars(v) eq 'vz_n') then vartarr=[vartarr,'mz_n','ro_n']
	if (vars(v) eq 'pr_n') then vartarr=[vartarr,'en_n','mx_n','my_n','mz_n','ro_n']
endfor
vres=vartarr[uniq(vartarr,sort(vartarr))]
;print,'to read: ', vres

  for v=0,n_elements(vres)-1 do begin
	     tmp=get_mpi_vars(vres[v],mpi_x,mpi_y,mpi_z,margin,ix_m,jx_m,kx_m,$
		              time_step,datapath)
	     command=vres(v)+"=tmp"
	     err=execute(command)
  endfor

;  vars[1]="VX"+strmid(vars[1],2)
;  vars[2]="VY"+strmid(vars[2],2)
;  vars[3]="VZ"+strmid(vars[3],2)
;  vars[4]="PR"+strmid(vars[4],2)
  var_names=[var_names,strupcase(vars)]  

for vloop=0,n_elements(vars)-1 do begin
;print,vars(vloop)

if vars(vloop) eq 'vx_n' then pv=create_struct(pv,vars(vloop),reform(mx_n/ro_n))
if vars(vloop) eq 'vy_n' then pv=create_struct(pv,vars(vloop),reform(my_n/ro_n))
if vars(vloop) eq 'vz_n' then pv=create_struct(pv,vars(vloop),reform(mz_n/ro_n))

if vars(vloop) eq 'vx_p' then pv=create_struct(pv,vars(vloop),reform(mx_p/ro_p))
if vars(vloop) eq 'vy_p' then pv=create_struct(pv,vars(vloop),reform(my_p/ro_p))
if vars(vloop) eq 'vz_p' then pv=create_struct(pv,vars(vloop),reform(mz_p/ro_p))

if vars(vloop) eq 'bx' then pv=create_struct(pv,vars(vloop),reform(bx))
if vars(vloop) eq 'by' then pv=create_struct(pv,vars(vloop),reform(by))
if vars(vloop) eq 'bz' then pv=create_struct(pv,vars(vloop),reform(bz))

if vars(vloop) eq 'ro_n' then pv=create_struct(pv,vars(vloop),reform(ro_n))
if vars(vloop) eq 'ro_p' then pv=create_struct(pv,vars(vloop),reform(ro_p))

if vars(vloop) eq 'pr_n' then pv=create_struct(pv,vars(vloop),(gm-1.0)*(reform(en_n)-0.5*reform(ro_n)*(reform(mx_n/ro_n)^2+reform(my_n/ro_n)^2+reform(mz_n/ro_n)^2)))
if vars(vloop) eq 'pr_p' then pv=create_struct(pv,vars(vloop),(gm-1.0)*(reform(en_p)-0.5*reform(ro_p)*(reform(mx_p/ro_p)^2+reform(my_p/ro_p)^2+reform(mz_p/ro_p)^2) -0.5*(reform(bx)^2+reform(by)^2+reform(bz)^2)))
endfor
;  if flag_mhd eq 1 then begin
;     vx_p=mx_p/ro_p
;     vy_p=my_p/ro_p
;     vz_p=mz_p/ro_p
;     pr_p=(gm-1.0)*(en_p-0.5*ro_p*(vx_p^2+vy_p^2+vz_p^2) $
;                    -0.5*(bx^2+by^2+bz^2))     

;     pv=create_struct(pv,vars,reform(ro_p),reform(vx_p),reform(vy_p), $
;                      reform(vz_p),reform(pr_p), $
;                      reform(bx),reform(by),reform(bz)) ;

;     if pv.info.flag_divb eq 1 and pv.info.ndim ge 2  then begin
;        tmp=get_mpi_vars("ps",mpi_x,mpi_y,mpi_z,margin,ix_m,jx_m,kx_m,$
;                         time_step,datapath)
;        pv=create_struct(pv,["ps"],reform(tmp))
;     endif
;  endif else begin
;     vx_n=mx_n/ro_n
;     vy_n=my_n/ro_n
;     vz_n=mz_n/ro_n
;     pr_n=(gm-1.0)*(en_n-0.5*ro_n*(vx_n^2+vy_n^2+vz_n^2))
;     pv=create_struct(pv,vars, $
;                      reform(ro_n),reform(vx_n),reform(vy_n),reform(vz_n),reform(pr_n))
;   endelse

end

