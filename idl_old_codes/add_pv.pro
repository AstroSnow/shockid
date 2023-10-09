pro add_pv,pv,vars,flag_mhd,var_names,mpi_x,mpi_y,mpi_z,margin,ix_m,jx_m,kx_m,$
                      time_step,datapath
  nv=n_elements(vars)
  gm=pv.info.gm
  if flag_mhd eq 1 then begin
     vnames=["ro_p","en_p","mx_p","my_p","mz_p","bx","by","bz"]
  endif else begin
     vnames=["ro_n","en_n","mx_n","my_n","mz_n"] 
  endelse



  for v=0,nv-1 do begin
     tmp=get_mpi_vars(vars[v],mpi_x,mpi_y,mpi_z,margin,ix_m,jx_m,kx_m,$
                      time_step,datapath)
     command=vnames[v]+"=tmp"
     err=execute(command)
  endfor

  vars[1]="VX"+strmid(vars[1],2)
  vars[2]="VY"+strmid(vars[2],2)
  vars[3]="VZ"+strmid(vars[3],2)
  vars[4]="PR"+strmid(vars[4],2)
  var_names=[var_names,strupcase(vars)]  

  if flag_mhd eq 1 then begin
     vx_p=mx_p/ro_p
     vy_p=my_p/ro_p
     vz_p=mz_p/ro_p
     pr_p=(gm-1.0)*(en_p-0.5*ro_p*(vx_p^2+vy_p^2+vz_p^2) $
                    -0.5*(bx^2+by^2+bz^2))     

     pv=create_struct(pv,vars,reform(ro_p),reform(vx_p),reform(vy_p), $
                      reform(vz_p),reform(pr_p), $
                      reform(bx),reform(by),reform(bz)) 

     if pv.info.flag_divb eq 1 and pv.info.ndim ge 2  then begin
        tmp=get_mpi_vars("ps",mpi_x,mpi_y,mpi_z,margin,ix_m,jx_m,kx_m,$
                         time_step,datapath)
        pv=create_struct(pv,["ps"],reform(tmp))
     endif
  endif else begin
     vx_n=mx_n/ro_n
     vy_n=my_n/ro_n
     vz_n=mz_n/ro_n
     pr_n=(gm-1.0)*(en_n-0.5*ro_n*(vx_n^2+vy_n^2+vz_n^2))
     pv=create_struct(pv,vars, $
                      reform(ro_n),reform(vx_n),reform(vy_n),reform(vz_n),reform(pr_n))
   endelse

end

