pro add_cv,pv,vars,flag_mhd,var_names,mpi_x,mpi_y,mpi_z,margin,ix_m,jx_m,kx_m,$
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
  
     pv=create_struct(pv,vnames,reform(ro_p),reform(mx_p),reform(my_p), $
                      reform(mz_p),reform(en_p), $
                      reform(bx),reform(by),reform(bz)) 

  endif else begin

     pv=create_struct(pv,vnames, $
                      reform(ro_n),reform(mx_n),reform(my_n),reform(mz_n),reform(en_n))
   endelse

end

