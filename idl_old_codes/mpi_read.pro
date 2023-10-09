pro mpi_read,output,files,mpi_x,mpi_y,mpi_z,$
             margin,ix_m,jx_m,kx_m,time_step=time_step
  if(n_elements(output) eq 0) then begin
     ix=(ix_m-2*margin[0])*mpi_x+2*margin[0]
     jx=(jx_m-2*margin[1])*mpi_y+2*margin[1]
     kx=(kx_m-2*margin[2])*mpi_z+2*margin[2]
     output=fltarr(ix,jx,kx)
  endif
  if(n_elements(time_step) eq 0) then time_step=-1
;   help,output
;   print,time_step

  for nf=0,n_elements(files)-1 do begin
     dacget3s,files[nf],var,narg=time_step
     ii=nf mod mpi_x
     jj=(nf/(mpi_x)) mod mpi_y
     kk=nf/(mpi_x*mpi_y)
     x0=ii*(ix_m-2*margin[0])
     y0=jj*(jx_m-2*margin[1])
     z0=kk*(kx_m-2*margin[2])
     output[x0:x0+ix_m-1,y0:y0+jx_m-1,z0:z0+kx_m-1,*]=(var)
  endfor

end
