function get_mpi_vars,key,mpi_x,mpi_y,mpi_z,margin,ix_m,jx_m,kx_m,$
                      time_step,datapath
  ix=(ix_m-2*margin[0])*mpi_x+2*margin[0]
  jx=(jx_m-2*margin[1])*mpi_y+2*margin[1]
  kx=(kx_m-2*margin[2])*mpi_z+2*margin[2]
  n_read=n_elements(time_step)
  output=fltarr(ix,jx,kx,n_read)  
  for tt=0,n_read-1 do begin     
     files=file_search(datapath+string(time_step[tt],form="(i4.4)")+key $
                       +".dac.*")
     tmp=!NULL
     mpi_read,tmp,files,mpi_x,mpi_y,mpi_z,margin,ix_m,jx_m,kx_m
     output[*,*,*,tt]=tmp     
  endfor
  return,output
end
