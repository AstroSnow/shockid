;PRO gridconvert

;Dependencies: all those from rdmpi plus
;add_cv.pro ;reads in the consered variables

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;TO DO
;1) write the config file correctly
;2) is a settings file needed?
;3) save the new data
;4) test everything
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

fnamein='/home/snow/Documents/Shock_stab/PIP/widthtests/ac1/'
tstepin=0

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;read in the initial data
datapath=fnamein
get_param2,datapath,info
gm=info.gm
cv=create_struct("info",info)
ts=info.nt
n_read=n_elements(tstepin)

margin=info.margin
eqs=info.eqs
flag_mhd=0
flag_pip=0
case eqs of 
    'HD':nvar=5
    'MHD':begin
    nvar=8
    flag_mhd=1
    end
    'PIP':begin
    nvar=13
    flag_mhd=1
    flag_pip=1
    end
end

  cv=create_struct(cv,["eqs","fl_mhd","fl_pip","fl_afr"], $
                   eqs,flag_mhd,flag_pip,flag_afr)

  tfile=file_search(datapath+"t.dac.*")
  dacget0s,tfile,t,narg=time_step
  if(eqs eq "AFR") then begin
     gridfile=file_search(datapath+"region.dac.*")
     dacget2s,gridfile,grid,narg=time_step
     cv=create_struct(cv,["grid"],grid)
  endif
  ix=info.ix
  jx=info.jx
  kx=info.kx


  if (where(tag_names(info) eq "MPI_DOMAINS"))[0] ne -1 then begin
     mpi_x=info.mpi_domains[0]
     mpi_y=info.mpi_domains[1]
     mpi_z=info.mpi_domains[2]
  endif else begin
     mpi_x=(mpi_y=(mpi_z=1))  
  endelse
  ix_m=(jx_m=(kx_m=1)) 
;  xfile=file_search(datapath+"x.dac.*")
  xfile=file_search(datapath+"x.dac."+string(indgen(mpi_x),form="(i4.4)"))
  ix_m=info.ix
  ix=margin[0]*2+mpi_x*(ix_m-margin[0]*2)
  x=findgen(ix)
  for n=0,mpi_x-1 do begin
     dacget1d,xfile[n],xc
     x0=n*(ix_m-2*margin[0])
     x[x0:x0+ix_m-1]=xc     
  endfor

  cv=create_struct(cv,"t",t,"x",x)
  ndim=cv.info.ndim
  if ndim ge 2 then begin     
     yfile=file_search(datapath+"y.dac."+string(indgen(mpi_y),form="(i4.4)"))
;       yfile=file_search(datapath+"y.dac.*")
     jx_m=info.jx
     jx=margin[1]*2+mpi_y*(jx_m-margin[1]*2)
     y=findgen(jx)
     for n=0,mpi_y-1 do begin
        dacget1d,yfile[n],yc
        y0=n*(jx_m-2*margin[1])
        y[y0:y0+jx_m-1]=yc     
     endfor
     cv=create_struct(cv,"y",y)
     if ndim ge 3 then begin
        zfile=file_search(datapath+"z.dac."+string(indgen(mpi_z),form="(i4.4)"))
        kx_m=info.kx
        kx=margin[2]*2+mpi_z*(kx_m-margin[2]*2)
        z=findgen(kx)
        for n=0,mpi_z-1 do begin
           dacget1d,zfile[n],zc
           z0=n*(kx_m-2*margin[2])
           z[z0:z0+kx_m-1]=zc     
        endfor
        cv=create_struct(cv,"z",z)
     endif
  endif
  cv=create_struct(cv,"ix",n_elements(x))
  if ndim ge 2 then jx=n_elements(y) else jx=1  
  if ndim ge 3 then kx=n_elements(z) else kx=1  
  cv=create_struct(cv,"jx",jx,"kx",kx)  

  n_cpu=mpi_x*mpi_y*mpi_z
  var_names=[]
  if time_step[0] eq -1 then begin
     time_step=indgen(n_read)
  endif
  output=dblarr(ix,jx,kx)

;Read in the conserved quantities
  if (flag_pip eq 1 or flag_mhd eq 0) then begin
     add_cv,cv,["ro_n","en_n","mx_n","my_n","mz_n"] ,0,var_names,$
            mpi_x,mpi_y,mpi_z,margin,ix_m,jx_m,kx_m,$
            time_step,datapath    
  endif  
  if (flag_mhd eq 1) then begin
     add_cv,cv,["ro_p","en_p","mx_p","my_p","mz_p","bx","by","bz"],$
            1,var_names,mpi_x,mpi_y,mpi_z,margin,ix_m,jx_m,kx_m,$
            time_step,datapath
  endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;Do something to the data
cv_new=cv
cv_new.bx=cv_new.bx+1.0

info_new=info

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;Save the config data
fnameout='./gridconverttest/'
tstepout=0

;Write a new config file
newconfigname=fnameout+'config.txt'

;NEED TO DO THIS PROPERLY
;$cp /home/snow/Documents/Shock_stab/PIP/widthtests/ac1/config.txt gridconverttest/config.txt

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;Save the new grid
 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;Save the new conservative variables


END
