pro get_param2,path,info
  result=path+"config.txt"
  openr,unit0,result,/get_lun
  ;read setting part----------------------
;  infile=path+"setting.txt"
 ; openr,unit,infile,/get_lun
 ; nline=file_lines(infile)
  eqs=''
  endflag=0
  tmp=""
  info={}
  while tmp ne "ENDSETTING" do begin
;  for i=0,nline-1 do begin
     readf,unit0,tmp
     tmp2=(strsplit(tmp,":",/extract))
     key=strcompress(tmp2[0],/remove_all)
     if(n_elements(tmp2)eq 1) then continue
;     if(n_elements(tmp2) gt 1) then begin
;        key=strcompress((strsplit(tmp2[1],";",/extract))[0],/remove_all)
;     endif else begin
;        continue
;     endelse     
     case key of 
        'flag_eqs':begin
           eqs=strcompress(tmp2[1],/remove_all)           
           info=create_struct(info,'eqs',eqs)
        end
        'fl_eqs':begin
           eqs=strcompress(tmp2[1],/remove_all)
           info=create_struct(info,'eqs',eqs)
        end
        'ndim':begin
           in=long(strcompress(tmp2[1],/remove_all))
           info=create_struct(info,key,in)
        end
        'safety':begin
           in=double(strcompress(tmp2[1],/remove_all))
           info=create_struct(info,key,in)
        end
        'flag_ini':begin
;           in=long(strcompress(tmp2[1],/remove_all))
           in=string(strcompress(tmp2[1],/remove_all))
           info=create_struct(info,key,in)
        end
        'flag_sch':begin
           in=long(strcompress(tmp2[1],/remove_all))
           info=create_struct(info,key,in)
        end
        'dtout':begin
           in=double(strcompress(tmp2[1],/remove_all))
           info=create_struct(info,key,in)
        end
        'tend':begin
           in=double(strcompress(tmp2[1],/remove_all))
           info=create_struct(info,key,in)
        end
        'beta':begin
           in=double(strcompress(tmp2[1],/remove_all))
           info=create_struct(info,key,in)
        end
        'flag_resi':begin
           in=long(strcompress(tmp2[1],/remove_all))
           info=create_struct(info,key,in)
        end
        'flag_IR':begin
           in=long(strcompress(tmp2[1],/remove_all))
           info=create_struct(info,key,in)
        end
        'flag_IR_type':begin
           in=long(strcompress(tmp2[1],/remove_all))
           info=create_struct(info,key,in)
        end
        'flag_grav':begin
           in=long(strcompress(tmp2[1],/remove_all))
           info=create_struct(info,key,in)
	end
        'flag_visc':begin
           in=long(strcompress(tmp2[1],/remove_all))
           info=create_struct(info,key,in)
        end
        'flag_amb':begin
           in=long(strcompress(tmp2[1],/remove_all))
           info=create_struct(info,key,in)
        end
        'eta_0':begin
           in=double(strcompress(tmp2[1],/remove_all))
           info=create_struct(info,key,in)
        end
        'flag_Hall':begin
           in=long(strcompress(tmp2[1],/remove_all))
           info=create_struct(info,key,in)
        end
        'et_Hall_0':begin
           in=double(strcompress(tmp2[1],/remove_all))
           info=create_struct(info,key,in)
        end
        'flag_divb':begin
           in=long(strcompress(tmp2[1],/remove_all))
           info=create_struct(info,key,in)
        end
        'n_fraction':begin
           in=double(strcompress(tmp2[1],/remove_all))
           info=create_struct(info,key,in)
        end
        'col':begin
           in=double(strcompress(tmp2[1],/remove_all))
           info=create_struct(info,key,in)
        end
        'flag_pip_imp':begin
           in=long(strcompress(tmp2[1],/remove_all))
           info=create_struct(info,key,in)
        end
        'flag_artvis':begin
           in=long(strcompress(tmp2[1],/remove_all))
           info=create_struct(info,key,in)
        end
        'db_clean':begin
           in=double(strcompress(tmp2[1],/remove_all))
           info=create_struct(info,key,in)
        end
        'theta':begin
           in=double(strcompress(tmp2[1],/remove_all))
           info=create_struct(info,key,in)
        end
        'debug_direction':begin
           in=long(strcompress(tmp2[1],/remove_all))
           info=create_struct(info,key,in)
        end
        's_order':begin
           in=long(strcompress(tmp2[1],/remove_all))
           info=create_struct(info,key,in)
        end
        'T_norm':begin
           in=long(strcompress(tmp2[1],/remove_all))
           info=create_struct(info,key,in)
        end
        else:
     endcase
  end
  




  ;---------------------------------------
  flag_mhd=(flag_pip=0)
  readf,unit0,flag_mhd,flag_pip
  nvar_m=(nvar_h=0)
  readf,unit0,nvar_h,nvar_m
  grids=lonarr(3)
  readf,unit0,grids
  ix=grids[0] &  jx=grids[1] &  kx=grids[2]
  margins=intarr(3)
  readf,unit0,margins
;  tmp=0
;  readf,unit0,tmp  & readf,unit0,tmp  
;  nt=tmp
  gm=0.0
  readf,unit0,gm
  flag_bnd=intarr(6)
  readf,unit0,flag_bnd

  spawn,"pwd",a_path
  outdir=a_path+path
  if flag_mhd eq 1 then begin
     files=timestamp_filter(file_search(path+"*bx.dac.0000"))
     nt=n_elements(files)
  endif else begin
     files=timestamp_filter(file_search(path+"*ro_n.dac.0000"))
     nt=n_elements(files)     
  endelse
  
  info=create_struct(info,"flag_mhd",flag_mhd)
  info=create_struct(info,"flag_pip",flag_pip)
  info=create_struct(info,"nvar_h",nvar_h)
  info=create_struct(info,"nvar_m",nvar_m)
  info=create_struct(info,"ix",ix)
  info=create_struct(info,"jx",jx)
  info=create_struct(info,"kx",kx)
  info=create_struct(info,"margin",margins)
  info=create_struct(info,"nt",nt)
  info=create_struct(info,"outdir",outdir)
  info=create_struct(info,"gm",gm)
  info=create_struct(info,"flag_bnd",flag_bnd)

;  info2={flag_mhd:flag_mhd,flag_pip:flag_pip,nvar_h:nvar_h,nvar_m:nvar_m, $
;  ix:ix,jx:jx,kx:kx,margin:margins,nt:nt,outdir:outdir,$
;       gm:gm,flag_bnd:flag_bnd}
  
  key="ixjxkx"
  in=string(ix-2*margins[0])+","+string(jx-2*margins[1]) $
     +","+string(kx-2*margins[2])
  info=create_struct(info,key,strcompress(in))
  

  if eof(unit0) ne 1 then begin
     mpi_domains=intarr(3)
     readf,unit0,mpi_domains     
     info=create_struct(info,"mpi_domains",mpi_domains)
;     ind_ix=(where(tag_names(info) eq "IX"))[0]
  endif

  free_lun,unit0


  


;  eqs=strcompress((strsplit(eqs,':',/extract))[0],/remove_all)

end




