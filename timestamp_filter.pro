function timestamp_filter,files
  spawn,'ls -lt '+strjoin(files," "),vars
  nv=n_elements(vars)
  bef=""
  for n=0,nv-1 do begin
     new=(strsplit(vars[n]," ",/extract))[-1]
     if(n gt 0 and new gt bef) then break
     bef=new
  endfor
  n=n<nv-1  
  return, files[0:n]  
end
