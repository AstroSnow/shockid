;; calculate az in 2D by using bx, by, dx, and dy data.
;; first, integrate bx along the y axis
;; dirct : 1 or -1

pro cal_az2d,bx,by,dx,dy,az,itype=itype,dirct=dirct,margin=margin,filename=filename

help,bx,by,dx,dy


;filename=''
print,'cal_az2d mod'

if not keyword_set(itype) then read,'itype (1 or 2): ',itype
if not keyword_set(margin) then read,'margin: ',margin
if not keyword_set(dirct) then read,'dirct: ',dirct
if not keyword_set(filename) then begin
   filename=''
   read,'filename: ',filename
endif
if n_elements(margin) ne 1 then begin
   print,'margin...'
   stop
endif

s=size(bx)
ix=s(1)
jx=s(2)
nx=s(3)

az=fltarr(ix,jx,nx)

t0=systime(1)

tempx = fltarr(ix)
tempy = fltarr(jx)

for n=0,nx-1 do begin
    print,'======== LOOP ',n,' ========'
    bxinteg=fltarr(ix,jx) 
    byinteg=fltarr(ix,jx)

    case itype of
       1: begin ;; integrate bx(margin,j) in y
          if(dirct eq 1) then begin ;; lower left
             byinteg[margin,*] = by[margin,*,n]*dx[margin]
             byinteg[margin-1,*] = byinteg[margin,*] $
                                   - by[margin-1,*,n]*dx[margin-1] ;; one step outside the boundary
             for i=margin+1,ix-1 do begin
                byinteg[i,*] = byinteg[i-1,*] + by[i,*,n]*dx[i]
             endfor
             bxinteg[*,margin] = bx[margin,margin,n]*dy[margin]
             bxinteg[*,margin-1] = bxinteg[*,margin] $
                                   - bx[margin,margin-1,n]*dy[margin-1] ;; one step outside the boundary
             for j=margin+1,jx-1 do begin
                tempx[*] = bx[margin,j,n]*dy[j]
                bxinteg[*,j] = bxinteg[*,j-1] + tempx[*]
             endfor
          endif else if(dirct eq -1) then begin ;; lower right
             byinteg[ix-margin-1,*] = by[ix-margin-1,*,n]*dx[ix-margin-1]
             byinteg[ix-margin,*] = byinteg[ix-margin-1,*] $
                                    + by[ix-margin,*,n]*dx[ix-margin] ;; one step outside the boundary
             for i=ix-margin-2,0,-1 do begin
                byinteg[i,*] = byinteg[i+1,*] - by[i,*,n]*dx[i]
             endfor
             bxinteg[*,margin] = bx[margin,margin,n]*dy[margin]
             bxinteg[*,margin-1] = bxinteg[*,margin] $
                                   - bx[margin,margin-1,n]*dy[margin-1] ;; one step outside the boundary             
             for j=margin+1,jx-1 do begin
                tempx[*] = bx[margin,j,n]*dy[j]
                bxinteg[*,j] = bxinteg[*,j-1] + tempx[*]
             endfor
          endif else begin
             print,'???'
             stop
          endelse
       end
       2: begin ;; integrate by(i,margin) in x
          if(dirct eq 1) then begin
             byinteg[margin,*] = by[margin,margin,n]*dx[margin]
             byinteg[margin-1,*] = byinteg[margin,*] $
                                   - by[margin-1,margin,n]*dx[margin-1] ;; one step outside the boundary
             for i=margin+1,ix-1 do begin
                tempy[*] = by[i,margin,n]*dx[i]
                byinteg[i,*] = byinteg[i-1,*] + tempy[*] 
             endfor
             bxinteg[*,margin] = bx[*,margin,n]*dy[margin]
             bxinteg[*,margin-1] = bxinteg[*,margin] $
                                   - bx[*,margin-1,n]*dy[margin-1]
             for j=margin+1,jx-1 do begin
                bxinteg[*,j] = bxinteg[*,j-1] + bx[*,j,n]*dy[j]
             endfor
          endif else if (dirct eq -1) then begin
             byinteg[ix-margin-1,*] = by[ix-margin-1,margin,n]*dx[ix-margin-1]
             byinteg[ix-margin,*] = byinteg[ix-margin-1,*] $
                                   + by[ix-margin,margin,n]*dx[ix-margin] ;; one step outside the boundary
             for i=ix-margin-2,0,-1 do begin
                tempy[*] = by[i,margin,n]*dx[i]
                byinteg[i,*] = byinteg[i+1,*] - tempy[*] 
             endfor
             bxinteg[*,margin] = bx[*,margin,n]*dy[margin]
             bxinteg[*,margin-1] = bxinteg[*,margin] $
                                   - bx[*,margin-1,n]*dy[margin-1]
             for j=margin+1,jx-1 do begin
                bxinteg[*,j] = bxinteg[*,j-1] + bx[*,j,n]*dy[j]
             endfor             
          endif else begin
             print,'???'
             stop
          endelse
       end
       endcase          

    az[*,*,n] = bxinteg - byinteg
;    t1=systime(1)
;    print,'TIME:',t1-t0
endfor;; n

print,'Saving ',filename
save,az,filename=filename

end
