pro recid_fun,cur,fname=fname,tread=tread

rdmpi,ds,datapath=fname,time_step=tread,/current

;E=ds.vx_p*ds.by-ds.vy_p*ds.bx

;c=contour(e,/fill,n_levels=101)

jz=ds.j_z

nx=n_elements(ds.x)
ny=n_elements(ds.y)
jz=jz(6:nx-7,6:ny-7)
;vx=(ds.vx_p(6:nx-7,6:ny-7))*sqrt(ds.ro_p(6:nx-7,6:ny-7))
;vy=(ds.vy_p(6:nx-7,6:ny-7))*sqrt(ds.ro_p(6:nx-7,6:ny-7))
;gridx=ds.x(6:nx-7)
;gridy=ds.y(6:ny-7)


jmin=max(abs(jz))/2.0

;Find all unsorted current
index=where(abs(jz) gt jmin)
s = SIZE(jz)
ncol = s(1)
col = index MOD ncol
row = index / ncol

;Sort the current by largest first
js=dblarr(n_elements(col))
for i=0,n_elements(col)-1 do begin
    js(i)=jz(col(i),row(i))
endfor
jso=sort(1.0/abs(js))
col=col(jso)
row=row(jso)


;Create initial structure
cur=create_struct('col',col)
cur=create_struct(cur,'row',row)

;c=image(jz,rgb_table=0)
;cb=colorbar(target=[c],orientation=1) 
;p=plot(col,row,/overplot,'.k')

i=0
;for i=0,n_elements(col)-1 do begin
while n_elements(col) ge 2 do begin

;Distance from test cell
d=double((col-col(0)))^2+double((row-row(0)))^2
index=where(d le 2.01) ;find adacent cells
oldl=-1
while (n_elements(index) ne oldl) do begin
    oldl=n_elements(index)    
    for j=0,n_elements(index)-1 do begin
        d=double((col-col(index(j))))^2+double((row-row(index(j))))^2
        index=[index,where(d le 2.01)] ;find adacent cells
        ;print,n_elements(index),j    
    endfor
;p=plot(col(index),row(index),'ob')
;stop
;index=index(uniq(index(sort(index)))) 
index=index(sort(index))
index=index(uniq(index))
;print,index
endwhile

if n_elements(index) gt 2 then cur=create_struct(cur,'a'+strtrim(i,1)+'x',col(index))
if n_elements(index) gt 2 then cur=create_struct(cur,'a'+strtrim(i,1)+'y',row(index))


;remove detected points from the array
col2=col(index)
row2=row(index)
;p=plot(col2,row2,/overplot,'ko-')
;wait,3
;p.hide=1
for k=0,n_elements(index)-1 do begin
    a=where((col eq col2(k)) and (row eq row2(k)))
    if n_elements(col) gt 1 then remove,a,col,row
endfor 

if n_elements(index) gt 2 then i=i+1
if n_elements(index) gt 2 then print,i,n_elements(col)

;p=plot(col,row,/overplot,'ko')
;wait,1
;p.hide=1

endwhile


cur=create_struct(cur,'count',i)

close,/all

END
