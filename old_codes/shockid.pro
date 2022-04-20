;Shock identification

rdmpi,ds,datapath='../MHD_test',time_step=100

margin=6

egx=N_elements(ds.x)-1
egy=N_elements(ds.y)-1

x=ds.x(margin:egx-margin)
y=ds.y(margin:egy-margin)
ro=ds.ro_p(margin:egx-margin,margin:egy-margin)
vx=ds.vx_p(margin:egx-margin,margin:egy-margin)
vy=ds.vy_p(margin:egx-margin,margin:egy-margin)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Div(v)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

divv=dblarr(n_elements(x),n_elements(y))
divvs=dblarr(n_elements(x),n_elements(y))

for i=0,n_elements(x)-1 do begin
for j=0,n_elements(y)-1 do begin
    divv(i,j)=(ds.vx_p(margin+i+1,margin+j)-ds.vx_p(margin+i-1,margin+j))/(ds.x(margin+i+1)-ds.x(margin+i-1)) $
+(ds.vy_p(margin+i,margin+j+1)-ds.vy_p(margin+i,margin+j+1))/(ds.y(margin+j+1)-ds.y(margin+j-1))

;if (divv(i,j) LT 0) then divvs(i,j)=divv(i,j)

;Filter out weak shocks
r=max([ds.ro_p(i+1,j)/ds.ro_p(i-1,j),ds.ro_p(i-1,j)/ds.ro_p(i+1,j),ds.ro_p(i,j+1)/ds.ro_p(i,j-1),ds.ro_p(i,j-1)/ds.ro_p(i,j+1)])
if ((divv(i,j) LT 0) and (r ge 1.03^2)) then divvs(i,j)=divv(i,j)

endfor
endfor

c=contour(alog10(-divvs),x,y,n_levels=101,/fill)
c.rgb_table=16
cb=colorbar(target=c)

c.save,'shockid_divvs.jpg'

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Find shocks and their centres
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

divvs2=divvs
snum=[]
sx=[]
sy=[]

sc=0

xpt=0
xmt=0
ypt=0
ymt=0

ncx=0
ncy=0

i3=0
j3=0

for i=2,n_elements(x)-3 do begin
for j=2,n_elements(y)-3 do begin
    if divvs(i,j) LT 0 then begin
        i2=i
        j2=j
        snum=[snum,sc]
        sx=[sx,i2]
        sy=[sy,j2]
        divvs2(i2,j2)=0 ;remove the old point

        ;define a test block
        tr=divvs(max([0,i2-1]):min([i2+1,n_elements(x)-2]),max([0,j2-1]):min([j2+1,n_elements(y)-2]))
        while (min(tr) lt 0) and (i2-1+i3+ncx lt n_elements(x)-1) and (j2-1+j3+ncy lt n_elements(y)-1) do begin
            mindv=0
            ncxt=0
            ncyt=0
            for i3=0,n_elements(tr(*,0))-1 do begin
                for j3=0,n_elements(tr(0,*))-1 do begin
                    if tr(i3,j3) LT 0 then begin
                            snum=[snum,sc]
                            sx=[sx,max([0,i2-1+i3+ncx])]
                            sy=[sy,max([0,j2-1+j3+ncy])]
                            if tr(i3,j3) lt mindv then begin
                                ncxt=i3
                                ncyt=j3
                            endif
                            divvs(max([0,i2-1+i3+ncx]):min([i2+1+i3+ncx,n_elements(x)-1])$
,max([0,j2-1+j3+ncy]):min([j2+1+j3+ncy,n_elements(y)-1]))=0
                    endif
                endfor
             endfor
            ncx=min([ncxt,n_elements(x)-2])
            ncy=min([ncyt,n_elements(y)-2])
tr=divvs(max([0,i2-1+ncx]):min([i2+1+ncx,n_elements(x)-1]),max([0,j2-1+ncy]):min([j2+1+ncy,n_elements(y)-1]))
         endwhile
        sc=sc+1
;        print,sc
    endif                   
                     
endfor
endfor

;This loop doesnt work. Supposed to remove small (lt 10 cells) shocks and plot but does't work
;i=0
;ei=0
;for s=0,max(snum) do begin
;    while snum(ei) eq s do begin
;         ei=ei+1
;    endwhile
;    if ei-i gt 10 then begin
;    sxt=sx(i:ei-1)/1024.0
;    syt=sy(i:ei-1)/1024.0
;    p=plot(sx(i:ei-1),sy(i:ei-1),/overplot) 
;    print,i,ei,s,max(snum)
;    endif
;    i=ei
;endfor

;p=plot(sx/1024.0,sy/1024.0,/overplot,'.')

c.save,'shockid_shockid.jpg'


END
