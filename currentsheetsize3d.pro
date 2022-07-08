;Get current sheet sizes

;fname='/media/snow/datadrive/OZdata/MHD_OZ_256/' & tread=10
fname='/media/snow/datadrive/OZdata/MHD_OZ_3D/' & tread=5
;fname='/media/snow/test2/OZdata/MHD_OZ_3D/'
;fname='/media/snow/datadrive/OZdata/MHD_OZ_256/' & tread=100
;fname='../../PIPoz/OZ_MHD_1024/' & tread=15
;fname='/media/snow/datadrive/OZdata/MHD_OZ_16384/Data' & tread=100


resolve_routine,'recid3d_fun'

recid3d_fun,cur,fname=fname,tread=tread

rdmpi,ds,datapath=fname,time_step=tread,/current,var=['bx','by','bz']

jz=ds.j_z
nx=n_elements(ds.x)
ny=n_elements(ds.y)
nz=n_elements(ds.z)
jz=jz(6:nx-7,6:ny-7,6:nz-7)

w=window(dim=[600,600])
;c=image(jz,rgb_table=0,/buffer)

;Arrays for widths and lengths
curw=[]
curl=[]

;i=2
for i=3,cur.count*3+1,3 do begin
print,i/3,cur.count
if n_elements(cur.(i)) GT 5 then begin
x=cur.(i)
y=cur.(i+1)
z=cur.(i+2)

;centre around zero
x0=x-mean(x)
y0=y-mean(y)
z0=z-mean(z)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;fit the data
;Polyfit isn't great for this.
;Need custom routine I think...
;a=poly_fit(x0,y0,1)
a1=indgen(101)/100.0
xp=max(x0)
xn=min(x0)
;xf=a1*(xp-xn)+xn ;fit point to the current sheet

;yf=a(0)*xf^0+a(1)*xf^1;+a(2)*xf^2+a(3)*xf^3+a(4)*xf^4

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

p=plot3d(x0+mean(x),y0+mean(y),z0+mean(z),'r.',/overplot)
;p=plot(xf+mean(x),yf+mean(y),'b',/overplot)

lp=0

if lp eq 1 then begin
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;Current sheet width
	;Based on direction of maximum gradient to jmax/2
	gradx=dblarr(n_elements(x),n_elements(y))
	grady=dblarr(n_elements(x),n_elements(y))
	gradz=dblarr(n_elements(x),n_elements(y))

	egx=N_elements(ds.x)-1
	egy=N_elements(ds.y)-1
	margin=6
	smthfac=3
	avecyl=20
	gradx=(smooth(ds.j_z(margin+1:egx-margin+1,margin:egy-margin),smthfac)-smooth(ds.j_z(margin-1:egx-margin-1,margin:egy-margin),smthfac))/(ds.x(10)-ds.x(8))
	grady=(smooth(ds.j_z(margin:egx-margin,margin+1:egy-margin+1),smthfac)-smooth(ds.j_z(margin:egx-margin,margin-1:egy-margin-1),smthfac))/(ds.y(10)-ds.y(8))
	gradz=0.0*dblarr(n_elements(gradx(*,0)),n_elements(gradx(0,*)))
	gradmag=sqrt(gradx^2+grady^2+gradz^2)

	ji=jz(cur.(i),cur.(i+1))
	jmax=max(ji,jml)


	;find shock normal based on current gradient
	normx=gradx(x(jml),y(jml))/gradmag(x(jml),y(jml))
	normy=grady(x(jml),y(jml))/gradmag(x(jml),y(jml))

	normlinex=normx*(indgen(2*avecyl+1)-avecyl)
	normliney=normy*(indgen(2*avecyl+1)-avecyl)  

	p=plot(normlinex+x(jml),normliney+y(jml),/overplot)

	tempx=indgen(2*avecyl+1)-avecyl+x(jml);mean(x)
	tempy=indgen(2*avecyl+1)-avecyl+y(jml);mean(y)
	tempx2=indgen(2*avecyl+1)-avecyl
	tempy2=indgen(2*avecyl+1)-avecyl

	;use periodic BC to fix negative values
	for ii=0,2*avecyl do begin
		if tempx(ii) lt 0 then tempx(ii)=N_elements(jz(*,0))+tempx(ii)
		if tempy(ii) lt 0 then tempy(ii)=N_elements(jz(0,*))+tempy(ii)
	endfor

	;use periodic BC to fix outside grid values
	for ii=0,2*avecyl do begin
		if tempx(ii) gt n_elements(jz(*,0))-1 then tempx(ii)=tempx(ii)-N_elements(jz(*,0))
		if tempy(ii) gt n_elements(jz(0,*))-1 then tempy(ii)=tempy(ii)-N_elements(jz(0,*))
	endfor

	jznorm=shocknormvals(jz,tempx,tempy,tempx2,tempy2,avecyl,normx,normy)

	lstep=sqrt((tempx(1)-tempx(0))^2+(tempy(1)-tempy(0))^2)
	lnorm=lstep*(indgen(2*avecyl+1)-avecyl)

	wmax=interpol(lnorm(avecyl:avecyl*2),jznorm(avecyl:avecyl*2),jmax/2.0)
	wmin=interpol(lnorm(avecyl:0:-1),jznorm(avecyl:0:-1),jmax/2.0)
	wl=wmax-wmin

	curw=[curw,wl]
endif

endif

endfor

p.save,'currentsheetwidth3dMHDt50.pdf'

END
