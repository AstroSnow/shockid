;Schematic picture

;Hau sonnerup sol
ax1=[]
ax2=[]

theta0=3.14/8.0
gamma=5.0/3.0
beta0=0.1

for i=0.0,3.5,0.001 do begin

	ax2=[ax2,i]

	a=i*((gamma-1.0)/gamma *((gamma+1.0)/(gamma-1.0)-tan(theta0)^2.0)*(i-1.0)^2.0+tan(theta0)^2.0 * ((gamma-1.0)/gamma *i -1.0)*(i-2.0))-beta0/(cos(theta0)^2)*(i-1.0)^2.0
	b=(gamma-1.0)/gamma *(i-1.0)^2.0/(cos(theta0)^2)-i*(tan(theta0)^2)*((gamma-1.0)/gamma *i -1.0)

	ax1=[ax1,a/b]	

endfor

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;Find switch-off shock
soy=1.0
sox=interpol(ax2,ax1,soy)

;Compressible ratio
sor=soy/sox

;Magnetic field
sobyu=0.8
sobyd=sobyu*sor*(soy-1.0)/(soy-sor)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;Intermediate shock
isy=1.1
isx=interpol(ax2,ax1,isy)

;Compressible ratio
isr=isy/isx

;Magnetic field
isbyu=0.8
isbyd=isbyu*isr*(isy-1.0)/(isy-isr)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;populate arrays
xarr=[-20:20:0.1]
yarr=0.5*(tanh(xarr)+1.0)

soa=(yarr)*(soy-sox)+sox
sob=(yarr)*(sobyu-sobyd)+sobyd

isa=(yarr)*(isy-isx)+isx
isb=(yarr)*(isbyu-isbyd)+isbyd

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;Plotting 

psoa=plot(xarr,soa,yr=[-0.2,1.4],thick=2,name='Switch-off',font_size=14, dimension=[970,600],position=[0.1,0.1,0.9,0.9],xtitle='x',ytitle='Squared Alfven Mach number')
pisa=plot(xarr,isa,/overplot,'--',thick=2,name='2-4 Intermediate')
pa=psoa.axes
pa[3].show=1
pa[3].color='b'
pa[3].title='By'
;pa[3].hide=1

psob=plot(xarr,sob,/overplot,'b',thick=2)
pisb=plot(xarr,isb,/overplot,'b--',thick=2)
pb=psoa.axes
;pb[1].hide=1
;pb[3].show=1

poly1=polygon([[-3,-3],[-3,2],[3,2],[3,-3]],/fill_background,fill_color="light steel blue",color='k',/Data)
poly1.pattern_spacing=10
poly1.pattern_orientation=45

t=text(-17,1.2,'Downstream',/Data,/overplot,font_size=14)
t=text(7,1.2,'Upstream',/Data,/overplot,font_size=14)

l=legend(target=[psoa,pisa])
l.font_size=14
l.position=[0.85,0.3]






END
