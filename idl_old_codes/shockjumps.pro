;Test the MHD shock jump solutions

ax1=[]
ax2=[]

;;;MHD solution
gamma=5.0/3.0
theta0=3.14/4.0
beta0=0.1

w=window(dim=[970,600])

carr=[]

for i=0.0,3.5,0.001 do begin

	ax2=[ax2,i]

	a=i*((gamma-1.0)/gamma *((gamma+1.0)/(gamma-1.0)-tan(theta0)^2.0)*(i-1.0)^2.0+tan(theta0)^2.0 * ((gamma-1.0)/gamma *i -1.0)*(i-2.0))-beta0/(cos(theta0)^2)*(i-1.0)^2.0
	b=(gamma-1.0)/gamma *(i-1.0)^2.0/(cos(theta0)^2)-i*(tan(theta0)^2)*((gamma-1.0)/gamma *i -1.0)

	ax1=[ax1,a/b]	

endfor

;carr=dblarr(3,N_elements(ax1))
;for i=0, N_elements(ax1) -1 do begin
;
;    if (ax1(i) LT ax2(i)) then carr(*,i)=[0,0,0]
;    if ((ax1(i) GT ax2(i)) and (ax1(i) LT 1.0 )) then carr(*,i)=[0,0,255]
;
;endfor

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;calculate the speeds
;Find downstream slow speed
usi=1
ax1dif=ax1-ax2
while ax1dif(usi) < 0 do begin
	usi=usi+1
endwhile
slowdown=ax2(usi)

;Find downstream fast speed
usi=N_elements(ax1)-10
ax1dif=ax1-ax2
while (ax1dif(usi) > 0) do begin
	usi=usi-1
endwhile
fastdown=ax2(usi)


;Find upstream slow speed
usi=1
ax1grad=deriv(ax2,ax1)
while ax1grad(usi-1)*ax1grad(usi+1) > 0 do begin
	usi=usi+1
endwhile
slowup=ax2(usi)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


p=plot(ax2,ax1,name='MHD',thick=4,/overplot,xr=[0,2],yr=[0,2],xtitle='$A^{d2}$', ytitle='$A^{u2}$',dim=[970,600],position=[0.1,0.1,0.9,0.9], font_size=20)
pl=plot([0:4],[0:4],/overplot)

axslow=ax1
axslow[where((ax1 GT 1.0) or (ax2 GT 1.0) or (ax2 GT ax1))]=1.0/0.0
pslow=plot(ax2,axslow,/overplot,thick=4,'r',name='Slow')

axfast=ax1
axfast[where((ax1 LT 1.0) or (ax2 GT ax1))]=1.0/0.0
pfast=plot(ax2,axfast,/overplot,thick=4,'b',name='Fast')

axint=ax1
axint[where((ax1 LT 1.0) or (ax2 GT 1.0))]=1.0/0.0 ;all intermediate shocks

axint1=axint
axint1[where(ax1 LT fastdown) ]=1.0/0.0
axint1[where(ax2 GT slowup) ]=1.0/0.0
pint1=plot(ax2,axint1,/overplot,thick=4,'m',name='Intermediate (1-4)')

axint2=axint
axint2[where(ax1 LT fastdown) ]=1.0/0.0
axint2[where(ax2 LT slowup) ]=1.0/0.0
pint2=plot(ax2,axint2,/overplot,thick=4,'y',name='Intermediate (1-3)')

axint3=axint
axint3[where(ax1 GT fastdown) ]=1.0/0.0
axint3[where(ax2 LT slowup) ]=1.0/0.0
pint3=plot(ax2,axint3,/overplot,thick=4,'c',name='Intermediate (2-3)')

axint4=axint
axint4[where(ax1 GT fastdown) ]=1.0/0.0
axint4[where(ax2 GT slowup) ]=1.0/0.0
pint4=plot(ax2,axint4,/overplot,thick=4,'g',name='Intermediate (2-4)')


pf=plot([0,4],[1,1]*fastdown,/overplot,name='Downstream fast','--')
ps=plot([0,4],[1,1]*slowdown,/overplot,name='Downstream slow','--')
ps=plot([1,1]*slowup,[0,4],/overplot,name='Upstream slow','--')
ps=plot([0,4],[1,1],/overplot,name='Alfven','--')
;p3=plot([1.0],[1.0],'o',thick=4,/overplot,name='Rotational Discontiunity')
;p3.sym_size=1.5
;p3.sym_thick=4.0

l=legend(target=[p,pslow,pfast,pint1,pint2,pint3,pint4])
l.font_size=14
l.position=[0.85,0.49]
stop

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

ax1=[]
ax2=[]

theta0=3.14/8.0
for i=0.0,3.5,0.001 do begin

	ax2=[ax2,i]

	a=i*((gamma-1.0)/gamma *((gamma+1.0)/(gamma-1.0)-tan(theta0)^2.0)*(i-1.0)^2.0+tan(theta0)^2.0 * ((gamma-1.0)/gamma *i -1.0)*(i-2.0))-beta0/(cos(theta0)^2)*(i-1.0)^2.0
	b=(gamma-1.0)/gamma *(i-1.0)^2.0/(cos(theta0)^2)-i*(tan(theta0)^2)*((gamma-1.0)/gamma *i -1.0)

	ax1=[ax1,a/b]	

endfor

p=plot(ax2,ax1,name='MHD',thick=4,/overplot,'--')

pl=plot([0:4],[0:4],/overplot)

end
