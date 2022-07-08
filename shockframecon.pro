function shockframecon,roup,vxup,vyup,bxup,byup,prup,rodown,vxdown,vydown,bxdown,bydown,prdown
;Estimate the shock frame by minimising all conservative variables

;Initial guess based on conserving mass
vxs=(rodown*vxdown-roup*vxup)/(roup-rodown)
;Use initial guess of vx and electric field to estimate vys
vys=(vxup+vxs)*byup/bxup-vyup


vxsa=[0.1*vxs:10.0*vxs:0.1*vxs]
vysa=[0.1*vys:10.0*vys:0.1*vys]

sol=dblarr(n_elements(vxsa),n_elements(vysa))
solmass=dblarr(n_elements(vxsa),n_elements(vysa))
solvx=dblarr(n_elements(vxsa),n_elements(vysa))
solvy=dblarr(n_elements(vxsa),n_elements(vysa))
solen=dblarr(n_elements(vxsa),n_elements(vysa))
solel=dblarr(n_elements(vxsa),n_elements(vysa))

;gas constant
gm=5.0/3.0

;This loop is far to slow and should be written in vector form when it is fixed
for i=0,n_elements(vxsa)-1 do begin
    for j=0,n_elements(vysa)-1 do begin
        sol(i,j)=abs(roup*(vxup+vxsa(i))-rodown*(vxdown+vxsa(i))) +$;Mass
            abs((roup*(vxup+vxsa(i))^2+prup+(byup^2)/2.0)-$
            (rodown*(vxdown+vxsa(i))^2+prdown+(bydown^2)/2.0)) +$;vx momentum
            abs((roup*(vxup+vxsa(i))*vyup+vysa(j)-bxup*byup)-$
            (rodown*(vxdown+vxsa(i))*vydown+vysa(j)-bxdown*bydown))+$;vy momentum
            abs(((vxup+vxsa(i))*(gm/(gm-1.0)*prup+0.5*roup*((vxup+vxsa(i))^2+(vyup+vysa(j))^2)))-$
            ((vxdown+vxsa(i))*(gm/(gm-1.0)*prdown+0.5*rodown*((vxdown+vxsa(i))^2+(vydown+vysa(j))^2))));+$;energy
            ;abs(((vxup+vxsa(i))*byup-(vyup+vysa(j))*bxup)-$
            ;((vxdown+vxsa(i))*bydown-(vydown+vysa(j))*bxdown))

;        solmass(i,j)=abs(roup*(vxup+vxsa(i))-rodown*(vxdown+vxsa(i)));Mass
;        solvx(i,j)=abs((roup*(vxup+vxsa(i))^2+prup+(byup^2)/2.0)-$
;            (rodown*(vxdown+vxsa(i))^2+prdown+(bydown^2)/2.0));vx momentum
;        solvy(i,j)=abs((roup*(vxup+vxsa(i))*vyup+vysa(j)-bxup*byup)-$
;            (rodown*(vxdown+vxsa(i))*vydown+vysa(j)-bxdown*bydown));vy momentum
;        solen(i,j)=abs(((vxup+vxsa(i))*(gm/(gm-1.0)*prup+0.5*roup*((vxup+vxsa(i))^2+(vyup+vysa(j))^2)))-$
;            ((vxdown+vxsa(i))*(gm/(gm-1.0)*prdown+0.5*rodown*((vxdown+vxsa(i))^2+(vydown+vysa(j))^2))));energy
;        solel(i,j)=abs(((vxup+vxsa(i))*byup-(vyup+vysa(j))*bxup)-$
;            ((vxdown+vxsa(i))*bydown-(vydown+vysa(j))*bxdown))
    endfor
endfor

a=min(sol,b)
col = b MOD n_elements(vxsa)
row = b / n_elements(vxsa)

vs=vxsa(col)

;contour,sol
;stop

return,vs

END
