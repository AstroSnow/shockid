function shockplane3d,normx,normy,normz
;Find and interpolate the data along the plane perpendicular to the shock norm
planearr=dblarr(3,3,3) 

const=normx^2+normy^2+normz^2

;2D perpendicular direction
;perpx=(-normy/normx)/sqrt(normy^2/normx^2+1.0)
;perpy=1.0/sqrt(normy^2/normx^2+1.0)

;
v1x=normy
v1y=-normx
v2x=-normz
v2z=normx

;if (abs(normy) le 1.0e-8) & (abs(normz) le 1.0e-8) then begin
;    v2z=1.0
;end

px=[-1:1]                                        
py=[-1:1]                              
for j=0,2 do begin                                
    for k=0,2 do begin                                
        planearr(0,j,k)=(v1x*px(j)+v2x*py(k))
        planearr(1,j,k)=v1y*px(j)
        planearr(2,j,k)=v2z*py(k)
    end
end

return,planearr

END
