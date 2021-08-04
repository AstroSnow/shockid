function shocknormvals3D, A, tempxyz,avecyl,normx,normy,normz

    tempa=dblarr(2*avecyl+1,2*avecyl+1,2*avecyl+1)

    tempx=tempxyz(0,*)
    tempy=tempxyz(1,*)
    tempz=tempxyz(2,*)
    tempx2=tempxyz(3,*)
    tempy2=tempxyz(4,*)
    tempz2=tempxyz(5,*)

    ;populate the array (density)
    for ii=0,2*avecyl do begin
        for jj=0,2*avecyl do begin
            for kk=0,2*avecyl do begin
                tempa(ii,jj,kk)=A(tempx(ii),tempy(jj),tempz(jj))
            endfor
        endfor
    endfor

    ;lind the normal line
    normlinex=normx*(indgen(2*avecyl+1)-avecyl)
    normliney=normy*(indgen(2*avecyl+1)-avecyl)
    normlinez=normz*(indgen(2*avecyl+1)-avecyl)   
print,normlinex
print,normliney
print,normlinez
stop
;    c=contour(tempa,tempx2,tempy2,/fill)
;    p=plot(normlinex,normliney,/overplot)
;stop
    ;Interpolate the values normal to the shock
    anorm=interp2d(tempa,tempx2,tempy2,normlinex,normliney)

    return,anorm
END
