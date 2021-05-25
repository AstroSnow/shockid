function shocknormvals, A, tempx,tempy,tempx2,tempy2,avecyl,normx,normy

    tempa=dblarr(2*avecyl+1,2*avecyl+1)

    ;populate the array (density)
    for ii=0,2*avecyl do begin
        for jj=0,2*avecyl do begin
            tempa(ii,jj)=A(tempx(ii),tempy(jj))
        endfor
    endfor

    ;lind the normal line
    normlinex=normx*(indgen(2*avecyl+1)-avecyl)
    normliney=normy*(indgen(2*avecyl+1)-avecyl)   

;    c=contour(tempa,tempx2,tempy2,/fill)
;    p=plot(normlinex,normliney,/overplot)
;stop
    ;Interpolate the values normal to the shock
    anorm=interp2d(tempa,tempx2,tempy2,normlinex,normliney)

    return,anorm
END
