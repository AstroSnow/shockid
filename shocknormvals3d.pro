function shocknormvals3d, A, tempxyz,avecyl,normx,normy,normz

    tempa=dblarr(2*avecyl+1,2*avecyl+1,2*avecyl+1)

    tempx=reform(tempxyz(0,*))
    tempy=reform(tempxyz(1,*))
    tempz=reform(tempxyz(2,*))
    tempx2=reform(tempxyz(3,*))
    tempy2=reform(tempxyz(4,*))
    tempz2=reform(tempxyz(5,*))

    ;populate the array (density)
    for ii=0,2*avecyl do begin
        for jj=0,2*avecyl do begin
            for kk=0,2*avecyl do begin
                tempa(ii,jj,kk)=A(tempx(ii),tempy(jj),tempz(kk))
            endfor
        endfor
    endfor

    ;lind the normal line
    normlinex=normx*(indgen(2*avecyl+1)-avecyl)+avecyl
    normliney=normy*(indgen(2*avecyl+1)-avecyl)+avecyl
    normlinez=normz*(indgen(2*avecyl+1)-avecyl)+avecyl
;print,normlinex,normliney,normlinez
;stop
;    c=contour(tempa(*,*,0),/fill)
;    p=plot(normlinex,normliney,/overplot)
;stop
    ;Interpolate the values normal to the shock
    ;anorm=interp2d(tempa,tempx2,tempy2,normlinex,normliney)
;print,tempa
    anorm=interpolate(tempa,normlinex,normliney,normlinez)	
;print,anorm
;stop
    return,anorm
END
