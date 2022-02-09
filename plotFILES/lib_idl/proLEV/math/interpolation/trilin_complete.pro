pro trilin_complete, xout, yout, zout, $
                     x1, x2, y1, y2, z1, z2, $
                     vala, valb, valc, vald, vale, valf, valg, valh, $
                     rada, radb, radc, radd, rade, radf, radg, radh, radp, $
                     expol, llogx, llogy, llogz, llogf, lr2, yinterp
;
if(lr2 eq 1) then begin
;applying interpolation of function values * r^2
   vala=vala*rada*rada
   valb=valb*radb*radb
   valc=valc*radc*radc
   vald=vald*radd*radd
   vale=vale*rade*rade
   valf=valf*radf*radf
   valg=valg*radg*radg
   valh=valh*radh*radh
endif
;
;perform only interpolation (expol=.false.)
if (expol eq 0) then begin
   trilin, xout, yout, zout, $
           x1, x2, y1, y2, z1, z2, $
           vala, valb, valc, vald, vale, valf, valg, valh, $
           llogx, llogy, llogz, llogf, yinterp
endif else begin
;set values to zero if extrapolation is needed (or perform extrapolation)
   yinterp=0.d0
   trilin, xout, yout, zout, $
           x1, x2, y1, y2, z1, z2, $
           vala, valb, valc, vald, vale, valf, valg, valh, $
           llogx, llogy, llogz, llogf, yinterp
endelse
;
if(lr2 eq 1) then begin
   yinterp=yinterp/radp/radp
endif
;
;
;
end
