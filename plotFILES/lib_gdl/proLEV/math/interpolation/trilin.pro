pro trilin, xout, yout, zout, $
            x1, x2, y1, y2, z1, z2, $
            vala, valb, valc, vald, vale, valf, valg, valh, $
            llogx, llogy, llogz, llogf, yinterp
;
if(llogf eq 1) then begin
;
;prepare input values for log-* interpolation
   dum_vala=alog10(vala)
   dum_valb=alog10(valb)
   dum_valc=alog10(valc)
   dum_vald=alog10(vald)
   dum_vale=alog10(vale)
   dum_valf=alog10(valf)
   dum_valg=alog10(valg)
   dum_valh=alog10(valh)
;
   if(llogz eq 1) then begin
;log-log interpolation in z-direction
      dum_z1=alog10(abs(z1))
      dum_z2=alog10(abs(z2))
      dum_zout=alog10(abs(zout))
   endif else begin
;log-lin interpolation in z-direction
      dum_z1=z1
      dum_z2=z2
      dum_zout=zout
   endelse
;
   if(llogx eq 1) then begin
;log-log-interpolation in x-direction
      dum_x1=alog10(abs(x1))
      dum_x2=alog10(abs(x2))
      dum_xout=alog10(abs(xout))
   endif else begin
;log-lin-interpolation in x-direction
      dum_x1=x1
      dum_x2=x2
      dum_xout=xout
   endelse
;
   if(llogy eq 1) then begin
;log-log-interpolation in y-direction
      dum_y1=alog10(abs(y1))
      dum_y2=alog10(abs(y2))
      dum_yout=alog10(abs(yout))
   endif else begin
;log-lin-interpolation in y-direction
      dum_y1=y1
      dum_y2=y2
      dum_yout=yout
   endelse
;
;perform interpolation in z-direction
   interpol_yp, dum_z2, dum_z1, dum_vald, dum_valb, dum_zout, yvalue_s0
   interpol_yp, dum_z2, dum_z1, dum_valc, dum_vala, dum_zout, yvalue_s1
   interpol_yp, dum_z2, dum_z1, dum_valh, dum_valf, dum_zout, yvalue_n0
   interpol_yp, dum_z2, dum_z1, dum_valg, dum_vale, dum_zout, yvalue_n1
;
;perform interpolation in x-direction
   interpol_yp, dum_x2, dum_x1, yvalue_n0, yvalue_n1, dum_xout, soln1
   interpol_yp, dum_x2, dum_x1, yvalue_s0, yvalue_s1, dum_xout, sols1
;
;perform interpolation in y-direction
   interpol_yp, dum_y2, dum_y1, soln1, sols1, dum_yout, yinterp
   yinterp = 10.d0^yinterp
;
endif else begin
;lin-lin-interpolation
;
   interpol_yp, z2, z1, vald, valb, zout, yvalue_s0
   interpol_yp, z2, z1, valc, vala, zout, yvalue_s1
   interpol_yp, z2, z1, valh, valf, zout, yvalue_n0
   interpol_yp, z2, z1, valg, vale, zout, yvalue_n1
;
   interpol_yp, x2, x1, yvalue_n0, yvalue_n1, xout, soln1
   interpol_yp, x2, x1, yvalue_s0, yvalue_s1, xout, sols1
;
   interpol_yp, y2, y1, soln1, sols1, yout, yinterp
;
endelse
;
;
end
