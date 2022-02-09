pro test_interpolation2d, oname=oname, windx=windx
;
;plots a dummy function as a function of radius,
;and interpolated values on upwind and downwind points
;
readcol, 'test_interpolation2da.dat', r0, f0
readcol, 'test_interpolation2db.dat', format='f,f,a,a,i', r, f, direction, face, mask
;save, /variables, filename='test_interpolation2d.idldat'
;direction specifies if upwind or downwind, face specifies which plane is hit
;(xy-plane, xz-plane, yz-plane)
;

;restore, 'test_interpolation2d.idldat'
;
ru=r(where(direction eq 'u'))
fu=f(where(direction eq 'u'))
rd=r(where(direction eq 'd'))
fd=f(where(direction eq 'd'))

r_s=r(where(face eq 's'))
f_s=f(where(face eq 's'))
r_xy=r(where(face eq 'xy'))
f_xy=f(where(face eq 'xy'))
r_xz=r(where(face eq 'xz'))
f_xz=f(where(face eq 'xz'))
r_yz=r(where(face eq 'yz'))
f_yz=f(where(face eq 'yz'))
;
;-----------------------------------------------------------------------
;
if(not keyword_set(windx)) then windx=0
;
if(keyword_set(oname)) then begin
   oname='ps_files/test_interpolation2da.ps'
   print, "writing output to: ", oname
   set_plot,'ps'
   device,file=oname, xsize=19., ysize=26.7, xoffset=1., yoffset=1. , $
    decomposed=0, color=1, bits_per_pixel=8
endif else begin
   window, windx
   device, decomposed=0
   windx=windx+1
endelse
;
loadct_lev, ci_black, ci_blue, ci_red, ci_cyan, ci_green, ci_magenta, ci_yellow, ci_white
;
plot, r0, f0, $
      /ylog, /ys,  $
      xrange=[0.9,1.1], $
      yrange=[0.90,1.01], $
      xtitle=textoidl('r'), $
      ytitle=textoidl('f(r)'), $
      charsize=2., $
      psym=3
oplot, ru, fu, $
      psym=3, $
      color=ci_red
oplot, rd, fd, $
      psym=3, $
      color=ci_green
;
lstr1='local points'
lstr2='upwind points'
lstr3='downwind points'
;
if(keyword_set(oname)) then begin
   legend, [lstr1, lstr2, lstr3], $
           psym=[3,3,3], $
           linestyle=[0,0,0], $
           color=[ci_black, ci_red, ci_green], $
           textcolor=ci_black, $
           /right_legend, $
           charsize=1.5
endif else begin
   legend, [lstr1, lstr2, lstr3], $
           psym=[3,3,3], $
           linestyle=[0,0,0], $
           color=[ci_white,ci_red,ci_green], $
           textcolor=ci_white, $
           /right_legend, $
           charsize=1.5
endelse



if keyword_set(oname) then begin
   device, /close
   set_plot,'x'
endif
;
;-----------------------------------------------------------------------
;
if(keyword_set(oname)) then begin
   oname='ps_files/test_interpolation2db.ps'
   print, "writing output to: ", oname
   set_plot,'ps'
   device,file=oname, xsize=19., ysize=26.7, xoffset=1., yoffset=1. , $
    decomposed=0, color=1, bits_per_pixel=8
endif else begin
   window, windx
   device, decomposed=0
   windx=windx+1
endelse
;
loadct_lev, ci_black, ci_blue, ci_red, ci_cyan, ci_green, ci_magenta, ci_yellow, ci_white
;
plot, r0, f0, $
      /ylog, /ys,  $
      xrange=[0.9,1.1], $
      yrange=[0.90,1.01], $
      xtitle=textoidl('r'), $
      ytitle=textoidl('f(r)'), $
      charsize=2., $
      psym=3
oplot, r_xy, f_xy, $
      psym=3, $
      color=ci_blue
stop
oplot, r_xz, f_xz, $
      psym=3, $
      color=ci_red
stop
oplot, r_yz, f_yz, $
      psym=3, $
      color=ci_green
stop
oplot, r_s, f_s, $
      psym=3, $
      color=ci_magenta


;
lstr1='local points'
lstr2='points in x-y plane'
lstr3='points in x-z plane'
lstr4='points in y-z plane'
lstr5='points on surface'
;
if(keyword_set(oname)) then begin
   legend, [lstr1, lstr2, lstr3, lstr4, lstr5], $
           psym=[3,3,3,3,3], $
           linestyle=[0,0,0,0,0], $
           color=[ci_black, ci_blue, ci_red, ci_green, ci_magenta], $
           textcolor=ci_black, $
           /right_legend, $
           charsize=1.5
endif else begin
   legend, [lstr1, lstr2, lstr3, lstr4, lstr5], $
           psym=[3,3,3,3,3], $
           linestyle=[0,0,0,0,0], $
           color=[ci_white,ci_blue,ci_red,ci_green, ci_magenta], $
           textcolor=ci_white, $
           /right_legend, $
           charsize=1.5
endelse



if keyword_set(oname) then begin
   device, /close
   set_plot,'x'
endif

end
