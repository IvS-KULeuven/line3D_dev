pro plotellipse, a, b, windx=windx, ylim=ylim, xlim=xlim, oplot=oplot, help=print_help
;
;+
; name:
;	plotellipse
;
; purpose:
;	this procedure plots an ellipse given by x^2/a^2 + y^2/b^2 = 1
;
; calling sequence:
;
;	plotellipse, a, b
;
; inputs:
;	a:   semi-major axis
;	b:   semi-minor axis
;	
; keyword parameters:
;       windx:  set this keyword to an integer defining the window, to which
;               output is plotted
;
;       ylim:   set this keyword to a 2-d array, which sets the yrange
;
;       xlim:   set this keyword to a 2-d array, which sets the xrange
;
;       oplot:  set this keyword (flag), to overplot in the current window
;               this keyword is not allowed, when oname is set
;
;      help:   set this keyword (flag) to show the documentation of the
;               procedure
;
; example:
;       plotellipse, 1.5, 1., windx=1. xlim=[0.,11.], ylim=[0.,11.]
;       plotellipse, 1.5, 1., /oplot
;-
;
;-----------------------if help is needed-------------------------------
;
if(keyword_set(print_help)) then begin
   doc_library, 'plotellipse'
   return
endif
;
if(n_elements(a) eq 0 or n_elements(b) eq 0) then begin
   doc_library, 'plotellipse'
   return
endif
;
;-----------------------------------------------------------------------
;
nphi=101
;
phi=findgen(nphi)*2.d0*!pi/float(nphi-1)
rad=a*b/sqrt(b^2*cos(phi)^2 + a^2*sin(phi)^2)
;
x=rad*cos(phi)
y=rad*sin(phi)
;
;---------------------calculate plot ranges-----------------------------
;
if(not keyword_set(ylim)) then begin
   ylim=[min(y), max(y)]
endif
;
if(not keyword_set(xlim)) then begin
   xlim=[min(x), max(x)]
endif
;
;------------------------titles-----------------------------------------
;
if(not keyword_set(windx)) then begin
      windx=0
endif

;check if overplot is used
if(not keyword_set(oplot)) then begin
   window, windx
   device, decomposed=0
   loadct, 0
endif
;
;-----------------------------------------------------------------------
;
if(keyword_set(oplot)) then begin
   oplot, x, y
endif else begin
   plot, x, y, $
      yrange=ylim, $
      xrange=xlim, $
      /xstyle, /ystyle
endelse
;
;
;
end

