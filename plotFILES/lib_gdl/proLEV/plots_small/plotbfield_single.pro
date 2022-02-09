pro plotBFIELD_single, windx=WINDX, nlines=NLINES, ralfven=RALFVEN, xlim=XLIM, ylim=YLIM, oplot=OPLOT, color=COLOR, help=print_help, mu_lim=mu_lim
;
;+
; NAME:
;	plotBFIELD_single
;
; PURPOSE:
;	This procedure plots magnetic dipole field lines
;
; CALLING SEQUENCE:
;
;	plotbfield
;
; INPUTS:
;
; KEYWORD PARAMETERS:
;       nlines: Set this keyword to an integer defining the number of B-field
;               lines that shall be plotted
;
;       ralfven: Set this keyword to the alfven-radius, for which B-field line
;                is plotted explicitly
;
;       windx:  Set this keyword to an integer defining the window, to which
;               output is plotted
;
;       ylim:   Set this keyword to a 2-d array, which sets the yrange
;
;       xlim:   Set this keyword to a 2-d array, which sets the xrange
;
;       oplot:  Set this keyword (flag), to overplot in the current window
;               This keyword is not allowed, when oname is set
;
;       mu_lim: Set this keyword to a maximum mu that shall be plotted for b-field-line
;
;       help:   Set this keyword (flag) to show the documentation of the
;               procedure
;
; EXAMPLE:
;       plotbfield, nlines=10, ralfven=5., windx=1. xlim=[0.,11.], ylim=[0.,11.]
;       plotbfield, /oplot
;-
;
;-----------------------IF HELP IS NEEDED-------------------------------
;
IF(KEYWORD_SET(print_help)) THEN BEGIN
   doc_library, 'plotBFIELD'
   return
ENDIF
;
;-----------------------------------------------------------------------
;
;number of field lines
if(not keyword_set(nlines)) then begin
   nlines=6
endif
;
;alfven radius (will always be plotted explicitly
if(not keyword_set(ralfven)) then begin
   ralfven=3.5
endif
;
;plot ranges
if(not keyword_set(xlim)) then begin
   xlim=[0.,5.]
endif
;
if(not keyword_set(ylim)) then begin
   ylim=[-2.5,2.5]
endif
;
;-----------------------------------------------------------------------
;
if(not keyword_set(windx)) then begin
   windx=0
endif
;
if(not keyword_set(oplot)) then begin
   ntheta=101
   theta=findgen(ntheta)*!pi/float(ntheta-1)
   x=sin(theta)
   y=cos(theta)
   window, windx
   device, decomposed=0
   loadct, 0
endif
;
if(not keyword_set(oplot)) then begin
   plot, x, y, $
         xrange=xlim, $
         yrange=ylim, $
         /ISOTROPIC
endif
;
;-----------------------------------------------------------------------
;
;overplot b-field line for alfven radius
if(keyword_set(mu_lim)) then begin
   nmu=101
   mu_alfven=sqrt(1.d0-1.d0/ralfven)
   mu_min=-1.d0*mu_alfven
   mu_max=-mu_lim
   mu=mu_min+findgen(nmu)*(mu_max-Mu_min)/(nmu-1)
   radius=(1.d0-mu*mu)/(1.d0-mu_alfven^2)
   x=sqrt(1.d0-mu*mu)*radius
   y=mu*radius   
   oplot, x, y, line=0, color=COLOR

   mu_min=mu_lim
   mu_max=mu_alfven
   mu=mu_min+findgen(nmu)*(mu_max-mu_min)/(nmu-1)
   radius=(1.d0-mu*mu)/(1.d0-mu_alfven^2)
   x=sqrt(1.d0-mu*mu)*radius
   y=mu*radius   
   oplot, x, y, line=0, color=COLOR  
endif else begin
   nmu=101
   mu_alfven=sqrt(1.d0-1.d0/ralfven)
   mu_min=-1.d0*mu_alfven
   mu_max=mu_alfven
   
   mu=mu_min+findgen(nmu)*(mu_max-Mu_min)/(nmu-1)
   radius=(1.d0-mu*mu)/(1.d0-mu_alfven^2)
   x=sqrt(1.d0-mu*mu)*radius
   y=mu*radius
;   
   oplot, x, y, line=0, color=COLOR
endelse

;
;-----------------------------------------------------------------------
;
;
;
end
