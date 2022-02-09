PRO plotCIRCLE, rad, windx=WINDX, ylim=YLIM, xlim=XLIM, oplot=OPLOT, linestyle=linestyle, help=print_help, isotropic=isotropic
;
;+
; NAME:
;	plotCIRCLE
;
; PURPOSE:
;	This procedure plots a circle with given radius:
;
; CALLING SEQUENCE:
;
;	plotcircle, rad
;
; INPUTS:
;	rad:	Radius of circle.
;	
; KEYWORD PARAMETERS:
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
;      help:   Set this keyword (flag) to show the documentation of the
;               procedure
;
; EXAMPLE:
;       plotcircle, 10., windx=1. xlim=[0.,11.], ylim=[0.,11.]
;       plotcircle, 5., /oplot
;-
;
;-----------------------IF HELP IS NEEDED-------------------------------
;
IF(KEYWORD_SET(print_help)) THEN BEGIN
   doc_library, 'plotCIRCLE'
   return
ENDIF
;
IF (N_ELEMENTS(rad) EQ 0) THEN BEGIN
   PRINT, 'SYNTAX: '
   PRINT, 'plotCIRCLE, rad'
   PRINT, 'with rad=float , e.g. plotCIRCLE, 10. '
   RETURN
ENDIF
;
;-----------------------------------------------------------------------
;
nphi=101
;
phi=findgen(nphi)*2.d0*!pi/float(nphi-1)
;
x=rad*sin(phi)
y=rad*cos(phi)
;
;---------------------CALCULATE PLOT RANGES-----------------------------
;
IF(NOT KEYWORD_SET(ylim)) THEN BEGIN
   ylim=[MIN(y), MAX(y)]
ENDIF
;
IF(NOT KEYWORD_SET(xlim)) THEN BEGIN
   xlim=[MIN(x), MAX(x)]
ENDIF
;
;------------------------TITLES-----------------------------------------
;
IF(NOT KEYWORD_SET(WINDX)) THEN BEGIN
      windx=0
ENDIF

;check if overplot is used
IF(NOT KEYWORD_SET(OPLOT)) THEN BEGIN
   window, windx
   device, decomposed=0
   loadct, 0
ENDIF
;
;-----------------------------------------------------------------------
;
IF(KEYWORD_SET(OPLOT)) THEN BEGIN
   OPLOT, x, y, linestyle=linestyle
ENDIF ELSE BEGIN
   PLOT, x, y, $
      linestyle=linestyle, $
      yrange=ylim, $
      xrange=xlim, $
         /XSTYLE, /YSTYLE, $
      isotropic=isotropic
ENDELSE
;
;
;
END

