PRO INTERPOL_YP, x1, x2, y1, y2, xp, yp, help=print_help
;
;+
; NAME:
;	INTERPOL_YP
;
; PURPOSE:
;	This procedure interpolates linearly between two points.
;       
;
; CALLING SEQUENCE:
;	INTERPOL_YP, x1, x2, y1, y2, xp, yp
;
; INPUTS:
;	x1, x2:	coordinates of given points
;	y1, y2:	function values at given points
;	xp:     coordinate of unknown point
;	
; KEYWORD PARAMETERS:
;	help:	Set this keyword (flag) tho show the documentation of this procedure
;
; OUTPUTS:
;	xp:     function value at unknown point
;
; EXAMPLE:
;	INTERPOL_YP, 1., 2., 1., 4., 1.5, yp
;-
;
;-------------------------OUTPUT IF HELP NEEDED-------------------------
;
if(keyword_set(print_help)) then begin
   doc_library, 'interpol_yp'
   return
endif
;
;-----------------------------------------------------------------------
;
yp = y2 + (y2-y1)/(x2-x1)*(xp-x2)

END
