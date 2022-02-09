function interpol_ypfct, xminus, xplus, yminus, yplus, xp, help=print_help
;
;-----------------------------------------------------------------------
;
;+
; NAME:
;	interpol_ypfct
;
; PURPOSE:
;	This function interpolates linearly between two points.
;       
;
; CALLING SEQUENCE:
;	fp=interpol_ypfct(x1,x2,y1,y2,xp)
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
;	fp:     function value at unknown point
;
; EXAMPLE:
;	interpol_ypfct(1.,2.,1.,4.,1.5)
;-
;
;-------------------------OUTPUT IF HELP NEEDED-------------------------
;
if(keyword_set(print_help)) then begin
   doc_library, 'interpol_ypfct'
   stop
endif
;
;-----------------------------------------------------------------------
;
interpol_yp = yplus + (yplus-yminus)/(xplus-xminus) * (xp-xplus)
;
return, interpol_yp
;
;
end function interpol_ypfct

