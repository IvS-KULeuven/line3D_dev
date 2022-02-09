function get_indx_obsdir, aindx, gindx, ngamma, help=print_help
;
;+
; NAME:
;       get_indx_obsdir
;
; PURPOSE:
;       This function calculate the index of observer-direction for
;       given index describing alpha, and total number of gamma-directions.
;
; CALLING SEQUENCE:
;       get_indx_obsdir, aindx, gindx, ngamma
;
; INPUTS:
;       aindx:  Integer describing observers angle alpha (w.r.t. z-axis)
;       gindx:  Integer describing observers angle gamma (w.r.t x-z-plane)
;       ngamma: Integer describing total number of given gamma-angles
;
; KEYWORD PARAMETERS:
;       help:   Set this keyword (flag) tho show the documentation of this procedure
;
; OUTPUTS:
;       Index for which wanted observers direction is stored.
;
; EXAMPLE:
;       get_indx_obsdir, 1, 10, 19
;
;-------------------------OUTPUT IF HELP NEEDED-------------------------
;
if(keyword_set(print_help)) then begin
   doc_library, 'contur_fluxem'
   return, 0
endif
;
;-----------------------------------------------------------------------
;
if(gindx gt ngamma) then begin
   print, 'error in get_indx_obsdir: gindx gt ngamma'
   stop
endif
;
;-----------------------------------------------------------------------
;
indx_obsdir=(aindx-1)*ngamma + gindx
return, indx_obsdir
;
;
;
end
