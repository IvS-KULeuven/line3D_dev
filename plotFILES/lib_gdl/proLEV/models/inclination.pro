function inclination, vrot, vsini, help=print_help
;
;+
; NAME:
;       inclination
;
; PURPOSE:
;       This function calculates inclination angle i for a given vrot and vsini
;
; CALLING SEQUENCE:
;       Result = inclination(vrot, vsini)
;
; INPUTS:
;       vrot:   rotational velocity at equator in arbitrary units
;       vsini:  vrot*sin(inclination) in same units
;
; KEYWORD PARAMETERS:
;       help:   Set this keyword (flag) tho show the documentation of this function
;
; EXAMPLE:
;       Result = inclination(400.d0, 200.d0)
;-
;
;-----------------------------------------------------------------------
;
if(keyword_set(print_help)) then begin
   doc_library, 'inclination'
   return, 0
endif
;
;-----------------------------------------------------------------------
;
sini = vsini/vrot
if(n_elements(sini) gt 1) then begin
   indx=where(sini gt 1.d0)
   if(n_elements(indx) eq 0) then begin
      if(indx ne -1) then sini=1.d0
   endif else begin
      sini(indx)=(indx+1.d0)/(indx+1.d0)
   endelse
endif else begin
   sini=min([1.d0,sini])
endelse
;
incl = asin(sini)
return, incl

end
