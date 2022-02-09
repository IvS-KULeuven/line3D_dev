function angle_obs_mag, incl, obl, phi, help=print_help
;
;+
; NAME:
;	angle_obs_mag
;
; PURPOSE:
;	This function calculates the angle between observers direction and
;       magnetic pole axis for an oblique rotator model at different phases
;
; CALLING SEQUENCE:
;	Result = angle_obs_mag(incl, obl, phi)
;
; INPUTS:
;	incl:  inclination in degree
;	obl:   obliquity in degree
;       phi:   phase angle in degree
;	
; KEYWORD PARAMETERS:
;	help:	Set this keyword (flag) tho show the documentation of this function
;
; EXAMPLE:
;	Result = angle_obs_mag(23.,73.,[0.,90.,180.,270.,360.])
;-
;
;-----------------------------------------------------------------------
;
if(keyword_set(print_help)) then begin
   doc_library, 'angle_obs_mag'
   return, 0
endif
;
if(n_params() lt 3) then begin
   print, 'Syntax - angle_obs_mag(inclination, obliquity, phi)'
   return, 0
endif
;
;-----------------------------------------------------------------------
;
;phi is phase angle, in terms of phase: phi=phi/2/pi
cosa = sin(incl*!pi/180.d0)*sin(obl*!pi/180.d0)*cos(phi*!pi/180.d0) + $
       cos(incl*!pi/180.d0)*cos(obl*!pi/180.d0)
;
return, acos(cosa)*180.d0/!pi

;
end
