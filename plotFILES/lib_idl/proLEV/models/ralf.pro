function ralf, beq, rstar, mdot, vinf, help=print_help
;
;+
; NAME:
;       ralf
;
; PURPOSE:
;       This function calculates the Alfven-radius for given input
;
; CALLING SEQUENCE:
;       Result = ralf(beq, rstar, mdot, vinf)
;
; INPUTS:
;       beq:    Equatorial magnetic field strength in Gauss (=B_pole/2)
;       rstar:  Radius in r_sun.
;       mdot:   Mass loss rate (for non-magnetic model) in M_sun/yr
;       vinf:   Terminal velocity in cm/s
;
; KEYWORD PARAMETERS:
;       help:   Set this keyword (flag) tho show the documentation of this function
;
; EXAMPLE:
;       Result = beq(2450., 14.5, 5.d-6, 2700.d5)
;-
;
;-----------------------------------------------------------------------
;
if(keyword_set(print_help)) then begin
   doc_library, 'ralf'
   return, 0
endif
;
;-----------------------------------------------------------------------
;
;convert units
mdot_cgs=mdot*!msu/!yr
rstar_cgs=rstar*!rsu

eta_star=beq^2 * rstar_cgs^2 / mdot_cgs /vinf
;
;low-field domain (eta_star << 1, udDoula2008)
ralf0=eta_star^0.25d0
;
;high-field domain (eta_star >> 1, udDoula2008)
ralf1=0.3d0+eta_star^0.25d0
;
;fit formula (udDoula2008)
ralf2=0.3d0+(eta_star+0.25d0)^0.25d0
;
;or find correct ralf by iteration (ralf^4-ralf^3)=eta_star for dipole)

print, ralf0, ralf1, ralf2
;
return, ralf2

end
