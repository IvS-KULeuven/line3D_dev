function lineforce, trad, help=print_help
;
;+
; NAME:
;       lineforce
;
; PURPOSE:
;       This function calculates the ratio of lineforce for a a different
;       radiation temperature only (the rest shall be equal), from the
;       tabulated values in Puls et al 2000
;
; CALLING SEQUENCE:
;       Result = lineforce(trad1, trad2)
;
; INPUTS:
;       trad1:   radiation temperature of model 1 in K
;       trad1:   radiation temperature of model 2 in K
;
; KEYWORD PARAMETERS:
;       help:   Set this keyword (flag) tho show the documentation of this function
;
; EXAMPLE:
;       Result = lineforce(35.d0, 30.d0)
;-
;
;-----------------------------------------------------------------------
;
if(keyword_set(print_help)) then begin
   doc_library, 'lineforce'
   return, 0
endif
;
;-----------------------------------------------------------------------
;
;define CAK-parameter for a range of effective temperatures (Puls et al 2000,
;table 2)
teff=[1.d4,2.d4,3.d4,4.d4,5.d4]
alpha=[0.44,0.58,0.64,0.67,0.66]
qbar=[915.d0,1497.d0,2498.d0,1954.d0,1939.d0]
qzero=[14505.d0,5171.d0,3630.d0,1778.d0,2260.d0]
;
;find CAK-parameter for the given radiation temperatures



rhoc = sqrt(rad-1.d0+mu^2)*sqrt(1.d0+3.d0*mu^2) / $
      sqrt(mu^2+delta^2/rad^2)/(4.d0*rad-3.d0+3.d0*mu^2)/rad^2

return, rhoc

end
