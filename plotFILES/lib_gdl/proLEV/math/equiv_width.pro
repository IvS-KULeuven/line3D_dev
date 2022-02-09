pro equiv_width, flux_tot, flux_cont, xarr, ew_lamb, ew_nue, ew_xobs, wave=WAVE, freq=FREQ, xobs=XOBS, lamb0=LAMB0, nue0=NUE0, vth_fiducial=VTH_FIDUCIAL, help=print_help
;
;+
; NAME:
;	equiv_width
;
; PURPOSE:
;	This procedure calculates the equivalent widths of a line-profile in
;	wavelength, frequency and xobs-regime.
;
;          W = integral of (f_tot/f_cont - 1) d(lambda,xobs,nue)
;
;          equivalent width > 0 when emission
;          equivalent width > 0 when absorption
;
;
; CALLING SEQUENCE:
;	CALL, equiv_width, flux_tot, flux_cont, xarr, ew_lamb, ew_nue, ew_xobs
;
; INPUTS:
;	flux_tot:  Total flux in the line
;	flux_cont: total continuum flux (in the absence of the line)
;       xarr:      wavelength in Angstrom, frequency in Hz or xobs in shift
;                  from line-center in doppler-units (depending on keywords)
;	
; KEYWORD PARAMETERS:
;       wave:   Set this keyword (flag) if xarr is given in
;               wavelength-regime.
;       freq:   Set this keyword (flag) if xarr is given in
;               frequency-regime.
;       xobs:   Set this keyword (flag) if xarr is given in
;               xobs-regime.
;       lamb0:  Set this keyword to the transition wavelength in Anstrom
;       freq0:  Set this keyword to the transition wavelength in Hz
;       vth_fiducial: Set this keyword to th fiducial Doppler velocity in cm/s
;	help:	Set this keyword (flag) tho show the documentation of this procedure
;
; OUTPUTS:
;       ew_lamb: Equivalent width in wavelength-regime
;       ew_freq: Equivalent width in frequency-regime
;       ew_xobs: Equivalent width in xobs-regime
;
; EXAMPLE:
;	equiv_width, ftot, fcont, lambda, ewl, ewn, ewx, /WAVE, lamb0=6562., $
;	vth_fiducial=100.
;-
;
;--------------------if help is needed----------------------------------
;
if(keyword_set(print_help)) then begin
   doc_library, 'equiv_width'
   return
endif
;
;------------------------CHECK KEYWORDS, ETC----------------------------
;
err_message=0
;
;check if arrays have correct input
if (n_elements(flux_tot) lt 2) then begin
   print, 'input error equiv_width: flux_tot needs at least two elements'
   err_message=1
endif
;
if(n_elements(flux_cont) lt 2) then begin
   print, 'input error equiv_width: flux_cont needs at least two elements'
   err_message=1
endif
;
if(n_elements(xarr) lt 2) then begin
   print, 'input error equiv_width: xarr needs at least two elements'
   err_message=1
endif
;
;check if output values are specified
if(n_elements(ew_lamb) eq 0) then begin
   print, 'input error equiv_width: ew_lamb not specified'
   err_message=1
endif
if(n_elements(ew_nue) eq 0) then begin
   print, 'input error equiv_width: ew_nue not specified'
   err_message=1
endif
if(n_elements(ew_xobs) eq 0) then begin
   print, 'input error equiv_width: ew_xobs not specified'
   err_message=1
endif
;
;check if keywords wave, freq, xobs are set correct
numb=0
if(keyword_set(wave)) then begin
   numb=numb+1
endif
if(keyword_set(freq)) then begin
   numb=numb+1
endif
if(keyword_set(xobs)) then begin
   numb=numb+1
endif
if(numb eq 0) then begin
   print, 'input error equiv_width: keyword, wave, freq or xobs is needed'
   err_message=1
endif
if(numb gt 1) then begin
   print, 'input error equiv_width: only one keyword wave, freq, xobs is allowed'
   err_message=1
endif
;
;check if keywords xnue0, lamb0 are set correct
if(not keyword_set(lamb0) and not keyword_set(nue0)) then begin
   print, 'input error equiv_width: keyword lamb0 or nue0 have to be specified'
   err_message=1
endif
;
if(not keyword_set(vth_fiducial)) then begin
   print, 'input error equiv_width: keyword vth_fiducial has to be specified'
   errr_message=1
endif
;
if(err_message eq 1) then begin
   print, 'syntax error in equiv_width'
   print, 'calling sequence:'
   print, '   equiv_width, flux_tot, flux_cont, xarr, ew_lamb, ew_nue, ew_xobs, wave=WAVE, freq=FREQ, xobs=XOBS, lamb0=LAMB0, nue0=NUE0, vth_fiducial=VTH_FIDUCIAL'
   stop
endif
;
if(keyword_set(nue0)) then begin
   lamb0=!cgs_clight/nue0/1.d-8
endif
if(keyword_set(lamb0)) then begin
   nue0=!cgs_clight/lamb0/1.d-8
endif
;
;-----------------------------------------------------------------------
;
nd1=n_elements(flux_tot)
nd2=n_elements(flux_cont)
nd3=n_elements(xarr)
if(nd1 ne nd2 or nd1 ne nd3 or nd2 ne nd3) then begin
   print, 'error equiv_width: different length of input arrays not allowed'
   stop
endif
nd=nd1
;
;integration with trapezoidal rule
equivalent_width=0.d0
for i=1, nd-1 do begin
   equivalent_width=equivalent_width + (0.5d0*(flux_tot(i-1)/flux_cont(i-1) + flux_tot(i)/flux_cont(i)) - 1.d0) * (xarr(i)-xarr(i-1))
endfor
;
;change sign if xarr is decreasing
if(xarr(1) lt xarr(0)) then begin
   equivalent_width=-1.d0*equivalent_width
endif
;
if(keyword_set(wave)) then begin
   ew_lamb = equivalent_width
   ew_nue = !cgs_clight/lamb0/lamb0/1.d-8 * ew_lamb
   ew_xobs = !cgs_clight/vth_fiducial/lamb0 * ew_lamb
endif
;
if(keyword_set(xobs)) then begin
   ew_xobs = equivalent_width
   ew_lamb = vth_fiducial*lamb0/!cgs_clight * ew_xobs
   ew_nue = !cgs_clight/lamb0/lamb0/1.d-8 * ew_lamb
endif
;
if(keyword_set(freq)) then begin
   ew_nue = equivalent_width
   ew_lamb = lamb0*lamb0*1.d-8/!cgs_clight * ew_nue
   ew_xobs = !cgs_clight/vht_fiducial/lamb0 * ew_lamb
endif
;
;print, '---------EQUIVALENT WIDTHS----------'
;print, 'ew_lambda = ', ew_lamb
;print, 'ew_nue    = ', ew_nue
;print, 'ew_x      = ', ew_xobs
;
return
;
end
