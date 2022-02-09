pro ionization_energy, T, h=h, he=he, o=o, c=c, nue, T, help=print_help
;
;+
; NAME:
;	ionization_energy
;
; PURPOSE:
;	This routine plots the planck-function for a given temperature T, 
;       and overplots estimated ionization energies from ground level
;       according to E_ion=911A*z^2/n^2
;
;
; CALLING SEQUENCE:
;	ionization_energy, T
;
; INPUTS:
;	T:	Temperature in K
;	
; KEYWORD PARAMETERS:
;	help:	Set this keyword (flag) tho show the documentation of this function
;       h:      Set this keyword to overplot hydrogen ionization
;       he:     Set this keyword to overplot he ionization
;       c:      Set this keyword to overplot carbon ionization
;       o:      Set this keyword to overplot oxygen ionization
;       si:     Set this keyword to overplot silicon ionization
;       fe:     Set this keyword to overplot iron ionization
;
; EXAMPLE:
;	ionization_energy, 20.d3, /he
;-
;
;-----------------------------------------------------------------------
;
if(keyword_set(print_help)) then begin
   doc_library, 'bnue'
   return, 0
endif
;
if(n_params() lt 2) then begin
   print, 'Syntax - BNUE(nue, T)'
   return, 0
endif
;
;-----------------------------------------------------------------------
;
;if(h=h) then begin
;z=1
;if(he=he) then begin
;z=2
;if(c=c) then begin
;z=6
;if(o=o) then begin
;z=8
;if(si=si) then begin
;z=14
;if(fe=fe) then begin
;z=26

;
end
