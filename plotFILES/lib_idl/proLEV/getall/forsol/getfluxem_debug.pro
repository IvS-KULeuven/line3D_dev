PRO getFLUXEM_DEBUG, dir, name_indx, FILE_EXIST, NXOBS, ALPHA, GAMMA, XOBS, FLUX_TOT, FLUX_CONT, FLUX_EMI, FLUX_ABS, file_name=file_name, help=print_help
;+
; NAME:
;       getFLUXEM_DEBUG
;
; PURPOSE:
;       This procedure reads in emergent flux profiles from output-files
;
; CALLING SEQUENCE:
;       getFLUXEM_DEBUG, dir, name_indx, file_exist, nxobs, alpha, gamma, $
;       xobs, flux_tot, flux_cont, flux_emi, flux_abs
;
; INPUTS:
;       dir:    String defining the directory where files are stored.
;       name_indx: Integer defining the file-name, which shall be read in
;
; KEYWORD PARAMETERS:
;       help:   Set this keyword (flag) to show the documentation of this
;       procedure
;       file_name: Set this keyword to a string, defining an explicit name of
;                  data-file (e.g. if the file is stored in unusual format)
;      
; OUTPUTS:
;       nxobs:   Integer defining the length of the arrays.
;       file_exist: Logical, which returns true if file exists.
;       alpha:   Floating point variable, returns inclination angle.
;       gamma:   FLoating point variable, returns rotational phase angle.
;       xobs:    Array of length nxobs. Frequency shift from line-center in
;       units of vth_fiducial (or wavelength in A, if conv-keyword is set)
;       flux_tot:  Array of length nxobs. Total flux in cgs
;       flux_cont: Array of length nxobs. Flux of neighbouring continuum in
;       cgs
;       flux_emi:  Array of length nxobs. Emission part of flux in cgs
;       flux_abs:  Array of length nxobs. Absorption part of flux in cgs
;
; EXAMPLE:
;       getFLUXEM_DEBUG, '.', 0, f_exist, 100, alpha, gamma, xobs, ftot,
;       fcont, femi, fabs
;       getFLUXEM_DEBUG, '.', 10, f_exist, 100, alpha, gamma, lamb, ftot,
;       fcont, femi, fabs, /CONV
;       getFLUXEM_DEBUG, '.', 10, f_exist, 100, alpha, gamma, xobs, ftot,
;       fcont, femi, fabs, file_name='unusual_output_name.dat'
;- 
;-------------------------OUTPUT IF HELP NEEDED-------------------------

if(keyword_set(print_help)) then begin 
 doc_library, 'getfluxem_debug'
 return
endif
;
;-----------------------DEFINE FILE-NAMES-------------------------------
;
IF(KEYWORD_SET (file_name)) THEN BEGIN
   fname=dir+'/'+file_name
ENDIF ELSE BEGIN
   dummy_name=STRING(name_indx, FORMAT='(I5.5)')
   fname=dir+'/FLUXEM_DEBUG_'+dummy_name+'.dat'
ENDELSE
;
;-----------------------CHECK IF FILES EXIST----------------------------
;
FILE_EXIST=FILE_TEST(fname)
IF(FILE_EXIST EQ 0) THEN BEGIN
   FILE_EXIST=0
   PRINT, 'FILE DOES NOT EXIST: ', fname
   RETURN
ENDIF
;
;------------------------READ IN DATA FROM FILES------------------------
;
HEADER=''
VAR0=0.D0
VAR1=0.D0
VAR2=0.D0
VAR3=0.D0
VAR4=0.D0
;
OPENR, 1, fname
   READF, 1, NXOBS
   READF, 1, FORMAT='(A6, E20.8)', HEADER, ALPHA
   READF, 1, FORMAT='(A6, E20.8)', HEADER, GAMMA
   READF, 1, HEADER
   XOBS=FLTARR(NXOBS)*1.D0
   flux_tot=FLTARR(NXOBS)*1.D0
   flux_cont=FLTARR(NXOBS)*1.D0
   flux_emi=FLTARR(NXOBS)*1.D0
   flux_abs=FLTARR(NXOBS)*1.D0
   FOR I=0, NXOBS-1 DO BEGIN
      READF, 1, VAR0, VAR1, VAR2, VAR3, VAR4
      XOBS(I)=VAR0
      flux_tot(I)=VAR1
      flux_cont(I)=VAR2
      flux_emi(I)=VAR3
      flux_abs(I)=VAR4
   ENDFOR
CLOSE, 1
;
RETURN
;
END
