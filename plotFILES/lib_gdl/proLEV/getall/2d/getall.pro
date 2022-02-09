PRO getALL, dir

; NAME:
;       getALL
; PURPOSE:
;       READ IN ALL PHYSICAL VARIABLES AND GRID PARAMETER
;
; CALLING SEQUENCE:
;       getALL, dir
;
; INPUTS:
;        dir: directory from where FORTRAN output files shall be read in

;

   IF (N_ELEMENTS(dir) EQ 0) THEN BEGIN
      PRINT, 'SYNTAX: '
      PRINT, 'getall, dir'
      PRINT, 'with dir=string, directory of model, e.g. "M1" '
      RETURN
   ENDIF

;
END
