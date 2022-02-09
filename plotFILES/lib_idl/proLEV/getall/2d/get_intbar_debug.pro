PRO get_intbar_debug, dir, FILE_EXIST, N1D, DIM, INTBAR
;
; NAME:
;       get_intbar_debug
; PURPOSE:
;       READ IN FREQUENCY INTEGRATED INTENSITY FOR EACH MU AND GRID-POINT
;
;-----------------------DEFINE FILE-NAMES-------------------------------
;
fname=dir+'/INTBAR_DEBUG.dat'
;
;-----------------------CHECK IF FILES EXIST----------------------------
;
FILE_EXIST=FILE_TEST(fname)
IF(FILE_EXIST EQ 0) THEN BEGIN
   PRINT, 'FILE DOES NOT EXIST: ', fname
   RETURN
ENDIF
;
;------------------------READ IN DATA FROM FILES------------------------
;
INTBAR=FLTARR(N1D, DIM)*1.D0
;
OPENR, 1, fname, /F77_UNFORMATTED
   READU, 1, INTBAR
CLOSE, 1
;
END
