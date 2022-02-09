PRO getITER_HISTO, dir, FILE_EXIST, N1D, ITMAX, epshistory
;
; NAME:
;       getITER_HISTO
; PURPOSE:
;       READ IN DEVIATION HISTORY OF ITERATION SCHEMEy
;
;-----------------------DEFINE FILE-NAMES-------------------------------
;
fname=dir+'/ithistory_IDL.dat'
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
epshistory=FLTARR(N1D, ITMAX)*1.D0
;
OPENR, 3, fname, /F77_UNFORMATTED
   READU, 3, epshistory
CLOSE, 3
;
RETURN
;
END
