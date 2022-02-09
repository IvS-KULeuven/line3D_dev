PRO getOPAL3D, dir, FILE_EXIST, NDXMAX, NDYMAX, NDZMAX, OPALBAR3D
;
; NAME:
;       getOPAL3D
; PURPOSE:
;       READ IN 3D- FREQUENCY INTEGRATED LINE OPACITY
;
;-----------------------DEFINE FILE-NAMES-------------------------------
;
fname=dir+'/MODEL_OPALBAR3D.dat'
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
OPALBAR3D=FLTARR(NDXMAX,NDYMAX,NDZMAX)*1.D0
;
OPENR, 1, fname, /F77_UNFORMATTED
   READU, 1, OPALBAR3D
CLOSE, 1
;
RETURN
;
END
