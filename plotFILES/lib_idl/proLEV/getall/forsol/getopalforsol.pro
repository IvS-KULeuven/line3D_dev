PRO getOPALFORSOL, dir, FILE_EXIST, np, nzeta, nz, opalbar_ray
;
; NAME:
;       getOPAL3D
; PURPOSE:
;       READ IN 3D- FREQUENCY INTEGRATED LINE OPACITY
;
;-----------------------DEFINE FILE-NAMES-------------------------------
;
fname=dir+'/forsol_opalbar.dat'
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
OPALBAR_RAY=FLTARR(NP, NZETA, NZ)*1.D0
;
OPENR, 1, fname, /F77_UNFORMATTED
   READU, 1, OPALBAR_RAY
CLOSE, 1
;
RETURN
;
END
