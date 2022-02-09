PRO getZETA, dir, FILE_EXIST, nzeta, zeta
;
;-----------------------DEFINE FILE-NAMES-------------------------------
;
fname=dir+'/forsol_zetagrid.dat'
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
zeta=FLTARR(nzeta)*1.D0
;
OPENR, 1, fname, /F77_UNFORMATTED
   READU, 1, zeta
CLOSE, 1
;
RETURN
;
END
