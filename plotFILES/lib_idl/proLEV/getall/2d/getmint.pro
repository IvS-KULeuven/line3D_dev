PRO getMINT, dir, FILE_EXIST, NDZMAX, NDXMAX, J2D
;
;-----------------------DEFINE FILE-NAMES-------------------------------
;
fname=dir+'/J2D_IDL.dat'
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
J2D=FLTARR(NDZMAX, NDXMAX)*1.D0
;
OPENR, 1, fname, /F77_UNFORMATTED
   READU, 1, J2D
CLOSE, 1
;
RETURN
;
END
