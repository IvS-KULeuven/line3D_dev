PRO getMINTBAR3D, dir, FILE_EXIST, NDXMAX, NDYMAX, NDZMAX, MINTBAR3D
;
;-----------------------DEFINE FILE-NAMES-------------------------------
;
fname=dir+'/SOLUTION_MINTBAR3D.dat'
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
MINTBAR3D=FLTARR(NDXMAX, NDYMAX, NDZMAX)*1.D0
;
OPENR, 1, fname, /F77_UNFORMATTED
   READU, 1, MINTBAR3D
CLOSE, 1
;
RETURN
;
END
