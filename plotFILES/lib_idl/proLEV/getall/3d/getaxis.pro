PRO getAXIS, dir, FILE_EXIST, NDXMAX, NDYMAX, NDZMAX, X, Y, Z
;
;-----------------------DEFINE FILE-NAMES-------------------------------
;
fnamex=dir+'/X_AXIS.dat'
fnamey=dir+'/Y_AXIS.dat'
fnamez=dir+'/Z_AXIS.dat'
;
;-----------------------CHECK IF FILES EXIST----------------------------
;
FILE_EXIST=FILE_TEST(fnamex)
IF(FILE_EXIST EQ 0) THEN BEGIN
   FILE_EXIST=0
   PRINT, 'FILE DOES NOT EXIST: ', fnamex
   RETURN
ENDIF
;
FILE_EXIST=FILE_TEST(fnamey)
IF(FILE_EXIST EQ 0) THEN BEGIN
   FILE_EXIST=0
   PRINT, 'FILE DOES NOT EXIST: ', fnamey
   RETURN
ENDIF
;
FILE_EXIST=FILE_TEST(fnamez)
IF(FILE_EXIST EQ 0) THEN BEGIN
   FILE_EXIST=0
   PRINT, 'FILE DOES NOT EXIST: ', fnamez
   RETURN
ENDIF
;
;------------------------READ IN AXES FROM FILES------------------------
;
X=FLTARR(NDXMAX)*1.D0
Y=FLTARR(NDYMAX)*1.D0
Z=FLTARR(NDZMAX)*1.D0
;
OPENR, 1, fnamex, /F77_UNFORMATTED
   READU, 1, X
CLOSE, 1
;
OPENR, 1, fnamey, /F77_UNFORMATTED
   READU, 1, Y
CLOSE, 1
;
OPENR, 1, fnamez, /F77_UNFORMATTED
   READU, 1, Z
CLOSE, 1
;
RETURN
;
END
