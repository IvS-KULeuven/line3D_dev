PRO getIMPACT, dir, FILE_EXIST, np, p
;
;-----------------------DEFINE FILE-NAMES-------------------------------
;
fname=dir+'/forsol_pgrid.dat'
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
p=FLTARR(np)*1.D0
;
OPENR, 1, fname, /F77_UNFORMATTED
   READU, 1, p
CLOSE, 1
;
RETURN
;
END
