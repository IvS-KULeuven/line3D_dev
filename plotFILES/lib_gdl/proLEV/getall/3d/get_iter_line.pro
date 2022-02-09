PRO get_iter_line, dir, FILE_EXIST, ITMAXL, DEVMAXL, EPSMAXL_ARR
;
;-----------------------DEFINE FILE-NAMES-------------------------------
;
fname1=dir+'/iterParam_LINE.dat'
fname2=dir+'/iterEPSMAX_LINE.dat'
;
;-----------------------CHECK IF FILES EXIST----------------------------
;
FILE_EXIST=FILE_TEST(fname1)
IF(FILE_EXIST EQ 0) THEN BEGIN
   PRINT, 'FILE DOES NOT EXIST: ', fname1
   RETURN
ENDIF
;
FILE_EXIST=FILE_TEST(fname2)
IF(FILE_EXIST EQ 0) THEN BEGIN
   PRINT, 'FILE DOES NOT EXIST: ', fname2
   RETURN
ENDIF
;
;------------------------READ IN DATA FROM FILES------------------------
;
HEADER=' '
OPENR, 1, fname1
   READF, 1, HEADER
   READF, 1, ITMAXL, DEVMAXL
CLOSE, 1
;
EPSMAXL_ARR=FLTARR(ITMAXL)*1.D0
;
OPENR, 1, fname2, /F77_UNFORMATTED
   READU, 1, EPSMAXL_ARR
CLOSE, 1
;
;
;
RETURN
;
END
