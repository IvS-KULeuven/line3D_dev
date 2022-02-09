PRO get_ithistory, dir, FILE_EXIST, ITMAX, DEVMAX, EPSMAX_ARR
;
;-----------------------DEFINE FILE-NAMES-------------------------------
;
fname1=dir+'/gridparamIDL.dat'
fname2=dir+'/iterParam.dat'
fname3=dir+'/ithistory_IDL.dat'
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
FILE_EXIST=FILE_TEST(fname3)
IF(FILE_EXIST EQ 0) THEN BEGIN
   PRINT, 'FILE DOES NOT EXIST: ', fname3
   RETURN
ENDIF
;
;------------------------READ IN DATA FROM FILES------------------------
;
HEADER=''
OPENR, 1, fname1
   READF, 1, HEADER
   READF, 1, N1D, NDZ, NCZ, NDX, NCX, NDZMAX, NDXMAX, DIM
CLOSE,1
;
OPENR, 2, fname2
   READF, 2, HEADER
   READF, 2, ITMAX, DEVMAX
CLOSE,2
;
epshistory=FLTARR(N1D, ITMAX)*1.D0
;
;
OPENR, 3, fname3, /F77_UNFORMATTED
   READU, 3, epshistory
CLOSE, 3
;
EPSMAX_ARR=FLTARR(ITMAX)*1.D0
DUMMY=FLTARR(N1D)*1.D0
;
FOR I=0, ITMAX-1 DO BEGIN
   DUMMY=ABS(epshistory(*,I))
   EPSMAX_ARR(I)=MAX(DUMMY)
ENDFOR

;
RETURN
;
END
