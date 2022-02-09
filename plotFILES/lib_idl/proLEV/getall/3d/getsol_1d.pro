PRO getSOL_1D, dir, FILE_EXIST, N1D, MINT1D, SCONT1D
;
;-----------------------CHECK IF FILES EXIST----------------------------
;
fname=dir+'/solution1D.dat'
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
HEADER=''
;
MINT1D=FLTARR(N1D)*1.D0
SCONT1D=FLTARR(N1D)*1.D0
;
OPENR, 1, fname
   READF, 1, HEADER
;
   FOR I=0, N1D-1 DO BEGIN
      READF, 1, VAR0, VAR1, VAR2
      MINT1D(I)=VAR1
      SCONT1D(I)=VAR2
   ENDFOR
CLOSE, 1
;
RETURN
;
END
