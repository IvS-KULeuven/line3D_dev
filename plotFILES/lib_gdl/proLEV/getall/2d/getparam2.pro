PRO getPARAM2, dir, FILE_EXIST, N1D, NCZ, NCX, NDZMAX, NDXMAX, DIM
;
; NAME:
;       getPARAM
; PURPOSE:
;       READ IN GRID PARAMETER (DIMENSIONS)
;
;-----------------------DEFINE FILE-NAMES-------------------------------
;
fname=dir+'/gridparamIDL.dat'
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
HEADER=''
;
OPENR, 1, fname
   READF, 1, HEADER
   READF, 1, N1D, NDZ, NCZ, NDX, NCX, NDZMAX, NDXMAX, DIM
CLOSE,1
;
RETURN
;
END
