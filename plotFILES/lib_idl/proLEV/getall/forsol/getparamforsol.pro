PRO getPARAMFORSOL, dir, FILE_EXIST, np, nzeta, nz
;
; NAME:
;       getPARAM3D
; PURPOSE:
;       READ IN GRID PARAMETER (DIMENSIONS)
;
;-----------------------DEFINE FILE-NAMES-------------------------------
;
fname=dir+'/forsol_dime.dat'
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
OPENR, 1, fname
   READF, 1, HEADER
   READF, 1, np, nzeta, nz
CLOSE,1
;
RETURN
;
END
