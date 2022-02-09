PRO getPARAM_JO_OLD, dir, FILE_EXIST, NDJO   
;                         
; NAME:                           
;       getPARAM_JO
; PURPOSE:
;       READ IN GRID PARAMETER (DIMENSIONS) FROM JO'S PROGRAM
;
;-----------------------DEFINE FILE-NAMES-------------------------------
;
fname=dir+'/param_JO.dat'
;
;-----------------------CHECK IF FILES EXIST----------------------------
;
FILE_EXIST=FILE_TEST(fname)
IF(FILE_EXIST EQ 0) THEN BEGIN
   PRINT, 'FILE DOES NOT EXIST: ', fname
   RETURN
ENDIF
;
NDJO=1
;
;------------------------READ IN DATA FROM FILES------------------------
;
HEADER=''
; 
OPENR, 1, fname 
   READF, 1, HEADER
   READF, 1, NDJO 
CLOSE,1
;
RETURN
;
END
