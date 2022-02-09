PRO getPARAM_JO, dir, FILE_EXIST, NDJO, XIC1=XIC1 
;                         
; NAME:                           
;       getPARAM_JO
; PURPOSE:
;       READ IN GRID PARAMETER (DIMENSIONS) FROM JO'S PROGRAM
;
;-----------------------DEFINE FILE-NAMES-------------------------------
;
fname=dir+'/gridPARAM_JO.dat'
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
XIC1=0.D0
; 
OPENR, 1, fname 
   READF, 1, HEADER
   READF, 1, NDJO, XIC1
CLOSE,1
;
RETURN
;
END
