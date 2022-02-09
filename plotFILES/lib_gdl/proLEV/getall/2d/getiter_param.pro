PRO getITER_PARAM, dir, FILE_EXIST, IT_INDX
;
; NAME:
;       getITER_PARAM
; PURPOSE:
;       READ IN PARAMETER FOR ITERATION SCHEME (ITMAX)
;
;-----------------------DEFINE FILE-NAMES-------------------------------
;
fname=dir+'/iteration_IDL.dat'
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
   READF, 1, IT_INDX
CLOSE,1
;
RETURN
;
END
