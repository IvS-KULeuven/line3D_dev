PRO getMU, dir, FILE_EXIST, DIM, ANGLENODES
;
; NAME:
;       getMU
; PURPOSE:
;       READ IN ANGULAR GRID (MU)
;
;-----------------------DEFINE FILE-NAMES-------------------------------
;
fname=dir+'/ANGLES_IDL.dat'
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
ANGLENODES=FLTARR(DIM)*1.D0
;
OPENR, 1, fname
   READF, 1, HEADER
   FOR I=0, DIM-1 DO BEGIN
      READF, 1, VAR0
      ANGLENODES(I)=VAR0
   ENDFOR
CLOSE, 1
;
RETURN
;
END
