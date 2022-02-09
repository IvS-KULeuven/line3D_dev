PRO getXOBS, dir, FILE_EXIST, NNUE, NUEOBS
;
; NAME:
;       getXOBS
; PURPOSE:
;       READ IN FREQUENCY GRID (NUEOBS)
;
;
;-----------------------DEFINE FILE-NAMES-------------------------------
;
fname=dir+'/NUEOBS_IDL.dat'
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
NUEOBS=FLTARR(NNUE)*1.D0
;
OPENR, 1, fname
   READF, 1, HEADER
   FOR I=0,NNUE-1 DO BEGIN
      READF, 1, VAR0
      NUEOBS(I)=VAR0
   ENDFOR
CLOSE,1
;
RETURN
;
END
