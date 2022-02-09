PRO getGRID, dir, FILE_EXIST, NDZMAX, NDXMAX, X, Z
;
; NAME:
;       getGRID
; PURPOSE:
;       READ IN SPATIAL GRID (X,Z-AXIS)
;
;-----------------------DEFINE FILE-NAMES-------------------------------
;
fnamex=dir+'/XAXIS_IDL.dat'
fnamez=dir+'/ZAXIS_IDL.dat'
;
;-----------------------CHECK IF FILES EXIST----------------------------
;
FILE_EXIST=FILE_TEST(fnamex)
IF(FILE_EXIST EQ 0) THEN BEGIN
   PRINT, 'FILE DOES NOT EXIST: ', fnamex
   RETURN
ENDIF
;
FILE_EXIST=FILE_TEST(fnamez)
IF(FILE_EXIST EQ 0) THEN BEGIN
   PRINT, 'FILE DOES NOT EXIST: ', fnamez
   RETURN
ENDIF
;
;------------------------READ IN DATA FROM FILES------------------------
;
HEADER=''
;
X=FLTARR((NDXMAX))*1.D0
Z=FLTARR((NDZMAX))*1.D0
;
;READ IN GRID
VAR0=1.D0
;
OPENR, 1, fnamex
   READF, 1, HEADER
   FOR I=0, NDXMAX-1 DO BEGIN
      READF, 1, VAR0
      X(I)=VAR0
   ENDFOR
CLOSE, 1
;
OPENR, 1, fnamez
   READF, 1, HEADER
   FOR I=0, NDZMAX-1 DO BEGIN
      READF, 1, VAR0
      Z(I)=VAR0
   ENDFOR
CLOSE, 1

END
