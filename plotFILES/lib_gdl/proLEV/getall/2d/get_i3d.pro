PRO get_i3d, dir, FILE_EXIST, DIM, NDZMAX, NDXMAX, INT3D
;
; NAME:
;       get_i3d
; PURPOSE:
;       READ IN INTENSITY AT EACH 2-D GRID POINTS FOR A GIVEN XOBS, MU
;
;-----------------------DEFINE FILE-NAMES-------------------------------
;
fname=dir+'/INT3D_IDL.dat'
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
INT3D=FLTARR(NDZMAX, NDXMAX, DIM)*1.D0
;
OPENR, 1, fname, /F77_UNFORMATTED
   READU, 1, INT3D
CLOSE, 1
;
END
