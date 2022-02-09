PRO get_i3d_debug, dir, FILE_EXIST, N1D, DIM, NNUE, INT3D
;
; NAME:
;       get_i3d_debug
; PURPOSE:
;       READ IN INTENSITY AT 1-D GRID POINTS FOR EACH ANGLE AND FREQUENCY
;
;-----------------------DEFINE FILE-NAMES-------------------------------
;
fname=dir+'/INT3D_DEBUG.dat'
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
INT3D=FLTARR(N1D, DIM, NNUE)*1.D0
;
OPENR, 1, fname, /F77_UNFORMATTED
   READU, 1, INT3D
CLOSE, 1
;
RETURN
;
END
