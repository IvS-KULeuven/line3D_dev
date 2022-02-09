PRO get_i3d_angdep, dir, FILE_EXIST, DIM_MU, DIM_PHI, NDXMAX, NDYMAX, NDZMAX, INT3D_ANGDEP
;
; NAME:
;       get_i3d_angdep
; PURPOSE:
;       READ IN INTENSITY AT EACH 3-D GRID POINTS FOR A GIVEN PHI, MU
;
;-----------------------DEFINE FILE-NAMES-------------------------------
;
fname=dir+'/INT3D_ANGDEP.dat'
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
INT3D_ANGDEP=FLTARR(NDXMAX, NDYMAX, NDZMAX, DIM_MU, DIM_PHI)*1.D0
;
OPENR, 1, fname, /F77_UNFORMATTED
   READU, 1, INT3D_ANGDEP
CLOSE, 1
;
RETURN
;
END
