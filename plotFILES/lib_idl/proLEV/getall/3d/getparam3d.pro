PRO getPARAM3D, dir, FILE_EXIST, N1D, NDX, NCX, NDXMAX, NDY, NCY, NDYMAX, NDZ, NCZ, NDZMAX, DIM_MU, DIM_PHI, NXOBS
;
; NAME:
;       getPARAM3D
; PURPOSE:
;       READ IN GRID PARAMETER (DIMENSIONS)
;
;-----------------------DEFINE FILE-NAMES-------------------------------
;
fname=dir+'/gridparam3D_IDL.dat'
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
   READF, 1, N1D, NDX, NCX, NDXMAX, NDY, NCY, NDYMAX, NDZ, NCZ, NDZMAX, DIM_MU, DIM_PHI, NXOBS
CLOSE,1
;
RETURN
;
END
