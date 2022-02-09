PRO get_options, dir, FILE_EXIST, OPT_AIT, OPT_NG, OPT_ALO, OPT_INVERT, SPATIAL_GRID, MU_GRID, PHI_GRID
;
; NAME:
;       getPARAM3D
; PURPOSE:
;       READ IN GRID PARAMETER (DIMENSIONS)
;
;-----------------------DEFINE FILE-NAMES-------------------------------
;
fname=dir+'/options.dat'
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
;------------------------READ IN DATA FROM FILE-------------------------
;
HEADER=''
;
OPENR, 1, fname
   READF, 1, HEADER
      totLen=STRLEN(HEADER)
      OPT_AIT_STR= STRMID(HEADER, totLen-1, 1)
   READF, 1, HEADER
      totLen=STRLEN(HEADER)
      OPT_NG_STR= STRMID(HEADER, totLen-1, 1)
   READF, 1, HEADER
      totLen=STRLEN(HEADER)
      OPT_ALO_STR= STRMID(HEADER, totLen-1, 1)
   READF, 1, HEADER
      totLen=STRLEN(HEADER)
      OPT_INVERT_STR= STRMID(HEADER, totLen-1, 1)
   READF, 1, HEADER
      totLen=STRLEN(HEADER)
      SPATIAL_GRID_STR= STRMID(HEADER, totLen-1, 1)
   READF, 1, HEADER
      totLen=STRLEN(HEADER)
      MU_GRID_STR= STRMID(HEADER, totLen-1, 1)
   READF, 1, HEADER
      totLen=STRLEN(HEADER)
      PHI_GRID_STR= STRMID(HEADER, totLen-1, 1)
CLOSE,1


IF(STRCMP('T', OPT_AIT_STR)) THEN BEGIN
   OPT_AIT=1
ENDIF ELSE BEGIN
   OPT_AIT=0
ENDELSE
;
IF(STRCMP('T', OPT_NG_STR)) THEN BEGIN
   OPT_NG=1
ENDIF ELSE BEGIN
   OPT_NG=0
ENDELSE

OPT_ALO=FIX(OPT_ALO_STR)
OPT_INVERT=FIX(OPT_INVERT_STR)
SPATIAL_GRID=FIX(SPATIAL_GRID_STR)
MU_GRID=FIX(MU_GRID_STR)
PHI_GRID=FIX(PHI_GRID_STR)
;
RETURN
;
END
