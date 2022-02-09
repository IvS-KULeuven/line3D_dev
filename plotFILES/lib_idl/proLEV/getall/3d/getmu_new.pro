PRO getMU_NEW, dir, FILE_EXIST, DIM_MU, NODES_MU
;
;-----------------------DEFINE FILE-NAMES-------------------------------
;
fname=dir+'/GRID_MU.dat'
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
NODES_MU=FLTARR(DIM_MU)*1.D0
;
OPENR, 1, fname, /F77_UNFORMATTED
   READU, 1, NODES_MU
CLOSE, 1
;
RETURN
;
END
