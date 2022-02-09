PRO getPHI, dir, FILE_EXIST, DIM_PHI, NODES_PHI
;
;-----------------------DEFINE FILE-NAMES-------------------------------
;
fname=dir+'/GRID_PHI.dat'
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
NODES_PHI=FLTARR(DIM_PHI)*1.D0
;
OPENR, 1, fname, /F77_UNFORMATTED
   READU, 1, NODES_PHI
CLOSE, 1
;
RETURN
;
END
