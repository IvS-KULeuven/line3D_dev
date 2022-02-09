PRO getPHI_PAIR, dir, FILE_EXIST, DIM_MU, DIM_PHI, NODES_PHI_PAIR
;
;-----------------------DEFINE FILE-NAMES-------------------------------
;
fname=dir+'/NODES_PHI_PAIR.dat'
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
NODES_PHI_PAIR=FLTARR(DIM_MU, DIM_PHI)*1.D0
;
OPENR, 1, fname, /F77_UNFORMATTED
   READU, 1, NODES_PHI_PAIR
CLOSE, 1
;
RETURN
;
END
