PRO getXOBS_3D, dir, FILE_EXIST, NXOBS, NODES_XOBS
;
;-----------------------DEFINE FILE-NAMES-------------------------------
;
fname=dir+'/GRID_XOBS.dat'
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
NODES_XOBS=FLTARR(NXOBS)*1.D0
;
OPENR, 1, fname, /F77_UNFORMATTED
   READU, 1, NODES_XOBS
CLOSE, 1
;
RETURN
;
END
