PRO getPHI_MASK, dir, FILE_EXIST, DIM_MU, DIM_PHI, PHI_MASK
;
;------------------------DEFINE FILE NAME-------------------------------
;
fname=dir+'/PHI_MASK.dat'
;
;---------------------CHECK IF FILE EXISTS------------------------------
;
FILE_EXIST=FILE_TEST(fname)

IF(FILE_EXIST EQ 0) THEN BEGIN
   PRINT, 'FILE DOES NOT EXIST: ', fname
   RETURN
ENDIF
;
;----------------------READ IN PHI_MASK---------------------------------
;
;NEED TO EXPLICITLY DEFINE INTEGER-ARRAY AS LONG
PHI_MASK=INDGEN(DIM_MU, DIM_PHI, /LONG)
;
OPENR, 1, fname, /F77_UNFORMATTED
   READU, 1, PHI_MASK
CLOSE, 1
;
RETURN
;
END
