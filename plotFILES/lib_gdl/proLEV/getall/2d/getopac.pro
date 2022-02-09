PRO getOPAC, dir, FILE_EXIST, N1D, OPAC1D, OPAC2D, NDZMAX, NDXMAX, NCZ
;
;-----------------------DEFINE FILE-NAMES-------------------------------
;
fname=dir+'/OPATH2D_IDL.dat'
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
OPAC2D=FLTARR(NDZMAX, NDXMAX)*1.D0
;
OPENR, 1, fname, /F77_UNFORMATTED
   READU, 1, OPAC2D
CLOSE, 1
;
OPAC1D=FLTARR(N1D)*1.D0
OPAC1D=OPAC2D(NDZMAX/2+NCZ : NDZMAX-1, NDXMAX/2)
;
RETURN
;
END
