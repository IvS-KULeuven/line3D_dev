PRO getGRID1D, dir, FILE_EXIST, N1D, R1D, VEL1D, OPATH1D, OPALBAR1D, T1D
;
;-----------------------DEFINE FILE-NAMES-------------------------------
;
fname=dir+'/grids1D.dat'
;
;-----------------------CHECK IF FILES EXIST----------------------------
;
FILE_EXIST=FILE_TEST(fname )
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
R1D=FLTARR(N1D)*1.D0
VEL1D=FLTARR(N1D)*1.D0
OPATH1D=FLTARR(N1D)*1.D0
OPALBAR1D=FLTARR(N1D)*1.D0
T1D=FLTARR(N1D)*1.D0
;
;
;
OPENR, 1, fname
   READF, 1, HEADER
   FOR I=0, N1D-1 DO BEGIN
      READF, 1, VAR0, VAR1, VAR2, VAR3, VAR4
      R1D(I)=VAR0
      VEL1D(I)=VAR1
      OPATH1D(I)=VAR2
      OPALBAR1D(I)=VAR3
      T1D(I)=VAR4
   ENDFOR
CLOSE, 1
;
RETURN
;
END
