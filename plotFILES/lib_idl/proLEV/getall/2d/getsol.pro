PRO getSOL, dir, FILE_EXIST, N1D, XIC1, SR, VMAX, R1D, VEL1D, SLINE1D, SSOBO1D, SCONT1D, SSOBOC1D
;
;-----------------------DEFINE FILE-NAMES-------------------------------
;
fname=dir+'/SOLUTION_1D_IDL.dat'
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
HEADER=''
;
OPENR, 1, fname
   READF, 1, HEADER
   READF, 1, N1D, XIC1, SR, VMAX
;
   R1D=FLTARR(N1D)*1.D0
   VEL1D=FLTARR(N1D)*1.D0
   SLINE1D=FLTARR(N1D)*1.D0
   SSOBO1D=FLTARR(N1D)*1.D0
   SCONT1D=FLTARR(N1D)*1.D0
   SSOBOC1D=FLTARR(N1D)*1.D0
;
   READF, 1, HEADER
   FOR I=0, N1D-1 DO BEGIN
      READF, 1, VAR0, VAR1, VAR2, VAR3, VAR4, VAR5
      R1D(I)=VAR0
      VEL1D(I)=VAR1
      SLINE1D(I)=VAR2
      SSOBO1D(I)=VAR3
      SCONT1D(I)=VAR4
      SSOBOC1D(I)=VAR5
   ENDFOR
CLOSE, 1
;
RETURN
;
END
