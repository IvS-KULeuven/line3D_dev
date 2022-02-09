PRO getFLUXEM_FORSOL, dir, FILE_EXIST, NXOBS, XOBS, FLUX_TOT, FLUX_CONT
;
;-----------------------DEFINE FILE-NAMES-------------------------------
;
fname=dir+'/FLUXEM.dat'
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
VAR0=0.D0
VAR1=0.D0
VAR2=0.D0
VAR3=0.D0
;
OPENR, 1, fname
   READF, 1, NXOBS
   READF, 1, HEADER
   XOBS=FLTARR(NXOBS)*1.D0
   flux_tot=FLTARR(NXOBS)*1.D0
   flux_cont=FLTARR(NXOBS)*1.D0
   FOR I=0, NXOBS-1 DO BEGIN
      READF, 1, VAR0, VAR1, VAR2, VAR3
      XOBS(I)=VAR0
      flux_tot(I)=VAR1
      flux_cont(I)=VAR2
   ENDFOR
CLOSE, 1
;
RETURN
;
END
