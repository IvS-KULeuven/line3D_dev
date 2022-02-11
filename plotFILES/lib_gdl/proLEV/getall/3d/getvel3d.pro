PRO getVEL3D, dir, FILE_EXIST, NDXMAX, NDYMAX, NDZMAX, VEL3D, VELX3D, VELY3D, VELZ3D
;
;-----------------------DEFINE FILE-NAMES-------------------------------
;
fname1=dir+'/MODEL_VEL3D.dat'
fname2=dir+'/MODEL_VELX3D.dat'
fname3=dir+'/MODEL_VELY3D.dat'
fname4=dir+'/MODEL_VELZ3D.dat'
;
;-----------------------CHECK IF FILES EXIST----------------------------
;
FILE_EXIST=FILE_TEST(fname1)
IF(FILE_EXIST EQ 0) THEN BEGIN
   FILE_EXIST=0
   PRINT, 'FILE DOES NOT EXIST: ', fname1
   RETURN
ENDIF
;
FILE_EXIST=FILE_TEST(fname2)
IF(FILE_EXIST EQ 0) THEN BEGIN
   FILE_EXIST=0
   PRINT, 'FILE DOES NOT EXIST: ', fname2
   RETURN
ENDIF
;
FILE_EXIST=FILE_TEST(fname3)
IF(FILE_EXIST EQ 0) THEN BEGIN
   FILE_EXIST=0
   PRINT, 'FILE DOES NOT EXIST: ', fname3
   RETURN
ENDIF
;
FILE_EXIST=FILE_TEST(fname4)
IF(FILE_EXIST EQ 0) THEN BEGIN
   FILE_EXIST=0
   PRINT, 'FILE DOES NOT EXIST: ', fname4
   RETURN
ENDIF
;
;------------------------READ IN DATA FROM FILES------------------------
;
VEL3D=FLTARR(NDXMAX, NDYMAX, NDZMAX)*1.D0
VELX3D=FLTARR(NDXMAX, NDYMAX, NDZMAX)*1.D0
VELY3D=FLTARR(NDXMAX, NDYMAX, NDZMAX)*1.D0
VELZ3D=FLTARR(NDXMAX, NDYMAX, NDZMAX)*1.D0

OPENR, 1, fname1, /F77_UNFORMATTED
   READU, 1, VEL3D
CLOSE, 1
;
OPENR, 1, fname2, /F77_UNFORMATTED
   READU, 1, VELX3D
CLOSE, 1
;
OPENR, 1, fname3, /F77_UNFORMATTED
   READU, 1, VELY3D
CLOSE, 1
;
OPENR, 1, fname4, /F77_UNFORMATTED
   READU, 1, VELZ3D
CLOSE, 1
;
RETURN
;
END