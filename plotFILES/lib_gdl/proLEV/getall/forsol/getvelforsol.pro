PRO getvelforsol, dir, FILE_EXIST, np, nzeta, nz, velx_ray, vely_ray, velz_ray
;
;-----------------------DEFINE FILE-NAMES-------------------------------
;
fname1=dir+'/forsol_velx.dat'
fname2=dir+'/forsol_vely.dat'
fname3=dir+'/forsol_velz.dat'
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
;------------------------READ IN DATA FROM FILES------------------------
;
velx_ray=FLTARR(np, nzeta, nz)*1.D0
vely_ray=FLTARR(np, nzeta, nz)*1.D0
velz_ray=FLTARR(np, nzeta, nz)*1.D0

OPENR, 1, fname1, /F77_UNFORMATTED
   READU, 1, velx_ray
CLOSE, 1
;
OPENR, 1, fname2, /F77_UNFORMATTED
   READU, 1, vely_ray
CLOSE, 1
;
OPENR, 1, fname3, /F77_UNFORMATTED
   READU, 1, velz_ray
CLOSE, 1
;
;
RETURN
;
END
