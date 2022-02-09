PRO get_raysol, dir, FILE_EXIST, dim_mu, dim_phi, nz_ray, int_ray, int_ray3d, z_ray, opath_ray
;
;-----------------------CHECK IF FILES EXIST----------------------------
;
fname1=dir+'/Z_RAYSOL.dat'
fname2=dir+'/INT1D_ANGDEP.dat'
fname3=dir+'/INT1D_3DANGDEP.dat'
;
;-----------------------CHECK IF FILES EXIST----------------------------
;
FILE_EXIST=FILE_TEST(fname1)
IF(FILE_EXIST EQ 0) THEN BEGIN
   FILE_EXIST=0
   PRINT, 'FILE DOES NOT EXIST: ', fname1
   RETURN
ENDIF
FILE_EXIST=FILE_TEST(fname2)
IF(FILE_EXIST EQ 0) THEN BEGIN
   FILE_EXIST=0
   PRINT, 'FILE DOES NOT EXIST: ', fname2
   RETURN
ENDIF
FILE_EXIST=FILE_TEST(fname3)
IF(FILE_EXIST EQ 0) THEN BEGIN
   FILE_EXIST=0
   PRINT, 'FILE DOES NOT EXIST: ', fname3
   RETURN
ENDIF
;
;------------------------READ IN DATA FROM FILES------------------------
;
OPENR, 1, fname1
   READF, 1, nz_ray

   z_ray=fltarr(nz_ray)*0.d0
   opath_ray=fltarr(nz_ray)*0.d0
;
   FOR I=0, nz_ray-1 DO BEGIN
      READF, 1, VAL0, VAL1
      z_ray(i) = VAL0
      opath_ray(i) = VAL1
   ENDFOR
CLOSE, 1
;
;
;
int_ray=fltarr(nz_ray, dim_mu, dim_phi)*0.d0
int_ray3d=fltarr(nz_ray, dim_mu, dim_phi)*0.d0
;
OPENR, 1, fname2, /F77_UNFORMATTED
   READU, 1, int_ray
CLOSE,1
;
OPENR, 1, fname3, /F77_UNFORMATTED
   READU, 1, int_ray3d
CLOSE,1

;
RETURN
;
END
