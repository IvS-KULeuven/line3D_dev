PRO get_grid_JO, dir, NDJO, RADIUS, VELO, OPAC, OPAL, TEMP, GRADV
;
;-----------------------DEFINE FILE-NAMES-------------------------------
;
fname=dir+'/PHYSICAL_GRID_JO.dat'
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
RADIUS=FLTARR(NDJO)*1.D0
VELO=FLTARR(NDJO)*1.D0
OPAC=FLTARR(NDJO)*1.D0
OPAL=FLTARR(NDJO)*1.D0
TEMP=FLTARR(NDJO)*1.D0
GRADV=FLTARR(NDJO)*1.D0
;
OPENR, 1, fname
   READF, 1, HEADER
   FOR I=0, NDJO-1 DO BEGIN
      READF, 1, VAR0, VAR1, VAR2, VAR3, VAR4, VAR5
      RADIUS(I)=VAR0
      VELO(I)=VAR1
      OPAC(I)=VAR2
      OPAL(I)=VAR3
      TEMP(I)=VAR4
      GRADV(I)=VAR5
   ENDFOR
CLOSE, 1
;
RETURN
;
END
