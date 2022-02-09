PRO getPHYSICAL_PARAM, dir, FILE_EXIST, VTH, DELDOP, XKLINE, XNUE0, EPSPARAM, VMAX, VMIN, BETA, SR, XIC1, TEFF
;
; NAME:
;       getPHYSICAL_PARAM
; PURPOSE:
;       READ IN PHYSICAL PARAMETER (THERMAL VELOCITY, ETC.)
;
;
;-----------------------DEFINE FILE-NAMES-------------------------------
;
fname=dir+'/physicalPARAM_IDL.dat'
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
   READF, 1, VTH, DELDOP, XKLINE, XNUE0, EPSPARAM, VMAX, VMIN, BETA, SR, XIC1, TEFF
CLOSE,1
;
END
