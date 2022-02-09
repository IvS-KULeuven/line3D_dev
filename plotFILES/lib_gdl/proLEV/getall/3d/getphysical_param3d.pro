PRO getPHYSICAL_PARAM3D, dir, FILE_EXIST, VTH, DELDOP, XKLINE, FAC, XNUE0, EPS_LINE, EPS_CONT, VMAX, VMIN, BETA, SR, XIC1, TEFF, $
                         rhoc_star=RHOC_STAR, rhow_star=RHOW_STAR, t_inf=T_INF, v_esc=V_ESC, obl=OBL, new=NEW, nnew=NNEW
;
; NAME:
;       getPHYSICAL_PARAM
; PURPOSE:
;       READ IN PHYSICAL PARAMETER (THERMAL VELOCITY, ETC.)
;
;-----------------------DEFINE FILE-NAMES-------------------------------
;
fname=dir+'/physicalPARAM_IDL.dat'
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
;
IF(KEYWORD_SET(NEW)) THEN BEGIN
   OPENR, 1, fname
      READF, 1, HEADER
      READF, 1, VTH, DELDOP, XKLINE, FAC, XNUE0, EPS_LINE, EPS_CONT, VMAX, VMIN, BETA, SR, XIC1, TEFF
   CLOSE,1
ENDIF ELSE BEGIN
   IF(KEYWORD_SET(NNEW)) THEN BEGIN
      OPENR, 1, fname
         READF, 1, HEADER
         READF, 1, VTH, DELDOP, XKLINE, FAC, XNUE0, EPS_LINE, EPS_CONT, VMAX, VMIN, BETA, SR, XIC1, TEFF, $
                   RHOC_STAR, RHOW_STAR, T_INF, V_ESC, OBL
      CLOSE, 1
   ENDIF ELSE BEGIN
      OPENR, 1, fname
         READF, 1, HEADER
         READF, 1, VTH, DELDOP, XKLINE, XNUE0, EPS_LINE, VMAX, VMIN, BETA, SR, XIC1, TEFF
      CLOSE,1
   ENDELSE
ENDELSE
;
;
;
RETURN
;
END
