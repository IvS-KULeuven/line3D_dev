PRO printPARAMS, dir

;+
; NAME:
;	printparams
; PURPOSE:
;	READ IN DEVIATION HISTORY OF ITERATION PROCEDURE
;
; CALLING SEQUENCE:
;	printparams, dir
;
; INPUTS:

   IF (N_ELEMENTS(dir) EQ 0) THEN BEGIN
      PRINT, 'SYNTAX: '
      PRINT, 'printparams, dir'
      PRINT, 'with dir=string , e.g. printparams, "M1" '
      RETURN
   ENDIF
;
;--------------------------GET GRID PARAMETER---------------------------
;
getPARAM3D, dir, FILE_EXIST1, N1D, NDX, NCX, NDXMAX, NDY, NCY, NDYMAX, NDZ, NCZ, NDZMAX, DIM_MU, DIM_PHI, NXOBS
;
;-----------------------GET ITERATION PARAMETER-------------------------
;
get_iter, dir, FILE_EXIST2, ITMAX, DEVMAX, EPSMAX
;
;-----------------------GET MODEL-ATMOSPHERE PARAMETER------------------
;
getPHYSICAL_PARAM, dir, FILE_EXIST3, VTH, DELDOP, XKLINE, XNUE0, EPSPARAM, VMAX, VMIN, BETA, SR, XIC1, TEFF
;
;--------------------------GET GLOBAL OPTIONS---------------------------
;
get_options, dir, FILE_EXIST4, OPT_AIT, OPT_NG, OPT_ALO, OPT_INVERT, SPATIAL_GRID, MU_GRID, PHI_GRID
;
;--------------------------GET GRID IN 1D-------------------------------
;
getGRID1D, dir, FILE_EXIST5, N1D, R1D, VEL1D, OPATH1D, OPALBAR1D, T1D
;
;--------------------------GET TIMING PROPERTIES------------------------
;
get_timing, dir, FILE_EXIST6, TTOT_SL, TTOT_ALO, TTOT_IT, T_TOT, DIM_MU, DIM_PHI, IT_TOT
;
;--------------------CHECK IF ALL FILES WHERE FOUND---------------------
;
IF(FILE_EXIST1*FILE_EXIST2*FILE_EXIST3*FILE_EXIST4*FILE_EXIST5*FILE_EXIST6 EQ 0) THEN BEGIN
   PRINT, 'ERROR IN OPENING FILES: RETURNING'
   RETURN
ENDIF
;
;-----------------------PRINT OUT ALL PARAMETER-------------------------
;
PRINT, '--------------------MODEL ATMOSPHERE INPUT-----------------'
PRINT, FORMAT='(A20, E20.8)', 'T_EFF [K]', TEFF
PRINT, FORMAT='(A20, A20)', 'MDOT ', 'SEE IN INDAT'
PRINT, FORMAT='(A20, E20.8)', 'VMIN [KM/S]', VMIN
PRINT, FORMAT='(A20, E20.8)', 'VMAX [KM/S]', VMAX
PRINT, FORMAT='(A20, E20.8)', 'VTH_FIDUCIAL [KM/S]', VTH/1.D5
PRINT, FORMAT='(A20, E20.8)', 'BETA', BETA
PRINT, FORMAT='(A20, E20.8)', 'R_MIN [R_SUN]', R1D(0)/(696.d8)
PRINT, FORMAT='(A20, E20.8)', 'R_MAX [R_STAR]', R1D(N1D-1)/SR
PRINT, ''
;
PRINT, '-----------------------OPTIONS-----------------------------'
PROPTIONS, OPT_AIT, OPT_NG, OPT_ALO, OPT_INVERT, SPATIAL_GRID, MU_GRID, PHI_GRID
PRINT, ''
;
PRINT, '-------------------ITERATION PARAMETER---------------------'
PRINT, FORMAT='(A30, E20.8)', 'Max number of iterations', ITMAX
PRINT, FORMAT='(A30, E20.8)', 'Convergence criterion DEVMAX', DEVMAX
PRINT, ''
;
PRINT, '-----------------SPATIAL GRID PARAMETER--------------------'
PRINT, FORMAT='(A32, I8)', 'No. radial grid points N1D     ', N1D
PRINT, FORMAT='(A32, I8)', 'No. grid points inside core NCX', NCX
PRINT, FORMAT='(A32, I8)', 'No. grid points inside core NCY', NCY
PRINT, FORMAT='(A32, I8)', 'No. grid points inside core NCZ', NCZ
PRINT, FORMAT='(A32, I8)', 'Total no. of grid points NDXMAX', NDXMAX
PRINT, FORMAT='(A32, I8)', 'Total no. of grid points NDYMAX', NDYMAX
PRINT, FORMAT='(A32, I8)', 'Total no. of grid points NDZMAX', NDZMAX
PRINT, ''
;
PRINT, '-----------------ANGULAR GRID PARAMETER--------------------'
PRINT, FORMAT='(A20, I8)', 'No. mu  grid points', DIM_MU
PRINT, FORMAT='(A20, I8)', 'No. phi grid points', DIM_PHI
PRINT, ''
;
PRINT, '-----------------OTHER PARAMETER---------------------------'
PRINT, FORMAT='(A20, E20.8)', 'Frequency point', XNUE0
PRINT, FORMAT='(A20, A20)', 'FAC         =  ', 'See in inMODEL-File'
PRINT, FORMAT='(A20, A20)', 'EPS_THOMSON =  ', 'See in inMODEL-File'
;
PRINT, '-----------------------TIMING------------------------------'
PRTIMING, TTOT_SL, TTOT_ALO, TTOT_IT, T_TOT, DIM_MU, DIM_PHI, IT_TOT
PRINT, ''
;
END
