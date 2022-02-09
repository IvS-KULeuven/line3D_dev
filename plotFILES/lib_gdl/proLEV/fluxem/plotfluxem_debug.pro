PRO plotFLUXEM_DEBUG, dir, fname=FNAME, obs_indx=OBS_INDX, oname=ONAME, windx=WINDX, ylim=YLIM, xlim=XLIM, oplot=OPLOT, cindx=CINDX, freq=FREQ, wave=WAVE, vel=VEL, lam0=LAM0, vth_fid=VTH_FID, xscalefac=xscalefac, isotropic=ISOTROPIC, charsize=charsize, help=print_help
;+
; NAME:
;	plotFLUXEM_DEBUG
;
; PURPOSE:
;	This procedure plots the emitted flux of a spectral line:
;
;          profile_tot = F_tot/F_cont
;          profile_emi = F_emi/F_cont   only emission profile
;          profile_abs = F_abs/F_cont   only absorption profile
;
;       where F_cont is the flux in the absence of the line
;
; CALLING SEQUENCE:
;
;	plotfluxem_debug, dir
;
; INPUTS:
;	dir:	Directory, where appropriate files are stored.
;	
; KEYWORD PARAMETERS:
;       fname:  Set this keyword to the name of a file inside directory dir, 
;               from which the input fluxes shall be read in
;
;       obs_indx: Set this keyword to an integer, to describe from which file
;                 inside directory dir, the input fluxes are read in:
;                 e.g. obs_indx=1 reads from 'dir/FLUXEM_00001.dat'
;
;       oname:  Set this keyword to the output-name of a ps-file, if output
;               shall be plotted to that file
;
;       windx:  Set this keyword to an integer defining the window, to which
;               output is plotted
;
;       ylim:   Set this keyword to a 2-d array, which sets the yrange
;
;       xlim:   Set this keyword to a 2-d array, which sets the xrange
;
;       oplot:  Set this keyword (flag), to overplot in the current window
;               This keyword is not allowed, when oname is set
;
;       cindx:  Set this keyword to an integer, describing the color in which
;               the emergent profile is plotted (affects only overplotted
;               profiles)
;
;       freq:   Set this keyword (flag) to plot profile in frequency-regime
;
;       wave:   Set this keyword (flag) to plot profile in wavelength-regime
;
;       vel:    Set this keyword (flag) to plot profile in velocity-regime
;
;       lam0:   Set this keyword to transition-wavelength (in Angstrom)
;
;       vth_fid: Set this keyword to fiducial thermal velocity (in km/s)
;
;       isotripic: Set this keyword to make an isotropic plot
;
;       xscalefac: Set this keyword to a floating point number, which
;       describes an up/downscaling of the x-axis (e.g. v_inf1/v_inf2)
;
;       help:   Set this keyword (flag) to show the documentation of the
;               procedure
;
; EXAMPLE:
;       plotfluxem, '.', obs_indx=1, ylim=[0.5,-0.5], xlim=[0.,2.5]
;       plotfluxem, '.', fname='FLUXEM_00002.dat', /oplot, cindx=200
;       plotfluxem, '.', obs_indx=1, /wave, lam0=4562.
;
;-
;
;-
;
;-----------------------IF HELP IS NEEDED-------------------------------
;
IF(KEYWORD_SET(print_help)) THEN BEGIN
   doc_library, 'plotfluxem_debug'
   return
ENDIF
;
IF (N_ELEMENTS(dir) EQ 0) THEN BEGIN
   PRINT, 'SYNTAX: '
   PRINT, 'plotFLUXEM_DEBUG, dir'
   PRINT, 'with dir=string , e.g. plotFLUXEM_DEBUG, "M1" '
   RETURN
ENDIF
IF (NOT KEYWORD_SET(obs_indx)) THEN BEGIN
   obs_indx=1
ENDIF
;
IF (KEYWORD_SET(freq) AND KEYWORD_SET(wave)) THEN BEGIN
   PRINT, 'ONLY ONE KEYWORD freq OR wave IS ALLOWED'
   RETURN
ENDIF
;
IF (KEYWORD_SET(freq) AND KEYWORD_SET(vel)) THEN BEGIN
   PRINT, 'ONLY ONE KEYWORD freq OR vel IS ALLOWED'
   RETURN
ENDIF
;
IF (KEYWORD_SET(vel) AND KEYWORD_SET(wave)) THEN BEGIN
   PRINT, 'ONLY ONE KEYWORD vel OR wave IS ALLOWED'
   RETURN
ENDIF
;
IF (KEYWORD_SET(vel)) THEN BEGIN
   IF(NOT KEYWORD_SET(vth_fid)) THEN BEGIN
      PRINT, 'WHEN PLOTS ARE PERFORMED IN VELOCITY-REGIME'
      PRINT, '   KEYWORD vth_fid HAS TO BE SPECIFIED'
      RETURN
   ENDIF
ENDIF
;
IF (KEYWORD_SET(freq) OR KEYWORD_SET(wave)) THEN BEGIN
   IF(NOT KEYWORD_SET(lam0)) THEN BEGIN
      PRINT, 'WHEN PLOTS ARE PERFORMED IN FREQUENCY OR WAVELENGTH-REGIME'
      PRINT, '   TRANSITION WAVELENGTH lam0 HAS TO BE SPECIFIED'
      RETURN
   ENDIF
   IF(NOT KEYWORD_SET(vth_fid)) THEN BEGIN
      PRINT, 'WHEN PLOTS ARE PERFORMED IN FREQUENCY OR WAVELENGTH-REGIME'
      PRINT, '   KEYWORD vth_fid HAS TO BE SPECIFIED'
      RETURN
   ENDIF
ENDIF
;
;-----------------------------------------------------------------------
;
getFLUXEM_DEBUG, dir, obs_indx, FILE_EXIST, NXOBS, ALPHA, GAMMA, XOBS, FLUX_TOT, FLUX_CONT, FLUX_EMI, FLUX_ABS
;
;-----------------------------------------------------------------------
;
IF(FILE_EXIST EQ 0) THEN BEGIN
   PRINT, 'ERROR IN OPENING FILES: RETURNING'
   RETURN
ENDIF
;
IF(KEYWORD_SET(XSCALEFAC)) THEN BEGIN
   XOBS=XOBS*XSCALEFAC
ENDIF
;
;---------------------CALCULATE PLOT RANGES-----------------------------
;
IF(NOT KEYWORD_SET(ylim)) THEN BEGIN
   ymin1=MIN(FLUX_TOT/FLUX_CONT)
   ymin2=MIN(FLUX_EMI/FLUX_CONT)
   ymin3=MIN(FLUX_ABS/FLUX_CONT)
   ymax1=MAX(FLUX_TOT/FLUX_CONT)
   ymax2=MAX(FLUX_EMI/FLUX_CONT)
   ymax3=MAX(FLUX_ABS/FLUX_CONT)

   ylim=[MIN([ymin1, ymin2, ymin3]), MAX([ymax1, ymax2, ymax3])]
;   ylim=[0.d0,2.5d0]
ENDIF
;
;
xarr=xobs
;
IF(KEYWORD_SET(freq)) THEN BEGIN
   nue0=lam0*1.d-8/!cgs_clight
   xarr=vth_fid*xobs*nue0/!cgs_clight + nue0
ENDIF
;
IF(KEYWORD_SET(wave)) THEN BEGIN
   xarr = lam0/(vth_fid*1.d5*xobs/!cgs_clight + 1.d0)
ENDIF
;
IF(KEYWORD_SET(vel)) THEN BEGIN
   xarr = vth_fid*xobs
ENDIF
;
IF(NOT KEYWORD_SET(xlim)) THEN BEGIN
   xlim=[MAX(xarr),MIN(xarr)]
   IF(KEYWORD_SET(wave)) THEN BEGIN
      xlim=[MIN(xarr),MAX(xarr)]
   ENDIF   
;   xlim=[-1.25d0, 1.25d0]
ENDIF
;
;------------------------TITLES-----------------------------------------
;
print, 'alpha', alpha
print, 'gamma', gamma
nobs_x=SIN(ALPHA)*COS(GAMMA)
nobs_y=SIN(ALPHA)*SIN(GAMMA)
nobs_z=COS(ALPHA)
;
;on default
xtitleStr=textoidl('x_{obs} in units of v_\infty')
;
IF(KEYWORD_SET(freq)) THEN BEGIN
   xtitleStr=textoidl('\nue in Hz')
ENDIF
;
IF(KEYWORD_SET(wave)) THEN BEGIN
   xtitleStr=textoidl('\lambda in A')
ENDIF
;
IF(KEYWORD_SET(vel)) THEN BEGIN
   xtitleStr=textoidl('v in km/s')
ENDIF
;
ytitleStr=textoidl('F / F_{c}')
titleStr=textoidl('observers direction n=') + STRING(nobs_x, FORMAT='(F9.5)') + ', ' $
                                            + STRING(nobs_y, FORMAT='(F9.5)') + ', ' $
                                            + STRING(nobs_z, FORMAT='(F9.5)')
;
;----------------------PLOT---------------------------------------------
;
IF KEYWORD_SET(ONAME) THEN BEGIN
;set output file
   PRINT, "WRITING OUTPUT TO: ", ONAME
   set_plot,'ps'
   device,file=ONAME
   loadct, 0
   IF(KEYWORD_SET(OPLOT)) THEN BEGIN
      PRINT, 'OVERPLOT NOT ALLOWED IN PS-MODE'
      RETURN
   ENDIF
ENDIF ELSE BEGIN
;set output to window
   IF(NOT KEYWORD_SET(WINDX)) THEN BEGIN
      windx=0
   ENDIF
;check if overplot is used
   IF(KEYWORD_SET(OPLOT)) THEN BEGIN
;check if color index is specified
      IF(KEYWORD_SET(CINDX)) THEN BEGIN
         loadct, 12
      ENDIF ELSE BEGIN
         CINDX=228
         loadct, 0
      ENDELSE
   ENDIF ELSE BEGIN
      window, windx
      device, decomposed=0
      loadct, 0
   ENDELSE
ENDELSE
;
;-----------------------------------------------------------------------
;
IF(KEYWORD_SET(OPLOT)) THEN BEGIN
   OPLOT, XARR, FLUX_TOT/FLUX_CONT, $
      line=0, $
      color=CINDX
   OPLOT, XARR, FLUX_EMI/FLUX_CONT, $
      line=1, $
      color=CINDX
   OPLOT, XARR, FLUX_ABS/FLUX_CONT, $
      line=2, $
      color=CINDX
ENDIF ELSE BEGIN
   PLOT, XARR, FLUX_TOT/FLUX_CONT, $
      line=0, $
      yrange=ylim, $
      xrange=xlim, $
      title=titleStr, $
      xtitle=xtitleStr, $
      isotropic=isotropic, $
      ytitle=ytitleStr, /XSTYLE, /YSTYLE, $
      charsize=charsize
   OPLOT, XARR, FLUX_EMI/FLUX_CONT, $
      line=1, $
      color=CINDX
   OPLOT, XARR, FLUX_ABS/FLUX_CONT, $
      line=2, $
      color=CINDX

ENDELSE
;

IF KEYWORD_SET(ONAME) THEN BEGIN
   device, /close
   set_plot,'x'
ENDIF
;
;-----------------------------------------------------------------------
;
END

