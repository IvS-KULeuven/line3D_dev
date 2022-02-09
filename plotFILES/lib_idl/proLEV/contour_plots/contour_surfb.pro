pro contour_surfb, fname=fname, windx=windx, xlim=xlim, ylim=ylim, clim=clim,  log=log, oname=oname, ogrid=ogrid, help=print_help
;
;+
; NAME:
;       contour_surfb
;
; PURPOSE:
;       This procedure plots the surface brightness at a certain wavelength
;
; CALLING SEQUENCE:
;
;       contour_surb, dir
;
; INPUTS:
;
; KEYWORD PARAMETERS:
;
;       windx:  Set this keyword to an integer defining the window, to which
;               output is plotted
;
;       fname:  Set this keyword to a string defining the file name
;
;       oname:  Set this keyword to a string, which defines output-file
;
;       ogrid:  Set this keyword (flag) to overplot the grid
;
;       ylim:   Set this keyword to a 2-d array, which sets the yrange
;
;       xlim:   Set this keyword to a 2-d array, which sets the xrange
;
;       clim:   Set this keyword to a 2-d array, whicht sets the contour-color-range
;
;       log:    Set this keyword (flag) to plot contours in log-space
;
;       help:   Set this keyword (flag) to show the documentation of the
;               procedure
;
; EXAMPLE:
;       plotfluxem, '.', ylim=[0.5,-0.5], xlim=[0.,2.5], clim=[0.,1.], /LOG
;
;-
;
;-----------------------IF HELP IS NEEDED-------------------------------
;
IF(KEYWORD_SET(print_help)) THEN BEGIN
   doc_library, 'contour_surfb'
   return
ENDIF
;
;-----------------------------------------------------------------------
;
if(not keyword_set(fname)) then begin
   doc_library, 'contour_surfb'
   return
endif
;
if(FILE_TEST(fname) ne 1) then begin
   print, 'file does not exist'
   return
endif
;
file_id = h5f_open(fname)
;
   group_id = H5G_OPEN(file_id, 'dimensions')
      att_id=H5A_OPEN_NAME(group_id, 'nzeta')
         nzeta=H5A_READ(att_id)
         nzeta=nzeta(0)
      h5a_close, att_id
      att_id=H5A_OPEN_NAME(group_id, 'np')
         np=H5A_READ(att_id)
         np=np(0)
      h5a_close, att_id
   h5g_close, group_id
;
   group_id = H5G_OPEN(file_id, 'parameter')
      att_id=H5A_OPEN_NAME(group_id, 'xic')
         xic1=H5A_READ(att_id)
         xic1=xic1(0)
      h5a_close, att_id
      att_id=H5A_OPEN_NAME(group_id, 'alpha')
         alpha=H5A_READ(att_id)
         alpha=alpha(0)
      h5a_close, att_id
      att_id=H5A_OPEN_NAME(group_id, 'gamma')
         gamma=H5A_READ(att_id)
         gamma=gamma(0)
      h5a_close, att_id
      att_id=H5A_OPEN_NAME(group_id, 'xobs')
         xobs=H5A_READ(att_id)
         xobs=xobs(0)
      h5a_close, att_id
   h5g_close, group_id
;
   group_id = h5g_open(file_id, 'coordinates')
      dset_id=h5d_open(group_id, 'p')
         p=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'zeta')
         zeta=h5d_read(dset_id)
      h5d_close, dset_id
   h5g_close, group_id
;
   group_id = h5g_open(file_id, 'surface')
      dset_id=h5d_open(group_id, 'intensity2d')
         Iem=h5d_read(dset_id)
      h5d_close, dset_id
;
;------new version: emission part and absorption part explicitly--------
;comment following lines for old version
      dset_id=h5d_open(group_id, 'intensity2d_emi')
         Iemi=h5d_read(dset_id)
      h5d_close, dset_id

      dset_id=h5d_open(group_id, 'intensity2d_abs')
         Iabs=h5d_read(dset_id)
      h5d_close, dset_id
;
   h5g_close, group_id
;
h5f_close, file_id
;
;--------------------DEFINE ALL PLOT-PROPERTIES-------------------------
;
;normalized intensity
Iem = Iem/XIC1
Iemi = Iemi/XIC1
Iabs = Iabs/XIC1
;
;DEFINE X,Y-AXIS
x=fltarr(np,nzeta)*0.d0
y=fltarr(np,nzeta)*0.d0
for i=0, np-1 do begin
   for j=0, nzeta-1 do begin
      x(i,j) = p(i)*cos(zeta(j))
      y(i,j) = p(i)*sin(zeta(j))
   endfor
endfor
;
;DEFINE RANGE
if(not keyword_set(xlim)) then begin
   xlim=[-MAX(p),MAX(p)]
endif
if(not keyword_set(ylim)) then begin
   ylim=[-MAX(p),MAX(p)]
endif
;
;DEFINE TITLE STRINGS
titleStr=textoidl('I_{em} (x_{obs}=')+string(xobs, format='(f10.4)')+textoidl(')')
xtitleStr=textoidl('x')
ytitleStr=textoidl('y')
;
;DEFINE CONTOUR-COLORS
if(keyword_set(log)) then begin
   min_color = min(Iem(where(iem gt 0.)))
   max_color = max(Iem)

   min_color = alog10(min_color)
   max_color = alog10(max_color)
   Iem=alog10(Iem)
   Iemi=alog10(Iemi)
   Iabs=alog10(Iabs)
endif else begin
   min_color = min(Iem)
   max_color = max(Iem)
endelse
;
if(keyword_set(clim)) then begin
   min_color = clim(0)
   max_color = clim(1)
   nx_colors=201
   ny_colors=201
   colors=fltarr(nx_colors, ny_colors)
   for i=0, nx_colors-1 do begin
      for j=0, ny_colors-1 do begin
         ran= min_color + RANDOMU(seed, /UNIFORM)*(max_color-min_color)
         colors(i,j) = ran
      endfor
   endfor

   get_contour_colors, colors, xobs, phase, ncolors, nlevels_iso, bottom, $
                       levels_final, levels_iso, c_colors_final,  $
                       cb_indx, cb_tickmark_name
;              
   ncolors_emi=ncolors
   nlevels_iso_emi=nlevels_iso
   bottom_emi=bottom
   levels_final_emi=levels_final
   levels_iso_emi=levels_iso
   c_colors_final_emi=c_colors_final
   cb_indx_emi=cb_indx
   cb_tickmark_name_emi=cb_tickmark_name
;
   ncolors_abs=ncolors
   nlevels_iso_abs=nlevels_iso
   bottom_abs=bottom
   levels_final_abs=levels_final
   levels_iso_abs=levels_iso
   c_colors_final_abs=c_colors_final
   cb_indx_abs=cb_indx
   cb_tickmark_name_abs=cb_tickmark_name
   
endif else begin
;
   get_contour_colors, Iem, p, zeta, ncolors, nlevels_iso, bottom, $
                       levels_final, levels_iso, c_colors_final,  $
                       cb_indx, cb_tickmark_name
;
   get_contour_colors, Iemi, p, zeta, ncolors_emi, nlevels_iso_emi, bottom_emi, $
                       levels_final_emi, levels_iso_emi, c_colors_final_emi,  $
                       cb_indx_emi, cb_tickmark_name_emi
;
   get_contour_colors, Iabs, p, zeta, ncolors_abs, nlevels_iso_abs, bottom_abs, $
                       levels_final_abs, levels_iso_abs, c_colors_final_abs,  $
                       cb_indx_abs, cb_tickmark_name_abs
;
endelse
;
;-------------------------MAKE THE PLOT---------------------------------
;
if(not keyword_set(windx)) then begin
   windx=0
endif

if(keyword_set(oname)) then begin
   print, 'writing output to: ', oname
   set_plot,'ps'
   device, file=ONAME, $ ;XSIZE=19., YSIZE=26.7, XOFFSET=1., YOFFSET=1. , $
   decomposed=0, color=1, BITS_PER_PIXEL=8
ENDIF ELSE BEGIN
   window, windx
   device, decomposed=0
   windx=windx+1
ENDELSE

;loadct, 0

CONTOURPLOTS, Iem, x, y, xlim, ylim, $
              titleStr, xtitleStr, ytitleStr, $
              ncolors, nlevels_iso, bottom, $
              levels_final, levels_iso, c_colors_final, $
              cb_indx, cb_tickmark_name, ctable=13, /ISOTROPIC

if(keyword_set(ogrid)) then begin
   oplot, x, y, psym=3
endif

IF KEYWORD_SET(ONAME) THEN BEGIN
   device, /close
   set_plot,'x'
ENDIF
;
;-----------------------------------------------------------------------
;
;emission and absorption part only on window
window, windx
device, decomposed=0
windx=windx+1

titleStr=textoidl('I_{emi} (x_{obs}=')+string(xobs, format='(f10.4)')+textoidl(')')

CONTOURPLOTS, Iemi, x, y, xlim, ylim, $
              titleStr, xtitleStr, ytitleStr, $
              ncolors_emi, nlevels_iso_emi, bottom_emi, $
              levels_final_emi, levels_iso_emi, c_colors_final_emi, $
              cb_indx_emi, cb_tickmark_name_emi, ctable=13, /ISOTROPIC
;
window, windx
device, decomposed=0
windx=windx+1

titleStr=textoidl('I_{abs} (x_{obs}=')+string(xobs, format='(f10.4)')+textoidl(')')

CONTOURPLOTS, Iabs, x, y, xlim, ylim, $
              titleStr, xtitleStr, ytitleStr, $
              ncolors_abs, nlevels_iso_abs, bottom_abs, $
              levels_final_abs, levels_iso_abs, c_colors_final_abs, $
              cb_indx_abs, cb_tickmark_name_abs, ctable=13, /ISOTROPIC
;
;-----------------------------------------------------------------------
;
end
