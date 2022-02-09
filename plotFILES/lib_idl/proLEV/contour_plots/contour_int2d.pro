pro contour_int2d, fname=fname, windx=windx, xlim=xlim, ylim=ylim, clim=clim,  log=log, oname=oname, ogrid=ogrid, new=NEW, help=print_help, emi=emi, abs=abs, v2d=v2d, t2d=t2d
;
;+
; NAME:
;       contour_int2d
;
; PURPOSE:
;       This procedure plots the intensity and optical depth along a ray (on a slice)
;
; CALLING SEQUENCE:
;
;       contour_int2d
;
; INPUTS:
;
; KEYWORD PARAMETERS:
;
;       windx:  Set this keyword to an integer defining the window, to which
;               output is plotted
;
;       fname:  Set this keyword to a string that defines the file name
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
;       new:    Set this keyword (flag) to read in new version of input-file:
;               including vn_2d
;
;       emi:    Set this keyword (flag) to plot only emitted intensity
;
;       abs:    Set this keyword (flag) to plot only absorbed intensity
;
;       v2d:    Set this keyword (flag) to plot only velocity
;
;       t2d:    Set this keyword (flag) to plot only optical depth
;
;       help:   Set this keyword (flag) to show the documentation of the
;               procedure
;
; EXAMPLE:
;       contour_int2d, '.', ylim=[0.5,-0.5], xlim=[0.,2.5], clim=[0.,1.], /LOG
;
;-
;
;-----------------------IF HELP IS NEEDED-------------------------------
;
IF(KEYWORD_SET(print_help)) THEN BEGIN
   doc_library, 'contour_int2d'
   return
ENDIF
;
;-----------------------------------------------------------------------
;
if(not keyword_set(fname)) then begin
   doc_library, 'contour_int2d'
   return
ENDIF
;
if(FILE_TEST(fname) ne 1) then begin
   print, 'file does not exist'
   return
endif
;
file_id = h5f_open(fname)
;
   group_id = H5G_OPEN(file_id, 'dimensions')
      att_id=H5A_OPEN_NAME(group_id, 'nz')
         nz=H5A_READ(att_id)
         nz=nz(0)
      h5a_close, att_id
      att_id=H5A_OPEN_NAME(group_id, 'nx')
         nx=H5A_READ(att_id)
         nx=nx(0)
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
      dset_id=h5d_open(group_id, 'x')
         x=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'z')
         z=h5d_read(dset_id)
      h5d_close, dset_id
   h5g_close, group_id
;
   group_id = h5g_open(file_id, 'along_ray')
      dset_id=h5d_open(group_id, 'int2d')
         int2d=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'iemi2d')
         iemi2d=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'iabs2d')
         iabs2d=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'tau2d')
         tau2d=h5d_read(dset_id)
      h5d_close, dset_id
      if(keyword_set(new)) then begin
         dset_id=h5d_open(group_id, 'vn2d')
            vn2d=h5d_read(dset_id)
         h5d_close, dset_id
      endif
;
   h5g_close, group_id
;
h5f_close, file_id
;
;------------------------plot triangles
;
;window, 0
;triangulate, x, z, tri
;ntri= (size(tri))[2]
;plot, x, z, psym=3, $
;   xrange=xlim, $
;   yrange=ylim

;for i=1, ntri-1 do begin
;   PLOTS, [x[tri[*,i]], x[tri[0,i]]], $
;         [z[tri[*,i]], z[tri[0,i]]]
;endfor

;stop

;nx = 2001                                    ;x−resolution of regular grid
;nz = 2001                                    ;y−resolution of regular grid

;int2d=int2d/xic1
;int2d=double(int2d)
;x=double(x)
;z=double(z)

;print, x(172,1)
;for i=0, nz-1 do begin
;   print, i, z(172,i), int2d(171,i), int2d(172,i), int2d(173,i)
;endfor

;iint2d = TRIGRID(x, z, int2d, tri, $                ;Interpolate to regular grid
;,NX    = nx, NY    = nz, $           ;Resolution of output grid
;XGRID = xx, YGRID = zz, $           ;Coordinates of output grid
;MISSING = !VALUES.F_NAN)            ;Points outside triangles are set to NaN
;x=xx
;z=zz
;int2d=iint2d
;stop
;
;--------------------DEFINE ALL PLOT-PROPERTIES-------------------------
;
;normalized intensity
int2d = int2d/XIC1
;
if(keyword_set(t2d)) then int2d=tau2d
if(keyword_set(v2d)) then int2d=vn2d
if(keyword_set(abs)) then int2d=iabs2d/xic1
if(keyword_set(emi)) then int2d=iemi2d/xic1
;
;-----------------------------------------------------------------------
;
;DEFINE RANGE
if(not keyword_set(xlim)) then begin
   xlim=[MIN(x),MAX(x)]
endif
if(not keyword_set(ylim)) then begin
   ylim=[MIN(z),MAX(z)]
endif
;
;DEFINE TITLE STRINGS
titleStr=textoidl('I_{em} (x_{obs}=')+string(xobs, format='(f10.4)')+textoidl(')')
xtitleStr=textoidl('x')
ytitleStr=textoidl('z')
;
;DEFINE CONTOUR-COLORS
if(keyword_set(log)) then begin
   min_color = min(int2d(where(int2d gt 0.)))
   max_color = max(int2d)
;
   min_color = alog10(min_color)
   max_color = alog10(max_color)
   int2d=alog10(int2d)
endif else begin
   min_color = min(int2d)
   max_color = max(int2d)
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
;
   get_contour_colors, colors, xobs, phase, ncolors, nlevels_iso, bottom, $
                       levels_final, levels_iso, c_colors_final,  $
                       cb_indx, cb_tickmark_name
;;              
;   ncolors_emi=ncolors
;   nlevels_iso_emi=nlevels_iso
;   bottom_emi=bottom
;   levels_final_emi=levels_final
;   levels_iso_emi=levels_iso
;   c_colors_final_emi=c_colors_final
;   cb_indx_emi=cb_indx
;   cb_tickmark_name_emi=cb_tickmark_name
;;
;   ncolors_abs=ncolors
;   nlevels_iso_abs=nlevels_iso
;   bottom_abs=bottom
;   levels_final_abs=levels_final
;   levels_iso_abs=levels_iso
;   c_colors_final_abs=c_colors_final
;   cb_indx_abs=cb_indx
;   cb_tickmark_name_abs=cb_tickmark_name
;   
endif else begin
;
   get_contour_colors, int2d, x, z, ncolors, nlevels_iso, bottom, $
                       levels_final, levels_iso, c_colors_final,  $
                       cb_indx, cb_tickmark_name
;;
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

CONTOURPLOTS, int2d, x, z, xlim, ylim, $
              titleStr, xtitleStr, ytitleStr, $
              ncolors, nlevels_iso, bottom, $
              levels_final, levels_iso, c_colors_final, $
              cb_indx, cb_tickmark_name, ctable=13, /ISOTROPIC, /IRREGULAR

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
end
