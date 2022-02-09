pro plot_surfb_vbin, fname=fname, windx=windx, oname=oname, xlim=xlim, ylim=ylim, clim=clim, logscale=logscale, help=print_help
;
;+
; NAME:
;	plot_surfb_vbin
;
; PURPOSE:
;	This procedure plots surface brightnes of spec_vbin.eo output
;
; CALLING SEQUENCE:
;
;	plot_surfb_vbin
;
; INPUTS:
;
; KEYWORD PARAMETERS:
;
;	fname:  Set this keyword to thee file-name of the model
;
;       help:   Set this keyword (flag) to show the documentation of the
;               procedure
;
; EXAMPLE:
;	plot_surfb, fname='output_surfb.h5'
;
;-
;
;-----------------------if help is needed-------------------------------
;
if(keyword_set(print_help)) then begin
   doc_library, 'plot_surfb_vbin'
   return
endif
;
if(n_elements(fname) eq 0) then begin
   print, 'error in plot_surfb_vbin: fname not specified'
   doc_library, 'plot_surfb_vbin'
   stop
endif
;
;------------------read all information from hdf5-file------------------
;
read_surfb_vbin, fname=fname, npoints=npoints, points_xcoord=points_xcoord, points_ycoord=points_ycoord, xic1=xic1, xic2=xic2, xobs=xobs, alpha=alpha, gamma=gamma, int2d_tot=int2d_tot, int2d_abs=int2d_abs, int2d_emi=int2d_emi
;
int2d_tot=int2d_tot/xic1
;
if(keyword_set(logscale)) then begin
   int2d_plot = int2d_tot
   indx1 = where(int2d_tot gt 0.)
   indx2 = where(int2d_tot le 0.)   
   int2d_plot(indx1) = alog10(int2d_tot(indx1))
   int2d_plot(indx2) = -15.d0
endif else begin
   int2d_plot = int2d_tot
endelse
;
;------------------------------------------------------------------------
;
if(not keyword_set(windx)) then windx=0
;
;------------------------------------------------------------------------
;
;if(keyword_set(oname)) then begin
;   oname='ps_files/spec_triangles.ps'
;   print, "writing output to: ", oname
;   set_plot,'ps'
;   device,file=oname, xsize=19., ysize=26.7, xoffset=1., yoffset=1. , $
;    decomposed=0, color=1, bits_per_pixel=8
;endif else begin
;   window, windx
;   device, decomposed=0
;   windx=windx+1
;endelse
;
;
;plot_triangles, '../outputFILES', obs_indx=1, oname=oname, xlim=xlim, ylim=ylim
;
;
;
;if keyword_set(oname) then begin
;   device, /close
;   set_plot,'x'
;endif
;
;-----------------------------------------------------------------------
;
if(not keyword_set(clim)) then begin
   cmin=min(int2d_plot)
   cmax=max(int2d_plot)
   clim=[cmin,cmax]
endif

;
get_contour_colors2, clim, ncolors, nlevels_iso, bottom, $
                     levels_final, levels_iso, c_colors_final,  $
                     cb_indx, cb_tickmark_name
;
;-----------------------------------------------------------------------
;
;
if(keyword_set(oname)) then begin
   oname='ps_files/contour_surfb.ps'
   print, "writing output to: ", oname
   set_plot,'ps'
   device,file=oname, color=1, bits_per_pixel=8
endif else begin
   window, windx
   device, decomposed=0
   windx=windx+1
endelse
;
;----------------------------TITLE-STRINGS------------------------------
;
if(keyword_set(logscale)) then begin
   ctitleStr=textoidl('\log(I/I_{c})')
endif else begin
   ctitleStr=textoidl('I/I_{c}')
endelse
xtitleStr=textoidl('x_p')
ytitleStr=textoidl('z_p')
;
;
contourplots_single, int2d_plot, points_xcoord, points_ycoord, $
              ncolors, nlevels_iso, bottom, $
              levels_final, levels_iso, c_colors_final, $
              cb_indx, cb_tickmark_name, $
              xlim=xlim, ylim=ylim, $
              titlestr=titlestr, xtitlestr=xtitlestr, ytitlestr=ytitlestr, $
              ctitlestr=ctitlestr, ctable=13, /background2, /isotropic, /irregular
;
;
if keyword_set(oname) then begin
   device, /close
   set_plot,'x'
endif
;
;
;
;
end
