pro plot_surfb, fname=fname, windx=windx, oname=oname, help=print_help
;
;+
; NAME:
;	plot_surfb
;
; PURPOSE:
;	This procedure plots surface brightnes of spec.eo output
;
; CALLING SEQUENCE:
;
;	plot_surfb
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
   doc_library, 'plot_surfb'
   return
endif
;
if(n_elements(fname) eq 0) then begin
   print, 'error in plot_surfb: fname not specified'
   doc_library, 'plot_surfb'
   stop
endif
;
;------------------read all information from hdf5-file------------------
;
read_surfb, fname=fname, np=np, p=p, zeta=zeta, nzeta=nzeta, xic1=xic1, xobs=xobs, alpha=alpha, gamma=gamma, int2d_tot=int2d_tot, int2d_abs=int2d_abs, int2d_emi=int2d_emi, int2d_cont=int2d_cont
;
int2d_tot=int2d_tot/xic1
;
;-----------------------------------------------------------------------
;
xlim=[-2.,2.]
ylim=[-2.,2.]
clim=[0.,1.]
;
get_contour_colors2, clim, ncolors, nlevels_iso, bottom, $
                     levels_final, levels_iso, c_colors_final,  $
                     cb_indx, cb_tickmark_name
;
;-----------------------------------------------------------------------
;
if(not keyword_set(windx)) then windx=0
;
if(keyword_set(oname)) then begin
   oname='ps_files/contour_surfb.ps'
   print, "writing output to: ", oname
   set_plot,'ps'
   device,file=oname, xsize=19., ysize=26.7, xoffset=1., yoffset=1. , $
    decomposed=0, color=1, bits_per_pixel=8
endif else begin
   window, windx, xsize=950, ysize=1200
   device, decomposed=0
   windx=windx+1
endelse
;

;
;----------------------------TITLE-STRINGS------------------------------
;
titleStr=''
ctitleStr=textoidl('I/I_{c}')
xtitleStr=textoidl('x')
ytitleStr=textoidl('z')
;
x=fltarr(np, nzeta)
z=fltarr(np, nzeta)
;
for i=0, np-1 do begin
   for j=0, nzeta-1 do begin
      x(i,j) = p(i)*cos(zeta(j))
      z(i,j) = p(i)*sin(zeta(j))
   endfor
endfor
;
loadct, 0
;
contourplots_single, int2d_tot, x, z, $
              ncolors, nlevels_iso, bottom, $
              levels_final, levels_iso, c_colors_final, $
              cb_indx, cb_tickmark_name, $
              xlim=xlim, ylim=ylim, $
              titlestr=titlestr, xtitlestr=xtitlestr, ytitlestr=ytitlestr, $
              ctitlestr=ctitlestr, ctable=13, /isotropic
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
