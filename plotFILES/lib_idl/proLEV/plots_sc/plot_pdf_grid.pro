pro plot_pdf_grid, dir=dir, ylim=ylim, xlim=xlim, oname=oname, windx=windx, help=print_help
;
;+
; NAME:
;       plot_pdf_grid
; PURPOSE:
;       plots probability density of created x,y,z grids
;
; CALLING SEQUENCE:
;       plot_pdf_grid
;
; INPUTS:
;
; KEYWORDS:
;       dir:     directory where grid-data is stored
;       oname:   output name (for output to ps-file)
;       windx:   window index
;       xlim:    plot range for x-axis
;       ylim:    plot range for y-axis
;-
;
;-----------------------------------------------------------------------
;
if(keyword_set(print_help)) then begin
   doc_library, 'plot_pdf_grid'
   return
endif
;
if(not keyword_set(dir)) then dir='.'
;
;-----------------------------------------------------------------------
;
if keyword_set(oname) then begin
   print, "writing output to: ", oname
   set_plot,'ps'
   device,file=oname, decomposed=0
endif else begin
   if(not keyword_set(windx)) then windx=0
   window, windx
   device, decomposed=0
   windx=windx+1
   set_plot, 'x'
endelse
;
readcol, dir+'/pdf_rgridf.dat', r_mid, pdf_r, p_r
readcol, dir+'/pdf_phigridf.dat', phi_mid, pdf_phi, p_phi
readcol, dir+'/pdf_xgridf_c.dat', xfine_mid_c, pdf_xfine_c, p_xfine_c
readcol, dir+'/pdf_xgridf_nc.dat', xfine_mid_nc, pdf_xfine_nc, p_xfine_nc
readcol, dir+'/pdf_xgridf.dat', xfine_mid, pdf_xfine, p_xfine
readcol, dir+'/pdf_xgrid1.dat', xfinal1_mid, pdf1_xfinal, p1_xfinal
readcol, dir+'/pdf_xgrid2.dat', xfinal2_mid, pdf2_xfinal, p2_xfinal
readcol, dir+'/pdf_xgrid3.dat', xfinal3_mid, pdf3_xfinal, p3_xfinal
;
;
loadct_lev, ci_black, ci_blue, ci_red, ci_cyan, ci_green, ci_magenta, ci_yellow, ci_white
;
titlestr=textoidl('probability density for grids')
xtitlestr=textoidl('x_{mid}, r_{mid}')
ytitlestr=textoidl('pdf(x), pdf(r)')
;

if(not keyword_set(xlim)) then xlim=[0.,max(xfine_mid)]
if(not keyword_set(ylim)) then ylim=[0.,max([pdf_xfine,pdf_xfine_nc,pdf_xfine_c])]

plot, [0.,0.], [0.,0.], $
      yrange=ylim, $
      xrange=xlim, $
      xtitle=xtitlestr, $
      ytitle=ytitlestr, $
      title=titlestr, $
      charsize=1.5
;
oplot, r_mid, pdf_r, line=0
oplot, xfine_mid_c, pdf_xfine_c, color=ci_blue
oplot, xfine_mid_nc, pdf_xfine_nc, color=ci_blue
oplot, xfine_mid, pdf_xfine, color=ci_red
oplot, xfinal1_mid, pdf1_xfinal, color=ci_green
oplot, xfinal2_mid, pdf2_xfinal, color=ci_magenta, psym=1
oplot, xfinal3_mid, pdf3_xfinal, color=ci_yellow
;
lstr1  = textoidl('pdf(r)')
lstr2  = textoidl('pdf(x_c), pdf(x_{nc})')
lstr3  = textoidl('pdf(x)')
lstr4  = textoidl('pdf(x_{final1})')
lstr5  = textoidl('pdf(x_{final2})')
lstr6  = textoidl('pdf(x_{final3})')

legend, [lstr1, lstr2, lstr3, lstr4, lstr5, lstr6], $
        psym=[0,0,0,0,0,0], $
        linestyle=[0,0,0,0,0,0], $
        color=[ci_white, ci_blue, ci_red, ci_green, ci_magenta, ci_yellow], $
        textcolor=ci_white, $
        charsize=1.5, $
        /right_legend
;
if keyword_set(oname) then begin
   device, /close
   set_plot,'x'
endif


end
;
;-----------------------------------------------------------------------
;-----------------------------------------------------------------------
;-----------------------------------------------------------------------
;
pro plotgrid_single, xarr, yarr, xtitlestr, ytitlestr, titlestr, xlimits, ylimits, rad=rad
;
ndx=n_elements(xarr)
ndy=n_elements(yarr)
;
ncircle=100
phi=findgen(ncircle)*2.*!pi/(ncircle-1)
zcircle=cos(phi)
xcircle=sin(phi)
;
if(keyword_set(rad)) then begin
   zcircle2=rad*cos(phi)
   xcircle2=rad*sin(phi)
endif
;
plot, [xarr(0),xarr(0)], [min(yarr), max(yarr)], $
      /isotropic, $
      xrange=xlimits, $
      yrange=ylimits, $
      line=0, $
      xtitle=xtitlestr, $
      ytitle=ytitlestr, $
      title=titlestr, $
      /xstyle, /ystyle, $
      charsize=2.
;
oplot, [min(xarr),max(xarr)], [yarr(0),yarr(0)], $
      line=0

for i=1, ndx-1 do begin
   oplot, [xarr(i),xarr(i)], [min(yarr), max(yarr)], $
   line=0
endfor
;
for i=1, ndy-1 do begin
   oplot, [min(xarr), max(xarr)], [yarr(i), yarr(i)], $
   line=0
endfor

oplot, zcircle,xcircle
;
if(keyword_set(rad)) then begin
   oplot, zcircle2, xcircle2
endif
;
end
