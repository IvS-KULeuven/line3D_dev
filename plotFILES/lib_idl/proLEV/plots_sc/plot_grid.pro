pro plot_grid, fname, oname=oname, windx=windx, xlim=xlim, ylim=ylim, zlim=zlim, rad=rad, help=print_help
;
;+
; NAME:
;       plot_grid
; PURPOSE:
;       plots grid
;
; CALLING SEQUENCE:
;       plot_grid, fname
;
; INPUTS:
;       fname:   file name where things are stored
;
; KEYWORDS:
;       oname:   output name (for output to ps-file)
;       windx:   window index
;       xlim, ylim, zlim:  plot ranges for x, y and z-axes
;       rad:     overplotting a circle with radius=rad
;
;-
;
;-----------------------------------------------------------------------
;
if(keyword_set(print_help)) then begin
   doc_library, 'plot_grid'
   return
endif
;
if(n_elements(fname) eq 0) then begin
   print, 'fname not specified'
   stop
endif
;
fexist=file_test(fname)
if(fexist eq 0) then begin
   print, 'file does not exist: ', fname
   stop
endif
;
get_axes, fname, x, y, z
;
;--------------------default ranges-------------------------------------
;
if(keyword_set(xlim)) then begin
   if(n_elements(xlim ne 2)) then begin
      print, 'keyword xlim needs to be 1d-array with two elements'
      return
   endif
endif else begin
   xlim=[min(x),max(x)]
endelse
;
if(keyword_set(ylim)) then begin
   if(n_elements(ylim ne 2)) then begin
      print, 'keyword ylim needs to be 1d-array with two elements'
      return
   endif
endif else begin
   ylim=[min(y),max(y)]
endelse
;
if(keyword_set(zlim)) then begin
   if(n_elements(zlim ne 2)) then begin
      print, 'keyword zlim needs to be 1d-array with two elements'
      return
   endif
endif else begin
   zlim=[min(z),max(z)]
endelse
;
;-----------------------------------------------------------------------
;
if keyword_set(oname) then begin
   print, "writing output to: ", oname
   set_plot,'ps'
   device,file=oname, decomposed=0
endif else begin
   if(not keyword_set(windx)) then begin
      windx=0
   endif
   window, windx, xsize=1200, ysize=400
   device, decomposed=0
   set_plot, 'x'
endelse
;
;
titlestr='grid'
;
;x-y-plane
!p.multi=[3,3,1]
xtitlestr=textoidl('x')
ytitlestr=textoidl('y')
plotgrid_single, x, y, xtitlestr, ytitlestr, titlestr, xlim, ylim, rad=rad

;x-z-plane
!p.multi=[2,3,1]
xtitlestr=textoidl('x')
ytitlestr=textoidl('z')
plotgrid_single, x, z, xtitlestr, ytitlestr, titlestr, xlim, zlim, rad=rad

;y-z-plane
!p.multi=[1,3,1]
xtitlestr=textoidl('y')
ytitlestr=textoidl('z')
plotgrid_single, y, z, xtitlestr, ytitlestr, titlestr, ylim, zlim, rad=rad



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
