pro plot_photprof, windx=windx, oname=oname, xlim=xlim,ylim=ylim, isotropic=istoropic
;
;
fname01 = '../outputFILES/photprof_star01.dat'
fname02 = '../outputFILES/photprof_star02.dat'  
;
readcol, fname01, xobs01, xic01, xicc01
readcol, fname02, xobs02, xic02, xicc02
;
loadct_lev, ci_black, ci_blue, ci_red, ci_cyan, ci_green, ci_magenta, ci_yellow, ci_white
;
if(not keyword_set(xlim)) then xlim=[max([xobs01,xobs02]),min([xobs01,xobs02])]
if(not keyword_set(ylim)) then ylim=[0.d0,1.1d0]  
if(not keyword_set(isotropic)) then isotropic=0
if(not keyword_set(windx)) then windx=0

if(keyword_set(oname)) then begin
   oname='ps_files/spec_photprof.ps'
   print, "writing output to: ", oname
   set_plot,'ps'
   device,file=oname, color=1, bits_per_pixel=8
endif else begin
   window, windx
   device, decomposed=0
   windx=windx+1
endelse
  
plot, [xlim(0),xlim(0)], [ylim(0),ylim(0)], psym=3, $
      isotropic=isotropic, $
      xtitle=textoidl('x_{obs}[v_{th}^*'), $
      ytitle=textoidl('F_L/F_c'), $
      xrange=xlim, /xs, $
      yrange=ylim, /ys
oplot, xobs01, xic01/xicc01, color=ci_blue
oplot, xobs02, xic02/xicc02, color=ci_red
;
if keyword_set(oname) then begin
   device, /close
   set_plot,'x'
endif
  
end
