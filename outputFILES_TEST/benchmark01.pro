pro benchmark01, oname=oname
;
;read model data and solution from fortran routines
;
read_model2d, ndxmax, ndzmax, x, z, x_u2d, z_u2d, x_d2d, z_d2d, int2d_sc, abs2d_sc, contr2d_sc, $
              int2d_fvm, abs2d_fvm, contr2d_fvm, $
              xu_opac2d, zu_opac2d, xd_opac2d, zd_opac2d, xu_int2d, zu_int2d, scont2d, scont_u2d, $
              scont_d2d, opac2d, opac_u2d, opac_d2d
read_model1d, nr, r, int1d, abs1d, contr1d, int1d_sc, abs1d_sc, contr1d_sc, $
              int1d_fvm, abs1d_fvm, contr1d_fvm, int1d_scray, int1d_fvmray, opac1d, scont1d
;
;
;-----------------------------plot the complete grid--------------------
;
;define filled circles as plot symbols
a = findgen(17) * (!pi*2/16.)
usersym, 2.*cos(a), 2.*sin(a), /fill ;radius 2.
;
windx=0
window, windx
device, decomposed=0
windx=windx+1
;
loadct, 0
plot_grid, x, z, ndxmax, ndzmax, xlim=[-1.5,1.5], ylim=[-1.5,1.5]
;
;oplot upwind and downwind points
loadct_lev, ci_black, ci_blue, ci_red, ci_cyan, ci_green, ci_magenta, ci_yellow, ci_white
for i=1, ndxmax-2 do begin
   for k=1, ndzmax-2 do begin

      oplot, [x(i),x_u2d(i,k)], [z(k),z_u2d(i,k)], color=ci_blue
      oplot, [x(i),x_d2d(i,k)], [z(k),z_d2d(i,k)], color=ci_red

;oplot points required for interpolation
      oplot, [xu_opac2d(i,k,0), xu_opac2d(i,k,1), xu_opac2d(i,k,2), xu_opac2d(i,k,3), xu_opac2d(i,k,4), xu_opac2d(i,k,5), xu_opac2d(i,k,6), xu_opac2d(i,k,7)], $
             [zu_opac2d(i,k,0), zu_opac2d(i,k,1), zu_opac2d(i,k,2), zu_opac2d(i,k,3), zu_opac2d(i,k,4), zu_opac2d(i,k,5), zu_opac2d(i,k,6), zu_opac2d(i,k,7)], $
             color=ci_blue, psym=8
      oplot, [xd_opac2d(i,k,0), xd_opac2d(i,k,1), xd_opac2d(i,k,2), xd_opac2d(i,k,3), xd_opac2d(i,k,4), xd_opac2d(i,k,5), xd_opac2d(i,k,6), xd_opac2d(i,k,7)], $
             [zd_opac2d(i,k,0), zd_opac2d(i,k,1), zd_opac2d(i,k,2), zd_opac2d(i,k,3), zd_opac2d(i,k,4), zd_opac2d(i,k,5), zd_opac2d(i,k,6), zd_opac2d(i,k,7)], $
             color=ci_yellow, psym=2
      oplot, [xu_int2d(i,k,0), xu_int2d(i,k,1), xu_int2d(i,k,2), xu_int2d(i,k,3), xu_int2d(i,k,4), xu_int2d(i,k,5), xu_int2d(i,k,6), xu_int2d(i,k,7)], $
             [zu_int2d(i,k,0), zu_int2d(i,k,1), zu_int2d(i,k,2), zu_int2d(i,k,3), zu_int2d(i,k,4), zu_int2d(i,k,5), zu_int2d(i,k,6), zu_int2d(i,k,7)], $
             color=ci_green, psym=2
;      read, test
   endfor
endfor

legstr1='interpolation opacity/sources upwind'
legstr2='interpolation opacity/sources downwind'
legstr3='interpolation intensity'
legend, [legstr1, legstr2, legstr3], $
        psym=[8,2,2], $
        color=[ci_blue,ci_yellow,ci_green]
;stop
;
;-----------to plot only special points for a better overview-----------
;
if(ndxmax gt 24 or ndzmax gt 24) then goto, jump
;
window, windx
device, decomposed=0
windx=windx+1
;
;
;is=max(where(x lt -1.5d0))
;ie=min(where(x gt 1.5d0))
;ks=max(where(z lt -1.5d0))
;ke=min(where(z gt 1.5d0))
is=1
ie=ndxmax-2
ks=1
ke=ndzmax-2
;
for i=is+1, ie-1 do begin
  for k=ks+1, ke-1 do begin

      loadct, 0
      plot_grid, x, z, ndxmax, ndzmax, xlim=[-3.,3.], ylim=[-3.,3.]
      loadct_lev, ci_black, ci_blue, ci_red, ci_cyan, ci_green, ci_magenta, ci_yellow, ci_white
      oplot, [x(i),x_u2d(i,k)], [z(k),z_u2d(i,k)], color=ci_blue
      oplot, [x(i),x_d2d(i,k)], [z(k),z_d2d(i,k)], color=ci_red

;oplot points required for interpolation
      oplot, [xu_opac2d(i,k,0), xu_opac2d(i,k,1), xu_opac2d(i,k,2), xu_opac2d(i,k,3), xu_opac2d(i,k,4), xu_opac2d(i,k,5), xu_opac2d(i,k,6), xu_opac2d(i,k,7)], $
             [zu_opac2d(i,k,0), zu_opac2d(i,k,1), zu_opac2d(i,k,2), zu_opac2d(i,k,3), zu_opac2d(i,k,4), zu_opac2d(i,k,5), zu_opac2d(i,k,6), zu_opac2d(i,k,7)], $
             color=ci_blue, psym=8
      oplot, [xd_opac2d(i,k,0), xd_opac2d(i,k,1), xd_opac2d(i,k,2), xd_opac2d(i,k,3), xd_opac2d(i,k,4), xd_opac2d(i,k,5), xd_opac2d(i,k,6), xd_opac2d(i,k,7)], $
             [zd_opac2d(i,k,0), zd_opac2d(i,k,1), zd_opac2d(i,k,2), zd_opac2d(i,k,3), zd_opac2d(i,k,4), zd_opac2d(i,k,5), zd_opac2d(i,k,6), zd_opac2d(i,k,7)], $
             color=ci_yellow, psym=2
      oplot, [xu_int2d(i,k,0), xu_int2d(i,k,1), xu_int2d(i,k,2), xu_int2d(i,k,3), xu_int2d(i,k,4), xu_int2d(i,k,5), xu_int2d(i,k,6), xu_int2d(i,k,7)], $
             [zu_int2d(i,k,0), zu_int2d(i,k,1), zu_int2d(i,k,2), zu_int2d(i,k,3), zu_int2d(i,k,4), zu_int2d(i,k,5), zu_int2d(i,k,6), zu_int2d(i,k,7)], $
             color=ci_green, psym=2

      legstr1='interpolation opacity/sources upwind'
      legstr2='interpolation opacity/sources downwind'
      legstr3='interpolation intensity'
      legend, [legstr1, legstr2, legstr3], $
              psym=[8,2,2], $
              color=[ci_blue,ci_yellow,ci_green]
;      print, i, k
      read, test
;
   endfor
endfor
;stop
jump: print, 'skip explicit grid-plots because ndxmax, ndzmax too large'
;
;-----plot the opacity and source model as a function of the radius-----
;
loadct, 0
window, windx, xsize=950, ysize=800
device, decomposed=0
windx=windx+1
;
xlim=[0.9d0,sqrt(x(ndxmax-1)^2+z(ndzmax-1)^2)]
;
;opacities
;NOTE: at r_upwind=1.d0 (on star), extrapolation of opacities
;      and therefore lower opacities than theoretical values
ymax=max([max(opac2d),max(opac_u2d),max(opac_d2d),max(opac1d)])
;normalize to maximum value
opac2d=opac2d/ymax
opac_u2d=opac_u2d/ymax
opac_d2d=opac_d2d/ymax
opac1d=opac1d/ymax
ymax=ymax/ymax
ymin=min([min(opac2d),min(opac_u2d),min(opac_d2d),min(opac1d)])
ylim=[ymin-0.1*(ymax-ymin), ymax+0.1*(ymax-ymin)]

!p.multi=[2,1,2]
plot, r, opac1d, $
      xrange=xlim, /xs, $
;      xrange=[0.9,1.1], /xs, $
      yrange=ylim, $
      xtitle='r', $
      ytitle=textoidl('\chi'), $
      line=2
;
loadct_lev, ci_black, ci_blue, ci_red, ci_cyan, ci_green, ci_magenta, ci_yellow, ci_white
;
for i=0, ndxmax-1 do begin
   for k=0, ndzmax-1 do begin
      rp=sqrt(x(i)^2 + z(k)^2)
      oplot, [rp,rp], [opac2d(i,k),opac2d(i,k)], psym=3, color=ci_blue

      ru=sqrt(x_u2d(i,k)^2 + z_u2d(i,k)^2)
      oplot, [ru,ru], [opac_u2d(i,k), opac_u2d(i,k)], psym=3, color=ci_yellow

      rd=sqrt(x_d2d(i,k)^2 + z_d2d(i,k)^2)
      oplot, [rd,rd], [opac_d2d(i,k), opac_d2d(i,k)], psym=3, color=ci_red
   endfor
endfor
legstr0=textoidl('opacity 1d')
legstr1=textoidl('opacity at r_p')
legstr2=textoidl('opacity at r_u')
legstr3=textoidl('opacity at r_d')
legend, [legstr1, legstr2, legstr3], $
         psym=[3,3,3], $
         color=[ci_blue,ci_yellow,ci_red], $
         /bottom
;
;
;source functions
ymax=max([max(scont2d),max(scont_u2d),max(scont_d2d),max(scont1d)])
;normalize to maximum value
if(ymax ne 0.) then begin
   scont2d=scont2d/ymax
   scont_u2d=scont_u2d/ymax
   scont_d2d=scont_d2d/ymax
   scont1d=scont1d/ymax
   ymax=ymax/ymax
endif
ymin=min([min(scont2d),min(scont_u2d),min(scont_d2d),min(scont1d)])
ylim=[ymin-0.1*(ymax-ymin), ymax+0.1*(ymax-ymin)]
;
!p.multi=[1,1,2]
loadct, 0
plot, r, scont1d, $
      xrange=xlim, /xs, $
      yrange=ylim, $
      xtitle='r', $
      ytitle=textoidl('source function'), $
      line=2
;
loadct_lev, ci_black, ci_blue, ci_red, ci_cyan, ci_green, ci_magenta, ci_yellow, ci_white
;
for i=0, ndxmax-1 do begin
   for k=0, ndzmax-1 do begin
      rp=sqrt(x(i)^2 + z(k)^2)
      oplot, [rp,rp], [scont2d(i,k),scont2d(i,k)], psym=3, color=ci_blue

      ru=sqrt(x_u2d(i,k)^2 + z_u2d(i,k)^2)
      oplot, [ru,ru], [scont_u2d(i,k), scont_u2d(i,k)], psym=3, color=ci_yellow

      rd=sqrt(x_d2d(i,k)^2 + z_d2d(i,k)^2)
      oplot, [rd,rd], [scont_d2d(i,k), scont_d2d(i,k)], psym=3, color=ci_red
   endfor
endfor
legstr0=textoidl('source 1d')
legstr1=textoidl('source at r_p')
legstr2=textoidl('source at r_u')
legstr3=textoidl('source at r_d')
legend, [legstr1, legstr2, legstr3], $
         psym=[1,1,1], $
         color=[ci_blue,ci_yellow,ci_red], $
         /bottom

!p.multi=0
loadct, 0
;
;----------------------------contour plots------------------------------
;
window, windx, xsize=950, ysize=600
device, decomposed=0
windx=windx+1
;
;intensity
clim=[0.,1.]
xlim=[min(x),max(x)]
ylim=[min(z),max(z)]
titleStr1='intensity sc (normalized to max)'
titleStr2='absorption part sc (normalized to max)'
titleStr3='source contribution sc (normalized to max)'
titleStr4='intensity fvm (normalized to max)'
titleStr5='absorption part fvm (normalized to max)'
titleStr6='source contribution fvm (normalized to max)'
xtitleStr='x'
ytitleStr='z'
;
loadct, 0
get_contour_colors2, clim, ncolors, nlevels_iso, bottom, $
                     levels_final, levels_iso, c_colors_final,  $
                     cb_indx, cb_tickmark_name
;
contourplots_2triple, int2d_sc/max(int2d_sc), abs2d_sc/max(abs2d_sc), contr2d_sc/max(contr2d_sc), $
                      int2d_fvm/max(int2d_fvm), abs2d_sc/max(abs2d_fvm), contr2d_fvm/max(contr2d_fvm), $
                      x, x, x, x, x, x, z, z, z, z, z, z, $
                      ncolors, nlevels_iso, bottom, $
                      levels_final, levels_iso, c_colors, $
                      cb_indx, cb_tickmark_name, $
                      xlim=xlim, $
                      ylim=ylim, $
                      titleStr1=titleStr1, $
                      titleStr2=titleStr2, $
                      titleStr3=titleStr3, $
                      titleStr4=titleStr4, $
                      titleStr5=titleStr5, $
                      titleStr6=titleStr6, $
                      xtitleStr=xtitleStr, $
                      ytitleStr=ytitleStr, $
                      ctitleStr=ctitleStr, $
                      ctable=13, $
                      /isotropic, /cb_top
!p.multi=0
;
;***************intensities with different contributions****************
;
window, windx, xsize=960, ysize=1200
device, decomposed=0
windx=windx+1
loadct, 0
;
;-----------------------intensities-------------------------------------
;
ymin=min([min(int1d), min(int1d_scray), min(int1d_fvmray), min(int1d_sc), min(int1d_fvm)])
ymax=max([max(int1d), max(int1d_scray), max(int1d_fvmray), max(int1d_sc), max(int1d_fvm)])
ymin=ymin-0.1*(ymax-ymin)
ymax=ymax+0.1*(ymax-ymin)
;
!p.multi=[3,1,3]
loadct_lev, ci_black, ci_blue, ci_red, ci_cyan, ci_green, ci_magenta, ci_yellow, ci_white
plot, r, int1d, $
      yrange=[ymin,ymax], $
      xtitle='r', $
      ytitle='intensity', $
      title='intensity', $
      charsize=2.
;
oplot, r, int1d_scray, $
      color=ci_blue
oplot, r, int1d_fvmray, $
      color=ci_cyan
oplot, r, int1d_sc, $
       line=2, $
       color=ci_red
oplot, r, int1d_fvm, $
       line=2, $
       color=ci_green
;
lstring1='exact'
lstring2='ray (sc2d)'
lstring3='ray (fvm2d)'
lstring4='sc1d'
lstring5='fvm1d'
lstring6=''
lstring7=''

legend, [lstring1, lstring2, lstring3, lstring4, lstring5], $
        /right_legend, $
        /bottom, $
        psym=[0,0,0,0,0], $
        linestyle=[0,0,0,2,2], $
        color= [255,ci_blue,ci_cyan,ci_red,ci_green], $
        textcolors=255
;
;-------------------absorption terms in each cell-----------------------
;
ymin=min([min(abs1d), min(abs1d_sc), min(abs1d_fvm)])
ymax=max([max(abs1d), max(abs1d_sc), max(abs1d_fvm)])
ymin=ymin-0.1*(ymax-ymin)
ymax=ymax+0.1*(ymax-ymin)
;
;
!p.multi=[2,1,3]
plot, r, abs1d, $
      yrange=[ymin,ymax], $
      title='absorbing terms (2d skipped since different dels)', $
      xtitle='r', $
      ytitle='absorbing terms', $
      charsize=2.
oplot, r, abs1d_sc, $
       line=2, $
       color=ci_red
oplot, r, abs1d_fvm, $
       line=2, $
       color=ci_green
;
lstring1='exact'
lstring2='sc1d'
lstring3='fvm1d'
;
legend, [lstring1, lstring2, lstring3], $
        /right_legend, $
        psym=[0,0,0], $
        linestyle=[0,2,2], $
        color= [255,ci_red,ci_green], $
        textcolors=255
;
;-----------------------source contribution terms-----------------------
;
ymin=min([min(contr1d), min(contr1d_sc), min(contr1d_fvm)])
ymax=max([max(contr1d), max(contr1d_sc), max(contr1d_fvm)])
ymin=ymin-0.1*(ymax-ymin)
ymax=ymax+0.1*(ymax-ymin)
;
!p.multi=[1,1,3]
plot, r, contr1d, $
      yrange=[ymin,ymax], $
      title='source contribution terms (2d skipped since different dels)', $
      xtitle='r', $
      ytitle='source contribution terms', $
      charsize=2.
loadct_lev, ci_black, ci_blue, ci_red, ci_cyan, ci_green, ci_magenta, ci_yellow, ci_white
oplot, r, contr1d_sc, $
       line=2, $
       color=ci_red
oplot, r, contr1d_fvm, $
       line=2, $
       color=ci_green
;
lstring1='exact'
lstring2='sc1d'
lstring3='fvm1d'
;
legend, [lstring1, lstring2, lstring3], $
        /right_legend, $
        psym=[0,0,0], $
        linestyle=[0,2,2], $
        color= [255,ci_red,ci_green], $
        textcolors=255
;
!p.multi=0
loadct, 0
;
;*************************output to ps file*****************************
;
if (not keyword_set(oname)) then return
;
;-----------------------contour plot of intensity-----------------------
;
oname='ps_files/benchmark01_a.ps'
print, 'writing output to: ', oname
;
aspect_ratio=1.6
set_plot,'ps'
xsize=17.780d0
ysize=xsize/aspect_ratio
device,file=oname, $
       decomposed=0, $
       color=1, $
       BITS_PER_PIXEL=8, $
       ysize=ysize, $
       xsize=xsize
;device, file=oname, decomposed=0, color=1, bits_per_pixel=8
;
clim=[0.,1.]
xlim=[-1.d0,min([5.d0,max(x)])]
ylim=[-1.d0,min([5.d0,max(z)])]
titleStr1='intensity'
xtitleStr='x'
ytitleStr='z'
;
loadct, 0
get_contour_colors2, clim, ncolors, nlevels_iso, bottom, $
                     levels_final, levels_iso, c_colors_final,  $
                     cb_indx, cb_tickmark_name
;
contourplots_single, int2d_sc/max(int2d_sc), $
                      x, z, $
                      ncolors, nlevels_iso, bottom, $
                      levels_final, levels_iso, c_colors, $
                      cb_indx, cb_tickmark_name, $
                      xlim=xlim, $
                      ylim=ylim, $
                      titleStr=titleStr1, $
                      xtitleStr=xtitleStr, $
                      ytitleStr=ytitleStr, $
                      ctitleStr=ctitleStr, $
                      ctable=13, $
                      /isotropic
if keyword_set(oname) then begin
   device, /close
   set_plot,'x'
endif
;
;---------------------------plots along ray-----------------------------
;
oname='ps_files/benchmark01_b.ps'
print, 'writing output to: ', oname
;
aspect_ratio=0.707
set_plot,'ps'
xsize=17.780d0
ysize=xsize/aspect_ratio
device,file=oname, $
       decomposed=0, $
       color=1, $
       BITS_PER_PIXEL=8, $
       ysize=ysize, $
       xsize=xsize
;device, file=oname, decomposed=0, color=1, bits_per_pixel=8
;
;intensities
!p.multi=[3,1,3]
;
ymin=min([min(int1d), min(int1d_scray), min(int1d_fvmray), min(int1d_sc), min(int1d_fvm)])
ymax=max([max(int1d), max(int1d_scray), min(int1d_fvmray), max(int1d_sc), min(int1d_fvm)])
ymin=ymin-0.1*(ymax-ymin)
ymax=ymax+0.1*(ymax-ymin)
;
loadct, 0
plot, r, int1d, $
      yrange=[ymin,ymax], $
      xtitle='r', $
      ytitle='intensity', $
      title='intensity', $
      charsize=1.2

;
oplot, r, int1d_scray, $
      color=ci_blue
oplot, r, int1d_fvmray, $
      color=ci_cyan
oplot, r, int1d_sc, $
       line=2, $
       color=ci_red
oplot, r, int1d_fvm, $
       line=2, $
       color=ci_green

loadct_lev, ci_black, ci_blue, ci_red, ci_cyan, ci_green, ci_magenta, ci_yellow, ci_white

oplot, r, int1d_scray, $
      color=ci_blue
oplot, r, int1d_scray, $
      psym=1, $
      color=ci_blue
oplot, r, int1d_fvmray, $
      color=ci_cyan
oplot, r, int1d_fvmray, $
      psym=1, $
      color=ci_cyan
oplot, r, int1d_sc, $
       line=2, $
       color=ci_red
oplot, r, int1d_fvm, $
       line=2, $
       color=ci_green
;
lstring1='exact'
lstring2='ray (sc2d)'
lstring3='ray (fvm2d)'
lstring4='sc1d'
lstring5='fvm1d'
;
legend, [lstring1, lstring2, lstring3, lstring4, lstring5], $
        /right_legend, $
        psym=[0,1,1,0,0], $
        linestyle=[0,1,1,2,2], $
        color= [0,ci_blue,ci_cyan,ci_red,ci_green]
;
;
;opacities
!p.multi=[2,1,3]
ymax=1.d0
ymin=min([min(opac2d),min(opac_u2d),min(opac_d2d),min(opac1d)])
ylim=[ymin-0.1*(ymax-ymin), ymax+0.1*(ymax-ymin)]
xlim=[0.9d0, sqrt(x(ndxmax-1)^2+z(ndzmax-1)^2)]
;
loadct, 0
plot, r, opac1d, $
      xrange=xlim, /xs, $
      yrange=ylim, $
      xtitle='r', $
      ytitle=textoidl('\chi'), $
      line=2, $
      charsize=1.2
;
loadct_lev, ci_black, ci_blue, ci_red, ci_cyan, ci_green, ci_magenta, ci_yellow, ci_white
;
for i=0, ndxmax-1 do begin
   for k=0, ndzmax-1 do begin
      rp=sqrt(x(i)^2 + z(k)^2)
      oplot, [rp,rp], [opac2d(i,k),opac2d(i,k)], psym=3, color=ci_blue

      ru=sqrt(x_u2d(i,k)^2 + z_u2d(i,k)^2)
      oplot, [ru,ru], [opac_u2d(i,k), opac_u2d(i,k)], psym=3, color=ci_red

      rd=sqrt(x_d2d(i,k)^2 + z_d2d(i,k)^2)
      oplot, [rd,rd], [opac_d2d(i,k), opac_d2d(i,k)], psym=3, color=ci_green
   endfor
endfor
legstr0=textoidl('opacity 1d')
legstr1=textoidl('opacity at r_p')
legstr2=textoidl('opacity at r_u')
legstr3=textoidl('opacity at r_d')
legend, [legstr0, legstr1, legstr2, legstr3], $
         linestyle=[2,0,0,0], $
         psym=[0,3,3,3], $
         color=[0,ci_blue,ci_red,ci_green], $
         /right_legend
;
;
;source functions
ymax=1.d0
ymin=min([min(scont2d),min(scont_u2d),min(scont_d2d),min(scont1d)])
ylim=[ymin-0.1*(ymax-ymin), ymax+0.1*(ymax-ymin)]
;
!p.multi=[1,1,3]
loadct, 0
plot, r, scont1d, $
      xrange=xlim, /xs, $
      yrange=ylim, $
      xtitle='r', $
      ytitle=textoidl('source function'), $
      line=2
;
loadct_lev, ci_black, ci_blue, ci_red, ci_cyan, ci_green, ci_magenta, ci_yellow, ci_white
;
for i=0, ndxmax-1 do begin
   for k=0, ndzmax-1 do begin
      rp=sqrt(x(i)^2 + z(k)^2)
      oplot, [rp,rp], [scont2d(i,k),scont2d(i,k)], psym=3, color=ci_blue

      ru=sqrt(x_u2d(i,k)^2 + z_u2d(i,k)^2)
      oplot, [ru,ru], [scont_u2d(i,k), scont_u2d(i,k)], psym=3, color=ci_red

      rd=sqrt(x_d2d(i,k)^2 + z_d2d(i,k)^2)
      oplot, [rd,rd], [scont_d2d(i,k), scont_d2d(i,k)], psym=3, color=ci_green
   endfor
endfor
legstr0=textoidl('source 1d')
legstr1=textoidl('source at r_p')
legstr2=textoidl('source at r_u')
legstr3=textoidl('source at r_d')
legend, [legstr0, legstr1, legstr2, legstr3], $
         psym=[0,3,3,3], $
         color=[0,ci_blue,ci_red,ci_green], $
         /bottom
;
if keyword_set(oname) then begin
   device, /close
   set_plot,'x'
endif


!p.multi=0
loadct, 0
;
end
;
;-----------------------------------------------------------------------
;-----------------------------------------------------------------------
;-----------------------------------------------------------------------
;
pro read_model2d, ndxmax, ndzmax, x, z, x_u2d, z_u2d, x_d2d, z_d2d, int2d_sc, abs2d_sc, contr2d_sc, $
                  int2d_fvm, abs2d_fvm, contr2d_fvm, xu_opac2d, zu_opac2d, $
                  xd_opac2d, zd_opac2d, xu_int2d, zu_int2d, scont2d, scont_u2d, scont_d2d, opac2d, opac_u2d, opac_d2d
;
fname='benchmark01.h5'
;
file_id = h5f_open(fname)
;
   group_id=h5g_open(file_id, 'dimensions')
      dset_id=h5a_open_name(group_id, 'ndxmax')
         ndxmax=h5a_read(dset_id)
         ndxmax=ndxmax(0)
      h5a_close, dset_id
      dset_id=h5a_open_name(group_id, 'ndzmax')
         ndzmax=h5a_read(dset_id)
         ndzmax=ndzmax(0)
      h5a_close, dset_id
   h5g_close, group_id
;
   x=fltarr(ndxmax)
   z=fltarr(ndzmax)
;
   group_id=h5g_open(file_id, 'coordinates')
      dset_id=h5d_open(group_id, 'x')
         x=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'z')
         z=h5d_read(dset_id)
      h5d_close, dset_id
   h5g_close, group_id
;
   int2d_sc=fltarr(ndxmax,ndzmax)
   abs2d_sc=fltarr(ndxmax,ndzmax)
   contr2d_sc=fltarr(ndxmax,ndzmax)
   int2d_fvm=fltarr(ndxmax,ndzmax)
   abs2d_fvm=fltarr(ndxmax,ndzmax)
   contr2d_fvm=fltarr(ndxmax,ndzmax)

   group_id=h5g_open(file_id, 'solution2d')
      dset_id=h5d_open(group_id, 'int2d_sc')
         int2d_sc=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'abs2d_sc')
         abs2d_sc=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'contr2d_sc')
         contr2d_sc=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'int2d_fvm')
         int2d_fvm=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'abs2d_fvm')
         abs2d_fvm=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'contr2d_fvm')
         contr2d_fvm=h5d_read(dset_id)
      h5d_close, dset_id
   h5g_close, group_id
;
   scont2d=fltarr(ndxmax,ndzmax)
   opac2d=fltarr(ndxmax,ndzmax)

   group_id=h5g_open(file_id, 'model2d')
      dset_id=h5d_open(group_id, 'scont2d')
         scont2d=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'opac2d')
         opac2d=h5d_read(dset_id)
      h5d_close, dset_id
   h5g_close, group_id
;
   x_u2d=fltarr(ndxmax,ndzmax)
   x_d2d=fltarr(ndxmax,ndzmax)
   z_u2d=fltarr(ndxmax,ndzmax)
   z_d2d=fltarr(ndxmax,ndzmax)
   xu_int2d=fltarr(ndxmax,ndzmax,8)
   zu_int2d=fltarr(ndxmax,ndzmax,8)
   xu_opac2d=fltarr(ndxmax,ndzmax,8)
   zu_opac2d=fltarr(ndxmax,ndzmax,8)
   xd_opac2d=fltarr(ndxmax,ndzmax,8)
   zd_opac2d=fltarr(ndxmax,ndzmax,8)
   scont_u2d=fltarr(ndxmax,ndzmax)
   opac_u2d=fltarr(ndxmax,ndzmax)
   scont_d2d=fltarr(ndxmax,ndzmax)
   opac_d2d=fltarr(ndxmax,ndzmax)

   group_id=h5g_open(file_id, 'debug')
      dset_id=h5d_open(group_id, 'x_upwind')
         x_u2d=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'z_upwind')
         z_u2d=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'x_downwind')
         x_d2d=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'z_downwind')
         z_d2d=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'xu_opac2d')
         xu_opac2d=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'zu_opac2d')
         zu_opac2d=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'xd_opac2d')
         xd_opac2d=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'zd_opac2d')
         zd_opac2d=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'xu_int2d')
         xu_int2d=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'zu_int2d')
         zu_int2d=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'scont_u2d')
         scont_u2d=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'scont_d2d')
         scont_d2d=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'opac_u2d')
         opac_u2d=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'opac_d2d')
         opac_d2d=h5d_read(dset_id)
      h5d_close, dset_id

   h5g_close, group_id
;
h5f_close, file_id
;
end
;
;-----------------------------------------------------------------------
;-----------------------------------------------------------------------
;-----------------------------------------------------------------------
;
pro read_model1d, nr, r, int1d, abs1d, contr1d, int1d_sc, abs1d_sc, contr1d_sc, $
                  int1d_fvm, abs1d_fvm, contr1d_fvm, int1d_scray, int1d_fvmray, opac1d, scont1d
;
fname='benchmark01.h5'
;
file_id = h5f_open(fname)
;
   group_id=h5g_open(file_id, 'dimensions')
      dset_id=h5a_open_name(group_id, 'nr')
         nr=h5a_read(dset_id)
         nr=nr(0)
      h5a_close, dset_id
   h5g_close, group_id
;
   r=fltarr(nr)
;
   group_id=h5g_open(file_id, 'coordinates')
      dset_id=h5d_open(group_id, 'r')
         r=h5d_read(dset_id)
      h5d_close, dset_id
   h5g_close, group_id
;
   int1d_sc=fltarr(nr)
   int1d=fltarr(nr)
   int1d_scray=fltarr(nr)
   int1d_fvmray=fltarr(nr)
   abs1d_sc=fltarr(nr)
   abs1d=fltarr(nr)
   contr1d_sc=fltarr(nr)
   contr1d=fltarr(nr)
   opac1d=fltarr(nr)
   scont1d=fltarr(nr)

   group_id=h5g_open(file_id, 'solution1d')
      dset_id=h5d_open(group_id, 'int1d_sc')
         int1d_sc=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'abs1d_sc')
         abs1d_sc=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'contr1d_sc')
         contr1d_sc=h5d_read(dset_id)
      h5d_close, dset_id

      dset_id=h5d_open(group_id, 'int1d_fvm')
         int1d_fvm=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'abs1d_fvm')
         abs1d_fvm=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'contr1d_fvm')
         contr1d_fvm=h5d_read(dset_id)
      h5d_close, dset_id

      dset_id=h5d_open(group_id, 'int1d_scray')
         int1d_scray=h5d_read(dset_id)
      h5d_close, dset_id

      dset_id=h5d_open(group_id, 'int1d_fvmray')
         int1d_fvmray=h5d_read(dset_id)
      h5d_close, dset_id

      dset_id=h5d_open(group_id, 'int1d')
         int1d=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'abs1d')
         abs1d=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'contr1d')
         contr1d=h5d_read(dset_id)
      h5d_close, dset_id

      dset_id=h5d_open(group_id, 'opac1d')
         opac1d=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'scont1d')
         scont1d=h5d_read(dset_id)
      h5d_close, dset_id
   h5g_close, group_id
;
h5f_close, file_id
;
end

;
;-----------------------------------------------------------------------
;-----------------------------------------------------------------------
;-----------------------------------------------------------------------
;
pro plot_grid, x, z, nx, nz, xlim=xlim, ylim=ylim
;
;plot grid, additional to a circle
if(not keyword_set(xlim)) then xlim=[min(x),max(x)]
if(not keyword_set(ylim)) then ylim=[min(z),max(z)]
;
plot, [x(0),x(1)], [z(0),z(0)], $
      xrange=xlim, $
      yrange=ylim, $
      /xs, /ys, /isotropic, xtitle='x', ytitle='z'

for i=0, nz-1 do begin
   oplot, [xlim(0),xlim(1)], [z(i),z(i)]
endfor
;
for i=0, nx-1 do begin
   oplot, [x(i),x(i)], [ylim(0),ylim(1)]
endfor
;
;oplot circle
nt=300
theta=findgen(nt)*2.d0*!pi/float(nt-1)
oplot, sin(theta), cos(theta)
;
end
