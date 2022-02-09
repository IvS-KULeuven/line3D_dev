pro benchmark05, dir=DIR, windx=WINDX, oname=ONAME, xlim=XLIM, ylim=YLIM
;
;plots contours of searchlight along a direction
;
;------------------read all information from hdf5-file------------------
;
if(not keyword_set(dir)) then dir='.'
fname=dir+'/benchmark05.h5'
;
file_id = h5f_open(fname)
;
   group_id = h5g_open(file_id, 'input_parameters')
      att_id=h5a_open_name(group_id, 'kcont')
         kcont=h5a_read(att_id)
         kcont=kcont(0)
      h5a_close, att_id
      att_id=h5a_open_name(group_id, 'eps_cont')
         eps_cont=h5a_read(att_id)
         eps_cont=eps_cont(0)
      h5a_close, att_id
      att_id=h5a_open_name(group_id, 'xic1')
         xic1=h5a_read(att_id)
         xic1=xic1(0)
      h5a_close, att_id
   h5g_close, group_id

   group_id = h5g_open(file_id, 'dimensions')
      att_id=h5a_open_name(group_id, 'nr')
         nr=h5a_read(att_id)
         nr=nr(0)
      h5a_close, att_id
   h5g_close, group_id
;
   group_id = h5g_open(file_id, 'coordinates')
      dset_id=h5d_open(group_id, 'r')
         r=h5d_read(dset_id)
      h5d_close, dset_id
   h5g_close, group_id
;
   group_id = h5g_open(file_id, 'convergence_behaviour')
      att_id=h5a_open_name(group_id, 'itmaxc')
         itmaxc=h5a_read(att_id)
         itmaxc=itmaxc(0)
      h5a_close, att_id
      dset_id=h5d_open(group_id, 'epsmaxc_sc')
         epsmaxc_sc=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'epsmaxc_fvm')
         epsmaxc_fvm=h5d_read(dset_id)
      h5d_close, dset_id
   h5g_close, group_id
;
   group_id = h5g_open(file_id, 'model_cr')
      dset_id=h5d_open(group_id, 't1d')
         t1d=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 't1d_jo')
         t1d_jo=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'opac1d')
         opac1d=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'opac1d_jo')
         opac1d_jo=h5d_read(dset_id)
      h5d_close, dset_id
   h5g_close, group_id
;
   group_id = h5g_open(file_id, 'solution_cr')
      dset_id=h5d_open(group_id, 'mint_sc')
         mint_sc=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'mint_fvm')
         mint_fvm=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'mint_joray')
         mint_joray=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'mint_jomom')
         mint_jomom=h5d_read(dset_id)
      h5d_close, dset_id
   h5g_close, group_id
;
h5f_close, file_id
;
;trad=44.5d3
;xnue0=2.191983810829996D15
;xic2=bnue(xnue0,1.*trad)
;print, xic2, xic1, xic2/xic1
;xic1=xic2
;stop
;normalize everything to xic1
mint_sc=mint_sc/xic1
mint_fvm=mint_fvm/xic1
mint_joray=mint_joray/xic1
mint_jomom=mint_jomom/xic1
;
;-----------------------------------------------------------------------
;-----------------------define range------------------------------------
;
if(not keyword_set(xlim)) then begin
   xmin=min(r)
   xmax=max(r)
   xlim=[xmin, xmax]
endif
;
if(not keyword_set(ylim)) then begin
   ymin=min([min(mint_sc*r^2), min(mint_fvm*r^2), min(mint_joray*r^2), min(mint_jomom*r^2)])
   ymax=max([max(mint_sc*r^2), max(mint_fvm*r^2), max(mint_joray*r^2), min(mint_jomom*r^2)])
   ylim=[ymin, ymax]
endif
;
;----------------------------title-strings------------------------------
;
;titlestr=textoidl('i (\theta=') + string(acos(mu), format='(f9.5)') + $
;         textoidl(', \phi=') + string(phi, format='(f9.5)') + ')'
titlestr=''
xtitlestr=textoidl('r')
ytitlestr=textoidl('J*r^2')
;
;-----------------------------------------------------------------------
;
if(not keyword_set(windx)) then windx=0
;
if(keyword_set(oname)) then begin
   oname='ps_files/benchmark05_a.ps'
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
loadct_lev, ci_black, ci_blue, ci_red, ci_cyan, ci_green, ci_magenta, ci_yellow, ci_white
;
!p.multi=[6,2,3]
;
;xlim2=[1.,1.5]
;ylim2=[0.8,1.2]
;w=0.5*(1.-sqrt(1.-1./r^2))
plot_io, r, mint_joray*r^2, $
      xrange=xlim, $
      yrange=ylim, $
      title=titleStr, $
      xtitle=xtitleStr, $
      ytitle=ytitleStr, $
         charsize=2., /xs, /ys

oplot, r, mint_jomom*r^2, $
      line=2

oplot, r, mint_sc*r^2, $
      color=ci_blue, $
      line=0

oplot, r, mint_sc*r^2, psym=1
;oplot, r, mint_sc*6.*r^2, color=ci_green

oplot, r, mint_fvm*r^2, $
      color=ci_red, $
      line=0
;
lstr1='1d (ray by ray)'
lstr2='1d (moments)'
lstr3='2d SC'
lstr4='2d FVM'
if(keyword_set(oname)) then begin
   legend, [lstr1, lstr2, lstr3, lstr4], $
           psym=[0,0,0,0], $
           linestyle=[0,2,0,0], $
           color=[ci_black, ci_black,ci_blue,ci_red], $
           textcolor=ci_black, $
           /right_legend
endif else begin
   legend, [lstr1, lstr2, lstr3, lstr4], $
           psym=[0,0,0,0], $
           linestyle=[0,2,0,0], $
           color=[ci_white,ci_white,ci_blue,ci_red], $
           textcolor=ci_white, $
           /right_legend
endelse
;

!p.multi=[5,2,3]
;
plot, 1.d0-1.d0/r, mint_joray*r^2, $
      xrange=[0.d0,1.d0], $
      yrange=ylim, $
      title=titleStr, $
      xtitle=textoidl('1-1/r'), $
      ytitle=ytitleStr, $
      charsize=2., /ylog

oplot, 1.d0-1.d0/r, mint_jomom*r^2, $
      line=2

oplot, 1.d0-1.d0/r, mint_sc*r^2, $
      color=ci_blue, $
      line=0

oplot, 1.d0-1.d0/r, mint_fvm*r^2, $
      color=ci_red, $
      line=0
;
lstr1='1d (ray by ray)'
lstr2='1d (moments)'
lstr3='2d SC'
lstr4='2d FVM'
if(keyword_set(oname)) then begin
   legend, [lstr1, lstr2, lstr3, lstr4], $
           psym=[0,0,0,0], $
           linestyle=[0,2,0,0], $
           color=[ci_black, ci_black,ci_blue,ci_red], $
           textcolor=ci_black, $
           /right_legend
endif else begin
   legend, [lstr1, lstr2, lstr3, lstr4], $
           psym=[0,0,0,0], $
           linestyle=[0,2,0,0], $
           color=[ci_white,ci_white,ci_blue,ci_red], $
           textcolor=ci_white, $
           /right_legend
endelse
;
;plot errors
ymin=min([min(mint_sc/mint_joray), min(mint_fvm/mint_joray)])
ymax=max([max(mint_sc/mint_joray), max(mint_fvm/mint_joray)])
ylim=[ymin, ymax]
xlim2=[0.5,2.5]
ylim2=[0.8,1.2]
!p.multi=[4,2,3]
plot, [0.,0.], [0.,0.], $
      xrange=xlim2, $
      yrange=ylim2, $
      title=titleStr, $
      xtitle=xtitleStr, $
      ytitle='J(2d)/J(1d)', $
      charsize=2.

oplot, r, mint_sc/mint_joray, $
      color=ci_blue, $
      line=0
;oplot, r, mint_sc*6./mint_joray, color=ci_green

oplot, r, mint_sc/mint_joray, psym=1
oplot, r, mint_fvm/mint_joray, $
      color=ci_red, $
      line=0

;print, mint_sc*6.*r^2
;stop
;
lstr1='2d SC'
lstr2='2d FVM'
if(keyword_set(oname)) then begin
   legend, [lstr1, lstr2], $
           psym=[0,0], $
           linestyle=[0,0], $
           color=[ci_blue,ci_red], $
           textcolor=ci_black, $
           /right_legend
endif else begin
   legend, [lstr1, lstr2], $
           psym=[0,0], $
           linestyle=[0,0], $
           color=[ci_blue,ci_red], $
           textcolor=ci_white, $
           /right_legend
endelse
;
!p.multi=[3,2,3]
plot, [0.,0.], [0.,0.], $
      xrange=[0.d0,1.d0], $
      yrange=ylim, $
      title=titleStr, $
      xtitle=textoidl('1-1/r'), $
      ytitle='J(2d)/J(1d)', $
      charsize=2.

oplot, 1.d0-1.d0/r, mint_sc/mint_joray, $
      color=ci_blue, $
      line=0
oplot, 1.-1./r, mint_sc/mint_joray, psym=1
oplot, 1.d0-1.d0/r, mint_fvm/mint_joray, $
      color=ci_red, $
      line=0
;
lstr1='2d SC'
lstr2='2d FVM'
if(keyword_set(oname)) then begin
   legend, [lstr1, lstr2], $
           psym=[0,0], $
           linestyle=[0,0], $
           color=[ci_blue,ci_red], $
           textcolor=ci_black, $
           /right_legend
endif else begin
   legend, [lstr1, lstr2], $
           psym=[0,0], $
           linestyle=[0,0], $
           color=[ci_blue,ci_red], $
           textcolor=ci_white, $
           /right_legend
endelse
;
;only to check if temperature model is consistent
!p.multi=[2,2,3]
plot, r, t1d_jo, $
      xtitle='r', $
      ytitle=textoidl('T'), $
      title='consistency check', $
      charsize=2., $
      yrange=[min(t1d_jo)-0.1*(max(t1d_jo)-min(t1d_jo)),max(t1d_jo)]
oplot, r, t1d, $
       psym=1
lstr1='temperature (1d)'
lstr2='temperature (2d)'
if(keyword_set(oname)) then begin
   legend, [lstr1, lstr2], $
         psym=[0,1], $
         linestyle=[0,0], $
         color=[ci_black,ci_black], $
         textcolor=ci_black, $
         /right_legend
endif else begin
   legend, [lstr1, lstr2], $
         psym=[0,1], $
         linestyle=[0,0], $
         color=[ci_white,ci_white], $
         textcolor=ci_white, $
         /right_legend
endelse
;
;
;only to check if opacity model is consistent
!p.multi=[1,2,3]
plot, r, opac1d_jo, $
      xtitle='r', $
      ytitle=textoidl('\chi'), $
      title='consistency check', $
      /ylog, $
      charsize=2.
oplot, r, opac1d, $
       psym=1
lstr1='opacity (1d)'
lstr2='opacity (2d)'
if(keyword_set(oname)) then begin
   legend, [lstr1, lstr2], $
         psym=[0,1], $
         linestyle=[0,0], $
         color=[ci_black,ci_black], $
         textcolor=ci_black, $
         /right_legend
endif else begin
   legend, [lstr1, lstr2], $
         psym=[0,1], $
         linestyle=[0,0], $
         color=[ci_white,ci_white], $
         textcolor=ci_white, $
         /right_legend
endelse
;
if keyword_set(oname) then begin
   device, /close
   set_plot,'x'
endif
;
!p.multi=0
;
;---------------------------convergence behaviour-----------------------
;
if(keyword_set(oname)) then begin
   oname='benchmark05_b.ps'
   print, "writing output to: ", oname
   set_plot,'ps'
   device,file=oname, xsize=19., ysize=26.7, xoffset=1., yoffset=1. , $
    decomposed=0, color=1, bits_per_pixel=8
endif else begin
   window, windx
   device, decomposed=0
   windx=windx+1
endelse
;
epsmaxc_sc=abs(epsmaxc_sc)
epsmaxc_fvm=abs(epsmaxc_fvm)
iternr=indgen(itmaxc)
;
xmax1=max(where(epsmaxc_sc gt 0))
xmax2=max(where(epsmaxc_fvm gt 0))
;
ymin1=epsmaxc_sc(xmax1)
ymin2=epsmaxc_fvm(xmax2)
;
ymax1=max(epsmaxc_sc)
ymax2=max(epsmaxc_fvm)
;
xmin=0
;
xlim=[xmin,max([xmax1,xmax2])]
ylim=[min([ymin1,ymin2]),max([ymax1,ymax2])]
;
ytitleStr=textoidl('((J^{(k-1)} - J^{(k)}) /J^{(k)})_{max}')
xtitleStr='# iterations'
titleStr='CONTINUUM'
;
plot, iternr, epsmaxc_sc, $
      /ylog, $
      line=0, $
      yrange=ylim, $
      xrange=xlim, $
      title=titlestr, $
      xtitle=xtitlestr, $
      ytitle=ytitlestr
oplot, iternr, epsmaxc_fvm, $
      line=2
lstr1='2d short characteristics'
lstr2='2d finite volume'
;
if(keyword_set(oname)) then begin
   legend, [lstr1, lstr2], $
         psym=[0,0], $
         linestyle=[0,2], $
         color=[0,0], $
         textcolor=0, $
         /right_legend
endif else begin
   legend, [lstr1, lstr2], $
         psym=[0,0], $
         linestyle=[0,2], $
         color=[255,255], $
         textcolor=255, $
        /right_legend
endelse

if keyword_set(oname) then begin
   device, /close
   set_plot,'x'
endif
;
end
