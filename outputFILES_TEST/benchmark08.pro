pro benchmark08, dir=DIR, windx=WINDX, oname=ONAME, xlim=XLIM, ylim=YLIM
;
;plot line source functions for spherically symmetric models
;
;------------------read all information from hdf5-file------------------
;
if(not keyword_set(dir)) then dir='.'
fname=dir+'/benchmark08.h5'
;
file_id = h5f_open(fname)
;
   group_id = h5g_open(file_id, 'input_parameters')
      att_id=h5a_open_name(group_id, 'opt_opal')
         opt_opal=h5a_read(att_id)
         opt_opal=opt_opal(0)
      h5a_close, att_id
      att_id=h5a_open_name(group_id, 'kline')
         kline=h5a_read(att_id)
         kline=kline(0)
      h5a_close, att_id
      att_id=h5a_open_name(group_id, 'alpha')
         alpha=h5a_read(att_id)
         alpha=alpha(0)
      h5a_close, att_id
      att_id=h5a_open_name(group_id, 'kappa0')
         kappa0=h5a_read(att_id)
         kappa0=kappa0(0)
      h5a_close, att_id
      att_id=h5a_open_name(group_id, 'eps_line')
         eps_line=h5a_read(att_id)
         eps_line=eps_line(0)
      h5a_close, att_id
      att_id=h5a_open_name(group_id, 'xic1')
         xic1=h5a_read(att_id)
         xic1=xic1(0)
      h5a_close, att_id
      att_id=h5a_open_name(group_id, 'vmax')
         vmax=h5a_read(att_id)
         vmax=vmax(0)*1.d5
      h5a_close, att_id
      att_id=h5a_open_name(group_id, 'rstar')
         rstar=h5a_read(att_id)
         rstar=rstar(0)
         sr=rstar*!rsu
      h5a_close, att_id
   h5g_close, group_id
;
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
      att_id=h5a_open_name(group_id, 'itmaxl')
         itmaxl=h5a_read(att_id)
         itmaxl=itmaxl(0)
      h5a_close, att_id
      dset_id=h5d_open(group_id, 'epsmaxl_sc')
         epsmaxl_sc=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'epsmaxl_fvm')
         epsmaxl_fvm=h5d_read(dset_id)
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
      dset_id=h5d_open(group_id, 'opalbar1d')
         opalbar1d=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'opalbar1d_jo')
         opalbar1d_jo=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'velr1d')
         velr=h5d_read(dset_id)
      h5d_close, dset_id
   h5g_close, group_id
;
   group_id = h5g_open(file_id, 'solution_cr')
      dset_id=h5d_open(group_id, 'sline_sc')
         sline_sc=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'sline_fvm')
         sline_fvm=h5d_read(dset_id)
      h5d_close, dset_id
     dset_id=h5d_open(group_id, 'ssobo1d_cr')
        ssobo_cr=h5d_read(dset_id)
     h5d_close, dset_id
     dset_id=h5d_open(group_id, 'sline_jo')
         sline_jo=h5d_read(dset_id)
      h5d_close, dset_id
     dset_id=h5d_open(group_id, 'ssobo1d_jo')
         ssobo_jo=h5d_read(dset_id)
      h5d_close, dset_id
   h5g_close, group_id
;
h5f_close, file_id
;
;normalize everything to xic1
sline_sc=sline_sc/xic1
sline_fvm=sline_fvm/xic1
ssobo_cr=ssobo_cr/xic1
sline_jo=sline_jo/xic1
ssobo_jo=ssobo_jo/xic1
;
velr=velr/vmax
;
;---------------------------calculate errors----------------------------
;
errors1d, r(0:nr-2), sline_jo(0:nr-2), r(0:nr-2), sline_sc(0:nr-2), erri, errm, err_max, devm
;
;-----------------------define range------------------------------------
;
if(not keyword_set(xlim)) then begin
   xmin=min(r)
   xmax=max(r)
   xlim=[xmin, xmax]
endif
;
if(not keyword_set(ylim)) then begin
   ymin=min([min(sline_sc*r^2), min(sline_fvm*r^2), min(sline_jo*r^2), min(ssobo_jo*r^2), min(ssobo_cr*r^2)])
   ymax=max([max(sline_sc*r^2), max(sline_fvm*r^2), max(sline_jo*r^2), min(ssobo_jo*r^2), max(ssobo_cr*r^2)])
   ylim=[ymin, ymax]
endif
;
;----------------------------title-strings------------------------------
;
if(opt_opal eq 0) then begin
   titlestr=textoidl('k_L=') + string(kline, format='(f9.5)')
endif else begin
   titlestr=textoidl('\kappa_0, \alpha=') + string(kappa0, format='(f9.5)') + $
            ', ' + string(alpha, format='(f9.5)')
endelse
;
xtitlestr=textoidl('r')
ytitlestr=textoidl('S_L*r^2/I_c')
;
;-----------------------------------------------------------------------
;
if(not keyword_set(windx)) then windx=0
;
if(keyword_set(oname)) then begin
   oname='ps_files/benchmark08_a.ps'
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
plot, r, sline_jo*r^2, $
      xrange=xlim, $
      yrange=ylim, $
      title=titleStr, $
      xtitle=xtitleStr, $
      ytitle=ytitleStr, $
      charsize=2.

oplot, r, ssobo_jo*r^2, $
      line=2

oplot, r, sline_sc*r^2, $
      color=ci_blue, $
      line=0
oplot, r, sline_fvm*r^2, $
      color=ci_red, $
      line=0
oplot, r, ssobo_cr*r^2, $
      color=ci_green, $
      line=0
;
lstr1='1d solution'
lstr2='1d (sobolev)'
lstr3='2d SC'
lstr4='2d FVM'
lstr5='3d (sobolev)'
;
if(keyword_set(oname)) then begin
   legend, [lstr1, lstr2, lstr3, lstr4, lstr5], $
           psym=[0,0,0,0,0], $
           linestyle=[0,2,0,0,0], $
           color=[ci_black, ci_black,ci_blue,ci_red,ci_green], $
           textcolor=ci_black, $
           /right_legend
endif else begin
   legend, [lstr1, lstr2, lstr3, lstr4, lstr5], $
           psym=[0,0,0,0,0], $
           linestyle=[0,2,0,0,0], $
           color=[ci_white,ci_white,ci_blue,ci_red,ci_green], $
           textcolor=ci_white, $
           /right_legend
endelse
;
!p.multi=[5,2,3]
;
plot, velr, sline_jo*r^2, $
      xrange=[0.,1.], $
      yrange=[0.,1.], $
      title=titleStr, $
      xtitle=textoidl('v'), $
      ytitle=ytitleStr, $
      charsize=2.

oplot, velr, ssobo_jo*r^2, $
      line=2

oplot, velr, sline_sc*r^2, $
      color=ci_blue, $
      line=0
oplot, velr, sline_fvm*r^2, $
      color=ci_red, $
      line=0
oplot, velr, ssobo_cr*r^2, $
      color=ci_green, $
      line=0
;
lstr1='1d solution'
lstr2='1d (sobolev)'
lstr3='2d SC'
lstr4='2d FVM'
lstr5='3d (sobolev)'
;
if(keyword_set(oname)) then begin
   legend, [lstr1, lstr2, lstr3, lstr4, lstr5], $
           psym=[0,0,0,0,0], $
           linestyle=[0,2,0,0,0], $
           color=[ci_black, ci_black,ci_blue,ci_red,ci_green], $
           textcolor=ci_black, $
           /right_legend
endif else begin
   legend, [lstr1, lstr2, lstr3, lstr4, lstr5], $
           psym=[0,0,0,0,0], $
           linestyle=[0,2,0,0,0], $
           color=[ci_white,ci_white,ci_blue,ci_red,ci_green], $
           textcolor=ci_white, $
           /right_legend
endelse
;
;plot error in r-space
ymin=min([min(sline_sc/sline_jo), min(sline_fvm/sline_jo)])
ymax=max([max(sline_sc/sline_jo), max(sline_fvm/sline_jo)])
ylim=[ymin, ymax]
;
!p.multi=[4,2,3]
plot, [0.,0.], [0.,0.], $
      xrange=xlim, $
      yrange=ylim, $
      title=titleStr, $
      xtitle=xtitleStr, $
      ytitle=textoidl('S_L(2d)/S_L(1d)'), $
      charsize=2.
oplot, r, sline_sc/sline_jo, $
      color=ci_blue, $
      line=0
oplot, r, sline_fvm/sline_jo, $
      color=ci_red, $
      line=0
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
;plot error in velocity-space
!p.multi=[3,2,3]
plot, [0.,0.], [0.,0.], $
      xrange=[0.,1.], $
      yrange=ylim, $
      title=titleStr, $
      xtitle=textoidl('v'), $
      ytitle=textoidl('S_L(2d)/S_L(1d)'), $
      charsize=2.
oplot, velr, sline_sc/sline_jo, $
      color=ci_blue, $
      line=0
oplot, velr, sline_fvm/sline_jo, $
      color=ci_red, $
      line=0
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
!p.multi=[6,3,6]
plot, r, t1d_jo, $
      xtitle='r', $
      ytitle=textoidl('T'), $
      title='consistency check', $
      charsize=2., $
      yrange=[min(t1d_jo)-0.1*(max(t1d_jo)-min(t1d_jo)),max(t1d_jo)]
oplot, r, t1d, $
       psym=1
lstr1='1d'
lstr2='2d'
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
;only to check if opacity model is consistent
!p.multi=[5,3,6]
plot, r, opalbar1d_jo, $
      xtitle='r', $
      ytitle=textoidl('\chi bar'), $
      title='consistency check', $
      /ylog, $
      charsize=2.
oplot, r, opalbar1d, $
       psym=1
lstr1='1d'
lstr2='2d'
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
;only to check if velocity model is consistent
!p.multi=[4,3,6]
plot, r, velr, $
      xtitle='r', $
      ytitle=textoidl('v_r / v_{max}'), $
      charsize=2.
oplot, r, velr, psym=1
;
;plot (frequency integrated) tau
tau=fltarr(nr)
tau(nr-1)=0.d0
dtau=fltarr(nr-1)
for i=nr-2, 0, -1 do begin
   dtau(i) = 0.5d0*(opalbar1d(i)+opalbar1d(i+1))*(r(i+1)-r(i))*sr
   tau(i) = tau(i+1)+dtau(i)
endfor
!p.multi=[3,3,6]
plot, r, tau, $
      xtitle=textoidl('r'), $
      ytitle=textoidl('\tau_r bar'), $
      charsize=2.
oplot, r, tau, psym=1
;
!p.multi=[2,3,6]
plot, r(0:nr-2), dtau, $
      xtitle=textoidl('r'), $
      ytitle=textoidl('\Delta tau_r bar'), $
      charsize=2.
oplot, r(0:nr-2), dtau, psym=1
;


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
   oname='ps_files/benchmark08_b.ps'
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
epsmaxl_sc=abs(epsmaxl_sc)
epsmaxl_fvm=abs(epsmaxl_fvm)
iternr=indgen(itmaxl)
;
xmax1=max(where(epsmaxl_sc gt 0))
xmax2=max(where(epsmaxl_fvm gt 0))
;
ymin1=epsmaxl_sc(xmax1)
ymin2=epsmaxl_fvm(xmax2)
;
ymax1=max(epsmaxl_sc)
ymax2=max(epsmaxl_fvm)
;
xmin=0
;
xlim=[xmin,max([xmax1,xmax2])]
ylim=[min([ymin1,ymin2]),max([ymax1,ymax2])]
;
ytitleStr=textoidl('((Jbar^{(k-1)} - Jbar^{(k)}) /Jbar^{(k)})_{max}')
xtitleStr='# iterations'
titleStr='LINE'
;
plot, iternr, epsmaxl_sc, $
      /ylog, $
      line=0, $
      yrange=ylim, $
      xrange=xlim, $
      title=titlestr, $
      xtitle=xtitlestr, $
      ytitle=ytitlestr
oplot, iternr, epsmaxl_fvm, $
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
