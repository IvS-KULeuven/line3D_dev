pro benchmark07, windx=WINDX, oname=ONAME, xlim=XLIM, ylim=YLIM
;
;plot line source functions for spherically symmetric models
;
;------------------read all information from hdf5-file------------------
;
fname='benchmark07.h5'
;
file_id = h5f_open(fname)
;
   group_id = h5g_open(file_id, 'input_parameters')
      att_id=h5a_open_name(group_id, 'kline')
         kline=h5a_read(att_id)
         kline=kline(0)
      h5a_close, att_id
      att_id=h5a_open_name(group_id, 'eps_line')
         eps_line=h5a_read(att_id)
         eps_line=eps_line(0)
      h5a_close, att_id
      att_id=h5a_open_name(group_id, 'xic1')
         xic1=h5a_read(att_id)
         xic1=xic1(0)
      h5a_close, att_id
   h5g_close, group_id
;
   group_id = h5g_open(file_id, 'dimensions')
      att_id=h5a_open_name(group_id, 'nr')
         nr=h5a_read(att_id)
         nr=nr(0)
      h5a_close, att_id
      att_id=h5a_open_name(group_id, 'ndxmax')
         ndxmax=h5a_read(att_id)
         ndxmax=ndxmax(0)
      h5a_close, att_id
      att_id=h5a_open_name(group_id, 'ndymax')
         ndymax=h5a_read(att_id)
         ndymax=ndymax(0)
      h5a_close, att_id
      att_id=h5a_open_name(group_id, 'ndzmax')
         ndzmax=h5a_read(att_id)
         ndzmax=ndzmax(0)
      h5a_close, att_id
   h5g_close, group_id
;
   group_id = h5g_open(file_id, 'coordinates')
      dset_id=h5d_open(group_id, 'r')
         r=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'x')
         x=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'y')
         y=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'z')
         z=h5d_read(dset_id)
      h5d_close, dset_id
   h5g_close, group_id
;
   group_id = h5g_open(file_id, 'model')
      dset_id=h5d_open(group_id, 'velr1d_cr')
         velr1d_cr=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'opalbar1d_cr')
         opalbar1d_cr=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 't1d_cr')
         t1d_cr=h5d_read(dset_id)
      h5d_close, dset_id
   h5g_close, group_id
;
   group_id = h5g_open(file_id, 'solution1d')
      dset_id=h5d_open(group_id, 'ssobo1d_crx')
         ssobo1d_crx=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'ssobo1d_cry')
         ssobo1d_cry=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'ssobo1d_crz')
         ssobo1d_crz=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'ssobo1d')
         ssobo1d=h5d_read(dset_id)
      h5d_close, dset_id
   h5g_close, group_id
;
   group_id = h5g_open(file_id, 'solution3d')
      dset_id=h5d_open(group_id, 'ssobo3d')
         ssobo3d=h5d_read(dset_id)
      h5d_close, dset_id
   h5g_close, group_id
;
h5f_close, file_id
;
;normalize everything to xic1
ssobo1d_crx=ssobo1d_crx/xic1
ssobo1d_cry=ssobo1d_cry/xic1
ssobo1d_crz=ssobo1d_crz/xic1
ssobo1d=ssobo1d/xic1
ssobo3d=ssobo3d/xic1
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
   ymin=min([min(ssobo1d_crx*r^2), min(ssobo1d_cry*r^2), min(ssobo1d_crz*r^2), min(ssobo1d*r^2)])
   ymax=max([max(ssobo1d_crx*r^2), max(ssobo1d_cry*r^2), max(ssobo1d_crz*r^2), max(ssobo1d*r^2)])
   ylim=[ymin, ymax]
   ylim=[0.,0.6]
endif
;
;----------------------------title-strings------------------------------
;
titlestr=textoidl('line-strength k_L=') + string(kline, format='(f9.5)')
;
xtitlestr=textoidl('r')
ytitlestr=textoidl('S_L*r^2/I_c')
;
;-----------------------------------------------------------------------
;
if(not keyword_set(windx)) then windx=0
;
if(keyword_set(oname)) then begin
   oname='ps_files/benchmark07.ps'
   print, "writing output to: ", oname
   set_plot,'ps'
   device,file=oname, $; xsize=19., ysize=26.7, xoffset=1., yoffset=1. , $
    decomposed=0, color=1, bits_per_pixel=8
endif else begin
   window, windx;, xsize=950, ysize=1200
   device, decomposed=0
   windx=windx+1
endelse
;
loadct_lev, ci_black, ci_blue, ci_red, ci_cyan, ci_green, ci_magenta, ci_yellow, ci_white
;
;

r3d=fltarr(ndxmax,ndymax,ndzmax)
for i=0, ndxmax-1 do begin
   for j=0, ndymax-1 do begin
      for k=0, ndzmax-1 do begin
         r3d(i,j,k)=sqrt(x(i)^2+y(j)^2+z(k)^2)
      endfor
   endfor
endfor
indx=where(r3d ge 1.d0 and r3d lt xmax)
;

plot, r3d(indx), ssobo3d(indx)*r3d(indx)^2, $
      psym=3, $
      xrange=xlim, $
      yrange=ylim, $
      title=titleStr, $
      xtitle=xtitleStr, $
      ytitle=ytitleStr, $
      charsize=2.
;
oplot, r, ssobo1d_crx*r^2, $
      color=ci_blue, $
      thick=3., $
      line=0
oplot, r, ssobo1d_cry*r^2, $
      color=ci_red, $
      thick=3., $
      line=0
oplot, r, ssobo1d_crz*r^2, $
      color=ci_green, $
      thick=3., $
      line=0
;
oplot, r, ssobo1d*r^2, $
      color=ci_magenta, $
      psym=1
      thick=3.

;
lstr1='1d sobolev'
lstr2='3d sobolev (x-axis)'
lstr3='3d sobolev (y-axis)'
lstr4='3d sobolev (z-axis)'
lstr5='3d sobolev'
;
if(keyword_set(oname)) then begin
   legend, [lstr1, lstr2, lstr3, lstr4, lstr5], $
           psym=[1,0,0,0,3], $
           linestyle=[0,0,0,0,0], $
           color=[ci_magenta, ci_blue,ci_red,ci_green,ci_black], $
           textcolor=ci_black, $
           /right_legend
endif else begin
   legend, [lstr1, lstr2, lstr3, lstr4, lstr5], $
           psym=[1,0,0,0,3], $
           linestyle=[0,0,0,0,0], $
           color=[ci_magenta,ci_blue,ci_red,ci_green,ci_white], $
           textcolor=ci_white, $
           /right_legend
endelse
;
;
if keyword_set(oname) then begin
   device, /close
   set_plot,'x'
endif
;
;
;
end
