pro test_interpolation, windx=WINDX, oname=ONAME
;
;plots contours of searchlight along a direction
;
;------------------read all information from hdf5-file------------------
;
fname='test_interpolation.h5'
;
file_id = h5f_open(fname)
;
   group_id = h5g_open(file_id, 'dimensions')
      att_id=h5a_open_name(group_id, 'n_data')
         n_data=h5a_read(att_id)
         n_data=n_data(0)
      h5a_close, att_id
      att_id=h5a_open_name(group_id, 'n_grid')
         n_grid=h5a_read(att_id)
         n_grid=n_grid(0)
      h5a_close, att_id
   h5g_close, group_id
;
   group_id = h5g_open(file_id, 'coordinates')
      dset_id=h5d_open(group_id, 'grid_data')
         s_data=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'grid_interp')
         s_interp=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'grid_theo')
         s_theo=h5d_read(dset_id)
      h5d_close, dset_id
   h5g_close, group_id
;
   group_id = h5g_open(file_id, 'function1')
      dset_id=h5d_open(group_id, 'y_interp11')
         y_interp11=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'y_interp12')
         y_interp12=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'y_interp13')
         y_interp13=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'y_interp14')
         y_interp14=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'y_interp15')
         y_interp15=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'y_interp16')
         y_interp16=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'y_interp17')
         y_interp17=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'y_interp18')
         y_interp18=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'y_data1')
         y_data1=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'y_theo1')
         y_theo1=h5d_read(dset_id)
      h5d_close, dset_id
   h5g_close, group_id
;
   group_id = h5g_open(file_id, 'function2')
      dset_id=h5d_open(group_id, 'y_interp21')
         y_interp21=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'y_interp22')
         y_interp22=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'y_interp23')
         y_interp23=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'y_interp24')
         y_interp24=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'y_interp25')
         y_interp25=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'y_interp26')
         y_interp26=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'y_interp27')
         y_interp27=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'y_interp28')
         y_interp28=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'y_data2')
         y_data2=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'y_theo2')
         y_theo2=h5d_read(dset_id)
      h5d_close, dset_id
   h5g_close, group_id
;
   group_id = h5g_open(file_id, 'function3')
      dset_id=h5d_open(group_id, 'y_interp31')
         y_interp31=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'y_interp32')
         y_interp32=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'y_interp33')
         y_interp33=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'y_interp34')
         y_interp34=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'y_interp35')
         y_interp35=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'y_interp36')
         y_interp36=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'y_interp37')
         y_interp37=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'y_interp38')
         y_interp38=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'y_data3')
         y_data3=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'y_theo3')
         y_theo3=h5d_read(dset_id)
      h5d_close, dset_id
   h5g_close, group_id
;
   group_id = h5g_open(file_id, 'function4')
      dset_id=h5d_open(group_id, 'y_interp41')
         y_interp41=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'y_interp42')
         y_interp42=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'y_interp43')
         y_interp43=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'y_interp44')
         y_interp44=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'y_interp45')
         y_interp45=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'y_interp46')
         y_interp46=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'y_interp47')
         y_interp47=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'y_interp48')
         y_interp48=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'y_data4')
         y_data4=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'y_theo4')
         y_theo4=h5d_read(dset_id)
      h5d_close, dset_id
   h5g_close, group_id
;
   group_id = h5g_open(file_id, 'function5')
      dset_id=h5d_open(group_id, 'y_interp51')
         y_interp51=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'y_interp52')
         y_interp52=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'y_interp53')
         y_interp53=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'y_interp54')
         y_interp54=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'y_interp55')
         y_interp55=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'y_interp56')
         y_interp56=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'y_interp57')
         y_interp57=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'y_interp58')
         y_interp58=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'y_data5')
         y_data5=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'y_theo5')
         y_theo5=h5d_read(dset_id)
      h5d_close, dset_id
   h5g_close, group_id
;
   group_id = h5g_open(file_id, 'function6')
      dset_id=h5d_open(group_id, 'y_interp61')
         y_interp61=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'y_interp62')
         y_interp62=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'y_interp63')
         y_interp63=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'y_interp64')
         y_interp64=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'y_interp65')
         y_interp65=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'y_interp66')
         y_interp66=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'y_interp67')
         y_interp67=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'y_interp68')
         y_interp68=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'y_data6')
         y_data6=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'y_theo6')
         y_theo6=h5d_read(dset_id)
      h5d_close, dset_id
   h5g_close, group_id
;
   group_id = h5g_open(file_id, 'function7')
      dset_id=h5d_open(group_id, 'y_interp71')
         y_interp71=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'y_interp72')
         y_interp72=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'y_interp73')
         y_interp73=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'y_interp74')
         y_interp74=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'y_interp75')
         y_interp75=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'y_interp76')
         y_interp76=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'y_interp77')
         y_interp77=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'y_interp78')
         y_interp78=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'y_data7')
         y_data7=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'y_theo7')
         y_theo7=h5d_read(dset_id)
      h5d_close, dset_id
   h5g_close, group_id
;
   group_id = h5g_open(file_id, 'function8')
      dset_id=h5d_open(group_id, 'y_interp81')
         y_interp81=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'y_interp82')
         y_interp82=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'y_interp83')
         y_interp83=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'y_interp84')
         y_interp84=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'y_interp85')
         y_interp85=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'y_interp86')
         y_interp86=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'y_interp87')
         y_interp87=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'y_interp88')
         y_interp88=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'y_data8')
         y_data8=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'y_theo8')
         y_theo8=h5d_read(dset_id)
      h5d_close, dset_id
   h5g_close, group_id
;
h5f_close, file_id
;
;-----------------------------------------------------------------------
;-----------------------define range------------------------------------
;
xlim=[min(s_theo),max(s_theo)]
;
ylim1=[min(y_theo1),max(y_theo1)]
ylim2=[min(y_theo2),max(y_theo2)]
ylim3=[min(y_theo3),max(y_theo3)]
ylim4=[min(y_theo4),max(y_theo4)]
ylim5=[min(y_theo5),max(y_theo5)]
ylim6=[min(y_theo6),max(y_theo6)]
ylim7=[min(y_theo7),max(y_theo7)]
ylim8=[min(y_theo8),max(y_theo8)]
;
;----------------------------title-strings------------------------------
;
titlestr1=textoidl('f(x)=m*x+t')
titlestr2=textoidl('f(x)=exp(-m*x+t)')
titlestr3=textoidl('step-function')
titlestr4=textoidl('step-function')
titlestr5=textoidl('step-function')
titlestr6=textoidl('f(x)=(x-xim1)^2')
titlestr7=textoidl('f(x)=(x-xi)^2')
titlestr8=textoidl('f(x)=(x-(xi-(xi-xim)/4))^2')
;
xtitlestr=textoidl('x')
ytitlestr=textoidl('f(x)')
;
lstr0='f(x)'
lstr1='linear interpolation'
lstr2='quadratic bezier spline'
lstr3='quadratic bezier spline (monotonic)'
lstr4='cubic catmull rom spline (4 points)'
lstr5='cubic spline (3 points)'
lstr6='cubic spline (4 points, monotonic Ibgui 2013)'
lstr7='cubic spline (4 points, monotonic Steffen 1990)'
lstr8='cubic spline (4 points, weighted mean derivatives)'
;
;-----------------------------------------------------------------------
;
if(not keyword_set(windx)) then windx=0
;
if(keyword_set(oname)) then begin
   oname='ps_files/test_interpolation.ps'
   print, "writing output to: ", oname
   set_plot,'ps'
   device,file=oname, $ ;xsize=19., ysize=26.7, xoffset=1., yoffset=1. , $
    decomposed=0, color=1, bits_per_pixel=8
endif else begin
   window, windx, xsize=1100, ysize=1200
   device, decomposed=0
   windx=windx+1
endelse
;
loadct_lev, ci_black, ci_blue, ci_red, ci_cyan, ci_green, ci_magenta, ci_yellow, ci_white
;

;
!p.multi=[9,3,3]
if(keyword_set(oname)) then begin
   plot, [0.,0.], [0.,0.], color=ci_white, charsize=2. ;blank plot
   legend, [lstr0,lstr1,lstr2,lstr3,lstr4,lstr5,lstr6,lstr7,lstr8], $
         psym=[0,0,0,0,0,0,0,0,0], $
         linestyle=[0,0,0,0,0,2,2,2,2], $
         color=[ci_black,ci_blue,ci_red,ci_cyan,ci_green, $
                         ci_blue,ci_red,ci_cyan,ci_green], $
         textcolor=ci_black, $
         /right_legend
endif else begin
   plot, [0.,0.], [0.,0.], color=ci_black, charsize=2.  ;blank plot
   legend, [lstr0,lstr1,lstr2,lstr3,lstr4,lstr5,lstr6,lstr7,lstr8], $
         psym=[0,0,0,0,0,0,0,0,0], $
         linestyle=[0,0,0,0,0,2,2,2,2], $
         color=[ci_white,ci_blue,ci_red,ci_cyan,ci_green, $
                         ci_blue,ci_red,ci_cyan,ci_green], $
         textcolor=ci_white, $
         /right_legend
endelse
;
!p.multi=[8,3,3]
plot, s_theo, y_theo1, $
      xrange=xlim, $
      yrange=ylim1, $
      title=titleStr1, $
      xtitle=xtitleStr, $
      ytitle=ytitleStr, $
      charsize=2.
oplot, s_data, y_data1, $
      psym=2, color=ci_red
oplot, s_interp, y_interp11, $
      color=ci_blue, $
      line=0
oplot, s_interp, y_interp12, $
      color=ci_red, $
      line=0
oplot, s_interp, y_interp13, $
      color=ci_cyan, $
      line=0
oplot, s_interp, y_interp14, $
      color=ci_green, $
      line=0
oplot, s_interp, y_interp15, $
      color=ci_blue, $
      line=2
oplot, s_interp, y_interp16, $
      color=ci_red, $
      line=2
oplot, s_interp, y_interp17, $
      color=ci_cyan, $
      line=2
oplot, s_interp, y_interp18, $
      color=ci_green, $
      line=2
;
;
;
!p.multi=[7,3,3]
plot, s_theo, y_theo2, $
      xrange=xlim, $
      yrange=ylim2, $
      title=titleStr2, $
      xtitle=xtitleStr, $
      ytitle=ytitleStr, $
      charsize=2.
oplot, s_data, y_data2, $
      psym=2, color=ci_red
oplot, s_interp, y_interp21, $
      color=ci_blue, $
      line=0
oplot, s_interp, y_interp22, $
      color=ci_red, $
      line=0
oplot, s_interp, y_interp23, $
      color=ci_cyan, $
      line=0
oplot, s_interp, y_interp24, $
      color=ci_green, $
      line=0
oplot, s_interp, y_interp25, $
      color=ci_blue, $
      line=2
oplot, s_interp, y_interp26, $
      color=ci_red, $
      line=2
oplot, s_interp, y_interp27, $
      color=ci_cyan, $
      line=2
oplot, s_interp, y_interp28, $
      color=ci_green, $
      line=2
;
;
;
!p.multi=[6,3,3]
plot, s_theo, y_theo3, $
      xrange=xlim, $
      yrange=ylim3, $
      title=titleStr3, $
      xtitle=xtitleStr, $
      ytitle=ytitleStr, $
      charsize=2.
oplot, s_data, y_data3, $
      psym=2, color=ci_red
oplot, s_interp, y_interp31, $
      color=ci_blue, $
      line=0
oplot, s_interp, y_interp32, $
      color=ci_red, $
      line=0
oplot, s_interp, y_interp33, $
      color=ci_cyan, $
      line=0
oplot, s_interp, y_interp34, $
      color=ci_green, $
      line=0
oplot, s_interp, y_interp35, $
      color=ci_blue, $
      line=2
oplot, s_interp, y_interp36, $
      color=ci_red, $
      line=2
oplot, s_interp, y_interp37, $
      color=ci_cyan, $
      line=2
oplot, s_interp, y_interp38, $
      color=ci_green, $
      line=2
;
;
;
!p.multi=[5,3,3]
plot, s_theo, y_theo4, $
      xrange=xlim, $
      yrange=ylim4, $
      title=titleStr4, $
      xtitle=xtitleStr, $
      ytitle=ytitleStr, $
      charsize=2.
oplot, s_data, y_data4, $
      psym=2, color=ci_red
oplot, s_interp, y_interp41, $
      color=ci_blue, $
      line=0
oplot, s_interp, y_interp42, $
      color=ci_red, $
      line=0
oplot, s_interp, y_interp43, $
      color=ci_cyan, $
      line=0
oplot, s_interp, y_interp44, $
      color=ci_green, $
      line=0
oplot, s_interp, y_interp45, $
      color=ci_blue, $
      line=2
oplot, s_interp, y_interp46, $
      color=ci_red, $
      line=2
oplot, s_interp, y_interp47, $
      color=ci_cyan, $
      line=2
oplot, s_interp, y_interp48, $
      color=ci_green, $
      line=2
;
;
;
!p.multi=[4,3,3]
plot, s_theo, y_theo5, $
      xrange=xlim, $
      yrange=ylim5, $
      title=titleStr5, $
      xtitle=xtitleStr, $
      ytitle=ytitleStr, $
      charsize=2.
oplot, s_data, y_data5, $
      psym=2, color=ci_red
oplot, s_interp, y_interp51, $
      color=ci_blue, $
      line=0
oplot, s_interp, y_interp52, $
      color=ci_red, $
      line=0
oplot, s_interp, y_interp53, $
      color=ci_cyan, $
      line=0
oplot, s_interp, y_interp54, $
      color=ci_green, $
      line=0
oplot, s_interp, y_interp55, $
      color=ci_blue, $
      line=2
oplot, s_interp, y_interp56, $
      color=ci_red, $
      line=2
oplot, s_interp, y_interp57, $
      color=ci_cyan, $
      line=2
oplot, s_interp, y_interp58, $
      color=ci_green, $
      line=2
;
;
;
!p.multi=[3,3,3]
plot, s_theo, y_theo6, $
      xrange=xlim, $
      yrange=ylim6, $
      title=titleStr6, $
      xtitle=xtitleStr, $
      ytitle=ytitleStr, $
      charsize=2.
oplot, s_data, y_data6, $
      psym=2, color=ci_red
oplot, s_interp, y_interp61, $
      color=ci_blue, $
      line=0
oplot, s_interp, y_interp62, $
      color=ci_red, $
      line=0
oplot, s_interp, y_interp63, $
      color=ci_cyan, $
      line=0
oplot, s_interp, y_interp64, $
      color=ci_green, $
      line=0
oplot, s_interp, y_interp65, $
      color=ci_blue, $
      line=2
oplot, s_interp, y_interp66, $
      color=ci_red, $
      line=2
oplot, s_interp, y_interp67, $
      color=ci_cyan, $
      line=2
oplot, s_interp, y_interp68, $
      color=ci_green, $
      line=2
;
;
;
!p.multi=[2,3,3]
plot, s_theo, y_theo7, $
      xrange=xlim, $
      yrange=ylim7, $
      title=titleStr7, $
      xtitle=xtitleStr, $
      ytitle=ytitleStr, $
      charsize=2.
oplot, s_data, y_data7, $
      psym=2, color=ci_red
oplot, s_interp, y_interp71, $
      color=ci_blue, $
      line=0
oplot, s_interp, y_interp72, $
      color=ci_red, $
      line=0
oplot, s_interp, y_interp73, $
      color=ci_cyan, $
      line=0
oplot, s_interp, y_interp74, $
      color=ci_green, $
      line=0
oplot, s_interp, y_interp75, $
      color=ci_blue, $
      line=2
oplot, s_interp, y_interp76, $
      color=ci_red, $
      line=2
oplot, s_interp, y_interp77, $
      color=ci_cyan, $
      line=2
oplot, s_interp, y_interp78, $
      color=ci_green, $
      line=2
;
;
;
!p.multi=[1,3,3]
plot, s_theo, y_theo8, $
      xrange=xlim, $
      yrange=ylim8, $
      title=titleStr8, $
      xtitle=xtitleStr, $
      ytitle=ytitleStr, $
      charsize=2.
oplot, s_data, y_data8, $
      psym=2, color=ci_red
oplot, s_interp, y_interp81, $
      color=ci_blue, $
      line=0
oplot, s_interp, y_interp82, $
      color=ci_red, $
      line=0
oplot, s_interp, y_interp83, $
      color=ci_cyan, $
      line=0
oplot, s_interp, y_interp84, $
      color=ci_green, $
      line=0
oplot, s_interp, y_interp85, $
      color=ci_blue, $
      line=2
oplot, s_interp, y_interp86, $
      color=ci_red, $
      line=2
oplot, s_interp, y_interp87, $
      color=ci_cyan, $
      line=2
oplot, s_interp, y_interp88, $
      color=ci_green, $
      line=2



;
if keyword_set(oname) then begin
   device, /close
   set_plot,'x'
endif
;
!p.multi=0
;
end
