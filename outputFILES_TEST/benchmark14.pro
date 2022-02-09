pro benchmark14, dir=dir, windx=windx, oname=oname, xlim=xlim, ylim=ylim
;
;plot line source functions for spherically symmetric models
;
;------------------read all information from hdf5-file------------------
;
if(not keyword_set(dir)) then dir='.'
fname=dir+'/benchmark14.h5'
;
file_id = h5f_open(fname)
;
   group_id = h5g_open(file_id, 'input_parameters')
      att_id=h5a_open_name(group_id, 'opt_opal')
         opt_opal=h5a_read(att_id)
         opt_opal=opt_opal(0)
      h5a_close, att_id
      att_id=h5a_open_name(group_id, 'kcont')
         kcont=h5a_read(att_id)
         kcont=kcont(0)
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
      att_id=h5a_open_name(group_id, 'eps_cont')
         eps_cont=h5a_read(att_id)
         eps_cont=eps_cont(0)
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
   group_id = h5g_open(file_id, 'convergence_behaviour')
      att_id=h5a_open_name(group_id, 'itmaxc')
         itmaxc=h5a_read(att_id)
         itmaxc=itmaxc(0)
      h5a_close, att_id
      att_id=h5a_open_name(group_id, 'itmaxl')
         itmaxl=h5a_read(att_id)
         itmaxl=itmaxl(0)
      h5a_close, att_id
      dset_id=h5d_open(group_id, 'epsmaxc_sc')
         epsmaxc_sc=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'epsmaxl_sc')
         epsmaxl_sc=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'epsmaxc_fvm')
         epsmaxc_fvm=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'epsmaxl_fvm')
         epsmaxl_fvm=h5d_read(dset_id)
      h5d_close, dset_id
   h5g_close, group_id
;
   group_id = h5g_open(file_id, 'model_cr')
      dset_id=h5d_open(group_id, 't1d_jo')
         t1d_jo=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'opac1d_jo')
         opac1d_jo=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'opalbar1d_jo')
         opalbar1d_jo=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'velr1d')
         velr1d=h5d_read(dset_id)
      h5d_close, dset_id
   h5g_close, group_id
;
   group_id = h5g_open(file_id, 'solution_cr')
      dset_id=h5d_open(group_id, 'scont1d_jo')
         scont1d_jo=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'sline1d_jo')
         sline1d_jo=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'ssobo1d_jo')
         ssobo1d_jo=h5d_read(dset_id)
      h5d_close, dset_id
   h5g_close, group_id
;
   group_id = h5g_open(file_id, 'model3d')
      dset_id=h5d_open(group_id, 'mask3d')
         mask3d=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 't3d')
         temp3d=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'opac3d')
         opac3d=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'opalbar3d')
         opalbar3d=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'velx3d')
         velx3d=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'vely3d')
         vely3d=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'velz3d')
         velz3d=h5d_read(dset_id)
      h5d_close, dset_id
   h5g_close, group_id
;
   group_id = h5g_open(file_id, 'solution3d')
      dset_id=h5d_open(group_id, 'mint3d_sc')
         mint3d_sc=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'mint3d_fvm')
         mint3d_fvm=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'mintbar3d_sc')
         mintbar3d_sc=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'mintbar3d_fvm')
         mintbar3d_fvm=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'scont3d_sc')
         scont3d_sc=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'scont3d_fvm')
         scont3d_fvm=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'sline3d_sc')
         sline3d_sc=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'sline3d_fvm')
         sline3d_fvm=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'ssobo3d')
         ssobo3d=h5d_read(dset_id)
      h5d_close, dset_id
   h5g_close, group_id
;
h5f_close, file_id
;
;normalize everything to xic1
scont3d_sc=scont3d_sc/xic1
scont3d_fvm=scont3d_fvm/xic1
sline3d_sc=sline3d_sc/xic1
sline3d_fvm=sline3d_fvm/xic1
ssobo3d=ssobo3d/xic1
scont1d_jo=scont1d_jo/xic1
sline1d_jo=sline1d_jo/xic1
ssobo1d_jo=ssobo1d_jo/xic1
;
;
;define radius
r3d=fltarr(ndxmax,ndymax,ndzmax)
velr3d=fltarr(ndxmax,ndymax,ndzmax)
for i=0, ndxmax-1 do begin
   for j=0, ndymax-1 do begin
      for k=0, ndzmax-1 do begin
         case mask3d(i,j,k) of
            1: begin
                  r3d(i,j,k)=sqrt(x(i)^2+y(j)^2+z(k)^2)
                  velr3d(i,j,k)=sqrt(velx3d(i,j,k)^2+vely3d(i,j,k)^2+velz3d(i,j,k)^2)
                  break
               end
            2: begin
                  r3d(i,j,k)=sqrt(x(i)^2+y(j)^2+z(k)^2)
                  velr3d(i,j,k)=sqrt(velx3d(i,j,k)^2+vely3d(i,j,k)^2+velz3d(i,j,k)^2)
                  break
               end
            3: begin
                  r3d(i,j,k)=sqrt(x(i)^2+y(j)^2+z(k)^2)
                  velr3d(i,j,k)=sqrt(velx3d(i,j,k)^2+vely3d(i,j,k)^2+velz3d(i,j,k)^2)
                  break
               end
            else: begin
                  r3d(i,j,k)=1.d0 ;dummy value
                  velr3d(i,j,k)=0.d0
               end
         endcase
      endfor
   endfor
endfor

;
;calculate sline1d_jo on 3d grid to obtain errors
sline3d_jo=fltarr(ndxmax,ndymax,ndzmax)*0.d0
scont3d_jo=fltarr(ndxmax,ndymax,ndzmax)*0.d0
for i=0, ndxmax-1 do begin
   for j=0, ndymax-1 do begin
      for k=0, ndzmax-1 do begin
         switch mask3d(i,j,k) of
            1: begin 
                  find_indx, r3d(i,j,k), r, nr, iim1, ii
                  sline3d_jo(i,j,k) = interpol_ypl(r(iim1),r(ii),sline1d_jo(iim1),sline1d_jo(ii),r3d(i,j,k))
                  scont3d_jo(i,j,k) = interpol_ypl(r(iim1),r(ii),scont1d_jo(iim1),scont1d_jo(ii),r3d(i,j,k))
                  break
               end
            2: begin 
                  find_indx, r3d(i,j,k), r, nr, iim1, ii
                  sline3d_jo(i,j,k) = interpol_ypl(r(iim1),r(ii),sline1d_jo(iim1),sline1d_jo(ii),r3d(i,j,k))
                  scont3d_jo(i,j,k) = interpol_ypl(r(iim1),r(ii),scont1d_jo(iim1),scont1d_jo(ii),r3d(i,j,k))
                  break
               end
            3: begin 
                  find_indx, r3d(i,j,k), r, nr, iim1, ii
                  sline3d_jo(i,j,k) = interpol_ypl(r(iim1),r(ii),sline1d_jo(iim1),sline1d_jo(ii),r3d(i,j,k))
                  scont3d_jo(i,j,k) = interpol_ypl(r(iim1),r(ii),scont1d_jo(iim1),scont1d_jo(ii),r3d(i,j,k))
                  break
               end
         endswitch
      endfor
   endfor
endfor
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
   indx=where(mask3d eq 1 or mask3d eq 2 or mask3d eq 3)
   ymin1=min(sline3d_fvm(indx)*r3d(indx)^2)
   ymax1=max(sline3d_sc(indx)*r3d(indx)^2)
   ymax2=max(sline3d_fvm(indx)*r3d(indx)^2)
   ymin=min([ymin1, min(sline1d_jo*r^2)])
   ymax=max([ymax1, ymax2, max(sline1d_jo*r^2)])
   ylim=[ymin, ymax]

   ymin1=min(scont3d_fvm(indx)*r3d(indx)^2)
   ymin2=min(scont3d_sc(indx)*r3d(indx)^2)
   ymin3=min(sline1d_jo*r^2)
   ymax1=max(scont3d_fvm(indx)*r3d(indx)^2)
   ymax2=max(scont3d_sc(indx)*r3d(indx)^2)
   ymax3=max(sline1d_jo*r^2)
   ylim_cont = [min([ymin1,ymin2,ymin3]),max([ymax1,ymax2,ymax3])]
endif
;
;----------------------------title-strings------------------------------
;
if(opt_opal eq 0) then begin
   titlestr=textoidl('k_L,k_c=') + string(kline, format='(f6.3)') + ', ' + string(kcont, format='(f6.3)')
endif else begin
   titlestr=textoidl('\kappa_0, \alpha, k_c=') + string(kappa0, format='(f6.3)') + $
            ', ' + string(alpha, format='(f6.3)') + ', ' + string(kcont, format='(f6.3)')
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
   oname='ps_files/benchmark13_a.ps'
   print, "writing output to: ", oname
   set_plot,'ps'
   device,file=oname, xsize=19., ysize=26.7, xoffset=1., yoffset=1. , $
    decomposed=0, color=1, bits_per_pixel=8
endif else begin
   window, windx, xsize=1100, ysize=1200
   device, decomposed=0
   windx=windx+1
endelse
;
loadct_lev, ci_black, ci_blue, ci_red, ci_cyan, ci_green, ci_magenta, ci_yellow, ci_white
;
!p.multi=[9,3,3]
;
plot, r, sline1d_jo*r^2, $
      xrange=xlim, $
      yrange=ylim, $
      title=titleStr, $
      xtitle=xtitleStr, $
      ytitle=ytitleStr, $
      charsize=2., /xs, /ys, $
      thick=2.

oplot, r, ssobo1d_jo*r^2, $
      line=2

oplot, r3d, sline3d_sc*r3d^2, $
      color=ci_blue, $
      psym=3
oplot, r3d, sline3d_fvm*r3d^2, $
      color=ci_red, $
      psym=3
oplot, r3d, ssobo3d*r3d^2, $
      color=ci_green, $
      psym=3
;
;
;oplot for positive and negative x axis
ix0=0
ix1=ndxmax/2
ix2=ndxmax-1
iy=ndymax/2
iz=ndzmax/2
oplot, abs(x(ix0:ix1)), sline3d_sc(ix0:ix1,iy,iz)*x(ix0:ix1)^2, color=ci_green, line=0
oplot, x(ix1:ix2), sline3d_sc(ix1:ix2,iy,iz)*x(ix1:ix2)^2, color=ci_green, psym=1
;
;oplot for positive and negative y axis
ix=ndxmax/2
iy0=0
iy1=ndymax/2
iy2=ndymax-1
iz=ndzmax/2
oplot, abs(y(iy0:iy1)), sline3d_sc(ix,iy0:iy1,iz)*y(iy0:iy1)^2, color=ci_magenta, line=2
oplot, y(iy1:iy2), sline3d_sc(ix,iy1:iy2,iz)*y(iy1:iy2)^2, color=ci_magenta, psym=1
;
;oplot for positive and negative z-directions
ix=ndxmax/2
iy=ndymax/2
iz0=0
iz1=ndzmax/2
iz2=ndzmax-1
oplot, abs(z(iz0:iz1)), sline3d_sc(ix,iy,iz0:iz1)*z(iz0:iz1)^2, color=ci_red, line=3
oplot, z(iz1:iz2), sline3d_sc(ix,iy,iz1:iz2)*z(iz1:iz2)^2, color=ci_red, psym=1

;
lstr1='1d (Jo)'
lstr2='1d (Sobolev)'
lstr3='3d SC'
lstr4='3d FVM'
lstr5='3d (Sobolev)'
;
if(keyword_set(oname)) then begin
   legend, [lstr1, lstr2, lstr3, lstr4, lstr5], $
           psym=[0,0,3,3,3], $
           linestyle=[0,2,0,0,0], $
           color=[ci_black, ci_black,ci_blue,ci_red,ci_green], $
           textcolor=ci_black, $
           /right_legend
endif else begin
   legend, [lstr1, lstr2, lstr3, lstr4, lstr5], $
           psym=[0,0,3,3,3], $
           linestyle=[0,2,0,0,0], $
           color=[ci_white,ci_white,ci_blue,ci_red,ci_green], $
           textcolor=ci_white, $
           /right_legend
endelse
;
;
;
!p.multi=[8,3,3]
;
velr1d=velr1d/vmax
;velr3d=velr3d*100.d0*1.d5/vmax ;delete this factor 100*1.d5 for later use
velr3d=velr3d/vmax
plot, velr1d, sline1d_jo*r^2, $
      xrange=[0.,1.], $
      yrange=[0.,1.], $
      title=titleStr, $
      xtitle=textoidl('v'), $
      ytitle=ytitleStr, $
      charsize=2.

oplot, velr1d, ssobo1d_jo*r^2, $
      line=2

oplot, velr3d, sline3d_sc*r3d^2, $
      color=ci_blue, $
      psym=3
oplot, velr3d, sline3d_fvm*r3d^2, $
      color=ci_red, $
      psym=3
oplot, velr3d, ssobo3d*r3d^2, $
      color=ci_green, $
      psym=3
;
lstr1='1d (Jo)'
lstr2='1d (Sobolev)'
lstr3='3d SC'
lstr4='3d FVM'
lstr5='3d (Sobolev)'
;
if(keyword_set(oname)) then begin
   legend, [lstr1, lstr2, lstr3, lstr4, lstr5], $
           psym=[0,0,3,3,3], $
           linestyle=[0,2,0,0,0], $
           color=[ci_black, ci_black,ci_blue,ci_red,ci_green], $
           textcolor=ci_black, $
           /right_legend
endif else begin
   legend, [lstr1, lstr2, lstr3, lstr4, lstr5], $
           psym=[0,0,0,0,0], $
           linestyle=[0,2,3,3,3], $
           color=[ci_white,ci_white,ci_blue,ci_red,ci_green], $
           textcolor=ci_white, $
           /right_legend
endelse
;
;plot error in r-space
indx=where(mask3d eq 1 or mask3d eq 2 or mask3d eq 3)
ymin=min([min(sline3d_sc(indx)/sline3d_jo(indx)), min(sline3d_fvm(indx)/sline3d_jo(indx))])
ymax=max([max(sline3d_sc(indx)/sline3d_jo(indx)), max(sline3d_fvm(indx)/sline3d_jo(indx))])
ylim=[ymin, ymax]
;
!p.multi=[6,3,3]
plot, [0.,0.], [0.,0.], $
      xrange=xlim, $
      yrange=ylim, $
      title=titleStr, $
      xtitle=xtitleStr, $
      ytitle=textoidl('S_L(3d)/S_L(1d)'), $
      charsize=2.
oplot, r3d, sline3d_sc/sline3d_jo, $
      color=ci_blue, $
      psym=3
oplot, r3d, sline3d_fvm/sline3d_jo, $
      color=ci_red, $
      psym=3
lstr1='3d SC'
lstr2='3d FVM'
if(keyword_set(oname)) then begin
   legend, [lstr1, lstr2], $
           psym=[3,3], $
           linestyle=[0,0], $
           color=[ci_blue,ci_red], $
           textcolor=ci_black, $
           /right_legend
endif else begin
   legend, [lstr1, lstr2], $
           psym=[3,3], $
           linestyle=[0,0], $
           color=[ci_blue,ci_red], $
           textcolor=ci_white, $
           /right_legend
endelse
;
;plot error in velocity-space
!p.multi=[5,3,3]
plot, [0.,0.], [0.,0.], $
      xrange=[0.,1.], $
      yrange=ylim, $
      title=titleStr, $
      xtitle=textoidl('v'), $
      ytitle=textoidl('S_L(3d)/S_L(1d)'), $
      charsize=2.
oplot, velr3d, sline3d_sc/sline3d_jo, $
      color=ci_blue, $
      psym=3
oplot, velr3d, sline3d_fvm/sline3d_jo, $
      color=ci_red, $
      psym=3
lstr1='3d SC'
lstr2='3d FVM'
if(keyword_set(oname)) then begin
   legend, [lstr1, lstr2], $
           psym=[3,3], $
           linestyle=[0,0], $
           color=[ci_blue,ci_red], $
           textcolor=ci_black, $
           /right_legend
endif else begin
   legend, [lstr1, lstr2], $
           psym=[3,3], $
           linestyle=[0,0], $
           color=[ci_blue,ci_red], $
           textcolor=ci_white, $
           /right_legend
endelse
;
;
;
;plot continuum source function
!p.multi=[7,3,3]
;
plot, r, scont1d_jo*r^2, $
      xrange=xlim_cont, $
      yrange=ylim_cont, $
      title=titleStr, $
      xtitle=textoidl('r'), $
      ytitle=textoidl('S_c*r^2/I_c'), $
      charsize=2., /xs, /ys, $
      thick=2.

oplot, r3d, scont3d_sc*r3d^2, $
      color=ci_blue, $
      psym=3
oplot, r3d, scont3d_fvm*r3d^2, $
      color=ci_red, $
      psym=3
;
lstr1='1d (Jo)'
lstr2='3d SC'
lstr3='3d FVM'
;
if(keyword_set(oname)) then begin
   legend, [lstr1, lstr2, lstr3], $
           psym=[0,3,3], $
           linestyle=[0,0,0], $
           color=[ci_black,ci_blue,ci_red], $
           textcolor=ci_black, $
           /right_legend
endif else begin
   legend, [lstr1, lstr2, lstr3], $
           psym=[0,3,3], $
           linestyle=[0,0,0], $
           color=[ci_white,ci_blue,ci_red], $
           textcolor=ci_white, $
           /right_legend
endelse
;
;
;
!p.multi=[4,3,3]
;
plot, velr1d, scont1d_jo*r^2, $
      xrange=[0.,1.], $
      yrange=[0.,1.], $
      title=titleStr, $
      xtitle=textoidl('v'), $
      ytitle=textoidl('S_c*r^2/I_c'), $
      charsize=2.
oplot, velr3d, scont3d_sc*r3d^2, $
      color=ci_blue, $
      psym=3
oplot, velr3d, scont3d_fvm*r3d^2, $
      color=ci_red, $
      psym=3
;
lstr1='1d (Jo)'
lstr3='3d SC'
lstr4='3d FVM'
;
;
if(keyword_set(oname)) then begin
   legend, [lstr1, lstr2, lstr3], $
           psym=[0,3,3], $
           linestyle=[0,0,0], $
           color=[ci_black,ci_blue,ci_red], $
           textcolor=ci_black, $
           /right_legend
endif else begin
   legend, [lstr1, lstr2, lstr3], $
           psym=[0,3,3], $
           linestyle=[0,0,0], $
           color=[ci_white,ci_blue,ci_red], $
           textcolor=ci_white, $
           /right_legend
endelse
;
;
;
;only to check if temperature model is consistent
!p.multi=[4,4,3]
plot, r, t1d_jo, $
      xtitle='r', $
      ytitle=textoidl('T'), $
      title='consistency check', $
      charsize=2., $
      yrange=[min(t1d_jo)-0.1*(max(t1d_jo)-min(t1d_jo)),max(t1d_jo)]
oplot, r3d, temp3d, $
       psym=3, $
       color=ci_blue
lstr1='temperature (1d)'
lstr2='temperature (3d)'
if(keyword_set(oname)) then begin
   legend, [lstr1, lstr2], $
         psym=[0,3], $
         linestyle=[0,0], $
         color=[ci_black,ci_blue], $
         textcolor=ci_black, $
         /right_legend
endif else begin
   legend, [lstr1, lstr2], $
         psym=[0,3], $
         linestyle=[0,0], $
         color=[ci_white,ci_blue], $
         textcolor=ci_white, $
         /right_legend
endelse
;
;only to check if line opacity model is consistent
!p.multi=[3,4,3]
plot, r, opalbar1d_jo, $
      xtitle='r', $
      ytitle=textoidl('\chi bar'), $
      title='consistency check', $
      /ylog, $
      charsize=2.
oplot, r3d, opalbar3d, $
       psym=3, $
       color=ci_blue
lstr1='1d'
lstr2='3d'
if(keyword_set(oname)) then begin
   legend, [lstr1, lstr2], $
         psym=[0,3], $
         linestyle=[0,0], $
         color=[ci_black,ci_blue], $
         textcolor=ci_black, $
         /right_legend
endif else begin
   legend, [lstr1, lstr2], $
         psym=[0,1], $
         linestyle=[0,0], $
         color=[ci_white,ci_blue], $
         textcolor=ci_white, $
         /right_legend
endelse
;
;only to check if velocity model is consistent
!p.multi=[2,4,3]
plot, r, velr1d, $
      xtitle='r', $
      ytitle=textoidl('v_r / v_{max}'), $
      title='consistency check', $
      charsize=2.
oplot, r3d, velr3d, psym=3, color=ci_blue
lstr1='1d'
lstr2='3d'
if(keyword_set(oname)) then begin
   legend, [lstr1, lstr2], $
         psym=[0,3], $
         linestyle=[0,0], $
         color=[ci_black,ci_blue], $
         textcolor=ci_black, $
         /right_legend, /bottom
endif else begin
   legend, [lstr1, lstr2], $
         psym=[0,1], $
         linestyle=[0,0], $
         color=[ci_white,ci_blue], $
         textcolor=ci_white, $
         /right_legend, /bottom
endelse
;
;
;
;only to check if continuum opacity model is consistent
!p.multi=[1,4,3]
plot, r, opac1d_jo, $
      xtitle='r', $
      ytitle=textoidl('\chi_c'), $
      title='consistency check', $
      /ylog, $
      charsize=2.
oplot, r3d, opac3d, $
       psym=3, $
       color=ci_blue
lstr1='1d'
lstr2='3d'
if(keyword_set(oname)) then begin
   legend, [lstr1, lstr2], $
         psym=[0,3], $
         linestyle=[0,0], $
         color=[ci_black,ci_blue], $
         textcolor=ci_black, $
         /right_legend
endif else begin
   legend, [lstr1, lstr2], $
         psym=[0,1], $
         linestyle=[0,0], $
         color=[ci_white,ci_blue], $
         textcolor=ci_white, $
         /right_legend
endelse
;
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
   oname='ps_files/benchmark13_b.ps'
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
iternrc=indgen(itmaxc)
epsmaxl_sc=abs(epsmaxl_sc)
epsmaxl_fvm=abs(epsmaxl_fvm)
iternrl=indgen(itmaxl)
;
xmaxl1=max(where(epsmaxl_sc gt 0))
xmaxl2=max(where(epsmaxl_fvm gt 0))
xmaxc1=max(where(epsmaxc_sc gt 0))
xmaxc2=max(where(epsmaxc_fvm gt 0))
;
yminl1=epsmaxl_sc(xmaxl1)
yminl2=epsmaxl_fvm(xmaxl2)
yminc1=epsmaxc_sc(xmaxc1)
yminc2=epsmaxc_fvm(xmaxc2)
;
ymaxl1=max(epsmaxc_sc)
ymaxl2=max(epsmaxc_fvm)
ymaxc1=max(epsmaxc_sc)
ymaxc2=max(epsmaxc_fvm)
;
xmin=0
;
xliml=[xmin,max([xmaxl1,xmaxl2])]
yliml=[min([yminl1,yminl2]),max([ymaxl1,ymaxl2])]
xlimc=[xmin,max([xmaxc1,xmaxc2])]
ylimc=[min([yminc1,yminc2]),max([ymaxc1,ymaxc2])]
;
ytitleStr=textoidl('((S_L^{(k-1)} - S_L^{(k)}) /S_L^{(k)})_{max}')
xtitleStr='# iterations'
titleStr='LINE'
;
!p.multi=[2,2,1]
plot, iternrl, epsmaxl_sc, $
      /ylog, $
      line=0, $
      yrange=yliml, $
      xrange=xliml, $
      title=titlestr, $
      xtitle=xtitlestr, $
      ytitle=ytitlestr
oplot, iternrl, epsmaxl_fvm, $
      line=2
lstr1='3d short characteristics'
lstr2='3d finite volume'
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
;
;
;
ytitleStr=textoidl('((S_c^{(k-1)} - S_c^{(k)}) /S_c^{(k)})_{max}')
xtitleStr='# iterations'
titleStr='CONTINUUM'
;
!p.multi=[1,2,1]
plot, iternrc, epsmaxc_sc, $
      /ylog, $
      line=0, $
      yrange=ylimc, $
      xrange=xlimc, $
      title=titlestr, $
      xtitle=xtitlestr, $
      ytitle=ytitlestr
oplot, iternrc, epsmaxc_fvm, $
      line=2
lstr1='3d short characteristics'
lstr2='3d finite volume'
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


!p.multi=0

if keyword_set(oname) then begin
   device, /close
   set_plot,'x'
endif
;
end
