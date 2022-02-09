pro benchmark12, dir=dir, windx=windx, oname=oname, xlim=xlim, ylim=ylim
;
;plots contours of searchlight along a direction
;
;------------------read all information from hdf5-file------------------
;
if(not keyword_set(dir)) then dir='.'
;
fcont_joray=1
fcont_jomom=1
fcontr3d_sc=1
fcontth3d_sc=1
fcontphi3d_sc=1
fcontr3d_fvm=1
fcontth3d_fvm=1
fcontphi3d_fvm=1
;
read_benchmark12, dir=dir, kcont=kcont, eps_cont=eps_cont, xic1=xic1, $
                  ndxmax=ndxmax, ndymax=ndymax, ndzmax=ndzmax, nr=nr, $
                  xarr=x, yarr=y, zarr=z, rarr=r, $
                  itmaxc=itmaxc, epsmaxc_sc=epsmaxc_sc, epsmaxc_fvm=epsmaxc_fvm, $
                  t1d_jo=t1d_jo, opac1d_jo=opac1d_jo, $
                  mask3d=mask3d, t3d=temp3d, opac3d=opac3d, $
                  mint1d_joray=mint_joray, mint1d_jomom=mint_jomom, $
                  fcont1d_joray=fcont_joray, fcont1d_jomom=fcont_jomom, $
                  mint3d_sc=mint3d_sc, mint3d_fvm=mint3d_fvm, $
                  fcontr3d_sc=fcontr3d_sc, fcontth3d_sc=fcontth3d_sc, fcontphi3d_sc=fcontphi3d_sc, $
                  fcontr3d_fvm=fcontr3d_fvm, fcontth3d_fvm=fcontth3d_fvm, fcontphi3d_fvm=fcontphi3d_fvm
;
;
;trad=44.5d3
;xnue0=2.191983810829996D15
;xic2=bnue(xnue0,1.*trad)
;print, xic2, xic1, xic2/xic1
;xic1=xic2
;stop
;normalize everything to xic1
mint3d_sc=mint3d_sc/xic1
mint3d_fvm=mint3d_fvm/xic1
mint_joray=mint_joray/xic1
mint_jomom=mint_jomom/xic1
fcont_joray=fcont_joray/xic1
fcont_jomom=fcont_jomom/xic1
fcontr3d_sc=fcontr3d_sc/xic1
fcontth3d_sc=fcontth3d_sc/xic1
fcontphi3d_sc=fcontphi3d_sc/xic1
fcontr3d_fvm=fcontr3d_fvm/xic1
fcontth3d_fvm=fcontth3d_fvm/xic1
fcontphi3d_fvm=fcontphi3d_fvm/xic1

;
;define radius
r3d=fltarr(ndxmax,ndymax,ndzmax)
for i=0, ndxmax-1 do begin
   for j=0, ndymax-1 do begin
      for k=0, ndzmax-1 do begin
         case mask3d(i,j,k) of
            1: r3d(i,j,k)=sqrt(x(i)^2+y(j)^2+z(k)^2)
            2: r3d(i,j,k)=sqrt(x(i)^2+y(j)^2+z(k)^2)
            3: r3d(i,j,k)=sqrt(x(i)^2+y(j)^2+z(k)^2)
            else: r3d(i,j,k)=1.d0 ;dummy value
         endcase
      endfor
   endfor
endfor
;
;calculate mint_joray on 3d grid to obtain errors
mint3d_joray=fltarr(ndxmax,ndymax,ndzmax)*0.d0
for i=0, ndxmax-1 do begin
   for j=0, ndymax-1 do begin
      for k=0, ndzmax-1 do begin
         case mask3d(i,j,k) of
            1: begin 
                  find_indx, r3d(i,j,k), r, nr, iim1, ii
                  mint3d_joray(i,j,k) = interpol_ypl(r(iim1),r(ii),mint_joray(iim1),mint_joray(ii),r3d(i,j,k))
               end
            2: begin 
                  find_indx, r3d(i,j,k), r, nr, iim1, ii
                  mint3d_joray(i,j,k) = interpol_ypl(r(iim1),r(ii),mint_joray(iim1),mint_joray(ii),r3d(i,j,k))
               end
            3: begin 
                  find_indx, r3d(i,j,k), r, nr, iim1, ii
                  mint3d_joray(i,j,k) = interpol_ypl(r(iim1),r(ii),mint_joray(iim1),mint_joray(ii),r3d(i,j,k))
               end
            else: ;do nothing
         endcase
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
   ymin1=min(mint3d_fvm(indx)*r3d(indx)^2)
   ymax1=max(mint3d_sc(indx)*r3d(indx)^2)
   ymax2=max(mint3d_fvm(indx)*r3d(indx)^2)
   ymin=min([ymin1, min(mint_joray*r^2), min(mint_jomom*r^2)])
   ymax=max([ymax1, ymax2, max(mint_joray*r^2), min(mint_jomom*r^2)])
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
   oname='ps_files/benchmark12_a.ps'
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
;
plot_io, r, mint_joray*r^2, $
;plot, r, mint_joray*r^2, $
      xrange=xlim, $
      yrange=ylim, $
;      yrange=[0.d0,1.d0], $
      title=titleStr, $
      xtitle=xtitleStr, $
      ytitle=ytitleStr, $
      charsize=2., /xs, /ys, $
      thick=2.

oplot, r, mint_jomom*r^2, $
      line=2

oplot, r3d, mint3d_sc*r3d^2, $
      color=ci_blue, $
      psym=3

oplot, r3d, mint3d_fvm*r3d^2, $
      color=ci_red, $
      psym=3
;
;
;oplot for positive and negative x axis
ix0=0
ix1=ndxmax/2
ix2=ndxmax-1
iy=ndymax/2
iz=ndzmax/2
oplot, abs(x(ix0:ix1)), mint3d_sc(ix0:ix1,iy,iz)*x(ix0:ix1)^2, color=ci_green, line=0
oplot, x(ix1:ix2), mint3d_sc(ix1:ix2,iy,iz)*x(ix1:ix2)^2, color=ci_green, psym=1
;
;oplot for positive and negative y axis
ix=ndxmax/2
iy0=0
iy1=ndymax/2
iy2=ndymax-1
iz=ndzmax/2
oplot, abs(y(iy0:iy1)), mint3d_sc(ix,iy0:iy1,iz)*y(iy0:iy1)^2, color=ci_magenta, line=2
oplot, y(iy1:iy2), mint3d_sc(ix,iy1:iy2,iz)*y(iy1:iy2)^2, color=ci_magenta, psym=1
;
;oplot for positive and negative z-directions
ix=ndxmax/2
iy=ndymax/2
iz0=0
iz1=ndzmax/2
iz2=ndzmax-1
oplot, abs(z(iz0:iz1)), mint3d_sc(ix,iy,iz0:iz1)*z(iz0:iz1)^2, color=ci_red, line=3
oplot, z(iz1:iz2), mint3d_sc(ix,iy,iz1:iz2)*z(iz1:iz2)^2, color=ci_red, psym=1
;
;
;test along diagonal axis
;rtest=fltarr(ndxmax/2+1)
;itest=fltarr(ndxmax/2+1)
;k=0
;for i=ndxmax/2,ndxmax-1 do begin
;   print, k, ndxmax/2
;   rtest(k)=r3d(i,i,i)
;   itest(k)=mint3d_sc(i,i,i)
;   k=k+1
;endfor
;oplot, rtest, itest*rtest^2, color=ci_green
;stop
;
lstr1='1d (ray by ray)'
lstr2='1d (moments)'
lstr3='3d SC'
lstr4='3d FVM'
if(keyword_set(oname)) then begin
   legend, [lstr1, lstr2, lstr3, lstr4], $
           psym=[0,0,3,3], $
           linestyle=[0,2,0,0], $
           color=[ci_black, ci_black,ci_blue,ci_red], $
           textcolor=ci_black, $
           /right_legend
endif else begin
   legend, [lstr1, lstr2, lstr3, lstr4], $
           psym=[0,0,3,3], $
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
      charsize=2., /ylog, $
      thick=2.

oplot, 1.d0-1.d0/r, mint_jomom*r^2, $
      line=2

oplot, 1.d0-1.d0/r3d, mint3d_sc*r3d^2, $
      color=ci_blue, $
      psym=3

oplot, 1.d0-1.d0/r3d, mint3d_fvm*r3d^2, $
      color=ci_red, $
      psym=3
;
lstr1='1d (ray by ray)'
lstr2='1d (moments)'
lstr3='3d SC'
lstr4='3d FVM'
if(keyword_set(oname)) then begin
   legend, [lstr1, lstr2, lstr3, lstr4], $
           psym=[0,0,3,3], $
           linestyle=[0,2,0,0], $
           color=[ci_black, ci_black,ci_blue,ci_red], $
           textcolor=ci_black, $
           /right_legend
endif else begin
   legend, [lstr1, lstr2, lstr3, lstr4], $
           psym=[0,0,0,0], $
           linestyle=[0,2,3,3], $
           color=[ci_white,ci_white,ci_blue,ci_red], $
           textcolor=ci_white, $
           /right_legend
endelse
;
indx=where(mask3d eq 1 or mask3d eq 2 or mask3d eq 3)
ymin=min([min(mint3d_sc(indx)/mint3d_joray(indx)), min(mint3d_fvm(indx)/mint3d_joray(indx))])
ymax=max([max(mint3d_sc(indx)/mint3d_joray(indx)), max(mint3d_fvm(indx)/mint3d_joray(indx))])
ylim=[ymin, ymax]
xlim2=[0.5,2.5]
ylim2=[0.8,1.2]
!p.multi=[4,2,3]
plot, [0.,0.], [0.,0.], $
      xrange=xlim2, $
      yrange=ylim2, $
      title=titleStr, $
      xtitle=xtitleStr, $
      ytitle='J(3d)/J(1d)', $
      charsize=2.

oplot, r3d(indx), mint3d_sc(indx)/mint3d_joray(indx), $
      color=ci_blue, $
      psym=3
oplot, r3d(indx), mint3d_fvm(indx)/mint3d_joray(indx), $
      color=ci_red, $
      psym=3

;print, mint_sc*6.*r^2
;stop
;
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
!p.multi=[3,2,3]
plot, [0.,0.], [0.,0.], $
      xrange=[0.d0,1.d0], $
      yrange=ylim, $
      title=titleStr, $
      xtitle=textoidl('1-1/r'), $
      ytitle='J(3d)/J(1d)', $
      charsize=2.

oplot, 1.d0-1.d0/r3d(indx), mint3d_sc(indx)/mint3d_joray(indx), $
      color=ci_blue, $
      psym=3
oplot, 1.d0-1.d0/r3d(indx), mint3d_fvm(indx)/mint3d_joray(indx), $
      color=ci_red, $
      psym=3
;
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
;only to check if temperature model is consistent
!p.multi=[2,2,3]
plot, r, t1d_jo, $
      xtitle='r', $
      ytitle=textoidl('T'), $
      title='consistency check', $
      charsize=2., $
      yrange=[min(t1d_jo)-0.1*(max(t1d_jo)-min(t1d_jo)),max(t1d_jo)], $
      thick=2.

oplot, r3d(indx), temp3d(indx), $
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
;
;only to check if opacity model is consistent
!p.multi=[1,2,3]
plot, r, opac1d_jo, $
      xtitle='r', $
      ytitle=textoidl('\chi'), $
      title='consistency check', $
      /ylog, $
      charsize=2.
oplot, r3d(indx), opac3d(indx), $
       psym=3, $
       color=ci_blue
lstr1='opacity (1d)'
lstr2='opacity (3d)'
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
   oname='benchmark12_b.ps'
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
;ylim=[2.905d-3,2.915d-3]
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
      ytitle=ytitlestr, /ys
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
;--------------------------flux components------------------------------
;
if(keyword_set(oname)) then begin
   oname='benchmark12_c.ps'
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

titlestr=''
xtitlestr=textoidl('r')
ytitlestr=textoidl('F_r*r^2/I_c')

xlim = [min(r),max(r)]
indx=where(mask3d eq 1 or mask3d eq 2 or mask3d eq 3)
ymin1=min(fcontr3d_fvm(indx)*r3d(indx)^2)
ymax1=max(fcontr3d_sc(indx)*r3d(indx)^2)
ymax2=max(fcontr3d_fvm(indx)*r3d(indx)^2)
ymin=min([ymin1, min(fcont_joray*r^2), min(fcont_jomom*r^2)])
ymax=max([ymax1, ymax2, max(fcont_joray*r^2), min(fcont_jomom*r^2)])
ylim=[ymin, ymax]
;
!p.multi=[6,2,3]
;
;
plot, r, fcont_joray*r^2, $
      xrange=xlim, $
      yrange=ylim, $
      title=titleStr, $
      xtitle=xtitleStr, $
      ytitle=ytitleStr, $
      charsize=2., /xs, /ys, $
      thick=2.

oplot, r, fcont_jomom*r^2, $
      line=2

oplot, r3d, fcontr3d_sc*r3d^2, $
      color=ci_blue, $
      psym=3

oplot, r3d, fcontr3d_fvm*r3d^2, $
      color=ci_red, $
      psym=3
;
lstr1='1d (ray by ray)'
lstr2='1d (moments)'
lstr3='3d SC'
lstr4='3d FVM'
if(keyword_set(oname)) then begin
   legend, [lstr1, lstr2, lstr3, lstr4], $
           psym=[0,0,3,3], $
           linestyle=[0,2,0,0], $
           color=[ci_black, ci_black,ci_blue,ci_red], $
           textcolor=ci_black, $
           /right_legend
endif else begin
   legend, [lstr1, lstr2, lstr3, lstr4], $
           psym=[0,0,3,3], $
           linestyle=[0,2,0,0], $
           color=[ci_white,ci_white,ci_blue,ci_red], $
           textcolor=ci_white, $
           /right_legend
endelse
;
!p.multi=[5,2,3]
;
plot, 1.d0-1.d0/r, fcont_joray*r^2, $
      xrange=[0.d0,1.d0], $
      yrange=ylim, $
      title=titleStr, $
      xtitle=textoidl('1-1/r'), $
      ytitle=ytitleStr, $
      charsize=2., /xs, /ys, $
      thick=2.

oplot, 1.d0-1.d0/r, fcont_jomom*r^2, $
      line=2

oplot, 1.d0-1.d0/r3d, fcontr3d_sc*r3d^2, $
      color=ci_blue, $
      psym=3

oplot, 1.d0-1.d0/r3d, fcontr3d_fvm*r3d^2, $
      color=ci_red, $
      psym=3
;
lstr1='1d (ray by ray)'
lstr2='1d (moments)'
lstr3='3d SC'
lstr4='3d FVM'
if(keyword_set(oname)) then begin
   legend, [lstr1, lstr2, lstr3, lstr4], $
           psym=[0,0,3,3], $
           linestyle=[0,2,0,0], $
           color=[ci_black, ci_black,ci_blue,ci_red], $
           textcolor=ci_black, $
           /right_legend
endif else begin
   legend, [lstr1, lstr2, lstr3, lstr4], $
           psym=[0,0,0,0], $
           linestyle=[0,2,3,3], $
           color=[ci_white,ci_white,ci_blue,ci_red], $
           textcolor=ci_white, $
           /right_legend
endelse
;
;-----------------------------------------------------------------------
;
indx=where(mask3d eq 1 or mask3d eq 2 or mask3d eq 3)
ymin=min([min(fcontth3d_sc(indx)), min(fcontth3d_fvm(indx))])
ymax=max([max(fcontth3d_sc(indx)), max(fcontth3d_fvm(indx))])
ylim=[ymin, ymax]
!p.multi=[4,2,3]
plot, [0.,0.], [0.,0.], $
      xrange=xlim, $
      yrange=ylim, $
      title=titleStr, $
      xtitle=xtitleStr, $
      ytitle=textoidl('F_\Theta/I_c'), $
      charsize=2.

oplot, r3d(indx), fcontth3d_sc(indx), $
      color=ci_blue, $
      psym=3
oplot, r3d(indx), fcontth3d_fvm(indx), $
      color=ci_red, $
      psym=3

;
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
!p.multi=[3,2,3]
plot, [0.,0.], [0.,0.], $
      xrange=[0.,1.], $
      yrange=ylim, $
      title=titleStr, $
      xtitle=xtitleStr, $
      ytitle=textoidl('F_\Theta/I_c'), $
      charsize=2.

oplot, 1.d0-1.d0/r3d(indx), fcontth3d_sc(indx), $
      color=ci_blue, $
      psym=3
oplot, 1.d0-1.d0/r3d(indx), fcontth3d_fvm(indx), $
      color=ci_red, $
      psym=3

;
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
;-----------------------------------------------------------------------
;
indx=where(mask3d eq 1 or mask3d eq 2 or mask3d eq 3)
ymin=min([min(fcontphi3d_sc(indx)), min(fcontphi3d_fvm(indx))])
ymax=max([max(fcontphi3d_sc(indx)), max(fcontphi3d_fvm(indx))])
ylim=[ymin, ymax]
!p.multi=[2,2,3]
plot, [0.,0.], [0.,0.], $
      xrange=xlim, $
      yrange=ylim, $
      title=titleStr, $
      xtitle=xtitleStr, $
      ytitle=textoidl('F_\Phi/I_c'), $
      charsize=2.

oplot, r3d(indx), fcontphi3d_sc(indx), $
      color=ci_blue, $
      psym=3
oplot, r3d(indx), fcontphi3d_fvm(indx), $
      color=ci_red, $
      psym=3

;
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
!p.multi=[1,2,3]
plot, [0.,0.], [0.,0.], $
      xrange=[0.,1.], $
      yrange=ylim, $
      title=titleStr, $
      xtitle=xtitleStr, $
      ytitle=textoidl('F_\Phi/I_c'), $
      charsize=2.

oplot, 1.d0-1.d0/r3d(indx), fcontphi3d_sc(indx), $
      color=ci_blue, $
      psym=3
oplot, 1.d0-1.d0/r3d(indx), fcontphi3d_fvm(indx), $
      color=ci_red, $
      psym=3

;
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
if keyword_set(oname) then begin
   device, /close
   set_plot,'x'
endif

!p.multi=0
;
end
