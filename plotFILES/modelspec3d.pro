pro modelspec3d, fname=fname, phic=phic, thetac=thetac, spc=spc, windx=windx, oname=oname, $
                 xlim=xlim, ylim=ylim;, xlim1d=xlim1d
;
;+
; NAME:
;	modelspec3d
;
; PURPOSE:
;	This procedure plots contours on a slice through the 3d grid from
;	output of modelspec.eo. The slice is defined by the angle phic
;       Contours are shown for:
;          3d model: rho, v_r, t
;
;	This procedure plots the radial stratification of the 3d grid from
;	output of modelspec.eo, at a given thetac, phic
;       Plots are shown for:
;          3d model: rho, v_r, v_th, v_phi, t
;
; CALLING SEQUENCE:
;	model3dspec, fname=fname
;
; INPUTS:
;	
; KEYWORD PARAMETERS:
;       thetac: theta coordinate for radial stratification  
;       phic:   phi coordinate for for slice and radial stratification
;       windx:  Set this keyword to an integer, defining the window-index
;       oname:  Set this keyword to a string, defining the output-name
;               (ps-file)
;       xlim:   Set this keyword to an array of length 2, defining the xrange
;       ylim:   Set this keyword to an array of length 2, defining the yrange
;       xlim1d: Set this keyword to an array of length 2, defining the xrange of the radial stratification
;       spc:    Set this keyword (flag) if input model in spherical coordinates
;	help:	Set this keyword (flag) tho show the documentation of this procedure
;
; OUTPUTS:
;	
; EXAMPLE:
;	modelspec3d, fname='modspec_model00.h5'
;       modelspec3d, fname='modspec_model00.h5', thetac=45, phic=45, xlim=[-3.,3.], xlim1d=[1.,10.
;-
;
;-------------------------OUTPUT IF HELP NEEDED-------------------------
;
pdefault=!p
;
if(keyword_set(print_help)) then begin
   doc_library, 'modelspec3d'
   return
endif
;
if(n_elements(fname) eq 0) then begin
   print, 'error in modelspec3d: fname not specified'
   doc_library, 'modelspec3d'
   stop
endif
;
if(not keyword_set(windx)) then windx=0

;thetac, phic specify the theta,phi angles for which radial stratification shall be plotted
if(not keyword_set(phic)) then phic=0.d0
if(not keyword_set(thetac)) then thetac=0.d0
phic=phic*!pi/180
thetac=thetac*!pi/180.
;
;-----------------------------------------------------------------------
;
if(not keyword_set(spc)) then begin
;input file with cartesian coordinates
   read_modelspec3d, fname, teff=teff, xnue0=xnue0, rstar=rstar, vth_fiducial=vth_fiducial, $
                     vmicro=vmicro, vmax=vmax, logg=logg, na=na, $
                     ndxmax=ndxmax, ndymax=ndymax, ndzmax=ndzmax, xcoord=xcoord, ycoord=ycoord, zcoord=zcoord, $
                     opac3d=opac3d, opalbar3d=opalbar3d, t3d=t3d, velx3d=velx3d, vely3d=vely3d, $
                     velz3d=velz3d, scont3d=scont3d, sline3d=sline3d
;
;opacity
   sr=rstar*!rsu
   if(keyword_set(oname)) then oname='modelspec3d_opalbar.ps'
   ctitlestr=textoidl('log(<\chi>)')
   clim=[-8.,2.7]
   opalbar3d=opalbar3d*sr
   opalbar3d(where(opalbar3d eq 0))=1.d-10
   contour_modelspec3d_cac, xcoord, ycoord, zcoord, opalbar3d, phic, xlim=xlim, ylim=ylim, clim=clim, ctitlestr=ctitlestr, windx=windx, oname=oname, /isotropic, /logscale
;
;radial velocity
   velr3d=fltarr(ndxmax,ndymax,ndzmax)
   for i=0, ndxmax-1 do begin
      for j=0, ndymax-1 do begin
         for k=0, ndzmax-1 do begin
            rad=sqrt(xcoord(i)^2+ycoord(j)^2+zcoord(k)^2)
            if(rad ge 1.d0) then begin
               velr3d(i,j,k) = velx3d(i,j,k)*xcoord(i)/rad + vely3d(i,j,k)*ycoord(j)/rad + velz3d(i,j,k)*zcoord(k)/rad
            endif else begin
               velr3d(i,j,k)=0.d0
            endelse
         endfor
      endfor
   endfor
   if(keyword_set(oname)) then oname='modelspec3d_velr.ps'
   ctitlestr=textoidl('v_r')
   contour_modelspec3d_cac, xcoord, ycoord, zcoord, velr3d, phic, xlim=xlim, ylim=ylim, ctitlestr=ctitlestr, windx=windx, oname=oname, /isotropic
;
;temperature
   if(keyword_set(oname)) then oname='modelspec3d_temperature.ps'
   ctitlestr=textoidl('T')
   contour_modelspec3d_cac, xcoord, ycoord, zcoord, t3d, phic, xlim=xlim, ylim=ylim, ctitlestr=ctitlestr, windx=windx, oname=oname, /isotropic, /logscale
;
;plots along radial stratification
   if(keyword_set(oname)) then oname='modelspec3d_radial.ps'
   plot_modelspec3d_cac, xcoord, ycoord, zcoord, thetac, phic, $
                         opac3d, opalbar3d, velx3d, vely3d, velz3d, velr3d, t3d, scont3d, sline3d, xlim=xlim1d, windx=windx, oname=oname
;
;-----------------------------------------------------------------------
;   
endif else begin  
;input file with spherical coordinates   
;
   read_modelspec3d_spc, fname, teff=teff, xnue0=xnue0, rstar=rstar, vth_fiducial=vth_fiducial, $
                         vmicro=vmicro, vmax=vmax, vmin=vmin, xlogg=xlogg, na=na, xic1=xic1, $
                         nr=nr, ntheta=ntheta, nphi=nphi, radius=radius, theta=theta, phi=phi, $
                         opac3d=opac3d, opalbar3d=opalbar3d, t3d=t3d, velx3d=velx3d, vely3d=vely3d, $
                         velz3d=velz3d, scont3d=scont3d, sline3d=sline3d
;
;opacity
   sr=rstar*!rsu
   if(keyword_set(oname)) then oname='modelspec3d_opalbar.ps'
   ctitlestr=textoidl('log(<\chi>)')
   clim=[-8.,2.7]
   opalbar3d=opalbar3d*sr
   contour_modelspec3d_spc, radius, theta, phi, opalbar3d, phic, $
                            xlim=xlim, ylim=ylim, clim=clim, ctitlestr=ctitlestr, $
                            windx=windx, oname=oname, /isotropic, /logscale
;radial velocity
   velr3d=fltarr(nr,ntheta,nphi)*0.d0
   velth3d=fltarr(nr,ntheta,nphi)*0.d0   
   for i=0, nr-1 do begin
      for j=0, ntheta-1 do begin
         for k=0, nphi-1 do begin
            velr3d(i,j,k) = velx3d(i,j,k)*sin(theta(j))*cos(phi(k)) + vely3d(i,j,k)*sin(theta(j))*sin(phi(k)) + velz3d(i,j,k)*cos(theta(j))
            velth3d(i,j,k) = velx3d(i,j,k)*cos(theta(j))*cos(phi(k)) + vely3d(i,j,k)*cos(theta(j))*sin(phi(k)) - velz3d(i,j,k)*sin(theta(j))
         endfor
      endfor
   endfor
   if(keyword_set(oname)) then oname='modelspec3d_velr.ps'
   ctitlestr=textoidl('v_r')
   contour_modelspec3d_spc, radius, theta, phi, velr3d, phic, $
                            xlim=xlim, ylim=ylim, ctitlestr=ctitlestr, $
                            windx=windx, oname=oname, /isotropic

   if(keyword_set(oname)) then oname='modelspec3d_velth.ps'
   ctitlestr=textoidl('v_\theta')
   contour_modelspec3d_spc, radius, theta, phi, velth3d, phic, $
                            xlim=xlim, ylim=ylim, ctitlestr=ctitlestr, $
                            windx=windx, oname=oname, /isotropic   

;temperature
;   if(keyword_set(oname)) then oname='modelspec3d_temperature.ps'
;   ctitlestr=textoidl('T')
;   contour_modelspec3d_spc, radius, theta, phi, t3d, phic, $
;                            xlim=xlim, ylim=ylim, ctitlestr=ctitlestr, $
;                            windx=windx, oname=oname, /isotropic, /logscale
;
;plots along radial stratification
   if(keyword_set(oname)) then oname='modelspec3d_radial.ps'
   opac3d=opac3d*sr
   plot_modelspec3d_spc, radius, theta, phi, thetac, phic, $
                         opac3d, opalbar3d, velx3d, vely3d, velz3d, velr3d, $
                         t3d, scont3d, sline3d, xlim=xlim, windx=windx, oname=oname
endelse
;
;-----------------------------------------------------------------------
;
!p=pdefault

end
;
;-----------------------------------------------------------------------
;-----------------------------------------------------------------------
;-----------------------------------------------------------------------
;
pro contour_modelspec3d_cac, x, y, z, arr3d, phic, xlim=xlim, ylim=ylim, clim=clim, ctitlestr=ctitlestr, windx=windx, oname=oname, isotropic=isotropic, downscale=downscale, logscale=logscale
;
ndxmax=n_elements(x)
ndymax=n_elements(y)
ndzmax=n_elements(z)
;
nr=501
rmin=1.d0
rmax=max([x(1:ndxmax-1),y(1:ndymax-1),z(1:ndzmax-1)])
grid_log, nr, radius, rmin, rmax
;
ntheta=91
grid_equi, ntheta, theta, 0.,!pi
;
arr2d=fltarr(nr,ntheta)*0.d0
;
if(keyword_set(logscale)) then begin
   arr3d_tmp=alog10(arr3d)
endif else begin
   arr3d_tmp=arr3d
endelse
;
for i=0, nr-1 do begin
   for j=0, ntheta-1 do begin
      xcoord = radius(i)*sin(theta(j))*cos(phic)
      ycoord = radius(i)*sin(theta(j))*sin(phic)
      zcoord = radius(i)*cos(theta(j))
   
      find_indx, xcoord, x, ndxmax, ii, iim1
      find_indx, ycoord, y, ndymax, jj, jjm1
      find_indx, zcoord, z, ndzmax, kk, kkm1

      coeff3d_8p_lin, x(iim1), x(ii), y(jjm1), y(jj), z(kkm1), z(kk), xcoord, ycoord, zcoord, $
                       acoeff, bcoeff, ccoeff, dcoeff, ecoeff, fcoeff, gcoeff, hcoeff
;
      arr2d(i,j) = acoeff*arr3d_tmp(iim1,jjm1,kkm1) + bcoeff*arr3d_tmp(ii,jjm1,kkm1) + $
                   ccoeff*arr3d_tmp(iim1,jj,kkm1) + dcoeff*arr3d_tmp(ii,jj,kkm1) + $
                   ecoeff*arr3d_tmp(iim1,jjm1,kk) + fcoeff*arr3d_tmp(ii,jjm1,kk) + $
                   gcoeff*arr3d_tmp(iim1,jj,kk) + hcoeff*arr3d_tmp(ii,jj,kk)
      if(keyword_set(logscale)) then arr2d(i,j)=10.d0^arr2d(i,j)
   endfor
endfor
;
if(keyword_set(downscale)) then begin
;downscale array to make plot smaller
   nr_new=81
   ntheta_new=82
   downsize_arr2d, radius, theta, arr2d, nr_new, ntheta_new, radius_new, theta_new, arr2d_new
endif else begin
   theta_new=theta
   radius_new=radius
   arr2d_new=arr2d
   nr_new=nr
   ntheta_new=ntheta
endelse
;
x_coord=fltarr(nr_new,ntheta_new)*0.d0
y_coord=fltarr(nr_new,ntheta_new)*0.d0
;
;-----------------------------------------------------------------------
;
if(keyword_set(logscale)) then begin
   indx=where(arr2d_new le 0.d0)
   if(n_elements(indx) ne 1) then begin
      arr2d_new(indx)=1.d10
      cmin=min(arr2d_new)
      arr2d_new(indx)=1.d-3*cmin
   endif else begin
      cmin=min(arr2d_new)
   endelse
   cmin=alog10(cmin)
   arr2d_new=alog10(arr2d_new)
endif   
if(not keyword_set(windx)) then windx=0
if(not keyword_set(xlim)) then xlim=[0.d0, 3.5d0]
if(not keyword_set(ylim)) then ylim=[-1.2d0,1.2d0]
if(not keyword_set(clim)) then begin
   if(not keyword_set(logscale)) then begin
      cmin=min(arr2d_new)
   endif
   cmax=max(arr2d_new)
   clim=[cmin,cmax]
endif
;
;define x and y coordinates
for i=0, nr_new-1 do begin
   for j=0, ntheta_new-1 do begin
      x_coord(i,j) = radius_new(i)*sin(theta_new(j))
      y_coord(i,j) = radius_new(i)*cos(theta_new(j))
   endfor
endfor
;
;set titles
xtitleStr=textoidl('x')
ytitleStr=textoidl('z')
titlestr=textoidl('at \phi= ') + string(phic*180./!pi, format='(f9.5)')
;
;get colors
get_contour_colors2, clim, ncolors, nlevels_iso, bottom, $
                    levels_final, levels_iso, c_colors, $
                     cb_indx, cb_tickmark_name
;
;-----------------------------------------------------------------------
;
if(keyword_set(oname)) then begin
   print, "writing output to: ", oname
   set_plot,'ps'
   device,file=oname, color=1, bits_per_pixel=8
endif else begin
   window, windx
   device, decomposed=0
   windx=windx+1
endelse
;
;make contour plot
contourplots_single, arr2d_new, x_coord, y_coord, $
              ncolors, nlevels_iso, bottom, levels_final, levels_iso, $
              c_colors, $
              cb_indx, cb_tickmark_name, $ 
              titleStr=titleStr, xtitleStr=xtitleStr, ytitleStr=ytitleStr, $
              xlim=xlim, ylim=ylim, $
                     ctable=13, ctitleStr=ctitleStr, /background2, isotropic=isotropic
;
;overplot the grid
loadct, 0
for i=0, ndxmax-1 do begin
   oplot, [x(i),x(i)], [ylim(0),ylim(1)]
endfor
for i=0, ndymax-1 do begin
   oplot, [xlim(0),xlim(1)], [y(i),y(i)]
endfor
;
if keyword_set(oname) then begin
   device, /close
   set_plot,'x'
endif
;
end
;
;-----------------------------------------------------------------------
;-----------------------------------------------------------------------
;-----------------------------------------------------------------------
;
pro plot_modelspec3d_cac, x, y, z, thetac, phic, $
                          opac3d, opalbar3d, velx3d, vely3d, velz3d, velr3d, t3d, scont3d, sline3d, xlim=xlim, windx=windx, oname=oname, vmax=vmax
;
if(not keyword_set(xlim)) then xlim=[0.9,10.]  
;-----------------------------------------------------------------------
;
ndxmax=n_elements(x)
ndymax=n_elements(y)
ndzmax=n_elements(z)
;
;create radial grid
nr=501
rmin=1.d0
rmax=max([x,y,z])
grid_log, nr, radius, rmin, rmax
;
;----------------------interpolate onto radial grid---------------------
;
opac1d=fltarr(nr)*0.d0
opalbar1d=fltarr(nr)*0.d0
t1d=fltarr(nr)*0.d0
velx1d=fltarr(nr)*0.d0
vely1d=fltarr(nr)*0.d0
velz1d=fltarr(nr)*0.d0
velr1d=fltarr(nr)*0.d0
scont1d=fltarr(nr)*0.d0
sline1d=fltarr(nr)*0.d0
;
opalbar3d_tmp = alog10(opalbar3d)
;
for i=0, nr-1 do begin
   xcoord = radius(i)*sin(thetac)*cos(phic)
   ycoord = radius(i)*sin(thetac)*sin(phic)
   zcoord = radius(i)*cos(thetac)
   
   find_indx, xcoord, x, ndxmax, ii, iim1
   find_indx, ycoord, y, ndymax, jj, jjm1
   find_indx, zcoord, z, ndzmax, kk, kkm1

   coeff3d_8p_lin, x(iim1), x(ii), y(jjm1), y(jj), z(kkm1), z(kk), xcoord, ycoord, zcoord, $
                    acoeff, bcoeff, ccoeff, dcoeff, ecoeff, fcoeff, gcoeff, hcoeff

   opac1d(i) = acoeff*opac3d(iim1,jjm1,kkm1) + bcoeff*opac3d(ii,jjm1,kkm1) + $
               ccoeff*opac3d(iim1,jj,kkm1) + dcoeff*opac3d(ii,jj,kkm1) + $
               ecoeff*opac3d(iim1,jjm1,kk) + fcoeff*opac3d(ii,jjm1,kk) + $
               gcoeff*opac3d(iim1,jj,kk) + hcoeff*opac3d(ii,jj,kk)

   opalbar1d(i) = acoeff*opalbar3d_tmp(iim1,jjm1,kkm1) + bcoeff*opalbar3d_tmp(ii,jjm1,kkm1) + $
               ccoeff*opalbar3d_tmp(iim1,jj,kkm1) + dcoeff*opalbar3d_tmp(ii,jj,kkm1) + $
               ecoeff*opalbar3d_tmp(iim1,jjm1,kk) + fcoeff*opalbar3d_tmp(ii,jjm1,kk) + $
                  gcoeff*opalbar3d_tmp(iim1,jj,kk) + hcoeff*opalbar3d_tmp(ii,jj,kk)
   opalbar1d(i)=10d0^opalbar1d(i)

   t1d(i) = acoeff*t3d(iim1,jjm1,kkm1) + bcoeff*t3d(ii,jjm1,kkm1) + $
               ccoeff*t3d(iim1,jj,kkm1) + dcoeff*t3d(ii,jj,kkm1) + $
               ecoeff*t3d(iim1,jjm1,kk) + fcoeff*t3d(ii,jjm1,kk) + $
               gcoeff*t3d(iim1,jj,kk) + hcoeff*t3d(ii,jj,kk)

   velx1d(i) = acoeff*velx3d(iim1,jjm1,kkm1) + bcoeff*velx3d(ii,jjm1,kkm1) + $
               ccoeff*velx3d(iim1,jj,kkm1) + dcoeff*velx3d(ii,jj,kkm1) + $
               ecoeff*velx3d(iim1,jjm1,kk) + fcoeff*velx3d(ii,jjm1,kk) + $
               gcoeff*velx3d(iim1,jj,kk) + hcoeff*velx3d(ii,jj,kk)

   vely1d(i) = acoeff*vely3d(iim1,jjm1,kkm1) + bcoeff*vely3d(ii,jjm1,kkm1) + $
               ccoeff*vely3d(iim1,jj,kkm1) + dcoeff*vely3d(ii,jj,kkm1) + $
               ecoeff*vely3d(iim1,jjm1,kk) + fcoeff*vely3d(ii,jjm1,kk) + $
               gcoeff*vely3d(iim1,jj,kk) + hcoeff*vely3d(ii,jj,kk)

   velz1d(i) = acoeff*velz3d(iim1,jjm1,kkm1) + bcoeff*velz3d(ii,jjm1,kkm1) + $
               ccoeff*velz3d(iim1,jj,kkm1) + dcoeff*velz3d(ii,jj,kkm1) + $
               ecoeff*velz3d(iim1,jjm1,kk) + fcoeff*velz3d(ii,jjm1,kk) + $
               gcoeff*velz3d(iim1,jj,kk) + hcoeff*velz3d(ii,jj,kk)

   velr1d(i) = acoeff*velr3d(iim1,jjm1,kkm1) + bcoeff*velr3d(ii,jjm1,kkm1) + $
               ccoeff*velr3d(iim1,jj,kkm1) + dcoeff*velr3d(ii,jj,kkm1) + $
               ecoeff*velr3d(iim1,jjm1,kk) + fcoeff*velr3d(ii,jjm1,kk) + $
               gcoeff*velr3d(iim1,jj,kk) + hcoeff*velr3d(ii,jj,kk)

   scont1d(i) = acoeff*scont3d(iim1,jjm1,kkm1) + bcoeff*scont3d(ii,jjm1,kkm1) + $
               ccoeff*scont3d(iim1,jj,kkm1) + dcoeff*scont3d(ii,jj,kkm1) + $
               ecoeff*scont3d(iim1,jjm1,kk) + fcoeff*scont3d(ii,jjm1,kk) + $
               gcoeff*scont3d(iim1,jj,kk) + hcoeff*scont3d(ii,jj,kk)

   sline1d(i) = acoeff*sline3d(iim1,jjm1,kkm1) + bcoeff*sline3d(ii,jjm1,kkm1) + $
               ccoeff*sline3d(iim1,jj,kkm1) + dcoeff*sline3d(ii,jj,kkm1) + $
               ecoeff*sline3d(iim1,jjm1,kk) + fcoeff*sline3d(ii,jjm1,kk) + $
                gcoeff*sline3d(iim1,jj,kk) + hcoeff*sline3d(ii,jj,kk)

endfor
;print, opac1d
;
;-----------------------------------------------------------------------
;
if(not keyword_set(windx)) then windx=0
window, windx, ysize=1200
windx=windx+1
device, decomposed=0
loadct_lev, ci_black, ci_blue, ci_red, ci_cyan, ci_green, ci_magenta, ci_yellow, ci_white
;
!p.multi=[6,2,3]
plot, radius, opac1d, /ylog, $
      xrange=xlim, /xs, $
      xtitle=textoidl('r'), $
      ytitle=textoidl('\chi_c'), $
      charsize=4.

!p.multi=[5,2,3]
plot, radius, opalbar1d, /ylog, $
      xrange=xlim, /xs, $
      xtitle=textoidl('r'), $
      ytitle=textoidl('<\chi_L>'), $
      charsize=4.

!p.multi=[4,2,3]
plot, radius, t1d, $
      xrange=xlim, /xs, $
      xtitle=textoidl('r'), $
      ytitle=textoidl('T_e'), $
      charsize=4.

!p.multi=[3,2,3]
if(not keyword_set(vmax)) then vmax=max([velx1d,vely1d,velz1d,velr1d])
plot, [0.,0.], [0.,0.], $
      xrange=xlim, /xs, $
      yrange=[0.,1.], /ys, $
      xtitle=textoidl('r'), $
      ytitle=textoidl('v_i/v_\infty'), $
      charsize=4.
oplot, radius, velr1d/vmax
oplot, radius, velx1d/vmax, color=ci_blue
oplot, radius, vely1d/vmax, color=ci_red
oplot, radius, velz1d/vmax, color=ci_green
legend, ['v_r', 'v_x','v_y', 'v_z'], $
         line=[0,0,0,0], $
         color=[ci_white,ci_blue,ci_red,ci_green], $
        charsize=1.5, $
        /bottom, $
        /right_legend

!p.multi=[2,2,3]
plot, radius, scont1d, $
      xrange=xlim, /xs, $
      xtitle=textoidl('r'), $
      ytitle=textoidl('S_c'), $
      charsize=4.

!p.multi=[1,2,3]
plot, radius, sline1d, /ylog, $
      xrange=xlim, /xs, $
      xtitle=textoidl('r'), $
      ytitle=textoidl('S_L'), $
      charsize=4.


end
;
;-----------------------------------------------------------------------
;-----------------------------------------------------------------------
;-----------------------------------------------------------------------
;
pro contour_modelspec3d_spc, radius, theta, phi, arr3d, phic, xlim=xlim, ylim=ylim, clim=clim, ctitlestr=ctitlestr, windx=windx, oname=oname, isotropic=isotropic, downscale=downscale, logscale=logscale
;
nr=n_elements(radius)
ntheta=n_elements(theta)
nphi=n_elements(phi)
;
arr2d=fltarr(nr,ntheta)*0.d0
;
if(keyword_set(logscale)) then begin
   arr3d_tmp=alog10(arr3d)
endif else begin
   arr3d_tmp=arr3d
endelse
;
find_indx, phic, phi, nphi, kk, kkm1
;
for i=0, nr-1 do begin
   for j=0, ntheta-1 do begin
      arr2d(i,j) = interpol1d_2p_lin(arr3d_tmp(i,j,kkm1),arr3d_tmp(i,j,kk),phi(kkm1),phi(kk),phic)
      if(keyword_set(logscale)) then arr2d(i,j)=10.d0^arr2d(i,j)
   endfor
endfor
;
if(keyword_set(downscale)) then begin
;downscale array to make plot smaller
   nr_new=81
   ntheta_new=82
   downsize_arr2d, radius, theta, arr2d, nr_new, ntheta_new, radius_new, theta_new, arr2d_new
endif else begin
   theta_new=theta
   radius_new=radius
   arr2d_new=arr2d
   nr_new=nr
   ntheta_new=ntheta
endelse
;
x_coord=fltarr(nr_new,ntheta_new)*0.d0
y_coord=fltarr(nr_new,ntheta_new)*0.d0
;
;-----------------------------------------------------------------------
;
if(keyword_set(logscale)) then begin
   indx=where(arr2d_new le 0.d0)
   if(n_elements(indx) ne 1) then begin
      arr2d_new(indx)=1.d10
      cmin=min(arr2d_new)
      arr2d_new(indx)=1.d-3*cmin
   endif else begin
      cmin=min(arr2d_new)
   endelse
   cmin=alog10(cmin)
   arr2d_new=alog10(arr2d_new)
endif   
if(not keyword_set(windx)) then windx=0
if(not keyword_set(xlim)) then xlim=[0.d0, 3.5d0]
if(not keyword_set(ylim)) then ylim=[-1.2d0,1.2d0]
if(not keyword_set(clim)) then begin
   if(not keyword_set(logscale)) then begin
      cmin=min(arr2d_new)
   endif
   cmax=max(arr2d_new)
   clim=[cmin,cmax]
endif
;
;define x and y coordinates
for i=0, nr_new-1 do begin
   for j=0, ntheta_new-1 do begin
      x_coord(i,j) = radius_new(i)*sin(theta_new(j))
      y_coord(i,j) = radius_new(i)*cos(theta_new(j))
   endfor
endfor
;
;set titles
xtitleStr=textoidl('x')
ytitleStr=textoidl('z')
titlestr=textoidl('at \phi= ') + string(phic*180./!pi, format='(f9.5)')
;
;get colors
get_contour_colors2, clim, ncolors, nlevels_iso, bottom, $
                    levels_final, levels_iso, c_colors, $
                     cb_indx, cb_tickmark_name
;
;-----------------------------------------------------------------------
;
if(keyword_set(oname)) then begin
   print, "writing output to: ", oname
   set_plot,'ps'
   device,file=oname, color=1, bits_per_pixel=8
endif else begin
   window, windx
   device, decomposed=0
   windx=windx+1
endelse
;
;make contour plot
contourplots_single, arr2d_new, x_coord, y_coord, $
              ncolors, nlevels_iso, bottom, levels_final, levels_iso, $
              c_colors, $
              cb_indx, cb_tickmark_name, $ 
              titleStr=titleStr, xtitleStr=xtitleStr, ytitleStr=ytitleStr, $
              xlim=xlim, ylim=ylim, $
                     ctable=13, ctitleStr=ctitleStr, /background2, isotropic=isotropic
;
;overplot the grid
loadct, 0
for i=0, ntheta_new-1 do begin
   oplot, [radius_new(0)*sin(theta_new(i)),radius_new(nr_new-1)*sin(theta_new(i))], $
          [radius_new(0)*cos(theta_new(i)),radius_new(nr_new-1)*cos(theta_new(i))]
endfor
;
if keyword_set(oname) then begin
   device, /close
   set_plot,'x'
endif
;
end
;
;-----------------------------------------------------------------------
;-----------------------------------------------------------------------
;-----------------------------------------------------------------------
;
pro plot_modelspec3d_spc, radius, theta, phi, thetac, phic, $
                          opac3d, opalbar3d, velx3d, vely3d, velz3d, velr3d, $
                          t3d, scont3d, sline3d, xlim=xlim, windx=windx, oname=oname, vmax=vmax
;
;----------------------interpolate onto radial grid---------------------
;
;
nr=n_elements(radius)
ntheta=n_elements(theta)
nphi=n_elements(phi)
;
opac1d=fltarr(nr)*0.d0
opalbar1d=fltarr(nr)*0.d0
t1d=fltarr(nr)*0.d0
velx1d=fltarr(nr)*0.d0
vely1d=fltarr(nr)*0.d0
velz1d=fltarr(nr)*0.d0
velr1d=fltarr(nr)*0.d0
scont1d=fltarr(nr)*0.d0
sline1d=fltarr(nr)*0.d0

find_indx, thetac, theta, ntheta, ii, iim1
find_indx, phic, phi, nphi, jj, jjm1

opalbar3d_tmp=alog10(opalbar3d)
;
for i=0, nr-1 do begin
   opac1d(i) = interpol2d_4p_lin(opac3d(i,iim1,jjm1), opac3d(i,ii,jjm1), opac3d(i,iim1,jj), opac3d(i,ii,jj), $
                                 theta(iim1),theta(ii),phi(jjm1),phi(jj),thetac,phic)
   opalbar1d(i) = interpol2d_4p_lin(opalbar3d_tmp(i,iim1,jjm1), opalbar3d_tmp(i,ii,jjm1), opalbar3d_tmp(i,iim1,jj), opalbar3d_tmp(i,ii,jj), $
                                    theta(iim1),theta(ii),phi(jjm1),phi(jj),thetac,phic)
   opalbar1d(i) = 10.d0^opalbar1d(i)
   t1d(i) = interpol2d_4p_lin(t3d(i,iim1,jjm1), t3d(i,ii,jjm1), t3d(i,iim1,jj), t3d(i,ii,jj), $
                                 theta(iim1),theta(ii),phi(jjm1),phi(jj),thetac,phic)
   velx1d(i) = interpol2d_4p_lin(velx3d(i,iim1,jjm1), velx3d(i,ii,jjm1), velx3d(i,iim1,jj), velx3d(i,ii,jj), $
                                 theta(iim1),theta(ii),phi(jjm1),phi(jj),thetac,phic)
   vely1d(i) = interpol2d_4p_lin(vely3d(i,iim1,jjm1), vely3d(i,ii,jjm1), vely3d(i,iim1,jj), vely3d(i,ii,jj), $
                                 theta(iim1),theta(ii),phi(jjm1),phi(jj),thetac,phic)   
   velz1d(i) = interpol2d_4p_lin(velz3d(i,iim1,jjm1), velz3d(i,ii,jjm1), velz3d(i,iim1,jj), velz3d(i,ii,jj), $
                                 theta(iim1),theta(ii),phi(jjm1),phi(jj),thetac,phic)
   velr1d(i) = interpol2d_4p_lin(velr3d(i,iim1,jjm1), velr3d(i,ii,jjm1), velr3d(i,iim1,jj), velr3d(i,ii,jj), $
                                 theta(iim1),theta(ii),phi(jjm1),phi(jj),thetac,phic)
   scont1d(i) = interpol2d_4p_lin(scont3d(i,iim1,jjm1), scont3d(i,ii,jjm1), scont3d(i,iim1,jj), scont3d(i,ii,jj), $
                                 theta(iim1),theta(ii),phi(jjm1),phi(jj),thetac,phic)
   sline1d(i) = interpol2d_4p_lin(sline3d(i,iim1,jjm1), sline3d(i,ii,jjm1), sline3d(i,iim1,jj), sline3d(i,ii,jj), $
                                 theta(iim1),theta(ii),phi(jjm1),phi(jj),thetac,phic)   
endfor
;
;-----------------------------------------------------------------------
;
if(not keyword_set(windx)) then windx=0
window, windx, ysize=1200
windx=windx+1
device, decomposed=0
loadct_lev, ci_black, ci_blue, ci_red, ci_cyan, ci_green, ci_magenta, ci_yellow, ci_white
;
!p.multi=[6,2,3]
plot, radius, opac1d, /ylog, $
      xrange=[0.9,10.], /xs, $
      xtitle=textoidl('r'), $
      ytitle=textoidl('\chi_c'), $
      charsize=4.

!p.multi=[5,2,3]
plot, radius, opalbar1d, /ylog, $
      xrange=[0.9,10.], /xs, $
      xtitle=textoidl('r'), $
      ytitle=textoidl('<\chi_L>'), $
      charsize=4.

!p.multi=[4,2,3]
plot, radius, t1d, $
      xrange=[0.9,10.], /xs, $
      xtitle=textoidl('r'), $
      ytitle=textoidl('T_e'), $
      charsize=4.

!p.multi=[3,2,3]
if(not keyword_set(vmax)) then vmax=max([velx1d,vely1d,velz1d,velr1d])
plot, [0.,0.], [0.,0.], $
      xrange=[0.9,10.], /xs, $
      yrange=[0.,1.], /ys, $
      xtitle=textoidl('r'), $
      ytitle=textoidl('v_i/v_\infty'), $
      charsize=4.
oplot, radius, velr1d/vmax
oplot, radius, velx1d/vmax, color=ci_blue
oplot, radius, vely1d/vmax, color=ci_red
oplot, radius, velz1d/vmax, color=ci_green
legend, ['v_r','v_x','v_y', 'v_z'], $
         line=[0,0,0,0], $
         color=[ci_white,ci_blue,ci_red,ci_green], $
        charsize=1.5, $
        /bottom, $
        /right_legend

!p.multi=[2,2,3]
plot, radius, scont1d, $
      xrange=[0.9,10.], /xs, $
      xtitle=textoidl('r'), $
      ytitle=textoidl('S_c'), $
      charsize=4.

!p.multi=[1,2,3]
plot, radius, sline1d, /ylog, $
      xrange=[0.9,10.], /xs, $
      xtitle=textoidl('r'), $
      ytitle=textoidl('S_L'), $
      charsize=4.


end
