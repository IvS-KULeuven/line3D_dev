pro model3d, fname=fname, thetac=thetac, phic=phic, windx=windx, oname=oname, $
             xlim=xlim, ylim=ylim;, xlim1d=xlim1d
;
;+
; NAME:
;	model3d
;
; PURPOSE:
;	This procedure plots contours on a slice through the 3d grid from
;	output of model.eo. The slice is defined by the angle phic
;       Contours are shown for:
;          3d model: rho, v_r, t
;
;	This procedure plots the radial stratification of the 3d grid from
;	output of model.eo, at a given thetac, phic
;       Plots are shown for:
;          3d model: rho, v_r, v_th, v_phi, t
;
; CALLING SEQUENCE:
;	model3d, fname=fname
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
;	help:	Set this keyword (flag) tho show the documentation of this procedure
;
; OUTPUTS:
;	
; EXAMPLE:
;	model3d, fname='model3d.h5'
;       model3d, fname='model3d.h5', thetac=45, phic=45, xlim=[-3.,3.], xlim1d=[1.,10.
;-
;
;-------------------------OUTPUT IF HELP NEEDED-------------------------
;
pdefault=!p
;
if(keyword_set(print_help)) then begin
   doc_library, 'model3d'
   return
endif
;
if(n_elements(fname) eq 0) then begin
   print, 'error in model3d: fname not specified'
   doc_library, 'model3d'
   stop
endif
;
if(not keyword_set(windx)) then windx=0
;
;thetac, phic specify the theta,phi angles for which radial stratification shall be plotted
if(not keyword_set(phic)) then phic=0.d0
if(not keyword_set(thetac)) then thetac=0.d0
phic=phic*!pi/180
thetac=thetac*!pi/180.
;
;-----------------------------------------------------------------------
;
read_model3d, fname, nr=nr, ntheta=ntheta, nphi=nphi, radius=radius, theta=theta, phi=phi, rho3d=rho3d, $
              velr3d=velr3d, velth3d=velth3d, velphi3d=velphi3d, t3d=t3d, $
              vth3d=vth3d
;
;r in r_star
radius=radius/radius(0)
;
;
;-----------------------------------------------------------------------
;
if(keyword_set(oname)) then oname='model3d_rho.ps'
ctitlestr=textoidl('log(\rho)')
indx=where(rho3d gt 0.)
rho_min = min(rho3d(indx))
indx=where(rho3d le 0.)
rho3d(indx) = rho_min
contour_model3d_spc, radius, theta, phi, alog10(rho3d), phic, xlim=xlim, ylim=ylim, ctitlestr=ctitlestr, windx=windx, oname=oname, /isotropic
;
if(keyword_set(oname)) then oname='model3d_velr.ps'
ctitlestr=textoidl('v_r')
contour_model3d_spc, radius, theta, phi, velr3d, phic, xlim=xlim, ylim=ylim, ctitlestr=ctitlestr, windx=windx, oname=oname
;
if(keyword_set(oname)) then oname='model3d_temperature.ps'
ctitlestr=textoidl('T')
contour_model3d_spc, radius, theta, phi, alog10(t3d), phic, xlim=xlim, ylim=ylim, ctitlestr=ctitlestr, windx=windx, oname=oname
;
;-----------------------------------------------------------------------
;;
if(keyword_set(oname)) then oname='model3d_radial.ps'
plot_model3d_spc, radius, theta, phi, thetac, phic, $
                  alog10(rho3d), velr3d, velth3d, velphi3d, t3d, xlim=xlim1d, windx=windx, oname=oname
;
!p=pdefault
end
;
;-----------------------------------------------------------------------
;-----------------------------------------------------------------------
;-----------------------------------------------------------------------
;
pro contour_model3d_spc, radius, theta, phi, arr3d, phic, xlim=xlim, ylim=ylim, clim=clim, ctitlestr=ctitlestr, windx=windx, oname=oname, isotropic=isotropic, downscale=downscale
;
;plot the 2d plane of a 3d array in spherical coordinates at phic
;
nr=n_elements(radius)
ntheta=n_elements(theta)
nphi=n_elements(phi)
;
arr2d=fltarr(nr,ntheta)*0.d0
;
find_indx, phic, phi, nphi, kk, kkm1
;
for i=0, nr-1 do begin
   for j=0, ntheta-1 do begin
      arr2d(i,j) = interpol1d_2p_lin(arr3d(i,j,kkm1),arr3d(i,j,kk),phi(kkm1),phi(kk),phic)
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
if(not keyword_set(windx)) then windx=0
if(not keyword_set(xlim)) then xlim=[0.d0, 3.5d0]
if(not keyword_set(ylim)) then ylim=[-1.2d0,1.2d0]
if(not keyword_set(clim)) then begin
   cmin=min(arr2d_new)
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
for i=0, ntheta_new-1 do begin
   oplot, [radius_new(0)*sin(theta_new(i)),radius_new(nr_new-1)*sin(theta_new(i))], $
          [radius_new(0)*cos(theta_new(i)),radius_new(nr_new-1)*cos(theta_new(i))]
endfor
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
pro plot_model3d_spc, radius, theta, phi, thetac, phic, $
                      rho3d, velr3d, velth3d, velphi3d, t3d, xlim=xlim, windx=windx, oname=oname
;
if(not keyword_set(xlim)) then xlim=[0.9,10.]  
;
;----------------------interpolate onto radial grid---------------------
;
nr=n_elements(radius)
ntheta=n_elements(theta)
nphi=n_elements(phi)
;
rho1d=fltarr(nr)*0.d0
velr1d=fltarr(nr)*0.d0
velth1d=fltarr(nr)*0.d0
velphi1d=fltarr(nr)*0.d0
t1d=fltarr(nr)*0.d0
;
find_indx, thetac, theta, ntheta, ii, iim1
find_indx, phic, phi, nphi, jj, jjm1
;
for i=0, nr-1 do begin
   rho1d(i) = interpol2d_4p_lin(rho3d(i,iim1,jjm1), rho3d(i,ii,jjm1), rho3d(i,iim1,jj), rho3d(i,ii,jj), $
                                 theta(iim1),theta(ii),phi(jjm1),phi(jj),thetac,phic)
   t1d(i) = interpol2d_4p_lin(t3d(i,iim1,jjm1), t3d(i,ii,jjm1), t3d(i,iim1,jj), t3d(i,ii,jj), $
                                  theta(iim1),theta(ii),phi(jjm1),phi(jj),thetac,phic)
   velr1d(i) = interpol2d_4p_lin(velr3d(i,iim1,jjm1), velr3d(i,ii,jjm1), velr3d(i,iim1,jj), velr3d(i,ii,jj), $
                                 theta(iim1),theta(ii),phi(jjm1),phi(jj),thetac,phic)
   velth1d(i) = interpol2d_4p_lin(velth3d(i,iim1,jjm1), velth3d(i,ii,jjm1), velth3d(i,iim1,jj), velth3d(i,ii,jj), $
                                 theta(iim1),theta(ii),phi(jjm1),phi(jj),thetac,phic)   
   velphi1d(i) = interpol2d_4p_lin(velphi3d(i,iim1,jjm1), velphi3d(i,ii,jjm1), velphi3d(i,iim1,jj), velphi3d(i,ii,jj), $
                                 theta(iim1),theta(ii),phi(jjm1),phi(jj),thetac,phic)
endfor
;
;-----------------------------------------------------------------------
;
titlestr=textoidl('at \phi, \theta= ') + string(phic*180./!pi, format='(f9.5)') + ', ' $
                                       + string(thetac*180./!pi, format='(f9.5)')
;
loadct, 0
if(not keyword_set(windx)) then windx=0
window, windx, ysize=1200
windx=windx+1
device, decomposed=0
loadct_lev, ci_black, ci_blue, ci_red, ci_cyan, ci_green, ci_magenta, ci_yellow, ci_white
;
!p.multi=[3,1,3]
plot, radius, rho1d, $
      xrange=xlim, /xs, $
      xtitle=textoidl('r'), $
      ytitle=textoidl('\rho'), $
      charsize=4.

!p.multi=[2,1,3]
plot, radius, t1d, /ylog, $
      xrange=xlim, /xs, $
      xtitle=textoidl('r'), $
      ytitle=textoidl('T'), $
      charsize=4.

vmax=max(sqrt(velr1d^2+velth1d^2+velphi1d^2))
ylim=[0.,1.]
!p.multi=[1,1,3]
plot, [0.,0.], [0.,0.], $
      xrange=xlim, /xs, $
      yrange=ylim, /ys, $
      xtitle=textoidl('r'), $
      ytitle=textoidl('v_i/v_\infty'), $
      charsize=4.
oplot, radius, velr1d/vmax, color=ci_blue
oplot, radius, velth1d/vmax, color=ci_red
oplot, radius, velphi1d/vmax, color=ci_green
legend, ['v_x','v_y', 'v_z'], $
         line=[0,0,0], $
         color=[ci_blue,ci_red,ci_green], $
        charsize=1.5, $
        /bottom, $
        /right_legend


end
