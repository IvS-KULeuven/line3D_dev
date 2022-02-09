pro model2d, dir=dir, windx=windx, oname=oname, clim=clim
;
;plots contours of searchlight along a direction
;
;------------------read all information from hdf5-file------------------
;
if(not keyword_set(dir)) then dir='.'
fname=dir+'/model2d.h5'
;
file_id = h5f_open(fname)
;
   group_id = h5g_open(file_id, 'dimensions')
      att_id=h5a_open_name(group_id, 'nr')
         nr=h5a_read(att_id)
         nr=nr(0)
      h5a_close, att_id
      att_id=h5a_open_name(group_id, 'ntheta')
         ntheta=h5a_read(att_id)
         ntheta=ntheta(0)
      h5a_close, att_id
   h5g_close, group_id
;
   group_id = h5g_open(file_id, 'coordinates')
      dset_id=h5d_open(group_id, 'r')
         r=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'theta')
         theta=h5d_read(dset_id)
      h5d_close, dset_id
   h5g_close, group_id
;
   group_id = h5g_open(file_id, 'model')
      dset_id=h5d_open(group_id, 'rho')
         rho2d=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'velr')
         velr2d=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'velth')
         velth2d=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'velphi')
         velphi2d=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'temperature')
         t2d=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'vth')
         vth2d=h5d_read(dset_id)
      h5d_close, dset_id
   h5g_close, group_id
;
h5f_close, file_id



;r in r_star
r=r/r(0)
;
;-----------------------define range------------------------------------
;
rlog = alog10(r-1.d0)
rlog(0) = -6.d0
if(not keyword_set(xlim)) then xlim=[min(rlog),max(rlog)]
if(not keyword_set(ylim)) then ylim=[0.,1.]
;
;----------------------------title-strings------------------------------
;
;titlestr=textoidl('i (\theta=') + string(acos(mu), format='(f9.5)') + $
;         textoidl(', \phi=') + string(phi, format='(f9.5)') + ')'
titlestr=''
xtitlestr=textoidl('r')
ytitlestr=textoidl('v')
;
;-----------------------------------------------------------------------
;
;plots along radial direction for several theta angles
;
theta0=0.d0
theta1=20.d0
theta2=40.d0
theta3=60.d0
theta4=80.d0
theta5=90.d0
;
find_indx, theta0*!pi/180.d0, theta, ntheta, indx0, ii
find_indx, theta1*!pi/180.d0, theta, ntheta, iim1, indx1
find_indx, theta2*!pi/180.d0, theta, ntheta, iim1, indx2
find_indx, theta3*!pi/180.d0, theta, ntheta, iim1, indx3
find_indx, theta4*!pi/180.d0, theta, ntheta, iim1, indx4
find_indx, theta5*!pi/180.d0, theta, ntheta, indx5, ii
;
if(not keyword_set(windx)) then windx=0
;
if(keyword_set(oname)) then begin
   oname='ps_files/model2d.ps'
   print, "writing output to: ", oname
   set_plot,'ps'
   device,file=oname, xsize=19., ysize=26.7, xoffset=1., yoffset=1. , $
          color=1, bits_per_pixel=8;, decomposed=0
endif else begin
   window, windx, xsize=950, ysize=1200
   device, decomposed=0
   windx=windx+1
endelse
;
loadct, 0
;
loadct_lev, ci_black, ci_blue, ci_red, ci_cyan, ci_green, ci_magenta, ci_yellow, ci_white
;
;--------------v_r(theta)/v_r(theta=0) at different theta---------------
;
!p.multi=[6,2,3]
;
ylim=[0.7,1.1]
;ylim=[0.9,1.1]
plot, rlog, velr2d(*,indx0)/velr2d(*,indx0), $
      yrange=ylim, /ys, $
      xrange=[-5.,2.], /xs, $
      xtitle=textoidl('log(R/R_*-1)'), $
      ytitle=textoidl('v_r(2d) / v_r(1d)'), $
      charsize=2., $
      linestyle=5
oplot, rlog, velr2d(*,indx1)/velr2d(*,indx0), linestyle=4
oplot, rlog, velr2d(*,indx2)/velr2d(*,indx0), linestyle=3
oplot, rlog, velr2d(*,indx3)/velr2d(*,indx0), linestyle=2
oplot, rlog, velr2d(*,indx4)/velr2d(*,indx0), linestyle=1
oplot, rlog, velr2d(*,indx5)/velr2d(*,indx0), linestyle=0

lstr_00 = 'THETA=' + string(theta(indx0)*180.d0/!pi, format='(i3)')
lstr_20 = 'THETA=' + string(theta(indx1)*180.d0/!pi, format='(i3)')
lstr_40 = 'THETA=' + string(theta(indx2)*180.d0/!pi, format='(i3)')
lstr_60 = 'THETA=' + string(theta(indx3)*180.d0/!pi, format='(i3)')
lstr_80 = 'THETA=' + string(theta(indx4)*180.d0/!pi, format='(i3)')
lstr_90 = 'THETA=' + string(theta(indx5)*180.d0/!pi, format='(i3)')
;
loadct_lev, ci_black, ci_blue, ci_red, ci_cyan, ci_green, ci_magenta, ci_yellow, ci_white
;
if(keyword_set(oname)) then begin
   legend, [lstr_90,lstr_80,lstr_60,lstr_40,lstr_20,lstr_00], $
           psym=[0,0,0,0,0,0], $
           linestyle=[0,1,2,3,4,5], $
           color=[ci_black, ci_black,ci_black,ci_black,ci_black,ci_black], $
           textcolor=ci_black, spacing=1.5, charsize=0.8, /bottom
endif else begin
   legend, [lstr_90,lstr_80,lstr_60,lstr_40,lstr_20,lstr_00], $
           psym=[0,0,0,0,0,0], $
           linestyle=[0,1,2,3,4,5], $
           color=[ci_white, ci_white,ci_white,ci_white,ci_white,ci_white], $
           textcolor=ci_white, spacing=1.5, charsize=0.8, /bottom
endelse
;
;------------v_phi(theta)/v_r(theta) for different theta----------------
;
!p.multi=[5,2,3]
;
ylim=[-5.,20.]
;ylim=[-2.,6.]
plot, rlog, velphi2d(*,indx0)/velr2d(*,indx0), $
      yrange=ylim, /ys, $
      xrange=[-5.,2.], /xs, $
      xtitle=textoidl('log(R/R_*-1)'), $
      ytitle=textoidl('v_\phi(2d) / v_r(2d)'), $
      charsize=2., $
      linestyle=5
oplot, rlog, velphi2d(*,indx1)/velr2d(*,indx0), linestyle=4
oplot, rlog, velphi2d(*,indx2)/velr2d(*,indx0), linestyle=3
oplot, rlog, velphi2d(*,indx3)/velr2d(*,indx0), linestyle=2
oplot, rlog, velphi2d(*,indx4)/velr2d(*,indx0), linestyle=1
oplot, rlog, velphi2d(*,indx5)/velr2d(*,indx0), linestyle=0

lstr_00 = 'THETA=' + string(theta(indx0)*180.d0/!pi, format='(i3)')
lstr_20 = 'THETA=' + string(theta(indx1)*180.d0/!pi, format='(i3)')
lstr_40 = 'THETA=' + string(theta(indx2)*180.d0/!pi, format='(i3)')
lstr_60 = 'THETA=' + string(theta(indx3)*180.d0/!pi, format='(i3)')
lstr_80 = 'THETA=' + string(theta(indx4)*180.d0/!pi, format='(i3)')
lstr_90 = 'THETA=' + string(theta(indx5)*180.d0/!pi, format='(i3)')
;
loadct_lev, ci_black, ci_blue, ci_red, ci_cyan, ci_green, ci_magenta, ci_yellow, ci_white
;
if(keyword_set(oname)) then begin
   legend, [lstr_90,lstr_80,lstr_60,lstr_40,lstr_20,lstr_00], $
           psym=[0,0,0,0,0,0], $
           linestyle=[0,1,2,3,4,5], $
           color=[ci_black, ci_black,ci_black,ci_black,ci_black,ci_black], $
           textcolor=ci_black, spacing=1.5, charsize=0.8, /right_legend
endif else begin
   legend, [lstr_90,lstr_80,lstr_60,lstr_40,lstr_20,lstr_00], $
           psym=[0,0,0,0,0,0], $
           linestyle=[0,1,2,3,4,5], $
           color=[ci_white, ci_white,ci_white,ci_white,ci_white,ci_white], $
           textcolor=ci_white, spacing=1.5, charsize=0.8, /right_legend
endelse
;
;------------v_theta(theta)/v_r(theta) for different theta--------------
;
!p.multi=[4,2,3]
;
ylim=[-.2,.8]
;ylim=[-.05,.2]
plot, rlog, velth2d(*,indx0)/velr2d(*,indx0), $
      yrange=ylim, /ys, $
      xrange=[-5.,2.], /xs, $
      xtitle=textoidl('log(R/R_*-1)'), $
      ytitle=textoidl('v_\theta(2d) / v_r(2d)'), $
      charsize=2., $
      linestyle=5
oplot, rlog, velth2d(*,indx1)/velr2d(*,indx0), linestyle=4
oplot, rlog, velth2d(*,indx2)/velr2d(*,indx0), linestyle=3
oplot, rlog, velth2d(*,indx3)/velr2d(*,indx0), linestyle=2
oplot, rlog, velth2d(*,indx4)/velr2d(*,indx0), linestyle=1
oplot, rlog, velth2d(*,indx5)/velr2d(*,indx0), linestyle=0

lstr_00 = 'THETA=' + string(theta(indx0)*180.d0/!pi, format='(i3)')
lstr_20 = 'THETA=' + string(theta(indx1)*180.d0/!pi, format='(i3)')
lstr_40 = 'THETA=' + string(theta(indx2)*180.d0/!pi, format='(i3)')
lstr_60 = 'THETA=' + string(theta(indx3)*180.d0/!pi, format='(i3)')
lstr_80 = 'THETA=' + string(theta(indx4)*180.d0/!pi, format='(i3)')
lstr_90 = 'THETA=' + string(theta(indx5)*180.d0/!pi, format='(i3)')
;
loadct_lev, ci_black, ci_blue, ci_red, ci_cyan, ci_green, ci_magenta, ci_yellow, ci_white
;
if(keyword_set(oname)) then begin
   legend, [lstr_90,lstr_80,lstr_60,lstr_40,lstr_20,lstr_00], $
           psym=[0,0,0,0,0,0], $
           linestyle=[0,1,2,3,4,5], $
           color=[ci_black, ci_black,ci_black,ci_black,ci_black,ci_black], $
           textcolor=ci_black, spacing=1.5, charsize=0.8, /left_legend
endif else begin
   legend, [lstr_90,lstr_80,lstr_60,lstr_40,lstr_20,lstr_00], $
           psym=[0,0,0,0,0,0], $
           linestyle=[0,1,2,3,4,5], $
           color=[ci_white, ci_white,ci_white,ci_white,ci_white,ci_white], $
           textcolor=ci_white, spacing=1.5, charsize=0.8, /left_legend
endelse
;
;------------v_theta(theta)/v_phi(theta) for different theta------------
;
!p.multi=[3,2,3]
;
ylim=[-.2,.8]
;ylim=[-.1,.3]
plot, rlog, velth2d(*,indx1)/velphi2d(*,indx1), $
      yrange=ylim, /ys, $
      xrange=[-5.,2.], /xs, $
      xtitle=textoidl('log(R/R_*-1)'), $
      ytitle=textoidl('v_\theta(2d) / v_\phi(2d)'), $
      charsize=2., $
      linestyle=4
oplot, rlog, velth2d(*,indx2)/velphi2d(*,indx2), linestyle=3
oplot, rlog, velth2d(*,indx3)/velphi2d(*,indx3), linestyle=2
oplot, rlog, velth2d(*,indx4)/velphi2d(*,indx4), linestyle=1
oplot, rlog, velth2d(*,indx5)/velphi2d(*,indx5), linestyle=0

lstr_20 = 'THETA=' + string(theta(indx1)*180.d0/!pi, format='(i3)')
lstr_40 = 'THETA=' + string(theta(indx2)*180.d0/!pi, format='(i3)')
lstr_60 = 'THETA=' + string(theta(indx3)*180.d0/!pi, format='(i3)')
lstr_80 = 'THETA=' + string(theta(indx4)*180.d0/!pi, format='(i3)')
lstr_90 = 'THETA=' + string(theta(indx5)*180.d0/!pi, format='(i3)')
;
loadct_lev, ci_black, ci_blue, ci_red, ci_cyan, ci_green, ci_magenta, ci_yellow, ci_white
;
if(keyword_set(oname)) then begin
   legend, [lstr_90,lstr_80,lstr_60,lstr_40,lstr_20], $
           psym=[0,0,0,0,0], $
           linestyle=[0,1,2,3,4], $
           color=[ci_black,ci_black,ci_black,ci_black,ci_black], $
           textcolor=ci_black, spacing=1.5, charsize=0.8, /left_legend
endif else begin
   legend, [lstr_90,lstr_80,lstr_60,lstr_40,lstr_20], $
           psym=[0,0,0,0,0], $
           linestyle=[0,1,2,3,4], $
           color=[ci_white,ci_white,ci_white,ci_white,ci_white], $
           textcolor=ci_white, spacing=1.5, charsize=0.8, /left_legend
endelse
;
;-----------------------density for different theta---------------------
;
!p.multi=[2,2,3]
;
;note: plot rho2d/rho1d needs to be scaled already in model (fortran) routine
;
ylim=[0.,3.]
;ylim=[.8,1.3]
plot, rlog, rho2d(*,indx0), $
      yrange=ylim, /ys, $
      xrange=[-3.,2.], /xs, $
      xtitle=textoidl('log(R/R_*-1)'), $
      ytitle=textoidl('\rho'), $
      charsize=2., $
      linestyle=5
oplot, rlog, rho2d(*,indx1), linestyle=4
oplot, rlog, rho2d(*,indx2), linestyle=3
oplot, rlog, rho2d(*,indx3), linestyle=2
oplot, rlog, rho2d(*,indx4), linestyle=1
oplot, rlog, rho2d(*,indx5), linestyle=0
;
lstr_00 = 'THETA=' + string(theta(indx0)*180.d0/!pi, format='(i3)')
lstr_20 = 'THETA=' + string(theta(indx1)*180.d0/!pi, format='(i3)')
lstr_40 = 'THETA=' + string(theta(indx2)*180.d0/!pi, format='(i3)')
lstr_60 = 'THETA=' + string(theta(indx3)*180.d0/!pi, format='(i3)')
lstr_80 = 'THETA=' + string(theta(indx4)*180.d0/!pi, format='(i3)')
lstr_90 = 'THETA=' + string(theta(indx5)*180.d0/!pi, format='(i3)')
;
loadct_lev, ci_black, ci_blue, ci_red, ci_cyan, ci_green, ci_magenta, ci_yellow, ci_white
;
if(keyword_set(oname)) then begin
   legend, [lstr_90,lstr_80,lstr_60,lstr_40,lstr_20,lstr_00], $
           psym=[0,0,0,0,0,0], $
           linestyle=[0,1,2,3,4,5], $
           color=[ci_black, ci_black,ci_black,ci_black,ci_black,ci_black], $
           textcolor=ci_black, spacing=1.5, charsize=0.8, /left_legend
endif else begin
   legend, [lstr_90,lstr_80,lstr_60,lstr_40,lstr_20,lstr_00], $
           psym=[0,0,0,0,0,0], $
           linestyle=[0,1,2,3,4,5], $
           color=[ci_white, ci_white,ci_white,ci_white,ci_white,ci_white], $
           textcolor=ci_white, spacing=1.5, charsize=0.8, /left_legend
endelse
;
if keyword_set(oname) then begin
   device, /close
   set_plot,'x'
endif
;
!p.multi=0
;
;***********************************************************************
;
;                           contour plots
;
;***********************************************************************
;
;downscale array to make plot smaller
;nr_new=81
;ntheta_new=82
;
;downsize_arr2d, r, theta, rho2d, nr_new, ntheta_new, radius_new, theta_new, rho_new
;downsize_arr2d, r, theta, velth2d, nr_new, ntheta_new, radius_new, theta_new, velth_new
;downsize_arr2d, r, theta, velr2d, nr_new, ntheta_new, radius_new, theta_new, velr_new
;
;theta=theta_new
;r=radius_new
;rho2d=rho_new
;velr2d=velr_new
;velth2d=velth_new
;nr=nr_new
;ntheta=ntheta_new
;
x_coord=fltarr(nr,ntheta)*0.d0
y_coord=fltarr(nr,ntheta)*0.d0
colors=fltarr(nr,ntheta)*0.d0
;
v_esc=1.d8
;
;-----------------------------------------------------------------------
;
aspect_ratio=1.6
if(keyword_set(oname)) then begin
   print, 'writing output to ', oname
   set_plot, 'ps'
   xsize=17.780d0
   ysize=xsize/aspect_ratio
   device,file=ONAME, $
;          decomposed=0, $
          color=1, $
          BITS_PER_PIXEL=8, $
          ysize=ysize, $
          xsize=xsize
;   device, file=oname, decomposed=0, color=1, bits_per_pixel=8
endif else begin
   window, windx
   device, decomposed=0
   xsize=!d.x_size
   ysize=xsize/aspect_ratio
   window, windx, xsize=xsize, ysize=ysize
   windx=windx+1
endelse
;
;-----------------------------------------------------------------------
;
xlim=[0.d0, 3.5d0]
ylim=[-1.2d0,1.2d0]

xlim=[0.,3.5]
ylim=[-1.4,1.4]

xlim=[-10.,10.]
ylim=[-10.,10.]

xlim=[0.,7.]
ylim=[-2.5,2.5]
ylim=[-7.,7.]
;
;--------------------------prepare velocity vectors---------------------
;
dim_vxslice=15
dim_vyslice=15

xmin=xlim(0)
xmax=0.9*xlim(1)
ymin=0.9*ylim(0)
ymax=0.9*ylim(1)

vx_slice=xmin+findgen(dim_vxslice)*(xmax-xmin)/float(dim_vxslice-1)
vy_slice=ymin+findgen(dim_vyslice)*(ymax-ymin)/float(dim_vyslice-1)

vx2d=fltarr(dim_vxslice,dim_vyslice)*0.d0
vy2d=fltarr(dim_vxslice,dim_vyslice)*0.d0

for i=0, dim_vxslice-1 do begin
   for j=0, dim_vyslice-1 do begin
      rad=sqrt(vx_slice(i)^2 + vy_slice(j)^2)
      th=acos(vy_slice(j)/rad)

      if(rad gt 1.) then begin
         interpol_bilinear, velth2d, r, theta, [rad,th], vtheta
         interpol_bilinear, velr2d, r, theta, [rad,th], vrad
         vx2d(i,j)=(vrad*sin(th)+vtheta*cos(th))/v_esc
         vy2d(i,j)=(vrad*cos(th)-vtheta*sin(th))/v_esc

;set maximum allowed v/v_esc to 2 (to have a nice plot)
         vlim=0.5d0
         vel=sqrt(vx2d(i,j)^2 + vy2d(i,j)^2)
         vnorm=vlim/vel
         if(vel gt vlim) then begin
            vx2d(i,j) = vx2d(i,j)*vnorm
            vy2d(i,j) = vy2d(i,j)*vnorm
;            print, vel, vnorm, vx2d(i,j), vy2d(i,j), sqrt(vx2d(i,j)^2 + vy2d(i,j)^2)
         endif


      endif

   endfor
endfor
;
;-------------------------plots-----------------------------------------
;
print, max(velr2d)/1.d5
;define colorbar for density: from [-1,1] and x, y-coordinates
if(not keyword_set(clim)) then begin
   clim=[-3.,0.]
   clim=[-14.,-11.]
   clim=[min(alog10(rho2d)),max(alog10(rho2d))]
;   clim=[min(velth2d),max(velth2d)]/2329016.8591245087d0
;   clim=[min(velr2d),max(velr2d)]/2329016.8591245087d0
endif
;
for i=0, nr-1 do begin
   for j=0, ntheta-1 do begin
      x_coord(i,j) = r(i)*sin(theta(j))
      y_coord(i,j) = r(i)*cos(theta(j))
   endfor
endfor
;
;set titles
;ctitleStr=textoidl('log(\rho/max(\rho))')
ctitleStr=textoidl('log(\rho)')
xtitleStr=textoidl('x')
ytitleStr=textoidl('z')
;
;get colors
get_contour_colors2, clim, ncolors, nlevels_iso, bottom, $
                    levels_final, levels_iso, c_colors, $
                    cb_indx, cb_tickmark_name
;
;make contour plot
contourplots_single, alog10(rho2d), x_coord, y_coord, $
;contourplots_single, velth2d/2329016.8591245087d0, x_coord, y_coord, $
;contourplots_single, velr2d/2329016.8591245087d0, x_coord, y_coord, $
              ncolors, nlevels_iso, bottom, levels_final, levels_iso, $
              c_colors, $
              cb_indx, cb_tickmark_name, $ 
              titleStr=titleStr, xtitleStr=xtitleStr, ytitleStr=ytitleStr, $
              xlim=xlim, ylim=ylim, $
              ctable=13, ctitleStr=ctitleStr, /background2, /isotropic
;
;overplot magnetic field
;loadct, 0
;plotbfield, nlines=8, ralfven=2.7, /oplot
;
;overplot velocity vectors
velovect, vx2d, vy2d, vx_slice, vy_slice, $
           xrange=xlim, $
           yrange=ylim, $
           length=.85, $
           /overplot
;
if keyword_set(oname) then begin
   device, xsize=17.780d0, ysize=12.7d0
   device, /close
   set_plot,'x'
endif

end
