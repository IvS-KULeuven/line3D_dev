pro contour_slice3d_vbin, cs1_r, cs1_theta, cs1_phi, cs1_arr3d, cs1_velx3d, cs1_vely3d, cs1_velz3d, $
                          cs2_r, cs2_theta, cs2_phi, cs2_arr3d, cs2_velx3d, cs2_vely3d, cs2_velz3d, $
                          x01, y01, z01, x02, y02, z02, vx01, vy01, vz01, vx02, vy02, vz02, $
                          unit_length0, unit_length1, unit_length2, $
                          offset=offset, alpha=alpha, gamma=gamma, $
                          logscale=logscale, xlim=xlim, ylim=ylim, clim=clim, $
                          ovel=ovel, titlestr1=titlestr1, titlestr2=titlestr2, ctitlestr=ctitlestr, isoc=isoc
;
;input: cs1_r, cs1_theta, cs1_phi    coordinate system of star 1 (r in unit_length1) 
;       cs1_arr3d                    array to be plotted on slice
;       cs1_velxyz3d                 velocity component of star 1 (in system 1)
;       x01, y01, z01                vector pointing towards origin of system 1 (in unit_length)
;       vx01, vy01, vz01             velocity components of coordinate system 1
;       cs2_r, cs1_theta, cs1_phi    coordinate system of star 2 (r in unit_length2) 
;       cs2_arr3d                    array to be plotted on slice
;       cs2_velxyz3d                 velocity component of star 2 (in system 2)
;       x02, y02, z02                vector pointing towards origin of system 2 (in unit_length)
;       vx02, vy02, vz02             velocity components of coordinate system 2
;
;       alpha, gamma, offset         defining normal vector of slice
;                                    to be plotted (in global coordinate system, offset in unit_length0)
;
;----------------------------define range-------------------------------
;
if(not keyword_set(offset)) then offset=0.d0
if(not keyword_set(alpha)) then alpha=0.d0
if(not keyword_set(gamma)) then gamma=0.d0
if(not keyword_set(clim))  then clim=[0.,0.]
;
nhat=[sin(alpha)*cos(gamma),sin(alpha)*sin(gamma),cos(alpha)]
;
;if(not keyword_set(xlim)) then xlim=[-max([cs1_r,cs2_r]),max([cs1_r,cs2_r])]
;if(not keyword_set(ylim)) then ylim=[-max([cs1_r,cs2_r]),max([cs1_r,cs2_r])]
;
cs1_nr=n_elements(cs1_r)
cs1_ntheta=n_elements(cs1_theta)
cs1_nphi=n_elements(cs1_phi)

cs2_nr=n_elements(cs2_r)
cs2_ntheta=n_elements(cs2_theta)
cs2_nphi=n_elements(cs2_phi)
;
rmin1=cs1_r(0)
rmax1=cs1_r(cs1_nr-1)

rmin2=cs2_r(0)
rmax2=cs2_r(cs2_nr-1)
;
;--------------calculate scalar array on slice for star 1---------------
;
;make the grid of the slice in polar coordinates
   dim_p=301
;note that dim_zeta = 4*uneven - 3 if 45 degree shall be plotted
   dim_zeta=4*27-3
;
;offset for coordinate system 1
   offset1 = offset - nhat(0)*x01 - nhat(1)*y01-nhat(2)*z01
   offset1 = offset1*unit_length0/unit_length1
   slice_cyc2, alpha, gamma, offset1, dim_p, dim_zeta, rmax1, rmin1, $
               cs1_r, cs1_theta, cs1_phi, cs1_nr, cs1_ntheta, cs1_nphi, cs1_arr3d, $
               cs1_xslice, cs1_yslice, cs1_arr2d
;
;---------------calculate velocity components on a slice----------------
;
;if(keyword_set(ovel)) then begin
;   dim_vxslice=11
;   dim_vyslice=11
;   slice_vel, alpha, gamma, offset, dim_vxslice, dim_vyslice, rmax, rmin, $
;              xmax, xmin, x, y, z, ndxmax, ndymax, ndzmax, velx3d, vely3d, velz3d, $
;              vx_slice, vy_slice, velx2d, vely2d
;endif
;
;--------------calculate scalar array on slice for star 2---------------
;
;make the grid of the slice in polar coordinates
   dim_p=301
;note that dim_zeta = 4*uneven - 3 if 45 degree shall be plotted
   dim_zeta=4*27-3
;
;offset for coordinate system 2
   offset2 = offset - nhat(0)*x02 - nhat(1)*y02-nhat(2)*z02
   offset2 = offset2*unit_length0/unit_length2
   slice_cyc2, alpha, gamma, offset2, dim_p, dim_zeta, rmax2, rmin2, $
               cs2_r, cs2_theta, cs2_phi, cs2_nr, cs2_ntheta, cs2_nphi, cs2_arr3d, $
               cs2_xslice, cs2_yslice, cs2_arr2d
;
;use different units for star 2
   cs2_xslice=cs2_xslice*unit_length2/unit_length0
   cs2_yslice=cs2_yslice*unit_length2/unit_length0   
;
;---------------calculate velocity components on a slice----------------
;
;if(keyword_set(ovel)) then begin
;   dim_vxslice=11
;   dim_vyslice=11
;   slice_vel, alpha, gamma, offset, dim_vxslice, dim_vyslice, rmax, rmin, $
;              xmax, xmin, x, y, z, ndxmax, ndymax, ndzmax, velx3d, vely3d, velz3d, $
;              vx_slice, vy_slice, velx2d, vely2d
;endif
   
;
;-------------------contour plots on slice------------------------------
;
;title-strings
xtitlestr=textoidl('x_{slice}')
ytitlestr=textoidl('y_{slice}')
;
;--------------------define color ranges--------------------------------
;
if(keyword_set(clim)) then begin
   clim2=clim
endif else begin
   min1_color = min(cs1_arr2d(where(cs1_arr2d gt 0.)))
   max1_color = max(cs1_arr2d)
   min2_color = min(cs2_arr2d(where(cs2_arr2d gt 0.)))
   max2_color = max(cs2_arr2d)   
   clim2=[min([min1_color,min2_color]),max([max1_color,max2_color])]
   if(keyword_set(logscale)) then clim2=alog10(clim2)
endelse

if(keyword_set(logscale)) then begin
   indx=where(cs1_arr2d gt 0.)
   cs1min=min(cs1_arr2d(indx))
   indx=where(cs1_arr2d eq 0.)
   cs1_arr2d(indx)=cs1min*1.d-8
   cs1_arr2d = alog10(cs1_arr2d)

   indx=where(cs2_arr2d gt 0.)
   cs2min=min(cs2_arr2d(indx))
   indx=where(cs2_arr2d eq 0.)
   cs2_arr2d(indx)=cs2min*1.d-8   
   cs2_arr2d = alog10(cs2_arr2d)
endif
;
;print, 'clim', clim
;print, 'clim2', clim2
;
;
;get colors
get_contour_colors2, clim2, ncolors, nlevels_iso, bottom, $
                     levels_final, levels_iso, c_colors, $
                     cb_indx, cb_tickmark_name
;
loadct, 0
;
;get new contour colors
contourplots_double, cs1_arr2d, cs2_arr2d, $
                     cs1_xslice, cs2_xslice, $
                     cs1_yslice, cs2_yslice, $
                     ncolors, nlevels_iso, bottom, levels_final, $
                     levels_iso, c_colors, cb_indx, cb_tickmark_name, $ 
                     titleStr1=titleStr1, titlestr2=titlestr2, $
                     xtitleStr=xtitleStr, ytitleStr=ytitleStr, $
                     xlim=xlim, ylim=ylim, isoc=isoc, $
                     ctable=13, ctitleStr=ctitleStr, /background1, /isotropic
;
;-----------------combine to global coordinate system-------------------
;
window, 1
device, decomposed=0
;
   calc_transmat2, alpha, gamma, transmat, transmat_inv=transmat_inv
;
   xslice=fltarr(2L*dim_p*dim_zeta)*0.d0
   yslice=fltarr(2L*dim_p*dim_zeta)*0.d0
   zslice=fltarr(2L*dim_p*dim_zeta)*0.d0
   arr2d=fltarr(2L*dim_p*dim_zeta)*0.d0
;
;transform system 1   
   ii=0
   for i=0, dim_p-1 do begin
      for j=0, dim_zeta-1 do begin
         vec_cyc=[cs1_xslice(i,j),cs1_yslice(i,j),offset1] ;in cylindrical system 1
         vec_cac=transmat#vec_cyc                           ;in carthesian system 1
         vec_cac=vec_cac*unit_length1/unit_length0 + [x01,y01,z01] ;in global system
         vec_cyc=transmat_inv#vec_cac                              ;in global cylindrica system
         xslice(ii)=vec_cyc(0)
         yslice(ii)=vec_cyc(1)
         arr2d(ii)=cs1_arr2d(i,j)
         ii=ii+1
      endfor
   endfor
;
   for i=0, dim_p-1 do begin
      for j=0, dim_zeta-1 do begin
         vec_cyc=[cs2_xslice(i,j),cs2_yslice(i,j),offset2]  ;in cylindrical system 2
         vec_cac=transmat#vec_cyc                           ;in carthesian system 2
         vec_cac=vec_cac*unit_length2/unit_length0 + [x02,y02,z02]     ;in global system
         vec_cyc=transmat_inv#vec_cac                                  ;in global cylindrica system
         xslice(ii)=vec_cyc(0)
         yslice(ii)=vec_cyc(1)
         arr2d(ii)=cs2_arr2d(i,j)
         ii=ii+1
      endfor
   endfor  
;
;-------------------contour plots on slice------------------------------
;
;title-strings
xtitlestr=textoidl('x_{slice}')
ytitlestr=textoidl('y_{slice}')
;
;--------------------define color ranges--------------------------------
;
;print, 'clim', clim
;print, 'clim2', clim2
;
;
;get colors
get_contour_colors2, clim2, ncolors, nlevels_iso, bottom, $
                     levels_final, levels_iso, c_colors, $
                     cb_indx, cb_tickmark_name
;
loadct, 0
;
;get new contour colors
contourplots_single, arr2d, xslice, yslice, $
                      ncolors, nlevels_iso, bottom, levels_final, $
                      levels_iso, c_colors, cb_indx, cb_tickmark_name, $ 
                      titleStr=titleStr, $
                      xtitleStr=xtitleStr, ytitleStr=ytitleStr, $
                      xlim=xlim, ylim=ylim, isoc=isoc, $
                      ctable=13, ctitleStr=ctitleStr, /background1, /isotropic, /irregular 



end
