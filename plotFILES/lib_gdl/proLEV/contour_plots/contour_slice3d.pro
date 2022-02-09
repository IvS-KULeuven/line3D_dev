pro contour_slice3d, x, y, z, arr3d, velx3d, vely3d, velz3d, offset=offset, alpha=alpha, gamma=gamma, dvec=dvec, $
                     logscale=logscale, xlim=xlim, ylim=ylim, polar=polar, clim=clim, ovel=ovel, titlestr=titlestr, $
                     ctitlestr=ctitlestr, isoc=isoc, np=np, nzeta=nzeta
;
;----------------------------define range-------------------------------
;
if(keyword_set(xlim)) then begin
   xmin=min(xlim)
   xmax=max(xlim)
   rangex=xlim
endif else begin
   xmin=min([min(x), min(y), min(z)])
   xmax=max([max(x), max(y), max(z)])
   rangex=[xmin, xmax]
endelse
;
if(keyword_set(ylim)) then begin
   ymin=min(ylim)
   ymax=max(ylim)
   rangey=ylim
endif else begin
   ymin=min([min(x), min(y), min(z)])
   ymax=max([max(x), max(y), max(z)])
   rangey=[ymin, ymax]
endelse
;
ndxmax=n_elements(x)
ndymax=n_elements(y)
ndzmax=n_elements(z)
;
rmin=1.d0
rmax=max([x(ndxmax-2), y(ndymax-2), z(ndzmax-2)])
;
;-------------------CALCULATE 3D SCALAR ARRAY ON A SLICE----------------
;
if(keyword_set(polar)) then begin
;make the grid of the slice in polar coordinates
   if(not keyword_set(np)) then np=301
   dim_p=np
;
;note that dim_zeta = 4*uneven - 3 if 45 degree shall be plotted
   if(not keyword_set(nzeta)) then nzeta=4*27-3
   dim_zeta=nzeta
;
   slice_cyc, alpha, gamma, offset, dim_p, dim_zeta, rmax, rmin, $
              x, y, z, ndxmax, ndymax, ndzmax, arr3d, $
              x_slice, y_slice, arr2d
;
endif else begin
;make the grid of the slice in carthesian coordinates
   dim_xslice=301
   dim_yslice=301
;
   slice_cac, alpha, gamma, offset, dim_xslice, dim_yslice, rmax, rmin, $
              x, y, z, ndxmax, ndymax, ndzmax, arr3d, $
              x_slice, y_slice, arr2d
endelse
;
;---------------calculate velocity components on a slice----------------
;
if(keyword_set(ovel)) then begin
   dim_vxslice=11
   dim_vyslice=11
   slice_vel, alpha, gamma, offset, dim_vxslice, dim_vyslice, rmax, rmin, $
              xmax, xmin, x, y, z, ndxmax, ndymax, ndzmax, velx3d, vely3d, velz3d, $
              vx_slice, vy_slice, velx2d, vely2d
endif
;
;-------------------contour plots on slice------------------------------
;
;title-strings
xtitlestr=textoidl('x_{slice}')
ytitlestr=textoidl('y_{slice}')
;
;--------------------define color ranges--------------------------------
;
if(keyword_set(logscale)) then begin
   if(not keyword_set(clim)) then begin
      min_color = min(arr2d(where(arr2d gt 0.)))
      max_color = max(arr2d)
      min_color = alog10(min_color)
      max_color = alog10(max_color)
      clim2=[min_color,max_color]
   endif else begin
      clim2=clim
   endelse
   arr2d = alog10(arr2d)
endif else begin
   if(not keyword_set(clim)) then begin
      min_color = min(arr2d)
      max_color = max(arr2d)
      clim2=[min_color,max_color]
   endif else begin
      clim2=clim
   endelse
endelse

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
contourplots_single2, arr2d, x_slice, y_slice, $
              ncolors, nlevels_iso, bottom, levels_final, levels_iso, $
              c_colors, $
              cb_indx, cb_tickmark_name, $ 
              titleStr=titleStr, xtitleStr=xtitleStr, ytitleStr=ytitleStr, $
              xlim=xlim, ylim=ylim, isoc=isoc, $
              ctable=13, ctitleStr=ctitleStr, /background2, /isotropic
;
theta=findgen(101) * 2.*!pi / 100.
;
xcircle=1.165*cos(theta)
ycircle=1.165*sin(theta)
;
xcircle1=3.4*cos(theta)
ycircle1=3.4*sin(theta)
;
xcircle2=6.1*cos(theta)
ycircle2=6.1*sin(theta)
;
loadct, 0
oplot, xcircle, ycircle
oplot, xcircle1, ycircle1
oplot, xcircle2, ycircle2
;
if(keyword_set(ovel)) then begin
;overplot velocity vectors
   velovect, velx2d, vely2d, vx_slice, vy_slice, $
              xrange=xlim, $
              yrange=ylim, $
              length=1., $
              /overplot
endif
;
;
;
end
