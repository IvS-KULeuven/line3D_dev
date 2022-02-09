pro contour_slice2d, x_slice, y_slice, arr2d, velx2d, vely2d, logscale=logscale, xlim=xlim, ylim=ylim, $
                     clim=clim, ovel=ovel, isoc=isoc, titlestr=titlestr
;
;----------------------------define range-------------------------------
;
if(keyword_set(xlim)) then begin
   xmin=min(xlim)
   xmax=max(xlim)
   rangex=xlim
endif else begin
   xmin=min(x_slice)
   xmax=max(x_slice)
   rangex=[xmin, xmax]
endelse
;
if(keyword_set(ylim)) then begin
   ymin=min(ylim)
   ymax=max(ylim)
   rangey=ylim
endif else begin
   ymin=min(y_slice)
   ymax=max(y_slice)
   rangey=[ymin, ymax]
endelse
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
   print, 'contour_slcie2d: keyword ovel to be implemented properly'
;   velovect, velx2d, vely2d, vx_slice, vy_slice, $
;              xrange=xlim, $
;              yrange=ylim, $
;              length=1., $
;              /overplot
endif
;
;
;
end
