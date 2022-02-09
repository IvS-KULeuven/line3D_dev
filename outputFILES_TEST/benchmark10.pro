pro benchmark10, dir=dir, windx=WINDX, oname=ONAME
;
;plots contours of searchlight along a direction
;
;------------------READ ALL INFORMATION FROM HDF5-FILE------------------
;
if(not keyword_set(dir)) then dir='.'
fname=dir+'/benchmark10.h5'
;
file_id = h5f_open(fname)
;
   group_id = H5G_OPEN(file_id, 'direction')
      att_id=H5A_OPEN_NAME(group_id, 'n_x')
         n_x=H5A_READ(att_id)
         n_x=n_x(0)
      h5a_close, att_id
      att_id=H5A_OPEN_NAME(group_id, 'n_y')
         n_y=H5A_READ(att_id)
         n_y=n_y(0)
      h5a_close, att_id
      att_id=H5A_OPEN_NAME(group_id, 'n_z')
         n_z=H5A_READ(att_id)
         n_z=n_z(0)
      h5a_close, att_id
   h5g_close, group_id

   group_id = H5G_OPEN(file_id, 'dimensions')
      att_id=H5A_OPEN_NAME(group_id, 'ndxmax')
         ndxmax=H5A_READ(att_id)
         ndxmax=ndxmax(0)
      h5a_close, att_id
      att_id=H5A_OPEN_NAME(group_id, 'ndymax')
         ndymax=H5A_READ(att_id)
         ndymax=ndymax(0)
      h5a_close, att_id
      att_id=H5A_OPEN_NAME(group_id, 'ndzmax')
         ndzmax=H5A_READ(att_id)
         ndzmax=ndzmax(0)
      h5a_close, att_id
   h5g_close, group_id
;
   group_id = h5g_open(file_id, 'coordinates')
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
   group_id = h5g_open(file_id, 'searchlight3d')
      dset_id=h5d_open(group_id, 'int3d_sc')
         int3d_sc=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'int3d_fvm')
         int3d_fvm=h5d_read(dset_id)
      h5d_close, dset_id
   h5g_close, group_id
;
h5f_close, file_id
;
;-----------------------------------------------------------------------
;
xlim=[-2.5d0,2.5d0]
ylim=[-2.5d0,2.5d0]
;
dist=5.d0;6.165d0
;
;------------CALCULATE 2D PROFILE (INTENSITY THROUGH A SURFACE)---------
;
;calculate coordinates of a plane perpendicular to ray direction
;   in 3d coordinate system
iperp2d, xlim=xlim, ylim=ylim, dist=dist, int3d_sc, int3d_fvm, x, y, z, n_x, n_y, n_z, $
             int2d_sc, int2d_fvm, xt, yt
;
;calculate the integral over the area
sumx_sc=fltarr(n_elements(yt))*0.d0
sumx_fvm=fltarr(n_elements(yt))*0.d0
for j=0, n_elements(yt)-1 do begin
   for i=1, n_elements(xt)-1 do begin
      sumx_sc(j)=sumx_sc(j)+(int2d_sc(i,j)+int2d_sc(i-1,j))*0.5d0*(xt(i)-xt(i-1))
      sumx_fvm(j)=sumx_fvm(j)+(int2d_fvm(i,j)+int2d_fvm(i-1,j))*0.5d0*(xt(i)-xt(i-1))
   endfor
endfor
sumy_sc=0.d0
sumy_fvm=0.d0
for j=1, n_elements(yt)-1 do begin
   sumy_sc=sumy_sc+(sumx_sc(j)+sumx_sc(j-1))*0.5d0*(yt(j)-yt(j-1))
   sumy_fvm=sumy_fvm+(sumx_fvm(j)+sumx_fvm(j-1))*0.5d0*(yt(j)-yt(j-1))
endfor
print, 'integrals'
print, sumy_sc, sumy_fvm, !pi, sumy_sc/!pi, sumy_fvm/!pi
;
;----------------CALCULATE 1D PROFILE ALONG XT AND YT-------------------
;
iperp1d, int2d_sc, int2d_fvm, xt, yt, int1dx_sc, int1dx_fvm, int1dx_theo, $
         int1dy_sc, int1dy_fvm, int1dy_theo
;
;------------CALCULATE INTENSITY ALONG DIRECTION (N_X,N_Y,N_Z)----------
;
itang1d, int3d_sc, int3d_fvm, x, y, z, n_x, n_y, n_z, coord_tang, int1dsc_tang, int1dfvm_tang, int1dtheo_tang
;
;************************2D CONTOUR PLOTS*******************************
;
;plot intensity through a surface perpendicular to ray direction
;
;-----------------------------------------------------------------------
;
;title strings
titlestr=textoidl('(n_x,n_y,n_z)=(') + string(n_x, format='(f4.2)') + ', ' + $
                                       string(n_y, format='(f4.2)') + ', ' + $
                                       string(n_z, format='(f4.2)') + '), ' + $
         textoidl('at r=') + string(dist, format='(f5.2)')
titlestr1=titlestr + ',   (SC)'
titlestr2=titlestr + ',   (FVM)'
;
ctitleStr=textoidl('I/I_{c}')
xtitleStr=textoidl('x_t')
ytitleStr=textoidl('y_t')
;
;-----------------------------------------------------------------------
;
;color range
min_color = 0.d0
max_color = 1.001d0   ;because otherwise, black line in contour when i>1 (for quadratic interpolations)
;
nx_colors=201
ny_colors=201
;
colors=fltarr(nx_colors, ny_colors)
for i=0, nx_colors-1 do begin
   for j=0, ny_colors-1 do begin
      ran= min_color + RANDOMU(seed, /UNIFORM)*(max_color-min_color)
      colors(i,j) = ran
   endfor
endfor
;
loadct, 0
get_contour_colors, colors, x, y, ncolors, nlevels_iso, bottom, $
                    levels_final, levels_iso, c_colors_final,  $
                    cb_indx, cb_tickmark_name
;
;-----------------------------------------------------------------------
;
;plot
;
if(not keyword_set(windx)) then windx=0
;
if(keyword_set(oname)) then begin
   oname='ps_files/benchmark10_a.ps'
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
loadct, 0
;
contourplots_double, int2d_sc, int2d_fvm, xt, xt, yt, yt, $
                     ncolors, nlevels_iso, bottom, $
                     levels_final, levels_iso, c_colors_final, $
                     cb_indx, cb_tickmark_name, $
                     xlim=xlim, $
                     ylim=ylim, $
                     titlestr1=titleStr1, $
                     titlestr2=titleStr2, $
                     xtitlestr=xtitleStr, $
                     ytitlestr=ytitlestr, $
                     ctitlestr=ctitlestr, $
                     /isoc, $
                     ctable=13, $
                     isotropic=isotropic, /cb_top, $
                     oplot_circle=1.d0
;
if keyword_set(oname) then begin
   device, /close
   set_plot,'x'
endif
;
;************************2D SURFACE PLOTS*******************************
;
;plot intensity through a surface perpendicular to ray direction
;
;-----------------------------------------------------------------------
;
;title strings
titlestr=textoidl('(n_x,n_y,n_z)=(') + string(n_x, format='(f4.2)') + ', ' + $
                                       string(n_y, format='(f4.2)') + ', ' + $
                                       string(n_z, format='(f4.2)') + '), ' + $
         textoidl('at r=') + string(dist, format='(f5.2)')
titlestr1=titlestr + ',   (SC)'
titlestr2=titlestr + ',   (FVM)'
;
ztitleStr=textoidl('I/I_{c}')
xtitleStr=textoidl('x_t')
ytitleStr=textoidl('y_t')
;
;-----------------------------------------------------------------------
;
;plot
;
if(keyword_set(oname)) then begin
   oname='ps_files/benchmark10_b.ps'
   print, "writing output to: ", oname
   set_plot,'ps'
   device,file=oname, $;, xsize=19., ysize=26.7, xoffset=1., yoffset=1. , $
    decomposed=0, color=1, bits_per_pixel=8
endif else begin
   window, windx, xsize=1400, ysize=600
   device, decomposed=0
   windx=windx+1
endelse
;
loadct, 0
;
!p.multi=[2,2,1]
surface, int2d_sc, xt, yt, skirt=0.d0, charsize=3., zrange=[0.d0,1.d0], xrange=xlim, yrange=ylim, $
         xtitle=xtitleStr, ytitle=ytitleStr, ztitle=ztitleStr, title='SC'

!p.multi=[1,2,1]
surface, int2d_fvm, xt, yt, skirt=0.d0, charsize=3., zrange=[0.d0,1.d0], xrange=xlim, yrange=ylim, $
         xtitle=xtitleStr, ytitle=ytitleStr, ztitle=ztitleStr, title='FVM'
;
!p.multi=0
;
if keyword_set(oname) then begin
   device, /close
   set_plot,'x'
endif
;
;******************1D PERPENDICULAR AND TANGENTIAL PLOTS****************
;
;intensity profiles along xt, yt
;
titlestr=textoidl('(n_x,n_y,n_z)=(') + string(n_x, format='(f4.2)') + ', ' + $
                                       string(n_y, format='(f4.2)') + ', ' + $
                                       string(n_z, format='(f4.2)') + '), ' + $
         textoidl('at r=') + string(dist, format='(f5.2)')
;
titlestr2=textoidl('(n_x,n_y,n_z)=(') + string(n_x, format='(f4.2)') + ', ' + $
                                       string(n_y, format='(f4.2)') + ', ' + $
                                       string(n_z, format='(f4.2)') + ')'

xtitleStr=textoidl('p')
ytitleStr=textoidl('I/I_{c}')
;
;-----------------------------------------------------------------------
;
if(keyword_set(oname)) then begin
   oname='ps_files/benchmark10_c.ps'
   print, "writing output to: ", oname
   set_plot,'ps'
   device,file=oname, $ ;xsize=19., ysize=26.7, xoffset=1., yoffset=1. , $
    decomposed=0, color=1, bits_per_pixel=8
endif else begin
   window, windx, ysize=800, xsize=800
   device, decomposed=0
   windx=windx+1
endelse
;
loadct_lev, ci_black, ci_blue, ci_red, ci_cyan, ci_green, ci_magenta, ci_yellow, ci_white
;
!p.multi=[2,1,2]
plot, xt, int1dx_theo, $
      title=titlestr, $
      xtitle=xtitlestr, $
      ytitle=ytitlestr, $
      xrange=xlim, $
      yrange=[0.,1.05], /ys, $
      charsize=1.2
;along xt direction
oplot, xt, int1dx_sc, color=ci_blue
oplot, xt, int1dx_fvm, color=ci_red
;along yt direction
oplot, yt, int1dy_sc, color=ci_blue, line=2
oplot, yt, int1dy_fvm, color=ci_red, line=2
;
legstr0='theoretical solution'
legstr1=textoidl('SC  along x_{t} at y_{t}=0')
legstr2=textoidl('FVM along x_{t} at y_{t}=0')
legstr3=textoidl('SC  along y_{t} at x_{t}=0')
legstr4=textoidl('FVM along y_{t} at x_{t}=0')
;
if(keyword_set(oname)) then begin
   legend, [legstr0, legstr1, legstr2, legstr3, legstr4], $
           line=[0,0,0,2,2], $
           color=[ci_black, ci_blue, ci_red, ci_blue, ci_red]
endif else begin
   legend, [legstr0, legstr1, legstr2, legstr3, legstr4], $
           line=[0,0,0,2,2], $
           color=[ci_white, ci_blue, ci_red, ci_blue, ci_red]
endelse
;
;-----------------------------------------------------------------------
;
!p.multi=[1,1,2]
;
loadct_lev, ci_black, ci_blue, ci_red, ci_cyan, ci_green, ci_magenta, ci_yellow, ci_white
;
plot, coord_tang, int1dtheo_tang, $
      title=titlestr2, $
      xtitle=textoidl('r'), $
      ytitle=ytitlestr, $
      xrange=[1.d0,10.d0], $
      yrange=[0.d0,1.05d0], /ys, $
      charsize=1.2
oplot, coord_tang, int1dsc_tang, color=ci_blue
oplot, coord_tang, int1dfvm_tang, color=ci_red
;
legstr0='theoretical solution'
legstr1=textoidl('short characteristics')
legstr2=textoidl('finite volume method')
;
if(keyword_set(oname)) then begin
   legend, [legstr0, legstr1, legstr2], $
           line=[0,0,0], $
           color=[ci_black, ci_blue, ci_red], $
           /bottom
endif else begin
   legend, [legstr0, legstr1, legstr2], $
           line=[0,0,0], $
           color=[ci_white, ci_blue, ci_red], $
           /bottom
endelse
;
!p.multi=0
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
pro calc_plane, dist=dist, xlim=xlim, ylim=ylim, n_dir, dim_xt, dim_yt, pos_3d, xt, yt
;
;input:   n_dir: normal vector of plane, for which physical value shall be
;                plotted (e.g. direction vector of intensity)
;         x,y,z: x, y, z-axes of original coordinate system
;         dim_xt, dim_yt: dimensions of new coordinate system (plane)
;
;keywords: dist:  (central) distance from origin, at which profile shall be
;                 calculated
;          xlim, ylim: range of 2d surface
;
;output:  pos_3d: positions of grid-points from new coordinate system (plane)
;                 in the old coordinate system
;         xt, yt: new coordinate system (plane)
;
;-------------------check if n_dir has correct dimension----------------
;
if(n_elements(n_dir) ne 3) then begin
   print, 'n_dir has wrong dimension'
   stop
endif
;
;--------------define distance where profile shall be taken-------------
;
if(not keyword_set(dist)) then dist=10.d0
;
;new origin:
origin=dist*n_dir
;
;---------calculate new unit_vectors of new coordinate system-----------
;
;unit vector in z-direction
e_z=n_dir
;
;unit vector in x-direction
e_x=unit_vec(e_z)
;
;unit vector in y-direction
e_y= cross_product(e_z, e_x)
;
;---------------calculate grid in new coordinate system-----------------
;
if(not keyword_set(xlim)) then xlim=[-3.d0,3.d0]
if(not keyword_set(ylim)) then ylim=[-3.d0,3.d0]
;
xmin=xlim(0)
xmax=xlim(1)
ymin=ylim(0)
ymax=ylim(1)
;
xt=xmin+findgen(dim_xt)*(xmax-xmin)/(dim_xt-1)
yt=ymin+findgen(dim_yt)*(ymax-ymin)/(dim_yt-1)
;
;---------calculate position vectors in old coordinate system-----------
;
pos_3d=fltarr(dim_xt,dim_yt,3)*0.d0
;
for i=0, dim_xt-1 do begin
   for j=0, dim_yt-1 do begin
      pos_3d(i,j,0)=origin(0)+xt(i)*e_x(0) + yt(j)*e_y(0)
      pos_3d(i,j,1)=origin(1)+xt(i)*e_x(1) + yt(j)*e_y(1)
      pos_3d(i,j,2)=origin(2)+xt(i)*e_x(2) + yt(j)*e_y(2)
   endfor
endfor
;
end
;
;-----------------------------------------------------------------------
;-----------------------------------------------------------------------
;-----------------------------------------------------------------------
;
function unit_vec, e_3
;
if(n_elements(e_3) ne 3) then begin
   print, 'wrong dimensions in function unit_vec'
   stop
endif
;   
xt1=1.d0
yt1=0.d0
zt1=-e_3(0)/e_3(2)
;
e_1=fltarr(3)*0.d0
e_1(0)=xt1
e_1(1)=yt1
e_1(2)=zt1
;
length=sqrt(xt1^2 + yt1^2 + zt1^2)
;
e_1=e_1/length
;
return, e_1
;
;
end
;
;-----------------------------------------------------------------------
;-----------------------------------------------------------------------
;-----------------------------------------------------------------------
;
function cross_product, vec1, vec2
;
if(n_elements(vec1) ne 3 or n_elements(vec2) ne 3) then begin
   print, 'wrong dimensions in function cross_product'
   stop
endif
;
xcoord=vec1(1)*vec2(2)-vec1(2)*vec2(1)
ycoord=vec1(2)*vec2(0)-vec1(0)*vec2(2)
zcoord=vec1(0)*vec2(1)-vec1(1)*vec2(0)
;
vec3=fltarr(3)*0.d0
vec3(0)=xcoord
vec3(1)=ycoord
vec3(2)=zcoord
;
return, vec3
;
;
end
;
;-----------------------------------------------------------------------
;-----------------------------------------------------------------------
;-----------------------------------------------------------------------
;
pro iperp2d, xlim=xlim, ylim=ylim, dist=dist, int3d_sc, int3d_fvm, x, y, z, n_x, n_y, n_z, $
             int2d_sc, int2d_fvm, xt, yt
;
;calculates intensity profile (through a surface perpendicular to direction
;   n=(n_x,n_y,n_z))
;
;surface limits
if(not keyword_set(xlim)) then xlim=[-2.d0,2.d0]
if(not keyword_set(ylim)) then ylim=[-2.d0,2.d0]
if(not keyword_set(dist)) then dist=2.5d0
;
nd_xt=151
nd_yt=151
;
;calculate perpendicular plane with corresponding 3d coordinates
calc_plane, dist=dist, xlim=xlim, ylim=ylim, [n_x,n_y,n_z], nd_xt, nd_yt, pos_3d, xt, yt
;
int2d_sc=fltarr(nd_xt,nd_yt)*0.d0
int2d_fvm=fltarr(nd_xt,nd_yt)*0.d0
;
;interpolate the intensity
pp=fltarr(3)*0.d0
for i=0, nd_xt-1 do begin
   for j=0, nd_yt-1 do begin
      pp(0)=pos_3d(i,j,0)
      pp(1)=pos_3d(i,j,1)
      pp(2)=pos_3d(i,j,2)
;
      interpol_trilinear, int3d_sc, x, y, z, pp, value_pp
;      interpol_triquad, int3d_sc, x, y, z, pp, value_pp
      int2d_sc(i,j)=value_pp
;
      interpol_trilinear, int3d_fvm, x, y, z, pp, value_pp
;      interpol_triquad, int3d_fvm, x, y, z, pp, value_pp
      int2d_fvm(i,j)=value_pp
   endfor
endfor
;
end
;
;-----------------------------------------------------------------------
;-----------------------------------------------------------------------
;-----------------------------------------------------------------------
;
pro iperp1d, int2d_sc, int2d_fvm, xt, yt, int1dx_sc, int1dx_fvm, int1dx_theo, $
             int1dy_sc, int1dy_fvm, int1dy_theo

nx=n_elements(xt)
ny=n_elements(yt)
;
int1dx_sc=fltarr(nx)*0.d0
int1dx_fvm=fltarr(nx)*0.d0
int1dx_theo=fltarr(nx)*0.d0
;
int1dy_sc=fltarr(ny)*0.d0
int1dy_fvm=fltarr(ny)*0.d0
int1dy_theo=fltarr(ny)*0.d0
;
for i=0, nx-1 do begin
   if(abs(xt(i)) le 1.d0) then int1dx_theo(i)=1.d0
   int1dx_sc(i) = int2d_sc(i,ny/2)
   int1dx_fvm(i) = int2d_fvm(i,ny/2)
endfor
;
for i=0, ny-1 do begin
   if(abs(yt(i)) le 1.d0) then int1dy_theo(i)=1.d0
   int1dy_sc(i) = int2d_sc(nx/2,i)
   int1dy_fvm(i) = int2d_fvm(nx/2,i)
endfor
;
end
;
;-----------------------------------------------------------------------
;-----------------------------------------------------------------------
;-----------------------------------------------------------------------
;
pro itang1d, int3d_sc, int3d_fvm, x, y, z, n_x, n_y, n_z, rcoord, int1d_sc, int1d_fvm, int1d_theo
;
;calculates intensity along a line parallel to ray direction
r1=1.d0
r2=10.d0
;
;create radial coordinate
nd=1001
del=alog10(r2/r1)/float(nd-1)
rcoord=fltarr(nd)
rcoord(0)=r1
for i=1, nd-1 do begin
   rcoord(i)=rcoord(i-1)*10.d0^del
endfor
;
;create coordinate along direction
coord_x=fltarr(nd)
coord_y=fltarr(nd)
coord_z=fltarr(nd)
int1d_sc=fltarr(nd)
int1d_fvm=fltarr(nd)
int1d_theo=fltarr(nd)
;
;calculate corresponding x and z coordinates
for i=0, nd-1 do begin
   coord_x(i)=rcoord(i)*n_x
   coord_y(i)=rcoord(i)*n_y
   coord_z(i)=rcoord(i)*n_z
endfor
;
;---------interpolate onto line and calculate theoretic solution--------
;
for i=0, nd-1 do begin
;   interpol_trilinear, int3d_sc, x, y, z, [coord_x(i),coord_y(i),coord_z(i)], valp
   interpol_triquad, int3d_sc, x, y, z, [coord_x(i),coord_y(i),coord_z(i)], valp

   int1d_sc(i)=valp
;   interpol_trilinear, int3d_fvm, x, y, z, [coord_x(i),coord_y(i),coord_z(i)], valp
   interpol_triquad, int3d_fvm, x, y, z, [coord_x(i),coord_y(i),coord_z(i)], valp
   int1d_fvm(i)=valp
   if(abs(rcoord(i)) ge 1.) then begin
      int1d_theo(i)=1.d0
   endif else begin
      int1d_theo(i)=0.d0
   endelse
endfor
;
;maniuplate int1d at the photosphere, because
;   interpolation scheme would interpolate to zero (value inside the core)
for i=0, nd-1 do begin
   rad=rcoord(i)
   if(rad lt 1.05d0) then begin
      int1d_sc(i)=1.d0
      int1d_fvm(i)=1.d0
   endif
endfor
;
end
