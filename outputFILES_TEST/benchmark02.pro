pro benchmark02, dir=dir, windx=windx, oname=oname, xlim=xlim, ylim=ylim, clim=clim, logscale=logscale
;
;plots contours of searchlight along a direction
;
;------------------READ ALL INFORMATION FROM HDF5-FILE------------------
;
if(not keyword_set(dir)) then dir='.'
fname=dir+'/benchmark02.h5'
;
file_id = h5f_open(fname)
;
   group_id = H5G_OPEN(file_id, 'input_parameters')
      att_id=H5A_OPEN_NAME(group_id, 'xic1')
         xic1=H5A_READ(att_id)
         xic1=xic1(0)
      h5a_close, att_id
   h5g_close, group_id

   group_id = H5G_OPEN(file_id, 'direction')
      att_id=H5A_OPEN_NAME(group_id, 'mu')
         mu=H5A_READ(att_id)
         mu=mu(0)
      h5a_close, att_id
   h5g_close, group_id

   group_id = H5G_OPEN(file_id, 'dimensions')
      att_id=H5A_OPEN_NAME(group_id, 'ndxmax')
         ndxmax=H5A_READ(att_id)
         ndxmax=ndxmax(0)
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
      dset_id=h5d_open(group_id, 'z')
         z=h5d_read(dset_id)
      h5d_close, dset_id
   h5g_close, group_id
;
   group_id = h5g_open(file_id, 'searchlight2d')
      dset_id=h5d_open(group_id, 'int2d_sc')
         int2d_sc=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'int2d_fvm')
         int2d_fvm=h5d_read(dset_id)
      h5d_close, dset_id
   h5g_close, group_id
;
h5f_close, file_id
;
;normalize intensities
int2d_sc=int2d_sc/xic1
int2d_fvm=int2d_fvm/xic1
;
;----------------------CALCULATE PROFILE ALONG A LINE-------------------
;
radius=3.d0
iperp, int2d_sc, x, z, mu, coord_line, xcoord_line, zcoord_line, nd, int1dsc_perp, int1d_theo, radius=radius
iperp, int2d_fvm, x, z, mu, coord_line2, xcoord_line2, zcoord_line2, nd2, int1dfvm_perp, int1d_theo2, radius=radius
;
;integrate along profile
sum_profilesc = integ_trapez(coord_line, int1dsc_perp)
sum_profilefvm = integ_trapez(coord_line2, int1dfvm_perp)
sum_profiletheo = integ_trapez(coord_line, int1d_theo)
;
print, 'integral over profile (sc)  ', sum_profilesc, int_tabulated(coord_line, int1dsc_perp)
print, 'integral over profile (fvm) ', sum_profilefvm, int_tabulated(coord_line2, int1dfvm_perp)
print, 'integral over profile (theo)', sum_profiletheo, int_tabulated(coord_line, int1d_theo)
;
;----------------------CALCULATE PROFILE ALONG DIRECTION----------------
;
itang, int2d_sc, x, z, mu, coord_tang, xcoord_tang, zcoord_tang, nd_tang, int1dsc_tang, int1d_tang_theo
itang, int2d_fvm, x, z, mu, coord_tang2, xcoord_tang2, zcoord_tang2, nd_tang2, int1dfvm_tang, int1d_tang_theo2
;
;-----------------------DEFINE RANGE------------------------------------
;
if(keyword_set(xlim)) then begin
   xmin=min(xlim)
   xmax=max(xlim)
endif else begin
   xmin=min([min(x), min(z)])
   xmax=max([max(x), max(z)])
   xlim=[xmin, xmax]
endelse
;
if(keyword_set(ylim)) then begin
   ymin=min(ylim)
   ymax=max(ylim)
endif else begin
   ymin=min([min(x), min(z)])
   ymax=max([max(x), max(z)])
   ylim=[ymin, ymax]
endelse
;
;----------------------------TITLE-STRINGS------------------------------
;
;titlestr=textoidl('i (\theta=') + string(acos(mu), format='(f9.5)') + $
;         textoidl(', \phi=') + string(phi, format='(f9.5)') + ')'
titleStr1='short characteristics'
titleStr2='finite volume method'
ctitleStr=textoidl('I/I_{c}')
xtitleStr=textoidl('x')
ytitleStr=textoidl('z')
;
;-----------------------------COLOR RANGE-------------------------------
;
if(keyword_set(logscale)) then begin
   min_color = min(int2d_sc(where(int2d_sc gt 0.)))
   max_color = max(int2d_sc)
   
   min_color = alog10(min_color)
   max_color = alog10(max_color)
   int2d_sc = alog10(int2d_sc)
endif else begin
   min_color = min(int2d_sc)
   max_color = max(int2d_sc)
endelse
;
;overwrite colors, if keyword clim is supplied
if(keyword_set(clim)) then begin
   min_color = clim(0)
   max_color = clim(1)
endif
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
get_contour_colors, colors, x, z, ncolors, nlevels_iso, bottom, $
                    levels_final, levels_iso, c_colors_final,  $
                    cb_indx, cb_tickmark_name
;
;-----------------------------------------------------------------------
;
if(not keyword_set(windx)) then windx=0
;
if(keyword_set(oname)) then begin
   oname='ps_files/benchmark02_a.ps'
   print, "writing output to: ", oname
   set_plot,'ps'
   device,file=oname, $ ;xsize=19., ysize=26.7, xoffset=1., yoffset=1. , $
    decomposed=0, color=1, bits_per_pixel=8
endif else begin
   window, windx
   device, decomposed=0
   windx=windx+1
endelse
;
loadct, 0
;
;contourplots, int2d_sc, x, z, xlim, ylim, $
;              titlestr, xtitlestr, ytitlestr, $
;              ncolors, nlevels_iso, bottom, $
;              levels_final, levels_iso, c_colors_final, $
;              cb_indx, cb_tickmark_name, ctable=13, /isoc, /isotropic, ctitlestr=ctitlestr

contourplots_double, int2d_sc, int2d_fvm, x, x, z, z, $
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
;                     /isoc, $
                     ctable=13, $
                     isotropic=isotropic, /cb_top, $
                     oplot_circle=1.d0, $
                     oplot1_x1=xcoord_line, oplot1_y1=zcoord_line, $
;                     oplot1_x2=xcoord_tang, oplot1_y2=zcoord_tang, $
                     oplot2_x1=xcoord_line2, oplot2_y1=zcoord_line2, $
                     /oplot_grid, $
                     oplot1_style=0, oplot2_style=0



;
if keyword_set(oname) then begin
   device, /close
   set_plot,'x'
endif
;
;-----------------------------------------------------------------------
;
if(keyword_set(oname)) then begin
   oname='ps_files/benchmark02_b.ps'
   print, "writing output to: ", oname
   set_plot,'ps'
   device,file=oname, $ ;xsize=19., ysize=26.7, xoffset=1., yoffset=1. , $
    decomposed=0, color=1, bits_per_pixel=8
endif else begin
   window, windx
   device, decomposed=0
   windx=windx+1
endelse
;
loadct, 0
loadct_lev, ci_black, ci_blue, ci_red, ci_cyan, ci_green, ci_magenta, ci_yellow, ci_white
!p.multi=[2,1,2]
plot, coord_line, int1d_theo, $
      title=textoidl('at radius') + string(radius, format='(f6.2)'), $
      xtitle=textoidl('q'), $
      ytitle=textoidl('I/I_{c}'), $
      yrange=[0.,1.05], /ys, $
      xrange=[-2.,2.]
oplot, coord_line, int1dsc_perp, line=0, color=ci_blue
oplot, coord_line2, int1dfvm_perp, line=0, color=ci_red
;
legstr1='theoretical solution'
legstr2='short characteristics'
legstr3='finite volume method'
if(keyword_set(oname)) then begin
   legend, [legstr1, legstr2, legstr3], $
           psym=[0,0,0], $
           linestyle=[0,0,0], $
           color=[ci_black,ci_blue,ci_red], $
           textcolor=ci_black
endif else begin
   legend, [legstr1, legstr2, legstr3], $
           psym=[0,0,0], $
           linestyle=[0,0,0], $
           color=[ci_white,ci_blue,ci_red], $
           textcolor=ci_white
endelse

!p.multi=[1,1,2]
plot, coord_tang, int1d_tang_theo, $
      xtitle=textoidl('r'), $
      ytitle=textoidl('I/I_{c}'), $
      yrange=[0.,1.05], /ys
oplot, coord_tang, int1dsc_tang, line=0, color=ci_blue
oplot, coord_tang2, int1dfvm_tang, line=0, color=ci_red
if(keyword_set(oname)) then begin
   legend, [legstr1, legstr2, legstr3], $
           psym=[0,0,0], $
           linestyle=[0,0,0], $
           color=[ci_black,ci_blue,ci_red], $
           textcolor=ci_black
endif else begin
   legend, [legstr1, legstr2, legstr3], $
           psym=[0,0,0], $
           linestyle=[0,0,0], $
           color=[ci_white,ci_blue,ci_red], $
           textcolor=ci_white
endelse

if keyword_set(oname) then begin
   device, /close
   set_plot,'x'
endif

!p.multi=0
;
end
;
;-----------------------------------------------------------------------
;-----------------------------------------------------------------------
;-----------------------------------------------------------------------
;
pro iperp, int2d, x, z, mu, q, xq, zq, nd, int1d, int1d_theo, radius=radius
;
;calculates intensity along a line perpendicular to ray direction
if(n_elements(radius) eq 1) then begin
   r1=radius
endif else begin
   r1=3.d0
endelse
x0=r1*sqrt(1.d0-mu^2)
z0=r1*mu
x1=0.d0
z1=r1/mu
;
;create coordinate along line
coord_min=-3.d0
coord_max=3.d0
nd=1001
coord_arr=coord_min+findgen(nd)*(coord_max-coord_min)/float(nd-1)
coord_x=fltarr(nd)
coord_z=fltarr(nd)
int1d=fltarr(nd)
int1d_theo=fltarr(nd)
;
;calculate corresponding x and z coordinates
for i=0, nd-1 do begin
   rlim=sqrt(r1^2 + coord_arr(i)^2)
   alpha=atan(coord_arr(i)/r1)
   alpha0=acos(mu)+alpha
   coord_x(i)=rlim*sin(alpha0)
   coord_z(i)=rlim*cos(alpha0)
endfor
;
;store in global array
q=coord_arr
xq=coord_x
zq=coord_z
;
;---------interpolate onto line and calculate theoretic solution--------
;
for i=0, nd-1 do begin
   interpol_bilinear, int2d, x, z, [xq(i),zq(i)], valp
   int1d(i)=valp
   if(abs(q(i)) le 1.) then begin
      int1d_theo(i)=1.d0
   endif else begin
      int1d_theo(i)=0.d0
   endelse
;   print, i, q(i), int1d(i), int1d_theo(i)
endfor
;
end
;
;-----------------------------------------------------------------------
;-----------------------------------------------------------------------
;-----------------------------------------------------------------------
;
pro itang, int2d, x, z, mu, q, xq, zq, nd, int1d, int1d_theo
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
coord_z=fltarr(nd)
int1d=fltarr(nd)
int1d_theo=fltarr(nd)
;
;calculate corresponding x and z coordinates
for i=0, nd-1 do begin
   coord_x(i)=rcoord(i)*sqrt(1.d0-mu^2)
   coord_z(i)=rcoord(i)*mu
endfor
;
;store in global array
q=rcoord
xq=coord_x
zq=coord_z
;
;---------interpolate onto line and calculate theoretic solution--------
;
for i=0, nd-1 do begin
   interpol_bilinear, int2d, x, z, [xq(i),zq(i)], valp
   int1d(i)=valp
   if(abs(q(i)) ge 1.) then begin
      int1d_theo(i)=1.d0
   endif else begin
      int1d_theo(i)=0.d0
   endelse
endfor
;
;maniuplate int1d at the photosphere, because
;   interpolation scheme would interpolate to zero (value inside the core)
for i=0, nd-1 do begin
   rad=q(i)
   if(rad lt 1.05d0) then begin
      int1d(i)=1.d0
   endif
endfor
;
;-----------------------------------------------------------------------
;even better: take grid points directly (only possible for direction 45
;degree)
;-----------------------------------------------------------------------
mu_crit=cos(45.d0*!pi/180.d0)
if(mu ne mu_crit) then return
;
for i=0, n_elements(x)-1 do begin
   if(x(i) ne z(i)) then begin
      print, 'x unequal z'
      stop
   endif
endfor
;
ndxmax=n_elements(x)
;
nx_origin=(n_elements(x)-1)/2
;
nr=0
for i=nx_origin, ndxmax-1 do begin
   rad=sqrt(x(i)^2+z(i)^2)
;   print, rad
   if(rad gt 1.d0) then begin
      nr=nr+1
   endif
endfor
;
int1d=fltarr(nr+1)
int1d_theo=fltarr(nr+1)
q=fltarr(nr+1)
xq=fltarr(nr+1)
zq=fltarr(nr+1)
nd=nr
;
for i=0, nr-1 do begin
   indx=ndxmax-nr+i
   q(i+1)=sqrt(x(indx)^2+z(indx)^2)
   xq(i+1)=x(indx)
   zq(i+1)=z(indx)
   int1d(i+1)=int2d(indx,indx)
   int1d_theo(i+1)=1.d0
endfor
;
q(0)=1.d0
xq(0)=sqrt(1.d0-mu^2)
zq(0)=mu
int1d(0)=1.d0
int1d_theo(0)=1.d0
;
;
end
