pro searchlight_test01, dir, muindx=muindx, phiindx=phiindx, windx=WINDX, oname=ONAME, xlim=XLIM, ylim=YLIM, clim=CLIM, logscale=LOGSCALE
;
;plots contours of searchlight along a direction
;
IF (N_ELEMENTS(dir) EQ 0) THEN BEGIN
   PRINT, 'SYNTAX: '
   PRINT, 'searchlight_test01, dir'
   PRINT, 'with dir=string, directory of model, e.g. "M1" '
   RETURN
ENDIF
;
;------------------READ ALL INFORMATION FROM HDF5-FILE------------------
;
fname=dir+'/searchlight_test01.h5'
;
file_id = h5f_open(fname)
;
   group_id = H5G_OPEN(file_id, 'dimensions')
      att_id=H5A_OPEN_NAME(group_id, 'dim_mu')
         dim_mu=H5A_READ(att_id)
         dim_mu=dim_mu(0)
      h5a_close, att_id
      att_id=H5A_OPEN_NAME(group_id, 'dim_phi')
         dim_phi=H5A_READ(att_id)
         dim_phi=dim_phi(0)
      h5a_close, att_id
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
   group_id = H5G_OPEN(file_id, 'options')
      att_id=H5A_OPEN_NAME(group_id, 'indx_mu')
         indx_mu=H5A_READ(att_id)
         indx_mu=indx_mu(0)
      h5a_close, att_id
      att_id=H5A_OPEN_NAME(group_id, 'indx_phi')
         indx_phi=H5A_READ(att_id)
         indx_phi=indx_phi(0)
      h5a_close, att_id
      att_id=H5A_OPEN_NAME(group_id, 'logical_all_angles')
         lang=H5A_READ(att_id)
         lang=lang(0)
      h5a_close, att_id
   h5g_close, group_id
;
   group_id = h5g_open(file_id, 'coordinates')
      dset_id=h5d_open(group_id, 'nodes_mu')
         nodes_mu=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'nodes_phi')
         nodes_phi=h5d_read(dset_id)
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
   group_id = h5g_open(file_id, 'searchlight')
      dset_id=h5d_open(group_id, 'int3d_angdep')
         int3d_angdep=h5d_read(dset_id)
      h5d_close, dset_id
   h5g_close, group_id
;
h5f_close, file_id
;
;-----------------------------------------------------------------------
;
;from fortran to idl indices
indx_mu=indx_mu-1
indx_phi=indx_phi-1
;
if(lang eq 1) then begin
   if(keyword_set(muindx)) then begin
      if(muindx gt dim_mu-1) then begin
         print, 'maximum allowed muindx:', dim_mu-1
         return
      endif
      indx_mu=muindx
   endif else begin
      indx_mu=0   ;default
   endelse
   if(keyword_set(phiindx)) then begin
      if(phiindx gt dim_phi-1) then begin
         print, 'maximum allowed phiindx:', dim_phi-1
         return
      endif
      indx_phi=phiindx
   endif else begin
      indx_phi=0   ;default
   endelse
endif
;
indx_ycoord=ndymax/2
mu=nodes_mu(indx_mu)
phi=nodes_phi(indx_phi)
;
if(lang eq 1) then begin
   int2d=int3d_angdep(*,indx_ycoord,*,indx_mu,indx_phi)
endif else begin
   int2d=fltarr(ndxmax,ndzmax)
   for i=0, ndxmax-1 do begin
      for j=0, ndzmax-1 do begin
         int2d(i,j)=int3d_angdep(i,indx_ycoord,j,0,0)
      endfor
   endfor
endelse
;
;----------------------CALCULATE PROFILE ALONG A LINE-------------------
;
iperp, int2d, x, z, mu, coord_line, xcoord_line, zcoord_line, nd, int1d, int1d_theo
;
;-----------------------DEFINE RANGE------------------------------------
;
IF(KEYWORD_SET(xlim)) THEN BEGIN
   xmin=MIN(xlim)
   xmax=MAX(xlim)
ENDIF ELSE BEGIN
   xmin=MIN([MIN(X), MIN(Z)])
   xmax=MAX([MAX(X), MAX(Z)])
   xlim=[xmin, xmax]
ENDELSE
;
IF(KEYWORD_SET(ylim)) THEN BEGIN
   ymin=MIN(ylim)
   ymax=MAX(ylim)
ENDIF ELSE BEGIN
   ymin=MIN([MIN(X), MIN(Z)])
   ymax=MAX([MAX(X), MAX(Z)])
   ylim=[ymin, ymax]
ENDELSE
;
;----------------------------TITLE-STRINGS------------------------------
;
;titleStr=textoidl('I (\Theta=') + STRING(acos(mu), FORMAT='(F9.5)') + $
;         textoidl(', \Phi=') + STRING(phi, FORMAT='(F9.5)') + ')'
titleStr=''
ctitleStr=textoidl('I/I_{c}')
xtitleStr=textoidl('x')
ytitleStr=textoidl('z')
;
;-----------------------------COLOR RANGE-------------------------------
;
IF(KEYWORD_SET(LOGSCALE)) THEN BEGIN
   min_color = min(int2d(where(int2d gt 0.)))
   max_color = max(int2d)
   
   min_color = alog10(min_color)
   max_color = alog10(max_color)
   int2d = alog10(int2d)
ENDIF ELSE BEGIN
   min_color = min(int2d)
   max_color = max(int2d)
ENDELSE
;
;OVERWRITE COLORS, IF KEYWORD CLIM IS SUPPLIED
IF(KEYWORD_SET(CLIM)) THEN BEGIN
   min_color = clim(0)
   max_color = clim(1)
ENDIF
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
get_contour_colors, colors, x, z, ncolors, nlevels_iso, bottom, $
                    levels_final, levels_iso, c_colors_final,  $
                    cb_indx, cb_tickmark_name
;
;-----------------------------------------------------------------------
;
IF(KEYWORD_SET(WINDX)) THEN BEGIN
   WINDOWINDX=WINDX
ENDIF ELSE BEGIN
   WINDOWINDX=0
ENDELSE
;
IF(KEYWORD_SET(ONAME)) THEN BEGIN
   PRINT, "WRITING OUTPUT TO: ", ONAME
   set_plot,'ps'
   device,file=ONAME, $ ;XSIZE=19., YSIZE=26.7, XOFFSET=1., YOFFSET=1. , $
    decomposed=0, color=1, BITS_PER_PIXEL=8
ENDIF ELSE BEGIN
   window, windowindx
   device, decomposed=0
   windowindx=windowindx+1
ENDELSE
;
loadct, 0
;
CONTOURPLOTS, INT2D, X, Z, xlim, ylim, $
              titleStr, xtitleStr, ytitleStr, $
              NCOLORS, NLEVELS_ISO, BOTTOM, $
              LEVELS_FINAL, LEVELS_ISO, C_COLORS_FINAL, $
              CB_INDX, CB_TICKMARK_NAME, ctable=13, /ISOC, /ISOTROPIC, ctitleStr=ctitleStr


;overplot grid
loadct, 0
for i=0, ndxmax-1 do begin
;   oplot, [x(i),x(i)], [ylim(0),ylim(1)]
endfor
for i=0, ndzmax-1 do begin
;   oplot, [xlim(0),xlim(1)], [z(i),z(i)]
endfor
;
;overplot photosphere
plotcircle, 1.d0, /oplot
;
;overplot line for which profile is calculated
oplot, [xcoord_line(0),xcoord_line(nd-1)], [zcoord_line(0),zcoord_line(nd-1)]

;
IF KEYWORD_SET(ONAME) THEN BEGIN
   device, /close
   set_plot,'x'
ENDIF
;
window, 1
plot, coord_line, int1d, $
      xtitle=textoidl('q'), $
      ytitle=textoidl('I/I_{c}'), $
      yrange=[0.,1.]
oplot, coord_line, int1d_theo, line=1

;
end

pro iperp, int2d, x, z, mu, q, xq, zq, nd, int1d, int1d_theo
;
;calculates intensity along a line perpendicular to ray direction
r1=3.d0
x0=r1*sqrt(1.d0-mu^2)
z0=r1*mu
x1=0.d0
z1=r1/mu
;
;create coordinate along line
coord_min=-2.d0
coord_max=2.d0
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
   print, i, q(i), int1d(i), int1d_theo(i)
endfor
;
end
