pro plot_sc3d, fname=fname, windx=windx, oname=oname, help=print_help
;
;+
; NAME:
;	plot_sc3d
;
; PURPOSE:
;	This procedure plots:
;          3d model as a function of radius of sc3d.eo output
;          3d solution (source functions) as function of radius of sc3d.eo output
;          convergence behaviour of sc3d.eo       
;
; CALLING SEQUENCE:
;
;	plot_sc3d
;
; INPUTS:
;
; KEYWORD PARAMETERS:
;
;	fname:  Set this keyword to thee file-name of the model
;
;       help:   Set this keyword (flag) to show the documentation of the
;               procedure
;
; EXAMPLE:
;	plot_sc3d, fname='output_model00.h5'
;
;-
;
;-----------------------if help is needed-------------------------------
;
if(keyword_set(print_help)) then begin
   doc_library, 'plot_sc3d'
   return
endif
;
if(n_elements(fname) eq 0) then begin
   print, 'error in plot_sc3d: fname not specified'
   doc_library, 'plot_sc3d'
   stop
endif
;
;------------------read all information from hdf5-file------------------
;
read_sc3d, fname, input_mod_dim=input_mod_dim, spatial_grid1d=spatial_grid1d, spatial_grid3d=spatial_grid3d, $
                  opt_method=opt_method, opt_opal=opt_opal, opt_angint_method=opt_angint_method, opt_sol2d=opt_sol2d, $
                  opt_incl_cont=opt_incl_cont, opt_start_cont=opt_start_cont, opt_ng_cont=opt_ng_cont, opt_ait_cont=opt_ait_cont, $
                  opt_incl_line=opt_incl_line, opt_start_line=opt_start_line, opt_ng_line=opt_ng_line, opt_ait_line=opt_ait_line, $
                  opt_alo_cont=opt_alo_cont, opt_alo_line=opt_alo_line, xic1=xic1, kcont=kcont, kline=kline, eps_cont=eps_cont, $
                  eps_line=eps_line, teff=teff, trad=trad, xnue0=xnue0, rstar=rstar, vth_fiducial=vth_fiducial, vmicro=vmicro, $
                  yhe=yhe, hei=hei, mdot=mdot, na=na, vmin=vmin, vmax=vmax, beta=beta, ndxmax=ndxmax, ndymax=ndymax, ndzmax=ndzmax, $
                  dim_mu=dim_mu, dim_phi=dim_phi, dim_omega=dim_omega, nxobs=nxobs, xarr=x, yarr=y, zarr=z, nodes_mu=nodes_mu, $
                  n_x=n_x, n_y=n_y, n_z=n_z, nodes_xobs=nodes_xobs, itmaxc=itmaxc, itmaxl=itmaxl, devmaxc=devmax, devmaxl=devmaxl, $
                  epsmaxc_arr=epsmaxc_arr, epsmaxl_arr=epsmaxl_arr, scont3d=scont3d, sline3d=sline3d, ssobo3d=ssobo3d, mint3d=mint3d, $
                  mintbar3d=mintbar3d, mask_totreg3d=mask_totreg3d, mask_innreg3d=mask_innreg3d, mask_bpoint3d=mask_bpoint3d, $
                  mask3d=mask3d, t3d=t3d, opac3d=opac3d, opalbar3d=opalbar3d, velx3d=velx3d, vely3d=vely3d, velz3d=velz3d
;
;------------------print out all options and parameters-----------------
;
print, '----------------------------------------'
print, 'input-file   ', fname
print, ''
print, '--------------options-------------------'
print, format='(a20,i20)', 'input_mod_dim', input_mod_dim
print, format='(a20,i20)', 'spatial_grid1d', spatial_grid1d
print, format='(a20,i20)', 'spatial_grid3d', spatial_grid3d
print, format='(a20,i20)', 'opt_method', opt_method
print, format='(a20,i20)', 'opt_opal', opt_opal
print, format='(a20,i20)', 'opt_angint_method', opt_angint_method
print, format='(a20,i20)', 'opt_sol2d', opt_sol2d
print, format='(a20,i20)', 'opt_incl_cont', opt_incl_cont
print, format='(a20,i20)', 'opt_start_cont', opt_start_cont
print, format='(a20,i20)', 'opt_ng_cont', opt_ng_cont
print, format='(a20,i20)', 'opt_ait_cont', opt_ait_cont
print, format='(a20,i20)', 'opt_incl_line', opt_incl_line
print, format='(a20,i20)', 'opt_start_line', opt_start_line
print, format='(a20,i20)', 'opt_ng_line', opt_ng_line
print, format='(a20,i20)', 'opt_ait_line', opt_ait_line
print, format='(a20,i20)', 'opt_alo_cont', opt_alo_cont
print, format='(a20,i20)', 'opt_alo_line', opt_alo_line
print, ''
print, '-------input/stellar parameters---------'
print, format='(a20,e20.8)', 'xic1', xic1
print, format='(a20,e20.8)', 'kcont', kcont
print, format='(a20,e20.8)', 'kline', kline
print, format='(a20,e20.8)', 'eps_cont', eps_cont
print, format='(a20,e20.8)', 'eps_line', eps_line
print, format='(a20,e20.8)', 'teff', teff
print, format='(a20,e20.8)', 'trad', trad
print, format='(a20,e20.8)', 'xnue0', xnue0
print, format='(a20,e20.8)', 'rstar', rstar
print, format='(a20,e20.8)', 'vth_fiducial', vth_fiducial
print, format='(a20,e20.8)', 'vmicro', vmicro
print, format='(a20,e20.8)', 'yhe', yhe
print, format='(a20,e20.8)', 'hei', hei
print, format='(a20,e20.8)', 'mdot', mdot
print, format='(a20,e20.8)', 'na', na
print, format='(a20,e20.8)', 'vmin', vmin
print, format='(a20,e20.8)', 'vmax', vmax
print, format='(a20,e20.8)', 'beta', beta
print, ''
print, '------------dimensions------------------'
print, format='(a20,i20)', 'ndxmax', ndxmax
print, format='(a20,i20)', 'ndymax', ndymax
print, format='(a20,i20)', 'ndzmax', ndzmax
print, format='(a20,i20)', 'dim_omega', dim_omega
print, format='(a20,i20)', 'nxobs', nxobs
print, '----------------------------------------'
;
;------------------prepare 3d radial array and normalization------------
;
;normalize everything to xic1
sline3d=sline3d/xic1
ssobo3d=ssobo3d/xic1
scont3d=scont3d/xic1
;
;calculate dilution factor
sline3d_dilfac=fltarr(ndxmax,ndymax,ndzmax)*0.d0
for i=0, ndxmax-1 do begin
   for j=0, ndymax-1 do begin
      for k=0, ndzmax-1 do begin
         rad=sqrt(x(i)^2+y(j)^2+z(k)^2)
         if(rad ge 1.d0) then sline3d_dilfac(i,j,k)=dilfac(1.d0,rad)
      endfor
   endfor
endfor
;
;define radius
r3d=fltarr(ndxmax,ndymax,ndzmax)
for i=0, ndxmax-1 do begin
   for j=0, ndymax-1 do begin
      for k=0, ndzmax-1 do begin
         case mask3d(i,j,k) of
            1: begin
                  r3d(i,j,k)=sqrt(x(i)^2+y(j)^2+z(k)^2)
                  break
               end
            2: begin
                  r3d(i,j,k)=sqrt(x(i)^2+y(j)^2+z(k)^2)
                  break
               end
            3: begin
                  r3d(i,j,k)=sqrt(x(i)^2+y(j)^2+z(k)^2)
                  break
               end
            else: begin
                  r3d(i,j,k)=1.d0 ;dummy value
               end
         endcase
      endfor
   endfor
endfor
;
;----------------------prepare convergence behaviour--------------------
;
epsmaxl_arr=abs(epsmaxl_arr)
epsmaxc_arr=abs(epsmaxc_arr)
iternrl_arr=indgen(itmaxl)
iternrc_arr=indgen(itmaxc)
;
;-----------------------define range------------------------------------
;
xmin=min(r3d)
xmax=max(r3d)
xlim=[xmin, xmax]
;
indx=where(mask3d eq 1 or mask3d eq 2 or mask3d eq 3)
;
;----------------------------title-strings------------------------------
;
if(opt_opal eq 0) then begin
   titlestr=textoidl('k_L=') + string(kline, format='(f10.2)') + ', ' + $
            textoidl('k_c=') + string(kcont, format='(f10.2)')
endif else begin
   titlestr=textoidl('\kappa_0, \alpha=') + string(kappa0, format='(f9.5)') + $
            ', ' + string(alpha, format='(f9.5)') + $
            ', ' + textoidl('k_c=') + string(kcont, format='(f9.5)')
endelse
;
xtitlestr=textoidl('r')
;
;-----------------------------------------------------------------------
;
if(not keyword_set(windx)) then windx=0
;
if(keyword_set(oname)) then begin
   oname='ps_files/benchmark13_a.ps'
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
;---------------------continuum source function-------------------------
;
!p.multi=[6,2,3]
;
ytitlestr=textoidl('S_c*r^2/I_c')
ymin=min(scont3d(indx)*r3d(indx)^2)
ymax=max(scont3d(indx)*r3d(indx)^2)
;
ylim=[ymin, ymax]
plot, [xmin,xmin], [ymin,ymin], $
      xrange=xlim, $
      yrange=ylim, $
      title=titleStr, $
      xtitle=xtitleStr, $
      ytitle=ytitleStr, $
      charsize=2., /xs, /ys
;
oplot, r3d, scont3d*r3d^2, $
      color=ci_blue, $
      psym=3
;
;oplot for positive and negative x axis
ix0=0
ix1=ndxmax/2
ix2=ndxmax-1
iy=ndymax/2
iz=ndzmax/2
oplot, abs(x(ix0:ix1)), scont3d(ix0:ix1,iy,iz)*x(ix0:ix1)^2, color=ci_green, line=0
oplot, x(ix1:ix2), scont3d(ix1:ix2,iy,iz)*x(ix1:ix2)^2, color=ci_green, psym=1
;
;oplot for positive and negative y axis
ix=ndxmax/2
iy0=0
iy1=ndymax/2
iy2=ndymax-1
iz=ndzmax/2
oplot, abs(y(iy0:iy1)), scont3d(ix,iy0:iy1,iz)*y(iy0:iy1)^2, color=ci_magenta, line=2
oplot, y(iy1:iy2), scont3d(ix,iy1:iy2,iz)*y(iy1:iy2)^2, color=ci_magenta, psym=1
;
;oplot for positive and negative z-directions
ix=ndxmax/2
iy=ndymax/2
iz0=0
iz1=ndzmax/2
iz2=ndzmax-1
oplot, abs(z(iz0:iz1)), scont3d(ix,iy,iz0:iz1)*z(iz0:iz1)^2, color=ci_red, line=3
oplot, z(iz1:iz2), scont3d(ix,iy,iz1:iz2)*z(iz1:iz2)^2, color=ci_red, psym=1
;
;------------------------line source function---------------------------
;
!p.multi=[5,2,3]
;
ytitlestr=textoidl('S_L*r^2/I_c')
ymin=min(sline3d(indx)*r3d(indx)^2)
ymax=max(sline3d(indx)*r3d(indx)^2)
ylim=[ymin, ymax]
plot, [xmin,xmin], [ymin,ymin], $
      xrange=xlim, $
      yrange=ylim, $
      title=titleStr, $
      xtitle=xtitleStr, $
      ytitle=ytitleStr, $
      charsize=2., /xs, /ys
;
oplot, r3d, sline3d*r3d^2, $
      color=ci_blue, $
      psym=3
;
;oplot for positive and negative x axis
ix0=0
ix1=ndxmax/2
ix2=ndxmax-1
iy=ndymax/2
iz=ndzmax/2
oplot, abs(x(ix0:ix1)), sline3d(ix0:ix1,iy,iz)*x(ix0:ix1)^2, color=ci_green, line=0
oplot, x(ix1:ix2), sline3d(ix1:ix2,iy,iz)*x(ix1:ix2)^2, color=ci_green, psym=1
;
;oplot for positive and negative y axis
ix=ndxmax/2
iy0=0
iy1=ndymax/2
iy2=ndymax-1
iz=ndzmax/2
oplot, abs(y(iy0:iy1)), sline3d(ix,iy0:iy1,iz)*y(iy0:iy1)^2, color=ci_magenta, line=2
oplot, y(iy1:iy2), sline3d(ix,iy1:iy2,iz)*y(iy1:iy2)^2, color=ci_magenta, psym=1
;
;oplot for positive and negative z-directions
ix=ndxmax/2
iy=ndymax/2
iz0=0
iz1=ndzmax/2
iz2=ndzmax-1
oplot, abs(z(iz0:iz1)), sline3d(ix,iy,iz0:iz1)*z(iz0:iz1)^2, color=ci_red, line=3
oplot, z(iz1:iz2), sline3d(ix,iy,iz1:iz2)*z(iz1:iz2)^2, color=ci_red, psym=1
;
;----------------------convergence behaviour continuum------------------
;
!p.multi=[4,2,3]
;
xmin=0
xmax=max(where(epsmaxc_arr gt 0))
if(xmax eq -1) then begin
   ymin=1.d-5
endif else begin
   ymin=epsmaxc_arr(xmax)
endelse
ymax=max(epsmaxc_arr)
;
ylim=[ymin,ymax]
;
ytitleStr=textoidl('((S_c^{(k-1)} - S_c^{(k)}) /S_c^{(k)})_{max}')
xtitleStr='# iterations'
;
plot, iternrc_arr, epsmaxc_arr, $
      /ylog, $
      line=0, $
      yrange=ylim, $
      xrange=[xmin,xmax], $
      title=titlestr, $
      xtitle=xtitlestr, $
      ytitle=ytitlestr, $
      charsize=2.
;
;----------------------convergence behaviour line-----------------------
;
!p.multi=[3,2,3]
;
xmin=0
xmax=max(where(epsmaxl_arr gt 0))
ymin=epsmaxl_arr(xmax)
ymax=max(epsmaxl_arr)
;
ylim=[ymin,ymax]
;
ytitleStr=textoidl('((S_L^{(k-1)} - S_L^{(k)}) /S_L^{(k)})_{max}')
xtitleStr='# iterations'
;
plot, iternrl_arr, epsmaxl_arr, $
      /ylog, $
      line=0, $
      yrange=ylim, $
      xrange=[xmin,xmax], $
      title=titlestr, $
      xtitle=xtitlestr, $
      ytitle=ytitlestr, $
      charsize=2.
;
;---------------------------temperature model---------------------------
;
ymin=min(t3d(indx))
ymax=max(t3d(indx))
ylim=[ymin-1.d3, ymax+1.d3]
!p.multi=[3,3,3]
plot, [0.,0.], [0.,0.], $
      xtitle=textoidl('r'), $
      ytitle=textoidl('T'), $
      title=titlestr, $
      charsize=2., $
      xrange=xlim, $
      yrange=ylim, /xs, /ys
oplot, r3d, t3d, $
       psym=3, $
       color=ci_blue
;
;---------------------continuum opacity model---------------------------
;
ymin=min(opac3d(indx))
ymax=max(opac3d(indx))
ylim=[ymin,ymax]
;
!p.multi=[2,3,3]
plot, [0.,0.], [0.,0.], $
      xtitle=textoidl('r'), $
      ytitle=textoidl('\chi_c'), $
      title=titlestr, $
      xrange=xlim, $
      yrange=ylim, $
      /xs, /ys, $
      /ylog, $
      charsize=2.
oplot, r3d, opac3d, $
       psym=3, $
       color=ci_blue
;
;---------------------------line opacity model--------------------------
;
ymin=min(opalbar3d(indx))
ymax=max(opalbar3d(indx))
ylim=[ymin,ymax]
;
!p.multi=[1,3,3]
plot, [0.,0.], [0.,0.], $
      xtitle=textoidl('r'), $
      ytitle=textoidl('\chi bar'), $
      title=titlestr, $
      xrange=xlim, $
      yrange=ylim, $
      /xs, /ys, $
      /ylog, $
      charsize=2.
oplot, r3d, opalbar3d, $
       psym=3, $
       color=ci_blue
;
;
;
if keyword_set(oname) then begin
   device, /close
   set_plot,'x'
endif
;
!p.multi=0
;
;------------------------mixed plots------------------------------------
;
window, windx
device, decomposed=0
;
ytitlestr=textoidl('S_L/WI_c')
ymin=min(sline3d(indx)/sline3d_dilfac(indx))
ymax=max(sline3d(indx)/sline3d_dilfac(indx))
ylim=[ymin, ymax]
plot, [xmin,xmin], [ymin,ymin], $
      xrange=xlim, $
      yrange=ylim, $
      title=titleStr, $
      xtitle=xtitleStr, $
      ytitle=ytitleStr, $
      charsize=2., /xs, /ys
;
oplot, r3d, sline3d/sline3d_dilfac, $
      color=ci_blue, $
      psym=3
;
;
end
