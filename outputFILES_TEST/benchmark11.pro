pro benchmark11, windx=WINDX, oname=ONAME, xlim=XLIM, ylim=YLIM, dir=dir, version=version
;
;plots contours of searchlight along a direction
;
loadct, 0
;
;------------------read all information from hdf5-file------------------
;
if(not keyword_set(dir)) then dir='.'
fname=dir+'/benchmark11.h5'

if(not keyword_set(version)) then version=2

if(version eq 1) then begin
;
   file_id = h5f_open(fname)
;
      group_id = h5g_open(file_id, 'options')
         att_id=h5a_open_name(group_id, 'opt_angint_method')
            opt_angint_method=h5a_read(att_id)
            opt_angint_method=opt_angint_method(0)
         h5a_close, att_id
      h5g_close, group_id
   
      group_id = h5g_open(file_id, 'dimensions')
         att_id=h5a_open_name(group_id, 'ndxmax')
            ndxmax=h5a_read(att_id)
            ndxmax=ndxmax(0)
         h5a_close, att_id
         att_id=h5a_open_name(group_id, 'ndymax')
            ndymax=h5a_read(att_id)
            ndymax=ndymax(0)
         h5a_close, att_id
         att_id=h5a_open_name(group_id, 'ndzmax')
            ndzmax=h5a_read(att_id)
            ndzmax=ndzmax(0)
         h5a_close, att_id
         att_id=h5a_open_name(group_id, 'dim_mu')
            dim_mu=h5a_read(att_id)
            dim_mu=dim_mu(0)
         h5a_close, att_id
         att_id=h5a_open_name(group_id, 'dim_phi')
            dim_phi=h5a_read(att_id)
            dim_phi=dim_phi(0)
         h5a_close, att_id
         att_id=h5a_open_name(group_id, 'n1d_angdep')
            n1d_angdep=h5a_read(att_id)
            n1d_angdep=n1d_angdep(0)
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
         dset_id=h5d_open(group_id, 'xcoord_angdep')
            xcoord_angdep=h5d_read(dset_id)
         h5d_close, dset_id
         dset_id=h5d_open(group_id, 'ycoord_angdep')
            ycoord_angdep=h5d_read(dset_id)
         h5d_close, dset_id
         dset_id=h5d_open(group_id, 'zcoord_angdep')
            zcoord_angdep=h5d_read(dset_id)
         h5d_close, dset_id
      h5g_close, group_id
;
      group_id = h5g_open(file_id, 'angles')
         dset_id=h5d_open(group_id, 'nodes_mu')
            nodes_mu=h5d_read(dset_id)
         h5d_close, dset_id
         dset_id=h5d_open(group_id, 'weight_mu')
            weight_mu=h5d_read(dset_id)
         h5d_close, dset_id
         dset_id=h5d_open(group_id, 'nodes_phi_pair')
            nodes_phi_pair=h5d_read(dset_id)
         h5d_close, dset_id
         dset_id=h5d_open(group_id, 'weight_phi_pair')
            weight_phi_pair=h5d_read(dset_id)
         h5d_close, dset_id
         dset_id=h5d_open(group_id, 'phi_mask')
            phi_mask=h5d_read(dset_id)
         h5d_close, dset_id
      h5g_close, group_id
   ;
      group_id = h5g_open(file_id, 'solution3d')
         dset_id=h5d_open(group_id, 'mint3d_sc')
            mint3d_sc=h5d_read(dset_id)
         h5d_close, dset_id
         dset_id=h5d_open(group_id, 'mint3d_fvm')
            mint3d_fvm=h5d_read(dset_id)
         h5d_close, dset_id
         dset_id=h5d_open(group_id, 'mint3d_theo')
            mint3d_theo=h5d_read(dset_id)
         h5d_close, dset_id
         dset_id=h5d_open(group_id, 'mask3d')
            mask3d=h5d_read(dset_id)
         h5d_close, dset_id
      h5g_close, group_id
   ;
      group_id = h5g_open(file_id, 'solution_angdep')
         dset_id=h5d_open(group_id, 'intsc_angdep')
            intsc_angdep=h5d_read(dset_id)
         h5d_close, dset_id
         dset_id=h5d_open(group_id, 'intfvm_angdep')
            intfvm_angdep=h5d_read(dset_id)
         h5d_close, dset_id
      h5g_close, group_id
;
   h5f_close, file_id
endif else begin
;
   file_id = h5f_open(fname)
;
      group_id = h5g_open(file_id, 'options')
         att_id=h5a_open_name(group_id, 'opt_angint_method')
            opt_angint_method=h5a_read(att_id)
            opt_angint_method=opt_angint_method(0)
         h5a_close, att_id
      h5g_close, group_id
   
      group_id = h5g_open(file_id, 'dimensions')
         att_id=h5a_open_name(group_id, 'ndxmax')
            ndxmax=h5a_read(att_id)
            ndxmax=ndxmax(0)
         h5a_close, att_id
         att_id=h5a_open_name(group_id, 'ndymax')
            ndymax=h5a_read(att_id)
            ndymax=ndymax(0)
         h5a_close, att_id
         att_id=h5a_open_name(group_id, 'ndzmax')
            ndzmax=h5a_read(att_id)
            ndzmax=ndzmax(0)
         h5a_close, att_id
         att_id=h5a_open_name(group_id, 'dim_omega')
            dim_omega=h5a_read(att_id)
            dim_omega=dim_omega(0)
         h5a_close, att_id
         att_id=h5a_open_name(group_id, 'n1d_angdep')
            n1d_angdep=h5a_read(att_id)
            n1d_angdep=n1d_angdep(0)
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
         dset_id=h5d_open(group_id, 'xcoord_angdep')
            xcoord_angdep=h5d_read(dset_id)
         h5d_close, dset_id
         dset_id=h5d_open(group_id, 'ycoord_angdep')
            ycoord_angdep=h5d_read(dset_id)
         h5d_close, dset_id
         dset_id=h5d_open(group_id, 'zcoord_angdep')
            zcoord_angdep=h5d_read(dset_id)
         h5d_close, dset_id
      h5g_close, group_id
;
      group_id = h5g_open(file_id, 'angles')
         dset_id=h5d_open(group_id, 'n_x')
            n_x=h5d_read(dset_id)
         h5d_close, dset_id
         dset_id=h5d_open(group_id, 'n_y')
            n_y=h5d_read(dset_id)
         h5d_close, dset_id
         dset_id=h5d_open(group_id, 'n_z')
            n_z=h5d_read(dset_id)
         h5d_close, dset_id
         dset_id=h5d_open(group_id, 'weight_omega')
            weight_omega=h5d_read(dset_id)
         h5d_close, dset_id
      h5g_close, group_id
   ;
      group_id = h5g_open(file_id, 'solution3d')
         dset_id=h5d_open(group_id, 'mint3d_sc')
            mint3d_sc=h5d_read(dset_id)
         h5d_close, dset_id
         dset_id=h5d_open(group_id, 'mint3d_fvm')
            mint3d_fvm=h5d_read(dset_id)
         h5d_close, dset_id
         dset_id=h5d_open(group_id, 'mint3d_theo')
            mint3d_theo=h5d_read(dset_id)
         h5d_close, dset_id
         dset_id=h5d_open(group_id, 'mask3d')
            mask3d=h5d_read(dset_id)
         h5d_close, dset_id
      h5g_close, group_id
   ;
      group_id = h5g_open(file_id, 'solution_angdep')
         dset_id=h5d_open(group_id, 'intsc_angdep')
            intsc_angdep=h5d_read(dset_id)
         h5d_close, dset_id
         dset_id=h5d_open(group_id, 'intfvm_angdep')
            intfvm_angdep=h5d_read(dset_id)
         h5d_close, dset_id
      h5g_close, group_id
   ;
   h5f_close, file_id
;
   nodes_mu=fltarr(dim_omega)
   nodes_phi=fltarr(dim_omega)
   for i=0, dim_omega-1 do begin
      get_angles_spc, n_x(i), n_y(i), n_z(i), theta, phi
      nodes_mu(i)=cos(theta)
      nodes_phi(i)=phi
   endfor
endelse


;
print, '------errors for sc method------'
errors3d, mint3d_sc, mint3d_theo, mask3d, erri, errm, err_max, devm
print, ''
print, '-----errors for fvm method------'
errors3d, mint3d_fvm, mint3d_theo, mask3d, erri, errm, err_max, devm
;
;-----------------------------------------------------------------------
;
;********************plot mean intensities******************************
;
;-----------------------define range------------------------------------
;
;define radius
r3d=fltarr(ndxmax,ndymax,ndzmax)
for i=0, ndxmax-1 do begin
   for j=0, ndymax-1 do begin
      for k=0, ndzmax-1 do begin
         case mask3d(i,j,k) of
            1: r3d(i,j,k)=sqrt(x(i)^2+y(j)^2+z(k)^2)
            2: r3d(i,j,k)=sqrt(x(i)^2+y(j)^2+z(k)^2)
            3: r3d(i,j,k)=sqrt(x(i)^2+y(j)^2+z(k)^2)
            else: r3d(i,j,k)=0.d0 ;dummy value
         endcase

;         dum=mint3d_sc(i,j,k)*r3d(i,j,k)^2
;         if(dum gt 0. and dum lt 0.1d0) then print, i, j, k, r3d(i,j,k), x(i), y(j), z(k)
      endfor
   endfor
endfor
;
;
if(not keyword_set(xlim)) then begin
   xlim=[1.d0, 12.d0]
endif
;
if(not keyword_set(ylim)) then begin
   ylim=[0.2d0,0.6d0]
   ylim=[-0.05d0,0.6d0]
endif
;
;----------------------------title-strings------------------------------
;
;titlestr=textoidl('i (\theta=') + string(acos(mu), format='(f9.5)') + $
;         textoidl(', \phi=') + string(phi, format='(f9.5)') + ')'
titlestr=''
xtitlestr=textoidl('r')
ytitlestr=textoidl('J*r^2')
;
;-----------------------------------------------------------------------
;
if(not keyword_set(windx)) then windx=0
;
if(keyword_set(oname)) then begin
   oname='ps_files/benchmark11_a.ps'
   print, "writing output to: ", oname
   set_plot,'ps'
   device,file=oname, $ ;xsize=19., ysize=26.7, xoffset=1., yoffset=1. , $
          color=1, bits_per_pixel=8
endif else begin
   window, windx
   device, decomposed=0
   windx=windx+1
endelse
;
loadct_lev, ci_black, ci_blue, ci_red, ci_cyan, ci_green, ci_magenta, ci_yellow, ci_white
;
!p.multi=[2,1,2]
;
plot, r3d, mint3d_theo*r3d^2, $
      xrange=xlim, $
      yrange=ylim, $
      title=titleStr, $
      xtitle=xtitleStr, $
      ytitle=ytitleStr, $
      psym=3

oplot, r3d, mint3d_sc*r3d^2, $
      color=ci_blue, $
      psym=3

oplot, r3d, mint3d_fvm*r3d^2, $
      color=ci_red, $
      psym=3

;plot, r3d, mint3d_theo, xrange=xlim, yrange=[1.d-4,1.d0], /ylog
;oplot, r3d, mint3d_sc, color=ci_blue, psym=3
;oplot, r3d, mint3d_fvm, color=ci_red, psym=3
;stop

;for i=ndxmax/2, ndxmax-1 do begin
;   oplot, r3d(i,ndymax/2,*), mint3d_sc(i,ndymax/2,*)*r3d(i,ndymax/2,*)^2, color=ci_green
;   oplot, r3d(i,ndymax/2,*), mint3d_sc(i,ndymax/2,*)*r3d(i,ndymax/2,*)^2, color=ci_green, psym=1
;   oplot, r3d(*,ndymax/2,i), mint3d_sc(*,ndymax/2,i)*r3d(*,ndymax/2,i)^2, color=ci_yellow
;   oplot, r3d(*,ndymax/2,i), mint3d_sc(*,ndymax/2,i)*r3d(*,ndymax/2,i)^2, color=ci_yellow, psym=1
;   oplot, r3d(ndxmax/2,*,i), mint3d_sc(ndxmax/2,*,i)*r3d(ndxmax/2,*,i)^2, color=ci_magenta
;   oplot, r3d(ndxmax/2,*,i), mint3d_sc(ndxmax/2,*,i)*r3d(ndxmax/2,*,i)^2, color=ci_magenta, psym=1
;   print, x(i), y(ndymax/2), z(i)
;   read, s
;   if(s eq 0) then break
;endfor

;
lstr1='theoretical (from dilution)'
lstr2='3d short characteristics'
lstr3='3d finite volume method'
if(keyword_set(oname)) then begin
   legend, [lstr1, lstr2, lstr3], $
         psym=[3,3,3], $
         linestyle=[0,0,0], $
         color=[ci_black,ci_blue,ci_red], $
         textcolor=ci_black, $
         /right_legend
endif else begin
   legend, [lstr1, lstr2, lstr3], $
         psym=[3,3,3], $
         linestyle=[0,0,0], $
         color=[ci_white,ci_blue,ci_red], $
         textcolor=ci_white, $
         /right_legend
endelse
;
!p.multi=[1,1,2]
plot, [0.,12.d0],[1.,1.], yrange=[0.,2.],xrange=[0.,12.d0], $
      line=1, $
      xtitle='r', $
      ytitle='J(num)/J(theo)'
oplot, r3d, mint3d_sc/mint3d_theo, $
       color=ci_blue, $
       psym=3
oplot, r3d, mint3d_fvm/mint3d_theo, $
       color=ci_red, $
       psym=3
lstr1='3d short characteristics'
lstr2='3d finite volume method'
if(keyword_set(oname)) then begin
   legend, [lstr1, lstr2], $
         psym=[3,3], $
         linestyle=[0,0], $
         color=[ci_blue,ci_red], $
         textcolor=ci_black, $
         /right_legend
endif else begin
   legend, [lstr1, lstr2], $
         psym=[3,3], $
         linestyle=[0,0], $
         color=[ci_blue,ci_red], $
         textcolor=ci_white, $
         /right_legend
endelse
;
if keyword_set(oname) then begin
   device, /close
   set_plot,'x'
endif
;
!p.multi=0
;
;********plot mean intensities at a certain radius on a sphere**********
;
if(keyword_set(oname)) then begin
   oname='ps_files/benchmark11_b.ps'
   print, "writing output to: ", oname
   set_plot,'ps'
   device,file=oname, color=1, bits_per_pixel=8
endif else begin
   window, windx
   device, decomposed=0
   windx=windx+1
endelse
;
plot_mint3d_map, ndxmax, ndymax, ndzmax, x, y, z, mint3d_sc, mint3d_fvm, radius=3.d0, clim=[0.95d0,1.5d0];clim=[0.95d0,1.15d0]
;
if keyword_set(oname) then begin
   device, /close
   set_plot,'x'
endif
;
;********************plot intensities as function of angle**************
;
;-----------------------define range------------------------------------
;
if(keyword_set(oname)) then begin
   oname='ps_files/benchmark11_c.ps'
   print, "writing output to: ", oname
   set_plot,'ps'
   device,file=oname, $ ;xsize=19., ysize=26.7, xoffset=1., yoffset=1. , $
   color=1, bits_per_pixel=8
endif else begin
   window, windx
   device, decomposed=0
   windx=windx+1
endelse
;
nline=ceil(sqrt(n1d_angdep))
ncols=floor(sqrt(n1d_angdep))
;
print, 'x   ', 'y   ', 'z   ', 'J_sc *r^2   ', 'J_theo * r^2   ', 'J_sc/J_theo   '
for i=0, n1d_angdep-1 do begin
;for i=0, 0 do begin

;
;perform integration
;   intsc2d=fltarr(dim_omega)*0.d0
;   inttheo2d=fltarr(dim_omega)*0.d0
;   jnum=0.d0
;   for k=0, dim_omega-1 do begin
;      muc, xcoord_angdep(i), ycoord_angdep(i), zcoord_angdep(i), nodes_phi(k), muc1, muc2
;      for j=0,dim_mu-1 do begin
;         intsc2d(j,k)=intsc_angdep(i,j,k)
;         if(nodes_mu(j) ge muc1 and nodes_mu(j) le muc2) then begin
;            inttheo2d(j,k)=1.d0
;            jnum = jnum+weight_mu(j)*weight_phi_pair(j,k)
;         endif
;      endfor
;   endfor
;   rad=sqrt(xcoord_angdep(i)^2+ycoord_angdep(i)^2+zcoord_angdep(i)^2)

;   find_indx, xcoord_angdep(i), x, ndxmax, iim1, ii
;   find_indx, ycoord_angdep(i), y, ndymax, jjm1, jj
;   find_indx, zcoord_angdep(i), z, ndzmax, kkm1, kk
;;   jnum = integ_trapez2d(nodes_mu, nodes_phi, intsc2d)/4.d0/!pi
;   jtheo = 0.5d0*(1.d0-sqrt(1.d0-1.d0/rad^2))
;   print, xcoord_angdep(i), ycoord_angdep(i), zcoord_angdep(i), jnum*rad^2, jtheo*rad^2, jnum/jtheo, mint3d_sc(ii,jj,kk)/jnum
;
;   plot_intangdep, xcoord_angdep(i), ycoord_angdep(i), zcoord_angdep(i),
;   nodes_mu, nodes_phi, intsc_angdep(i,*,*), intfvm_angdep(i,*,*), phi_mask,
;   oname=oname, opt_angint_method=opt_angint_method
   plot_intangdep2, xcoord_angdep(i), ycoord_angdep(i), zcoord_angdep(i), nodes_mu, nodes_phi, intsc_angdep(i,*), intfvm_angdep(i,*), oname=oname
   read, enter
endfor
;
if keyword_set(oname) then begin
   device, /close
   set_plot,'x'
endif



;
!p.multi=0


end
;
;-----------------------------------------------------------------------
;-----------------------------------------------------------------------
;-----------------------------------------------------------------------
;
pro plot_intangdep, xp, yp, zp, nodes_mu, nodes_phi, intsc_angdep, intfvm_angdep, phi_mask, oname=oname, opt_angint_method=opt_angint_method
;
;-------------------define mu_crit on each phi-level--------------------
;
np=1000
mu_crit1=fltarr(np)*0.d0
mu_crit2=fltarr(np)*0.d0
;
phi=findgen(np)*2.d0*!pi/(np-1)
;
for i=0, np-1 do begin
   muc, xp, yp, zp, phi(i), muc1, muc2
   mu_crit1(i) = muc1
   mu_crit2(i) = muc2
endfor
;
;-----------------------define range------------------------------------
;
xlim=[-1.d0,1.d0]
ylim=[0.d0,2.d0*!pi]
;
zmin=min([min(intsc_angdep),min(intfvm_angdep)])
zmax=max([max(intsc_angdep),max(intfvm_angdep)])
;
dz=zmax-zmin
if(dz eq 0.) then begin
   zmin=zmin-0.1d0
   zmax=zmax+0.1d0
endif else begin
   zmin=zmin-0.1d0*dz
   zmax=zmax+0.1d0*dz
endelse
zlim=[zmin,zmax]
;
;------------------if triangles are used: irregular arrays--------------
;
dim_mu=n_elements(nodes_mu)
dim_phi=n_elements(nodes_phi)
;
ntot=0L
for i=0, dim_mu-1 do begin
   for j=0, dim_phi-1 do begin
      if(phi_mask(i,j) eq 1) then ntot=ntot+1
   endfor
endfor
;
nodes_mu2=fltarr(ntot)
nodes_phi2=fltarr(ntot)
intsc_angdep2=fltarr(ntot)
intfvm_angdep2=fltarr(ntot)
;
k=0
for i=0, dim_mu-1 do begin
   for j=0, dim_phi-1 do begin
      if(phi_mask(i,j) eq 1) then begin
         nodes_mu2(k)=nodes_mu(i)
         nodes_phi2(k)=nodes_phi(j)
         intsc_angdep2(k)=intsc_angdep(0,i,j)
         intfvm_angdep2(k)=intfvm_angdep(0,i,j)
         k=k+1
      endif
   endfor
endfor

;for debugging
;window, 0
;plot, [0.,0.],[0.,0.], xrange=[-1.d0,1.d0], yrange=[0.d0,2.d0*!pi]
;for i=0, dim_mu-1 do begin
;   for j=0, dim_phi-1 do begin
;      if(phi_mask(i,j) eq 1) then begin
;         oplot, [nodes_mu(i),nodes_mu(i)], [nodes_phi(j),nodes_phi(j)], psym=1
;      endif
;   endfor
;endfor
;
;nodes_mu2=[-1.d0, -0.707106781186547d0, 0.d0, 0.707106781186548d0, 1.d0]
;nodes_phi2=[0.d0, 0.785398163397448d0, 1.57079632679490d0, 2.35619449019234d0, $
;            3.14159265358979d0, 3.92699081698724d0, 4.71238898038469, 5.49778714378214, 6.28318530717959]
;for i=0, 4 do begin
;   oplot, [nodes_mu2(i),nodes_mu2(i)], [0.d0, 2.d0*!pi]
;endfor
;for i=0, 8 do begin
;   oplot, [-1.d0, 1.d0], [nodes_phi2(i),nodes_phi2(i)]
;endfor
;stop


;
;----------------------------title-strings------------------------------
;
;titlestr=textoidl('i (\theta=') + string(acos(mu), format='(f9.5)') + $
;         textoidl(', \phi=') + string(phi, format='(f9.5)') + ')'
titlestr='at point (x,y,z)=' + string(xp, format='(f6.3)') + ',' $
                             + string(yp, format='(f6.3)') + ',' $
                             + string(zp, format='(f6.3)')
titlestr1='SC ' + titlestr
titlestr2='FVM ' + titlestr

xtitlestr=textoidl('\mu')
ytitlestr=textoidl('\phi')
ztitlestr=textoidl('I(\mu,\phi)')
;
;---------------------for surface plot----------------------------------
;
;loadct_lev, ci_black, ci_blue, ci_red, ci_cyan, ci_green, ci_magenta, ci_yellow, ci_white
;
;surface, intfvm_angdep, nodes_mu, nodes_phi, skirt=0.d0, charsize=3., zrange=zlim, xrange=xlim, yrange=ylim, $
;         xtitle=xtitleStr, ytitle=ytitleStr, ztitle=ztitleStr, title=titlestr
;
;----------------------for contour plot---------------------------------
;
;color range
min_color = 0.01d0
max_color = 1.d0
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
get_contour_colors, colors, nodes_mu, nodes_phi, ncolors, nlevels_iso, bottom, $
                    levels_final, levels_iso, c_colors_final,  $
                    cb_indx, cb_tickmark_name
;
if(opt_angint_method eq 6 or opt_angint_method eq 7) then begin
;irregular grid (triangles)
   contourplots_double, intsc_angdep2, intfvm_angdep2, /irregular, $
                        nodes_mu2, nodes_mu2, nodes_phi2, nodes_phi2, $
                        ncolors, nlevels_iso, bottom, $
                        levels_final, levels_iso, c_colors_final, $
                        cb_indx, cb_tickmark_name, $
                        xlim=xlim, $
                        ylim=ylim, $
                        titlestr1=titlestr1, $
                        titlestr2=titlestr2, $
                        xtitlestr=xtitleStr, $
                        ytitlestr=ytitlestr, $
                        ctitlestr=ztitlestr, $
                        /oplot_grid, $
                        ctable=13, $
                        oplot1_x1=mu_crit1, oplot1_y1=phi, $
                        oplot1_x2=mu_crit2, oplot1_y2=phi, $
                        oplot2_x1=mu_crit1, oplot2_y1=phi, $
                        oplot2_x2=mu_crit2, oplot2_y2=phi, $
                        oplot1_style=1, oplot2_style=1, $
                        /background1
endif else begin
;
   contourplots_double, intsc_angdep, intfvm_angdep, $
                        nodes_mu, nodes_mu, nodes_phi, nodes_phi, $
                        ncolors, nlevels_iso, bottom, $
                        levels_final, levels_iso, c_colors_final, $
                        cb_indx, cb_tickmark_name, $
                        xlim=xlim, $
                        ylim=ylim, $
                        titlestr1=titlestr1, $
                        titlestr2=titlestr2, $
                        xtitlestr=xtitleStr, $
                        ytitlestr=ytitlestr, $
                        ctitlestr=ztitlestr, $
                        /oplot_grid, $
                        ctable=13, $
                        oplot1_x1=mu_crit1, oplot1_y1=phi, $
                        oplot1_x2=mu_crit2, oplot1_y2=phi, $
                        oplot2_x1=mu_crit1, oplot2_y1=phi, $
                        oplot2_x2=mu_crit2, oplot2_y2=phi, $
                        oplot1_style=1, oplot2_style=1, $
                        /background1
;
endelse
;
loadct_lev, ci_black, ci_blue, ci_red, ci_cyan, ci_green, ci_magenta, ci_yellow, ci_white
oplot, mu_crit1, phi, psym=1, thick=4., color=ci_red
oplot, mu_crit2, phi, psym=1, thick=4., color=ci_red
;
;
;
end
;
;-----------------------------------------------------------------------
;-----------------------------------------------------------------------
;-----------------------------------------------------------------------
;
pro plot_intangdep2, xp, yp, zp, nodes_mu, nodes_phi, intsc_angdep, intfvm_angdep, oname=oname
;
;for irregular mu,phi spacing
;
;-------------------define mu_crit on each phi-level--------------------
;
np=1000
mu_crit1=fltarr(np)*0.d0
mu_crit2=fltarr(np)*0.d0
;
phi=findgen(np)*2.d0*!pi/(np-1)
;
for i=0, np-1 do begin
   muc, xp, yp, zp, phi(i), muc1, muc2
   mu_crit1(i) = muc1
   mu_crit2(i) = muc2
endfor
;
;-----------------------define range------------------------------------
;
xlim=[-1.d0,1.d0]
ylim=[0.d0,2.d0*!pi]
;
zmin=min([min(intsc_angdep),min(intfvm_angdep)])
zmax=max([max(intsc_angdep),max(intfvm_angdep)])
;
dz=zmax-zmin
if(dz eq 0.) then begin
   zmin=zmin-0.1d0
   zmax=zmax+0.1d0
endif else begin
   zmin=zmin-0.1d0*dz
   zmax=zmax+0.1d0*dz
endelse
zlim=[zmin,zmax]
;
;----------------------------title-strings------------------------------
;
;titlestr=textoidl('i (\theta=') + string(acos(mu), format='(f9.5)') + $
;         textoidl(', \phi=') + string(phi, format='(f9.5)') + ')'
titlestr='at point (x,y,z)=' + string(xp, format='(f6.3)') + ',' $
                             + string(yp, format='(f6.3)') + ',' $
                             + string(zp, format='(f6.3)')
titlestr1='SC ' + titlestr
titlestr2='FVM ' + titlestr
titlestr1=''
titlestr2=''

xtitlestr=textoidl('sin(\theta)')
ytitlestr=textoidl('\phi')
ztitlestr=textoidl('I(\mu,\phi)')
;
;---------------------for surface plot----------------------------------
;
;loadct_lev, ci_black, ci_blue, ci_red, ci_cyan, ci_green, ci_magenta, ci_yellow, ci_white
;
;surface, intfvm_angdep, nodes_mu, nodes_phi, skirt=0.d0, charsize=3., zrange=zlim, xrange=xlim, yrange=ylim, $
;         xtitle=xtitleStr, ytitle=ytitleStr, ztitle=ztitleStr, title=titlestr
;
;----------------------for contour plot---------------------------------
;
;color range
min_color = 0.01d0
max_color = 1.d0
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
get_contour_colors, colors, nodes_mu, nodes_phi, ncolors, nlevels_iso, bottom, $
                    levels_final, levels_iso, c_colors_final,  $
                    cb_indx, cb_tickmark_name
;
contourplots_double, intsc_angdep, intfvm_angdep, /irregular, $
                     nodes_mu, nodes_mu, nodes_phi, nodes_phi, $
                     ncolors, nlevels_iso, bottom, $
                     levels_final, levels_iso, c_colors_final, $
                     cb_indx, cb_tickmark_name, $
                     xlim=xlim, $
                     ylim=ylim, $
                     titlestr1=titlestr1, $
                     titlestr2=titlestr2, $
                     xtitlestr=xtitleStr, $
                     ytitlestr=ytitlestr, $
                     ctitlestr=ztitlestr, $
                     /oplot_grid, $
                     ctable=13, $
;                     oplot1_x1=mu_crit1, oplot1_y1=phi, $
;                     oplot1_x2=mu_crit2, oplot1_y2=phi, $
;                     oplot2_x1=mu_crit1, oplot2_y1=phi, $
;                     oplot2_x2=mu_crit2, oplot2_y2=phi, $
;                     oplot1_style=1, oplot2_style=1, $
                     /background1
;
loadct_lev, ci_black, ci_blue, ci_red, ci_cyan, ci_green, ci_magenta, ci_yellow, ci_white
;oplot, mu_crit1, phi, psym=1, thick=4., color=ci_red
;oplot, mu_crit2, phi, psym=1, thick=4., color=ci_red
;
;
;
end
;
;-----------------------------------------------------------------------
;-----------------------------------------------------------------------
;-----------------------------------------------------------------------
;
pro plot_mint3d_map, ndxmax, ndymax, ndzmax, x, y, z, mint3d_sc, mint3d_fvm, radius=radius, clim=clim
;
;-----------------------define angular grids----------------------------
;
if(not keyword_set(radius)) then radius=1.5d0
;
n_theta=31
n_phi=2*n_theta-1
;
tht_min=0.d0
tht_max=!pi
;
phi_min=0.d0
phi_max=2.d0*!pi
;
mu_grid= tht_min+findgen(n_theta)*(tht_max-tht_min)/(n_theta-1)
phi_grid=phi_min+findgen(n_phi)*(phi_max-phi_min)/(n_phi-1)
;
;-----------------------------------------------------------------------
;
mint2d_sc=fltarr(n_phi, n_theta)*1.d0
mint2d_fvm=fltarr(n_phi, n_theta)*1.d0
;
;------interpolate mint3d onto sphere with radius and angle-grids------
;
pp=fltarr(3)*0.d0
for i=0, n_theta - 1 do begin
   for j=0, n_phi-1 do begin
      pp(0)=radius*sin(mu_grid(i))*cos(phi_grid(j))
      pp(1)=radius*sin(mu_grid(i))*sin(phi_grid(j))
      pp(2)=radius*cos(mu_grid(i))
;
;find indices of cube-vertices for interpolation
      get_xyz_indx, pp(0), pp(1), pp(2), x, y, z, ndxmax, ndymax, ndzmax, indx_x1, indx_x2, indx_y1, indx_y2, indx_z1, indx_z2, expol
;store all spatial values on cube-vertices for interpolation and check if
;interpolation in logspace is allowed
      get_xyz_values1, pp(0), pp(1), pp(2), x, y, z, ndxmax, ndymax, ndzmax, $
                       indx_x1, indx_x2, indx_y1, indx_y2, indx_z1, indx_z2, $
                       x1, x2, y1, y2, z1, z2, $
                       rada, radb, radc, radd, rade, radf, radg, radh, radp, $
                       llogx, llogy, llogz
;
;------------------------for short characteristics----------------------
;
;store all physical values on cube-vertices for interpolation
      get_xyz_values2, ndxmax, ndymax, ndzmax, mint3d_sc, $
                       indx_x1, indx_x2, indx_y1, indx_y2, indx_z1, indx_z2, $
                       vala, valb, valc, vald, vale, valf, valg, valh, llogf

;following statements can switch to pure linear interpolation
;      llogx=0
;      llogy=0
;      llogz=0
;      llogf=0
;here, can decide if values shall be interpolated by function*r^2
      lr2=1
;actual interpolation
     trilin_complete, pp(0), pp(1), pp(2), $
                      x1, x2, y1, y2, z1, z2, $
                      vala, valb, valc, vald, vale, valf, valg, valh, $
                      rada, radb, radc, radd, rade, radf, radg, radh, radp, $
                      expol, llogx, llogy, llogz, llogf, lr2, valp

      mint2d_sc(j,i)=valp
;
;--------------------------------for fvm--------------------------------
;
;store all physical values on cube-vertices for interpolation
      get_xyz_values2, ndxmax, ndymax, ndzmax, mint3d_fvm, $
                       indx_x1, indx_x2, indx_y1, indx_y2, indx_z1, indx_z2, $
                       vala, valb, valc, vald, vale, valf, valg, valh, llogf

;following statements can switch to pure linear interpolation
;      llogx=0
;      llogy=0
;      llogz=0
;      llogf=0
;here, can decide if values shall be interpolated by function*r^2
      lr2=1
;actual interpolation
     trilin_complete, pp(0), pp(1), pp(2), $
                      x1, x2, y1, y2, z1, z2, $
                      vala, valb, valc, vald, vale, valf, valg, valh, $
                      rada, radb, radc, radd, rade, radf, radg, radh, radp, $
                      expol, llogx, llogy, llogz, llogf, lr2, valp

;      interpol_trilinear, mint3d, x, y, z, pp, valp

      mint2d_fvm(j,i)=valp
   endfor
endfor
;
dilfac=1.d0-sqrt(1.d0-1.d0/radius/radius)
mint_theo=0.5*dilfac
;
mint2d_sc=mint2d_sc/mint_theo
mint2d_fvm=mint2d_fvm/mint_theo
;
;-----------------------------------------------------------------------
;
;define range
rangex=[min(phi_grid),max(phi_grid)]
rangey=[min(mu_grid),max(mu_grid)]
;
;---------------PLOTS FOR POSITIVE X_OFFSET-----------------------------
;
;TITLE-STRINGS
titleStr1=textoidl('SC at r= ') + STRING(RADIUS, FORMAT='(F3.1)')
titleStr2=textoidl('FVM at r= ') + STRING(RADIUS, FORMAT='(F3.1)')
ctitleStr=textoidl('J / J_{theo}')
xtitleStr=textoidl('\phi')
ytitleStr=textoidl('\theta')
;
;*************************NOW: MAP PROJECTION***************************
;
;
;transform colatitude to latitude
mu_grid=mu_grid*180.d0/!pi - 90.d0
;
;calculate azimuth in degree
phi_grid=phi_grid*180.d0/!pi
;
;define center of map (in degree)
center_lat=45.d0
center_lon=45.d0
;
;define limits to plot
;limits=[latmin, lonmin, latmax, lonmax]
plimit=[-90.d0,0.d0,90.d0,360.d0]
;plimit=[30.d0, 30.d0, 60.d0, 60.d0]
;
;
;--------------------DEFINE COLOR RANGES--------------------------------
;
if(not keyword_set(clim)) then begin
   min_color = min(mint2d(where(mint2d gt 0.)))
   max_color = max(mint2d)
   clim=[min_color, max_color]
endif else begin
   min_color=clim(0)
   max_color=clim(1)
endelse
;
;
get_contour_colors2, clim, ncolors, nlevels_iso, bottom, $
                     levels_final, levels_iso, c_colors_final,  $
                     cb_indx, cb_tickmark_name

contourplots_map_double, mint2d_sc, mint2d_fvm, phi_grid, phi_grid, mu_grid, mu_grid, $
                         ncolors, nlevels_iso, bottom, $
                         levels_final, levels_iso, c_colors, $
                         cb_indx, cb_tickmark_name, $
                         ctable=13, $
                         center_lat=center_lat, center_lon=center_lon, $
                         plimit=plimit, titlestr1=titlestr1, titlestr2=titlestr2, $
                         ctitlestr=ctitlestr, /isotropic, /cb_top
;
end
;
;-----------------------------------------------------------------------
;-----------------------------------------------------------------------
;-----------------------------------------------------------------------
;
pro muc, xp, yp, zp, phi, muc1, muc2
;
;calculate critical mu (tangential ray at surface) at point (xp,yp,zp)
;for a given phi-angle
;
a=cos(phi)*xp+sin(phi)*yp
b=zp
c=-sqrt(xp^2+yp^2+zp^2-1.d0)
;
disc=a^2+b^2-c^2
if(disc ge 0.) then begin
   muc1 = (-b*c-a*sqrt(disc))/(a^2+b^2)
   muc2 = (-b*c+a*sqrt(disc))/(a^2+b^2)

   f1=a*sqrt(1.d0-muc1^2)+b*muc1+c
   f2=a*sqrt(1.d0-muc2^2)+b*muc2+c
   if(abs(f1) gt 1.d-6) then muc1=10.d0
   if(abs(f2) gt 1.d-6) then muc2=10.d0

endif else begin
;return dummy value if no tangential ray exists at all
   muc1 = 10.d0
   muc2 = 10.d0
endelse
;
mucmin=min([muc1,muc2])
mucmax=max([muc1,muc2])
muc1=mucmin
muc2=mucmax
;
end
