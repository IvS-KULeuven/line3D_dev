pro modelspec1d_vbin, fname=fname, cs1_thetac=cs1_thetac, cs1_phic=cs1_phic, cs2_thetac=cs2_thetac, cs2_phic=cs2_phic, windx=windx, oname=oname, $
                      xlim=xlim, print_help=print_help
;
;+
; NAME:
;	modelspec1d_vbin
;
; PURPOSE:
;	This procedure plots the radial stratification through a 3d
;	grid specified by angles cs1_thetac, cs1_phic (for coordinate
;	system 1) and cs1_thetac, cs1_phic (for coordinate system 2) and 
;
; CALLING SEQUENCE:
;	model1dspec_vbin, fname=fname
;
; INPUTS:
;	
; KEYWORD PARAMETERS:
;       cs1_thetac:  latitude in system1
;       cs2_thetac:  latitude in system2
;       cs1_phic:  azimuth in system1
;       cs2_phic:  azimuth in system2  
;       windx:  Set this keyword to an integer, defining the window-index
;       oname:  Set this keyword to a string, defining the output-name
;               (ps-file)
;       xlim:   Set this keyword to an array of length 2, defining the xrange
;	help:	Set this keyword (flag) tho show the documentation of this procedure
;
; OUTPUTS:
;	
; EXAMPLE:
;	modelspec1d_vbin, fname='modspec_model00.h5'
;       modelspec1d_vbin, fname='modspec_model00.h5', cs1_thetac=90., cs1_phic=90.
;-
;
;-------------------------OUTPUT IF HELP NEEDED-------------------------
;
pdefault=!p
;
if(keyword_set(print_help)) then begin
   doc_library, 'modelspec1d_vbin'
   return
endif

fname='/home/levin/Postdoc/line3D/outputFILES/modspec_model00.h5'
;
if(n_elements(fname) eq 0) then begin
   print, 'error in modelspec1d_vbin: fname not specified'
   doc_library, 'modelspec1d_vbin'
   stop
endif
;
if(not keyword_set(windx)) then windx=0
;
if(not keyword_set(cs1_thetac)) then cs1_thetac=0.d0
if(not keyword_set(cs1_phic)) then cs1_phic=0.d0
if(not keyword_set(cs2_thetac)) then cs2_thetac=cs1_thetac
if(not keyword_set(cs2_phic)) then cs2_phic=cs1_phic

cs1_thetac=cs1_thetac*!pi/180
cs1_phic=cs1_phic*!pi/180.
cs2_thetac=cs2_thetac*!pi/180
cs2_phic=cs2_phic*!pi/180.
;
;-----------------------------------------------------------------------
;
read_modelspec3d_vbin, fname, $
         x01=x01, y01=y01, z01=z01, x02=x02, y02=y02, z02=z02, $
         vx01=vx01, vy01=vy01, vz01=vz01, vx02=vx02, vy02=vy02, vz02=vz02, $
         rstar1=unit_length1, teff1=teff1, logg1=logg1, lstar1=lstar1, yhe1=yhe1, vrot1=vrot1, $
         rstar2=unit_length2, teff2=teff2, logg2=logg2, lstar2=lstar2, yhe2=yhe2, vrot2=vrot2, $
         vth_fiducial=vth_fiducial, vmicro=vmicro, vmax=vmax, xnue0=xnue0, na=na, $
         unit_length=unit_length0, $
         cs1_nr=cs1_nr, cs1_ntheta=cs1_ntheta, cs1_nphi=cs1_nphi, $
         cs2_nr=cs2_nr, cs2_ntheta=cs2_ntheta, cs2_nphi=cs2_nphi, $                           
         cs1_radius=cs1_radius, cs1_theta=cs1_theta, cs1_phi=cs1_phi, $
         cs2_radius=cs2_radius, cs2_theta=cs2_theta, cs2_phi=cs2_phi, $
         cs1_scont3d=cs1_scont3d, cs1_sline3d=cs1_sline3d, cs1_t3d=cs1_t3d, cs1_rho3d=cs1_rho3d, cs1_opac3d=cs1_opac3d, $
         cs1_opalbar3d=cs1_opalbar3d, cs1_velx3d=cs1_velx3d, cs1_vely3d=cs1_vely3d, cs1_velz3d=cs1_velz3d, $
         cs2_scont3d=cs2_scont3d, cs2_sline3d=cs2_sline3d, cs2_t3d=cs2_t3d, cs2_rho3d=cs2_rho3d, cs2_opac3d=cs2_opac3d, $
         cs2_opalbar3d=cs2_opalbar3d, cs2_velx3d=cs2_velx3d, cs2_vely3d=cs2_vely3d, cs2_velz3d=cs2_velz3d            

;opacity units in 1/rstar
sr1 = unit_length1*!rsu
cs1_opac3d=cs1_opac3d*sr1
cs1_opalbar3d=cs1_opalbar3d*sr1
;
sr2 = unit_length2*!rsu
cs2_opac3d=cs2_opac3d*sr2
cs2_opalbar3d=cs2_opalbar3d*sr2

;
;------------------------------density----------------------------------
;
if(keyword_set(oname)) then oname='modelspec1d_vbin_cs1.ps'
;
clim=[-14.d0,-3.d0]
logscale=1
;
titlestr1='star 1'
titlestr2='star 2'
ctitlestr=textoidl('log(<\chi>)')


if(keyword_set(oname)) then begin
   print, "writing output to: ", oname
   set_plot,'ps'
   device,file=oname, $
    decomposed=0, color=1, bits_per_pixel=8
endif else begin
   window, windx
   device, decomposed=0
   windx=windx+1
endelse
;
plot_modelspec3d_spc, cs1_radius, cs1_theta, cs1_phi, cs1_thetac, cs1_phic, $
                      cs1_rho3d, cs1_opac3d, cs1_opalbar3d, cs1_velx3d, cs1_vely3d, cs1_velz3d, $
                      cs1_t3d, cs1_scont3d, cs1_sline3d, xlim=[0.9,10.], windx=windx, oname=oname, vmax=1.d5
;
!p=pdefault

end
;
;-----------------------------------------------------------------------
;-----------------------------------------------------------------------
;-----------------------------------------------------------------------
;
pro plot_modelspec3d_spc, radius, theta, phi, thetac, phic, $
                          rho3d, opac3d, opalbar3d, velx3d, vely3d, velz3d, $
                          t3d, scont3d, sline3d, xlim=xlim, windx=windx, oname=oname, vmax=vmax
;
;----------------------interpolate onto radial grid---------------------
;
;
nr=n_elements(radius)
ntheta=n_elements(theta)
nphi=n_elements(phi)
;
opac1d=fltarr(nr)*0.d0
opalbar1d=fltarr(nr)*0.d0
rho1d=fltarr(nr)*0.d0
t1d=fltarr(nr)*0.d0
velx1d=fltarr(nr)*0.d0
vely1d=fltarr(nr)*0.d0
velz1d=fltarr(nr)*0.d0
velr1d=fltarr(nr)*0.d0
velphi1d=fltarr(nr)*0.d0
scont1d=fltarr(nr)*0.d0
sline1d=fltarr(nr)*0.d0

find_indx, thetac, theta, ntheta, ii, iim1
find_indx, phic, phi, nphi, jj, jjm1
;
sint = sin(thetac)
cost = cos(thetac)
sinp = sin(phic)
cosp = cos(phic)
;
;
for i=0, nr-1 do begin
   opac1d(i) = interpol2d_4p_lin(opac3d(i,iim1,jjm1), opac3d(i,ii,jjm1), opac3d(i,iim1,jj), opac3d(i,ii,jj), $
                                 theta(iim1),theta(ii),phi(jjm1),phi(jj),thetac,phic)
   opalbar1d(i) = interpol2d_4p_lin(opalbar3d(i,iim1,jjm1), opalbar3d(i,ii,jjm1), opalbar3d(i,iim1,jj), opalbar3d(i,ii,jj), $
                                    theta(iim1),theta(ii),phi(jjm1),phi(jj),thetac,phic)
   rho1d(i) = interpol2d_4p_lin(rho3d(i,iim1,jjm1), rho3d(i,ii,jjm1), rho3d(i,iim1,jj), rho3d(i,ii,jj), $
                                  theta(iim1),theta(ii),phi(jjm1),phi(jj),thetac,phic)   
   t1d(i) = interpol2d_4p_lin(t3d(i,iim1,jjm1), t3d(i,ii,jjm1), t3d(i,iim1,jj), t3d(i,ii,jj), $
                                 theta(iim1),theta(ii),phi(jjm1),phi(jj),thetac,phic)
   velx1d(i) = interpol2d_4p_lin(velx3d(i,iim1,jjm1), velx3d(i,ii,jjm1), velx3d(i,iim1,jj), velx3d(i,ii,jj), $
                                 theta(iim1),theta(ii),phi(jjm1),phi(jj),thetac,phic)
   vely1d(i) = interpol2d_4p_lin(vely3d(i,iim1,jjm1), vely3d(i,ii,jjm1), vely3d(i,iim1,jj), vely3d(i,ii,jj), $
                                 theta(iim1),theta(ii),phi(jjm1),phi(jj),thetac,phic)   
   velz1d(i) = interpol2d_4p_lin(velz3d(i,iim1,jjm1), velz3d(i,ii,jjm1), velz3d(i,iim1,jj), velz3d(i,ii,jj), $
                                 theta(iim1),theta(ii),phi(jjm1),phi(jj),thetac,phic)
   scont1d(i) = interpol2d_4p_lin(scont3d(i,iim1,jjm1), scont3d(i,ii,jjm1), scont3d(i,iim1,jj), scont3d(i,ii,jj), $
                                 theta(iim1),theta(ii),phi(jjm1),phi(jj),thetac,phic)
   sline1d(i) = interpol2d_4p_lin(sline3d(i,iim1,jjm1), sline3d(i,ii,jjm1), sline3d(i,iim1,jj), sline3d(i,ii,jj), $
                                 theta(iim1),theta(ii),phi(jjm1),phi(jj),thetac,phic)   

   velr1d(i) = velx1d(i)*sint*cosp + vely1d(i)*sint*sinp + velz1d(i)*cost
   velphi1d(i) = -velx1d(i)*sinp + vely1d(i)*cosp
   
endfor
;
;-----------------------------------------------------------------------
;
if(not keyword_set(windx)) then windx=0

if(keyword_set(oname)) then begin
   print, "writing output to: ", oname
   set_plot,'ps'
   device,file=oname, $
    decomposed=0, color=1, bits_per_pixel=8
endif else begin
   window, windx, ysize=1200
   device, decomposed=0
   windx=windx+1
endelse

loadct_lev, ci_black, ci_blue, ci_red, ci_cyan, ci_green, ci_magenta, ci_yellow, ci_white, ci_dgreen=ci_dgreen
;
!p.multi=[6,2,3]
plot, radius, opac1d, /ylog, $
      xrange=[0.9,10.], /xs, $
      xtitle=textoidl('r'), $
      ytitle=textoidl('\chi_c [1/R_\ast]'), $
      charsize=4.

!p.multi=[5,2,3]
plot, radius, opalbar1d, /ylog, $
      xrange=[0.9,10.], /xs, $
      xtitle=textoidl('r'), $
      ytitle=textoidl('<\chi_L> [1/s/R_\ast]'), $
      charsize=4.

!p.multi=[4,2,3]
plot, radius, rho1d, /ylog, $
      xrange=[0.9,10.], /xs, $
      xtitle=textoidl('r'), $
      ytitle=textoidl('\rho [g/cm^3]'), $
      charsize=4.

!p.multi=[3,2,3]
if(not keyword_set(vmax)) then vmax=max([velx1d,vely1d,velz1d,velr1d])
plot, [0.,0.], [0.,0.], $
      xrange=[0.9,10.], /xs, $
      yrange=[min([velr1d/vmax,velphi1d/vmax,velx1d/vmax,vely1d/vmax,velz1d/vmax]), $
              max([velr1d/vmax,velphi1d/vmax,velx1d/vmax,vely1d/vmax,velz1d/vmax])], $
      /ys, $
      xtitle=textoidl('r'), $
      ytitle=textoidl('v_i'), $
      charsize=4.
oplot, radius, velx1d/vmax, color=ci_blue
oplot, radius, vely1d/vmax, color=ci_red
oplot, radius, velz1d/vmax, color=ci_green
oplot, radius, velphi1d/vmax, colo=ci_dgreen
oplot, radius, velr1d/vmax
legend, [textoidl('v_r'),$
         textoidl('v_\phi'), $
         textoidl('v_x'), $
         textoidl('v_y'), $
         textoidl('v_z')], $
         line=[0,0,0,0,0], $
         color=[ci_white,ci_dgreen,ci_blue,ci_red,ci_green], $
        charsize=1.5, $
;        /bottom, $
        /right_legend

!p.multi=[2,2,3]
plot, radius, scont1d, $
      xrange=[0.9,10.], /xs, $
      xtitle=textoidl('r'), $
      ytitle=textoidl('S_c'), $
      charsize=4.

!p.multi=[1,2,3]
plot, radius, sline1d, /ylog, $
      xrange=[0.9,10.], /xs, $
      xtitle=textoidl('r'), $
      ytitle=textoidl('S_L'), $
      charsize=4.


if keyword_set(oname) then begin
   device, /close
   set_plot,'x'
endif



end
