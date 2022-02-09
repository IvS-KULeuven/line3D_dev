pro surface_rot1, r_pole, m_star, v_rot, r_eq, v_crit, omega, theta, r_surface, x_surface, z_surface, l_star=l_star, teff_surface=teff_surface, beta_z=beta_z, help=print_help
;
;+
; NAME:
;       surface_rot1
;
; PURPOSE:
;       This procedure calculates the surface of a rotating star from
;       a roche model (e.g. Cranmer/Owocki 1995)
;
; CALLING SEQUENCE:
;       surface_rot1, r_pole, m_star, v_rot, r_eq, v_crit, omega, theta, r_surface, x_surface, z_surface
;
; INPUTS:
;       r_pole: polar radius in r_sun
;       m_star: stellar mass in m_sun (without eddington factor!!!)
;       v_rot:  rotational velocity at equator in km/s
;
; OUTPUTS:
;       r_eq:   equatorial radius in units of r_pole
;       v_crit: critical velocity in km/s
;       omega:  omega = omega_eq / omega_crit   (omega_eq, omega_crit the
;               angular velocities at equator and critical)
;       theta:       co-latitude (for which radius is calculated)
;       r_surface:   radius(theta) of stellar surface in units of r_pole
;       x_surface:   x-coordinate of stellar surface in units of r_pole
;       z_surface:   z-coordinate of stellar surface in units of r_pole
; 
; KEYWORD PARAMETERS:
;       help:   Set this keyword (flag) tho show the documentation of this
;               function
;       teff_surface:  Set this keyword to calculate the gravity darkening
;       l_star: stellar luminosity in l_sun
;       beta_z: von-Zeipel parameter (default: 1/4)
;
; EXAMPLE:
;       surace_rot1, 19., 52.5, 379.d0, r_eq, v_crit, omega, theta, r_surface, $
;                    x_surface, z_surface, teff_surface=teff_surface
;-
;
;-----------------------------------------------------------------------
;
if(keyword_set(print_help)) then begin
   doc_library, 'surface_rot1'
   return
endif
;
;-----------------------------------------------------------------------
;
r_pole_cgs = r_pole*!rsu
m_star_cgs = m_star*!msu
v_rot_cgs = v_rot*1.d5
;
w0 = v_rot_cgs^2*r_pole_cgs/2.d0/!cgs_grav/m_star_cgs
;
r_eq_cgs = r_pole_cgs/(1.d0-w0)
r_eq = r_eq_cgs/r_pole_cgs
;
v_crit_cgs = sqrt(2.d0*!cgs_grav*m_star_cgs/3.d0/r_pole_cgs)
v_crit = v_crit_cgs/1.d5
;
omega = v_rot_cgs*1.5d0*r_pole_cgs/r_eq_cgs/v_crit_cgs
;
;
;calculate surface (in units of r_pole)
nd=101
theta_min=1.d-4
theta_max=!pi-1.d-4
;
theta=theta_min+findgen(nd)*(theta_max-theta_min)/(nd-1)
;
r_surface=3.d0/omega/sin(theta)*cos((!pi+acos(omega*sin(theta)))/3.d0)
;
x_surface=r_surface*sin(theta)
z_surface=r_surface*cos(theta)
;
;
;----------------------calculate gravity darkening law------------------
;-------------following petrenz/puls 1996, cranmer/owocki 1995----------
;---------------------(neglecting eddington factor!!!)------------------
;
nd_fine=101
theta_min=1.d-4
theta_max=!pi/2.d0-1.d-4
;
if(not keyword_set(beta_z)) then beta_z = 1.d0/4.d0
;
theta_fine=theta_min+findgen(nd_fine)*(theta_max-theta_min)/(nd_fine-1)
;
if(not keyword_set(l_star)) then l_star=1.d5
l_star_cgs = l_star*!lsu
;
w0 = v_rot_cgs^2*r_pole_cgs/2.d0/!cgs_grav/m_star_cgs
omega = sqrt(27.d0/4.d0*w0*(1.d0-w0)^2)
;
;
gr_fine=fltarr(nd_fine)*0.d0
gtheta_fine=fltarr(nd_fine)*0.d0
gperp_fine=fltarr(nd_fine)*0.d0
integrand_fine=fltarr(nd_fine)*0.d0
rsurface_fine=fltarr(nd_fine)*0.d0
;
for i=0, nd_fine-1 do begin
;
;radius as function of theta in units of rstar
   rsurface_fine(i)=3.d0/omega/sin(theta_fine(i))*cos((!pi+acos(omega*sin(theta_fine(i))))/3.d0)
;gravity components
   gr_fine(i) = !cgs_grav*m_star_cgs/r_pole_cgs^2 * (-1.d0/rsurface_fine(i)^2 + 8.d0/27.d0*rsurface_fine(i)*omega^2*sin(theta_fine(i))^2)
   gtheta_fine(i) = !cgs_grav*m_star_cgs/r_pole_cgs^2 * 8.d0/27.d0*rsurface_fine(i)*omega^2*sin(theta_fine(i))*cos(theta_fine(i))
   gperp_fine(i) = sqrt(gr_fine(i)^2+gtheta_fine(i)^2)
;
;integrand to determine sigma
   integrand_fine(i) = 4.d0*!pi*(gperp_fine(i))^(4.d0*beta_z) * (rsurface_fine(i)*r_pole_cgs)^2*sin(theta_fine(i)) / (-gr_fine(i)/gperp_fine(i))
;   print, theta_fine(i), rsurface_fine(i), gr_fine(i), gtheta_fine(i), gperp_fine(i), integrand_fine(i)
endfor

;
;calculate sigma
sigma=0.d0
for i=1, nd-1 do begin
   sigma = sigma + (integrand_fine(i)+integrand_fine(i-1))*(theta_fine(i)-theta_fine(i-1))*0.5d0
endfor
sigma_fit = 4.d0*!pi*!cgs_grav*m_star_cgs*(1.d0-0.1969*omega^2 - 0.094292*omega^4 + $
            0.33812*omega^6 - 1.3066*omega^8 + 1.8286*omega^10 - 0.92714*omega^12)
;
;calculate effective temperature
teff_surface_fine = (l_star_cgs*gperp_fine^(4.d0*beta_z)/!cgs_sb/sigma)^0.25d0
;
;print, 'sigma', sigma, sigma_fit
;print, 'r_pole', r_pole_cgs/!rsu
;print, 'm_star', m_star_cgs/!msu
;print, 'v_rot', v_rot_cgs/1.d5
;print, 'l_star', l_star_cgs/!lsu
;print, 'w0', w0
;print, 'omega', omega
;for i=0, nd_fine-1 do begin
;   print, theta_fine(i)*180.d0/!pi, rsurface_fine(i), alog10(gperp_fine(i)), integrand_fine(i), teff_surface_fine(i)
;endfor
;stop
;
;
teff_surface=fltarr(nd)
;
for i=0, nd-1 do begin
   find_indx, r_surface(i), rsurface_fine, nd_fine, iim1, ii
   teff_surface(i) = interpol_ypfct(rsurface_fine(iim1), rsurface_fine(ii), teff_surface_fine(iim1), teff_surface_fine(ii), r_surface(i))
endfor
;
end
