pro surface_rot2, r_pole, m_star, v_rot, r_eq, v_crit, omega, theta, r_surface, x_surface, z_surface, help=print_help
;
;+
; NAME:
;       surface_rot2
;
; PURPOSE:
;       This procedure calculates the surface of a rotating star from
;       an ellipsoid
;
; CALLING SEQUENCE:
;       surface_rot2, r_pole, m_star, v_rot, r_eq, v_crit, omega, theta, r_surface, x_surface, z_surface
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
;       help:   Set this keyword (flag) tho show the documentation of this function
;
; EXAMPLE:
;       surace_rot2, 19., 52.5, 379.d0, r_eq, v_crit, omega, theta, r_surface, $
;                    x_surface, z_surface
;-
;
;-----------------------------------------------------------------------
;
if(keyword_set(print_help)) then begin
   doc_library, 'surface_rot2'
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
a=r_eq
b=1.d0
r_surface = a*b/sqrt(b^2*sin(theta)*sin(theta)+a^2*cos(theta)*cos(theta))
;
x_surface=r_surface*sin(theta)
z_surface=r_surface*cos(theta)
;
;
;
end
