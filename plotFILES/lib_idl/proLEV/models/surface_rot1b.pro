
;
;-----------------------------------------------------------------------
;-----------------------------------------------------------------------
;-----------------------------------------------------------------------
;
function g_fct, x, const
;
;
g_fct = cos(x) + alog(tan(x/2.d0)) - const
return, g_fct
;
end
;
;-----------------------------------------------------------------------
;-----------------------------------------------------------------------
;-----------------------------------------------------------------------
;
pro regula_falsi, xl, xu, x_root, const
;
;-----------------------------------------------------------------------
;
;FINDS THE ROOT OF A FUNCTION G_FCT WITH THE REGULA-FALSI METHOD
;
;ON INPUT: xl, xu: FIRST GUESS IN WHICH INTERVAL THE ROOT MAY BE
;                  BOTH ARE SET IN CORRECT ORDER, SUCH THAT xl<xu
;          FUNCTION FOR WHICH ROOT SHALL BE FOUND (SPECIFIED BELOW AT 'local functions')
;
;ON OUTPUT: x_root: ROOT OF A FUNCTION f(x)=0
;
nit=10000000L
eps=1.d-14
;
;----first: enlarge interval if y_lower and y_upper have equal sign-----
;---second: if they still have equal sign, then reduce the interval-----
;
;swap lower and upper values, if necessary
if(xl gt xu) then begin
   swap=xl
   xl=xu
   xu=swap
endif
;
;enlarge interval if necessary
x_lower=xl
x_upper=xu
;
y_lower = g_fct(x_lower, const)
y_upper = g_fct(x_upper, const)
;
if(y_upper ne y_upper) then begin
   for i=0, 100-1 do begin
      x_mid = (x_lower+x_upper)/2.d0
      y_mid = g_fct(x_mid,const)
;      print, 'mid', x_mid, y_mid, x_upper, x_lower
      if(y_mid ne y_mid) then begin
         x_upper = x_mid
         y_upper = y_mid
      endif else begin
         if(y_lower*y_mid gt 0.d0) then begin
            x_lower = x_mid
            y_lower = y_mid
         endif else begin
            x_upper = x_mid
            y_upper = y_mid
            break
         endelse
      endelse
   endfor
endif
;
;
enlint, x_lower, x_upper, const
;
;
;
y_lower = g_fct(x_lower, const)
y_upper = g_fct(x_upper, const)
;
;reduce right interval-boundary (from initial start values) if necessary
if(y_lower*y_upper gt 0.d0) then begin
   x_lower=xl
   x_upper=xu
   redint, x_lower, x_upper, const
   y_lower=g_fct(x_lower, const)
   y_upper=g_fct(x_upper, const)
endif
;
if(y_lower*y_upper gt 0.d0) then begin
;error if no interval is found
   print,  'error in regula_falsi: no interval for which f(xl) < f(xu) is found'
   stop
endif
;
;--now, when start values give unequal signs, make regula falsi method--
;
for i=1, nit do begin
;
   x_new = x_lower - y_lower * (x_upper-x_lower)/(y_upper-y_lower)
   y_new = g_fct(x_new, const)

;   print, i
;   print, x_lower, x_new, x_upper
;   print, y_lower, y_new, y_upper
;   print, ''
;
   if(abs(y_new) le eps) then break
;
   if(y_new*y_lower lt 0.d0 ) then begin
      x_upper=x_new
      y_upper=y_new
   endif else begin
      x_lower=x_new
      y_lower=y_new
   endelse
;
endfor
;
x_root=x_new
;
;
;
end
;
;-----------------------------------------------------------------------
;-----------------------------------------------------------------------
;-----------------------------------------------------------------------
;
pro enlint, x_lower, x_upper, const
;
;THIS ROUTINE ENLARGES INTERVAL [x_lower, x_upper] UNTIL f(x_lower)*f(x_upper)<0.
;INTERVAL IS ENLARGED WITH INCREASING STEPS
;
;(SEE ALSO NUMERICAL RECIPES IN F77, p. 345)
;
nit=100
fac=1.5d0
;
y_lower=g_fct(x_lower, const)
y_upper=g_fct(x_upper, const)
;
for i=1, nit do begin
;
   if(y_lower*y_upper gt 0.d0) then begin
;calculate new range
      if(abs(y_lower) le abs(y_upper)) then begin
;lower the left boundary
         x_lower=x_lower-fac*(x_upper-x_lower)
         y_lower=g_fct(x_lower, const)
      endif else begin
;enlarge right boundary
         x_upper=x_upper+fac*(x_upper-x_lower)
         y_upper=g_fct(x_upper, const)
      endelse
   endif else begin
      break
   endelse
;
endfor
;
end
;
;-----------------------------------------------------------------------
;-----------------------------------------------------------------------
;-----------------------------------------------------------------------
;
pro redint, x_lower, x_upper, const
;
;THIS ROUTINE REDUCES INTERVAL [x_lower, x_upper] UNTIL f(x_lower)*f(x_upper)<0.
;INTERVAL IS REDUCED WITH CONSTANT STEPS
;NOTE: ONLY x_upper is reduced down to x_lower
;
nit=100
;
;note: division by nit, not by (nit-1), such that x_upper ne x_lower at nit=100
del=(x_upper-x_lower)/nit
;
y_lower=g_fct(x_lower, const)
y_upper=g_fct(x_upper, const)
;
for i=1, nit do begin
;
   if(y_lower*y_upper gt 0.d0) then begin
;reduce right boundary
      x_upper = x_upper - del
      y_upper = g_fct(x_upper, const)
   endif else begin 
      break
   endelse
endfor
;
;  
;
end
;
;-----------------------------------------------------------------------
;-----------------------------------------------------------------------
;-----------------------------------------------------------------------
;
pro surface_rot1b, r_pole, m_star, v_rot, r_eq, v_crit, omega, theta, r_surface, x_surface, z_surface, l_star=l_star, teff_surface=teff_surface, help=print_help
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
;                      according to omega-model (Espinosa 2011)
;       l_star: stellar luminosity in l_sun
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
;-----------------following omega-model from espinosa 2011--------------
;
omega_1 = v_rot_cgs/r_eq_cgs
omega_k = sqrt(!cgs_grav*m_star_cgs/r_eq_cgs^3)
omega_f = omega_1/omega_k

;print, omega_f, omega
;stop
;
nd_fine=101
theta_min=1.d-4
theta_max=!pi/2.d0-1.d-4
;
theta_fine=theta_min+findgen(nd_fine)*(theta_max-theta_min)/(nd_fine-1)
;
if(not keyword_set(l_star)) then l_star=1.d5
l_star_cgs = l_star*!lsu
;
w0 = v_rot_cgs^2*r_pole_cgs/2.d0/!cgs_grav/m_star_cgs
omega = sqrt(27.d0/4.d0*w0*(1.d0-w0)^2)
;
gr_fine=fltarr(nd_fine)*0.d0
gtheta_fine=fltarr(nd_fine)*0.d0
gperp_fine=fltarr(nd_fine)*0.d0
integrand_fine=fltarr(nd_fine)*0.d0
rsurface_fine=fltarr(nd_fine)*0.d0
teff_surface_fine=fltarr(nd_fine)*0.d0
;
;
for i=0, nd_fine-1 do begin
;
;radius as function of theta in units of r_pole
   rsurface_fine(i)=3.d0/omega/sin(theta_fine(i))*cos((!pi+acos(omega*sin(theta_fine(i))))/3.d0)
;radius as function of theta in units of r_eq
   rsurface_eq = rsurface_fine(i)/r_eq
;
;calculate function psi
   const = (omega_f^2*rsurface_eq^3*cos(theta_fine(i))^3)/3.d0 + cos(theta_fine(i)) + alog(tan(theta_fine(i)/2.d0))
   psi_min=1.d-10
   psi_max=10.d0
   regula_falsi, psi_min, psi_max, psi, const

;effective temperature
   fdum1 = (l_star_cgs/4.d0/!pi/!cgs_sb/r_eq_cgs^2)^0.25d0
   fdum2 = (1.d0/rsurface_eq^4 + omega_f^4*rsurface_eq^2*sin(theta_fine(i))^2 - 2.d0*omega_f^2*sin(theta_fine(i))^2/rsurface_eq)^(1.d0/8.d0)
   fdum3 = sqrt(tan(psi)/tan(theta_fine(i)))
   teff_surface_fine(i) = fdum1*fdum2*fdum3
   print, fdum1, fdum2, fdum3, psi, theta_fine(i), tan(psi), tan(theta_fine(i)), cos(psi)+alog(tan(psi/2.d0))-const
endfor

;at the equator: analytic solution
teff_surface_fine(nd_fine-1) = sqrt(2.d0/(2.d0+omega_f^2))*(1.d0-omega_f^2)^(1.d0/12.d0) * exp(-4/3.d0*omega_f^2/(2.d0+omega_f^2)^3) * teff_surface_fine(0)
;print, teff_surface_fine
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
