import numpy as np
from const_cgs import *
#
def grid_equi(xmin,xmax,nx):
#
#--------------calculates equidistant grid between xmin,xmax------------
#
   xarr = xmin + np.arange(0,nx,1)*(xmax-xmin)/(nx-1)
   return xarr

########################################################################

def grid_log(xmin,xmax,nx):

   if xmax/xmin <= 0.:
      print('error in grid_log: xmax/xmin le 0.')
      exit()

   delta = np.log10(xmax/xmin)/np.float(nx-1)
   xarr = np.zeros(nx)
   xarr[0] = xmin
   for i in np.arange(1,nx):
      xarr[i] = xarr[i-1]*10.**delta

   return xarr
#
########################################################################  
#
def get_angles_spc(x, y, z):
#
#--------------calculates angles in spherical coordinates---------------
#---------------for given carthesian coordinates x, y, z----------------
#----------------------for right-handed system--------------------------
#
#on input: coordinates x, y, z
#
#on output: angles theta, phi
#
#(from jon's formal solver)
#
#-----------------------calculation of theta----------------------------
#
   rad = np.sqrt(x**2+y**2+z**2)
#
   if rad != 0.:
      theta = np.arccos(z/rad)
   else:
      theta = np.double(0.)
#
#------------------------calculation of phi-----------------------------
#
   if z == rad:
      phi=np.double(0.)
      return theta,phi
#
   if x == 0.:
      if y == 0.:
         phi = np.double(0.)
      else:
         if y > 0.:
            phi = np.pi/2.
         else:
            phi = 3.*np.pi/2.
   else:
      if y > 0. and x > 0.:
#first quadrant
         phi = np.arctan(y/x)
      else:
         if x > 0.:
#fourth quadrant
            phi = np.double(2.)*np.pi + np.arctan(y/x)
         else:
#second and third quadrant
            phi = np.pi + np.arctan(y/x)

#   theta=round(theta,6)
#   phi=round(phi,6)
#
#-------------------------test if everything worked---------------------
#
   x_test=rad*np.sin(theta)*np.cos(phi)
   y_test=rad*np.sin(theta)*np.sin(phi)
   z_test=rad*np.cos(theta)
#
   tol=1.e-6
   if np.abs(x_test-x) > tol:
      print('error in get_angles_spc', x_test, x)
      exit()
   if np.abs(y_test-y) > tol:
      print('error in get_angles_spc', y_test, y)
      exit()
   if np.abs(z_test-z) > tol:
      print('error in get_angles_spc', z_test, z)
      exit()
#
   return theta, phi
#
#-----------------------------------------------------------------------
#--------------transform vector from cartesian components to------------
#-------------------------spherical components--------------------------
#-----------------------------------------------------------------------
#
def transvec_cac_spc(x,y,z,vec_x,vec_y,vec_z):
#
#input: coordinates in cartesian system, x,y,z
#       vector in cartesian system vec_x, vec_y, vec_z
   
   rad = x**2 + y**2 + z**2
   [theta,phi] = get_angles_spc(x, y, z)
   sint = np.sin(theta)
   cost = np.cos(theta)
   sinp = np.sin(phi)
   cosp = np.cos(phi)

   vec_r = vec_x*sint*cosp + vec_y*sint*sinp + vec_z*cost
   vec_theta = vec_x*cost*cosp + vec_y*cost*sinp - vec_z*sint
   vec_phi = -vec_x*sinp + vec_y*cosp

   return vec_r, vec_theta, vec_phi
#
#-----------------------------------------------------------------------
#--------------transform 3x3 tensfor from cartesian components to-------
#-------------------------spherical components--------------------------
#-----------------------------------------------------------------------
#
def transten_cac_spc(x,y,z,t_xx,t_yy,t_zz, t_xy,t_xz, t_yz):
#
#input: coordinates in cartesian system, x,y,z
#       (symmetric) tensort in cartesian system t_xx, t_xy, t_xz
#                                               t_xy, t_yy, t_yz
#                                               t_xz, t_yz, t_zz
#output: tensor components in spherical system
   
   rad = x**2 + y**2 + z**2
   [theta,phi] = get_angles_spc(x, y, z)
   sint = np.sin(theta)
   cost = np.cos(theta)
   sinp = np.sin(phi)
   cosp = np.cos(phi)

   a11 = t_xx*sint*cosp + t_xy*sint*sinp + t_xz*cost
   a12 = t_xx*cost*cosp + t_xy*cost*sinp - t_xz*sint
   a13 = -t_xx*sinp + t_xy*cosp
   a21 = t_xy*sint*cosp + t_yy*sint*sinp + t_yz*cost
   a22 = t_xy*cost*cosp + t_yy*cost*sinp - t_yz*sint
   a23 = -t_xy*sinp + t_yy*cosp
   a31 = t_xz*sint*cosp + t_yz*sint*sinp + t_zz*cost
   a32 = t_xz*cost*cosp + t_yz*cost*sinp - t_zz*sint
   a33 = -t_xz*sinp + t_yz*cosp
   t_rr = a11*sint*cosp + a21*sint*sinp + a31*cost
   t_rth = a12*sint*cosp + a22*sint*sinp + a32*cost
   t_rphi = a13*sint*cosp + a23*sint*sinp + a33*cost
   t_thth = a12*cost*cosp + a22*cost*sinp - a32*sint
   t_thphi = a13*cost*cosp + a23*cost*sinp - a33*sint
   t_phiphi = -a13*sinp + a23*cosp

   return t_rr, t_thth, t_phiphi, t_rth, t_rphi, t_thphi
#
#-----------------------------------------------------------------------
#--------------read columns from a file---------------------------------
#-----------------------------------------------------------------------
#
def readcol(fname,ncols=0,nskip=0):

   if ncols == 0:
      print('error in readcol: ncols = 0 not allowed')
      exit()   
   if ncols > 23:
      print('error in readcol: ncols > 23 not implemented yet')
      exit()
   
#dirty, but working
   f = open(fname, 'r')

   for i in range(0,nskip):
      print('iskip', i)
      header = f.readline()

   nd = 0
   for line in f.readlines():
      nd=nd+1
   f.close()

#default: only working for float columns
   f = open(fname, 'r')

   for i in range(0,nskip):
      header = f.readline()


   data2d = np.zeros(shape=(nd,ncols))
   i=0
   for line in f.readlines():
      line = line.strip()
      columns = line.split()
      for jj in np.arange(0,ncols):
#         print(jj, columns[jj])
         data2d[i][jj] = columns[jj]
      i=i+1

   f.close()
#
#
#
   if ncols == 1:
      return data2d[:,0]

   if ncols == 2:
      return data2d[:,0], data2d[:,1]

   if ncols == 3:
      return data2d[:,0], data2d[:,1], data2d[:,2]

   if ncols == 4:
      return data2d[:,0], data2d[:,1], data2d[:,2], data2d[:,3]

   if ncols == 5:
      return data2d[:,0], data2d[:,1], data2d[:,2], data2d[:,3], data2d[:,4]

   if ncols == 6:
      return data2d[:,0], data2d[:,1], data2d[:,2], data2d[:,3], data2d[:,4], data2d[:,5]
   
   if ncols == 7:
      return data2d[:,0], data2d[:,1], data2d[:,2], data2d[:,3], data2d[:,4], data2d[:,5], \
             data2d[:,6]
   
   if ncols == 8:
      return data2d[:,0], data2d[:,1], data2d[:,2], data2d[:,3], data2d[:,4], data2d[:,5], \
             data2d[:,6], data2d[:,7]
   
   if ncols == 9:
      return data2d[:,0], data2d[:,1], data2d[:,2], data2d[:,3], data2d[:,4], data2d[:,5], \
             data2d[:,6], data2d[:,7], data2d[:,8]
   
   if ncols == 10:
      return data2d[:,0], data2d[:,1], data2d[:,2], data2d[:,3], data2d[:,4], data2d[:,5], \
             data2d[:,6], data2d[:,7], data2d[:,8], data2d[:,9]
   
   if ncols == 11:
      return data2d[:,0], data2d[:,1], data2d[:,2], data2d[:,3], data2d[:,4], data2d[:,5], \
             data2d[:,6], data2d[:,7], data2d[:,8], data2d[:,9], data2d[:,10]
   
   if ncols == 12:
      return data2d[:,0], data2d[:,1], data2d[:,2], data2d[:,3], data2d[:,4], data2d[:,5], \
             data2d[:,6], data2d[:,7], data2d[:,8], data2d[:,9], data2d[:,10], data2d[:,11]
   
   if ncols == 13:
      return data2d[:,0], data2d[:,1], data2d[:,2], data2d[:,3], data2d[:,4], data2d[:,5], \
             data2d[:,6], data2d[:,7], data2d[:,8], data2d[:,9], data2d[:,10], data2d[:,11], \
             data2d[:,12]
   
   if ncols == 14:
      return data2d[:,0], data2d[:,1], data2d[:,2], data2d[:,3], data2d[:,4], data2d[:,5], \
             data2d[:,6], data2d[:,7], data2d[:,8], data2d[:,9], data2d[:,10], data2d[:,11], \
             data2d[:,12], data2d[:,13]
   
   if ncols == 15:
      return data2d[:,0], data2d[:,1], data2d[:,2], data2d[:,3], data2d[:,4], data2d[:,5], \
             data2d[:,6], data2d[:,7], data2d[:,8], data2d[:,9], data2d[:,10], data2d[:,11], \
             data2d[:,12], data2d[:,13], data2d[:,14]
   
   if ncols == 16:
      return data2d[:,0], data2d[:,1], data2d[:,2], data2d[:,3], data2d[:,4], data2d[:,5], \
             data2d[:,6], data2d[:,7], data2d[:,8], data2d[:,9], data2d[:,10], data2d[:,11], \
             data2d[:,12], data2d[:,13], data2d[:,14], data2d[:,15]

   if ncols == 17:
      return data2d[:,0], data2d[:,1], data2d[:,2], data2d[:,3], data2d[:,4], data2d[:,5], \
             data2d[:,6], data2d[:,7], data2d[:,8], data2d[:,9], data2d[:,10], data2d[:,11], \
             data2d[:,12], data2d[:,13], data2d[:,14], data2d[:,15], data2d[:,16]

   if ncols == 18:
      return data2d[:,0], data2d[:,1], data2d[:,2], data2d[:,3], data2d[:,4], data2d[:,5], \
             data2d[:,6], data2d[:,7], data2d[:,8], data2d[:,9], data2d[:,10], data2d[:,11], \
             data2d[:,12], data2d[:,13], data2d[:,14], data2d[:,15], data2d[:,16], data2d[:,17]

   if ncols == 19:
      return data2d[:,0], data2d[:,1], data2d[:,2], data2d[:,3], data2d[:,4], data2d[:,5], \
             data2d[:,6], data2d[:,7], data2d[:,8], data2d[:,9], data2d[:,10], data2d[:,11], \
             data2d[:,12], data2d[:,13], data2d[:,14], data2d[:,15], data2d[:,16], data2d[:,17], \
             data2d[:,18]

   if ncols == 20:
      return data2d[:,0], data2d[:,1], data2d[:,2], data2d[:,3], data2d[:,4], data2d[:,5], \
             data2d[:,6], data2d[:,7], data2d[:,8], data2d[:,9], data2d[:,10], data2d[:,11], \
             data2d[:,12], data2d[:,13], data2d[:,14], data2d[:,15], data2d[:,16], data2d[:,17], \
             data2d[:,18], data2d[:,19]

   if ncols == 21:
      return data2d[:,0], data2d[:,1], data2d[:,2], data2d[:,3], data2d[:,4], data2d[:,5], \
             data2d[:,6], data2d[:,7], data2d[:,8], data2d[:,9], data2d[:,10], data2d[:,11], \
             data2d[:,12], data2d[:,13], data2d[:,14], data2d[:,15], data2d[:,16], data2d[:,17], \
             data2d[:,18], data2d[:,19], data2d[:,20]

   if ncols == 22:
      return data2d[:,0], data2d[:,1], data2d[:,2], data2d[:,3], data2d[:,4], data2d[:,5], \
             data2d[:,6], data2d[:,7], data2d[:,8], data2d[:,9], data2d[:,10], data2d[:,11], \
             data2d[:,12], data2d[:,13], data2d[:,14], data2d[:,15], data2d[:,16], data2d[:,17], \
             data2d[:,18], data2d[:,19], data2d[:,20], data2d[:,21]

   if ncols == 23:
      return data2d[:,0], data2d[:,1], data2d[:,2], data2d[:,3], data2d[:,4], data2d[:,5], \
             data2d[:,6], data2d[:,7], data2d[:,8], data2d[:,9], data2d[:,10], data2d[:,11], \
             data2d[:,12], data2d[:,13], data2d[:,14], data2d[:,15], data2d[:,16], data2d[:,17], \
             data2d[:,18], data2d[:,19], data2d[:,20], data2d[:,21], data2d[:,22]   
#
########################################################################
def smooth_gauss(x,y,sigma=1.,ngauss=101):
#
#smooth out a distribution y using a gaussian filter of width sigma
#over ngauss points
#note: works only for equidistant steps
#
   nd = np.size(x)
#
   ysmooth=np.copy(y)
#
   for i in range(0,nd):
      jstart = np.max(np.array([0, i-np.floor(ngauss/2.)]))
      jend = np.min(np.array([nd-1, i+np.floor(ngauss/2.)]))
      jstart=int(jstart)
      jend=int(jend)
#
      sums=0.
      norm=0.
      for j in range(jstart,jend):
         gauss = np.exp(-((x[i]-x[j])/sigma)**2)
         sums = sums + y[j]*gauss
         norm = norm + gauss

#      print(sums,norm)
      ysmooth[i] = sums/norm
#      print('t1', sum*(x(1)-x(0)), norm*(x(1)-x(0)))

   return ysmooth
#
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
#
def calc_transmat2(alpha,gamma):
#
#+
# NAME:
#	CALC_TRANSMAT2
#
# PURPOSE:
#	This procedure calculates a transformation matrix to transform
#       a cylindrical coordinate system (e_x, e_y, e_z) to 
#       a carthesian coordinate system (ee_x, ee_y, ee_z).
#       The cylindrical coordinate system is specfied by viewing angles.
#       
#
# CALLING SEQUENCE:
#	CALC_TRANSMAT2, alpha, gamma, transmat
#
# INPUTS:
#	alpha:	viewing angle with respect to the z-axis of carthesian
#     	        coordinates in radiant
#	gamma:	viewing angle with respect to the x-axis of carthesian   
# 	        coordinates (in x-y-plane) in radiant
#	
# KEYWORD PARAMETERS:
#       transmat_inv:   Set this keyword if also the inverse of the
#                       transformation matrix shall be returned
#	help:	Set this keyword (flag) tho show the documentation of this procedure
#
# OUTPUTS:
#	transformation-matrix for coordinate transformation, transmat
#
# EXAMPLE:
#	CALC_TRANSMAT2, !PI/2., 0., transmat
#-
#
#
   if alpha > np.pi:
      print('error in calc_transmat2: alpha needs to be in radiant')
      exit()
#
   if gamma > 2.*np.pi:
      print('error in calc_transmat2: gamma needs to be in radiant')
      exit()
#
#---------------CALCULATION OF TRANSFORMATION MATRIX--------------------
#--------FOR SYSTEM (ex, ey, ez) TO SYSTEM (eex, eey, eez)--------------
#
   ex=np.zeros(3)
   ey=np.zeros(3)
   ez=np.zeros(3)
#
   eex=np.zeros(3)
   eey=np.zeros(3)
   eez=np.zeros(3)
#
   transmat=np.zeros(shape=(3,3))
#
   print('----------calculate transformation matrix--------------')
   print('alpha= ', alpha)
   print('gamma= ', gamma)
#
#calculate unit vectors of new coordinate system:  ez=e_r, ex=e_phi, ey=e_theta
#or, to get correct azimuthal velocity directions: ez=e_r, ex_e_theta, ey=e_phi
   ex=np.array([ -np.sin(gamma), np.cos(gamma), 0. ])
   ey=np.array([np.cos(alpha)*np.cos(gamma), np.cos(alpha)*np.sin(gamma), -np.sin(alpha) ])
   ez=np.array([np.sin(alpha)*np.cos(gamma), np.sin(alpha)*np.sin(gamma), np.cos(alpha)])
#
   if np.abs(ex[0]) < 1e-8: ex[0]=0.
   if np.abs(ex[1]) < 1e-8: ex[1]=0.
   if np.abs(ex[2]) < 1e-8: ex[2]=0.
   if np.abs(ey[0]) < 1e-8: ey[0]=0.
   if np.abs(ey[1]) < 1e-8: ey[1]=0.
   if np.abs(ey[2]) < 1e-8: ey[2]=0.
   if np.abs(ez[0]) < 1e-8: ez[0]=0.
   if np.abs(ez[1]) < 1e-8: ez[1]=0.
   if np.abs(ez[2]) < 1e-8: ez[2]=0.
#
#check for orthogonality
   check1=np.dot(ex,ey)
   check2=np.dot(ex,ez)
   check3=np.dot(ey,ez)
#
   if np.abs(check1) > 1.e-7:
      print(check1)
      print('error in calc_transmat2: ex, ey not orthogonal')
      exit()
   if np.abs(check2) > 1.e-7:
      print(check2)
      print('error in calc_transmat2: ex, ez not orthogonal')
      exit()
   if np.abs(check3) > 1.e-7:
      print(check3)
      print('error in calc_transmat2: ey, ez not orthogonal')
      exit()      
#
#----transformation matrix from system (ex,ey,ez) to (eex, eey, eez)----
#
   eex = np.array([ 1., 0., 0. ])
   eey = np.array([ 0., 1., 0. ])
   eez = np.array([ 0., 0., 1. ])
#
   transmat = np.array([[np.dot(ex,eex), np.dot(ex,eey), np.dot(ex,eez)],
                          [np.dot(ey,eex), np.dot(ey,eey), np.dot(ey,eez)], 
                          [np.dot(ez,eex), np.dot(ez,eey), np.dot(ez,eez)]])
   transmat_inv = np.transpose(transmat)
#
   print('unit vectors:')
   print('e_x_slice', ex)
   print('e_y_slice', ey)
   print('e_z_slice', ez)
   print('')
   print('transformation matrix:')
   print(transmat)
   print('')
   print('inverse transformation matrix:')
   print(transmat_inv)
   print('')

   return transmat, transmat_inv
#
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
#
def angles_spc(x, y, z):
#
#+
# NAME:
#	angles_spc
#
# PURPOSE:
#	This procedure calculates the latitude and polar angle of a vector
#       given in carthesian coordinates.
#
# CALLING SEQUENCE:
#	angles_spc, x, y, z, theta, phi
#
# INPUTS:
#	x:	x-coordinate of vector
#	y:	y-coordinate of vector
#	z:	z-coordinate of vector
#	
# KEYWORD PARAMETERS:
#	help:	Set this keyword (flag) tho show the documentation of this proce#ure
#
# OUTPUTS:
#	theta:   Latitude-angle of vector
#       phi:     Azimuth-angle of vector
#
# EXAMPLE:
#	angles_spc, 1., 1., 1., theta, phi
#-
#
#
#-----------------------CALCULATION OF THETA----------------------------
#
   rad = np.sqrt(x*x+y*y+z*z)
#
   if rad != 0.:
      theta = np.arccos(z/rad)
   else:
      theta = 0.
#
#------------------------CALCULATION OF PHI-----------------------------
#
   if z == rad:
      phi=0.
      return theta, phi
#
#
   if x == 0.:
      if y == 0.:
         phi = 0.0
      elif y > 0.:
         phi = np.pi/2.0
      else:
         phi = 3.0*np.pi/2.0
   elif y > 0 and x > 0.:
      phi = np.arctan(y/x) #first quadrant
   elif x > 0.:
      phi = 2.0*np.pi + np.arctan(y/x) #fourth quadrant
   else:
      phi = np.pi + np.arctan(y/x) #second and third quadrant
#
#-------------------------test if everything worked---------------------
#
   x_test=rad*np.sin(theta)*np.cos(phi)
   y_test=rad*np.sin(theta)*np.sin(phi)
   z_test=rad*np.cos(theta)
#
   if np.abs(x_test-x) > 1.e-8:
      print('error in get_angles_spc', x_test, x)
      exit()
   if np.abs(y_test-y) > 1.e-8:
      print('error in get_angles_spc', x_test, x)
      exit()
   if np.abs(z_test-z) > 1.e-8:
      print('error in get_angles_spc', x_test, x)
      exit()

   return theta, phi
#
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
#
def trad(xnue0,intensity):
   term1 = cgs_planck * xnue0 / cgs_kb
   term2 = 2.0 * cgs_planck * xnue0**3.0 / cgs_clight / cgs_clight / intensity + 1.0
   term2 = np.log(term2)
   trad = term1 / term2
   return trad
#
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
#
def dilfac(rstar, r):

   if r < rstar:
      print('error in dilfac: radius less than stellar radius')
      exit()

   x=(rstar/r)**2
   dilfac=0.5*(1.-np.sqrt(1.-x))
   return dilfac
#
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
#
def equivalent_width(xobs,ftot,fcont,xmin,xmax):
   ew = 0.
   nd=np.size(xobs)
   for i in range(1,nd):
      if xobs[i] >= xmin and xobs[i] <= xmax:
         ew = ew + 0.5*(ftot[i]/fcont[i] + ftot[i-1]/fcont[i-1] - 2.)*(xobs[i]-xobs[i-1])
   return ew
#
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
#
def conv_indx_1d_to_2d(indx_1d, ndx):

   indx_x = indx_1d % ndx
   indx_y = ((indx_1d - indx_x)/ndx)

   return int(indx_x), int(indx_y)
  
########################################################################

def phase_to_time(phase, eccentricity):
   from scipy import integrate
   #
   #transforms a phase-angle of an eccentric orbit
   #to the phase in time

   cc = (1.-eccentricity**2.)**1.5 / 2./np.pi

   nd = 101
   phases = np.linspace(0.,phase,nd)
   fct = 1./(1.+eccentricity*np.cos(phases))**2

   integral = integrate.simps(fct, phases)

   return cc*integral
      
