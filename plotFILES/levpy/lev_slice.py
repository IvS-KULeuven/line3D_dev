import numpy as np
from lev_misc import *
from lev_interpol3d import *

def slice_cyc(x, y, z, arr3d, npr=101 ,nzeta=81, n_vec=[0.,0.,1.], offset=0., pmin=0.01, pmax=10., rmin=1., rmax=10.):
#calculates a circular slice with normal-vector n_vec

   zeta_arr = grid_equi(0.,2.*np.pi,nzeta)
   p_arr = grid_log(pmin,pmax,npr)

   arr2d = np.zeros(shape=(nzeta,npr))
   xcoord = np.zeros(shape=(nzeta,npr))
   ycoord = np.zeros(shape=(nzeta,npr))

#transformation matrix
   alpha, gamma = get_angles_spc(n_vec[0], n_vec[1], n_vec[2])
   transmat, transmat_inv = calc_transmat2(alpha, gamma)

   vec_slice = np.zeros(3)
   vec_cac = np.zeros(3)

   for i in range(0,npr):
      for j in range(0,nzeta):
         xcoord[j][i] = p_arr[i]*np.sin(zeta_arr[j])
         ycoord[j][i] = p_arr[i]*np.cos(zeta_arr[j])

         vec_slice[0] = xcoord[j][i]
         vec_slice[1] = ycoord[j][i]
         vec_slice[2] = offset
         vec_cac = np.dot(transmat,vec_slice)
         radp=np.sqrt(vec_cac[0]*vec_cac[0] + vec_cac[1]*vec_cac[1] + vec_cac[2]*vec_cac[2])
#need to make sure that values on axes are treated the same
         if np.abs(vec_cac[0]) < 1.e-5: vec_cac[0]=0.
         if np.abs(vec_cac[1]) < 1.e-5: vec_cac[1]=0.
         if np.abs(vec_cac[2]) < 1.e-5: vec_cac[2]=0.

         if radp < rmin:
            arr2d[j][i]=-1.e10
         elif radp > rmax:
            arr2d[j][i]=-1.e10
         else:
#find indices of cube-vertices for interpolation
            ix1, ix2, iy1, iy2, iz1, iz2, expol = get_xyz_indx(vec_cac[0], vec_cac[1], vec_cac[2], x, y, z, np.size(x), np.size(y), np.size(z))
#
#----------------------interpolate 3d arrays onto slice-----------------
#
#density
            x1, x2, y1, y2, z1, z2, rada, radb, radc, radd, rade, radf, radg, radh, radp, llogx, llogy, llogz = \
                           get_xyz_values1(vec_cac[0],vec_cac[1],vec_cac[2], x, y, z, np.size(x), np.size(y), np.size(z), ix1, ix2, iy1, iy2, iz1, iz2)
            vala, valb, valc, vald, vale, valf, valg, valh, llogf = get_xyz_values2(np.size(x), np.size(y), np.size(z), arr3d, ix1, ix2, iy1, iy2, iz1, iz2)
            arr2d[j][i] = trilin_complete(vec_cac[0], vec_cac[1], vec_cac[2], x1, x2, y1, y2, z1, z2, \
                                       vala, valb, valc, vald, vale, valf, valg, valh, \
                                       rada, radb, radc, radd, rade, radf, radg, radh, radp, \
                                       expol, llogx=llogx, llogy=llogy, llogz=llogz, llogf=llogf, lr2=1)

   return arr2d, xcoord, ycoord
