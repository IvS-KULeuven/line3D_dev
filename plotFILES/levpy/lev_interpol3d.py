import numpy as np
import lev_interpol1d
#
def get_xyz_indx(xin, yin, zin, xarr, yarr, zarr, ndxmax, ndymax, ndzmax):
#
#-----------------------------------------------------------------------
#
#  finds indices of grid for a cube surrounding the point (xin,yin,zin)
#
#  input:  coordinates of point:             xin, yin, zin
#          dimension of grid:                ndxmax, ndymax, ndzmax
#          x,y,z-grid:                       xarr, yarr, zarr
#  output: indices of x-grid:                indx_x1, indx_x2
#          indices of y-grid:                indx_y1, indx_y2
#          indices of z-grid:                indx_z1, indx_z2
#          flag if extrapolation is needed:  expol
#
#-----------------------------------------------------------------------
#
   rmax=np.max(xarr)

   if xin >= 0.:
      startx=ndxmax-2
      endx=0
      alpha=-1
   else:
      startx=1
      endx=ndxmax-1
      alpha=1
#
   if yin >= 0.:
      starty=ndymax-2
      endy=0
      beta=-1
   else:
      starty=1
      endy=ndymax-1
      beta=1
#
   if zin >= 0.:
      startz=ndzmax-2
      endz=0
      gamma=-1
   else:
      startz=1
      endz=ndzmax-1
      gamma=1
#
   indx_x1=startx
   indx_x2=startx-alpha
   for i in range(startx,endx,alpha):
      if alpha*xarr[i] >= alpha*xin:
         indx_x1=i
         indx_x2=i-alpha
         break
#
   indx_y1=starty
   indx_y2=starty-beta
   for i in range(starty, endy, beta):
      if beta*yarr[i] >= beta*yin:
         indx_y1=i
         indx_y2=i-beta
         break
#
   indx_z1=startz
   indx_z2=startz-gamma
   for i in range(startz, endz, gamma):
      if gamma*zarr[i] >= gamma*zin:
         indx_z1=i
         indx_z2=i-gamma
         break

#
#--------------extrapolation for grid-points near photosphere-----------
#
#check if extrapolation is needed at inner part of the star
   rad=np.sqrt(xarr[indx_x1]*xarr[indx_x1]+yarr[indx_y1]*yarr[indx_y1]+zarr[indx_z1]*zarr[indx_z1])
#
   if rad < 1.:
      expol=1
      for i in range(1,10):
         indx_x1=indx_x1-alpha
         indx_x2=indx_x2-alpha
         indx_y1=indx_y1-beta
         indx_y2=indx_y2-beta
         indx_z1=indx_z1-gamma
         indx_z2=indx_z2-gamma
         rad=np.sqrt(xarr[indx_x1]*xarr[indx_x1]+yarr[indx_y1]*yarr[indx_y1]+zarr[indx_z1]*zarr[indx_z1])
         if rad >= 1.:
            break
      if rad < 1.:
         print('error in get_xyz_indx: rad lt 1.d0 => extrapolation over more than 10 grid points')
         exit()

      return indx_x1, indx_x2, indx_y1, indx_y2, indx_z1, indx_z2, expol
   else:
      expol=0
#
#-------------extrapolation for grid points larger than rmax------------
#
#check if extrapolation is needed at outer part of the star
   rad=np.sqrt(xarr[indx_x2]*xarr[indx_x2]+yarr[indx_y2]*yarr[indx_y2]+zarr[indx_z2]*zarr[indx_z2])
#   
   if rad > rmax:
#by default: extrapolation
      expol1=1
      expol2=1
      expol3=1
      if np.abs(xin) < 1.e-6  and np.abs(yin) < 1.e-6:
         print('point is on z-axis => no extrapolation needed')
         expol1=0
      if np.abs(xin) < 1.e-6 and np.abs(zin) < 1.e-6:
         print('point is on y-axis => no extrapolation needed')
         expol2=0
      if np.abs(yin) < 1.e-6 and np.abs(zin) < 1.e-6:
         print('point is on x-axis => no extrapolation needed')
         expol3=0
      if expol1 == 0 or expol2 == 0 or expol3 == 0:
         expol=0
      else:
         expol=1
         for i in range(1,10):
            indx_x1=indx_x1+alpha
            indx_x2=indx_x2+alpha
            indx_y1=indx_y1+beta
            indx_y2=indx_y2+beta
            indx_z1=indx_z1+gamma
            indx_z2=indx_z2+gamma
            rad=np.sqrt(xarr[indx_x2]*xarr[indx_x2]+yarr[indx_y2]*yarr[indx_y2]+zarr[indx_z2]*zarr[indx_z2])
            if rad <= rmax:
               break
         if rad > rmax:
            print('error in get_xyz_indx: linfo_max eq false => extrapolation over more than 10 grid points')
   else:
      expol=0

   return indx_x1, indx_x2, indx_y1, indx_y2, indx_z1, indx_z2, expol
#
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
#
def get_rtp_indx(rin, thetain, phiin, rarr, tharr, phiarr, nr, ntheta, nphi):
#
#-----------------------------------------------------------------------
#
#  finds indices of grid for a cube surrounding the point (rin,thetain,phiin)
#
#  input:  coordinates of point:             xin, yin, zin
#          dimension of grid:                ndxmax, ndymax, ndzmax
#          x,y,z-grid:                       xarr, yarr, zarr
#  output: indices of x-grid:                indx_x1, indx_x2
#          indices of y-grid:                indx_y1, indx_y2
#          indices of z-grid:                indx_z1, indx_z2
#          flag if extrapolation is needed:  expol
#
#-----------------------------------------------------------------------
#
   expol=0
#
   rmin=np.min(rarr)
   rmax=np.max(rarr)
#
   indx_r1=0
   indx_r2=1
   for i in range(1,nr):
      if rarr[i] >= rin:
         indx_r1=i-1
         indx_r2=i
         break
#
   if rin > rmax:
      indx_r1=nr-2
      indx_r2=nr-1
      expol=1
#
   if rin < rmin:
      indx_r1=0
      indx_r2=1
      expol=1
#
#;
#
   thmin=np.min(tharr)
   thmax=np.max(tharr)
   indx_th1=0
   indx_th2=1
   for i in range(1,ntheta):
      if tharr[i] >= thetain:
         indx_th1=i-1
         indx_th2=i
         break
#
   if thetain > thmax:
      indx_th1=ntheta-2
      indx_th2=ntheta-1
      expol=1
#
   if thetain < thmin:
      indx_th1=0
      indx_th2=1
      expol=1
#
#
   phimin=np.min(phiarr)
   phimax=np.max(phiarr)
   indx_phi1=0
   indx_phi2=1
   for i in range(1,nphi):
      if phiarr[i] >= phiin:
         indx_phi1=i-1
         indx_phi2=i
         break
#
   if phiin > phimax:
      indx_phi1=nphi-2
      indx_phi2=nphi-1
      expol=1
#
   if phiin < phimin:
      indx_phi1=0
      indx_phi2=1
      expol=1

   return indx_r1, indx_r2, indx_th1, indx_th2, indx_phi1, indx_phi2, expol
#
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
#
def get_xyz_values1(x_coord, y_coord, z_coord, xarr, yarr, zarr, ndxmax, ndymax, ndzmax,
                     indx_x1, indx_x2, indx_y1, indx_y2, indx_z1, indx_z2):
#
#-----------------------------------------------------------------------
#
#   get coordinates and radii of a cube
#   with vertices given by indx_x1, ... indx_z2
#
#input: dimension of arrays x,y,z: ndxmax, ndymax, ndzmax
#       arrays x,y,z
#       coordinates of a point inside cube: x_coord, y_coord, z_coord
#       indices of cube-vertices: indx_x1, ... indx_z2
#
#output: grid value: x1, x2, y1, y2, z1, z2
#        radii of vertices and of point: rada, ... radh, radp
#        flags to decide if logarithmic interpolation is allowed:
#               llogx, llogy, llogz
#
#-----------------------------------------------------------------------
#
   x1=xarr[indx_x1]
   x2=xarr[indx_x2]
   y1=yarr[indx_y1]
   y2=yarr[indx_y2]
   z1=zarr[indx_z1]
   z2=zarr[indx_z2]
#
   rada = np.sqrt(x1*x1 + y1*y1 + z1*z1)
   radb = np.sqrt(x2*x2 + y1*y1 + z1*z1)
   radc = np.sqrt(x1*x1 + y1*y1 + z2*z2)
   radd = np.sqrt(x2*x2 + y1*y1 + z2*z2)
   rade = np.sqrt(x1*x1 + y2*y2 + z1*z1)
   radf = np.sqrt(x2*x2 + y2*y2 + z1*z1)
   radg = np.sqrt(x1*x1 + y2*y2 + z2*z2)
   radh = np.sqrt(x2*x2 + y2*y2 + z2*z2)
#
   radp = np.sqrt(x_coord*x_coord + y_coord*y_coord + z_coord*z_coord)
#
#-------------check if logarithmic interpolation is allowed-------------
#
   if x1*x2 <= 0.:
      llogx=0.
   else:
      if x1*x_coord <= 0.:
         llogx=0
      else:
         llogx=1
   
   if y1*y2 <= 0.:
      llogy=0
   else:
      if y1*y_coord <= 0.:
         llogy=0
      else:
         llogy=1
   
   if z1*z2 <= 0.:
      llogz=0
   else:
      if z1*z_coord <= 0.:
         llogz=0
      else:
         llogz=1
   
   return x1, x2, y1, y2, z1, z2, rada, radb, radc, radd, rade, radf, radg, radh, radp, llogx, llogy, llogz
#
#-----------------------------------------------------------------------
#
def get_xyz_values2(ndxmax, ndymax, ndzmax, yvalue3d,
                    indx_x1, indx_x2, indx_y1, indx_y2, indx_z1, indx_z2):
#
#-----------------------------------------------------------------------
#
#   get physical quantities of a cube 
#   with vertices given by indx_x1, ... indx_z2
#
#input: dimension of 3-d array: ndxmax, ndymax, ndzmax
#       physical value on grid: yvalue3d
#       indices of cube-vertices: indx_x1, ... indx_z2
#
#output: physical value on cube vertices: vala ... valh
#        flag to decide if logarithmic interpolation is allowed:
#               llogf
#
#-----------------------------------------------------------------------
#
#   vala = yvalue3d[indx_x1,indx_y1,indx_z1]
#   valb = yvalue3d[indx_x2,indx_y1,indx_z1]
#   valc = yvalue3d[indx_x1,indx_y1,indx_z2]
#   vald = yvalue3d[indx_x2,indx_y1,indx_z2]
#   vale = yvalue3d[indx_x1,indx_y2,indx_z1]
#   valf = yvalue3d[indx_x2,indx_y2,indx_z1]
#   valg = yvalue3d[indx_x1,indx_y2,indx_z2]
#   valh = yvalue3d[indx_x2,indx_y2,indx_z2]
   
   vala = yvalue3d[indx_z1][indx_y1][indx_x1]
   valb = yvalue3d[indx_z1][indx_y1][indx_x2]
   valc = yvalue3d[indx_z2][indx_y1][indx_x1]
   vald = yvalue3d[indx_z2][indx_y1][indx_x2]
   vale = yvalue3d[indx_z1][indx_y2][indx_x1]
   valf = yvalue3d[indx_z1][indx_y2][indx_x2]
   valg = yvalue3d[indx_z2][indx_y2][indx_x1]
   valh = yvalue3d[indx_z2][indx_y2][indx_x2]   
#
   if vala<=0. or valb<=0. or valc<=0. or vald<=0. or vale<=0. or valf<=0. or valg<=0. or valh<=0.:
      llogf=0
   else:
      llogf=1

   return vala, valb, valc, vald, vale, valf, valg, valh, llogf
#
#-----------------------------------------------------------------------
#
def trilin(xout, yout, zout, 
           x1, x2, y1, y2, z1, z2, 
           vala, valb, valc, vald, vale, valf, valg, valh, 
           llogx, llogy, llogz, llogf):

   if llogf == 1:
#
#prepare input values for log-* interpolation
      dum_vala=np.log10(vala)
      dum_valb=np.log10(valb)
      dum_valc=np.log10(valc)
      dum_vald=np.log10(vald)
      dum_vale=np.log10(vale)
      dum_valf=np.log10(valf)
      dum_valg=np.log10(valg)
      dum_valh=np.log10(valh)
#
      if llogz == 1:
#log-log interpolation in z-direction
         dum_z1=np.log10(np.abs(z1))
         dum_z2=np.log10(np.abs(z2))
         dum_zout=np.log10(np.abs(zout))
      else:
#log-lin interpolation in z-direction
         dum_z1=z1
         dum_z2=z2
         dum_zout=zout
#
      if llogx == 1:
#log-log-interpolation in x-direction
         dum_x1=np.log10(np.abs(x1))
         dum_x2=np.log10(np.abs(x2))
         dum_xout=np.log10(np.abs(xout))
      else:
#log-lin-interpolation in x-direction
         dum_x1=x1
         dum_x2=x2
         dum_xout=xout
#
      if llogy == 1:
#log-log-interpolation in y-direction
         dum_y1=np.log10(np.abs(y1))
         dum_y2=np.log10(np.abs(y2))
         dum_yout=np.log10(np.abs(yout))
      else:
#log-lin-interpolation in y-direction
         dum_y1=y1
         dum_y2=y2
         dum_yout=yout
#
#perform interpolation in z-direction
      yvalue_s0=lev_interpol1d.interpol_yp(dum_z2, dum_z1, dum_vald, dum_valb, dum_zout)
      yvalue_s1=lev_interpol1d.interpol_yp(dum_z2, dum_z1, dum_valc, dum_vala, dum_zout)
      yvalue_n0=lev_interpol1d.interpol_yp(dum_z2, dum_z1, dum_valh, dum_valf, dum_zout)
      yvalue_n1=lev_interpol1d.interpol_yp(dum_z2, dum_z1, dum_valg, dum_vale, dum_zout)
#
#perform interpolation in x-direction
      soln1=lev_interpol1d.interpol_yp(dum_x2, dum_x1, yvalue_n0, yvalue_n1, dum_xout)
      sols1=lev_interpol1d.interpol_yp(dum_x2, dum_x1, yvalue_s0, yvalue_s1, dum_xout)
#
#perform interpolation in y-direction
      yinterp=lev_interpol1d.interpol_yp(dum_y2, dum_y1, soln1, sols1, dum_yout)
      yinterp = 10.**yinterp
#
   else:
#lin-lin-interpolation

      yvalue_s0=lev_interpol1d.interpol_yp(z2, z1, vald, valb, zout)
      yvalue_s1=lev_interpol1d.interpol_yp(z2, z1, valc, vala, zout)
      yvalue_n0=lev_interpol1d.interpol_yp(z2, z1, valh, valf, zout)
      yvalue_n1=lev_interpol1d.interpol_yp(z2, z1, valg, vale, zout)
#
      soln1=lev_interpol1d.interpol_yp(x2, x1, yvalue_n0, yvalue_n1, xout)
      sols1=lev_interpol1d.interpol_yp(x2, x1, yvalue_s0, yvalue_s1, xout)
#
      yinterp=lev_interpol1d.interpol_yp(y2, y1, soln1, sols1, yout)

   return yinterp
#
#-----------------------------------------------------------------------
#
def trilin_complete(xout, yout, zout, 
                    x1, x2, y1, y2, z1, z2,
                    vala, valb, valc, vald, vale, valf, valg, valh,
                    rada, radb, radc, radd, rade, radf, radg, radh, radp,
                    expol, llogx=0, llogy=0, llogz=0, llogf=0, lr2=0):
#
   if lr2 == 1:
#applying interpolation of function values * r^2
      vala=vala*rada*rada
      valb=valb*radb*radb
      valc=valc*radc*radc
      vald=vald*radd*radd
      vale=vale*rade*rade
      valf=valf*radf*radf
      valg=valg*radg*radg
      valh=valh*radh*radh

#perform only interpolation (expol=.false.)
   if expol == 0:
      yinterp=trilin(xout, yout, zout, \
                     x1, x2, y1, y2, z1, z2, \
                     vala, valb, valc, vald, vale, valf, valg, valh, \
                     llogx, llogy, llogz, llogf)
   else:
#set values to zero if extrapolation is needed (or perform extrapolation)
      yinterp=0.
      yinterp = trilin(xout, yout, zout, \
                       x1, x2, y1, y2, z1, z2, \
                       vala, valb, valc, vald, vale, valf, valg, valh, \
                       llogx, llogy, llogz, llogf)

   if lr2 == 1:
      yinterp=yinterp/radp/radp

   return yinterp
#
#-----------------------------------------------------------------------
#
def interpol3d_8p_lin(f_im1jm1km1, f_ijm1km1, f_im1jkm1, f_ijkm1, 
                      f_im1jm1k, f_ijm1k, f_im1jk, f_ijk, 
                      x_im1, x_i, y_jm1, y_j, z_km1, z_k, x_p, y_p, z_p):
#
#         interpolates values given on a 3d grid onto point x_p, y_p, z_p
#                      using triliniear interpolation
#
#   f(x,y,z) = f_km1 + (f_k-f_km1)*(z_p-z_km1)/(z_k-z_km1)
#
#      with:   f_km1 = f_jm1km1 + (f_jkm1-f_jm1km1)*(y_p-y_im1)/(y_j-y_jm1)
#              f_k   = f_jm1    + (f_jk  -f_jm1k  )*(y_p-y_im1)/(y_j-y_jm1)
#              f_jk     = f_im1jkm1   + (f_ijkm1-f_im1jkm1) *   (x_p-x_im1)/(x_i-x_im1)
#              f_jm1km1 = f_im1jm1km1 + (f_ijm1km1-f_im1jm1km1)*(x_p-x_im1)/(x_i-x_im1)
#
#on input: 
#
#          f_im1jk-------------f_ijk          z_k
#            /|                 /|             |  y_j 
#           / |                / |             |   /
#          /  |               /  |             |  /
#         /   |              /   |             | /
#        /    |             /    |             |/
#       / f_im1jm1---------/--f_ijkm1        x_im1--------------x_i
#      /     /  x(x,y,z)  /     /            y_jm1
# f_im1jm1k-----------f_ijm1k  /             z_km1
#     |    /             |    /
#     |   /              |   /
#     |  /               |  /
#     | /                | /
#     |/                 |/
#f_im1jm1km1---------f_ijm1km1
#                    
#        x_p, y_p, z_p: coordinates of point onto which shall be interpolated
#
#on output: interpolated value at x_p, y_p, z_p
#
#define deltax, deltay
   dxi=x_i-x_im1
   dx=x_p-x_im1
   dyj=y_j-y_jm1
   dy=y_p-y_jm1
   dzk=z_k-z_km1
   dz=z_p-z_km1
#
#define deltax, deltay-ratios
   rdx=dx/dxi
   rdy=dy/dyj
   rdz=dz/dzk
   rdxdy=rdx*rdy
#
#trilinear interpolation
   fdum1=1.0-rdx-rdy+rdxdy
   fdum2=rdx-rdxdy
   fdum3=rdy-rdxdy
   fdum4=rdxdy
#
   ecoeff=fdum1*rdz
   fcoeff=fdum2*rdz
   gcoeff=fdum3*rdz
   hcoeff=fdum4*rdz
   acoeff=fdum1-ecoeff
   bcoeff=fdum2-fcoeff
   ccoeff=fdum3-gcoeff
   dcoeff=fdum4-hcoeff
#
   return acoeff*f_im1jm1km1 + bcoeff*f_ijm1km1 + ccoeff*f_im1jkm1 + dcoeff*f_ijkm1 + \
          ecoeff*f_im1jm1k   + fcoeff*f_ijm1k   + gcoeff*f_im1jk   + hcoeff*f_ijk
