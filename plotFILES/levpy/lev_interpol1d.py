import numpy as np

def interpol_yp(x1, x2, y1, y2, xp):
#
# name:
#	interpol_yp
#
# purpose:
#	this procedure interpolates linearly between two points.
#       
# inputs:
#	x1, x2:	coordinates of given points
#	y1, y2:	function values at given points
#	xp:     coordinate of unknown point
#	
#
# outputs:
#	yp:     function value at unknown point
#
#
#
#-----------------------------------------------------------------------
#
   yp = y2 + (y2-y1)/(x2-x1)*(xp-x2)
   return yp
#
########################################################################

def find_indx(xin,xarr,nd):
#
# NAME:
#       find_indx
# PURPOSE:
#       finds index for interpolation of two points
#
#       for interpolation: xarr(iim1)----xin-----xarr(ii)
#       for extrapolation: xin----xarr(iim1)-----xarr(ii)
#                          xarr(iim)----xarr(ii)-----xin
#
#
# INPUTS:
#       coordinate of point:             xin
#       dimension of grid:               nd
#       grid:                            xarr
#
# OUTPUTS: 
#       indices:   iim1, ii
#
#-----------------------------------------------------------------------
#
   if xin >= xarr[nd-1]:
      iim1=nd-2
      ii=nd-1
      return iim1, ii
#
   iim1=0
   ii=1
   for i in range(1,nd):
      if xarr[i] >= xin:
         ii=i
         iim1=i-1
         break

   return iim1, ii

########################################################################

def find_eindx(xin,xarr,nd):
#
# NAME:
#       find_eindx
# PURPOSE:
#       finds index for interpolation of two points for an equidistant grid
#
#       for interpolation: xarr(iim1)----xin-----xarr(ii)
#       for extrapolation: xin----xarr(iim1)-----xarr(ii)
#                          xarr(iim)----xarr(ii)-----xin
#
#
# INPUTS:
#       coordinate of point:             xin
#       dimension of grid:               nd
#       grid:                            xarr
#
# OUTPUTS: 
#       indices:   iim1, ii
#
#-----------------------------------------------------------------------
#
   if xin >= xarr[nd-1]:
      iim1=nd-2
      ii=nd-1
      return iim1, ii
#
   if xin <= xarr[0]:
      iim1=0
      ii=1
      return iim1, ii

   dx = xarr[1]-xarr[0]
   xmin = xarr[0]
   ii = np.int(np.ceil((xin-xmin)/dx))   
   iim1=ii-1

   return iim1, ii
#
########################################################################

def interpol1d_2p_lin(f_im1, f_i, x_im1, x_i, x_p):
#
#         interpolates values given on a 1d grid onto point x_p
#                 using linear interpolation
#
#   f(x) = f_im1 + (f_i-f_im1)*(x_p-x_im1)/(x_j-x_im1)
#
#on input: 
#
#          f_im1--------------f_i
#          x_im1--------------x_i
#                           
#        x_p: coordinates of point onto which shall be interpolated
#
#on output: interpolated value at x_p
#
#define deltax
   dxi=x_i-x_im1
   dx=x_p-x_im1
#
#define deltax-ratios
   rdx=dx/dxi
#
#linear interpolation
   acoeff=1.0-rdx
   bcoeff=rdx
#
   return acoeff*f_im1 + bcoeff*f_i

#
########################################################################
#
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
#
def coeff_typ_quadc(f_im1, f_i, f_ip1, x_im1, x_i, x_ip1):
#
#         interpolation coefficients on a 1d-axis
#                 using quadratic function
#
#   f(x) = a*(x-xi)^2 + b*(x-xi) + c*(x-xi)
#
#on input: 
#
#                    f_im1----.-----f_i--------f_ip1
#                    x_im1---x_p----x_i--x_p---x_ip1
#
#        x_p: coordinate of point onto which shall be interpolated
#
#on output: coefficients acoeff, bcoeff, ccoeff
#
#

   delx_i=x_i-x_im1
   delx_ip1=x_ip1-x_i
   delf_i=f_i-f_im1
   delf_ip1=f_ip1-f_i
   
#calculate coefficients (numerically not the fastests method, but for the moment...)
   acoeff=(delf_ip1/delx_ip1 - delf_i/delx_i)/(delx_ip1+delx_i)
   bcoeff=(delf_ip1*delx_i/delx_ip1 + delf_i*delx_ip1/delx_i)/(delx_ip1+delx_i)
   ccoeff=f_i

   return acoeff, bcoeff, ccoeff
