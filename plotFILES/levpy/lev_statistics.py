import numpy as np
from const_cgs import *
from lev_misc import *
from lev_interpol1d import *
#
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
#
def calc_mean(x):
   mean = np.sum(x)/np.float(np.size(x))
   
   return mean
#
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
#
def calc_variance(x, xmean):
#calculate unbiased sample variance given the mean
   nd=np.size(x)
   variance = 0.   
   for i in np.arange(nd):
       variance = variance + (x[i]-xmean)**2
       
   variance = variance/np.float(nd-1)
   
   return variance
#
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
#
def get_value_snr(mean, snr):
#
#calculates a value around the mean with given signal-to-noise
#assuming the noise is gaussian distributed with pdf
# f(x) = 1/sqrt(2.*pi*sigma^2) * exp(- (x-xmean)^2/2/sigma^2)
#   
   sigma = mean/snr
   xvalue = np.random.normal(mean, sigma)

   return xvalue


#
########################################################################
######################## wavelet transforms ############################
########################################################################
#
#
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
#
#
def calc_fwtransf(x1d,y1d, alim=[-6000.,6000.], blim=[1.,1000.],na=61,nb=101):
#
#calculate wavelet transform for mexican hat wavelet (forward transform)
#
#input: x1d     x-coordinate of function to be transformed
#       y1d     function values at the x-coordinates
#       alim    limits of the position-parameter of the wavelet transform
#       blim    limits of the scale-parameter of the wavelet transform
#       na, nb  number of data points to sample the scale-position plane
#
#output: a1d    the coordinates of the position-parameter
#        b1d    the coordinates of the scale-parameter
#        wtransform2d   the wavelet transform at each a1d,b1d
#
#calculate the grid
   amin = alim[0]
   amax = alim[1]
   a1d = np.linspace(amin,amax,na)

   bmin = blim[0]
   bmax = blim[1]
#   b1d = np.linspace(bmin,bmax,nb)
   b1d = grid_log(bmin,bmax,nb)

#nd: number of data points to sample each wavelet
   nd=101

#nx: number of data points used for the input function
   nx=np.size(x1d)
#
   wtransform2d=np.zeros(shape=(nb,na))

   for i in np.arange(na):   
      for j in np.arange(nb):
#calculate wavelet
         wmin = np.max([a1d[i]-10.*b1d[j],np.min(x1d)])
         wmax = np.min([a1d[i]+10.*b1d[j],np.max(x1d)])
         wx1d = np.linspace(wmin,wmax,nd)
         fdum = (wx1d-a1d[i])/b1d[j]
#         print(wmin,wmax)         
         wavelet1d = (1.-fdum**2)*np.exp(-0.5*fdum**2)/b1d[j]
#perform integration
         iim1, ii = find_indx(wx1d[0], x1d, nx)
         ykkm1 = interpol1d_2p_lin(y1d[iim1], y1d[ii], x1d[iim1], x1d[ii], wx1d[0])
         fkkm1 = wavelet1d[0]*ykkm1
         sum=0.
         for k in np.arange(1,nd):
            iim1, ii = find_indx(wx1d[k], x1d, nx)
            ykk = interpol1d_2p_lin(y1d[iim1], y1d[ii], x1d[iim1], x1d[ii], wx1d[k])
            fkk = wavelet1d[k]*ykk
            dx = wx1d[k]-wx1d[k-1]
            sum = sum + 0.5*(fkkm1+fkk)*dx
            fkkm1 = fkk
         wtransform2d[j][i] = sum*np.pi/8.
#         print(norm)

#
   return a1d, b1d, wtransform2d
#
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
#
def calc_bwtransf(a1d, b1d, wtransform2d):
#
#calculate backward wavelet transform for mexican hat wavelet for reproduction
#    of the signal from the wavelet transform
#
#can be used as benchmarking only if nb and range of b1d is very large
#
#input: a1d    the coordinates of the position-parame
#       b1d    the coordinates of the scale-parameter
#       wtransform2d   the wavelet transform at each a1d,b1d
#
#
#output: y1d    signal at each a1d
#
# 
   na=np.size(a1d)
   nb=np.size(b1d)
   y1d=np.zeros(na)

   for i in np.arange(na):
#      print(a1d[i])
#perfrom integration over b
      sum = 0
      fim1 = wtransform2d[0][i]/b1d[0]
#      print(b1d[0], fim1)
      for j in np.arange(1,nb):
         fii = wtransform2d[j][i]/b1d[j]
         dx = b1d[j]-b1d[j-1]
         sum = sum + 0.5*(fim1+fii)*dx
#         print(b1d[j], fii, dx, sum)         
         fim1 = fii
      y1d[i] = sum

#   exit()
#
   return y1d
#
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
#
def calc_power_scale(a1d, b1d, wtransform2d, blim=[0.1,1.]):
#
#calculate the power of the signal for a given range of scales blim
#
#input: a1d    the coordinates of the position-parame
#       b1d    the coordinates of the scale-parameter
#       wtransform2d   the wavelet transform at each a1d,b1d
#
#output: y1d   the power at each a1d
   
#
#calculate the integration grid
   na=np.size(a1d)
   nb=np.size(b1d)
   y1d=np.zeros(na)

   bmin = blim[0]
   bmax = blim[1]

   iim1l, iil = find_indx(bmin, b1d, nb)
   iim1u, iiu = find_indx(bmax, b1d, nb)   

   for i in np.arange(na):
#perfrom integration over b
      fiim1 = interpol1d_2p_lin(wtransform2d[iim1l][i], wtransform2d[iil][i], b1d[iim1l], b1d[iil], bmin)**2
      fii = wtransform2d[iil][i]**2
      dx = b1d[iil]-bmin
      sum = 0.5*(fiim1+fii)*dx
      for j in np.arange(iil+1,iim1u+1):
         fiim1 = fii
         fii = wtransform2d[j][i]**2
         dx = b1d[j]-b1d[j-1]
         sum = sum + 0.5*(fiim1+fii)*dx
      fiim1 = wtransform2d[iim1u][i]**2
      fii = interpol1d_2p_lin(wtransform2d[iim1u][i], wtransform2d[iiu][i], b1d[iim1u], b1d[iiu], bmax)**2
      dx = bmax - b1d[iim1u]
      sum = sum + 0.5*(fiim1+fii)*dx
#      
      y1d[i] = sum/(bmax-bmin)

   return y1d
#
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
#
def calc_power_loc(a1d, b1d, wtransform2d, alim=[-2.,-1.]):
#
#calculate the power of the signal for a given range of locations alim
#
#input: a1d    the coordinates of the position-parame
#       b1d    the coordinates of the scale-parameter
#       wtransform2d   the wavelet transform at each a1d,b1d
#
#output: y1d   the power at each a1d
#
#calculate the integration grid
   na=np.size(a1d)
   nb=np.size(b1d)
   y1d=np.zeros(nb)

   amin = alim[0]
   amax = alim[1]

   iim1l, iil = find_indx(amin, a1d, na)
   iim1u, iiu = find_indx(amax, a1d, na)   

   for j in np.arange(nb):
#perfrom integration over a
      fiim1 = interpol1d_2p_lin(wtransform2d[j][iim1l], wtransform2d[j][iil], a1d[iim1l], a1d[iil], amin)**2
      fii = wtransform2d[j][iil]**2
      dx = a1d[iil]-amin
      sum = 0.5*(fiim1+fii)*dx
      for i in np.arange(iil+1,iim1u+1):
         fiim1 = fii
         fii = wtransform2d[j][i]**2
         dx = a1d[i]-a1d[i-1]
         sum = sum + 0.5*(fiim1+fii)*dx
      fiim1 = wtransform2d[j][iim1u]**2
      fii = interpol1d_2p_lin(wtransform2d[j][iim1u], wtransform2d[j][iiu], a1d[iim1u], a1d[iiu], amax)**2
      dx = amax - a1d[iim1u]
      sum = sum + 0.5*(fiim1+fii)*dx
      
      y1d[j] = sum/(amax-amin)

#   for j in np.arange(nb):
##perfrom integration over a
#      fiim1 = wtransform2d[j][0]**2
#      sum=0.
#      for i in np.arange(1,na):
#         fii = wtransform2d[j][i]**2
#         dx = a1d[i]-a1d[i-1]
#         sum = sum + 0.5*(fiim1+fii)*dx
#         fiim1 = fii
#
#      y1d[j] = sum/(amax-amin)

   return y1d
