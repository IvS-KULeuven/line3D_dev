from lev_interpol1d import *
import numpy as np

def get_linecenter(xobs,fnorm,flim=1.,xlim=5.):
#
#
#get the center of an emission line with the barycenter method
#input: xobs in vth_fiducial
#       fnorm: normalized flux
#
#flim: only consider points left and right of a certain height
#      specified by flim (in units of max(fnorm)-1.)
#xlim: only consider points with xobs in [xobs(flim),xmax]
#
#
   dfmax = np.max(fnorm)-1.0
   dfmin = 1.0-np.min(fnorm)
#
   if dfmax >= dfmin:
      labs=0 
      lemi=1  #emission profile
   else:
      labs=1  #absorption profile       
      lemi=0

   if lemi == 1:
      fp = 1.0+(np.max(fnorm)-1.0)*flim
   else:
      fp = 1.0-(1.0-min(fnorm))*flim
#
   nd = np.size(xobs)
#
#--------------find indices for left part of the profile----------------
#
#starting index
   indxl2=0
   for i in range(0,nd):
#      print, i, fnorm(i), fmax, fmin
      if fnorm[i] >= fp and lemi == 1:
         indxl2=i
         break
      if fnorm[i] <= fp and labs == 1:
         indxl2=i
         break
#
   if indxl2 == 0:
      print('error in get_linecenter: indxl2 = 0')
      exit()
#
#interpolate xobs to find limit
   xobsl2 = interpol1d_2p_lin(xobs[indxl2-1],xobs[indxl2],fnorm[indxl2-1],fnorm[indxl2], fp)
#
   xobsl1 = xobsl2 - xlim
#
   if xobsl1 >= xobsl2:
      print('error in get linecenter: xobsl1 > xobsl2')
      exit()
#
   if xobsl1 < np.min(xobs):
      print('error in get linecenter: xobsl1 < min(xobs)')
      exit
#
#define flux profile with doubled resolution
   indx = np.where(xobs >= xobsl1)
   xobs_dum = xobs[indx]
   indx = np.where(xobs_dum <= xobsl2)
   ndl = np.size(indx)*2
   xobsl = np.linspace(xobsl1,xobsl2,ndl)
   fnorml = np.zeros(ndl)
   for i in range(0,ndl):
      iim1, ii = find_indx(xobsl[i], xobs, nd)
      fnorml[i] = interpol1d_2p_lin(fnorm[iim1],fnorm[ii],xobs[iim1],xobs[ii],xobsl[i])
#
#--------------find indices for right part of the profile---------------
#
   if lemi == 1:
      fp = 1.0+(np.max(fnorm)-1.0)*flim
   else:
      fp = 1.0-(1.0-np.min(fnorm))*flim
##
##starting index
   indxr1=nd-1
   for i in range(nd-1,-1,-1):
      if fnorm[i] >= fp and lemi == 1:
         indxr1=i
         break
      if fnorm[i] <= fp and labs == 1:
         indxr1=i
         break
#
   if indxr1 == nd-1:
      print('error in get_linecenter: indxr1 = nd-1')
      exit()
#
#interpolate xobs to find limit
   xobsr1 = interpol1d_2p_lin(xobs[indxr1],xobs[indxr1+1],fnorm[indxr1],fnorm[indxr1+1],fp)
#
   xobsr2 = xobsr1 + xlim
#
   if xobsr2 <= xobsr1:
      print('error in get linecenter: xobsr2 < xobsr1')
      exit()
#
#print, max(xobs)
#print, xobsr2
#print, ''
   if xobsr2 > np.max(xobs):
      print('error in get linecenter: xobsr2 > max(xobs)')
      exit()
#
#define flux profile with doubled resolution
   indx = np.where(xobs >= xobsr1)
   xobs_dum = xobs[indx]
   indx = np.where(xobs <= xobsr2)
   ndr = np.size(indx)*2
   xobsr = np.linspace(xobsr1,xobsr2,ndr)
   fnormr = np.zeros(ndr)
   for i in range(0,ndr):
      iim1,ii = find_indx(xobsr[i], xobs, nd)
      fnormr[i] = interpol1d_2p_lin(fnorm[iim1],fnorm[ii],xobs[iim1],xobs[ii],xobsr[i])
#
#-------------------------get the barycenter----------------------------
#
#left part
   suml1=0.0
   suml2=0.0
   for i in range(1,ndl):
      suml1 = suml1 + 0.5*(xobsl[i-1]*fnorml[i-1]+xobsl[i]*fnorml[i])*(xobsl[i]-xobsl[i-1])
      suml2 = suml2 + 0.5*(fnorml[i-1]+fnorml[i])*(xobsl[i]-xobsl[i-1])
#
#right part
   sumr1=0.0
   sumr2=0.0
   for i in range(1,ndr):
      sumr1 = sumr1 + 0.5*(xobsr[i-1]*fnormr[i-1]+xobsr[i]*fnormr[i])*(xobsr[i]-xobsr[i-1])
      sumr2 = sumr2 + 0.5*(fnormr[i-1]+fnormr[i])*(xobsr[i]-xobsr[i-1])
#central part
   sumc1=0.5*(fnorml[ndl-1]*xobsl[ndl-1]+fnormr[0]*xobsr[0])*(xobsr[0]-xobsl[ndl-1])
   sumc2=0.5*(fnorml[ndl-1]+fnormr[0])*(xobsr[0]-xobsl[ndl-1])
#sumc1=0.d0
#sumc2=0.d0
#
#print, suml1, sumr1, sumc1
#print, suml2, sumr2, sumc2
   xcenter=(suml1+sumr1+sumc1)/(suml2+sumr2+sumc2)


   return xcenter, xobsl, xobsr, fnorml, fnormr
