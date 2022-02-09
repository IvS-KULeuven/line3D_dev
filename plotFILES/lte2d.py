import matplotlib.pyplot as plt
import numpy as np
from const_cgs import *
from lev_contour import *
from lev_misc import conv_indx_1d_to_2d, readcol

from matplotlib.colors import BoundaryNorm

def read_lte2d(fname=None):

#dirty, but working
   f = open(fname, 'r')
   
   header = f.readline()
   na = f.readline()
   na = na.strip()
   na = na.split()
   na = float(na[0])
   
   nlines = f.readline()
   nlines = nlines.strip()
   nlines = nlines.split()
   nlines = int(nlines[0])
   
   ntemp = f.readline()
   ntemp = ntemp.strip()
   ntemp = ntemp.split()
   ntemp = int(ntemp[0])

   nrho = f.readline()
   nrho = nrho.strip()
   nrho = nrho.split()  
   nrho = int(nrho[0])

   glower = f.readline()
   glower = glower.strip()
   glower = glower.split()
   glower = float(glower[0])

   line = f.readline()
   line = f.readline()
   for i in np.arange(nlines):
      line = f.readline()

   line = f.readline()
   line = f.readline()
   line = f.readline()

   #temperature
   line = f.readline()
   line = line.strip()
   columns = line.split()

   t1d=np.zeros(ntemp)
   for i in np.arange(1,ntemp+1):
      t1d[i-1] = float(columns[i])

   rho1d=np.zeros(ntemp)
   nlower2d=np.zeros(shape=(nrho,ntemp))
   for i in np.arange(nrho):
      line = f.readline()
      line = line.strip()
      columns = line.split()
      rho1d[i] = float(columns[0])
      for j in np.arange(ntemp):
         nlower2d[i][j] = columns[j+1]

   nlower2d=np.transpose(nlower2d)

   return rho1d, t1d, nlower2d
#
def main(fnames='../lte_tables/Y09800/06_03_09.txt', clim=[-10.,12.], windx=0):
#
#
   r_model, temp_model, rho_model, opalbar_model, ltemp_model, lrho_model, lopalbar_model = readcol('../TRASH/for_luka3.dat', ncols=7,nskip=0)

   fig = plt.figure(windx, constrained_layout=True, figsize=(8,8))
#
   plt.ion()
   plt.show()
#
   ax = fig.subplots(4,4)

   clevels, ticks = get_clevels(clim=clim)

   i=0
   j=0
   k=0
   for fname in fnames:
      rho1d, t1d, nlower2d = read_lte2d(fname)

      i,j = conv_indx_1d_to_2d(k, 4)
      k=k+1

      ax[j,i].set_title(fname[-12:-4])
      ax[j,i].set_xlabel(r'$\log(\rho)$')
      ax[j,i].set_ylabel(r'$\log(T)$')


#     clevels, ticks = get_clevels(clim=[np.min(nlower2d),np.max(nlower2d)])
      cmap_lev = get_cmap('jet', upper='black', lower='grey')
      contourplot = ax[j,i].pcolormesh(rho1d, t1d, nlower2d, vmin=clim[0], vmax=clim[1], cmap=cmap_lev, shading='auto')
      ax[j,i].plot(lrho_model,ltemp_model, color='black', lw=1)

   cbar = fig.colorbar(contourplot,ax=ax[0:4,3],orientation='vertical')
   cbar.set_label(r'$\log(n_l)$')

   
#   oname1 = oname2d+'.png'
#   oname2 = oname2d+'.ps'
#   fig2.savefig(oname1, bbox_inches='tight')
#   fig2.savefig(oname2, bbox_inches='tight')

#############################################################################
#---------------------plot everything for a certain model--------------------
#############################################################################

dir_base = '../lte_tables/Y09800'

fnames_hei=[dir_base+'/02_01_01.txt',
           dir_base+'/02_01_02.txt',
           dir_base+'/02_01_03.txt',
           dir_base+'/02_01_04.txt',
           dir_base+'/02_01_05.txt',
           dir_base+'/02_01_06.txt',
           dir_base+'/02_01_07.txt',
           dir_base+'/02_01_08.txt',
           dir_base+'/02_01_09.txt',
           dir_base+'/02_01_10.txt',
           dir_base+'/02_01_11.txt',
           dir_base+'/02_01_12.txt',
           dir_base+'/02_01_13.txt',
           dir_base+'/02_01_14.txt']

fnames_heii=[dir_base+'/02_02_01.txt',
           dir_base+'/02_02_02.txt',
           dir_base+'/02_02_03.txt',
           dir_base+'/02_02_04.txt',
           dir_base+'/02_02_05.txt',
           dir_base+'/02_02_06.txt',
           dir_base+'/02_02_07.txt',
           dir_base+'/02_02_08.txt',
           dir_base+'/02_02_09.txt',
           dir_base+'/02_02_10.txt',
           dir_base+'/02_02_11.txt',
           dir_base+'/02_02_12.txt',
           dir_base+'/02_02_13.txt',
           dir_base+'/02_02_14.txt']


fnames_ciii=[dir_base+'/06_03_01.txt',
             dir_base+'/06_03_02.txt',
             dir_base+'/06_03_03.txt',
             dir_base+'/06_03_04.txt',
             dir_base+'/06_03_05.txt',
             dir_base+'/06_03_06.txt',
             dir_base+'/06_03_07.txt',
             dir_base+'/06_03_08.txt',
             dir_base+'/06_03_09.txt',
             dir_base+'/06_03_10.txt',
             dir_base+'/06_03_11.txt',
             dir_base+'/06_03_12.txt',
             dir_base+'/06_03_13.txt',
             dir_base+'/06_03_14.txt',
             dir_base+'/06_03_15.txt',
             dir_base+'/06_03_16.txt']


fnames_civ=[dir_base+'/06_04_01.txt',
             dir_base+'/06_04_02.txt',
             dir_base+'/06_04_03.txt',
             dir_base+'/06_04_04.txt',
             dir_base+'/06_04_05.txt',
             dir_base+'/06_04_06.txt',
             dir_base+'/06_04_07.txt',
             dir_base+'/06_04_08.txt',
             dir_base+'/06_04_09.txt',
             dir_base+'/06_04_10.txt',
             dir_base+'/06_04_11.txt',
             dir_base+'/06_04_12.txt',
             dir_base+'/06_04_13.txt',
             dir_base+'/06_04_14.txt',
             dir_base+'/06_04_15.txt',
             dir_base+'/06_04_16.txt']




fnames_niii=[dir_base+'/07_03_01.txt',
             dir_base+'/07_03_02.txt',
             dir_base+'/07_03_03.txt',
             dir_base+'/07_03_04.txt',
             dir_base+'/07_03_05.txt',
             dir_base+'/07_03_06.txt',
             dir_base+'/07_03_07.txt',
             dir_base+'/07_03_08.txt',
             dir_base+'/07_03_09.txt',
             dir_base+'/07_03_10.txt',
             dir_base+'/07_03_11.txt',
             dir_base+'/07_03_12.txt',
             dir_base+'/07_03_13.txt',
             dir_base+'/07_03_14.txt',
             dir_base+'/07_03_15.txt',
             dir_base+'/07_03_16.txt']


fnames_niv=[dir_base+'/07_04_01.txt',
             dir_base+'/07_04_02.txt',
             dir_base+'/07_04_03.txt',
             dir_base+'/07_04_04.txt',
             dir_base+'/07_04_05.txt',
             dir_base+'/07_04_06.txt',
             dir_base+'/07_04_07.txt',
             dir_base+'/07_04_08.txt',
             dir_base+'/07_04_09.txt',
             dir_base+'/07_04_10.txt',
             dir_base+'/07_04_11.txt',
             dir_base+'/07_04_12.txt',
             dir_base+'/07_04_13.txt',
             dir_base+'/07_04_14.txt',
             dir_base+'/07_04_15.txt',
             dir_base+'/07_04_16.txt']  

main(fnames=fnames_hei, windx=1)
main(fnames=fnames_heii, windx=2)
main(fnames=fnames_ciii, windx=3)
main(fnames=fnames_civ, windx=4)
main(fnames=fnames_niii, windx=5)
main(fnames=fnames_niv, windx=6)


sdum = input("Press [q] to exit.")
if sdum == 'q':
    exit()
