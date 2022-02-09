import numpy as np
import matplotlib.cm as cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import copy
#
def get_clevels(clim=[0.,1.],nticks=5):
   cmin=clim[0]
   cmax=clim[1]
   if cmax == cmin:
      cmin = cmin - 0.1*cmin
      cmax = cmax + 0.1*cmax

   dcol=(cmax-cmin)/99.
   clevels=np.arange(cmin,cmax+dcol,dcol)
   ticks=np.linspace(cmin,cmax,nticks)
   return clevels, ticks

def get_logvals(x):   
   indx = np.where(x > 0.)
   xmin = np.min(x[indx])
   logx = np.copy(x)
   indx = np.where(logx <= 0.)
   logx[indx] = xmin
   logx = np.log10(logx)
   return logx

def get_cmap(cmap, upper=None, lower=None):

   if cmap == 'blues_lev':
      cdict = {'red':   [[0.0,  0.0, 0.0],
                         [0.4,  0.0, 0.0],
                         [1.0,  0.9, 0.9]],
               'green': [[0.0,  0.0, 0.0],
                         [0.4,  0.0, 0.0],
                         [1.0,  0.9, 0.9]],
               'blue':  [[0.0,  0.0, 0.0],
                         [0.5,  0.7, 0.7],
                         [1.0,  1.0, 1.0]]}

      cmap_lev = LinearSegmentedColormap('cmap_lev', segmentdata=cdict, N=256)
   
   elif cmap == 'jet':
      cmap_lev = copy.copy(cm.get_cmap('jet'))
   elif cmap == 'Blues':
      cmap_lev = copy.copy(cm.get_cmap('Blues'))   
   elif cmap == 'viridis':
      cmap_lev = copy.copy(cm.get_cmap('viridis'))
   elif cmap == 'seismic':
      cmap_lev = copy.copy(cm.get_cmap('seismic'))         
   elif cmap == 'Greys':
      cmap_lev = copy.copy(cm.get_cmap('Greys'))
   elif cmap == 'Greys_r':
      cmap_lev = copy.copy(cm.get_cmap('Greys_r'))
   else:
      print('error in get_cmap: cmap-string not included in get_cmap')
      cmap_lev=None


   if upper == 'black':
      cmap_lev.set_over('black')
   elif upper == 'grey':
      cmap_lev.set_over('grey')      

   if lower == 'black':
      cmap_lev.set_under('black')
   elif lower == 'grey':
      cmap_lev.set_under('grey')
      
      
   return cmap_lev
   
