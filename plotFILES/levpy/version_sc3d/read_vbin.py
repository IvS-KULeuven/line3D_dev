import numpy as np
import h5py
#
def read_surfb_vbin(fname='spec_surface.h5'):
#
# NAME:
#	read_surfb_vbin
#
# PURPOSE:
#	This procedure reads surface-brightnes *.h5 file from spec_vbin.eo
#
# CALLING SEQUENCE:
#
#	read_surfb_vbin
#
# INPUTS:
#
#       all inputs correspond to the variables from spec_vbin.f90
#	
# KEYWORD PARAMETERS:
#
#       help:   Set this keyword (flag) to show the documentation of the
#               procedure
#
#       fname:  Set this keyword to the filename
#
#       Set keyword to read in corresponding parameter/array
#       Available keywords are:
#
#       dimensions:
#          npoints
#
#       coordinates:
#          points_xcoord, points_ycoord
#
#       parameters:
#          xic1, xic2, xobs, alpha, gamma
#
#       surface brigthness:
#          int2d, int2d_abs, int2d_emi
#
# EXAMPLE:
#
#
#
#
   try:
       f = open(fname)
   except IOError:
       print('error in read_surfb_vbin: file not found')
   finally:
       f.close()
#
#------------------read all information from hdf5-file------------------
#

   file_id = h5py.File(fname, 'r')

   group_id = file_id['input_parameters']
   xic1 = group_id.attrs['xic1']
   xic2 = group_id.attrs['xic2']
   alpha = group_id.attrs['alpha']
   gamma = group_id.attrs['gamma']
   xobs = group_id.attrs['xobs']      
   xic1 = xic1[0]
   xic2 = xic2[0]
   alpha = alpha[0]
   gamma = gamma[0]
   xobs = xobs[0]

   group_id = file_id['dimensions']
   npoints = group_id.attrs['npoints']
   npoints = npoints[0]
 

   group_id = file_id['coordinates']
   points_xcoord = np.array(group_id['points_xcoord'])
   points_ycoord = np.array(group_id['points_ycoord'])   
   
   group_id = file_id['surfaceb']
   int2d_tot = np.array(group_id['iem_surface'])
   int2d_abs = np.array(group_id['iabs_surface'])
   int2d_emi = np.array(group_id['iemi_surface'])   
   

   file_id.close()

   return npoints, points_xcoord, points_ycoord, xic1, xic2, xobs, alpha, gamma, int2d_tot, int2d_abs, int2d_emi


#################################################################

def get_triangles (fname):

#+
# NAME:
#       getTRIANGLES
#
# PURPOSE:
#       This procedure reads in triangles
#
# CALLING SEQUENCE:
#
#
# INPUTS:
#
# KEYWORD PARAMETERS:
#      
# OUTPUTS:
#
# EXAMPLE:
#- 
#
   try:
       f = open(fname)
   except:
       print('error in get_triangles: file not found')
       exit()
#
#------------------------READ IN DATA FROM FILES------------------------
#
#dirty, but working
   f = open(fname, 'r')
   header = f.readline()
   npoints = 0
   for line in f.readlines():
      npoints=npoints+1
   f.close()

   x1 = np.zeros(npoints)
   y1 = np.zeros(npoints)
   x2 = np.zeros(npoints)
   y2 = np.zeros(npoints)
   x3 = np.zeros(npoints)
   y3 = np.zeros(npoints)
   f = open(fname, 'r')
   header = f.readline()
   i=0
   for line in f.readlines():
      line = line.strip()
      columns = line.split()
      x1[i] = float(columns[3])
      y1[i] = float(columns[4])
      x2[i] = float(columns[5])
      y2[i] = float(columns[6])
      x3[i] = float(columns[7])
      y3[i] = float(columns[8])
      i=i+1

   return x1, y1, x2, y2, x3, y3


#################################################################

def get_fluxem (fname):

#+
# NAME:
#       get_fluxem
#
# PURPOSE:
#       This procedure reads in emergent fluxes
#
# CALLING SEQUENCE:
#
#
# INPUTS:
#
# KEYWORD PARAMETERS:
#      
# OUTPUTS:
#
# EXAMPLE:
#- 
#
   try:
       f = open(fname)
   except:
       print('error in get_fluxem: file not found', fname)
       exit()
#
#------------------------READ IN DATA FROM FILES------------------------
#
#dirty, but working
   f = open(fname, 'r')
   
   nxobs = f.readline()
   nxobs = int(nxobs)

   line = f.readline()
   columns = line.split()
   alpha = float(columns[1])
   
   line = f.readline()
   columns = line.split()
   gamma = float(columns[1])   
   
   header = f.readline()

   xobs = np.zeros(nxobs)
   ftot = np.zeros(nxobs)
   fcont = np.zeros(nxobs)

   i=0
   for line in f.readlines():
      line = line.strip()
      columns = line.split()
      xobs[i] = float(columns[0])
      ftot[i] = float(columns[1])
      fcont[i] = float(columns[2])
      i=i+1

   return nxobs, alpha, gamma, xobs, ftot, fcont


#################################################################

def get_fluxem_debug (fname):

#+
# NAME:
#       get_fluxem_debug
#
# PURPOSE:
#       This procedure reads in emergent fluxes
#
# CALLING SEQUENCE:
#
#
# INPUTS:
#
# KEYWORD PARAMETERS:
#      
# OUTPUTS:
#
# EXAMPLE:
#- 
#
   try:
       f = open(fname)
   except:
       print('error in get_fluxem_debug: file not found')
       exit()
#
#------------------------READ IN DATA FROM FILES------------------------
#
#dirty, but working
   f = open(fname, 'r')
   
   nxobs = f.readline()
   nxobs = int(nxobs)

   line = f.readline()
   columns = line.split()
   alpha = float(columns[1])
   
   line = f.readline()
   columns = line.split()
   gamma = float(columns[1])   
   
   header = f.readline()

   xobs = np.zeros(nxobs)
   ftot = np.zeros(nxobs)
   fcont = np.zeros(nxobs)
   femi = np.zeros(nxobs)
   fabs = np.zeros(nxobs)      

   i=0
   for line in f.readlines():
      line = line.strip()
      columns = line.split()
      xobs[i] = float(columns[0])
      ftot[i] = float(columns[1])
      fcont[i] = float(columns[2])
      femi[i] = float(columns[3])      
      fabs[i] = float(columns[4])
      i=i+1

   return nxobs, alpha, gamma, xobs, ftot, fcont, femi, fabs
