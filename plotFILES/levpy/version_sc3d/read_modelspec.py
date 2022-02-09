import numpy as np
import h5py
#
def read_modelspec3d_spc(fname):


   file_id = h5py.File(fname, 'r')

   group_id = file_id['dimensions']
   nr = group_id.attrs['nr']
   nr = nr[0]
   ntheta = group_id.attrs['ntheta']
   ntheta = ntheta[0]
   nphi = group_id.attrs['nphi']
   nphi = nphi[0]


   group_id = file_id['input_parameters']
   teff = group_id.attrs['teff']
   teff = teff[0]   
   logg = group_id.attrs['logg']
   logg = logg[0]   
   yhe = group_id.attrs['yhe']
   yhe = yhe[0]   
   trad = group_id.attrs['trad']
   trad = trad[0]   
   xnue0 = group_id.attrs['xnue0']
   xnue0 = xnue0[0]   
   rstar = group_id.attrs['rstar']
   rstar = rstar[0]   
   lstar = group_id.attrs['lstar']
   lstar = lstar[0]   
   vth_fiducial = group_id.attrs['vth_fiducial']
   vth_fiducial = vth_fiducial[0]   
   vmicro = group_id.attrs['vmicro']
   vmicro = vmicro[0]   
   vmax = group_id.attrs['vmax']
   vmax = vmax[0]   
   vrot = group_id.attrs['vrot']
   vrot = vrot[0]   
   vmin = group_id.attrs['vmin']
   vmin = vmin[0]   
   na = group_id.attrs['na']
   na = na[0]


   group_id = file_id['bcondition']
   xic1 = group_id.attrs['xic1']
   xic1 = xic1[0]   
   
   group_id = file_id['coordinates']
   radius = np.array(group_id['r'])
   theta = np.array(group_id['theta'])
   phi = np.array(group_id['phi'])

   group_id = file_id['model3d']
   t3d = np.array(group_id['t3d'])
   opac3d = np.array(group_id['opac3d'])
   opalbar3d = np.array(group_id['opalbar3d'])
   velx3d = np.array(group_id['velx3d'])   
   vely3d = np.array(group_id['vely3d'])   
   velz3d = np.array(group_id['velz3d'])

   group_id = file_id['solution3d']
   sline3d = np.array(group_id['sline3d'])
   scont3d = np.array(group_id['scont3d'])

   file_id.close()

   return teff, trad, xnue0, xic1, rstar, vth_fiducial, vmicro, vmax, logg, na, \
           nr, ntheta, nphi, radius, theta, phi, t3d, opac3d, opalbar3d, velx3d, vely3d, velz3d, \
           sline3d, scont3d   

########################################################################


def read_modelspec3d_vbin(fname,version='v00'):
#
#+
# NAME:
#	read_modelspec3d_vbin(fname)
#
# PURPOSE:
#	This procedure reads modspec.h5 file, that has been used as input for
#	spec_vbin.eo (spherical coordinates), binary version
#
# CALLING SEQUENCE:
#
#	cs1_nr, cs1_ntheta, cs1_nphi, cs1_radius, cs1_theta,
#          cs1_phi, cs1_velx3d, cs1_vely3d, cs1_velz3d, cs1_t3d,
#          cs1_opac3d, cs1_opalbar3d, cs1_scont3d, cs1_sline3d
#          cs2_nr, cs2_ntheta, cs2_nphi, cs2_radius, cs2_theta,
#          cs2_phi, cs2_velx3d, cs2_vely3d, cs2_velz3d, cs2_t3d,
#          cs2_opac3d, cs2_opalbar3d, cs2_scont3d, cs2_sline3d = read_modelspec3d_vbin(fname)
#
# INPUTS:
#	fname:  file name, that will be read in.
#
# OUTPUTS:
#       all outputs correspond to the variables from modelspec_vbin.eo   
#
#          cs1_nr, cs1_ntheta, cs1_nphi, cs1_radius, cs1_theta,
#          cs1_phi, cs1_velx3d, cs1_vely3d, cs1_velz3d, cs1_t3d,
#          cs1_opac3d, cs1_opalbar3d, cs1_scont3d, cs1_sline3d
#          cs2_nr, cs2_ntheta, cs2_nphi, cs2_radius, cs2_theta,
#          cs2_phi, cs2_velx3d, cs2_vely3d, cs2_velz3d, cs2_t3d,
#          cs2_opac3d, cs2_opalbar3d, cs2_scont3d, cs2_sline3d
#
#
   file_id = h5py.File(fname, 'r')


   group_id = file_id['input_parameters']
   x01 = group_id.attrs['x01']
   x01 = x01[0]
   y01 = group_id.attrs['y01']
   y01 = y01[0]
   z01 = group_id.attrs['z01']
   z01 = z01[0]
   vx01 = group_id.attrs['vx01']
   vx01 = vx01[0]
   vy01 = group_id.attrs['vy01']
   vy01 = vy01[0]
   vz01 = group_id.attrs['vz01']
   vz01 = vz01[0]

   x02 = group_id.attrs['x02']
   x02 = x02[0]
   y02 = group_id.attrs['y02']
   y02 = y02[0]
   z02 = group_id.attrs['z02']
   z02 = z02[0]
   vx02 = group_id.attrs['vx02']
   vx02 = vx02[0]
   vy02 = group_id.attrs['vy02']
   vy02 = vy02[0]
   vz02 = group_id.attrs['vz02']
   vz02 = vz02[0]         
   
   rstar1 = group_id.attrs['rstar1']
   rstar1 = rstar1[0]
   rstar2 = group_id.attrs['rstar2']
   rstar2 = rstar2[0]

   teff1 = group_id.attrs['teff1']
   teff1 = teff1[0]
   teff2 = group_id.attrs['teff2']
   teff2 = teff2[0]      

   logg1 = group_id.attrs['logg1']
   logg1 = logg1[0]
   logg2 = group_id.attrs['logg2']
   logg2 = logg2[0]

   lstar1 = group_id.attrs['lstar1']
   lstar1 = lstar1[0]
   lstar2 = group_id.attrs['lstar2']
   lstar2 = lstar2[0]

   yhe1 = group_id.attrs['yhe1']
   yhe1 = yhe1[0]
   yhe2 = group_id.attrs['yhe2']
   yhe2 = yhe2[0]
   
   vrot1 = group_id.attrs['vrot1']
   vrot1 = vrot1[0]
   vrot2 = group_id.attrs['vrot2']
   vrot2 = vrot2[0]
   
   xnue0 = group_id.attrs['xnue0']
   xnue0 = xnue0[0]
   
   vth_fiducial = group_id.attrs['vth_fiducial']
   vth_fiducial = vth_fiducial[0]

   if version == 'v00':
      vmicro = group_id.attrs['vmicro']
      vmicro = vmicro[0]
      vmicro1 = vmicro
      vmicro2 = vmicro
   if version == 'v01':
      vmicro1 = group_id.attrs['vmicro1']
      vmicro1 = vmicro1[0]
      vmicro2 = group_id.attrs['vmicro2']
      vmicro2 = vmicro2[0]      
      
   vmax = group_id.attrs['vmax']
   vmax = vmax[0]

   na = group_id.attrs['na']
   na = na[0]
   
   iline = group_id.attrs['iline']
   iline = iline[0]   

   unit_length = group_id.attrs['unit_length']
   unit_length = unit_length[0]   
#
#----------------------------star1-------------------------------------
#
   group1_id = file_id['star1']
   group2_id = group1_id['dimensions']
   cs1_nr = group2_id.attrs['nr']
   cs1_nr = cs1_nr[0]
   cs1_ntheta = group2_id.attrs['ntheta']
   cs1_ntheta = cs1_ntheta[0]   
   cs1_nphi = group2_id.attrs['nphi']
   cs1_nphi = cs1_nphi[0]
#
   group2_id = group1_id['coordinates']
   cs1_radius = np.array(group2_id['r'])
   cs1_theta = np.array(group2_id['theta'])
   cs1_phi = np.array(group2_id['phi'])
#
   group2_id = group1_id['model3d']
   cs1_velx3d = np.array(group2_id['velx3d'])
   cs1_vely3d = np.array(group2_id['vely3d'])
   cs1_velz3d = np.array(group2_id['velz3d'])   
   cs1_t3d = np.array(group2_id['t3d'])
   cs1_rho3d = np.array(group2_id['rho3d'])      
   cs1_opac3d = np.array(group2_id['opac3d'])   
   cs1_opalbar3d = np.array(group2_id['opalbar3d'])      
#
   group2_id = group1_id['solution3d']
   cs1_scont3d = np.array(group2_id['scont3d'])
   cs1_sline3d = np.array(group2_id['sline3d'])      
#
#----------------------------star2-------------------------------------
#
   group1_id = file_id['star2']
   group2_id = group1_id['dimensions']
   cs2_nr = group2_id.attrs['nr']
   cs2_nr = cs2_nr[0]
   cs2_ntheta = group2_id.attrs['ntheta']
   cs2_ntheta = cs2_ntheta[0]   
   cs2_nphi = group2_id.attrs['nphi']
   cs2_nphi = cs2_nphi[0]
#
   group2_id = group1_id['coordinates']
   cs2_radius = np.array(group2_id['r'])
   cs2_theta = np.array(group2_id['theta'])
   cs2_phi = np.array(group2_id['phi'])
#
   group2_id = group1_id['model3d']
   cs2_velx3d = np.array(group2_id['velx3d'])
   cs2_vely3d = np.array(group2_id['vely3d'])
   cs2_velz3d = np.array(group2_id['velz3d'])   
   cs2_t3d = np.array(group2_id['t3d'])
   cs2_rho3d = np.array(group2_id['rho3d'])      
   cs2_opac3d = np.array(group2_id['opac3d'])   
   cs2_opalbar3d = np.array(group2_id['opalbar3d'])      
#
   group2_id = group1_id['solution3d']
   cs2_scont3d = np.array(group2_id['scont3d'])
   cs2_sline3d = np.array(group2_id['sline3d'])   
#
   file_id.close()

   return iline, xnue0, vth_fiducial, vmax, unit_length, \
      x01, y01, z01, vx01, vy01, vz01, rstar1, teff1, logg1, lstar1, yhe1, vrot1, vmicro1, \
      x02, y02, z02, vx02, vy02, vz02, rstar2, teff2, logg2, lstar2, yhe2, vrot2, vmicro2, \
      cs1_nr,  cs1_ntheta, cs1_nphi, cs1_radius, cs1_theta, cs1_phi, \
      cs1_velx3d, cs1_vely3d, cs1_velz3d, cs1_t3d, cs1_rho3d, cs1_opac3d, cs1_opalbar3d, cs1_scont3d, cs1_sline3d, \
      cs2_nr,  cs2_ntheta, cs2_nphi, cs2_radius, cs2_theta, cs2_phi, \
      cs2_velx3d, cs2_vely3d, cs2_velz3d, cs2_t3d, cs2_rho3d, cs2_opac3d, cs2_opalbar3d, cs2_scont3d, cs2_sline3d      
