import numpy as np
import h5py
import lev_misc
#
def benchmark11(fname='./benchmark11.h5'):


   file_id = h5py.File(fname, 'r')

   group_id = file_id['options']
   opt_angint_method = group_id.attrs['opt_angint_method']
   opt_angint_method = opt_angint_method[0]

   group_id = file_id['dimensions']
   ndxmax = group_id.attrs['ndxmax']
   ndymax = group_id.attrs['ndymax']
   ndzmax = group_id.attrs['ndzmax']
   dim_omega = group_id.attrs['dim_omega']
   n1d_angdep = group_id.attrs['n1d_angdep']
   ndxmax = ndxmax[0]
   ndymax = ndymax[0]
   ndzmax = ndzmax[0]
   dim_omega = dim_omega[0]
   n1d_angdep = n1d_angdep[0]

   group_id = file_id['coordinates']
   x = np.array(group_id['x'])
   y = np.array(group_id['y'])
   z = np.array(group_id['z'])
   xcoord_angdep = np.array(group_id['xcoord_angdep'])
   ycoord_angdep = np.array(group_id['ycoord_angdep'])
   zcoord_angdep = np.array(group_id['zcoord_angdep'])
   
   group_id = file_id['angles']
   n_x = np.array(group_id['n_x'])
   n_y = np.array(group_id['n_y'])
   n_z = np.array(group_id['n_z'])
   weight_omega = np.array(group_id['weight_omega'])
   
   group_id = file_id['solution3d']
   mint3d_sc = np.array(group_id['mint3d_sc'])
   mint3d_fvm = np.array(group_id['mint3d_fvm'])
   mint3d_theo = np.array(group_id['mint3d_theo'])
   mask3d = np.array(group_id['mask3d'])
   
   group_id = file_id['solution_angdep']
   intsc_angdep = np.array(group_id['intsc_angdep'])
   intfvm_angdep = np.array(group_id['intfvm_angdep'])

   file_id.close()

   nodes_mu=np.zeros(dim_omega)
   nodes_phi=np.zeros(dim_omega)

   for i in range(0,dim_omega):
      fdum=lev_misc.get_angles_spc(n_x[i],n_y[i],n_z[i])
      nodes_mu[i] = fdum[0]
      nodes_phi[i] = fdum[1]   

   return opt_angint_method, ndxmax, ndymax, ndzmax, \
       dim_omega, n1d_angdep, x, y, z, \
       xcoord_angdep, ycoord_angdep, zcoord_angdep, \
       n_x, n_y, n_z, weight_omega, nodes_mu, nodes_phi, \
       mint3d_sc, mint3d_fvm, mint3d_theo, mask3d, \
       intsc_angdep, intfvm_angdep
