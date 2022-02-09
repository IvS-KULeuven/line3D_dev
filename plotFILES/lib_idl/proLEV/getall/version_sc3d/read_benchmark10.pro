pro read_benchmark10, dir=dir, $
                      n_x=n_x, n_y=n_y, n_z=n_z, $
                      ndxmax=ndxmax, ndymax=ndymax, ndzmax=ndzmax, $
                      xarr=xarr, yarr=yarr, zarr=zarr, $
                      int3d_sc=int3d_sc, int3d_fvm=int3d_fvm, $
                      help=print_help


;+
; NAME:
;	read_benchmark10
;
; PURPOSE:
;	This procedure reads benchmark10.h5 file from sc3d.eo output
;
; CALLING SEQUENCE:
;
;	read_benchmark10
;
; INPUTS:
;
;       all inputs correspond to the variables from sc3d.f90
;	
; KEYWORD PARAMETERS:
;
;       help:   Set this keyword (flag) to show the documentation of the
;               procedure
;
;       dir:    Set this keyword to the directory, where benchmark13.h5 shall
;               be read
;
;       Set keyword to read in corresponding parameter/array
;       Available keywords are:
;
;       options:
;
;       input parameters / stellar parameters:
;
;       dimensions:
;          ndxmax, ndymax, ndzmax
;
;       spatial/angular/frequency grids:
;          n_x, n_y, n_z, xarr, yarr, zarr
;
;       convergence behaviour:
;
;       model and solution on central ray:
;
;       model and solution on 3d grid:
;         int3d_sc, int3d_fvm
;
; EXAMPLE:
;	read_benchmark10, dir='.', ndxmax=ndxmax, ndymax=ndymax, ndzmax=ndzmax, $
;                  xarr=xarr, yarr=yarr, zarr=zarr, int3d_sc=int3d_sc
;
;-
;
;-----------------------if help is needed-------------------------------
;
if(keyword_set(print_help)) then begin
   doc_library, 'read_benchmark10'
   return
endif
;
if(not keyword_set(dir)) then dir='.'
fname=dir+'/benchmark10.h5'
;
if(file_test(fname) ne 1) then begin
   print, 'error in read_benchmark10: file does not exist'
   return
endif
;
;
;------------------read all information from hdf5-file------------------
;
file_id = h5f_open(fname)
;
   group_id = H5G_OPEN(file_id, 'direction')
      att_id=H5A_OPEN_NAME(group_id, 'n_x')
         n_x=H5A_READ(att_id)
         n_x=n_x(0)
      h5a_close, att_id
      att_id=H5A_OPEN_NAME(group_id, 'n_y')
         n_y=H5A_READ(att_id)
         n_y=n_y(0)
      h5a_close, att_id
      att_id=H5A_OPEN_NAME(group_id, 'n_z')
         n_z=H5A_READ(att_id)
         n_z=n_z(0)
      h5a_close, att_id
   h5g_close, group_id

   group_id = H5G_OPEN(file_id, 'dimensions')
      att_id=H5A_OPEN_NAME(group_id, 'ndxmax')
         ndxmax=H5A_READ(att_id)
         ndxmax=ndxmax(0)
      h5a_close, att_id
      att_id=H5A_OPEN_NAME(group_id, 'ndymax')
         ndymax=H5A_READ(att_id)
         ndymax=ndymax(0)
      h5a_close, att_id
      att_id=H5A_OPEN_NAME(group_id, 'ndzmax')
         ndzmax=H5A_READ(att_id)
         ndzmax=ndzmax(0)
      h5a_close, att_id
   h5g_close, group_id
;
   group_id = h5g_open(file_id, 'coordinates')
      dset_id=h5d_open(group_id, 'x')
         xarr=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'y')
         yarr=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'z')
         zarr=h5d_read(dset_id)
      h5d_close, dset_id
   h5g_close, group_id
;
   group_id = h5g_open(file_id, 'searchlight3d')
      dset_id=h5d_open(group_id, 'int3d_sc')
         int3d_sc=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'int3d_fvm')
         int3d_fvm=h5d_read(dset_id)
      h5d_close, dset_id
   h5g_close, group_id
;
h5f_close, file_id
;
;
end
