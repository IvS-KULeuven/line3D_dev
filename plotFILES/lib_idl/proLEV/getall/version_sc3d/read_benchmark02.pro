pro read_benchmark02, dir=dir, xic1=xic1, mu=mu, $
                      ndxmax=ndxmax, ndzmax=ndzmax, $
                      xarr=xarr, zarr=zarr, $
                      int2d_sc=int2d_sc, int2d_fvm=int2d_fvm, $
                      help=print_help


;+
; NAME:
;	read_benchmark02
;
; PURPOSE:
;	This procedure reads benchmark02.h5 file from sc3d.eo output
;
; CALLING SEQUENCE:
;
;	read_benchmark02
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
;          xic1, mu
;
;       dimensions:
;          ndxmax, ndzmax
;
;       spatial/angular/frequency grids:
;          xarr, zarr
;
;       convergence behaviour:
;
;       model and solution on central ray:
;         
;       model and solution on 2d grid:
;         int2d_sc, int2d_fvm
;
;       intensity for certain directions:
;
; EXAMPLE:
;	read_benchmark02, dir='.', ndxmax=ndxmax, ndzmax=ndzmax, $
;                  xarr=xarr, zarr=zarr, int2d_sc=int2d_sc
;
;-
;
;-----------------------if help is needed-------------------------------
;
if(keyword_set(print_help)) then begin
   doc_library, 'read_benchmark02'
   return
endif
;
if(not keyword_set(dir)) then dir='.'
fname=dir+'/benchmark02.h5'
;
if(file_test(fname) ne 1) then begin
   print, 'error in read_benchmark02: file does not exist'
   return
endif
;
;
;------------------read all information from hdf5-file------------------
;
file_id = h5f_open(fname)
;
   group_id = H5G_OPEN(file_id, 'input_parameters')
      att_id=H5A_OPEN_NAME(group_id, 'xic1')
         xic1=H5A_READ(att_id)
         xic1=xic1(0)
      h5a_close, att_id
   h5g_close, group_id

   group_id = H5G_OPEN(file_id, 'direction')
      att_id=H5A_OPEN_NAME(group_id, 'mu')
         mu=H5A_READ(att_id)
         mu=mu(0)
      h5a_close, att_id
   h5g_close, group_id

   group_id = H5G_OPEN(file_id, 'dimensions')
      att_id=H5A_OPEN_NAME(group_id, 'ndxmax')
         ndxmax=H5A_READ(att_id)
         ndxmax=ndxmax(0)
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
      dset_id=h5d_open(group_id, 'z')
         zarr=h5d_read(dset_id)
      h5d_close, dset_id
   h5g_close, group_id
;
   group_id = h5g_open(file_id, 'searchlight2d')
      dset_id=h5d_open(group_id, 'int2d_sc')
         int2d_sc=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'int2d_fvm')
         int2d_fvm=h5d_read(dset_id)
      h5d_close, dset_id
   h5g_close, group_id
;
h5f_close, file_id
;
;
end
