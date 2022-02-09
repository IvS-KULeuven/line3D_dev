pro read_searchlight3d, fname, ndxmax=ndxmax, ndymax=ndymax, ndzmax=ndzmax, $
                      xarr=xarr, yarr=yarr, zarr=zarr, int3d=int3d, inttheo3d=inttheo3d, $
                      mask3d=mask3d, opac3d=opac3d, nn_x=nn_x, nn_y=nn_y, nn_z=nn_z, nhelp=print_help
;+
; NAME:
;
;-
;
;-----------------------if help is needed-------------------------------
;
if(keyword_set(print_help)) then begin
   doc_library, 'read_searchlight3d'
   return
endif
;
if(n_elements(fname) eq 0) then begin
   print, 'error in read_searchlight3d: fname not specified'
   doc_library, 'read_searchlight3d'
   stop
endif
;
if(file_test(fname) ne 1) then begin
   print, 'error in read_searchlight3d: file does not exist'
   return
endif
;
;------------------read all information from hdf5-file------------------
;
file_id = h5f_open(fname)
;
;-----------------------dimensions--------------------------------------
;
   group_id = h5g_open(file_id, 'dimensions')
      att_id=h5a_open_name(group_id, 'ndxmax')
         ndxmax=h5a_read(att_id)
         ndxmax=floor(ndxmax(0))
         h5a_close, att_id
      att_id=h5a_open_name(group_id, 'ndymax')
         ndymax=h5a_read(att_id)
         ndymax=floor(ndymax(0))
         h5a_close, att_id
      att_id=h5a_open_name(group_id, 'ndzmax')
         ndzmax=h5a_read(att_id)
         ndzmax=floor(ndzmax(0))
      h5a_close, att_id
   h5g_close, group_id

;
;----------------------spatial coordinates------------------------------
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
;----------------------angular grids------------------------------------
;
   group_id = h5g_open(file_id, 'angles')
      dset_id=h5d_open(group_id, 'n_x')
         nn_x=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'n_y')
         nn_y=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'n_z')
         nn_z=h5d_read(dset_id)
      h5d_close, dset_id
   h5g_close, group_id
   nn_x=nn_x(0)
   nn_y=nn_y(0)
   nn_z=nn_z(0)

;
;-------------------------3d solution-----------------------------------
;
   group_id = h5g_open(file_id, 'solution3d')
      dset_id=h5d_open(group_id, 'int3d')
         int3d=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'int3d_theo')
         inttheo3d=h5d_read(dset_id)
      h5d_close, dset_id
   h5g_close, group_id
;
;-------------------------3d model--------------------------------------
;
   group_id = h5g_open(file_id, 'model3d')
      dset_id=h5d_open(group_id, 'mask3d')
         mask3d=h5d_read(dset_id)
         mask3d=floor(mask3d)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'opac3d')
         opac3d=h5d_read(dset_id)
      h5d_close, dset_id
   h5g_close, group_id
;
h5f_close, file_id
;
;
end
