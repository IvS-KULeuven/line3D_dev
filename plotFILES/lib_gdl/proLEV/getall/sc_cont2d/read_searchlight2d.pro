pro read_searchlight2d, fname, ndxmax=ndxmax, ndzmax=ndzmax, $
                      xarr=xarr, zarr=zarr, int2d=int2d, inttheo2d=inttheo2d, $
                      mask2d=mask2d, opac2d=opac2d, epsmaxi_arr=epsmaxi_arr, $
                      nn_z=nn_z, nhelp=print_help
;+
; NAME:
;
;-
;
;-----------------------if help is needed-------------------------------
;
if(keyword_set(print_help)) then begin
   doc_library, 'read_searchlight2d'
   return
endif
;
if(n_elements(fname) eq 0) then begin
   print, 'error in read_searchlight2d: fname not specified'
   doc_library, 'read_searchlight2d'
   stop
endif
;
if(file_test(fname) ne 1) then begin
   print, 'error in read_searchlight2d: file does not exist'
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
      dset_id=h5d_open(group_id, 'z')
         zarr=h5d_read(dset_id)
      h5d_close, dset_id
   h5g_close, group_id
;
;----------------------angular grids------------------------------------
;
   group_id = h5g_open(file_id, 'angles')
      dset_id=h5d_open(group_id, 'nodes_mu')
         nodes_mu=h5d_read(dset_id)
      h5d_close, dset_id
   h5g_close, group_id
   nn_z=nodes_mu(0)   
;
;-------------------------convergence behaviour-------------------------
;
   group_id = h5g_open(file_id, 'convergence_behaviour')
      att_id=h5a_open_name(group_id, 'itmaxi')
         itmaxi=h5a_read(att_id)
         itmaxi=itmaxi(0)
      h5a_close, att_id
      att_id=h5a_open_name(group_id, 'devmaxi')
         devmaxi=h5a_read(att_id)
         devmaxi=devmaxi(0)
      h5a_close, att_id
      att_id=h5a_open_name(group_id, 'nconvi')
         nconvi=h5a_read(att_id)
         nconvi=nconvi(0)
      h5a_close, att_id
      dset_id=h5d_open(group_id, 'epsmaxi_arr')
         epsmaxi_arr=h5d_read(dset_id)
         epsmaxi_arr=epsmaxi_arr(where(abs(epsmaxi_arr) gt 0.))
      h5d_close, dset_id
   h5g_close, group_id
;
;-------------------------2d solution-----------------------------------
;
   group_id = h5g_open(file_id, 'solution2d')
      dset_id=h5d_open(group_id, 'int2d')
         int2d=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'int2d_theo')
         inttheo2d=h5d_read(dset_id)
      h5d_close, dset_id
   h5g_close, group_id
;
;-------------------------3d model--------------------------------------
;
   group_id = h5g_open(file_id, 'model2d')
      dset_id=h5d_open(group_id, 'mask2d')
         mask2d=h5d_read(dset_id)
         mask2d=floor(mask2d)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'opac2d')
         opac2d=h5d_read(dset_id)
      h5d_close, dset_id
   h5g_close, group_id
;
h5f_close, file_id
;
;
end
