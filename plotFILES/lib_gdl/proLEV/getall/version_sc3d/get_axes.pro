pro get_axes, fname, x, y, z, help=print_help
;
;+
; NAME:
;	get_axes, fname, x, y, z
;
; PURPOSE:
;	This procedure reads only x, y and z-axes from *.h5 file of sc3d.eo output
;
; CALLING SEQUENCE:
;
;	get_axes, fname, x, y, z
;
; INPUTS:
;	fname:  file name, that will be read in.
;	
; KEYWORD PARAMETERS:
;
;       help:   Set this keyword (flag) to show the documentation of the
;               procedure
;
; EXAMPLE:
;	get_axes, fname, x, y, z
;
;-
;
;-----------------------if help is needed-------------------------------
;
if(keyword_set(print_help)) then begin
   doc_library, 'get_axes'
   return
endif
;
if(n_elements(fname) eq 0) then begin
   print, 'error in get_axse: fname not specified'
   doc_library, 'get_axes'
   stop
endif
;
;------------------read all information from hdf5-file------------------
;
file_id = h5f_open(fname)
;
   group_id = h5g_open(file_id, 'coordinates')
      dset_id=h5d_open(group_id, 'x')
         x=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'y')
         y=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'z')
         z=h5d_read(dset_id)
      h5d_close, dset_id

   h5g_close, group_id
;
h5f_close, file_id
;
;
end
