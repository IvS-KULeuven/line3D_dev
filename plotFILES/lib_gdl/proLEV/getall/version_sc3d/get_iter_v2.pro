pro get_iter_v2, fname, itmaxc, itmaxl, devmaxc, devmaxl, epsmaxc_arr, epsmaxl_arr, $
             help=print_help
;
;+
; NAME:
;	get_iter2, itmaxc, itmaxl, devmaxc, devmaxl, epsmaxc_arr, epsmaxl_arr
;
; PURPOSE:
;	This procedure reads only iteration history from *.h5 file of sc3d.eo output
;
; CALLING SEQUENCE:
;
;	get_iter2, fname, itmaxc, itmaxl, devmaxc, devmaxl, epsmaxc_arr, epsmaxl_arr
;
; INPUTS:
;	fname:  file name, that will be read in.
;
;       all other inputs correspond to the variables from sc3d.f90
;	
; KEYWORD PARAMETERS:
;
;       help:   Set this keyword (flag) to show the documentation of the
;               procedure
;
; EXAMPLE:
;	get_iter2, fname, itmaxc, itmaxl, devmaxc, devmaxl, epsmaxc_arr, epsmaxl_arr
;
;-
;
;-----------------------if help is needed-------------------------------
;
if(keyword_set(print_help)) then begin
   doc_library, 'get_iter_v2'
   return
endif
;
if(n_elements(fname) eq 0) then begin
   print, 'error in get_all: fname not specified'
   doc_library, 'get_iter_v2'
   stop
endif
;
;------------------read all information from hdf5-file------------------
;
file_id = h5f_open(fname)
;
;----------------------convergence behaviour----------------------------
;
   group_id = h5g_open(file_id, 'convergence_behaviour')
      att_id=h5a_open_name(group_id, 'itmaxc')
         itmaxc=h5a_read(att_id)
         itmaxc=itmaxc(0)
      h5a_close, att_id
      att_id=h5a_open_name(group_id, 'itmaxl')
         itmaxl=h5a_read(att_id)
         itmaxl=itmaxl(0)
      h5a_close, att_id
      att_id=h5a_open_name(group_id, 'devmaxc')
         devmaxc=h5a_read(att_id)
         devmaxc=devmaxc(0)
      h5a_close, att_id
      att_id=h5a_open_name(group_id, 'devmaxl')
         devmaxl=h5a_read(att_id)
         devmaxl=devmaxl(0)
      h5a_close, att_id
      dset_id=h5d_open(group_id, 'epsmaxc_arr')
         epsmaxc_arr=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'epsmaxl_arr')
         epsmaxl_arr=h5d_read(dset_id)
      h5d_close, dset_id
   h5g_close, group_id
;
h5f_close, file_id
;
;
end
