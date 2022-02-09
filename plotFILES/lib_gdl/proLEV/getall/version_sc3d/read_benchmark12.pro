pro read_benchmark12, dir=dir, xic1=xic1, kcont=kcont, $
                      eps_cont=eps_cont, rstar=rstar, $
                      ndxmax=ndxmax, ndymax=ndymax, ndzmax=ndzmax, nr=nr, $
                      xarr=xarr, yarr=yarr, zarr=zarr, rarr=rarr, $
                      itmaxc=itmaxc, epsmaxc_sc=epsmaxc_sc, epsmaxc_fvm=epsmaxc_fvm, $
                      t1d_jo=t1d_jo, opac1d_jo=opac1d_jo, mint1d_joray=mint1d_joray, mint1d_jomom=mint1d_jomom, $
                      fcont1d_joray=fcont1d_joray, fcont1d_jomom=fcont1d_jomom, $
                      mask3d=mask3d, t3d=t3d, opac3d=opac3d, mint3d_sc=mint3d_sc, mint3d_fvm=mint3d_fvm, $
                      fcontr3d_sc=fcontr3d_sc, fcontth3d_sc=fcontth3d_sc, fcontphi3d_sc=fcontphi3d_sc, $
                      fcontr3d_fvm=fcontr3d_fvm, fcontth3d_fvm=fcontth3d_fvm, fcontphi3d_fvm=fcontphi3d_fvm, $
                      help=print_help, version=version


;+
; NAME:
;	read_benchmark12
;
; PURPOSE:
;	This procedure reads benchmark12.h5 file from sc3d.eo output
;
; CALLING SEQUENCE:
;
;	read_benchmark12
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
;       dir:    Set this keyword to the directory, where benchmark12.h5 shall
;               be read
;
;       version: Set this keyword to version of h5-file
;
;       Set keyword to read in corresponding parameter/array
;       Available keywords are:
;
;       options:
;
;       input parameters / stellar parameters:
;          kcont, eps_cont
;
;       dimensions:
;          ndxmax, ndymax, ndzmax, nr
;
;       spatial/angular/frequency grids:
;          xarr, yarr, zarr, rarr
;
;       convergence behaviour:
;          itmaxc, epsmaxc_sc, epsmaxc_fvm
;
;       model and solution on central ray:
;         t1d_jo, opac1d_jo, mint1d_joray, mint1d_jomom, fcont1d_joray, fcont1d_jomom
;
;       model and solution on 3d grid:
;         mask3d, t3d, opac3d, mint3d_sc, mint3d_fvm, 
;         fcontr3d_sc, fcontth3d_sc, fcontphi3d_sc, 
;         fcontr3d_fvm, fcontth3d_fvm, fcontphi3d_fvm
;
; EXAMPLE:
;	read_benchmark12, dir='.', ndxmax=ndxmax, ndymax=ndymax, ndzmax=ndzmax, $
;                  xarr=xarr, yarr=yarr, zarr=zarr, mint3d_sc=mint3d_sc
;
;-
;
;-----------------------if help is needed-------------------------------
;
if(keyword_set(print_help)) then begin
   doc_library, 'read_benchmark12'
   return
endif
;
if(not keyword_set(dir)) then dir='.'
fname=dir+'/benchmark12.h5'
;
if(file_test(fname) ne 1) then begin
   print, 'error in read_benchmark12: file does not exist'
   return
endif
;
if(not keyword_set(version)) then version=0
;
;------------------read all information from hdf5-file------------------
;
file_id = h5f_open(fname)
;
   group_id = h5g_open(file_id, 'input_parameters')
      att_id=h5a_open_name(group_id, 'kcont')
         kcont=h5a_read(att_id)
         kcont=kcont(0)
      h5a_close, att_id
      att_id=h5a_open_name(group_id, 'eps_cont')
         eps_cont=h5a_read(att_id)
         eps_cont=eps_cont(0)
      h5a_close, att_id
      att_id=h5a_open_name(group_id, 'xic1')
         xic1=h5a_read(att_id)
         xic1=xic1(0)
      h5a_close, att_id
   h5g_close, group_id

   group_id = h5g_open(file_id, 'dimensions')
      att_id=h5a_open_name(group_id, 'ndxmax')
         ndxmax=h5a_read(att_id)
         ndxmax=ndxmax(0)
      h5a_close, att_id
      att_id=h5a_open_name(group_id, 'ndymax')
         ndymax=h5a_read(att_id)
         ndymax=ndymax(0)
      h5a_close, att_id
      att_id=h5a_open_name(group_id, 'ndzmax')
         ndzmax=h5a_read(att_id)
         ndzmax=ndzmax(0)
      h5a_close, att_id
      att_id=h5a_open_name(group_id, 'nr')
         nr=h5a_read(att_id)
         nr=nr(0)
      h5a_close, att_id
   h5g_close, group_id
;
   group_id = h5g_open(file_id, 'coordinates')
      dset_id=h5d_open(group_id, 'r')
         rarr=h5d_read(dset_id)
      h5d_close, dset_id
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
   group_id = h5g_open(file_id, 'convergence_behaviour')
      att_id=h5a_open_name(group_id, 'itmaxc')
         itmaxc=h5a_read(att_id)
         itmaxc=itmaxc(0)
      h5a_close, att_id
      dset_id=h5d_open(group_id, 'epsmaxc_sc')
         epsmaxc_sc=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'epsmaxc_fvm')
         epsmaxc_fvm=h5d_read(dset_id)
      h5d_close, dset_id
   h5g_close, group_id
;
   group_id = h5g_open(file_id, 'model_cr')
      dset_id=h5d_open(group_id, 't1d_jo')
         t1d_jo=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'opac1d_jo')
         opac1d_jo=h5d_read(dset_id)
      h5d_close, dset_id
   h5g_close, group_id
;
   group_id = h5g_open(file_id, 'model3d')
      dset_id=h5d_open(group_id, 'mask3d')
         mask3d=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 't3d')
         t3d=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'opac3d')
         opac3d=h5d_read(dset_id)
      h5d_close, dset_id
   h5g_close, group_id
;
   group_id = h5g_open(file_id, 'solution_cr')
      dset_id=h5d_open(group_id, 'mint_joray')
         mint1d_joray=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'mint_jomom')
         mint1d_jomom=h5d_read(dset_id)
      h5d_close, dset_id
      if(version ne 0) then begin
         dset_id=h5d_open(group_id, 'fcont_joray')
            fcont1d_joray=h5d_read(dset_id)
         h5d_close, dset_id
      endif
      if(version ne 0) then begin
         dset_id=h5d_open(group_id, 'fcont_jomom')
            fcont1d_jomom=h5d_read(dset_id)
         h5d_close, dset_id
      endif
   h5g_close, group_id
;
   group_id = h5g_open(file_id, 'solution3d')
      dset_id=h5d_open(group_id, 'mint3d_sc')
         mint3d_sc=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'mint3d_fvm')
         mint3d_fvm=h5d_read(dset_id)
      h5d_close, dset_id
      if(version ne 0) then begin
         dset_id=h5d_open(group_id, 'fcontr3d_sc')
            fcontr3d_sc=h5d_read(dset_id)
         h5d_close, dset_id
      endif
      if(version ne 0) then begin
         dset_id=h5d_open(group_id, 'fcontth3d_sc')
            fcontth3d_sc=h5d_read(dset_id)
         h5d_close, dset_id
      endif
      if(version ne 0) then begin
         dset_id=h5d_open(group_id, 'fcontphi3d_sc')
            fcontphi3d_sc=h5d_read(dset_id)
         h5d_close, dset_id
      endif
      if(version ne 0) then begin
         dset_id=h5d_open(group_id, 'fcontr3d_fvm')
            fcontr3d_fvm=h5d_read(dset_id)
         h5d_close, dset_id
      endif
      if(version ne 0) then begin
         dset_id=h5d_open(group_id, 'fcontth3d_fvm')
            fcontth3d_fvm=h5d_read(dset_id)
         h5d_close, dset_id
      endif
      if(version ne 0) then begin
         dset_id=h5d_open(group_id, 'fcontphi3d_fvm')
            fcontphi3d_fvm=h5d_read(dset_id)
         h5d_close, dset_id
      endif
       
   h5g_close, group_id
;
h5f_close, file_id
;
;
end
