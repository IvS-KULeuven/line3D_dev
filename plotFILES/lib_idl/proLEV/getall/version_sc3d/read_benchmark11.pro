pro read_benchmark11, dir=dir, opt_angint_method=opt_angint_method, $
                      ndxmax=ndxmax, ndymax=ndymax, ndzmax=ndzmax, dim_omega=dim_omega, $
                      n1d_angdep=n1d_angdep, xarr=xarr, yarr=yarr, zarr=zarr, $
                      xcoord_angdep=xcoord_angdep, ycoord_angdep=ycoord_angdep, zcoord_angdep=zcoord_angdep, $
                      n_x=n_x, n_y=n_y, n_z=n_z, weight_omega=weight_omega, $
                      mask3d=mask3d, mint3d_sc=mint3d_sc, mint3d_fvm=mint3d_fvm, mint3d_theo=mint3d_theo, $
                      fcontr3d_theo=fcontr3d_theo, fcontr3d_sc=fcontr3d_sc, fcontr3d_fvm=fcontr3d_fvm, $
                      fcontth3d_theo=fcontth3d_theo, fcontth3d_sc=fcontth3d_sc, fcontth3d_fvm=fcontth3d_fvm, $
                      fcontphi3d_theo=fcontphi3d_theo, fcontphi3d_sc=fcontphi3d_sc, fcontphi3d_fvm=fcontphi3d_fvm, $
                      intsc_angdep=intsc_angdep, intfvm_angdep=intfvm_angdep, $
                      help=print_help


;+
; NAME:
;	read_benchmark11
;
; PURPOSE:
;	This procedure reads benchmark11.h5 file from sc3d.eo output
;
; CALLING SEQUENCE:
;
;	read_benchmark11
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
;          opt_angint_method
;
;       input parameters / stellar parameters:
;
;       dimensions:
;          ndxmax, ndymax, ndzmax, dim_omega, n1d_angdep
;
;       spatial/angular/frequency grids:
;          xarr, yarr, zarr, xcoord_angdep, ycoord_angdep, zcoord_angdep,
;          n_x, n_y, n_z, weight_omega
;
;       convergence behaviour:
;
;       model and solution on central ray:
;         
;       model and solution on 3d grid:
;         mask3d, mint3d_sc, mint3d_fvm, mint3d_theo
;
;       intensity for certain directions:
;          intsc_angdep, intfvm_angdep
;
; EXAMPLE:
;	read_benchmark11, dir='.', ndxmax=ndxmax, ndymax=ndymax, ndzmax=ndzmax, $
;                  xarr=xarr, yarr=yarr, zarr=zarr, mint3d_sc=mint3d_sc, mint3d_theo=mint3d_theo
;
;-
;
;-----------------------if help is needed-------------------------------
;
if(keyword_set(print_help)) then begin
   doc_library, 'read_benchmark11'
   return
endif
;
if(not keyword_set(dir)) then dir='.'
fname=dir+'/benchmark11.h5'
;
if(file_test(fname) ne 1) then begin
   print, 'error in read_benchmark11: file does not exist'
   return
endif
;
;
;------------------read all information from hdf5-file------------------
;
   file_id = h5f_open(fname)
;
      group_id = h5g_open(file_id, 'options')
         att_id=h5a_open_name(group_id, 'opt_angint_method')
            opt_angint_method=h5a_read(att_id)
            opt_angint_method=opt_angint_method(0)
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
         att_id=h5a_open_name(group_id, 'dim_omega')
            dim_omega=h5a_read(att_id)
            dim_omega=dim_omega(0)
         h5a_close, att_id
         att_id=h5a_open_name(group_id, 'n1d_angdep')
            n1d_angdep=h5a_read(att_id)
            n1d_angdep=n1d_angdep(0)
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
         dset_id=h5d_open(group_id, 'xcoord_angdep')
            xcoord_angdep=h5d_read(dset_id)
         h5d_close, dset_id
         dset_id=h5d_open(group_id, 'ycoord_angdep')
            ycoord_angdep=h5d_read(dset_id)
         h5d_close, dset_id
         dset_id=h5d_open(group_id, 'zcoord_angdep')
            zcoord_angdep=h5d_read(dset_id)
         h5d_close, dset_id
      h5g_close, group_id
;
      group_id = h5g_open(file_id, 'angles')
         dset_id=h5d_open(group_id, 'n_x')
            n_x=h5d_read(dset_id)
         h5d_close, dset_id
         dset_id=h5d_open(group_id, 'n_y')
            n_y=h5d_read(dset_id)
         h5d_close, dset_id
         dset_id=h5d_open(group_id, 'n_z')
            n_z=h5d_read(dset_id)
         h5d_close, dset_id
         dset_id=h5d_open(group_id, 'weight_omega')
            weight_omega=h5d_read(dset_id)
         h5d_close, dset_id
      h5g_close, group_id
   ;
      group_id = h5g_open(file_id, 'solution3d')
         dset_id=h5d_open(group_id, 'mint3d_sc')
            mint3d_sc=h5d_read(dset_id)
         h5d_close, dset_id
         dset_id=h5d_open(group_id, 'mint3d_fvm')
            mint3d_fvm=h5d_read(dset_id)
         h5d_close, dset_id
         dset_id=h5d_open(group_id, 'mint3d_theo')
            mint3d_theo=h5d_read(dset_id)
         h5d_close, dset_id
         dset_id=h5d_open(group_id, 'fcontr3d_sc')
            fcontr3d_sc=h5d_read(dset_id)
         h5d_close, dset_id
         dset_id=h5d_open(group_id, 'fcontr3d_fvm')
            fcontr3d_fvm=h5d_read(dset_id)
         h5d_close, dset_id
         dset_id=h5d_open(group_id, 'fcontr3d_theo')
            fcontr3d_theo=h5d_read(dset_id)
         h5d_close, dset_id
         dset_id=h5d_open(group_id, 'fcontth3d_sc')
            fcontth3d_sc=h5d_read(dset_id)
         h5d_close, dset_id
         dset_id=h5d_open(group_id, 'fcontth3d_fvm')
            fcontth3d_fvm=h5d_read(dset_id)
         h5d_close, dset_id
         dset_id=h5d_open(group_id, 'fcontth3d_theo')
            fcontth3d_theo=h5d_read(dset_id)
         h5d_close, dset_id
         dset_id=h5d_open(group_id, 'fcontphi3d_sc')
            fcontphi3d_sc=h5d_read(dset_id)
         h5d_close, dset_id
         dset_id=h5d_open(group_id, 'fcontphi3d_fvm')
            fcontphi3d_fvm=h5d_read(dset_id)
         h5d_close, dset_id
         dset_id=h5d_open(group_id, 'fcontphi3d_theo')
            fcontphi3d_theo=h5d_read(dset_id)
         h5d_close, dset_id
         dset_id=h5d_open(group_id, 'mask3d')
            mask3d=h5d_read(dset_id)
         h5d_close, dset_id
      h5g_close, group_id
   ;
      group_id = h5g_open(file_id, 'solution_angdep')
         dset_id=h5d_open(group_id, 'intsc_angdep')
            intsc_angdep=h5d_read(dset_id)
         h5d_close, dset_id
         dset_id=h5d_open(group_id, 'intfvm_angdep')
            intfvm_angdep=h5d_read(dset_id)
         h5d_close, dset_id
      h5g_close, group_id
   ;
   h5f_close, file_id
;
;
end
