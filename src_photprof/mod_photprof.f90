!
!-----------------------------------------------------------------------
!----------------photospheric profile from external---------------------
!-----------------------------------------------------------------------
!
module photprof_ext
!
   use prog_type
!
   integer(i4b):: nxobs_ext
!
   character(len=10) :: fname1, fname2
   character(len=6) :: fnamep_iijjkk, fnamec_iijjkk
   character(len=25), parameter :: dirphotprofp_v0='phot_flux/dir_prof_v0'
   character(len=30), parameter :: dirphotprofp_herrero1='phot_flux/sym/dir_prof_herrero'
   character(len=30), parameter :: dirphotprofc_herrero1='phot_flux/sym/dir_cont_herrero'
   character(len=30), parameter :: dirphotprofp_herrero2='phot_flux/sym/dir_prof_herrero'
   character(len=30), parameter :: dirphotprofc_herrero2='phot_flux/sym/dir_cont_herrero'
   character(len=500) :: dirphotprofp_herrero, dirphotprofc_herrero
!dirphotprof = 'phot_flux/blend'     if photospheric profile is used including he-blend
!dirphotprof = 'phot_flux/sym'       if photospheric profile is symmetrized in order to
!                                  remove he-component


   character(len=18), parameter :: dirphotprofp_fastwind='phot_flux/fastwind'
!
   character(len=20), parameter :: dirphotprofp_coelho05='phot_flux/s_coelho05'
   character(len=26), parameter :: dirphotprofp_coelho14='phot_flux/s_coelho14_hrplc'
!nxobs_ext: number of frequency points of external photspheric profile
!xnue_ext: frequency points of external photospheric profile
!xic_ext: intensity for given frequency points of external photospheric profile
!
!note: two different profiles are read in, and interpolation is performed if needed
!
end module photprof_ext
