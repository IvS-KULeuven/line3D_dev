!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!--------MAIN PROGRAM TO CALCULATE ATMOSPHERIC STRUCTURE----------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!purpose: calculate atmospheric structure and store everything in
!         an external hdf5-file
!
!v0: including model 0-5:
!     0: spherical wind (1D)
!     1: ablation with snapshot from dylan kee (2D)
!     2: ablation with initial conditions from dylan kee (2D)
!     3: adm (2D)
!     4: adm with beta-vel-law in wind outflow region (2D)
!     5: mhd snapshot from asif ud-doula (3D)
!     6: pseudo 2d rotation: radial beta-velocity + azimuthal velocities
!        (with 1d density structure)
!     7: bjorkman&casinelli(1993) model
!     8: 2d VH-1 model
!     9: input files from Cis (ku leuven 1d ldi simulations), output to 1D smooth model
!    10: input files from Cis (ku leuven 1d ldi simulations), output to 3D random model
!    11: input files from Luka (ku leuven simulations), output to 1D model
!    12: input files from Florian (ku leuven simulations), output to 3D model (azimuthal symmetry)
!    13: input files from Florian (ku leuven 3d mhd simulations)
!   all files are stored as h5-files, which can easily be read in into
!   main radiative transfer program
!    14: analytical Be-star disc
!    15: 1d test (beta-velocity wind)
!    16: christi's mhd models
!    17: nicos 3d slab with map-projections
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
program main
!
use prog_type
use params_input, only: input_mod
!
implicit none
!
! ... local scalars
!
! ... local arrays
!
! ... local logicals
!
! ... local functions
!
!---------------------read input parameters-----------------------------
!
call read_input
!
!---------------------print input parameters----------------------------
!
call print_input
!
!-------------------calculating model atmosphere------------------------
!
write(*,*) '-------------------------calculate model atmosphere----------------------------'
write(*,*)
!
select case(input_mod)
!
   case(0)
      write(*,*) '------calculating standard model 0--------'
      write(*,*)
      call calc_mod1d
      call output_mod1d
!
   case(1)
      write(*,*) '------calculating DYLANs model 1----------'
      write(*,*)
      call calc_mod2d_abla
      call output_mod2d
!
   case(2)
      write(*,*) '------calculating DYLANs model 2----------'
      write(*,*)
      call calc_mod2d_ablb
      call output_mod2d
!
   case(3)
      write(*,*) '--------calculating ADM model-------------'
      write(*,*)
      call calc_mod2d_adma
      call output_mod2d
!
   case(4)
      write(*,*) '---calculating ADM model (beta-vel-law)---'
      write(*,*)
      call calc_mod2d_admb
      call output_mod2d
!
   case(5)
      write(*,*) '--------calculating 3D MHD model----------'
      write(*,*)
      call calc_mod3d_mhd
      call output_mod3d
!
   case(6)
      write(*,*) '---calculating pseudo 2d rotation model---'
      write(*,*)
      call calc_mod2d_rot
      call output_mod2d
!
   case(7)
      write(*,*) '---------calculating 2D BC model----------'
      write(*,*)
      call calc_mod2d_bc
      call output_mod2d
!
   case(8)
      write(*,*) '---------calculating 2D VH1 model---------'
      write(*,*)
      call calc_mod2d_vh1
      call output_mod2d
!
   case(9)
      write(*,*) '---calculating 1D LDI model (from Cis)----'
      write(*,*)
      call calc_mod1d_cis
      call output_mod1d
!
   case(10)
      write(*,*) '---calculating 3D LDI model (from Cis)----'
      write(*,*)
      call calc_mod3d_cis
      call output_mod3d
!
   case(11)
      write(*,*) '----calculating 1D model (from Luka)------'
      write(*,*)
      call calc_mod1d_luka
      call output_mod1d
!
   case(12)
      write(*,*) '--calculating 2D MHD model (from Florian)-'
      write(*,*)
      call calc_mod2d_florian
      call output_mod2d
!
   case(13)
      write(*,*) '--calculating 3D MHD model (from Florian)-'
      write(*,*)
      call calc_mod3d_florian
      call output_mod3d
!
   case(14)
      write(*,*) '---calculating analytical Be star disc----'
      write(*,*)
      call calc_mod2d_be
      call output_mod2d
!
   case(15)
      write(*,*) '---calculating analytical 1d test model---'
      write(*,*)
      call calc_mod1d_test
      call output_mod1d
!
   case(16)
      write(*,*) '--calculating 3D MHD model (from Christi)-'
      write(*,*)
      call calc_mod3d_christi
      call output_mod3d
!
   case(17)
      write(*,*) '---calculating 3D RHD model (from Nico)---'
      write(*,*)
      call calc_mod3d_nicowr3d
      call output_mod3d      
!
    case default
       stop 'error: unvalid input-model specified, check input_mod'
!
end select
!
!
end program main
!
!***********************************************************************
!***********************************************************************
!***********************************************************************
!
!---------------------SUBROUTINES AND FUNCTIONS-------------------------
!
!***********************************************************************
!***********************************************************************
!***********************************************************************
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine read_input
!
!-----------------------------------------------------------------------
!-------------------read input parameters-------------------------------
!-----------------------------------------------------------------------
!
use prog_type
use mod_directories, only: indat_file, model_dir
use params_input, only: trad, teff, tmin, xlogg, rstar, lstar, rmax, rlim, tmin, xmloss, vmin, &
                        vmax, beta, yhe, hei, input_mod, na, vmicro, rmin, opt_incl_line, &
                       vrot, eps_cont, kcont
use params_stellar, only: sr
use fund_const, only: rsu, pi, cgs_grav, xmsu
!
!!
implicit none
!
! ... local scalars
real(dp) :: xnue0, vth_fiducial
integer(i4b) :: input_mod_dim, spatial_grid1d, spatial_grid3d, opt_opac, opt_opal, opt_ltec, &
                opt_angint_method, opt_method, opt_alo_cont, opt_alo_line
logical :: opt_sol2d, opt_incl_cont, opt_start_cont, opt_start_line, opt_ng_cont, &
           opt_ng_line, opt_ait_cont, opt_ait_line, opt_incl_gdark, opt_incl_sdist
!
! ... local characters
character(len=100) :: output_file
!
! ... local functions
!
! ... namelist
namelist / input_options / model_dir, output_file, input_mod, input_mod_dim, spatial_grid1d, spatial_grid3d, &
                           opt_opac, opt_opal, opt_angint_method, opt_method, opt_sol2d, opt_ltec, opt_incl_cont, &
                           opt_start_cont, opt_ng_cont, opt_ait_cont, opt_incl_line, opt_start_line, &
                           opt_ng_line, opt_ait_line, opt_alo_cont, opt_alo_line, opt_incl_gdark, opt_incl_sdist
namelist / input_mod_1d / teff, trad, xlogg, rstar, lstar, rmax, tmin, xmloss, vmin, vmax, &
                          vmicro, vth_fiducial, vrot, beta, yhe, hei, xnue0, na

namelist / input_infreg / rmin, rlim
namelist / input_cont / eps_cont, kcont
!
!-----------------------------------------------------------------------
!
write(*,*) '-----------------------------read input----------------------------------------'
write(*,*)
write(*,*) 'input file name (*.nml) to define model-atmosphere'
read(*,*) indat_file
write(*,*) 'reading input from file: ', trim(indat_file)
write(*,*)
!
!-----------------------------------------------------------------------
!
open(1, file=trim(indat_file), status='old', form='formatted')

!read options, i.e. which model shall be used
   rewind 1
   read(1, nml=input_options)
!
!read 1d model parameters
   rewind 1
   read(1, nml=input_mod_1d)
!
   sr=rstar*rsu
   tmin=tmin*teff
!
   rewind 1
   read(1, nml=input_cont)
!
!-----------------------------------------------------------------------
!
   rewind 1
   read(1, nml=input_infreg)
!
!
close(1)
!
!
!
end subroutine read_input
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine print_input
!
use prog_type
use fund_const
use params_input
use params_stellar
!
IMPLICIT NONE
!
write(*,*) '-------------------------summary of input parameter----------------------------'
write(*,*)
!
write(*,'(a20, i20)') 'input model', input_mod
write(*,*)
write(*,*) 'GENERAL INPUT'
write(*,'(a20, 2es20.8)') 'r_star [cm], [r_sun]', sr, rstar
write(*,'(a20, es20.8)') 'v_min [km/s]', vmin
write(*,'(a20, es20.8)') 'v_max [km/s]', vmax
write(*,'(a20, es20.8)') 'v_micro [km/s]', vmicro
write(*,'(a20, es20.8)') 'vrot [km/s]', vrot
write(*,'(a20, es20.8)') 'beta', beta
write(*,'(a20, es20.8)') 'xlogg', xlogg
write(*,'(a20, es20.8)') 'mdot', xmloss
write(*,'(a20, es20.8)') 'teff', teff
write(*,'(a20, es20.8)') 'yhe', yhe
write(*,'(a20, es20.8)') 'hei', hei
write(*,'(a20, i20)') 'na', na
write(*,'(a20, es20.8)') 't_min [K]', tmin
write(*,'(a20, es20.8)') 't_eff [K]', teff
write(*,*)
write(*,*)
!
!
!
end subroutine print_input
