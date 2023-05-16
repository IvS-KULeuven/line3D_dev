!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
module options_modspec
!
   use prog_type
!
   implicit none
!
   character(len=500) :: indat_file, input_file, input_file2, output_file
!
   integer(i4b) :: input_mod
!input_mod=0  1D model: output is read from input-file (from JO's program),
!                       and converted to standard-output that can be used by spec.eo
!input_mod=1  1D model: velocities from beta-velocity law
!                       opacities from line-strength parameterization
!                       source function from sobolev approach
!input_mod=2  1D model: velocities from beta-velocity law
!                       opacities from hamann-parameterization
!                       source function from sobolev approach
!input_mod=3  1D model: velocities from benchmark08-file
!                       opacities from benchmark08-file
!                       source function from benchmark08-file (2d SC solution)
!input_mod=4  1D model: velocities from benchmark08-file
!                       opacities from benchmark08-file
!                       source function from benchmark08-file (2d FVM solution)
!input_mod=5  1D model: velocities from benchmark08-file
!                       opacities from benchmark08-file
!                       source function from benchmark08-file (from JO's program)
!input_mod=6  1D model: velocities from beta velocity law,
!                       opacities and source functions are calculated according to
!                          petrenz&puls 1995 (h-alpha)
!input_mod=7  3D model: velocities from input h5-file
!                       opacities from input h5-file
!                       source function from input h5-file (3d SC solution)
!input_mod=8  3D model: velocities from input h5-file
!                       opacities from input h5-file
!                       source function from input h5-file (3d FVM solution)
!input_mod=9  3D model: velocities from input h5-file
!                       opacities from input h5-file
!                       source function from input h5-file (from JO's program on 3d grid)
!input_mod=10  2D model: velocities and densities from input h5-file
!                        opacities and source functions are calculated according to
!                          petrenz&puls 1995 (h-alpha)
!input_mod=11  3d model: standard ouput from sc3c.eo (3d cartesian model)
!input_mod=12  3d model: standard ouput from sc3c.eo (3d cartesian model) interpolated onto spherical grid
!
!
!
end module options_modspec
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
module dime_modspec
!
   use prog_type
!
   implicit none
!
!--------------------dimension of 1d spherical grid---------------------
!
   integer(i4b) :: n1d
   real(dp), dimension(:), allocatable :: r1d, velr1d, opalbar1d, opac1d, &
      t1d, vth1d, sline1d, scont1d
!
!---------------------dimensions of cartesian 3d grid-------------------
!-----------------and-dimensions of spherical 3d grid-------------------
!
   integer(i4b) :: ndxmax, ndymax, ndzmax
   integer(i4b) :: nr, ntheta, nphi
   real(dp), dimension(:), allocatable :: x, y, z, r, theta, phi
   integer, dimension(:,:,:), allocatable :: imask3d
   real(dp), dimension(:,:,:), allocatable :: velx3d, vely3d, velz3d, opac3d, opalbar3d, &
      sline3d, scont3d, t3d, trad3d
!
end module dime_modspec
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
module params_modspec
!
!
   use prog_type
!
   implicit none
!
   real(dp) :: unit_length
!
   real(dp) :: xlogg, lstar, rstar, mstar, rmin, rmax, xmloss, yhe, yhe_mass, hei, sr
!yhe = n_He/n_H, yhe_mass = m_He/m_tot
!
   real(dp) :: teff, tmin, trad, xic1, xic2
!
   real(dp) :: vrot, vmin, vmax, vmicro, beta, vth_fiducial
!
   real(dp) :: eps_line
!
   integer(i4b) :: opt_opal
!option which opacity law has been used
!
!
end module params_modspec
