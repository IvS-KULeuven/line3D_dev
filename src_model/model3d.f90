subroutine calc_mod3d_mhd
!
use prog_type
use hdf5
use dime_modext, only: nr_modext, ntheta_modext, nphi_modext
use model3d, only: r_modext3d, theta_modext3d, phi_modext3d, rho_modext3d, t_modext3d, trad_modext3d, &
                   velr_modext3d, velth_modext3d, velphi_modext3d, vth_modext3d, eps_cont_modext3d
use params_input, only: teff, vmicro, na, eps_cont
!
implicit none
!
! ... local scalars
integer(hid_t) :: file_id, group_id, dset_id
integer :: i, err
real(dp) :: vth
!
! ... local characters
character(len=28), parameter :: model3d_mhd='models/mhd/3d/model_mhd3d.h5'
!character(len=26), parameter :: model3d_mhd='models/mhd/3d/model_060.h5'
!character(len=22), parameter :: model3d_mhd='models/mhd/2d/model.h5'
!model3d_mhd: input-file for Asif's MHD model
!
character(len=10), parameter :: gname_dimensions='dimensions'
character(len=11), parameter :: gname_coords='coordinates'
character(len=5), parameter :: gname_model='model'
character(len=2), parameter :: aname_nr='nr'
character(len=6), parameter :: aname_nth='ntheta'
character(len=4), parameter :: aname_nphi='nphi'
character(len=3), parameter :: dname_density='rho'
character(len=11), parameter :: dname_temp='temperature'
character(len=6), parameter :: dname_radius='radius'
character(len=5), parameter :: dname_theta='theta'
character(len=3), parameter :: dname_phi='phi'
character(len=5), parameter :: dname_velr='vel_r'
character(len=6), parameter :: dname_velth='vel_th'
character(len=7), parameter :: dname_velphi='vel_phi'
!
! ... local arrays
integer(hsize_t), dimension(1), parameter :: dims0 = (/ 1 /)
integer(hsize_t), dimension(1) :: dims_radius, dims_theta, dims_phi
integer(hsize_t), dimension(3) :: dims_3d
!
! ... local functions
real(dp) :: vthermal
!
!

write(*,*) '---------------------read model from 3d MHD atmosphere-------------------------'
write(*,*) 'file name: ', model3d_mhd
write(*,*)
!
call h5open_f(err)
call h5fopen_f(model3d_mhd, h5f_acc_rdonly_f, file_id, err)
!
!--------------------------read dimensions------------------------------
!
call h5gopen_f(file_id, gname_dimensions, group_id, err)
!
   call h5aopen_f(group_id, aname_nr, dset_id, err)
      call h5aread_f(dset_id, H5T_NATIVE_INTEGER, nr_modext, dims0, err)
   call h5aclose_f(dset_id, err)
!
   call h5aopen_f(group_id, aname_nth, dset_id, err)
      call h5aread_f(dset_id, H5T_NATIVE_INTEGER, ntheta_modext, dims0, err)
   call h5aclose_f(dset_id, err)
!
   call h5aopen_f(group_id, aname_nphi, dset_id, err)
      call h5aread_f(dset_id, H5T_NATIVE_INTEGER, nphi_modext, dims0, err)
   call h5aclose_f(dset_id, err)
!
call h5gclose_f(group_id, err)
!
!
!set dimension-arrays
dims_radius=(/ nr_modext /)
dims_theta=(/ ntheta_modext /)
dims_phi=(/ nphi_modext /)
dims_3d=(/nr_modext, ntheta_modext, nphi_modext /)
!
!-----------------------coordinates-------------------------------------
!
!allocate arrays
allocate(r_modext3d(nr_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod3d_mhd'
allocate(theta_modext3d(ntheta_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod3d_mhd'
allocate(phi_modext3d(nphi_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod3d_mhd'
!
call h5gopen_f(file_id, gname_coords, group_id, err)
!
   call h5dopen_f(group_id, dname_radius, dset_id, err)
      call h5dread_f(dset_id, h5t_native_double, r_modext3d, dims_radius, err)
   call h5dclose_f(dset_id, err)
!
   call h5dopen_f(group_id, dname_theta, dset_id, err)
      call h5dread_f(dset_id, h5t_native_double, theta_modext3d, dims_theta, err)
   call h5dclose_f(dset_id, err)
!
   call h5dopen_f(group_id, dname_phi, dset_id, err)
      call h5dread_f(dset_id, h5t_native_double, phi_modext3d, dims_phi, err)
   call h5dclose_f(dset_id, err)
!
call h5gclose_f(group_id, err)
!
!--------------------------model----------------------------------------
!   
allocate(rho_modext3d(nr_modext, ntheta_modext, nphi_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod3d_mhd'
allocate(t_modext3d(nr_modext, ntheta_modext, nphi_modext), stat=err)
if(err.ne.0) stop 'error: allocation in calc_mod3d_mhd'
allocate(trad_modext3d(nr_modext, ntheta_modext, nphi_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod3d_mhd'
allocate(velr_modext3d(nr_modext, ntheta_modext, nphi_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod3d_mhd'
allocate(velth_modext3d(nr_modext, ntheta_modext, nphi_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod3d_mhd'
allocate(velphi_modext3d(nr_modext, ntheta_modext, nphi_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod3d_mhd'
!
call h5gopen_f(file_id, gname_model, group_id, err)
!
   call h5dopen_f(group_id, dname_density, dset_id, err)
      call h5dread_f(dset_id, h5t_native_double, rho_modext3d, dims_3d, err)
   call h5dclose_f(dset_id, err)
!
   call h5dopen_f(group_id, dname_temp, dset_id, err)
      call h5dread_f(dset_id, h5t_native_double, t_modext3d, dims_3d, err)
   call h5dclose_f(dset_id, err)
   trad_modext3d=t_modext3d   
!
   call h5dopen_f(group_id, dname_velr, dset_id, err)
      call h5dread_f(dset_id, h5t_native_double, velr_modext3d, dims_3d, err)
   call h5dclose_f(dset_id, err)
!
   call h5dopen_f(group_id, dname_velth, dset_id, err)
      call h5dread_f(dset_id, h5t_native_double, velth_modext3d, dims_3d, err)
   call h5dclose_f(dset_id, err)
!
   call h5dopen_f(group_id, dname_velphi, dset_id, err)
      call h5dread_f(dset_id, h5t_native_double, velphi_modext3d, dims_3d, err)
   call h5dclose_f(dset_id, err)
!
call h5gclose_f(group_id, err)
!
call h5fclose_f(file_id, err)
call h5close_f(err)
!
!
!need to set thermal velocities
allocate(vth_modext3d(nr_modext, ntheta_modext, nphi_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod3d_mhd'
allocate(eps_cont_modext3d(nr_modext, ntheta_modext, nphi_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod3d_mhd'
vth = vthermal(vmicro*1.d5, teff, na)
vth_modext3d = vth
eps_cont_modext3d = eps_cont
!
!
!
end subroutine calc_mod3d_mhd
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine calc_mod3d_christi
!
use prog_type
use fund_const
use hdf5
use dime_modext, only: nr_modext, ntheta_modext, nphi_modext
use model3d, only: r_modext3d, theta_modext3d, phi_modext3d, rho_modext3d, t_modext3d, trad_modext3d, &
                   velr_modext3d, velth_modext3d, velphi_modext3d, vth_modext3d, eps_cont_modext3d
use params_input, only: teff, vmicro, na, eps_cont
use params_stellar, only: sr
!
implicit none
!
! ... local scalars
integer(i4b) :: i, j, k, nx, ny, nz, err
integer(i4b) :: iim2, iim1, ii, iip1, jjm2, jjm1, jj, jjp1, kkm2, kkm1, kk, kkp1
real(dp) :: xp, yp, zp, sint, cost, sinp, cosp
real(dp) :: velx, vely, velz, rho, temp, rho_min, t_min
real(dp) :: vth, rmin, rmax, bx, by, bz, fdum
!
! ... local characters
character(len=26), parameter :: model3d_mhd='models/christi/hdfaa200.h5'
character(len=50) :: cdum
!model3d_mhd: input-file for Christis MHD model
!
! ... local arrays
real(dp), dimension(:), allocatable :: x, y, z
real(dp), dimension(:,:,:), allocatable :: velx3d, vely3d, velz3d, rho3d, t3d
!
! ... local functions
real(dp) :: vthermal, interpol3d_8p_lin
!
! ... for hdf5
integer(hid_t) :: file_id, dset_id
integer(hsize_t), dimension(1), parameter :: dims0 = (/ 1 /)
integer(hsize_t), dimension(1) :: dims_radius, dims_theta, dims_phi
integer(hsize_t), dimension(3) :: dims_3d
!

write(*,*) '---------------------read model from 3d MHD atmosphere-------------------------'
write(*,*) 'file name: ', model3d_mhd
write(*,*)
!
!
!read in all data
!
!
call h5open_f(err)
call h5fopen_f(model3d_mhd, h5f_acc_rdonly_f, file_id, err)
!
!--------------------------read dimensions------------------------------
!
call h5dopen_f(file_id, 'nr', dset_id, err)
   call h5dread_f(dset_id, h5t_native_integer, nr_modext, dims0, err)
call h5dclose_f(dset_id, err)

call h5dopen_f(file_id, 'ntheta', dset_id, err)
   call h5dread_f(dset_id, h5t_native_integer, ntheta_modext, dims0, err)
call h5dclose_f(dset_id, err)

call h5dopen_f(file_id, 'nphi', dset_id, err)
   call h5dread_f(dset_id, h5t_native_integer, nphi_modext, dims0, err)
call h5dclose_f(dset_id, err)
!
!-----------------------coordinates-------------------------------------
!
!set dimension-arrays
dims_radius=(/ nr_modext /)
dims_theta=(/ ntheta_modext /)
dims_phi=(/ nphi_modext /)
dims_3d=(/nr_modext, ntheta_modext, nphi_modext /)   
!allocate arrays
allocate(r_modext3d(nr_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod3d_christi'
allocate(theta_modext3d(ntheta_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod3d_christi'
allocate(phi_modext3d(nphi_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod3d_christi'
!

call h5dopen_f(file_id, 'radius', dset_id, err)
   call h5dread_f(dset_id, h5t_native_double, r_modext3d, dims_radius, err)
call h5dclose_f(dset_id, err)

call h5dopen_f(file_id, 'theta', dset_id, err)
   call h5dread_f(dset_id, h5t_native_double, theta_modext3d, dims_theta, err)
call h5dclose_f(dset_id, err)

call h5dopen_f(file_id, 'phi', dset_id, err)
   call h5dread_f(dset_id, h5t_native_double, phi_modext3d, dims_phi, err)
call h5dclose_f(dset_id, err)
!
!manipulate the angular grid to have complete boundaries
theta_modext3d(1)=0.d0
theta_modext3d(ntheta_modext)=pi
phi_modext3d(1)=0.d0
phi_modext3d(nphi_modext)=2.d0*pi
!
!--------------------------model----------------------------------------
!   
allocate(rho_modext3d(nr_modext, ntheta_modext, nphi_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod3d_mhd'
allocate(t_modext3d(nr_modext, ntheta_modext, nphi_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod3d_mhd'
allocate(trad_modext3d(nr_modext, ntheta_modext, nphi_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod3d_mhd'
allocate(velr_modext3d(nr_modext, ntheta_modext, nphi_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod3d_mhd'
allocate(velth_modext3d(nr_modext, ntheta_modext, nphi_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod3d_mhd'
allocate(velphi_modext3d(nr_modext, ntheta_modext, nphi_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod3d_mhd'
!
!
call h5dopen_f(file_id, 'rho3d', dset_id, err)
   call h5dread_f(dset_id, h5t_native_double, rho_modext3d, dims_3d, err)
call h5dclose_f(dset_id, err)
!
call h5dopen_f(file_id, 't3d', dset_id, err)
   call h5dread_f(dset_id, h5t_native_double, t_modext3d, dims_3d, err)
call h5dclose_f(dset_id, err)
trad_modext3d=t_modext3d   
!
call h5dopen_f(file_id, 'velr3d', dset_id, err)
   call h5dread_f(dset_id, h5t_native_double, velr_modext3d, dims_3d, err)
call h5dclose_f(dset_id, err)
!
call h5dopen_f(file_id, 'velth3d', dset_id, err)
   call h5dread_f(dset_id, h5t_native_double, velth_modext3d, dims_3d, err)
call h5dclose_f(dset_id, err)
!
call h5dopen_f(file_id, 'velphi3d', dset_id, err)
   call h5dread_f(dset_id, h5t_native_double, velphi_modext3d, dims_3d, err)
call h5dclose_f(dset_id, err)
!
!
call h5fclose_f(file_id, err)
call h5close_f(err)
!
!manipulate temperatures
t_modext3d=teff
trad_modext3d=t_modext3d
!
!need to set thermal velocities
allocate(vth_modext3d(nr_modext, ntheta_modext, nphi_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod3d_mhd'
allocate(eps_cont_modext3d(nr_modext, ntheta_modext, nphi_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod3d_mhd'
vth = vthermal(vmicro*1.d5, teff, na)
vth_modext3d = vth
eps_cont_modext3d = eps_cont
!!
!!
!j=ntheta_modext/2
!k=1
!do i=1, nr_modext
!   write(*,*) r_modext3d(i)/sr, rho_modext3d(i,j,k), velr_modext3d(i,j,k)/1.d5, velth_modext3d(i,j,k)/1.d5, velphi_modext3d(i,j,k)/1.d5
!enddo
!stop
!
!
!
end subroutine calc_mod3d_christi
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine calc_mod3d_cis
!
use prog_type
use fund_const, only: cgs_mp, cgs_kb, pi
use dime_modext, only: nr_modext, ntheta_modext, nphi_modext
use model3d, only: r_modext3d, theta_modext3d, phi_modext3d, rho_modext3d, t_modext3d, trad_modext3d, &
                   velr_modext3d, velth_modext3d, velphi_modext3d, vth_modext3d, eps_cont_modext3d
use params_input, only: teff, vmicro, na, yhe, hei, eps_cont
use params_stellar, only: sr
use mod_grid, only: grid_equi
!
implicit none
!
! ... local scalars
!for b dwarf
!integer(i4b), parameter :: nr_1d_cis=2000, tstep_min=32, tstep_max=432
!for o dwarf
integer(i4b), parameter :: nr_1d_cis=2000, tstep_min=32, tstep_max=332
!for o star
!integer(i4b), parameter :: nr_1d_cis=2000, tstep_min=22, tstep_max=122
integer(i4b) :: i,j,k
integer(i4b) :: err, nt_cis, tstep, idum
real(dp) :: cconst1, vth
real(dp) :: fdum
!
! ... local characters
!for b dwarf
!character(len=41), parameter :: model1d_cis='models/cis/output_bstar_cooling_Feldmeier'
!for o dwarf (no cooling correction)
!character(len=35), parameter :: model1d_cis='models/cis/output_ostar_100xQmax_rc'
!for o dwarf (with cooling correction)
character(len=45), parameter :: model1d_cis='models/cis/output_ostar_100xQmax_rc_feldmeier'
!for o star
!character(len=34), parameter :: model1d_cis='models/cis/output_Ostar_rc_15xTeff'
character(len=500) :: fname
character(len=4) :: chdum
!
! ... local logicals
logical :: check1
!
! ... local arrays
real(dp), dimension(:), allocatable :: p_modext1d  !pressure
real(dp), dimension(:,:), allocatable :: r2d, rho2d, velr2d, t2d, p2d
integer(i4b), dimension(64), parameter :: seed = (/ &
               676865, 350848,  32171, 881468, 615122, 479556, 455219,  89650, &
               835154, 575764, 830269, 479197, 746278, 918314,  62049, 352171, &
               489558, 465378, 903709, 768839, 104036, 606566, 706051, 168587, &
               972615, 980707, 564674, 112480, 559149, 991281, 744474, 619200, &
               684163, 145540,   8131, 896623, 257203, 761081,  44394, 861176, &
               601220, 145688, 366560, 843924, 941293, 376982,   5005, 405884, &
               170702, 848387, 432408, 742485, 239070, 782087, 573526,   7287, &
               430405, 352846, 747140, 620387, 818214, 290087, 335369, 644222 /)
!
! ... local functions
real(dp) :: vthermal
!
!
!-----------------------allocate arrays---------------------------------
!
nr_modext=nr_1d_cis
nt_cis=tstep_max-tstep_min+1
!
!standard clump size
ntheta_modext=51
!large and few clumps
!ntheta_modext=25
!small and many clumps
!ntheta_modext=101
nphi_modext=2*ntheta_modext-1
!
!arrays with timestamps
allocate(r2d(nr_modext,nt_cis), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod1d_cis'
allocate(rho2d(nr_modext,nt_cis), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod1d_cis'
allocate(velr2d(nr_modext,nt_cis), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod1d_cis'
allocate(t2d(nr_modext,nt_cis), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod1d_cis'
allocate(p2d(nr_modext,nt_cis), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod1d_cis'
!
!coordinates
allocate(r_modext3d(nr_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod3d_cis'
allocate(theta_modext3d(ntheta_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod3d_cis'
allocate(phi_modext3d(nphi_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod3d_cis'
!
!--------------------------model----------------------------------------
!   
allocate(rho_modext3d(nr_modext, ntheta_modext, nphi_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod3d_cis'
allocate(t_modext3d(nr_modext, ntheta_modext, nphi_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod3d_cis'
allocate(trad_modext3d(nr_modext, ntheta_modext, nphi_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod3d_cis'
allocate(velr_modext3d(nr_modext, ntheta_modext, nphi_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod3d_cis'
allocate(velth_modext3d(nr_modext, ntheta_modext, nphi_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod3d_cis'
allocate(velphi_modext3d(nr_modext, ntheta_modext, nphi_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod3d_cis'
allocate(vth_modext3d(nr_modext, ntheta_modext, nphi_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod3d_cis'
allocate(eps_cont_modext3d(nr_modext, ntheta_modext, nphi_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod3d_cis'
!
!------------------calculate/define required constants------------------
!
!mean molecular weight
cconst1=(1.d0+4.d0*yhe)/(2.d0+(1.d0+hei)*yhe)
!
!-------------------calculate/define angular grids----------------------
!
call grid_equi(0.d0,2.d0*pi,nphi_modext,phi_modext3d)
!call grid_equi(1.d0,-1.d0,ntheta_modext,theta_modext3d)
!theta_modext3d=acos(theta_modext3d)
call grid_equi(0.d0,pi,ntheta_modext,theta_modext3d)
!
!---------------read in all timesteps of 1d atmosheric structure--------
!
do j=1, nt_cis
   tstep = tstep_min+j-1
!
!get the filename for tstep
   write(chdum,'(i4.4)') tstep
   write(fname,*) './'//model1d_cis//'/LDI_'//chdum//'.blk'
!delete spaces on left side of string (why are those occurring?)   
   fname=adjustl(fname)
   inquire(file=trim(fname), exist=check1)
   if(.not.check1) then
!      write(*,*) trim(fname)
      write(*,*) 'error in calc_mod1d_cis: file "', trim(fname), '" does not exist'
      stop
   endif
   !
   open(1, file=trim(fname), form='formatted')
      read(1,*) chdum
      read(1,*) chdum
      read(1,*) chdum
      do i=1, nr_modext
         read(1,'(7e14.6)') r2d(i,j), rho2d(i,j), velr2d(i,j), p2d(i,j), fdum, fdum, fdum
         t2d(i,j) = cconst1*p2d(i,j)*cgs_mp/rho2d(i,j)/cgs_kb
!replace negative velocities by (almost) zero
         if(velr2d(i,j).le.0.) then
            velr2d(i,j) = 1.d-10
         endif
      enddo
!      if(minval(velr2d).lt.0.) then
!         k=k+1
!         write(*,*) j, fname, minval(velr2d)/1.d5
!      endif
   close(1)
   !
enddo

!write(*,*) k, nt_cis
!stop 'go on in here'
!
!rescaling to match r=1 (only small correction due to numerical rounding errors,
!r=1 is actually set via input parameter rstar)
r_modext3d=r2d(:,1)/r2d(1,1)*sr
!
!-----------------------------------------------------------------------
!
!always use same seed in order that reproducable results
call random_seed(put=seed)
do i=1, ntheta_modext
!   write(*,*) theta_modext3d(i)
   do j=1, nphi_modext
      call random_number(fdum)
      idum=nint(1. + fdum*(nt_cis-1))
      if(idum.lt.1) stop 'error in calc_mod3d_cis: random indx lt 1'
      if(idum.gt.nt_cis) stop 'error in calc_mod3d_cis: random indx gt maximum allowed'
!
!set density etc. at (theta,phi) to 1D structure at random time step
      rho_modext3d(:,i,j)=rho2d(:,idum)
      velr_modext3d(:,i,j)=velr2d(:,idum)
      t_modext3d(:,i,j)=t2d(:,idum)
      trad_modext3d(:,i,j)=t2d(:,idum)
   enddo
enddo
!
velth_modext3d=0.d0
velphi_modext3d=0.d0
!
!constant thermal velocity
vth = vthermal(vmicro*1.d5, teff, na)
vth_modext3d = vth
!
!constant thermalization parameter
eps_cont_modext3d = eps_cont
!
!
end subroutine calc_mod3d_cis
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine calc_mod3d_florian
!
! reads in Florian's 2D MHD/LDI simulations from is_min to is_max using random snapshots
!
use prog_type
use fund_const, only: pi, zero, one, two
use dime_modext, only: nr_modext, ntheta_modext, nphi_modext
use model3d, only: r_modext3d, theta_modext3d, phi_modext3d, t_modext3d, trad_modext3d, rho_modext3d, velr_modext3d, &
                   velth_modext3d, velphi_modext3d, vth_modext3d, eps_cont_modext3d
use params_input, only: teff, vmicro, na, eps_cont
use mod_grid, only: grid_equi
use mod_directories, only: indat_file
use mod_sort, only: arr_equal
use hdf5
!
implicit none
!
! ... local scalars
integer(i4b) :: i, j, k
integer(i4b) :: idum, err, nr_test, ntheta_test
integer(i4b) :: is_min, is_max, nphi, opt_bvel
real(dp) :: fdum, fsum
!
! ... local logicals
logical :: lcheck
! ... local arrays
real(dp), dimension(:), allocatable :: r_modext2d, theta_modext2d
real(dp), dimension(:), allocatable :: t_modext1d, rho_modext1d, velr_modext1d, &
                                         velth_modext1d, velphi_modext1d, vth_modext1d
real(dp), dimension(:,:), allocatable :: t_modext2d, rho_modext2d, velr_modext2d, &
                                         velth_modext2d, velphi_modext2d, vth_modext2d
!
! ... local characters!
!character(len=32), parameter :: model3d_florian='models/florian/100G_zetapup3d.h5'
!character(len=33), parameter :: model3d_florian='models/florian/1000G_zetapup3d.h5'
!character(len=31), parameter :: model3d_florian='models/florian/cak_zetapup3d.h5'
!character(len=30), parameter :: model3d_florian='models/florian/0G_zetapup3d.h5'
character(len=500)  :: fname_model, fname
!
! ... namelists
namelist / input_usr / fname_model, is_min, is_max, nphi, opt_bvel
!
! ... local functions
real(dp) :: vthermal
!
! ... for hdf5-files
integer(hid_t) :: file_id, group_id, dset_id, attr_id
integer(hsize_t), dimension(1) :: dims_scalars
integer(hsize_t), dimension(1) :: dims_rad
integer(hsize_t), dimension(1) :: dims_lat
integer(hsize_t), dimension(1) :: dims_az
integer(hsize_t), dimension(2) :: dims2d
integer(hsize_t), dimension(3) :: dims3d
!
! ... for random numbers
integer(i4b), dimension(64), parameter :: seed = (/ &
               676865, 350848,  32171, 881468, 615122, 479556, 455219,  89650, &
               835154, 575764, 830269, 479197, 746278, 918314,  62049, 352171, &
               489558, 465378, 903709, 768839, 104036, 606566, 706051, 168587, &
               972615, 980707, 564674, 112480, 559149, 991281, 744474, 619200, &
               684163, 145540,   8131, 896623, 257203, 761081,  44394, 861176, &
               601220, 145688, 366560, 843924, 941293, 376982,   5005, 405884, &
               170702, 848387, 432408, 742485, 239070, 782087, 573526,   7287, &
               430405, 352846, 747140, 620387, 818214, 290087, 335369, 644222 /)
!

!
write(*,*) '---------------------read 2d random models from florian------------------------'
write(*,*) 'reading input from file: ', trim(indat_file)
write(*,*)
!
open(1, file=trim(indat_file), status='old', form='formatted')
   rewind 1
   read(1, nml=input_usr)
close(1)
!
nphi_modext=nphi
!
write(*,*) 'file base name: ', trim(fname_model)
write(*,*) 'is_min, is_max:', is_min, is_max
write(*,*) 'nphi', nphi_modext
write(*,*) 'opt_bvel', opt_bvel
write(*,*) 'for opt_bvel = 1 ', 'using random snapshots from is_min to is_max (3D model)'
write(*,*) 'for opt_bvel = 2 ', 'averaging all random snapshots from is_min to is_max (2D model)'
write(*,*) 'for opt_bvel = 3 ', 'averaging all random snapshots from is_min to is_max (1D model)'
write(*,*)
!
!read in dimensions and set coordinates for first snapshot
idum = is_min
!
!---------------------read in dimensions--------------------------------
!
write(fname,'(a,i4.4,a)') trim(fname_model), idum, '.h5'
!
fname=adjustl(fname)
inquire(file=trim(fname), exist=lcheck)
if(.not.lcheck) then
   write(*,*) 'error in calc_mod3d_nicowr3d: file "', trim(fname), '" does not exist'
   stop
endif
!
write(*,*) 'file name: ', trim(fname)
write(*,*)
!
call h5open_f(err)
call h5fopen_f(trim(fname), h5f_acc_rdonly_f, file_id, err)
!
call h5gopen_f(file_id, 'dimensions', group_id, err)
   call h5aopen_f(group_id, 'nr', attr_id, err)
      call h5aread_f(attr_id, h5t_native_integer, nr_modext, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5aopen_f(group_id, 'ntheta', attr_id, err)
      call h5aread_f(attr_id, h5t_native_integer, ntheta_modext, dims_scalars, err)
   call h5aclose_f(attr_id, err)
call h5gclose_f(group_id, err)
!
!-------------------------get coordiantes-------------------------------
!
allocate(r_modext3d(nr_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod3d_florian'
allocate(theta_modext3d(ntheta_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod3d_florian'
allocate(phi_modext3d(nphi_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod3d_florian'
!
call h5gopen_f(file_id, 'coordinates', group_id, err)
   call h5dopen_f(group_id, 'r', dset_id, err)
      call h5dread_f(dset_id, h5t_native_double, r_modext3d, dims_rad, err)
   call h5dclose_f(dset_id, err)
   call h5dopen_f(group_id, 'theta', dset_id, err)
      call h5dread_f(dset_id, h5t_native_double, theta_modext3d, dims_lat, err)
   call h5dclose_f(dset_id, err)
call h5gclose_f(group_id, err)

call h5fclose_f(file_id, err)
call h5close_f(err)

call grid_equi(zero,two*pi,nphi_modext,phi_modext3d)   
!
!-----------------------allocate arrays---------------------------------
!
allocate(velr_modext3d(nr_modext, ntheta_modext, nphi_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod3d_florian'
allocate(velth_modext3d(nr_modext, ntheta_modext, nphi_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod3d_florian'
allocate(velphi_modext3d(nr_modext, ntheta_modext, nphi_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod3d_florian'
allocate(rho_modext3d(nr_modext, ntheta_modext, nphi_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod3d_florian'
allocate(t_modext3d(nr_modext, ntheta_modext, nphi_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod3d_florian'
allocate(trad_modext3d(nr_modext, ntheta_modext, nphi_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod3d_florian'
allocate(vth_modext3d(nr_modext, ntheta_modext, nphi_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod3d_florian'
allocate(eps_cont_modext3d(nr_modext, ntheta_modext, nphi_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod3d_florian'

allocate(r_modext2d(nr_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod3d_florian'
allocate(theta_modext2d(ntheta_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod3d_florian'
allocate(velr_modext2d(nr_modext, ntheta_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod3d_florian'
allocate(velth_modext2d(nr_modext, ntheta_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod3d_florian'
allocate(velphi_modext2d(nr_modext, ntheta_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod3d_florian'
allocate(rho_modext2d(nr_modext, ntheta_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod3d_florian'
allocate(t_modext2d(nr_modext, ntheta_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod3d_florian'
allocate(vth_modext2d(nr_modext, ntheta_modext), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod3d_florian'
!
dims_rad= (/ nr_modext /)
dims_lat= (/ ntheta_modext /)
dims2d=(/ nr_modext, ntheta_modext /)
!
!------------------read in random snapshots-----------------------------
!
if (opt_bvel.eq.1) then
!
!always use same seed in order that reproducable results   
   call random_seed(put=seed)
   !
   do k=1, nphi

      call random_number(fdum)
      idum=nint(is_min + fdum*(is_max-is_min))
      if(idum.lt.is_min) stop 'error in calc_mod3d_nicowr3d: random indx lt is_min'
      if(idum.gt.is_max) stop 'error in calc_mod3d_nicowr3d: random indx gt is_max'


      write(fname,'(a,i4.4,a)') trim(fname_model), idum, '.h5'
      !
      fname=adjustl(fname)
      inquire(file=trim(fname), exist=lcheck)
      if(.not.lcheck) then
         write(*,*) 'error in calc_mod3d_nicowr3d: file "', trim(fname), '" does not exist'
         stop
      endif
      !
      write(*,*) 'phi indx, nphi, file name: ', k, nphi, idum, trim(fname)
      !
      !
      call h5open_f(err)
      call h5fopen_f(trim(fname), h5f_acc_rdonly_f, file_id, err)

      !read dimensions
      call h5gopen_f(file_id, 'dimensions', group_id, err)
         call h5aopen_f(group_id, 'nr', attr_id, err)
            call h5aread_f(attr_id, h5t_native_integer, nr_test, dims_scalars, err)
         call h5aclose_f(attr_id, err)
         call h5aopen_f(group_id, 'ntheta', attr_id, err)
            call h5aread_f(attr_id, h5t_native_integer, ntheta_test, dims_scalars, err)
         call h5aclose_f(attr_id, err)
      call h5gclose_f(group_id, err)

      if(nr_test.ne.nr_modext) stop 'error in calc_mod3d_florian: nr not matching'
      if(ntheta_test.ne.ntheta_modext) stop 'error in calc_mod3d_florian: ntheta not matching'
      
      call h5gopen_f(file_id, 'coordinates', group_id, err)
         call h5dopen_f(group_id, 'r', dset_id, err)
            call h5dread_f(dset_id, h5t_native_double, r_modext2d, dims_rad, err)
         call h5dclose_f(dset_id, err)
         call h5dopen_f(group_id, 'theta', dset_id, err)
            call h5dread_f(dset_id, h5t_native_double, theta_modext2d, dims_rad, err)
         call h5dclose_f(dset_id, err)
      call h5gclose_f(group_id, err)
      lcheck = arr_equal(r_modext2d, r_modext3d, nr_modext)
      if(.not.lcheck) stop 'error in calc_mod3d_florian: r-coordinates not equal'
      lcheck = arr_equal(theta_modext2d, theta_modext3d, ntheta_modext)      
      if(.not.lcheck) stop 'error in calc_mod3d_florian: theta-coordinates not equal'      

      call h5gopen_f(file_id, 'model', group_id, err)
         call h5dopen_f(group_id, 'rho', dset_id, err)
            call h5dread_f(dset_id, h5t_native_double, rho_modext2d, dims2d, err)
         call h5dclose_f(dset_id, err)
         call h5dopen_f(group_id, 'temperature', dset_id, err)
            call h5dread_f(dset_id, h5t_native_double, t_modext2d, dims2d, err)
         call h5dclose_f(dset_id, err)
         call h5dopen_f(group_id, 'velr', dset_id, err)
            call h5dread_f(dset_id, h5t_native_double, velr_modext2d, dims2d, err)
         call h5dclose_f(dset_id, err)
         call h5dopen_f(group_id, 'velth', dset_id, err)
            call h5dread_f(dset_id, h5t_native_double, velth_modext2d, dims2d, err)
         call h5dclose_f(dset_id, err)
         call h5dopen_f(group_id, 'velphi', dset_id, err)
            call h5dread_f(dset_id, h5t_native_double, velphi_modext2d, dims2d, err)
         call h5dclose_f(dset_id, err)
         call h5dopen_f(group_id, 'vth', dset_id, err)
            call h5dread_f(dset_id, h5t_native_double, vth_modext2d, dims2d, err)
         call h5dclose_f(dset_id, err)
      
      call h5gclose_f(group_id, err)
      call h5fclose_f(file_id, err)
      call h5close_f(err)

      !store in 3d array
      rho_modext3d(:,:,k) = rho_modext2d
      t_modext3d(:,:,k) = t_modext2d
      velr_modext3d(:,:,k) = velr_modext2d
      velth_modext3d(:,:,k) = velth_modext2d
      velphi_modext3d(:,:,k) = velphi_modext2d
      vth_modext3d(:,:,k) = vth_modext2d
   enddo
!
!------------------calculate average from all snapshots-----------------
!--------------------------(2D model)-----------------------------------   
!
elseif(opt_bvel.eq.2) then
   
   rho_modext3d = zero
   t_modext3d = zero
   velr_modext3d = zero
   velth_modext3d = zero
   velphi_modext3d = zero
   vth_modext3d = zero
   fsum = zero
   !
   do k=is_min, is_max
      idum = k

      write(fname,'(a,i4.4,a)') trim(fname_model), idum, '.h5'
      !
      fname=adjustl(fname)
      inquire(file=trim(fname), exist=lcheck)
      if(.not.lcheck) then
         write(*,*) 'error in calc_mod3d_nicowr3d: file "', trim(fname), '" does not exist'
         stop
      endif
      !
      write(*,*) 'phi indx, nphi, file name: ', k-is_min+1, is_max-is_min+1, idum, trim(fname)
      !
      !
      call h5open_f(err)
      call h5fopen_f(trim(fname), h5f_acc_rdonly_f, file_id, err)

      !read dimensions
      call h5gopen_f(file_id, 'dimensions', group_id, err)
         call h5aopen_f(group_id, 'nr', attr_id, err)
            call h5aread_f(attr_id, h5t_native_integer, nr_test, dims_scalars, err)
         call h5aclose_f(attr_id, err)
         call h5aopen_f(group_id, 'ntheta', attr_id, err)
            call h5aread_f(attr_id, h5t_native_integer, ntheta_test, dims_scalars, err)
         call h5aclose_f(attr_id, err)
      call h5gclose_f(group_id, err)

      if(nr_test.ne.nr_modext) stop 'error in calc_mod3d_florian: nr not matching'
      if(ntheta_test.ne.ntheta_modext) stop 'error in calc_mod3d_florian: ntheta not matching'
      
      call h5gopen_f(file_id, 'coordinates', group_id, err)
         call h5dopen_f(group_id, 'r', dset_id, err)
            call h5dread_f(dset_id, h5t_native_double, r_modext2d, dims_rad, err)
         call h5dclose_f(dset_id, err)
         call h5dopen_f(group_id, 'theta', dset_id, err)
            call h5dread_f(dset_id, h5t_native_double, theta_modext2d, dims_rad, err)
         call h5dclose_f(dset_id, err)
      call h5gclose_f(group_id, err)
      lcheck = arr_equal(r_modext2d, r_modext3d, nr_modext)
      if(.not.lcheck) stop 'error in calc_mod3d_florian: r-coordinates not equal'
      lcheck = arr_equal(theta_modext2d, theta_modext3d, ntheta_modext)      
      if(.not.lcheck) stop 'error in calc_mod3d_florian: theta-coordinates not equal'      

      call h5gopen_f(file_id, 'model', group_id, err)
         call h5dopen_f(group_id, 'rho', dset_id, err)
            call h5dread_f(dset_id, h5t_native_double, rho_modext2d, dims2d, err)
         call h5dclose_f(dset_id, err)
         call h5dopen_f(group_id, 'temperature', dset_id, err)
            call h5dread_f(dset_id, h5t_native_double, t_modext2d, dims2d, err)
         call h5dclose_f(dset_id, err)
         call h5dopen_f(group_id, 'velr', dset_id, err)
            call h5dread_f(dset_id, h5t_native_double, velr_modext2d, dims2d, err)
         call h5dclose_f(dset_id, err)
         call h5dopen_f(group_id, 'velth', dset_id, err)
            call h5dread_f(dset_id, h5t_native_double, velth_modext2d, dims2d, err)
         call h5dclose_f(dset_id, err)
         call h5dopen_f(group_id, 'velphi', dset_id, err)
            call h5dread_f(dset_id, h5t_native_double, velphi_modext2d, dims2d, err)
         call h5dclose_f(dset_id, err)
         call h5dopen_f(group_id, 'vth', dset_id, err)
            call h5dread_f(dset_id, h5t_native_double, vth_modext2d, dims2d, err)
         call h5dclose_f(dset_id, err)
      
      call h5gclose_f(group_id, err)
      call h5fclose_f(file_id, err)
      call h5close_f(err)

      !store in 3d array
      rho_modext3d(:,:,1) = rho_modext3d(:,:,1) + rho_modext2d
      t_modext3d(:,:,1) = t_modext3d(:,:,1) + t_modext2d
      velr_modext3d(:,:,1) = velr_modext3d(:,:,1) + velr_modext2d
      velth_modext3d(:,:,1) = velth_modext3d(:,:,1) + velth_modext2d
      velphi_modext3d(:,:,1) = velphi_modext3d(:,:,1) + velphi_modext2d
      vth_modext3d(:,:,1) = vth_modext3d(:,:,1) + vth_modext2d
      fsum = fsum + one
      !
   enddo

   rho_modext3d(:,:,1) = rho_modext3d(:,:,1)/fsum
   t_modext3d(:,:,1) = t_modext3d(:,:,1)/fsum
   velr_modext3d(:,:,1) = velr_modext3d(:,:,1)/fsum
   velth_modext3d(:,:,1) = velth_modext3d(:,:,1)/fsum
   velphi_modext3d(:,:,1) = velphi_modext3d(:,:,1)/fsum
   vth_modext3d(:,:,1) = vth_modext3d(:,:,1)/fsum

   do k=2, nphi_modext
      rho_modext3d(:,:,k) = rho_modext3d(:,:,1)
      t_modext3d(:,:,k) = t_modext3d(:,:,1)      
      velr_modext3d(:,:,k) = velr_modext3d(:,:,1)
      velth_modext3d(:,:,k) = velth_modext3d(:,:,1)
      velphi_modext3d(:,:,k) = velphi_modext3d(:,:,1)
      vth_modext3d(:,:,k) = vth_modext3d(:,:,1)
   enddo
!
!------------------calculate average from all snapshots-----------------
!--------------------------(1D model)-----------------------------------   
!
elseif(opt_bvel.eq.3) then
   !
   allocate(velr_modext1d(nr_modext), stat=err)
      if(err.ne.0) stop 'error: allocation in calc_mod3d_florian'
   allocate(velth_modext1d(nr_modext), stat=err)
      if(err.ne.0) stop 'error: allocation in calc_mod3d_florian'
   allocate(velphi_modext1d(nr_modext), stat=err)
      if(err.ne.0) stop 'error: allocation in calc_mod3d_florian'
   allocate(rho_modext1d(nr_modext), stat=err)
      if(err.ne.0) stop 'error: allocation in calc_mod3d_florian'
   allocate(t_modext1d(nr_modext), stat=err)
      if(err.ne.0) stop 'error: allocation in calc_mod3d_florian'
   allocate(vth_modext1d(nr_modext), stat=err)
      if(err.ne.0) stop 'error: allocation in calc_mod3d_florian'   
   
   rho_modext1d = zero
   t_modext1d = zero
   velr_modext1d = zero
   velth_modext1d = zero
   velphi_modext1d = zero
   vth_modext1d = zero
   fsum = zero
   !
   do k=is_min, is_max
      idum = k

      write(fname,'(a,i4.4,a)') trim(fname_model), idum, '.h5'
      !
      fname=adjustl(fname)
      inquire(file=trim(fname), exist=lcheck)
      if(.not.lcheck) then
         write(*,*) 'error in calc_mod3d_nicowr3d: file "', trim(fname), '" does not exist'
         stop
      endif
      !
      write(*,*) 'phi indx, nphi, file name: ', k-is_min+1, is_max-is_min+1, idum, trim(fname)      
      !
      !
      call h5open_f(err)
      call h5fopen_f(trim(fname), h5f_acc_rdonly_f, file_id, err)

      !read dimensions
      call h5gopen_f(file_id, 'dimensions', group_id, err)
         call h5aopen_f(group_id, 'nr', attr_id, err)
            call h5aread_f(attr_id, h5t_native_integer, nr_test, dims_scalars, err)
         call h5aclose_f(attr_id, err)
         call h5aopen_f(group_id, 'ntheta', attr_id, err)
            call h5aread_f(attr_id, h5t_native_integer, ntheta_test, dims_scalars, err)
         call h5aclose_f(attr_id, err)
      call h5gclose_f(group_id, err)

      if(nr_test.ne.nr_modext) stop 'error in calc_mod3d_florian: nr not matching'
      if(ntheta_test.ne.ntheta_modext) stop 'error in calc_mod3d_florian: ntheta not matching'
      
      call h5gopen_f(file_id, 'coordinates', group_id, err)
         call h5dopen_f(group_id, 'r', dset_id, err)
            call h5dread_f(dset_id, h5t_native_double, r_modext2d, dims_rad, err)
         call h5dclose_f(dset_id, err)
         call h5dopen_f(group_id, 'theta', dset_id, err)
            call h5dread_f(dset_id, h5t_native_double, theta_modext2d, dims_rad, err)
         call h5dclose_f(dset_id, err)
      call h5gclose_f(group_id, err)
      lcheck = arr_equal(r_modext2d, r_modext3d, nr_modext)
      if(.not.lcheck) stop 'error in calc_mod3d_florian: r-coordinates not equal'
      lcheck = arr_equal(theta_modext2d, theta_modext3d, ntheta_modext)      
      if(.not.lcheck) stop 'error in calc_mod3d_florian: theta-coordinates not equal'      

      call h5gopen_f(file_id, 'model', group_id, err)
         call h5dopen_f(group_id, 'rho', dset_id, err)
            call h5dread_f(dset_id, h5t_native_double, rho_modext2d, dims2d, err)
         call h5dclose_f(dset_id, err)
         call h5dopen_f(group_id, 'temperature', dset_id, err)
            call h5dread_f(dset_id, h5t_native_double, t_modext2d, dims2d, err)
         call h5dclose_f(dset_id, err)
         call h5dopen_f(group_id, 'velr', dset_id, err)
            call h5dread_f(dset_id, h5t_native_double, velr_modext2d, dims2d, err)
         call h5dclose_f(dset_id, err)
         call h5dopen_f(group_id, 'velth', dset_id, err)
            call h5dread_f(dset_id, h5t_native_double, velth_modext2d, dims2d, err)
         call h5dclose_f(dset_id, err)
         call h5dopen_f(group_id, 'velphi', dset_id, err)
            call h5dread_f(dset_id, h5t_native_double, velphi_modext2d, dims2d, err)
         call h5dclose_f(dset_id, err)
         call h5dopen_f(group_id, 'vth', dset_id, err)
            call h5dread_f(dset_id, h5t_native_double, vth_modext2d, dims2d, err)
         call h5dclose_f(dset_id, err)
      
      call h5gclose_f(group_id, err)
      call h5fclose_f(file_id, err)
      call h5close_f(err)

      !average over all thetas
      do j=1, ntheta_modext
         rho_modext1d = rho_modext1d + rho_modext2d(:,j)
         t_modext1d = t_modext1d + t_modext2d(:,j)
         velr_modext1d = velr_modext1d + velr_modext2d(:,j)
         velth_modext1d = velth_modext1d + velth_modext2d(:,j)
         velphi_modext1d = velphi_modext1d + velphi_modext2d(:,j)
         vth_modext1d = vth_modext1d + vth_modext2d(:,j)
         fsum = fsum + one
      enddo
      !
   enddo
   !
   do j=1, ntheta_modext
      do k=1, nphi_modext
         rho_modext3d(:,j,k) = rho_modext1d/fsum
         t_modext3d(:,j,k) = t_modext1d/fsum
         velr_modext3d(:,j,k) = velr_modext1d/fsum
         velth_modext3d(:,j,k) = velth_modext1d/fsum
         velphi_modext3d(:,j,k) = velphi_modext1d/fsum
         vth_modext3d(:,j,k) = vth_modext1d/fsum
      enddo
   enddo
   !
else
   stop 'error in calc_mod3d_florian: opt_bvel not properly specified'
endif
!
!--------------------------calculate new thermal velocity---------------
!
do i=1, nr_modext
   do j=1, ntheta_modext
      do k=1, nphi_modext
         vth_modext3d(i,j,k) = vthermal(vmicro*1.d5, teff, na)
         eps_cont_modext3d(i,j,k) = eps_cont
         trad_modext3d(i,j,k) = t_modext3d(i,j,k)
      enddo
   enddo
enddo
!
!finally, update theta-array at the boundaries (has to be equal to 0 and 180 to avoid extrapolation in main program)
theta_modext3d(1) = zero
theta_modext3d(ntheta_modext) = pi
write(*,*)
!
!
end subroutine calc_mod3d_florian
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine calc_mod3d_nicowr3d
!
!using tangential spherical cube to map 3D slab simulations
!onto a sphere
!  
!----------------------------------------------------------------------
!
use prog_type
use fund_const
use model3d, only: r_modext3d, theta_modext3d, phi_modext3d, &
                   rho_modext3d, t_modext3d, trad_modext3d, eps_cont_modext3d, vth_modext3d, velr_modext3d, &
                   velth_modext3d, velphi_modext3d
use mod_directories, only: indat_file
use dime_modext, only: nr_modext, ntheta_modext, nphi_modext
use params_input, only: eps_cont, vmicro, na
use params_stellar, only: sr
use mod_amrvac_reader
use mod_interp2d, only: coeff2d_4p_lin
!
implicit none
!
! ... local scalars
integer(i4b) :: nx, ny, nz, nx2, ny2, nz2, is_min, is_max, nx_snaps, ny_snaps, n_snaps, max_refinement
integer(i4b) :: isnap, i, j, k, err, idum
integer(i4b) :: indx_x, indx_y, il, iu, jl, ju
integer(i4b) :: iim2, iim1, ii, iip1, &
                jjm2, jjm1, jj, jjp1
integer(i4b) :: irho, ivelx, ively, ivelz, itgas, itrad, nw
real(dp) :: bconst, zshift
real(dp) :: velr, rad, dx, dy
real(dp) :: phi_min, phi_max, theta_min, theta_max, rmin, rmax
real(dp) :: xp, yp, zp, radp, sint, cost, tant, sinp, cosp, tanp, costt, sintt
real(dp) :: s1, s2, s3, s4, s5, s6, fdum, fdum1, fdum2, fdum3, dels
real(dp) :: acoeff, bcoeff, ccoeff, dcoeff
real(dp) :: rho, tgas, trad, velx, vely, velz, rho_avg, tgas_avg, trad_avg, velr_avg, mdot_avg, mdot
real(dp) :: unit_length, unit_velocity, unit_temperature, unit_density, beta, vmin, vinf
!
! ... local characters
character(len=500)  :: fname_model, fname
!
! ... local logicals
logical :: check1, lxfirst1, lxfirst2, lxfirst3, lxfirst4, lyfirst1, lyfirst2, lyfirst3, lyfirst4
!
! ... namelists
integer(i4b) :: opt_bvel
integer(i4b), dimension(:), allocatable :: nd
namelist / input_usr / fname_model, unit_length, unit_velocity, unit_temperature, unit_density, is_min, is_max, opt_bvel, beta, vmin, vinf, mdot, max_refinement, nd
!
! ... local arrays
real(dp), dimension(:), allocatable :: x_cac3d, y_cac3d, z_cac3d, x_test3d, y_test3d, z_test3d, stretching
real(dp), dimension(:,:,:), allocatable :: rho_cac3d, tgas_cac3d, trad_cac3d, velx_cac3d, vely_cac3d, velz_cac3d
real(dp), dimension(:), allocatable :: x2_cac3d, y2_cac3d, z2_cac3d
real(dp), dimension(:,:,:), allocatable :: rho2_cac3d, tgas2_cac3d, trad2_cac3d, velx2_cac3d, vely2_cac3d, velz2_cac3d
real(dp), dimension(3,3) :: transmat1, transmat2, transmat3, transmat4, transmat5, transmat6, transmat, &
                            transmatt1, transmatt2, transmatt3, transmatt4, transmatt5, transmatt6, transmatt
real(dp), dimension(3) :: vec_cube, vec_plane, eex, eey, eez
!
!
! ... local functions
real(dp) :: vthermal, bvel
!
! ... for random numbers
integer(i4b), dimension(64), parameter :: seed = (/ &
               676865, 350848,  32171, 881468, 615122, 479556, 455219,  89650, &
               835154, 575764, 830269, 479197, 746278, 918314,  62049, 352171, &
               489558, 465378, 903709, 768839, 104036, 606566, 706051, 168587, &
               972615, 980707, 564674, 112480, 559149, 991281, 744474, 619200, &
               684163, 145540,   8131, 896623, 257203, 761081,  44394, 861176, &
               601220, 145688, 366560, 843924, 941293, 376982,   5005, 405884, &
               170702, 848387, 432408, 742485, 239070, 782087, 573526,   7287, &
               430405, 352846, 747140, 620387, 818214, 290087, 335369, 644222 /)
!
! ... local derived types
type(alldata) :: my_data
type(grid) :: my_grid
!
!
allocate(nd(3), stat=err)
!
write(*,*) '---------------------read model from Nicos 3D WR models------------------------'
write(*,*) 'reading input from file: ', trim(indat_file)
write(*,*)
!
open(1, file=trim(indat_file), status='old', form='formatted')
   rewind 1
   read(1, nml=input_usr)
close(1)
!
!always use same seed in order that reproducable results
call random_seed(put=seed)
call random_number(fdum)
idum=nint(is_min + fdum*(is_max-is_min))
if(idum.lt.is_min) stop 'error in calc_mod3d_nicowr3d: random indx lt is_min'
if(idum.gt.is_max) stop 'error in calc_mod3d_nicowr3d: random indx gt is_max'
!
write(fname,'(a,a,i4.4,a)') trim(fname_model), '_', idum, '.dat'
!
fname=adjustl(fname)
inquire(file=trim(fname), exist=check1)
if(.not.check1) then
   write(*,*) 'error in calc_mod3d_nicowr3d: file "', trim(fname), '" does not exist'
   stop
endif
!
!--------------------read dimensions and allocate arrays----------------
!
write(*,*) 'file name: ', trim(fname)
write(*,*)
!
!


allocate(stretching(3), stat=err)
stretching=(/ one, one, one /)
!
!get the data and the grid for this snapshot
my_data = get_data(fname, levmax_usr=max_refinement, stretching=stretching, nd=nd)
my_grid = my_data%data%mesh
!
!
!
nx = my_data%data%data_shape(3)
ny = my_data%data%data_shape(2)
nz = my_data%data%data_shape(1)
nw = my_data%data%data_shape(4)

!write(*,*) nx, ny, nz
!stop 'go on in nico_wr3d'
!
!allocate the data arrays
allocate(x_cac3d(nx), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod3d_nicowr'
allocate(y_cac3d(ny), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod3d_nicowr'
allocate(z_cac3d(nz), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod3d_nicowr'

allocate(x_test3d(nx), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod3d_nicowr'
allocate(y_test3d(ny), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod3d_nicowr'
allocate(z_test3d(nz), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod3d_nicowr'

allocate(rho_cac3d(nx,ny,nz), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod3d_nicowr'
allocate(velx_cac3d(nx,ny,nz), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod3d_nicowr'
allocate(vely_cac3d(nx,ny,nz), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod3d_nicowr'
allocate(velz_cac3d(nx,ny,nz), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod3d_nicowr'
allocate(tgas_cac3d(nx,ny,nz), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod3d_nicowr'
allocate(trad_cac3d(nx,ny,nz), stat=err)
if(err.ne.0) stop 'error: allocation in calc_mod3d_nicowr'
!
!find indices where data is stored
do i=1, nw
!   write(*,*) fdata%w_names(i)
   if(trim(my_data%data%w_names(i)).eq.'rho') irho=i
   if(trim(my_data%data%w_names(i)).eq.'v1') ivelz=i
   if(trim(my_data%data%w_names(i)).eq.'v2') ively=i
   if(trim(my_data%data%w_names(i)).eq.'v3') ivelx=i
   if(trim(my_data%data%w_names(i)).eq.'Trad') itrad=i
   if(trim(my_data%data%w_names(i)).eq.'Tgas') itgas=i
enddo
!
!-------------calculate number of snapshots to be read in---------------
!
x_cac3d = my_data%data%mesh%zgrid%coord
y_cac3d = my_data%data%mesh%ygrid%coord
z_cac3d = my_data%data%mesh%xgrid%coord
!
!
!
dx = x_cac3d(nx)-x_cac3d(1)
dy = y_cac3d(ny)-y_cac3d(1)
!
!number of snapshots required for x and y direction in order to obtain
!complete x-y range from [-1,1]x[-1,1]
!
nx_snaps = nint(2.d0/dx)
ny_snaps = nint(2.d0/dy)
!
!total number of snap-shots required
n_snaps = nx_snaps*ny_snaps
!
if(nx_snaps*dx.ne.2.d0) stop 'error in calc_mod3d_nicowr3d: x-range not matching'
if(ny_snaps*dy.ne.2.d0) stop 'error in calc_mod3d_nicowr3d: y-range not matching'
!
!allocate large cube from [-1,1]x[-1,1]
nx2 = nx_snaps*nx-(nx_snaps-1)
ny2 = ny_snaps*ny-(ny_snaps-1)
nz2 = nz
!
allocate(x2_cac3d(nx2), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod3d_nicowr'
allocate(y2_cac3d(ny2), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod3d_nicowr'
allocate(z2_cac3d(nz2), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod3d_nicowr'
!
allocate(rho2_cac3d(nx2,ny2,nz2), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod3d_nicowr'
allocate(velx2_cac3d(nx2,ny2,nz2), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod3d_nicowr'
allocate(vely2_cac3d(nx2,ny2,nz2), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod3d_nicowr'
allocate(velz2_cac3d(nx2,ny2,nz2), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod3d_nicowr'
allocate(tgas2_cac3d(nx2,ny2,nz2), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod3d_nicowr'
allocate(trad2_cac3d(nx2,ny2,nz2), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_mod3d_nicowr'
!
!-----------------------------------------------------------------------
!
lxfirst1 = .true.
lxfirst2 = .true.
lxfirst3 = .true.
lxfirst4 = .true.

lyfirst1 = .true.
lyfirst2 = .true.
lyfirst3 = .true.
lyfirst4 = .true.
!
!
!
do isnap=1, n_snaps

   call random_number(fdum)
   idum=nint(is_min + fdum*(is_max-is_min))
   if(idum.lt.is_min) stop 'error in calc_mod3d_nicowr3d: random indx lt is_min'
   if(idum.gt.is_max) stop 'error in calc_mod3d_nicowr3d: random indx gt is_max'
!
   write(fname,'(a,a,i4.4,a)') trim(fname_model), '_', idum, '.dat'   
!
   fname=adjustl(fname)
   inquire(file=trim(fname), exist=check1)
   if(.not.check1) then
      write(*,*) 'error in calc_mod3d_nicowr3d: file "', trim(fname), '" does not exist'
      stop
   endif
!
   write(*,*) 'file name: ', trim(fname)
!
!--------------------------read in data--------------------------------
!
   my_data = get_data(fname, levmax_usr=max_refinement, nd=nd, stretching=stretching, grid_out=my_grid)
!   
   x_test3d = my_data%data%mesh%zgrid%coord
   y_test3d = my_data%data%mesh%ygrid%coord
   z_test3d = my_data%data%mesh%xgrid%coord
!
   if(sum(abs(x_test3d-x_cac3d)).gt.1.d-10) stop 'error in calc_mod3d_nicowr3d: x-coordinates not matching'
   if(sum(abs(y_test3d-y_cac3d)).gt.1.d-10) stop 'error in calc_mod3d_nicowr3d: y-coordinates not matching'
   if(sum(abs(z_test3d-z_cac3d)).gt.1.d-10) stop 'error in calc_mod3d_nicowr3d: z-coordinates not matching'   
!
   do i=1, nx
      do j=1, ny
         do k=1, nz
            rho_cac3d(i,j,k) = my_data%data%data(k,j,i,irho)
            velx_cac3d(i,j,k) = my_data%data%data(k,j,i,ivelx)
            vely_cac3d(i,j,k) = my_data%data%data(k,j,i,ively)
            velz_cac3d(i,j,k) = my_data%data%data(k,j,i,ivelz)
            tgas_cac3d(i,j,k) = my_data%data%data(k,j,i,itgas)
            trad_cac3d(i,j,k) = my_data%data%data(k,j,i,itrad)                  
         enddo
      enddo
   enddo
!
!---------------calculate indices where to store data-------------------
!
   call conv_indx_1d_to_2d (isnap, nx_snaps, indx_x, indx_y)
!
!calculate data ranges
   select case(indx_x)
      case(1)
         il = 1
         iu = nx
         if(lxfirst1) then
            x2_cac3d(il:iu) = -one - x_cac3d(1) + x_cac3d
            lxfirst1=.false.
         endif
      case(2)
         il = nx
         iu = 2*nx-1
         if(lxfirst2) then
            x2_cac3d(il:iu) = x2_cac3d(il) - x_cac3d(1) + x_cac3d
            lxfirst2=.false.
         endif
      case(3)
         il = 2*nx-1
         iu = 3*nx-2
         if(lxfirst3) then
            x2_cac3d(il:iu) = x2_cac3d(il) - x_cac3d(1) + x_cac3d
            lxfirst3=.false.
         endif
      case(4)
         il = 3*nx-2
         iu = 4*nx-3
         if(lxfirst4) then
            x2_cac3d(il:iu) = x2_cac3d(il) - x_cac3d(1) + x_cac3d
            lxfirst4=.false.
         endif
      case default
         stop 'error in calc_mod3d_nicowr3d: indx_x out of range'
   end select
!
   select case(indx_y)
      case(1)
         jl = 1
         ju = ny
         if(lyfirst1) then
            y2_cac3d(jl:ju) = -one - y_cac3d(1) + y_cac3d
            lyfirst1=.false.
         endif
      case(2)
         jl = ny
         ju = 2*ny-1
         if(lyfirst2) then
            y2_cac3d(jl:ju) = y2_cac3d(jl) - y_cac3d(1) + y_cac3d
            lyfirst2=.false.
         endif
      case(3)
         jl = 2*ny-1
         ju = 3*ny-2
         if(lyfirst3) then         
            y2_cac3d(jl:ju) = y2_cac3d(jl) - y_cac3d(1) + y_cac3d
            lyfirst3=.false.
         endif
      case(4)
         jl = 3*ny-2
         ju = 4*ny-3
         if(lyfirst4) then
            y2_cac3d(jl:ju) = y2_cac3d(jl) - y_cac3d(1) + y_cac3d
            lyfirst4=.false.
         endif
      case default
         stop 'error in calc_mod3d_nicowr3d: indx_x out of range'
   end select
!
!------------------------store all the data-----------------------------
!
   rho2_cac3d(il:iu,jl:ju,:) = rho_cac3d
   velx2_cac3d(il:iu,jl:ju,:) = velx_cac3d
   vely2_cac3d(il:iu,jl:ju,:) = vely_cac3d
   velz2_cac3d(il:iu,jl:ju,:) = velz_cac3d
   tgas2_cac3d(il:iu,jl:ju,:) = tgas_cac3d
   trad2_cac3d(il:iu,jl:ju,:) = trad_cac3d   
 
enddo
!
!-----------------------deallocate snapshots----------------------------
!
deallocate(rho_cac3d)
deallocate(velx_cac3d)
deallocate(vely_cac3d)
deallocate(velz_cac3d)
deallocate(tgas_cac3d)
deallocate(trad_cac3d)
!
!transform to cgs
rho2_cac3d = rho2_cac3d*unit_density
tgas2_cac3d = tgas2_cac3d*unit_temperature
trad2_cac3d = trad2_cac3d*unit_temperature
velx2_cac3d = velx2_cac3d*unit_velocity
vely2_cac3d = vely2_cac3d*unit_velocity
velz2_cac3d = velz2_cac3d*unit_velocity
!
!----------------create the spherical coordinate grid-------------
!
nr_modext = nz2
rmin = minval(z_cac3d)
rmax = maxval(z_cac3d)
!
ntheta_modext = 2*nx2
theta_min = 0.d0
theta_max = pi
!
nphi_modext = 4*ny2
phi_min = 0.d0
phi_max = 2.d0*pi
!
!
!write(*,*) nx2, ny2, nz2, nr_modext, ntheta_modext, nphi_modext, my_data%data_shape
!write(*,*) 'test1'
!stop 'go on in nico_wr3d'
!
allocate(r_modext3d(nr_modext), stat=err)
allocate(theta_modext3d(ntheta_modext), stat=err)
allocate(phi_modext3d(nphi_modext), stat=err)
!
r_modext3d = z_cac3d
call grid_equi(theta_min, theta_max, ntheta_modext, theta_modext3d)
call grid_equi(phi_min, phi_max, nphi_modext, phi_modext3d)
!
!allocate 3d arrays
allocate(velr_modext3d(nr_modext, ntheta_modext, nphi_modext), stat=err)
allocate(velth_modext3d(nr_modext, ntheta_modext, nphi_modext), stat=err)
allocate(velphi_modext3d(nr_modext, ntheta_modext, nphi_modext), stat=err)
allocate(rho_modext3d(nr_modext, ntheta_modext, nphi_modext), stat=err)
allocate(t_modext3d(nr_modext, ntheta_modext, nphi_modext), stat=err)
allocate(trad_modext3d(nr_modext, ntheta_modext, nphi_modext), stat=err)
allocate(vth_modext3d(nr_modext, ntheta_modext, nphi_modext), stat=err)
allocate(eps_cont_modext3d(nr_modext, ntheta_modext, nphi_modext), stat=err)
!
!---------------define the cube transformations-------------------------
!
!transformation matrix from plane->cube and from cube -> plane for plane 1
transmatt1 = reshape((/ zero, one, zero, &
                       zero, zero, one, &
                       one, zero, zero /), shape(transmat1))
transmat1 = transpose(transmatt1)
!
!transformation matrix from plane->cube and from cube -> plane for plane 2
transmatt2 = reshape((/ -one, zero, zero, &
                       zero, zero, one, &
                       zero, one, zero /), shape(transmat2))
transmat2 = transpose(transmatt2)
!
!transformation matrix from plane->cube and from cube -> plane for plane 3
transmatt3 = reshape((/ zero, -one, zero, &
                       zero, zero, one, &
                       -one, zero, zero /), shape(transmat3))
transmat3 = transpose(transmatt3)
!
!transformation matrix from plane->cube and from cube -> plane for plane 4
transmatt4 = reshape((/ one, zero, zero, &
                       zero, zero, one, &
                       zero, -one, zero /), shape(transmat4))
transmat4 = transpose(transmatt4)
!
!transformation matrix from plane->cube and from cube -> plane for plane 5
transmatt5 = reshape((/ zero, -one, zero, &
                       -one, zero, zero, &
                       zero, zero, one /), shape(transmat5))
transmat5 = transpose(transmatt5)
!
!transformation matrix from plane->cube and from cube -> plane for plane 6
transmatt6 = reshape((/ zero, one, zero, &
                       one, zero, zero, &
                       zero, zero, -one /), shape(transmat6))
transmat6 = transpose(transmatt6)
!
!-----------------------------------------------------------------------
!
write(*,*)
write(*,*) 'transforming data from cube to sphere'
write(*,*)

!vec_cube = (/ 1., 2., 3. /)
!write(*,*) matmul(transmatt1,vec_cube)
do j=1, ntheta_modext!ntheta_modext/4, ntheta_modext
!do j=200, 200
   sint = sin(theta_modext3d(j))
   cost = cos(theta_modext3d(j))
   tant = tan(theta_modext3d(j))

   costt = cos(pi/two-theta_modext3d(j))
   sintt = sin(pi/two-theta_modext3d(j))

   do k=1, nphi_modext
!   do k=85, 85
      sinp = sin(phi_modext3d(k))
      cosp = cos(phi_modext3d(k))
      tanp = tan(phi_modext3d(k))
!
!calculate coordinate on the cube      
!decide which plane is to be used by calculating distance to each plane
      fdum1 = sint*cosp
      if(fdum1.eq.zero) then
         s1 = 1.d10
         s3 = 1.d10
      elseif(fdum1.lt.zero) then
         s1 = 1.d10
         s3 = -one/fdum1
      else
         s1 = one/fdum1
         s3 = 1.d10
      endif

      fdum2 = sint*sinp
      if(fdum2.eq.zero) then
         s2 = 1.d10
         s4 = 1.d10
      elseif(fdum2.lt.zero) then
         s2 = 1.d10
         s4 = -one/fdum2
      else
         s2 = one/fdum2
         s4 = 1.d10
      endif 

      fdum3 = cost
      if(fdum3.eq.zero) then
         s5 = 1.d10
         s6 = 1.d10
      elseif(fdum3.lt.zero) then
         s5 = 1.d10
         s6 = -one/fdum3
      else
         s5 = one/fdum3
         s6 = 1.d10
      endif

      dels = min(s1,s2,s3,s4,s5,s6)
!
!calculate position on the cube and set transformation-matrices
      if(dels.eq.s1) then
         vec_cube = (/ one, tanp, one/tant/cosp /)
         transmat = transmat1
         transmatt = transmatt1
      elseif(dels.eq.s2) then
         vec_cube = (/ one/tanp, one, one/tant/sinp /)
         transmat = transmat2
         transmatt = transmatt2
      elseif(dels.eq.s3) then
         vec_cube = (/ -one, -tanp, -one/tant/cosp /)
         transmat = transmat3
         transmatt = transmatt3
      elseif(dels.eq.s4) then
         vec_cube = (/ -one/tanp, -one, -one/tant/sinp /)
         transmat = transmat4
         transmatt = transmatt4
      elseif(dels.eq.s5) then
         vec_cube = (/ tant*cosp, tant*sinp, one /)
         transmat = transmat5
         transmatt = transmatt5
      elseif(dels.eq.s6) then
         vec_cube = (/ -tant*cosp, -tant*sinp, -one /)
         transmat = transmat6
         transmatt = transmatt6
      else
         stop 'error in calc_mod3d_nicowr3d: dels not valid'
      endif
!calculate coordinates on the slab
      vec_plane = matmul(transmat,vec_cube)
!calculate interpolation coefficients
      call find_index(vec_plane(1), x2_cac3d, nx2, iim2, iim1, ii, iip1)
      call find_index(vec_plane(2), y2_cac3d, ny2, jjm2, jjm1, jj, jjp1)
      call coeff2d_4p_lin(x2_cac3d(iim1), x2_cac3d(ii), y2_cac3d(jjm1), y2_cac3d(jj), &
                          vec_plane(1), vec_plane(2), &
                          acoeff, bcoeff, ccoeff, dcoeff)         
!
!perform interpolation on each radial point
      do i=1, nz
          rho = acoeff*rho2_cac3d(iim1,jjm1,i) + bcoeff*rho2_cac3d(ii,jjm1,i) + &
                ccoeff*rho2_cac3d(iim1,jj,i)   + dcoeff*rho2_cac3d(ii,jj,i)
          tgas = acoeff*tgas2_cac3d(iim1,jjm1,i) + bcoeff*tgas2_cac3d(ii,jjm1,i) + &
                 ccoeff*tgas2_cac3d(iim1,jj,i)   + dcoeff*tgas2_cac3d(ii,jj,i)
          trad = acoeff*trad2_cac3d(iim1,jjm1,i) + bcoeff*trad2_cac3d(ii,jjm1,i) + &
                ccoeff*trad2_cac3d(iim1,jj,i)   + dcoeff*trad2_cac3d(ii,jj,i)          
          trad = acoeff*trad2_cac3d(iim1,jjm1,i) + bcoeff*trad2_cac3d(ii,jjm1,i) + &
                ccoeff*trad2_cac3d(iim1,jj,i)   + dcoeff*trad2_cac3d(ii,jj,i)
          velx = acoeff*velx2_cac3d(iim1,jjm1,i) + bcoeff*velx2_cac3d(ii,jjm1,i) + &
                ccoeff*velx2_cac3d(iim1,jj,i)   + dcoeff*velx2_cac3d(ii,jj,i)
          vely = acoeff*vely2_cac3d(iim1,jjm1,i) + bcoeff*vely2_cac3d(ii,jjm1,i) + &
                ccoeff*vely2_cac3d(iim1,jj,i)   + dcoeff*vely2_cac3d(ii,jj,i)
          velz = acoeff*velz2_cac3d(iim1,jjm1,i) + bcoeff*velz2_cac3d(ii,jjm1,i) + &
               ccoeff*velz2_cac3d(iim1,jj,i)   + dcoeff*velz2_cac3d(ii,jj,i)
!
!transform velocity vector on slab to velocity vector on cube
          vec_plane = (/ velx, vely, velz /)
!          vec_cube = matmul(transmatt,vec_plane)
!          velx = vec_cube(1)*eex(1) + vec_cube(2)*eey(1) + vec_cube(3)*eez(1)
!          vely = vec_cube(1)*eex(2) + vec_cube(2)*eey(2) + vec_cube(3)*eez(2)
!          velz = vec_cube(1)*eex(3) + vec_cube(2)*eey(3) + vec_cube(3)*eez(3)
!
!store everything in spherical coordinates
          rho_modext3d(i,j,k) = rho
          t_modext3d(i,j,k) = tgas
          trad_modext3d(i,j,k) = trad
          velr_modext3d(i,j,k) = velz !sint*cosp*velx + sint*sinp*vely + cost*velz
          velth_modext3d(i,j,k) = velx*cosp + vely*sinp !cost*cosp*velx + cost*sinp*vely - sint*velz
          velphi_modext3d(i,j,k) = -velx*sinp + vely*cosp !-sinp*velx + cosp*vely
!          write(*,*) velx, vely, velz
          vth_modext3d(i,j,k) = vthermal(vmicro*1.d5,tgas,na)
          eps_cont_modext3d(i,j,k) = eps_cont

!          write(*,*) asin(sint)*180./pi, asin(sinp)*180./pi , velx/1.d5, vely/1.d5, velz/1.d5, velr_modext3d(i,j,k)/1.d5
       enddo
   enddo
enddo

!stop 'go on in wr3d'
!
!
!
!-----overwrite everything with beta-velocity-law if option set---------
!
if(opt_bvel.eq.1.or.opt_bvel.eq.2) then
   write(*,*) 'calculating beta-velocity law'
   write(*,*)
!
!calculate lateral and azimuthal average velocities, average temperatures, average densities
   do i=1, nr_modext
      trad_avg = 0.d0
      tgas_avg = 0.d0
      rho_avg = 0.d0
      velr_avg = 0.d0
!
      idum = 0
      do j=1, ntheta_modext
         do k=1, nphi_modext
            trad_avg = trad_avg + trad_modext3d(i,j,k)
            tgas_avg = tgas_avg + t_modext3d(i,j,k)
            velr_avg = velr_avg + velr_modext3d(i,j,k)
            rho_avg = rho_avg + rho_modext3d(i,j,k)
            idum = idum+1
         enddo
      enddo


      rho_modext3d(i,:,:) = rho_avg/idum
      t_modext3d(i,:,:) = tgas_avg/idum
      trad_modext3d(i,:,:) = trad_avg/idum
      velr_modext3d(i,:,:) = velr_avg/idum
      velth_modext3d(i,:,:) = zero
      velphi_modext3d(i,:,:) = zero
      vth_modext3d(i,:,:) = vthermal(vmicro*1.d5,tgas_avg,na)
      eps_cont_modext3d(i,:,:) = eps_cont
   enddo
!
!calculate radial average radial velociies and mass-loss rates (only outermost 40 grid points)
   idum=0
   mdot_avg = 0.d0
   velr_avg = 0.d0
   do i=nr_modext-40, nr_modext
      velr_avg = velr_avg + velr_modext3d(i,1,1)
      mdot_avg = mdot_avg + four*pi*(r_modext3d(i)*sr)**2 * rho_modext3d(i,1,1) * velr_modext3d(i,1,1)
      idum=idum+1
   enddo
   velr_avg = velr_avg/idum
   mdot_avg = mdot_avg/idum

!   beta=1.d0
!   vmin=10.d5

   if(opt_bvel.eq.1) then
      !calculate v_inf such that vel_r(beta-law) = vel_r(model) at r_max
      vmin = vmin*1.d5 !in cm/s      
      fdum = (vmin/velr_avg)**(one/beta)
      bconst = (one - fdum)/(one - fdum/r_modext3d(nr_modext))
      vinf = vmin/(one-bconst)**beta
      mdot = mdot_avg
   elseif(opt_bvel.eq.2) then
      vmin=vmin*1.d5 !in cm/s
      vinf=vinf*1.d5 !in cm/s
      mdot=mdot*xmsu/yr !in g/s
      bconst = one - (vmin/vinf)**(one/beta)
   endif

   write(*,*) 'using'
   write(*,*) 'beta', beta
   write(*,*) 'vmin [km/s]', vmin/1.d5
   write(*,*) 'v_r(r_max) [km/s]', velr_avg/1.d5   
   write(*,*) 'vinf [km/s]', vinf/1.d5
   write(*,*) 'mdot [Msun/yr]', mdot*yr/xmsu
   write(*,*) 'bconst', bconst   
   
   do i=1, nr_modext
      velr_modext3d(i,:,:) = bvel(r_modext3d(i), vinf, bconst, beta)
      rho_modext3d(i,:,:) = mdot/four/pi/velr_modext3d(i,1,1)/(r_modext3d(i)*sr)**2
   enddo
   

   

endif
!stop 'go on in nico3d'
!
!
end subroutine calc_mod3d_nicowr3d
