!
!-----------------------------------------------------------------------
!
!   calculates or reads 1d/2d/3d model atmospheres that can be read in
!               spec_vbin.eo to calculate line profiles
!
!v02: including jet from dylan bollen
!     updated indat files
!     including coordinate bases in global system (allowing for tilts of complete coordinate systems)
!     including rotation axes of individual systems in global system (allowing for stellar rotations)
!     including sunyaev shakura disc
!     including jet from olivier verhamme
!-----------------------------------------------------------------------
!
program modspec
!
use prog_type
use options_modspec, only: input_mod
!
implicit none
!
! ... local scalars
!
call read_input
!
select case(input_mod)
   case(0)
      call calc_model3da
      call output3d_spc
   case(1)
!Star 1 with Sobolev solution
      call calc_model3db
      call output3d_spc
   case(2)
!Star 2 with Sobolev solution
      call calc_model3dc
      call output3d_spc
   case(3)
!Star 1 with Be-disc from Dylans initial condition
      call calc_model3dd
      call output3d_spc
   case(4)
!Star 1 with Be-disc from Ileyks hydro simulations
      call calc_model3de
      call output3d_spc
   case(5)
!Star 1 and star 2 with disc from Dylan
      call calc_model3df
      call output3d_spc
   case(6)
!Star 1 with disc from dylan, star 2 with wind
      call calc_model3dg
      call output3d_spc 
   case(7)
!Star 1 as post AGB star, star 2 with jet
      call calc_model3d_jeta
      call output3d_spc
   case(8)
!Star 1, and star 2 with alpha-disc
      call calc_model3d_adisc
      call output3d_spc      
   case(9)
!Star 1 as post AGB star, star 2 with jet from Oliviers model
      call calc_model3d_jetb
      call output3d_spc      
   case(1000)
!only for some testing
      call calc_model3d_test
   case default
      stop 'error in modspec: input_mod not specified'
end select
!

end program modspec
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine read_input
!
!-----------------------------------------------------------------------
!-------------------read in all input parameters------------------------
!-----------------------------------------------------------------------
!
use prog_type
use fund_const, only: rsu, xlsu, four, pi, cgs_sb, one
use options_modspec, only: indat_file, input_file, input_file2, output_file, input_mod
use params_modspec, only: x01, y01, z01, vx01, vy01, vz01, rstar1, rmin1, rmax1, sr1, teff1, trad1, &
                          logg1, lstar1, yhe1, yhe1_mass, fehe1, aenh1, vrot1, vmicro1, &
                          ex01, ey01, ez01, rot_axis01, &
                          x02, y02, z02, vx02, vy02, vz02, rstar2, rmin2, rmax2, sr2, teff2, trad2, &
                          logg2, lstar2, yhe2, yhe2_mass, fehe2, aenh2, vrot2, vmicro2, &
                          ex02, ey02, ez02, rot_axis02, &
                          eps_line, unit_length, vth_fiducial
use mod_iline, only: iline, na, xnue0, kline, get_iline
use mod_lte, only: get_lte_table
!
!
implicit none
!
! ... local scalars
real(dp) :: fcheck1, fcheck2, fcheck3
!
! ... local characters
!
! ... local functions
real(dp), dimension(3) :: p_object01, v_object01, p_object02, v_object02
!
! ... namelist
namelist / input_options / input_file, input_file2, output_file, input_mod
namelist / input_model1 / rstar1, rmin1, rmax1, teff1, trad1, logg1, yhe1, fehe1, aenh1, vrot1, vmicro1, p_object01, v_object01, ex01, ey01, ez01, rot_axis01
namelist / input_model2 / rstar2, rmin2, rmax2, teff2, trad2, logg2, yhe2, fehe2, aenh2, vrot2, vmicro2, p_object02, v_object02, ex02, ey02, ez02, rot_axis02
namelist / input_units / unit_length, vth_fiducial
namelist / input_line / iline, eps_line, kline
!
!-----------------------------------------------------------------------
!
write(*,*) '----------------------------read input-----------------------------------------'
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

!read model parameters of star 1
   rewind 1
   read(1, nml=input_model1)
   sr1=rstar1*rsu
   lstar1 = four*pi*cgs_sb * sr1**2 * teff1**4 / xlsu
   !
!positions as scalars
   x01 = p_object01(1)
   y01 = p_object01(2)
   z01 = p_object01(3)   
!velocities as scalars and in cm/s
   vx01 = v_object01(1)*1.d5
   vy01 = v_object01(2)*1.d5
   vz01 = v_object01(3)*1.d5  
   vrot1 = vrot1*1.d5
   vmicro1 = vmicro1*1.d5
!
!normalize vectors describing the coordinate system
   ex01 = ex01/sqrt(ex01(1)**2 + ex01(2)**2 + ex01(3)**2)
   ey01 = ey01/sqrt(ey01(1)**2 + ey01(2)**2 + ey01(3)**2)      
   ez01 = ez01/sqrt(ez01(1)**2 + ez01(2)**2 + ez01(3)**2)

   fcheck1 = dot_product(ex01,ey01)
   fcheck2 = dot_product(ex01,ez01)
   fcheck3 = dot_product(ey01,ez01)
   if(abs(fcheck1).gt.1.d-14) then
      write(*,*) 'error in read_input: ex01, ey01 not perpendicular'
      write(*,*) 'ex01*ey01', fcheck1
      stop
   endif
   if(abs(fcheck2).gt.1.d-14) then
      write(*,*) 'error in read_input: ex01, ez01 not perpendicular'
      write(*,*) 'ex01*ez01', fcheck2
      stop
   endif
   if(abs(fcheck3).gt.1.d-14) then
      write(*,*) 'error in read_input: ey01, ez01 not perpendicular'
      write(*,*) 'ey01*ez01', fcheck3
      stop
   endif   

!normalize rotation axis of star 1
   rot_axis01 = rot_axis01/sqrt(rot_axis01(1)**2 + rot_axis01(2)**2 + rot_axis01(3)**2)
!
!
!read model parameters of star 2
   rewind 1
   read(1, nml=input_model2)
   sr2=rstar2*rsu
   lstar2 = four*pi*cgs_sb * sr2**2 * teff2**4 / xlsu   

!positions as scalars
   x02 = p_object02(1)
   y02 = p_object02(2)
   z02 = p_object02(3)   
!velocities as scalars and in cm/s
   vx02 = v_object02(1)*1.d5
   vy02 = v_object02(2)*1.d5
   vz02 = v_object02(3)*1.d5  
   vrot2 = vrot2*1.d5
   vmicro2 = vmicro2*1.d5
!
!normalize vectors describing the coordinate system
   ex02 = ex02/sqrt(ex02(1)**2 + ex02(2)**2 + ex02(3)**2)
   ey02 = ey02/sqrt(ey02(1)**2 + ey02(2)**2 + ey02(3)**2)      
   ez02 = ez02/sqrt(ez02(1)**2 + ez02(2)**2 + ez02(3)**2)

   fcheck1 = dot_product(ex02,ey02)
   fcheck2 = dot_product(ex02,ez02)
   fcheck3 = dot_product(ey02,ez02)      
   if(fcheck1.ne.0.) stop 'error in read_input: ex02, ey02 not perpendicular'
   if(fcheck2.ne.0.) stop 'error in read_input: ex02, ez02 not perpendicular'
   if(fcheck3.ne.0.) stop 'error in read_input: ey02, ez02 not perpendicular'
   
!normalize rotation axis of star 2
   rot_axis02 = rot_axis02/sqrt(rot_axis02(1)**2 + rot_axis02(2)**2 + rot_axis02(3)**2)
!
!read line strength parameters etc
   rewind 1
   read(1, nml=input_line)
!
!read units
   rewind 1
   read(1, nml=input_units)
   vth_fiducial=vth_fiducial*1.d5
!
close(1)
!
!
call get_iline(iline)
!
!
!
end subroutine read_input
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine calc_model3da
!
!blank two coordinate system model
!
use prog_type
use fund_const, only: zero, two, pi
use options_modspec, only: input_file, input_file2, input_mod
use dime_modspec, only: cs1_nr, cs1_ntheta, cs1_nphi, cs1_r, cs1_theta, cs1_phi, &
                        cs1_velx3d, cs1_vely3d, cs1_velz3d, cs1_opac3d, &
                        cs1_opalbar3d, cs1_sline3d, cs1_scont3d, cs1_imask3d, cs1_t3d, cs1_rho3d, &
                        cs2_nr, cs2_ntheta, cs2_nphi, cs2_r, cs2_theta, cs2_phi, &
                        cs2_velx3d, cs2_vely3d, cs2_velz3d, cs2_opac3d, &
                        cs2_opalbar3d, cs2_sline3d, cs2_scont3d, cs2_imask3d, cs2_t3d, cs2_rho3d
use params_modspec, only: vx01, vy01, vz01, vx02, vy02, vz02, x01, y01, z01, x02, y02, z02, &
                          vmin1, vmax1, beta1, vrot1, rstar1, rmin1, rmax1, sr1, teff1, yhe1, hei1, logg1, lstar1, mdot1, &
                          vmin2, vmax2, beta2, vrot2, rstar2, rmin2, rmax2, sr2, teff2, yhe2, hei2, logg2, lstar2, mdot2, &
                          vmax, vth_fiducial, vmicro1, vmicro2
use hdf5
use mod_grid, only: grid_log, grid_equi
!
implicit none
!
! ... local scalars
integer(i4b) :: i, j, k, err
real(dp) :: vel, opac, opalbar, rho, sline, scont
!
! ... local arrays
!
! ... local characters
!
! ... local functions
!
! ... for hdf5 file
!
!definitions for the line-transition
!velocity in cm/s (note that we are allowed to change vth_fiducial in spec_vbin.eo)
vth_fiducial=1.d8
!
write(*,*) '----------------calc_model3da: calculating 3d model for star 1-----------------'
write(*,*)
!
!need to be defined (if gravity darkening and/or photospheric profile shall be accounted for in spec_vbin.eo)
vmicro1=1.d8
logg1=3.6d0
lstar1=1.d6
vrot1=0.d8
yhe1=0.1d0
!
!
cs1_nr=51
cs1_ntheta=51
cs1_nphi=2*cs1_ntheta-1
!
!allocate arrays
allocate(cs1_r(cs1_nr), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3da'
allocate(cs1_theta(cs1_ntheta), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3da'
allocate(cs1_phi(cs1_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3da'
allocate(cs1_opac3d(cs1_nr,cs1_ntheta,cs1_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3da'
allocate(cs1_opalbar3d(cs1_nr,cs1_ntheta,cs1_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3da'
allocate(cs1_scont3d(cs1_nr,cs1_ntheta,cs1_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3da'
allocate(cs1_sline3d(cs1_nr,cs1_ntheta,cs1_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3da'
allocate(cs1_velx3d(cs1_nr,cs1_ntheta,cs1_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3da'
allocate(cs1_vely3d(cs1_nr,cs1_ntheta,cs1_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3da'
allocate(cs1_velz3d(cs1_nr,cs1_ntheta,cs1_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3da'
allocate(cs1_t3d(cs1_nr,cs1_ntheta,cs1_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3da'
allocate(cs1_rho3d(cs1_nr,cs1_ntheta,cs1_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3da'
!
!radius in units of sr1   
call grid_log(rmin1, rmax1, cs1_nr,  cs1_r)
!
!theta-grid equidistant
call grid_equi(zero, pi, cs1_ntheta, cs1_theta)
!
!phi-grid equidistant
call grid_equi(zero, two*pi, cs1_nphi, cs1_phi)
!
!
!
do i=1, cs1_nr
   do j=1, cs1_ntheta
      do k=1, cs1_nphi
!        
!----------------------line source function-----------------------------
!
         sline = zero
         cs1_sline3d(i,j,k)=sline
!
!----------------------continuum source function------------------------
!
         scont = zero
         cs1_scont3d(i,j,k)=scont
!
!------------------------velocity components----------------------------
!
         cs1_velx3d(i,j,k)=zero
         cs1_vely3d(i,j,k)=zero
         cs1_velz3d(i,j,k)=zero
!
!---------------------------temperature---------------------------------
!
         cs1_t3d(i,j,k) = teff1
!
!-----------------------------opacity-----------------------------------
!
         rho = zero
         opalbar=zero
         opac=zero
!
         cs1_rho3d(i,j,k) = rho
         cs1_opalbar3d(i,j,k) = opalbar
         cs1_opac3d(i,j,k) = opac
!
      enddo
   enddo
enddo
!
!
!
write(*,*) '----------------calc_model3da: calculating 3d model for star 2-----------------'
write(*,*)
!
!need to be defined (if gravity darkening and/or photospheric profile shall be accounted for in spec_vbin.eo)
vmicro2=1.d8
logg2=3.6d0
lstar2=1.d6
vrot2=0.d8
yhe2=0.1d0
!
cs2_nr=31
cs2_ntheta=31
cs2_nphi=2*cs2_ntheta-1
!
!allocate arrays
allocate(cs2_r(cs2_nr), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3da'
allocate(cs2_theta(cs2_ntheta), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3da'
allocate(cs2_phi(cs2_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3da'
allocate(cs2_opac3d(cs2_nr,cs2_ntheta,cs2_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3da'
allocate(cs2_opalbar3d(cs2_nr,cs2_ntheta,cs2_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3da'
allocate(cs2_scont3d(cs2_nr,cs2_ntheta,cs2_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3da'
allocate(cs2_sline3d(cs2_nr,cs2_ntheta,cs2_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3da'
allocate(cs2_velx3d(cs2_nr,cs2_ntheta,cs2_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3da'
allocate(cs2_vely3d(cs2_nr,cs2_ntheta,cs2_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3da'
allocate(cs2_velz3d(cs2_nr,cs2_ntheta,cs2_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3da'
allocate(cs2_t3d(cs2_nr,cs2_ntheta,cs2_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3da'
allocate(cs2_rho3d(cs2_nr,cs2_ntheta,cs2_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3da'
!
!radius in units of sr1   
call grid_log(rmin2, rmax2, cs2_nr,  cs2_r)
!
!theta-grid equidistant
call grid_equi(zero, pi, cs2_ntheta, cs2_theta)
!
!phi-grid equidistant
call grid_equi(zero, two*pi, cs2_nphi, cs2_phi)
!
!
!
do i=1, cs2_nr
   do j=1, cs2_ntheta
      do k=1, cs2_nphi
!        
!----------------------line source function-----------------------------
!
         sline = zero
         cs2_sline3d(i,j,k)=sline
!
!----------------------continuum source function------------------------
!
         scont = zero
         cs2_scont3d(i,j,k)=scont
!
!------------------------velocity components----------------------------
!
         cs2_velx3d(i,j,k)=zero
         cs2_vely3d(i,j,k)=zero
         cs2_velz3d(i,j,k)=zero
!
!---------------------------temperature---------------------------------
!
         cs2_t3d(i,j,k) = teff2
!
!-----------------------------opacity-----------------------------------
!
         rho = zero
         opalbar=zero
         opac=zero
!
         cs2_rho3d(i,j,k) = rho
         cs2_opalbar3d(i,j,k) = opalbar
         cs2_opac3d(i,j,k) = opac
!
      enddo
   enddo
enddo
!
!
!-----------------------------set velocities------------------------------
!
!maximum velocity in global system
!
vmax=zero
do i=1, cs1_nr
   do j=1, cs1_ntheta
      do k=1, cs1_nphi
         vel=sqrt((vx01+cs1_velx3d(i,j,k))**2+(vy01+cs1_vely3d(i,j,k))**2+(vz01+cs1_velz3d(i,j,k))**2)
         if(vel.gt.vmax) vmax=vel
      enddo
   enddo
enddo

do i=1, cs2_nr
   do j=1, cs2_ntheta
      do k=1, cs2_nphi
         vel=sqrt((vx02+cs2_velx3d(i,j,k))**2+(vy02+cs2_vely3d(i,j,k))**2+(vz02+cs2_velz3d(i,j,k))**2)
         if(vel.gt.vmax) vmax=vel
      enddo
   enddo
enddo
!
!
!
!
end subroutine calc_model3da
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine calc_model3db
!
!
use prog_type
use fund_const, only: one, zero, two, pi, xmsu
use options_modspec, only: input_file, input_file2, input_mod
use dime_modspec, only: cs1_nr, cs1_ntheta, cs1_nphi, cs1_r, cs1_theta, cs1_phi, &
                        cs1_velx3d, cs1_vely3d, cs1_velz3d, cs1_opac3d, &
                        cs1_opalbar3d, cs1_sline3d, cs1_scont3d, cs1_imask3d, cs1_t3d, cs1_rho3d, &
                        cs2_nr, cs2_ntheta, cs2_nphi, cs2_r, cs2_theta, cs2_phi, &
                        cs2_velx3d, cs2_vely3d, cs2_velz3d, cs2_opac3d, &
                        cs2_opalbar3d, cs2_sline3d, cs2_scont3d, cs2_imask3d, cs2_t3d, cs2_rho3d
use params_modspec, only: vx01, vy01, vz01, vx02, vy02, vz02, x01, y01, z01, x02, y02, z02, &
                          vmin1, vmax1, beta1, vrot1, rstar1, rmin1, rmax1, sr1, teff1, yhe1, hei1, logg1, lstar1, mdot1, vmicro1, &
                          vmin2, vmax2, beta2, vrot2, rstar2, rmin2, rmax2, sr2, teff2, yhe2, hei2, logg2, lstar2, mdot2, vmicro2, &
                          vmax, vth_fiducial, trad1, unit_length
use hdf5
use mod_grid, only: grid_log, grid_equi
use mod_opacities, only: get_opalbar
use mod_iline, only: iline, na, xnue0, kline
!
implicit none
!
! ... local scalars
integer(i4b) :: i, j, k, err
real(dp) :: vel, opac, opalbar, rho, sline, scont
real(dp) :: beta, vmin, vinf, mdot, mdot_cgs, bconst, xic1, eps_line, gradv, temp, velr, vth
real(dp) :: sint, cost, sinp, cosp, theta1, phi1, xcoord0, ycoord0, zcoord0, xcoord1, ycoord1, zcoord1, &
            velx0, vely0, velz0, velx1, vely1, velz1, rad, rad1, rlim1
!
! ... local arrays
!
! ... local characters
!
! ... local functions
real(dp) :: bnue, bvel, vthermal, sobo1d
!
! ... for hdf5 file
!
write(*,*) '----------------calc_model3db 3d Sobolev beta-vel-wind for star 1-------------'
write(*,*)
!
!define beta velocity wind
beta=1.d0
vmin=10.d5
vinf=2000.d5
mdot = 1.d-6
mdot_cgs=mdot*xmsu/(365.25d0*24.d0*3600.d0)
xic1=bnue(xnue0, trad1)
bconst=1.d0-(vmin/vinf)**(1.d0/beta)
!
vmicro1=1.d7
vth_fiducial=1.d7
logg1=3.6d0
lstar1=1.d6
vrot1=0.d8
yhe1=0.1d0
hei1=2.d0
eps_line=1.d-6
kline=1.d3
rlim1=11.d0
!
!
!
cs1_nr=101
cs1_ntheta=51
cs1_nphi=2*cs1_ntheta-1
!
!allocate arrays
allocate(cs1_r(cs1_nr), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3db'
allocate(cs1_theta(cs1_ntheta), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3db'
allocate(cs1_phi(cs1_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3db'
allocate(cs1_opac3d(cs1_nr,cs1_ntheta,cs1_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3db'
allocate(cs1_opalbar3d(cs1_nr,cs1_ntheta,cs1_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3db'
allocate(cs1_scont3d(cs1_nr,cs1_ntheta,cs1_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3db'
allocate(cs1_sline3d(cs1_nr,cs1_ntheta,cs1_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3db'
allocate(cs1_velx3d(cs1_nr,cs1_ntheta,cs1_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3db'
allocate(cs1_vely3d(cs1_nr,cs1_ntheta,cs1_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3db'
allocate(cs1_velz3d(cs1_nr,cs1_ntheta,cs1_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3db'
allocate(cs1_t3d(cs1_nr,cs1_ntheta,cs1_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3db'
allocate(cs1_rho3d(cs1_nr,cs1_ntheta,cs1_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3db'
!
!radius in units of sr1   
call grid_log(rmin1, rmax1, cs1_nr,  cs1_r)
!
!theta-grid equidistant
call grid_equi(zero, pi, cs1_ntheta, cs1_theta)
!
!phi-grid equidistant
call grid_equi(zero, two*pi, cs1_nphi, cs1_phi)
!
!
if(cs1_r(cs1_nr-1).le.rlim1) then
   write(*,*) 'error: we require always two zero-opacity ghost points at outer boundary'
   stop
endif
!
cs1_sline3d=zero
cs1_scont3d=zero
cs1_velx3d=zero
cs1_vely3d=zero
cs1_velz3d=zero
cs1_t3d=zero
cs1_rho3d=zero
cs1_opalbar3d=zero
cs1_opac3d=zero
!
do i=1, cs1_nr
   velr = bvel(cs1_r(i), vinf, bconst, beta)
   temp = trad1
   vth = vthermal(vmicro1, temp, na)
   rho = mdot_cgs/(4.d0*pi*cs1_r(i)**2 * sr1**2 * velr)
   opalbar = get_opalbar(iline, kline, sr1, yhe1, hei1, temp, vth_fiducial, xnue0, one, one, rho)/sr1  !in cgs
   gradv = velr * bconst*beta/cs1_r(i)**2/(1.d0-bconst/cs1_r(i))
   sline = sobo1d(cs1_r(i), velr/vth, gradv/vth, opalbar*sr1, temp, xic1, xnue0, eps_line)
   scont = zero
   opac = zero
   
   do j=1, cs1_ntheta
      do k=1, cs1_nphi
!        
!----------------------line source function-----------------------------
!         
         if(cs1_r(i).le.rlim1) cs1_sline3d(i,j,k)=sline
!
!----------------------continuum source function------------------------
!
         if(cs1_r(i).le.rlim1) cs1_scont3d(i,j,k)=scont
!
!------------------------velocity components----------------------------
!
         sint=sin(cs1_theta(j))
         cost=cos(cs1_theta(j))
         sinp=sin(cs1_phi(k))
         cosp=cos(cs1_phi(k))

         cs1_velx3d(i,j,k)=velr*sint*cosp
         cs1_vely3d(i,j,k)=velr*sint*sinp
         cs1_velz3d(i,j,k)=velr*cost
!
!---------------------------temperature---------------------------------
!
         cs1_t3d(i,j,k) = temp
!
!-----------------------------opacity-----------------------------------
!
         if(cs1_r(i).le.rlim1) cs1_rho3d(i,j,k) = rho
         if(cs1_r(i).le.rlim1) cs1_opalbar3d(i,j,k) = opalbar*sr1   !in units of sr1
         if(cs1_r(i).le.rlim1) cs1_opac3d(i,j,k) = opac*sr1
!
      enddo
   enddo
enddo
!
!write(*,*) iline, kline, rstar1, sr1, yhe1, hei1, trad1, vth_fiducial, xnue0
!stop 'go on in sobo_vbin'
!
!
write(*,*) '----------------calc_model3db: calculating 3d model for star 2-----------------'
write(*,*)
!
!need to be defined (if gravity darkening and/or photospheric profile shall be accounted for in spec_vbin.eo)
vmicro2=1.d7
logg2=3.6d0
lstar2=1.d6
vrot2=0.d7
yhe2=0.1d0
!
cs2_nr=31
cs2_ntheta=31
cs2_nphi=2*cs2_ntheta-1
!
!allocate arrays
allocate(cs2_r(cs2_nr), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3db'
allocate(cs2_theta(cs2_ntheta), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3db'
allocate(cs2_phi(cs2_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3db'
allocate(cs2_opac3d(cs2_nr,cs2_ntheta,cs2_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3db'
allocate(cs2_opalbar3d(cs2_nr,cs2_ntheta,cs2_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3db'
allocate(cs2_scont3d(cs2_nr,cs2_ntheta,cs2_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3db'
allocate(cs2_sline3d(cs2_nr,cs2_ntheta,cs2_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3db'
allocate(cs2_velx3d(cs2_nr,cs2_ntheta,cs2_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3db'
allocate(cs2_vely3d(cs2_nr,cs2_ntheta,cs2_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3db'
allocate(cs2_velz3d(cs2_nr,cs2_ntheta,cs2_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3db'
allocate(cs2_t3d(cs2_nr,cs2_ntheta,cs2_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3db'
allocate(cs2_rho3d(cs2_nr,cs2_ntheta,cs2_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3db'
!
!radius in units of sr1   
call grid_log(rmin2, rmax2, cs2_nr,  cs2_r)
!
!theta-grid equidistant
call grid_equi(zero, pi, cs2_ntheta, cs2_theta)
!
!phi-grid equidistant
call grid_equi(zero, two*pi, cs2_nphi, cs2_phi)
!
!
!
do i=1, cs2_nr
   do j=1, cs2_ntheta
      do k=1, cs2_nphi
!
         velr=zero
         temp=zero
         
         rad=cs2_r(i)         
         sint=sin(cs2_theta(j))
         cost=cos(cs2_theta(j))
         sinp=sin(cs2_phi(k))
         cosp=cos(cs2_phi(k))

!position in global coordinate system
         xcoord0 = x02 + rad*sint*cosp*rstar2/unit_length
         ycoord0 = y02 + rad*sint*sinp*rstar2/unit_length
         zcoord0 = z02 + rad*cost*rstar2/unit_length
!postion in system of star 1
         xcoord1 = (xcoord0-x01)*unit_length/rstar1
         ycoord1 = (ycoord0-y01)*unit_length/rstar1
         zcoord1 = (zcoord0-z01)*unit_length/rstar1
         rad1=sqrt(xcoord1**2+ycoord1**2+zcoord1**2)
         call get_angles_spc(xcoord1, ycoord1, zcoord1, theta1, phi1)         
!
!atmosphere for star 1
         if(rad1.lt.rmin1) then
            velr = zero
            temp = zero
            vth = zero
            opalbar = zero
            sline = zero
            scont = zero
            opac = zero         
         elseif(rad1.lt.rlim1) then
            velr = bvel(rad1, vinf, bconst, beta)
            temp = trad1
            vth = vthermal(vmicro1, temp, na)
            rho = mdot_cgs/(4.d0*pi*rad1**2 * sr1**2 * velr)
            opalbar = get_opalbar(iline, kline, sr1, yhe1, hei1, temp, vth_fiducial, xnue0, one, one, rho)/sr1  !in cgs
            gradv = velr * bconst*beta/rad1**2/(1.d0-bconst/rad1)
            sline = sobo1d(rad1, velr/vth, gradv/vth, opalbar*sr1, temp, xic1, xnue0, eps_line)
            scont = zero
            opac = zero
         endif         
!        
!----------------------line source function-----------------------------
!
         cs2_sline3d(i,j,k)=sline
!
!----------------------continuum source function------------------------
!
         cs2_scont3d(i,j,k)=scont
!
!------------------------velocity components----------------------------
!
!velocity components in system of star 1
         velx1 = velr*sin(theta1)*cos(phi1)
         vely1 = velr*sin(theta1)*sin(phi1)
         velz1 = velr*cos(theta1)
!velocity components in global system
         velx0 = velx1 + vx01
         vely0 = vely1 + vy01
         velz0 = velz1 + vz01
!velocity components in system of star 2         
         cs2_velx3d(i,j,k)=velx0 - vx02
         cs2_vely3d(i,j,k)=vely0 - vy02
         cs2_velz3d(i,j,k)=velz0 - vz02
!
!---------------------------temperature---------------------------------
!
         cs2_t3d(i,j,k) = temp
!
!-----------------------------opacity-----------------------------------
!
         cs2_rho3d(i,j,k) = zero !rho
         cs2_opalbar3d(i,j,k) = zero !opalbar*sr2  !in units of sr2
         cs2_opac3d(i,j,k) = zero!opac*sr2

         if(rad1.lt.rlim1) then
            cs2_opalbar3d(i,j,k)=opalbar*sr2
            cs2_sline3d(i,j,k)=sline
         endif
!
      enddo
   enddo
enddo
!
!
!-----------------------------set velocities------------------------------
!
!maximum velocity in global system
!
vmax=zero
do i=1, cs1_nr
   do j=1, cs1_ntheta
      do k=1, cs1_nphi
         vel=sqrt((vx01+cs1_velx3d(i,j,k))**2+(vy01+cs1_vely3d(i,j,k))**2+(vz01+cs1_velz3d(i,j,k))**2)
         if(vel.gt.vmax) vmax=vel
      enddo
   enddo
enddo

do i=1, cs2_nr
   do j=1, cs2_ntheta
      do k=1, cs2_nphi
         vel=sqrt((vx02+cs2_velx3d(i,j,k))**2+(vy02+cs2_vely3d(i,j,k))**2+(vz02+cs2_velz3d(i,j,k))**2)
         if(vel.gt.vmax) vmax=vel
      enddo
   enddo
enddo
!
!
!
!
end subroutine calc_model3db
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine calc_model3dc
!
!
use prog_type
use fund_const, only: zero, two, pi, xmsu, one
use options_modspec, only: input_file, input_file2, input_mod
use dime_modspec, only: cs1_nr, cs1_ntheta, cs1_nphi, cs1_r, cs1_theta, cs1_phi, &
                        cs1_velx3d, cs1_vely3d, cs1_velz3d, cs1_opac3d, &
                        cs1_opalbar3d, cs1_sline3d, cs1_scont3d, cs1_imask3d, cs1_t3d, cs1_rho3d, &
                        cs2_nr, cs2_ntheta, cs2_nphi, cs2_r, cs2_theta, cs2_phi, &
                        cs2_velx3d, cs2_vely3d, cs2_velz3d, cs2_opac3d, &
                        cs2_opalbar3d, cs2_sline3d, cs2_scont3d, cs2_imask3d, cs2_t3d, cs2_rho3d
use params_modspec, only: vx01, vy01, vz01, vx02, vy02, vz02, x01, y01, z01, x02, y02, z02, &
                          vmin1, vmax1, beta1, vrot1, rstar1, rmin1, rmax1, sr1, teff1, yhe1, hei1, logg1, lstar1, mdot1, vmicro1, &
                          vmin2, vmax2, beta2, vrot2, rstar2, rmin2, rmax2, sr2, teff2, yhe2, hei2, logg2, lstar2, mdot2, vmicro2, &
                          vmax, vth_fiducial, trad2, unit_length
use hdf5
use mod_grid, only: grid_log, grid_equi
use mod_opacities, only: get_opalbar
use mod_iline, only: iline, na, xnue0, kline
!
implicit none
!
! ... local scalars
integer(i4b) :: i, j, k, err
real(dp) :: vel, opac, opalbar, rho, sline, scont
real(dp) :: beta, vmin, vinf, mdot, mdot_cgs, bconst, xic2, eps_line, gradv, temp, velr, vth
real(dp) :: sint, cost, sinp, cosp, theta2, phi2, xcoord0, ycoord0, zcoord0, xcoord2, ycoord2, zcoord2, &
            velx0, vely0, velz0, velx2, vely2, velz2, rad, rad2
!
! ... local arrays
!
! ... local characters
!
! ... local functions
real(dp) :: bnue, bvel, vthermal, sobo1d
!
! ... for hdf5 file
!
write(*,*) '----------------calc_model3dc 3d Sobolev beta-vel-wind for star 2-------------'
write(*,*)
!
!define beta velocity wind
beta=1.d0
vmin=10.d5
vinf=2000.d5
mdot = 1.d-6
mdot_cgs=mdot*xmsu/(365.25d0*24.d0*3600.d0)
xic2=bnue(xnue0, trad2)
bconst=1.d0-(vmin/vinf)**(1.d0/beta)
!
vmicro2=1.d7
vth_fiducial=1.d7
logg2=3.6d0
lstar2=1.d6
vrot2=0.d8
yhe2=0.1d0
hei2=2.d0
eps_line=1.d-6
kline=1.d3
!
!
!
cs2_nr=101
cs2_ntheta=51
cs2_nphi=2*cs2_ntheta-1
!
!allocate arrays
allocate(cs2_r(cs2_nr), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3dc'
allocate(cs2_theta(cs2_ntheta), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3dc'
allocate(cs2_phi(cs2_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3dc'
allocate(cs2_opac3d(cs2_nr,cs2_ntheta,cs2_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3dc'
allocate(cs2_opalbar3d(cs2_nr,cs2_ntheta,cs2_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3dc'
allocate(cs2_scont3d(cs2_nr,cs2_ntheta,cs2_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3dc'
allocate(cs2_sline3d(cs2_nr,cs2_ntheta,cs2_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3dc'
allocate(cs2_velx3d(cs2_nr,cs2_ntheta,cs2_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3dc'
allocate(cs2_vely3d(cs2_nr,cs2_ntheta,cs2_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3dc'
allocate(cs2_velz3d(cs2_nr,cs2_ntheta,cs2_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3dc'
allocate(cs2_t3d(cs2_nr,cs2_ntheta,cs2_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3dc'
allocate(cs2_rho3d(cs2_nr,cs2_ntheta,cs2_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3dc'
!
!radius in units of sr2   
call grid_log(rmin2, rmax2, cs2_nr,  cs2_r)
!
!theta-grid equidistant
call grid_equi(zero, pi, cs2_ntheta, cs2_theta)
!
!phi-grid equidistant
call grid_equi(zero, two*pi, cs2_nphi, cs2_phi)
!
!
!
do i=1, cs2_nr
   velr = bvel(cs2_r(i), vinf, bconst, beta)
   temp = trad2
   vth = vthermal(vmicro2, temp, na)
   rho = mdot_cgs/(4.d0*pi*cs2_r(i)**2 * sr2**2 * velr)
   opalbar = get_opalbar(iline, kline, sr2, yhe2, hei2, temp, vth_fiducial, xnue0, one, one, rho)/sr2   !in cgs
   gradv = velr * bconst*beta/cs2_r(i)**2/(1.d0-bconst/cs2_r(i))
   sline = sobo1d(cs2_r(i), velr/vth, gradv/vth, opalbar*sr2, temp, xic2, xnue0, eps_line)
   scont = zero
   opac = zero
   
   do j=1, cs2_ntheta
      do k=1, cs2_nphi
!        
!----------------------line source function-----------------------------
!
         cs2_sline3d(i,j,k)=sline
!
!----------------------continuum source function------------------------
!
         cs2_scont3d(i,j,k)=scont
!
!------------------------velocity components----------------------------
!
         sint=sin(cs2_theta(j))
         cost=cos(cs2_theta(j))
         sinp=sin(cs2_phi(k))
         cosp=cos(cs2_phi(k))

         cs2_velx3d(i,j,k)=velr*sint*cosp
         cs2_vely3d(i,j,k)=velr*sint*sinp
         cs2_velz3d(i,j,k)=velr*cost
!
!---------------------------temperature---------------------------------
!
         cs2_t3d(i,j,k) = temp
!
!-----------------------------opacity-----------------------------------
!
         cs2_rho3d(i,j,k) = rho
         cs2_opalbar3d(i,j,k) = opalbar*sr2   !in units of sr2
         cs2_opac3d(i,j,k) = opac*sr2
!
      enddo
   enddo
enddo
!
!
!
write(*,*) '----------------calc_model3dc: calculating 3d model for star 1-----------------'
write(*,*)
!
!need to be defined (if gravity darkening and/or photospheric profile shall be accounted for in spec_vbin.eo)
vmicro1=1.d7
logg1=3.6d0
lstar1=1.d6
vrot1=0.d7
yhe1=0.1d0
!
cs1_nr=31
cs1_ntheta=31
cs1_nphi=2*cs1_ntheta-1
!
!allocate arrays
allocate(cs1_r(cs1_nr), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3dc'
allocate(cs1_theta(cs1_ntheta), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3dc'
allocate(cs1_phi(cs1_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3dc'
allocate(cs1_opac3d(cs1_nr,cs1_ntheta,cs1_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3dc'
allocate(cs1_opalbar3d(cs1_nr,cs1_ntheta,cs1_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3dc'
allocate(cs1_scont3d(cs1_nr,cs1_ntheta,cs1_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3dc'
allocate(cs1_sline3d(cs1_nr,cs1_ntheta,cs1_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3dc'
allocate(cs1_velx3d(cs1_nr,cs1_ntheta,cs1_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3dc'
allocate(cs1_vely3d(cs1_nr,cs1_ntheta,cs1_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3dc'
allocate(cs1_velz3d(cs1_nr,cs1_ntheta,cs1_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3dc'
allocate(cs1_t3d(cs1_nr,cs1_ntheta,cs1_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3dc'
allocate(cs1_rho3d(cs1_nr,cs1_ntheta,cs1_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3dc'
!
!radius in units of sr1   
call grid_log(rmin1, rmax1, cs1_nr,  cs1_r)
!
!theta-grid equidistant
call grid_equi(zero, pi, cs1_ntheta, cs1_theta)
!
!phi-grid equidistant
call grid_equi(zero, two*pi, cs1_nphi, cs1_phi)
!
!
!
do i=1, cs1_nr
   do j=1, cs1_ntheta
      do k=1, cs1_nphi
!
         rad=cs1_r(i)         
         sint=sin(cs1_theta(j))
         cost=cos(cs1_theta(j))
         sinp=sin(cs1_phi(k))
         cosp=cos(cs1_phi(k))

!position in global coordinate system
         xcoord0 = x01 + rad*sint*cosp*rstar1/unit_length
         ycoord0 = y01 + rad*sint*sinp*rstar1/unit_length
         zcoord0 = z01 + rad*cost*rstar1/unit_length
!postion in system of star 2
         xcoord2 = (xcoord0-x02)*unit_length/rstar2
         ycoord2 = (ycoord0-y02)*unit_length/rstar2
         zcoord2 = (zcoord0-z02)*unit_length/rstar2
         rad2=sqrt(xcoord2**2+ycoord2**2+zcoord2**2)
         call get_angles_spc(xcoord2, ycoord2, zcoord2, theta2, phi2)         
!
!atmosphere for star 2
         if(rad2.lt.rmin2) then
            velr = zero
            temp = zero
            vth = zero
            opalbar = zero
            sline = zero
            scont = zero
            opac = zero         
         else
            velr = bvel(rad2, vinf, bconst, beta)
            temp = trad2
            vth = vthermal(vmicro2, temp, na)
            rho = mdot_cgs/(4.d0*pi*rad2**2 * sr2**2 * velr)
            opalbar = get_opalbar(iline, kline, sr2, yhe2, hei2, temp, vth_fiducial, xnue0, one, one, rho)/sr2   !in cgs            
            gradv = velr * bconst*beta/rad2**2/(1.d0-bconst/rad2)
            sline = sobo1d(rad2, velr/vth, gradv/vth, opalbar*sr2, temp, xic2, xnue0, eps_line)
            scont = zero
            opac = zero
         endif         
!        
!----------------------line source function-----------------------------
!
         cs1_sline3d(i,j,k)=sline
!
!----------------------continuum source function------------------------
!
         cs1_scont3d(i,j,k)=scont
!
!------------------------velocity components----------------------------
!
!velocity components in system of star 2
         velx2 = velr*sin(theta2)*cos(phi2)
         vely2 = velr*sin(theta2)*sin(phi2)
         velz2 = velr*cos(theta2)
!velocity components in global system
         velx0 = velx2 + vx02
         vely0 = vely2 + vy02
         velz0 = velz2 + vz02
!velocity components in system of star 1         
         cs1_velx3d(i,j,k)=velx0 - vx01
         cs1_vely3d(i,j,k)=vely0 - vy01
         cs1_velz3d(i,j,k)=velz0 - vz01
!
!---------------------------temperature---------------------------------
!
         cs1_t3d(i,j,k) = temp
!
!-----------------------------opacity-----------------------------------
!
         cs1_rho3d(i,j,k) = rho
         cs1_opalbar3d(i,j,k) = opalbar*sr1  !in units of sr1
         cs1_opac3d(i,j,k) = opac*sr1
!
      enddo
   enddo
enddo
!
!
!-----------------------------set velocities------------------------------
!
!maximum velocity in global system
!
vmax=zero
do i=1, cs1_nr
   do j=1, cs1_ntheta
      do k=1, cs1_nphi
         vel=sqrt((vx01+cs1_velx3d(i,j,k))**2+(vy01+cs1_vely3d(i,j,k))**2+(vz01+cs1_velz3d(i,j,k))**2)
         if(vel.gt.vmax) vmax=vel
      enddo
   enddo
enddo

do i=1, cs2_nr
   do j=1, cs2_ntheta
      do k=1, cs2_nphi
         vel=sqrt((vx02+cs2_velx3d(i,j,k))**2+(vy02+cs2_vely3d(i,j,k))**2+(vz02+cs2_velz3d(i,j,k))**2)
         if(vel.gt.vmax) vmax=vel
      enddo
   enddo
enddo
!
!
!
!
end subroutine calc_model3dc
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine calc_model3dd
!
!star 1 has a Be-type disc from Dylans initial conditions  
!
use prog_type
use fund_const
use options_modspec, only: input_file, input_file2, input_mod
use dime_modspec, only: cs1_nr, cs1_ntheta, cs1_nphi, cs1_r, cs1_theta, cs1_phi, &
                        cs1_velx3d, cs1_vely3d, cs1_velz3d, cs1_opac3d, &
                        cs1_opalbar3d, cs1_sline3d, cs1_scont3d, cs1_imask3d, cs1_t3d, cs1_rho3d, &
                        cs2_nr, cs2_ntheta, cs2_nphi, cs2_r, cs2_theta, cs2_phi, &
                        cs2_velx3d, cs2_vely3d, cs2_velz3d, cs2_opac3d, &
                        cs2_opalbar3d, cs2_sline3d, cs2_scont3d, cs2_imask3d, cs2_t3d, cs2_rho3d
use params_modspec, only: vx01, vy01, vz01, vx02, vy02, vz02, x01, y01, z01, x02, y02, z02, &
                          vmin1, vmax1, beta1, vrot1, rstar1, rmin1, rmax1, sr1, teff1, yhe1, hei1, logg1, lstar1, mdot1, vmicro1, &
                          vmin2, vmax2, beta2, vrot2, rstar2, rmin2, rmax2, sr2, teff2, yhe2, hei2, logg2, lstar2, mdot2, vmicro2, &
                          vmax, vth_fiducial, trad1, unit_length
use hdf5
use mod_grid, only: grid_log, grid_equi
use mod_opacities, only: opac_thomson, get_opalbar
use mod_iline, only: iline, na, xnue0, kline
!
implicit none
!
! ... local scalars
integer(i4b) :: i, j, k, err, idum1, idum2, ii
real(dp) :: vel, opac, opalbar, rho, sline, scont
real(dp) :: beta, vmin, vinf, mdot, mdot_cgs, bconst, xic1, eps_line, gradv, temp, velr, vth
real(dp) :: sint1, cost1, sinp1, cosp1, theta1, phi1, rad1, rad2, sint2, cost2, sinp2, cosp2, &
            xcoord0, ycoord0, zcoord0, xcoord1, ycoord1, zcoord1, &
            velx0, vely0, velz0, velx1, vely1, velz1, xcoord2, ycoord2, zcoord2, kcont

real(dp) :: hi1, mmw, mstar1_cgs, tdisc, mdisc, dtheta_disc, theta_disc0, theta_disc1, csound, velphi, rho_disc0, b2, b3, rdisc_max, slope
!
! ... local arrays
real(dp), dimension(:), allocatable :: fdum1_arr, fdum2_arr
!
! ... local characters
!
! ... local functions
real(dp) :: bnue, bvel, vthermal, vsound, mean_molecular_weight, sline_depcoeff
!
! ... for hdf5 file
!
kcont=1.d0
kline=1.d0
!
write(*,*) '----------------calc_model3dd: Be-type disc for star 1 (Dylan)----------------'
write(*,*)
!
hei1=2.d0
!
!define local parameters (and maybe overwrite global ones)
hi1 = 1.d0   !number free electrons for each hydrogen atom
mmw = mean_molecular_weight(hi1,hei1,yhe1)  !mean molecular weight
mstar1_cgs = sr1**2 * ten**logg1/cgs_grav
vrot1 = sqrt(cgs_grav*mstar1_cgs/sr1)    !breakup velocity (consistent with v_phi of disc model)
tdisc = 10.d3            !isothermal temperature of the disc
mdisc = 5.d-10 * xmsu     !*1.d-10 set to almost zero
dtheta_disc = 45.d0*pi/180.d0    !opening angle of the disc (plus/minus 5 degree here)
theta_disc0 = pi/two - dtheta_disc
theta_disc1 = pi/two + dtheta_disc
!maximum radius of the disc (in rstar1)
rdisc_max=13.5d0
!rdisc_max=26.d0
!rdisc_max=40.d0
!slope of the disc
slope=3.5d0
!
!
cs1_nr=101
idum1=12
idum2=20
cs1_ntheta=2*idum1+2*idum2+1
cs1_nphi=2*cs1_ntheta-1
!
!allocate arrays
allocate(cs1_r(cs1_nr), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3dd'
allocate(cs1_theta(cs1_ntheta), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3dd'
allocate(cs1_phi(cs1_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3dd'
allocate(cs1_opac3d(cs1_nr,cs1_ntheta,cs1_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3dd'
allocate(cs1_opalbar3d(cs1_nr,cs1_ntheta,cs1_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3dd'
allocate(cs1_scont3d(cs1_nr,cs1_ntheta,cs1_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3dd'
allocate(cs1_sline3d(cs1_nr,cs1_ntheta,cs1_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3dd'
allocate(cs1_velx3d(cs1_nr,cs1_ntheta,cs1_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3dd'
allocate(cs1_vely3d(cs1_nr,cs1_ntheta,cs1_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3dd'
allocate(cs1_velz3d(cs1_nr,cs1_ntheta,cs1_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3dd'
allocate(cs1_t3d(cs1_nr,cs1_ntheta,cs1_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3dd'
allocate(cs1_rho3d(cs1_nr,cs1_ntheta,cs1_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3dd'
!
!radius in units of sr1   
call grid_log(rmin1, rmax1, cs1_nr,  cs1_r)
!
if(cs1_r(cs1_nr-1).le.rdisc_max) then
   write(*,*) 'error: we require always two zero-opacity ghost points at outer boundary'
   stop
endif
!
!theta-grid with small steps at equator (always ntheta/2 for the disc)
allocate(fdum1_arr(idum1))
call grid_log(pi/2.d0+5.d-3, pi/2.d0+dtheta_disc-1.d-3, idum1, fdum1_arr)
allocate(fdum2_arr(idum2))
call grid_log(2.d0*fdum1_arr(idum1)-fdum1_arr(idum1-1), pi, idum2, fdum2_arr)
ii=1
do i=1, idum2
   cs1_theta(ii)=pi-fdum2_arr(idum2+1-i)
   ii=ii+1
enddo
do i=1, idum1
   cs1_theta(ii)=pi-fdum1_arr(idum1+1-i)
   ii=ii+1
enddo
cs1_theta(ii)=pi/two
ii=ii+1
do i=1, idum1
   cs1_theta(ii)=fdum1_arr(i)
   ii=ii+1
enddo
do i=1, idum2
   cs1_theta(ii)=fdum2_arr(i)
   ii=ii+1
enddo
cs1_theta(1)=zero
cs1_theta(cs1_ntheta)=pi
!
!
!call grid_equi(zero, pi, cs1_ntheta, cs1_theta)
!
!phi-grid equidistant
call grid_equi(zero, two*pi, cs1_nphi, cs1_phi)
!
!
!
do i=1, cs1_nr
   do j=1, cs1_ntheta
      do k=1, cs1_nphi
!   do j=33, 33!cs1_ntheta
!      do k=1, 1!cs1_nphi
!
!zero values if inside secondary         
         rad1=cs1_r(i)
         sint1=sin(cs1_theta(j))
         cost1=cos(cs1_theta(j))
         sinp1=sin(cs1_phi(k))
         cosp1=cos(cs1_phi(k))
!position in global coordinate system
         xcoord0 = x01 + rad1*sint1*cosp1*rstar1/unit_length
         ycoord0 = y01 + rad1*sint1*sinp1*rstar1/unit_length
         zcoord0 = z01 + rad1*cost1*rstar1/unit_length
!postion in system of star 2
         xcoord2 = (xcoord0-x02)*unit_length/rstar2
         ycoord2 = (ycoord0-y02)*unit_length/rstar2
         zcoord2 = (zcoord0-z02)*unit_length/rstar2
         rad2=sqrt(xcoord2**2+ycoord2**2+zcoord2**2)      
!
         if(rad2.lt.rmin2) then
            velphi = zero
            temp = zero
            vth = zero
            opalbar = zero
            sline = zero
            scont = zero
            opac = zero
         elseif(rad1.gt.rdisc_max) then
            velphi = zero
            temp = zero
            vth = zero
            opalbar = zero
            sline = zero
            scont = zero
            opac = zero
         elseif(cs1_theta(j).gt.theta_disc0.and.&
                cs1_theta(j).lt.theta_disc1) then
!
            sint1=sin(cs1_theta(j))
            velphi = sqrt(cgs_grav*mstar1_cgs/cs1_r(i)/sr1/sint1)
            temp = tdisc
            vth = vthermal(vmicro1, temp, na)
            csound = vsound(temp,mmw)
            rho_disc0 = mdisc*sqrt(cgs_grav*mstar1_cgs)*log(ten)/two/sqrt(pi**3)/csound/sr1**3.5d0
            rho = rho_disc0*(cs1_r(i)*sint1)**(-slope) * exp(cgs_grav*mstar1_cgs/csound**2 * (one/sr1/cs1_r(i)-one/sr1/cs1_r(i)/sint1))
!
!departure coeffcients in LTE
            b2=one
            b3=one
            opalbar = get_opalbar(iline, kline, sr1, yhe1, hei1, temp, vth_fiducial, xnue0, b2, b3, rho)/sr1   !in cgs
            sline = sline_depcoeff(xnue0, temp, b2, b3)
            scont = bnue(xnue0,tdisc)
            opac = opac_thomson(yhe1, hei1, rho, kcont)
!            write(*,*) cs1_r(i), temp, b2, b3, rho, opalbar*sr1
         else
            velphi = zero
            temp = trad1
            vth = vthermal(vmicro1, temp, na)
            rho = zero
            opalbar = zero
            sline = zero
            scont = zero
            opac = zero         
         endif   
!        
!----------------------line source function-----------------------------
!
         cs1_sline3d(i,j,k)=sline
!
!----------------------continuum source function------------------------
!
         cs1_scont3d(i,j,k)=scont
!
!------------------------velocity components----------------------------
!
         sint1=sin(cs1_theta(j))
         cost1=cos(cs1_theta(j))
         sinp1=sin(cs1_phi(k))
         cosp1=cos(cs1_phi(k))

         cs1_velx3d(i,j,k) = -velphi*sinp1
         cs1_vely3d(i,j,k) = velphi*cosp1
         cs1_velz3d(i,j,k) = zero
!
!---------------------------temperature---------------------------------
!
         cs1_t3d(i,j,k) = temp
!
!-----------------------------opacity-----------------------------------
!
         cs1_rho3d(i,j,k) = rho
         cs1_opalbar3d(i,j,k) = opalbar*sr1   !in units of sr1
         cs1_opac3d(i,j,k) = opac*sr1
!
      enddo
   enddo
enddo
!
!stop 'go on in modelspec'

!j=30
!write(*,*) cs1_theta(j)*180./pi
!k=1
!do i=1, cs1_nr
!   write(*,*) i, cs1_r(i), cs1_opalbar3d(i,j,k)
!enddo
!stop

!
!
write(*,*) '----------------calc_model3dd: calculating 3d model for star 2-----------------'
write(*,*)
!
!need to be defined (if gravity darkening and/or photospheric profile shall be accounted for in spec_vbin.eo)
!
cs2_nr=31
cs2_ntheta=31
cs2_nphi=2*cs2_ntheta-1
!
!allocate arrays
allocate(cs2_r(cs2_nr), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3dd'
allocate(cs2_theta(cs2_ntheta), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3dd'
allocate(cs2_phi(cs2_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3dd'
allocate(cs2_opac3d(cs2_nr,cs2_ntheta,cs2_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3dd'
allocate(cs2_opalbar3d(cs2_nr,cs2_ntheta,cs2_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3dd'
allocate(cs2_scont3d(cs2_nr,cs2_ntheta,cs2_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3dd'
allocate(cs2_sline3d(cs2_nr,cs2_ntheta,cs2_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3dd'
allocate(cs2_velx3d(cs2_nr,cs2_ntheta,cs2_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3dd'
allocate(cs2_vely3d(cs2_nr,cs2_ntheta,cs2_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3dd'
allocate(cs2_velz3d(cs2_nr,cs2_ntheta,cs2_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3dd'
allocate(cs2_t3d(cs2_nr,cs2_ntheta,cs2_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3dd'
allocate(cs2_rho3d(cs2_nr,cs2_ntheta,cs2_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3dd'
!
!radius in units of sr1   
call grid_log(rmin2, rmax2, cs2_nr,  cs2_r)
!
!theta-grid equidistant
call grid_equi(zero, pi, cs2_ntheta, cs2_theta)
!
!phi-grid equidistant
call grid_equi(zero, two*pi, cs2_nphi, cs2_phi)
!
!
!
do i=1, cs2_nr
   do j=1, cs2_ntheta
      do k=1, cs2_nphi
!
         rad2=cs2_r(i)         
         sint2=sin(cs2_theta(j))
         cost2=cos(cs2_theta(j))
         sinp2=sin(cs2_phi(k))
         cosp2=cos(cs2_phi(k))

!position in global coordinate system
         xcoord0 = x02 + rad2*sint2*cosp2*rstar2/unit_length
         ycoord0 = y02 + rad2*sint2*sinp2*rstar2/unit_length
         zcoord0 = z02 + rad2*cost2*rstar2/unit_length
!postion in system of star 1
         xcoord1 = (xcoord0-x01)*unit_length/rstar1
         ycoord1 = (ycoord0-y01)*unit_length/rstar1
         zcoord1 = (zcoord0-z01)*unit_length/rstar1
         rad1=sqrt(xcoord1**2+ycoord1**2+zcoord1**2)
         call get_angles_spc(xcoord1, ycoord1, zcoord1, theta1, phi1)         
!
!atmosphere for star 1
         if(rad1.lt.rmin1.or.rad1.gt.rdisc_max) then
            velphi = zero
            temp = zero
            vth = zero
            opalbar = zero
            sline = zero
            scont = zero
            opac = zero         
         else

            if(theta1.gt.theta_disc0.and.&
               theta1.lt.theta_disc1) then
!              
               sint1=sin(theta1)
               velphi = sqrt(cgs_grav*mstar1_cgs/rad1/sr1/sint1)
               temp = tdisc
               vth = vthermal(vmicro1, temp, na)
               csound = vsound(temp,mmw)
               rho_disc0 = mdisc*sqrt(cgs_grav*mstar1_cgs)*log(ten)/two/sqrt(pi**3)/csound/sr1**3.5d0
               rho = rho_disc0*(rad1*sint1)**(-slope) * exp(cgs_grav*mstar1_cgs/csound**2 * (one/sr1/rad1-one/sr1/rad1/sint1))
!
               b2=one
               b3=one
         
               opalbar = get_opalbar(iline, kline, sr1, yhe1, hei1, temp, vth_fiducial, xnue0, b2, b3, rho)/sr1   !in cgs                        
               sline = sline_depcoeff(xnue0, temp, b2, b3)
               scont = bnue(xnue0,tdisc)
               opac = opac_thomson(yhe1, hei1, rho, kcont)


!               write(*,*) rad2, rad1, velphi, theta1*180.d0/pi, rho_disc0, rho, opalbar
            else
               velphi = zero
               temp = trad1
               vth = vthermal(vmicro1, temp, na)
               rho = zero
               opalbar = zero
               sline = zero
               scont = zero
               opac = zero         
            endif
         endif

!         write(*,*) rad1, rad2, opalbar*sr2, theta1, phi1, cs2_theta(j), cs2_phi(k)
!        
!----------------------line source function-----------------------------
!
         cs2_sline3d(i,j,k)=sline
!
!----------------------continuum source function------------------------
!
         cs2_scont3d(i,j,k)=scont
!
!------------------------velocity components----------------------------
!
!velocity components in system of star 1
         velx1 = -velphi*sin(phi1)
         vely1 = velphi*cos(phi1)
         velz1 = zero
!velocity components in global system
         velx0 = velx1 + vx01
         vely0 = vely1 + vy01
         velz0 = velz1 + vz01
!velocity components in system of star 2         
         cs2_velx3d(i,j,k)=velx0 - vx02
         cs2_vely3d(i,j,k)=vely0 - vy02
         cs2_velz3d(i,j,k)=velz0 - vz02
!
!---------------------------temperature---------------------------------
!
         cs2_t3d(i,j,k) = temp
!
!-----------------------------opacity-----------------------------------
!
         cs2_rho3d(i,j,k) = rho
         cs2_opalbar3d(i,j,k) = opalbar*sr2  !in units of sr2
         cs2_opac3d(i,j,k) = opac*sr2
!
      enddo
   enddo
enddo
!
!stop
!
!-----------------------------set velocities------------------------------
!
!maximum velocity in global system
!
vmax=zero
do i=1, cs1_nr
   do j=1, cs1_ntheta
      do k=1, cs1_nphi
         vel=sqrt((vx01+cs1_velx3d(i,j,k))**2+(vy01+cs1_vely3d(i,j,k))**2+(vz01+cs1_velz3d(i,j,k))**2)
         if(vel.gt.vmax) vmax=vel
      enddo
   enddo
enddo

do i=1, cs2_nr
   do j=1, cs2_ntheta
      do k=1, cs2_nphi
         vel=sqrt((vx02+cs2_velx3d(i,j,k))**2+(vy02+cs2_vely3d(i,j,k))**2+(vz02+cs2_velz3d(i,j,k))**2)
         if(vel.gt.vmax) vmax=vel
      enddo
   enddo
enddo
!
!
!test a rotational velocity
!vrot1=300.d5
!vrot2=300.d5
!
!for now, set zero continuum
!cs1_opac3d=0.d0
!cs1_scont3d=0d0
!cs2_opac3d=0.d0
!cs2_scont3d=0.d0
!
!
end subroutine calc_model3dd
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine calc_model3df
!
!star 1 has a Be-type disc from Dylans initial conditions  
!star 2 has a similar disc
!
use prog_type
use fund_const
use options_modspec, only: input_file, input_file2, input_mod, indat_file
use dime_modspec, only: cs1_nr, cs1_ntheta, cs1_nphi, cs1_r, cs1_theta, cs1_phi, &
                        cs1_velx3d, cs1_vely3d, cs1_velz3d, cs1_opac3d, &
                        cs1_opalbar3d, cs1_sline3d, cs1_scont3d, cs1_imask3d, cs1_t3d, cs1_rho3d, &
                        cs2_nr, cs2_ntheta, cs2_nphi, cs2_r, cs2_theta, cs2_phi, &
                        cs2_velx3d, cs2_vely3d, cs2_velz3d, cs2_opac3d, &
                        cs2_opalbar3d, cs2_sline3d, cs2_scont3d, cs2_imask3d, cs2_t3d, cs2_rho3d
use params_modspec, only: vx01, vy01, vz01, vx02, vy02, vz02, x01, y01, z01, x02, y02, z02, &
                          vmin1, vmax1, beta1, vrot1, rstar1, rmin1, rmax1, sr1, teff1, yhe1, yhe1_mass, hei1, logg1, lstar1, mdot1, vmicro1, &
                          vmin2, vmax2, beta2, vrot2, rstar2, rmin2, rmax2, sr2, teff2, yhe2, yhe2_mass, hei2, logg2, lstar2, mdot2, vmicro2, &
                          vmax, vth_fiducial, trad1, unit_length
use hdf5
use mod_grid, only: grid_log, grid_equi
use mod_opacities, only: opac_thomson, opac_thomson3, get_opalbar, get_opalbar2, get_opalbar3
use mod_iline, only: iline, na, xnue0, kline
use mod_lte, only: get_lte_table
!
implicit none
!
! ... local scalars
integer(i4b) :: opt_temperature1, opt_temperature2, opt_vrot1, opt_vrot2
integer(i4b) :: i, j, k, err, idum1, idum2, ii
real(dp) :: vel, opac, opalbar, rho, sline, scont
real(dp) :: beta, vmin, vinf, mdot, mdot_cgs, bconst, xic1, eps_line, gradv, temp, velr, vth
real(dp) :: sint1, cost1, sinp1, cosp1, theta1, phi1, rad1, theta2, phi2, rad2, sint2, cost2, sinp2, cosp2, &
            xcoord0, ycoord0, zcoord0, xcoord1, ycoord1, zcoord1, &
            velx0, vely0, velz0, velx1, vely1, velz1, velx2, vely2, velz2, xcoord2, ycoord2, zcoord2, kcont, fdum1, fdum2, fdum3, fdum4

real(dp) :: hi1, mmw1, mstar1_cgs, tdisc1, mdisc1, dtheta_disc1, theta_disc1l, theta_disc1u, rho0_disc1, rdisc1_max, slope1
real(dp) :: hi2, mmw2, mstar2_cgs, tdisc2, mdisc2, dtheta_disc2, theta_disc2l, theta_disc2u, rho0_disc2, rdisc2_max, slope2
real(dp) :: csound, velphi, b2, b3, vrot2fac, tdisc1a, tdisc1b, tacoeff1, tbcoeff1, tdisc2a, tdisc2b, tacoeff2, tbcoeff2
real(dp) :: vrot1_local, vrot2_local
!
! ... local arrays
real(dp), dimension(:), allocatable :: fdum1_arr, fdum2_arr
!
! ... local characters
!
! ... local functions
real(dp) :: bnue, bvel, vthermal, vsound, mean_molecular_weight, sline_depcoeff, dilfac
!
! ... for temperature structure
integer(i4b) :: cs1_np, cs2_np
integer(i4b), parameter :: nz=301
real(dp) :: zmin, zmax, tauz, dtau, rad
real(dp), dimension(:), allocatable :: cs1_p1d, cs1_tauz1d, cs2_p1d, cs2_tauz1d
real(dp), dimension(nz) :: opac1d, z1d
!

namelist / input_usr/ opt_temperature1, opt_temperature2, opt_vrot1, opt_vrot2, kcont, kline, vrot2fac, vth_fiducial, &
     vmicro1, vrot1, logg1, yhe1, hi1, hei1, lstar1, tdisc1, tdisc1a, tdisc1b, mdisc1, rdisc1_max, dtheta_disc1, slope1, &
     vmicro2, vrot2, logg2, yhe2, hi2, hei2, lstar2, tdisc2, tdisc2a, tdisc2b, mdisc2, rdisc2_max, dtheta_disc2, slope2
!
! ... for hdf5 file
!
write(*,*) '-----------------------------reading indat file--------------------------------'
write(*,*) 'reading usr input from file: ', trim(indat_file)
write(*,*)
!
!opt_temperature
!=0 for constant temperature at tdisc
!=1 for linear temperature law
!=2 for 1/r temperature law
!=3 for lucy optically thin temperature law
!=4 for carciofi temperature law
!=5 for temperature from dilution of b-star
!
open(1, file=trim(indat_file), status='old', form='formatted')
   rewind 1
   read(1, nml=input_usr)
close(1)
!
!
vmicro1=vmicro1*1.d5
vmicro2=vmicro2*1.d5
vth_fiducial=vth_fiducial*1.d5
!
mdisc1 = mdisc1 * xmsu
mdisc2 = mdisc2 * xmsu
!
!maximum radius of the disc (in rstar1)
rdisc1_max = rdisc1_max*unit_length/rstar1  !10.d0*3.7/rstar1
rdisc2_max = rdisc2_max*unit_length/rstar2
!
dtheta_disc1=dtheta_disc1*pi/180.d0
dtheta_disc2=dtheta_disc2*pi/180.d0
!
yhe1_mass = one / (one + one/four/yhe1)
yhe2_mass = one / (one + one/four/yhe2)
!
!get the lte table
call get_lte_table(yhe1_mass)
!
!----------------------define disc 1--------------------------
!
mmw1 = mean_molecular_weight(hi1,hei1,yhe1)  !mean molecular weight
mstar1_cgs = sr1**2 * ten**logg1/cgs_grav
vrot1_local = sqrt(cgs_grav*mstar1_cgs/sr1)    !breakup velocity (consistent with v_phi of disc model)
!
!opening angle of the disc (plus/minus 45 degree here)
theta_disc1l = pi/two - dtheta_disc1
theta_disc1u = pi/two + dtheta_disc1
!
!temperature coefficients
tacoeff1 = tdisc1a-(tdisc1b-tdisc1a)*one/(rdisc1_max-one)
tbcoeff1 = (tdisc1b-tdisc1a)/(rdisc1_max-one)
!
!-------------------------define disc 2-----------------------
!
!define disc 2
mmw2 = mean_molecular_weight(hi2,hei2,yhe2)  !mean molecular weight
mstar2_cgs = sr2**2 * ten**logg2/cgs_grav
vrot2_local = sqrt(cgs_grav*mstar2_cgs/sr2)    !breakup velocity (consistent with v_phi of disc model)
!
theta_disc2l = pi/two - dtheta_disc2
theta_disc2u = pi/two + dtheta_disc2
!
tacoeff2 = tdisc2a-(tdisc2b-tdisc2a)*one/(rdisc2_max-one)
tbcoeff2 = (tdisc2b-tdisc2a)/(rdisc2_max-one)
!
!
!
write(*,*) '----------------calc_model3df: Be-type disc for star 1 (Dylan)----------------'
write(*,*)
!
!
cs1_nr=101
idum1=12
idum2=20
cs1_ntheta=2*idum1+2*idum2+1
cs1_nphi=2*cs1_ntheta-1
!
!allocate arrays
allocate(cs1_r(cs1_nr), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3df'
allocate(cs1_theta(cs1_ntheta), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3df'
allocate(cs1_phi(cs1_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3df'
allocate(cs1_opac3d(cs1_nr,cs1_ntheta,cs1_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3df'
allocate(cs1_opalbar3d(cs1_nr,cs1_ntheta,cs1_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3df'
allocate(cs1_scont3d(cs1_nr,cs1_ntheta,cs1_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3df'
allocate(cs1_sline3d(cs1_nr,cs1_ntheta,cs1_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3df'
allocate(cs1_velx3d(cs1_nr,cs1_ntheta,cs1_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3df'
allocate(cs1_vely3d(cs1_nr,cs1_ntheta,cs1_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3df'
allocate(cs1_velz3d(cs1_nr,cs1_ntheta,cs1_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3df'
allocate(cs1_t3d(cs1_nr,cs1_ntheta,cs1_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3df'
allocate(cs1_rho3d(cs1_nr,cs1_ntheta,cs1_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3df'
!
!radius in units of sr1   
call grid_log(rmin1, rmax1, cs1_nr,  cs1_r)
!
if(cs1_r(cs1_nr-1).le.rdisc1_max) then
   write(*,*) 'error: we require always two zero-opacity ghost points at outer boundary'
   stop
endif
!
!theta-grid with small steps at equator (always ntheta/2 for the disc)
allocate(fdum1_arr(idum1))
call grid_log(pi/2.d0+5.d-3, pi/2.d0+dtheta_disc1-1.d-3, idum1, fdum1_arr)
allocate(fdum2_arr(idum2))
call grid_log(2.d0*fdum1_arr(idum1)-fdum1_arr(idum1-1), pi, idum2, fdum2_arr)
ii=1
do i=1, idum2
   cs1_theta(ii)=pi-fdum2_arr(idum2+1-i)
   ii=ii+1
enddo
do i=1, idum1
   cs1_theta(ii)=pi-fdum1_arr(idum1+1-i)
   ii=ii+1
enddo
cs1_theta(ii)=pi/two
ii=ii+1
do i=1, idum1
   cs1_theta(ii)=fdum1_arr(i)
   ii=ii+1
enddo
do i=1, idum2
   cs1_theta(ii)=fdum2_arr(i)
   ii=ii+1
enddo
cs1_theta(1)=zero
cs1_theta(cs1_ntheta)=pi
!
!
!call grid_equi(zero, pi, cs1_ntheta, cs1_theta)
!
!phi-grid equidistant
call grid_equi(zero, two*pi, cs1_nphi, cs1_phi)
!

cs2_nr=101
idum1=12
idum2=20
cs2_ntheta=2*idum1+2*idum2+1
cs2_nphi=2*cs1_ntheta-1
!
!allocate arrays
allocate(cs2_r(cs2_nr), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3df'
allocate(cs2_theta(cs2_ntheta), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3df'
allocate(cs2_phi(cs2_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3df'
allocate(cs2_opac3d(cs2_nr,cs2_ntheta,cs2_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3df'
allocate(cs2_opalbar3d(cs2_nr,cs2_ntheta,cs2_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3df'
allocate(cs2_scont3d(cs2_nr,cs2_ntheta,cs2_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3df'
allocate(cs2_sline3d(cs2_nr,cs2_ntheta,cs2_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3df'
allocate(cs2_velx3d(cs2_nr,cs2_ntheta,cs2_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3df'
allocate(cs2_vely3d(cs2_nr,cs2_ntheta,cs2_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3df'
allocate(cs2_velz3d(cs2_nr,cs2_ntheta,cs2_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3df'
allocate(cs2_t3d(cs2_nr,cs2_ntheta,cs2_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3df'
allocate(cs2_rho3d(cs2_nr,cs2_ntheta,cs2_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3df'
!
!radius in units of sr2
call grid_log(rmin2, rmax2, cs2_nr,  cs2_r)
!
if(cs2_r(cs2_nr-1).le.rdisc2_max) then
   write(*,*) 'error: we require always two zero-opacity ghost points at outer boundary'
   stop
endif
!
!theta-grid with small steps at equator (always ntheta/2 for the disc)
deallocate(fdum1_arr)
deallocate(fdum2_arr)
allocate(fdum1_arr(idum1))
call grid_log(pi/2.d0+5.d-3, pi/2.d0+dtheta_disc2-1.d-3, idum1, fdum1_arr)
allocate(fdum2_arr(idum2))
call grid_log(2.d0*fdum1_arr(idum1)-fdum1_arr(idum1-1), pi, idum2, fdum2_arr)
ii=1
do i=1, idum2
   cs2_theta(ii)=pi-fdum2_arr(idum2+1-i)
   ii=ii+1
enddo
do i=1, idum1
   cs2_theta(ii)=pi-fdum1_arr(idum1+1-i)
   ii=ii+1
enddo
cs2_theta(ii)=pi/two
ii=ii+1
do i=1, idum1
   cs2_theta(ii)=fdum1_arr(i)
   ii=ii+1
enddo
do i=1, idum2
   cs2_theta(ii)=fdum2_arr(i)
   ii=ii+1
enddo
cs2_theta(1)=zero
cs2_theta(cs2_ntheta)=pi
!
!
!call grid_equi(zero, pi, cs1_ntheta, cs1_theta)
!
!phi-grid equidistant
call grid_equi(zero, two*pi, cs2_nphi, cs2_phi)
!
!-------calculate temperature structure as function of r*sin(theta)-----
!
cs1_np=cs1_nr
allocate(cs1_p1d(cs1_np))
allocate(cs1_tauz1d(cs1_np))
cs1_p1d=cs1_r
!
csound = vsound(teff1,mmw1)
rho0_disc1 = mdisc1*sqrt(cgs_grav*mstar1_cgs)*log(ten)/two/sqrt(pi**3)/csound/sr1**3.5d0
!
do i=1, cs1_np
   zmin = cs1_p1d(i)*cos(theta_disc1u)   
   zmax = cs1_p1d(i)*cos(theta_disc1l)
   call grid_equi(zmin, zmax, nz, z1d)
   do j=1, nz
      rad = sqrt(z1d(j)**2 + cs1_p1d(i)**2)
      rho = rho0_disc1*(cs1_p1d(i))**(-slope1) * exp(cgs_grav*mstar1_cgs/csound**2 * (one/sr1/rad-one/sr1/cs1_p1d(i)))
      opac1d(j) = opac_thomson(yhe1, hei1, rho, kcont)*sr1
   enddo
   tauz = 0.d0
   do j=2, nz
      dtau = 0.5d0*(opac1d(j-1)+opac1d(j))*(z1d(j)-z1d(j-1))
      tauz = tauz + dtau
   enddo
   cs1_tauz1d(i) = tauz
!   write(*,*) cs1_p1d(i), zmin, zmax, tauz
enddo
!
!
!
cs2_np=cs2_nr
allocate(cs2_p1d(cs2_np))
allocate(cs2_tauz1d(cs2_np))
cs2_p1d=cs2_r
!
csound = vsound(teff2,mmw2)
rho0_disc2 = mdisc2*sqrt(cgs_grav*mstar2_cgs)*log(ten)/two/sqrt(pi**3)/csound/sr2**3.5d0
!write(*,*) rho0_disc1, mstar1_cgs/xmsu, vrot1_local/1.d8
!write(*,*) rho0_disc2, mstar2_cgs/xmsu, vrot2_local/1.d8
!stop
do i=1, cs2_np
   zmin = cs2_p1d(i)*cos(theta_disc2u)   
   zmax = cs2_p1d(i)*cos(theta_disc2l)
   call grid_equi(zmin, zmax, nz, z1d)
   do j=1, nz
      rad = sqrt(z1d(j)**2 + cs2_p1d(i)**2)
      rho = rho0_disc2*(cs2_p1d(i))**(-slope2) * exp(cgs_grav*mstar2_cgs/csound**2 * (one/sr2/rad-one/sr2/cs2_p1d(i)))
      opac1d(j) = opac_thomson(yhe2, hei2, rho, kcont)*sr2
   enddo
   tauz = 0.d0
   do j=2, nz
      dtau = 0.5d0*(opac1d(j-1)+opac1d(j))*(z1d(j)-z1d(j-1))
      tauz = tauz + dtau
   enddo
   cs2_tauz1d(i) = tauz
!  write(*,*) cs2_p1d(i), zmin, zmax, tauz
enddo
!stop
!
!-----------------------------------------------------------------------
!
!calculate corresponding mass-accretion-rate if alpha-disc model
!
csound = vsound(tdisc1,mmw1)
rho0_disc1 = mdisc1*sqrt(cgs_grav*mstar1_cgs)*log(ten)/two/sqrt(pi**3)/csound/sr1**3.5d0
!
!stellar radius in 10^10 cm
fdum1 = (sr1/1.d10)**(15.d0/8.d0)
!
!mass on solar masses
fdum2 = (mstar1_cgs/xmsu)**(-5.d0/8.d0)
!
!alpha parameter
fdum3 = (0.01d0)**(7.d0/10.d0)
!
!mass accretion rate in 10^16 g/s
fdum4 = rho0_disc1*fdum1*fdum2*fdum3/3.1d-8
fdum4 = fdum4**(20.d0/11.d0)
!
!mass accretion rate in msun/yr
fdum4 = fdum4*1.d16 * yr/xmsu
!
write(*,*) rho0_disc1, sr1/rsu, mstar1_cgs/xmsu, mstar2_cgs/xmsu, fdum4
!stop 'go on here'
!-----------------------------------------------------------------------
!
do i=1, cs1_nr
   do j=1, cs1_ntheta
      do k=1, cs1_nphi
!   do j=33, 33!cs1_ntheta
!      do k=1, 1!cs1_nphi
!
!zero values if inside secondary         
         rad1=cs1_r(i)
         sint1=sin(cs1_theta(j))
         cost1=cos(cs1_theta(j))
         sinp1=sin(cs1_phi(k))
         cosp1=cos(cs1_phi(k))
!position in global coordinate system
         xcoord0 = x01 + rad1*sint1*cosp1*rstar1/unit_length
         ycoord0 = y01 + rad1*sint1*sinp1*rstar1/unit_length
         zcoord0 = z01 + rad1*cost1*rstar1/unit_length
!postion in system of star 2
         xcoord2 = (xcoord0-x02)*unit_length/rstar2
         ycoord2 = (ycoord0-y02)*unit_length/rstar2
         zcoord2 = (zcoord0-z02)*unit_length/rstar2
         rad2=sqrt(xcoord2**2+ycoord2**2+zcoord2**2)
         call get_angles_spc(xcoord2, ycoord2, zcoord2, theta2, phi2)
!
         if(rad2.lt.rmin2) then
            velphi = zero
            temp = zero
            vth = zero
            opalbar = zero
            sline = zero
            scont = zero
            opac = zero
         elseif(rad1.gt.rdisc1_max) then
            velphi = zero
            temp = zero
            vth = zero
            opalbar = zero
            sline = zero
            scont = zero
            opac = zero
            if(rad2.lt.rdisc2_max.and. &
               theta2.gt.theta_disc2l.and. &
               theta2.lt.theta_disc2u) then   !set disc 2

               sint2=sin(theta2)
               velphi = sqrt(cgs_grav*mstar2_cgs/rad2/sr2/sint2) * vrot2fac
!velocity components in system of star 2
               velx2 = -velphi*sin(phi2)
               vely2 = velphi*cos(phi2)
               velz2 = zero
!velocity components in global system
               velx0 = velx2 + vx02
               vely0 = vely2 + vy02
               velz0 = velz2 + vz02
!velocity components in system of star 1
               velx1 = velx0 - vx01
               velx1 = vely0 - vy01
               velz1 = velz0 - vz01
               velphi = -sinp1*velx1 + cosp1*vely1

               if(opt_temperature2.eq.0) then
                  temp=tdisc2
               elseif(opt_temperature2.eq.1) then
                  temp = tacoeff2 + tbcoeff2*rad2
               elseif(opt_temperature2.eq.2) then
                  temp = max(teff2/rad2,5d3)
               elseif(opt_temperature2.eq.3) then
                  temp = teff2 * dilfac(rad2)**0.25d0
               elseif(opt_temperature2.eq.4) then
!neglect  sint term
!                  temp = teff2/pi**0.25d0 * (asin(one/rad2/sint2) - sqrt(one-one/(rad2*sint2)**2)/rad2/sint2)**0.25d0
                  temp = teff2/pi**0.25d0 * (asin(one/rad2) - sqrt(one-one/(rad2)**2)/rad2)**0.25d0
                  temp = max(temp,0.6d0*teff2)  !poor mans temperature algorithm (instead of finding correct tauz)
!                  if(cs2_tauz1d(**).lt.0.1d0) temp = 0.6d0*teff
               else
                  stop 'opt_temperature not properly specified'
               endif
               
               vth = vthermal(vmicro2, temp, na)
               csound = vsound(temp,mmw1)
               rho0_disc2 = mdisc2*sqrt(cgs_grav*mstar2_cgs)*log(ten)/two/sqrt(pi**3)/csound/sr2**3.5d0
               rho = rho0_disc2*(rad2*sint2)**(-slope2) * exp(cgs_grav*mstar2_cgs/csound**2 * (one/sr2/rad2-one/sr2/rad2/sint2))
!departure coeffcients in LTE (everything calculated for mmw and yhe of star 1)
               b2=one
               b3=one
               opalbar = get_opalbar(iline, kline, sr2, yhe1, hei1, temp, vth_fiducial, xnue0, b2, b3, rho)/sr2   !in cgs
               sline = sline_depcoeff(xnue0, temp, b2, b3)
               scont = bnue(xnue0,tdisc2)
               opac = opac_thomson(yhe1, hei1, rho, kcont)
            endif
               
         elseif(cs1_theta(j).gt.theta_disc1l.and.&
                cs1_theta(j).lt.theta_disc1u) then
!
            sint1=sin(cs1_theta(j))
            velphi = sqrt(cgs_grav*mstar1_cgs/cs1_r(i)/sr1/sint1)

            if(opt_temperature1.eq.0) then
               temp=tdisc1
            elseif(opt_temperature1.eq.1) then
               temp = tacoeff1 + tbcoeff1*rad1
            elseif(opt_temperature1.eq.2) then
               temp = max(teff1/rad1,5d3)
            elseif(opt_temperature1.eq.3) then
               temp = teff1 * dilfac(rad1)**0.25d0
            elseif(opt_temperature1.eq.4) then
               temp = teff1/pi**0.25d0 * (asin(one/rad1) - sqrt(one-one/(rad1)**2)/rad1)**0.25d0
               if(cs1_tauz1d(i).lt.0.1d0) temp = max(temp, 0.6d0*teff1)
            elseif(opt_temperature1.eq.5) then
               temp = dilfac(rad2)**0.25d0 * teff2
               temp = max(temp,5.d3)
!               write(*,*) rad1, rad2, dilfac(rad2), temp               
            else
               stop 'opt_temperature not properly specified'
            endif
!            
            vth = vthermal(vmicro1, temp, na)
            csound = vsound(temp,mmw1)
            rho0_disc1 = mdisc1*sqrt(cgs_grav*mstar1_cgs)*log(ten)/two/sqrt(pi**3)/csound/sr1**3.5d0
            rho = rho0_disc1*(cs1_r(i)*sint1)**(-slope1) * exp(cgs_grav*mstar1_cgs/csound**2 * (one/sr1/cs1_r(i)-one/sr1/cs1_r(i)/sint1))
            if(rho.ne.rho) then
               write(*,*) rho, rho0_disc1, cs1_r(i), sint1, slope1, mstar1_cgs, csound, sr1, rho0_disc1*(cs1_r(i)*sint1)**(-slope1)
               stop
            endif
!
!departure coeffcients in LTE
            b2=one
            b3=one
            if(temp.lt.6.1d3) then
               opalbar = get_opalbar2(iline, kline, sr1, yhe1, 1.d0, temp, vth_fiducial, xnue0, b2, b3, rho)/sr1  !in cgs               
            else
               opalbar = get_opalbar(iline, kline, sr1, yhe1, hei1, temp, vth_fiducial, xnue0, b2, b3, rho)/sr1   !in cgs
            endif
            sline = sline_depcoeff(xnue0, temp, b2, b3)
            scont = bnue(xnue0,temp)
            if(temp.lt.6.1d3) then            
               opac = opac_thomson3(yhe1, 1.d0, 0.d0, rho, kcont) !helium only partially ionized, no electrons from hydrogen
            else            
               opac = opac_thomson(yhe1, hei1, rho, kcont)
            endif
               
!            write(*,*) rho, opalbar
!            write(*,*) cs1_r(i), temp, b2, b3, rho, opalbar*sr1
         else
            velphi = zero
            temp = 0.8d0*(teff1+teff2)/2.d0
            vth = vthermal(vmicro1, temp, na)
            rho = zero
            opalbar = zero
            sline = zero
            scont = zero
            opac = zero         
         endif   
!        
!----------------------line source function-----------------------------
!
         cs1_sline3d(i,j,k)=sline
!
!----------------------continuum source function------------------------
!
         cs1_scont3d(i,j,k)=scont
!
!------------------------velocity components----------------------------
!
         sint1=sin(cs1_theta(j))
         cost1=cos(cs1_theta(j))
         sinp1=sin(cs1_phi(k))
         cosp1=cos(cs1_phi(k))

         cs1_velx3d(i,j,k) = -velphi*sinp1
         cs1_vely3d(i,j,k) = velphi*cosp1
         cs1_velz3d(i,j,k) = zero
!
!---------------------------temperature---------------------------------
!
         cs1_t3d(i,j,k) = temp
!
!-----------------------------opacity-----------------------------------
!
         cs1_rho3d(i,j,k) = rho
         cs1_opalbar3d(i,j,k) = opalbar*sr1   !in units of sr1
         cs1_opac3d(i,j,k) = opac*sr1
!
      enddo
   enddo
enddo
!
!stop 'go on in modelspec'

!j=30
!write(*,*) cs1_theta(j)*180./pi
!k=1
!do i=1, cs1_nr
!   write(*,*) i, cs1_r(i), cs1_opalbar3d(i,j,k)
!enddo
!stop

!
!
write(*,*) '----------------calc_model3df: Be-type disc for star 2 (Dylan)-----------------'
write(*,*)
!
!
!rho0_disc2 = mdisc2*sqrt(cgs_grav*mstar2_cgs)*log(ten)/two/sqrt(pi**3)/csound/sr2**3.5d0
!rho0_disc1 = mdisc1*sqrt(cgs_grav*mstar1_cgs)*log(ten)/two/sqrt(pi**3)/csound/sr1**3.5d0
!write(*,*) rho0_disc1, rho0_disc2
!write(*,*) sr1, sr2
!stop
!
do i=1, cs2_nr
   do j=1, cs2_ntheta
      do k=1, cs2_nphi
!   do j=16, 16!cs2_ntheta
!      do k=1, 1!cs2_nphi
!
         rad2=cs2_r(i)
         theta2=cs2_theta(j)
         sint2=sin(cs2_theta(j))
         cost2=cos(cs2_theta(j))
         sinp2=sin(cs2_phi(k))
         cosp2=cos(cs2_phi(k))

!position in global coordinate system
         xcoord0 = x02 + rad2*sint2*cosp2*rstar2/unit_length
         ycoord0 = y02 + rad2*sint2*sinp2*rstar2/unit_length
         zcoord0 = z02 + rad2*cost2*rstar2/unit_length
!postion in system of star 1
         xcoord1 = (xcoord0-x01)*unit_length/rstar1
         ycoord1 = (ycoord0-y01)*unit_length/rstar1
         zcoord1 = (zcoord0-z01)*unit_length/rstar1
         rad1=sqrt(xcoord1**2+ycoord1**2+zcoord1**2)
         call get_angles_spc(xcoord1, ycoord1, zcoord1, theta1, phi1)         
!
!atmosphere for star 1
         if(rad1.lt.rmin1) then
            velphi = zero
            temp = zero
            vth = zero
            opalbar = zero
            sline = zero
            scont = zero
            opac = zero
         elseif(rad1.gt.rdisc1_max) then
            velphi = zero
            temp = zero
            vth = zero
            opalbar = zero
            sline = zero
            scont = zero
            opac = zero
            if(rad2.lt.rdisc2_max.and. &
               theta2.gt.theta_disc2l.and. &
               theta2.lt.theta_disc2u) then   !set disc 2

               sint2=sin(theta2)
               velphi = sqrt(cgs_grav*mstar2_cgs/rad2/sr2/sint2) * vrot2fac

               if(opt_temperature2.eq.0) then
                  temp=tdisc2
               elseif(opt_temperature2.eq.1) then
                  temp = tacoeff2 + tbcoeff2*rad2
               elseif(opt_temperature2.eq.2) then
                  temp = max(teff2/rad2,5d3)
               elseif(opt_temperature2.eq.3) then
                  temp = teff2 * dilfac(rad2)**0.25d0
               elseif(opt_temperature2.eq.4) then
                  temp = teff2/pi**0.25d0 * (asin(one/rad2) - sqrt(one-one/(rad2)**2)/rad2)**0.25d0
                  if(cs2_tauz1d(i).lt.0.15d0) temp = max(temp, 0.6d0*teff2)
               else
                  stop 'opt_temperature not properly specified'
               endif
               
               vth = vthermal(vmicro2, temp, na)
               csound = vsound(temp,mmw1)
               rho0_disc2 = mdisc2*sqrt(cgs_grav*mstar2_cgs)*log(ten)/two/sqrt(pi**3)/csound/sr2**3.5d0
               rho = rho0_disc2*(rad2*sint2)**(-slope2) * exp(cgs_grav*mstar2_cgs/csound**2 * (one/sr2/rad2-one/sr2/rad2/sint2))
!departure coeffcients in LTE (everything calculated for mmw and yhe of star 1)
               b2=one
               b3=one
!               write(*,*) theta2*180./pi, theta_disc2l*180./pi,  theta_disc2u*180./pi, temp, rad2, rho0_disc2
               opalbar = get_opalbar(iline, kline, sr2, yhe2, hei2, temp, vth_fiducial, xnue0, b2, b3, rho)/sr2   !in cgs
!               write(*,*) 't2'
               sline = sline_depcoeff(xnue0, temp, b2, b3)
               scont = bnue(xnue0,temp)
               opac = opac_thomson(yhe2, hei2, rho, kcont)
!               write(*,*) theta2*180.d0/pi, theta_disc2l*180.d0/pi, theta_disc2u*180.d0/pi
!               write(*,*) rho, opalbar
            endif            
         elseif(theta1.gt.theta_disc1l.and.&
                theta1.lt.theta_disc1u) then
!              
               sint1=sin(theta1)
               velphi = sqrt(cgs_grav*mstar1_cgs/rad1/sr1/sint1)
!velocity components in system of star 1
               velx1 = -velphi*sin(phi1)
               vely1 = velphi*cos(phi1)
               velz1 = zero
!velocity components in global system
               velx0 = velx1 + vx01
               vely0 = vely1 + vy01
               velz0 = velz1 + vz01
!velocity components in system of star 2         
               velx2 = velx0 - vx02
               vely2 = vely0 - vy02
               velz2 = velz0 - vz02
               velphi = -sinp2*velx2 + cosp2*vely2
               
               if(opt_temperature1.eq.0) then
                  temp=tdisc1
               elseif(opt_temperature1.eq.1) then
                  temp = tacoeff1 + tbcoeff1*rad1
               elseif(opt_temperature1.eq.2) then
                  temp = max(teff1/rad1,5d3)
               elseif(opt_temperature1.eq.3) then
                  temp = teff1 * dilfac(rad1)**0.25d0
               elseif(opt_temperature1.eq.4) then
                  temp = teff1/pi**0.25d0 * (asin(one/rad1) - sqrt(one-one/(rad1)**2)/rad1)**0.25d0
                  temp = max(temp,0.6d0*teff1)  !poor mans temperature algorithm (instead of finding correct tauz)
               else
                  stop 'opt_temperature not properly specified'
               endif

               
               vth = vthermal(vmicro1, temp, na)
               csound = vsound(tdisc1,mmw1)
               rho0_disc1 = mdisc1*sqrt(cgs_grav*mstar1_cgs)*log(ten)/two/sqrt(pi**3)/csound/sr1**3.5d0
               rho = rho0_disc1*(rad1*sint1)**(-slope1) * exp(cgs_grav*mstar1_cgs/csound**2 * (one/sr1/rad1-one/sr1/rad1/sint1))
!
               b2=one
               b3=one
         
               opalbar = get_opalbar(iline, kline, sr1, yhe1, hei1, temp, vth_fiducial, xnue0, b2, b3, rho)/sr1   !in cgs                        
               sline = sline_depcoeff(xnue0, temp, b2, b3)
               scont = bnue(xnue0,temp)
               opac = opac_thomson(yhe1, hei1, rho, kcont)


!               write(*,*) rad2, rad1, velphi, theta1*180.d0/pi, rho_disc0, rho, opalbar
         else
            velphi = zero
            temp = 0.8d0*(teff1+teff2)/2.d0
            vth = vthermal(vmicro1, temp, na)
            rho = zero
            opalbar = zero
            sline = zero
            scont = zero
            opac = zero         
         endif

!         write(*,*) rad1, rad2, opalbar*sr2, theta1, phi1, cs2_theta(j), cs2_phi(k)
!        
!----------------------line source function-----------------------------
!
         cs2_sline3d(i,j,k)=sline
!
!----------------------continuum source function------------------------
!
         cs2_scont3d(i,j,k)=scont
!
!------------------------velocity components----------------------------
!
!velocity components in system of star 2
         sint2=sin(cs2_theta(j))
         cost2=cos(cs2_theta(j))
         sinp2=sin(cs2_phi(k))
         cosp2=cos(cs2_phi(k))

         cs2_velx3d(i,j,k) = -velphi*sinp2
         cs2_vely3d(i,j,k) = velphi*cosp2
         cs2_velz3d(i,j,k) = zero         
!
!---------------------------temperature---------------------------------
!
         cs2_t3d(i,j,k) = temp
!
!-----------------------------opacity-----------------------------------
!
         cs2_rho3d(i,j,k) = rho
         cs2_opalbar3d(i,j,k) = opalbar*sr2  !in units of sr2
         cs2_opac3d(i,j,k) = opac*sr2
!
      enddo
   enddo
enddo
!
!stop 'go on in modelspec'
!
!-----------------------------set velocities------------------------------
!
!maximum velocity in global system
!
vmax=zero
do i=1, cs1_nr
   do j=1, cs1_ntheta
      do k=1, cs1_nphi
         vel=sqrt((vx01+cs1_velx3d(i,j,k))**2+(vy01+cs1_vely3d(i,j,k))**2+(vz01+cs1_velz3d(i,j,k))**2)
         if(vel.gt.vmax) vmax=vel
      enddo
   enddo
enddo

do i=1, cs2_nr
   do j=1, cs2_ntheta
      do k=1, cs2_nphi
         vel=sqrt((vx02+cs2_velx3d(i,j,k))**2+(vy02+cs2_vely3d(i,j,k))**2+(vz02+cs2_velz3d(i,j,k))**2)
         if(vel.gt.vmax) vmax=vel
      enddo
   enddo
enddo
!
!
!
!zero rotation of primary
select case(opt_vrot1)
   case(1)
      vrot1=vrot1_local
   case default
      vrot1=vrot1*1.d5
end select
!
select case(opt_vrot2)
   case(1)
      vrot2=vrot2_local
   case default
      vrot2=vrot2*1.d5
end select
!
!for now, set zero continuum
!cs1_opac3d=0.d0
!cs1_scont3d=0d0
!cs2_opac3d=0.d0
!cs2_scont3d=0.d0
!
end subroutine calc_model3df
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine calc_model3d_adisc
!
!star 1 has a alpha-disc model
!star 2 has no disc
!
use prog_type
use fund_const
use options_modspec, only: input_file, input_file2, input_mod, indat_file
use dime_modspec, only: cs1_nr, cs1_ntheta, cs1_nphi, cs1_r, cs1_theta, cs1_phi, &
                        cs1_velx3d, cs1_vely3d, cs1_velz3d, cs1_opac3d, &
                        cs1_opalbar3d, cs1_sline3d, cs1_scont3d, cs1_imask3d, cs1_t3d, cs1_rho3d, &
                        cs2_nr, cs2_ntheta, cs2_nphi, cs2_r, cs2_theta, cs2_phi, &
                        cs2_velx3d, cs2_vely3d, cs2_velz3d, cs2_opac3d, &
                        cs2_opalbar3d, cs2_sline3d, cs2_scont3d, cs2_imask3d, cs2_t3d, cs2_rho3d
use params_modspec, only: vx01, vy01, vz01, vx02, vy02, vz02, x01, y01, z01, x02, y02, z02, &
                          vmin1, vmax1, beta1, vrot1, rstar1, rmin1, rmax1, sr1, teff1, yhe1, hei1, logg1, lstar1, mdot1, vmicro1, &
                          vmin2, vmax2, beta2, vrot2, rstar2, rmin2, rmax2, sr2, teff2, yhe2, hei2, logg2, lstar2, mdot2, vmicro2, &
                          vmax, vth_fiducial, trad1, unit_length
use hdf5
use mod_grid, only: grid_log, grid_equi
use mod_opacities, only: opac_thomson, get_opalbar, get_opalbar2, get_opalbar3
use mod_iline, only: iline, na, xnue0, kline
!
implicit none
!
! ... local scalars
integer(i4b) :: opt_temperature1, opt_temperature2
integer(i4b) :: i, j, k, err, idum1, idum2, rdum1, rdum2, ii
real(dp) :: vel, opac, opalbar, rho, sline, scont
real(dp) :: beta, vmin, vinf, mdot, mdot_cgs, bconst, xic1, eps_line, gradv, temp, velr, vth
real(dp) :: sint1, cost1, sinp1, cosp1, theta1, phi1, rad1, theta2, phi2, rad2, sint2, cost2, sinp2, cosp2, &
            xcoord0, ycoord0, zcoord0, xcoord1, ycoord1, zcoord1, &
            velx0, vely0, velz0, velx1, vely1, velz1, velx2, vely2, velz2, xcoord2, ycoord2, zcoord2, kcont

real(dp) :: hi1, mmw1, mstar1_cgs, mdot1_cgs, mdot1_16, rschwarz1, rschwarz1_cgs, tdisc1, mdisc1, dtheta_disc1, theta_disc1l, theta_disc1u, rho0_disc1, rdisc1_max, rdisc1_min, slope1
real(dp) :: hi2, mmw2, mstar2_cgs, tdisc2, mdisc2, dtheta_disc2, theta_disc2l, theta_disc2u, rho0_disc2, rdisc2_max, slope2
real(dp) :: rcyc, rcyc_10, zcyc, csound, velphi, b2, b3, vrot1fac, vrot2fac, tdisc1a, tdisc1b, tacoeff1, tbcoeff1, tdisc2a, tdisc2b, tacoeff2, tbcoeff2

real(dp) :: fdum
!
! ... local arrays
real(dp), dimension(:), allocatable :: fdum1_arr, fdum2_arr
!
! ... local characters
!
! ... local functions
real(dp) :: bnue, bvel, vthermal, vsound, mean_molecular_weight, sline_depcoeff, dilfac
!
! ... for temperature structure
real(dp), dimension(:), allocatable :: cs1_rdum
!
! ... namelist
real(dp) :: alpha1, mstar1, mstar2, tmin1
namelist / input_usr/ mstar1, mstar2, tmin1, vmicro1, vmicro2, vrot1, vrot2, yhe1, yhe2, hei1, hei2, hi1, hi2, mdot1, alpha1, kcont, kline, vth_fiducial, rdisc1_min, rdisc1_max
!
! ... for hdf5 file
!
write(*,*) '-----------------------------reading indat file--------------------------------'
write(*,*) 'reading usr input from file: ', trim(indat_file)
write(*,*)
!
open(1, file=trim(indat_file), status='old', form='formatted')
   rewind 1
   read(1, nml=input_usr)
close(1)
!
!in cgs
mdot1_cgs = mdot1*xmsu/yr
mstar1_cgs = mstar1*xmsu
mstar2_cgs = mstar2*xmsu
logg1 = log10(cgs_grav*mstar1_cgs/sr1**2)
logg2 = log10(cgs_grav*mstar2_cgs/sr2**2)
!
!calculate schwarzschild-radius
rschwarz1_cgs = two*cgs_grav*mstar1_cgs/cgs_clight**2
rschwarz1 = rschwarz1_cgs/rsu
!
!minimum radius of disc at 6*schwarzschild-radius
rmin1 = 6.d0*rschwarz1
!
!velocities in cgs
vmicro1=vmicro1*1.d5
vmicro2=vmicro2*1.d5
vrot1 = vrot1*1.d5
vrot2 = vrot2*1.d5
vth_fiducial = vth_fiducial*1.d5
!
if(abs(rmin1-rstar1).gt.small_number) then
   write(*,*) 'error in calc_model3d_adisc: minimum radius of disc ne input rstar1'
   write(*,*) 'set rstar1 to', rmin1
   stop
endif
rmin1 = one
!
!temperature at the innermost point of the disc
teff1 = (three*cgs_grav*mstar1_cgs*mdot1_cgs/eight/pi/cgs_sb/sr1**3)**0.25d0
!
!
!define parameters that are always required, however not important here
lstar1=1.d6
lstar2=1.d6
!
!calculate mean molecular weight
mmw1 = mean_molecular_weight(hi1,hei1,yhe1)
mmw2 = mean_molecular_weight(hi2,hei2,yhe2)
!
!calculate new rotational velocities
vrot1 = sqrt(cgs_grav*mstar1_cgs/sr1)
vrot2 = sqrt(cgs_grav*mstar2_cgs/sr2)

!opening angle of disc (plus/minus 45 degree here)
dtheta_disc1 = 45.d0*pi/180.d0
theta_disc1l = pi/two - dtheta_disc1
theta_disc1u = pi/two + dtheta_disc1
!
dtheta_disc2 = 45.d0*pi/180.d0
theta_disc2l = pi/two - dtheta_disc2
theta_disc2u = pi/two + dtheta_disc2
!
!
!minimum and maximum disc radii
rdisc2_max=0.d0    !in rstar2
!
!scale-factors for rotational velocities
vrot1fac=1.d0
vrot2fac=1.d0
!
!define exponents for mid-plane structure of alpha-disc model
!
write(*,*) '-------------------calc_model3d_adisc: alpha-disc model-----------------------'
write(*,*)
!
!
rdum1 = 97
rdum2 = 4
cs1_nr = rdum1+rdum2
idum1=12
idum2=20
cs1_ntheta=2*idum1+2*idum2+1
cs1_nphi=2*cs1_ntheta-1
!
!allocate arrays
allocate(cs1_r(cs1_nr), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3d_adisc'
allocate(cs1_theta(cs1_ntheta), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3d_adisc'
allocate(cs1_phi(cs1_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3d_adisc'
allocate(cs1_opac3d(cs1_nr,cs1_ntheta,cs1_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3d_adisc'
allocate(cs1_opalbar3d(cs1_nr,cs1_ntheta,cs1_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3d_adisc'
allocate(cs1_scont3d(cs1_nr,cs1_ntheta,cs1_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3d_adisc'
allocate(cs1_sline3d(cs1_nr,cs1_ntheta,cs1_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3d_adisc'
allocate(cs1_velx3d(cs1_nr,cs1_ntheta,cs1_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3d_adisc'
allocate(cs1_vely3d(cs1_nr,cs1_ntheta,cs1_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3d_adisc'
allocate(cs1_velz3d(cs1_nr,cs1_ntheta,cs1_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3d_adisc'
allocate(cs1_t3d(cs1_nr,cs1_ntheta,cs1_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3d_adisc'
allocate(cs1_rho3d(cs1_nr,cs1_ntheta,cs1_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3d_adisc'
!
!radius in units of sr1
!call grid_log(rmin1, rmax1, cs1_nr,  cs1_r)
allocate(cs1_rdum(rdum1))
call grid_log(rmin1, 0.999d0*rdisc1_max, rdum1, cs1_rdum)
cs1_r(1:rdum1) = cs1_rdum
deallocate(cs1_rdum)
allocate(cs1_rdum(rdum2))
call grid_equi(1.0001*rdisc1_max, rmax1, rdum2, cs1_rdum)
cs1_r(rdum1+1:cs1_nr) = cs1_rdum
deallocate(cs1_rdum)
!
if(cs1_r(cs1_nr-1).le.rdisc1_max) then
   write(*,*) 'error: we require always two zero-opacity ghost points at outer boundary'
   stop
endif
!
!theta-grid with small steps at equator (always ntheta/2 for the disc)
allocate(fdum1_arr(idum1))
call grid_log(pi/2.d0+5.d-3, pi/2.d0+dtheta_disc1-1.d-3, idum1, fdum1_arr)
allocate(fdum2_arr(idum2))
call grid_log(2.d0*fdum1_arr(idum1)-fdum1_arr(idum1-1), pi, idum2, fdum2_arr)
ii=1
do i=1, idum2
   cs1_theta(ii)=pi-fdum2_arr(idum2+1-i)
   ii=ii+1
enddo
do i=1, idum1
   cs1_theta(ii)=pi-fdum1_arr(idum1+1-i)
   ii=ii+1
enddo
cs1_theta(ii)=pi/two
ii=ii+1
do i=1, idum1
   cs1_theta(ii)=fdum1_arr(i)
   ii=ii+1
enddo
do i=1, idum2
   cs1_theta(ii)=fdum2_arr(i)
   ii=ii+1
enddo
cs1_theta(1)=zero
cs1_theta(cs1_ntheta)=pi
!
!
!call grid_equi(zero, pi, cs1_ntheta, cs1_theta)
!
!phi-grid equidistant
call grid_equi(zero, two*pi, cs1_nphi, cs1_phi)
!
!
!
cs2_nr=101
idum1=12
idum2=20
cs2_ntheta=2*idum1+2*idum2+1
cs2_nphi=2*cs1_ntheta-1
!
!allocate arrays
allocate(cs2_r(cs2_nr), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3d_adisc'
allocate(cs2_theta(cs2_ntheta), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3d_adisc'
allocate(cs2_phi(cs2_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3d_adisc'
allocate(cs2_opac3d(cs2_nr,cs2_ntheta,cs2_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3d_adisc'
allocate(cs2_opalbar3d(cs2_nr,cs2_ntheta,cs2_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3d_adisc'
allocate(cs2_scont3d(cs2_nr,cs2_ntheta,cs2_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3d_adisc'
allocate(cs2_sline3d(cs2_nr,cs2_ntheta,cs2_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3d_adisc'
allocate(cs2_velx3d(cs2_nr,cs2_ntheta,cs2_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3d_adisc'
allocate(cs2_vely3d(cs2_nr,cs2_ntheta,cs2_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3d_adisc'
allocate(cs2_velz3d(cs2_nr,cs2_ntheta,cs2_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3d_adisc'
allocate(cs2_t3d(cs2_nr,cs2_ntheta,cs2_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3d_adisc'
allocate(cs2_rho3d(cs2_nr,cs2_ntheta,cs2_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3d_adisc'
!
!radius in units of sr2
call grid_log(rmin2, rmax2, cs2_nr,  cs2_r)
!
if(cs2_r(cs2_nr-1).le.rdisc2_max) then
   write(*,*) 'error: we require always two zero-opacity ghost points at outer boundary'
   stop
endif
!
!theta-grid with small steps at equator (always ntheta/2 for the disc)
deallocate(fdum1_arr)
deallocate(fdum2_arr)
allocate(fdum1_arr(idum1))
call grid_log(pi/2.d0+5.d-3, pi/2.d0+dtheta_disc2-1.d-3, idum1, fdum1_arr)
allocate(fdum2_arr(idum2))
call grid_log(2.d0*fdum1_arr(idum1)-fdum1_arr(idum1-1), pi, idum2, fdum2_arr)
ii=1
do i=1, idum2
   cs2_theta(ii)=pi-fdum2_arr(idum2+1-i)
   ii=ii+1
enddo
do i=1, idum1
   cs2_theta(ii)=pi-fdum1_arr(idum1+1-i)
   ii=ii+1
enddo
cs2_theta(ii)=pi/two
ii=ii+1
do i=1, idum1
   cs2_theta(ii)=fdum1_arr(i)
   ii=ii+1
enddo
do i=1, idum2
   cs2_theta(ii)=fdum2_arr(i)
   ii=ii+1
enddo
cs2_theta(1)=zero
cs2_theta(cs2_ntheta)=pi
!
!
!call grid_equi(zero, pi, cs1_ntheta, cs1_theta)
!
!phi-grid equidistant
call grid_equi(zero, two*pi, cs2_nphi, cs2_phi)
!
!-----------------------------------------------------------------------
!
do i=1, cs1_nr
   do j=1, cs1_ntheta
      do k=1, cs1_nphi
!   do j=33, 33!cs1_ntheta
!      do k=1, 1!cs1_nphi
!
!zero values if inside secondary
         rad1=cs1_r(i)
         sint1=sin(cs1_theta(j))
         cost1=cos(cs1_theta(j))
         sinp1=sin(cs1_phi(k))
         cosp1=cos(cs1_phi(k))
!position in global coordinate system
         xcoord0 = x01 + rad1*sint1*cosp1*rstar1/unit_length
         ycoord0 = y01 + rad1*sint1*sinp1*rstar1/unit_length
         zcoord0 = z01 + rad1*cost1*rstar1/unit_length
!postion in system of star 2
         xcoord2 = (xcoord0-x02)*unit_length/rstar2
         ycoord2 = (ycoord0-y02)*unit_length/rstar2
         zcoord2 = (zcoord0-z02)*unit_length/rstar2
         rad2=sqrt(xcoord2**2+ycoord2**2+zcoord2**2)
         call get_angles_spc(xcoord2, ycoord2, zcoord2, theta2, phi2)
!
         if(rad2.lt.rmin2) then
            velphi = zero
            temp = zero
            vth = zero
            opalbar = zero
            sline = zero
            scont = zero
            opac = zero
         elseif(rad1.gt.rdisc1_max) then
            velphi = zero
            temp = zero
            vth = zero
            opalbar = zero
            sline = zero
            scont = zero
            opac = zero
!            if(rad2.lt.rdisc2_max.and. &
!               theta2.gt.theta_disc2l.and. &
!               theta2.lt.theta_disc2u) then   !set disc 2
!
!               sint2=sin(theta2)
!               velphi = sqrt(cgs_grav*mstar2_cgs/rad2/sr2/sint2) * vrot2fac
!!velocity components in system of star 2
!               velx2 = -velphi*sin(phi2)
!               vely2 = velphi*cos(phi2)
!               velz2 = zero
!!velocity components in global system
!               velx0 = velx2 + vx02
!               vely0 = vely2 + vy02
!               velz0 = velz2 + vz02
!!velocity components in system of star 1
!               velx1 = velx0 - vx01
!               velx1 = vely0 - vy01
!               velz1 = velz0 - vz01
!               velphi = -sinp1*velx1 + cosp1*vely1
!
!               if(opt_temperature2.eq.0) then
!                  temp=tdisc2
!               elseif(opt_temperature2.eq.1) then
!                  temp = tacoeff2 + tbcoeff2*rad2
!               elseif(opt_temperature2.eq.2) then
!                  temp = max(teff2/rad2,5d3)
!               elseif(opt_temperature2.eq.3) then
!                  temp = teff2 * dilfac(rad2)**0.25d0
!               elseif(opt_temperature2.eq.4) then
!!neglect  sint term
!!                  temp = teff2/pi**0.25d0 * (asin(one/rad2/sint2) - sqrt(one-one/(rad2*sint2)**2)/rad2/sint2)**0.25d0
!                  temp = teff2/pi**0.25d0 * (asin(one/rad2) - sqrt(one-one/(rad2)**2)/rad2)**0.25d0
!                  temp = max(temp,0.6d0*teff2)  !poor mans temperature algorithm (instead of finding correct tauz)
!!                  if(cs2_tauz1d(**).lt.0.1d0) temp = 0.6d0*teff
!               else
!                  stop 'opt_temperature not properly specified'
!               endif
!               
!               vth = vthermal(vmicro2, temp, na)
!               csound = vsound(temp,mmw1)
!               rho0_disc2 = mdisc2*sqrt(cgs_grav*mstar2_cgs)*log(ten)/two/sqrt(pi**3)/csound/sr2**3.5d0
!               rho = rho0_disc2*(rad2*sint2)**(-slope2) * exp(cgs_grav*mstar2_cgs/csound**2 * (one/sr2/rad2-one/sr2/rad2/sint2))
!!departure coeffcients in LTE (everything calculated for mmw and yhe of star 1)
!               b2=one
!               b3=one
!               opalbar = get_opalbar(iline, kline, sr2, yhe1, hei1, temp, vth_fiducial, xnue0, b2, b3, rho)/sr2   !in cgs
!               sline = sline_depcoeff(xnue0, temp, b2, b3)
!               scont = bnue(xnue0,tdisc2)
!               opac = opac_thomson(yhe1, hei1, rho, kcont)
!            endif
               
         elseif(cs1_theta(j).gt.theta_disc1l.and.&
                cs1_theta(j).lt.theta_disc1u.and.rad1.ge.rdisc1_min) then
!
            sint1=sin(cs1_theta(j))
            cost1=cos(cs1_theta(j))
!calculate corresponding rcyc and zcyc in cylindrical coordinates (and use minimum of 1.d0)
            rcyc=max(rad1*sint1,one)
            zcyc=abs(rad1*cost1)

            call model_adisc(rcyc, zcyc, rstar1, mstar1, alpha1, mdot1, mmw1, temp, rho, velr, velphi)
!            temp = max(temp,tmin1)
!
!departure coeffcients in LTE
            b2=one
            b3=one
            if(temp.lt.6.1d3) then
               !assume helium not ionized, hydrogen from ionization balance
               opalbar = get_opalbar3(iline, kline, sr1, yhe1, 0.d0, temp, vth_fiducial, xnue0, b2, b3, rho)/sr1   !in cgs
               opac = 0.d0
            else
               !assume helium ionized, hydrogen from ionization balance
               opalbar = get_opalbar3(iline, kline, sr1, yhe1, hei1, temp, vth_fiducial, xnue0, b2, b3, rho)/sr1   !in cgs
               opac = opac_thomson(yhe1, hei1, rho, kcont)               
            endif
!
            sline = sline_depcoeff(xnue0, temp, b2, b3)
            scont = bnue(xnue0,temp)            

!            if(opalbar.ne.opalbar) then
!              write(*,*) opalbar, rho, temp, rcyc, zcyc
!              stop 'error'
!           endif
!            write(*,*) rho, opalbar
!            write(*,*) cs1_r(i), temp, b2, b3, rho, opalbar*sr1
         else
            velphi = zero
            temp = tmin1!0.8d0*(teff1+teff2)/2.d0
            rho = zero
            opalbar = zero
            sline = zero
            scont = zero
            opac = zero         
         endif   
!        
!----------------------line source function-----------------------------
!
         cs1_sline3d(i,j,k)=sline
!
!----------------------continuum source function------------------------
!
         cs1_scont3d(i,j,k)=scont
!
!------------------------velocity components----------------------------
!
         sint1=sin(cs1_theta(j))
         cost1=cos(cs1_theta(j))
         sinp1=sin(cs1_phi(k))
         cosp1=cos(cs1_phi(k))

         cs1_velx3d(i,j,k) = -velphi*sinp1
         cs1_vely3d(i,j,k) = velphi*cosp1
         cs1_velz3d(i,j,k) = zero
!
!---------------------------temperature---------------------------------
!
         cs1_t3d(i,j,k) = temp
!
!-----------------------------opacity-----------------------------------
!
         cs1_rho3d(i,j,k) = rho
         cs1_opalbar3d(i,j,k) = opalbar*sr1   !in units of sr1
         cs1_opac3d(i,j,k) = opac*sr1



         !         write(*,*) cs1_r(i)*rstar1, rho, temp, velr/1.d5, velphi/1.d5, rstar1
!
      enddo
   enddo
enddo
!
!stop 'go on in modelspec'

!j=30
!write(*,*) cs1_theta(j)*180./pi
!k=1
!do i=1, cs1_nr
!   write(*,*) i, cs1_r(i), cs1_opalbar3d(i,j,k)
!enddo
!stop

!
!
write(*,*) '----------------calc_model3df: Be-type disc for star 2 (Dylan)-----------------'
write(*,*)
!
!
!rho0_disc2 = mdisc2*sqrt(cgs_grav*mstar2_cgs)*log(ten)/two/sqrt(pi**3)/csound/sr2**3.5d0
!rho0_disc1 = mdisc1*sqrt(cgs_grav*mstar1_cgs)*log(ten)/two/sqrt(pi**3)/csound/sr1**3.5d0
!write(*,*) rho0_disc1, rho0_disc2
!write(*,*) sr1, sr2
!stop
!
do i=1, cs2_nr
   do j=1, cs2_ntheta
      do k=1, cs2_nphi
!   do j=16, 16!cs2_ntheta
!      do k=1, 1!cs2_nphi
!
         rad2=cs2_r(i)
         theta2=cs2_theta(j)
         sint2=sin(cs2_theta(j))
         cost2=cos(cs2_theta(j))
         sinp2=sin(cs2_phi(k))
         cosp2=cos(cs2_phi(k))

!position in global coordinate system
         xcoord0 = x02 + rad2*sint2*cosp2*rstar2/unit_length
         ycoord0 = y02 + rad2*sint2*sinp2*rstar2/unit_length
         zcoord0 = z02 + rad2*cost2*rstar2/unit_length
!postion in system of star 1
         xcoord1 = (xcoord0-x01)*unit_length/rstar1
         ycoord1 = (ycoord0-y01)*unit_length/rstar1
         zcoord1 = (zcoord0-z01)*unit_length/rstar1
         rad1=sqrt(xcoord1**2+ycoord1**2+zcoord1**2)
         call get_angles_spc(xcoord1, ycoord1, zcoord1, theta1, phi1)         
!
!atmosphere for star 1
         if(rad1.lt.rmin1) then
            velphi = zero
            temp = zero
            vth = zero
            opalbar = zero
            sline = zero
            scont = zero
            opac = zero
         elseif(rad1.gt.rdisc1_max) then
            velphi = zero
            temp = zero
            vth = zero
            opalbar = zero
            sline = zero
            scont = zero
            opac = zero
!            if(rad2.lt.rdisc2_max.and. &
!               theta2.gt.theta_disc2l.and. &
!               theta2.lt.theta_disc2u) then   !set disc 2
!
!               sint2=sin(theta2)
!               velphi = sqrt(cgs_grav*mstar2_cgs/rad2/sr2/sint2) * vrot2fac
!
!               if(opt_temperature2.eq.0) then
!                  temp=tdisc2
!               elseif(opt_temperature2.eq.1) then
!                  temp = tacoeff2 + tbcoeff2*rad2
!               elseif(opt_temperature2.eq.2) then
!                  temp = max(teff2/rad2,5d3)
!               elseif(opt_temperature2.eq.3) then
!                  temp = teff2 * dilfac(rad2)**0.25d0
!               elseif(opt_temperature2.eq.4) then
!                  temp = teff2/pi**0.25d0 * (asin(one/rad2) - sqrt(one-one/(rad2)**2)/rad2)**0.25d0
!                  if(cs2_tauz1d(i).lt.0.15d0) temp = max(temp, 0.6d0*teff2)
!               else
!                  stop 'opt_temperature not properly specified'
!               endif
!               
!               vth = vthermal(vmicro2, temp, na)
!               csound = vsound(temp,mmw1)
!               rho0_disc2 = mdisc2*sqrt(cgs_grav*mstar2_cgs)*log(ten)/two/sqrt(pi**3)/csound/sr2**3.5d0
!               rho = rho0_disc2*(rad2*sint2)**(-slope2) * exp(cgs_grav*mstar2_cgs/csound**2 * (one/sr2/rad2-one/sr2/rad2/sint2))
!!departure coeffcients in LTE (everything calculated for mmw and yhe of star 1)
!               b2=one
!               b3=one
!               opalbar = get_opalbar(iline, kline, sr2, yhe1, hei1, temp, vth_fiducial, xnue0, b2, b3, rho)/sr2   !in cgs
!               sline = sline_depcoeff(xnue0, temp, b2, b3)
!               scont = bnue(xnue0,temp)
!               opac = opac_thomson(yhe1, hei1, rho, kcont)
!!               write(*,*) theta2*180.d0/pi, theta_disc2l*180.d0/pi, theta_disc2u*180.d0/pi
!!               write(*,*) rho, opalbar
!            endif            
         elseif(theta1.gt.theta_disc1l.and.&
              theta1.lt.theta_disc1u) then
!
            sint1 = sin(theta1)
            cost1=cos(theta1)
!calculate corresponding rcyc and zcyc in cylindrical coordinates (and use minimum of 1.d0)
            rcyc=max(rad1*sint1,one)
            zcyc=abs(rad1*cost1)

            call model_adisc(rcyc, zcyc, rstar1, mstar1, alpha1, mdot1, mmw1, temp, rho, velr, velphi)
!            temp = max(temp,tmin1)
!
!departure coeffcients in LTE
            b2=one
            b3=one
            if(temp.lt.6.d3) then
               !assume helium not ionized, hydrogen from ionization balance
               opalbar = get_opalbar3(iline, kline, sr1, yhe1, 0.d0, temp, vth_fiducial, xnue0, b2, b3, rho)/sr1   !in cgs
               opac = 0.d0
            else
               !assume helium ionized, hydrogen from ionization balance
               opalbar = get_opalbar3(iline, kline, sr1, yhe1, hei1, temp, vth_fiducial, xnue0, b2, b3, rho)/sr1   !in cgs
               opac = opac_thomson(yhe1, hei1, rho, kcont)               
            endif
!
            sline = sline_depcoeff(xnue0, temp, b2, b3)
            scont = bnue(xnue0,temp)            
!
!velocity components in system of star 1
            velx1 = -velphi*sin(phi1)
            vely1 = velphi*cos(phi1)
            velz1 = zero
!velocity components in global system
            velx0 = velx1 + vx01
            vely0 = vely1 + vy01
            velz0 = velz1 + vz01
!velocity components in system of star 2         
            velx2 = velx0 - vx02
            vely2 = vely0 - vy02
            velz2 = velz0 - vz02
            velphi = -sinp2*velx2 + cosp2*vely2
!            
         else
            velphi = zero
            temp = tmin1
            rho = zero
            opalbar = zero
            sline = zero
            scont = zero
            opac = zero         
         endif

!         write(*,*) rad1, rad2, opalbar*sr2, theta1, phi1, cs2_theta(j), cs2_phi(k)
!        
!----------------------line source function-----------------------------
!
         cs2_sline3d(i,j,k)=sline
!
!----------------------continuum source function------------------------
!
         cs2_scont3d(i,j,k)=scont
!
!------------------------velocity components----------------------------
!
!velocity components in system of star 2
         sint2=sin(cs2_theta(j))
         cost2=cos(cs2_theta(j))
         sinp2=sin(cs2_phi(k))
         cosp2=cos(cs2_phi(k))

         cs2_velx3d(i,j,k) = -velphi*sinp2
         cs2_vely3d(i,j,k) = velphi*cosp2
         cs2_velz3d(i,j,k) = zero         
!
!---------------------------temperature---------------------------------
!
         cs2_t3d(i,j,k) = temp
!
!-----------------------------opacity-----------------------------------
!
         cs2_rho3d(i,j,k) = rho
         cs2_opalbar3d(i,j,k) = opalbar*sr2  !in units of sr2
         cs2_opac3d(i,j,k) = opac*sr2
!
      enddo
   enddo
enddo
!
!stop 'go on in modelspec'
!
!-----------------------------set velocities------------------------------
!
!maximum velocity in global system
!
vmax=zero
do i=1, cs1_nr
   do j=1, cs1_ntheta
      do k=1, cs1_nphi
         vel=sqrt((vx01+cs1_velx3d(i,j,k))**2+(vy01+cs1_vely3d(i,j,k))**2+(vz01+cs1_velz3d(i,j,k))**2)
         if(vel.gt.vmax) vmax=vel
      enddo
   enddo
enddo

do i=1, cs2_nr
   do j=1, cs2_ntheta
      do k=1, cs2_nphi
         vel=sqrt((vx02+cs2_velx3d(i,j,k))**2+(vy02+cs2_vely3d(i,j,k))**2+(vz02+cs2_velz3d(i,j,k))**2)
         if(vel.gt.vmax) vmax=vel
      enddo
   enddo
enddo
!
!

!zero rotation of primary
vrot2=0.d0
!
!for now, set zero continuum
!cs1_opac3d=0.d0
!cs1_scont3d=0d0
!cs2_opac3d=0.d0
!cs2_scont3d=0.d0
!
end subroutine calc_model3d_adisc
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine calc_model3d_jeta
!
!star 1 is post AGB star
!star 2 has a jet
!
use prog_type
use fund_const
use options_modspec, only: indat_file, input_file, input_file2, input_mod
use dime_modspec, only: cs1_nr, cs1_ntheta, cs1_nphi, cs1_r, cs1_theta, cs1_phi, &
                        cs1_velx3d, cs1_vely3d, cs1_velz3d, cs1_opac3d, &
                        cs1_opalbar3d, cs1_sline3d, cs1_scont3d, cs1_imask3d, cs1_t3d, cs1_rho3d, &
                        cs2_nr, cs2_ntheta, cs2_nphi, cs2_r, cs2_theta, cs2_phi, &
                        cs2_velx3d, cs2_vely3d, cs2_velz3d, cs2_opac3d, &
                        cs2_opalbar3d, cs2_sline3d, cs2_scont3d, cs2_imask3d, cs2_t3d, cs2_rho3d
use params_modspec, only: vx01, vy01, vz01, vx02, vy02, vz02, x01, y01, z01, x02, y02, z02, ex02, ey02, ez02, &
                          vmin1, vmax1, beta1, vrot1, rstar1, rmin1, rmax1, sr1, teff1, yhe1, hei1, logg1, lstar1, mdot1, vmicro1, &
                          vmin2, vmax2, beta2, vrot2, rstar2, rmin2, rmax2, sr2, teff2, yhe2, hei2, logg2, lstar2, mdot2, vmicro2, &
                          vmax, vth_fiducial, trad1, unit_length
use hdf5
use mod_grid, only: grid_log, grid_equi
use mod_opacities, only: opac_thomson, get_opalbar2
use mod_iline, only: iline, na, xnue0, kline
!
implicit none
!
! ... input scalars
real(dp) :: p_rho, logc_n, p_vel, theta_in, theta_out, theta_cav, vel_in, vel_out, tjet, phi_tilt, rmax_jet
namelist / input_usr / p_rho, logc_n, p_vel, theta_in, theta_out, theta_cav, vel_in, vel_out, tjet, phi_tilt, rmax_jet
!
! ... local scalars
integer(i4b) :: opt_temperature1, opt_temperature2
integer(i4b) :: i, j, k, err, idum1, idum2, ii
real(dp) :: vel, opac, opalbar, rho, sline, scont, temp, velphi, velr, velth, b2, b3, n_out
real(dp) :: fac1, fac2, rho0, theta_cav2, theta_out2, theta2_dum
real(dp) :: sint1, cost1, sinp1, cosp1, theta1, phi1, rad1, theta2, phi2, rad2, sint2, cost2, sinp2, cosp2, &
            xcoord0, ycoord0, zcoord0, xcoord1, ycoord1, zcoord1, &
            velx0, vely0, velz0, velx1, vely1, velz1, velx2, vely2, velz2, xcoord2, ycoord2, zcoord2
!
! ... local arrays
real(dp), dimension(:), allocatable :: fdum1_arr, fdum2_arr
!
! ... local characters
!
! ... local functions
real(dp) :: bnue, bvel, sline_depcoeff, vthermal
!
! ... for temperature structure
integer(i4b) :: cs1_np, cs2_np
integer(i4b), parameter :: nz=301
real(dp) :: zmin, zmax, tauz, dtau, rad
real(dp), dimension(:), allocatable :: cs1_p1d, cs1_tauz1d, cs2_p1d, cs2_tauz1d
real(dp), dimension(nz) :: opac1d, z1d
!
!
! ... for hdf5 file
!
write(*,*) '-------------calc_model3d_jeta: AGB-star + jet at secondary (Dylan B.)--------'
write(*,*)
write(*,*) 'reading input from file: ', trim(indat_file)
write(*,*)
!
!-----------------------------------------------------------------------
!
open(1, file=trim(indat_file), status='old', form='formatted')
!read only user input
   rewind 1
   read(1, nml=input_usr)
close(1)
!
!transform velocities to cm/s
vel_in = vel_in*1.d5
vel_out = vel_out*1.d5
!
!transform angles to rad
theta_in = theta_in * pi/180.d0
theta_out = theta_out * pi/180.d0
theta_cav = theta_cav * pi/180.d0
phi_tilt = phi_tilt * pi/180.d0
!
!get angles for opposite hemisphere
theta_cav2 = pi-theta_cav
theta_out2 = pi-theta_out
!
!transform base density to cgs
!number density
n_out = ten**logc_n / 1.d6
rho0 = cgs_mp * n_out !(one + four * yhe2) * n_out
!
!-----------------------------------------------------------------------
!
!coordinate system transformation to incorporate phi_tilt
!warning: for a given tilt, need to calculate line profiles at different phases
!         by input of p_object01 and p_object02 (instead of rotating observer around the system!!!)
ex02 = (/ cos(phi_tilt), zero, -sin(phi_tilt) /)
ez02 = (/ sin(phi_tilt), zero, cos(phi_tilt) /)
ey02 = (/ ex02(2)*ez02(3)-ex02(3)*ez02(2), &
          ex02(1)*ez02(3)-ex02(3)*ez02(1), &
          ex02(1)*ez02(2)-ex02(2)*ez02(1) /)
!
!-----------------------------------------------------------------------
!
!helium ionization
hei1 = 0.d0
hei2 = 0.d0
!
!-----------------------------------------------------------------------
!
!b2 = one
!b3 = one
!temp = tjet
!write(*,*) rho0, iline, kline, sr2, yhe2, hei2, temp, vth_fiducial, xnue0, b2, b3
!vth_fiducial=vthermal(0.d0, temp, 1)
!yhe2=0.d0
!opalbar = get_opalbar2(iline, kline, sr2, yhe2, hei2, temp, vth_fiducial, xnue0, b2, b3, rho0)/sr2   !in cgs     
!write(*,*) opalbar/sqrt(pi), opalbar*xnue0*vth_fiducial/cgs_clight
!stop
!
write(*,*) '-----calc_model3d_jeta: AGB-star + jet at secondary (Dylan B.) for star 1-----'
write(*,*)
!
!
!no need for high resolution here
cs1_nr=5
cs1_ntheta=5
cs1_nphi=2*cs1_ntheta-1
!
!allocate arrays
allocate(cs1_r(cs1_nr), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3df'
allocate(cs1_theta(cs1_ntheta), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3df'
allocate(cs1_phi(cs1_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3df'
allocate(cs1_opac3d(cs1_nr,cs1_ntheta,cs1_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3df'
allocate(cs1_opalbar3d(cs1_nr,cs1_ntheta,cs1_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3df'
allocate(cs1_scont3d(cs1_nr,cs1_ntheta,cs1_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3df'
allocate(cs1_sline3d(cs1_nr,cs1_ntheta,cs1_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3df'
allocate(cs1_velx3d(cs1_nr,cs1_ntheta,cs1_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3df'
allocate(cs1_vely3d(cs1_nr,cs1_ntheta,cs1_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3df'
allocate(cs1_velz3d(cs1_nr,cs1_ntheta,cs1_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3df'
allocate(cs1_t3d(cs1_nr,cs1_ntheta,cs1_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3df'
allocate(cs1_rho3d(cs1_nr,cs1_ntheta,cs1_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3df'
!
!radius in units of sr1
call grid_log(rmin1, rmax1, cs1_nr,  cs1_r)
!
!equidistant theta grid
call grid_equi(zero, pi, cs1_ntheta, cs1_theta)
!
!phi-grid equidistant
call grid_equi(zero, two*pi, cs1_nphi, cs1_phi)
!
!
!high resolution grid for jet
cs2_nr=101
idum1=41
idum2=21
cs2_ntheta = 2*idum1+idum2
cs2_nphi=2*cs1_ntheta-1
!
!allocate arrays
allocate(cs2_r(cs2_nr), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3df'
allocate(cs2_theta(cs2_ntheta), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3df'
allocate(cs2_phi(cs2_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3df'
allocate(cs2_opac3d(cs2_nr,cs2_ntheta,cs2_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3df'
allocate(cs2_opalbar3d(cs2_nr,cs2_ntheta,cs2_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3df'
allocate(cs2_scont3d(cs2_nr,cs2_ntheta,cs2_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3df'
allocate(cs2_sline3d(cs2_nr,cs2_ntheta,cs2_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3df'
allocate(cs2_velx3d(cs2_nr,cs2_ntheta,cs2_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3df'
allocate(cs2_vely3d(cs2_nr,cs2_ntheta,cs2_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3df'
allocate(cs2_velz3d(cs2_nr,cs2_ntheta,cs2_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3df'
allocate(cs2_t3d(cs2_nr,cs2_ntheta,cs2_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3df'
allocate(cs2_rho3d(cs2_nr,cs2_ntheta,cs2_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3df'
!
!radius in units of sr2
!call grid_log(rmin2, rmax2, cs2_nr,  cs2_r)
call grid_equi(rmin2, rmax2, cs2_nr,  cs2_r)
!
if(cs2_r(cs2_nr-1).le.rmax_jet) then
   write(*,*) 'error in jeta: we require always two zero-opacity ghost points at outer boundary'
   stop
endif
!
!theta-grid with small steps near the poles (always ntheta/2 for the disc)
allocate(fdum1_arr(idum1))
call grid_equi(zero, theta_out, idum1, fdum1_arr)
allocate(fdum2_arr(idum2))
call grid_equi(two*theta_out-fdum1_arr(idum1-1), pi - two*theta_out + fdum1_arr(idum1-1), idum2, fdum2_arr)
ii=1
do i=1, idum1
   cs2_theta(ii)=fdum1_arr(i)
   ii=ii+1
enddo
do i=1, idum2
   cs2_theta(ii)=fdum2_arr(i)
   ii=ii+1
enddo
do i=1, idum1
   cs2_theta(ii)= pi - fdum1_arr(idum1+1-i)
   ii=ii+1
enddo
cs2_theta(1)=zero
cs2_theta(cs2_ntheta)=pi
!
!call grid_equi(zero, pi, cs1_ntheta, cs1_theta)
!
!phi-grid equidistant
call grid_equi(zero, two*pi, cs2_nphi, cs2_phi)
!
!-----------------------------------------------------------------------
!
do i=1, cs1_nr
   do j=1, cs1_ntheta
      do k=1, cs1_nphi
!   do j=33, 33!cs1_ntheta
!      do k=1, 1!cs1_nphi

!zero values if inside secondary         
         rad1=cs1_r(i)
         sint1=sin(cs1_theta(j))
         cost1=cos(cs1_theta(j))
         sinp1=sin(cs1_phi(k))
         cosp1=cos(cs1_phi(k))
!position in global coordinate system
         xcoord0 = x01 + rad1*sint1*cosp1*rstar1/unit_length
         ycoord0 = y01 + rad1*sint1*sinp1*rstar1/unit_length
         zcoord0 = z01 + rad1*cost1*rstar1/unit_length
!postion in system of star 2
         xcoord2 = (xcoord0-x02)*unit_length/rstar2
         ycoord2 = (ycoord0-y02)*unit_length/rstar2
         zcoord2 = (zcoord0-z02)*unit_length/rstar2
         rad2=sqrt(xcoord2**2+ycoord2**2+zcoord2**2)
         call get_angles_spc(xcoord2, ycoord2, zcoord2, theta2, phi2)
!
!default values
         temp = zero
         velr = zero
         velth = zero
         velphi = zero
         velx1 = zero
         vely1 = zero
         velz1 = zero
         rho = zero
         opac = zero            
         opalbar = zero
         scont = zero            
         sline = zero
!
         if(rad2.lt.rmin2) then
!if inside core of star 2
            velr = zero
            velth = zero
            velphi = zero
            temp = zero
            rho = zero
            opac = zero            
            opalbar = zero
            scont = zero            
            sline = zero
         elseif(rad2.lt.rmax_jet) then
!if inside computational domain of star 2            
            if(theta2.le.theta_out) then
!if inside of the jet
!densities
               rho = rho0 * (theta2/theta_out)**p_rho / (zcoord2*sr2/cgs_au)**2
!velocity components in system of star 2
               fac2 = ((theta2-theta_cav)/(theta_out-theta_cav))**2
               fac1 = (exp(-p_vel*fac2)-exp(-p_vel))/(one-exp(-p_vel))
               velr = vel_out + (vel_in-vel_out)*fac1
               velth = zero
               velphi = zero
               velx2 = velr*sint2*cosp2 + velth*cost2*cosp2 - velphi*sinp2
               vely2 = velr*sint2*sinp2 + velth*cost2*sinp2 + velphi*cosp2
               velz2 = velr*cost2 - velth*sint2
!velocity components in global system
               velx0 = velx2 + vx02
               vely0 = vely2 + vy02
               velz0 = velz2 + vz02
!velocity components in system of star 1
               velx1 = velx0 - vx01
               velx1 = vely0 - vy01
               velz1 = velz0 - vz01
!temperatures and turbulent velocities
               temp = tjet
!departure coefficients (here: LTE)
               b2 = one
               b3 = one
!opacities and source functions
               opac = zero
               opalbar = get_opalbar2(iline, kline, sr2, yhe2, hei2, temp, vth_fiducial, xnue0, b2, b3, rho)/sr2   !in cgs     
               scont = zero
               sline = zero !sline_depcoeff(xnue0, temp, b2, b3)
            endif
         endif
!        
!----------------------line source function-----------------------------
!
         cs1_sline3d(i,j,k)=sline
!
!----------------------continuum source function------------------------
!
         cs1_scont3d(i,j,k)=scont
!
!------------------------velocity components----------------------------
!
         cs1_velx3d(i,j,k) = velx1
         cs1_vely3d(i,j,k) = vely1
         cs1_velz3d(i,j,k) = velz1
!
!---------------------------temperature---------------------------------
!
         cs1_t3d(i,j,k) = temp
!
!-----------------------------opacity-----------------------------------
!
         cs1_rho3d(i,j,k) = rho
         cs1_opalbar3d(i,j,k) = opalbar*sr1   !in units of sr1
         cs1_opac3d(i,j,k) = opac*sr1
!
      enddo
   enddo
enddo
!
!
!
write(*,*) '-----calc_model3d_jeta: AGB-star + jet at secondary (Dylan B.) for star 2-----'
write(*,*)
!
!
!
do i=1, cs2_nr
   do j=1, cs2_ntheta
      do k=1, cs2_nphi
!
         rad2=cs2_r(i)
         theta2=cs2_theta(j)
         sint2=sin(cs2_theta(j))
         cost2=cos(cs2_theta(j))
         sinp2=sin(cs2_phi(k))
         cosp2=cos(cs2_phi(k))
!
!cartesian coordinates in system 2
         xcoord2 = rad2*sint2*cosp2
         ycoord2 = rad2*sint2*sinp2
         zcoord2 = rad2*cost2
!position in global coordinate system
         xcoord0 = x02 + rad2*sint2*cosp2*rstar2/unit_length
         ycoord0 = y02 + rad2*sint2*sinp2*rstar2/unit_length
         zcoord0 = z02 + rad2*cost2*rstar2/unit_length
!postion in system of star 1
         xcoord1 = (xcoord0-x01)*unit_length/rstar1
         ycoord1 = (ycoord0-y01)*unit_length/rstar1
         zcoord1 = (zcoord0-z01)*unit_length/rstar1
         rad1=sqrt(xcoord1**2+ycoord1**2+zcoord1**2)
         call get_angles_spc(xcoord1, ycoord1, zcoord1, theta1, phi1)         
!
!default values
         temp = zero
         velr = zero
         velth = zero
         velphi = zero
         rho = zero
         opac = zero            
         opalbar = zero
         scont = zero            
         sline = zero
!
         if(rad1.lt.rmin1) then
!if inside core of star 1
            velr = zero
            velth = zero
            velphi = zero
            temp = zero
            rho = zero
            opac = zero            
            opalbar = zero
            scont = zero            
            sline = zero
!         elseif(rad1.lt.rmin1) then
!if inside computational domain of star 1            
         elseif(rad2.lt.rmax_jet) then
!if inside computational domain of star 2            
            if(theta2.le.theta_out) then
!if inside of the jet upper hemisphere
!densities
               rho = rho0 * (theta2/theta_out)**p_rho / (zcoord2*sr2/cgs_au)**2
!velocities
               fac2 = ((theta2-theta_cav)/(theta_out-theta_cav))**2
               fac1 = (exp(-p_vel*fac2)-exp(-p_vel))/(one-exp(-p_vel))
               velr = vel_out + (vel_in-vel_out)*fac1
               velth = zero
               velphi = zero
!temperatures and turbulent velocities
               temp = tjet
!departure coefficients (here: LTE)
               b2 = one
               b3 = one
!opacities and source functions
               opac = zero
               opalbar = get_opalbar2(iline, kline, sr2, yhe2, hei2, temp, vth_fiducial, xnue0, b2, b3, rho)/sr2   !in cgs     
               scont = zero
               sline = zero!sline_depcoeff(xnue0, temp, b2, b3)
            endif           
            if(theta2.ge.theta_out2) then
!if inside of the jet lower hemisphere
               theta2_dum = pi - theta2               
!densities
               rho = rho0 * (theta2_dum/theta_out)**p_rho / (zcoord2*sr2/cgs_au)**2
!velocities
               fac2 = ((theta2_dum-theta_cav)/(theta_out-theta_cav))**2
               fac1 = (exp(-p_vel*fac2)-exp(-p_vel))/(one-exp(-p_vel))
               velr = vel_out + (vel_in-vel_out)*fac1
               velth = zero
               velphi = zero
!temperatures and turbulent velocities
               temp = tjet
!departure coefficients (here: LTE)
               b2 = one
               b3 = one
!opacities and source functions
               opac = zero
               opalbar = get_opalbar2(iline, kline, sr2, yhe2, hei2, temp, vth_fiducial, xnue0, b2, b3, rho)/sr2   !in cgs     
               scont = zero
               sline = zero!sline_depcoeff(xnue0, temp, b2, b3)
            endif
         endif
!        
!----------------------line source function-----------------------------
!
         cs2_sline3d(i,j,k)=sline
!
!----------------------continuum source function------------------------
!
         cs2_scont3d(i,j,k)=scont
!
!------------------------velocity components----------------------------
!
!velocity components in system of star 2
         cs2_velx3d(i,j,k) = velr*sint2*cosp2 + velth*cost2*cosp2 - velphi*sinp2
         cs2_vely3d(i,j,k) = velr*sint2*sinp2 + velth*cost2*sinp2 + velphi*cosp2
         cs2_velz3d(i,j,k) = velr*cost2 - velth*sint2
!
!---------------------------temperature---------------------------------
!
         cs2_t3d(i,j,k) = temp
!
!-----------------------------opacity-----------------------------------
!
         cs2_rho3d(i,j,k) = rho
         cs2_opalbar3d(i,j,k) = opalbar*sr2  !in units of sr2
         cs2_opac3d(i,j,k) = opac*sr2
!
      enddo
   enddo
enddo
!
!k=1
!j=41
!do i=1, cs2_nr
!   write(*,*) i, cs2_theta(j)*180./pi, cs2_r(i)*sr2/cgs_au, cs2_rho3d(i,j,k), rho0/(cs2_r(i)*cos(cs2_theta(j))*sr2/cgs_au)**2, rho0
!enddo!
!
!write(*,*) 'go on in jeta'
!stop
!
!-----------------------------set velocities------------------------------
!
!maximum velocity in global system
!
vmax=zero
do i=1, cs1_nr
   do j=1, cs1_ntheta
      do k=1, cs1_nphi
         vel=sqrt((vx01+cs1_velx3d(i,j,k))**2+(vy01+cs1_vely3d(i,j,k))**2+(vz01+cs1_velz3d(i,j,k))**2)
         if(vel.gt.vmax) vmax=vel
      enddo
   enddo
enddo

do i=1, cs2_nr
   do j=1, cs2_ntheta
      do k=1, cs2_nphi
         vel=sqrt((vx02+cs2_velx3d(i,j,k))**2+(vy02+cs2_vely3d(i,j,k))**2+(vz02+cs2_velz3d(i,j,k))**2)
         if(vel.gt.vmax) vmax=vel
      enddo
   enddo
enddo
!
!
end subroutine calc_model3d_jeta
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine calc_model3d_jetb
!
!star 1 is post AGB star
!star 2 has a jet from Oliviers model
!
use prog_type
use fund_const
use options_modspec, only: indat_file, input_file, input_file2, input_mod
use dime_modspec, only: cs1_nr, cs1_ntheta, cs1_nphi, cs1_r, cs1_theta, cs1_phi, &
                        cs1_velx3d, cs1_vely3d, cs1_velz3d, cs1_opac3d, &
                        cs1_opalbar3d, cs1_sline3d, cs1_scont3d, cs1_imask3d, cs1_t3d, cs1_rho3d, &
                        cs2_nr, cs2_ntheta, cs2_nphi, cs2_r, cs2_theta, cs2_phi, &
                        cs2_velx3d, cs2_vely3d, cs2_velz3d, cs2_opac3d, &
                        cs2_opalbar3d, cs2_sline3d, cs2_scont3d, cs2_imask3d, cs2_t3d, cs2_rho3d
use params_modspec, only: vx01, vy01, vz01, vx02, vy02, vz02, x01, y01, z01, x02, y02, z02, ex02, ey02, ez02, &
                          vmin1, vmax1, beta1, vrot1, rstar1, rmin1, rmax1, sr1, teff1, yhe1, hei1, logg1, lstar1, mdot1, vmicro1, &
                          vmin2, vmax2, beta2, vrot2, rstar2, rmin2, rmax2, sr2, teff2, yhe2, hei2, logg2, lstar2, mdot2, vmicro2, &
                          vmax, vth_fiducial, trad1, unit_length
use hdf5
use mod_grid, only: grid_log, grid_equi
use mod_opacities, only: opac_thomson, get_opalbar2
use mod_iline, only: iline, na, xnue0, kline
use mod_interp1d, only: find_index
use mod_interp2d, only: interpol2d_4p_lin
!
implicit none
!
! ... input scalars
!
! ... local scalars
integer(i4b) :: opt_temperature1, opt_temperature2
integer(i4b) :: i, j, k, err, idum1, idum2
real(dp) :: vel, opac, opalbar, rho, sline, scont, temp, velphi, velr, velth, b2, b3, n_out
real(dp) :: fac1, fac2, rho0, theta_cav2, theta_out2, theta2_dum
real(dp) :: sint1, cost1, sinp1, cosp1, theta1, phi1, rad1, theta2, phi2, rad2, sint2, cost2, sinp2, cosp2, &
            xcoord0, ycoord0, zcoord0, xcoord1, ycoord1, zcoord1, &
            velx0, vely0, velz0, velx1, vely1, velz1, velx2, vely2, velz2, xcoord2, ycoord2, zcoord2
!
! ... local arrays
real(dp), dimension(:), allocatable :: fdum1_arr, fdum2_arr
!
! ... local characters
!
! ... local functions
real(dp) :: bnue, bvel, sline_depcoeff, vthermal
!
! ... for interpolation
integer(i4b) :: iim2, iim1, ii, iip1, kkm2, kkm1, kk, kkp1
real(dp) :: rcoord_cyc, zcoord_cyc, velr_cyc, velz_cyc, velphi_cyc
!

! ... for the model
integer(i4b) :: nr_cyc, nz_cyc
logical :: lcheck
real(dp), dimension(:), allocatable :: r_cyc, z_cyc
real(dp), dimension(:,:), allocatable :: rho2d_cyc, velr2d_cyc, velz2d_cyc, velphi2d_cyc, t2d_cyc

! ... for hdf5-files
integer(hid_t) :: file_id, group_id, dset_id, attr_id
integer(hsize_t), dimension(1) :: dims_scalars
integer(hsize_t), dimension(1) :: dims_r
integer(hsize_t), dimension(1) :: dims_z
integer(hsize_t), dimension(2) :: dims2d
!
! ... namelist
character(len=500)  :: fname_model
real(dp) :: model_unitl, model_unitd, model_unitv, model_unitt, tjet
namelist / input_usr / fname_model, model_unitl, model_unitd, model_unitv, model_unitt, cs2_nphi, tjet
!
! ... for hdf5 file
!
write(*,*) '-------------calc_model3d_jetb: AGB-star + jet at secondary (Olivier V.)------'
write(*,*)
write(*,*) 'reading input from file: ', trim(indat_file)
write(*,*)
!
!-----------------------------------------------------------------------
!
open(1, file=trim(indat_file), status='old', form='formatted')
!read only user input
   rewind 1
   read(1, nml=input_usr)
close(1)
!
!---------------------read in data--------------------------------------
!
fname_model = adjustl(fname_model)
inquire(file=trim(fname_model), exist=lcheck)
if(.not.lcheck) then
   write(*,*) 'error in calc_model3d_jetb: file "', trim(fname_model), '" does not exist'
   stop
endif
!
!
!
write(*,*) 'file name: ', trim(fname_model)
write(*,*)
!
call h5open_f(err)
call h5fopen_f(trim(fname_model), h5f_acc_rdonly_f, file_id, err)
!
call h5gopen_f(file_id, 'dimensions', group_id, err)
   call h5aopen_f(group_id, 'nx', attr_id, err)
      call h5aread_f(attr_id, h5t_native_integer, nr_cyc, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5aopen_f(group_id, 'nz', attr_id, err)
      call h5aread_f(attr_id, h5t_native_integer, nz_cyc, dims_scalars, err)
   call h5aclose_f(attr_id, err)
call h5gclose_f(group_id, err)
!
!-------------------------get coordiantes-------------------------------
!
allocate(r_cyc(nr_cyc), stat=err)
allocate(z_cyc(nz_cyc), stat=err)
allocate(velr2d_cyc(nr_cyc, nz_cyc), stat=err)
allocate(velz2d_cyc(nr_cyc, nz_cyc), stat=err)
allocate(velphi2d_cyc(nr_cyc, nz_cyc), stat=err)
allocate(rho2d_cyc(nr_cyc, nz_cyc), stat=err)
allocate(t2d_cyc(nr_cyc, nz_cyc), stat=err)
!

dims_r= (/ nr_cyc /)
dims_z= (/ nz_cyc /)
dims2d=(/ nr_cyc, nz_cyc /)
!
call h5gopen_f(file_id, 'coordinates', group_id, err)
   call h5dopen_f(group_id, 'x', dset_id, err)
   call h5dread_f(dset_id, h5t_native_double, r_cyc, dims_r, err)
   call h5dclose_f(dset_id, err)
   call h5dopen_f(group_id, 'z', dset_id, err)
      call h5dread_f(dset_id, h5t_native_double, z_cyc, dims_z, err)
   call h5dclose_f(dset_id, err)
call h5gclose_f(group_id, err)
!
call h5gopen_f(file_id, 'model', group_id, err)
   call h5dopen_f(group_id, 'rho', dset_id, err)
      call h5dread_f(dset_id, h5t_native_double, rho2d_cyc, dims2d, err)
   call h5dclose_f(dset_id, err)
   call h5dopen_f(group_id, 'temperature', dset_id, err)
      call h5dread_f(dset_id, h5t_native_double, t2d_cyc, dims2d, err)
   call h5dclose_f(dset_id, err)
   call h5dopen_f(group_id, 'velr', dset_id, err)
      call h5dread_f(dset_id, h5t_native_double, velr2d_cyc, dims2d, err)
   call h5dclose_f(dset_id, err)
   call h5dopen_f(group_id, 'velz', dset_id, err)
      call h5dread_f(dset_id, h5t_native_double, velz2d_cyc, dims2d, err)
   call h5dclose_f(dset_id, err)
   call h5dopen_f(group_id, 'velphi', dset_id, err)
      call h5dread_f(dset_id, h5t_native_double, velphi2d_cyc, dims2d, err)
   call h5dclose_f(dset_id, err)

call h5gclose_f(group_id, err)
call h5fclose_f(file_id, err)
call h5close_f(err)
!
!
!transform everything to cgs
velr2d_cyc = velr2d_cyc * model_unitv
velz2d_cyc = velz2d_cyc * model_unitv
velphi2d_cyc = velphi2d_cyc * model_unitv
rho2d_cyc = rho2d_cyc * model_unitd
t2d_cyc = tjet * model_unitt !t2d_cyc * model_unitt
!
!transform length scales to sr2
r_cyc = r_cyc * model_unitl/sr2
z_cyc = z_cyc * model_unitl/sr2
!
!-----------------------------------------------------------------------
!
!helium ionization
hei1 = 0.d0
hei2 = 0.d0
!
!-----------------------------------------------------------------------
!
write(*,*) '-----calc_model3d_jeta: AGB-star + jet at secondary (Olivier.) for star 1-----'
write(*,*)
!
!
!no need for high resolution here
cs1_nr=5
cs1_ntheta=5
cs1_nphi=2*cs1_ntheta-1
!
!allocate arrays
allocate(cs1_r(cs1_nr), stat=err)
   if(err.ne.0) stop 'allocation error calc_model_jetb'
allocate(cs1_theta(cs1_ntheta), stat=err)
   if(err.ne.0) stop 'allocation error calc_model_jetb'
allocate(cs1_phi(cs1_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model_jetb'
allocate(cs1_opac3d(cs1_nr,cs1_ntheta,cs1_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model_jetb'
allocate(cs1_opalbar3d(cs1_nr,cs1_ntheta,cs1_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model_jetb'
allocate(cs1_scont3d(cs1_nr,cs1_ntheta,cs1_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model_jetb'
allocate(cs1_sline3d(cs1_nr,cs1_ntheta,cs1_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model_jetb'
allocate(cs1_velx3d(cs1_nr,cs1_ntheta,cs1_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model_jetb'
allocate(cs1_vely3d(cs1_nr,cs1_ntheta,cs1_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model_jetb'
allocate(cs1_velz3d(cs1_nr,cs1_ntheta,cs1_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model_jetb'
allocate(cs1_t3d(cs1_nr,cs1_ntheta,cs1_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model_jetb'
allocate(cs1_rho3d(cs1_nr,cs1_ntheta,cs1_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model_jetb'
!
!radius in units of sr1
call grid_log(rmin1, rmax1, cs1_nr,  cs1_r)
   
!
!equidistant theta grid
call grid_equi(zero, pi, cs1_ntheta, cs1_theta)
!
!phi-grid equidistant
call grid_equi(zero, two*pi, cs1_nphi, cs1_phi)
!
!
!high resolution for jet
cs2_nr = 401  !nr_cyc
cs2_ntheta = 201
!cs2_nphi = 2*cs2_ntheta-1
!
!allocate arrays
allocate(cs2_r(cs2_nr), stat=err)
   if(err.ne.0) stop 'allocation error calc_model_jetb'
allocate(cs2_theta(cs2_ntheta), stat=err)
   if(err.ne.0) stop 'allocation error calc_model_jetb'
allocate(cs2_phi(cs2_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model_jetb'
allocate(cs2_opac3d(cs2_nr,cs2_ntheta,cs2_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model_jetb'
allocate(cs2_opalbar3d(cs2_nr,cs2_ntheta,cs2_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model_jetb'
allocate(cs2_scont3d(cs2_nr,cs2_ntheta,cs2_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model_jetb'
allocate(cs2_sline3d(cs2_nr,cs2_ntheta,cs2_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model_jetb'
allocate(cs2_velx3d(cs2_nr,cs2_ntheta,cs2_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model_jetb'
allocate(cs2_vely3d(cs2_nr,cs2_ntheta,cs2_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model_jetb'
allocate(cs2_velz3d(cs2_nr,cs2_ntheta,cs2_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model_jetb'
allocate(cs2_t3d(cs2_nr,cs2_ntheta,cs2_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model_jetb'
allocate(cs2_rho3d(cs2_nr,cs2_ntheta,cs2_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model_jetb'
!
!radius in units of sr2
!call grid_log(rmin2, rmax2, cs2_nr,  cs2_r)
call grid_equi(rmin2, rmax2, cs2_nr,  cs2_r)
!
if(cs2_r(cs2_nr-1).le.sqrt(r_cyc(nr_cyc)**2 + z_cyc(nz_cyc)**2)) then
   write(*,*) 'error in jetb: we require always two zero-opacity ghost points at outer boundary'
   stop
endif
!
!phi and theta grid equidistant
call grid_equi(zero, pi, cs2_ntheta, cs2_theta)
call grid_equi(zero, two*pi, cs2_nphi, cs2_phi)
!
!-----------------------------------------------------------------------
!
do i=1, cs1_nr
   do j=1, cs1_ntheta
      do k=1, cs1_nphi
!
         rad1=cs1_r(i)
         sint1=sin(cs1_theta(j))
         cost1=cos(cs1_theta(j))
         sinp1=sin(cs1_phi(k))
         cosp1=cos(cs1_phi(k))
!position in global coordinate system
         xcoord0 = x01 + rad1*sint1*cosp1*rstar1/unit_length
         ycoord0 = y01 + rad1*sint1*sinp1*rstar1/unit_length
         zcoord0 = z01 + rad1*cost1*rstar1/unit_length
!postion in system of star 2
         xcoord2 = (xcoord0-x02)*unit_length/rstar2
         ycoord2 = (ycoord0-y02)*unit_length/rstar2
         zcoord2 = (zcoord0-z02)*unit_length/rstar2
         rad2=sqrt(xcoord2**2+ycoord2**2+zcoord2**2)
         call get_angles_spc(xcoord2, ycoord2, zcoord2, theta2, phi2)
!
!default values
         temp = zero
         velr = zero
         velth = zero
         velphi = zero
         velx1 = zero
         vely1 = zero
         velz1 = zero
         rho = zero
         opac = zero            
         opalbar = zero
         scont = zero            
         sline = zero
!
!         if(rad2.lt.rmin2) then
!if inside core of star 2
!            velr = zero
!            velth = zero
!            velphi = zero
!            temp = zero
!            rho = zero
!            opac = zero            
!            opalbar = zero
!            scont = zero            
!            sline = zero
!         elseif(rad2.lt.sqrt(r_cyc(nr_cyc)**2 + z_cyc(nz_cyc)**2) then
!if inside computational domain of star 2            
!         endif
!        
!----------------------line source function-----------------------------
!
         cs1_sline3d(i,j,k)=sline
!
!----------------------continuum source function------------------------
!
         cs1_scont3d(i,j,k)=scont
!
!------------------------velocity components----------------------------
!
         cs1_velx3d(i,j,k) = velx1
         cs1_vely3d(i,j,k) = vely1
         cs1_velz3d(i,j,k) = velz1
!
!---------------------------temperature---------------------------------
!
         cs1_t3d(i,j,k) = temp
!
!-----------------------------opacity-----------------------------------
!
         cs1_rho3d(i,j,k) = rho
         cs1_opalbar3d(i,j,k) = opalbar*sr1   !in units of sr1
         cs1_opac3d(i,j,k) = opac*sr1
!
      enddo
   enddo
enddo
!
!
!
write(*,*) '-----calc_model3d_jeta: AGB-star + jet at secondary (Olivier.) for star 2-----'
write(*,*)
!
!
!
do i=1, cs2_nr
   do j=1, cs2_ntheta
      do k=1, cs2_nphi
!
         rad2=cs2_r(i)
         theta2=cs2_theta(j)
         sint2=sin(cs2_theta(j))
         cost2=cos(cs2_theta(j))
         sinp2=sin(cs2_phi(k))
         cosp2=cos(cs2_phi(k))
!
!cartesian coordinates in system 2
         xcoord2 = rad2*sint2*cosp2
         ycoord2 = rad2*sint2*sinp2
         zcoord2 = rad2*cost2
!
!coordinates in cylindrical system of input grid
         rcoord_cyc = rad2*sint2
         zcoord_cyc = rad2*cost2         

!position in global coordinate system
         xcoord0 = x02 + rad2*sint2*cosp2*rstar2/unit_length
         ycoord0 = y02 + rad2*sint2*sinp2*rstar2/unit_length
         zcoord0 = z02 + rad2*cost2*rstar2/unit_length
!postion in system of star 1
         xcoord1 = (xcoord0-x01)*unit_length/rstar1
         ycoord1 = (ycoord0-y01)*unit_length/rstar1
         zcoord1 = (zcoord0-z01)*unit_length/rstar1
         rad1=sqrt(xcoord1**2+ycoord1**2+zcoord1**2)
         call get_angles_spc(xcoord1, ycoord1, zcoord1, theta1, phi1)         
!
!default values
         temp = zero
         velr = zero
         velth = zero
         velphi = zero
         rho = zero
         opac = zero            
         opalbar = zero
         scont = zero            
         sline = zero
!
         if(rad1.lt.rmin1) then
!if inside core of star 1
            velr = zero
            velth = zero
            velphi = zero
            temp = zero
            rho = zero
            opac = zero            
            opalbar = zero
            scont = zero            
            sline = zero
         elseif(rcoord_cyc.lt.r_cyc(nr_cyc).and. &
!                zcoord_cyc.gt.z_cyc(1).and.&
                zcoord_cyc.ge.zero.and.&              
                zcoord_cyc.lt.z_cyc(nz_cyc)) then
            !if inside computational domain of star 2
            !
            !find indices for 2d interpolation
            call find_index(rcoord_cyc, r_cyc, nr_cyc, iim2, iim1, ii, iip1)
            call find_index(zcoord_cyc, z_cyc, nz_cyc, kkm2, kkm1, kk, kkp1)

!if interpolated in log-space, take care about zeros
            rho = interpol2d_4p_lin(rho2d_cyc(iim1,kkm1), rho2d_cyc(ii,kkm1), &
                                    rho2d_cyc(iim1, kk), rho2d_cyc(ii,kk), &
                                    r_cyc(iim1), r_cyc(ii), z_cyc(kkm1), z_cyc(kk), &
                                    rcoord_cyc, zcoord_cyc)

            temp = interpol2d_4p_lin(t2d_cyc(iim1,kkm1), t2d_cyc(ii,kkm1), &
                                    t2d_cyc(iim1, kk), t2d_cyc(ii,kk), &
                                    r_cyc(iim1), r_cyc(ii), z_cyc(kkm1), z_cyc(kk), &
                                    rcoord_cyc, zcoord_cyc)

            velr_cyc = interpol2d_4p_lin(velr2d_cyc(iim1,kkm1), velr2d_cyc(ii,kkm1), &
                                    velr2d_cyc(iim1, kk), velr2d_cyc(ii,kk), &
                                    r_cyc(iim1), r_cyc(ii), z_cyc(kkm1), z_cyc(kk), &
                                    rcoord_cyc, zcoord_cyc)

            velz_cyc = interpol2d_4p_lin(velz2d_cyc(iim1,kkm1), velz2d_cyc(ii,kkm1), &
                                    velz2d_cyc(iim1, kk), velz2d_cyc(ii,kk), &
                                    r_cyc(iim1), r_cyc(ii), z_cyc(kkm1), z_cyc(kk), &
                                         rcoord_cyc, zcoord_cyc)

            velphi = interpol2d_4p_lin(velphi2d_cyc(iim1,kkm1), velphi2d_cyc(ii,kkm1), &
                                       velphi2d_cyc(iim1, kk), velphi2d_cyc(ii,kk), &
                                       r_cyc(iim1), r_cyc(ii), z_cyc(kkm1), z_cyc(kk), &
                                       rcoord_cyc, zcoord_cyc)

            velr = velr_cyc*sint2 + velz_cyc*cost2
            velth = velr_cyc*cost2 - velz_cyc*sint2
            velphi = velphi  !try counter-rotation

            !write(*,*) velr_cyc**2 + velz_cyc**2 + velphi**2, velr**2+velth**2+velphi**2            
!
!departure coefficients (here: LTE)
               b2 = one
               b3 = one
!opacities and source functions
               opac = zero
               opalbar = get_opalbar2(iline, kline, sr2, yhe2, hei2, temp, vth_fiducial, xnue0, b2, b3, rho)/sr2   !in cgs     
               scont = zero
               sline = zero !sline_depcoeff(xnue0, temp, b2, b3)
               
!
         endif
!        
!----------------------line source function-----------------------------
!
         cs2_sline3d(i,j,k)=sline
!
!----------------------continuum source function------------------------
!
         cs2_scont3d(i,j,k)=scont
!
!------------------------velocity components----------------------------
!
!velocity components in system of star 2
         cs2_velx3d(i,j,k) = velr*sint2*cosp2 + velth*cost2*cosp2 - velphi*sinp2
         cs2_vely3d(i,j,k) = velr*sint2*sinp2 + velth*cost2*sinp2 + velphi*cosp2
         cs2_velz3d(i,j,k) = velr*cost2 - velth*sint2
!
!---------------------------temperature---------------------------------
!
         cs2_t3d(i,j,k) = temp
!
!-----------------------------opacity-----------------------------------
!
         cs2_rho3d(i,j,k) = rho
         cs2_opalbar3d(i,j,k) = opalbar*sr2  !in units of sr2
         cs2_opac3d(i,j,k) = opac*sr2

!         if(i.eq.23.and.(k.eq.1.or.k.eq.87)) then
!            write(*,*) k, j, cs2_theta(j)*180./pi, cs2_rho3d(i,j,k), velr/1.d5, velth/1.d5, velphi/1.d5
!         endif
!
      enddo
   enddo
enddo
!
!do j=1, cs2_ntheta
!   write(*,*) j, cs2_theta(j)*180./pi
!enddo
!
!do k=1, cs2_nphi
!   write(*,*) k, cs2_phi(k)*180./pi
!enddo
!stop
!j=51
!k=51
!do i=1, cs2_nr
!  write(*,*) i, cs2_theta(j)*180./pi, cs2_r(i), log10(cs2_rho3d(i,j,k)), cs2_opalbar3d(i,j,k)
!enddo!
!stop
!i=23
!do j=1, cs2_ntheta
!   write(*,*) j, cs2_theta(j)*180./pi, cs2_rho3d(i,j,k), cs2_velr3d(i,j,k), cs2_velth3d(i,j,k), cs2_velphi3d(i,j,k)
!enddo
!
!write(*,*) 'go on in jetb'
!stop
!
!-----------------------------set velocities------------------------------
!
!maximum velocity in global system
!
vmax=zero
do i=1, cs1_nr
   do j=1, cs1_ntheta
      do k=1, cs1_nphi
         vel=sqrt((vx01+cs1_velx3d(i,j,k))**2+(vy01+cs1_vely3d(i,j,k))**2+(vz01+cs1_velz3d(i,j,k))**2)
         if(vel.gt.vmax) vmax=vel
      enddo
   enddo
enddo

do i=1, cs2_nr
   do j=1, cs2_ntheta
      do k=1, cs2_nphi
         vel=sqrt((vx02+cs2_velx3d(i,j,k))**2+(vy02+cs2_vely3d(i,j,k))**2+(vz02+cs2_velz3d(i,j,k))**2)
         if(vel.gt.vmax) vmax=vel
      enddo
   enddo
enddo
!
!
end subroutine calc_model3d_jetb
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine calc_model3d_test
!
!some testing
!
use prog_type
use fund_const
use options_modspec, only: indat_file, input_file, input_file2, input_mod
use dime_modspec, only: cs1_nr, cs1_ntheta, cs1_nphi, cs1_r, cs1_theta, cs1_phi, &
                        cs1_velx3d, cs1_vely3d, cs1_velz3d, cs1_opac3d, &
                        cs1_opalbar3d, cs1_sline3d, cs1_scont3d, cs1_imask3d, cs1_t3d, cs1_rho3d, &
                        cs2_nr, cs2_ntheta, cs2_nphi, cs2_r, cs2_theta, cs2_phi, &
                        cs2_velx3d, cs2_vely3d, cs2_velz3d, cs2_opac3d, &
                        cs2_opalbar3d, cs2_sline3d, cs2_scont3d, cs2_imask3d, cs2_t3d, cs2_rho3d
use params_modspec, only: vx01, vy01, vz01, vx02, vy02, vz02, x01, y01, z01, x02, y02, z02, ex02, ey02, ez02, &
                          vmin1, vmax1, beta1, vrot1, rstar1, rmin1, rmax1, sr1, teff1, yhe1, hei1, logg1, lstar1, mdot1, vmicro1, &
                          vmin2, vmax2, beta2, vrot2, rstar2, rmin2, rmax2, sr2, teff2, yhe2, hei2, logg2, lstar2, mdot2, vmicro2, &
                          vmax, vth_fiducial, trad1, unit_length
use hdf5
use mod_opacities, only: get_opalbar, get_opalbar2, get_opalbar3
use mod_iline, only: iline, na, xnue0, kline
!
implicit none
!
! ... input scalars
!
! ... local scalars
integer(i4b) :: i
integer(i4b), parameter :: nt=100
real(dp) :: rho0, b2, b3
!
! ... local arrays
real(dp), dimension(nt) :: temp1d, opalbar1d_a, opalbar1d_b, opalbar1d_c
!
! ... local characters
!
! ... local functions
!
! ... for temperature structure
!
! ... for hdf5 file
!
write(*,*) '------------------calc_model3d_test: only some testing------------------------'
write(*,*)
!
!-----------------------------------------------------------------------
!
!helium ionization
hei1 = 0.d0
hei2 =0.d0
!
yhe1 = 0.1d0
yhe2 = 0.1d0
!
sr1 = rsu
sr2 = rsu
!
iline=1
kline=1.d0
!
vth_fiducial=1.d7
!
xnue0= 4.5680294d14
!
rho0=1.d-8
rho0=1.d-10
rho0=1.d-12
rho0=1.d-14
!
b2 = one
b3 = one
!
!-----------------------------------------------------------------------
!
do i=1, nt
   temp1d(i) = 1.d3 + (i-1)*(20.d3-1.d3)/(nt-1)
   opalbar1d_a(i) = get_opalbar(iline, kline, sr1, yhe1, hei1, temp1d(i), vth_fiducial, xnue0, b2, b3, rho0)/sr1   !in cgs
   opalbar1d_b(i) = get_opalbar2(iline, kline, sr1, yhe1, hei1, temp1d(i), vth_fiducial, xnue0, b2, b3, rho0)/sr1   !in cgs
   opalbar1d_c(i) = get_opalbar3(iline, kline, sr1, yhe1, hei1, temp1d(i), vth_fiducial, xnue0, b2, b3, rho0)/sr1   !in cgs           
   write(*,*) temp1d(i), opalbar1d_a(i), opalbar1d_b(i), opalbar1d_c(i)   
enddo
!
!
!
end subroutine calc_model3d_test
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine calc_model3dg
!
!star 1 has a Be-type disc from Dylans initial conditions  
!star 2 has a beta-velocity wind
!
use prog_type
use fund_const
use options_modspec, only: input_file, input_file2, input_mod
use dime_modspec, only: cs1_nr, cs1_ntheta, cs1_nphi, cs1_r, cs1_theta, cs1_phi, &
                        cs1_velx3d, cs1_vely3d, cs1_velz3d, cs1_opac3d, &
                        cs1_opalbar3d, cs1_sline3d, cs1_scont3d, cs1_imask3d, cs1_t3d, cs1_rho3d, &
                        cs2_nr, cs2_ntheta, cs2_nphi, cs2_r, cs2_theta, cs2_phi, &
                        cs2_velx3d, cs2_vely3d, cs2_velz3d, cs2_opac3d, &
                        cs2_opalbar3d, cs2_sline3d, cs2_scont3d, cs2_imask3d, cs2_t3d, cs2_rho3d
use params_modspec, only: vx01, vy01, vz01, vx02, vy02, vz02, x01, y01, z01, x02, y02, z02, &
                          vmin1, vmax1, beta1, vrot1, rstar1, rmin1, rmax1, sr1, teff1, yhe1, hei1, logg1, lstar1, mdot1, vmicro1, &
                          vmin2, vmax2, beta2, vrot2, rstar2, rmin2, rmax2, sr2, teff2, yhe2, hei2, logg2, lstar2, mdot2, vmicro2, &
                          vmax, vth_fiducial, trad1, unit_length
use hdf5
use mod_grid, only: grid_log, grid_equi
use mod_opacities, only: opac_thomson, get_opalbar
use mod_iline, only: iline, na, xnue0, kline
!
implicit none
!
! ... local scalars
integer(i4b) :: i, j, k, err, idum1, idum2, ii
real(dp) :: vel, opac, opalbar, rho, sline, scont
real(dp) :: beta, vmin, vinf, mdot, mdot_cgs, xic1, eps_line, gradv, temp, velr, vth
real(dp) :: sint1, cost1, sinp1, cosp1, theta1, phi1, rad1, theta2, phi2, rad2, sint2, cost2, sinp2, cosp2, &
            xcoord0, ycoord0, zcoord0, xcoord1, ycoord1, zcoord1, &
            velx0, vely0, velz0, velx1, vely1, velz1, velx2, vely2, velz2, xcoord2, ycoord2, zcoord2, kcont

real(dp) :: hi1, mmw1, mstar1_cgs, tdisc1, mdisc1, dtheta_disc1, theta_disc1l, theta_disc1u, rho0_disc1, rdisc1_max, slope1
real(dp) :: hi2, mmw2, mstar2_cgs, twind2, mdisc2, dtheta_disc2, theta_disc2l, theta_disc2u, rho0_disc2, rdisc2_max, slope2
real(dp) :: csound, velphi, b2, b3, bconst2, tdisc1a, tdisc1b, tacoeff, tbcoeff
!
! ... local arrays
real(dp), dimension(:), allocatable :: fdum1_arr, fdum2_arr
!
! ... local characters
!
! ... local functions
real(dp) :: bnue, bvel, vthermal, vsound, mean_molecular_weight, sline_depcoeff, dilfac
!
! ... for temperature structure
integer(i4b) :: opt_temperature
integer(i4b) :: cs1_np
integer(i4b), parameter :: nz=301
real(dp) :: zmin, zmax, tauz, dtau, rad, tfloor
real(dp), dimension(:), allocatable :: cs1_p1d, cs1_tauz1d
real(dp), dimension(nz) :: opac1d, z1d
!
! ... for hdf5 file
!
tfloor=0.8d0*(teff1+teff2)/2.d0
!write(*,*) tfloor
!stop
opt_temperature = 4
!=0 for constant temperature at tdisc
!=1 for linear temperature law
!=2 for 1/r temperature law
!=3 for lucy optically thin temperature law
!=4 for carciofi temperature law
!
kcont=1.d0
kline=1.d0
!
!define parameters that are always required
vth_fiducial=1.d7
vmicro1=5.d6
logg1=4.15d0
lstar1=1.d6
vrot1=0.d8
yhe1=0.1d0
hei1=2.d0
!
!need to be defined (if gravity darkening and/or photospheric profile shall be accounted for in spec_vbin.eo)
vmicro2=5.d6
logg2=3.1d0
lstar2=1.d6
vrot2=0.d7
yhe2=0.2d0
hei2=2.d0
!
!define disc 1
hi1 = 1.d0   !number free electrons for each hydrogen atom
mmw1 = mean_molecular_weight(hi1,hei1,yhe1)  !mean molecular weight
mstar1_cgs = sr1**2 * ten**logg1/cgs_grav
vrot1 = sqrt(cgs_grav*mstar1_cgs/sr1)    !breakup velocity (consistent with v_phi of disc model)
tdisc1 = 10.d3            !isothermal temperature of the disc
tdisc1a=18.d3
tdisc1b=5.d3
mdisc1 = 2.d-11 * xmsu     !*1.d-10 set to almost zero
dtheta_disc1 = 45.d0*pi/180.d0    !opening angle of the disc (plus/minus 5 degree here)
theta_disc1l = pi/two - dtheta_disc1
theta_disc1u = pi/two + dtheta_disc1
!maximum radius of the disc (in rstar1)
rdisc1_max=13.5d0
tacoeff = tdisc1a-(tdisc1b-tdisc1a)*one/(rdisc1_max-one)
tbcoeff = (tdisc1b-tdisc1a)/(rdisc1_max-one)
!slope of the disc
slope1=3.5d0

!define wind for star 2
hi2 = 1.d0   !number free electrons for each hydrogen atom
vrot2 = 0.d0
twind2 = 0.8d0*teff2
mdot2 = 4.d-8   !in msun/yr
vmin2 = 10.d0   !in km/s
vmax2 = 150.d0
beta2 = 1.d0
bconst2 = (1.d0-vmin2/vmax2)**(1.d0/beta2)

!in cgs
mdot2 = mdot2*xmsu/365.25d0/24.d0/3600.d0
vmin2 = vmin2*1.d5
vmax2 = vmax2*1.d5
!
!
!
write(*,*) '----------------calc_model3dg: Be-type disc for star 1 (Dylan)----------------'
write(*,*)
!
!
cs1_nr=101
idum1=12
idum2=20
cs1_ntheta=2*idum1+2*idum2+1
cs1_nphi=2*cs1_ntheta-1
!
!allocate arrays
allocate(cs1_r(cs1_nr), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3df'
allocate(cs1_theta(cs1_ntheta), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3df'
allocate(cs1_phi(cs1_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3df'
allocate(cs1_opac3d(cs1_nr,cs1_ntheta,cs1_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3df'
allocate(cs1_opalbar3d(cs1_nr,cs1_ntheta,cs1_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3df'
allocate(cs1_scont3d(cs1_nr,cs1_ntheta,cs1_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3df'
allocate(cs1_sline3d(cs1_nr,cs1_ntheta,cs1_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3df'
allocate(cs1_velx3d(cs1_nr,cs1_ntheta,cs1_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3df'
allocate(cs1_vely3d(cs1_nr,cs1_ntheta,cs1_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3df'
allocate(cs1_velz3d(cs1_nr,cs1_ntheta,cs1_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3df'
allocate(cs1_t3d(cs1_nr,cs1_ntheta,cs1_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3df'
allocate(cs1_rho3d(cs1_nr,cs1_ntheta,cs1_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3df'
!
!radius in units of sr1   
call grid_log(rmin1, rmax1, cs1_nr,  cs1_r)
!
if(cs1_r(cs1_nr-1).le.rdisc1_max) then
   write(*,*) 'error: we require always two zero-opacity ghost points at outer boundary'
   stop
endif
!
!theta-grid with small steps at equator (always ntheta/2 for the disc)
allocate(fdum1_arr(idum1))
call grid_log(pi/2.d0+5.d-3, pi/2.d0+dtheta_disc1-1.d-3, idum1, fdum1_arr)
allocate(fdum2_arr(idum2))
call grid_log(2.d0*fdum1_arr(idum1)-fdum1_arr(idum1-1), pi, idum2, fdum2_arr)
ii=1
do i=1, idum2
   cs1_theta(ii)=pi-fdum2_arr(idum2+1-i)
   ii=ii+1
enddo
do i=1, idum1
   cs1_theta(ii)=pi-fdum1_arr(idum1+1-i)
   ii=ii+1
enddo
cs1_theta(ii)=pi/two
ii=ii+1
do i=1, idum1
   cs1_theta(ii)=fdum1_arr(i)
   ii=ii+1
enddo
do i=1, idum2
   cs1_theta(ii)=fdum2_arr(i)
   ii=ii+1
enddo
cs1_theta(1)=zero
cs1_theta(cs1_ntheta)=pi
!
!
!call grid_equi(zero, pi, cs1_ntheta, cs1_theta)
!
!phi-grid equidistant
call grid_equi(zero, two*pi, cs1_nphi, cs1_phi)
!
!-------calculate temperature structure as function of r*sin(theta)-----
!
cs1_np=cs1_nr
allocate(cs1_p1d(cs1_np))
allocate(cs1_tauz1d(cs1_np))
cs1_p1d=cs1_r
!
csound = vsound(teff1,mmw1)
rho0_disc1 = mdisc1*sqrt(cgs_grav*mstar1_cgs)*log(ten)/two/sqrt(pi**3)/csound/sr1**3.5d0
do i=1, cs1_np
   zmin = cs1_p1d(i)*cos(theta_disc1u)   
   zmax = cs1_p1d(i)*cos(theta_disc1l)
   call grid_equi(zmin, zmax, nz, z1d)
   do j=1, nz
      rad = sqrt(z1d(j)**2 + cs1_p1d(i)**2)
      rho = rho0_disc1*(cs1_p1d(i))**(-slope1) * exp(cgs_grav*mstar1_cgs/csound**2 * (one/sr1/rad-one/sr1/cs1_p1d(i)))
      opac1d(j) = opac_thomson(yhe1, hei1, rho, kcont)*sr1
   enddo
   tauz = 0.d0
   do j=2, nz
      dtau = 0.5d0*(opac1d(j-1)+opac1d(j))*(z1d(j)-z1d(j-1))
      tauz = tauz + dtau
   enddo
   cs1_tauz1d(i) = tauz

!   rad = cs1_p1d(i)
!   rho = rho0_disc1*(cs1_p1d(i))**(-slope1) * exp(cgs_grav*mstar1_cgs/csound**2 * (one/sr1/rad-one/sr1/cs1_p1d(i)))   
!   write(*,*) cs1_p1d(i), zmin, zmax, tauz, rho
enddo

!write(*,*) rho0_disc1
!stop
!
!-----------------------------------------------------------------------
!
!
do i=1, cs1_nr
   do j=1, cs1_ntheta
      do k=1, cs1_nphi
!   do j=33, 33!cs1_ntheta
!      do k=1, 1!cs1_nphi
!
!zero values if inside secondary         
         rad1=cs1_r(i)
         sint1=sin(cs1_theta(j))
         cost1=cos(cs1_theta(j))
         sinp1=sin(cs1_phi(k))
         cosp1=cos(cs1_phi(k))
!position in global coordinate system
         xcoord0 = x01 + rad1*sint1*cosp1*rstar1/unit_length
         ycoord0 = y01 + rad1*sint1*sinp1*rstar1/unit_length
         zcoord0 = z01 + rad1*cost1*rstar1/unit_length
!postion in system of star 2
         xcoord2 = (xcoord0-x02)*unit_length/rstar2
         ycoord2 = (ycoord0-y02)*unit_length/rstar2
         zcoord2 = (zcoord0-z02)*unit_length/rstar2
         rad2=sqrt(xcoord2**2+ycoord2**2+zcoord2**2)
         call get_angles_spc(xcoord2, ycoord2, zcoord2, theta2, phi2)
!
         if(rad2.lt.rmin2) then
            rho = zero
            temp = zero
            vth = zero
            opalbar = zero
            sline = zero
            scont = zero
            opac = zero
            velx1 = zero
            vely1 = zero
            velz1 = zero
         elseif(rad1.gt.rdisc1_max) then
            rho = zero
            temp = zero
            vth = zero
            opalbar = zero
            sline = zero
            scont = zero
            opac = zero
            velx1 = zero
            vely1 = zero
            velz1 = zero
            if(rad2.lt.rmax2) then   !set wind of star 2
               velr = bvel(rad2, vmax2, bconst2, beta2)
               sint2=sin(theta2)
!velocity components in system of star 2
               velx2 = velr*sin(theta2)*cos(phi2)
               vely2 = velr*sin(theta2)*sin(phi2)
               velz2 = velr*cos(theta2)
!velocity components in global system
               velx0 = velx2 + vx02
               vely0 = vely2 + vy02
               velz0 = velz2 + vz02
!velocity components in system of star 1
               velx1 = velx0 - vx01
               velx1 = vely0 - vy01
               velz1 = velz0 - vz01
               
               temp = twind2
               vth = vthermal(vmicro2, temp, na)
               rho = mdot2/4.d0/pi/(rad2*sr2)**2/velr
!departure coeffcients in LTE
               b2=one
               b3=one
               opalbar = get_opalbar(iline, kline, sr2, yhe2, hei2, temp, vth_fiducial, xnue0, b2, b3, rho)/sr2   !in cgs
               sline = sline_depcoeff(xnue0, temp, b2, b3)
               scont = bnue(xnue0,twind2)
               opac = opac_thomson(yhe2, hei2, rho, kcont)
            endif
               
         elseif(cs1_theta(j).gt.theta_disc1l.and.&
                cs1_theta(j).lt.theta_disc1u) then
!
            sint1=sin(cs1_theta(j))
            velphi = sqrt(cgs_grav*mstar1_cgs/cs1_r(i)/sr1/sint1)
            velx1 = -velphi*sinp1
            vely1 = velphi*cosp1
            velz1 = zero

            if(opt_temperature.eq.0) then
               temp=tdisc1
            elseif(opt_temperature.eq.1) then
               temp = tacoeff + tbcoeff*rad1
            elseif(opt_temperature.eq.2) then
               temp = max(teff1/rad1,5d3)
            elseif(opt_temperature.eq.3) then
               temp = teff1 * dilfac(rad1)**0.25d0
            elseif(opt_temperature.eq.4) then
               temp = teff1/pi**0.25d0 * (asin(one/rad1) - sqrt(one-one/(rad1)**2)/rad1)**0.25d0
               if(cs1_tauz1d(i).lt.0.1d0) temp = max(temp, 0.6d0*teff1)
            else
               stop 'opt_temperature not properly specified'
            endif

            vth = vthermal(vmicro1, temp, na)
            csound = vsound(tdisc1,mmw1)
            rho0_disc1 = mdisc1*sqrt(cgs_grav*mstar1_cgs)*log(ten)/two/sqrt(pi**3)/csound/sr1**3.5d0
            rho = rho0_disc1*(cs1_r(i)*sint1)**(-slope1) * exp(cgs_grav*mstar1_cgs/csound**2 * (one/sr1/cs1_r(i)-one/sr1/cs1_r(i)/sint1))
!
!departure coeffcients in LTE
            b2=one
            b3=one
            opalbar = get_opalbar(iline, kline, sr1, yhe1, hei1, temp, vth_fiducial, xnue0, b2, b3, rho)/sr1   !in cgs
            sline = sline_depcoeff(xnue0, temp, b2, b3)
            scont = bnue(xnue0,temp)
            opac = opac_thomson(yhe1, hei1, rho, kcont)
!            write(*,*) rho, opalbar
!            write(*,*) cs1_r(i), temp, b2, b3, rho, opalbar*sr1
         else
            velx1 = zero
            vely1 = zero
            velz1 = zero
            temp = 0.8d0*(teff1+teff2)/2.d0
            vth = vthermal(vmicro1, temp, na)
            rho = zero
            opalbar = zero
            sline = zero
            scont = zero
            opac = zero         
         endif   
!        
!----------------------line source function-----------------------------
!
         cs1_sline3d(i,j,k)=sline
!
!----------------------continuum source function------------------------
!
         cs1_scont3d(i,j,k)=scont
!
!------------------------velocity components----------------------------
!
         cs1_velx3d(i,j,k) = velx1
         cs1_vely3d(i,j,k) = vely1
         cs1_velz3d(i,j,k) = velz1
!
!---------------------------temperature---------------------------------
!
         cs1_t3d(i,j,k) = temp
!
!-----------------------------opacity-----------------------------------
!
         cs1_rho3d(i,j,k) = rho
         cs1_opalbar3d(i,j,k) = opalbar*sr1   !in units of sr1
         cs1_opac3d(i,j,k) = opac*sr1
!
      enddo
   enddo
enddo
!
!stop 'go on in modelspec'

!j=30
!write(*,*) cs1_theta(j)*180./pi
!k=1
!do i=1, cs1_nr
!   write(*,*) i, cs1_r(i), cs1_opalbar3d(i,j,k)
!enddo
!stop

!
!
write(*,*) '----------------calc_model3dg: calculating 3d model for star 2-----------------'
write(*,*)
!
!
cs2_nr=101
cs2_ntheta=65
cs2_nphi=2*cs1_ntheta-1
!
!allocate arrays
allocate(cs2_r(cs2_nr), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3dg'
allocate(cs2_theta(cs2_ntheta), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3dg'
allocate(cs2_phi(cs2_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3dg'
allocate(cs2_opac3d(cs2_nr,cs2_ntheta,cs2_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3dg'
allocate(cs2_opalbar3d(cs2_nr,cs2_ntheta,cs2_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3dg'
allocate(cs2_scont3d(cs2_nr,cs2_ntheta,cs2_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3dg'
allocate(cs2_sline3d(cs2_nr,cs2_ntheta,cs2_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3dg'
allocate(cs2_velx3d(cs2_nr,cs2_ntheta,cs2_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3dg'
allocate(cs2_vely3d(cs2_nr,cs2_ntheta,cs2_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3dg'
allocate(cs2_velz3d(cs2_nr,cs2_ntheta,cs2_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3dg'
allocate(cs2_t3d(cs2_nr,cs2_ntheta,cs2_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3dg'
allocate(cs2_rho3d(cs2_nr,cs2_ntheta,cs2_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3dg'
!
!radius in units of sr2
call grid_log(rmin2, rmax2, cs2_nr,  cs2_r)
!
if(cs2_r(cs2_nr-1).le.rdisc2_max) then
   write(*,*) 'error: we require always two zero-opacity ghost points at outer boundary'
   stop
endif
!
!
call grid_equi(zero, pi, cs2_ntheta, cs2_theta)
!
!phi-grid equidistant
call grid_equi(zero, two*pi, cs2_nphi, cs2_phi)
!
!
!
do i=1, cs2_nr
   do j=1, cs2_ntheta
      do k=1, cs2_nphi
!   do j=16, 16!cs2_ntheta
!      do k=1, 1!cs2_nphi
!
         rad2=cs2_r(i)
         theta2=cs2_theta(j)
         sint2=sin(cs2_theta(j))
         cost2=cos(cs2_theta(j))
         sinp2=sin(cs2_phi(k))
         cosp2=cos(cs2_phi(k))

!position in global coordinate system
         xcoord0 = x02 + rad2*sint2*cosp2*rstar2/unit_length
         ycoord0 = y02 + rad2*sint2*sinp2*rstar2/unit_length
         zcoord0 = z02 + rad2*cost2*rstar2/unit_length
!postion in system of star 1
         xcoord1 = (xcoord0-x01)*unit_length/rstar1
         ycoord1 = (ycoord0-y01)*unit_length/rstar1
         zcoord1 = (zcoord0-z01)*unit_length/rstar1
         rad1=sqrt(xcoord1**2+ycoord1**2+zcoord1**2)
         call get_angles_spc(xcoord1, ycoord1, zcoord1, theta1, phi1)         
!
!atmosphere for star 1
         if(rad1.lt.rmin1) then
            rho = zero
            velx2 = zero
            vely2 = zero
            velz2 = zero
            temp = zero
            vth = zero
            opalbar = zero
            sline = zero
            scont = zero
            opac = zero
         elseif(rad1.gt.rdisc1_max) then
            rho = zero
            velx2 = zero
            vely2 = zero
            velz2 = zero
            temp = zero
            vth = zero
            opalbar = zero
            sline = zero
            scont = zero
            opac = zero
            if(rad2.lt.rmax2-1.d-8) then   !set wind of star 2
               velr = bvel(rad2, vmax2, bconst2, beta2)
!               write(*,*) rad2, rmax2
               velx2 = velr*sint2*cosp2
               vely2 = velr*sint2*sinp2
               velz2 = velr*cost2
               temp = twind2
               rho = mdot2/4.d0/pi/(rad2*sr2)**2/velr
!departure coeffcients in LTE
               b2=one
               b3=one
               opalbar = get_opalbar(iline, kline, sr2, yhe2, hei2, temp, vth_fiducial, xnue0, b2, b3, rho)/sr2   !in cgs
               sline = sline_depcoeff(xnue0, temp, b2, b3)
               scont = bnue(xnue0,twind2)
               opac = opac_thomson(yhe2, hei2, rho, kcont)
            endif            
         elseif(theta1.gt.theta_disc1l.and.&
                theta1.lt.theta_disc1u) then
!              
               sint1=sin(theta1)
               velphi = sqrt(cgs_grav*mstar1_cgs/rad1/sr1/sint1)
!velocity components in system of star 1
               velx1 = -velphi*sin(phi1)
               vely1 = velphi*cos(phi1)
               velz1 = zero
!velocity components in global system
               velx0 = velx1 + vx01
               vely0 = vely1 + vy01
               velz0 = velz1 + vz01
!velocity components in system of star 2         
               velx2 = velx0 - vx02
               vely2 = vely0 - vy02
               velz2 = velz0 - vz02
!
               temp = tacoeff + tbcoeff*rad1 !tdisc1
               vth = vthermal(vmicro1, temp, na)
               csound = vsound(tdisc1,mmw1)
               rho0_disc1 = mdisc1*sqrt(cgs_grav*mstar1_cgs)*log(ten)/two/sqrt(pi**3)/csound/sr1**3.5d0
               rho = rho0_disc1*(rad1*sint1)**(-slope1) * exp(cgs_grav*mstar1_cgs/csound**2 * (one/sr1/rad1-one/sr1/rad1/sint1))
!
               b2=one
               b3=one
         
               opalbar = get_opalbar(iline, kline, sr1, yhe1, hei1, temp, vth_fiducial, xnue0, b2, b3, rho)/sr1   !in cgs                        
               sline = sline_depcoeff(xnue0, temp, b2, b3)
               scont = bnue(xnue0,temp)
               opac = opac_thomson(yhe1, hei1, rho, kcont)


!               write(*,*) rad2, rad1, velphi, theta1*180.d0/pi, rho_disc0, rho, opalbar
         else
            rho = zero
            velx2 = zero
            vely2 = zero
            velz2 = zero
            temp = trad1
            vth = vthermal(vmicro1, temp, na)
            rho = zero
            opalbar = zero
            sline = zero
            scont = zero
            opac = zero         
         endif

!         write(*,*) rad1, rad2, opalbar*sr2, theta1, phi1, cs2_theta(j), cs2_phi(k)
!        
!----------------------line source function-----------------------------
!
         cs2_sline3d(i,j,k)=sline
!
!----------------------continuum source function------------------------
!
         cs2_scont3d(i,j,k)=scont
!
!------------------------velocity components----------------------------
!
         cs2_velx3d(i,j,k) = velx2
         cs2_vely3d(i,j,k) = vely2
         cs2_velz3d(i,j,k) = velz2
!
!---------------------------temperature---------------------------------
!
         cs2_t3d(i,j,k) = temp
!
!-----------------------------opacity-----------------------------------
!
         cs2_rho3d(i,j,k) = rho
         cs2_opalbar3d(i,j,k) = opalbar*sr2  !in units of sr2
         cs2_opac3d(i,j,k) = opac*sr2
!         write(*,*) rad2, aopac
         
!
      enddo
   enddo
enddo
!
!stop 'go on in modelspec'
!
!-----------------------------set velocities------------------------------
!
!maximum velocity in global system
!
vmax=zero
do i=1, cs1_nr
   do j=1, cs1_ntheta
      do k=1, cs1_nphi
         vel=sqrt((vx01+cs1_velx3d(i,j,k))**2+(vy01+cs1_vely3d(i,j,k))**2+(vz01+cs1_velz3d(i,j,k))**2)
         if(vel.gt.vmax) vmax=vel
      enddo
   enddo
enddo

do i=1, cs2_nr
   do j=1, cs2_ntheta
      do k=1, cs2_nphi
         vel=sqrt((vx02+cs2_velx3d(i,j,k))**2+(vy02+cs2_vely3d(i,j,k))**2+(vz02+cs2_velz3d(i,j,k))**2)
         if(vel.gt.vmax) vmax=vel
      enddo
   enddo
enddo
!
!

!zero rotation of primary
vrot2=0.d0
!
!for now, set zero continuum
!cs1_opac3d=0.d0
!cs1_scont3d=0d0
!cs2_opac3d=0.d0
!cs2_scont3d=0.d0
!

!
end subroutine calc_model3dg
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine calc_model3de
!
!star 1 has a Be-type disc from Ileyks simulations
!
use prog_type
use fund_const
use options_modspec, only: input_file, input_file2, input_mod
use dime_modspec, only: cs1_nr, cs1_ntheta, cs1_nphi, cs1_r, cs1_theta, cs1_phi, &
                        cs1_velx3d, cs1_vely3d, cs1_velz3d, cs1_opac3d, &
                        cs1_opalbar3d, cs1_sline3d, cs1_scont3d, cs1_imask3d, cs1_t3d, cs1_rho3d, &
                        cs2_nr, cs2_ntheta, cs2_nphi, cs2_r, cs2_theta, cs2_phi, &
                        cs2_velx3d, cs2_vely3d, cs2_velz3d, cs2_opac3d, &
                        cs2_opalbar3d, cs2_sline3d, cs2_scont3d, cs2_imask3d, cs2_t3d, cs2_rho3d
use params_modspec, only: vx01, vy01, vz01, vx02, vy02, vz02, x01, y01, z01, x02, y02, z02, &
                          vmin1, vmax1, beta1, vrot1, rstar1, rmin1, rmax1, sr1, teff1, yhe1, hei1, logg1, lstar1, mdot1, vmicro1, &
                          vmin2, vmax2, beta2, vrot2, rstar2, rmin2, rmax2, sr2, teff2, yhe2, hei2, logg2, lstar2, mdot2, vmicro2, &
                          vmax, vth_fiducial, trad1, unit_length
use hdf5
use mod_interp3d, only: get_rtp_indx, coeff3d_8p_lin
use mod_grid, only: grid_log, grid_equi
use mod_opacities, only: opac_thomson, get_opalbar
use mod_iline, only: iline, na, xnue0, kline
!
implicit none
!
! ... local scalars
integer(i4b) :: i, j, k, err, idum1, idum2, iim1, ii, jjm1, jj, kkm1, kk
real(dp) :: vel, opac, opalbar, rho, sline, scont
real(dp) :: beta, vmin, vinf, mdot, mdot_cgs, bconst, xic1, eps_line, gradv, temp, vth
real(dp) :: sint1, cost1, sinp1, cosp1, theta1, phi1, rad1, rad2, sint2, cost2, sinp2, cosp2, &
            xcoord0, ycoord0, zcoord0, xcoord1, ycoord1, zcoord1, &
            velx0, vely0, velz0, velx1, vely1, velz1, xcoord2, ycoord2, zcoord2, kcont, &
            acoeff, bcoeff, ccoeff, dcoeff, ecoeff, fcoeff, gcoeff, hcoeff
real(dp), parameter :: rho_scalefac=.1d0

real(dp) :: hi1, mmw, mstar1_cgs, csound, velphi, velr, velth, b2, b3, tdisc1, tdisc2, tacoeff, tbcoeff, rdisc_max
!
! ... local arrays
real(dp), dimension(:), allocatable :: fdum1_arr, fdum2_arr
real(dp), dimension(:,:,:), allocatable :: cs1_velr3d, cs1_velth3d, cs1_velphi3d
!
! ... local characters
character(len=19), parameter :: cs1_fname='models/ileyk/f10.h5'
!
! ... local logicals
logical :: expol
!
! ... local functions
real(dp) :: bnue, bvel, vthermal, vsound, mean_molecular_weight, sline_depcoeff
!
! ... for hdf5 file
integer(hid_t) :: file_id, dset_id, attr_id, group_id
integer(hsize_t), dimension(1) :: dims_scalars = (/1/)
integer(hsize_t), dimension(1) :: dims_r, dims_theta, dims_phi
integer(hsize_t), dimension(3) :: dims_3d
!
kcont=1.d0
kline=1.d0


tdisc1 = 5.d3!18.d3
tdisc2 = 5.d3
rdisc_max=10.d0
tacoeff = tdisc1-(tdisc2-tdisc1)*one/(rdisc_max-one)
tbcoeff = (tdisc2-tdisc1)/(rdisc_max-one)
!
write(*,*) '----------------calc_model3de: Be-type disc from Ileyk------------------------'
write(*,*)
!
!define parameters that are always required
vth_fiducial=1.d7
vmicro1=1.5d6
logg1=4.d0
lstar1=1.d6
vrot1=0.d8
yhe1=0.1d0
hei1=2.d0
!
!define local parameters (and maybe overwrite global ones)
hi1 = 1.d0   !number free electrons for each hydrogen atom
mmw = mean_molecular_weight(hi1,hei1,yhe1)  !mean molecular weight
mstar1_cgs = sr1**2 * ten**logg1/cgs_grav
vrot1 = sqrt(cgs_grav*mstar1_cgs/sr1)
!
!
call h5open_f (err)
   call h5fopen_f(cs1_fname, h5f_acc_rdonly_f, file_id, err)
!
      call h5gopen_f(file_id, 'dimensions', group_id, err)
         call h5aopen_f(group_id, 'nr', attr_id, err)
            call h5aread_f(attr_id, h5t_native_integer, cs1_nr, dims_scalars, err)
         call h5aclose_f(attr_id, err)
         call h5aopen_f(group_id, 'ntheta', attr_id, err)
            call h5aread_f(attr_id, h5t_native_integer, cs1_ntheta, dims_scalars, err)
         call h5aclose_f(attr_id, err)
         call h5aopen_f(group_id, 'nphi', attr_id, err)
            call h5aread_f(attr_id, h5t_native_integer, cs1_nphi, dims_scalars, err)
         call h5aclose_f(attr_id, err)
      call h5gclose_f(group_id, err)

      dims_r = (/ cs1_nr /)
      dims_theta = (/ cs1_ntheta /)
      dims_phi = (/ cs1_nphi /)      
      dims_3d = (/ cs1_nr, cs1_ntheta, cs1_nphi /)

!allocate arrays
      allocate(cs1_r(cs1_nr), stat=err)
      allocate(cs1_theta(cs1_ntheta), stat=err)
      allocate(cs1_phi(cs1_nphi), stat=err)
      allocate(cs1_rho3d(cs1_nr,cs1_ntheta,cs1_nphi), stat=err)
      allocate(cs1_velr3d(cs1_nr,cs1_ntheta,cs1_nphi), stat=err)
      allocate(cs1_velth3d(cs1_nr,cs1_ntheta,cs1_nphi), stat=err)
      allocate(cs1_velphi3d(cs1_nr,cs1_ntheta,cs1_nphi), stat=err)      
      allocate(cs1_t3d(cs1_nr,cs1_ntheta,cs1_nphi), stat=err)

!
      call h5gopen_f(file_id, 'coordinates', group_id, err)
         call h5dopen_f(group_id, 'r', dset_id, err)
            call h5dread_f(dset_id, h5t_native_double, cs1_r, dims_r, err)
         call h5dclose_f(dset_id, err)
         call h5dopen_f(group_id, 'theta', dset_id, err)
            call h5dread_f(dset_id, h5t_native_double, cs1_theta, dims_theta, err)
         call h5dclose_f(dset_id, err)
         call h5dopen_f(group_id, 'phi', dset_id, err)
            call h5dread_f(dset_id, h5t_native_double, cs1_phi, dims_phi, err)
         call h5dclose_f(dset_id, err)
      call h5gclose_f(group_id, err)
!
      call h5gopen_f(file_id, 'model3d', group_id, err)
         call h5dopen_f(group_id, 't3d', dset_id, err)
            call h5dread_f(dset_id, h5t_native_double, cs1_t3d, dims_3d, err)
         call h5dclose_f(dset_id, err)
         call h5dopen_f(group_id, 'rho3d', dset_id, err)
            call h5dread_f(dset_id, h5t_native_double, cs1_rho3d, dims_3d, err)
         call h5dclose_f(dset_id, err)
         cs1_rho3d=cs1_rho3d*rho_scalefac   
         call h5dopen_f(group_id, 'velr3d', dset_id, err)
            call h5dread_f(dset_id, h5t_native_double, cs1_velr3d, dims_3d, err)
         call h5dclose_f(dset_id, err)
         call h5dopen_f(group_id, 'velth3d', dset_id, err)
            call h5dread_f(dset_id, h5t_native_double, cs1_velth3d, dims_3d, err)
         call h5dclose_f(dset_id, err)
         call h5dopen_f(group_id, 'velphi3d', dset_id, err)
            call h5dread_f(dset_id, h5t_native_double, cs1_velphi3d, dims_3d, err)
         call h5dclose_f(dset_id, err)
      call h5gclose_f(group_id, err)
!
   call h5fclose_f(file_id, err)
call h5close_f(err)
!
!
!scale to rstar1
cs1_r=cs1_r/cs1_r(1)*sr1

!radius in units of sr1
cs1_r=cs1_r/sr1
!
if(abs(rmin1-minval(cs1_r)).gt.small_number) then
   write(*,*) 'error in calc_model3de: rmin1 not matching'
   write(*,*) 'set rmin1 to', minval(cs1_r)
   stop
endif
if(abs(rmax1-maxval(cs1_r)).gt.small_number) then
   write(*,*) 'error in calc_model3de: rmax1 not matching'
   write(*,*) 'set rmax1 to', maxval(cs1_r)
   stop
endif
!
!
!
allocate(cs1_opac3d(cs1_nr,cs1_ntheta,cs1_nphi), stat=err)
allocate(cs1_opalbar3d(cs1_nr,cs1_ntheta,cs1_nphi), stat=err)
allocate(cs1_scont3d(cs1_nr,cs1_ntheta,cs1_nphi), stat=err)
allocate(cs1_sline3d(cs1_nr,cs1_ntheta,cs1_nphi), stat=err)
allocate(cs1_velx3d(cs1_nr,cs1_ntheta,cs1_nphi), stat=err)
allocate(cs1_vely3d(cs1_nr,cs1_ntheta,cs1_nphi), stat=err)
allocate(cs1_velz3d(cs1_nr,cs1_ntheta,cs1_nphi), stat=err)
!
!
!
do i=1, cs1_nr
   do j=1, cs1_ntheta
      do k=1, cs1_nphi
!
!manipulate temperature
         cs1_t3d(i,j,k) = tacoeff + tbcoeff*cs1_r(i)

!
!zero values if inside secondary         
         rad1=cs1_r(i)
         sint1=sin(cs1_theta(j))
         cost1=cos(cs1_theta(j))
         sinp1=sin(cs1_phi(k))
         cosp1=cos(cs1_phi(k))
!position in global coordinate system
         xcoord0 = x01 + rad1*sint1*cosp1*rstar1/unit_length
         ycoord0 = y01 + rad1*sint1*sinp1*rstar1/unit_length
         zcoord0 = z01 + rad1*cost1*rstar1/unit_length
!postion in system of star 2
         xcoord2 = (xcoord0-x02)*unit_length/rstar2
         ycoord2 = (ycoord0-y02)*unit_length/rstar2
         zcoord2 = (zcoord0-z02)*unit_length/rstar2
         rad2=sqrt(xcoord2**2+ycoord2**2+zcoord2**2)      
!
         if(rad2.lt.rmin2) then
            velr = zero
            velth = zero
            velphi = zero
            temp = zero
            vth = zero
            opalbar = zero
            sline = zero
            scont = zero
            opac = zero
            rho = zero
         else
            velr = cs1_velr3d(i,j,k)
            velth = cs1_velth3d(i,j,k)
            velphi = cs1_velphi3d(i,j,k)
            temp = cs1_t3d(i,j,k)
            rho = cs1_rho3d(i,j,k)
!departure coeffcients in LTE
            b2=one
            b3=one
            opalbar = get_opalbar(iline, kline, sr1, yhe1, hei1, temp, vth_fiducial, xnue0, b2, b3, rho)/sr1   !in cgs                     
            sline = sline_depcoeff(xnue0, temp, b2, b3)
            scont = bnue(xnue0,temp)
            opac = opac_thomson(yhe1, hei1, rho, kcont)
         endif   
!        
!----------------------line source function-----------------------------
!
         cs1_sline3d(i,j,k)=sline
!
!----------------------continuum source function------------------------
!
         cs1_scont3d(i,j,k)=scont
!
!------------------------velocity components----------------------------
!
         sint1=sin(cs1_theta(j))
         cost1=cos(cs1_theta(j))
         sinp1=sin(cs1_phi(k))
         cosp1=cos(cs1_phi(k))
         
         cs1_velx3d(i,j,k) = velr*sint1*cosp1 + velth*cost1*cosp1 - velphi*sinp1
         cs1_vely3d(i,j,k) = velr*sint1*sinp1 + velth*cost1*sinp1 + velphi*cosp1
         cs1_velz3d(i,j,k) = velr*cost1 - velth*sint1
!
!---------------------------temperature---------------------------------
!
         cs1_t3d(i,j,k) = temp
!
!-----------------------------opacity-----------------------------------
!
         cs1_rho3d(i,j,k) = rho
         cs1_opalbar3d(i,j,k) = opalbar*sr1   !in units of sr1
         cs1_opac3d(i,j,k) = opac*sr1
!
      enddo
   enddo
enddo
!
!
!
write(*,*) '----------------calc_model3dd: calculating 3d model for star 2-----------------'
write(*,*)
!
!need to be defined (if gravity darkening and/or photospheric profile shall be accounted for in spec_vbin.eo)
vmicro2=1.5d6
logg2=3.1d0
lstar2=1.d6
vrot2=0.d7
yhe2=0.2d0
!
cs2_nr=31
cs2_ntheta=31
cs2_nphi=2*cs2_ntheta-1
!
!allocate arrays
allocate(cs2_r(cs2_nr), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3de'
allocate(cs2_theta(cs2_ntheta), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3de'
allocate(cs2_phi(cs2_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3de'
allocate(cs2_opac3d(cs2_nr,cs2_ntheta,cs2_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3de'
allocate(cs2_opalbar3d(cs2_nr,cs2_ntheta,cs2_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3de'
allocate(cs2_scont3d(cs2_nr,cs2_ntheta,cs2_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3de'
allocate(cs2_sline3d(cs2_nr,cs2_ntheta,cs2_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3de'
allocate(cs2_velx3d(cs2_nr,cs2_ntheta,cs2_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3de'
allocate(cs2_vely3d(cs2_nr,cs2_ntheta,cs2_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3de'
allocate(cs2_velz3d(cs2_nr,cs2_ntheta,cs2_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3de'
allocate(cs2_t3d(cs2_nr,cs2_ntheta,cs2_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3de'
allocate(cs2_rho3d(cs2_nr,cs2_ntheta,cs2_nphi), stat=err)
   if(err.ne.0) stop 'allocation error calc_model3de'
!
!radius in units of sr1   
call grid_log(rmin2, rmax2, cs2_nr,  cs2_r)
!
!theta-grid equidistant
call grid_equi(zero, pi, cs2_ntheta, cs2_theta)
!
!phi-grid equidistant
call grid_equi(zero, two*pi, cs2_nphi, cs2_phi)
!
!
!
do i=1, cs2_nr
   do j=1, cs2_ntheta
      do k=1, cs2_nphi
!
         rad2=cs2_r(i)         
         sint2=sin(cs2_theta(j))
         cost2=cos(cs2_theta(j))
         sinp2=sin(cs2_phi(k))
         cosp2=cos(cs2_phi(k))

!position in global coordinate system
         xcoord0 = x02 + rad2*sint2*cosp2*rstar2/unit_length
         ycoord0 = y02 + rad2*sint2*sinp2*rstar2/unit_length
         zcoord0 = z02 + rad2*cost2*rstar2/unit_length
!postion in system of star 1
         xcoord1 = (xcoord0-x01)*unit_length/rstar1
         ycoord1 = (ycoord0-y01)*unit_length/rstar1
         zcoord1 = (zcoord0-z01)*unit_length/rstar1
         rad1=sqrt(xcoord1**2+ycoord1**2+zcoord1**2)
         call get_angles_spc(xcoord1, ycoord1, zcoord1, theta1, phi1)         
!
!atmosphere for star 1
         if(rad1.lt.rmin1.or.rad1.gt.rmax1) then
            velr = zero
            velth = zero
            velphi = zero
            temp = zero
            opalbar = zero
            sline = zero
            scont = zero
            opac = zero         
         else
!interpolate from system 1
            call get_rtp_indx(rad1, theta1, phi1, cs1_r, cs1_theta, cs1_phi, cs1_nr, cs1_ntheta, cs1_nphi, &
                              iim1, ii, jjm1, jj, kkm1, kk, expol, rmin1, rmax1)
!
            call coeff3d_8p_lin(cs1_r(iim1), cs1_r(ii), cs1_theta(jjm1), cs1_theta(jj), &
                                cs1_phi(kkm1), cs1_phi(kk), rad1, theta1, phi1, &
                                acoeff, bcoeff, ccoeff, dcoeff, ecoeff, fcoeff, gcoeff, hcoeff)

            velr = acoeff*cs1_velr3d(iim1,jjm1,kkm1) + bcoeff*cs1_velr3d(ii,jjm1,kkm1) + &
                   ccoeff*cs1_velr3d(iim1,jj,kkm1) + dcoeff*cs1_velr3d(ii,jj,kkm1) + &
                   ecoeff*cs1_velr3d(iim1,jjm1,kk) + fcoeff*cs1_velr3d(ii,jjm1,kk) + &
                   gcoeff*cs1_velr3d(iim1,jj,kk) + hcoeff*cs1_velr3d(ii,jj,kk)

            velth = acoeff*cs1_velth3d(iim1,jjm1,kkm1) + bcoeff*cs1_velth3d(ii,jjm1,kkm1) + &
                    ccoeff*cs1_velth3d(iim1,jj,kkm1) + dcoeff*cs1_velth3d(ii,jj,kkm1) + &
                    ecoeff*cs1_velth3d(iim1,jjm1,kk) + fcoeff*cs1_velth3d(ii,jjm1,kk) + &
                    gcoeff*cs1_velth3d(iim1,jj,kk) + hcoeff*cs1_velth3d(ii,jj,kk)
         
            velphi = acoeff*cs1_velphi3d(iim1,jjm1,kkm1) + bcoeff*cs1_velphi3d(ii,jjm1,kkm1) + &
                     ccoeff*cs1_velphi3d(iim1,jj,kkm1) + dcoeff*cs1_velphi3d(ii,jj,kkm1) + &
                     ecoeff*cs1_velphi3d(iim1,jjm1,kk) + fcoeff*cs1_velphi3d(ii,jjm1,kk) + &
                     gcoeff*cs1_velphi3d(iim1,jj,kk) + hcoeff*cs1_velphi3d(ii,jj,kk)

            temp = acoeff*cs1_t3d(iim1,jjm1,kkm1) + bcoeff*cs1_t3d(ii,jjm1,kkm1) + &
                   ccoeff*cs1_t3d(iim1,jj,kkm1) + dcoeff*cs1_t3d(ii,jj,kkm1) + &
                   ecoeff*cs1_t3d(iim1,jjm1,kk) + fcoeff*cs1_t3d(ii,jjm1,kk) + &
                   gcoeff*cs1_t3d(iim1,jj,kk) + hcoeff*cs1_t3d(ii,jj,kk)
            
            sline = acoeff*cs1_sline3d(iim1,jjm1,kkm1) + bcoeff*cs1_sline3d(ii,jjm1,kkm1) + &
                    ccoeff*cs1_sline3d(iim1,jj,kkm1)   + dcoeff*cs1_sline3d(ii,jj,kkm1) + &
                    ecoeff*cs1_sline3d(iim1,jjm1,kk) + fcoeff*cs1_sline3d(ii,jjm1,kk) + &
                    gcoeff*cs1_sline3d(iim1,jj,kk) + hcoeff*cs1_sline3d(ii,jj,kk)

            scont = acoeff*cs1_scont3d(iim1,jjm1,kkm1) + bcoeff*cs1_scont3d(ii,jjm1,kkm1) + &
                    ccoeff*cs1_scont3d(iim1,jj,kkm1)   + dcoeff*cs1_scont3d(ii,jj,kkm1) + &
                    ecoeff*cs1_scont3d(iim1,jjm1,kk) + fcoeff*cs1_scont3d(ii,jjm1,kk) + &
                    gcoeff*cs1_scont3d(iim1,jj,kk) + hcoeff*cs1_scont3d(ii,jj,kk)

            opac = acoeff*cs1_opac3d(iim1,jjm1,kkm1) + bcoeff*cs1_opac3d(ii,jjm1,kkm1) + &
                   ccoeff*cs1_opac3d(iim1,jj,kkm1)   + dcoeff*cs1_opac3d(ii,jj,kkm1) + &
                   ecoeff*cs1_opac3d(iim1,jjm1,kk) + fcoeff*cs1_opac3d(ii,jjm1,kk) + &
                   gcoeff*cs1_opac3d(iim1,jj,kk) + hcoeff*cs1_opac3d(ii,jj,kk)

            opalbar = acoeff*cs1_opalbar3d(iim1,jjm1,kkm1) + bcoeff*cs1_opalbar3d(ii,jjm1,kkm1) + &
                      ccoeff*cs1_opalbar3d(iim1,jj,kkm1)   + dcoeff*cs1_opalbar3d(ii,jj,kkm1) + &
                      ecoeff*cs1_opalbar3d(iim1,jjm1,kk) + fcoeff*cs1_opalbar3d(ii,jjm1,kk) + &
                      gcoeff*cs1_opalbar3d(iim1,jj,kk) + hcoeff*cs1_opalbar3d(ii,jj,kk)

            rho = acoeff*cs1_rho3d(iim1,jjm1,kkm1) + bcoeff*cs1_rho3d(ii,jjm1,kkm1) + &
                      ccoeff*cs1_rho3d(iim1,jj,kkm1)   + dcoeff*cs1_rho3d(ii,jj,kkm1) + &
                      ecoeff*cs1_rho3d(iim1,jjm1,kk) + fcoeff*cs1_rho3d(ii,jjm1,kk) + &
                      gcoeff*cs1_rho3d(iim1,jj,kk) + hcoeff*cs1_rho3d(ii,jj,kk)

         endif
!        
!----------------------line source function-----------------------------
!
         cs2_sline3d(i,j,k)=sline
!
!----------------------continuum source function------------------------
!
         cs2_scont3d(i,j,k)=scont
!
!------------------------velocity components----------------------------
!
!velocity components in system of star 1
         sint1=sin(theta1)
         cost1=cos(theta1)
         sinp1=sin(phi1)
         cosp1=cos(phi1)
         velx1 = velr*sint1*cosp1 + velth*cost1*cosp1 - velphi*sinp1
         vely1 = velr*sint1*sinp1 + velth*cost1*sinp1 + velphi*cosp1
         velz1 = velr*cost1 - velth*sint1
!velocity components in global system
         velx0 = velx1 + vx01
         vely0 = vely1 + vy01
         velz0 = velz1 + vz01
!velocity components in system of star 2         
         cs2_velx3d(i,j,k)=velx0 - vx02
         cs2_vely3d(i,j,k)=vely0 - vy02
         cs2_velz3d(i,j,k)=velz0 - vz02

!         if(j.eq.16.and.k.eq.1) then
!            write(*,*) cs2_r(i), cs2_vely3d(i,j,k), vely0, vy02, vely0-vy02
!         endif
!
!---------------------------temperature---------------------------------
!
         cs2_t3d(i,j,k) = temp
!
!-----------------------------opacity-----------------------------------
!
         cs2_rho3d(i,j,k) = rho
         cs2_opalbar3d(i,j,k) = opalbar*sr2/sr1  !in units of sr2
         cs2_opac3d(i,j,k) = opac*sr2/sr1
!
      enddo
   enddo
enddo
!
!stop
!
!-----------------------------set velocities------------------------------
!
!maximum velocity in global system
!
vmax=zero
do i=1, cs1_nr
   do j=1, cs1_ntheta
      do k=1, cs1_nphi
         vel=sqrt((vx01+cs1_velx3d(i,j,k))**2+(vy01+cs1_vely3d(i,j,k))**2+(vz01+cs1_velz3d(i,j,k))**2)
         if(vel.gt.vmax) vmax=vel
      enddo
   enddo
enddo

do i=1, cs2_nr
   do j=1, cs2_ntheta
      do k=1, cs2_nphi
         vel=sqrt((vx02+cs2_velx3d(i,j,k))**2+(vy02+cs2_vely3d(i,j,k))**2+(vz02+cs2_velz3d(i,j,k))**2)
         if(vel.gt.vmax) vmax=vel
      enddo
   enddo
enddo
!
!
!test a rotational velocity
!vrot1=300.d5
!vrot2=300.d5
!
!for now, set zero continuum
!cs1_opac3d=0.d0
!cs1_scont3d=0d0
!cs2_opac3d=0.d0
!cs2_scont3d=0.d0
!
!
end subroutine calc_model3de
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine output3d_spc
!
use prog_type
use fund_const
use dime_modspec, only: cs1_nr, cs1_ntheta, cs1_nphi, cs1_r, cs1_theta, cs1_phi, cs1_t3d, cs1_opac3d, cs1_opalbar3d, &
                        cs1_velx3d, cs1_vely3d, cs1_velz3d, cs1_sline3d, cs1_scont3d, cs1_rho3d, &
                        cs2_nr, cs2_ntheta, cs2_nphi, cs2_r, cs2_theta, cs2_phi, cs2_t3d, cs2_opac3d, cs2_opalbar3d, &
                        cs2_velx3d, cs2_vely3d, cs2_velz3d, cs2_sline3d, cs2_scont3d, cs2_rho3d
use params_modspec
use options_modspec, only: output_file
use hdf5
use mod_iline, only: iline, na, xnue0
!
!
implicit none
!
! ... arguments
!
! ... local scalars
integer(i4b) :: i, j, k
real(dp), parameter :: rho_scalefac=0.1d0
!
! ... local arrays
!
! ... for output to hdf5
integer(i4b) :: err
integer(hid_t) :: file_id, dspace_id, dset_id, aspace_id, attr_id, group_id, group_id0
integer(hsize_t), dimension(1) :: dims_scalars = (/1/), dims_vectors = (/ 3 /)
integer(hsize_t), dimension(1) :: dims_r, dims_theta, dims_phi
integer(hsize_t), dimension(3) :: dims_3d
!
!-----------------------------------------------------------------------
!
!output stored here required for spec_vbin.eo postprocessing
!x01, y01, z01:    coordinates of star 1 in global coordinate system in corresponding units
!vx01, vy01, vz01: velocity components of star 1 coordinate system in global coordinate system (cm/s)
!rstar1:   stellar radius of star 1 in rsun
!rmin1:    minimum radius of star 1 in rstar1
!rmax1:    maximum radius of star 1 in rstar1
!teff1:    effective temperature of star 1 [K]
!trad1:    radiation temperature of star 1 [K]
!vrot1:    rotational velocity of star 1 [cm/s]
!logg1:    log(g) of star 1
!lstar1:   luminosity of star 1 [lsun]
!yhe1:     helium abundance of star 1


!
write(*,*) '-------------------------output to directory-----------------------------------'
write(*,*) 'output to file ', trim(output_file)
write(*,*)
!
!-----------------------------------------------------------------------
! 
call h5open_f(err)
   call h5fcreate_f(trim(output_file), h5f_acc_trunc_f, file_id, err)

!----------------------------input parameters---------------------------
!
      call h5gcreate_f(file_id, 'input_parameters', group_id, err)
         call h5screate_simple_f(1, dims_scalars, aspace_id, err)
!x01, y01, z01, x02, y02, z02 in unit_length
            call h5acreate_f(group_id, 'x01', h5t_native_double, aspace_id, attr_id, err)
               call h5awrite_f(attr_id, h5t_native_double, x01, dims_scalars, err)
            call h5aclose_f(attr_id, err)
            call h5acreate_f(group_id, 'y01', h5t_native_double, aspace_id, attr_id, err)
               call h5awrite_f(attr_id, h5t_native_double, y01, dims_scalars, err)
            call h5aclose_f(attr_id, err)               
            call h5acreate_f(group_id, 'z01', h5t_native_double, aspace_id, attr_id, err)
               call h5awrite_f(attr_id, h5t_native_double, z01, dims_scalars, err)
            call h5aclose_f(attr_id, err)
            call h5acreate_f(group_id, 'x02', h5t_native_double, aspace_id, attr_id, err)
               call h5awrite_f(attr_id, h5t_native_double, x02, dims_scalars, err)
            call h5aclose_f(attr_id, err)
            call h5acreate_f(group_id, 'y02', h5t_native_double, aspace_id, attr_id, err)
               call h5awrite_f(attr_id, h5t_native_double, y02, dims_scalars, err)
            call h5aclose_f(attr_id, err)               
            call h5acreate_f(group_id, 'z02', h5t_native_double, aspace_id, attr_id, err)
               call h5awrite_f(attr_id, h5t_native_double, z02, dims_scalars, err)
            call h5aclose_f(attr_id, err)

!vx01, vy01, vz01, vx02, vy02, vz02 (cm/s)
            call h5acreate_f(group_id, 'vx01', h5t_native_double, aspace_id, attr_id, err)
               call h5awrite_f(attr_id, h5t_native_double, vx01, dims_scalars, err)
            call h5aclose_f(attr_id, err)
            call h5acreate_f(group_id, 'vy01', h5t_native_double, aspace_id, attr_id, err)
               call h5awrite_f(attr_id, h5t_native_double, vy01, dims_scalars, err)
            call h5aclose_f(attr_id, err)               
            call h5acreate_f(group_id, 'vz01', h5t_native_double, aspace_id, attr_id, err)
               call h5awrite_f(attr_id, h5t_native_double, vz01, dims_scalars, err)
               call h5aclose_f(attr_id, err)
            call h5acreate_f(group_id, 'vx02', h5t_native_double, aspace_id, attr_id, err)
               call h5awrite_f(attr_id, h5t_native_double, vx02, dims_scalars, err)
            call h5aclose_f(attr_id, err)
            call h5acreate_f(group_id, 'vy02', h5t_native_double, aspace_id, attr_id, err)
               call h5awrite_f(attr_id, h5t_native_double, vy02, dims_scalars, err)
            call h5aclose_f(attr_id, err)               
            call h5acreate_f(group_id, 'vz02', h5t_native_double, aspace_id, attr_id, err)
               call h5awrite_f(attr_id, h5t_native_double, vz02, dims_scalars, err)
            call h5aclose_f(attr_id, err)               

!rstar1, rstar2 in r_sun
!rmin1, rmax1 in rstar1
!rmin2, rmax2 in rstar2            
            call h5acreate_f(group_id, 'rstar1', h5t_native_double, aspace_id, attr_id, err)
               call h5awrite_f(attr_id, h5t_native_double, rstar1, dims_scalars, err)
            call h5aclose_f(attr_id, err)
            call h5acreate_f(group_id, 'rstar2', h5t_native_double, aspace_id, attr_id, err)
               call h5awrite_f(attr_id, h5t_native_double, rstar2, dims_scalars, err)
            call h5aclose_f(attr_id, err)               
            call h5acreate_f(group_id, 'rmin1', h5t_native_double, aspace_id, attr_id, err)
               call h5awrite_f(attr_id, h5t_native_double, rmin1, dims_scalars, err)
            call h5aclose_f(attr_id, err)
            call h5acreate_f(group_id, 'rmin2', h5t_native_double, aspace_id, attr_id, err)
               call h5awrite_f(attr_id, h5t_native_double, rmin2, dims_scalars, err)
            call h5aclose_f(attr_id, err)               
            call h5acreate_f(group_id, 'rmax1', h5t_native_double, aspace_id, attr_id, err)
               call h5awrite_f(attr_id, h5t_native_double, rmax1, dims_scalars, err)
            call h5aclose_f(attr_id, err)
            call h5acreate_f(group_id, 'rmax2', h5t_native_double, aspace_id, attr_id, err)
               call h5awrite_f(attr_id, h5t_native_double, rmax2, dims_scalars, err)
            call h5aclose_f(attr_id, err)            

!temperatures in kelvin
            call h5acreate_f(group_id, 'trad1', h5t_native_double, aspace_id, &
                             attr_id, err)
               call h5awrite_f(attr_id, h5t_native_double, trad1, dims_scalars, err)
            call h5aclose_f(attr_id, err)
            call h5acreate_f(group_id, 'trad2', h5t_native_double, aspace_id, &
                             attr_id, err)
               call h5awrite_f(attr_id, h5t_native_double, trad2, dims_scalars, err)
            call h5aclose_f(attr_id, err)
            call h5acreate_f(group_id, 'teff1', h5t_native_double, aspace_id, &
                             attr_id, err)
               call h5awrite_f(attr_id, h5t_native_double, teff1, dims_scalars, err)
            call h5aclose_f(attr_id, err)
            call h5acreate_f(group_id, 'teff2', h5t_native_double, aspace_id, &
                             attr_id, err)
               call h5awrite_f(attr_id, h5t_native_double, teff2, dims_scalars, err)
            call h5aclose_f(attr_id, err)
            call h5acreate_f(group_id, 'logg1', h5t_native_double, aspace_id, &
                             attr_id, err)
               call h5awrite_f(attr_id, h5t_native_double, logg1, dims_scalars, err)
            call h5aclose_f(attr_id, err)
            call h5acreate_f(group_id, 'logg2', h5t_native_double, aspace_id, &
                             attr_id, err)
               call h5awrite_f(attr_id, h5t_native_double, logg2, dims_scalars, err)
            call h5aclose_f(attr_id, err)
            call h5acreate_f(group_id, 'lstar1', h5t_native_double, aspace_id, &
                             attr_id, err)
               call h5awrite_f(attr_id, h5t_native_double, lstar1, dims_scalars, err)
            call h5aclose_f(attr_id, err)
            call h5acreate_f(group_id, 'lstar2', h5t_native_double, aspace_id, &
                             attr_id, err)
               call h5awrite_f(attr_id, h5t_native_double, lstar2, dims_scalars, err)
            call h5aclose_f(attr_id, err)
            call h5acreate_f(group_id, 'yhe1', h5t_native_double, aspace_id, &
                             attr_id, err)
               call h5awrite_f(attr_id, h5t_native_double, yhe1, dims_scalars, err)
               call h5aclose_f(attr_id, err)
            call h5acreate_f(group_id, 'yhe2', h5t_native_double, aspace_id, &
                             attr_id, err)
               call h5awrite_f(attr_id, h5t_native_double, yhe2, dims_scalars, err)
               call h5aclose_f(attr_id, err)
            call h5acreate_f(group_id, 'fehe1', h5t_native_double, aspace_id, &
                             attr_id, err)
               call h5awrite_f(attr_id, h5t_native_double, fehe1, dims_scalars, err)
            call h5aclose_f(attr_id, err)
            call h5acreate_f(group_id, 'fehe2', h5t_native_double, aspace_id, &
                             attr_id, err)
               call h5awrite_f(attr_id, h5t_native_double, fehe2, dims_scalars, err)
            call h5aclose_f(attr_id, err)
            call h5acreate_f(group_id, 'aenh1', h5t_native_double, aspace_id, &
                             attr_id, err)
               call h5awrite_f(attr_id, h5t_native_double, aenh1, dims_scalars, err)
            call h5aclose_f(attr_id, err)
            call h5acreate_f(group_id, 'aenh2', h5t_native_double, aspace_id, &
                             attr_id, err)
               call h5awrite_f(attr_id, h5t_native_double, aenh2, dims_scalars, err)
            call h5aclose_f(attr_id, err)                                             

!all velocities in cm/s
            call h5acreate_f(group_id, 'vrot1', h5t_native_double, aspace_id, &
                             attr_id, err)
               call h5awrite_f(attr_id, h5t_native_double, vrot1, dims_scalars, err)
            call h5aclose_f(attr_id, err)
            call h5acreate_f(group_id, 'vrot2', h5t_native_double, aspace_id, &
                             attr_id, err)
               call h5awrite_f(attr_id, h5t_native_double, vrot2, dims_scalars, err)
            call h5aclose_f(attr_id, err)
            call h5acreate_f(group_id, 'vth_fiducial', h5t_native_double, aspace_id, &
                             attr_id, err)
               call h5awrite_f(attr_id, h5t_native_double, vth_fiducial, dims_scalars, err)
            call h5aclose_f(attr_id, err)
            call h5acreate_f(group_id, 'vmicro1', h5t_native_double, aspace_id, &
                             attr_id, err)
               call h5awrite_f(attr_id, h5t_native_double, vmicro1, dims_scalars, err)
            call h5aclose_f(attr_id, err)
            call h5acreate_f(group_id, 'vmicro2', h5t_native_double, aspace_id, &
                             attr_id, err)
               call h5awrite_f(attr_id, h5t_native_double, vmicro2, dims_scalars, err)
            call h5aclose_f(attr_id, err)               
            call h5acreate_f(group_id, 'vmax', h5t_native_double, aspace_id, &
                             attr_id, err)
               call h5awrite_f(attr_id, h5t_native_double, vmax, dims_scalars, err)
            call h5aclose_f(attr_id, err)

!unit length of global coordinate system [rsun]
            call h5acreate_f(group_id, 'unit_length', h5t_native_double, aspace_id, &
                             attr_id, err)
               call h5awrite_f(attr_id, h5t_native_double, unit_length, dims_scalars, err)
            call h5aclose_f(attr_id, err)
               
!transition frequency
            call h5acreate_f(group_id, 'xnue0', h5t_native_double, aspace_id, &
                             attr_id, err)
               call h5awrite_f(attr_id, h5t_native_double, xnue0, dims_scalars, err)
            call h5aclose_f(attr_id, err)
               
            call h5acreate_f(group_id, 'na', h5t_native_integer, aspace_id, &
                             attr_id, err)
               call h5awrite_f(attr_id, h5t_native_integer, na, dims_scalars, err)
            call h5aclose_f(attr_id, err)
            call h5acreate_f(group_id, 'iline', h5t_native_integer, aspace_id, &
                             attr_id, err)
               call h5awrite_f(attr_id, h5t_native_integer, iline, dims_scalars, err)
            call h5aclose_f(attr_id, err)               
         call h5sclose_f(aspace_id, err)

         call h5screate_simple_f(1, dims_vectors, dspace_id, err)
            call h5dcreate_f(group_id, 'ex01', h5t_native_double, dspace_id, dset_id, err)
               call h5dwrite_f(dset_id, h5t_native_double, ex01, dims_vectors, err)
            call h5dclose_f(dset_id, err)
            call h5dcreate_f(group_id, 'ey01', h5t_native_double, dspace_id, dset_id, err)
               call h5dwrite_f(dset_id, h5t_native_double, ey01, dims_vectors, err)
            call h5dclose_f(dset_id, err)
            call h5dcreate_f(group_id, 'ez01', h5t_native_double, dspace_id, dset_id, err)
               call h5dwrite_f(dset_id, h5t_native_double, ez01, dims_vectors, err)
            call h5dclose_f(dset_id, err)
            call h5dcreate_f(group_id, 'ex02', h5t_native_double, dspace_id, dset_id, err)
               call h5dwrite_f(dset_id, h5t_native_double, ex02, dims_vectors, err)
            call h5dclose_f(dset_id, err)
            call h5dcreate_f(group_id, 'ey02', h5t_native_double, dspace_id, dset_id, err)
               call h5dwrite_f(dset_id, h5t_native_double, ey02, dims_vectors, err)
            call h5dclose_f(dset_id, err)
            call h5dcreate_f(group_id, 'ez02', h5t_native_double, dspace_id, dset_id, err)
               call h5dwrite_f(dset_id, h5t_native_double, ez02, dims_vectors, err)
            call h5dclose_f(dset_id, err)
            call h5dcreate_f(group_id, 'rot_axis01', h5t_native_double, dspace_id, dset_id, err)
               call h5dwrite_f(dset_id, h5t_native_double, rot_axis01, dims_vectors, err)
            call h5dclose_f(dset_id, err)
            call h5dcreate_f(group_id, 'rot_axis02', h5t_native_double, dspace_id, dset_id, err)
               call h5dwrite_f(dset_id, h5t_native_double, rot_axis02, dims_vectors, err)
            call h5dclose_f(dset_id, err)                              
         call h5sclose_f(dspace_id, err)            
     call h5gclose_f(group_id, err)
!
!-----------------------star 1------------------------------------------   
!
      dims_r = (/ cs1_nr /)
      dims_theta = (/ cs1_ntheta /)
      dims_phi = (/ cs1_nphi /)
      dims_3d = (/ cs1_nr, cs1_ntheta, cs1_nphi /)
      
      call h5gcreate_f(file_id, 'star1', group_id0, err)
!
!----------------------------dimensions---------------------------------
!
         call h5gcreate_f(group_id0, 'dimensions', group_id, err)
            call h5screate_simple_f(1, dims_scalars, aspace_id, err)
               call h5acreate_f(group_id, 'nr', h5t_native_integer, aspace_id, &
                                attr_id, err)
                  call h5awrite_f(attr_id, h5t_native_integer, cs1_nr, dims_scalars, err)
               call h5aclose_f(attr_id, err)
               call h5acreate_f(group_id, 'ntheta', h5t_native_integer, aspace_id, &
                                attr_id, err)
                  call h5awrite_f(attr_id, h5t_native_integer, cs1_ntheta, dims_scalars, err)
               call h5aclose_f(attr_id, err)
               call h5acreate_f(group_id, 'nphi', h5t_native_integer, aspace_id, &
                                attr_id, err)
                  call h5awrite_f(attr_id, h5t_native_integer, cs1_nphi, dims_scalars, err)
               call h5aclose_f(attr_id, err)
            call h5sclose_f(aspace_id, err)
         call h5gclose_f(group_id, err)
!
!----------------------------coordinates--------------------------------
!
!radius in rstar 1         
         call h5gcreate_f(group_id0, 'coordinates', group_id, err)
            call h5screate_simple_f(1, dims_r, dspace_id, err)
               call h5dcreate_f(group_id, 'r', h5t_native_double, dspace_id, dset_id, err)
                  call h5dwrite_f(dset_id, h5t_native_double, cs1_r, dims_r, err)
               call h5dclose_f(dset_id, err)
            call h5sclose_f(dspace_id, err)
            call h5screate_simple_f(1, dims_theta, dspace_id, err)
               call h5dcreate_f(group_id, 'theta', h5t_native_double, dspace_id, dset_id, err)
                  call h5dwrite_f(dset_id, h5t_native_double, cs1_theta, dims_theta, err)
               call h5dclose_f(dset_id, err)
            call h5sclose_f(dspace_id, err)
            call h5screate_simple_f(1, dims_phi, dspace_id, err)
               call h5dcreate_f(group_id, 'phi', h5t_native_double, dspace_id, dset_id, err)
                  call h5dwrite_f(dset_id, h5t_native_double, cs1_phi, dims_phi, err)
               call h5dclose_f(dset_id, err)
            call h5sclose_f(dspace_id, err)
        call h5gclose_f(group_id, err)
!
!----------------------------3d solution--------------------------------
!
        call h5gcreate_f(group_id0, 'solution3d', group_id, err)
            call h5screate_simple_f(3, dims_3d, dspace_id, err)
               call h5dcreate_f(group_id, 'sline3d', h5t_native_double, dspace_id, dset_id, err)
                  call h5dwrite_f(dset_id, h5t_native_double, cs1_sline3d, dims_3d, err)
               call h5dclose_f(dset_id, err)
               call h5dcreate_f(group_id, 'scont3d', h5t_native_double, dspace_id, dset_id, err)
                  call h5dwrite_f(dset_id, h5t_native_double, cs1_scont3d, dims_3d, err)
               call h5dclose_f(dset_id, err)
            call h5sclose_f(dspace_id, err)
        call h5gclose_f(group_id, err)
!
!--------------------------------3d model-------------------------------
!
!temperature in kelvin
!opalbar in cgs
!opac in cgs
         cs1_opac3d=cs1_opac3d/sr1
         cs1_opalbar3d=cs1_opalbar3d/sr1
         call h5gcreate_f(group_id0, 'model3d', group_id, err)
            call h5screate_simple_f(3, dims_3d, dspace_id, err)
               call h5dcreate_f(group_id, 'rho3d', h5t_native_double, dspace_id, dset_id, err)
                  call h5dwrite_f(dset_id, h5t_native_double, cs1_rho3d, dims_3d, err)
               call h5dclose_f(dset_id, err)
               call h5dcreate_f(group_id, 't3d', h5t_native_double, dspace_id, dset_id, err)
                  call h5dwrite_f(dset_id, h5t_native_double, cs1_t3d, dims_3d, err)
               call h5dclose_f(dset_id, err)
               call h5dcreate_f(group_id, 'opac3d', h5t_native_double, dspace_id, dset_id, err)
                  call h5dwrite_f(dset_id, h5t_native_double, cs1_opac3d, dims_3d, err)
               call h5dclose_f(dset_id, err)
               call h5dcreate_f(group_id, 'opalbar3d', h5t_native_double, dspace_id, dset_id, err)
                  call h5dwrite_f(dset_id, h5t_native_double, cs1_opalbar3d, dims_3d, err)
               call h5dclose_f(dset_id, err)
               call h5dcreate_f(group_id, 'velx3d', h5t_native_double, dspace_id, dset_id, err)
                  call h5dwrite_f(dset_id, h5t_native_double, cs1_velx3d, dims_3d, err)
               call h5dclose_f(dset_id, err)
               call h5dcreate_f(group_id, 'vely3d', h5t_native_double, dspace_id, dset_id, err)
                  call h5dwrite_f(dset_id, h5t_native_double, cs1_vely3d, dims_3d, err)
               call h5dclose_f(dset_id, err)
               call h5dcreate_f(group_id, 'velz3d', h5t_native_double, dspace_id, dset_id, err)
                  call h5dwrite_f(dset_id, h5t_native_double, cs1_velz3d, dims_3d, err)
               call h5dclose_f(dset_id, err)
            call h5sclose_f(dspace_id, err)
         call h5gclose_f(group_id, err)
!restore units            
         cs1_opac3d=cs1_opac3d*sr1
         cs1_opalbar3d=cs1_opalbar3d*sr1
!
      call h5gclose_f(group_id0, err)
!
!-----------------------star 2------------------------------------------   
!
      dims_r = (/ cs2_nr /)
      dims_theta = (/ cs2_ntheta /)
      dims_phi = (/ cs2_nphi /)
      dims_3d = (/ cs2_nr, cs2_ntheta, cs2_nphi /)
      
      call h5gcreate_f(file_id, 'star2', group_id0, err)
!
!----------------------------dimensions---------------------------------
!
         call h5gcreate_f(group_id0, 'dimensions', group_id, err)
            call h5screate_simple_f(1, dims_scalars, aspace_id, err)
               call h5acreate_f(group_id, 'nr', h5t_native_integer, aspace_id, &
                                attr_id, err)
                  call h5awrite_f(attr_id, h5t_native_integer, cs2_nr, dims_scalars, err)
               call h5aclose_f(attr_id, err)
               call h5acreate_f(group_id, 'ntheta', h5t_native_integer, aspace_id, &
                                attr_id, err)
                  call h5awrite_f(attr_id, h5t_native_integer, cs2_ntheta, dims_scalars, err)
               call h5aclose_f(attr_id, err)
               call h5acreate_f(group_id, 'nphi', h5t_native_integer, aspace_id, &
                                attr_id, err)
                  call h5awrite_f(attr_id, h5t_native_integer, cs2_nphi, dims_scalars, err)
               call h5aclose_f(attr_id, err)
            call h5sclose_f(aspace_id, err)
         call h5gclose_f(group_id, err)
!
!----------------------------coordinates--------------------------------
!
!radius in rstar 12     
         call h5gcreate_f(group_id0, 'coordinates', group_id, err)
            call h5screate_simple_f(1, dims_r, dspace_id, err)
               call h5dcreate_f(group_id, 'r', h5t_native_double, dspace_id, dset_id, err)
                  call h5dwrite_f(dset_id, h5t_native_double, cs2_r, dims_r, err)
               call h5dclose_f(dset_id, err)
            call h5sclose_f(dspace_id, err)
            call h5screate_simple_f(1, dims_theta, dspace_id, err)
               call h5dcreate_f(group_id, 'theta', h5t_native_double, dspace_id, dset_id, err)
                  call h5dwrite_f(dset_id, h5t_native_double, cs2_theta, dims_theta, err)
               call h5dclose_f(dset_id, err)
            call h5sclose_f(dspace_id, err)
            call h5screate_simple_f(1, dims_phi, dspace_id, err)
               call h5dcreate_f(group_id, 'phi', h5t_native_double, dspace_id, dset_id, err)
                  call h5dwrite_f(dset_id, h5t_native_double, cs2_phi, dims_phi, err)
               call h5dclose_f(dset_id, err)
            call h5sclose_f(dspace_id, err)
        call h5gclose_f(group_id, err)
!
!----------------------------3d solution--------------------------------
!
        call h5gcreate_f(group_id0, 'solution3d', group_id, err)
            call h5screate_simple_f(3, dims_3d, dspace_id, err)
               call h5dcreate_f(group_id, 'sline3d', h5t_native_double, dspace_id, dset_id, err)
                  call h5dwrite_f(dset_id, h5t_native_double, cs2_sline3d, dims_3d, err)
               call h5dclose_f(dset_id, err)
               call h5dcreate_f(group_id, 'scont3d', h5t_native_double, dspace_id, dset_id, err)
                  call h5dwrite_f(dset_id, h5t_native_double, cs2_scont3d, dims_3d, err)
               call h5dclose_f(dset_id, err)
            call h5sclose_f(dspace_id, err)
        call h5gclose_f(group_id, err)
!
!--------------------------------3d model-------------------------------
!
!temperature in kelvin
!opalbar in cgs
!opac in cgs
!velocity in cm/s        
         cs2_opac3d=cs2_opac3d/sr2
         cs2_opalbar3d=cs2_opalbar3d/sr2
         call h5gcreate_f(group_id0, 'model3d', group_id, err)
            call h5screate_simple_f(3, dims_3d, dspace_id, err)
               call h5dcreate_f(group_id, 'rho3d', h5t_native_double, dspace_id, dset_id, err)
                  call h5dwrite_f(dset_id, h5t_native_double, cs2_rho3d, dims_3d, err)
               call h5dclose_f(dset_id, err)
               call h5dcreate_f(group_id, 't3d', h5t_native_double, dspace_id, dset_id, err)
                  call h5dwrite_f(dset_id, h5t_native_double, cs2_t3d, dims_3d, err)
               call h5dclose_f(dset_id, err)
               call h5dcreate_f(group_id, 'opac3d', h5t_native_double, dspace_id, dset_id, err)
                  call h5dwrite_f(dset_id, h5t_native_double, cs2_opac3d, dims_3d, err)
               call h5dclose_f(dset_id, err)
               call h5dcreate_f(group_id, 'opalbar3d', h5t_native_double, dspace_id, dset_id, err)
                  call h5dwrite_f(dset_id, h5t_native_double, cs2_opalbar3d, dims_3d, err)
               call h5dclose_f(dset_id, err)
               call h5dcreate_f(group_id, 'velx3d', h5t_native_double, dspace_id, dset_id, err)
                  call h5dwrite_f(dset_id, h5t_native_double, cs2_velx3d, dims_3d, err)
               call h5dclose_f(dset_id, err)
               call h5dcreate_f(group_id, 'vely3d', h5t_native_double, dspace_id, dset_id, err)
                  call h5dwrite_f(dset_id, h5t_native_double, cs2_vely3d, dims_3d, err)
               call h5dclose_f(dset_id, err)
               call h5dcreate_f(group_id, 'velz3d', h5t_native_double, dspace_id, dset_id, err)
                  call h5dwrite_f(dset_id, h5t_native_double, cs2_velz3d, dims_3d, err)
               call h5dclose_f(dset_id, err)
            call h5sclose_f(dspace_id, err)
         call h5gclose_f(group_id, err)
!restore units            
         cs2_opac3d=cs2_opac3d*sr2
         cs2_opalbar3d=cs2_opalbar3d*sr2
!
      call h5gclose_f(group_id0, err)
!
!-----------------------------------------------------------------------
!
   call h5fclose_f(file_id, err)
call h5close_f(err)
!
!-----------------------------------------------------------------------
!
write(*,'(a20,es20.8)') 'xnue0', xnue0
write(*,'(a20,es20.8)') 'vth_fiducial', vth_fiducial
write(*,'(a20,i20)') 'na', na
write(*,*)
!
write(*,*) '-------------------star 1---------------------'
write(*,'(a20,i20)') 'nr', cs1_nr
write(*,'(a20,i20)') 'ntheta', cs1_ntheta
write(*,'(a20,i20)') 'nphi', cs1_nphi
write(*,'(a20,es20.8)') 'teff', teff1
write(*,'(a20,es20.8)') 'trad', trad1
write(*,'(a20,es20.8)') 'rstar', rstar1
write(*,'(a20,es20.8)') 'logg', logg1
write(*,'(a20,es20.8)') 'lstar', lstar1
write(*,'(a20,es20.8)') 'yhe', yhe1
write(*,'(a20,es20.8)') 'vrot', vrot1
write(*,'(a20,es20.8)') 'vmicro', vmicro1
write(*,'(a20,es20.8)') 'rmin', rmin1
write(*,'(a20,es20.8)') 'rmax', rmax1
write(*,'(a20,es20.8)') 'x01', x01
write(*,'(a20,es20.8)') 'y01', y01
write(*,'(a20,es20.8)') 'z01', z01
write(*,'(a20,es20.8)') 'vx01', vx01
write(*,'(a20,es20.8)') 'vy01', vy01
write(*,'(a20,es20.8)') 'vz01', vz01
write(*,*)
!
j=(cs1_ntheta-1)/2+1
k=1
write(*,'(a20,2f8.4)') 'at (theta,phi)= ', cs1_theta(j)*180./pi, cs1_phi(k)*180./pi
write(*,'(10a20)') 'r', 'T(r)', 'rho(r)', 'opac(r)', 'opalbar(r)', &
                  'velx(r)', 'vely(r)', 'velz(r)', 'sline(r)', 'scont(r)'
write(*,'(10a20)') '[rstar1]', '[K]', '[g/cm^3]', '[1/sr1]', '[1/s/sr1]', &
                  '[cm/s]', '[cm/s]', '[cm/s]', '[erg/s/cm^2/Hz]', '[erg/s/cm^2/Hz]'
do i=1, cs1_nr
   write(*,'(10es20.8)') cs1_r(i), cs1_t3d(i,j,k), cs1_rho3d(i,j,k), cs1_opac3d(i,j,k), cs1_opalbar3d(i,j,k), cs1_velx3d(i,j,k), &
                        cs1_vely3d(i,j,k), cs1_velz3d(i,j,k), cs1_sline3d(i,j,k), cs1_scont3d(i,j,k)
enddo

j=1
k=1
write(*,'(a20,2f8.4)') 'at (theta,phi)= ', cs1_theta(j)*180./pi, cs1_phi(k)*180./pi
write(*,'(10a20)') 'r', 'T(r)', 'rho(r)', 'opac(r)', 'opalbar(r)', &
                  'velx(r)', 'vely(r)', 'velz(r)', 'sline(r)', 'scont(r)'
write(*,'(10a20)') '[rstar1]', '[K]', '[g/cm^3]', '[1/sr1]', '[1/s/sr1]', &
                  '[cm/s]', '[cm/s]', '[cm/s]', '[erg/s/cm^2/Hz]', '[erg/s/cm^2/Hz]'
do i=1, cs1_nr
   write(*,'(10es20.8)') cs1_r(i), cs1_t3d(i,j,k), cs1_rho3d(i,j,k), cs1_opac3d(i,j,k), cs1_opalbar3d(i,j,k), cs1_velx3d(i,j,k), &
                        cs1_vely3d(i,j,k), cs1_velz3d(i,j,k), cs1_sline3d(i,j,k), cs1_scont3d(i,j,k)
enddo
!
write(*,*)
write(*,*)
!
write(*,*) '-------------------star 2---------------------'
write(*,'(a20,i20)') 'nr', cs2_nr
write(*,'(a20,i20)') 'ntheta', cs2_ntheta
write(*,'(a20,i20)') 'nphi', cs2_nphi
write(*,'(a20,es20.8)') 'teff', teff2
write(*,'(a20,es20.8)') 'trad', trad2
write(*,'(a20,es20.8)') 'rstar', rstar2
write(*,'(a20,es20.8)') 'logg', logg2
write(*,'(a20,es20.8)') 'lstar', lstar2
write(*,'(a20,es20.8)') 'yhe', yhe2
write(*,'(a20,es20.8)') 'vrot', vrot2
write(*,'(a20,es20.8)') 'vmicro', vmicro2
write(*,'(a20,es20.8)') 'rmin', rmin2
write(*,'(a20,es20.8)') 'rmax', rmax2
write(*,'(a20,es20.8)') 'x02', x02
write(*,'(a20,es20.8)') 'y02', y02
write(*,'(a20,es20.8)') 'z02', z02
write(*,'(a20,es20.8)') 'vx02', vx02
write(*,'(a20,es20.8)') 'vy02', vy02
write(*,'(a20,es20.8)') 'vz02', vz02
write(*,*)
!
j=(cs2_ntheta-1)/2+1
k=1!cs2_nphi/2+1
write(*,'(a20,2f8.4)') 'at (theta,phi)= ', cs2_theta(j)*180./pi, cs2_phi(k)*180./pi
write(*,'(10a20)') 'r', 'T(r)', 'rho(r)', 'opac(r)', 'opalbar(r)', &
                  'velx(r)', 'vely(r)', 'velz(r)', 'sline(r)', 'scont(r)'
write(*,'(10a20)') '[rstar2]', '[K]', '[g/cm^3]', '[1/sr2]', '[1/s/sr2]', &
                  '[cm/s]', '[cm/s]', '[cm/s]', '[erg/s/cm^2/Hz]', '[erg/s/cm^2/Hz]'

do i=1, cs2_nr
   write(*,'(10es20.8)') cs2_r(i), cs2_t3d(i,j,k), cs2_rho3d(i,j,k), cs2_opac3d(i,j,k), cs2_opalbar3d(i,j,k), cs2_velx3d(i,j,k), &
                        cs2_vely3d(i,j,k), cs2_velz3d(i,j,k), cs2_sline3d(i,j,k), cs2_scont3d(i,j,k)
enddo

j=1
k=1!cs2_nphi/2+1
write(*,'(a20,2f8.4)') 'at (theta,phi)= ', cs2_theta(j)*180./pi, cs2_phi(k)*180./pi
write(*,'(10a20)') 'r', 'T(r)', 'rho(r)', 'opac(r)', 'opalbar(r)', &
                  'velx(r)', 'vely(r)', 'velz(r)', 'sline(r)', 'scont(r)'
write(*,'(10a20)') '[rstar2]', '[K]', '[g/cm^3]', '[1/sr2]', '[1/s/sr2]', &
                  '[cm/s]', '[cm/s]', '[cm/s]', '[erg/s/cm^2/Hz]', '[erg/s/cm^2/Hz]'

do i=1, cs2_nr
   write(*,'(10es20.8)') cs2_r(i), cs2_t3d(i,j,k), cs2_rho3d(i,j,k), cs2_opac3d(i,j,k), cs2_opalbar3d(i,j,k), cs2_velx3d(i,j,k), &
                        cs2_vely3d(i,j,k), cs2_velz3d(i,j,k), cs2_sline3d(i,j,k), cs2_scont3d(i,j,k)
enddo
!
!
!
end subroutine output3d_spc
