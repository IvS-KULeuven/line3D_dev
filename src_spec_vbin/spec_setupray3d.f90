!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine setup_ray3d_cs1spc(indx_point, cs1_x, cs1_y, lcore)
!
!------sets up a p-ray along z in coordinate system of star 1-----------
!-----------------------for a given xobs--------------------------------
!--------coordinate system of star 1 in spherical coordiantes-----------
!
!to calculate for star1:
!   1. z-ray         from z_ray=sqrt(r^2-p^2)
!   2. z-ray-coordinates in cartesian coordinate system
!   3. opalbar_ray   trilinear interpolation from given 3d-grid
!   4. sline_ray     trilinear interpolation from given 3d-grid
!   5. velx(y,z)_ray trilinear interpolation from given 3d-grid
!   5. velz_ray      projected velocity along the ray
!
!input: cs0_p         radial distance of point in cylindrical system
!       cs0_x, cs0_y  x,y coordinates of point in cylindrical system
  
!
!output: lcore    logical to describe if core or non-core ray
!        all physical quantities along ray
!        cs1_p         radial distance of point cylindrical system of star 1
!        cs1_x, cs1_y  x,y coordinates of point in cylindrical system of star 1
!
use prog_type
use fund_const  
use dime_spec, only:  cs1_nr
use dime_model3d, only: cs1_nr_spc, cs1_ntheta_spc, cs1_nphi_spc, cs1_r_spc, cs1_theta_spc, cs1_phi_spc, &
                        cs1_opalbar3d, cs1_sline3d, cs1_scont3d, cs1_opac3d, cs1_velx3d, cs1_vely3d, cs1_velz3d, cs1_vth3d
use mod_spectrum, only: cs1_nz_ray, cs1_z_ray, cs1_opalbar_ray, cs1_opac_ray, cs1_scont_ray, cs1_sline_ray, &
                        cs1_vth_ray, cs1_temp_ray, cs1_velz_ray, cs1_imask_ray, cs1_r, &
                        unit_length, transmat, transmat_inv, translvec1, rotmat1, rotmat1_inv, nhat, &
                        del_vel2, vphot1_proj
use params_spec, only: vth_fiducial, vth_min, xic1, vmax, tmin, na, vmicro1, rmax0, rmin1, rmax1, rstar1, vrot1
use timing_spec, only: t_setup1
use mod_triangles, only: points_xcoord, points_ycoord
use omp_lib
use mod_interp3d, only: get_rtp_indx, get_rtp_values1, get_rtp_values2, trilin_complete
!
implicit none
!
! ... arguments
integer(i4b), intent(in) :: indx_point
real(dp), intent(out) :: cs1_x, cs1_y
logical, intent(out) :: lcore
!
! ... local scalars
integer(i4b) :: i,j,k
integer(i4b) :: err
integer(i4b) :: indx_r1, indx_r2, indx_t1, indx_t2, indx_p1, indx_p2, indx
integer(i4b) :: iz, iz_dum, nadd
integer(i4b), parameter :: nz_max = 20000
real(dp) :: delz, zmax, zdum, dum_vel1, dum_vel2, dum_delv
real(dp) :: dum_velx, dum_vely, dum_velz, dum_vel, dum_gradv, dum_vmicro, dum_opac, dum_scont 
real(dp) :: r1, r2, theta1, theta2, phi1, phi2
real(dp) :: vala, valb, valc, vald, vale, valf, valg, valh, valp
real(dp) :: velx, vely, velz, velr, velphi, vel_abs
real(dp) :: radp, thetap, phip, fdum
real(dp) :: cs0_x, cs0_y, cs1_p, cs1_z, cs0_z
real(dp) :: tstart, tend

!
! ... local arrays
real(dp), dimension(3) :: vec_cac, vec_cyc, vec_vel
integer(i1b), dimension(nz_max) :: imaskdum_ray, imaskdum2_ray
real(dp), dimension(nz_max) :: zdum_ray, veldum_ray, vthdum_ray, &
                               opacdum_ray, opalbardum_ray, scontdum_ray, slinedum_ray, tempdum_ray
real(dp), dimension(nz_max) :: zdum2_ray, veldum2_ray, vthdum2_ray, opacdum2_ray, opalbardum2_ray, scontdum2_ray, &
                               slinedum2_ray, tempdum2_ray
!
! ... local logicals
logical :: expol, linfo_phot, linfo_max, llogr, llogt, llogp, llogf, lr2, ldum
!
! ... local functions
real(dp) :: interpol_yp
real(dp) :: interpol_yp_spline
real(dp) :: calc_vmicro
real(dp) :: vthermal
logical :: boundary
!
tstart = omp_get_wtime()
!
!------------------transform given cs0_x, cs0_y-------------------------
!----------------to coordinate system of star 1-------------------------
!
cs0_x = points_xcoord(indx_point)
cs0_y = points_ycoord(indx_point)
cs0_z = zero
!
call transform_cs0_cs1(cs0_x, cs0_y, cs0_z, unit_length, rstar1, translvec1, transmat, transmat_inv, cs1_x, cs1_y, cs1_z)

cs1_p=sqrt(cs1_x**2+cs1_y**2)
!
!-----------------------------------------------------------------------
!
if(cs1_p.ge.rmax1) then
!if outside of region where information is stored, set to dummy values
   cs1_nz_ray=2
   lcore=.false.
!
!get bounds of global coordinate system in units of local coordinate system
   zmax = rmax0*unit_length/rstar1
!   
!store dummy values
   call allocate_fs1dcs1
   cs1_z_ray(1) = -zmax
   cs1_z_ray(2) = zmax
   cs1_imask_ray(1) = 0
   cs1_imask_ray(2) = 2


   cs1_opalbar_ray = zero
   cs1_opac_ray = zero
   cs1_scont_ray = zero   
   cs1_sline_ray = zero
   cs1_vth_ray = 1.d-8
   cs1_temp_ray = zero   
   cs1_velz_ray = zero
!
else
!
!--------------get the first guess of the z-vector----------------------
!----------------in cylindrical system of star 1------------------------

   iz=0
!
   do i=cs1_nr, 1, -1
      zdum=cs1_r(i)*cs1_r(i)-cs1_p*cs1_p
      if(zdum.gt.zero) then
         iz=iz+1
         zdum_ray(iz)=sqrt(zdum)
      else
         iz=iz+1
         zdum_ray(iz)=zero
         exit
      endif
   enddo
!
!check if core ray
   lcore=.false.
!
!inner-most point in carthesian coordinates of star 1
   vec_cyc(1)=cs1_x
   vec_cyc(2)=cs1_y
   vec_cyc(3)=zdum_ray(iz)
   vec_cac=matmul(transmat, vec_cyc)
!transform to tilted coordinate system of star 1
   vec_cac=matmul(rotmat1_inv, vec_cac)
!
!note: function can in general calculate arbitrary boundary shapes, here: spherical
   lcore=boundary(vec_cac)
!
!
   if(.not.lcore) then
!set non core rays
      iz_dum=2*iz-1
      do i=1, iz-1
         zdum_ray(iz+i)=-zdum_ray(iz-i)
      enddo
      iz=iz_dum
      vphot1_proj=zero
   else
!calculation of photospheric velocity (only rotational velocity) projected onto ray
!for the moment: assume rotation axis always in global z-direction (no tilts)
      vec_cac=matmul(transmat, vec_cyc)      
      radp=sqrt(vec_cac(1)**2 + vec_cac(2)**2 + vec_cac(3)**2)
      call vel_phot(vec_cac(1), vec_cac(2), vec_cac(3), radp, velx, vely, velz, vrot1)
      vphot1_proj=nhat(1)*velx + nhat(2)*vely + nhat(3)*velz
!      write(*,*) vrot1, velx, vely, velz, vphot1_proj      
   endif
!
!-----------------------------------------------------------------------
!
   do i=1, iz
!
!calculate z_ray in carthesian coordinates of star 1
      vec_cyc(1)=cs1_x
      vec_cyc(2)=cs1_y
      vec_cyc(3)=zdum_ray(i)
      vec_cac=matmul(transmat, vec_cyc)
!transform to tilted coordinate system
      vec_cac=matmul(rotmat1_inv, vec_cac)
!
!check if point of ray lies within region where information is stored
!note: info_region2 is used, because vec_cac on inner boundary is not exactly r=1
      call info_region2(vec_cac(1), vec_cac(2), vec_cac(3), rmin1, rmax1, linfo_phot, linfo_max, ldum)
!
      if(i.eq.iz) then
         if(linfo_phot.eqv..false.)  then
            write(*,*) vec_cac(1)**2+vec_cac(2)**2+vec_cac(3)**2, linfo_phot, linfo_max
            stop 'error in setup_ray3d_cs1spc: linfo_phot is false'
         endif
      endif
!interpolation only if point lies within region where information is stored
      if(linfo_phot.and.linfo_max) then
!
!calculate corresponding points in spherical coordinate system
         radp=sqrt(vec_cac(1)**2 + vec_cac(2)**2 + vec_cac(3)**2)
         call get_angles_spc(vec_cac(1), vec_cac(2), vec_cac(3), thetap, phip)
!
!search for indices of a the surrounding grid-cell (or neighbouring grid-cells for extrapolation)
         call get_rtp_indx(radp, thetap, phip, cs1_r_spc, cs1_theta_spc, cs1_phi_spc, &
                           cs1_nr_spc, cs1_ntheta_spc, cs1_nphi_spc, &
                           indx_r1, indx_r2, indx_t1, indx_t2, indx_p1, indx_p2, expol, rmin1, rmax1)
!
!get coordinates and radii of a cube surrounding point p of interest
!llogx, llogy, llogz are flags to interpolate in logspace, if allowed
         call get_rtp_values1(radp, thetap, phip, cs1_r_spc, cs1_theta_spc, cs1_phi_spc, cs1_nr_spc, cs1_ntheta_spc, cs1_nphi_spc, &
                             indx_r1, indx_r2, indx_t1, indx_t2, indx_p1, indx_p2, &
                             r1, r2, theta1, theta2, phi1, phi2, &
                             llogr, llogt, llogp)
!
!----------------------interpolation of velocity components-------------
!
!get velx on cell vertices
         call get_rtp_values2(cs1_nr_spc, cs1_ntheta_spc, cs1_nphi_spc, cs1_velx3d, &
                              indx_r1, indx_r2, indx_t1, indx_t2, indx_p1, indx_p2, &
                              vala, valb, valc, vald, vale, valf, valg, valh, llogf)
         llogr=.false.
         llogt=.false.
         llogp=.false.
         llogf=.false.
         lr2=.false.
!actual interpolation
         call trilin_complete(radp, thetap, phip, &
                              r1, r2, theta1, theta2, phi1, phi2, &
                              vala, valb, valc, vald, vale, valf, valg, valh, &
                              r1, r2, r1, r2, r1, r2, r1, r2, radp, &
                              expol, .true., llogr, llogt, llogp, llogf, lr2, velx)
!
!get vely on cell vertices
         call get_rtp_values2(cs1_nr_spc, cs1_ntheta_spc, cs1_nphi_spc, cs1_vely3d, &
                              indx_r1, indx_r2, indx_t1, indx_t2, indx_p1, indx_p2, &
                              vala, valb, valc, vald, vale, valf, valg, valh, llogf)
         llogr=.false.
         llogt=.false.
         llogp=.false.
         llogf=.false.
         lr2=.false.
!actual interpolation
         call trilin_complete(radp, thetap, phip, &
                              r1, r2, theta1, theta2, phi1, phi2, &
                              vala, valb, valc, vald, vale, valf, valg, valh, &
                              r1, r2, r1, r2, r1, r2, r1, r2, radp, &
                              expol, .true., llogr, llogt, llogp, llogf, lr2, vely)
!
!get velz on cell vertices
         call get_rtp_values2(cs1_nr_spc, cs1_ntheta_spc, cs1_nphi_spc, cs1_velz3d, &
                              indx_r1, indx_r2, indx_t1, indx_t2, indx_p1, indx_p2, &
                              vala, valb, valc, vald, vale, valf, valg, valh, llogf)
         llogr=.false.
         llogt=.false.
         llogp=.false.
         llogf=.false.
         lr2=.false.
!actual interpolation
         call trilin_complete(radp, thetap, phip, &
                              r1, r2, theta1, theta2, phi1, phi2, &
                              vala, valb, valc, vald, vale, valf, valg, valh, &
                              r1, r2, r1, r2, r1, r2, r1, r2, radp, &
                              expol, .true., llogr, llogt, llogp, llogf, lr2, velz)
!
!for debugging
!      fdum=0.d0
!      call vel_t3(radp*sin(thetap)*cos(phip), radp*sin(thetap)*sin(phip), radp*cos(thetap), radp, velx, vely, velz, fdum)

!transform to formal-solver carthesian coordinates
         vec_vel = matmul(rotmat1, (/ velx, vely, velz /))  !(/ velx, vely, velz /) !
!
!calculation of velocity projected onto ray
         veldum_ray(i)=nhat(1)*vec_vel(1) + nhat(2)*vec_vel(2) + nhat(3)*vec_vel(3)
!
!----------------------interpolation of line opacity--------------------
!
!get line-opacities on cube vertices
         call get_rtp_values2(cs1_nr_spc, cs1_ntheta_spc, cs1_nphi_spc, cs1_opalbar3d, &
                              indx_r1, indx_r2, indx_t1, indx_t2, indx_p1, indx_p2, &
                              vala, valb, valc, vald, vale, valf, valg, valh, llogf)
!         llogr=.false.
!         llogt=.false.
!         llogp=.false.
!         llogf=.false.
         lr2=.true.
!actual interpolation
         call trilin_complete(radp, thetap, phip, &
                              r1, r2, theta1, theta2, phi1, phi2, &
                              vala, valb, valc, vald, vale, valf, valg, valh, &
                              r1, r2, r1, r2, r1, r2, r1, r2, radp, &
                              expol, .true., llogr, llogt, llogp, llogf, lr2, valp)
!line opacity in units of 1/unit_length
         opalbardum_ray(i) = valp
!
!------------------------line source function---------------------------
!
!get sline on cube vertices
         call get_rtp_values2(cs1_nr_spc, cs1_ntheta_spc, cs1_nphi_spc, cs1_sline3d, &
                              indx_r1, indx_r2, indx_t1, indx_t2, indx_p1, indx_p2, &
                              vala, valb, valc, vald, vale, valf, valg, valh, llogf)
!         llogr=.false.
!         llogt=.false.
!         llogp=.false.
!         llogf=.false.
         lr2=.false.
!actual interpolation
         call trilin_complete(radp, thetap, phip, &
                              r1, r2, theta1, theta2, phi1, phi2, &
                              vala, valb, valc, vald, vale, valf, valg, valh, &
                              r1, r2, r1, r2, r1, r2, r1, r2, radp, &
                              expol, .true., llogr, llogt, llogp, llogf, lr2, valp)
         slinedum_ray(i)=valp
!
!------------------------continuum source function----------------------
!
!get scont on cube vertices
         call get_rtp_values2(cs1_nr_spc, cs1_ntheta_spc, cs1_nphi_spc, cs1_scont3d, &
                              indx_r1, indx_r2, indx_t1, indx_t2, indx_p1, indx_p2, &
                              vala, valb, valc, vald, vale, valf, valg, valh, llogf)
!         llogr=.false.
!         llogt=.false.
!         llogp=.false.
!         llogf=.false.
         lr2=.false.
!actual interpolation
         call trilin_complete(radp, thetap, phip, &
                              r1, r2, theta1, theta2, phi1, phi2, &
                              vala, valb, valc, vald, vale, valf, valg, valh, &
                              r1, r2, r1, r2, r1, r2, r1, r2, radp, &
                              expol, .true., llogr, llogt, llogp, llogf, lr2, valp)
         scontdum_ray(i)=valp
!
!------------------------continuum opacity------------------------------
!
!get opac on cube vertices
         call get_rtp_values2(cs1_nr_spc, cs1_ntheta_spc, cs1_nphi_spc, cs1_opac3d, &
                              indx_r1, indx_r2, indx_t1, indx_t2, indx_p1, indx_p2, &
                              vala, valb, valc, vald, vale, valf, valg, valh, llogf)
!         llogr=.false.
!         llogt=.false.
!         llogp=.false.
!         llogf=.false.
         lr2=.false.
!actual interpolation
         call trilin_complete(radp, thetap, phip, &
                              r1, r2, theta1, theta2, phi1, phi2, &
                              vala, valb, valc, vald, vale, valf, valg, valh, &
                              r1, r2, r1, r2, r1, r2, r1, r2, radp, &
                              expol, .true., llogr, llogt, llogp, llogf, lr2, valp)
         opacdum_ray(i)=valp
!
!-----------------------------------------------------------------------
!
!----temporary solution: set v_micro, temperature, continuum to zero----
!
!calculation of absolute value of velocity (needed to calculate vmicro)
         vel_abs=sqrt(velx**2 + vely**2 + velz**2)*vth_fiducial
!
!note: minimum microturbulent velocity set to input-vmicro
         dum_vmicro = calc_vmicro(vmicro1, vmax*1.d5, vel_abs)
         dum_vmicro = vmicro1
!
!calculate corresponding thermal velocity (for a fixed tmin)
         vthdum_ray(i) = vthermal(dum_vmicro, tmin, na)
!         vthdum_ray(i) = 5.d5   !test low thermal velocities

!for now: set temperature to zero (not needed, only later for correct thermal velocity calculation)
         tempdum_ray(i)=zero
!
!------set all quantities to zero, if outside information region--------
!
      else
         opalbardum_ray(i) = zero
         slinedum_ray(i) = zero
         veldum_ray(i) = zero
         vthdum_ray(i) = 1.d-8
         tempdum_ray(i) = zero
         scontdum_ray(i) = zero
         opacdum_ray(i) = zero
!
      endif

!      write(*,'(5es20.8)') zdum_ray(i), slinedum_ray(i), opacdum_ray(i), &
!                           opalbardum_ray(i), veldum_ray(i)
                                                
!
   enddo
!
!------------------store everything in global arrays--------------------
!
   cs1_nz_ray=iz
!
   call allocate_fs1dcs1
!
   do i=1, cs1_nz_ray
      cs1_z_ray(i) = zdum_ray(cs1_nz_ray+1-i)
      cs1_opalbar_ray(i) = opalbardum_ray(cs1_nz_ray+1-i)
      cs1_opac_ray(i) = opacdum_ray(cs1_nz_ray+1-i)
      cs1_scont_ray(i) = scontdum_ray(cs1_nz_ray+1-i)
      cs1_sline_ray(i) = slinedum_ray(cs1_nz_ray+1-i)
      cs1_vth_ray(i) = vthdum_ray(cs1_nz_ray+1-i)
      cs1_temp_ray(i) = tempdum_ray(cs1_nz_ray+1-i)
      cs1_velz_ray(i) = veldum_ray(cs1_nz_ray+1-i)
      cs1_imask_ray(i) = 1
   enddo
!
!
endif
!
!
!
tend = omp_get_wtime()
t_setup1 = t_setup1 + tend-tstart
!
!
end subroutine setup_ray3d_cs1spc
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine setup_ray3d_cs2spc(indx_point, cs2_x, cs2_y, lcore)
!
!------sets up a p-ray along z in coordinate system of star 2-----------
!-----------------------for a given xobs--------------------------------
!--------coordinate system of star 2 in spherical coordiantes-----------
!
!to calculate for star1:
!   1. z-ray         from z_ray=sqrt(r^2-p^2)
!   2. z-ray-coordinates in cartesian coordinate system
!   3. opalbar_ray   trilinear interpolation from given 3d-grid
!   4. sline_ray     trilinear interpolation from given 3d-grid
!   5. velx(y,z)_ray trilinear interpolation from given 3d-grid
!   5. velz_ray      projected velocity along the ray
!
!input: cs0_p         radial distance of point in cylindrical system
!       cs0_x, cs0_y  x,y coordinates of point in cylindrical system
!        cs2_p         radial distance of point cylindrical system of star 2
!        cs2_x, cs2_y  x,y coordinates of point in cylindrical system of star 2
  
!
!output: lcore    logical to describe if core or non-core ray
!        all physical quantities along ray
!
use prog_type
use fund_const  
use dime_spec, only:  cs2_nr
use dime_model3d, only: cs2_nr_spc, cs2_ntheta_spc, cs2_nphi_spc, cs2_r_spc, cs2_theta_spc, cs2_phi_spc, &
                        cs2_opalbar3d, cs2_sline3d, cs2_scont3d, cs2_opac3d, cs2_velx3d, cs2_vely3d, cs2_velz3d, cs2_vth3d
use mod_spectrum, only: cs2_nz_ray, cs2_z_ray, cs2_opalbar_ray, cs2_opac_ray, cs2_scont_ray, cs2_sline_ray, &
                        cs2_vth_ray, cs2_temp_ray, cs2_velz_ray, cs2_imask_ray, cs2_r, &
                        unit_length, transmat, transmat_inv, translvec2, rotmat2, rotmat2_inv, nhat, &
                        del_vel2, vphot2_proj
use params_spec, only: vth_fiducial, vth_min, xic2, vmax, tmin, na, vmicro2, rmax0, rmin2, rmax2, rstar2, vrot2
use timing_spec, only: t_setup2
use mod_triangles, only: points_xcoord, points_ycoord
use omp_lib
use mod_interp3d, only: get_rtp_indx, get_rtp_values1, get_rtp_values2, trilin_complete
!
implicit none
!
! ... arguments
integer(i4b), intent(in) :: indx_point
real(dp), intent(out) :: cs2_x, cs2_y
logical, intent(out) :: lcore
!
! ... local scalars
integer(i4b) :: i,j,k
integer(i4b) :: err
integer(i4b) :: indx_r1, indx_r2, indx_t1, indx_t2, indx_p1, indx_p2, indx
integer(i4b) :: iz, iz_dum, nadd
integer(i4b), parameter :: nz_max = 20000
real(dp) :: delz, zmax, zdum, dum_vel1, dum_vel2, dum_delv
real(dp) :: dum_velx, dum_vely, dum_velz, dum_vel, dum_gradv, dum_vmicro, dum_opac, dum_scont 
real(dp) :: r1, r2, theta1, theta2, phi1, phi2
real(dp) :: vala, valb, valc, vald, vale, valf, valg, valh, valp
real(dp) :: velx, vely, velz, velr, velphi, vel_abs
real(dp) :: radp, thetap, phip, fdum
real(dp) :: cs2_p, cs0_x, cs0_y, cs0_z, cs2_z
real(dp) :: tstart, tend
!
! ... local arrays
real(dp), dimension(3) :: vec_cac, vec_cyc, vec_vel
real(dp), dimension(nz_max) :: zdum_ray, veldum_ray, vthdum_ray, &
                               opacdum_ray, opalbardum_ray, scontdum_ray, slinedum_ray, tempdum_ray
real(dp), dimension(nz_max) :: zdum2_ray, veldum2_ray, vthdum2_ray, opacdum2_ray, opalbardum2_ray, scontdum2_ray, &
                               slinedum2_ray, tempdum2_ray
!
! ... local logicals
logical :: expol, linfo_phot, linfo_max, llogr, llogt, llogp, llogf, lr2, ldum
!
! ... local functions
real(dp) :: interpol_yp
real(dp) :: interpol_yp_spline
real(dp) :: calc_vmicro
real(dp) :: vthermal
logical :: boundary
!
tstart = omp_get_wtime()
!
!-------------transform given cs0_p, cs0_x, cs0_y-----------------------
!----------------to coordinate system of star 2-------------------------
!
cs0_x = points_xcoord(indx_point)
cs0_y = points_ycoord(indx_point)
cs0_z = zero
!
call transform_cs0_cs1(cs0_x, cs0_y, cs0_z, unit_length, rstar2, translvec2, transmat, transmat_inv, cs2_x, cs2_y, cs2_z)

cs2_p=sqrt(cs2_x**2+cs2_y**2)
!
!-----------------------------------------------------------------------
!
if(cs2_p.ge.rmax2) then
!if outside of region where information is stored, set to dummy values
   cs2_nz_ray=2
   lcore=.false.
!
!get bounds of global coordinate system in units of local coordinate system
   zmax = rmax0*unit_length/rstar2
!   
!store dummy values
   call allocate_fs1dcs2
   cs2_z_ray(1) = -zmax
   cs2_z_ray(2) = zmax
   cs2_imask_ray(1) = 0
   cs2_imask_ray(2) = 2

   cs2_opalbar_ray = zero
   cs2_opac_ray = zero
   cs2_scont_ray = zero   
   cs2_sline_ray = zero
   cs2_vth_ray = 1.d-8
   cs2_temp_ray = zero   
   cs2_velz_ray = zero
!
   else
!
!--------------get the first guess of the z-vector----------------------
!----------------in cylindrical system of star 1------------------------
!
   iz=0
!
   do i=cs2_nr, 1, -1
      zdum=cs2_r(i)*cs2_r(i)-cs2_p*cs2_p
      if(zdum.gt.zero) then
         iz=iz+1
         zdum_ray(iz)=sqrt(zdum)
      else
         iz=iz+1
         zdum_ray(iz)=zero
         exit
      endif
   enddo
!
!check if core ray
   lcore=.false.
!
!inner-most point in carthesian coordinates of star 2
   vec_cyc(1)=cs2_x
   vec_cyc(2)=cs2_y
   vec_cyc(3)=zdum_ray(iz)
   vec_cac=matmul(transmat, vec_cyc)
!transform to tilted coordinate system of star 2
   vec_cac=matmul(rotmat2_inv, vec_cac)   
!
!note: function can in general calculate arbitrary boundary shapes, here: spherical
   lcore=boundary(vec_cac)
!
!
   if(.not.lcore) then
!set non core rays
      iz_dum=2*iz-1
      do i=1, iz-1
         zdum_ray(iz+i)=-zdum_ray(iz-i)
      enddo
      iz=iz_dum
      vphot2_proj=zero
   else
!calculation of photospheric velocity (only rotational velocity) projected onto ray
!for the moment: assume rotation axis aligned with global z-axis (no tilts)
      vec_cac=matmul(transmat, vec_cyc)      
      radp=sqrt(vec_cac(1)**2 + vec_cac(2)**2 + vec_cac(3)**2)
      call vel_phot(vec_cac(1), vec_cac(2), vec_cac(3), radp, velx, vely, velz, vrot2)
      vphot2_proj=nhat(1)*velx + nhat(2)*vely + nhat(3)*velz
   endif
!
!-----------------------------------------------------------------------
!
   do i=1, iz
!
!calculate z_ray in carthesian coordinates of star 2
      vec_cyc(1)=cs2_x
      vec_cyc(2)=cs2_y
      vec_cyc(3)=zdum_ray(i)
      vec_cac=matmul(transmat, vec_cyc)
!transform to tilted coordinate system of star 2
      vec_cac=matmul(rotmat2_inv, vec_cac)      
!
!check if point of ray lies within region where information is stored
!note: info_region2 is used, because vec_cac on inner boundary is not exactly r=1
      call info_region2(vec_cac(1), vec_cac(2), vec_cac(3), rmin2, rmax2, linfo_phot, linfo_max, ldum)
!
      if(i.eq.iz) then
         if(linfo_phot.eqv..false.) then
            write(*,*) vec_cac(1)**2+vec_cac(2)**2+vec_cac(3)**2, linfo_phot, linfo_max
            stop 'error'
         endif
      endif
!interpolation only if point lies within region where information is stored
      if(linfo_phot.and.linfo_max) then
!
!calculate corresponding points in spherical coordinate system
         radp=sqrt(vec_cac(1)**2 + vec_cac(2)**2 + vec_cac(3)**2)
         call get_angles_spc(vec_cac(1), vec_cac(2), vec_cac(3), thetap, phip)
!
!search for indices of a the surrounding grid-cell (or neighbouring grid-cells for extrapolation)
         call get_rtp_indx(radp, thetap, phip, cs2_r_spc, cs2_theta_spc, cs2_phi_spc, cs2_nr_spc, cs2_ntheta_spc, cs2_nphi_spc, &
                           indx_r1, indx_r2, indx_t1, indx_t2, indx_p1, indx_p2, expol, rmin2, rmax2)
!
!get coordinates and radii of a cube surrounding point p of interest
!llogx, llogy, llogz are flags to interpolate in logspace, if allowed
         call get_rtp_values1(radp, thetap, phip, cs2_r_spc, cs2_theta_spc, cs2_phi_spc, cs2_nr_spc, cs2_ntheta_spc, cs2_nphi_spc, &
                             indx_r1, indx_r2, indx_t1, indx_t2, indx_p1, indx_p2, &
                             r1, r2, theta1, theta2, phi1, phi2, &
                             llogr, llogt, llogp)
!
!----------------------interpolation of velocity components-------------
!
!get velx on cell vertices
         call get_rtp_values2(cs2_nr_spc, cs2_ntheta_spc, cs2_nphi_spc, cs2_velx3d, &
                              indx_r1, indx_r2, indx_t1, indx_t2, indx_p1, indx_p2, &
                              vala, valb, valc, vald, vale, valf, valg, valh, llogf)
         llogr=.false.
         llogt=.false.
         llogp=.false.
         llogf=.false.
         lr2=.false.
!actual interpolation
         call trilin_complete(radp, thetap, phip, &
                              r1, r2, theta1, theta2, phi1, phi2, &
                              vala, valb, valc, vald, vale, valf, valg, valh, &
                              r1, r2, r1, r2, r1, r2, r1, r2, radp, &
                              expol, .true., llogr, llogt, llogp, llogf, lr2, velx)
!
!get vely on cell vertices
         call get_rtp_values2(cs2_nr_spc, cs2_ntheta_spc, cs2_nphi_spc, cs2_vely3d, &
                              indx_r1, indx_r2, indx_t1, indx_t2, indx_p1, indx_p2, &
                              vala, valb, valc, vald, vale, valf, valg, valh, llogf)
         llogr=.false.
         llogt=.false.
         llogp=.false.
         llogf=.false.
         lr2=.false.
!actual interpolation
         call trilin_complete(radp, thetap, phip, &
                              r1, r2, theta1, theta2, phi1, phi2, &
                              vala, valb, valc, vald, vale, valf, valg, valh, &
                              r1, r2, r1, r2, r1, r2, r1, r2, radp, &
                              expol, .true., llogr, llogt, llogp, llogf, lr2, vely)
!
!get velz on cell vertices
         call get_rtp_values2(cs2_nr_spc, cs2_ntheta_spc, cs2_nphi_spc, cs2_velz3d, &
                              indx_r1, indx_r2, indx_t1, indx_t2, indx_p1, indx_p2, &
                              vala, valb, valc, vald, vale, valf, valg, valh, llogf)
         llogr=.false.
         llogt=.false.
         llogp=.false.
         llogf=.false.
         lr2=.false.
!actual interpolation
         call trilin_complete(radp, thetap, phip, &
                              r1, r2, theta1, theta2, phi1, phi2, &
                              vala, valb, valc, vald, vale, valf, valg, valh, &
                              r1, r2, r1, r2, r1, r2, r1, r2, radp, &
                              expol, .true., llogr, llogt, llogp, llogf, lr2, velz)
!
!for debugging
!      fdum=0.d0
!      call vel_t3(radp*sin(thetap)*cos(phip), radp*sin(thetap)*sin(phip), radp*cos(thetap), radp, velx, vely, velz, fdum)

!transform to formal-solver carthesian coordinates
         vec_vel = matmul(rotmat2, (/ velx, vely, velz /)) !(/ velx, vely, velz /) !
!
!calculation of velocity projected onto ray
         veldum_ray(i)=nhat(1)*vec_vel(1) + nhat(2)*vec_vel(2) + nhat(3)*vec_vel(3)
!
!----------------------interpolation of line opacity--------------------
!
!get line-opacities on cube vertices
         call get_rtp_values2(cs2_nr_spc, cs2_ntheta_spc, cs2_nphi_spc, cs2_opalbar3d, &
                              indx_r1, indx_r2, indx_t1, indx_t2, indx_p1, indx_p2, &
                              vala, valb, valc, vald, vale, valf, valg, valh, llogf)
!         llogr=.false.
!         llogt=.false.
!         llogp=.false.
!         llogf=.false.
         lr2=.true.
!actual interpolation
         call trilin_complete(radp, thetap, phip, &
                              r1, r2, theta1, theta2, phi1, phi2, &
                              vala, valb, valc, vald, vale, valf, valg, valh, &
                              r1, r2, r1, r2, r1, r2, r1, r2, radp, &
                              expol, .true., llogr, llogt, llogp, llogf, lr2, valp)
!line opacity in units of  1/unit_length
         opalbardum_ray(i) = valp
!         write(*,*) r1, r2, indx_t1, indx_t2, cs2_theta_spc(indx_t1), cs2_theta_spc(indx_t2)!, thetap*180./pi, theta1*180./pi, theta2*180./pi, phi1*180./pi, phi2*180./pi, llogr
!
!------------------------line source function---------------------------
!
!get sline on cube vertices
         call get_rtp_values2(cs2_nr_spc, cs2_ntheta_spc, cs2_nphi_spc, cs2_sline3d, &
                              indx_r1, indx_r2, indx_t1, indx_t2, indx_p1, indx_p2, &
                              vala, valb, valc, vald, vale, valf, valg, valh, llogf)
!         llogr=.false.
!         llogt=.false.
!         llogp=.false.
!         llogf=.false.
         lr2=.false.
!actual interpolation
         call trilin_complete(radp, thetap, phip, &
                              r1, r2, theta1, theta2, phi1, phi2, &
                              vala, valb, valc, vald, vale, valf, valg, valh, &
                              r1, r2, r1, r2, r1, r2, r1, r2, radp, &
                              expol, .true., llogr, llogt, llogp, llogf, lr2, valp)
         slinedum_ray(i)=valp
!
!----------------------interpolation of continuum opacity---------------
!
!get line-opacities on cube vertices
         call get_rtp_values2(cs2_nr_spc, cs2_ntheta_spc, cs2_nphi_spc, cs2_opac3d, &
                              indx_r1, indx_r2, indx_t1, indx_t2, indx_p1, indx_p2, &
                              vala, valb, valc, vald, vale, valf, valg, valh, llogf)
!         llogr=.false.
!         llogt=.false.
!         llogp=.false.
!         llogf=.false.
         lr2=.true.
!actual interpolation
         call trilin_complete(radp, thetap, phip, &
                              r1, r2, theta1, theta2, phi1, phi2, &
                              vala, valb, valc, vald, vale, valf, valg, valh, &
                              r1, r2, r1, r2, r1, r2, r1, r2, radp, &
                              expol, .true., llogr, llogt, llogp, llogf, lr2, valp)
!line opacity in units of  1/unit_length
         opacdum_ray(i) = valp
!
!------------------------line source function---------------------------
!
!get sline on cube vertices
         call get_rtp_values2(cs2_nr_spc, cs2_ntheta_spc, cs2_nphi_spc, cs2_scont3d, &
                              indx_r1, indx_r2, indx_t1, indx_t2, indx_p1, indx_p2, &
                              vala, valb, valc, vald, vale, valf, valg, valh, llogf)
!         llogr=.false.
!         llogt=.false.
!         llogp=.false.
!         llogf=.false.
         lr2=.false.
!actual interpolation
         call trilin_complete(radp, thetap, phip, &
                              r1, r2, theta1, theta2, phi1, phi2, &
                              vala, valb, valc, vald, vale, valf, valg, valh, &
                              r1, r2, r1, r2, r1, r2, r1, r2, radp, &
                              expol, .true., llogr, llogt, llogp, llogf, lr2, valp)
         scontdum_ray(i)=valp
!
!-----------------------------------------------------------------------
!
!----temporary solution: set v_micro, temperature, continuum to zero----
!
!calculation of absolute value of velocity (needed to calculate vmicro)
         vel_abs=sqrt(velx**2 + vely**2 + velz**2)*vth_fiducial
!
!note: minimum microturbulent velocity set to input-vmicro
         dum_vmicro = calc_vmicro(vmicro2, vmax*1.d5, vel_abs)
         dum_vmicro = vmicro2
!
!calculate corresponding thermal velocity (for a fixed tmin)
         vthdum_ray(i) = vthermal(dum_vmicro, tmin, na)
!      vthdum_ray(i) = 5.d5   !test low thermal velocities

!for now: set temperature to zero (not needed, only later for correct thermal velocity calculation)
         tempdum_ray(i)=zero
!
!------set all quantities to zero, if outside information region--------
!
      else
         opalbardum_ray(i) = zero
         slinedum_ray(i) = zero
         veldum_ray(i) = zero
         vthdum_ray(i) = 1.d-8
         tempdum_ray(i) = zero
         scontdum_ray(i) = zero
         opacdum_ray(i) = zero
!
      endif

!      write(*,'(5es20.8)') zdum_ray(i), slinedum_ray(i), opacdum_ray(i), opalbardum_ray(i), &
!                           veldum_ray(i)
                                                
!
   enddo
!
!------------------store everything in global arrays--------------------
!
   cs2_nz_ray=iz
!
   call allocate_fs1dcs2
!
!
   do i=1, cs2_nz_ray
!      write(*,*) i, cs2_nz_ray, cs2_nz_ray+1-i
      cs2_z_ray(i) = zdum_ray(cs2_nz_ray+1-i)
      cs2_opalbar_ray(i) = opalbardum_ray(cs2_nz_ray+1-i)
      cs2_opac_ray(i) = opacdum_ray(cs2_nz_ray+1-i)
      cs2_scont_ray(i) = scontdum_ray(cs2_nz_ray+1-i)
      cs2_sline_ray(i) = slinedum_ray(cs2_nz_ray+1-i)
      cs2_vth_ray(i) = vthdum_ray(cs2_nz_ray+1-i)
      cs2_temp_ray(i) = tempdum_ray(cs2_nz_ray+1-i)
      cs2_velz_ray(i) = veldum_ray(cs2_nz_ray+1-i)
      cs2_imask_ray(i) = 1
   enddo
!
!
!
endif
!
tend = omp_get_wtime()
t_setup2 = t_setup2 + tend-tstart


!stop 'go on here'
!
end subroutine setup_ray3d_cs2spc
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine setup_ray3d_cs0spc(indx_point, cs1_x, cs1_y, cs2_x, cs2_y, lcore1, lcore2)
!
!----sets up a p-ray along z in global cylindrical coordinate system----
!-----------------------for a given xobs--------------------------------
!--------coordinate system of star 2 in spherical coordiantes-----------
!
!to calculate for star1:
!   1. z-ray         from z_ray=sqrt(r^2-p^2)
!   2. z-ray-coordinates in cartesian coordinate system
!   3. opalbar_ray   trilinear interpolation from given 3d-grid
!   4. sline_ray     trilinear interpolation from given 3d-grid
!   5. velx(y,z)_ray trilinear interpolation from given 3d-grid
!   5. velz_ray      projected velocity along the ray
!
!input: 
!       cs0_x, cs0_y  x,y coordinates of point in cylindrical system
!       cs1_x, cs1_y  x,y coordinates of point in cylindrical system  of star 1
!       cs2_x, cs2_y  x,y coordinates of point in cylindrical system  of star 2
!       lcore1    logical to describe if core or non-core ray of star1
!       lcore2    logical to describe if core or non-core ray of star2
!
!output: 
!        all physical quantities along ray in global cylindrical system
!
use prog_type
use fund_const  
use mod_spectrum, only: nz_ray, z_ray, opalbar_ray, opac_ray, scont_ray, sline_ray, &
                        vth_ray, temp_ray, velz_ray, &
                        cs1_nz_ray, cs1_z_ray, cs1_opalbar_ray, cs1_opac_ray, cs1_scont_ray, cs1_sline_ray, &
                        cs1_vth_ray, cs1_temp_ray, cs1_velz_ray, cs1_imask_ray, &
                        cs2_nz_ray, cs2_z_ray, cs2_opalbar_ray, cs2_opac_ray, cs2_scont_ray, cs2_sline_ray, &
                        cs2_vth_ray, cs2_temp_ray, cs2_velz_ray,  cs2_imask_ray, &
                        unit_length, transmat, transmat_inv, translvec1, translvec2, &
                        del_vel2, vphot_proj, vphot1_proj, vphot2_proj, nhat
use params_spec, only: rmax1, rmax2, rstar1, rstar2, vx01, vy01, vz01, vx02, vy02, vz02
use timing_spec, only: t_setup3
use omp_lib
use mod_interp1d, only: interpol_yp, interpol_yp_spline
!
implicit none
!
! ... arguments
integer(i4b), intent(in) :: indx_point
real(dp), intent(in) :: cs1_x, cs1_y, cs2_x, cs2_y
logical, intent(inout) :: lcore1, lcore2 
!
! ... local scalars
integer(i4b) :: i,j,k,kstart
integer(i4b) :: indx, err
integer(i4b) :: i1, i2, iz, iz_dum, nadd
integer(i4b), parameter :: nz_max = 20000
real(dp) :: x1, y1, z1, x2, y2, z2, z1_max, z2_max
real(dp) :: delz, zmax, zdum, dum_vel1, dum_vel2, dum_delv
real(dp) :: cs0_x, cs0_y
real(dp) :: ziim1,  zii, ziip1, znew
real(dp) :: v01_proj, v02_proj
real(dp) :: tstart, tend

!
! ... local arrays
real(dp), dimension(3) :: vec_cac, vec_cyc
integer(i1b), dimension(nz_max) :: imaskdum_ray, imaskdum2_ray
real(dp), dimension(nz_max) :: zdum_ray, veldum_ray, vthdum_ray, &
                               opacdum_ray, opalbardum_ray, scontdum_ray, slinedum_ray, tempdum_ray
real(dp), dimension(nz_max) :: zdum2_ray, veldum2_ray, vthdum2_ray, &
                               opacdum2_ray, opalbardum2_ray, scontdum2_ray, slinedum2_ray, tempdum2_ray
!
! ... local logicals
logical :: llcore1, llcore2
!
! ... local functions
real(dp) :: calc_vmicro
real(dp) :: vthermal
logical :: boundary
!
tstart = omp_get_wtime()
!
!-----------------since projected velocities along z--------------------
!--------------only given in local coordinate systems-------------------
!
v01_proj = nhat(1)*vx01 + nhat(2)*vy01 + nhat(3)*vz01
v02_proj = nhat(1)*vx02 + nhat(2)*vy02 + nhat(3)*vz02
!
vphot1_proj = vphot1_proj + v01_proj
vphot2_proj = vphot2_proj + v02_proj
!
vphot_proj=zero
!
!-------------------------get the total complete z-ray------------------
!
!store total of cs1_nz_ray*cs2_nz_ray points
!
!debug start
!do i=1, cs1_nz_ray
!   call transform_cs1_cs0(cs1_x, cs1_y, cs1_z_ray(i), rstar1, unit_length, translvec1, transmat, transmat_inv, x1, y1, z1)   
!   write(*,*) cs1_z_ray(i), z1, cs1_opalbar_ray(i), cs1_opac_ray(i), cs1_scont_ray(i), cs1_sline_ray(i), cs1_velz_ray(i)
!enddo
!write(*,*)
!do i=1, cs2_nz_ray
!   call transform_cs1_cs0(cs2_x, cs2_y, cs2_z_ray(i), rstar2, unit_length, translvec2, transmat, transmat_inv, x2, y2, z2)      
!   write(*,*) cs2_z_ray(i), z2, cs2_opalbar_ray(i), cs2_opac_ray(i), cs2_scont_ray(i), cs2_sline_ray(i), cs2_velz_ray(i)
!enddo
!write(*,*)
!debug end
!stop
!
!
!set maximum possible z-values of both grids
call transform_cs1_cs0(cs1_x, cs1_y, cs1_z_ray(cs1_nz_ray), rstar1, unit_length, translvec1, transmat, transmat_inv, x1, y1, z1_max)
call transform_cs1_cs0(cs2_x, cs2_y, cs2_z_ray(cs2_nz_ray), rstar2, unit_length, translvec2, transmat, transmat_inv, x2, y2, z2_max)
!
i1=1
i2=1
kstart=1
llcore1=.false.
llcore2=.false.
iz=0
!
!transform coordinates of star 1 and star 2
call transform_cs1_cs0(cs1_x, cs1_y, cs1_z_ray(i1), rstar1, unit_length, translvec1, transmat, transmat_inv, x1, y1, z1)
call transform_cs1_cs0(cs2_x, cs2_y, cs2_z_ray(i2), rstar2, unit_length, translvec2, transmat, transmat_inv, x2, y2, z2)
!
!
!
do i=1, cs1_nz_ray+cs2_nz_ray
!
!-------------store the z-coordinates in increasing order---------------   
!
   if(z1.le.z2) then
!note: due to less equals, next advance z2 will be set with same z-coordiante (will be accounted for later)
!
!---------------set coordinates of point in system 1---------------------
!      
      zdum_ray(i) = z1
      veldum_ray(i) = cs1_velz_ray(i1)+v01_proj
      vthdum_ray(i) = cs1_vth_ray(i1)
      opacdum_ray(i) = cs1_opac_ray(i1)
      opalbardum_ray(i) = cs1_opalbar_ray(i1)
      slinedum_ray(i) = cs1_sline_ray(i1)
      scontdum_ray(i) = cs1_scont_ray(i1)
      tempdum_ray(i) = cs1_temp_ray(i1)
      imaskdum_ray(i) = cs1_imask_ray(i1)
!     if(i1.eq.1) then
!        write(*,*) zdum_ray(i), cs1_opac_ray(i1), opacdum_ray(i)
!     endif
!advance z1-grid
      i1=i1+1
      if(i1.eq.cs1_nz_ray+1) then
!set large dummy value if out of bounds (such that only z2-grid will be advanced in all following steps)
         z1 = z2_max + one
      else
         call transform_cs1_cs0(cs1_x, cs1_y, cs1_z_ray(i1), rstar1, unit_length, translvec1, transmat, transmat_inv, x1, y1, z1)
      endif
!
!advance the z0-grid      
      if(lcore1.and.i1.eq.2) then
!set actual starting point of z-grid if boundary is hit
!and reset counter for total number of considered points
!and set corresponding photospheric velocity         
         llcore1=.true.
         llcore2=.false.
         vphot_proj=vphot1_proj
         kstart=i
         iz=1
      else
         iz=iz+1
      endif

   elseif(z1.gt.z2) then
!
!---------------set coordinates of point in system 2---------------------
!
      zdum_ray(i) = z2
      veldum_ray(i) = cs2_velz_ray(i2)+v02_proj
      vthdum_ray(i) = cs2_vth_ray(i2)
      opacdum_ray(i) = cs2_opac_ray(i2)
      opalbardum_ray(i) = cs2_opalbar_ray(i2)
      slinedum_ray(i) = cs2_sline_ray(i2)
      scontdum_ray(i) = cs2_scont_ray(i2)
      tempdum_ray(i) = cs2_temp_ray(i2)
      imaskdum_ray(i) = cs2_imask_ray(i2)

!      write(*,*) i, iz, z1, z2, i2, zdum_ray(i), imaskdum_ray(i), cs1_nz_ray, cs2_nz_ray
!advance z2-grid
      i2=i2+1
      if(i2.eq.cs2_nz_ray+1) then
!set large dummy value if out of bounds (such that only z1-grid will be advanced in all following steps)
         z2 = z1_max + one
      else
         call transform_cs1_cs0(cs2_x, cs2_y, cs2_z_ray(i2), rstar2, unit_length, translvec2, transmat, transmat_inv, x2, y2, z2)
      endif
!
!advance the z0-grid
      if(lcore2.and.i2.eq.2) then
!set actual starting point of z-grid if boundary is hit
!and reset counter for total number of considered points
         vphot_proj=vphot2_proj
         llcore2=.true.
         llcore1=.false.
         kstart=i
         iz=1      
      else
         iz=iz+1
      endif
!
      
   endif
!
enddo
!
!do i=1, iz
!   write(*,*) i, zdum_ray(i), veldum_ray(i), v02_proj, v01_proj
!enddo
!write(*,*) lcore1, lcore2
!stop 'go on 1'
!
lcore1=llcore1
lcore2=llcore2
!
!-----------------------------------------------------------------------
!
!copy everything to new grid, in order that begin of ray is always at surface (if hit)
zdum2_ray = zero
veldum2_ray = zero
vthdum2_ray = zero
opacdum2_ray = zero
opalbardum2_ray = zero
slinedum2_ray = zero
scontdum2_ray = zero
tempdum2_ray = zero
imaskdum2_ray = 10

indx=1
do i=kstart, kstart+iz-1
   zdum2_ray(indx) = zdum_ray(i)   
   veldum2_ray(indx) = veldum_ray(i)
   vthdum2_ray(indx) = vthdum_ray(i)
   opacdum2_ray(indx) = opacdum_ray(i)
   opalbardum2_ray(indx) = opalbardum_ray(i)
   slinedum2_ray(indx) = slinedum_ray(i)
   scontdum2_ray(indx) = scontdum_ray(i)
   tempdum2_ray(indx) = tempdum_ray(i)
   imaskdum2_ray(indx) = imaskdum_ray(i)
!  write(*,*) kstart, kstart+iz-1, iz, i, indx, zdum2_ray(indx), opacdum_ray(i)
   indx=indx+1
enddo
!do i=1, iz
!   write(*,*) i, zdum2_ray(i)
!enddo
!stop 'go on 2'
!
!--------------set correct terminal velocities--------------------------
!-------------if outside of calculation volume--------------------------
!----------to avoid spurious refinement scheme--------------------------
!
!on left side of the ray
if(imaskdum2_ray(1).eq.0.and.imaskdum2_ray(2).eq.0) then
   veldum2_ray(1)=veldum2_ray(3)
   veldum2_ray(2)=veldum2_ray(3)
elseif(imaskdum2_ray(1).eq.0) then
   veldum2_ray(1)=veldum2_ray(2)
elseif(imaskdum2_ray(2).eq.0) then
   stop 'error in setup_ray3d_cs0spc: mask for outer boundary erroneous, left'
endif

!on right side of the ray
if(imaskdum2_ray(iz).eq.2.and.imaskdum2_ray(iz-1).eq.2) then
   veldum2_ray(iz)=veldum2_ray(iz-2)
   veldum2_ray(iz-1)=veldum2_ray(iz-2)
elseif(imaskdum2_ray(iz).eq.2) then
   veldum2_ray(iz)=veldum2_ray(iz-1)
elseif(imaskdum2_ray(iz-1).eq.2) then
   do i=1, iz
      write(*,*) i, zdum2_ray(i), imaskdum2_ray(i)
   enddo
   stop 'error in setup_ray3d_cs0spc: mask for outer boundary erroneous, right'
endif
!
!-------------------------account for same z-coordinates----------------
!        
do i=2, iz-1
!
   ziim1 = zdum2_ray(i-1)
   zii = zdum2_ray(i)
   ziip1 = zdum2_ray(i+1)   
   if(abs(zii-ziim1).lt.small_number) then
!shift point i towards right and interpolate
      znew = (ziip1+zii)/two
      zdum2_ray(i) = znew
      veldum2_ray(i) = interpol_yp(zii, ziip1, veldum2_ray(i), veldum2_ray(i+1), znew)
      opacdum2_ray(i) = interpol_yp(zii, ziip1, opacdum2_ray(i), opacdum2_ray(i+1), znew)
      opalbardum2_ray(i) = interpol_yp(zii, ziip1, opalbardum2_ray(i), opalbardum2_ray(i+1), znew)
      scontdum2_ray(i) = interpol_yp(zii, ziip1, scontdum2_ray(i), scontdum2_ray(i+1), znew)
      slinedum2_ray(i) = interpol_yp(zii, ziip1, slinedum2_ray(i), slinedum2_ray(i+1), znew)
      tempdum2_ray(i) = interpol_yp(zii, ziip1, tempdum2_ray(i), tempdum2_ray(i+1), znew)
      vthdum2_ray(i) = interpol_yp(zii, ziip1, vthdum2_ray(i), vthdum2_ray(i+1), znew)
   endif   
   if(abs(ziip1-zii).lt.small_number) then
!shift point i towards right and interpolate
      znew = (ziim1+zii)/two
      zdum2_ray(i) = znew
      veldum2_ray(i) = interpol_yp(ziim1, zii, veldum2_ray(i-1), veldum2_ray(i), znew)
      opacdum2_ray(i) = interpol_yp(ziim1, zii, opacdum2_ray(i-1), opacdum2_ray(i), znew)
      opalbardum2_ray(i) = interpol_yp(ziim1, zii, opalbardum2_ray(i-1), opalbardum2_ray(i), znew)
      scontdum2_ray(i) = interpol_yp(ziim1, zii, scontdum2_ray(i-1), scontdum2_ray(i), znew)
      slinedum2_ray(i) = interpol_yp(ziim1, zii, slinedum2_ray(i-1), slinedum2_ray(i), znew)
      tempdum2_ray(i) = interpol_yp(ziim1, zii, tempdum2_ray(i-1), tempdum2_ray(i), znew)
      vthdum2_ray(i) = interpol_yp(ziim1, zii, vthdum2_ray(i-1), vthdum2_ray(i), znew)
   endif
!
enddo   
!   
!----------------check if resonance zones are resolved------------------
!-----------if not resolved: add grid points and interpolate------------
!(will be performed only on very final ray)
!
!velocity at beginning point
dum_vel1=veldum2_ray(1)
!
veldum_ray(1)=veldum2_ray(1)
zdum_ray(1)=zdum2_ray(1)
opacdum_ray(1)=opacdum2_ray(1)
opalbardum_ray(1)=opalbardum2_ray(1)
scontdum_ray(1)=scontdum2_ray(1)
slinedum_ray(1)=slinedum2_ray(1)
tempdum_ray(1)=tempdum2_ray(1)
vthdum_ray(1)=vthdum2_ray(1)
!
k=0
iz_dum=iz
!
do i=2, iz
!
   dum_vel2=veldum2_ray(i)
!
!!calculation of velocity steps
   dum_delv=dum_vel2-dum_vel1
!
!  if(dum_vel1*dum_vel2.lt.0.d0) then
!     write(*,'(a96)') 'error in setup_ray_scp: velocity jump along ray from positive velocities <-> negative velocities'
!     write(*,'(a82)') '   => 1. trust the interpolation scheme for shear velocites, etc (not recommended)'
!     write(*,'(a68)') '         at least, you should check against an analytic velocity law'
!     write(*,'(a82)') '   => 2. use an analytic description of the velocity law, instead of interpolation'
!     stop
!  endif
!
   if(abs(dum_delv).gt.del_vel2) then
!always round up. note: ceiling is only f90
      nadd=ceiling(abs(dum_delv)/del_vel2)
!
      do j=1, nadd-1
!add nadd-1 additional points in velocity space
         veldum_ray(i+k)=veldum2_ray(i-1) + j*(veldum2_ray(i)-veldum2_ray(i-1))/nadd
!interpolation of all variables
         zdum_ray(i+k)=interpol_yp(veldum2_ray(i-1), veldum2_ray(i), zdum2_ray(i-1), zdum2_ray(i), veldum_ray(i+k))
         opacdum_ray(i+k)=interpol_yp(veldum2_ray(i-1), veldum2_ray(i), opacdum2_ray(i-1), opacdum2_ray(i), veldum_ray(i+k))
         opalbardum_ray(i+k)=interpol_yp(veldum2_ray(i-1), veldum2_ray(i), opalbardum2_ray(i-1), opalbardum2_ray(i), veldum_ray(i+k))
         scontdum_ray(i+k)=interpol_yp(veldum2_ray(i-1), veldum2_ray(i), scontdum2_ray(i-1), scontdum2_ray(i), veldum_ray(i+k))
         slinedum_ray(i+k)=interpol_yp(veldum2_ray(i-1), veldum2_ray(i), slinedum2_ray(i-1), slinedum2_ray(i), veldum_ray(i+k))
         tempdum_ray(i+k)=interpol_yp(veldum2_ray(i-1), veldum2_ray(i), tempdum2_ray(i-1), tempdum2_ray(i), veldum_ray(i+k))
         vthdum_ray(i+k)=interpol_yp(veldum2_ray(i-1), veldum2_ray(i), vthdum2_ray(i-1), vthdum2_ray(i), veldum_ray(i+k))
!
         k=k+1
         iz_dum=iz_dum+1
      enddo
   endif

   veldum_ray(i+k)=veldum2_ray(i)
   zdum_ray(i+k)=zdum2_ray(i)
   opacdum_ray(i+k)=opacdum2_ray(i)
   opalbardum_ray(i+k)=opalbardum2_ray(i)
   scontdum_ray(i+k)=scontdum2_ray(i)
   slinedum_ray(i+k)=slinedum2_ray(i)
   tempdum_ray(i+k)=tempdum2_ray(i)
   vthdum_ray(i+k)=vthdum2_ray(i)
!
   dum_vel1=dum_vel2
! 
enddo
!
if(iz_dum.gt.nz_max) stop 'error in setup_ray3d_cs0spc: increase nz_max'
!
!------------------store everything in global arrays--------------------
!
nz_ray=iz_dum
!
call allocate_fs1d
!
do i=1, nz_ray
   z_ray(i) = zdum_ray(i)
   opalbar_ray(i) = opalbardum_ray(i)
   opac_ray(i) = opacdum_ray(i)
   scont_ray(i) = scontdum_ray(i)
   sline_ray(i) = slinedum_ray(i)
   vth_ray(i) = vthdum_ray(i)
   temp_ray(i) = tempdum_ray(i)
   velz_ray(i) = veldum_ray(i)
enddo
!
!stop
!
!do i=1, nz_ray
!   if(vth_ray(i).le.30.d5) then
!      write(*,*) z_ray(i), opalbar_ray(i), velz_ray(i), vth_ray(i)/1.d5
!   endif
!enddo
!write(*,*)
!stop
!
tend = omp_get_wtime()
t_setup3 = t_setup3 + tend-tstart
!
!
!
end subroutine setup_ray3d_cs0spc
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
function boundary(vec)
!
use fund_const  
use prog_type
!
implicit none
!
! ... arguments
real(dp), dimension(3), intent(in) :: vec
logical :: boundary
!
! ... local scalars
real(dp) :: rad
!
rad=sqrt(vec(1)**2+vec(2)**2+vec(3)**2)
!
!note: possible rounding errors
if(rad.lt.one+small_number) then
   boundary=.true.
else
   boundary=.false.
endif
!
end function boundary
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine vel_phot(xcoord, ycoord, zcoord, rad, velx, vely, velz, vrot)
!
!-----------calculate velocity components for velocity-law--------------
!----------------v_phi(r,theta) = vrot * sin(theta)/r-------------------
!
!input: xcoord, ycoord, zcoord
!       vrot in units of vth_fiducial  
!
use prog_type
use params_spec, only: vth_fiducial
!
implicit none
!
! ... arguments
real(dp), intent(in) :: xcoord, ycoord, zcoord, rad, vrot
real(dp) :: velx, vely, velz
!
! ... local scalars
real(dp) :: dum_phi, dum_r
!
! ... local functions
!
!in units of vthfiducial
!dum_r = vmin/rad/vth_fiducial   !if radial velocity is included for photosphere
dum_r = 0.d0                    !if radial velocity is neglected for photosphere
dum_phi = vrot/rad/rad!/vth_fiducial
!
velx = dum_r*xcoord - dum_phi*ycoord
vely = dum_r*ycoord + dum_phi*xcoord
velz = dum_r*zcoord
!
!
end subroutine vel_phot
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine vel_t3(xcoord, ycoord, zcoord, rad, velx, vely, velz, velr, vrot)
!
!-----------calculate velocity components for velocity-law--------------
!----------------v_phi(r,theta) = vrot * sin(theta)/r-------------------
!----------------v_r(r) = beta-velocity-law-----------------------------
!
!note: xcoord = r*sin(theta)*cos(phi)
!      vel_x = velr*sin(theta)*cos(phi) - velphi*sin(phi)
!      vel_y = velr*sin(theta)*sin(phi) + velphi*cos(phi)
!      vel_z = velr*cos(theta)
!   can be reduced to velx = velr/r * xcoord - vrot/r^2 * ycoord
!
!
use prog_type
use params_spec, only: vth_fiducial
!
implicit none
!
! ... arguments
real(dp), intent(in) :: xcoord, ycoord, zcoord, rad, velr, vrot
real(dp) :: velx, vely, velz
!
! ... local scalars
real(dp) :: dum_r, dum_phi, dum_velr
!
! ... local functions
!
!in units of vthfiducial
dum_r = velr / rad
dum_phi = vrot*1.d5/rad/rad/vth_fiducial
!
!for rigid body rotation:
!dum_phi = vrot*1.d5/vthfiducial
!
velx = dum_r*xcoord - dum_phi*ycoord
vely = dum_r*ycoord + dum_phi*xcoord
velz = dum_r*zcoord
!
!write(*,*) rad, sqrt(velx**2 + vely**2 + velz**2)
!
end subroutine vel_t3
