!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine info_region(xp, yp, zp, rmin, rlim, linfo_phot, linfo_max, linfo_boundary)
!
!-----------------------------------------------------------------------
!
!      check if a point (xp, yp, zp) lies within region
!                where information is stored
!
!   input: coordinates:                   xcoord, ycoord, zcoord
!          minimum and maximum radius     rmin, rlim
!   output: logical linfo_phot: information region for minimum radius
!                   linfo_max:  information region for maximum radius
!                   linfo_boundary: if boundary point is exactly hit
!
!-----------------------------------------------------------------------
!
use prog_type
!
implicit none
!
! ... arguments
real(dp), intent(in) :: xp, yp, zp, rmin, rlim
logical :: linfo_phot, linfo_max, linfo_boundary
!
! ... local scalars
real(dp) :: rad
!
! ... local functions
!
rad=sqrt(xp**2+yp**2+zp**2)
!
!be careful here: used in main program as well as in formal solver!!!
linfo_phot=.true.
if(rad.lt.rmin) linfo_phot=.false.
!
linfo_max=.true.
if(rad.gt.rlim) linfo_max=.false.
!
linfo_boundary=.false.
if(rad.eq.rmin) linfo_boundary=.true.
!
!
end subroutine info_region
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine info_region_rot(xp, yp, zp, sma, smb, smc, rlim, linfo_phot, linfo_max, linfo_boundary)
!
!-----------------------------------------------------------------------
!
!      check if a point (xp, yp, zp) lies within region
!                where information is stored
! (accounting for stellar surface distortion, where the stellar surface 
!    is approximated by ellipsoid with semi-major-axes sma, smb, smc)
!          
!
!   input: coordinates:       xcoord, ycoord, zcoord
!          maximum radius:    rlim
!          semi-major-axes:   sma, smb, smc
!
!   output: logical linfo_phot: information region for minimum radius
!                   linfo_max:  information region for maximum radius
!                   linfo_boundary: if boundary point is exactly hit
!
!-----------------------------------------------------------------------
!
use prog_type
!
implicit none
!
! ... arguments
real(dp), intent(in) :: xp, yp, zp, sma, smb, smc, rlim
logical :: linfo_phot, linfo_max, linfo_boundary
!
! ... local scalars
real(dp) :: rad, fdum
!
! ... local functions
!
fdum = sqrt((xp/sma)**2 + (yp/smb)**2 + (zp/smc)**2)
rad = sqrt(xp**2 + yp**2 + zp**2)
!
linfo_phot=.true.
if(fdum.lt.1.d0-1.d-14) linfo_phot=.false.
!write(*,'(a10,es30.20)') 'fdum', fdum
!
linfo_max=.true.
if(rad.gt.rlim) linfo_max=.false.
!
linfo_boundary=.false.
if(abs(fdum-1.d0).lt.1.d-14) linfo_boundary=.true.
!
!
end subroutine info_region_rot
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine info_region2(xp, yp, zp, rmin, rlim, linfo_phot, linfo_max, linfo_boundary)
!
!-----------------------------------------------------------------------
!
!      check if a point (xp, yp, zp) lies within region
!                where information is stored
!
!   input: coordinates:                   xcoord, ycoord, zcoord
!          minimum and maximum radius     rmin, rlim
!   output: logical linfo_phot: information region for minimum radius
!                   linfo_max:  information region for maximum radius
!                   linfo_boundary: if boundary point is exactly hit
!
!
!-----------------------------------------------------------------------
!
use prog_type
!
implicit none
!
! ... arguments
real(dp), intent(in) :: xp, yp, zp, rmin, rlim
logical :: linfo_phot, linfo_max, linfo_boundary
!
! ... local scalars
real(dp) :: rad
!
! ... local functions
!
rad=sqrt(xp**2+yp**2+zp**2)
!
!be careful here: used in main program as well as in formal solver!!!
linfo_phot=.true.
if(rad-rmin.lt.-1.d-14) linfo_phot=.false.
!
linfo_max=.true.
if(rad-rlim.gt.1.d-13) linfo_max=.false.
!
linfo_boundary=.false.
if(rad.eq.rmin) linfo_boundary=.true.
!
!
end subroutine info_region2
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine info_region_spc(rad, rmin, rlim, linfo_phot, linfo_max)
!
!-----------------------------------------------------------------------
!
!      check if a point (with given radius rad) lies within region
!                    where information is stored
!
!   input: radial coordinate:             rad
!          minimum and maximum radius     rmin, rlim
!   output: logical linfo_phot: information region on photosphere
!                   linfo_max:  information region for maximum radius
!
!-----------------------------------------------------------------------
!
use prog_type
!
implicit none
!
! ... arguments
real(dp), intent(in) :: rad, rmin, rlim
logical :: linfo_phot, linfo_max
!
! ... local scalars
!
! ... local functions
!
!be careful here: used in main program as well as in formal solver!!!
if(rad+1.d-14.lt.rmin) then
   linfo_phot=.false.
else
   linfo_phot=.true.
endif
!
!
if(rad-rlim.gt.1.d-14) then
!outside of information region
   linfo_max=.false.
else
   linfo_max=.true.
endif
!
!
!
end subroutine info_region_spc
