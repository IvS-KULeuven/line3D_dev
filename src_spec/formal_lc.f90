!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine formal_ray(z_ray, opalbar_ray, opac_ray, sline_ray, scont_ray, profile_ray, nz_ray, iin, iin_c, iem, iem_c, iemi, iabs)
!
!----------------calculates intensity along a given ray-----------------
!input: z_ray          spatial coordinate along ray
!       opalbar_ray    frequency integrated line opacity along ray
!       opac_ray       continuum opacity along ray
!       sline_ray      line source function along ray
!       scont_ray      continuum source function along ray
!       profile_ray    profile function along ray
!       nz_ray         number of data points of ray
!       iin            boundary intensity
!
!output: iem           total emergent intensity (from the line and continuum)
!        iem_c         emergent intensity (from the continuum)
!        iemi          only emission part of profile
!        iabs          only absorption part of profile
!-----------------------------------------------------------------------
!
use prog_type
use fund_const, only: zero  
!
implicit none
!
! ... arguments
integer(i4b), intent(in) :: nz_ray
real(dp), dimension(nz_ray), intent(in) :: z_ray, opalbar_ray, opac_ray, &
                                           scont_ray, sline_ray, profile_ray
real(dp), intent(in) :: iin, iin_c
real(dp) :: iem, iem_c, iemi, iabs
!
! ... local scalars
integer(i4b) :: i
real(dp) :: idum, idum1, idum2, idumc, idumc1, idumc2, opalbar1, opalbar2, &
            opac1, opac2, opatot1, opatot2, dtau, tau, dtauc, dtaul, tauc, dumc1, &
            dumc2, duml1, duml2, scont1, scont2, sline1, sline2
real(dp) :: idum1_abs, idum1_emi, idum2_abs, idum2_emi
!
!-----------------------------------------------------------------------
!
!first step
idum=0.d0
idum1=iin
idum2=0.d0
idumc=0.d0
idumc1=iin_c
idumc2=0.d0
tau=0.d0
tauc=0.d0
!
!**debug
idum1_abs=iin
idum2_abs=0.d0
idum1_emi=0.d0
idum2_emi=0.d0
!**
!
opalbar1 = opalbar_ray(1) * profile_ray(1)
opac1  = opac_ray(1)
opatot1 = opac1 + opalbar1
!
scont1 = scont_ray(1)
sline1 = sline_ray(1)
!
!open(1, file='trash/iray.dat', form='formatted')
!
do i=2, nz_ray
   opalbar2 = opalbar_ray(i) * profile_ray(i)
   opac2 = opac_ray(i)
   opatot2 = opac2 + opalbar2
!
   scont2 = scont_ray(i)
   sline2 = sline_ray(i)
!
!without continuum
!    dtau= (opalbar2+opalbar1)*(z_ray(i-1)-z_ray(i))/2.d0
!    idum = idum + (1.d0-exp(dtau)) * (sline_ray(i-1)+sline_ray(i)) / 2.d0
!
!with continuum
   dtau= (opatot2+opatot1)*(z_ray(i-1)-z_ray(i))/2.d0
   dtauc = (opac2+opac1)*(z_ray(i-1)-z_ray(i))/2.d0
   dtaul = (opalbar2+opalbar1)*(z_ray(i-1)-z_ray(i))/2.d0
!note: need to distinguish four cases in order that division by zero never occurs!!!
!      (e.g. if outside information region...)
!   if(opac_ray(i-1).eq.0.d0) then
!      dumc1=0.d0
!   else
!      dumc1=(opac_ray(i-1) * scont_ray(i-1)) / (opac_ray(i-1) + opalbar_ray(i-1))
!   endif
!!
!   if(opac_ray(i).eq.0.d0) then
!      dumc2=0.d0
!   else
!      dumc2=(opac_ray(i)   * scont_ray(i))   / (opac_ray(i)   + opalbar_ray(i))
!   endif
!!
!   if(opalbar_ray(i-1).eq.0.d0) then
!      duml1=0.d0
!   else
!      duml1=(opalbar_ray(i-1) * sline_ray(i-1)) / (opac_ray(i-1) + opalbar_ray(i-1))
!   endif
!!
!   if(opalbar_ray(i).eq.0.d0) then
!      duml2=0.d0
!   else
!      duml2=(opalbar_ray(i)   * sline_ray(i))   / (opac_ray(i)   + opalbar_ray(i))   
!   endif

   if(opac1.eq.zero) then
      dumc1 = zero
   else
      dumc1 = scont1*opac1/opatot1
   endif
!
   if(opac2.eq.zero) then
      dumc2 = zero
   else
      dumc2 = scont2*opac2/opatot2
   endif
!
   if(opalbar1.eq.zero) then
      duml1 = zero
   else
      duml1 = sline1*opalbar1/opatot1
   endif
!
   if(opalbar2.eq.zero) then
      duml2 = zero
   else
      duml2 = sline2*opalbar2/opatot2
   endif

   
   idum2 = idum1 * exp(dtau) + (1.d0-exp(dtau)) * (dumc1+dumc2+duml1+duml2)/2.d0
   idumc2 = idumc1 * exp(dtauc) + (1.d0-exp(dtauc)) * (scont_ray(i-1)+scont_ray(i)) / 2.d0

!   write(*,'(10es20.8)') z_ray(i-1)-z_ray(i), sline_ray(i), opalbar_ray(i), profile_ray(i), dtau, opatot2, opatot1, idum2
!   if(idum2.lt.0.) stop 'go in in formal_ray'   
!
!only source contribution
!   idum = idum + (1.d0-exp(dtau)) * (dumc1+dumc2+duml1+duml2)/2.d0
!   idumc = idumc + (1.d0-exp(dtauc)) * (scont_ray(i-1)+scont_ray(i)) / 2.d0
!
!only absorption part, and only emission part (for line+continuum)
   idum2_abs = idum1_abs*exp(dtau)
   idum2_emi = idum1_emi*exp(dtau) +  (1.d0-exp(dtau)) * (dumc1+dumc2+duml1+duml2)/2.d0
   idum1_emi=idum2_emi
   idum1_abs=idum2_abs


!
   tau=tau-dtau
   tauc=tauc-dtauc
   opalbar1 = opalbar2
   opac1 = opac2
   opatot1 = opatot2
   sline1 = sline2
   scont1 = scont2
   idum1=idum2
   idumc1=idumc2

!   write(*,'(6es20.8)') z_ray(i-1), z_ray(i), sline_ray(i), opalbar_ray(i), profile_ray(i), dtau
!   if(idum1.lt.0.) stop 'go in in formal_ray'
!   write(*,*) z_ray(i), idum2, profile_ray(i)
!

!
enddo
!
!close(1)
!
!additionally contribution from incident intensity
!note: iin, iin_c is set to zero automatically in setup_ray if noncore-ray
!iem = iin*exp(-tau) + idum
!iem_c = iin*exp(-tauc) + idumc
iem=idum2
iem_c=idumc2
!
iemi=idum2_emi
iabs=idum2_abs
!
!
end subroutine formal_ray
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine formal_ray2(z_ray, opalbar_ray, opac_ray, sline_ray, scont_ray, profile_ray, nz_ray, iin, iin_c, int1d, iemi1d, iabs1d, tau1d)
!
!----------------calculates intensity along a given ray-----------------
!input: z_ray          spatial coordinate along ray
!       opalbar_ray    frequency integrated line opacity along ray
!       opac_ray       continuum opacity along ray
!       sline_ray      line source function along ray
!       scont_ray      continuum source function along ray
!       profile_ray    profile function along ray
!       nz_ray         number of data points of ray
!       iin            boundary intensity
!
!output: int1d         intensity along ray
!        iemi1d        only emission part
!        iabs1d        only absorption part
!        tau1d         opatical depth
!-----------------------------------------------------------------------
!
use prog_type
!
implicit none
!
! ... arguments
integer(i4b), intent(in) :: nz_ray
real(dp), dimension(nz_ray), intent(in) :: z_ray, opalbar_ray, opac_ray, &
                                           scont_ray, sline_ray, profile_ray
real(dp), intent(in) :: iin, iin_c
real(dp), dimension(nz_ray) :: int1d, iemi1d, iabs1d, tau1d
!
! ... local scalars
integer(i4b) :: i
real(dp) :: idum, idum1, idum2, idumc, idumc1, idumc2, opalbar1, opalbar2, &
            opac1, opac2, opatot1, opatot2, dtau, tau, dtauc, tauc, dumc1, &
            dumc2, duml1, duml2
real(dp) :: idum1_abs, idum1_emi, idum2_abs, idum2_emi
!
!-----------------------------------------------------------------------
!
!first step
idum=0.d0
idum1=iin
idum2=0.d0
idumc=0.d0
idumc1=iin_c
idumc2=0.d0
tau=0.d0
tauc=0.d0
!
!**debug
idum1_abs=iin
idum2_abs=0.d0
idum1_emi=0.d0
idum2_emi=0.d0
!**
!
opalbar1 = opalbar_ray(1) * profile_ray(1)
opac1  = opac_ray(1)
opatot1 = opac1 + opalbar1
!
int1d(1)=iin
iemi1d(1)=0.d0
iabs1d(1)=iin
tau1d(1)=0.d0
!
!
do i=2, nz_ray
   opalbar2 = opalbar_ray(i) * profile_ray(i)
   opac2 = opac_ray(i)
   opatot2 = opac2 + opalbar2
!
!without continuum
!    dtau= (opalbar2+opalbar1)*(z_ray(i-1)-z_ray(i))/2.d0
!    idum = idum + (1.d0-exp(dtau)) * (sline_ray(i-1)+sline_ray(i)) / 2.d0
!
!with continuum
   dtau= (opatot2+opatot1)*(z_ray(i-1)-z_ray(i))/2.d0
   dtauc= (opac2+opac1)*(z_ray(i-1)-z_ray(i))/2.d0

!note: need to distinguish four cases in order that division by zero never occurs!!!
!      (e.g. if outside information region...)
   if(opac_ray(i-1).eq.0.d0) then
      dumc1=0.d0
   else
      dumc1=(opac_ray(i-1) * scont_ray(i-1)) / (opac_ray(i-1) + opalbar_ray(i-1))
   endif
!
   if(opac_ray(i).eq.0.d0) then
      dumc2=0.d0
   else
      dumc2=(opac_ray(i)   * scont_ray(i))   / (opac_ray(i)   + opalbar_ray(i))
   endif
!
   if(opalbar_ray(i-1).eq.0.d0) then
      duml1=0.d0
   else
      duml1=(opalbar_ray(i-1) * sline_ray(i-1)) / (opac_ray(i-1) + opalbar_ray(i-1))
   endif
!
   if(opalbar_ray(i).eq.0.d0) then
      duml2=0.d0
   else
      duml2=(opalbar_ray(i)   * sline_ray(i))   / (opac_ray(i)   + opalbar_ray(i))   
   endif
!
   idum2 = idum1 * exp(dtau) + (1.d0-exp(dtau)) * (dumc1+dumc2+duml1+duml2)/2.d0
   idumc2 = idumc1 * exp(dtauc) + (1.d0-exp(dtauc)) * (scont_ray(i-1)+scont_ray(i)) / 2.d0
!
!only source contribution
!   idum = idum + (1.d0-exp(dtau)) * (dumc1+dumc2+duml1+duml2)/2.d0
!   idumc = idumc + (1.d0-exp(dtauc)) * (scont_ray(i-1)+scont_ray(i)) / 2.d0
!
!**debug
!only absorption part, and only emission part
   idum2_abs = idum1_abs*exp(dtau)
   idum2_emi = idum1_emi*exp(dtau) +  (1.d0-exp(dtau)) * (dumc1+dumc2+duml1+duml2)/2.d0
   idum1_emi=idum2_emi
   idum1_abs=idum2_abs
!**
!
   tau=tau-dtau
   tauc=tauc-dtauc
   opalbar1 = opalbar2
   opac1 = opac2
   opatot1 = opatot2
   idum1=idum2
   idumc1=idumc2

   tau1d(i)=tau
   int1d(i)=idum2
   iemi1d(i)=idum2_emi
   iabs1d(i)=idum2_abs

!   write(1,'(6es20.8)') z_ray(i), sline_ray(i), opalbar2, idum2, idumc2, profile_ray(i)
!   write(*,*) z_ray(i), idum2, profile_ray(i)
!

!
enddo
!
!reverse array, such that tau is measured from outside to inside
do i=1, nz_ray
   tau1d(i) = tau1d(nz_ray) - tau1d(i)
enddo
!
!
!
end subroutine formal_ray2
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine print_err(iem, iem_theo, ierr)
!
use prog_type
!
implicit none
!
! ... arguments
real(dp), intent(in) :: iem, iem_theo
real(dp) :: ierr
!
ierr=(iem-iem_theo)/iem_theo
!
write(*,'(3(a20))') 'iem', 'iem (theo)', 'relative error'
write(*,'(3(e20.8))') iem, iem_theo, ierr
write(*,*)
!
end subroutine print_err
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine allocate_fs1d
!
!----------------set up one arbitrary ray (for test purposes)-----------
!
use prog_type
use mod_spectrum, only: nz_ray, z_ray, opalbar_ray, opac_ray, sline_ray, scont_ray, &
                  velz_ray, profile_ray, vth_ray, temp_ray
!
implicit none
!
! ... arguments
!
! ... local scalars
integer(i4b) :: err
!
!-----------------------------------------------------------------------
!
if(allocated(z_ray)) deallocate(z_ray)
if(allocated(opalbar_ray)) deallocate(opalbar_ray)
if(allocated(opac_ray)) deallocate(opac_ray)
if(allocated(sline_ray)) deallocate(sline_ray)
if(allocated(scont_ray)) deallocate(scont_ray)
if(allocated(velz_ray)) deallocate(velz_ray)
if(allocated(vth_ray)) deallocate(vth_ray)
if(allocated(temp_ray)) deallocate(temp_ray)
if(allocated(profile_ray)) deallocate(profile_ray)
!
!-----------------------------------------------------------------------
!
allocate(z_ray(nz_ray), stat=err)
   if(err.ne.0) stop 'allocation error allocate_fs1d: z_ray'
allocate(opalbar_ray(nz_ray), stat=err)
   if(err.ne.0) stop 'allocation error allocate_fs1d: opalbar_ray'
allocate(opac_ray(nz_ray), stat=err)
   if(err.ne.0) stop 'allocation error allocate_fs1d: opac_ray'
allocate(sline_ray(nz_ray), stat=err)
   if(err.ne.0) stop 'allocation error allocate_fs1d: sline_ray'
allocate(scont_ray(nz_ray), stat=err)
   if(err.ne.0) stop 'allocation error allocate_fs1d: scont_ray'
allocate(velz_ray(nz_ray), stat=err)
   if(err.ne.0) stop 'allocation error allocate_fs1d: velz_ray'
allocate(vth_ray(nz_ray), stat=err)
   if(err.ne.0) stop 'allocation error allocate_fs1d: vth_ray'
allocate(temp_ray(nz_ray), stat=err)
   if(err.ne.0) stop 'allocation error allocate_fs1d: temp_ray'
allocate(profile_ray(nz_ray), stat=err)
   if(err.ne.0) stop 'allocation error allocate_fs1d: profile_ray'


end subroutine allocate_fs1d

!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine deallocate_fs1d
!
!----------------set up one arbitrary ray (for test purposes)-----------
!
use prog_type
use mod_spectrum, only: nz_ray, z_ray, opalbar_ray, opac_ray, sline_ray, scont_ray, &
                  velz_ray, profile_ray, vth_ray, temp_ray
!
implicit none
!
! ... arguments
!
! ... local scalars
!
!-----------------------------------------------------------------------
!
deallocate(z_ray)
deallocate(opalbar_ray)
deallocate(opac_ray)
deallocate(sline_ray)
deallocate(scont_ray)
deallocate(velz_ray)
deallocate(vth_ray)
deallocate(temp_ray)
deallocate(profile_ray)
!
!
!
end subroutine deallocate_fs1d
