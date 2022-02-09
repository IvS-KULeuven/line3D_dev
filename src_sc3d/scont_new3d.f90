subroutine scont_new3d
!
!-----------------------------------------------------------------------
!---------calculates new iterate of continuum source function-----------
!   different options: classical lambda-iteration
!                      diagonal of lambda-matrix
!                      nearest neighbours
!-----------------------------------------------------------------------
!
use prog_type
use options, only: opt_alo_cont
!
implicit none
!
! ... local scalars
integer(i4b) :: err
!
write(*,*) '--------------calculating new iterate for source function (alo)----------------'
write(*,*)
!
!-----------------------------------------------------------------------
!
select case(opt_alo_cont)
   case(0)
      call scont_new3d_classic
   case(1)
      call scont_new3d_diag
   case(2)
      call scont_new3d_dn
   case(3)
      call scont_new3d_nn
   case default
      stop 'set option opt_alo_cont'
end select
!
!
end subroutine scont_new3d
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine scont_new3d_classic
!
!-----------------------------------------------------------------------
!--------calculates new iterate of continuum source function------------
!-------------------classical lambda-iteration--------------------------
!-----------------------------------------------------------------------
!
use prog_type
use fund_const
use dime3d, only: ndxmax, ndymax, ndzmax, scont3d, mint3d, imask3d, t3d, eps_cont3d
use freq, only: xnue0, lfreqint
!
implicit none
!
! ... local scalars
integer(i4b) :: i, j, k
!
! ... local functions
real(dp) :: bnue2
!
!-----------------------------------------------------------------------
!
do i=1, ndxmax
   do j=1, ndymax
      do k=1, ndzmax
         select case(imask3d(i,j,k))
            case(1,2,3)
               scont3d(i,j,k) = (1.d0-eps_cont3d(i,j,k)) * mint3d(i,j,k) + eps_cont3d(i,j,k)*bnue2(xnue0, t3d(i,j,k), lfreqint)
            case default
         endselect
      enddo
   enddo
enddo
!
end subroutine scont_new3d_classic
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine scont_new3d_diag
!
!-----------------------------------------------------------------------
!--------calculates new iterate of continuum source function------------
!---------approximate lambda iteration using only diagonal--------------
!-----------------------------------------------------------------------
!
use prog_type
use fund_const
use dime3d, only: ndxmax, ndymax, ndzmax, scont3d, alocont_nn3d, mint3d, imask3d, t3d, eps_cont3d
use freq, only: xnue0, lfreqint
!
implicit none
!
! ... local scalars
integer(i4b) :: i, j, k
integer(i4b) :: indx_1d
real(dp) :: dummy1, dummy2
!
! ... local functions
real(dp) :: bnue2
!
!----------------calculating snew directly for 3-d arrays---------------
!
!do i=1, ndzmax
!   write(*,*) scont3d(ndxmax/2+1,ndymax/2+1,i)
!enddo
!
do i=1, ndxmax
   do j=1, ndymax
      do k=1, ndzmax
         select case(imask3d(i,j,k))
            case(1,2,3)
            dummy2 = one - eps_cont3d(i,j,k)
            dummy1 = one - (one - eps_cont3d(i,j,k))*alocont_nn3d(i,j,k,14)
            scont3d(i,j,k) = (dummy2/dummy1) * mint3d(i,j,k) - &
                             (dummy2/dummy1) * alocont_nn3d(i,j,k,14) * scont3d(i,j,k) + &
                             (eps_cont3d(i,j,k)/dummy1)*bnue2(xnue0, t3d(i,j,k),lfreqint)
            case default
         endselect
      enddo
   enddo
enddo
!


!do i=1, ndzmax
!   write(*,*) scont3d(ndxmax/2+1,ndymax/2+1,i)
!enddo
!stop
!
!
end subroutine scont_new3d_diag

!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine scont_new3d_dn
!
!
!-----------------------------------------------------------------------
!------------calculates new iterate of source function------------------
!---approximate lambda iteration using nearest neighbour contribution---
!----------alo has to be stored in sparse-matrix formats----------------
!-----------------------------------------------------------------------
!
use prog_type
use fund_const
use dime3d, only: ndxmax, ndymax, ndzmax, scont3d, mint3d, imask3d, t3d, eps_cont3d, &
                  alocont_nn3d, alocont_rowindx, alocont_colindx, alocont_data, alocont_data_diag
use freq, only: xnue0, lfreqint
use timing, only: ts_alo, te_alo, ttot_alo
!
implicit none
!
! ... local scalars
integer(i4b) :: i, j, k
integer(i4b) :: indx_1d, indx_x, indx_y, indx_z, err, indx_rspec
real(dp) :: dummy1, dummy2, rspec
!
! ... local arrays
real(dp), dimension(:), allocatable :: mint_vec, scont_vec, bnue_vec, dummy_vec, sum_vec, epsc_vec
logical, dimension(:), allocatable :: mask_vec
!
! ... local functions
real(dp) :: bnue2
!
!-----------------------------------------------------------------------
!---------------------allocation of local arrays------------------------
!-----------------------------------------------------------------------
!
allocate(bnue_vec(ndxmax*ndymax*ndzmax), stat=err)
   if(err.ne.0) stop 'allocation error scont_new3d_dn: bnue_vec'
!
allocate(scont_vec(ndxmax*ndymax*ndzmax), stat=err)
   if(err.ne.0) stop 'allocation error scont_new3d_dn: scont_vec'
!
allocate(mint_vec(ndxmax*ndymax*ndzmax), stat=err)
   if(err.ne.0) stop 'allocation error scont_new3d_dn: mint_vec'
!
allocate(epsc_vec(ndxmax*ndymax*ndzmax), stat=err)
   if(err.ne.0) stop 'allocation error scont_new3d_dn: epsc_vec'
!
allocate(dummy_vec(ndxmax*ndymax*ndzmax), stat=err)
   if(err.ne.0) stop 'allocation error scont_new3d_dn: dummy_vec'
!
allocate(mask_vec(ndxmax*ndymax*ndzmax), stat=err)
   if(err.ne.0) stop 'allocation error scont_new3d_dn: mask_vec'
!
allocate(sum_vec(ndxmax*ndymax*ndzmax), stat=err)
   if(err.ne.0) stop 'allocation error scont_new3d_dn: sum_vec'
!
!----------------transform 3d-arrays to 1-d array-----------------------
!
mint_vec=0.d0
scont_vec=0.d0
bnue_vec=0.d0
epsc_vec=0.d0
dummy_vec=0.d0
mask_vec=.false.
!
do i=1, ndxmax
   do j=1, ndymax
      do k=1, ndzmax
         select case(imask3d(i,j,k))
            case(1,2,3)
               call conv_indx_3d_to_1d(i,j,k,ndxmax,ndymax,indx_1d)
               mint_vec(indx_1d)=mint3d(i,j,k)
               scont_vec(indx_1d)=scont3d(i,j,k)
               epsc_vec(indx_1d)=eps_cont3d(i,j,k)               
               mask_vec(indx_1d)=.true.
               bnue_vec(indx_1d)=bnue2(xnue0, t3d(i,j,k),lfreqint)
            case default
         endselect
      enddo
   enddo
enddo
!
call calc_alocont_dn_coo
!
!------------set up linear system that needs to be solved---------------
!-------------------in order to obtain s_new----------------------------
!
call matmul_coo(alocont_data, alocont_colindx, alocont_rowindx, scont_vec, dummy_vec, ndxmax*ndymax*ndzmax, 7*ndxmax*ndymax*ndzmax)
!
scont_vec=(one-epsc_vec)*dummy_vec
!
bnue_vec=epsc_vec*bnue_vec
!
mint_vec=(one-epsc_vec)*mint_vec
!
dummy_vec=mint_vec-scont_vec+bnue_vec
!
!alocont_data=(one-epsc_vec)*alocont_data
alocont_data_diag = one - (one - epsc_vec)*alocont_data_diag
!
do i=1, 7*ndxmax*ndymax*ndzmax
   if(alocont_colindx(i).eq.alocont_rowindx(i)) then
      alocont_data(i) = one - (one-epsc_vec(alocont_rowindx(i)))*alocont_data(i)
   else
      alocont_data(i) = -(one-epsc_vec(alocont_rowindx(i)))*alocont_data(i)
   endif
enddo
!
!--------------check if matrix is diagonal dominant---------------------
!-----------------and estimate spectral radius--------------------------
!
sum_vec=0.d0
!
do i=1, 7*ndxmax*ndymax*ndzmax
   if(alocont_rowindx(i).ne.alocont_colindx(i)) then
      sum_vec(alocont_rowindx(i)) = sum_vec(alocont_rowindx(i)) + abs(alocont_data(i))
   endif
enddo
!
rspec=maxval(sum_vec/alocont_data_diag)
indx_rspec=maxloc(sum_vec/alocont_data_diag, 1)
call conv_indx_1d_to_3d(indx_rspec, ndxmax, ndymax, indx_x, indx_y, indx_z)
write(*,'(a30,es20.8)') 'estimate of spectral radius', rspec
write(*,'(a30,i10,3i5)') 'at 1d, 3d indices', indx_rspec, indx_x, indx_y, indx_z
write(*,*)
!
!--------------------linear system now reads:---------------------------
!              alocont_matrix * s_new = scont_vec
!
write(*,*) '-----------------------inverting nearest neighbour alo-------------------------'
call cpu_time(ts_alo)
!
dummy_vec=-dummy_vec
!
call output_alocont_coo(alocont_data, alocont_colindx, alocont_rowindx, alocont_data_diag, dummy_vec, ndxmax*ndymax*ndzmax, 7*ndxmax*ndymax*ndzmax)
!
call jsor_coo(alocont_data, alocont_colindx, alocont_rowindx, alocont_data_diag, dummy_vec, 1.d0, ndxmax*ndymax*ndzmax, 7*ndxmax*ndymax*ndzmax, .false., scont_vec)
!
call cpu_time(te_alo)
ttot_alo=ttot_alo + te_alo-ts_alo
!
!write(*,*) 'total time needed to invert alo: ', te_alo-ts_alo
!
!---------back-transformation of source function on 3-d grid------------
!
do i=1, ndxmax*ndymax*ndzmax
   if(mask_vec(i)) then
      call conv_indx_1d_to_3d(i, ndxmax, ndymax, indx_x, indx_y, indx_z)
      scont3d(indx_x, indx_y, indx_z) = scont_vec(i)
   endif
enddo
!
write(*,*) '-------------------------------------------------------------------------------'
write(*,*)
!
!-----------------------------------------------------------------------
!
end subroutine scont_new3d_dn
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine scont_new3d_nn
!
!
!-----------------------------------------------------------------------
!------------calculates new iterate of source function------------------
!---approximate lambda iteration using nearest neighbour contribution---
!----------alo has to be stored in sparse-matrix formats----------------
!-----------------------------------------------------------------------
!
use prog_type
use fund_const
use dime3d, only: ndxmax, ndymax, ndzmax, scont3d, mint3d, imask3d, t3d, eps_cont3d, x, y, z, &
                  alocont_nn3d, alocont_rowindx, alocont_colindx, alocont_data, alocont_data_diag
use freq, only: xnue0, lfreqint
use timing, only: ts_alo, te_alo, ttot_alo
!
implicit none
!
! ... local scalars
integer(i4b) :: i, j, k
integer(i4b) :: indx_1d, indx_x, indx_y, indx_z, err, indx_rspec
real(dp) :: dummy1, dummy2, rspec
!
! ... local arrays
real(dp), dimension(:), allocatable :: mint_vec, scont_vec, bnue_vec, dummy_vec, sum_vec, epsc_vec
logical, dimension(:), allocatable :: mask_vec
!
! ... local functions
real(dp) :: bnue2
!
!-----------------------------------------------------------------------
!---------------------allocation of local arrays------------------------
!-----------------------------------------------------------------------
!
allocate(bnue_vec(ndxmax*ndymax*ndzmax), stat=err)
   if(err.ne.0) stop 'allocation error scont_new3d_nn: bnue_vec'
!
allocate(scont_vec(ndxmax*ndymax*ndzmax), stat=err)
   if(err.ne.0) stop 'allocation error scont_new3d_nn: scont_vec'
!
allocate(mint_vec(ndxmax*ndymax*ndzmax), stat=err)
   if(err.ne.0) stop 'allocation error scont_new3d_nn: mint_vec'
!
allocate(epsc_vec(ndxmax*ndymax*ndzmax), stat=err)
   if(err.ne.0) stop 'allocation error scont_new3d_nn: epsc_vec'
!
allocate(dummy_vec(ndxmax*ndymax*ndzmax), stat=err)
   if(err.ne.0) stop 'allocation error scont_new3d_nn: dummy_vec'
!
allocate(mask_vec(ndxmax*ndymax*ndzmax), stat=err)
   if(err.ne.0) stop 'allocation error scont_new3d_nn: mask_vec'
!
allocate(sum_vec(ndxmax*ndymax*ndzmax), stat=err)
   if(err.ne.0) stop 'allocation error  scont_new3d_nn: sum_vec'
!
!----------------transform 3d-arrays to 1-d array-----------------------
!
mint_vec=0.d0
scont_vec=0.d0
bnue_vec=0.d0
dummy_vec=0.d0
epsc_vec = 1.d0
mask_vec=.false.
!
do i=1, ndxmax
   do j=1, ndymax
      do k=1, ndzmax
         select case(imask3d(i,j,k))
            case(1,2,3)
               call conv_indx_3d_to_1d(i,j,k,ndxmax,ndymax,indx_1d)
               mint_vec(indx_1d)=mint3d(i,j,k)
               scont_vec(indx_1d)=scont3d(i,j,k)
               epsc_vec(indx_1d)=eps_cont3d(i,j,k)
               mask_vec(indx_1d)=.true.
               bnue_vec(indx_1d)=bnue2(xnue0, t3d(i,j,k),lfreqint)
            case default
         endselect
      enddo
   enddo
enddo
!
call calc_alocont_nn_coo
!
!------------set up linear system that needs to be solved---------------
!-------------------in order to obtain s_new----------------------------
!
call matmul_coo(alocont_data, alocont_colindx, alocont_rowindx, scont_vec, dummy_vec, ndxmax*ndymax*ndzmax, 27*ndxmax*ndymax*ndzmax)
!
scont_vec = (one-epsc_vec)*dummy_vec
!
bnue_vec = epsc_vec*bnue_vec
!
mint_vec = (one-epsc_vec)*mint_vec
!
dummy_vec=mint_vec-scont_vec+bnue_vec
!
!alocont_data=(1.d0-eps_cont)*alocont_data
alocont_data_diag = one-(one-epsc_vec)*alocont_data_diag
!
do i=1, 27*ndxmax*ndymax*ndzmax
   if(alocont_colindx(i).eq.alocont_rowindx(i)) then
      alocont_data(i) = one - (one-epsc_vec(alocont_rowindx(i)))*alocont_data(i)
   else
      alocont_data(i) = -(one-epsc_vec(alocont_rowindx(i)))*alocont_data(i)
   endif
!   write(*,*) alocont_data(i)
enddo
!stop
!
!--------------check if matrix is diagonal dominant---------------------
!-----------------and estimate spectral radius--------------------------
!
sum_vec=0.d0
!
do i=1, 27*ndxmax*ndymax*ndzmax
   if(alocont_rowindx(i).ne.alocont_colindx(i)) then
      sum_vec(alocont_rowindx(i)) = sum_vec(alocont_rowindx(i)) + abs(alocont_data(i))
   endif
enddo
!
rspec=maxval(sum_vec/alocont_data_diag)
indx_rspec=maxloc(sum_vec/alocont_data_diag, 1)
call conv_indx_1d_to_3d(indx_rspec, ndxmax, ndymax, indx_x, indx_y, indx_z)
write(*,'(a30,es20.8)') 'estimate of spectral radius', rspec
write(*,'(a30,i10,3i5)') 'at 1d, 3d indices', indx_rspec, indx_x, indx_y, indx_z
write(*,*)
!
!--------------------linear system now reads:---------------------------
!              alocont_matrix * s_new = scont_vec
!
write(*,*) '-----------------------inverting nearest neighbour alo-------------------------'
call cpu_time(ts_alo)
!
dummy_vec=-dummy_vec
!
call output_alocont_coo(alocont_data, alocont_colindx, alocont_rowindx, alocont_data_diag, dummy_vec, ndxmax*ndymax*ndzmax, 27*ndxmax*ndymax*ndzmax)
!
call jsor_coo(alocont_data, alocont_colindx, alocont_rowindx, alocont_data_diag, dummy_vec, 1.d0, ndxmax*ndymax*ndzmax, 27*ndxmax*ndymax*ndzmax, .false., scont_vec)
!
call cpu_time(te_alo)
ttot_alo=ttot_alo + te_alo-ts_alo
!
!write(*,*) 'total time needed to invert alo: ', te_alo-ts_alo
!
!---------back-transformation of source function on 3-d grid------------
!
do i=1, ndxmax*ndymax*ndzmax
   if(mask_vec(i)) then
      call conv_indx_1d_to_3d(i, ndxmax, ndymax, indx_x, indx_y, indx_z)
      scont3d(indx_x, indx_y, indx_z) = scont_vec(i)
   endif
enddo
!
write(*,*) '-------------------------------------------------------------------------------'
write(*,*)
!
!-----------------------------------------------------------------------
!
end subroutine scont_new3d_nn
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine calc_alocont_dn_coo
!
!-----------------------------------------------------------------------
!------------calculates approximate lambda operator---------------------
!-------using only direct neighbours (6 neighbours + local point)-------
!----------------------coo storage format-------------------------------
!-----------------------------------------------------------------------
!
use prog_type
use fund_const
use dime3d, only: ndxmax, ndymax, ndzmax, alocont_nn3d, &
                  alocont_data, alocont_data_diag, alocont_colindx, alocont_rowindx
!
implicit none
!
! ... local scalars
integer(i4b) :: i, j, k, err
integer(i4b) :: iindx, indx_1d, indx_1d_col, indx_x, indx_y, indx_z
integer :: ndiags_ijk, ndiags_max, ndiags_min
real(dp) :: rspec, rspec_max
!
! ... local arrays
!
! ... local functions
!
! ... local characters
character(len=50) :: enter
!
!-------------------allocation of alo-matrix----------------------------
!
if(.not.allocated(alocont_data)) then
   allocate(alocont_data(7*ndxmax*ndymax*ndzmax), stat=err)
      if(err.ne.0) stop 'allocation error calc_alocont_dn_coo: alocont_data'
endif

if(.not.allocated(alocont_data_diag)) then
   allocate(alocont_data_diag(ndxmax*ndymax*ndzmax), stat=err)
      if(err.ne.0) stop 'allocation error calc_alocont_dn_coo: alocont_data_diag'
endif
!
if(.not.allocated(alocont_colindx)) then
   allocate(alocont_colindx(7*ndxmax*ndymax*ndzmax), stat=err)
      if(err.ne.0) stop 'allocation error calc_alocont_dn_coo: alocont_col_indx'
endif
!
if(.not.allocated(alocont_rowindx)) then
   allocate(alocont_rowindx(7*ndxmax*ndymax*ndzmax), stat=err)
      if(err.ne.0) stop 'allocation error calc_alocont_dn_coo: alocont_row_indx'
endif
!
!-----------------calculate alo in coo storage format-------------------
!
!define maximum allowed spectral radius
rspec_max=1.d0
!
alocont_data=0.d0
alocont_data_diag=0.d0
alocont_rowindx=1
alocont_colindx=1
!
iindx=1
ndiags_max = 0
ndiags_min = 7
!
do i=2, ndxmax-1
   do j=2, ndymax-1
      do k=2, ndzmax-1
!
         rspec=0.d0
!
!include diagonal part
         rspec=rspec+abs(alocont_nn3d(i,j,k,14))
         call conv_indx_3d_to_1d(i,j,k, ndxmax, ndymax, indx_1d_col)
         alocont_data(iindx) = alocont_nn3d(i,j,k,14)
         alocont_colindx(iindx)=indx_1d_col
         alocont_rowindx(iindx)=indx_1d_col
         alocont_data_diag(indx_1d_col) = alocont_nn3d(i,j,k,14)
         ndiags_ijk=1
!
!include direct neighbour in positive x- direction
         rspec=rspec+abs(alocont_nn3d(i+1,j,k,13))
         if(rspec.lt.rspec_max) then
            call conv_indx_3d_to_1d(i+1,j,k,ndxmax,ndymax,indx_1d)
            alocont_data(iindx+1) = alocont_nn3d(i+1,j,k,13)
            alocont_colindx(iindx+1)=indx_1d
            alocont_rowindx(iindx+1)=indx_1d_col
            ndiags_ijk=ndiags_ijk+1
            ndiags_max = max(ndiags_ijk,ndiags_max)
         else
            ndiags_min = min(ndiags_ijk,ndiags_min)
            iindx=iindx+7
            cycle
         endif
!
!include direct neighbour in negative x-direction
         rspec=rspec+abs(alocont_nn3d(i-1,j,k,15))
         if(rspec.lt.rspec_max) then
            call conv_indx_3d_to_1d(i-1,j,k,ndxmax,ndymax,indx_1d)
            alocont_data(iindx+2) = alocont_nn3d(i-1,j,k,15)
            alocont_colindx(iindx+2)=indx_1d
            alocont_rowindx(iindx+2)=indx_1d_col
            ndiags_ijk=ndiags_ijk+1
            ndiags_max = max(ndiags_ijk,ndiags_max)
         else
            ndiags_min = min(ndiags_ijk,ndiags_min)
            iindx=iindx+7
            cycle
         endif
!
!include direct neighbour in positive y-direction
         rspec=rspec+abs(alocont_nn3d(i,j+1,k,11))
         if(rspec.lt.rspec_max) then
            call conv_indx_3d_to_1d(i,j+1,k,ndxmax,ndymax,indx_1d)
            alocont_data(iindx+3) = alocont_nn3d(i,j+1,k,11)
            alocont_colindx(iindx+3)=indx_1d
            alocont_rowindx(iindx+3)=indx_1d_col
            ndiags_ijk=ndiags_ijk+1
            ndiags_max = max(ndiags_ijk,ndiags_max)
         else
            ndiags_min = min(ndiags_ijk,ndiags_min)
            iindx=iindx+7
            cycle
         endif         
!
!include direct neighbour in negative y-direction
         rspec=rspec+abs(alocont_nn3d(i,j-1,k,17))
         if(rspec.lt.rspec_max) then
            call conv_indx_3d_to_1d(i,j-1,k,ndxmax,ndymax,indx_1d)
            alocont_data(iindx+4) = alocont_nn3d(i,j-1,k,17)
            alocont_colindx(iindx+4)=indx_1d
            alocont_rowindx(iindx+4)=indx_1d_col
            ndiags_ijk=ndiags_ijk+1
            ndiags_max = max(ndiags_ijk,ndiags_max)
         else
            ndiags_min = min(ndiags_ijk,ndiags_min)
            iindx=iindx+7
            cycle
         endif
!
!include direct neighbour in positive z-direction
         rspec=rspec+abs(alocont_nn3d(i,j,k+1,5))
         if(rspec.lt.rspec_max) then
            call conv_indx_3d_to_1d(i,j,k+1,ndxmax,ndymax,indx_1d)
            alocont_data(iindx+5) = alocont_nn3d(i,j,k+1,5)
            alocont_colindx(iindx+5)=indx_1d
            alocont_rowindx(iindx+5)=indx_1d_col
            ndiags_ijk=ndiags_ijk+1
            ndiags_max = max(ndiags_ijk,ndiags_max)
         else
            ndiags_min = min(ndiags_ijk,ndiags_min)
            iindx=iindx+7
            cycle
         endif
!
!include direct neighbour in negative z-direction
         rspec=rspec+abs(alocont_nn3d(i,j,k-1,23))
         if(rspec.lt.rspec_max) then
            call conv_indx_3d_to_1d(i,j,k-1,ndxmax,ndymax,indx_1d)
            alocont_data(iindx+6) = alocont_nn3d(i,j,k-1,23)
            alocont_colindx(iindx+6)=indx_1d
            alocont_rowindx(iindx+6)=indx_1d_col
            ndiags_ijk=ndiags_ijk+1
            ndiags_max = max(ndiags_ijk,ndiags_max)
         else
            ndiags_min = min(ndiags_ijk,ndiags_min)
            iindx=iindx+7
            cycle
         endif
!
!-----------------------------------------------------------------------
!
         iindx=iindx+7
!
      enddo
   enddo
enddo
!
write(*,*)
write(*,*) 'maximum number of used neighbours for ALO calculation', ndiags_max-1
write(*,*) 'minimum number of used neighbours for ALO calculation', ndiags_min-1
write(*,*)
return
!
!***********************************************************************
!*******************old version: per column*****************************
!***********************************************************************
!
!-----------------calculate alo in coo storage format-------------------
!
!alocont_data=0.d0
!alocont_rowindx=1
!alocont_colindx=1
!!
!iindx=1
!!
!do i=1, ndxmax
!   do j=1, ndymax
!      do k=1, ndzmax
!!
!!-------------------------diagonal part---------------------------------
!!
!         call conv_indx_3d_to_1d(i,j,k, ndxmax, ndymax, indx_1d_col)
!         alocont_data(iindx) = alocont_nn3d(i,j,k,14)
!         alocont_colindx(iindx)=indx_1d_col
!         alocont_rowindx(iindx)=indx_1d_col
!         alocont_data_diag(indx_1d_col) = alocont_nn3d(i,j,k,14)
!!
!!-----------nearest neighbour in positive x-direction-------------------
!!
!         if(i.eq.ndxmax) then
!!set row_indx to over-next neighbour
!!(for nn-alo, entries are zero anyhow, otherwise: spurious results)
!            call conv_indx_3d_to_1d(i-2,j,k,ndxmax,ndymax,indx_1d)
!            alocont_data(iindx+1) = 0.d0
!            alocont_colindx(iindx+1) = indx_1d_col
!            alocont_rowindx(iindx+1) = indx_1d
!         else
!            call conv_indx_3d_to_1d(i+1,j,k,ndxmax,ndymax,indx_1d)
!            alocont_data(iindx+1) = alocont_nn3d(i,j,k,15)
!            alocont_colindx(iindx+1)=indx_1d_col
!            alocont_rowindx(iindx+1)=indx_1d
!         endif
!!
!!-----------nearest neighbour in negative x-direction-------------------
!
!         if(i.eq.1) then
!!set row_indx to over-next neighbour
!!(for nn-alo, entries are zero anyhow, otherwise: spurious results)
!            call conv_indx_3d_to_1d(i+2,j,k,ndxmax,ndymax,indx_1d)
!            alocont_data(iindx+2) = 0.d0
!            alocont_colindx(iindx+2) = indx_1d_col
!            alocont_rowindx(iindx+2) = indx_1d
!         else
!            call conv_indx_3d_to_1d(i-1,j,k,ndxmax,ndymax,indx_1d)
!            alocont_data(iindx+2) = alocont_nn3d(i,j,k,13)
!            alocont_colindx(iindx+2)=indx_1d_col
!            alocont_rowindx(iindx+2)=indx_1d
!         endif
!!
!!-----------nearest neighbour in positive y-direction-------------------
!!
!         if(j.eq.ndymax) then
!!set row_indx to over-next neighbour
!!(for nn-alo, entries are zero anyhow, otherwise: spurious results)
!            call conv_indx_3d_to_1d(i,j-2,k,ndxmax,ndymax,indx_1d)
!            alocont_data(iindx+3) = 0.d0
!            alocont_colindx(iindx+3) = indx_1d_col
!            alocont_rowindx(iindx+3) = indx_1d
!         else
!            call conv_indx_3d_to_1d(i,j+1,k,ndxmax,ndymax,indx_1d)
!            alocont_data(iindx+3) = alocont_nn3d(i,j,k,17)
!            alocont_colindx(iindx+3)=indx_1d_col
!            alocont_rowindx(iindx+3)=indx_1d
!         endif
!!
!!-----------nearest neighbour in negative y-direction-------------------
!!
!         if(j.eq.1) then
!!set row_indx to over-next neighbour
!!(for nn-alo, entries are zero anyhow, otherwise: spurious results)
!            call conv_indx_3d_to_1d(i,j+2,k,ndxmax,ndymax,indx_1d)
!            alocont_data(iindx+4) = 0.d0
!            alocont_colindx(iindx+4) = indx_1d_col
!            alocont_rowindx(iindx+4) = indx_1d
!         else
!            call conv_indx_3d_to_1d(i,j-1,k,ndxmax,ndymax,indx_1d)
!            alocont_data(iindx+4) = alocont_nn3d(i,j,k,11)
!            alocont_colindx(iindx+4)=indx_1d_col
!            alocont_rowindx(iindx+4)=indx_1d
!         endif
!!
!!-----------nearest neighbour in positive z-direction-------------------
!!
!         if(k.eq.ndzmax) then
!!set row_indx to over-next neighbour
!!(for nn-alo, entries are zero anyhow, otherwise: spurious results)
!            call conv_indx_3d_to_1d(i,j,k-2,ndxmax,ndymax,indx_1d)
!            alocont_data(iindx+5) = 0.d0
!            alocont_colindx(iindx+5) = indx_1d_col
!            alocont_rowindx(iindx+5) = indx_1d
!         else
!            call conv_indx_3d_to_1d(i,j,k+1,ndxmax,ndymax,indx_1d)
!            alocont_data(iindx+5) = alocont_nn3d(i,j,k,23)
!            alocont_colindx(iindx+5)=indx_1d_col
!            alocont_rowindx(iindx+5)=indx_1d
!         endif
!!
!!-----------nearest neighbour in negative x-direction-------------------
!!
!         if(k.eq.1) then
!!set row_indx to over-next neighbour
!!(for nn-alo, entries are zero anyhow, otherwise: spurious results)
!            call conv_indx_3d_to_1d(i,j,k+2,ndxmax,ndymax,indx_1d)
!            alocont_data(iindx+6) = 0.d0
!            alocont_colindx(iindx+6) = indx_1d_col
!            alocont_rowindx(iindx+6) = indx_1d
!         else
!            call conv_indx_3d_to_1d(i,j,k-1,ndxmax,ndymax,indx_1d)
!            alocont_data(iindx+6) = alocont_nn3d(i,j,k,5)
!            alocont_colindx(iindx+6)=indx_1d_col
!            alocont_rowindx(iindx+6)=indx_1d
!         endif
!!
!         iindx=iindx+7
!      enddo
!   enddo
!enddo
!
!
!
end subroutine calc_alocont_dn_coo
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine calc_alocont_nn_coo
!
!-----------------------------------------------------------------------
!------------calculates approximate lambda operator---------------------
!-------using only nearest neighbours (26 neighbours + local point)-----
!----------------------coo storage format-------------------------------
!-----------------------------------------------------------------------
!
use prog_type
use fund_const
use dime3d, only: ndxmax, ndymax, ndzmax, alocont_nn3d, &
                  alocont_data, alocont_data_diag, alocont_colindx, alocont_rowindx, &
                  imask3d, x, y, z
!
implicit none
!
! ... local scalars
integer(i4b) :: i, j, k, err
integer(i4b) :: iindx, indx_1d, indx_1d_col, indx_x, indx_y, indx_z
integer :: ndiags_ijk, ndiags_max, ndiags_min
real(dp) :: sum, eps, rspec, rspec_max
!
! ... local arrays
!
! ... local functions
!
! ... local characters
character(len=50) :: enter
!
! ... for debug
real(dp), dimension(:), allocatable :: alocont_data2, alocont_data_diag2
integer(i4b), dimension(:), allocatable :: alocont_rowindx2, alocont_colindx2
!
!-------------------allocation of alo-matrix----------------------------
!
if(.not.allocated(alocont_data)) then
   allocate(alocont_data(27*ndxmax*ndymax*ndzmax), stat=err)
      if(err.ne.0) stop 'allocation error calc_alocont_nn_coo: alocont_data'
endif

if(.not.allocated(alocont_data_diag)) then
   allocate(alocont_data_diag(ndxmax*ndymax*ndzmax), stat=err)
      if(err.ne.0) stop 'allocation error calc_alocont_nn_coo: alocont_data_diag'
endif
!
if(.not.allocated(alocont_colindx)) then
   allocate(alocont_colindx(27*ndxmax*ndymax*ndzmax), stat=err)
      if(err.ne.0) stop 'allocation error calc_alocont_nn_coo: alocont_col_indx'
endif
!
if(.not.allocated(alocont_rowindx)) then
   allocate(alocont_rowindx(27*ndxmax*ndymax*ndzmax), stat=err)
      if(err.ne.0) stop 'allocation error calc_alocont_nn_coo: alocont_row_indx'
endif
!
!-----------------calculate alo in coo storage format-------------------
!
!define maximum allowed spectral radius
rspec_max=1.d0
!
alocont_data=0.d0
alocont_data_diag=0.d0
alocont_rowindx=1
alocont_colindx=1
!
iindx=1
ndiags_max = 0
ndiags_min = 27
!
do i=2, ndxmax-1
   do j=2, ndymax-1
      do k=2, ndzmax-1

!do i=17, 17
!   do j=31, 31
!      do k=31, 31

!        select case(imask3d(i,j,k))
!           case(0,4)
!              write(*,*) x(i), y(j), z(k)
!              write(*,'(9es20.8)') alocont_nn3d(i,j,k,1), alocont_nn3d(i,j,k,2), alocont_nn3d(i,j,k,3), &
!                                   alocont_nn3d(i,j,k,4), alocont_nn3d(i,j,k,5), alocont_nn3d(i,j,k,6), &
!                                   alocont_nn3d(i,j,k,7), alocont_nn3d(i,j,k,8), alocont_nn3d(i,j,k,9), &
!                                   alocont_nn3d(i,j,k,10), alocont_nn3d(i,j,k,11), alocont_nn3d(i,j,k,12), &
!                                   alocont_nn3d(i,j,k,13), alocont_nn3d(i,j,k,14), alocont_nn3d(i,j,k,15), &
!                                   alocont_nn3d(i,j,k,16), alocont_nn3d(i,j,k,17), alocont_nn3d(i,j,k,18), &
!                                   alocont_nn3d(i,j,k,19), alocont_nn3d(i,j,k,20), alocont_nn3d(i,j,k,21), &
!                                   alocont_nn3d(i,j,k,22), alocont_nn3d(i,j,k,23), alocont_nn3d(i,j,k,24), &
!                                   alocont_nn3d(i,j,k,25), alocont_nn3d(i,j,k,26), alocont_nn3d(i,j,k,27)
!              write(*,*)
!        end select

         rspec=0.d0
!
!
!include diagonal part
         rspec=rspec+abs(alocont_nn3d(i,j,k,14))
         call conv_indx_3d_to_1d(i,j,k, ndxmax, ndymax, indx_1d_col)
!!         write(*,*) 'including diagonal', indx_1d_col, rspec
         alocont_data(iindx) = alocont_nn3d(i,j,k,14)
         alocont_colindx(iindx)=indx_1d_col
         alocont_rowindx(iindx)=indx_1d_col
         alocont_data_diag(indx_1d_col) = alocont_nn3d(i,j,k,14)
         ndiags_ijk=1
!
!
!!include direct neighbour in positive x- direction
         rspec=rspec+abs(alocont_nn3d(i+1,j,k,13))
         if(rspec.lt.rspec_max) then
            call conv_indx_3d_to_1d(i+1,j,k,ndxmax,ndymax,indx_1d)
!            write(*,*) 'including x+ diagonal', indx_1d, rspec
            alocont_data(iindx+1) = alocont_nn3d(i+1,j,k,13)
            alocont_colindx(iindx+1)=indx_1d
            alocont_rowindx(iindx+1)=indx_1d_col
            ndiags_ijk=ndiags_ijk+1
            ndiags_max = max(ndiags_ijk,ndiags_max)
         else
            ndiags_min = min(ndiags_ijk,ndiags_min)
            iindx=iindx+27
            cycle
         endif
!
!!include direct neighbour in negative x-direction
         rspec=rspec+abs(alocont_nn3d(i-1,j,k,15))
         if(rspec.lt.rspec_max) then
            call conv_indx_3d_to_1d(i-1,j,k,ndxmax,ndymax,indx_1d)
!            write(*,*) 'including x- diagonal', indx_1d, rspec
            alocont_data(iindx+2) = alocont_nn3d(i-1,j,k,15)
            alocont_colindx(iindx+2)=indx_1d
            alocont_rowindx(iindx+2)=indx_1d_col
            ndiags_ijk=ndiags_ijk+1
            ndiags_max = max(ndiags_ijk,ndiags_max)
         else
            ndiags_min = min(ndiags_ijk,ndiags_min)
            iindx=iindx+27
            cycle
         endif
!
!!include direct neighbour in positive y-direction
         rspec=rspec+abs(alocont_nn3d(i,j+1,k,11))
         if(rspec.lt.rspec_max) then
            call conv_indx_3d_to_1d(i,j+1,k,ndxmax,ndymax,indx_1d)
!            write(*,*) 'including y+ diagonal', indx_1d, rspec
            alocont_data(iindx+3) = alocont_nn3d(i,j+1,k,11)
            alocont_colindx(iindx+3)=indx_1d
            alocont_rowindx(iindx+3)=indx_1d_col
            ndiags_ijk=ndiags_ijk+1
            ndiags_max = max(ndiags_ijk,ndiags_max)
         else
            ndiags_min = min(ndiags_ijk,ndiags_min)
            iindx=iindx+27
            cycle
         endif         
!
!!include direct neighbour in negative y-direction
         rspec=rspec+abs(alocont_nn3d(i,j-1,k,17))
         if(rspec.lt.rspec_max) then
            call conv_indx_3d_to_1d(i,j-1,k,ndxmax,ndymax,indx_1d)
!            write(*,*) 'including y- diagonal', indx_1d, rspec
            alocont_data(iindx+4) = alocont_nn3d(i,j-1,k,17)
            alocont_colindx(iindx+4)=indx_1d
            alocont_rowindx(iindx+4)=indx_1d_col
            ndiags_ijk=ndiags_ijk+1
            ndiags_max = max(ndiags_ijk,ndiags_max)
         else
            ndiags_min = min(ndiags_ijk,ndiags_min)
            iindx=iindx+27
            cycle
         endif
!
!include direct neighbour in positive z-direction
         rspec=rspec+abs(alocont_nn3d(i,j,k+1,5))
         if(rspec.lt.rspec_max) then
            call conv_indx_3d_to_1d(i,j,k+1,ndxmax,ndymax,indx_1d)
!            write(*,*) 'including z+ diagonal', indx_1d, rspec
            alocont_data(iindx+5) = alocont_nn3d(i,j,k+1,5)
            alocont_colindx(iindx+5)=indx_1d
            alocont_rowindx(iindx+5)=indx_1d_col
            ndiags_ijk=ndiags_ijk+1
            ndiags_max = max(ndiags_ijk,ndiags_max)
         else
            ndiags_min = min(ndiags_ijk,ndiags_min)
            iindx=iindx+27
            cycle
         endif
!
!include direct neighbour in negative z-direction
         rspec=rspec+abs(alocont_nn3d(i,j,k-1,23))
         if(rspec.lt.rspec_max) then
            call conv_indx_3d_to_1d(i,j,k-1,ndxmax,ndymax,indx_1d)
!            write(*,*) 'inculding z- diagonal', indx_1d, rspec
            alocont_data(iindx+6) = alocont_nn3d(i,j,k-1,23)
            alocont_colindx(iindx+6)=indx_1d
            alocont_rowindx(iindx+6)=indx_1d_col
            ndiags_ijk=ndiags_ijk+1
            ndiags_max = max(ndiags_ijk,ndiags_max)
         else
            ndiags_min = min(ndiags_ijk,ndiags_min)
            iindx=iindx+27
            cycle
         endif
!
!include neighbour in positive x-y-direction
         rspec=rspec+abs(alocont_nn3d(i+1,j+1,k,10))
         if(rspec.lt.rspec_max) then
            call conv_indx_3d_to_1d(i+1,j+1,k,ndxmax,ndymax,indx_1d)
!            write(*,*) 'including x+y+ diagonal', indx_1d, rspec
            alocont_data(iindx+7) = alocont_nn3d(i+1,j+1,k,10)
            alocont_colindx(iindx+7)=indx_1d
            alocont_rowindx(iindx+7)=indx_1d_col
            ndiags_ijk=ndiags_ijk+1
            ndiags_max = max(ndiags_ijk,ndiags_max)
         else
            ndiags_min = min(ndiags_ijk,ndiags_min)
            iindx=iindx+27
            cycle
         endif

!include neighbour in negative x- and negative y-direction
         rspec=rspec+abs(alocont_nn3d(i-1,j-1,k,18))
         if(rspec.lt.rspec_max) then
            call conv_indx_3d_to_1d(i-1,j-1,k,ndxmax,ndymax,indx_1d)
!            write(*,*) 'including x-y- diagonal', indx_1d, rspec
            alocont_data(iindx+8) = alocont_nn3d(i-1,j-1,k,18)
            alocont_colindx(iindx+8)=indx_1d
            alocont_rowindx(iindx+8)=indx_1d_col
            ndiags_ijk=ndiags_ijk+1
            ndiags_max = max(ndiags_ijk,ndiags_max)
         else
            ndiags_min = min(ndiags_ijk,ndiags_min)
            iindx=iindx+27
            cycle
         endif

!include neighbour in positive x- and negative y-direction
         rspec=rspec+abs(alocont_nn3d(i+1,j-1,k,16))
         if(rspec.lt.rspec_max) then
            call conv_indx_3d_to_1d(i+1,j-1,k,ndxmax,ndymax,indx_1d)
!            write(*,*) 'including x+y- diagonal', indx_1d, rspec
            alocont_data(iindx+9) = alocont_nn3d(i+1,j-1,k,16)
            alocont_colindx(iindx+9)=indx_1d
            alocont_rowindx(iindx+9)=indx_1d_col
            ndiags_ijk=ndiags_ijk+1
            ndiags_max = max(ndiags_ijk,ndiags_max)
         else
            ndiags_min = min(ndiags_ijk,ndiags_min)
            iindx=iindx+27
            cycle
         endif

!include neighbour in negative x- and positive y-direction
         rspec=rspec+abs(alocont_nn3d(i-1,j+1,k,12))
         if(rspec.lt.rspec_max) then
            call conv_indx_3d_to_1d(i-1,j+1,k,ndxmax,ndymax,indx_1d)
!            write(*,*) 'including x-y+ diagonal', indx_1d, rspec
            alocont_data(iindx+10) = alocont_nn3d(i-1,j+1,k,12)
            alocont_colindx(iindx+10)=indx_1d
            alocont_rowindx(iindx+10)=indx_1d_col
            ndiags_ijk=ndiags_ijk+1
            ndiags_max = max(ndiags_ijk,ndiags_max)
         else
            ndiags_min = min(ndiags_ijk,ndiags_min)
            iindx=iindx+27
            cycle
         endif

!including neighbour in positive x- and positive z-direction
         rspec=rspec+abs(alocont_nn3d(i+1,j,k+1,4))
         if(rspec.lt.rspec_max) then
            call conv_indx_3d_to_1d(i+1,j,k+1,ndxmax,ndymax,indx_1d)
!            write(*,*) 'including x+z+ diagonal', indx_1d, rspec
            alocont_data(iindx+11) = alocont_nn3d(i+1,j,k+1,4)
            alocont_colindx(iindx+11)=indx_1d
            alocont_rowindx(iindx+11)=indx_1d_col
            ndiags_ijk=ndiags_ijk+1
            ndiags_max = max(ndiags_ijk,ndiags_max)
         else
            ndiags_min = min(ndiags_ijk,ndiags_min)
            iindx=iindx+27
            cycle
         endif

!including neighbour in negative x- and negative z-direction
         rspec=rspec+abs(alocont_nn3d(i-1,j,k-1,24))
         if(rspec.lt.rspec_max) then
            call conv_indx_3d_to_1d(i-1,j,k-1,ndxmax,ndymax,indx_1d)
!            write(*,*) 'including x-z- diagonal', indx_1d, rspec
            alocont_data(iindx+12) = alocont_nn3d(i-1,j,k-1,24)
            alocont_colindx(iindx+12)=indx_1d
            alocont_rowindx(iindx+12)=indx_1d_col
            ndiags_ijk=ndiags_ijk+1
            ndiags_max = max(ndiags_ijk,ndiags_max)
         else
            ndiags_min = min(ndiags_ijk,ndiags_min)
            iindx=iindx+27
            cycle
         endif

!including neighbour in positive x- and negative z-direction
         rspec=rspec+abs(alocont_nn3d(i+1,j,k-1,22))
         if(rspec.lt.rspec_max) then
            call conv_indx_3d_to_1d(i+1,j,k-1,ndxmax,ndymax,indx_1d)
!            write(*,*) 'including x+z- diagonal', indx_1d, rspec
            alocont_data(iindx+13) = alocont_nn3d(i+1,j,k-1,22)
            alocont_colindx(iindx+13)=indx_1d
            alocont_rowindx(iindx+13)=indx_1d_col
            ndiags_ijk=ndiags_ijk+1
            ndiags_max = max(ndiags_ijk,ndiags_max)
         else
            ndiags_min = min(ndiags_ijk,ndiags_min)
            iindx=iindx+27
            cycle
         endif

!including neighbour in negative x- and positive z-direction
         rspec=rspec+abs(alocont_nn3d(i-1,j,k+1,6))
         if(rspec.lt.rspec_max) then
            call conv_indx_3d_to_1d(i-1,j,k+1,ndxmax,ndymax,indx_1d)
!            write(*,*) 'including x-z+ diagonal', indx_1d, rspec
            alocont_data(iindx+14) = alocont_nn3d(i-1,j,k+1,6)
            alocont_colindx(iindx+14)=indx_1d
            alocont_rowindx(iindx+14)=indx_1d_col
            ndiags_ijk=ndiags_ijk+1
            ndiags_max = max(ndiags_ijk,ndiags_max)
         else
            ndiags_min = min(ndiags_ijk,ndiags_min)
            iindx=iindx+27
            cycle
         endif

!include neighbour in positive y- and positive z-direction
         rspec=rspec+abs(alocont_nn3d(i,j+1,k+1,2))
         if(rspec.lt.rspec_max) then
         call conv_indx_3d_to_1d(i,j+1,k+1,ndxmax,ndymax,indx_1d)
!            write(*,*) 'including y+z+ diagonal', indx_1d, rspec
            alocont_data(iindx+15) = alocont_nn3d(i,j+1,k+1,2)
            alocont_colindx(iindx+15)=indx_1d
            alocont_rowindx(iindx+15)=indx_1d_col
            ndiags_ijk=ndiags_ijk+1
            ndiags_max = max(ndiags_ijk,ndiags_max)
         else
            ndiags_min = min(ndiags_ijk,ndiags_min)
            iindx=iindx+27
            cycle
         endif

!include neighbour in negative y- and negative z-direction
         rspec=rspec+abs(alocont_nn3d(i,j-1,k-1,26))
         if(rspec.lt.rspec_max) then
            call conv_indx_3d_to_1d(i,j-1,k-1,ndxmax,ndymax,indx_1d)
!            write(*,*) 'including y-z- diagonal', indx_1d, rspec
            alocont_data(iindx+16) = alocont_nn3d(i,j-1,k-1,26)
            alocont_colindx(iindx+16)=indx_1d
            alocont_rowindx(iindx+16)=indx_1d_col
            ndiags_ijk=ndiags_ijk+1
            ndiags_max = max(ndiags_ijk,ndiags_max)
         else
            ndiags_min = min(ndiags_ijk,ndiags_min)
            iindx=iindx+27
            cycle
         endif
!
!include neighbour in positive y- and negative z-direction
         rspec=rspec+abs(alocont_nn3d(i,j+1,k-1,20))
         if(rspec.lt.rspec_max) then
            call conv_indx_3d_to_1d(i,j+1,k-1,ndxmax,ndymax,indx_1d)
!            write(*,*) 'including y+z-diagonal', indx_1d, rspec
            alocont_data(iindx+17) = alocont_nn3d(i,j+1,k-1,20)
            alocont_colindx(iindx+17)=indx_1d
            alocont_rowindx(iindx+17)=indx_1d_col
            ndiags_ijk=ndiags_ijk+1
            ndiags_max = max(ndiags_ijk,ndiags_max)
         else
            ndiags_min = min(ndiags_ijk,ndiags_min)
            iindx=iindx+27
            cycle
         endif

!include neighbour in negative y- and positive z-direction
         rspec=rspec+abs(alocont_nn3d(i,j-1,k+1,8))
         if(rspec.lt.rspec_max) then
            call conv_indx_3d_to_1d(i,j-1,k+1,ndxmax,ndymax,indx_1d)
!            write(*,*) 'including y-z+ diagonal', indx_1d, rspec
            alocont_data(iindx+18) = alocont_nn3d(i,j-1,k+1,8)
            alocont_colindx(iindx+18)=indx_1d
            alocont_rowindx(iindx+18)=indx_1d_col
            ndiags_ijk=ndiags_ijk+1
            ndiags_max = max(ndiags_ijk,ndiags_max)
         else
            ndiags_min = min(ndiags_ijk,ndiags_min)
            iindx=iindx+27
            cycle
         endif

!include neighbour in positive x-, positive y- and positive z-direction
         rspec=rspec+abs(alocont_nn3d(i+1,j+1,k+1,1))
         if(rspec.lt.rspec_max) then
            call conv_indx_3d_to_1d(i+1,j+1,k+1,ndxmax,ndymax,indx_1d)
!            write(*,*) 'including x+y+z+ diagonal', indx_1d, rspec
            alocont_data(iindx+19) = alocont_nn3d(i+1,j+1,k+1,1)
            alocont_colindx(iindx+19)=indx_1d
            alocont_rowindx(iindx+19)=indx_1d_col
            ndiags_ijk=ndiags_ijk+1
            ndiags_max = max(ndiags_ijk,ndiags_max)
         else
            ndiags_min = min(ndiags_ijk,ndiags_min)
            iindx=iindx+27
            cycle
         endif

!include neighbour in negative x-, negative y- and negative z-direction
         rspec=rspec+abs(alocont_nn3d(i-1,j-1,k-1,27))
         if(rspec.lt.rspec_max) then
            call conv_indx_3d_to_1d(i-1,j-1,k-1,ndxmax,ndymax,indx_1d)
!            write(*,*) 'including x-y-z- diagonal', indx_1d, rspec
            alocont_data(iindx+20) = alocont_nn3d(i-1,j-1,k-1,27)
            alocont_colindx(iindx+20)=indx_1d
            alocont_rowindx(iindx+20)=indx_1d_col
            ndiags_ijk=ndiags_ijk+1
            ndiags_max = max(ndiags_ijk,ndiags_max)
         else
            ndiags_min = min(ndiags_ijk,ndiags_min)
            iindx=iindx+27
            cycle
         endif
!
!include neighbour in positive x-, negative y- and positive z-direction
         rspec=rspec+abs(alocont_nn3d(i+1,j-1,k+1,7))
         if(rspec.lt.rspec_max) then
            call conv_indx_3d_to_1d(i+1,j-1,k+1,ndxmax,ndymax,indx_1d)
!            write(*,*) 'including x-y-z+ diagonal', indx_1d, rspec
            alocont_data(iindx+21) = alocont_nn3d(i+1,j-1,k+1,7)
            alocont_colindx(iindx+21)=indx_1d
            alocont_rowindx(iindx+21)=indx_1d_col
            ndiags_ijk=ndiags_ijk+1
            ndiags_max = max(ndiags_ijk,ndiags_max)
         else
            ndiags_min = min(ndiags_ijk,ndiags_min)
            iindx=iindx+27
            cycle
         endif
!
!include neighbour in negative x-, positive y- and negative z-direction
         rspec=rspec+abs(alocont_nn3d(i-1,j+1,k-1,21))
         if(rspec.lt.rspec_max) then
            call conv_indx_3d_to_1d(i-1,j+1,k-1,ndxmax,ndymax,indx_1d)
!            write(*,*) 'including x-y+z- diagonal', indx_1d, rspec
            alocont_data(iindx+22) = alocont_nn3d(i-1,j+1,k-1,21)
            alocont_colindx(iindx+22)=indx_1d
            alocont_rowindx(iindx+22)=indx_1d_col
            ndiags_ijk=ndiags_ijk+1
            ndiags_max = max(ndiags_ijk,ndiags_max)
         else
            ndiags_min = min(ndiags_ijk,ndiags_min)
            iindx=iindx+27
            cycle
         endif
!
!include neighbour in positive x-, positive y- and negative z-direction
         rspec=rspec+abs(alocont_nn3d(i+1,j+1,k-1,19))
         if(rspec.lt.rspec_max) then
            call conv_indx_3d_to_1d(i+1,j+1,k-1,ndxmax,ndymax,indx_1d)
!            write(*,*) 'including x+y+z- diagonal', indx_1d, rspec
            alocont_data(iindx+23) = alocont_nn3d(i+1,j+1,k-1,19)
            alocont_colindx(iindx+23)=indx_1d
            alocont_rowindx(iindx+23)=indx_1d_col
            ndiags_ijk=ndiags_ijk+1
            ndiags_max = max(ndiags_ijk,ndiags_max)
         else
            ndiags_min = min(ndiags_ijk,ndiags_min)
            iindx=iindx+27
            cycle
         endif
!
!include neighbour in negative x-, negative y- and positive z-direction
         rspec=rspec+abs(alocont_nn3d(i-1,j-1,k+1,9))
         if(rspec.lt.rspec_max) then
            call conv_indx_3d_to_1d(i-1,j-1,k+1,ndxmax,ndymax,indx_1d)
!            write(*,*) 'including x-y-z+ diagonal', indx_1d, rspec
            alocont_data(iindx+24) = alocont_nn3d(i-1,j-1,k+1,9)
            alocont_colindx(iindx+24)=indx_1d
            alocont_rowindx(iindx+24)=indx_1d_col
            ndiags_ijk=ndiags_ijk+1
            ndiags_max = max(ndiags_ijk,ndiags_max)
         else
            ndiags_min = min(ndiags_ijk,ndiags_min)
            iindx=iindx+27
            cycle
         endif

!include neighbour in positive x-, negative y- and negative z-direction
         rspec=rspec+abs(alocont_nn3d(i+1,j-1,k-1,25))
         if(rspec.lt.rspec_max) then
            call conv_indx_3d_to_1d(i+1,j-1,k-1,ndxmax,ndymax,indx_1d)
!            write(*,*) 'including x+y-z- diagonal', indx_1d, rspec
            alocont_data(iindx+25) = alocont_nn3d(i+1,j-1,k-1,25)
            alocont_colindx(iindx+25)=indx_1d
            alocont_rowindx(iindx+25)=indx_1d_col
            ndiags_ijk=ndiags_ijk+1
            ndiags_max = max(ndiags_ijk,ndiags_max)
         else
            ndiags_min = min(ndiags_ijk,ndiags_min)
            iindx=iindx+27
            cycle
         endif

!include neighbour in negative x-, positive y- and positive z-direction
         rspec=rspec+abs(alocont_nn3d(i-1,j+1,k+1,3))
         if(rspec.lt.rspec_max) then
            call conv_indx_3d_to_1d(i-1,j+1,k+1,ndxmax,ndymax,indx_1d)
!            write(*,*) 'including x-y+z+ diagonal', indx_1d, rspec
            alocont_data(iindx+26) = alocont_nn3d(i-1,j+1,k+1,3)
            alocont_colindx(iindx+26)=indx_1d
            alocont_rowindx(iindx+26)=indx_1d_col
            ndiags_ijk=ndiags_ijk+1
            ndiags_max = max(ndiags_ijk,ndiags_max)
!            write(*,*) ndiags_ijk
         else
            ndiags_min = min(ndiags_ijk,ndiags_min)
            iindx=iindx+27
            cycle
         endif
!
!-----------------------------------------------------------------------
!
         iindx=iindx+27
!
      enddo
   enddo
enddo
!
write(*,*)
write(*,*) 'maximum number of used neighbours for ALO calculation', ndiags_max-1
write(*,*) 'minimum number of used neighbours for ALO calculation', ndiags_min-1
write(*,*)
return
!
!***********************************************************************
!***********************old version: per column*************************
!***********************************************************************
!
!-----------------calculate alo in coo storage format-------------------
!
!alocont_data=0.d0
!alocont_data_diag=0.d0
!alocont_rowindx=1
!alocont_colindx=1
!
!iindx=1
!
!do i=2, ndxmax-1
!   do j=2, ndymax-1
!      do k=2, ndzmax-1
!
!!-----------------------neighbours on level k-1-------------------------
!!
!         call conv_indx_3d_to_1d(i,j,k,ndxmax,ndymax,indx_1d_col)
!         alocont_data(iindx) = alocont_nn3d(i,j,k,14)
!         alocont_colindx(iindx)=indx_1d_col
!         alocont_rowindx(iindx)=indx_1d_col
!         alocont_data_diag(indx_1d_col) = alocont_nn3d(i,j,k,14)
!!
!         call conv_indx_3d_to_1d(i-1,j-1,k-1,ndxmax,ndymax,indx_1d)
!         alocont_data(iindx+1) = alocont_nn3d(i,j,k,1)
!         alocont_colindx(iindx+1)=indx_1d_col
!         alocont_rowindx(iindx+1)=indx_1d
!!
!         call conv_indx_3d_to_1d(i,j-1,k-1,ndxmax,ndymax,indx_1d)
!         alocont_data(iindx+2) = alocont_nn3d(i,j,k,2)
!         alocont_colindx(iindx+2)=indx_1d_col
!         alocont_rowindx(iindx+2)=indx_1d
!!
!         call conv_indx_3d_to_1d(i+1,j-1,k-1,ndxmax,ndymax,indx_1d)
!         alocont_data(iindx+3) = alocont_nn3d(i,j,k,3)
!         alocont_colindx(iindx+3)=indx_1d_col
!         alocont_rowindx(iindx+3)=indx_1d
!!
!         call conv_indx_3d_to_1d(i-1,j,k-1,ndxmax,ndymax,indx_1d)
!         alocont_data(iindx+4) = alocont_nn3d(i,j,k,4)
!         alocont_colindx(iindx+4)=indx_1d_col
!         alocont_rowindx(iindx+4)=indx_1d
!!
!         call conv_indx_3d_to_1d(i,j,k-1,ndxmax,ndymax,indx_1d)
!         alocont_data(iindx+5) = alocont_nn3d(i,j,k,5)
!         alocont_colindx(iindx+5)=indx_1d_col
!         alocont_rowindx(iindx+5)=indx_1d
!!
!         call conv_indx_3d_to_1d(i+1,j,k-1,ndxmax,ndymax,indx_1d)
!         alocont_data(iindx+6) = alocont_nn3d(i,j,k,6)
!         alocont_colindx(iindx+6)=indx_1d_col
!         alocont_rowindx(iindx+6)=indx_1d
!!
!         call conv_indx_3d_to_1d(i-1,j+1,k-1,ndxmax,ndymax,indx_1d)
!         alocont_data(iindx+7) = alocont_nn3d(i,j,k,7)
!         alocont_colindx(iindx+7)=indx_1d_col
!         alocont_rowindx(iindx+7)=indx_1d
!!
!         call conv_indx_3d_to_1d(i,j+1,k-1,ndxmax,ndymax,indx_1d)
!         alocont_data(iindx+8) = alocont_nn3d(i,j,k,8)
!         alocont_colindx(iindx+8)=indx_1d_col
!         alocont_rowindx(iindx+8)=indx_1d
!!
!         call conv_indx_3d_to_1d(i+1,j+1,k-1,ndxmax,ndymax,indx_1d)
!         alocont_data(iindx+9) = alocont_nn3d(i,j,k,9)
!         alocont_colindx(iindx+9)=indx_1d_col
!         alocont_rowindx(iindx+9)=indx_1d
!!
!!------------------------neighbours on level k--------------------------
!!
!         call conv_indx_3d_to_1d(i-1,j-1,k,ndxmax,ndymax,indx_1d)
!         alocont_data(iindx+10) = alocont_nn3d(i,j,k,10)
!         alocont_colindx(iindx+10)=indx_1d_col
!         alocont_rowindx(iindx+10)=indx_1d
!!
!         call conv_indx_3d_to_1d(i,j-1,k,ndxmax,ndymax,indx_1d)
!         alocont_data(iindx+11) = alocont_nn3d(i,j,k,11)
!         alocont_colindx(iindx+11)=indx_1d_col
!         alocont_rowindx(iindx+11)=indx_1d
!!
!         call conv_indx_3d_to_1d(i+1,j-1,k,ndxmax,ndymax,indx_1d)
!         alocont_data(iindx+12) = alocont_nn3d(i,j,k,12)
!         alocont_colindx(iindx+12)=indx_1d_col
!         alocont_rowindx(iindx+12)=indx_1d
!!
!         call conv_indx_3d_to_1d(i-1,j,k,ndxmax,ndymax,indx_1d)
!         alocont_data(iindx+13) = alocont_nn3d(i,j,k,13)
!         alocont_colindx(iindx+13)=indx_1d_col
!         alocont_rowindx(iindx+13)=indx_1d
!!
!         call conv_indx_3d_to_1d(i+1,j,k,ndxmax,ndymax,indx_1d)
!         alocont_data(iindx+14) = alocont_nn3d(i,j,k,15)
!         alocont_colindx(iindx+14)=indx_1d_col
!         alocont_rowindx(iindx+14)=indx_1d
!!
!         call conv_indx_3d_to_1d(i-1,j+1,k,ndxmax,ndymax,indx_1d)
!         alocont_data(iindx+15) = alocont_nn3d(i,j,k,16)
!         alocont_colindx(iindx+15)=indx_1d_col
!         alocont_rowindx(iindx+15)=indx_1d
!!
!         call conv_indx_3d_to_1d(i,j+1,k,ndxmax,ndymax,indx_1d)
!         alocont_data(iindx+16) = alocont_nn3d(i,j,k,17)
!         alocont_colindx(iindx+16)=indx_1d_col
!         alocont_rowindx(iindx+16)=indx_1d
!!
!         call conv_indx_3d_to_1d(i+1,j+1,k,ndxmax,ndymax,indx_1d)
!         alocont_data(iindx+17) = alocont_nn3d(i,j,k,18)
!         alocont_colindx(iindx+17)=indx_1d_col
!         alocont_rowindx(iindx+17)=indx_1d
!
!!-----------------------neighbours on level k+1-------------------------
!
!         call conv_indx_3d_to_1d(i-1,j-1,k+1,ndxmax,ndymax,indx_1d)
!         alocont_data(iindx+18) = alocont_nn3d(i,j,k,19)
!         alocont_colindx(iindx+18)=indx_1d_col
!         alocont_rowindx(iindx+18)=indx_1d
!
!         call conv_indx_3d_to_1d(i,j-1,k+1,ndxmax,ndymax,indx_1d)
!         alocont_data(iindx+19) = alocont_nn3d(i,j,k,20)
!         alocont_colindx(iindx+19)=indx_1d_col
!         alocont_rowindx(iindx+19)=indx_1d
!
!         call conv_indx_3d_to_1d(i+1,j-1,k+1,ndxmax,ndymax,indx_1d)
!         alocont_data(iindx+20) = alocont_nn3d(i,j,k,21)
!         alocont_colindx(iindx+20)=indx_1d_col
!         alocont_rowindx(iindx+20)=indx_1d
!
!         call conv_indx_3d_to_1d(i-1,j,k+1,ndxmax,ndymax,indx_1d)
!         alocont_data(iindx+21) = alocont_nn3d(i,j,k,22)
!         alocont_colindx(iindx+21)=indx_1d_col
!         alocont_rowindx(iindx+21)=indx_1d
!
!         call conv_indx_3d_to_1d(i,j,k+1,ndxmax,ndymax,indx_1d)
!         alocont_data(iindx+22) = alocont_nn3d(i,j,k,23)
!         alocont_colindx(iindx+22)=indx_1d_col
!         alocont_rowindx(iindx+22)=indx_1d
!
!         call conv_indx_3d_to_1d(i+1,j,k+1,ndxmax,ndymax,indx_1d)
!         alocont_data(iindx+23) = alocont_nn3d(i,j,k,24)
!         alocont_colindx(iindx+23)=indx_1d_col
!         alocont_rowindx(iindx+23)=indx_1d
!
!         call conv_indx_3d_to_1d(i-1,j+1,k+1,ndxmax,ndymax,indx_1d)
!         alocont_data(iindx+24) = alocont_nn3d(i,j,k,25)
!         alocont_colindx(iindx+24)=indx_1d_col
!         alocont_rowindx(iindx+24)=indx_1d
!
!         call conv_indx_3d_to_1d(i,j+1,k+1,ndxmax,ndymax,indx_1d)
!         alocont_data(iindx+25) = alocont_nn3d(i,j,k,26)
!         alocont_colindx(iindx+25)=indx_1d_col
!         alocont_rowindx(iindx+25)=indx_1d
!
!         call conv_indx_3d_to_1d(i+1,j+1,k+1,ndxmax,ndymax,indx_1d)
!         alocont_data(iindx+26) = alocont_nn3d(i,j,k,27)
!         alocont_colindx(iindx+26)=indx_1d_col
!         alocont_rowindx(iindx+26)=indx_1d
!
!         iindx=iindx+27
!
!      enddo
!   enddo
!enddo
!
!
end subroutine calc_alocont_nn_coo
!

!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine check_alocont_o(oindx)
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
use prog_type
use fund_const
use dime3d, only: ndxmax, ndymax, ndzmax, alocont_o_nn3d, imask3d, x, y, z
use angles, only: n_x, n_y, n_z
!
implicit none
!
! ... arguments
integer(i4b), intent(in) :: oindx
!
! ... local scalars
integer(i4b) :: i, j, k
real(dp) :: sum
!
! ... local arrays
!
! ... local functions
!
! ... local characters
character(len=50) :: enter
!
return
!
do i=2, ndxmax-1
   do j=2, ndymax-1
      do k=2, ndzmax-1

         select case(imask3d(i,j,k))
            case(1,2,3)
               sum = alocont_o_nn3d(i-1,j-1,k-1,27) + alocont_o_nn3d(i,j-1,k-1,26) + alocont_o_nn3d(i+1,j-1,k-1,25) + &
                     alocont_o_nn3d(i-1,j,k-1,24) + alocont_o_nn3d(i,j,k-1,23) + alocont_o_nn3d(i+1,j,k-1,22) + &
                     alocont_o_nn3d(i-1,j+1,k-1,21) + alocont_o_nn3d(i,j+1,k-1,20) + alocont_o_nn3d(i+1,j+1,k-1,19) + &
                     alocont_o_nn3d(i-1,j-1,k,18) + alocont_o_nn3d(i,j-1,k,17) + alocont_o_nn3d(i+1,j-1,k,16) + &
                     alocont_o_nn3d(i-1,j,k,15) + alocont_o_nn3d(i+1,j,k,13) + alocont_o_nn3d(i,j,k,14) + &
                     alocont_o_nn3d(i-1,j+1,k,12) + alocont_o_nn3d(i,j+1,k,11) + alocont_o_nn3d(i+1,j+1,k,10) + &
                     alocont_o_nn3d(i-1,j-1,k+1,9) + alocont_o_nn3d(i,j-1,k+1,8) + alocont_o_nn3d(i+1,j-1,k+1,7) + &
                     alocont_o_nn3d(i-1,j,k+1,6) + alocont_o_nn3d(i,j,k+1,5) + alocont_o_nn3d(i+1,j,k+1,4) + &
                     alocont_o_nn3d(i-1,j+1,k+1,3) + alocont_o_nn3d(i,j+1,k+1,2) + alocont_o_nn3d(i+1,j+1,k+1,1)
               if(sum.gt.1.d0+1.d-8) then
                  write(*,*) alocont_o_nn3d(i-1,j-1,k-1,27) , alocont_o_nn3d(i,j-1,k-1,26) , alocont_o_nn3d(i+1,j-1,k-1,25) , &
                             alocont_o_nn3d(i-1,j,k-1,24) , alocont_o_nn3d(i,j,k-1,23) , alocont_o_nn3d(i+1,j,k-1,22) , &
                             alocont_o_nn3d(i-1,j+1,k-1,21) , alocont_o_nn3d(i,j+1,k-1,20) , alocont_o_nn3d(i+1,j+1,k-1,19)
                  write(*,*) 
                  write(*,*) alocont_o_nn3d(i-1,j-1,k,18) , alocont_o_nn3d(i,j-1,k,17) , alocont_o_nn3d(i+1,j-1,k,16), &
                             alocont_o_nn3d(i-1,j,k,15) , alocont_o_nn3d(i,j,k,14) , alocont_o_nn3d(i+1,j,k,13) , &
                             alocont_o_nn3d(i-1,j+1,k,12) , alocont_o_nn3d(i,j+1,k,11) , alocont_o_nn3d(i+1,j+1,k,10)
                  write(*,*)
                  write(*,*) alocont_o_nn3d(i-1,j-1,k+1,9) , alocont_o_nn3d(i,j-1,k+1,8) , alocont_o_nn3d(i+1,j-1,k+1,7) , &
                             alocont_o_nn3d(i-1,j,k+1,6) , alocont_o_nn3d(i,j,k+1,5) , alocont_o_nn3d(i+1,j,k+1,4) , &
                             alocont_o_nn3d(i-1,j+1,k+1,3) , alocont_o_nn3d(i,j+1,k+1,2) , alocont_o_nn3d(i+1,j+1,k+1,1)
                  write(*,*)
                  write(*,*) i, j, k
                  write(*,*) x(i), y(j), z(k)
                  write(*,*)
                  write(*,*) n_x(oindx), n_y(oindx), n_z(oindx)
                  write(*,*)
                  write(*,*) sum
                  stop 'error in check_alocont_o: spectral radius too large'
               endif
         end select

      enddo
   enddo
enddo
!
end subroutine check_alocont_o
