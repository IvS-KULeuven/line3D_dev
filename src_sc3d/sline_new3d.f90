subroutine sline_new3d
!
!-----------------------------------------------------------------------
!-----------calculates new iterate of line source function--------------
!   different options: classical lambda-iteration
!                      diagonal of lambda-matrix
!                      nearest neighbours
!-----------------------------------------------------------------------
!
use prog_type
use options, only: opt_alo_line
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
select case(opt_alo_line)
   case(0)
      call sline_new3d_classic
   case(1)
      call sline_new3d_diag
   case(2)
      call sline_new3d_dn
   case(3)
      call sline_new3d_nn
   case default
      stop 'set option opt_alo_line'
end select
!
!
end subroutine sline_new3d
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine sline_new3d_classic
!
!-----------------------------------------------------------------------
!-----------calculates new iterate of line source function--------------
!-------------------classical lambda-iteration--------------------------
!-----------------------------------------------------------------------
!
use prog_type
use fund_const
use dime3d, only: ndxmax, ndymax, ndzmax, sline3d, mintbar3d, imask3d, t3d
use params_input, only: eps_line
use freq, only: xnue0
!
implicit none
!
! ... local scalars
integer(i4b) :: i, j, k
!
! ... local functions
real(dp) :: bnue
!
!-----------------------------------------------------------------------
!
do i=1, ndxmax
   do j=1, ndymax
      do k=1, ndzmax
         select case(imask3d(i,j,k))
            case(1,2,3)
               sline3d(i,j,k) = (1.d0-eps_line) * mintbar3d(i,j,k) + eps_line*bnue(xnue0, t3d(i,j,k))
            case default
         endselect
      enddo
   enddo
enddo
!
end subroutine sline_new3d_classic
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine sline_new3d_diag
!
!-----------------------------------------------------------------------
!-----------calculates new iterate of line source function--------------
!---------approximate lambda iteration using only diagonal--------------
!-----------------------------------------------------------------------
!
use prog_type
use fund_const
use dime3d, only: ndxmax, ndymax, ndzmax, sline3d, ssobo3d, aloline_nn3d, mintbar3d, imask3d, t3d, x, y, z
use params_input, only: eps_line
use freq, only: xnue0
!
implicit none
!
! ... local scalars
integer(i4b) :: i, j, k
integer(i4b) :: indx_1d
real(dp) :: dummy1, dummy2
!
! ... local functions
real(dp) :: bnue
!
! ... local logicals
logical :: lnegative
!
!----------------calculating snew directly for 3-d arrays---------------
!
dummy2=1.d0-eps_line

!do i=1, ndzmax
!   write(*,*) sline3d(ndxmax/2+1,ndymax/2+1,i)
!enddo
!
do i=1, ndxmax
   do j=1, ndymax
      do k=1, ndzmax
         select case(imask3d(i,j,k))
            case(1,2,3)
            dummy1=1.d0-(1.d0-eps_line)*aloline_nn3d(i,j,k,14)
            sline3d(i,j,k) = (dummy2/dummy1) * mintbar3d(i,j,k) - &
                             (dummy2/dummy1) * aloline_nn3d(i,j,k,14) * sline3d(i,j,k) + &
                             (eps_line/dummy1)*bnue(xnue0, t3d(i,j,k))
            case default
         endselect
      enddo
   enddo
enddo
!
!
!
!check if negative source functions occurr
lnegative=.false.
do i=1, ndxmax
   do j=1, ndymax
      do k=1, ndzmax
         if(sline3d(i,j,k).lt.0.d0) then
            write(*,'(3i5,7es20.8,i5)') i, j, k, x(i), y(j), z(k), sline3d(i,j,k), ssobo3d(i,j,k), mintbar3d(i,j,k), aloline_nn3d(i,j,k,14), imask3d(i,j,k)
            lnegative=.true.
         endif
      enddo
   enddo
enddo
if(lnegative) then
   write(*,*) 'error in sline_new3d_diag: new iterate gives negative source functions'
   write(*,*) 'try: linear interpolations or lower band alo'
   stop
endif
!
write(*,*) '-------------------------------------------------------------------------------'
write(*,*)
!
!
end subroutine sline_new3d_diag

!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine sline_new3d_dn
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
use dime3d, only: ndxmax, ndymax, ndzmax, sline3d, ssobo3d, mintbar3d, imask3d, x, y, z, t3d, &
                  aloline_nn3d, aloline_rowindx, aloline_colindx, aloline_data, aloline_data_diag
use params_input, only: eps_line
use freq, only: xnue0
use timing, only: ts_alo, te_alo, ttot_alo
use mod_interp2d, only: wpa_interp2d, wpb_interp2d, wp_interp2d, wpa_interp1d, &
                        wpb_interp1d, wp_interp1d, wpa_integ1d, wpb_integ1d, &
                        wp_integ1d, lng_expol
!
implicit none
!
! ... local scalars
integer(i4b) :: i, j, k
integer(i4b) :: indx_1d, indx_x, indx_y, indx_z, err, indx_rspec
real(dp) :: dummy1, dummy2, rspec
!
! ... local arrays
real(dp), dimension(:), allocatable :: mintbar_vec, sline_vec, bnue_vec, dummy_vec, sum_vec
logical, dimension(:), allocatable :: mask_vec
!
! ... local logicals
logical :: lnegative
!
! ... local functions
real(dp) :: bnue
!
!-----------------------------------------------------------------------
!---------------------allocation of local arrays------------------------
!-----------------------------------------------------------------------
!
allocate(bnue_vec(ndxmax*ndymax*ndzmax), stat=err)
   if(err.ne.0) stop 'allocation error sline_new3d_dn: bnue_vec'
!
allocate(sline_vec(ndxmax*ndymax*ndzmax), stat=err)
   if(err.ne.0) stop 'allocation error sline_new3d_dn: sline_vec'
!
allocate(mintbar_vec(ndxmax*ndymax*ndzmax), stat=err)
   if(err.ne.0) stop 'allocation error sline_new3d_dn: mintbar_vec'
!
allocate(dummy_vec(ndxmax*ndymax*ndzmax), stat=err)
   if(err.ne.0) stop 'allocation error sline_new3d_dn: dummy_vec'
!
allocate(mask_vec(ndxmax*ndymax*ndzmax), stat=err)
   if(err.ne.0) stop 'allocation error sline_new3d_dn: mask_vec'
!
allocate(sum_vec(ndxmax*ndymax*ndzmax), stat=err)
   if(err.ne.0) stop 'allocation error sline_new3d_dn: sum_vec'
!
!----------------transform 3d-arrays to 1-d array-----------------------
!
mintbar_vec=0.d0
sline_vec=0.d0
bnue_vec=0.d0
dummy_vec=0.d0
mask_vec=.false.
!
do i=1, ndxmax
   do j=1, ndymax
      do k=1, ndzmax
         select case(imask3d(i,j,k))
            case(1,2,3)
               call conv_indx_3d_to_1d(i,j,k,ndxmax,ndymax,indx_1d)
               mintbar_vec(indx_1d)=mintbar3d(i,j,k)
               sline_vec(indx_1d)=sline3d(i,j,k)
               mask_vec(indx_1d)=.true.
               bnue_vec(indx_1d)=bnue(xnue0, t3d(i,j,k))
            case default
         endselect
      enddo
   enddo
enddo
!
call calc_aloline_dn_coo
!
!------------set up linear system that needs to be solved---------------
!-------------------in order to obtain s_new----------------------------
!
call matmul_coo(aloline_data, aloline_colindx, aloline_rowindx, sline_vec, dummy_vec, ndxmax*ndymax*ndzmax, 7*ndxmax*ndymax*ndzmax)
!
sline_vec=(1.d0-eps_line)*dummy_vec
!
bnue_vec=eps_line*bnue_vec
!
mintbar_vec=(1.d0-eps_line)*mintbar_vec
!
dummy_vec=mintbar_vec-sline_vec+bnue_vec
!
aloline_data=(1.d0-eps_line)*aloline_data
aloline_data_diag=1.d0-(1.d0-eps_line)*aloline_data_diag
!
do i=1, 7*ndxmax*ndymax*ndzmax
   if(aloline_colindx(i).eq.aloline_rowindx(i)) then
      aloline_data(i)=1.d0-aloline_data(i)
   else
      aloline_data(i)=-1.d0*aloline_data(i)
   endif
enddo
!
!--------------check if matrix is diagonal dominant---------------------
!-----------------and estimate spectral radius--------------------------
!
sum_vec=0.d0
!
do i=1, 7*ndxmax*ndymax*ndzmax
   if(aloline_rowindx(i).ne.aloline_colindx(i)) then
      sum_vec(aloline_rowindx(i)) = sum_vec(aloline_rowindx(i)) + abs(aloline_data(i))
   endif
enddo
!
rspec=maxval(sum_vec/aloline_data_diag)
indx_rspec=maxloc(sum_vec/aloline_data_diag, 1)
call conv_indx_1d_to_3d(indx_rspec, ndxmax, ndymax, indx_x, indx_y, indx_z)
write(*,'(a30,es20.8)') 'estimate of spectral radius', rspec
write(*,'(a30,i10,3i5)') 'at 1d, 3d indices', indx_rspec, indx_x, indx_y, indx_z
write(*,*)
!
!--------------------linear system now reads:---------------------------
!              aloline_matrix * s_new = sline_vec
!
write(*,*) '-----------------------inverting nearest neighbour alo-------------------------'
call cpu_time(ts_alo)
!
dummy_vec=-dummy_vec
!
call jsor_coo(aloline_data, aloline_colindx, aloline_rowindx, aloline_data_diag, dummy_vec, 1.d0, ndxmax*ndymax*ndzmax, 7*ndxmax*ndymax*ndzmax, .false., sline_vec)
!
call cpu_time(te_alo)
ttot_alo=ttot_alo + te_alo-ts_alo
!
!write(*,*) 'total time needed to invert alo: ', te_alo-ts_alo
!
!---------back-transformation of source function on 3-d grid------------
!
!check if negative source function occurred
lng_expol=.false.
do i=1, ndxmax*ndymax*ndzmax
   if(mask_vec(i)) then
      if(sline_vec(i).lt.0.d0) then
         write(*,*) 'warning in sline_new3d_dn: new iterate gives negative source functions'
         wpa_interp2d=wpa_interp2d+1.d0
         wpb_interp2d=wpb_interp2d+1.d0
         wp_interp2d=wpa_interp2d/wpb_interp2d
         wpa_interp1d=wpa_interp1d+1.d0
         wpb_interp1d=wpb_interp1d+1.d0
         wp_interp1d=wpa_interp1d/wpb_interp1d
         wpa_integ1d=wpa_integ1d+1.d0
         wpb_integ1d=wpb_integ1d+1.d0
         wp_integ1d=wpa_integ1d/wpb_integ1d
         write(*,*) 'setting derivative weights for 2d bezier interpolation from, to', (wpa_interp2d-1.d0)/(wpb_interp2d-1.d0), wp_interp2d
         write(*,*) 'setting derivative weights for 1d bezier interpolation from, to', (wpa_interp1d-1.d0)/(wpb_interp1d-1.d0), wp_interp1d
         write(*,*) 'setting derivative weights for 1d bezier integration   from, to', (wpa_integ1d-1.d0)/(wpb_integ1d-1.d0), wp_integ1d
         write(*,*) 'using old source function'
         lng_expol=.true.
         return
      endif
   endif
enddo
!
!
!if no negative source function occured, go on
do i=1, ndxmax*ndymax*ndzmax
   if(mask_vec(i)) then
      call conv_indx_1d_to_3d(i, ndxmax, ndymax, indx_x, indx_y, indx_z)
      sline3d(indx_x, indx_y, indx_z) = sline_vec(i)
   endif
enddo
!
!check if negative source functions occurred (again? -> just to ensure that no shit has happened)
!                                                    -> costs no time
lnegative=.false.
do i=1, ndxmax
   do j=1, ndymax
      do k=1, ndzmax
         if(sline3d(i,j,k).lt.0.d0) then
            write(*,'(3i5,6es20.8,i5)') i, j, k, x(i), y(j), z(k), sline3d(i,j,k), ssobo3d(i,j,k), mintbar3d(i,j,k), imask3d(i,j,k)
            lnegative=.true.
         endif
      enddo
   enddo
enddo
if(lnegative) then
   write(*,*) 'error in sline_new3d_dn: new iterate gives negative source functions'
   write(*,*) 'try: linear interpolations or lower band alo'
   stop
endif
!
write(*,*) '-------------------------------------------------------------------------------'
write(*,*)
!
!-----------------------------------------------------------------------
!
end subroutine sline_new3d_dn
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine sline_new3d_nn
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
use dime3d, only: ndxmax, ndymax, ndzmax, sline3d, ssobo3d, mintbar3d, imask3d, t3d, x, y, z, &
                  aloline_nn3d, aloline_rowindx, aloline_colindx, aloline_data, aloline_data_diag
use params_input, only: eps_line
use freq, only: xnue0
use timing, only: ts_alo, te_alo, ttot_alo
use mod_interp2d, only: wpa_interp2d, wpb_interp2d, wp_interp2d, wpa_interp1d, &
                        wpb_interp1d, wp_interp1d, wpa_integ1d, wpb_integ1d, &
                        wp_integ1d, lng_expol
!
implicit none
!
! ... local scalars
integer(i4b) :: i, j, k
integer(i4b) :: indx_1d, indx_x, indx_y, indx_z, err, indx_rspec
real(dp) :: dummy1, dummy2, rspec
!
! ... local arrays
real(dp), dimension(:), allocatable :: mintbar_vec, sline_vec, bnue_vec, dummy_vec, sum_vec
logical, dimension(:), allocatable :: mask_vec
!
! ... local logicals
logical :: lnegative
!
! ... local functions
real(dp) :: bnue
!
!-----------------------------------------------------------------------
!---------------------allocation of local arrays------------------------
!-----------------------------------------------------------------------
!
allocate(bnue_vec(ndxmax*ndymax*ndzmax), stat=err)
   if(err.ne.0) stop 'allocation error sline_new3d_nn: bnue_vec'
!
allocate(sline_vec(ndxmax*ndymax*ndzmax), stat=err)
   if(err.ne.0) stop 'allocation error sline_new3d_nn: sline_vec'
!
allocate(mintbar_vec(ndxmax*ndymax*ndzmax), stat=err)
   if(err.ne.0) stop 'allocation error sline_new3d_nn: mintbar_vec'
!
allocate(dummy_vec(ndxmax*ndymax*ndzmax), stat=err)
   if(err.ne.0) stop 'allocation error sline_new3d_nn: dummy_vec'
!
allocate(mask_vec(ndxmax*ndymax*ndzmax), stat=err)
   if(err.ne.0) stop 'allocation error sline_new3d_nn: mask_vec'
!
allocate(sum_vec(ndxmax*ndymax*ndzmax), stat=err)
   if(err.ne.0) stop 'allocation error  sline_new3d_nn: sum_vec'
!
!----------------transform 3d-arrays to 1-d array-----------------------
!
mintbar_vec=0.d0
sline_vec=0.d0
bnue_vec=0.d0
dummy_vec=0.d0
mask_vec=.false.
!
do i=1, ndxmax
   do j=1, ndymax
      do k=1, ndzmax
         select case(imask3d(i,j,k))
            case(1,2,3)
               call conv_indx_3d_to_1d(i,j,k,ndxmax,ndymax,indx_1d)
               mintbar_vec(indx_1d)=mintbar3d(i,j,k)
               sline_vec(indx_1d)=sline3d(i,j,k)
               mask_vec(indx_1d)=.true.
               bnue_vec(indx_1d)=bnue(xnue0, t3d(i,j,k))
            case default
         endselect
      enddo
   enddo
enddo
!
call calc_aloline_nn_coo
!
!------------set up linear system that needs to be solved---------------
!-------------------in order to obtain s_new----------------------------
!
call matmul_coo(aloline_data, aloline_colindx, aloline_rowindx, sline_vec, dummy_vec, ndxmax*ndymax*ndzmax, 27*ndxmax*ndymax*ndzmax)
!
sline_vec=(1.d0-eps_line)*dummy_vec
!
bnue_vec=eps_line*bnue_vec
!
mintbar_vec=(1.d0-eps_line)*mintbar_vec
!
dummy_vec=mintbar_vec-sline_vec+bnue_vec
!
aloline_data=(1.d0-eps_line)*aloline_data
aloline_data_diag=1.d0-(1.d0-eps_line)*aloline_data_diag
!
do i=1, 27*ndxmax*ndymax*ndzmax
   if(aloline_colindx(i).eq.aloline_rowindx(i)) then
      aloline_data(i)=1.d0-aloline_data(i)
   else
      aloline_data(i)=-1.d0*aloline_data(i)
   endif
enddo
!
!--------------check if matrix is diagonal dominant---------------------
!-----------------and estimate spectral radius--------------------------
!
sum_vec=0.d0
!
do i=1, 27*ndxmax*ndymax*ndzmax
   if(aloline_rowindx(i).ne.aloline_colindx(i)) then
      sum_vec(aloline_rowindx(i)) = sum_vec(aloline_rowindx(i)) + abs(aloline_data(i))
   endif
enddo
!
rspec=maxval(sum_vec/aloline_data_diag)
indx_rspec=maxloc(sum_vec/aloline_data_diag, 1)
call conv_indx_1d_to_3d(indx_rspec, ndxmax, ndymax, indx_x, indx_y, indx_z)
write(*,'(a30,es20.8)') 'estimate of spectral radius', rspec
write(*,'(a30,i10,3i5)') 'at 1d, 3d indices', indx_rspec, indx_x, indx_y, indx_z
write(*,*)
!
!--------------------linear system now reads:---------------------------
!              aloline_matrix * s_new = sline_vec
!
write(*,*) '-----------------------inverting nearest neighbour alo-------------------------'
call cpu_time(ts_alo)
!
dummy_vec=-dummy_vec
!
call jsor_coo(aloline_data, aloline_colindx, aloline_rowindx, aloline_data_diag, dummy_vec, 1.d0, ndxmax*ndymax*ndzmax, 27*ndxmax*ndymax*ndzmax, .false., sline_vec)
!
call cpu_time(te_alo)
ttot_alo=ttot_alo + te_alo-ts_alo
!
!write(*,*) 'total time needed to invert alo: ', te_alo-ts_alo
!
!---------back-transformation of source function on 3-d grid------------
!
!check if negative source function occurred and ensure that ng_extrapolation is started
!from beginning if such situation occurr
lng_expol=.false.
do i=1, ndxmax*ndymax*ndzmax
   if(mask_vec(i)) then
      if(sline_vec(i).lt.0.d0) then
         write(*,*) 'warning in sline_new3d_nn: new iterate gives negative source functions'
         wpa_interp2d=wpa_interp2d+1.d0
         wpb_interp2d=wpb_interp2d+1.d0
         wp_interp2d=wpa_interp2d/wpb_interp2d
         wpa_interp1d=wpa_interp1d+1.d0
         wpb_interp1d=wpb_interp1d+1.d0
         wp_interp1d=wpa_interp1d/wpb_interp1d
         wpa_integ1d=wpa_integ1d+1.d0
         wpb_integ1d=wpb_integ1d+1.d0
         wp_integ1d=wpa_integ1d/wpb_integ1d
         write(*,*) 'setting derivative weights for 2d bezier interpolation from, to', (wpa_interp2d-1.d0)/(wpb_interp2d-1.d0), wp_interp2d
         write(*,*) 'setting derivative weights for 1d bezier interpolation from, to', (wpa_interp1d-1.d0)/(wpb_interp1d-1.d0), wp_interp1d
         write(*,*) 'setting derivative weights for 1d bezier integration   from, to', (wpa_integ1d-1.d0)/(wpb_integ1d-1.d0), wp_integ1d
         write(*,*) 'using old source function'
         lng_expol=.true.
         return
      endif
   endif
enddo
!
do i=1, ndxmax*ndymax*ndzmax
   if(mask_vec(i)) then
      call conv_indx_1d_to_3d(i, ndxmax, ndymax, indx_x, indx_y, indx_z)
      sline3d(indx_x, indx_y, indx_z) = sline_vec(i)
   endif
enddo
!
!check if negative source functions occurr (again? -> see sline_new3d_dn
lnegative=.false.
do i=1, ndxmax
   do j=1, ndymax
      do k=1, ndzmax
         if(sline3d(i,j,k).lt.0.d0) then
            write(*,'(3i5,6es20.8,i5)') i, j, k, x(i), y(j), z(k), sline3d(i,j,k), ssobo3d(i,j,k), mintbar3d(i,j,k), imask3d(i,j,k)
            lnegative=.true.
         endif
      enddo
   enddo
enddo
if(lnegative) then
   write(*,*) 'error in sline_new3d_nn: new iterate gives negative source functions'
   write(*,*) 'try: linear interpolations or lower band alo'
   stop
endif
!
write(*,*) '-------------------------------------------------------------------------------'
write(*,*)
!
!-----------------------------------------------------------------------
!
end subroutine sline_new3d_nn
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine calc_aloline_dn_coo
!
!-----------------------------------------------------------------------
!------------calculates approximate lambda operator---------------------
!-------using only direct neighbours (6 neighbours + local point)-------
!----------------------coo storage format-------------------------------
!-----------------------------------------------------------------------
!
use prog_type
use fund_const
use dime3d, only: ndxmax, ndymax, ndzmax, aloline_nn3d, &
                  aloline_data, aloline_data_diag, aloline_colindx, aloline_rowindx
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
if(.not.allocated(aloline_data)) then
   allocate(aloline_data(7*ndxmax*ndymax*ndzmax), stat=err)
      if(err.ne.0) stop 'allocation error calc_aloline_dn_coo: aloline_data'
endif

if(.not.allocated(aloline_data_diag)) then
   allocate(aloline_data_diag(ndxmax*ndymax*ndzmax), stat=err)
      if(err.ne.0) stop 'allocation error calc_aloline_dn_coo: aloline_data_diag'
endif
!
if(.not.allocated(aloline_colindx)) then
   allocate(aloline_colindx(7*ndxmax*ndymax*ndzmax), stat=err)
      if(err.ne.0) stop 'allocation error calc_aloline_dn_coo: aloline_col_indx'
endif
!
if(.not.allocated(aloline_rowindx)) then
   allocate(aloline_rowindx(7*ndxmax*ndymax*ndzmax), stat=err)
      if(err.ne.0) stop 'allocation error calc_aloline_dn_coo: aloline_rowindx'
endif
!
!-----------------calculate alo in coo storage format-------------------
!
!define maximum allowed spectral radius
rspec_max=1.d0
!
aloline_data=0.d0
aloline_data_diag=0.d0
aloline_rowindx=1
aloline_colindx=1
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
         rspec=rspec+abs(aloline_nn3d(i,j,k,14))
         call conv_indx_3d_to_1d(i,j,k, ndxmax, ndymax, indx_1d_col)
         aloline_data(iindx) = aloline_nn3d(i,j,k,14)
         aloline_colindx(iindx)=indx_1d_col
         aloline_rowindx(iindx)=indx_1d_col
         aloline_data_diag(indx_1d_col) = aloline_nn3d(i,j,k,14)
         ndiags_ijk=1
!
!include direct neighbour in positive x- direction
         rspec=rspec+abs(aloline_nn3d(i+1,j,k,13))
         if(rspec.lt.rspec_max) then
            call conv_indx_3d_to_1d(i+1,j,k,ndxmax,ndymax,indx_1d)
            aloline_data(iindx+1) = aloline_nn3d(i+1,j,k,13)
            aloline_colindx(iindx+1)=indx_1d
            aloline_rowindx(iindx+1)=indx_1d_col
            ndiags_ijk=ndiags_ijk+1
            ndiags_max = max(ndiags_ijk,ndiags_max)
         else
            ndiags_min = min(ndiags_ijk,ndiags_min)
            iindx=iindx+7
            cycle
         endif
!
!include direct neighbour in negative x-direction
         rspec=rspec+abs(aloline_nn3d(i-1,j,k,15))
         if(rspec.lt.rspec_max) then
            call conv_indx_3d_to_1d(i-1,j,k,ndxmax,ndymax,indx_1d)
            aloline_data(iindx+2) = aloline_nn3d(i-1,j,k,15)
            aloline_colindx(iindx+2)=indx_1d
            aloline_rowindx(iindx+2)=indx_1d_col
            ndiags_ijk=ndiags_ijk+1
            ndiags_max = max(ndiags_ijk,ndiags_max)
         else
            ndiags_min = min(ndiags_ijk,ndiags_min)
            iindx=iindx+7
            cycle
         endif
!
!include direct neighbour in positive y-direction
         rspec=rspec+abs(aloline_nn3d(i,j+1,k,11))
         if(rspec.lt.rspec_max) then
            call conv_indx_3d_to_1d(i,j+1,k,ndxmax,ndymax,indx_1d)
            aloline_data(iindx+3) = aloline_nn3d(i,j+1,k,11)
            aloline_colindx(iindx+3)=indx_1d
            aloline_rowindx(iindx+3)=indx_1d_col
            ndiags_ijk=ndiags_ijk+1
            ndiags_max = max(ndiags_ijk,ndiags_max)
         else
            ndiags_min = min(ndiags_ijk,ndiags_min)
            iindx=iindx+7
            cycle
         endif         
!
!include direct neighbour in negative y-direction
         rspec=rspec+abs(aloline_nn3d(i,j-1,k,17))
         if(rspec.lt.rspec_max) then
            call conv_indx_3d_to_1d(i,j-1,k,ndxmax,ndymax,indx_1d)
            aloline_data(iindx+4) = aloline_nn3d(i,j-1,k,17)
            aloline_colindx(iindx+4)=indx_1d
            aloline_rowindx(iindx+4)=indx_1d_col
            ndiags_ijk=ndiags_ijk+1
            ndiags_max = max(ndiags_ijk,ndiags_max)
         else
            ndiags_min = min(ndiags_ijk,ndiags_min)
            iindx=iindx+7
            cycle
         endif
!
!include direct neighbour in positive z-direction
         rspec=rspec+abs(aloline_nn3d(i,j,k+1,5))
         if(rspec.lt.rspec_max) then
            call conv_indx_3d_to_1d(i,j,k+1,ndxmax,ndymax,indx_1d)
            aloline_data(iindx+5) = aloline_nn3d(i,j,k+1,5)
            aloline_colindx(iindx+5)=indx_1d
            aloline_rowindx(iindx+5)=indx_1d_col
            ndiags_ijk=ndiags_ijk+1
            ndiags_max = max(ndiags_ijk,ndiags_max)
         else
            ndiags_min = min(ndiags_ijk,ndiags_min)
            iindx=iindx+7
            cycle
         endif
!
!include direct neighbour in negative z-direction
         rspec=rspec+abs(aloline_nn3d(i,j,k-1,23))
         if(rspec.lt.rspec_max) then
            call conv_indx_3d_to_1d(i,j,k-1,ndxmax,ndymax,indx_1d)
            aloline_data(iindx+6) = aloline_nn3d(i,j,k-1,23)
            aloline_colindx(iindx+6)=indx_1d
            aloline_rowindx(iindx+6)=indx_1d_col
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
!
!
end subroutine calc_aloline_dn_coo
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine calc_aloline_nn_coo
!
!-----------------------------------------------------------------------
!------------calculates approximate lambda operator---------------------
!-------using only nearest neighbours (26 neighbours + local point)-----
!----------------------coo storage format-------------------------------
!-----------------------------------------------------------------------
!
use prog_type
use fund_const
use dime3d, only: ndxmax, ndymax, ndzmax, aloline_nn3d, &
                  aloline_data, aloline_data_diag, aloline_colindx, aloline_rowindx, &
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
real(dp), dimension(:), allocatable :: aloline_data2, aloline_data_diag2
integer(i4b), dimension(:), allocatable :: aloline_rowindx2, aloline_colindx2
!
!-------------------allocation of alo-matrix----------------------------
!
if(.not.allocated(aloline_data)) then
   allocate(aloline_data(27*ndxmax*ndymax*ndzmax), stat=err)
      if(err.ne.0) stop 'allocation error calc_aloline_nn_coo: aloline_data'
endif

if(.not.allocated(aloline_data_diag)) then
   allocate(aloline_data_diag(ndxmax*ndymax*ndzmax), stat=err)
      if(err.ne.0) stop 'allocation error calc_aloline_nn_coo: aloline_data_diag'
endif
!
if(.not.allocated(aloline_colindx)) then
   allocate(aloline_colindx(27*ndxmax*ndymax*ndzmax), stat=err)
      if(err.ne.0) stop 'allocation error calc_aloline_nn_coo: aloline_col_indx'
endif
!
if(.not.allocated(aloline_rowindx)) then
   allocate(aloline_rowindx(27*ndxmax*ndymax*ndzmax), stat=err)
      if(err.ne.0) stop 'allocation error calc_aloline_nn_coo: aloline_rowindx'
endif
!
!-----------------calculate alo in coo storage format-------------------
!
!define maximum allowed spectral radius
rspec_max=1.d0
!
aloline_data=0.d0
aloline_data_diag=0.d0
aloline_rowindx=1
aloline_colindx=1
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
!              write(*,'(9es20.8)') aloline_nn3d(i,j,k,1), aloline_nn3d(i,j,k,2), aloline_nn3d(i,j,k,3), &
!                                   aloline_nn3d(i,j,k,4), aloline_nn3d(i,j,k,5), aloline_nn3d(i,j,k,6), &
!                                   aloline_nn3d(i,j,k,7), aloline_nn3d(i,j,k,8), aloline_nn3d(i,j,k,9), &
!                                   aloline_nn3d(i,j,k,10), aloline_nn3d(i,j,k,11), aloline_nn3d(i,j,k,12), &
!                                   aloline_nn3d(i,j,k,13), aloline_nn3d(i,j,k,14), aloline_nn3d(i,j,k,15), &
!                                   aloline_nn3d(i,j,k,16), aloline_nn3d(i,j,k,17), aloline_nn3d(i,j,k,18), &
!                                   aloline_nn3d(i,j,k,19), aloline_nn3d(i,j,k,20), aloline_nn3d(i,j,k,21), &
!                                   aloline_nn3d(i,j,k,22), aloline_nn3d(i,j,k,23), aloline_nn3d(i,j,k,24), &
!                                   aloline_nn3d(i,j,k,25), aloline_nn3d(i,j,k,26), aloline_nn3d(i,j,k,27)
!              write(*,*)
!        end select

         rspec=0.d0
!
!
!include diagonal part
         rspec=rspec+abs(aloline_nn3d(i,j,k,14))
         call conv_indx_3d_to_1d(i,j,k, ndxmax, ndymax, indx_1d_col)
!!         write(*,*) 'including diagonal', indx_1d_col, rspec
         aloline_data(iindx) = aloline_nn3d(i,j,k,14)
         aloline_colindx(iindx)=indx_1d_col
         aloline_rowindx(iindx)=indx_1d_col
         aloline_data_diag(indx_1d_col) = aloline_nn3d(i,j,k,14)
         ndiags_ijk=1
!
!
!!include direct neighbour in positive x- direction
         rspec=rspec+abs(aloline_nn3d(i+1,j,k,13))
         if(rspec.lt.rspec_max) then
            call conv_indx_3d_to_1d(i+1,j,k,ndxmax,ndymax,indx_1d)
!            write(*,*) 'including x+ diagonal', indx_1d, rspec
            aloline_data(iindx+1) = aloline_nn3d(i+1,j,k,13)
            aloline_colindx(iindx+1)=indx_1d
            aloline_rowindx(iindx+1)=indx_1d_col
            ndiags_ijk=ndiags_ijk+1
            ndiags_max = max(ndiags_ijk,ndiags_max)
         else
            ndiags_min = min(ndiags_ijk,ndiags_min)
            iindx=iindx+27
            cycle
         endif
!
!!include direct neighbour in negative x-direction
         rspec=rspec+abs(aloline_nn3d(i-1,j,k,15))
         if(rspec.lt.rspec_max) then
            call conv_indx_3d_to_1d(i-1,j,k,ndxmax,ndymax,indx_1d)
!            write(*,*) 'including x- diagonal', indx_1d, rspec
            aloline_data(iindx+2) = aloline_nn3d(i-1,j,k,15)
            aloline_colindx(iindx+2)=indx_1d
            aloline_rowindx(iindx+2)=indx_1d_col
            ndiags_ijk=ndiags_ijk+1
            ndiags_max = max(ndiags_ijk,ndiags_max)
         else
            ndiags_min = min(ndiags_ijk,ndiags_min)
            iindx=iindx+27
            cycle
         endif
!
!!include direct neighbour in positive y-direction
         rspec=rspec+abs(aloline_nn3d(i,j+1,k,11))
         if(rspec.lt.rspec_max) then
            call conv_indx_3d_to_1d(i,j+1,k,ndxmax,ndymax,indx_1d)
!            write(*,*) 'including y+ diagonal', indx_1d, rspec
            aloline_data(iindx+3) = aloline_nn3d(i,j+1,k,11)
            aloline_colindx(iindx+3)=indx_1d
            aloline_rowindx(iindx+3)=indx_1d_col
            ndiags_ijk=ndiags_ijk+1
            ndiags_max = max(ndiags_ijk,ndiags_max)
         else
            ndiags_min = min(ndiags_ijk,ndiags_min)
            iindx=iindx+27
            cycle
         endif         
!
!!include direct neighbour in negative y-direction
         rspec=rspec+abs(aloline_nn3d(i,j-1,k,17))
         if(rspec.lt.rspec_max) then
            call conv_indx_3d_to_1d(i,j-1,k,ndxmax,ndymax,indx_1d)
!            write(*,*) 'including y- diagonal', indx_1d, rspec
            aloline_data(iindx+4) = aloline_nn3d(i,j-1,k,17)
            aloline_colindx(iindx+4)=indx_1d
            aloline_rowindx(iindx+4)=indx_1d_col
            ndiags_ijk=ndiags_ijk+1
            ndiags_max = max(ndiags_ijk,ndiags_max)
         else
            ndiags_min = min(ndiags_ijk,ndiags_min)
            iindx=iindx+27
            cycle
         endif
!
!include direct neighbour in positive z-direction
         rspec=rspec+abs(aloline_nn3d(i,j,k+1,5))
         if(rspec.lt.rspec_max) then
            call conv_indx_3d_to_1d(i,j,k+1,ndxmax,ndymax,indx_1d)
!            write(*,*) 'including z+ diagonal', indx_1d, rspec
            aloline_data(iindx+5) = aloline_nn3d(i,j,k+1,5)
            aloline_colindx(iindx+5)=indx_1d
            aloline_rowindx(iindx+5)=indx_1d_col
            ndiags_ijk=ndiags_ijk+1
            ndiags_max = max(ndiags_ijk,ndiags_max)
         else
            ndiags_min = min(ndiags_ijk,ndiags_min)
            iindx=iindx+27
            cycle
         endif
!
!include direct neighbour in negative z-direction
         rspec=rspec+abs(aloline_nn3d(i,j,k-1,23))
         if(rspec.lt.rspec_max) then
            call conv_indx_3d_to_1d(i,j,k-1,ndxmax,ndymax,indx_1d)
!            write(*,*) 'inculding z- diagonal', indx_1d, rspec
            aloline_data(iindx+6) = aloline_nn3d(i,j,k-1,23)
            aloline_colindx(iindx+6)=indx_1d
            aloline_rowindx(iindx+6)=indx_1d_col
            ndiags_ijk=ndiags_ijk+1
            ndiags_max = max(ndiags_ijk,ndiags_max)
         else
            ndiags_min = min(ndiags_ijk,ndiags_min)
            iindx=iindx+27
            cycle
         endif
!
!include neighbour in positive x-y-direction
         rspec=rspec+abs(aloline_nn3d(i+1,j+1,k,10))
         if(rspec.lt.rspec_max) then
            call conv_indx_3d_to_1d(i+1,j+1,k,ndxmax,ndymax,indx_1d)
!            write(*,*) 'including x+y+ diagonal', indx_1d, rspec
            aloline_data(iindx+7) = aloline_nn3d(i+1,j+1,k,10)
            aloline_colindx(iindx+7)=indx_1d
            aloline_rowindx(iindx+7)=indx_1d_col
            ndiags_ijk=ndiags_ijk+1
            ndiags_max = max(ndiags_ijk,ndiags_max)
         else
            ndiags_min = min(ndiags_ijk,ndiags_min)
            iindx=iindx+27
            cycle
         endif

!include neighbour in negative x- and negative y-direction
         rspec=rspec+abs(aloline_nn3d(i-1,j-1,k,18))
         if(rspec.lt.rspec_max) then
            call conv_indx_3d_to_1d(i-1,j-1,k,ndxmax,ndymax,indx_1d)
!            write(*,*) 'including x-y- diagonal', indx_1d, rspec
            aloline_data(iindx+8) = aloline_nn3d(i-1,j-1,k,18)
            aloline_colindx(iindx+8)=indx_1d
            aloline_rowindx(iindx+8)=indx_1d_col
            ndiags_ijk=ndiags_ijk+1
            ndiags_max = max(ndiags_ijk,ndiags_max)
         else
            ndiags_min = min(ndiags_ijk,ndiags_min)
            iindx=iindx+27
            cycle
         endif

!include neighbour in positive x- and negative y-direction
         rspec=rspec+abs(aloline_nn3d(i+1,j-1,k,16))
         if(rspec.lt.rspec_max) then
            call conv_indx_3d_to_1d(i+1,j-1,k,ndxmax,ndymax,indx_1d)
!            write(*,*) 'including x+y- diagonal', indx_1d, rspec
            aloline_data(iindx+9) = aloline_nn3d(i+1,j-1,k,16)
            aloline_colindx(iindx+9)=indx_1d
            aloline_rowindx(iindx+9)=indx_1d_col
            ndiags_ijk=ndiags_ijk+1
            ndiags_max = max(ndiags_ijk,ndiags_max)
         else
            ndiags_min = min(ndiags_ijk,ndiags_min)
            iindx=iindx+27
            cycle
         endif

!include neighbour in negative x- and positive y-direction
         rspec=rspec+abs(aloline_nn3d(i-1,j+1,k,12))
         if(rspec.lt.rspec_max) then
            call conv_indx_3d_to_1d(i-1,j+1,k,ndxmax,ndymax,indx_1d)
!            write(*,*) 'including x-y+ diagonal', indx_1d, rspec
            aloline_data(iindx+10) = aloline_nn3d(i-1,j+1,k,12)
            aloline_colindx(iindx+10)=indx_1d
            aloline_rowindx(iindx+10)=indx_1d_col
            ndiags_ijk=ndiags_ijk+1
            ndiags_max = max(ndiags_ijk,ndiags_max)
         else
            ndiags_min = min(ndiags_ijk,ndiags_min)
            iindx=iindx+27
            cycle
         endif

!including neighbour in positive x- and positive z-direction
         rspec=rspec+abs(aloline_nn3d(i+1,j,k+1,4))
         if(rspec.lt.rspec_max) then
            call conv_indx_3d_to_1d(i+1,j,k+1,ndxmax,ndymax,indx_1d)
!            write(*,*) 'including x+z+ diagonal', indx_1d, rspec
            aloline_data(iindx+11) = aloline_nn3d(i+1,j,k+1,4)
            aloline_colindx(iindx+11)=indx_1d
            aloline_rowindx(iindx+11)=indx_1d_col
            ndiags_ijk=ndiags_ijk+1
            ndiags_max = max(ndiags_ijk,ndiags_max)
         else
            ndiags_min = min(ndiags_ijk,ndiags_min)
            iindx=iindx+27
            cycle
         endif

!including neighbour in negative x- and negative z-direction
         rspec=rspec+abs(aloline_nn3d(i-1,j,k-1,24))
         if(rspec.lt.rspec_max) then
            call conv_indx_3d_to_1d(i-1,j,k-1,ndxmax,ndymax,indx_1d)
!            write(*,*) 'including x-z- diagonal', indx_1d, rspec
            aloline_data(iindx+12) = aloline_nn3d(i-1,j,k-1,24)
            aloline_colindx(iindx+12)=indx_1d
            aloline_rowindx(iindx+12)=indx_1d_col
            ndiags_ijk=ndiags_ijk+1
            ndiags_max = max(ndiags_ijk,ndiags_max)
         else
            ndiags_min = min(ndiags_ijk,ndiags_min)
            iindx=iindx+27
            cycle
         endif

!including neighbour in positive x- and negative z-direction
         rspec=rspec+abs(aloline_nn3d(i+1,j,k-1,22))
         if(rspec.lt.rspec_max) then
            call conv_indx_3d_to_1d(i+1,j,k-1,ndxmax,ndymax,indx_1d)
!            write(*,*) 'including x+z- diagonal', indx_1d, rspec
            aloline_data(iindx+13) = aloline_nn3d(i+1,j,k-1,22)
            aloline_colindx(iindx+13)=indx_1d
            aloline_rowindx(iindx+13)=indx_1d_col
            ndiags_ijk=ndiags_ijk+1
            ndiags_max = max(ndiags_ijk,ndiags_max)
         else
            ndiags_min = min(ndiags_ijk,ndiags_min)
            iindx=iindx+27
            cycle
         endif

!including neighbour in negative x- and positive z-direction
         rspec=rspec+abs(aloline_nn3d(i-1,j,k+1,6))
         if(rspec.lt.rspec_max) then
            call conv_indx_3d_to_1d(i-1,j,k+1,ndxmax,ndymax,indx_1d)
!            write(*,*) 'including x-z+ diagonal', indx_1d, rspec
            aloline_data(iindx+14) = aloline_nn3d(i-1,j,k+1,6)
            aloline_colindx(iindx+14)=indx_1d
            aloline_rowindx(iindx+14)=indx_1d_col
            ndiags_ijk=ndiags_ijk+1
            ndiags_max = max(ndiags_ijk,ndiags_max)
         else
            ndiags_min = min(ndiags_ijk,ndiags_min)
            iindx=iindx+27
            cycle
         endif

!include neighbour in positive y- and positive z-direction
         rspec=rspec+abs(aloline_nn3d(i,j+1,k+1,2))
         if(rspec.lt.rspec_max) then
         call conv_indx_3d_to_1d(i,j+1,k+1,ndxmax,ndymax,indx_1d)
!            write(*,*) 'including y+z+ diagonal', indx_1d, rspec
            aloline_data(iindx+15) = aloline_nn3d(i,j+1,k+1,2)
            aloline_colindx(iindx+15)=indx_1d
            aloline_rowindx(iindx+15)=indx_1d_col
            ndiags_ijk=ndiags_ijk+1
            ndiags_max = max(ndiags_ijk,ndiags_max)
         else
            ndiags_min = min(ndiags_ijk,ndiags_min)
            iindx=iindx+27
            cycle
         endif

!include neighbour in negative y- and negative z-direction
         rspec=rspec+abs(aloline_nn3d(i,j-1,k-1,26))
         if(rspec.lt.rspec_max) then
            call conv_indx_3d_to_1d(i,j-1,k-1,ndxmax,ndymax,indx_1d)
!            write(*,*) 'including y-z- diagonal', indx_1d, rspec
            aloline_data(iindx+16) = aloline_nn3d(i,j-1,k-1,26)
            aloline_colindx(iindx+16)=indx_1d
            aloline_rowindx(iindx+16)=indx_1d_col
            ndiags_ijk=ndiags_ijk+1
            ndiags_max = max(ndiags_ijk,ndiags_max)
         else
            ndiags_min = min(ndiags_ijk,ndiags_min)
            iindx=iindx+27
            cycle
         endif
!
!include neighbour in positive y- and negative z-direction
         rspec=rspec+abs(aloline_nn3d(i,j+1,k-1,20))
         if(rspec.lt.rspec_max) then
            call conv_indx_3d_to_1d(i,j+1,k-1,ndxmax,ndymax,indx_1d)
!            write(*,*) 'including y+z-diagonal', indx_1d, rspec
            aloline_data(iindx+17) = aloline_nn3d(i,j+1,k-1,20)
            aloline_colindx(iindx+17)=indx_1d
            aloline_rowindx(iindx+17)=indx_1d_col
            ndiags_ijk=ndiags_ijk+1
            ndiags_max = max(ndiags_ijk,ndiags_max)
         else
            ndiags_min = min(ndiags_ijk,ndiags_min)
            iindx=iindx+27
            cycle
         endif

!include neighbour in negative y- and positive z-direction
         rspec=rspec+abs(aloline_nn3d(i,j-1,k+1,8))
         if(rspec.lt.rspec_max) then
            call conv_indx_3d_to_1d(i,j-1,k+1,ndxmax,ndymax,indx_1d)
!            write(*,*) 'including y-z+ diagonal', indx_1d, rspec
            aloline_data(iindx+18) = aloline_nn3d(i,j-1,k+1,8)
            aloline_colindx(iindx+18)=indx_1d
            aloline_rowindx(iindx+18)=indx_1d_col
            ndiags_ijk=ndiags_ijk+1
            ndiags_max = max(ndiags_ijk,ndiags_max)
         else
            ndiags_min = min(ndiags_ijk,ndiags_min)
            iindx=iindx+27
            cycle
         endif

!include neighbour in positive x-, positive y- and positive z-direction
         rspec=rspec+abs(aloline_nn3d(i+1,j+1,k+1,1))
         if(rspec.lt.rspec_max) then
            call conv_indx_3d_to_1d(i+1,j+1,k+1,ndxmax,ndymax,indx_1d)
!            write(*,*) 'including x+y+z+ diagonal', indx_1d, rspec
            aloline_data(iindx+19) = aloline_nn3d(i+1,j+1,k+1,1)
            aloline_colindx(iindx+19)=indx_1d
            aloline_rowindx(iindx+19)=indx_1d_col
            ndiags_ijk=ndiags_ijk+1
            ndiags_max = max(ndiags_ijk,ndiags_max)
         else
            ndiags_min = min(ndiags_ijk,ndiags_min)
            iindx=iindx+27
            cycle
         endif

!include neighbour in negative x-, negative y- and negative z-direction
         rspec=rspec+abs(aloline_nn3d(i-1,j-1,k-1,27))
         if(rspec.lt.rspec_max) then
            call conv_indx_3d_to_1d(i-1,j-1,k-1,ndxmax,ndymax,indx_1d)
!            write(*,*) 'including x-y-z- diagonal', indx_1d, rspec
            aloline_data(iindx+20) = aloline_nn3d(i-1,j-1,k-1,27)
            aloline_colindx(iindx+20)=indx_1d
            aloline_rowindx(iindx+20)=indx_1d_col
            ndiags_ijk=ndiags_ijk+1
            ndiags_max = max(ndiags_ijk,ndiags_max)
         else
            ndiags_min = min(ndiags_ijk,ndiags_min)
            iindx=iindx+27
            cycle
         endif
!
!include neighbour in positive x-, negative y- and positive z-direction
         rspec=rspec+abs(aloline_nn3d(i+1,j-1,k+1,7))
         if(rspec.lt.rspec_max) then
            call conv_indx_3d_to_1d(i+1,j-1,k+1,ndxmax,ndymax,indx_1d)
!            write(*,*) 'including x-y-z+ diagonal', indx_1d, rspec
            aloline_data(iindx+21) = aloline_nn3d(i+1,j-1,k+1,7)
            aloline_colindx(iindx+21)=indx_1d
            aloline_rowindx(iindx+21)=indx_1d_col
            ndiags_ijk=ndiags_ijk+1
            ndiags_max = max(ndiags_ijk,ndiags_max)
         else
            ndiags_min = min(ndiags_ijk,ndiags_min)
            iindx=iindx+27
            cycle
         endif
!
!include neighbour in negative x-, positive y- and negative z-direction
         rspec=rspec+abs(aloline_nn3d(i-1,j+1,k-1,21))
         if(rspec.lt.rspec_max) then
            call conv_indx_3d_to_1d(i-1,j+1,k-1,ndxmax,ndymax,indx_1d)
!            write(*,*) 'including x-y+z- diagonal', indx_1d, rspec
            aloline_data(iindx+22) = aloline_nn3d(i-1,j+1,k-1,21)
            aloline_colindx(iindx+22)=indx_1d
            aloline_rowindx(iindx+22)=indx_1d_col
            ndiags_ijk=ndiags_ijk+1
            ndiags_max = max(ndiags_ijk,ndiags_max)
         else
            ndiags_min = min(ndiags_ijk,ndiags_min)
            iindx=iindx+27
            cycle
         endif
!
!include neighbour in positive x-, positive y- and negative z-direction
         rspec=rspec+abs(aloline_nn3d(i+1,j+1,k-1,19))
         if(rspec.lt.rspec_max) then
            call conv_indx_3d_to_1d(i+1,j+1,k-1,ndxmax,ndymax,indx_1d)
!            write(*,*) 'including x+y+z- diagonal', indx_1d, rspec
            aloline_data(iindx+23) = aloline_nn3d(i+1,j+1,k-1,19)
            aloline_colindx(iindx+23)=indx_1d
            aloline_rowindx(iindx+23)=indx_1d_col
            ndiags_ijk=ndiags_ijk+1
            ndiags_max = max(ndiags_ijk,ndiags_max)
         else
            ndiags_min = min(ndiags_ijk,ndiags_min)
            iindx=iindx+27
            cycle
         endif
!
!include neighbour in negative x-, negative y- and positive z-direction
         rspec=rspec+abs(aloline_nn3d(i-1,j-1,k+1,9))
         if(rspec.lt.rspec_max) then
            call conv_indx_3d_to_1d(i-1,j-1,k+1,ndxmax,ndymax,indx_1d)
!            write(*,*) 'including x-y-z+ diagonal', indx_1d, rspec
            aloline_data(iindx+24) = aloline_nn3d(i-1,j-1,k+1,9)
            aloline_colindx(iindx+24)=indx_1d
            aloline_rowindx(iindx+24)=indx_1d_col
            ndiags_ijk=ndiags_ijk+1
            ndiags_max = max(ndiags_ijk,ndiags_max)
         else
            ndiags_min = min(ndiags_ijk,ndiags_min)
            iindx=iindx+27
            cycle
         endif

!include neighbour in positive x-, negative y- and negative z-direction
         rspec=rspec+abs(aloline_nn3d(i+1,j-1,k-1,25))
         if(rspec.lt.rspec_max) then
            call conv_indx_3d_to_1d(i+1,j-1,k-1,ndxmax,ndymax,indx_1d)
!            write(*,*) 'including x+y-z- diagonal', indx_1d, rspec
            aloline_data(iindx+25) = aloline_nn3d(i+1,j-1,k-1,25)
            aloline_colindx(iindx+25)=indx_1d
            aloline_rowindx(iindx+25)=indx_1d_col
            ndiags_ijk=ndiags_ijk+1
            ndiags_max = max(ndiags_ijk,ndiags_max)
         else
            ndiags_min = min(ndiags_ijk,ndiags_min)
            iindx=iindx+27
            cycle
         endif

!include neighbour in negative x-, positive y- and positive z-direction
         rspec=rspec+abs(aloline_nn3d(i-1,j+1,k+1,3))
         if(rspec.lt.rspec_max) then
            call conv_indx_3d_to_1d(i-1,j+1,k+1,ndxmax,ndymax,indx_1d)
!            write(*,*) 'including x-y+z+ diagonal', indx_1d, rspec
            aloline_data(iindx+26) = aloline_nn3d(i-1,j+1,k+1,3)
            aloline_colindx(iindx+26)=indx_1d
            aloline_rowindx(iindx+26)=indx_1d_col
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
!
end subroutine calc_aloline_nn_coo
