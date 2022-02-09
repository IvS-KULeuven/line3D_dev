!***********************************************************************
!***********************************************************************
!
!   benchmark01 performas following test:
!
!   -calculate opacity and source function model (power-law in 1/r)
!   -calculate intensity along a predefined ray in 2d (short characteristics and fvm)
!   -calculate intensity along a predefined ray in 1d (extra grid calculated; short characteristics and fvm)
!   -calculate intensity along a predefined ray analytically
!   -output everything (plots in plotFILES/benchmark01.pro)
!
!***********************************************************************
!***********************************************************************
!
subroutine benchmark01_grid
!
use prog_type
use inf_reg, only: rmin
use params_input, only: rmax
use mod_benchmark, only: nr, r1d, scont1d, opac1d, tau1d, tau_min, tau_max, &
                         source_min, source_max, im_source, im_opacity, &
                         x_u2d, z_u2d, x_d2d, z_d2d, abs2d_sc, abs2d_fvm, contr2d_sc, contr2d_fvm, &
                         scont_u2d, opac_u2d, scont_d2d, opac_d2d, int_u2d, &
                         xu_interp_int2d, xu_interp_opac2d, zu_interp_int2d, zu_interp_opac2d, &
                         xd_interp_opac2d, zd_interp_opac2d, opac2d, scont2d, int2d_sc, int2d_fvm, n_z
use dime3d, only: imask_innreg3d, x, z, opac3d, scont3d, ndxmax, ndymax, ndzmax
use bcondition, only: xic1, xic2
use angles, only: nodes_mu
use mod_interp1d, only: find_index, interpol_fyp_cube3
!
implicit none
!
! ... local scalars
integer(i4b) :: i, j, k, iim2, iim1, ii, iip1, err
real(dp) :: c, tau_local, rad
real(dp) :: rmax_dum
!
! ... local functions
real(dp) :: benchmark01_source
!
!-----------------boundary condition and considered direction-----------
!
allocate(opac2d(ndxmax,ndzmax), stat=err)
   if(err.ne.0) stop 'allocation error in benchmark01_grid'
allocate(scont2d(ndxmax,ndzmax), stat=err)
   if(err.ne.0) stop 'allocation error in benchmark01_grid'
!
xic1=1.d0
xic2=0.d0
!
nodes_mu=n_z
!
!---------------------------spatial grid--------------------------------
!
rmax_dum=sqrt(2.d0)*x(ndxmax)
do i=2, nr-1
   r1d(i)=rmin+(i-2)*(rmax_dum-rmin)/(nr-3)
enddo
!
!include ghospoints
r1d(1)=r1d(2)-(r1d(3)-r1d(2))
r1d(nr)=r1d(nr-1)+(r1d(nr-1)-r1d(nr-2))
!
!---------------------------opacity model-------------------------------
!oapcity as powerlaw, such that tau_max is always hit:
!   chi=c/s^im_opacity
if(im_opacity.lt.0) then
   stop 'error in grid1d_model: n less than 0 not allowed'
elseif(im_opacity.eq.0) then
   c=(tau_max-tau_min)/(rmax-rmin)
elseif(im_opacity.eq.1) then
   c=(tau_max-tau_min)/log(rmax/rmin)
else
   c=(1-im_opacity)*(tau_max-tau_min)/(rmax**(1-im_opacity)-rmin**(1-im_opacity))
endif
!
opac1d=0.d0
do i=2, nr-1
   opac1d(i)=c/r1d(i)**im_opacity
enddo
!
!----------------calculate corresponding tau-grid-----------------------
!
tau1d=0.d0
if(im_opacity.eq.0) then
   do i=2, nr
      tau1d(i) = tau_min + c*(r1d(i)-rmin)
   enddo
elseif(im_opacity.eq.1) then
   do i=2, nr
      tau1d(i) = tau_min + c*log(r1d(i)/rmin)
   enddo
else
   do i=2, nr
      tau1d(i) = tau_min + c*(r1d(i)**(1-im_opacity) - rmin**(1-im_opacity))/(1-im_opacity)
   enddo
endif
!
!------------source function model--------------------------------------
!
!different source function models
do i=1, nr
   scont1d(i) = benchmark01_source(tau1d(i), tau1d(2), tau1d(nr), source_min, source_max, im_source)
enddo
!
!******************************3D grid**********************************
!
opac2d=0.d0
scont2d=0.d0
opac3d=0.d0
scont3d=0.d0
!
j=ndymax/2+1
do i=1, ndxmax
   do k=1, ndzmax
      if(imask_innreg3d(i,j,k).ne.1) then
         rad=sqrt(x(i)**2 + z(k)**2)
!opacity
         opac2d(i,k)=c/rad**im_opacity
         opac3d(i,j,k)=opac2d(i,k)
       
!source function
         call find_index(rad, r1d, nr, iim2, iim1, ii, iip1)
         tau_local = interpol_fyp_cube3(tau1d(iim2), tau1d(iim1), tau1d(ii), tau1d(iip1), &
                                       r1d(iim2), r1d(iim1), r1d(ii), r1d(iip1), &
                                       rad, ii, nr)
         scont2d(i,k) = benchmark01_source(tau_local, tau1d(2), tau1d(nr), source_min, source_max, im_source)
         scont3d(i,j,k) = scont2d(i,k)
      endif
   enddo
enddo
!
!upwind and downwind coordinates
allocate(x_u2d(ndxmax,ndzmax), stat=err)
   if(err.ne.0) stop 'allocation error in benchmark01_grid'
allocate(z_u2d(ndxmax,ndzmax), stat=err)
   if(err.ne.0) stop 'allocation error in benchmark01_grid'
allocate(x_d2d(ndxmax,ndzmax), stat=err)
   if(err.ne.0) stop 'allocation error in benchmark01_grid'
allocate(z_d2d(ndxmax,ndzmax), stat=err)
   if(err.ne.0) stop 'allocation error in benchmark01_grid'
!absorption and source-contribution part
allocate(int2d_sc(ndxmax,ndzmax), stat=err)
   if(err.ne.0) stop 'allocation error in benchmark01_grid'
allocate(int2d_fvm(ndxmax,ndzmax), stat=err)
   if(err.ne.0) stop 'allocation error in benchmark01_grid'
allocate(abs2d_sc(ndxmax,ndzmax), stat=err)
   if(err.ne.0) stop 'allocation error in benchmark01_grid'
allocate(abs2d_fvm(ndxmax,ndzmax), stat=err)
   if(err.ne.0) stop 'allocation error in benchmark01_grid'
allocate(contr2d_sc(ndxmax,ndzmax), stat=err)
   if(err.ne.0) stop 'allocation error in benchmark01_grid'
allocate(contr2d_fvm(ndxmax,ndzmax), stat=err)
   if(err.ne.0) stop 'allocation error in benchmark01_grid'
!intensities, opacities and source functions at upwind and downwind points
allocate(int_u2d(ndxmax,ndzmax), stat=err)
   if(err.ne.0) stop 'allocation error in benchmark01_grid'
allocate(scont_u2d(ndxmax,ndzmax), stat=err)
   if(err.ne.0) stop 'allocation error in benchmark01_grid'
allocate(scont_d2d(ndxmax,ndzmax), stat=err)
   if(err.ne.0) stop 'allocation error in benchmark01_grid'
allocate(opac_u2d(ndxmax,ndzmax), stat=err)
   if(err.ne.0) stop 'allocation error in benchmark01_grid'
allocate(opac_d2d(ndxmax,ndzmax), stat=err)
   if(err.ne.0) stop 'allocation error in benchmark01_grid'
!coordinates from which everything has been interpolated
allocate(xu_interp_int2d(ndxmax,ndzmax,8), stat=err)
   if(err.ne.0) stop 'allocation error in benchmark01_grid'
allocate(xu_interp_opac2d(ndxmax,ndzmax,8), stat=err)
   if(err.ne.0) stop 'allocation error in benchmark01_grid'
allocate(zu_interp_int2d(ndxmax,ndzmax,8), stat=err)
   if(err.ne.0) stop 'allocation error in benchmark01_grid'
allocate(zu_interp_opac2d(ndxmax,ndzmax,8), stat=err)
   if(err.ne.0) stop 'allocation error in benchmark01_grid'
allocate(xd_interp_opac2d(ndxmax,ndzmax,8), stat=err)
   if(err.ne.0) stop 'allocation error in benchmark01_grid'
allocate(zd_interp_opac2d(ndxmax,ndzmax,8), stat=err)
   if(err.ne.0) stop 'allocation error in benchmark01_grid'
!
end subroutine benchmark01_grid
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
function benchmark01_source(tau, tau_min, tau_max, source_min, source_max, indx_source)
!
use prog_type
!
implicit none
!
! ... arguments
integer(i4b), intent(in) :: indx_source
real(dp), intent(in) :: tau, tau_min, tau_max, source_min, source_max
real(dp) :: benchmark01_source
!
select case(indx_source)
   case(0)
!constant in tau
      benchmark01_source=source_min
   case(1)
!constant in tau
      benchmark01_source=source_max
   case(2)
!linearly increasing in tau
      benchmark01_source=source_min + (source_max-source_min)*(tau-tau_min)/(tau_max-tau_min)
   case(3)
!linearly decreasing in tau
      benchmark01_source = source_max + (source_min-source_max)*(tau-tau_min)/(tau_max-tau_min)
   case(4)
!exponentially increasing in tau
      benchmark01_source=source_min*exp(tau-tau_min)
   case(5)
!exponentially decreasing in tau
      benchmark01_source=source_max*exp(-tau_max - tau)
   case default
      stop 'error in funct: undefined index'
end select
!
end function benchmark01_source

!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
function benchmark01_sourcecontr(tau_i, tau_im1, s_i, s_im1, indx_source)
!
!calculate source contribtion for each cell given a source-fct model
!
use prog_type
!
implicit none
!
! ... arguments
integer(i4b), intent(in) :: indx_source
real(dp), intent(in) :: tau_i, tau_im1, s_i, s_im1
real(dp) :: benchmark01_sourcecontr
!
! ... local scalars
real(dp) :: del_tau
!
del_tau=tau_i-tau_im1
!
select case(indx_source)
   case(0)
!constant in tau
      benchmark01_sourcecontr = s_i * (1.d0-exp(-del_tau))
   case(1)
!constant in tau
      benchmark01_sourcecontr = s_i * (1.d0-exp(-del_tau))
   case(2)
!linearly increasing in tau
      benchmark01_sourcecontr = (1.d0-exp(-del_tau))*s_im1 + (del_tau+exp(-del_tau)-1.d0)*(s_i-s_im1)/(tau_i-tau_im1)
   case(3)
!linearly decreasing in tau
      benchmark01_sourcecontr = (1.d0-exp(-del_tau))*s_im1 + (del_tau+exp(-del_tau)-1.d0)*(s_i-s_im1)/(tau_i-tau_im1)
   case(4)
!exponentially increasing in tau
      benchmark01_sourcecontr = 0.5d0*s_im1*(exp(del_tau)-exp(-del_tau))
   case(5)
!exponentially decreasing in tau
      benchmark01_sourcecontr = s_im1*del_tau*exp(-del_tau)
   case default
      stop 'error in funct: undefined index'
end select
!
end function benchmark01_sourcecontr
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine benchmark01_solution
!
use prog_type
use bcondition, only: xic1
use mod_benchmark, only: nr, tau1d, opac1d, scont1d, int1d, int1d_sc, int1d_fvm, &
                         abs1d, abs1d_sc, abs1d_fvm, contr1d, contr1d_sc, contr1d_fvm, r1d, im_source
!
implicit none
!
! ... local scalars
integer(i4b) :: i
real(dp) :: int_u, opac_u, opac_p, opac_d, scont_u, scont_p, scont_d, &
            dels_u, dels_d, abs_sc, contr_sc, int_sc, del_tau, alo_mm, alo_u, alo_d
!
! ... local functions
real(dp) :: benchmark01_sourcecontr
!
!
int1d(1) = xic1
int1d(2) = xic1
int1d_sc(1) = xic1
int1d_sc(2) = xic1
int1d_fvm(1) = xic1
int1d_fvm(2) = xic1
!
do i=3, nr-1
!
!------------------------theoretical solution---------------------------
!
   abs1d(i) = exp(tau1d(i-1)-tau1d(i))
   contr1d(i) =  benchmark01_sourcecontr(tau1d(i), tau1d(i-1), scont1d(i), scont1d(i-1), im_source)
   int1d(i) = int1d(i-1)*exp(tau1d(i-1)-tau1d(i)) + contr1d(i)
!
!------------------------1d short characteristics-----------------------
!
   int_u = int1d_sc(i-1)
!
   opac_u=opac1d(i-1)
   opac_p=opac1d(i)
   opac_d=opac1d(i+1)
!
   scont_u=scont1d(i-1)
   scont_p=scont1d(i)
   scont_d=scont1d(i+1)
!
   dels_u=r1d(i)-r1d(i-1)
   dels_d=r1d(i+1)-r1d(i)
!
   call fsc_cont(int_u, opac_u, opac_p, opac_d, scont_u, scont_p, scont_d, &
                  dels_u, dels_d, abs_sc, contr_sc, int_sc, alo_mm, alo_u, alo_d)

   abs1d_sc(i) = abs_sc
   contr1d_sc(i) =  contr_sc
   int1d_sc(i) = int_sc
!
!-----------------------------1d finite volume--------------------------
!
   del_tau=0.5d0*opac1d(i)*(r1d(i+1)-r1d(i-1))
   int1d_fvm(i) = del_tau/(1.d0+del_tau)*scont1d(i) + int1d_fvm(i-1)/(1.d0+del_tau)
   abs1d_fvm(i) = 1.d0/(1.d0+del_tau)
   contr1d_fvm(i) = del_tau/(1.d0+del_tau)*scont1d(i)
!
enddo
!
!
!
abs1d(1)=abs1d(3)
abs1d(2)=abs1d(3)
abs1d_sc(2)=abs1d_sc(3)
abs1d_sc(1)=abs1d_sc(3)
abs1d_fvm(2)=abs1d_fvm(3)
abs1d_fvm(1)=abs1d_fvm(3)
!
!
!
end subroutine benchmark01_solution
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine benchmark01_ray
!
use prog_type
use mod_benchmark, only: nr, r1d, int1d_scray, int1d_fvmray, n_z, int2d_sc, int2d_fvm
use dime3d, only: ndxmax, ndymax, ndzmax, x, z
use inf_reg, only: rmin
use params_input, only: rmax
use mod_interp2d, only: get_xy_indx, get_xy_values1, get_xy_values2, bilin
!
implicit none
!
! ... local scalars
integer(i4b) :: i, err
integer(i4b) :: ii_x, iim1_x, ii_z, iim1_z
real(dp) :: xp, zp, vala, valb, valc, vald, valp, x_iim1, x_ii, z_iim1, z_ii
real(dp) :: n_x, rlim
!
! ... local arrays
!
! ... local logicals
logical :: llogx, llogz, llogf, lexpol
!
n_x=sqrt(1.d0-n_z**2)
if(n_x.gt.n_z) then
   rlim=z(ndzmax-1)/n_z
elseif(n_x.lt.n_z) then
   rlim=x(ndxmax-1)/n_x
else
   rlim=sqrt(2.d0)*rmax
endif
!
!
do i=1, nr
!calculate x and z coordinates for given radius and direction
   xp=r1d(i)*n_x
   zp=r1d(i)*n_z
   if(r1d(i).le.rlim) then
!find index on x and z axes
      call get_xy_indx(xp, zp, x, z, ndxmax, ndzmax, iim1_x, ii_x, iim1_z, ii_z, lexpol, rmin, rlim)

!get values on neighbouring points
      call get_xy_values1(xp, zp, x, z, ndxmax, ndzmax, &
                          iim1_x, ii_x, iim1_z, ii_z, &
                          x_iim1, x_ii, z_iim1, z_ii, llogx, llogz)
!
!--------------------interpolating intensities--------------------------
!
!for short characteristics method
      call get_xy_values2(ndxmax, ndzmax, int2d_sc, &
                          iim1_x, ii_x, iim1_z, ii_z, &
                          vala, valb, valc, vald, llogf)
      call bilin(xp, zp, x_iim1, x_ii, z_iim1, z_ii, &
                 vala, valb, valc, vald, llogx, llogz, llogf, valp)
      int1d_scray(i)=valp

!for finite volume method
      call get_xy_values2(ndxmax, ndzmax, int2d_fvm, &
                          iim1_x, ii_x, iim1_z, ii_z, &
                          vala, valb, valc, vald, llogf)
      call bilin(xp, zp, x_iim1, x_ii, z_iim1, z_ii, &
                 vala, valb, valc, vald, llogx, llogz, llogf, valp)
      int1d_fvmray(i)=valp

      if(lexpol) then
         int1d_scray(i)=0.d0
         int1d_fvmray(i)=0.d0
      endif
!
   endif
enddo

end subroutine benchmark01_ray
