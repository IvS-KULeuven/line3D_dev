!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine calcnodes_omega
!
   use prog_type
   use angles, only: dim_mu, dim_phi, weight_omega, dim_omega, q_alo, n_x, n_y, n_z
   use options, only: opt_angint_method
!
   implicit none
!
! ... local scalars
   integer(i4b) :: i, err
!
!
   write(*,*) '----------------calculating number of mu and phi grid points-------------------'
   write(*,*)
!
   select case (opt_angint_method)
!
    case(0)
!trapezoidal rule
      call calcnodes_omega0
    case(1)
!trapezoidal rule with direcitonal distribution from lobell&blomme
      call calcnodes_omega1
    case(2)
!simpsons rule
      call calcnodes_omega2
    case(3)
!boole's rule
      call calcnodes_omega3
    case(4)
!hermite spline
      call calcnodes_omega4
    case(5)
!gauss-legendre integration
      call calcnodes_omega5
    case(6)
!gauss-chebyshev itnegration
      call calcnodes_omega6
    case(7)
!triangulation with triangles
      call calcnodes_omega7
    case(8)
!triangulation with triangles (pseudo gauss)
      call calcnodes_omega8
    case(9)
!lebedev quadrature
      call calcnodes_omega9
    case default
      stop 'error in dime_omega_nodes: opt_angint_method not specified'
!
   end select
!
!check if omega-nodes are normalized
   if(abs(1.d0-sum(weight_omega)).gt.1.d-5) then
      write(*,*) 'error in calcnodes_omega: integration weights not normalized'
      write(*,*) 'sum(weight_omega)', sum(weight_omega)
!   stop
   endif
!
!check if only positive weights occurr
   do i=1, dim_omega
      if(weight_omega(i).lt.0.d0) then
         write(*,*) 'error in calcnodes_omega: negative weights are not allowed due to stability of solution scheme'
         write(*,*) '   -> choose differen ntheta or different opt_angint_method'
         stop
      endif
   enddo
!
!set indices for nearest neighbour alo-calculations (direction dependent)
   allocate(q_alo(dim_omega,27), stat=err)
   if(err.ne.0) stop 'allocation error in calcnodes_omega'
   q_alo=0
   do i=1, dim_omega
      if(n_x(i).gt.0..and.n_y(i).gt.0..and.n_z(i).gt.0.) then
         q_alo(i,:) = (/ 1, 2, 3, 4, 5, 6, 7, 8, 9, &
            10, 11, 12, 13, 14, 15, 16, 17, 18, &
            19, 20, 21, 22, 23, 24, 25, 26, 27 /)
      elseif(n_x(i).gt.0..and.n_y(i).gt.0..and.n_z(i).lt.0.) then
         q_alo(i,:) = (/ 19, 20, 21, 22, 23, 24, 25, 26, 27, &
            10, 11, 12, 13, 14, 15, 16, 17, 18, &
            1, 2, 3, 4, 5, 6, 7, 8, 9 /)
      elseif(n_x(i).gt.0..and.n_y(i).lt.0..and.n_z(i).gt.0.) then
         q_alo(i,:) = (/ 7, 8, 9, 4, 5, 6, 1, 2, 3, &
            16, 17, 18, 13, 14, 15, 10, 11, 12, &
            25, 26, 27, 22, 23, 24, 19, 20, 21 /)
      elseif(n_x(i).lt.0..and.n_y(i).gt.0..and.n_z(i).gt.0.) then
         q_alo(i,:) = (/ 3, 2, 1, 6, 5, 4, 9, 8, 7, &
            12, 11, 10, 15, 14, 13, 18, 17, 16, &
            21, 20, 19, 24, 23, 22, 27, 26, 25 /)
      elseif(n_x(i).gt.0..and.n_y(i).lt.0..and.n_z(i).lt.0.) then
         q_alo(i,:) = (/ 25, 26, 27, 22, 23, 24, 19, 20, 21, &
            16, 17, 18, 13, 14, 15, 10, 11, 12, &
            7, 8, 9, 4, 5, 6, 1, 2, 3 /)
      elseif(n_x(i).lt.0..and.n_y(i).gt.0..and.n_z(i).lt.0.) then
         q_alo(i,:) = (/ 21, 20, 19, 24, 23, 22, 27, 26, 25, &
            12, 11, 10, 15, 14, 13, 18, 17, 16, &
            3, 2, 1, 6, 5, 4, 9, 8, 7 /)
      elseif(n_x(i).lt.0..and.n_y(i).lt.0..and.n_z(i).gt.0.) then
         q_alo(i,:) = (/ 9, 8, 7, 6, 5, 4, 3, 2, 1, &
            18, 17, 16, 15, 14, 13, 12, 11, 10, &
            27, 26, 25, 24, 23, 22, 21, 20, 19 /)
      elseif(n_x(i).lt.0..and.n_y(i).lt.0..and.n_z(i).lt.0.) then
         q_alo(i,:) = (/ 27, 26, 25, 24, 23, 22, 21, 20, 19, &
            18, 17, 16, 15, 14, 13, 12, 11, 10, &
            9, 8, 7, 6, 5, 4, 3, 2, 1 /)
      endif
   enddo
!
!
end subroutine calcnodes_omega
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine calcnodes_omega0
!
!-----------------------------------------------------------------------
!   equidistant theta-grid
!   equidistant phi-grid, with delta_phi=delta-theta
!   trapezoidal rule
!-----------------------------------------------------------------------
!
   use prog_type
   use fund_const, only: pi
   use angles, only: n_theta, dim_mu, dim_phi, dim_omega, nodes_mu, weight_mu, &
      n_x, n_y, n_z, weight_omega
   use mod_integ1d, only: precalc_oweight_trapez
!
   implicit none
!
! ... local scalars
   integer(i4b) :: i, j, k, err
   real(dp), parameter :: theta_min=1.d-6, theta_max=pi-1.d-5   !open interval, since n_x=0 or n_z=0 not allowed,
!                                                             and slightly asymmetric, to avoid mu=0.
   real(dp), parameter :: phi_min=1.d-6, phi_max=2.d0*pi-1.d-5   !open interval, n_x=0, n_y=0, n_z=0 not allowed
!                                                             !and slightly asymmetric therefore as well
!
! ... local arrays
   real(dp), dimension(:), allocatable :: nodes_phi, weight_phi

!
   dim_mu=2*n_theta-1
   dim_phi=2*dim_mu-1
   dim_omega=dim_mu*dim_phi
!
!----------------------------mu grid------------------------------------
!
   allocate(nodes_mu(dim_mu), stat=err)
   if(err.ne.0) stop 'allocation error in calcnodes_omega0'
   allocate(weight_mu(dim_mu), stat=err)
   if(err.ne.0) stop 'allocation error in calcnodes_omega0'
!
   do i=1, dim_mu
      nodes_mu(i) = theta_max - (i-1)*(theta_max-theta_min)/(dim_mu-1)
      nodes_mu(i) = cos(nodes_mu(i))
      if(abs(nodes_mu(i)).lt.1.d-15) then
         nodes_mu(i)=0.d0
      endif
   enddo
!
   call precalc_oweight_trapez(nodes_mu, dim_mu, -1.d0, 1.d0, weight_mu)
!
!normalization
   weight_mu=weight_mu/2.d0
!
!----------------------------phi grid-----------------------------------
!
   allocate(nodes_phi(dim_phi), stat=err)
   if(err.ne.0) stop 'allocation error in calcnodes_omega0'
   allocate(weight_phi(dim_phi), stat=err)
   if(err.ne.0) stop 'allocation error in calcnodes_omega0'
!
   do i=1, dim_phi
      nodes_phi(i) = phi_min + (i-1)*(phi_max-phi_min)/(dim_phi-1)
   enddo
!
   call precalc_oweight_trapez(nodes_phi, dim_phi, 0.d0, 2.d0*pi, weight_phi)
!
!normalization
   weight_phi=weight_phi/2.d0/pi
!
!------------------------all directions and weights---------------------
!
   allocate(n_x(dim_omega), stat=err)
   if(err.ne.0) stop 'allocation error in calcnodes_omega0'
   allocate(n_y(dim_omega), stat=err)
   if(err.ne.0) stop 'allocation error in calcnodes_omega0'
   allocate(n_z(dim_omega), stat=err)
   if(err.ne.0) stop 'allocation error in calcnodes_omega0'
   allocate(weight_omega(dim_omega), stat=err)
   if(err.ne.0) stop 'allocation error in calcnodes_omega0'
!
   k=1
   do i=1, dim_mu
      do j=1, dim_phi
         n_x(k)=sqrt(1.d0-nodes_mu(i)*nodes_mu(i))*cos(nodes_phi(j))
         n_y(k)=sqrt(1.d0-nodes_mu(i)*nodes_mu(i))*sin(nodes_phi(j))
         n_z(k)=nodes_mu(i)
!
         if(abs(n_x(k)).lt.1.d-14) n_x(k)=0.d0
         if(abs(n_y(k)).lt.1.d-14) n_y(k)=0.d0
         if(abs(n_z(k)).lt.1.d-14) n_z(k)=0.d0
!
         weight_omega(k) = weight_mu(i)*weight_phi(j)
         k=k+1
      enddo
   enddo
!
!
!
end subroutine calcnodes_omega0
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine calcnodes_omega1
!
!-----------------------------------------------------------------------
!   equidistant theta-grid
!   equidistant phi-grid, calculated for each theta-level
!      such that domega=const (see Lobell&Blomme 2008)
!   trapezoidal rule
!-----------------------------------------------------------------------
!
   use prog_type
   use fund_const, only: pi
   use angles, only: n_theta, dim_mu, dim_phi, dim_omega, nodes_mu, weight_mu, &
      n_x, n_y, n_z, weight_omega
   use mod_integ1d, only: precalc_oweight_trapez
!
   implicit none
!
! ... local scalars
   integer(i4b) :: i, j, k, err, n_phi
   real(dp), parameter :: theta_min=1.d-6, theta_max=pi-1.d-5   !open interval, since n_x=0 or n_z=0 not allowed,
!                                                             and slightly asymmetric, to avoid mu=0.

   real(dp), parameter :: phi_min=1.d-6, phi_max=2.d0*pi-1.d-5   !open interval, n_x=0, n_y=0, n_z=0 not allowed
!                                                             !and slightly asymmetric therefore as well

   real(dp) :: del_theta, del_phi, del_omega
!
! ... local arrays
   real(dp), dimension(:), allocatable :: nodes_phi, weight_phi, nodes_theta
   integer(i4b), dimension(:), allocatable :: dim_phi_arr
   real(dp), dimension(:,:), allocatable :: nodes_phi_pair, weight_phi_pair
   logical, dimension(:,:), allocatable :: phi_mask
!
! ... local functions
   integer(i4b) :: n_phi_lobl

!
   dim_mu=2*n_theta-1
   dim_phi=2*dim_mu-1
!dim_omega=dim_mu*dim_phi
!
!----------------------------mu grid------------------------------------
!
   allocate(nodes_mu(dim_mu), stat=err)
   if(err.ne.0) stop 'allocation error in calcnodes_omega1'
   allocate(weight_mu(dim_mu), stat=err)
   if(err.ne.0) stop 'allocation error in calcnodes_omega1'
!
   do i=1, dim_mu
      nodes_mu(i) = theta_max - (i-1)*(theta_max-theta_min)/(dim_mu-1)
      nodes_mu(i) = cos(nodes_mu(i))
      if(abs(nodes_mu(i)).lt.1.d-15) then
         nodes_mu(i)=0.d0
      endif
   enddo
!
   call precalc_oweight_trapez(nodes_mu, dim_mu, -1.d0, 1.d0, weight_mu)
!
!normalization
   weight_mu=weight_mu/2.d0
!
!----------------calculate theta-array from (given) nodes_mu-------------
!
   allocate(nodes_theta(dim_mu), stat=err)
   if(err.ne.0) stop 'allocation error calcnodes_omega1'
!
!also reverse array
   do i=1, dim_mu
      nodes_theta(i) = acos(nodes_mu(dim_mu+1-i))
   enddo
!
!-------------------calculate equidistant phi array---------------------
!----------------with different n_phi on each theta-level---------------
!
   allocate(nodes_phi_pair(dim_mu,dim_phi), stat=err)
   if(err.ne.0) stop 'allocation error in calcnodes_omega1'
   allocate(weight_phi_pair(dim_mu,dim_phi), stat=err)
   if(err.ne.0) stop 'allocation error in calcnodes_omega1'
   allocate(phi_mask(dim_mu,dim_phi), stat=err)
   if(err.ne.0) stop 'allocation error in calcnodes_omega1'
   allocate(dim_phi_arr(dim_mu), stat=err)
   if(err.ne.0) stop 'allocation error in calcnodes_omega1'
!
   phi_mask=.false.
   nodes_phi_pair=0.d0
   weight_phi_pair=0.d0
!
!del_theta, del_phi and del_omega in equatorial plane
   del_theta=nodes_theta(n_theta)-nodes_theta(n_theta-1)
   del_phi=2.d0*pi/(dim_phi-1)
   del_omega=del_theta*del_phi
!
   dim_omega=0
!
   do i=1, dim_mu
!
!-----------------calculate phi-grid for given theta--------------------
!-------new version: del_theta is allowed to be latitude dependent------
!
!calculate del_theta (of theta)
      if(i.eq.1) then
         del_theta=nodes_theta(i+1)-nodes_theta(i)
      else if(i.eq.dim_mu) then
         del_theta=nodes_theta(i)-nodes_theta(i-1)
      else
!use central differences in order to provide symmetry
         del_theta=(nodes_theta(i+1)-nodes_theta(i-1))/2.d0
      endif
!
!n_phi minimum to be 2 and has to be odd
      n_phi=n_phi_lobl(del_omega, del_theta, nodes_theta(i))
      n_phi=maxval((/n_phi,2/))
      n_phi=n_phi/2
      n_phi=2*n_phi+1
!
      if(n_phi.gt.dim_phi) stop 'error calcnodes_omega1: n_phi gt dim_phi'
!
      dim_phi_arr(i)=n_phi
!
!calculate phi grid
      do j=1, n_phi
         nodes_phi_pair(i,j)=phi_min + float(j-1)*(phi_max-phi_min)/(n_phi-1)
         phi_mask(i,j) = .true.
         dim_omega=dim_omega+1
      enddo
!
!calculate weights
      allocate(nodes_phi(n_phi), stat=err)
      if(err.ne.0) stop 'allocation error in calcnodes_omega1'
      allocate(weight_phi(n_phi), stat=err)
      if(err.ne.0) stop 'allocation error in calcnodes_omega1'
      nodes_phi=nodes_phi_pair(i,1:n_phi)
      call precalc_oweight_trapez(nodes_phi, n_phi, 0.d0, 2.d0*pi, weight_phi)
      weight_phi_pair(i,1:n_phi)=weight_phi
      deallocate(nodes_phi)
      deallocate(weight_phi)
!
   enddo
!
   weight_phi_pair=weight_phi_pair/2.d0/pi
!
!------------------------all directions and weights---------------------
!
   allocate(n_x(dim_omega), stat=err)
   if(err.ne.0) stop 'allocation error in calcnodes_omega1'
   allocate(n_y(dim_omega), stat=err)
   if(err.ne.0) stop 'allocation error in calcnodes_omega1'
   allocate(n_z(dim_omega), stat=err)
   if(err.ne.0) stop 'allocation error in calcnodes_omega1'
   allocate(weight_omega(dim_omega), stat=err)
   if(err.ne.0) stop 'allocation error in calcnodes_omega1'
!
   k=1
   do i=1, dim_mu
      do j=1, dim_phi
         if(phi_mask(i,j)) then
            n_x(k)=sqrt(1.d0-nodes_mu(i)*nodes_mu(i))*cos(nodes_phi_pair(i,j))
            n_y(k)=sqrt(1.d0-nodes_mu(i)*nodes_mu(i))*sin(nodes_phi_pair(i,j))
            n_z(k)=nodes_mu(i)

!         write(*,*) nodes_phi_pair(i,j), n_x(k)
!
            if(abs(n_x(k)).lt.1.d-14) n_x(k)=0.d0
            if(abs(n_y(k)).lt.1.d-14) n_y(k)=0.d0
            if(abs(n_z(k)).lt.1.d-14) n_z(k)=0.d0
!
            weight_omega(k) = weight_mu(i)*weight_phi_pair(i,j)
            k=k+1
         endif
      enddo
   enddo
!
!
!
end subroutine calcnodes_omega1
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine calcnodes_omega2
!
!-----------------------------------------------------------------------
!   equidistant theta-grid
!   equidistant phi-grid, with delta_phi=delta-theta
!   simpsons rule
!-----------------------------------------------------------------------
!
   use prog_type
   use fund_const, only: pi
   use angles, only: n_theta, dim_mu, dim_phi, dim_omega, nodes_mu, weight_mu, &
      n_x, n_y, n_z, weight_omega
   use mod_integ1d, only: precalc_oweight_simps
!
   implicit none
!
! ... local scalars
   integer(i4b) :: i, j, k, err
   real(dp), parameter :: theta_min=1.d-6, theta_max=pi-1.d-5   !open interval, since n_x=0 or n_z=0 not allowed,
!                                                             and slightly asymmetric, to avoid mu=0.
   real(dp), parameter :: phi_min=1.d-6, phi_max=2.d0*pi-1.d-5   !open interval, n_x=0, n_y=0, n_z=0 not allowed
!                                                             !and slightly asymmetric therefore as well
!
! ... local arrays
   real(dp), dimension(:), allocatable :: nodes_phi, weight_phi

!
   dim_mu=2*n_theta-1
   dim_phi=2*dim_mu-1
   dim_omega=dim_mu*dim_phi
!
!----------------------------mu grid------------------------------------
!
   allocate(nodes_mu(dim_mu), stat=err)
   if(err.ne.0) stop 'allocation error in calcnodes_omega2'
   allocate(weight_mu(dim_mu), stat=err)
   if(err.ne.0) stop 'allocation error in calcnodes_omega2'
!
   do i=1, dim_mu
      nodes_mu(i) = theta_max - (i-1)*(theta_max-theta_min)/(dim_mu-1)
      nodes_mu(i) = cos(nodes_mu(i))
      if(abs(nodes_mu(i)).lt.1.d-15) then
         nodes_mu(i)=0.d0
      endif
   enddo
!
   call precalc_oweight_simps(nodes_mu, dim_mu, -1.d0, 1.d0, weight_mu)
!
!normalization
   weight_mu=weight_mu/2.d0
!
!----------------------------phi grid-----------------------------------
!
   allocate(nodes_phi(dim_phi), stat=err)
   if(err.ne.0) stop 'allocation error in calcnodes_omega2'
   allocate(weight_phi(dim_phi), stat=err)
   if(err.ne.0) stop 'allocation error in calcnodes_omega2'
!
   do i=1, dim_phi
      nodes_phi(i) = phi_min + (i-1)*(phi_max-phi_min)/(dim_phi-1)
   enddo
!
   call precalc_oweight_simps(nodes_phi, dim_phi, 0.d0, 2.d0*pi, weight_phi)
!
!normalization
   weight_phi=weight_phi/2.d0/pi
!
!------------------------all directions and weights---------------------
!
   allocate(n_x(dim_omega), stat=err)
   if(err.ne.0) stop 'allocation error in calcnodes_omega2'
   allocate(n_y(dim_omega), stat=err)
   if(err.ne.0) stop 'allocation error in calcnodes_omega2'
   allocate(n_z(dim_omega), stat=err)
   if(err.ne.0) stop 'allocation error in calcnodes_omega2'
   allocate(weight_omega(dim_omega), stat=err)
   if(err.ne.0) stop 'allocation error in calcnodes_omega2'
!
   k=1
   do i=1, dim_mu
      do j=1, dim_phi
         n_x(k)=sqrt(1.d0-nodes_mu(i)*nodes_mu(i))*cos(nodes_phi(j))
         n_y(k)=sqrt(1.d0-nodes_mu(i)*nodes_mu(i))*sin(nodes_phi(j))
         n_z(k)=nodes_mu(i)
!
         if(abs(n_x(k)).lt.1.d-14) n_x(k)=0.d0
         if(abs(n_y(k)).lt.1.d-14) n_y(k)=0.d0
         if(abs(n_z(k)).lt.1.d-14) n_z(k)=0.d0
!
         weight_omega(k) = weight_mu(i)*weight_phi(j)
         k=k+1
      enddo
   enddo
!
!
!
end subroutine calcnodes_omega2
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine calcnodes_omega3
!
!-----------------------------------------------------------------------
!   equidistant theta-grid
!   equidistant phi-grid, with delta_phi=delta-theta
!   boole's rule
!-----------------------------------------------------------------------
!
   use prog_type
   use fund_const, only: pi
   use angles, only: n_theta, dim_mu, dim_phi, dim_omega, nodes_mu, weight_mu, &
      n_x, n_y, n_z, weight_omega
   use mod_integ1d, only: precalc_oweight_boole
!
   implicit none
!
! ... local scalars
   integer(i4b) :: i, j, k, err
   real(dp), parameter :: theta_min=1.d-6, theta_max=pi-1.d-5   !open interval, since n_x=0 or n_z=0 not allowed,
!                                                             and slightly asymmetric, to avoid mu=0.
   real(dp), parameter :: phi_min=1.d-6, phi_max=2.d0*pi-1.d-5   !open interval, n_x=0, n_y=0, n_z=0 not allowed
!                                                             !and slightly asymmetric therefore as well
!
! ... local arrays
   real(dp), dimension(:), allocatable :: nodes_phi, weight_phi

!
   dim_mu=2*n_theta-1
   dim_phi=2*dim_mu-1
   dim_omega=dim_mu*dim_phi
!
!----------------------------mu grid------------------------------------
!
   allocate(nodes_mu(dim_mu), stat=err)
   if(err.ne.0) stop 'allocation error in calcnodes_omega3'
   allocate(weight_mu(dim_mu), stat=err)
   if(err.ne.0) stop 'allocation error in calcnodes_omega3'
!
   do i=1, dim_mu
      nodes_mu(i) = theta_max - (i-1)*(theta_max-theta_min)/(dim_mu-1)
      nodes_mu(i) = cos(nodes_mu(i))
      if(abs(nodes_mu(i)).lt.1.d-15) then
         nodes_mu(i)=0.d0
      endif
   enddo
!
   call precalc_oweight_boole(nodes_mu, dim_mu, -1.d0, 1.d0, weight_mu)
!
!normalization
   weight_mu=weight_mu/2.d0
!
!----------------------------phi grid-----------------------------------
!
   allocate(nodes_phi(dim_phi), stat=err)
   if(err.ne.0) stop 'allocation error in calcnodes_omega3'
   allocate(weight_phi(dim_phi), stat=err)
   if(err.ne.0) stop 'allocation error in calcnodes_omega3'
!
   do i=1, dim_phi
      nodes_phi(i) = phi_min + (i-1)*(phi_max-phi_min)/(dim_phi-1)
   enddo
!
   call precalc_oweight_boole(nodes_phi, dim_phi, 0.d0, 2.d0*pi, weight_phi)
!
!normalization
   weight_phi=weight_phi/2.d0/pi
!
!------------------------all directions and weights---------------------
!
   allocate(n_x(dim_omega), stat=err)
   if(err.ne.0) stop 'allocation error in calcnodes_omega3'
   allocate(n_y(dim_omega), stat=err)
   if(err.ne.0) stop 'allocation error in calcnodes_omega3'
   allocate(n_z(dim_omega), stat=err)
   if(err.ne.0) stop 'allocation error in calcnodes_omega3'
   allocate(weight_omega(dim_omega), stat=err)
   if(err.ne.0) stop 'allocation error in calcnodes_omega3'
!
   k=1
   do i=1, dim_mu
      do j=1, dim_phi
         n_x(k)=sqrt(1.d0-nodes_mu(i)*nodes_mu(i))*cos(nodes_phi(j))
         n_y(k)=sqrt(1.d0-nodes_mu(i)*nodes_mu(i))*sin(nodes_phi(j))
         n_z(k)=nodes_mu(i)
!
         if(abs(n_x(k)).lt.1.d-14) n_x(k)=0.d0
         if(abs(n_y(k)).lt.1.d-14) n_y(k)=0.d0
         if(abs(n_z(k)).lt.1.d-14) n_z(k)=0.d0
!
         weight_omega(k) = weight_mu(i)*weight_phi(j)
         k=k+1
      enddo
   enddo
!
!
!
end subroutine calcnodes_omega3
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine calcnodes_omega4
!
!-----------------------------------------------------------------------
!   equidistant theta-grid
!   equidistant phi-grid, with delta_phi=delta-theta
!   hermite spline integration
!-----------------------------------------------------------------------
!
   use prog_type
   use fund_const, only: pi
   use angles, only: n_theta, dim_mu, dim_phi, dim_omega, nodes_mu, weight_mu, &
      n_x, n_y, n_z, weight_omega
!
   implicit none
!
! ... local scalars
   integer(i4b) :: i, j, k, err
   real(dp), parameter :: theta_min=1.d-6, theta_max=pi-1.d-5   !open interval, since n_x=0 or n_z=0 not allowed,
!                                                             and slightly asymmetric, to avoid mu=0.
   real(dp), parameter :: phi_min=1.d-6, phi_max=2.d0*pi-1.d-5   !open interval, n_x=0, n_y=0, n_z=0 not allowed
!                                                             !and slightly asymmetric therefore as well
!
! ... local arrays
   real(dp), dimension(:), allocatable :: nodes_phi, weight_phi
!
!
   dim_mu=2*n_theta-1
   dim_phi=2*dim_mu-1
   dim_omega=dim_mu*dim_phi
!
!----------------------------mu grid------------------------------------
!
   allocate(nodes_mu(dim_mu), stat=err)
   if(err.ne.0) stop 'allocation error in calcnodes_omega4'
   allocate(weight_mu(dim_mu), stat=err)
   if(err.ne.0) stop 'allocation error in calcnodes_omega4'
!
   do i=1, dim_mu
      nodes_mu(i) = theta_max - (i-1)*(theta_max-theta_min)/(dim_mu-1)
      nodes_mu(i) = cos(nodes_mu(i))
      if(abs(nodes_mu(i)).lt.1.d-15) then
         nodes_mu(i)=0.d0
      endif
   enddo
!
   stop 'precalc_weight_spline needs to be rewritten for open boundaries'
   call precalc_weight_spline(nodes_mu, dim_mu, weight_mu, .false.)
!
!normalization
   weight_mu=weight_mu/2.d0
!
!----------------------------phi grid-----------------------------------
!
   allocate(nodes_phi(dim_phi), stat=err)
   if(err.ne.0) stop 'allocation error in calcnodes_omega4'
   allocate(weight_phi(dim_phi), stat=err)
   if(err.ne.0) stop 'allocation error in calcnodes_omega4'
!
   do i=1, dim_phi
      nodes_phi(i) = phi_min + (i-1)*(phi_max-phi_min)/(dim_phi-1)
   enddo
!
   stop 'precalc_weight_spline needs to be rewritten for open boundaries'
   call precalc_weight_spline(nodes_phi, dim_phi, weight_phi, .true.)
!
!normalization
   weight_phi=weight_phi/2.d0/pi
!
!------------------------all directions and weights---------------------
!
   allocate(n_x(dim_omega), stat=err)
   if(err.ne.0) stop 'allocation error in calcnodes_omega4'
   allocate(n_y(dim_omega), stat=err)
   if(err.ne.0) stop 'allocation error in calcnodes_omega4'
   allocate(n_z(dim_omega), stat=err)
   if(err.ne.0) stop 'allocation error in calcnodes_omega4'
   allocate(weight_omega(dim_omega), stat=err)
   if(err.ne.0) stop 'allocation error in calcnodes_omega4'
!
   k=1
   do i=1, dim_mu
      do j=1, dim_phi
         n_x(k)=sqrt(1.d0-nodes_mu(i)*nodes_mu(i))*cos(nodes_phi(j))
         n_y(k)=sqrt(1.d0-nodes_mu(i)*nodes_mu(i))*sin(nodes_phi(j))
         n_z(k)=nodes_mu(i)
!
         if(abs(n_x(k)).lt.1.d-14) n_x(k)=0.d0
         if(abs(n_y(k)).lt.1.d-14) n_y(k)=0.d0
         if(abs(n_z(k)).lt.1.d-14) n_z(k)=0.d0
!
         weight_omega(k) = weight_mu(i)*weight_phi(j)
         k=k+1
      enddo
   enddo
!
!
!
end subroutine calcnodes_omega4
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine calcnodes_omega5
!
!-----------------------------------------------------------------------
!   gauss legendre integration with n_theta points in [0.d0,pi/2.d0]
!      (i.e., same integration in each octant)
!-----------------------------------------------------------------------
!
   use prog_type
   use fund_const, only: pi
   use angles, only: n_theta, dim_mu, dim_phi, dim_omega, nodes_mu, weight_mu, &
      n_x, n_y, n_z, weight_omega
   use mod_integ1d, only: precalc_weight_legendre
!
   implicit none
!
! ... local scalars
   integer(i4b) :: i, j, k, err
   real(dp) :: theta_min, theta_max, phi_min, phi_max, sum
!
! ... local arrays
   real(dp), dimension(:), allocatable :: nodes_dum, weight_dum, nodes_phi, weight_phi, &
      nodes_rev, weight_rev
!
   dim_mu=2*n_theta
   dim_phi=4*n_theta
   dim_omega=dim_mu*dim_phi
!
!---------------------standard grid in [-1,d0,1.d0]---------------------
!
   allocate(nodes_dum(n_theta), stat=err)
   if(err.ne.0) stop 'allocation error in calcnodes_omega5'
   allocate(weight_dum(n_theta), stat=err)
   if(err.ne.0) stop 'allocation error in calcnodes_omega5'
!
   call precalc_weight_legendre(n_theta, -1.d0, 1.d0, nodes_dum, weight_dum)
!
!-----------------------------mu grid-----------------------------------
!
   allocate(nodes_mu(dim_mu), stat=err)
   if(err.ne.0) stop 'allocation error in calcnodes_omega5'
   allocate(weight_mu(dim_mu), stat=err)
   if(err.ne.0) stop 'allocation error in calcnodes_omega5'
!
   theta_min=0.d0
   theta_max=pi/2.d0
   nodes_mu(1:n_theta)=(theta_max-theta_min)/2.d0 * nodes_dum + (theta_min+theta_max)/2.d0
   weight_mu(1:n_theta)=weight_dum*(theta_max-theta_min)/2.d0*sin(nodes_mu(1:n_theta))
!
   theta_min=pi/2.d0
   theta_max=pi
   nodes_mu(n_theta+1:dim_mu) = (theta_max-theta_min)/2.d0 * nodes_dum + (theta_min+theta_max)/2.d0
   weight_mu(n_theta+1:dim_mu)=weight_dum*sin(nodes_mu(n_theta+1:dim_mu))*(theta_max-theta_min)/2.d0
!
   nodes_mu=cos(nodes_mu)
!
   weight_mu = weight_mu/2.d0
!
   allocate(nodes_rev(dim_mu), stat=err)
   if(err.ne.0) stop 'allocation error calcnodes_omega5'
   allocate(weight_rev(dim_mu), stat=err)
   if(err.ne.0) stop 'allocation error calcnodes_omega5'
!
   do i=1, dim_mu
      nodes_rev(i) = nodes_mu(dim_mu+1-i)
      weight_rev(i) = weight_mu(dim_mu+1-i)
   enddo
!
   nodes_mu=nodes_rev
   weight_mu=weight_rev
!
!-------------------------------phi grid--------------------------------
!
   allocate(nodes_phi(dim_phi), stat=err)
   if(err.ne.0) stop 'allocation error in calcnodes_omega5'
   allocate(weight_phi(dim_phi), stat=err)
   if(err.ne.0) stop 'allocation error in calcnodes_omega5'
!
   phi_min=0.d0
   phi_max=pi/2.d0
   nodes_phi(1:n_theta) = (phi_max-phi_min)/2.d0 * nodes_dum + (phi_min+phi_max)/2.d0
   weight_phi(1:n_theta) = weight_dum*(phi_max-phi_min)/2.d0
!
   phi_min=pi/2.d0
   phi_max=pi
   nodes_phi(n_theta+1:2*n_theta) = (phi_max-phi_min)/2.d0 * nodes_dum + (phi_min+phi_max)/2.d0
   weight_phi(n_theta+1:2*n_theta) = weight_dum*(phi_max-phi_min)/2.d0
!
   phi_min=pi
   phi_max=3.d0*pi/2.d0
   nodes_phi(2*n_theta+1:3*n_theta) = (phi_max-phi_min)/2.d0 * nodes_dum + (phi_min+phi_max)/2.d0
   weight_phi(2*n_theta+1:3*n_theta) = weight_dum*(phi_max-phi_min)/2.d0
!
   phi_min=3.d0*pi/2.d0
   phi_max=2.d0*pi
   nodes_phi(3*n_theta+1:4*n_theta) = (phi_max-phi_min)/2.d0 * nodes_dum + (phi_min+phi_max)/2.d0
   weight_phi(3*n_theta+1:4*n_theta) = weight_dum*(phi_max-phi_min)/2.d0
!
   weight_phi=weight_phi/2.d0/pi
!
!------------------------all directions and weights---------------------
!
   allocate(n_x(dim_omega), stat=err)
   if(err.ne.0) stop 'allocation error in calcnodes_omega5'
   allocate(n_y(dim_omega), stat=err)
   if(err.ne.0) stop 'allocation error in calcnodes_omega5'
   allocate(n_z(dim_omega), stat=err)
   if(err.ne.0) stop 'allocation error in calcnodes_omega5'
   allocate(weight_omega(dim_omega), stat=err)
   if(err.ne.0) stop 'allocation error in calcnodes_omega5'
!
   k=1
   do i=1, dim_mu
      do j=1, dim_phi
         n_x(k)=sqrt(1.d0-nodes_mu(i)*nodes_mu(i))*cos(nodes_phi(j))
         n_y(k)=sqrt(1.d0-nodes_mu(i)*nodes_mu(i))*sin(nodes_phi(j))
         n_z(k)=nodes_mu(i)
!
         if(abs(n_x(k)).lt.1.d-14) n_x(k)=0.d0
         if(abs(n_y(k)).lt.1.d-14) n_y(k)=0.d0
         if(abs(n_z(k)).lt.1.d-14) n_z(k)=0.d0
!
         weight_omega(k) = weight_mu(i)*weight_phi(j)
         k=k+1
      enddo
   enddo
!
!
!
end subroutine calcnodes_omega5
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine calcnodes_omega6
!
!-----------------------------------------------------------------------
!   gauss chebyshev integration with n_theta points in [0.d0,pi/2.d0]
!      (i.e., same integration in each octant)
!
! NOTE: EVENTUALLY TO BE DEBUGGED BECAUSE WEIGHTS DO NOT SUM TO 1!!!!
!-----------------------------------------------------------------------
!
   use prog_type
   use fund_const, only: pi
   use angles, only: n_theta, dim_mu, dim_phi, dim_omega, nodes_mu, weight_mu, &
      n_x, n_y, n_z, weight_omega
   use mod_integ1d, only: precalc_weight_chebyshev
!
   implicit none
!
! ... local scalars
   integer(i4b) :: i, j, k, err
   real(dp) :: theta_min, theta_max, phi_min, phi_max, sum
!
! ... local arrays
   real(dp), dimension(:), allocatable :: nodes_dum, weight_dum, nodes_phi, weight_phi, &
      nodes_rev, weight_rev
!
   dim_mu=2*n_theta
   dim_phi=4*n_theta
   dim_omega=dim_mu*dim_phi
!
!---------------------standard grid in [-1,d0,1.d0]---------------------
!
   allocate(nodes_dum(n_theta), stat=err)
   if(err.ne.0) stop 'allocation error in calcnodes_omega6'
   allocate(weight_dum(n_theta), stat=err)
   if(err.ne.0) stop 'allocation error in calcnodes_omega6'
!
   call precalc_weight_chebyshev(n_theta, -1.d0, 1.d0, nodes_dum, weight_dum)
!
!-----------------------------mu grid-----------------------------------
!
   allocate(nodes_mu(dim_mu), stat=err)
   if(err.ne.0) stop 'allocation error in calcnodes_omega6'
   allocate(weight_mu(dim_mu), stat=err)
   if(err.ne.0) stop 'allocation error in calcnodes_omega6'
!
   theta_min=0.d0
   theta_max=pi/2.d0
   nodes_mu(1:n_theta)=(theta_max-theta_min)/2.d0 * nodes_dum + (theta_min+theta_max)/2.d0
   weight_mu(1:n_theta)=weight_dum*(theta_max-theta_min)/2.d0*sin(nodes_mu(1:n_theta))
!
   theta_min=pi/2.d0
   theta_max=pi
   nodes_mu(n_theta+1:dim_mu) = (theta_max-theta_min)/2.d0 * nodes_dum + (theta_min+theta_max)/2.d0
   weight_mu(n_theta+1:dim_mu)=weight_dum*sin(nodes_mu(n_theta+1:dim_mu))*(theta_max-theta_min)/2.d0
!
   nodes_mu=cos(nodes_mu)
!
   weight_mu = weight_mu/2.d0
!
   allocate(nodes_rev(dim_mu), stat=err)
   if(err.ne.0) stop 'allocation error calcnodes_omega6'
   allocate(weight_rev(dim_mu), stat=err)
   if(err.ne.0) stop 'allocation error calcnodes_omega6'
!
   do i=1, dim_mu
      nodes_rev(i) = nodes_mu(dim_mu+1-i)
      weight_rev(i) = weight_mu(dim_mu+1-i)
   enddo
!
   nodes_mu=nodes_rev
   weight_mu=weight_rev
!
!-------------------------------phi grid--------------------------------
!
   allocate(nodes_phi(dim_phi), stat=err)
   if(err.ne.0) stop 'allocation error in calcnodes_omega6'
   allocate(weight_phi(dim_phi), stat=err)
   if(err.ne.0) stop 'allocation error in calcnodes_omega6'
!
   phi_min=0.d0
   phi_max=pi/2.d0
   nodes_phi(1:n_theta) = (phi_max-phi_min)/2.d0 * nodes_dum + (phi_min+phi_max)/2.d0
   weight_phi(1:n_theta) = weight_dum*(phi_max-phi_min)/2.d0
!
   phi_min=pi/2.d0
   phi_max=pi
   nodes_phi(n_theta+1:2*n_theta) = (phi_max-phi_min)/2.d0 * nodes_dum + (phi_min+phi_max)/2.d0
   weight_phi(n_theta+1:2*n_theta) = weight_dum*(phi_max-phi_min)/2.d0
!
   phi_min=pi
   phi_max=3.d0*pi/2.d0
   nodes_phi(2*n_theta+1:3*n_theta) = (phi_max-phi_min)/2.d0 * nodes_dum + (phi_min+phi_max)/2.d0
   weight_phi(2*n_theta+1:3*n_theta) = weight_dum*(phi_max-phi_min)/2.d0
!
   phi_min=3.d0*pi/2.d0
   phi_max=2.d0*pi
   nodes_phi(3*n_theta+1:4*n_theta) = (phi_max-phi_min)/2.d0 * nodes_dum + (phi_min+phi_max)/2.d0
   weight_phi(3*n_theta+1:4*n_theta) = weight_dum*(phi_max-phi_min)/2.d0
!
   weight_phi=weight_phi/2.d0/pi
!
!------------------------all directions and weights---------------------
!
   allocate(n_x(dim_omega), stat=err)
   if(err.ne.0) stop 'allocation error in calcnodes_omega6'
   allocate(n_y(dim_omega), stat=err)
   if(err.ne.0) stop 'allocation error in calcnodes_omega6'
   allocate(n_z(dim_omega), stat=err)
   if(err.ne.0) stop 'allocation error in calcnodes_omega6'
   allocate(weight_omega(dim_omega), stat=err)
   if(err.ne.0) stop 'allocation error in calcnodes_omega6'
!
   k=1
   do i=1, dim_mu
      do j=1, dim_phi
         n_x(k)=sqrt(1.d0-nodes_mu(i)*nodes_mu(i))*cos(nodes_phi(j))
         n_y(k)=sqrt(1.d0-nodes_mu(i)*nodes_mu(i))*sin(nodes_phi(j))
         n_z(k)=nodes_mu(i)
!
         if(abs(n_x(k)).lt.1.d-14) n_x(k)=0.d0
         if(abs(n_y(k)).lt.1.d-14) n_y(k)=0.d0
         if(abs(n_z(k)).lt.1.d-14) n_z(k)=0.d0
!
         weight_omega(k) = weight_mu(i)*weight_phi(j)
         k=k+1
      enddo
   enddo
!
!
end subroutine calcnodes_omega6
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine calcnodes_omega7
!
!
!-----------------------------------------------------------------------
!   split 2d integration domain in triangles
!   integration with barycentric interpolation in each triangle
!-----------------------------------------------------------------------
!
   use prog_type
   use fund_const, only: pi
   use angles, only: dim_mu, dim_phi, dim_omega, n_theta, nodes_mu, weight_omega, n_x, n_y, n_z
!
   implicit none
!
! ... local scalars
   integer(i4b) :: i, j, k, err
   real(dp) :: del_im, del_ip, phi1, phi2, mu1, mu2, w1, w2, w3, sum
   real(dp), parameter :: phi_min=1.d-6, phi_max=2.d0*pi-1.d-5   !open interval, n_x=0, n_y=0, n_z=0 not allowed
!                                                             !and slightly asymmetric therefore as well
   real(dp), parameter :: theta_min=1.d-6, theta_max=pi-1.d-5   !open interval, since n_x=0 or n_z=0 not allowed,
!                                                             and slightly asymmetric, to avoid mu=0.
!
! ... local arrays
   real(dp), dimension(:), allocatable :: nodes_phi
!
!
   dim_mu=2*n_theta-1
   dim_phi=2*dim_mu-1
!
   if(dim_mu.lt.3) stop 'error in calcnodes_omega7: dim_mu has to be ge 3'
!if(mod(dim_mu,2).eq.0) stop 'error in calcnodes_omega7: dim_mu has to be odd'
   if(dim_phi.lt.3) stop 'error in calcnodes_omega7:: dim_phi has to be ge 3'
!if(mod(dim_phi,2).eq.0) stop 'error in calcnodes_omega7: dim_phi has to be odd'
!
!----------------------------mu grid------------------------------------
!
   allocate(nodes_mu(dim_mu), stat=err)
   if(err.ne.0) stop 'allocation error in calcnodes_omega7'
!
   do i=1, dim_mu
      nodes_mu(i) = theta_max - (i-1)*(theta_max-theta_min)/(dim_mu-1)
   enddo
   nodes_mu=cos(nodes_mu)
!
!make mu-grid equidistant for three subsequent nodes
   do i=2, dim_mu, 2
      nodes_mu(i) = (nodes_mu(i+1)+nodes_mu(i-1))/2.d0
   enddo
!
!check if three subsequent x-points are really equidistant
   do i=2, dim_mu-1, 2
      del_im=nodes_mu(i)-nodes_mu(i-1)
      del_ip=nodes_mu(i+1)-nodes_mu(i)
      if(abs(del_im-del_ip).gt.1.d-15) then
         stop 'error in calcnodes_omega7: mu-grid not equidistant'
      endif
   enddo
!
!----------------------------phi grid-----------------------------------
!
   allocate(nodes_phi(dim_phi), stat=err)
   if(err.ne.0) stop 'allocation error in calcnodes_omega7'
!
   do i=1, dim_phi
      nodes_phi(i) = phi_min + (i-1)*(phi_max-phi_min)/(dim_phi-1)
   enddo
!
   dim_omega=ceiling(dim_mu/2.)*ceiling(dim_phi/2.) + floor(dim_mu/2.)*floor(dim_phi/2.)
!
!------------------------all directions and weights---------------------
!
   allocate(n_x(dim_omega), stat=err)
   if(err.ne.0) stop 'allocation error in calcnodes_omega7'
   allocate(n_y(dim_omega), stat=err)
   if(err.ne.0) stop 'allocation error in calcnodes_omega7'
   allocate(n_z(dim_omega), stat=err)
   if(err.ne.0) stop 'allocation error in calcnodes_omega7'
   allocate(weight_omega(dim_omega), stat=err)
   if(err.ne.0) stop 'allocation error in calcnodes_omega7'
!
   weight_omega=0.d0
!
   phi1=0.d0
   phi2=2.d0*pi
   mu1=-1.d0
   mu2=1.d0
!
   k=1
   do i=1, dim_mu
      do j=1, dim_phi
         if(i.eq.1.and.j.eq.1) then   !weight at (1,1)
            w1=((nodes_mu(i+1)-nodes_mu(i))*(nodes_phi(j+2)-nodes_phi(j)) + &
               (nodes_mu(i+2)-nodes_mu(i))*(nodes_phi(j+1)-nodes_phi(j)))/6.d0
            w2=((nodes_mu(i+2)-nodes_mu(i))*(nodes_phi(j)-phi1) + &
               (nodes_mu(i)-mu1)*(nodes_phi(j+2)-nodes_phi(j)))/2.d0
            w3=(nodes_mu(i)-mu1)*(nodes_phi(j)-phi1)
            weight_omega(k)=w1+w2+w3
            n_x(k)=sqrt(1.d0-nodes_mu(i)*nodes_mu(i))*cos(nodes_phi(j))
            n_y(k)=sqrt(1.d0-nodes_mu(i)*nodes_mu(i))*sin(nodes_phi(j))
            n_z(k)=nodes_mu(i)
            k=k+1
         elseif(i.eq.dim_mu.and.j.eq.1) then   !weight at (dim_mu,1)
            w1=((nodes_mu(i)-nodes_mu(i-2))*(nodes_phi(j+1)-nodes_phi(j)) + &
               (nodes_mu(i)-nodes_mu(i-1))*(nodes_phi(j+2)-nodes_phi(j)))/6.d0
            w2=((nodes_mu(i)-nodes_mu(i-2))*(nodes_phi(j)-phi1) + &
               (mu2-nodes_mu(i))*(nodes_phi(j+2)-nodes_phi(j)))/2.d0
            w3=(mu2-nodes_mu(i))*(nodes_phi(j)-phi1)
            weight_omega(k)=w1+w2+w3
            n_x(k)=sqrt(1.d0-nodes_mu(i)*nodes_mu(i))*cos(nodes_phi(j))
            n_y(k)=sqrt(1.d0-nodes_mu(i)*nodes_mu(i))*sin(nodes_phi(j))
            n_z(k)=nodes_mu(i)
            k=k+1
         elseif(i.eq.1.and.j.eq.dim_phi) then   !weight at (1,dim_phi)
            w1=((nodes_mu(i+1)-nodes_mu(i))*(nodes_phi(j)-nodes_phi(j-2)) + &
               (nodes_mu(i+2)-nodes_mu(i))*(nodes_phi(j)-nodes_phi(j-1)))/6.d0
            w2=((nodes_mu(i+2)-nodes_mu(i))*(phi2-nodes_phi(j)) + &
               (nodes_mu(i)-mu1)*(nodes_phi(j)-nodes_phi(j-2)))/2.d0
            w3=(nodes_mu(i)-mu1)*(phi2-nodes_phi(j))
            weight_omega(k)=w1+w2+w3
            n_x(k)=sqrt(1.d0-nodes_mu(i)*nodes_mu(i))*cos(nodes_phi(j))
            n_y(k)=sqrt(1.d0-nodes_mu(i)*nodes_mu(i))*sin(nodes_phi(j))
            n_z(k)=nodes_mu(i)
            k=k+1
         elseif(i.eq.dim_mu.and.j.eq.dim_phi) then   !weight at (dim_mu,dim_phi)
            w1=((nodes_mu(i)-nodes_mu(i-2))*(nodes_phi(j)-nodes_phi(j-1)) + &
               (nodes_mu(i)-nodes_mu(i-1))*(nodes_phi(j)-nodes_phi(j-2)))/6.d0
            w2=((nodes_mu(i)-nodes_mu(i-2))*(phi2-nodes_phi(j)) + &
               (mu2-nodes_mu(i))*(nodes_phi(j)-nodes_phi(j-2)))/2.d0
            w3=(mu2-nodes_mu(i))*(phi2-nodes_phi(j))
            weight_omega(k)=w1+w2+w3
            n_x(k)=sqrt(1.d0-nodes_mu(i)*nodes_mu(i))*cos(nodes_phi(j))
            n_y(k)=sqrt(1.d0-nodes_mu(i)*nodes_mu(i))*sin(nodes_phi(j))
            n_z(k)=nodes_mu(i)
            k=k+1
         elseif(j.eq.1.and.mod(i,2).ne.0) then   !weights along edge at j=1 (only odd indices contribute)
            w1=((nodes_mu(i+2)-nodes_mu(i-2))*(nodes_phi(j+1)-nodes_phi(j)) + &
               (nodes_mu(i+1)-nodes_mu(i-1))*(nodes_phi(j+2)-nodes_phi(j)))/6.d0
            w2=(nodes_mu(i+2)-nodes_mu(i-2))*(nodes_phi(j)-phi1)/2.d0
            weight_omega(k)=w1+w2
            n_x(k)=sqrt(1.d0-nodes_mu(i)*nodes_mu(i))*cos(nodes_phi(j))
            n_y(k)=sqrt(1.d0-nodes_mu(i)*nodes_mu(i))*sin(nodes_phi(j))
            n_z(k)=nodes_mu(i)
            k=k+1
         elseif(j.eq.dim_phi.and.mod(i,2).ne.0) then   !weights along edge at j=dim_phi (only odd indices contribute)
            w1=((nodes_mu(i+2)-nodes_mu(i-2))*(nodes_phi(j)-nodes_phi(j-1)) + &
               (nodes_mu(i+1)-nodes_mu(i-1))*(nodes_phi(j)-nodes_phi(j-2)))/6.d0
            w2=(nodes_mu(i+2)-nodes_mu(i-2))*(phi2-nodes_phi(j))/2.d0
            weight_omega(k)=w1+w2
            n_x(k)=sqrt(1.d0-nodes_mu(i)*nodes_mu(i))*cos(nodes_phi(j))
            n_y(k)=sqrt(1.d0-nodes_mu(i)*nodes_mu(i))*sin(nodes_phi(j))
            n_z(k)=nodes_mu(i)
            k=k+1
         elseif(i.eq.1.and.mod(j,2).ne.0) then   !weights along edge at i=1 (only odd indices contribute)
            w1=((nodes_mu(i+1)-nodes_mu(i))*(nodes_phi(j+2)-nodes_phi(j-2)) + &
               (nodes_mu(i+2)-nodes_mu(i))*(nodes_phi(j+1)-nodes_phi(j-1)))/6.d0
            w2=(nodes_mu(i)-mu1)*(nodes_phi(j+2)-nodes_phi(j-2))/2.d0
            weight_omega(k)=w1+w2
            n_x(k)=sqrt(1.d0-nodes_mu(i)*nodes_mu(i))*cos(nodes_phi(j))
            n_y(k)=sqrt(1.d0-nodes_mu(i)*nodes_mu(i))*sin(nodes_phi(j))
            n_z(k)=nodes_mu(i)
            k=k+1
         elseif(i.eq.dim_mu.and.mod(j,2).ne.0) then   !weights along edge at i=dim_mu (only odd indices contribute)
            w1=((nodes_mu(i)-nodes_mu(i-1))*(nodes_phi(j+2)-nodes_phi(j-2)) + &
               (nodes_mu(i)-nodes_mu(i-2))*(nodes_phi(j+1)-nodes_phi(j-1)))/6.d0
            w2=(mu2-nodes_mu(i))*(nodes_phi(j+2)-nodes_phi(j-2))/2.d0
            weight_omega(k)=w1+w2
            n_x(k)=sqrt(1.d0-nodes_mu(i)*nodes_mu(i))*cos(nodes_phi(j))
            n_y(k)=sqrt(1.d0-nodes_mu(i)*nodes_mu(i))*sin(nodes_phi(j))
            n_z(k)=nodes_mu(i)
            k=k+1
         elseif(mod(i,2).eq.0.and.mod(j,2).eq.0) then   !weights at each rectangle center (only even indices)
            weight_omega(k)=(nodes_mu(i+1)-nodes_mu(i-1))*(nodes_phi(j+1)-nodes_phi(j-1))/3.d0
            n_x(k)=sqrt(1.d0-nodes_mu(i)*nodes_mu(i))*cos(nodes_phi(j))
            n_y(k)=sqrt(1.d0-nodes_mu(i)*nodes_mu(i))*sin(nodes_phi(j))
            n_z(k)=nodes_mu(i)
            k=k+1
         elseif(mod(i,2).ne.0.and.mod(j,2).ne.0) then   !weights at each rectangle corner (only odd indices)
            w1=(nodes_mu(i+2)-nodes_mu(i-2))*(nodes_phi(j+1)-nodes_phi(j-1))
            w2=(nodes_mu(i+1)-nodes_mu(i-1))*(nodes_phi(j+2)-nodes_phi(j-2))
            weight_omega(k)=(w1+w2)/6.d0
            n_x(k)=sqrt(1.d0-nodes_mu(i)*nodes_mu(i))*cos(nodes_phi(j))
            n_y(k)=sqrt(1.d0-nodes_mu(i)*nodes_mu(i))*sin(nodes_phi(j))
            n_z(k)=nodes_mu(i)
            k=k+1
         endif
      enddo
   enddo
!
!normalize weights
   weight_omega=weight_omega/4.d0/pi
!
!
end subroutine calcnodes_omega7
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine calcnodes_omega8
!
!-----------------------------------------------------------------------
!   split 2d integration domain in triangles
!   integration with 7 point 'pseudo gauss' in each triangle
!-----------------------------------------------------------------------
!
   use prog_type
   use fund_const, only: pi
   use angles, only: dim_mu, dim_phi, dim_omega, n_theta, nodes_mu, weight_omega, n_x, n_y, n_z
!
   implicit none
!
! ... local scalars
   integer(i4b) :: i, j, k, kx, ky, err
   integer(i4b) :: dim_mu2, dim_phi2
   real(dp) :: jaca, jacb, jacc, jacd, w1, w2
!
! ... local arrays
   real(dp), dimension(:), allocatable :: nodes_mu2, nodes_phi2, nodes_phi
   real(dp), dimension(:,:), allocatable :: weight_omega2
   logical, dimension(:,:), allocatable :: phi_mask
!
   dim_mu=7*n_theta
   dim_phi=2*dim_mu
!
   if(dim_mu.lt.7) stop 'error in calcnodes_omega8: dim_mu has to be ge 7'
   if(dim_phi.lt.7) stop 'error in calcnodes_omega8: dim_phi has to be ge 7'
!
   dim_mu2=dim_mu/7 + 1
   dim_phi2=dim_phi/7 + 1
!
!actual grid, for which gauss points are distributed
   allocate(nodes_mu2(dim_mu2), stat=err)
   if(err.ne.0) stop 'allocation error in calcnodes_omega8'
   allocate(nodes_phi2(dim_phi2), stat=err)
   if(err.ne.0) stop 'allocation error in calcnodes_omega8'
!
!-------------------------equidistant theta grid-------------------------
!
   do i=1, dim_mu2
      nodes_mu2(i) = pi - (i-1)*pi/(dim_mu2-1)
   enddo
   nodes_mu2=cos(nodes_mu2)
!
!-------------------------equidistant phi grid--------------------------
!
   do i=1, dim_phi2
      nodes_phi2(i) = (i-1)*2.d0*pi/(dim_phi2-1)
   enddo
!
!-----------------------distribute gauss points-------------------------
!
   allocate(nodes_mu(dim_mu), stat=err)
   if(err.ne.0) stop 'allocation error in calcnodes_omega8'
   allocate(nodes_phi(dim_phi), stat=err)
   if(err.ne.0) stop 'allocation error in calcnodes_omega8'
!
   nodes_mu=0.d0
   k=1
   do i=2, dim_mu2
      nodes_mu(k)   = (9*nodes_mu2(i-1)+nodes_mu2(i))/10.d0
      nodes_mu(k+1) = (5*nodes_mu2(i-1)+nodes_mu2(i))/6.d0
      nodes_mu(k+2) = (7*nodes_mu2(i-1)+3*nodes_mu2(i))/10.d0
      nodes_mu(k+3) = (nodes_mu2(i-1)+nodes_mu2(i))/2.d0
      nodes_mu(k+4) = (7*nodes_mu2(i)+3*nodes_mu2(i-1))/10.d0
      nodes_mu(k+5) = (5*nodes_mu2(i)+nodes_mu2(i-1))/6.d0
      nodes_mu(k+6) = (9*nodes_mu2(i)+nodes_mu2(i-1))/10.d0
      k=k+7
   enddo
!
   k=1
   do i=2, dim_phi2
      nodes_phi(k)   = (9*nodes_phi2(i-1)+nodes_phi2(i))/10.d0
      nodes_phi(k+1) = (5*nodes_phi2(i-1)+nodes_phi2(i))/6.d0
      nodes_phi(k+2) = (7*nodes_phi2(i-1)+3*nodes_phi2(i))/10.d0
      nodes_phi(k+3) = (nodes_phi2(i-1)+nodes_phi2(i))/2.d0
      nodes_phi(k+4) = (7*nodes_phi2(i)+3*nodes_phi2(i-1))/10.d0
      nodes_phi(k+5) = (5*nodes_phi2(i)+nodes_phi2(i-1))/6.d0
      nodes_phi(k+6) = (9*nodes_phi2(i)+nodes_phi2(i-1))/10.d0
      k=k+7
   enddo
!
!--------------------calculate weights and mask-------------------------
!
   allocate(phi_mask(dim_mu,dim_phi), stat=err)
   if(err.ne.0) stop 'allocation error in calcnodes_omega8'
   allocate(weight_omega2(dim_mu,dim_phi), stat=err)
   if(err.ne.0) stop 'allocation error in calcnodes_omega8'
!
   phi_mask=.false.
   weight_omega2=0.d0
!
!standard weights (on unit triangle)
   w1=25.d0/48.d0/2.d0
   w2=-27.d0/48.d0/2.d0
!
   dim_omega=0
!
   ky=1
   do j=2, dim_phi2
      kx=1
      do i=2, dim_mu2
!jacobians for each triangle
         jaca=(nodes_mu2(i)-nodes_mu2(i-1))*(nodes_phi2(j)-nodes_phi2(j-1))/2.d0
         jacb=(nodes_mu2(i)-nodes_mu2(i-1))*(nodes_phi2(j)-nodes_phi2(j-1))/2.d0
         jacc=(nodes_mu2(i)-nodes_mu2(i-1))*(nodes_phi2(j)-nodes_phi2(j-1))/2.d0
         jacd=(nodes_mu2(i)-nodes_mu2(i-1))*(nodes_phi2(j)-nodes_phi2(j-1))/2.d0
!triangle a
         weight_omega2(kx+2,ky)=w1*jaca
         weight_omega2(kx+4,ky)=w1*jaca
         weight_omega2(kx+3,ky+1)=w2*jaca
         weight_omega2(kx+3,ky+2)=w1*jaca
         phi_mask(kx+2,ky)=.true.
         phi_mask(kx+4,ky)=.true.
         phi_mask(kx+3,ky+1)=.true.
         phi_mask(kx+3,ky+2)=.true.
!triangle b
         weight_omega2(kx+6,ky+2)=w1*jacb
         weight_omega2(kx+4,ky+3)=w1*jacb
         weight_omega2(kx+5,ky+3)=w2*jacb
         weight_omega2(kx+6,ky+4)=w1*jacb
         phi_mask(kx+6,ky+2)=.true.
         phi_mask(kx+4,ky+3)=.true.
         phi_mask(kx+5,ky+3)=.true.
         phi_mask(kx+6,ky+4)=.true.
!triangle c
         weight_omega2(kx+3,ky+4)=w1*jacc
         weight_omega2(kx+3,ky+5)=w2*jacc
         weight_omega2(kx+2,ky+6)=w1*jacc
         weight_omega2(kx+4,ky+6)=w1*jacc
         phi_mask(kx+3,ky+4)=.true.
         phi_mask(kx+3,ky+5)=.true.
         phi_mask(kx+2,ky+6)=.true.
         phi_mask(kx+4,ky+6)=.true.
!triangle c
         weight_omega2(kx,ky+2)=w1*jacd
         weight_omega2(kx+1,ky+3)=w2*jacd
         weight_omega2(kx+2,ky+3)=w1*jacd
         weight_omega2(kx,ky+4)=w1*jacd
         phi_mask(kx,ky+2)=.true.
         phi_mask(kx+1,ky+3)=.true.
         phi_mask(kx+2,ky+3)=.true.
         phi_mask(kx,ky+4)=.true.
         dim_omega=dim_omega+16
         kx=kx+7
      enddo
      ky=ky+7
   enddo

!
!------------------------all directions and weights---------------------
!
   allocate(n_x(dim_omega), stat=err)
   if(err.ne.0) stop 'allocation error in calcnodes_omega0'
   allocate(n_y(dim_omega), stat=err)
   if(err.ne.0) stop 'allocation error in calcnodes_omega0'
   allocate(n_z(dim_omega), stat=err)
   if(err.ne.0) stop 'allocation error in calcnodes_omega0'
   allocate(weight_omega(dim_omega), stat=err)
   if(err.ne.0) stop 'allocation error in calcnodes_omega0'
!
   k=1
   do i=1, dim_mu
      do j=1, dim_phi
         if(phi_mask(i,j)) then
            n_x(k)=sqrt(1.d0-nodes_mu(i)*nodes_mu(i))*cos(nodes_phi(j))
            n_y(k)=sqrt(1.d0-nodes_mu(i)*nodes_mu(i))*sin(nodes_phi(j))
            n_z(k)=nodes_mu(i)
!
            if(abs(n_x(k)).lt.1.d-14) n_x(k)=0.d0
            if(abs(n_y(k)).lt.1.d-14) n_y(k)=0.d0
            if(abs(n_z(k)).lt.1.d-14) n_z(k)=0.d0
!
            weight_omega(k) = weight_omega2(i,j)
            k=k+1
         endif
      enddo
   enddo
!
!normalize weights
   weight_omega=weight_omega/4.d0/pi
!
!write(*,'(10es20.8)') weight_omega
!stop
!
end subroutine calcnodes_omega8
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine calcnodes_omega9
!
!-----------------------------------------------------------------------
!   lebedev quadrature with routines from
!   https://people.sc.fsu.edu/~jburkardt/f_src/sphere_lebedev_rule/sphere_lebedev_rule.html
!-----------------------------------------------------------------------
!
   use prog_type
   use fund_const, only: pi
   use angles, only: n_theta, dim_mu, dim_phi, dim_omega, nodes_mu, weight_mu, &
      n_x, n_y, n_z, weight_omega
!
   implicit none
!
! ... local scalars
   integer(i4b) :: i, err, available, n_omega
   real(dp) :: cb, sb, ca, sa, cc, sc, n_x2, n_y2, n_z2
!
! ... local arrays
!
! ... local functions
   integer(i4b) :: order_table, available_table
!
!only rough estimate for number of integration points
   dim_mu=2*n_theta-1
   dim_phi=2*dim_mu-1
   n_omega=dim_mu*dim_phi
!
!allocate nodes_mu and weight mu (although not required, but for output-routines...)
   allocate(nodes_mu(dim_mu), stat=err)
   if(err.ne.0) stop 'allocation error in calcnodes_omega9'
   allocate(weight_mu(dim_mu), stat=err)
   if(err.ne.0) stop 'allocation error in calcnodes_omega9'
!
!
!check if rule is available and calculate the order of integration scheme
   do i=1, 200
      available = available_table (i)
      dim_omega=order_table(i)
      if(available.eq.1.and.dim_omega.ge.n_omega) exit
   enddo
!
!calculate lebedev nodes and weights
   allocate(n_x(dim_omega), stat=err)
   if(err.ne.0) stop 'allocation error in calcnodes_omega9'
   allocate(n_y(dim_omega), stat=err)
   if(err.ne.0) stop 'allocation error in calcnodes_omega9'
   allocate(n_z(dim_omega), stat=err)
   if(err.ne.0) stop 'allocation error in calcnodes_omega9'
   allocate(weight_omega(dim_omega), stat=err)
   if(err.ne.0) stop 'allocation error in calcnodes_omega9'
!
!calculate nodes and weights
   call ld_by_order (dim_omega, n_x, n_y, n_z, weight_omega)
!
!!exclude zero directions (by rotation of all directions around x, y, and z-axis)
   ca=0.999999999d0
   cb=ca
   cc=ca
   sa=sqrt(1.d0-ca**2)
   sb=sa
   sc=sa
!
   do i=1, dim_omega
      n_x2=cc*cb*n_x(i)+(cc*sb*sa-sc*ca)*n_y(i)+(cc*sb*ca+sc*sa)*n_z(i)
      n_y2=sc*cb*n_x(i)+(sc*sb*sa+cc*ca)*n_y(i)+(sc*sb*ca-cc*sa)*n_z(i)
      n_z2=-sb*n_x(i)+cb*sa*n_y(i)+cb*ca*n_z(i)
      n_x(i)=n_x2
      n_y(i)=n_y2
      n_z(i)=n_z2
   enddo
!
!
!
end subroutine calcnodes_omega9
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
function n_phi_lobl(del_omega, del_theta, theta)
!
!--------calculates number of phi integration nodes for a given---------
!--(constant) del_omega at mu-level theta with corresponding del_theta--
!
   use prog_type
   use fund_const
!
   implicit none
!
! ... arguments
   integer(i4b) :: n_phi_lobl
   real(dp), intent(in) :: del_omega, del_theta, theta
!
! ... local scalars
   real(dp) :: n_phi_real
!
   n_phi_real=2.d0*pi*sin(theta)*del_theta / del_omega + 1.
   n_phi_lobl=nint(n_phi_real)
!
end function n_phi_lobl
