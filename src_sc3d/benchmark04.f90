subroutine benchmark04_solution
!
!--------------setting opacity and source-functions to zero-------------
!--------------calculating intensities for different angles-------------
!-----------integrating intensities over complete solid angle-----------
!-----------------------for 3d radiative transfer-----------------------
!------------------------(see subroutine mint_3d)-----------------------
!
use prog_type
use dimecr, only: n1d_cr, r1d_cr, norm1d_cr
use dime3d, only: mint3d, int3d, opac3d, scont3d, ndxmax, ndymax, ndzmax, z
use angles, only: dim_mu, dim_phi, weight_mu, nodes_mu, n_x, n_y, n_z
use bcondition, only: xic1, xic2
use mod_benchmark, only: mint1d_theo, mint1d_sc, mint1d_fvm, r1d_angdep, n1d_angdep, intsc_angdep, intfvm_angdep
use fund_const, only: pi
use params_input, only: rmax
use mod_interp1d, only: find_index
!
implicit none
!
! ... local scalars
integer(i4b) :: i, j, k, muindx, err
integer(i4b) :: nr_dum, iim2, iim1, ii, iip1
real(dp) :: dilfac, muc, sum
!
! ... local arrays
integer(i4b), dimension(:), allocatable :: indx1d_angdep
real(dp), dimension(:,:), allocatable :: weight_mu2d
!
! ... local functions
real(dp) :: w_otrapez_step, w_otrapez_step2
!
!-----------------------------------------------------------------------
!
allocate(mint1d_theo(n1d_cr), stat=err)
   if(err.ne.0) stop 'allocation error in benchmark04'
allocate(mint1d_sc(n1d_cr), stat=err)
   if(err.ne.0) stop 'allocation error in benchmark04'
allocate(mint1d_fvm(n1d_cr), stat=err)
   if(err.ne.0) stop 'allocation error in benchmark04'
!
!for angular dependence of intensity at 6 different radii
n1d_angdep=6
!
allocate(r1d_angdep(n1d_angdep), stat=err)
   if(err.ne.0) stop 'allocation error in benchmark04'
allocate(indx1d_angdep(n1d_angdep), stat=err)
   if(err.ne.0) stop 'allocation error in benchmark04'
allocate(intsc_angdep(n1d_angdep, dim_mu), stat=err)
   if(err.ne.0) stop 'allocation error in benchmark04'
allocate(intfvm_angdep(n1d_angdep, dim_mu), stat=err)
   if(err.ne.0) stop 'allocation error in benchmark04'

r1d_angdep=(/ 4.0d0, 4.53d0, 5.21d0, 5.79d0, 6.5d0, 7.78d0 /)
!r1d_angdep=(/ 1.1d0, 1.2d0, 1.5d0, 2.5d0, 3.d0, 7.78d0 /)
r1d_angdep=(/ 2.2d0, 3.1d0, 4.d0, 4.53d0, 5.4d0, 9.d0 /)
do i=1, n1d_angdep
!   r1d_angdep(n1d_angdep+1-i) = 1.d0+(rmax-1.d0)/(2.d0**i)
   call find_index(r1d_angdep(n1d_angdep+1-i), r1d_cr, n1d_cr, iim2, iim1, ii, iip1)
   r1d_angdep(n1d_angdep+1-i) = r1d_cr(ii)
   indx1d_angdep(n1d_angdep+1-i) = ii
enddo
!
!-----------------------------------------------------------------------
!
xic1=1.d0
xic2=0.d0
!
opac3d=1.d-10
scont3d=1.d-10
!
!-------------calculating mu-weights at each radial grid point----------
!
!nodes_mu(1)=-0.99999999d0
!nodes_mu(dim_mu)=0.9999999d0
!do i=2, dim_mu-1
!   nodes_mu(i) = nodes_mu(1) + i*(nodes_mu(dim_mu)-nodes_mu(1))/(dim_mu-1)
!enddo
!call precalc_weight_boole(nodes_mu, dim_mu, weight_mu)
!weight_mu=weight_mu/2.d0

allocate(weight_mu2d(n1d_cr,dim_mu), stat=err)
   if(err.ne.0) stop 'allocation error in subroutine benchmark04_solution'
weight_mu2d=0.d0
!
do i=1, n1d_cr
!i=10
   muc=sqrt(1.d0-1.d0/z(ndzmax-n1d_cr-1+i)**2)
!   write(*,*) 'at', z(ndzmax-n1d_cr-1+i), muc
!
   sum=0.d0
   weight_mu2d(i,1)=(nodes_mu(1)+nodes_mu(2)+2.d0)/2.d0
   if(muc.le.nodes_mu(2).and.muc.ge.nodes_mu(1)) weight_mu2d(i,1)=1.d0+muc
   sum=sum+weight_mu2d(i,1)
!
!from left
   do j=2, dim_mu-1
!      weight_mu2d(i,j)=(nodes_mu(j+1)-nodes_mu(j-1))/2.d0
      weight_mu2d(i,j)=2.d0*weight_mu(j)
      if(muc.lt.nodes_mu(j+1)) then
         weight_mu2d(i,j)=1.d0+muc-sum
         exit
      endif
      sum=sum+weight_mu2d(i,j)
   enddo
!
!from right
   sum=0.d0
   weight_mu2d(i,dim_mu)=(2.d0-nodes_mu(dim_mu)-nodes_mu(dim_mu-1))/2.d0
   if(muc.le.nodes_mu(dim_mu).and.muc.ge.nodes_mu(dim_mu-1)) weight_mu2d(i,dim_mu)=1.d0-muc
   sum=sum+weight_mu2d(i,dim_mu)
   do j=dim_mu-1, 1, -1
!      weight_mu2d(i,j)=(nodes_mu(j+1)-nodes_mu(j-1))/2.d0
      weight_mu2d(i,j)=2.d0*weight_mu(j)
      if(muc.gt.nodes_mu(j-1)) then
         weight_mu2d(i,j)=1.d0-muc-sum
         exit
      endif
      sum=sum+weight_mu2d(i,j)
   enddo
!
!   do j=1, dim_mu
!      write(*,'(i5,3es20.8)') j, nodes_mu(j), weight_mu2d(i,j)/2.d0, weight_mu(j)
!   enddo
!   sum=0.d0
!   do j=26, dim_mu
!      sum=sum+weight_mu2d(i,j)/2.d0
!   enddo
!   write(*,*) sum, (1.d0-muc)/2.d0
enddo
!
!normalization
weight_mu2d=weight_mu2d/2.d0
!
!------------calculating intensities for given nodes--------------------
!
mint1d_sc=0.d0
mint1d_fvm=0.d0
norm1d_cr=0.d0
!
write(*,*) '--------------------calculating intensities for all mu-------------------------'
write(*,*)
!
do muindx=1, dim_mu
!
   write(*,*) 'calculating muindx of dim_mu', muindx, dim_mu
!
!-------------------------short characteristics solution----------------
!
   call fsc_cont2d(muindx)

   do j=1, n1d_cr
!on positive z-axis

!***debug start
      muc=sqrt(1.d0-1.d0/z(ndzmax-n1d_cr-1+j)**2)
!      weight_mu(muindx)=w_otrapez_step2(-1.d0,1.d0,muc,muindx,dim_mu,nodes_mu)/2.d0
!overwrite intensities
!      if(muc.le.nodes_mu(muindx)) then
!         int3d(ndxmax/2+1,ndymax/2+1,ndzmax-n1d_cr-1+j)=1.d0
!      else
!         int3d(ndxmax/2+1,ndymax/2+1,ndzmax-n1d_cr-1+j)=0.d0
!      endif
!***debug end

!      mint1d_sc(j)=mint1d_sc(j) + int3d(ndxmax/2+1,ndymax/2+1,ndzmax-n1d_cr-1+j)*weight_mu(muindx)
!      norm1d_cr(j)=norm1d_cr(j) + weight_mu(muindx)
      mint1d_sc(j)=mint1d_sc(j) + int3d(ndxmax/2+1,ndymax/2+1,ndzmax-n1d_cr-1+j)*weight_mu2d(j,muindx)
      norm1d_cr(j)=norm1d_cr(j) + weight_mu2d(j,muindx)
!store angular dependence as function of intensity
      do k=1, n1d_angdep
         if(j.eq.indx1d_angdep(k)) then
            intsc_angdep(k,muindx) = int3d(ndxmax/2+1,ndymax/2+1,ndzmax-n1d_cr-1+j)
         endif
      enddo
!on negative z-axis
!      mint1d_sc(j)=mint1d_fvm(j) + int3d(ndxmax/2+1, ndymax/2+1, n1d_cr+2-j)*weight_mu(muindx)
!      norm1d_cr(j)=norm1d_cr(j) + weight_mu(muindx)
   enddo
!
!-------------------------finite volume method solution-----------------
!
   call ffvm_cont2d(muindx)
   do j=1, n1d_cr
!on positive z-axis
!
!***debug start
      muc=sqrt(1.d0-1.d0/z(ndzmax-n1d_cr-1+j)**2)
      weight_mu(muindx)=w_otrapez_step2(-1.d0,1.d0,muc,muindx,dim_mu,nodes_mu)/2.d0
!***debug end
!
      mint1d_fvm(j)=mint1d_fvm(j) + int3d(ndxmax/2+1,ndymax/2+1,ndzmax-n1d_cr-1+j)*weight_mu(muindx)
      norm1d_cr(j)=norm1d_cr(j) + weight_mu(muindx)
!store angular dependence as function of intensity
      do k=1, n1d_angdep
         if(j.eq.indx1d_angdep(k)) then
            intfvm_angdep(k,muindx) = int3d(ndxmax/2+1,ndymax/2+1,ndzmax-n1d_cr-1+j)
         endif
      enddo
!on negative z-axis
!      mint1d_fvm(j)=mint1d_fvm(j) + int3d(ndxmax/2+1, ndymax/2+1, n1d_cr+2-j)*weight_mu(muindx)
!      norm1d_cr(j)=norm1d_cr(j) + weight_mu(muindx)
   enddo
!
enddo
!
!-----------------------------------------------------------------------
!
!theoretical solution
write(*,'(4a20)') 'radius', 'J(theo)', 'J(2d-sc)', 'J(2d-fvm)'
do i=1, n1d_cr
   dilfac=1.d0-sqrt(1.d0-1.d0/r1d_cr(i)**2)
   mint1d_theo(i) = 0.5d0*xic1*dilfac
   write(*,'(4es20.8)') r1d_cr(i), mint1d_theo(i), mint1d_sc(i), mint1d_fvm(i)
enddo
write(*,*)
!
end subroutine benchmark04_solution
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
function w_otrapez_step(a, b, xc, ii, n, x)
!
!   calculates weight for an integral between [a,b] with open boundaries at a, b
!                for a step-fct approach at a position ix
!
!   I = w1*f1 + w2*f2 + ... wi*fi + .... wn*fn
!
! on input:   integration bounds:   a,b
!             xc:                   x-value of step
!             x:                    x-array
!             n:                   data points of x,f array
!             ii:                   current index for which weight shall be calculated
!
use prog_type
use mod_interp1d, only: find_index  
!
implicit none
!
! ... arguments
real(dp) :: w_otrapez_step
real(dp), intent(in) :: a, b, xc
integer(i4b), intent(in) :: ii, n
real(dp), dimension(n), intent(in) :: x
!
! ... local scalars
integer(i4b) :: jjm2, jjm1, jj, jjp1
real(dp) :: dxi, dxip1, dxc, dxc1, dxc2
!
call find_index(xc, x, n, jjm2, jjm1, jj, jjp1)
!
if(ii.eq.1) then
   if(ii.eq.jj-2) then
      dxip1=x(ii+1)-x(ii)
      dxc=xc-x(ii+1)
      w_otrapez_step = x(ii)-a+dxip1/2.d0-dxc**2/2.d0/dxip1
   elseif(ii.eq.jj-1) then
      dxip1=x(ii+1)-x(ii)
      dxc=xc-x(ii)
      w_otrapez_step = xc-a-dxc**2/2.d0/dxip1
   else
      dxip1=x(2)-x(1)
      w_otrapez_step = x(1)-a + dxip1/2.d0
   endif
elseif(ii.eq.2) then
   if(ii.eq.jj-2) then
      dxi=x(ii)-x(ii-1)
      dxip1=x(ii+1)-x(ii)
      dxc=xc-x(ii+1)
      w_otrapez_step = (dxi+dxip1)/2.d0 - dxc**2/2.d0/dxip1
   elseif(ii.eq.jj-1) then
      dxi=x(ii)-x(ii-1)
      dxc=xc-x(ii)
      w_otrapez_step = dxi/2.d0 +dxc + dxc**2/2.d0/dxi
   elseif(ii.eq.jj) then
      dxi=x(ii)-x(ii-1)
      dxip1=x(ii+1)-x(ii)
      dxc1=x(ii)-xc
      dxc2=xc-x(ii-1)
      w_otrapez_step = dxc1 + dxc2**2/2.d0/dxi + dxc1**2/2.d0/dxip1 + dxip1/2.d0
   else
      dxi=x(ii)-x(ii-1)
      dxip1=x(ii+1)-x(ii)
      w_otrapez_step = (dxi+dxip1)/2.d0
   endif
elseif(ii.eq.n) then
   if(ii.eq.jj) then
      dxi=x(ii)-x(ii-1)
      dxc=xc-x(ii)
      w_otrapez_step = b-xc - dxc**2/2.d0/dxi
   elseif(ii.eq.jj+1) then
      dxi=x(ii)-x(ii-1)
      dxc=xc-x(ii-1)
      w_otrapez_step = b-x(ii) + dxi/2.d0 - dxc**2/2.d0/dxi
   else
      dxi=x(ii)-x(ii-1)
      w_otrapez_step = b-x(ii)+dxi/2.d0
   endif
elseif(ii.eq.n-1) then
   if(ii.eq.jj-1) then
      dxi=x(ii)-x(ii-1)
      dxip1=x(ii+1)-x(ii)
      dxc1=xc-x(ii)
      dxc2=xc-x(ii+1)
      w_otrapez_step = dxc1 + dxi/2.d0 + dxc1**2/2.d0/dxi + dxc2**2/2.d0/dxip1
   elseif(ii.eq.jj) then
      dxip1=x(ii+1)-x(ii)
      dxc=x(ii)-xc
      w_otrapez_step = dxc + dxip1/2.d0 + dxc**2/2.d0/dxip1
   elseif(ii.eq.jj+1) then
      dxi=x(ii)-x(ii-1)
      dxip1=x(ii+1)-x(ii)
      dxc=xc-x(ii-1)
      w_otrapez_step = (dxi+dxip1)/2.d0 - dxc**2/2.d0/dxi
   else
      dxi=x(ii)-x(ii-1)
      dxip1=x(ii+1)-x(ii)
      w_otrapez_step = (dxi+dxip1)/2.d0
   endif
elseif(ii.eq.jj-2) then
   dxi=x(ii)-x(ii-1)
   dxip1=x(ii+1)-x(ii)
   dxc=xc-x(ii+1)
   w_otrapez_step = (dxi+dxip1)/2.d0-dxc**2/2.d0/dxip1
elseif(ii.eq.jj-1) then
   dxi=x(ii)-x(ii-1)
   dxc=xc-x(ii)
   w_otrapez_step = dxi/2.d0 + dxc + dxc**2/2.d0/dxi
elseif(ii.eq.jj) then
   dxip1=x(ii+1)-x(ii)
   dxc=x(ii)-xc
   w_otrapez_step = dxip1/2.d0 + dxc + dxc**2/2.d0/dxip1
elseif(ii.eq.jj+1) then
   dxi=x(ii)-x(ii-1)
   dxip1=x(ii+1)-x(ii)
   dxc=xc-x(ii-1)
   w_otrapez_step = (dxi+dxip1)/2.d0-dxc**2/2.d0/dxi
else
   dxi=x(ii)-x(ii-1)
   dxip1=x(ii+1)-x(ii)
   w_otrapez_step = (dxi+dxip1)/2.d0
endif  
!
!
!
end function w_otrapez_step
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
function w_otrapez_step2(a, b, xc, ii, n, x)
!
!   calculates weight for an integral between [a,b] with open boundaries at a, b
!                for a step-fct approach at a position ix
!                  assuming constant 'extrapolation'
!
!   I = w1*f1 + w2*f2 + ... wi*fi + .... wn*fn
!
! on input:   integration bounds:   a,b
!             xc:                   x-value of step
!             x:                    x-array
!             n:                    data points of x,f array
!             ii:                   current index for which weight shall be calculated
!
use prog_type
use mod_interp1d, only: find_index  
!
implicit none
!
! ... arguments
real(dp) :: w_otrapez_step2
real(dp), intent(in) :: a, b, xc
integer(i4b), intent(in) :: ii, n
real(dp), dimension(n), intent(in) :: x
!
! ... local scalars
integer(i4b) :: jjm2, jjm1, jj, jjp1
real(dp) :: dxi, dxip1, dxc, dxc1, dxc2
!
call find_index(xc, x, n, jjm2, jjm1, jj, jjp1)
!
if(ii.eq.1) then
   if(ii.eq.jj-1) then
      w_otrapez_step2 = xc-a
   else
      w_otrapez_step2 = (x(1)+x(2)-2.d0*a)/2.d0
   endif
elseif(ii.eq.n) then
   if(ii.eq.jjp1) then
      w_otrapez_step2 = (2.d0*b-x(ii)-x(ii-1))/2.d0
   elseif(ii.eq.jj) then
      w_otrapez_step2 = b-xc
   else
      w_otrapez_step2 = (2.d0*b-x(ii)-x(ii-1))/2.d0
   endif
else
   if(ii.eq.jj-1) then
      w_otrapez_step2 = (2.d0*xc-x(ii)-x(ii-1))/2.d0
   elseif(ii.eq.jj) then
      w_otrapez_step2 = (x(ii+1)+x(ii)-2.d0*xc)/2.d0
   else
      w_otrapez_step2 = (x(ii+1)-x(ii-1))/2.d0
   endif
endif
!
!
!
end function w_otrapez_step2
