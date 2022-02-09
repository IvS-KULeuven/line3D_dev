subroutine gridxyz
!
!-----------------------------------------------------------------------
!--------------sets up x, y, z-grid (from radial grid)------------------
!------------------(with equidistant core points)-----------------------
!-----------------------------------------------------------------------
!
use prog_type
use fund_const
use dime3d, only: ndx, ndy, ndz, ncx, ncy, ncz, ndxmax, ndymax, ndzmax, x, y, z
use dime1d, only: n1d, r1d, r1d_dum
use inf_reg, only: rmin
!
implicit none
!
! ... local scalars
integer(i4b) :: i, err
real(dp) :: del, delcx, delcy, delcz
!
!set up radial grid
call gridr
!
!----------------------set up x, y, z grid----------------------------
!
allocate(x(ndxmax), stat=err)
   if(err.ne.0) stop 'error in gridxyz: allocation'
allocate(y(ndymax), stat=err)
   if(err.ne.0) stop 'error in gridxyz: allocation'
allocate(z(ndzmax), stat=err)
   if(err.ne.0) stop 'error in gridxyz: allocation'
!
!equidistant grid steps inside core
delcx=rmin/float(ncx)
delcy=rmin/float(ncy)
delcz=rmin/float(ncz)
!
!fill gridpoints from -rmax to photosphere
!both outermost grid points shall have same coordinates
!   (will be overwritten with phantom point later on)
z(1)=-r1d(n1d)
do i=2, ndz-ncz
   z(i)=-r1d(n1d-i+2)
enddo
!
!fill ncx grifpoints from photosphere to zero
do i=ndz-ncz+1, ndz
   z(i) = z(i-1) + delcz
enddo
!
!set origint to exactly zero
z(ndz)=zero
!
!mirror grid on forward hemisphere
do i=1, ndz-1
   z(ndzmax-i+1)=abs(z(i))
enddo
!
!check if dimensions of x and z grid are compatible
if(ndzmax.ne.ndxmax.or.ncx.ne.ncz) stop 'error in gridxyz: dimensions of x and z not compatible'
!set x equal to z-grid
x=z
!
if(ndymax.ne.ndzmax.or.ncy.ne.ncz) stop 'error in gridxyz: dimensions of y and z not compatible'
y=z
!
!deallocate 1d arrays, since not needed anymore
deallocate(r1d)
deallocate(r1d_dum)
!
!
end subroutine gridxyz
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine gridxyz_equi
!
!-----------------------------------------------------------------------
!----------------------sets up x, y, z -grid----------------------------
!----------with all points distributed equidistantly--------------------
!-----------------------------------------------------------------------
!
use prog_type
use fund_const
use dime3d, only: ndxmax, ndymax, ndzmax, x, y, z
use params_input, only: rmax
!
implicit none
!
! ... local scalars
integer(i4b) :: i, err
reAl(dp) :: del, delcx, delcy, delcz
!
!----------------------set up x, y, And z grid--------------------------
!
if(mod(ndxmax, 2).eq.0) stop 'error in gridxyz_equi: ndxmax needs to be odd!'
if(mod(ndymax, 2).eq.0) stop 'error in gridxyz_equi: ndymax needs to be odd!'
if(mod(ndzmax, 2).eq.0) stop 'error in gridxyz_equi: ndzmax needs to be odd!'
!
allocate(x(ndxmax), stat=err)
   if(err.ne.0) stop 'error in gridxyz_equi: allocation'
allocate(y(ndymax), stat=err)
   if(err.ne.0) stop 'error in gridxyz_equi: allocation'
allocate(z(ndzmax), stat=err)
   if(err.ne.0) stop 'error in gridxyz_equi: allocation'
!
do i=2, ndxmax-1
   x(i) = -rmax + (i-2) * two*rmax/(ndxmax-3)
enddo
x(1) = -rmax
x(ndxmax)=rmax

do i=2, ndymax-1
   y(i) = -rmax + (i-2) * two*rmax/(ndymax-3)
enddo
y(1)=-rmax
y(ndymax)=rmax

do i=2, ndzmax-1
   z(i) = -rmax + (i-2) * two*rmax/(ndzmax-3)
enddo
z(1)=-rmax
z(ndzmax)=rmax
!
!
!include grid point at radius 1
do i=1, ndxmax
   if(x(i).eq.one) exit
   if(x(i).gt.one) then
      x(i) = one
      exit
   endif
enddo
do i=1, ndxmax
   if(x(i).eq.-one) exit
   if(x(i).gt.-one) then
      x(i-1) = -one
      exit
   endif
enddo
do i=1, ndymax
   if(y(i).eq.one) exit
   if(y(i).gt.one) then
      y(i) = one
      exit
   endif
enddo
do i=1, ndymax
   if(y(i).eq.-one) exit
   if(y(i).gt.-one) then
      y(i-1) = -one
      exit
   endif
enddo
do i=1, ndzmax
   if(z(i).eq.one) exit
   if(z(i).gt.one) then
      z(i) = one
      exit
   endif
enddo
do i=1, ndzmax
   if(z(i).eq.-one) exit
   if(z(i).gt.-one) then
      z(i-1) = -one
      exit
   endif
enddo
!
!
!
end subroutine gridxyz_equi
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine gridxyz_opt
!
!
use prog_type
use fund_const
use dime1d, only: n1d_dum, r1d_dum, r1d
use dime3d, only: ndx, ndy, ndz, ndxmax, ndymax, ndzmax, x, y, z
use mod_sort, only: bubblesort
use mod_grid, only: recalc_grid1d
!use params_stellar, only: sr
!
implicit none
!
! ... local scalars
integer(i4b) :: i, j, k, lx, ly, lz
integer(i4b) :: err, nuniq_x, nuniq_y, nuniq_z, nav, nrest_uniq, nrest_grid, indx, nphi, ntheta, ndx_dum
integer(i4b), parameter :: nphi_max=45
real(dp), parameter :: deltau_max=one/3.d0
real(dp) :: val1, val2, xav
real(dp) :: delx, delx1, delx2, max_delx, maxx, dum
!
! ... local errors
real(dp), dimension(:), allocatable :: xcoords, ycoords, zcoords, xcoords_dum, ycoords_dum, zcoords_dum
real(dp), dimension(:), allocatable :: phi, theta
real(dp), dimension(:), allocatable :: xdum
!
!----------------------set up radial grid-------------------------------
!
call gridr
!
!---------------create dummy phi and theta-grid-------------------------
!phi is in interval [0,pi/2]
!theta is in interval [0,pi/2]
!   => first octant
!
!calculate optimum number of nphi
!dum=acos(deltau_max/opath1d(1)/sr) * two/pi
!nphi=int((dum-two)/(dum-one))
!nphi=min(nphi,nphi_max)
!
!
!nphi=9 (nc=65 with nd=120)
nphi=45
ntheta=5
!
allocate(phi(nphi), stat=err)
   if(err.ne.0) stop 'error in gridxyz_opt: allocation'
allocate(theta(ntheta), stat=err)
   if(err.ne.0) stop 'error in gridxyz_opt: allocation'
!
!
do i=1, nphi
   phi(i) = float(i-1) * pi / float(nphi-1) / two
enddo
!
do i=1, ntheta
   theta(i) = float(i-1) * pi / float(ntheta-1) / two
enddo
!
!-----------------------------------------------------------------------
!
!*******************old (2d) version************************************
!-----------------------------------------------------------------------
!
!number of unique coordinates of all radial vectors
nuniq_x=n1d_dum*(nphi-1) + 1
!
allocate(xcoords(nuniq_x), stat=err)
   if(err.ne.0) stop 'allocation error gridxyz_opt: xcoords'
xcoords=ten
!
!calculate all x-coordinates (in stellar radii)
k=1
do i=1, nphi-1
   dum=cos(phi(i))
   do j=1, n1d_dum
      xcoords(k) = r1d_dum(j) * dum
      k=k+1
   enddo
enddo
!
xcoords(k)=zero
!
!---------------------create grid from zero zo rmax---------------------
!
!xdum only up to ndx-1 since phantom point will be calculated later on
ndx_dum=ndx-1
allocate(xdum(ndx_dum), stat=err)
   if(err.ne.0) stop 'allocation error gridxyz_opt: xdum'
!
!calculate new grid from xcoords
call recalc_grid1d(xcoords, nuniq_x, ndx_dum, 3, xdum)
!
!want special fixed points (zero, one, rmax)
xdum(ndx_dum-2)=one
xdum(ndx_dum-1)=zero
xdum(ndx_dum)=xcoords(nuniq_x)
!
call bubblesort(xdum, ndx_dum)

!
!debug
!write(*,*) nphi, n1d_dum
!write(*,*) r1d_dum
!write(*,*)
!write(*,*) xdum
!open(1,file='TRASH/grid1.dat')
!   do i=1, ndx_dum
!      write(1,*) i, xdum(i)
!   enddo
!close(1)
!stop 'go on in grid_spat'
!
!--------------------check if grid-values occur twice-------------------
!
do i=2, ndx_dum-1
   val1=xdum(i-1)
   val2=xdum(i)
   if(abs(val2-val1).lt.1.d-14) then
!insert value on half range (to upper grid point)
      xdum(i)=(xdum(i+1)+xdum(i))/two
   endif
enddo
!
!same on outer boundary
val1=xdum(ndx_dum)
val2=xdum(ndx_dum-1)
if(abs(val2-val1).lt.1.d-14) then
!insert value on half range (to lower grid point)
   xdum(ndx_dum-1)=(xdum(ndx_dum-2)+xdum(ndx_dum-1))/two
endif
!
!-----make outer part of grid equidistant in order that delta(xdum)-----
!---------------------is not decreasing---------------------------------
!note: maybe logarithmic increment is much better
!
do i=2, ndx_dum
!apply algorithm only if xdum gt 1
   if(xdum(i).gt.one) then
      delx1=xdum(i)-xdum(i-1)
      max_delx=delx1*(ndx_dum-i)
      maxx=xdum(i)+max_delx
      if(maxx.lt.xcoords(nuniq_x)) then
         indx=i
      endif
   endif
enddo
!
if(ndx_dum.ne.indx) then
   delx=(xcoords(nuniq_x)-xdum(indx))/(ndx_dum-indx)
   do i=indx+1, ndx_dum
      xdum(i)=xdum(i-1)+delx
   enddo
endif
!
!-------------------create complete grid--------------------------------
!
allocate(x(ndxmax), stat=err)
   if(err.ne.0) stop 'error in gridxyz_opt: allocation'
allocate(y(ndymax), stat=err)
   if(err.ne.0) stop 'error in gridxyz_opt: allocation'
allocate(z(ndzmax), stat=err)
   if(err.ne.0) stop 'error in gridxyz_opt: allocation'
!
x(1)=-xdum(ndx_dum)
do i=2, ndx_dum
   x(i)=-xdum(ndx_dum+2-i)
enddo
!
do i=1, ndx_dum
   x(ndx_dum + i) = xdum(i)
enddo
!
x(ndxmax)=xdum(ndx_dum)
!
!check if dimensions of x and z grid are compatible
if(ndzmax.ne.ndxmax.or.ndz.ne.ndx) stop 'dimensions of z and x not compatible'
!set x equal to z-grid
z=x
!
!check if dimensions of y and z grid are compatible
if(ndymax.ne.ndxmax.or.ndy.ne.ndx) stop 'dimensions of y and x not compatible'
!set y equal to z-grid
y=x
!
!deallocate 1d arrays, since not needed anymore
deallocate(r1d)
deallocate(r1d_dum)

!write(*,*) x
!stop 'go on in gridxyz_opt'
!
!
end subroutine gridxyz_opt
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine gridxyz_optb
!
!calculates a radial grid from input parameters
!calculates a (not necessarily) uniform angle-grid
!
!calculates the probability density of the radial grid,
!and corresponding probability density of the x,y,z-coordinates
!
!calculates the actualy x,y,z-coordinates following this distrubution
!   (with input nc_final, nnc_final: number of points inside core and outside core)
!
!
use prog_type
use fund_const
use mod_directories, only: output_dir
use dime1d, only: n1d, n1d_dum, r1d_dum, r1d
use dime3d, only: ncx, ndx, ncz, ndz, ndy, ndxmax, ndymax, ndzmax, x, y, z
use params_stellar, only: sr
use mod_interp1d, only: find_index
use mod_grid, only: calc_pdf, grid_log
!
implicit none
!
! ... local scalars
integer(i4b) :: i, j, k, err
integer(i4b) :: iim2, iim1, ii, iip1
integer(i4b) :: nphi, nx, nx_c, nx_nc, nz, nz_c, nz_nc, &
                nxfine_c, nxfine_nc, nxfine, ndx_dum, nncx_final, nx_final, ncx_final, &
                nzfine_c, nzfine_nc, nzfine, ndz_dum, nncz_final, nz_final, ncz_final
real(dp) :: integral, cc, grad, xinterc, zinterc, del, rmin, rmax, xmin, xmax, zmin, zmax, &
            x_jm1, x_j, z_jm1, z_j, r_jm1, r_j, phi_jm1, phi_j, pdfr_jm1, pdfr_j, &
            pdfp_jm1, pdfp_j, h_jm1, h_j, prob_nc, prob_c, pc, pc_i, xi, zi, delp, &
            delx1, maxx, max_delx
!
! ... local arrays
!real(dp), dimension(:), allocatable :: xcoords, ycoords, zcoords, xcoords_dum, ycoords_dum, zcoords_dum
real(dp), dimension(:), allocatable :: phi, phidum, phi_mid, pdf_phi, p_phi, &
                                       r1d_mid, pdf_r1d, p_r1d, &
                                       xfine_c, xfine_mid_c, pdf_xfine_c, p_xfine_c, &
                                       zfine_c, zfine_mid_c, pdf_zfine_c, p_zfine_c, &
                                       xfine_nc, xfine_mid_nc, pdf_xfine_nc, p_xfine_nc, &
                                       zfine_nc, zfine_mid_nc, pdf_zfine_nc, p_zfine_nc, &
                                       xfine, xfine_mid, pdf_xfine, p_xfine, &
                                       zfine, zfine_mid, pdf_zfine, p_zfine
real(dp), dimension(:), allocatable :: xdum, xfinal, xfinal_c, xfinal_nc, xfinal_mid, pdf_xfinal, p_xfinal, &
                                       zdum, zfinal, zfinal_c, zfinal_nc, zfinal_mid, pdf_zfinal, p_zfinal
!
!----------------------set up radial grid-------------------------------
!
call gridr
!
rmin=r1d_dum(1)
rmax=r1d_dum(n1d_dum)
!
allocate(r1d_mid(n1d_dum-1), stat=err)
   if(err.ne.0) stop 'error in gridxyz_optb: allocation'
allocate(pdf_r1d(n1d_dum-1), stat=err)
   if(err.ne.0) stop 'error in gridxyz_optb: allocation'
allocate(p_r1d(n1d_dum-1), stat=err)
   if(err.ne.0) stop 'error in gridxyz_optb: allocation'
!
!calculate probability density
call calc_pdf(n1d_dum, r1d_dum, r1d_mid, pdf_r1d, p_r1d)
!
!---------------create dummy phi grid (measured from x to z-axis)-------
!phi is in interval [0,pi/2]
!
nphi=45
!
allocate(phi(nphi), stat=err)
   if(err.ne.0) stop 'error in gridxyz_optb: allocation'
allocate(phidum(nphi), stat=err)
   if(err.ne.0) stop 'error in gridxyz_optb: allocation'
allocate(phi_mid(nphi-1), stat=err)
   if(err.ne.0) stop 'error in gridxyz_optb: allocation'
allocate(pdf_phi(nphi-1), stat=err)
   if(err.ne.0) stop 'error in gridxyz_optb: allocation'
allocate(p_phi(nphi-1), stat=err)
   if(err.ne.0) stop 'error in gridxyz_optb: allocation'
!
!phi grid equidistant
do i=1, nphi
   phi(i) = float(i-1) * pi / float(nphi-1) / two
enddo
!
!phi grid with high resolution at small angles
!call grid_log(1.d-3,pi/two,nphi,phi)
!!
!!phi grid with high resolution at large angles
!phidum(nphi)=pi/two
!k=2
!do i=nphi-1, 1, -1
!   phidum(i) = phidum(i+1)-(phi(k)-phi(k-1))
!   k=k+1
!enddo
!phi=phidum
!
!calculate probability density
call calc_pdf(nphi, phi, phi_mid, pdf_phi, p_phi)
!
!-------------calculate pdf for x-coordinates on a fine grid------------
!---------------1. pdf for core points (normalized to core)-------------
!---------------2. pdf for non-core points (normlized to non-core)------
!---------------3. complete pdf (normalized to complete range)----------
!
!----------------------------core pdf-----------------------------------
!
nxfine_c=100
allocate(xfine_c(nxfine_c), stat=err)
allocate(xfine_mid_c(nxfine_c-1), stat=err)
allocate(pdf_xfine_c(nxfine_c-1), stat=err)
allocate(p_xfine_c(nxfine_c-1), stat=err)
!
!core grid confined towards rmin
cc=half
grad = rmin/log10((cc+one)/cc)
zinterc = -rmin*log10(cc)/log10((cc+one)/cc)
do i=1, nxfine_c
   xfine_c(i) = grad * log10(float(i)/(nxfine_c)+cc) + zinterc
enddo
xfine_c(1)=zero
xfine_c(nxfine_c)=rmin
!
!
!create mid points
do i=1, nxfine_c-1
   xfine_mid_c(i) = (xfine_c(i+1)+xfine_c(i))/two
enddo
!
!
!calculate pdf for core coordinates
nz_c=200
nz_nc=400
nz=nz_c+nz_nc
allocate(zdum(nz), stat=err)
!
do i=1, nxfine_c-1
!
!create z-coordinates for z-integration (obtain marignal distribution by integrating over z)
   zmin = sqrt(rmin**2-xfine_mid_c(i)**2)
   zmax = sqrt(rmax**2-xfine_mid_c(i)**2)
!
   cc=half
   grad = (rmin-zmin)/log10((cc+one)/cc)
   zinterc = (zmin*log10(cc+one)-rmin*log10(cc))/log10((cc+one)/cc)
   k=1
   do j=1, nz_c
      zdum(k) = grad * log10(float(j-1)/(nz_c-1)+cc) + zinterc
      k=k+1
   enddo
   del=log10(zmax/rmin)/float(nz_nc)
   do j=1, nz_nc
      zdum(k) = zdum(k-1)*ten**del
      k=k+1
   enddo
!
!calculate pdf for the given x
   integral=zero
   z_jm1 = zdum(1)
!
!calculate pdf for r and phi distribtion of the given (x,z)
   r_jm1 = sqrt(xfine_mid_c(i)**2+zdum(1)**2)
   phi_jm1 = atan(zdum(1)/xfine_mid_c(i))
   call find_index(r_jm1, r1d_dum, n1d_dum, iim2, iim1, ii, iip1)
   pdfr_jm1 = pdf_r1d(iim1)
   call find_index(phi_jm1, phi, nphi, iim2, iim1, ii, iip1)
   pdfp_jm1 = pdf_phi(iim1)
!   pdfp_jm1 = two/pi
!calculate joint pdf
   h_jm1 = pdfr_jm1*pdfp_jm1/sqrt(xfine_mid_c(i)**2+z_jm1**2)
!
   do j=2, nz
      z_j = zdum(j)
!calculate pdf for r and phi distribtion of the given (x,z)
      r_j = sqrt(xfine_mid_c(i)**2+zdum(j)**2)
      phi_j = atan(zdum(j)/xfine_mid_c(i))
      call find_index(r_j, r1d_dum, n1d_dum, iim2, iim1, ii, iip1)
      pdfr_j = pdf_r1d(iim1)
      call find_index(phi_j, phi, nphi, iim2, iim1, ii, iip1)
      pdfp_j = pdf_phi(iim1)
!      pdfp_j = two/!pi
!calculate joint pdf
      h_j = pdfr_j*pdfp_j/sqrt(xfine_mid_c(i)**2+z_j**2)
!
!perform integration over z
      integral = integral + half*(h_jm1+h_j)*(z_j-z_jm1)
!
      z_jm1 = z_j
      h_jm1 = h_j
!
   enddo
!
   pdf_xfine_c(i) = integral
!
enddo
!
!
!
!calculate probability of finding x in corresponding bins
do i=1, nxfine_c-1
   p_xfine_c(i) = pdf_xfine_c(i)*(xfine_c(i+1)-xfine_c(i))
!   write(*,*) xfine_mid_c(i), pdf_xfine_c(i)*four, p_xfine_c(i)
enddo
!
prob_c = sum(p_xfine_c)
!
!normalize pdf to core interval
pdf_xfine_c=pdf_xfine_c/prob_c
!
!--------------------------non core pdf---------------------------------
!
nxfine_nc=200
allocate(xfine_nc(nxfine_nc), stat=err)
allocate(xfine_mid_nc(nxfine_nc-1), stat=err)
allocate(pdf_xfine_nc(nxfine_nc-1), stat=err)
allocate(p_xfine_nc(nxfine_nc-1), stat=err)
!
!non-core grid logarithmic
call grid_log(rmin, rmax, nxfine_nc, xfine_nc)
!
!create mid points
do i=1, nxfine_nc-1
   xfine_mid_nc(i) = (xfine_nc(i+1)+xfine_nc(i))/two
enddo
!
!
!
deallocate(zdum)
nz_c=100
nz_nc=100
nz=nz_c+nz_nc
allocate(zdum(nz), stat=err)
!
!calculate pdf for non-core x-coordinates (on half grid)
do i=1, nxfine_nc-1
!
!create z-coordinates for z-integration (obtain marignal distribution by integrating over z)
   zmin = zero
   zmax = sqrt(rmax**2-xfine_mid_nc(i)**2)
!
   cc=half
   grad = (rmin-zmin)/log10((cc+one)/cc)
   zinterc = (zmin*log10(cc+one)-rmin*log10(cc))/log10((cc+one)/cc)
   zdum(1) = zmin
   k=2
   do j=2, nz_c
      zdum(k) = grad * log10(float(j-1)/(nz_c-1)+cc) + zinterc
      k=k+1
   enddo
   del=log10(zmax/rmin)/float(nz_nc)
   do j=1, nz_nc
      zdum(k) = zdum(k-1)*ten**del
      k=k+1
   enddo
!
!
!calculate pdf for the given x
   integral=zero
   z_jm1 = zdum(1)
!
!calculate pdf for r and phi distribtion of the given (x,z)
   r_jm1 = sqrt(xfine_mid_nc(i)**2+zdum(1)**2)
   phi_jm1 = atan(zdum(1)/xfine_mid_nc(i))
   call find_index(r_jm1, r1d_dum, n1d_dum, iim2, iim1, ii, iip1)
   pdfr_jm1 = pdf_r1d(iim1)
   call find_index(phi_jm1, phi, nphi, iim2, iim1, ii, iip1)
   pdfp_jm1 = pdf_phi(iim1)
!   pdfp_jm1 = two/pi
!calculate joint pdf
   h_jm1 = pdfr_jm1*pdfp_jm1/sqrt(xfine_mid_nc(i)**2+z_jm1**2)
!
   do j=2, nz
      z_j = zdum(j)
!calculate pdf for r and phi distribtion of the given (x,z)
      r_j = sqrt(xfine_mid_nc(i)**2+zdum(j)**2)
      phi_j = atan(zdum(j)/xfine_mid_nc(i))
      call find_index(r_j, r1d_dum, n1d_dum, iim2, iim1, ii, iip1)
      pdfr_j = pdf_r1d(iim1)
      call find_index(phi_j, phi, nphi, iim2, iim1, ii, iip1)
      pdfp_j = pdf_phi(iim1)
!      pdfp_j = two/!pi
!calculate joint pdf
      h_j = pdfr_j*pdfp_j/sqrt(xfine_mid_nc(i)**2+z_j**2)
!
!perform integration over z
      integral = integral + half*(h_jm1+h_j)*(z_j-z_jm1)
!
      z_jm1 = z_j
      h_jm1 = h_j
!
   enddo
!
   pdf_xfine_nc(i) = integral
!
enddo
!
!
!
!calculate probability of finding x in corresponding bins
do i=1, nxfine_nc-1
   p_xfine_nc(i) = pdf_xfine_nc(i)*(xfine_nc(i+1)-xfine_nc(i))
enddo
!
!
prob_nc=sum(p_xfine_nc)
!
!normalize pdf to non-core interval
pdf_xfine_nc=pdf_xfine_nc/prob_nc
!
!------------complete pdf normalized to complete range [0.,rmax]-------------
!-----------------------(only required for checks)---------------------------
!
nxfine = nxfine_c + nxfine_nc -1 
allocate(xfine(nxfine), stat=err)
allocate(xfine_mid(nxfine-1), stat=err)
allocate(pdf_xfine(nxfine-1), stat=err)
allocate(p_xfine(nxfine-1), stat=err)
!
!
k=1
do i=1, nxfine_c
   xfine(k) = xfine_c(i)
   k=k+1
enddo
do i=2, nxfine_nc
   xfine(k) = xfine_nc(i)
   k=k+1
enddo
!
k=1
do i=1, nxfine_c-1
   xfine_mid(k) = xfine_mid_c(i)
   pdf_xfine(k) = pdf_xfine_c(i)*prob_c
   p_xfine(k) = p_xfine_c(i)
   k=k+1
enddo
do i=1, nxfine_nc-1
   xfine_mid(k) = xfine_mid_nc(i)
   pdf_xfine(k) = pdf_xfine_nc(i)*prob_nc
   p_xfine(k) = p_xfine_nc(i)
   k=k+1
enddo
!
pdf_xfine=pdf_xfine/sum(p_xfine)
!
!---------------------create final grid from zero to rmax---------------
!
!ndx_dum=ndx-1
!nx_final=ndx_dum
!nncx_final=ndx_dum-ncx_final+1   !+1 because rmin point occurrs twice
!write(*,*) nx_final, ncx_final, nncx_final
!
ndx_dum=ndx-1
nx_final=ndx_dum
ncx_final = ncx
nncx_final=ndx_dum-ncx_final+1   !+1 because rmin point occurrs twice
!
write(*,*) 'probability of finding x-coordinate inside core ', prob_c
write(*,*) 'probability of finding x-coordinate outside core', prob_nc
write(*,*) ' => optimum number of core points to match distribution    ', nint(ndx_dum*prob_c)
write(*,*) ' => optimum number of non-core points to match distribution', ndx_dum-nint(ndx_dum*prob_c)+1
write(*,*)
write(*,*) 'actually used number of core points    ', ncx_final
write(*,*) 'actually used number of non-core points', nncx_final
write(*,*) '# core points / # points    ', float(ncx_final)/float(nncx_final+ncx_final)
write(*,*) '# non-core points / # points', float(nncx_final)/float(nncx_final+ncx_final)
write(*,*)
!
!xdum only up to ndx-1 since phantom point will be calculated later on
allocate(xfinal_c(ncx_final), stat=err)
allocate(xfinal_nc(nncx_final), stat=err)
allocate(xfinal(nx_final), stat=err)
!
!-------------------distribute nxc_final grid points in core------------
!-------------------------(following the core pdf)----------------------
!
!probability that has to be matched for each interval
pc = one/float(ncx_final-1)
xfinal_c(1)=zero
!
pc_i = zero
do i=2, ncx_final
!probability that has to be matched
   pc_i = pc_i + pc
!
   delp = zero
   do j=1, nxfine_c-1
      xi =  (pc_i-delp)/pdf_xfine_c(j)  + xfine_c(j)
      xfinal_c(i) = xi
      if(xi.le.xfine_c(j+1)) then
!correct position found: go to next x_final-coordinate
         exit
      else
!go to next interval
         delp = delp + pdf_xfine_c(j)*(xfine_c(j+1)-xfine_c(j))
      endif
   enddo
enddo
!
!---------------distribute nncx_final grid points outside core----------
!-----------------------(following the non-core pdf)--------------------
!
!probability that has to be matched for each interval
pc = one/float(nncx_final-1)
xfinal_nc(1)=xfine_nc(1)
!
pc_i = zero
do i=2, nncx_final
!probability that has to be matched
   pc_i = pc_i + pc
!
   delp = zero
   do j=1, nxfine_nc-1
      xi =  (pc_i-delp)/pdf_xfine_nc(j)  + xfine_nc(j)
      xfinal_nc(i) = xi
      if(xi.le. xfine_nc(j+1)) then
!;correct position found: go to next x_final-coordinate
         exit
      else
!go to next interval
         delp = delp + pdf_xfine_nc(j)*(xfine_nc(j+1)-xfine_nc(j))
      endif
   enddo
enddo
!
!----------------------combined distribution----------------------------
!
k=1
do i=1, ncx_final
   xfinal(k) = xfinal_c(i)
   k=k+1
enddo
xfinal_c(ncx_final)=rmin
do i=1, nncx_final-1
   xfinal(k) = xfinal_nc(i+1)
   k=k+1
enddo
xfinal(nx_final)=rmax
!
!calculate final pdf (only for checks)
allocate(xfinal_mid(nx_final-1))
allocate(pdf_xfinal(nx_final-1))
allocate(p_xfinal(nx_final-1))
!
call calc_pdf(nx_final, xfinal, xfinal_mid, pdf_xfinal, p_xfinal)
!
!--------------------------print out everything-------------------------
!
open(1, file=output_dir//'/pdf_rgridf.dat')
   write(1,*) 'radial grid (input)'
   write(1,'(3a20)') 'r_mid', 'pdf(r)', 'p(r)'
   do i=1, n1d_dum-1 
      write(1,'(3es20.8)') r1d_mid(i), pdf_r1d(i), p_r1d(i)
   enddo
close(1)
!
open(1, file=output_dir//'/pdf_phigridf.dat')
   write(1,*) 'phi grid (input)'
   write(1,'(3a20)') 'phi_mid', 'pdf(phi)', 'p(phi)'
   do i=1, nphi-1 
      write(1,'(3es20.8)') phi_mid(i), pdf_phi(i), p_phi(i)
   enddo
close(1)

!stop 'go on in grid_spatial'
!
open(1, file=output_dir//'/pdf_xgridf_c.dat')
   write(1,*) 'optimum core-distribution of x-coordinates (w.r.t. radial grid)'
   write(1,'(3a20)') 'x_mid', 'pdf(x)', 'p(x)'
   do i=1, nxfine_c-1 
      write(1,'(3es20.8)') xfine_mid_c(i), pdf_xfine_c(i), p_xfine_c(i)
   enddo
close(1)
!
open(1, file=output_dir//'/pdf_xgridf_nc.dat')
   write(1,*) 'optimum non-core-distribution of x-coordinates (w.r.t. radial grid)'
   write(1,'(3a20)') 'x_mid', 'pdf(x)', 'p(x)'
   do i=1, nxfine_nc-1 
      write(1,'(3es20.8)') xfine_mid_nc(i), pdf_xfine_nc(i), p_xfine_nc(i)
   enddo
close(1)
!
open(1, file=output_dir//'/pdf_xgridf.dat')
   write(1,*) 'optimum (complete) distribution of x-coordinates w.r.t. radial grid'
   write(1,'(3a20)') 'x_mid', 'pdf(x)', 'p(x)'
   do i=1, nxfine-1 
      write(1,'(3es20.8)') xfine_mid(i), pdf_xfine(i), p_xfine(i)
   enddo
close(1)
!
open(1, file=output_dir//'/pdf_xgrid1.dat')
   write(1,*) 'actually used distribution at first step: gridxyz_optb'
   write(1,'(3a20)') 'x_mid', 'pdf(x)', 'p(x)'
   do i=1, nx_final-1 
      write(1,'(3es20.8)') xfinal_mid(i), pdf_xfinal(i), p_xfinal(i)
   enddo
close(1)
!
!-------------------create complete grid--------------------------------
!
allocate(x(ndxmax), stat=err)
   if(err.ne.0) stop 'error in gridxyz_opt: allocation'
!
x(1)=-xfinal(nx_final)
do i=2, nx_final
   x(i)=-xfinal(nx_final+2-i)
enddo
!
do i=1, nx_final
   x(nx_final + i) = xfinal(i)
enddo
!
x(ndxmax)=xfinal(nx_final)
!
!-------------calculate pdf for z-coordinates on a fine grid------------
!---------------1. pdf for core points (normalized to core)-------------
!---------------2. pdf for non-core points (normlized to non-core)------
!---------------3. complete pdf (normalized to complete range)----------
!
!----------------------------core pdf-----------------------------------
!
nzfine_c=100
allocate(zfine_c(nzfine_c), stat=err)
allocate(zfine_mid_c(nzfine_c-1), stat=err)
allocate(pdf_zfine_c(nzfine_c-1), stat=err)
allocate(p_zfine_c(nzfine_c-1), stat=err)
!
!core grid confined towards rmin
cc=half
grad = rmin/log10((cc+one)/cc)
xinterc = -rmin*log10(cc)/log10((cc+one)/cc)
do i=1, nzfine_c
   zfine_c(i) = grad * log10(float(i)/(nzfine_c)+cc) + xinterc
enddo
zfine_c(1)=zero
zfine_c(nzfine_c)=rmin
!
!create mid points
do i=1, nzfine_c-1
   zfine_mid_c(i) = (zfine_c(i+1)+zfine_c(i))/two
enddo
!
!
!calculate pdf for core coordinates
nx_c=200
nx_nc=400
nx=nx_c+nx_nc
allocate(xdum(nx), stat=err)
!
do i=1, nzfine_c-1
!
!create x-coordinates for x-integration (obtain marignal distribution by integrating over x)
!   write(*,*) zfine_mid_c(i)
   xmin = sqrt(rmin**2-zfine_mid_c(i)**2)
   xmax = sqrt(rmax**2-zfine_mid_c(i)**2)
!   write(*,*) xmin, xmax
!
   cc=half
   grad = (rmin-xmin)/log10((cc+one)/cc)
   xinterc = (xmin*log10(cc+one)-rmin*log10(cc))/log10((cc+one)/cc)
   k=1
   do j=1, nx_c
      xdum(k) = grad * log10(float(j-1)/(nx_c-1)+cc) + xinterc
      k=k+1
   enddo
   del=log10(xmax/rmin)/float(nx_nc)
   do j=1, nx_nc
      xdum(k) = xdum(k-1)*ten**del
      k=k+1
   enddo
!   write(*,*) xdum
!
!calculate pdf for the given x
   integral=zero
   x_jm1 = xdum(1)
!
!calculate pdf for r and phi distribtion of the given (x,z)
   r_jm1 = sqrt(zfine_mid_c(i)**2+xdum(1)**2)
   phi_jm1 = atan(zfine_mid_c(i)/xdum(1))
!   write(*,*) r_jm1, phi_jm1
   call find_index(r_jm1, r1d_dum, n1d_dum, iim2, iim1, ii, iip1)
   pdfr_jm1 = pdf_r1d(iim1)
   call find_index(phi_jm1, phi, nphi, iim2, iim1, ii, iip1)
   pdfp_jm1 = pdf_phi(iim1)
!   write(*,*) pdfr_jm1, pdfp_jm1
!   pdfp_jm1 = two/pi
!calculate joint pdf
   h_jm1 = pdfr_jm1*pdfp_jm1/sqrt(zfine_mid_c(i)**2+x_jm1**2)
!   write(*,*) h_jm1
!
   do j=2, nx
      x_j = xdum(j)
!calculate pdf for r and phi distribtion of the given (x,z)
      r_j = sqrt(zfine_mid_c(i)**2+xdum(j)**2)
      phi_j = atan(zfine_mid_c(i)/xdum(j))
      call find_index(r_j, r1d_dum, n1d_dum, iim2, iim1, ii, iip1)
      pdfr_j = pdf_r1d(iim1)
      call find_index(phi_j, phi, nphi, iim2, iim1, ii, iip1)
      pdfp_j = pdf_phi(iim1)
!      pdfp_j = two/!pi
!calculate joint pdf
      h_j = pdfr_j*pdfp_j/sqrt(zfine_mid_c(i)**2+x_j**2)
!
!perform integration over z
      integral = integral + half*(h_jm1+h_j)*(x_j-x_jm1)
!
      x_jm1 = x_j
      h_jm1 = h_j 
!
   enddo
!
   pdf_zfine_c(i) = integral
!
enddo
!
!
!
!calculate probability of finding x in corresponding bins
do i=1, nzfine_c-1
   p_zfine_c(i) = pdf_zfine_c(i)*(zfine_c(i+1)-zfine_c(i))
!   write(*,*) xfine_mid_c(i), pdf_xfine_c(i)*four, p_xfine_c(i)
enddo
!
prob_c = sum(p_zfine_c)
!
!normalize pdf to core interval
pdf_zfine_c=pdf_zfine_c/prob_c
!
!--------------------------non core pdf---------------------------------
!
nzfine_nc=200
allocate(zfine_nc(nzfine_nc), stat=err)
allocate(zfine_mid_nc(nzfine_nc-1), stat=err)
allocate(pdf_zfine_nc(nzfine_nc-1), stat=err)
allocate(p_zfine_nc(nzfine_nc-1), stat=err)
!
!non-core grid logarithmic
call grid_log(rmin, rmax, nzfine_nc, zfine_nc)
!
!create mid points
do i=1, nzfine_nc-1
   zfine_mid_nc(i) = (zfine_nc(i+1)+zfine_nc(i))/two
enddo
!
!
!
deallocate(xdum)
nx_c=100
nx_nc=100
nx=nx_c+nx_nc
allocate(xdum(nx), stat=err)
!
!calculate pdf for non-core x-coordinates (on half grid)
do i=1, nzfine_nc-1
!
!create x-coordinates for x-integration (obtain marignal distribution by integrating over x)
   xmin = zero
   xmax = sqrt(rmax**2-zfine_mid_nc(i)**2)
!
   cc=half
   grad = (rmin-xmin)/log10((cc+one)/cc)
   xinterc = (xmin*log10(cc+one)-rmin*log10(cc))/log10((cc+one)/cc)
   xdum(1) = xmin
   k=2
   do j=2, nx_c
      xdum(k) = grad * log10(float(j-1)/(nx_c-1)+cc) + xinterc
      k=k+1
   enddo
   del=log10(xmax/rmin)/float(nx_nc)
   do j=1, nx_nc
      xdum(k) = xdum(k-1)*ten**del
      k=k+1
   enddo
!
!
!
!calculate pdf for the given x
   integral=zero
   x_jm1 = xdum(1)
!
!calculate pdf for r and phi distribtion of the given (x,z)
   r_jm1 = sqrt(zfine_mid_nc(i)**2+xdum(1)**2)
!   phi_jm1 = atan(zfine_mid_nc(i)/xdum(1))
   phi_jm1 = pi/two !avoid division by zero
   call find_index(r_jm1, r1d_dum, n1d_dum, iim2, iim1, ii, iip1)
   pdfr_jm1 = pdf_r1d(iim1)
   call find_index(phi_jm1, phi, nphi, iim2, iim1, ii, iip1)
   pdfp_jm1 = pdf_phi(iim1)
!   pdfp_jm1 = two/pi
!calculate joint pdf
   h_jm1 = pdfr_jm1*pdfp_jm1/sqrt(zfine_mid_nc(i)**2+x_jm1**2)
!
   do j=2, nx
      x_j = xdum(j)
!calculate pdf for r and phi distribtion of the given (x,z)
      r_j = sqrt(zfine_mid_nc(i)**2+xdum(j)**2)
      phi_j = atan(zfine_mid_nc(i)/xdum(j))
      call find_index(r_j, r1d_dum, n1d_dum, iim2, iim1, ii, iip1)
      pdfr_j = pdf_r1d(iim1)
      call find_index(phi_j, phi, nphi, iim2, iim1, ii, iip1)
      pdfp_j = pdf_phi(iim1)
!      pdfp_j = two/!pi
!calculate joint pdf
      h_j = pdfr_j*pdfp_j/sqrt(zfine_mid_nc(i)**2+x_j**2)
!
!perform integration over z
      integral = integral + half*(h_jm1+h_j)*(x_j-x_jm1)
!
      x_jm1 = x_j
      h_jm1 = h_j
!
   enddo
!
   pdf_zfine_nc(i) = integral
!
enddo
!
!
!
!calculate probability of finding x in corresponding bins
do i=1, nzfine_nc-1
   p_zfine_nc(i) = pdf_zfine_nc(i)*(zfine_nc(i+1)-zfine_nc(i))
enddo
!
!
prob_nc=sum(p_zfine_nc)
!
!normalize pdf to non-core interval
pdf_zfine_nc=pdf_zfine_nc/prob_nc
!
!------------complete pdf normalized to complete range [0.,rmax]-------------
!-----------------------(only required for checks)---------------------------
!
nzfine = nzfine_c + nzfine_nc -1 
allocate(zfine(nzfine), stat=err)
allocate(zfine_mid(nzfine-1), stat=err)
allocate(pdf_zfine(nzfine-1), stat=err)
allocate(p_zfine(nzfine-1), stat=err)
!
!
k=1
do i=1, nzfine_c
   zfine(k) = zfine_c(i)
   k=k+1
enddo
do i=2, nzfine_nc
   zfine(k) = zfine_nc(i)
   k=k+1
enddo
!
k=1
do i=1, nzfine_c-1
   zfine_mid(k) = zfine_mid_c(i)
   pdf_zfine(k) = pdf_zfine_c(i)*prob_c
   p_zfine(k) = p_zfine_c(i)
   k=k+1
enddo
do i=1, nzfine_nc-1
   zfine_mid(k) = zfine_mid_nc(i)
   pdf_zfine(k) = pdf_zfine_nc(i)*prob_nc
   p_zfine(k) = p_zfine_nc(i)
   k=k+1
enddo
!
pdf_zfine=pdf_zfine/sum(p_zfine)
!
!---------------------create final grid from zero to rmax---------------
!
ndz_dum=ndz-1
nz_final=ndz_dum
ncz_final = ncz
nncz_final=ndz_dum-ncz_final+1   !+1 because rmin point occurrs twice
!
write(*,*) 'probability of finding z-coordinate inside core ', prob_c
write(*,*) 'probability of finding z-coordinate outside core', prob_nc
write(*,*) ' => optimum number of core points to match distribution    ', nint(ndz_dum*prob_c)
write(*,*) ' => optimum number of non-core points to match distribution', ndz_dum-nint(ndz_dum*prob_c)+1
write(*,*)
write(*,*) 'actually used number of core points    ', ncz_final
write(*,*) 'actually used number of non-core points', nncz_final
write(*,*) '# core points / # points    ', float(ncz_final)/float(nncz_final+ncz_final)
write(*,*) '# non-core points / # points', float(nncz_final)/float(nncz_final+ncz_final)
write(*,*)
!
!zdum only up to ndx-1 since phantom point will be calculated later on
allocate(zfinal_c(ncz_final), stat=err)
allocate(zfinal_nc(nncz_final), stat=err)
allocate(zfinal(nz_final), stat=err)
!
!-------------------distribute nxz_final grid points in core------------
!-------------------------(following the core pdf)----------------------
!
!probability that has to be matched for each interval
pc = one/float(ncz_final-1)
zfinal_c(1)=zero
!
pc_i = zero
do i=2, ncz_final
!probability that has to be matched
   pc_i = pc_i + pc
!
   delp = zero
   do j=1, nzfine_c-1
      zi =  (pc_i-delp)/pdf_zfine_c(j)  + zfine_c(j)
      zfinal_c(i) = zi
      if(zi.le.zfine_c(j+1)) then
!correct position found: go to next z_final-coordinate
         exit
      else
!go to next interval
         delp = delp + pdf_zfine_c(j)*(zfine_c(j+1)-zfine_c(j))
      endif
   enddo
enddo
!
!---------------distribute nncz_final grid points outside core----------
!-----------------------(following the non-core pdf)--------------------
!
!probability that has to be matched for each interval
pc = one/float(nncz_final-1)
zfinal_nc(1)=zfine_nc(1)
!
pc_i = zero
do i=2, nncz_final
!probability that has to be matched
   pc_i = pc_i + pc
!
   delp = zero
   do j=1, nzfine_nc-1
      zi =  (pc_i-delp)/pdf_zfine_nc(j)  + zfine_nc(j)
      zfinal_nc(i) = zi
      if(zi.le. zfine_nc(j+1)) then
!;correct position found: go to next z_final-coordinate
         exit
      else
!go to next interval
         delp = delp + pdf_zfine_nc(j)*(zfine_nc(j+1)-zfine_nc(j))
      endif
   enddo
enddo
!
!----------------------combined distribution----------------------------
!
k=1
do i=1, ncz_final
   zfinal(k) = zfinal_c(i)
   k=k+1
enddo
zfinal_c(ncz_final)=rmin
do i=1, nncz_final-1
   zfinal(k) = zfinal_nc(i+1)
   k=k+1
enddo
zfinal(nz_final)=rmax
!
!calculate final pdf (only for checks)
allocate(zfinal_mid(nz_final-1))
allocate(pdf_zfinal(nz_final-1))
allocate(p_zfinal(nz_final-1))
!
call calc_pdf(nz_final, zfinal, zfinal_mid, pdf_zfinal, p_zfinal)
!
!--------------------------print out everything-------------------------
!
open(1, file=output_dir//'/pdf_zgridf_c.dat')
   write(1,*) 'optimum core-distribution of z-coordinates (w.r.t. radial grid)'
   write(1,'(3a20)') 'z_mid', 'pdf(z)', 'p(z)'
   do i=1, nzfine_c-1 
      write(1,'(3es20.8)') zfine_mid_c(i), pdf_zfine_c(i), p_zfine_c(i)
   enddo
close(1)
!
open(1, file=output_dir//'/pdf_zgridf_nc.dat')
   write(1,*) 'optimum non-core-distribution of z-coordinates (w.r.t. radial grid)'
   write(1,'(3a20)') 'z_mid', 'pdf(z)', 'p(z)'
   do i=1, nzfine_nc-1 
      write(1,'(3es20.8)') zfine_mid_nc(i), pdf_zfine_nc(i), p_zfine_nc(i)
   enddo
close(1)
!
open(1, file=output_dir//'/pdf_zgridf.dat')
   write(1,*) 'optimum (complete) distribution of z-coordinates w.r.t. radial grid'
   write(1,'(3a20)') 'z_mid', 'pdf(z)', 'p(z)'
   do i=1, nzfine-1 
      write(1,'(3es20.8)') zfine_mid(i), pdf_zfine(i), p_zfine(i)
   enddo
close(1)
!
open(1, file=output_dir//'/pdf_zgrid1.dat')
   write(1,*) 'actually used distribution at first step: gridxyz_optb'
   write(1,'(3a20)') 'z_mid', 'pdf(z)', 'p(z)'
   do i=1, nz_final-1 
      write(1,'(3es20.8)') zfinal_mid(i), pdf_zfinal(i), p_zfinal(i)
   enddo
close(1)
!
!-------------------create complete z grid------------------------------
!
allocate(z(ndzmax), stat=err)
   if(err.ne.0) stop 'error in gridxyz_opt: allocation'
!
z(1)=-zfinal(nz_final)
do i=2, nz_final
   z(i)=-zfinal(nz_final+2-i)
enddo
!
do i=1, nz_final
   z(nz_final + i) = zfinal(i)
enddo
!
z(ndzmax)=zfinal(nz_final)
!
!-------------------------------create y grid---------------------------
!
!!check if dimensions of x and z grid are compatible
!if(ndzmax.ne.ndxmax.or.ndz.ne.ndx) stop 'dimensions of z and x not compatible'
!!set x equal to z-grid
!z=x
!
allocate(y(ndymax), stat=err)
   if(err.ne.0) stop 'error in gridxyz_opt: allocation'
!check if dimensions of y and x grid are compatible
if(ndymax.ne.ndxmax.or.ndy.ne.ndx) stop 'dimensions of y and x not compatible'
!set y equal to x-grid
y=x
!
!deallocate 1d arrays, since not needed anymore
deallocate(r1d)
deallocate(r1d_dum)
!
!write(*,*) x
!stop 'go on in gridxyz_optb'
!
!
end subroutine gridxyz_optb
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine gridxyz_opt1d
!
!reads in radial grid from input model
!calculates a (not necessarily) uniform angle-grid
!
!calculates the probability density of the radial grid,
!and corresponding probability density of the x,y,z-coordinates
!
!calculates the actualy x,y,z-coordinates following this distrubution
!   (with input nc_final, nnc_final: number of points inside core and outside core)
!
!
use prog_type
use fund_const
use mod_directories, only: output_dir
use dime1d, only: n1d_dum, r1d_dum
use dime3d, only: ncx, ndx, ncz, ndz, ndy, ndxmax, ndymax, ndzmax, x, y, z
use params_stellar, only: sr
use params_input, only: rmax
use modext, only: nr_modext, r_modext
use mod_interp1d, only: find_index
use mod_grid, only: calc_pdf, grid_log
!
implicit none
!
! ... local scalars
integer(i4b) :: i, j, k, err
integer(i4b) :: iim2, iim1, ii, iip1
integer(i4b) :: nphi, nx, nx_c, nx_nc, nz, nz_c, nz_nc, &
                nxfine_c, nxfine_nc, nxfine, ndx_dum, nncx_final, nx_final, ncx_final, &
                nzfine_c, nzfine_nc, nzfine, ndz_dum, nncz_final, nz_final, ncz_final
real(dp) :: integral, cc, grad, xinterc, zinterc, del, rmin, xmin, xmax, zmin, zmax, &
            x_jm1, x_j, z_jm1, z_j, r_jm1, r_j, phi_jm1, phi_j, pdfr_jm1, pdfr_j, &
            pdfp_jm1, pdfp_j, h_jm1, h_j, prob_nc, prob_c, pc, pc_i, xi, zi, delp, &
            delx1, maxx, max_delx
!
! ... local arrays
!real(dp), dimension(:), allocatable :: xcoords, ycoords, zcoords, xcoords_dum, ycoords_dum, zcoords_dum
real(dp), dimension(:), allocatable :: phi, phidum, phi_mid, pdf_phi, p_phi, &
                                       r1d_mid, pdf_r1d, p_r1d, &
                                       xfine_c, xfine_mid_c, pdf_xfine_c, p_xfine_c, &
                                       zfine_c, zfine_mid_c, pdf_zfine_c, p_zfine_c, &
                                       xfine_nc, xfine_mid_nc, pdf_xfine_nc, p_xfine_nc, &
                                       zfine_nc, zfine_mid_nc, pdf_zfine_nc, p_zfine_nc, &
                                       xfine, xfine_mid, pdf_xfine, p_xfine, &
                                       zfine, zfine_mid, pdf_zfine, p_zfine
real(dp), dimension(:), allocatable :: xdum, xfinal, xfinal_c, xfinal_nc, xfinal_mid, pdf_xfinal, p_xfinal, &
                                       zdum, zfinal, zfinal_c, zfinal_nc, zfinal_mid, pdf_zfinal, p_zfinal
!
!----------------------set up radial grid-------------------------------
!
call read_mod1d
!
n1d_dum=nr_modext
allocate(r1d_dum(n1d_dum), stat=err)
   if(err.ne.0) stop 'error in gridxyz_opt1d: allocation'
r1d_dum=r_modext/sr
!
rmin=r1d_dum(1)
!rmax=r1d_dum(n1d_dum)
!
allocate(r1d_mid(n1d_dum-1), stat=err)
   if(err.ne.0) stop 'error in gridxyz_opt1d: allocation'
allocate(pdf_r1d(n1d_dum-1), stat=err)
   if(err.ne.0) stop 'error in gridxyz_opt1d: allocation'
allocate(p_r1d(n1d_dum-1), stat=err)
   if(err.ne.0) stop 'error in gridxyz_opt1d: allocation'
!
!calculate probability density
call calc_pdf(n1d_dum, r1d_dum, r1d_mid, pdf_r1d, p_r1d)
!
!---------------create dummy phi grid (measured from x to z-axis)-------
!phi is in interval [0,pi/2]
!
nphi=45
!
allocate(phi(nphi), stat=err)
   if(err.ne.0) stop 'error in gridxyz_opt1d: allocation'
allocate(phidum(nphi), stat=err)
   if(err.ne.0) stop 'error in gridxyz_opt1d: allocation'
allocate(phi_mid(nphi-1), stat=err)
   if(err.ne.0) stop 'error in gridxyz_opt1d: allocation'
allocate(pdf_phi(nphi-1), stat=err)
   if(err.ne.0) stop 'error in gridxyz_opt1d: allocation'
allocate(p_phi(nphi-1), stat=err)
   if(err.ne.0) stop 'error in gridxyz_opt1d: allocation'
!
!phi grid equidistant
do i=1, nphi
   phi(i) = float(i-1) * pi / float(nphi-1) / two
enddo
!
!phi grid with high resolution at small angles
!call grid_log(1.d-3,pi/two,nphi,phi)
!!
!!phi grid with high resolution at large angles
!phidum(nphi)=pi/two
!k=2
!do i=nphi-1, 1, -1
!   phidum(i) = phidum(i+1)-(phi(k)-phi(k-1))
!   k=k+1
!enddo
!phi=phidum
!
!calculate probability density
call calc_pdf(nphi, phi, phi_mid, pdf_phi, p_phi)
!
!-------------calculate pdf for x-coordinates on a fine grid------------
!---------------1. pdf for core points (normalized to core)-------------
!---------------2. pdf for non-core points (normlized to non-core)------
!---------------3. complete pdf (normalized to complete range)----------
!
!----------------------------core pdf-----------------------------------
!
nxfine_c=100
allocate(xfine_c(nxfine_c), stat=err)
allocate(xfine_mid_c(nxfine_c-1), stat=err)
allocate(pdf_xfine_c(nxfine_c-1), stat=err)
allocate(p_xfine_c(nxfine_c-1), stat=err)
!
!core grid confined towards rmin
cc=half
grad = rmin/log10((cc+one)/cc)
zinterc = -rmin*log10(cc)/log10((cc+one)/cc)
do i=1, nxfine_c
   xfine_c(i) = grad * log10(float(i)/(nxfine_c)+cc) + zinterc
enddo
xfine_c(1)=zero
xfine_c(nxfine_c)=rmin
!
!
!create mid points
do i=1, nxfine_c-1
   xfine_mid_c(i) = (xfine_c(i+1)+xfine_c(i))/two
enddo
!
!
!calculate pdf for core coordinates
nz_c=200
nz_nc=400
nz=nz_c+nz_nc
allocate(zdum(nz), stat=err)
!
do i=1, nxfine_c-1
!
!create z-coordinates for z-integration (obtain marignal distribution by integrating over z)
   zmin = sqrt(rmin**2-xfine_mid_c(i)**2)
   zmax = sqrt(rmax**2-xfine_mid_c(i)**2)
!
   cc=half
   grad = (rmin-zmin)/log10((cc+one)/cc)
   zinterc = (zmin*log10(cc+one)-rmin*log10(cc))/log10((cc+one)/cc)
   k=1
   do j=1, nz_c
      zdum(k) = grad * log10(float(j-1)/(nz_c-1)+cc) + zinterc
      k=k+1
   enddo
   del=log10(zmax/rmin)/float(nz_nc)
   do j=1, nz_nc
      zdum(k) = zdum(k-1)*ten**del
      k=k+1
   enddo
!
!calculate pdf for the given x
   integral=zero
   z_jm1 = zdum(1)
!
!calculate pdf for r and phi distribtion of the given (x,z)
   r_jm1 = sqrt(xfine_mid_c(i)**2+zdum(1)**2)
   phi_jm1 = atan(zdum(1)/xfine_mid_c(i))
   call find_index(r_jm1, r1d_dum, n1d_dum, iim2, iim1, ii, iip1)
   pdfr_jm1 = pdf_r1d(iim1)
   call find_index(phi_jm1, phi, nphi, iim2, iim1, ii, iip1)
   pdfp_jm1 = pdf_phi(iim1)
!   pdfp_jm1 = two/pi
!calculate joint pdf
   h_jm1 = pdfr_jm1*pdfp_jm1/sqrt(xfine_mid_c(i)**2+z_jm1**2)
!
   do j=2, nz
      z_j = zdum(j)
!calculate pdf for r and phi distribtion of the given (x,z)
      r_j = sqrt(xfine_mid_c(i)**2+zdum(j)**2)
      phi_j = atan(zdum(j)/xfine_mid_c(i))
      call find_index(r_j, r1d_dum, n1d_dum, iim2, iim1, ii, iip1)
      pdfr_j = pdf_r1d(iim1)
      call find_index(phi_j, phi, nphi, iim2, iim1, ii, iip1)
      pdfp_j = pdf_phi(iim1)
!      pdfp_j = two/!pi
!calculate joint pdf
      h_j = pdfr_j*pdfp_j/sqrt(xfine_mid_c(i)**2+z_j**2)
!
!perform integration over z
      integral = integral + half*(h_jm1+h_j)*(z_j-z_jm1)
!
      z_jm1 = z_j
      h_jm1 = h_j
!
   enddo
!
   pdf_xfine_c(i) = integral
!
enddo
!
!
!
!calculate probability of finding x in corresponding bins
do i=1, nxfine_c-1
   p_xfine_c(i) = pdf_xfine_c(i)*(xfine_c(i+1)-xfine_c(i))
!   write(*,*) xfine_mid_c(i), pdf_xfine_c(i)*four, p_xfine_c(i)
enddo
!
prob_c = sum(p_xfine_c)
!
!normalize pdf to core interval
pdf_xfine_c=pdf_xfine_c/prob_c
!
!--------------------------non core pdf---------------------------------
!
nxfine_nc=200
allocate(xfine_nc(nxfine_nc), stat=err)
allocate(xfine_mid_nc(nxfine_nc-1), stat=err)
allocate(pdf_xfine_nc(nxfine_nc-1), stat=err)
allocate(p_xfine_nc(nxfine_nc-1), stat=err)
!
!non-core grid logarithmic
call grid_log(rmin, rmax, nxfine_nc, xfine_nc)
!
!create mid points
do i=1, nxfine_nc-1
   xfine_mid_nc(i) = (xfine_nc(i+1)+xfine_nc(i))/two
enddo
!
!
!
deallocate(zdum)
nz_c=100
nz_nc=100
nz=nz_c+nz_nc
allocate(zdum(nz), stat=err)
!
!calculate pdf for non-core x-coordinates (on half grid)
do i=1, nxfine_nc-1
!
!create z-coordinates for z-integration (obtain marignal distribution by integrating over z)
   zmin = zero
   zmax = sqrt(rmax**2-xfine_mid_nc(i)**2)
!
   cc=half
   grad = (rmin-zmin)/log10((cc+one)/cc)
   zinterc = (zmin*log10(cc+one)-rmin*log10(cc))/log10((cc+one)/cc)
   zdum(1) = zmin
   k=2
   do j=2, nz_c
      zdum(k) = grad * log10(float(j-1)/(nz_c-1)+cc) + zinterc
      k=k+1
   enddo
   del=log10(zmax/rmin)/float(nz_nc)
   do j=1, nz_nc
      zdum(k) = zdum(k-1)*ten**del
      k=k+1
   enddo
!
!
!calculate pdf for the given x
   integral=zero
   z_jm1 = zdum(1)
!
!calculate pdf for r and phi distribtion of the given (x,z)
   r_jm1 = sqrt(xfine_mid_nc(i)**2+zdum(1)**2)
   phi_jm1 = atan(zdum(1)/xfine_mid_nc(i))
   call find_index(r_jm1, r1d_dum, n1d_dum, iim2, iim1, ii, iip1)
   pdfr_jm1 = pdf_r1d(iim1)
   call find_index(phi_jm1, phi, nphi, iim2, iim1, ii, iip1)
   pdfp_jm1 = pdf_phi(iim1)
!   pdfp_jm1 = two/pi
!calculate joint pdf
   h_jm1 = pdfr_jm1*pdfp_jm1/sqrt(xfine_mid_nc(i)**2+z_jm1**2)
!
   do j=2, nz
      z_j = zdum(j)
!calculate pdf for r and phi distribtion of the given (x,z)
      r_j = sqrt(xfine_mid_nc(i)**2+zdum(j)**2)
      phi_j = atan(zdum(j)/xfine_mid_nc(i))
      call find_index(r_j, r1d_dum, n1d_dum, iim2, iim1, ii, iip1)
      pdfr_j = pdf_r1d(iim1)
      call find_index(phi_j, phi, nphi, iim2, iim1, ii, iip1)
      pdfp_j = pdf_phi(iim1)
!      pdfp_j = two/!pi
!calculate joint pdf
      h_j = pdfr_j*pdfp_j/sqrt(xfine_mid_nc(i)**2+z_j**2)
!
!perform integration over z
      integral = integral + half*(h_jm1+h_j)*(z_j-z_jm1)
!
      z_jm1 = z_j
      h_jm1 = h_j
!
   enddo
!
   pdf_xfine_nc(i) = integral
!
enddo
!
!
!
!calculate probability of finding x in corresponding bins
do i=1, nxfine_nc-1
   p_xfine_nc(i) = pdf_xfine_nc(i)*(xfine_nc(i+1)-xfine_nc(i))
enddo
!
!
prob_nc=sum(p_xfine_nc)
!
!normalize pdf to non-core interval
pdf_xfine_nc=pdf_xfine_nc/prob_nc
!
!------------complete pdf normalized to complete range [0.,rmax]-------------
!-----------------------(only required for checks)---------------------------
!
nxfine = nxfine_c + nxfine_nc -1 
allocate(xfine(nxfine), stat=err)
allocate(xfine_mid(nxfine-1), stat=err)
allocate(pdf_xfine(nxfine-1), stat=err)
allocate(p_xfine(nxfine-1), stat=err)
!
!
k=1
do i=1, nxfine_c
   xfine(k) = xfine_c(i)
   k=k+1
enddo
do i=2, nxfine_nc
   xfine(k) = xfine_nc(i)
   k=k+1
enddo
!
k=1
do i=1, nxfine_c-1
   xfine_mid(k) = xfine_mid_c(i)
   pdf_xfine(k) = pdf_xfine_c(i)*prob_c
   p_xfine(k) = p_xfine_c(i)
   k=k+1
enddo
do i=1, nxfine_nc-1
   xfine_mid(k) = xfine_mid_nc(i)
   pdf_xfine(k) = pdf_xfine_nc(i)*prob_nc
   p_xfine(k) = p_xfine_nc(i)
   k=k+1
enddo
!
pdf_xfine=pdf_xfine/sum(p_xfine)
!
!---------------------create final grid from zero to rmax---------------
!
!ndx_dum=ndx-1
!nx_final=ndx_dum
!nncx_final=ndx_dum-ncx_final+1   !+1 because rmin point occurrs twice
!write(*,*) nx_final, ncx_final, nncx_final
!
ndx_dum=ndx-1
nx_final=ndx_dum
ncx_final = ncx
nncx_final=ndx_dum-ncx_final+1   !+1 because rmin point occurrs twice
!
write(*,*) 'probability of finding x-coordinate inside core ', prob_c
write(*,*) 'probability of finding x-coordinate outside core', prob_nc
write(*,*) ' => optimum number of core points to match distribution    ', nint(ndx_dum*prob_c)
write(*,*) ' => optimum number of non-core points to match distribution', ndx_dum-nint(ndx_dum*prob_c)+1
write(*,*)
write(*,*) 'actually used number of core points    ', ncx_final
write(*,*) 'actually used number of non-core points', nncx_final
write(*,*) '# core points / # points    ', float(ncx_final)/float(nncx_final+ncx_final)
write(*,*) '# non-core points / # points', float(nncx_final)/float(nncx_final+ncx_final)
write(*,*)
!
!xdum only up to ndx-1 since phantom point will be calculated later on
allocate(xfinal_c(ncx_final), stat=err)
allocate(xfinal_nc(nncx_final), stat=err)
allocate(xfinal(nx_final), stat=err)
!
!-------------------distribute nxc_final grid points in core------------
!-------------------------(following the core pdf)----------------------
!
!probability that has to be matched for each interval
pc = one/float(ncx_final-1)
xfinal_c(1)=zero
!
pc_i = zero
do i=2, ncx_final
!probability that has to be matched
   pc_i = pc_i + pc
!
   delp = zero
   do j=1, nxfine_c-1
      xi =  (pc_i-delp)/pdf_xfine_c(j)  + xfine_c(j)
      xfinal_c(i) = xi
      if(xi.le.xfine_c(j+1)) then
!correct position found: go to next x_final-coordinate
         exit
      else
!go to next interval
         delp = delp + pdf_xfine_c(j)*(xfine_c(j+1)-xfine_c(j))
      endif
   enddo
enddo
!
!---------------distribute nncx_final grid points outside core----------
!-----------------------(following the non-core pdf)--------------------
!
!probability that has to be matched for each interval
pc = one/float(nncx_final-1)
xfinal_nc(1)=xfine_nc(1)
!
pc_i = zero
do i=2, nncx_final
!probability that has to be matched
   pc_i = pc_i + pc
!
   delp = zero
   do j=1, nxfine_nc-1
      xi =  (pc_i-delp)/pdf_xfine_nc(j)  + xfine_nc(j)
      xfinal_nc(i) = xi
      if(xi.le. xfine_nc(j+1)) then
!;correct position found: go to next x_final-coordinate
         exit
      else
!go to next interval
         delp = delp + pdf_xfine_nc(j)*(xfine_nc(j+1)-xfine_nc(j))
      endif
   enddo
enddo
!
!----------------------combined distribution----------------------------
!
k=1
do i=1, ncx_final
   xfinal(k) = xfinal_c(i)
   k=k+1
enddo
xfinal_c(ncx_final)=rmin
do i=1, nncx_final-1
   xfinal(k) = xfinal_nc(i+1)
   k=k+1
enddo
xfinal(nx_final)=rmax
!
!calculate final pdf (only for checks)
allocate(xfinal_mid(nx_final-1))
allocate(pdf_xfinal(nx_final-1))
allocate(p_xfinal(nx_final-1))
!
call calc_pdf(nx_final, xfinal, xfinal_mid, pdf_xfinal, p_xfinal)
!
!--------------------------print out everything-------------------------
!
open(1, file=output_dir//'/pdf_rgridf.dat')
   write(1,*) 'radial grid (input)'
   write(1,'(3a20)') 'r_mid', 'pdf(r)', 'p(r)'
   do i=1, n1d_dum-1 
      write(1,'(3es20.8)') r1d_mid(i), pdf_r1d(i), p_r1d(i)
   enddo
close(1)
!
open(1, file=output_dir//'/pdf_phigridf.dat')
   write(1,*) 'phi grid (input)'
   write(1,'(3a20)') 'phi_mid', 'pdf(phi)', 'p(phi)'
   do i=1, nphi-1 
      write(1,'(3es20.8)') phi_mid(i), pdf_phi(i), p_phi(i)
   enddo
close(1)

!stop 'go on in grid_spatial'
!
open(1, file=output_dir//'/pdf_xgridf_c.dat')
   write(1,*) 'optimum core-distribution of x-coordinates (w.r.t. radial grid)'
   write(1,'(3a20)') 'x_mid', 'pdf(x)', 'p(x)'
   do i=1, nxfine_c-1 
      write(1,'(3es20.8)') xfine_mid_c(i), pdf_xfine_c(i), p_xfine_c(i)
   enddo
close(1)
!
open(1, file=output_dir//'/pdf_xgridf_nc.dat')
   write(1,*) 'optimum non-core-distribution of x-coordinates (w.r.t. radial grid)'
   write(1,'(3a20)') 'x_mid', 'pdf(x)', 'p(x)'
   do i=1, nxfine_nc-1 
      write(1,'(3es20.8)') xfine_mid_nc(i), pdf_xfine_nc(i), p_xfine_nc(i)
   enddo
close(1)
!
open(1, file=output_dir//'/pdf_xgridf.dat')
   write(1,*) 'optimum (complete) distribution of x-coordinates w.r.t. radial grid'
   write(1,'(3a20)') 'x_mid', 'pdf(x)', 'p(x)'
   do i=1, nxfine-1 
      write(1,'(3es20.8)') xfine_mid(i), pdf_xfine(i), p_xfine(i)
   enddo
close(1)
!
open(1, file=output_dir//'/pdf_xgrid1.dat')
   write(1,*) 'actually used distribution at first step: gridxyz_optb'
   write(1,'(3a20)') 'x_mid', 'pdf(x)', 'p(x)'
   do i=1, nx_final-1 
      write(1,'(3es20.8)') xfinal_mid(i), pdf_xfinal(i), p_xfinal(i)
   enddo
close(1)
!
!-------------------create complete grid--------------------------------
!
allocate(x(ndxmax), stat=err)
   if(err.ne.0) stop 'error in gridxyz_opt: allocation'
!
x(1)=-xfinal(nx_final)
do i=2, nx_final
   x(i)=-xfinal(nx_final+2-i)
enddo
!
do i=1, nx_final
   x(nx_final + i) = xfinal(i)
enddo
!
x(ndxmax)=xfinal(nx_final)
!
!-------------calculate pdf for z-coordinates on a fine grid------------
!---------------1. pdf for core points (normalized to core)-------------
!---------------2. pdf for non-core points (normlized to non-core)------
!---------------3. complete pdf (normalized to complete range)----------
!
!----------------------------core pdf-----------------------------------
!
nzfine_c=100
allocate(zfine_c(nzfine_c), stat=err)
allocate(zfine_mid_c(nzfine_c-1), stat=err)
allocate(pdf_zfine_c(nzfine_c-1), stat=err)
allocate(p_zfine_c(nzfine_c-1), stat=err)
!
!core grid confined towards rmin
cc=half
grad = rmin/log10((cc+one)/cc)
xinterc = -rmin*log10(cc)/log10((cc+one)/cc)
do i=1, nzfine_c
   zfine_c(i) = grad * log10(float(i)/(nzfine_c)+cc) + xinterc
enddo
zfine_c(1)=zero
zfine_c(nzfine_c)=rmin
!
!create mid points
do i=1, nzfine_c-1
   zfine_mid_c(i) = (zfine_c(i+1)+zfine_c(i))/two
enddo
!
!
!calculate pdf for core coordinates
nx_c=200
nx_nc=400
nx=nx_c+nx_nc
allocate(xdum(nx), stat=err)
!
do i=1, nzfine_c-1
!
!create x-coordinates for x-integration (obtain marignal distribution by integrating over x)
   xmin = sqrt(rmin**2-zfine_mid_c(i)**2)
   xmax = sqrt(rmax**2-zfine_mid_c(i)**2)
!
   cc=half
   grad = (rmin-xmin)/log10((cc+one)/cc)
   xinterc = (xmin*log10(cc+one)-rmin*log10(cc))/log10((cc+one)/cc)
   k=1
   do j=1, nx_c
      xdum(k) = grad * log10(float(j-1)/(nx_c-1)+cc) + xinterc
      k=k+1
   enddo
   del=log10(xmax/rmin)/float(nx_nc)
   do j=1, nx_nc
      xdum(k) = xdum(k-1)*ten**del
      k=k+1
   enddo
!
!calculate pdf for the given x
   integral=zero
   x_jm1 = xdum(1)
!
!calculate pdf for r and phi distribtion of the given (x,z)
   r_jm1 = sqrt(zfine_mid_c(i)**2+xdum(1)**2)
   phi_jm1 = atan(zfine_mid_c(i)/xdum(1))
!   write(*,*) r_jm1, phi_jm1
   call find_index(r_jm1, r1d_dum, n1d_dum, iim2, iim1, ii, iip1)
   pdfr_jm1 = pdf_r1d(iim1)
   call find_index(phi_jm1, phi, nphi, iim2, iim1, ii, iip1)
   pdfp_jm1 = pdf_phi(iim1)
!   write(*,*) pdfr_jm1, pdfp_jm1
!   pdfp_jm1 = two/pi
!calculate joint pdf
   h_jm1 = pdfr_jm1*pdfp_jm1/sqrt(zfine_mid_c(i)**2+x_jm1**2)
!   write(*,*) h_jm1
!
   do j=2, nx
      x_j = xdum(j)
!calculate pdf for r and phi distribtion of the given (x,z)
      r_j = sqrt(zfine_mid_c(i)**2+xdum(j)**2)
      phi_j = atan(zfine_mid_c(i)/xdum(j))
      call find_index(r_j, r1d_dum, n1d_dum, iim2, iim1, ii, iip1)
      pdfr_j = pdf_r1d(iim1)
      call find_index(phi_j, phi, nphi, iim2, iim1, ii, iip1)
      pdfp_j = pdf_phi(iim1)
!      pdfp_j = two/!pi
!calculate joint pdf
      h_j = pdfr_j*pdfp_j/sqrt(zfine_mid_c(i)**2+x_j**2)
!
!perform integration over z
      integral = integral + half*(h_jm1+h_j)*(x_j-x_jm1)
!
      x_jm1 = x_j
      h_jm1 = h_j 
!
   enddo
!
   pdf_zfine_c(i) = integral
!
enddo
!
!
!
!calculate probability of finding x in corresponding bins
do i=1, nzfine_c-1
   p_zfine_c(i) = pdf_zfine_c(i)*(zfine_c(i+1)-zfine_c(i))
!   write(*,*) xfine_mid_c(i), pdf_xfine_c(i)*four, p_xfine_c(i)
enddo
!
prob_c = sum(p_zfine_c)
!
!normalize pdf to core interval
pdf_zfine_c=pdf_zfine_c/prob_c
!
!--------------------------non core pdf---------------------------------
!
nzfine_nc=200
allocate(zfine_nc(nzfine_nc), stat=err)
allocate(zfine_mid_nc(nzfine_nc-1), stat=err)
allocate(pdf_zfine_nc(nzfine_nc-1), stat=err)
allocate(p_zfine_nc(nzfine_nc-1), stat=err)
!
!non-core grid logarithmic
call grid_log(rmin, rmax, nzfine_nc, zfine_nc)
!
!create mid points
do i=1, nzfine_nc-1
   zfine_mid_nc(i) = (zfine_nc(i+1)+zfine_nc(i))/two
enddo
!
!
!
deallocate(xdum)
nx_c=100
nx_nc=100
nx=nx_c+nx_nc
allocate(xdum(nx), stat=err)
!
!calculate pdf for non-core x-coordinates (on half grid)
do i=1, nzfine_nc-1
!
!create x-coordinates for x-integration (obtain marignal distribution by integrating over x)
   xmin = zero
   xmax = sqrt(rmax**2-zfine_mid_nc(i)**2)
!
   cc=half
   grad = (rmin-xmin)/log10((cc+one)/cc)
   xinterc = (xmin*log10(cc+one)-rmin*log10(cc))/log10((cc+one)/cc)
   xdum(1) = xmin
   k=2
   do j=2, nx_c
      xdum(k) = grad * log10(float(j-1)/(nx_c-1)+cc) + xinterc
      k=k+1
   enddo
   del=log10(xmax/rmin)/float(nx_nc)
   do j=1, nx_nc
      xdum(k) = xdum(k-1)*ten**del
      k=k+1
   enddo
!
!
!
!calculate pdf for the given x
   integral=zero
   x_jm1 = xdum(1)
!
!calculate pdf for r and phi distribtion of the given (x,z)
   r_jm1 = sqrt(zfine_mid_nc(i)**2+xdum(1)**2)
!   phi_jm1 = atan(zfine_mid_nc(i)/xdum(1))
   phi_jm1 = pi/two !avoid division by zero
   call find_index(r_jm1, r1d_dum, n1d_dum, iim2, iim1, ii, iip1)
   pdfr_jm1 = pdf_r1d(iim1)
   call find_index(phi_jm1, phi, nphi, iim2, iim1, ii, iip1)
   pdfp_jm1 = pdf_phi(iim1)
!   pdfp_jm1 = two/pi
!calculate joint pdf
   h_jm1 = pdfr_jm1*pdfp_jm1/sqrt(zfine_mid_nc(i)**2+x_jm1**2)
!
   do j=2, nx
      x_j = xdum(j)
!calculate pdf for r and phi distribtion of the given (x,z)
      r_j = sqrt(zfine_mid_nc(i)**2+xdum(j)**2)
      phi_j = atan(zfine_mid_nc(i)/xdum(j))
      call find_index(r_j, r1d_dum, n1d_dum, iim2, iim1, ii, iip1)
      pdfr_j = pdf_r1d(iim1)
      call find_index(phi_j, phi, nphi, iim2, iim1, ii, iip1)
      pdfp_j = pdf_phi(iim1)
!      pdfp_j = two/!pi
!calculate joint pdf
      h_j = pdfr_j*pdfp_j/sqrt(zfine_mid_nc(i)**2+x_j**2)
!
!perform integration over z
      integral = integral + half*(h_jm1+h_j)*(x_j-x_jm1)
!
      x_jm1 = x_j
      h_jm1 = h_j
!
   enddo
!
   pdf_zfine_nc(i) = integral
!
enddo
!
!
!
!calculate probability of finding x in corresponding bins
do i=1, nzfine_nc-1
   p_zfine_nc(i) = pdf_zfine_nc(i)*(zfine_nc(i+1)-zfine_nc(i))
enddo
!
!
prob_nc=sum(p_zfine_nc)
!
!normalize pdf to non-core interval
pdf_zfine_nc=pdf_zfine_nc/prob_nc
!
!------------complete pdf normalized to complete range [0.,rmax]-------------
!-----------------------(only required for checks)---------------------------
!
nzfine = nzfine_c + nzfine_nc -1 
allocate(zfine(nzfine), stat=err)
allocate(zfine_mid(nzfine-1), stat=err)
allocate(pdf_zfine(nzfine-1), stat=err)
allocate(p_zfine(nzfine-1), stat=err)
!
!
k=1
do i=1, nzfine_c
   zfine(k) = zfine_c(i)
   k=k+1
enddo
do i=2, nzfine_nc
   zfine(k) = zfine_nc(i)
   k=k+1
enddo
!
k=1
do i=1, nzfine_c-1
   zfine_mid(k) = zfine_mid_c(i)
   pdf_zfine(k) = pdf_zfine_c(i)*prob_c
   p_zfine(k) = p_zfine_c(i)
   k=k+1
enddo
do i=1, nzfine_nc-1
   zfine_mid(k) = zfine_mid_nc(i)
   pdf_zfine(k) = pdf_zfine_nc(i)*prob_nc
   p_zfine(k) = p_zfine_nc(i)
   k=k+1
enddo
!
pdf_zfine=pdf_zfine/sum(p_zfine)
!
!---------------------create final grid from zero to rmax---------------
!
ndz_dum=ndz-1
nz_final=ndz_dum
ncz_final = ncz
nncz_final=ndz_dum-ncz_final+1   !+1 because rmin point occurrs twice
!
write(*,*) 'probability of finding z-coordinate inside core ', prob_c
write(*,*) 'probability of finding z-coordinate outside core', prob_nc
write(*,*) ' => optimum number of core points to match distribution    ', nint(ndz_dum*prob_c)
write(*,*) ' => optimum number of non-core points to match distribution', ndz_dum-nint(ndz_dum*prob_c)+1
write(*,*)
write(*,*) 'actually used number of core points    ', ncz_final
write(*,*) 'actually used number of non-core points', nncz_final
write(*,*) '# core points / # points    ', float(ncz_final)/float(nncz_final+ncz_final)
write(*,*) '# non-core points / # points', float(nncz_final)/float(nncz_final+ncz_final)
write(*,*)
!
!zdum only up to ndx-1 since phantom point will be calculated later on
allocate(zfinal_c(ncz_final), stat=err)
allocate(zfinal_nc(nncz_final), stat=err)
allocate(zfinal(nz_final), stat=err)
!
!-------------------distribute nxz_final grid points in core------------
!-------------------------(following the core pdf)----------------------
!
!probability that has to be matched for each interval
pc = one/float(ncz_final-1)
zfinal_c(1)=zero
!
pc_i = zero
do i=2, ncz_final
!probability that has to be matched
   pc_i = pc_i + pc
!
   delp = zero
   do j=1, nzfine_c-1
      zi =  (pc_i-delp)/pdf_zfine_c(j)  + zfine_c(j)
      zfinal_c(i) = zi
      if(zi.le.zfine_c(j+1)) then
!correct position found: go to next z_final-coordinate
         exit
      else
!go to next interval
         delp = delp + pdf_zfine_c(j)*(zfine_c(j+1)-zfine_c(j))
      endif
   enddo
enddo
!
!---------------distribute nncz_final grid points outside core----------
!-----------------------(following the non-core pdf)--------------------
!
!probability that has to be matched for each interval
pc = one/float(nncz_final-1)
zfinal_nc(1)=zfine_nc(1)
!
pc_i = zero
do i=2, nncz_final
!probability that has to be matched
   pc_i = pc_i + pc
!
   delp = zero
   do j=1, nzfine_nc-1
      zi =  (pc_i-delp)/pdf_zfine_nc(j)  + zfine_nc(j)
      zfinal_nc(i) = zi
      if(zi.le. zfine_nc(j+1)) then
!;correct position found: go to next z_final-coordinate
         exit
      else
!go to next interval
         delp = delp + pdf_zfine_nc(j)*(zfine_nc(j+1)-zfine_nc(j))
      endif
   enddo
enddo
!
!----------------------combined distribution----------------------------
!
k=1
do i=1, ncz_final
   zfinal(k) = zfinal_c(i)
   k=k+1
enddo
zfinal_c(ncz_final)=rmin
do i=1, nncz_final-1
   zfinal(k) = zfinal_nc(i+1)
   k=k+1
enddo
zfinal(nz_final)=rmax
!
!calculate final pdf (only for checks)
allocate(zfinal_mid(nz_final-1))
allocate(pdf_zfinal(nz_final-1))
allocate(p_zfinal(nz_final-1))
!
call calc_pdf(nz_final, zfinal, zfinal_mid, pdf_zfinal, p_zfinal)
!
!--------------------------print out everything-------------------------
!
open(1, file=output_dir//'/pdf_zgridf_c.dat')
   write(1,*) 'optimum core-distribution of z-coordinates (w.r.t. radial grid)'
   write(1,'(3a20)') 'z_mid', 'pdf(z)', 'p(z)'
   do i=1, nzfine_c-1 
      write(1,'(3es20.8)') zfine_mid_c(i), pdf_zfine_c(i), p_zfine_c(i)
   enddo
close(1)
!
open(1, file=output_dir//'/pdf_zgridf_nc.dat')
   write(1,*) 'optimum non-core-distribution of z-coordinates (w.r.t. radial grid)'
   write(1,'(3a20)') 'z_mid', 'pdf(z)', 'p(z)'
   do i=1, nzfine_nc-1 
      write(1,'(3es20.8)') zfine_mid_nc(i), pdf_zfine_nc(i), p_zfine_nc(i)
   enddo
close(1)
!
open(1, file=output_dir//'/pdf_zgridf.dat')
   write(1,*) 'optimum (complete) distribution of z-coordinates w.r.t. radial grid'
   write(1,'(3a20)') 'z_mid', 'pdf(z)', 'p(z)'
   do i=1, nzfine-1 
      write(1,'(3es20.8)') zfine_mid(i), pdf_zfine(i), p_zfine(i)
   enddo
close(1)
!
open(1, file=output_dir//'/pdf_zgrid1.dat')
   write(1,*) 'actually used distribution at first step: gridxyz_optb'
   write(1,'(3a20)') 'z_mid', 'pdf(z)', 'p(z)'
   do i=1, nz_final-1 
      write(1,'(3es20.8)') zfinal_mid(i), pdf_zfinal(i), p_zfinal(i)
   enddo
close(1)
!
!-------------------create complete z grid------------------------------
!
allocate(z(ndzmax), stat=err)
   if(err.ne.0) stop 'error in gridxyz_opt: allocation'
!
z(1)=-zfinal(nz_final)
do i=2, nz_final
   z(i)=-zfinal(nz_final+2-i)
enddo
!
do i=1, nz_final
   z(nz_final + i) = zfinal(i)
enddo
!
z(ndzmax)=zfinal(nz_final)
!
!-------------------------------create y grid---------------------------
!
!!check if dimensions of x and z grid are compatible
!if(ndzmax.ne.ndxmax.or.ndz.ne.ndx) stop 'dimensions of z and x not compatible'
!!set x equal to z-grid
!z=x
!
allocate(y(ndymax), stat=err)
   if(err.ne.0) stop 'error in gridxyz_opt: allocation'
!check if dimensions of y and x grid are compatible
if(ndymax.ne.ndxmax.or.ndy.ne.ndx) stop 'dimensions of y and x not compatible'
!set y equal to x-grid
y=x
!
!deallocate 1d arrays, since not needed anymore
deallocate(r1d_dum)
!
!write(*,*) x
!write(*,*) 
!write(*,*) y
!write(*,*) 
!write(*,*) z
!write(*,*)
!write(*,*) x-z
!stop 'go on in gridxyz_opt1d'
!
!
end subroutine gridxyz_opt1d
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine gridxyz_opt2d
!
!
use prog_type
use fund_const
use dime3d, only: ndx, ndy, ndz, ndxmax, ndymax, ndzmax, x, y, z
use params_input, only: rmax
use params_stellar, only: sr
use modext, only: nr_modext, ntheta_modext, r_modext, theta_modext
use mod_directories, only: model2d_file, model_dir
use hdf5
use mod_sort, only: bubblesort, quicksort, uniq_arr, uniq_elements
use mod_grid, only: recalc_grid1d
!
implicit none
!
! ... local scalars
integer(i4b) :: i, j, k, lx, ly, lz
integer(i4b) :: err, nuniq_x, nuniq_y, nuniq_z, nav, nrest_uniq, nrest_grid, indx, nphi, ntheta
integer(i4b) :: ndx_dum, ndy_dum, ndz_dum
integer(i4b) :: indx_dum, nr_dum
integer(i4b), parameter :: nphi_max=45
real(dp), parameter :: deltau_max=one/3.d0
real(dp) :: val1, val2, xav
real(dp) :: delx, delx1, delx2, max_delx, maxx, dum
real(dp) :: dely, dely1, dely2, max_dely, maxy
real(dp) :: delz, delz1, delz2, max_delz, maxz
real(dp) :: ts, te
!
! ... local errors
real(dp), dimension(:), allocatable :: xcoords, ycoords, zcoords, xcoords_dum, ycoords_dum, zcoords_dum
real(dp), dimension(:), allocatable :: phi, theta
real(dp), dimension(:), allocatable :: xdum, ydum, zdum
real(dp), dimension(:), allocatable :: radius_dum
!
! ... for hdf5-file
integer(hid_t) :: file_id, group_id, attr_id
integer(hsize_t), dimension(1) :: dims_rad, dims_theta, dims_scalars
!
!----------------read in radius and theta grid from 2d file-------------
!
write(*,*) '-------------------read spatial grid from 2d input file------------------------'
write(*,*) 'file name: ', trim(model_dir)//'/'//model2d_file
write(*,*)
!
dims_scalars=(/ 1 /)
!
call h5open_f(err)
call h5fopen_f(trim(model_dir)//'/'//model2d_file, h5f_acc_rdonly_f, file_id, err)
!
!read dimensions
call h5gopen_f(file_id, 'dimensions', group_id, err)
   call h5aopen_f(group_id, 'nr', attr_id, err)
      call h5aread_f(attr_id, h5t_native_integer, nr_modext, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5aopen_f(group_id, 'ntheta', attr_id, err)
      call h5aread_f(attr_id, h5t_native_integer, ntheta_modext, dims_scalars, err)
   call h5aclose_f(attr_id, err)
call h5gclose_f(group_id, err)
!
!read coordinates
dims_rad=(/nr_modext/)
dims_theta=(/ntheta_modext/)
!
allocate(r_modext(nr_modext), stat=err)
   if(err.ne.0) stop 'error in gridxyz_opt2d: allocation'
allocate(theta_modext(ntheta_modext), stat=err)
   if(err.ne.0) stop 'error in gridxyz_opt2d: allocation'
!
call h5gopen_f(file_id, 'coordinates', group_id, err)
   call h5dopen_f(group_id, 'r', attr_id, err)
      call h5dread_f(attr_id, h5t_native_double, r_modext, dims_rad, err)
   call h5dclose_f(attr_id, err)
   call h5dopen_f(group_id, 'theta', attr_id, err)
      call h5dread_f(attr_id, h5t_native_double, theta_modext, dims_theta, err)
   call h5dclose_f(attr_id, err)
call h5gclose_f(group_id, err)
!
call h5fclose_f(file_id, err)
call h5close_f(err)
!
!-----------take only those grid points which are less than rmax--------
!
indx_dum=0
do i=1, nr_modext
   if(r_modext(i)/sr.le.rmax) then
      indx_dum=indx_dum+1
   endif
enddo
!
nr_dum=indx_dum
allocate(radius_dum(nr_dum), stat=err)
   if(err.ne.0) stop 'error in gridxyz_opt2d: radius_dum'
!
do i=1, nr_dum
   radius_dum(i) = r_modext(i)/sr
enddo
!
!-------------create dummy azimuthal (phi) grid for x-y-plane-----------
!phi is in interval [0,pi/2]
!theta is in interval [0,pi/2]
!   => first octant
!
nphi=45
!
allocate(phi(nphi), stat=err)
   if(err.ne.0) stop 'error in gridxyz_opt2d: allocation'
!
do i=1, nphi
   phi(i) = float(i-1) * pi / float(nphi-1) / two
enddo
!
!********************new (3d) version***********************************
!
!allocate coordinate vectors for all possible combinations
allocate(xcoords_dum(nr_dum * nphi * ntheta_modext), stat=err)
   if(err.ne.0) stop 'error in gridxyz_opt2d: allocation'
allocate(ycoords_dum(nr_dum * nphi * ntheta_modext), stat=err)
   if(err.ne.0) stop 'error in gridxyz_opt2d: allocation'
allocate(zcoords_dum(nr_dum * nphi * ntheta_modext), stat=err)
   if(err.ne.0) stop 'error in gridxyz_opt2d: allocation'
!
lx=1
!
do i=1, ntheta_modext
   do j=1, nphi
      do k=1, nr_dum
            xcoords_dum(lx)=radius_dum(k)*cos(phi(j))*sin(theta_modext(i))
            ycoords_dum(lx)=radius_dum(k)*sin(phi(j))*sin(theta_modext(i))
            zcoords_dum(lx)=radius_dum(k)*cos(theta_modext(i))
            lx=lx+1
      enddo
   enddo
enddo
!
!round coordinates to 8th digit
xcoords_dum=nint(xcoords_dum*1.d8)/1.d8
ycoords_dum=nint(ycoords_dum*1.d8)/1.d8
zcoords_dum=nint(zcoords_dum*1.d8)/1.d8
!
!sort the arrays
call cpu_time(ts)
call quicksort(xcoords_dum, nr_dum*ntheta_modext*nphi, 1, nr_dum*ntheta_modext*nphi)
call quicksort(ycoords_dum, nr_dum*ntheta_modext*nphi, 1, nr_dum*ntheta_modext*nphi)
call quicksort(zcoords_dum, nr_dum*ntheta_modext*nphi, 1, nr_dum*ntheta_modext*nphi)
call cpu_time(te)
!write(*,*) 'quicksort done', te-ts
!
!find number of unique elements
call cpu_time(ts)
call uniq_arr(xcoords_dum, nr_dum*nphi*ntheta_modext, nuniq_x)
call uniq_arr(ycoords_dum, nr_dum*nphi*ntheta_modext, nuniq_y)
call uniq_arr(zcoords_dum, nr_dum*nphi*ntheta_modext, nuniq_z)
call cpu_time(te)
!write(*,*) 'number of unique elements found', te-ts
!
!store unique elements
allocate(xcoords(nuniq_x), stat=err)
   if(err.ne.0) stop 'error in gridxyz_opt2d: allocation'
allocate(ycoords(nuniq_y), stat=err)
   if(err.ne.0) stop 'error in gridxyz_opt2d: allocation'
allocate(zcoords(nuniq_z), stat=err)
   if(err.ne.0) stop 'error in gridxyz_opt2d: allocation'
!
call cpu_time(ts)
call uniq_elements(nr_dum*ntheta_modext*nphi, nuniq_x, xcoords_dum, xcoords)
call uniq_elements(nr_dum*ntheta_modext*nphi, nuniq_y, ycoords_dum, ycoords)
call uniq_elements(nr_dum*ntheta_modext*nphi, nuniq_z, zcoords_dum, zcoords)
call cpu_time(te)
!write(*,*) 'unique elements done', te-ts
!
!deallocate dummy arrays, which are not needed anymore
deallocate(xcoords_dum, stat=err)
   if(err.ne.0) stop 'error in gridxyz_opt2d: deallocation'
deallocate(ycoords_dum, stat=err)
   if(err.ne.0) stop 'error in gridxyz_opt2d: deallocation'
deallocate(zcoords_dum, stat=err)
   if(err.ne.0) stop 'error in gridxyz_opt2d: deallocation'
!
!
!recalculate grids from large arrays
call cpu_time(ts)
ndx_dum=ndx-1
allocate(xdum(ndx_dum), stat=err)
   if(err.ne.0) stop 'allocation error gridxyz_kee: xdum'
call recalc_grid1d(xcoords, nuniq_x, ndx_dum, 3, xdum)
xdum(ndx_dum-2)=one
xdum(ndx_dum-1)=zero
xdum(ndx_dum)=rmax
xcoords(nuniq_x)=rmax
!
ndy_dum=ndy-1
allocate(ydum(ndy_dum), stat=err)
   if(err.ne.0) stop 'allocation error gridxyz_kee: ydum'
call recalc_grid1d(ycoords, nuniq_y, ndy_dum, 3, ydum)
ydum(ndy_dum-2)=one
ydum(ndy_dum-1)=zero
ydum(ndy_dum)=rmax
ycoords(nuniq_y)=rmax
!
ndz_dum=ndz-1
allocate(zdum(ndz_dum), stat=err)
   if(err.ne.0) stop 'allocation error gridxyz_kee: zdum'
call recalc_grid1d(zcoords, nuniq_z, ndz_dum, 3, zdum)
zdum(ndz_dum-2)=one
zdum(ndz_dum-1)=zero
zdum(ndz_dum)=rmax
zcoords(nuniq_z)=rmax
!
call cpu_time(te)
!write(*,*) 'recalculating grids done', te-ts
!
!sort the final arrays
call bubblesort(xdum, ndx_dum)
call bubblesort(ydum, ndy_dum)
call bubblesort(zdum, ndz_dum)
!
!--------------------check if grid-values occur twice-------------------
!
!for x-axis
!
do i=2, ndx_dum-1
   val1=xdum(i-1)
   val2=xdum(i)
   if(abs(val2-val1).lt.1.d-14) then
!insert value on half range (to upper grid point)
      xdum(i)=(xdum(i+1)+xdum(i))/two
   endif
enddo
!
!same on outer boundary
val1=xdum(ndx_dum)
val2=xdum(ndx_dum-1)
if(abs(val2-val1).lt.1.d-14) then
!insert value on half range (to lower grid point)
   xdum(ndx_dum-1)=(xdum(ndx_dum-2)+xdum(ndx_dum-1))/two
endif
!
!
!for y-axis
do i=2, ndy_dum-1
   val1=ydum(i-1)
   val2=ydum(i)
   if(abs(val2-val1).lt.1.d-14) then
!insert value on half range (to upper grid point)
      ydum(i)=(ydum(i+1)+ydum(i))/two
   endif
enddo
!
!same on outer boundary
val1=ydum(ndy_dum)
val2=ydum(ndy_dum-1)
if(abs(val2-val1).lt.1.d-14) then
!insert value on half range (to lower grid point)
   ydum(ndy_dum-1)=(ydum(ndy_dum-2)+ydum(ndy_dum-1))/two
endif
!
!
!for z-axis
do i=2, ndz_dum-1
   val1=zdum(i-1)
   val2=zdum(i)
   if(abs(val2-val1).lt.1.d-14) then
!insert value on half range (to upper grid point)
      zdum(i)=(zdum(i+1)+zdum(i))/two
   endif
enddo
!
!same on outer boundary
val1=zdum(ndz_dum)
val2=zdum(ndz_dum-1)
if(abs(val2-val1).lt.1.d-14) then
!insert value on half range (to lower grid point)
   zdum(ndz_dum-1)=(zdum(ndz_dum-2)+zdum(ndz_dum-1))/two
endif
!
!-----make outer part of grid equidistant in order that delta(xdum)-----
!---------------------is not decreasing---------------------------------
!note: maybe logarithmic increment is much better
!
!for x-axis
!
do i=2, ndx_dum
!apply algorithm only if xdum gt 1
   if(xdum(i).gt.one) then
      delx1=xdum(i)-xdum(i-1)
      max_delx=delx1*(ndx_dum-i)
      maxx=xdum(i)+max_delx
      if(maxx.lt.xcoords(nuniq_x)) then
         indx=i
      endif
   endif
enddo
!
if(ndx_dum.ne.indx) then
   delx=(xcoords(nuniq_x)-xdum(indx))/(ndx_dum-indx)
   do i=indx+1, ndx_dum
      xdum(i)=xdum(i-1)+delx
   enddo
endif
!
!
!for y-axis
!
do i=2, ndy_dum
!apply algorithm only if xdum gt 1
   if(ydum(i).gt.one) then
      dely1=ydum(i)-ydum(i-1)
      max_dely=dely1*(ndy_dum-i)
      maxy=ydum(i)+max_dely
      if(maxy.lt.ycoords(nuniq_y)) then
         indx=i
      endif
   endif
enddo
!
if(ndy_dum.ne.indx) then
   dely=(ycoords(nuniq_x)-ydum(indx))/(ndy_dum-indx)
   do i=indx+1, ndy_dum
      ydum(i)=ydum(i-1)+dely
   enddo
endif
!
!
!for z-axis
!
do i=2, ndz_dum
!apply algorithm only if xdum gt 1
   if(zdum(i).gt.one) then
      delz1=zdum(i)-zdum(i-1)
      max_delz=delz1*(ndz_dum-i)
      maxz=zdum(i)+max_delz
      if(maxz.lt.zcoords(nuniq_z)) then
         indx=i
      endif
   endif
enddo
!
if(ndz_dum.ne.indx) then
   delz=(zcoords(nuniq_z)-zdum(indx))/(ndz_dum-indx)
   do i=indx+1, ndz_dum
      zdum(i)=zdum(i-1)+delz
   enddo
endif
!
!-------------------create complete grid--------------------------------
!
allocate(x(ndxmax), stat=err)
   if(err.ne.0) stop 'error in gridxyz_opt2d: allocation'
allocate(y(ndymax), stat=err)
   if(err.ne.0) stop 'error in gridxyz_opt2d: allocation'
allocate(z(ndzmax), stat=err)
   if(err.ne.0) stop 'error in gridxyz_opt2d: allocation'
!
!x-grid
x(1)=-xdum(ndx_dum)
do i=2, ndx_dum
   x(i)=-xdum(ndx_dum+2-i)
enddo
do i=1, ndx_dum
   x(ndx_dum + i) = xdum(i)
enddo
x(ndxmax)=xdum(ndx_dum)
!
!y-grid
y(1)=-ydum(ndy_dum)
do i=2, ndy_dum
   y(i)=-ydum(ndy_dum+2-i)
enddo
do i=1, ndy_dum
   y(ndy_dum + i) = ydum(i)
enddo
y(ndymax)=ydum(ndy_dum)
!
!z-grid
z(1)=-zdum(ndz_dum)
do i=2, ndz_dum
   z(i)=-zdum(ndz_dum+2-i)
enddo
do i=1, ndz_dum
   z(ndz_dum + i) = zdum(i)
enddo
z(ndzmax)=zdum(ndz_dum)
!
do i=1, ndzmax
   write(*,*) z(i)
enddo
stop
!
!deallocate r and theta array from model, since will be allocated later on
deallocate(r_modext)
deallocate(theta_modext)
!
end subroutine gridxyz_opt2d
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine gridxyz_opt2db
!
!reads in radial grid from input model
!reads in azimuthal angle grid  from input model
!
!calculates the probability density of the radial grid,
!and corresponding probability density of the x,y,z-coordinates
!
!calculates the actualy x,y,z-coordinates following this distrubution
!   (with input nc_final, nnc_final: number of points inside core and outside core)
!
!
use prog_type
use fund_const
use options, only: input_mod_dim
use mod_directories, only: output_dir
use dime1d, only: n1d_dum, r1d_dum
use dime3d, only: ncx, ndx, ncz, ndz, ndy, ndxmax, ndymax, ndzmax, x, y, z
use params_stellar, only: sr
use modext, only: nr_modext, r_modext, ntheta_modext, theta_modext
use params_input, only: rmax
use mod_interp1d, only: find_index
use mod_grid, only: calc_pdf, grid_log
!
implicit none
!
! ... local scalars
integer(i4b) :: i, j, k, err
integer(i4b) :: iim2, iim1, ii, iip1
integer(i4b) :: nphi, nx, nx_c, nx_nc, nz, nz_c, nz_nc, &
                nxfine_c, nxfine_nc, nxfine, ndx_dum, nncx_final, nx_final, ncx_final, &
                nzfine_c, nzfine_nc, nzfine, ndz_dum, nncz_final, nz_final, ncz_final
real(dp) :: integral, cc, grad, xinterc, zinterc, del, rmin, xmin, xmax, zmin, zmax, &
            x_jm1, x_j, z_jm1, z_j, r_jm1, r_j, phi_jm1, phi_j, pdfr_jm1, pdfr_j, &
            pdfp_jm1, pdfp_j, h_jm1, h_j, prob_nc, prob_c, pc, pc_i, xi, zi, delp, &
            delx1, maxx, max_delx, fdum
!
! ... local arrays
!real(dp), dimension(:), allocatable :: xcoords, ycoords, zcoords, xcoords_dum, ycoords_dum, zcoords_dum
real(dp), dimension(:), allocatable :: phi, phidum, phi_mid, pdf_phi, p_phi, &
                                       r1d_mid, pdf_r1d, p_r1d, &
                                       xfine_c, xfine_mid_c, pdf_xfine_c, p_xfine_c, &
                                       zfine_c, zfine_mid_c, pdf_zfine_c, p_zfine_c, &
                                       xfine_nc, xfine_mid_nc, pdf_xfine_nc, p_xfine_nc, &
                                       zfine_nc, zfine_mid_nc, pdf_zfine_nc, p_zfine_nc, &
                                       xfine, xfine_mid, pdf_xfine, p_xfine, &
                                       zfine, zfine_mid, pdf_zfine, p_zfine
real(dp), dimension(:), allocatable :: xdum, xfinal, xfinal_c, xfinal_nc, xfinal_mid, pdf_xfinal, p_xfinal, &
                                       zdum, zfinal, zfinal_c, zfinal_nc, zfinal_mid, pdf_zfinal, p_zfinal
!
!----------------------set up radial grid-------------------------------
!
if(input_mod_dim.ne.2) then
   stop 'error in gridxyz_opt2db: input_mod_dim ne 2'
endif
!
call read_mod2d
!
n1d_dum=nr_modext
allocate(r1d_dum(n1d_dum), stat=err)
   if(err.ne.0) stop 'error in gridxyz_opt2d: allocation'
r1d_dum=r_modext/sr
!
rmin=r1d_dum(1)
!rmax=r1d_dum(n1d_dum)
!
allocate(r1d_mid(n1d_dum-1), stat=err)
   if(err.ne.0) stop 'error in gridxyz_opt2d: allocation'
allocate(pdf_r1d(n1d_dum-1), stat=err)
   if(err.ne.0) stop 'error in gridxyz_opt2d: allocation'
allocate(p_r1d(n1d_dum-1), stat=err)
   if(err.ne.0) stop 'error in gridxyz_opt2d: allocation'
!
!calculate probability density
call calc_pdf(n1d_dum, r1d_dum, r1d_mid, pdf_r1d, p_r1d)
!
!---------------create dummy phi grid (measured from x to z-axis)-------
!phi is in interval [0,pi/2]
!
nphi=45
fdum=ntheta_modext/two
nphi=floor(fdum)+1
!
allocate(phi(nphi), stat=err)
   if(err.ne.0) stop 'error in gridxyz_opt2d: allocation'
allocate(phidum(nphi), stat=err)
   if(err.ne.0) stop 'error in gridxyz_opt2d: allocation'
allocate(phi_mid(nphi-1), stat=err)
   if(err.ne.0) stop 'error in gridxyz_opt2d: allocation'
allocate(pdf_phi(nphi-1), stat=err)
   if(err.ne.0) stop 'error in gridxyz_opt2d: allocation'
allocate(p_phi(nphi-1), stat=err)
   if(err.ne.0) stop 'error in gridxyz_opt2d: allocation'
!
!phi grid equidistant
phi(1)=zero
do i=2, nphi-1
!   phi(i) = float(i-1) * pi / float(nphi-1) / two
   phi(i) = theta_modext(i)
!   write(*,*) i, phi(i), theta_modext(i)*180./pi
enddo
phi(nphi)=pi/two
!
!
if(abs(phi(1)).gt.small_number) stop 'error in gridxyz_opt2d: angular grid needs to start at 0 (needs to be checked)'
if(abs(phi(nphi)-pi/two).gt.small_number) stop 'error in gridxyz_opt2d: angular grid needs to end at pi/2 (needs to be checked)'
do i=2, nphi
   if(phi(i).lt.phi(i-1)) stop 'error in gridxyz_opt2d: angular grid needs to be monotonic'
enddo
!
!phi grid with high resolution at small angles
!call grid_log(1.d-3,pi/two,nphi,phi)
!!
!!phi grid with high resolution at large angles
!phidum(nphi)=pi/two
!k=2
!do i=nphi-1, 1, -1
!   phidum(i) = phidum(i+1)-(phi(k)-phi(k-1))
!   k=k+1
!enddo
!phi=phidum
!
!calculate probability density
call calc_pdf(nphi, phi, phi_mid, pdf_phi, p_phi)
!
!-------------calculate pdf for x-coordinates on a fine grid------------
!---------------1. pdf for core points (normalized to core)-------------
!---------------2. pdf for non-core points (normlized to non-core)------
!---------------3. complete pdf (normalized to complete range)----------
!
!----------------------------core pdf-----------------------------------
!
nxfine_c=100
allocate(xfine_c(nxfine_c), stat=err)
allocate(xfine_mid_c(nxfine_c-1), stat=err)
allocate(pdf_xfine_c(nxfine_c-1), stat=err)
allocate(p_xfine_c(nxfine_c-1), stat=err)
!
!core grid confined towards rmin
cc=half
grad = rmin/log10((cc+one)/cc)
zinterc = -rmin*log10(cc)/log10((cc+one)/cc)
do i=1, nxfine_c
   xfine_c(i) = grad * log10(float(i)/(nxfine_c)+cc) + zinterc
enddo
xfine_c(1)=zero
xfine_c(nxfine_c)=rmin
!
!
!create mid points
do i=1, nxfine_c-1
   xfine_mid_c(i) = (xfine_c(i+1)+xfine_c(i))/two
enddo
!
!
!calculate pdf for core coordinates
nz_c=200
nz_nc=400
nz=nz_c+nz_nc
allocate(zdum(nz), stat=err)
!
do i=1, nxfine_c-1
!
!create z-coordinates for z-integration (obtain marignal distribution by integrating over z)
   zmin = sqrt(rmin**2-xfine_mid_c(i)**2)
   zmax = sqrt(rmax**2-xfine_mid_c(i)**2)
!
   cc=half
   grad = (rmin-zmin)/log10((cc+one)/cc)
   zinterc = (zmin*log10(cc+one)-rmin*log10(cc))/log10((cc+one)/cc)
   k=1
   do j=1, nz_c
      zdum(k) = grad * log10(float(j-1)/(nz_c-1)+cc) + zinterc
      k=k+1
   enddo
   del=log10(zmax/rmin)/float(nz_nc)
   do j=1, nz_nc
      zdum(k) = zdum(k-1)*ten**del
      k=k+1
   enddo
!
!calculate pdf for the given x
   integral=zero
   z_jm1 = zdum(1)
!
!calculate pdf for r and phi distribtion of the given (x,z)
   r_jm1 = sqrt(xfine_mid_c(i)**2+zdum(1)**2)
   phi_jm1 = atan(zdum(1)/xfine_mid_c(i))
   call find_index(r_jm1, r1d_dum, n1d_dum, iim2, iim1, ii, iip1)
   pdfr_jm1 = pdf_r1d(iim1)
   call find_index(phi_jm1, phi, nphi, iim2, iim1, ii, iip1)
   pdfp_jm1 = pdf_phi(iim1)
!   pdfp_jm1 = two/pi
!calculate joint pdf
   h_jm1 = pdfr_jm1*pdfp_jm1/sqrt(xfine_mid_c(i)**2+z_jm1**2)
!
   do j=2, nz
      z_j = zdum(j)
!calculate pdf for r and phi distribtion of the given (x,z)
      r_j = sqrt(xfine_mid_c(i)**2+zdum(j)**2)
      phi_j = atan(zdum(j)/xfine_mid_c(i))
      call find_index(r_j, r1d_dum, n1d_dum, iim2, iim1, ii, iip1)
      pdfr_j = pdf_r1d(iim1)
      call find_index(phi_j, phi, nphi, iim2, iim1, ii, iip1)
      pdfp_j = pdf_phi(iim1)
!      pdfp_j = two/!pi
!calculate joint pdf
      h_j = pdfr_j*pdfp_j/sqrt(xfine_mid_c(i)**2+z_j**2)
!
!perform integration over z
      integral = integral + half*(h_jm1+h_j)*(z_j-z_jm1)
!
      z_jm1 = z_j
      h_jm1 = h_j
!
   enddo
!
   pdf_xfine_c(i) = integral
!
enddo
!
!
!
!calculate probability of finding x in corresponding bins
do i=1, nxfine_c-1
   p_xfine_c(i) = pdf_xfine_c(i)*(xfine_c(i+1)-xfine_c(i))
!   write(*,*) xfine_mid_c(i), pdf_xfine_c(i)*four, p_xfine_c(i)
enddo
!
prob_c = sum(p_xfine_c)
!
!normalize pdf to core interval
pdf_xfine_c=pdf_xfine_c/prob_c
!
!--------------------------non core pdf---------------------------------
!
nxfine_nc=200
allocate(xfine_nc(nxfine_nc), stat=err)
allocate(xfine_mid_nc(nxfine_nc-1), stat=err)
allocate(pdf_xfine_nc(nxfine_nc-1), stat=err)
allocate(p_xfine_nc(nxfine_nc-1), stat=err)
!
!non-core grid logarithmic
call grid_log(rmin, rmax, nxfine_nc, xfine_nc)
!
!create mid points
do i=1, nxfine_nc-1
   xfine_mid_nc(i) = (xfine_nc(i+1)+xfine_nc(i))/two
enddo
!
!
!
deallocate(zdum)
nz_c=100
nz_nc=100
nz=nz_c+nz_nc
allocate(zdum(nz), stat=err)
!
!calculate pdf for non-core x-coordinates (on half grid)
do i=1, nxfine_nc-1
!
!create z-coordinates for z-integration (obtain marignal distribution by integrating over z)
   zmin = zero
   zmax = sqrt(rmax**2-xfine_mid_nc(i)**2)
!
   cc=half
   grad = (rmin-zmin)/log10((cc+one)/cc)
   zinterc = (zmin*log10(cc+one)-rmin*log10(cc))/log10((cc+one)/cc)
   zdum(1) = zmin
   k=2
   do j=2, nz_c
      zdum(k) = grad * log10(float(j-1)/(nz_c-1)+cc) + zinterc
      k=k+1
   enddo
   del=log10(zmax/rmin)/float(nz_nc)
   do j=1, nz_nc
      zdum(k) = zdum(k-1)*ten**del
      k=k+1
   enddo
!
!
!calculate pdf for the given x
   integral=zero
   z_jm1 = zdum(1)
!
!calculate pdf for r and phi distribtion of the given (x,z)
   r_jm1 = sqrt(xfine_mid_nc(i)**2+zdum(1)**2)
   phi_jm1 = atan(zdum(1)/xfine_mid_nc(i))
   call find_index(r_jm1, r1d_dum, n1d_dum, iim2, iim1, ii, iip1)
   pdfr_jm1 = pdf_r1d(iim1)
   call find_index(phi_jm1, phi, nphi, iim2, iim1, ii, iip1)
   pdfp_jm1 = pdf_phi(iim1)
!   pdfp_jm1 = two/pi
!calculate joint pdf
   h_jm1 = pdfr_jm1*pdfp_jm1/sqrt(xfine_mid_nc(i)**2+z_jm1**2)
!
   do j=2, nz
      z_j = zdum(j)
!calculate pdf for r and phi distribtion of the given (x,z)
      r_j = sqrt(xfine_mid_nc(i)**2+zdum(j)**2)
      phi_j = atan(zdum(j)/xfine_mid_nc(i))
      call find_index(r_j, r1d_dum, n1d_dum, iim2, iim1, ii, iip1)
      pdfr_j = pdf_r1d(iim1)
      call find_index(phi_j, phi, nphi, iim2, iim1, ii, iip1)
      pdfp_j = pdf_phi(iim1)
!      pdfp_j = two/!pi
!calculate joint pdf
      h_j = pdfr_j*pdfp_j/sqrt(xfine_mid_nc(i)**2+z_j**2)
!
!perform integration over z
      integral = integral + half*(h_jm1+h_j)*(z_j-z_jm1)
!
      z_jm1 = z_j
      h_jm1 = h_j
!
   enddo
!
   pdf_xfine_nc(i) = integral
!
enddo
!
!
!
!calculate probability of finding x in corresponding bins
do i=1, nxfine_nc-1
   p_xfine_nc(i) = pdf_xfine_nc(i)*(xfine_nc(i+1)-xfine_nc(i))
enddo
!
!
prob_nc=sum(p_xfine_nc)
!
!normalize pdf to non-core interval
pdf_xfine_nc=pdf_xfine_nc/prob_nc
!
!------------complete pdf normalized to complete range [0.,rmax]-------------
!-----------------------(only required for checks)---------------------------
!
nxfine = nxfine_c + nxfine_nc -1 
allocate(xfine(nxfine), stat=err)
allocate(xfine_mid(nxfine-1), stat=err)
allocate(pdf_xfine(nxfine-1), stat=err)
allocate(p_xfine(nxfine-1), stat=err)
!
!
k=1
do i=1, nxfine_c
   xfine(k) = xfine_c(i)
   k=k+1
enddo
do i=2, nxfine_nc
   xfine(k) = xfine_nc(i)
   k=k+1
enddo
!
k=1
do i=1, nxfine_c-1
   xfine_mid(k) = xfine_mid_c(i)
   pdf_xfine(k) = pdf_xfine_c(i)*prob_c
   p_xfine(k) = p_xfine_c(i)
   k=k+1
enddo
do i=1, nxfine_nc-1
   xfine_mid(k) = xfine_mid_nc(i)
   pdf_xfine(k) = pdf_xfine_nc(i)*prob_nc
   p_xfine(k) = p_xfine_nc(i)
   k=k+1
enddo
!
pdf_xfine=pdf_xfine/sum(p_xfine)
!
!---------------------create final grid from zero to rmax---------------
!
!ndx_dum=ndx-1
!nx_final=ndx_dum
!nncx_final=ndx_dum-ncx_final+1   !+1 because rmin point occurrs twice
!write(*,*) nx_final, ncx_final, nncx_final
!
ndx_dum=ndx-1
nx_final=ndx_dum
ncx_final = ncx
nncx_final=ndx_dum-ncx_final+1   !+1 because rmin point occurrs twice
!
write(*,*) 'probability of finding x-coordinate inside core ', prob_c
write(*,*) 'probability of finding x-coordinate outside core', prob_nc
write(*,*) ' => optimum number of core points to match distribution    ', nint(ndx_dum*prob_c)
write(*,*) ' => optimum number of non-core points to match distribution', ndx_dum-nint(ndx_dum*prob_c)+1
write(*,*)
write(*,*) 'actually used number of core points    ', ncx_final
write(*,*) 'actually used number of non-core points', nncx_final
write(*,*) '# core points / # points    ', float(ncx_final)/float(nncx_final+ncx_final)
write(*,*) '# non-core points / # points', float(nncx_final)/float(nncx_final+ncx_final)
write(*,*)
!
!xdum only up to ndx-1 since phantom point will be calculated later on
allocate(xfinal_c(ncx_final), stat=err)
allocate(xfinal_nc(nncx_final), stat=err)
allocate(xfinal(nx_final), stat=err)
!
!-------------------distribute nxc_final grid points in core------------
!-------------------------(following the core pdf)----------------------
!
!probability that has to be matched for each interval
pc = one/float(ncx_final-1)
xfinal_c(1)=zero
!
pc_i = zero
do i=2, ncx_final
!probability that has to be matched
   pc_i = pc_i + pc
!
   delp = zero
   do j=1, nxfine_c-1
      xi =  (pc_i-delp)/pdf_xfine_c(j)  + xfine_c(j)
      xfinal_c(i) = xi
      if(xi.le.xfine_c(j+1)) then
!correct position found: go to next x_final-coordinate
         exit
      else
!go to next interval
         delp = delp + pdf_xfine_c(j)*(xfine_c(j+1)-xfine_c(j))
      endif
   enddo
enddo
!
!---------------distribute nncx_final grid points outside core----------
!-----------------------(following the non-core pdf)--------------------
!
!probability that has to be matched for each interval
pc = one/float(nncx_final-1)
xfinal_nc(1)=xfine_nc(1)
!
pc_i = zero
do i=2, nncx_final
!probability that has to be matched
   pc_i = pc_i + pc
!
   delp = zero
   do j=1, nxfine_nc-1
      xi =  (pc_i-delp)/pdf_xfine_nc(j)  + xfine_nc(j)
      xfinal_nc(i) = xi
      if(xi.le. xfine_nc(j+1)) then
!;correct position found: go to next x_final-coordinate
         exit
      else
!go to next interval
         delp = delp + pdf_xfine_nc(j)*(xfine_nc(j+1)-xfine_nc(j))
      endif
   enddo
enddo
!
!----------------------combined distribution----------------------------
!
k=1
do i=1, ncx_final
   xfinal(k) = xfinal_c(i)
   k=k+1
enddo
xfinal_c(ncx_final)=rmin
do i=1, nncx_final-1
   xfinal(k) = xfinal_nc(i+1)
   k=k+1
enddo
xfinal(nx_final)=rmax
!
!calculate final pdf (only for checks)
allocate(xfinal_mid(nx_final-1))
allocate(pdf_xfinal(nx_final-1))
allocate(p_xfinal(nx_final-1))
!
call calc_pdf(nx_final, xfinal, xfinal_mid, pdf_xfinal, p_xfinal)
!
!--------------------------print out everything-------------------------
!
open(1, file=output_dir//'/pdf_rgridf.dat')
   write(1,*) 'radial grid (input)'
   write(1,'(3a20)') 'r_mid', 'pdf(r)', 'p(r)'
   do i=1, n1d_dum-1 
      write(1,'(3es20.8)') r1d_mid(i), pdf_r1d(i), p_r1d(i)
   enddo
close(1)
!
open(1, file=output_dir//'/pdf_phigridf.dat')
   write(1,*) 'phi grid (input)'
   write(1,'(3a20)') 'phi_mid', 'pdf(phi)', 'p(phi)'
   do i=1, nphi-1 
      write(1,'(3es20.8)') phi_mid(i), pdf_phi(i), p_phi(i)
   enddo
close(1)

!stop 'go on in grid_spatial'
!
open(1, file=output_dir//'/pdf_xgridf_c.dat')
   write(1,*) 'optimum core-distribution of x-coordinates (w.r.t. radial grid)'
   write(1,'(3a20)') 'x_mid', 'pdf(x)', 'p(x)'
   do i=1, nxfine_c-1 
      write(1,'(3es20.8)') xfine_mid_c(i), pdf_xfine_c(i), p_xfine_c(i)
   enddo
close(1)
!
open(1, file=output_dir//'/pdf_xgridf_nc.dat')
   write(1,*) 'optimum non-core-distribution of x-coordinates (w.r.t. radial grid)'
   write(1,'(3a20)') 'x_mid', 'pdf(x)', 'p(x)'
   do i=1, nxfine_nc-1 
      write(1,'(3es20.8)') xfine_mid_nc(i), pdf_xfine_nc(i), p_xfine_nc(i)
   enddo
close(1)
!
open(1, file=output_dir//'/pdf_xgridf.dat')
   write(1,*) 'optimum (complete) distribution of x-coordinates w.r.t. radial grid'
   write(1,'(3a20)') 'x_mid', 'pdf(x)', 'p(x)'
   do i=1, nxfine-1 
      write(1,'(3es20.8)') xfine_mid(i), pdf_xfine(i), p_xfine(i)
   enddo
close(1)
!
open(1, file=output_dir//'/pdf_xgrid1.dat')
   write(1,*) 'actually used distribution at first step: gridxyz_optb'
   write(1,'(3a20)') 'x_mid', 'pdf(x)', 'p(x)'
   do i=1, nx_final-1 
      write(1,'(3es20.8)') xfinal_mid(i), pdf_xfinal(i), p_xfinal(i)
   enddo
close(1)
!
!-------------------create complete grid--------------------------------
!
allocate(x(ndxmax), stat=err)
   if(err.ne.0) stop 'error in gridxyz_opt3d: allocation'
!
x(1)=-xfinal(nx_final)
do i=2, nx_final
   x(i)=-xfinal(nx_final+2-i)
enddo
!
do i=1, nx_final
   x(nx_final + i) = xfinal(i)
enddo
!
x(ndxmax)=xfinal(nx_final)
!
do i=1, ndxmax
   if(abs(x(i)-one).lt.1.d-8) x(i)=one
enddo
!
!-------------calculate pdf for z-coordinates on a fine grid------------
!---------------1. pdf for core points (normalized to core)-------------
!---------------2. pdf for non-core points (normlized to non-core)------
!---------------3. complete pdf (normalized to complete range)----------
!
!----------------------------core pdf-----------------------------------
!
nzfine_c=100
allocate(zfine_c(nzfine_c), stat=err)
allocate(zfine_mid_c(nzfine_c-1), stat=err)
allocate(pdf_zfine_c(nzfine_c-1), stat=err)
allocate(p_zfine_c(nzfine_c-1), stat=err)
!
!core grid confined towards rmin
cc=half
grad = rmin/log10((cc+one)/cc)
xinterc = -rmin*log10(cc)/log10((cc+one)/cc)
do i=1, nzfine_c
   zfine_c(i) = grad * log10(float(i)/(nzfine_c)+cc) + xinterc
enddo
zfine_c(1)=zero
zfine_c(nzfine_c)=rmin
!
!create mid points
do i=1, nzfine_c-1
   zfine_mid_c(i) = (zfine_c(i+1)+zfine_c(i))/two
enddo
!
!
!calculate pdf for core coordinates
nx_c=200
nx_nc=400
nx=nx_c+nx_nc
allocate(xdum(nx), stat=err)
!
do i=1, nzfine_c-1
!
!create x-coordinates for x-integration (obtain marignal distribution by integrating over x)
!   write(*,*) zfine_mid_c(i)
   xmin = sqrt(rmin**2-zfine_mid_c(i)**2)
   xmax = sqrt(rmax**2-zfine_mid_c(i)**2)
!   write(*,*) xmin, xmax
!
   cc=half
   grad = (rmin-xmin)/log10((cc+one)/cc)
   xinterc = (xmin*log10(cc+one)-rmin*log10(cc))/log10((cc+one)/cc)
   k=1
   do j=1, nx_c
      xdum(k) = grad * log10(float(j-1)/(nx_c-1)+cc) + xinterc
      k=k+1
   enddo
   del=log10(xmax/rmin)/float(nx_nc)
   do j=1, nx_nc
      xdum(k) = xdum(k-1)*ten**del
      k=k+1
   enddo
!   write(*,*) xdum
!
!calculate pdf for the given x
   integral=zero
   x_jm1 = xdum(1)
!
!calculate pdf for r and phi distribtion of the given (x,z)
   r_jm1 = sqrt(zfine_mid_c(i)**2+xdum(1)**2)
   phi_jm1 = atan(zfine_mid_c(i)/xdum(1))
!   write(*,*) r_jm1, phi_jm1
   call find_index(r_jm1, r1d_dum, n1d_dum, iim2, iim1, ii, iip1)
   pdfr_jm1 = pdf_r1d(iim1)
   call find_index(phi_jm1, phi, nphi, iim2, iim1, ii, iip1)
   pdfp_jm1 = pdf_phi(iim1)
!   write(*,*) pdfr_jm1, pdfp_jm1
!   pdfp_jm1 = two/pi
!calculate joint pdf
   h_jm1 = pdfr_jm1*pdfp_jm1/sqrt(zfine_mid_c(i)**2+x_jm1**2)
!   write(*,*) h_jm1
!
   do j=2, nx
      x_j = xdum(j)
!calculate pdf for r and phi distribtion of the given (x,z)
      r_j = sqrt(zfine_mid_c(i)**2+xdum(j)**2)
      phi_j = atan(zfine_mid_c(i)/xdum(j))
      call find_index(r_j, r1d_dum, n1d_dum, iim2, iim1, ii, iip1)
      pdfr_j = pdf_r1d(iim1)
      call find_index(phi_j, phi, nphi, iim2, iim1, ii, iip1)
      pdfp_j = pdf_phi(iim1)
!      pdfp_j = two/!pi
!calculate joint pdf
      h_j = pdfr_j*pdfp_j/sqrt(zfine_mid_c(i)**2+x_j**2)
!
!perform integration over z
      integral = integral + half*(h_jm1+h_j)*(x_j-x_jm1)
!
      x_jm1 = x_j
      h_jm1 = h_j 
!
   enddo
!
   pdf_zfine_c(i) = integral
!
enddo
!
!
!
!calculate probability of finding x in corresponding bins
do i=1, nzfine_c-1
   p_zfine_c(i) = pdf_zfine_c(i)*(zfine_c(i+1)-zfine_c(i))
!   write(*,*) xfine_mid_c(i), pdf_xfine_c(i)*four, p_xfine_c(i)
enddo
!
prob_c = sum(p_zfine_c)
!
!normalize pdf to core interval
pdf_zfine_c=pdf_zfine_c/prob_c
!
!--------------------------non core pdf---------------------------------
!
nzfine_nc=200
allocate(zfine_nc(nzfine_nc), stat=err)
allocate(zfine_mid_nc(nzfine_nc-1), stat=err)
allocate(pdf_zfine_nc(nzfine_nc-1), stat=err)
allocate(p_zfine_nc(nzfine_nc-1), stat=err)
!
!non-core grid logarithmic
call grid_log(rmin, rmax, nzfine_nc, zfine_nc)
!
!create mid points
do i=1, nzfine_nc-1
   zfine_mid_nc(i) = (zfine_nc(i+1)+zfine_nc(i))/two
enddo
!
!
!
deallocate(xdum)
nx_c=100
nx_nc=100
nx=nx_c+nx_nc
allocate(xdum(nx), stat=err)
!
!calculate pdf for non-core x-coordinates (on half grid)
do i=1, nzfine_nc-1
!
!create x-coordinates for x-integration (obtain marignal distribution by integrating over x)
   xmin = zero
   xmax = sqrt(rmax**2-zfine_mid_nc(i)**2)
!
   cc=half
   grad = (rmin-xmin)/log10((cc+one)/cc)
   xinterc = (xmin*log10(cc+one)-rmin*log10(cc))/log10((cc+one)/cc)
   xdum(1) = xmin
   k=2
   do j=2, nx_c
      xdum(k) = grad * log10(float(j-1)/(nx_c-1)+cc) + xinterc
      k=k+1
   enddo
   del=log10(xmax/rmin)/float(nx_nc)
   do j=1, nx_nc
      xdum(k) = xdum(k-1)*ten**del
      k=k+1
   enddo
!
!
!
!calculate pdf for the given x
   integral=zero
   x_jm1 = xdum(1)
!
!calculate pdf for r and phi distribtion of the given (x,z)
   r_jm1 = sqrt(zfine_mid_nc(i)**2+xdum(1)**2)
!   phi_jm1 = atan(zfine_mid_nc(i)/xdum(1))
   phi_jm1 = pi/two !avoid division by zero
   call find_index(r_jm1, r1d_dum, n1d_dum, iim2, iim1, ii, iip1)
   pdfr_jm1 = pdf_r1d(iim1)
   call find_index(phi_jm1, phi, nphi, iim2, iim1, ii, iip1)
   pdfp_jm1 = pdf_phi(iim1)
!   pdfp_jm1 = two/pi
!calculate joint pdf
   h_jm1 = pdfr_jm1*pdfp_jm1/sqrt(zfine_mid_nc(i)**2+x_jm1**2)
!
   do j=2, nx
      x_j = xdum(j)
!calculate pdf for r and phi distribtion of the given (x,z)
      r_j = sqrt(zfine_mid_nc(i)**2+xdum(j)**2)
      phi_j = atan(zfine_mid_nc(i)/xdum(j))
      call find_index(r_j, r1d_dum, n1d_dum, iim2, iim1, ii, iip1)
      pdfr_j = pdf_r1d(iim1)
      call find_index(phi_j, phi, nphi, iim2, iim1, ii, iip1)
      pdfp_j = pdf_phi(iim1)
!      pdfp_j = two/!pi
!calculate joint pdf
      h_j = pdfr_j*pdfp_j/sqrt(zfine_mid_nc(i)**2+x_j**2)
!
!perform integration over z
      integral = integral + half*(h_jm1+h_j)*(x_j-x_jm1)
!
      x_jm1 = x_j
      h_jm1 = h_j
!
   enddo
!
   pdf_zfine_nc(i) = integral
!
enddo
!
!
!
!calculate probability of finding x in corresponding bins
do i=1, nzfine_nc-1
   p_zfine_nc(i) = pdf_zfine_nc(i)*(zfine_nc(i+1)-zfine_nc(i))
enddo
!
!
prob_nc=sum(p_zfine_nc)
!
!normalize pdf to non-core interval
pdf_zfine_nc=pdf_zfine_nc/prob_nc
!
!------------complete pdf normalized to complete range [0.,rmax]-------------
!-----------------------(only required for checks)---------------------------
!
nzfine = nzfine_c + nzfine_nc -1 
allocate(zfine(nzfine), stat=err)
allocate(zfine_mid(nzfine-1), stat=err)
allocate(pdf_zfine(nzfine-1), stat=err)
allocate(p_zfine(nzfine-1), stat=err)
!
!
k=1
do i=1, nzfine_c
   zfine(k) = zfine_c(i)
   k=k+1
enddo
do i=2, nzfine_nc
   zfine(k) = zfine_nc(i)
   k=k+1
enddo
!
k=1
do i=1, nzfine_c-1
   zfine_mid(k) = zfine_mid_c(i)
   pdf_zfine(k) = pdf_zfine_c(i)*prob_c
   p_zfine(k) = p_zfine_c(i)
   k=k+1
enddo
do i=1, nzfine_nc-1
   zfine_mid(k) = zfine_mid_nc(i)
   pdf_zfine(k) = pdf_zfine_nc(i)*prob_nc
   p_zfine(k) = p_zfine_nc(i)
   k=k+1
enddo
!
pdf_zfine=pdf_zfine/sum(p_zfine)
!
!---------------------create final grid from zero to rmax---------------
!
ndz_dum=ndz-1
nz_final=ndz_dum
ncz_final = ncz
nncz_final=ndz_dum-ncz_final+1   !+1 because rmin point occurrs twice
!
write(*,*) 'probability of finding z-coordinate inside core ', prob_c
write(*,*) 'probability of finding z-coordinate outside core', prob_nc
write(*,*) ' => optimum number of core points to match distribution    ', nint(ndz_dum*prob_c)
write(*,*) ' => optimum number of non-core points to match distribution', ndz_dum-nint(ndz_dum*prob_c)+1
write(*,*)
write(*,*) 'actually used number of core points    ', ncz_final
write(*,*) 'actually used number of non-core points', nncz_final
write(*,*) '# core points / # points    ', float(ncz_final)/float(nncz_final+ncz_final)
write(*,*) '# non-core points / # points', float(nncz_final)/float(nncz_final+ncz_final)
write(*,*)
!
!zdum only up to ndx-1 since phantom point will be calculated later on
allocate(zfinal_c(ncz_final), stat=err)
allocate(zfinal_nc(nncz_final), stat=err)
allocate(zfinal(nz_final), stat=err)
!
!-------------------distribute nxz_final grid points in core------------
!-------------------------(following the core pdf)----------------------
!
!probability that has to be matched for each interval
pc = one/float(ncz_final-1)
zfinal_c(1)=zero
!
pc_i = zero
do i=2, ncz_final
!probability that has to be matched
   pc_i = pc_i + pc
!
   delp = zero
   do j=1, nzfine_c-1
      zi =  (pc_i-delp)/pdf_zfine_c(j)  + zfine_c(j)
      zfinal_c(i) = zi
      if(zi.le.zfine_c(j+1)) then
!correct position found: go to next z_final-coordinate
         exit
      else
!go to next interval
         delp = delp + pdf_zfine_c(j)*(zfine_c(j+1)-zfine_c(j))
      endif
   enddo
enddo
!
!---------------distribute nncz_final grid points outside core----------
!-----------------------(following the non-core pdf)--------------------
!
!probability that has to be matched for each interval
pc = one/float(nncz_final-1)
zfinal_nc(1)=zfine_nc(1)
!
pc_i = zero
do i=2, nncz_final
!probability that has to be matched
   pc_i = pc_i + pc
!
   delp = zero
   do j=1, nzfine_nc-1
      zi =  (pc_i-delp)/pdf_zfine_nc(j)  + zfine_nc(j)
      zfinal_nc(i) = zi
      if(zi.le. zfine_nc(j+1)) then
!;correct position found: go to next z_final-coordinate
         exit
      else
!go to next interval
         delp = delp + pdf_zfine_nc(j)*(zfine_nc(j+1)-zfine_nc(j))
      endif
   enddo
enddo
!
!----------------------combined distribution----------------------------
!
k=1
do i=1, ncz_final
   zfinal(k) = zfinal_c(i)
   k=k+1
enddo
zfinal_c(ncz_final)=rmin
do i=1, nncz_final-1
   zfinal(k) = zfinal_nc(i+1)
   k=k+1
enddo
zfinal(nz_final)=rmax
!
!calculate final pdf (only for checks)
allocate(zfinal_mid(nz_final-1))
allocate(pdf_zfinal(nz_final-1))
allocate(p_zfinal(nz_final-1))
!
call calc_pdf(nz_final, zfinal, zfinal_mid, pdf_zfinal, p_zfinal)
!
!--------------------------print out everything-------------------------
!
open(1, file=output_dir//'/pdf_zgridf_c.dat')
   write(1,*) 'optimum core-distribution of z-coordinates (w.r.t. radial grid)'
   write(1,'(3a20)') 'z_mid', 'pdf(z)', 'p(z)'
   do i=1, nzfine_c-1 
      write(1,'(3es20.8)') zfine_mid_c(i), pdf_zfine_c(i), p_zfine_c(i)
   enddo
close(1)
!
open(1, file=output_dir//'/pdf_zgridf_nc.dat')
   write(1,*) 'optimum non-core-distribution of z-coordinates (w.r.t. radial grid)'
   write(1,'(3a20)') 'z_mid', 'pdf(z)', 'p(z)'
   do i=1, nzfine_nc-1 
      write(1,'(3es20.8)') zfine_mid_nc(i), pdf_zfine_nc(i), p_zfine_nc(i)
   enddo
close(1)
!
open(1, file=output_dir//'/pdf_zgridf.dat')
   write(1,*) 'optimum (complete) distribution of z-coordinates w.r.t. radial grid'
   write(1,'(3a20)') 'z_mid', 'pdf(z)', 'p(z)'
   do i=1, nzfine-1 
      write(1,'(3es20.8)') zfine_mid(i), pdf_zfine(i), p_zfine(i)
   enddo
close(1)
!
open(1, file=output_dir//'/pdf_zgrid1.dat')
   write(1,*) 'actually used distribution at first step: gridxyz_optb'
   write(1,'(3a20)') 'z_mid', 'pdf(z)', 'p(z)'
   do i=1, nz_final-1 
      write(1,'(3es20.8)') zfinal_mid(i), pdf_zfinal(i), p_zfinal(i)
   enddo
close(1)
!
!-------------------create complete z grid------------------------------
!
allocate(z(ndzmax), stat=err)
   if(err.ne.0) stop 'error in gridxyz_opt3d: allocation'
!
z(1)=-zfinal(nz_final)
do i=2, nz_final
   z(i)=-zfinal(nz_final+2-i)
enddo
!
do i=1, nz_final
   z(nz_final + i) = zfinal(i)
enddo
!
z(ndzmax)=zfinal(nz_final)
!
do i=1, ndzmax
   if(abs(z(i)-one).lt.1.d-8) z(i)=one
enddo
!
!-------------------------------create y grid---------------------------
!
!!check if dimensions of x and z grid are compatible
!if(ndzmax.ne.ndxmax.or.ndz.ne.ndx) stop 'dimensions of z and x not compatible'
!!set x equal to z-grid
!z=x
!
allocate(y(ndymax), stat=err)
   if(err.ne.0) stop 'error in gridxyz_opt3d: allocation'
!check if dimensions of y and x grid are compatible
if(ndymax.ne.ndxmax.or.ndy.ne.ndx) stop 'dimensions of y and x not compatible'
!set y equal to x-grid
y=x
!
!deallocate 1d arrays, since not needed anymore
deallocate(r1d_dum)
!
!write(*,*) x
!write(*,*) 
!write(*,*) y
!write(*,*) 
!write(*,*) z
!stop 'go on in gridxyz_opt3d'
!
!
end subroutine gridxyz_opt2db
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine gridxyz_opt3d
!
!reads in radial grid from input model
!reads in azimuthal angle grid  from input model
!
!calculates the probability density of the radial grid,
!and corresponding probability density of the x,y,z-coordinates
!
!calculates the actualy x,y,z-coordinates following this distrubution
!   (with input nc_final, nnc_final: number of points inside core and outside core)
!
!
use prog_type
use fund_const
use options, only: input_mod_dim
use mod_directories, only: output_dir
use dime1d, only: n1d_dum, r1d_dum
use dime3d, only: ncx, ndx, ncz, ndz, ndy, ndxmax, ndymax, ndzmax, x, y, z
use params_stellar, only: sr
use modext, only: nr_modext, r_modext, ntheta_modext, theta_modext
use params_input, only: rmax
use mod_interp1d, only: find_index
use mod_grid, only: calc_pdf, grid_log
!
implicit none
!
! ... local scalars
integer(i4b) :: i, j, k, err
integer(i4b) :: iim2, iim1, ii, iip1
integer(i4b) :: nphi, nx, nx_c, nx_nc, nz, nz_c, nz_nc, &
                nxfine_c, nxfine_nc, nxfine, ndx_dum, nncx_final, nx_final, ncx_final, &
                nzfine_c, nzfine_nc, nzfine, ndz_dum, nncz_final, nz_final, ncz_final
real(dp) :: integral, cc, grad, xinterc, zinterc, del, rmin, xmin, xmax, zmin, zmax, &
            x_jm1, x_j, z_jm1, z_j, r_jm1, r_j, phi_jm1, phi_j, pdfr_jm1, pdfr_j, &
            pdfp_jm1, pdfp_j, h_jm1, h_j, prob_nc, prob_c, pc, pc_i, xi, zi, delp, &
            delx1, maxx, max_delx, fdum
!
! ... local arrays
!real(dp), dimension(:), allocatable :: xcoords, ycoords, zcoords, xcoords_dum, ycoords_dum, zcoords_dum
real(dp), dimension(:), allocatable :: phi, phidum, phi_mid, pdf_phi, p_phi, &
                                       r1d_mid, pdf_r1d, p_r1d, &
                                       xfine_c, xfine_mid_c, pdf_xfine_c, p_xfine_c, &
                                       zfine_c, zfine_mid_c, pdf_zfine_c, p_zfine_c, &
                                       xfine_nc, xfine_mid_nc, pdf_xfine_nc, p_xfine_nc, &
                                       zfine_nc, zfine_mid_nc, pdf_zfine_nc, p_zfine_nc, &
                                       xfine, xfine_mid, pdf_xfine, p_xfine, &
                                       zfine, zfine_mid, pdf_zfine, p_zfine
real(dp), dimension(:), allocatable :: xdum, xfinal, xfinal_c, xfinal_nc, xfinal_mid, pdf_xfinal, p_xfinal, &
                                       zdum, zfinal, zfinal_c, zfinal_nc, zfinal_mid, pdf_zfinal, p_zfinal
!
!----------------------set up radial grid-------------------------------
!
if(input_mod_dim.ne.3) then
   stop 'error in gridxyz_opt3d: input_mod_dim ne 3'
endif
!
call read_mod3d
!
n1d_dum=nr_modext
allocate(r1d_dum(n1d_dum), stat=err)
   if(err.ne.0) stop 'error in gridxyz_opt3d: allocation'
r1d_dum=r_modext/sr
!
rmin=r1d_dum(1)
!rmax=r1d_dum(n1d_dum)
!
allocate(r1d_mid(n1d_dum-1), stat=err)
   if(err.ne.0) stop 'error in gridxyz_opt3d: allocation'
allocate(pdf_r1d(n1d_dum-1), stat=err)
   if(err.ne.0) stop 'error in gridxyz_opt3d: allocation'
allocate(p_r1d(n1d_dum-1), stat=err)
   if(err.ne.0) stop 'error in gridxyz_opt3d: allocation'
!
!calculate probability density
call calc_pdf(n1d_dum, r1d_dum, r1d_mid, pdf_r1d, p_r1d)
!
!---------------create dummy phi grid (measured from x to z-axis)-------
!phi is in interval [0,pi/2]
!
nphi=45
fdum=ntheta_modext/two
nphi=floor(fdum)+1
!
allocate(phi(nphi), stat=err)
   if(err.ne.0) stop 'error in gridxyz_opt3d: allocation'
allocate(phidum(nphi), stat=err)
   if(err.ne.0) stop 'error in gridxyz_opt3d: allocation'
allocate(phi_mid(nphi-1), stat=err)
   if(err.ne.0) stop 'error in gridxyz_opt3d: allocation'
allocate(pdf_phi(nphi-1), stat=err)
   if(err.ne.0) stop 'error in gridxyz_opt3d: allocation'
allocate(p_phi(nphi-1), stat=err)
   if(err.ne.0) stop 'error in gridxyz_opt3d: allocation'
!
!phi grid equidistant
do i=1, nphi
   phi(i) = float(i-1) * pi / float(nphi-1) / two
!   phi(i) = theta_modext(i)
!   write(*,*) i, phi(i), theta_modext(i)*180./pi
enddo
!
if(abs(phi(1)).gt.1.d-10) stop 'error in gridxyz_opt3d: angular grid needs to start at 0 (needs to be checked)'
if(abs(phi(nphi)-pi/two).gt.1.d-10) stop 'error in gridxyz_opt3d: angular grid needs to end at pi/2 (needs to be checked)'
!
!phi grid with high resolution at small angles
!call grid_log(1.d-3,pi/two,nphi,phi)
!!
!!phi grid with high resolution at large angles
!phidum(nphi)=pi/two
!k=2
!do i=nphi-1, 1, -1
!   phidum(i) = phidum(i+1)-(phi(k)-phi(k-1))
!   k=k+1
!enddo
!phi=phidum
!
!calculate probability density
call calc_pdf(nphi, phi, phi_mid, pdf_phi, p_phi)
!
!-------------calculate pdf for x-coordinates on a fine grid------------
!---------------1. pdf for core points (normalized to core)-------------
!---------------2. pdf for non-core points (normlized to non-core)------
!---------------3. complete pdf (normalized to complete range)----------
!
!----------------------------core pdf-----------------------------------
!
nxfine_c=100
allocate(xfine_c(nxfine_c), stat=err)
allocate(xfine_mid_c(nxfine_c-1), stat=err)
allocate(pdf_xfine_c(nxfine_c-1), stat=err)
allocate(p_xfine_c(nxfine_c-1), stat=err)
!
!core grid confined towards rmin
cc=half
grad = rmin/log10((cc+one)/cc)
zinterc = -rmin*log10(cc)/log10((cc+one)/cc)
do i=1, nxfine_c
   xfine_c(i) = grad * log10(float(i)/(nxfine_c)+cc) + zinterc
enddo
xfine_c(1)=zero
xfine_c(nxfine_c)=rmin
!
!
!create mid points
do i=1, nxfine_c-1
   xfine_mid_c(i) = (xfine_c(i+1)+xfine_c(i))/two
enddo
!
!
!calculate pdf for core coordinates
nz_c=200
nz_nc=400
nz=nz_c+nz_nc
allocate(zdum(nz), stat=err)
!
do i=1, nxfine_c-1
!
!create z-coordinates for z-integration (obtain marignal distribution by integrating over z)
   zmin = sqrt(rmin**2-xfine_mid_c(i)**2)
   zmax = sqrt(rmax**2-xfine_mid_c(i)**2)
!
   cc=half
   grad = (rmin-zmin)/log10((cc+one)/cc)
   zinterc = (zmin*log10(cc+one)-rmin*log10(cc))/log10((cc+one)/cc)
   k=1
   do j=1, nz_c
      zdum(k) = grad * log10(float(j-1)/(nz_c-1)+cc) + zinterc
      k=k+1
   enddo
   del=log10(zmax/rmin)/float(nz_nc)
   do j=1, nz_nc
      zdum(k) = zdum(k-1)*ten**del
      k=k+1
   enddo
!
!calculate pdf for the given x
   integral=zero
   z_jm1 = zdum(1)
!
!calculate pdf for r and phi distribtion of the given (x,z)
   r_jm1 = sqrt(xfine_mid_c(i)**2+zdum(1)**2)
   phi_jm1 = atan(zdum(1)/xfine_mid_c(i))
   call find_index(r_jm1, r1d_dum, n1d_dum, iim2, iim1, ii, iip1)
   pdfr_jm1 = pdf_r1d(iim1)
   call find_index(phi_jm1, phi, nphi, iim2, iim1, ii, iip1)
   pdfp_jm1 = pdf_phi(iim1)
!   pdfp_jm1 = two/pi
!calculate joint pdf
   h_jm1 = pdfr_jm1*pdfp_jm1/sqrt(xfine_mid_c(i)**2+z_jm1**2)
!
   do j=2, nz
      z_j = zdum(j)
!calculate pdf for r and phi distribtion of the given (x,z)
      r_j = sqrt(xfine_mid_c(i)**2+zdum(j)**2)
      phi_j = atan(zdum(j)/xfine_mid_c(i))
      call find_index(r_j, r1d_dum, n1d_dum, iim2, iim1, ii, iip1)
      pdfr_j = pdf_r1d(iim1)
      call find_index(phi_j, phi, nphi, iim2, iim1, ii, iip1)
      pdfp_j = pdf_phi(iim1)
!      pdfp_j = two/!pi
!calculate joint pdf
      h_j = pdfr_j*pdfp_j/sqrt(xfine_mid_c(i)**2+z_j**2)
!
!perform integration over z
      integral = integral + half*(h_jm1+h_j)*(z_j-z_jm1)
!
      z_jm1 = z_j
      h_jm1 = h_j
!
   enddo
!
   pdf_xfine_c(i) = integral
!
enddo
!
!
!
!calculate probability of finding x in corresponding bins
do i=1, nxfine_c-1
   p_xfine_c(i) = pdf_xfine_c(i)*(xfine_c(i+1)-xfine_c(i))
!   write(*,*) xfine_mid_c(i), pdf_xfine_c(i)*four, p_xfine_c(i)
enddo
!
prob_c = sum(p_xfine_c)
!
!normalize pdf to core interval
pdf_xfine_c=pdf_xfine_c/prob_c
!
!--------------------------non core pdf---------------------------------
!
nxfine_nc=200
allocate(xfine_nc(nxfine_nc), stat=err)
allocate(xfine_mid_nc(nxfine_nc-1), stat=err)
allocate(pdf_xfine_nc(nxfine_nc-1), stat=err)
allocate(p_xfine_nc(nxfine_nc-1), stat=err)
!
!non-core grid logarithmic
call grid_log(rmin, rmax, nxfine_nc, xfine_nc)
!
!create mid points
do i=1, nxfine_nc-1
   xfine_mid_nc(i) = (xfine_nc(i+1)+xfine_nc(i))/two
enddo
!
!
!
deallocate(zdum)
nz_c=100
nz_nc=100
nz=nz_c+nz_nc
allocate(zdum(nz), stat=err)
!
!calculate pdf for non-core x-coordinates (on half grid)
do i=1, nxfine_nc-1
!
!create z-coordinates for z-integration (obtain marignal distribution by integrating over z)
   zmin = zero
   zmax = sqrt(rmax**2-xfine_mid_nc(i)**2)
!
   cc=half
   grad = (rmin-zmin)/log10((cc+one)/cc)
   zinterc = (zmin*log10(cc+one)-rmin*log10(cc))/log10((cc+one)/cc)
   zdum(1) = zmin
   k=2
   do j=2, nz_c
      zdum(k) = grad * log10(float(j-1)/(nz_c-1)+cc) + zinterc
      k=k+1
   enddo
   del=log10(zmax/rmin)/float(nz_nc)
   do j=1, nz_nc
      zdum(k) = zdum(k-1)*ten**del
      k=k+1
   enddo
!
!
!calculate pdf for the given x
   integral=zero
   z_jm1 = zdum(1)
!
!calculate pdf for r and phi distribtion of the given (x,z)
   r_jm1 = sqrt(xfine_mid_nc(i)**2+zdum(1)**2)
   phi_jm1 = atan(zdum(1)/xfine_mid_nc(i))
   call find_index(r_jm1, r1d_dum, n1d_dum, iim2, iim1, ii, iip1)
   pdfr_jm1 = pdf_r1d(iim1)
   call find_index(phi_jm1, phi, nphi, iim2, iim1, ii, iip1)
   pdfp_jm1 = pdf_phi(iim1)
!   pdfp_jm1 = two/pi
!calculate joint pdf
   h_jm1 = pdfr_jm1*pdfp_jm1/sqrt(xfine_mid_nc(i)**2+z_jm1**2)
!
   do j=2, nz
      z_j = zdum(j)
!calculate pdf for r and phi distribtion of the given (x,z)
      r_j = sqrt(xfine_mid_nc(i)**2+zdum(j)**2)
      phi_j = atan(zdum(j)/xfine_mid_nc(i))
      call find_index(r_j, r1d_dum, n1d_dum, iim2, iim1, ii, iip1)
      pdfr_j = pdf_r1d(iim1)
      call find_index(phi_j, phi, nphi, iim2, iim1, ii, iip1)
      pdfp_j = pdf_phi(iim1)
!      pdfp_j = two/!pi
!calculate joint pdf
      h_j = pdfr_j*pdfp_j/sqrt(xfine_mid_nc(i)**2+z_j**2)
!
!perform integration over z
      integral = integral + half*(h_jm1+h_j)*(z_j-z_jm1)
!
      z_jm1 = z_j
      h_jm1 = h_j
!
   enddo
!
   pdf_xfine_nc(i) = integral
!
enddo
!
!
!
!calculate probability of finding x in corresponding bins
do i=1, nxfine_nc-1
   p_xfine_nc(i) = pdf_xfine_nc(i)*(xfine_nc(i+1)-xfine_nc(i))
enddo
!
!
prob_nc=sum(p_xfine_nc)
!
!normalize pdf to non-core interval
pdf_xfine_nc=pdf_xfine_nc/prob_nc
!
!------------complete pdf normalized to complete range [0.,rmax]-------------
!-----------------------(only required for checks)---------------------------
!
nxfine = nxfine_c + nxfine_nc -1 
allocate(xfine(nxfine), stat=err)
allocate(xfine_mid(nxfine-1), stat=err)
allocate(pdf_xfine(nxfine-1), stat=err)
allocate(p_xfine(nxfine-1), stat=err)
!
!
k=1
do i=1, nxfine_c
   xfine(k) = xfine_c(i)
   k=k+1
enddo
do i=2, nxfine_nc
   xfine(k) = xfine_nc(i)
   k=k+1
enddo
!
k=1
do i=1, nxfine_c-1
   xfine_mid(k) = xfine_mid_c(i)
   pdf_xfine(k) = pdf_xfine_c(i)*prob_c
   p_xfine(k) = p_xfine_c(i)
   k=k+1
enddo
do i=1, nxfine_nc-1
   xfine_mid(k) = xfine_mid_nc(i)
   pdf_xfine(k) = pdf_xfine_nc(i)*prob_nc
   p_xfine(k) = p_xfine_nc(i)
   k=k+1
enddo
!
pdf_xfine=pdf_xfine/sum(p_xfine)
!
!---------------------create final grid from zero to rmax---------------
!
!ndx_dum=ndx-1
!nx_final=ndx_dum
!nncx_final=ndx_dum-ncx_final+1   !+1 because rmin point occurrs twice
!write(*,*) nx_final, ncx_final, nncx_final
!
ndx_dum=ndx-1
nx_final=ndx_dum
ncx_final = ncx
nncx_final=ndx_dum-ncx_final+1   !+1 because rmin point occurrs twice
!
write(*,*) 'probability of finding x-coordinate inside core ', prob_c
write(*,*) 'probability of finding x-coordinate outside core', prob_nc
write(*,*) ' => optimum number of core points to match distribution    ', nint(ndx_dum*prob_c)
write(*,*) ' => optimum number of non-core points to match distribution', ndx_dum-nint(ndx_dum*prob_c)+1
write(*,*)
write(*,*) 'actually used number of core points    ', ncx_final
write(*,*) 'actually used number of non-core points', nncx_final
write(*,*) '# core points / # points    ', float(ncx_final)/float(nncx_final+ncx_final)
write(*,*) '# non-core points / # points', float(nncx_final)/float(nncx_final+ncx_final)
write(*,*)
!
!xdum only up to ndx-1 since phantom point will be calculated later on
allocate(xfinal_c(ncx_final), stat=err)
allocate(xfinal_nc(nncx_final), stat=err)
allocate(xfinal(nx_final), stat=err)
!
!-------------------distribute nxc_final grid points in core------------
!-------------------------(following the core pdf)----------------------
!
!probability that has to be matched for each interval
pc = one/float(ncx_final-1)
xfinal_c(1)=zero
!
pc_i = zero
do i=2, ncx_final
!probability that has to be matched
   pc_i = pc_i + pc
!
   delp = zero
   do j=1, nxfine_c-1
      xi =  (pc_i-delp)/pdf_xfine_c(j)  + xfine_c(j)
      xfinal_c(i) = xi
      if(xi.le.xfine_c(j+1)) then
!correct position found: go to next x_final-coordinate
         exit
      else
!go to next interval
         delp = delp + pdf_xfine_c(j)*(xfine_c(j+1)-xfine_c(j))
      endif
   enddo
enddo
!
!---------------distribute nncx_final grid points outside core----------
!-----------------------(following the non-core pdf)--------------------
!
!probability that has to be matched for each interval
pc = one/float(nncx_final-1)
xfinal_nc(1)=xfine_nc(1)
!
pc_i = zero
do i=2, nncx_final
!probability that has to be matched
   pc_i = pc_i + pc
!
   delp = zero
   do j=1, nxfine_nc-1
      xi =  (pc_i-delp)/pdf_xfine_nc(j)  + xfine_nc(j)
      xfinal_nc(i) = xi
      if(xi.le. xfine_nc(j+1)) then
!;correct position found: go to next x_final-coordinate
         exit
      else
!go to next interval
         delp = delp + pdf_xfine_nc(j)*(xfine_nc(j+1)-xfine_nc(j))
      endif
   enddo
enddo
!
!----------------------combined distribution----------------------------
!
k=1
do i=1, ncx_final
   xfinal(k) = xfinal_c(i)
   k=k+1
enddo
xfinal_c(ncx_final)=rmin
do i=1, nncx_final-1
   xfinal(k) = xfinal_nc(i+1)
   k=k+1
enddo
xfinal(nx_final)=rmax
!
!calculate final pdf (only for checks)
allocate(xfinal_mid(nx_final-1))
allocate(pdf_xfinal(nx_final-1))
allocate(p_xfinal(nx_final-1))
!
call calc_pdf(nx_final, xfinal, xfinal_mid, pdf_xfinal, p_xfinal)
!
!--------------------------print out everything-------------------------
!
open(1, file=output_dir//'/pdf_rgridf.dat')
   write(1,*) 'radial grid (input)'
   write(1,'(3a20)') 'r_mid', 'pdf(r)', 'p(r)'
   do i=1, n1d_dum-1 
      write(1,'(3es20.8)') r1d_mid(i), pdf_r1d(i), p_r1d(i)
   enddo
close(1)
!
open(1, file=output_dir//'/pdf_phigridf.dat')
   write(1,*) 'phi grid (input)'
   write(1,'(3a20)') 'phi_mid', 'pdf(phi)', 'p(phi)'
   do i=1, nphi-1 
      write(1,'(3es20.8)') phi_mid(i), pdf_phi(i), p_phi(i)
   enddo
close(1)

!stop 'go on in grid_spatial'
!
open(1, file=output_dir//'/pdf_xgridf_c.dat')
   write(1,*) 'optimum core-distribution of x-coordinates (w.r.t. radial grid)'
   write(1,'(3a20)') 'x_mid', 'pdf(x)', 'p(x)'
   do i=1, nxfine_c-1 
      write(1,'(3es20.8)') xfine_mid_c(i), pdf_xfine_c(i), p_xfine_c(i)
   enddo
close(1)
!
open(1, file=output_dir//'/pdf_xgridf_nc.dat')
   write(1,*) 'optimum non-core-distribution of x-coordinates (w.r.t. radial grid)'
   write(1,'(3a20)') 'x_mid', 'pdf(x)', 'p(x)'
   do i=1, nxfine_nc-1 
      write(1,'(3es20.8)') xfine_mid_nc(i), pdf_xfine_nc(i), p_xfine_nc(i)
   enddo
close(1)
!
open(1, file=output_dir//'/pdf_xgridf.dat')
   write(1,*) 'optimum (complete) distribution of x-coordinates w.r.t. radial grid'
   write(1,'(3a20)') 'x_mid', 'pdf(x)', 'p(x)'
   do i=1, nxfine-1 
      write(1,'(3es20.8)') xfine_mid(i), pdf_xfine(i), p_xfine(i)
   enddo
close(1)
!
open(1, file=output_dir//'/pdf_xgrid1.dat')
   write(1,*) 'actually used distribution at first step: gridxyz_optb'
   write(1,'(3a20)') 'x_mid', 'pdf(x)', 'p(x)'
   do i=1, nx_final-1 
      write(1,'(3es20.8)') xfinal_mid(i), pdf_xfinal(i), p_xfinal(i)
   enddo
close(1)
!
!-------------------create complete grid--------------------------------
!
allocate(x(ndxmax), stat=err)
   if(err.ne.0) stop 'error in gridxyz_opt3d: allocation'
!
x(1)=-xfinal(nx_final)
do i=2, nx_final
   x(i)=-xfinal(nx_final+2-i)
enddo
!
do i=1, nx_final
   x(nx_final + i) = xfinal(i)
enddo
!
x(ndxmax)=xfinal(nx_final)
!
do i=1, ndxmax
   if(abs(x(i)-one).lt.1.d-8) x(i)=one
enddo
!
!-------------calculate pdf for z-coordinates on a fine grid------------
!---------------1. pdf for core points (normalized to core)-------------
!---------------2. pdf for non-core points (normlized to non-core)------
!---------------3. complete pdf (normalized to complete range)----------
!
!----------------------------core pdf-----------------------------------
!
nzfine_c=100
allocate(zfine_c(nzfine_c), stat=err)
allocate(zfine_mid_c(nzfine_c-1), stat=err)
allocate(pdf_zfine_c(nzfine_c-1), stat=err)
allocate(p_zfine_c(nzfine_c-1), stat=err)
!
!core grid confined towards rmin
cc=half
grad = rmin/log10((cc+one)/cc)
xinterc = -rmin*log10(cc)/log10((cc+one)/cc)
do i=1, nzfine_c
   zfine_c(i) = grad * log10(float(i)/(nzfine_c)+cc) + xinterc
enddo
zfine_c(1)=zero
zfine_c(nzfine_c)=rmin
!
!create mid points
do i=1, nzfine_c-1
   zfine_mid_c(i) = (zfine_c(i+1)+zfine_c(i))/two
enddo
!
!
!calculate pdf for core coordinates
nx_c=200
nx_nc=400
nx=nx_c+nx_nc
allocate(xdum(nx), stat=err)
!
do i=1, nzfine_c-1
!
!create x-coordinates for x-integration (obtain marignal distribution by integrating over x)
!   write(*,*) zfine_mid_c(i)
   xmin = sqrt(rmin**2-zfine_mid_c(i)**2)
   xmax = sqrt(rmax**2-zfine_mid_c(i)**2)
!   write(*,*) xmin, xmax
!
   cc=half
   grad = (rmin-xmin)/log10((cc+one)/cc)
   xinterc = (xmin*log10(cc+one)-rmin*log10(cc))/log10((cc+one)/cc)
   k=1
   do j=1, nx_c
      xdum(k) = grad * log10(float(j-1)/(nx_c-1)+cc) + xinterc
      k=k+1
   enddo
   del=log10(xmax/rmin)/float(nx_nc)
   do j=1, nx_nc
      xdum(k) = xdum(k-1)*ten**del
      k=k+1
   enddo
!   write(*,*) xdum
!
!calculate pdf for the given x
   integral=zero
   x_jm1 = xdum(1)
!
!calculate pdf for r and phi distribtion of the given (x,z)
   r_jm1 = sqrt(zfine_mid_c(i)**2+xdum(1)**2)
   phi_jm1 = atan(zfine_mid_c(i)/xdum(1))
!   write(*,*) r_jm1, phi_jm1
   call find_index(r_jm1, r1d_dum, n1d_dum, iim2, iim1, ii, iip1)
   pdfr_jm1 = pdf_r1d(iim1)
   call find_index(phi_jm1, phi, nphi, iim2, iim1, ii, iip1)
   pdfp_jm1 = pdf_phi(iim1)
!   write(*,*) pdfr_jm1, pdfp_jm1
!   pdfp_jm1 = two/pi
!calculate joint pdf
   h_jm1 = pdfr_jm1*pdfp_jm1/sqrt(zfine_mid_c(i)**2+x_jm1**2)
!   write(*,*) h_jm1
!
   do j=2, nx
      x_j = xdum(j)
!calculate pdf for r and phi distribtion of the given (x,z)
      r_j = sqrt(zfine_mid_c(i)**2+xdum(j)**2)
      phi_j = atan(zfine_mid_c(i)/xdum(j))
      call find_index(r_j, r1d_dum, n1d_dum, iim2, iim1, ii, iip1)
      pdfr_j = pdf_r1d(iim1)
      call find_index(phi_j, phi, nphi, iim2, iim1, ii, iip1)
      pdfp_j = pdf_phi(iim1)
!      pdfp_j = two/!pi
!calculate joint pdf
      h_j = pdfr_j*pdfp_j/sqrt(zfine_mid_c(i)**2+x_j**2)
!
!perform integration over z
      integral = integral + half*(h_jm1+h_j)*(x_j-x_jm1)
!
      x_jm1 = x_j
      h_jm1 = h_j 
!
   enddo
!
   pdf_zfine_c(i) = integral
!
enddo
!
!
!
!calculate probability of finding x in corresponding bins
do i=1, nzfine_c-1
   p_zfine_c(i) = pdf_zfine_c(i)*(zfine_c(i+1)-zfine_c(i))
!   write(*,*) xfine_mid_c(i), pdf_xfine_c(i)*four, p_xfine_c(i)
enddo
!
prob_c = sum(p_zfine_c)
!
!normalize pdf to core interval
pdf_zfine_c=pdf_zfine_c/prob_c
!
!--------------------------non core pdf---------------------------------
!
nzfine_nc=200
allocate(zfine_nc(nzfine_nc), stat=err)
allocate(zfine_mid_nc(nzfine_nc-1), stat=err)
allocate(pdf_zfine_nc(nzfine_nc-1), stat=err)
allocate(p_zfine_nc(nzfine_nc-1), stat=err)
!
!non-core grid logarithmic
call grid_log(rmin, rmax, nzfine_nc, zfine_nc)
!
!create mid points
do i=1, nzfine_nc-1
   zfine_mid_nc(i) = (zfine_nc(i+1)+zfine_nc(i))/two
enddo
!
!
!
deallocate(xdum)
nx_c=100
nx_nc=100
nx=nx_c+nx_nc
allocate(xdum(nx), stat=err)
!
!calculate pdf for non-core x-coordinates (on half grid)
do i=1, nzfine_nc-1
!
!create x-coordinates for x-integration (obtain marignal distribution by integrating over x)
   xmin = zero
   xmax = sqrt(rmax**2-zfine_mid_nc(i)**2)
!
   cc=half
   grad = (rmin-xmin)/log10((cc+one)/cc)
   xinterc = (xmin*log10(cc+one)-rmin*log10(cc))/log10((cc+one)/cc)
   xdum(1) = xmin
   k=2
   do j=2, nx_c
      xdum(k) = grad * log10(float(j-1)/(nx_c-1)+cc) + xinterc
      k=k+1
   enddo
   del=log10(xmax/rmin)/float(nx_nc)
   do j=1, nx_nc
      xdum(k) = xdum(k-1)*ten**del
      k=k+1
   enddo
!
!
!
!calculate pdf for the given x
   integral=zero
   x_jm1 = xdum(1)
!
!calculate pdf for r and phi distribtion of the given (x,z)
   r_jm1 = sqrt(zfine_mid_nc(i)**2+xdum(1)**2)
!   phi_jm1 = atan(zfine_mid_nc(i)/xdum(1))
   phi_jm1 = pi/two !avoid division by zero
   call find_index(r_jm1, r1d_dum, n1d_dum, iim2, iim1, ii, iip1)
   pdfr_jm1 = pdf_r1d(iim1)
   call find_index(phi_jm1, phi, nphi, iim2, iim1, ii, iip1)
   pdfp_jm1 = pdf_phi(iim1)
!   pdfp_jm1 = two/pi
!calculate joint pdf
   h_jm1 = pdfr_jm1*pdfp_jm1/sqrt(zfine_mid_nc(i)**2+x_jm1**2)
!
   do j=2, nx
      x_j = xdum(j)
!calculate pdf for r and phi distribtion of the given (x,z)
      r_j = sqrt(zfine_mid_nc(i)**2+xdum(j)**2)
      phi_j = atan(zfine_mid_nc(i)/xdum(j))
      call find_index(r_j, r1d_dum, n1d_dum, iim2, iim1, ii, iip1)
      pdfr_j = pdf_r1d(iim1)
      call find_index(phi_j, phi, nphi, iim2, iim1, ii, iip1)
      pdfp_j = pdf_phi(iim1)
!      pdfp_j = two/!pi
!calculate joint pdf
      h_j = pdfr_j*pdfp_j/sqrt(zfine_mid_nc(i)**2+x_j**2)
!
!perform integration over z
      integral = integral + half*(h_jm1+h_j)*(x_j-x_jm1)
!
      x_jm1 = x_j
      h_jm1 = h_j
!
   enddo
!
   pdf_zfine_nc(i) = integral
!
enddo
!
!
!
!calculate probability of finding x in corresponding bins
do i=1, nzfine_nc-1
   p_zfine_nc(i) = pdf_zfine_nc(i)*(zfine_nc(i+1)-zfine_nc(i))
enddo
!
!
prob_nc=sum(p_zfine_nc)
!
!normalize pdf to non-core interval
pdf_zfine_nc=pdf_zfine_nc/prob_nc
!
!------------complete pdf normalized to complete range [0.,rmax]-------------
!-----------------------(only required for checks)---------------------------
!
nzfine = nzfine_c + nzfine_nc -1 
allocate(zfine(nzfine), stat=err)
allocate(zfine_mid(nzfine-1), stat=err)
allocate(pdf_zfine(nzfine-1), stat=err)
allocate(p_zfine(nzfine-1), stat=err)
!
!
k=1
do i=1, nzfine_c
   zfine(k) = zfine_c(i)
   k=k+1
enddo
do i=2, nzfine_nc
   zfine(k) = zfine_nc(i)
   k=k+1
enddo
!
k=1
do i=1, nzfine_c-1
   zfine_mid(k) = zfine_mid_c(i)
   pdf_zfine(k) = pdf_zfine_c(i)*prob_c
   p_zfine(k) = p_zfine_c(i)
   k=k+1
enddo
do i=1, nzfine_nc-1
   zfine_mid(k) = zfine_mid_nc(i)
   pdf_zfine(k) = pdf_zfine_nc(i)*prob_nc
   p_zfine(k) = p_zfine_nc(i)
   k=k+1
enddo
!
pdf_zfine=pdf_zfine/sum(p_zfine)
!
!---------------------create final grid from zero to rmax---------------
!
ndz_dum=ndz-1
nz_final=ndz_dum
ncz_final = ncz
nncz_final=ndz_dum-ncz_final+1   !+1 because rmin point occurrs twice
!
write(*,*) 'probability of finding z-coordinate inside core ', prob_c
write(*,*) 'probability of finding z-coordinate outside core', prob_nc
write(*,*) ' => optimum number of core points to match distribution    ', nint(ndz_dum*prob_c)
write(*,*) ' => optimum number of non-core points to match distribution', ndz_dum-nint(ndz_dum*prob_c)+1
write(*,*)
write(*,*) 'actually used number of core points    ', ncz_final
write(*,*) 'actually used number of non-core points', nncz_final
write(*,*) '# core points / # points    ', float(ncz_final)/float(nncz_final+ncz_final)
write(*,*) '# non-core points / # points', float(nncz_final)/float(nncz_final+ncz_final)
write(*,*)
!
!zdum only up to ndx-1 since phantom point will be calculated later on
allocate(zfinal_c(ncz_final), stat=err)
allocate(zfinal_nc(nncz_final), stat=err)
allocate(zfinal(nz_final), stat=err)
!
!-------------------distribute nxz_final grid points in core------------
!-------------------------(following the core pdf)----------------------
!
!probability that has to be matched for each interval
pc = one/float(ncz_final-1)
zfinal_c(1)=zero
!
pc_i = zero
do i=2, ncz_final
!probability that has to be matched
   pc_i = pc_i + pc
!
   delp = zero
   do j=1, nzfine_c-1
      zi =  (pc_i-delp)/pdf_zfine_c(j)  + zfine_c(j)
      zfinal_c(i) = zi
      if(zi.le.zfine_c(j+1)) then
!correct position found: go to next z_final-coordinate
         exit
      else
!go to next interval
         delp = delp + pdf_zfine_c(j)*(zfine_c(j+1)-zfine_c(j))
      endif
   enddo
enddo
!
!---------------distribute nncz_final grid points outside core----------
!-----------------------(following the non-core pdf)--------------------
!
!probability that has to be matched for each interval
pc = one/float(nncz_final-1)
zfinal_nc(1)=zfine_nc(1)
!
pc_i = zero
do i=2, nncz_final
!probability that has to be matched
   pc_i = pc_i + pc
!
   delp = zero
   do j=1, nzfine_nc-1
      zi =  (pc_i-delp)/pdf_zfine_nc(j)  + zfine_nc(j)
      zfinal_nc(i) = zi
      if(zi.le. zfine_nc(j+1)) then
!;correct position found: go to next z_final-coordinate
         exit
      else
!go to next interval
         delp = delp + pdf_zfine_nc(j)*(zfine_nc(j+1)-zfine_nc(j))
      endif
   enddo
enddo
!
!----------------------combined distribution----------------------------
!
k=1
do i=1, ncz_final
   zfinal(k) = zfinal_c(i)
   k=k+1
enddo
zfinal_c(ncz_final)=rmin
do i=1, nncz_final-1
   zfinal(k) = zfinal_nc(i+1)
   k=k+1
enddo
zfinal(nz_final)=rmax
!
!calculate final pdf (only for checks)
allocate(zfinal_mid(nz_final-1))
allocate(pdf_zfinal(nz_final-1))
allocate(p_zfinal(nz_final-1))
!
call calc_pdf(nz_final, zfinal, zfinal_mid, pdf_zfinal, p_zfinal)
!
!--------------------------print out everything-------------------------
!
open(1, file=output_dir//'/pdf_zgridf_c.dat')
   write(1,*) 'optimum core-distribution of z-coordinates (w.r.t. radial grid)'
   write(1,'(3a20)') 'z_mid', 'pdf(z)', 'p(z)'
   do i=1, nzfine_c-1 
      write(1,'(3es20.8)') zfine_mid_c(i), pdf_zfine_c(i), p_zfine_c(i)
   enddo
close(1)
!
open(1, file=output_dir//'/pdf_zgridf_nc.dat')
   write(1,*) 'optimum non-core-distribution of z-coordinates (w.r.t. radial grid)'
   write(1,'(3a20)') 'z_mid', 'pdf(z)', 'p(z)'
   do i=1, nzfine_nc-1 
      write(1,'(3es20.8)') zfine_mid_nc(i), pdf_zfine_nc(i), p_zfine_nc(i)
   enddo
close(1)
!
open(1, file=output_dir//'/pdf_zgridf.dat')
   write(1,*) 'optimum (complete) distribution of z-coordinates w.r.t. radial grid'
   write(1,'(3a20)') 'z_mid', 'pdf(z)', 'p(z)'
   do i=1, nzfine-1 
      write(1,'(3es20.8)') zfine_mid(i), pdf_zfine(i), p_zfine(i)
   enddo
close(1)
!
open(1, file=output_dir//'/pdf_zgrid1.dat')
   write(1,*) 'actually used distribution at first step: gridxyz_optb'
   write(1,'(3a20)') 'z_mid', 'pdf(z)', 'p(z)'
   do i=1, nz_final-1 
      write(1,'(3es20.8)') zfinal_mid(i), pdf_zfinal(i), p_zfinal(i)
   enddo
close(1)
!
!-------------------create complete z grid------------------------------
!
allocate(z(ndzmax), stat=err)
   if(err.ne.0) stop 'error in gridxyz_opt3d: allocation'
!
z(1)=-zfinal(nz_final)
do i=2, nz_final
   z(i)=-zfinal(nz_final+2-i)
enddo
!
do i=1, nz_final
   z(nz_final + i) = zfinal(i)
enddo
!
z(ndzmax)=zfinal(nz_final)
!
do i=1, ndzmax
   if(abs(z(i)-one).lt.1.d-8) z(i)=one
enddo
!
!-------------------------------create y grid---------------------------
!
!!check if dimensions of x and z grid are compatible
!if(ndzmax.ne.ndxmax.or.ndz.ne.ndx) stop 'dimensions of z and x not compatible'
!!set x equal to z-grid
!z=x
!
allocate(y(ndymax), stat=err)
   if(err.ne.0) stop 'error in gridxyz_opt3d: allocation'
!check if dimensions of y and x grid are compatible
if(ndymax.ne.ndxmax.or.ndy.ne.ndx) stop 'dimensions of y and x not compatible'
!set y equal to x-grid
y=x
!
!deallocate 1d arrays, since not needed anymore
deallocate(r1d_dum)
!
!write(*,*) x
!write(*,*) 
!write(*,*) y
!write(*,*) 
!write(*,*) z
!stop 'go on in gridxyz_opt3d'
!
!
end subroutine gridxyz_opt3d
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine recalc_gridxyz
!
!---------------set equidistant steps on coordinate axes----------------
!-----if delx, dely, delz is larger than maximum allowed increments-----
!
!to ensure that the original grid properties are conserved, 
!   all other grid points are scaled up by a certain constant
!   (actually, only the delta-x,y,z are scaled up)
!
use prog_type
use fund_const
use mod_directories, only: output_dir
use dime3d, only: ndx, ndy, ndz, ndxmax, ndymax, ndzmax, delx_max, &
                  dely_max, delz_max, x, y, z
use params_input, only: rmax
use mod_grid, only: calc_pdf
!
implicit none
!
! ... local scalars
integer(i4b) :: i, j, err
integer(i4b) :: indx, indx_nc
integer(i4b) :: ndxmax_dum, ndymax_dum, ndzmax_dum
real(dp) :: delx, dely, delz
real(dp) :: const
real(dp) :: x_max, xnew_im1, y_max, ynew_im1, z_max, znew_im1
!
! ... local arrays
real(dp), dimension(:), allocatable :: x_dum, y_dum, z_dum, delx_dum, dely_dum, delz_dum, &
                                       xdum_mid, pdf_xdum, p_xdum
!
!------------calculate minimum possible radial increment----------------
!-----------------(if equidistant grid was used)------------------------
!
delx=two*rmax/(ndxmax-1)
dely=two*rmax/(ndymax-1)
delz=two*rmax/(ndzmax-1)
if(delx.gt.delx_max) stop 'error in recalc_gridxyz: ndxmax needs to be increased'
if(dely.gt.dely_max) stop 'error in recalc_gridxyz: ndymax needs to be increased'
if(delz.gt.delz_max) stop 'error in recalc_gridxyz: ndzmax needs to be increased'
!
!--------work only on positive xyz-axis (will be mirrored later on)-----
!
ndxmax_dum=ndxmax/2+1
ndymax_dum=ndymax/2+1
ndzmax_dum=ndzmax/2+1
!
allocate(x_dum(ndxmax_dum), stat=err)
   if(err.ne.0) stop 'allocation error recalc_gridxyz: x_dum'
allocate(y_dum(ndymax_dum), stat=err)
   if(err.ne.0) stop 'allocation error recalc_gridxyz: y_dum'
allocate(z_dum(ndzmax_dum), stat=err)
   if(err.ne.0) stop 'allocation error recalc_gridxyz: z_dum'
allocate(delx_dum(ndxmax_dum), stat=err)
   if(err.ne.0) stop 'allocation error recalc_gridxyz: delx_dum'
allocate(dely_dum(ndymax_dum), stat=err)
   if(err.ne.0) stop 'allocation error recalc_gridxyz: dely_dum'
allocate(delz_dum(ndzmax_dum), stat=err)
   if(err.ne.0) stop 'allocation error recalc_gridxyz: delz_dum'
!
x_dum=x(ndxmax/2+1:ndxmax)
y_dum=y(ndymax/2+1:ndymax)
z_dum=z(ndzmax/2+1:ndzmax)
!
!-----------------------reset x-grid------------------------------------
!
!loop as long as needed to ensure that grid points which are inside the up-scaling region
!   have delx/dely/delz lt delx_max also after up-scaling!
!
x_max=x_dum(ndxmax_dum-1)
!
outerx: do j=1, 100

!find index from where increments are too large
   indx=0
   delx=delx_max
   do i=2, ndxmax_dum
      delx=x_dum(i)-x_dum(i-1)
      if(delx.gt.delx_max) then
         indx=i
         exit
      endif
   enddo
!
!exit outer loop if no delx exists which is gt delx_max
   if(indx.eq.0) exit outerx
!exit outer loop if delx .gt. delx_max only due to precision errors
   if(abs(delx-delx_max).lt.1.d-14) exit outerx
!
!---------------free one point inside the core if j is odd--------------
!
  if(mod(j,2).ne.0) then
!find number of core points
      do i=1, ndxmax_dum
         if(x_dum(i).ge.one) then
            indx_nc=i-1
            exit
         endif
      enddo
!
      if(indx_nc.le.2) stop 'error in recalc_gridxyz: not enough core points' 
!
!scale up all core points such that x_dum(indx_nc)=one
      const=x_dum(indx_nc+1)/x_dum(indx_nc)
      x_dum(1:indx_nc)=const*x_dum(1:indx_nc)
!
!shift all points from indx_nc+1 to indx-1 one to the left
      do i=indx_nc+1, indx-1 
         x_dum(i-1)=x_dum(i)
      enddo
!
!make equidistant steps from indx-1 to ndxmax_dum-1 with steps delx_max
      do i=indx-1, ndxmax_dum-1
         x_dum(i) = x_dum(i-1) + delx_max
      enddo
!
!scale up/down all non-core points in order to have x_dum(ndxmax_dum-1)=xmax
      delx_dum(1)=zero
      do i=2, ndxmax_dum
         delx_dum(i)=x_dum(i)-x_dum(i-1)
      enddo
      const=(x_max-x_dum(indx_nc))/sum(delx_dum(indx_nc+1:ndxmax_dum-1))
!
      do i=indx_nc+1, ndxmax_dum-1
         x_dum(i) = x_dum(i-1) + const*delx_dum(i)
      enddo

!      do i=1, ndxmax_dum
!         write(*,*) i, x_dum(i)
!      enddo
!
  endif
!
!------------------scale up outer points if j is even-------------------
!
  if(mod(j,2).eq.0) then
      xnew_im1 = x_max - delx_max*(ndxmax_dum-1-indx+1)
      delx_dum(1)=zero
      do i=2, ndxmax_dum
         delx_dum(i)=x_dum(i)-x_dum(i-1)
      enddo

      const=(xnew_im1-x_dum(indx_nc))/sum(delx_dum(indx_nc+1:indx-1))
!
!scale down points from indx_nc+1 to indx-2
      do i=indx_nc+1, indx-2
         x_dum(i)=x_dum(i-1)+const*delx_dum(i)
      enddo
!
!from that on, make equidistant grid (and forget about outermost point...)
      do i=indx-1, ndxmax_dum-1
         x_dum(i) = xnew_im1 + (i-indx+1)*delx_max
      enddo
!
   endif
!
enddo outerx
!
!
!
if(j.eq.101) then
   write(*,*) 'error in recalc_gridxyz: upscaling of x-axis did not work:'
   write(*,*) '#', 'x_dum', 'delx'
   do i=2, ndxmax_dum
      write(*,*) i, x_dum(i), x_dum(i)-x_dum(i-1)
   enddo
   write(*,*)
   stop 
endif
!
!now, overwrite acutally used x-coordinate from -x_dum to x_dum
do i=1, ndxmax_dum
   x(i) = -x_dum(ndxmax_dum+1-i)
   x(ndxmax_dum+i-1) = x_dum(i)
enddo
!
!
!calculate pdf (only for checks)
allocate(xdum_mid(ndxmax_dum-2))
allocate(pdf_xdum(ndxmax_dum-2))
allocate(p_xdum(ndxmax_dum-2))
!
call calc_pdf(ndxmax_dum-1, x_dum(1:ndxmax_dum-1), xdum_mid, pdf_xdum, p_xdum)
!
open(1, file=output_dir//'/pdf_xgrid2.dat')
   write(1,*) 'actually used distribution at second step: recalc_gridxyz'
   write(1,'(3a20)') 'x_mid', 'pdf(x)', 'p(x)'
   do i=1, ndxmax_dum-2
      write(1,'(3es20.8)') xdum_mid(i), pdf_xdum(i), p_xdum(i)
   enddo
   close(1)
!
!-----------------------reset y-grid------------------------------------
!
!loop as long as needed to ensure that grid points which are inside the up-scaling region
!   have delx/dely/delz lt dely_max also after up-scaling!
!
y_max=y_dum(ndymax_dum-1)
!
outery: do j=1, 100

!find index from where increments are too large
   indx=0
   dely=dely_max
   do i=2, ndymax_dum
      dely=y_dum(i)-y_dum(i-1)
      if(dely.gt.dely_max) then
         indx=i
         exit
      endif
   enddo
!
!exit outer loop if no dely exists which is gt dely_max
   if(indx.eq.0) exit outery
!exit outer loop if dely .gt. dely_max only due to precision errors
   if(abs(dely-dely_max).lt.1.d-14) exit outery
!
!---------------free one point inside the core if j is odd--------------
!
   if(mod(j,2).ne.0) then
!find number of core points
      do i=1, ndymax_dum
         if(y_dum(i).ge.one) then
            indx_nc=i-1
            exit
         endif
      enddo
!
!scale up all core points such that y_dum(indx_nc)=one
      const=y_dum(indx_nc+1)/y_dum(indx_nc)
      y_dum(1:indx_nc)=const*y_dum(1:indx_nc)
!
!shift all points from indx_nc+1 to indx-1 one to the left
      do i=indx_nc+1, indx-1 
         y_dum(i-1)=y_dum(i)
      enddo
!
!make equidistant steps from indx-2 to ndymax_dum-1 with steps dely_max
      do i=indx-1, ndymax_dum-1
         y_dum(i) = y_dum(i-1) + dely_max
      enddo
!
!scale up/down all non-core points in order to have y_dum(ndymax_dum-1)=ymax
      dely_dum(1)=zero
      do i=2, ndymax_dum
         dely_dum(i)=y_dum(i)-y_dum(i-1)
      enddo
      const=(y_max-y_dum(indx_nc))/sum(dely_dum(indx_nc+1:ndymax_dum-1))
!
      do i=indx_nc+1, ndymax_dum-1
         y_dum(i) = y_dum(i-1) + const*dely_dum(i)
      enddo
!
   endif
!
!------------------scale up outer points if j is even-------------------
!
   if(mod(j,2).eq.0) then
      ynew_im1 = y_max - dely_max*(ndymax_dum-1-indx+1)
      dely_dum(1)=zero
      do i=2, ndymax_dum
         dely_dum(i)=y_dum(i)-y_dum(i-1)
      enddo

      const=(ynew_im1-y_dum(indx_nc))/sum(dely_dum(indx_nc+1:indx-1))
!
!scale down points from indx_nc+1 to indx-2
      do i=indx_nc+1, indx-2
         y_dum(i)=y_dum(i-1)+const*dely_dum(i)
      enddo
!
!from that on, make equidistant grid (and forget about outermost point...)
      do i=indx-1, ndymax_dum-1
         y_dum(i) = ynew_im1 + (i-indx+1)*dely_max
      enddo
!
   endif
!
enddo outery
!
if(j.eq.101) then
   write(*,*) 'error in recalc_gridxyz: upscaling of y-axis did not work:'
   write(*,'(a5, 2a20)') '#', 'y_dum', 'dely'
   do i=2, ndymax_dum
      write(*,'(i5, 2e20.8)') i, y_dum(i), y_dum(i)-y_dum(i-1)
   enddo
   write(*,*)
   stop 
endif
!
!now, overwrite acutally used y-coordinate from -y_dum to y_dum
do i=1, ndymax_dum
   y(i) = -y_dum(ndymax_dum+1-i)
   y(ndymax_dum+i-1) = y_dum(i)
enddo
!
!-----------------------reset z-grid------------------------------------
!
!loop as long as needed to ensure that grid points which are inside the up-scaling region
!   have delx/dely/delz lt delz_max also after up-scaling!
!
z_max=z_dum(ndzmax_dum-1)
!
outerz: do j=1, 100

!find index from where increments are too large
   indx=0
   delz=delz_max
   do i=2, ndzmax_dum
      delz=z_dum(i)-z_dum(i-1)
      if(delz.gt.delz_max) then
         indx=i
         exit
      endif
   enddo
!
!exit outer loop if no delz exists which is gt delz_max
   if(indx.eq.0) exit outerz
!exit outer loop if delz .gt. delz_max only due to precision errors
   if(abs(delz-delz_max).lt.1.d-14) exit outerz
!
!---------------free one point inside the core if j is odd--------------
!
   if(mod(j,2).ne.0) then
!find number of core points
      do i=1, ndzmax_dum
         if(z_dum(i).ge.one) then
            indx_nc=i-1
            exit
         endif
      enddo
!
!scale up all core points such that z_dum(indx_nc)=one
      const=z_dum(indx_nc+1)/z_dum(indx_nc)
      z_dum(1:indx_nc)=const*z_dum(1:indx_nc)
!
!shift all points from indx_nc+1 to indx-1 one to the left
      do i=indx_nc+1, indx-1 
         z_dum(i-1)=z_dum(i)
      enddo
!
!make equidistant steps from indx-2 to ndzmax_dum-1 with steps delz_max
      do i=indx-1, ndzmax_dum-1
         z_dum(i) = z_dum(i-1) + delz_max
      enddo
!
!scale up/down all non-core points in order to have z_dum(ndzmax_dum-1)=zmax
      delz_dum(1)=zero
      do i=2, ndzmax_dum
         delz_dum(i)=z_dum(i)-z_dum(i-1)
      enddo
      const=(z_max-z_dum(indx_nc))/sum(delz_dum(indx_nc+1:ndzmax_dum-1))
!
      do i=indx_nc+1, ndzmax_dum-1
         z_dum(i) = z_dum(i-1) + const*delz_dum(i)
      enddo
!
   endif
!
!------------------scale up outer points if j is even-------------------
!
   if(mod(j,2).eq.0) then
      znew_im1 = z_max - delz_max*(ndzmax_dum-1-indx+1)
      delz_dum(1)=zero
      do i=2, ndzmax_dum
         delz_dum(i)=z_dum(i)-z_dum(i-1)
      enddo
!
      const=(znew_im1-z_dum(indx_nc))/sum(delz_dum(indx_nc+1:indx-1))
!
!scale down points from indx_nc+1 to indx-2
      do i=indx_nc+1, indx-2
         z_dum(i)=z_dum(i-1)+const*delz_dum(i)
      enddo
!
!from that on, make equidistant grid (and forget about outermost point...)
      do i=indx-1, ndzmax_dum-1
         z_dum(i) = znew_im1 + (i-indx+1)*delz_max
      enddo
!
   endif
!
enddo outerz
!
if(j.eq.101) then
   write(*,*) 'error in recalc_gridxyz: upscaling of z-axis did not work:'
   write(*,'(a5, 2a20)') '#', 'z_dum', 'delz'
   do i=2, ndzmax_dum
      write(*,'(i5, 2e20.8)') i, z_dum(i), z_dum(i)-z_dum(i-1)
   enddo
   write(*,*)
   stop 
endif
!
!now, overwrite acutally used z-coordinate from -z_dum to z_dum
do i=1, ndzmax_dum
   z(i) = -z_dum(ndzmax_dum+1-i)
   z(ndzmax_dum+i-1) = z_dum(i)
enddo
!
!
!write(*,*) x
!write(*,*)
!x=(/ -12.797851d0, -12.000000d0, -11.300000d0, -10.600000d0, -9.9000000d0, -9.2000000d0, &
!-8.5000000d0, -7.8499020d0, -7.2586310d0, -6.7176547d0, -6.2136517d0, -5.7462186d0, -5.3033940d0, -4.8886315d0, &
!-4.4991898d0, -4.1246273d0, -3.7750043d0, -3.4481124d0, -3.1557507d0, -2.8917742d0, -2.6571460d0, -2.4499398d0, &
!-2.2658623d0, -2.1017464d0, -1.9572221d0, -1.8287458d0, -1.7128406d0, -1.6100105d0, -1.5179043d0, -1.4336432d0, &
!-1.3572476d0, -1.2877396d0, -1.2257713d0, -1.1692105d0, -1.1191722d0, -1.0783439d0, -1.0446928d0, &
!-1.0217302d0, -1.0055393d0, -1.0000000d0, -0.98450830d0, -0.96424246d0, -0.94137250d0, -0.91612835d0, -0.88866442d0, &
!-0.85948188d0, -0.82881014d0, -0.79667655d0, -0.76310841d0, -0.72820222d0, -0.69198234d0, -0.65415142d0, -0.61607401d0, -0.57573725d0, &
!-0.53487804d0, -0.49352992d0, -0.44990084d0, -0.40786581d0, -0.36673137d0, -0.32365376d0, -0.27946078d0, -0.23535626d0, &
!-0.19116861d0, -0.14363409d0, -0.096092413d0, -0.049473501d0, zero, 0.049473501d0, 0.096092413d0, 0.14363409d0, &
!0.19116861d0, 0.23535626d0, 0.27946078d0, 0.32365376d0, 0.36673137d0, 0.40786581d0, 0.44990084d0, 0.49352992d0, &
!0.53487804d0, 0.57573725d0, 0.61607401d0, 0.65415142d0, 0.69198234d0, 0.72820222d0, 0.76310841d0, 0.79667655d0, &
!0.82881014d0, 0.85948188d0, 0.88866442d0, 0.91612835d0, 0.94137250d0, 0.96424246d0, 0.98450830d0, 1.0000000d0, &
!1.0055393d0, 1.0217302d0, 1.0446928d0, 1.0783439d0, 1.1191722d0, 1.1692105d0, 1.2257713d0, 1.2877396d0, &
!1.3572476d0, 1.4336432d0, 1.5179043d0, 1.6100105d0, 1.7128406d0, 1.8287458d0, 1.9572221d0, 2.1017464d0, &
!2.2658623d0, 2.4499398d0, 2.6571460d0, 2.8917742d0, 3.1557507d0, 3.4481124d0, &
!3.7750043d0, 4.1246273d0, 4.4991898d0, 4.8886315d0, 5.3033940d0, 5.7462186d0, 6.2136517d0, 6.7176547d0, &
!7.2586310d0, 7.8499020d0, 8.5000000d0, 9.2000000d0, 9.9000000d0, 10.600000d0, 11.300000d0, 12.000000d0, 12.797851d0 /)
!y=x
!z=x
!write(*,*) x
!write(*,*) 
!write(*,*) z
!
!stop 'go on in recalc_gridxyz'



!
end subroutine recalc_gridxyz
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine recalc_vgridxyz
!
!------------set equidistant velocity steps on coordinate axes----------
!----------------------from stellar surface to rmax,--------------------
!----------------given a previous velocity stratification---------------
!
use prog_type
use fund_const
use dime3d, only: ndxmax, ndymax, ndzmax, x, y, z, velx3d, vely3d, velz3d
use dime1d, only: delv
use inf_reg, only: rmin
use params_input, only: vth_fiducial
use mod_interp1d, only: find_index, cube_mono_coeff, interpol_yp_spline, spline_coeff
use mod_grid, only: calc_pdf
!
implicit none
!
! ... local scalars
integer(i4b), parameter :: itmax=100000
integer(i4b) :: i, j, k, iit, err
integer(i4b) :: iim2, iim1, ii, iip1, ir, nr, jj, kk
real(dp) :: delv_x, cc, delx_new, delx_old, delv_old, delv_new, dv1, dv2
!
! ... local logicals
logical :: lx, lgood
!
! ... local arrays
real(dp), dimension(:), allocatable :: r_orig, r_old, r_new1, r_new2, r_new3, &
                                       vel_orig, vel_old, vel_new, &
                                       acoeff_vel, bcoeff_vel, ccoeff_vel, dcoeff_vel, &
                                       x_new, y_new, z_new
!
! ... local functions
!
!-------------------------------x-grid----------------------------------
!
ir=0
do i=1, ndxmax
!   if(abs(x(i)-rmin).lt.1.d-14) ir = i
   if(x(i).ge.rmin) then
      ir=i
      exit
   endif
enddo
!
if(ir.eq.0) stop 'error in recalc_vgridxyz: x(i)=rmin not found'
!
nr=ndxmax-ir
!
ii=ndxmax-nr
j=ndxmax/2+1
k=ndzmax/2+1
do i=2, nr-1
   dv1 = velx3d(ii,j,k)-velx3d(ii-1,j,k)
   dv2 = velx3d(ii+1,j,k)-velx3d(ii,j,k)
   if(dv1*dv2.lt.zero) then
      write(*,*) 'recalc_vgridxyz not performed for non-monotonic velocity field in x'
      return
   endif
   ii=ii+1
enddo
!
allocate(r_orig(nr), stat=err)
   if(err.ne.0) stop 'allocation error in recalc_vgridxyz'
allocate(r_old(nr), stat=err)
   if(err.ne.0) stop 'allocation error in recalc_vgridxyz'
allocate(r_new1(nr), stat=err)
   if(err.ne.0) stop 'allocation error in recalc_vgridxyz'
allocate(r_new2(nr), stat=err)
   if(err.ne.0) stop 'allocation error in recalc_vgridxyz'
allocate(r_new3(nr), stat=err)
   if(err.ne.0) stop 'allocation error in recalc_vgridxyz'
allocate(vel_orig(nr), stat=err)
   if(err.ne.0) stop 'allocation error in recalc_vgridxyz'
allocate(vel_old(nr), stat=err)
   if(err.ne.0) stop 'allocation error in recalc_vgridxyz'
allocate(vel_new(nr), stat=err)
   if(err.ne.0) stop 'allocation error in recalc_vgridxyz'
allocate(acoeff_vel(nr), stat=err)
   if(err.ne.0) stop 'allocation error in recalc_vgridxyz'
allocate(bcoeff_vel(nr), stat=err)
   if(err.ne.0) stop 'allocation error in recalc_vgridxyz'
allocate(ccoeff_vel(nr), stat=err)
   if(err.ne.0) stop 'allocation error in recalc_vgridxyz'
allocate(dcoeff_vel(nr), stat=err)
   if(err.ne.0) stop 'allocation error in recalc_vgridxyz'
allocate(x_new(ndxmax), stat=err)
   if(err.ne.0) stop 'allocation error in recalc_vgridxyz'
!
!set radial stratification
j=1
do i=ndxmax-nr, ndxmax-1
   r_old(j) = x(i)
   vel_old(j) = velx3d(i,ndymax/2+1,ndzmax/2+1)
   j=j+1
!   write(*,*) velx3d(i,ndymax/2+1,ndzmax/2+1), y(ndymax/2+1), z(ndzmax/2+1)
enddo
!
r_orig=r_old
vel_orig=vel_old
!
!spline-interpolation coefficients
!call spline_coeff(nr, r_orig, vel_orig, acoeff_vel, bcoeff_vel, ccoeff_vel, dcoeff_vel)
call cube_mono_coeff(nr, r_orig, vel_orig, acoeff_vel, bcoeff_vel, ccoeff_vel, dcoeff_vel)
!call lin_coeff(nr, r_orig, vel_orig, ccoeff_vel, dcoeff_vel)
!acoeff_vel=zero
!bcoeff_vel=zero

!
!check if radial stratification has too large velocity-steps
lx = .false.
do i=2, nr
   delv_x = abs(vel_old(i)-vel_old(i-1))
   if(delv_x.gt.delv) then
      lx=.true.
      exit
   endif
enddo
!
if(.not.lx) then
!input x-stratification is okay
   x_new = x
else
!calculate new x-stratification
!
   delv_x = abs(vel_old(nr)-vel_old(1))/(nr-1)
   delv_x = max(delv,delv_x)
!
   iterationx: do iit=1, itmax
      r_new1 = r_old
      vel_new = vel_old
      lgood = .true.
      do i=2, nr
         delv_old = vel_old(i)-vel_old(i-1)
         if(abs(delv_old)/delv_x.gt.1.01d0) then
!calculate new coordinates at (i-1) and at i
            delx_old = r_old(i) - r_old(i-1)
            delx_new = delv_x*delx_old/delv_old
            if(i.eq.2) then
               r_new1(i) = r_old(i) - delx_new
            elseif(i.eq.nr) then
               r_new1(i-1) = r_old(i-1) + delx_new
            else
               r_new1(i-1) = r_old(i-1) + (delx_old-delx_new)/two
               r_new1(i) = r_old(i) - (delx_old-delx_new)/two
            endif
!
!scale up rest of the grid towards left by keeping delx-ratios from original grid
            r_new2=r_new1
            do j=i-2, 1, -1
               r_new2(j) = r_new2(j+1) - abs((r_new2(j+2)-r_new2(j+1)))*(r_old(j+1)-r_old(j))/(r_old(j+2)-r_old(j+1))
            enddo
!scale up rest of the grid towards right by keeping delx-ratios from original grid
            do j=i+1, nr
               r_new2(j) = r_new2(j-1) + abs((r_new2(j-1)-r_new2(j-2)))*(r_old(j)-r_old(j-1))/(r_old(j-1)-r_old(j-2))
            enddo
!
!scale the grid to the left to ensure that left boundary is correct
            r_new3=r_new2
            cc = (r_new2(i-1)-r_old(1))/(r_new2(i-1)-r_new2(1))
            do j=i-2, 1, -1
            r_new3(j) = r_new3(j+1)-cc*(r_new2(j+1)-r_new2(j))
            enddo
!scale the grid to the left to ensure that left boundary is correct
            cc = (r_old(nr)-r_new2(i))/(r_new2(nr)-r_new2(i))
            do j=i+1, nr
               r_new3(j) = r_new3(j-1)+cc*(r_new2(j)-r_new2(j-1))
            enddo
!
!finally: interpolate velocities
            do j=1, nr
               call find_index(r_new3(j), r_orig, nr, iim2, iim1, ii, iip1)
!               vel_new(j) = vel_old(iim1) + (r_new3(j)-r_old(iim1))*(vel_old(ii)-vel_old(iim1))/(r_old(ii)-r_old(iim1))
                vel_new(j)=interpol_yp_spline(acoeff_vel(ii), bcoeff_vel(ii), ccoeff_vel(ii), dcoeff_vel(ii), r_orig(ii), r_new3(j))
            enddo
!
!            write(*,'(2i5,10es20.8)') iit, i, delv_old, delv_x
!            do j=1, nr
!               write(*,'(10es20.8)') r_orig(j), r_old(j), r_new1(j), r_new2(j), r_new3(j), vel_orig(j), vel_new(j)
!            enddo
!            write(*,*)
!
!overwrite original array and begin from start
            r_old=r_new3
            vel_old=vel_new
            lgood=.false.
!            if(abs(r_old(1)-one).gt.1.d-8) stop
            exit
         endif
!grid is okay and stop iteration
      enddo
!      write(*,*) iit, lgood
      if(lgood) exit iterationx
   enddo iterationx
!
!check if everything has worked
   if(.not.lgood) then
      write(*,*) 'error in recalc_vgridxyz: no solution found for nice x velocity-grid properties'
      write(*,*) '   either: increase number of iterations, or error in algorithm (e.g. for non-monotonic velocity fields?)'
      stop
   endif
   if(abs(r_old(1)-r_orig(1)).gt.1.d-10) stop 'error in recalc_vgridxyz: x minimum radius not conserved'
   if(abs(r_old(nr)-r_orig(nr)).gt.1.d-10) stop 'error in recalc_vgridxyz: x maximum radius not conserved'
   do i=2, nr
      if(r_old(i).le.r_old(i-1)) then
         write(*,*) 'error in recalc_vgridxyz: new radial x grid is not monotonic'
         write(*,*) '   either: trust original grid, and comment this procedure in main program'
         write(*,*) '   or find bug in algorithm (e.g. for non-monotonic velocity fields?'
         stop
      endif
   enddo

!   do i=2, nr
!      write(*,'(10es20.8)') r_old(i), vel_old(i), 1.d8/vth_fiducial*(one-0.9999d0/r_old(i))**0.5, vel_old(i)-vel_old(i-1), r_orig(i), vel_orig(i), 1.d8/vth_fiducial*(one-0.9999d0/r_orig(i))**0.5, vth_fiducial
!   enddo
!   write(*,*)
!   open(1, file='TRASH/test_vel.dat')
!   do i=1, 300
!      cc = one + (i-1)*0.1d0/299.d0
!      call find_index(cc, r_orig, nr, iim2, iim1, ii, iip1)
!      write(1,'(10es20.8)') cc, interpol_yp_spline(acoeff_vel(ii), bcoeff_vel(ii), ccoeff_vel(ii), dcoeff_vel(ii), r_orig(ii), cc), 1.d8/vth_fiducial*(one-0.9999d0/cc)**0.5, &
!                 interpol_yp(r_orig(iim1), r_orig(ii), vel_orig(iim1), vel_orig(ii), cc)
!   enddo
!   close(1)
!  stop 'go on in vgrid 2'
!
!if everything has worked: reset minimum and maximum radius to exact value, and reset x-coordinates
   r_old(1) = r_orig(1)
   r_old(nr) = r_orig(nr)
   x_new=x
   do i=2, nr+1
      x_new(i) = -r_old(nr-i+2)
   enddo
   j=1
   do i=ndxmax-nr, ndxmax-1
      x_new(i) = r_old(j)
      j=j+1
   enddo
!
!   do j=1, ndxmax
!      write(*,*) j, x(j)
!   enddo
!   stop
!   write(*,*)
!   open(1, file='TRASH/grid.dat')
!      do i=1, nr
!         write(1,'(4es20.8)') r_old(i), vel_new(i), r_orig(i), vel_orig(i)
!         write(*,'(4es20.8)') r_old(i), vel_new(i), r_orig(i), vel_orig(i)
!      enddo
!   close(1)
!
endif
!
deallocate(r_orig)
deallocate(r_old)
deallocate(r_new1)
deallocate(r_new2)
deallocate(r_new3)
deallocate(vel_orig)
deallocate(vel_old)
deallocate(vel_new)
deallocate(acoeff_vel)
deallocate(bcoeff_vel)
deallocate(ccoeff_vel)
deallocate(dcoeff_vel)
!
!x-coord done
!
!-------------------------------y-grid----------------------------------
!
ir=0
do i=1, ndymax
!   if(abs(y(i)-rmin).lt.1.d-14) ir = i
   if(y(i).ge.rmin) then
      ir=i
      exit
   endif
enddo
!
if(ir.eq.0) stop 'error in recalc_vgridxyz: y(i)=rmin not found'
!
nr=ndymax-ir
!
i=ndxmax/2+1
jj=ndymax-nr
k=ndzmax/2+1
do j=2, nr-1
   dv1 = vely3d(i,jj,k)-vely3d(i,jj-1,k)
   dv2 = vely3d(i,jj+1,k)-vely3d(i,jj,k)
   if(dv1*dv2.lt.zero) then
      write(*,*) 'recalc_vgridxyz not performed for non-monotonic velocity field in y'
      return
   endif
   jj=jj+1
enddo
!
allocate(r_orig(nr), stat=err)
   if(err.ne.0) stop 'allocation error in recalc_vgridxyz'
allocate(r_old(nr), stat=err)
   if(err.ne.0) stop 'allocation error in recalc_vgridxyz'
allocate(r_new1(nr), stat=err)
   if(err.ne.0) stop 'allocation error in recalc_vgridxyz'
allocate(r_new2(nr), stat=err)
   if(err.ne.0) stop 'allocation error in recalc_vgridxyz'
allocate(r_new3(nr), stat=err)
   if(err.ne.0) stop 'allocation error in recalc_vgridxyz'
allocate(vel_orig(nr), stat=err)
   if(err.ne.0) stop 'allocation error in recalc_vgridxyz'
allocate(vel_old(nr), stat=err)
   if(err.ne.0) stop 'allocation error in recalc_vgridxyz'
allocate(vel_new(nr), stat=err)
   if(err.ne.0) stop 'allocation error in recalc_vgridxyz'
allocate(acoeff_vel(nr), stat=err)
   if(err.ne.0) stop 'allocation error in recalc_vgridxyz'
allocate(bcoeff_vel(nr), stat=err)
   if(err.ne.0) stop 'allocation error in recalc_vgridxyz'
allocate(ccoeff_vel(nr), stat=err)
   if(err.ne.0) stop 'allocation error in recalc_vgridxyz'
allocate(dcoeff_vel(nr), stat=err)
   if(err.ne.0) stop 'allocation error in recalc_vgridxyz'
allocate(y_new(ndymax), stat=err)
   if(err.ne.0) stop 'allocation error in recalc_vgridxyz'
!
!set radial stratification
j=1
do i=ndymax-nr, ndymax-1
   r_old(j) = y(i)
   vel_old(j) = vely3d(ndxmax/2+1,i,ndzmax/2+1)
   j=j+1
enddo
!
r_orig=r_old
vel_orig=vel_old
!
!spline-interpolation coefficients
call spline_coeff(nr, r_orig, vel_orig, acoeff_vel, bcoeff_vel, ccoeff_vel, dcoeff_vel)
!
!check if radial stratification has too large velocity-steps
lx = .false.
do i=2, nr
   delv_x = abs(vel_old(i)-vel_old(i-1))
   if(delv_x.gt.delv) then
      lx=.true.
      exit
   endif
enddo
!
if(.not.lx) then
!input x-stratification is okay
   y_new = y
else
!calculate new x-stratification
!
   delv_x = (vel_old(nr)-vel_old(1))/(nr-1)
   delv_x = max(delv,delv_x)
!

   iterationy: do iit=1, itmax
      r_new1 = r_old
      vel_new = vel_old
      lgood = .true.
      do i=2, nr
         delv_old = vel_old(i)-vel_old(i-1)
         if(abs(delv_old)/delv_x.gt.1.01) then
!calculate new coordinates at (i-1) and at i
            delx_old = r_old(i) - r_old(i-1)
            delx_new = delv_x*delx_old/delv_old
            if(i.eq.2) then
               r_new1(i) = r_old(i) - delx_new
            elseif(i.eq.nr) then
               r_new1(i-1) = r_old(i-1) + delx_new
            else
               r_new1(i-1) = r_old(i-1) + (delx_old-delx_new)/two
               r_new1(i) = r_old(i) - (delx_old-delx_new)/two
            endif
!
!scale up rest of the grid towards left by keeping delx-ratios from original grid
            r_new2=r_new1
            do j=i-2, 1, -1
               r_new2(j) = r_new2(j+1) - abs((r_new2(j+2)-r_new2(j+1)))*(r_old(j+1)-r_old(j))/(r_old(j+2)-r_old(j+1))
            enddo
!scale up rest of the grid towards right by keeping delx-ratios from original grid
            do j=i+1, nr
               r_new2(j) = r_new2(j-1) + abs((r_new2(j-1)-r_new2(j-2)))*(r_old(j)-r_old(j-1))/(r_old(j-1)-r_old(j-2))
            enddo
!
!scale the grid to the left to ensure that left boundary is correct
            r_new3=r_new2
            cc = (r_new2(i-1)-r_old(1))/(r_new2(i-1)-r_new2(1))
            do j=i-2, 1, -1
            r_new3(j) = r_new3(j+1)-cc*(r_new2(j+1)-r_new2(j))
            enddo
!scale the grid to the left to ensure that left boundary is correct
            cc = (r_old(nr)-r_new2(i))/(r_new2(nr)-r_new2(i))
            do j=i+1, nr
               r_new3(j) = r_new3(j-1)+cc*(r_new2(j)-r_new2(j-1))
            enddo
!
!finally: interpolate velocities
            do j=1, nr
               call find_index(r_new3(j), r_orig, nr, iim2, iim1, ii, iip1)
!               vel_new(j) = vel_old(iim1) + (r_new3(j)-r_old(iim1))*(vel_old(ii)-vel_old(iim1))/(r_old(ii)-r_old(iim1))
                vel_new(j)=interpol_yp_spline(acoeff_vel(ii), bcoeff_vel(ii), ccoeff_vel(ii), dcoeff_vel(ii), r_orig(ii), r_new3(j))
            enddo
!
!            write(*,*) i
!            do j=1, nr
!               write(*,'(10es20.8)') r_orig(j), r_old(j), r_new1(j), r_new2(j), r_new3(j), vel_orig(j), vel_new(j)
!            enddo
!            write(*,*)
!
!overwrite original array and begin from start
            r_old=r_new3
            vel_old=vel_new
            lgood=.false.
!            if(abs(r_old(1)-one).gt.1.d-8) stop
            exit
         endif
!grid is okay and stop iteration
      enddo
!      write(*,*) iit, lgood
      if(lgood) exit iterationy
   enddo iterationy
!
!check if everything has worked
   if(.not.lgood) then
      write(*,*) 'error in recalc_vgridxyz: no solution found for nice y velocity-grid properties'
      write(*,*) '   either: increase number of iterations, or error in algorithm (e.g. for non-monotonic velocity fields?)'
      stop
   endif
   if(abs(r_old(1)-r_orig(1)).gt.1.d-10) stop 'error in recalc_vgridxyz: y minimum radius not conserved'
   if(abs(r_old(nr)-r_orig(nr)).gt.1.d-10) stop 'error in recalc_vgridxyz: y maximum radius not conserved'
   do i=2, nr
      if(r_old(i).le.r_old(i-1)) then
         write(*,*) 'error in recalc_vgridxyz: new radial y grid is not monotonic'
         write(*,*) '   either: trust original grid, and comment this procedure in main program'
         write(*,*) '   or find bug in algorithm (e.g. for non-monotonic velocity fields?'
         stop
      endif
   enddo
!
!if everything has worked: reset minimum and maximum radius to exact value, and reset x-coordinates
   r_old(1) = r_orig(1)
   r_old(nr) = r_orig(nr)
   y_new=y
   do i=2, nr+1
      y_new(i) = -r_old(nr-i+2)
   enddo
   j=1
   do i=ndymax-nr, ndymax-1
      y_new(i) = r_old(j)
      j=j+1
   enddo
!
endif
!
deallocate(r_orig)
deallocate(r_old)
deallocate(r_new1)
deallocate(r_new2)
deallocate(r_new3)
deallocate(vel_orig)
deallocate(vel_old)
deallocate(vel_new)
deallocate(acoeff_vel)
deallocate(bcoeff_vel)
deallocate(ccoeff_vel)
deallocate(dcoeff_vel)
!
!y-grid done
!
!-------------------------------z-grid----------------------------------
!
ir=0
do i=1, ndzmax
!   if(abs(y(i)-rmin).lt.1.d-14) ir = i
   if(z(i).ge.rmin) then
      ir = i
      exit
   endif
enddo
!
if(ir.eq.0) stop 'error in recalc_vgridxyz: z(i)=rmin not found'
!
nr=ndzmax-ir

i=ndxmax/2+1
j=ndymax/2+1
kk=ndzmax-nr
do j=2, nr-1
   dv1 = velz3d(i,j,kk)-velz3d(i,j,kk-1)
   dv2 = velz3d(i,j,kk+1)-velz3d(i,j,kk)
   if(dv1*dv2.lt.zero) then
      write(*,*) 'recalc_vgridxyz not performed for non-monotonic velocity field in z'
      return
   endif
   kk=kk+1
enddo
!
allocate(r_orig(nr), stat=err)
   if(err.ne.0) stop 'allocation error in recalc_vgridxyz'
allocate(r_old(nr), stat=err)
   if(err.ne.0) stop 'allocation error in recalc_vgridxyz'
allocate(r_new1(nr), stat=err)
   if(err.ne.0) stop 'allocation error in recalc_vgridxyz'
allocate(r_new2(nr), stat=err)
   if(err.ne.0) stop 'allocation error in recalc_vgridxyz'
allocate(r_new3(nr), stat=err)
   if(err.ne.0) stop 'allocation error in recalc_vgridxyz'
allocate(vel_orig(nr), stat=err)
   if(err.ne.0) stop 'allocation error in recalc_vgridxyz'
allocate(vel_old(nr), stat=err)
   if(err.ne.0) stop 'allocation error in recalc_vgridxyz'
allocate(vel_new(nr), stat=err)
   if(err.ne.0) stop 'allocation error in recalc_vgridxyz'
allocate(acoeff_vel(nr), stat=err)
   if(err.ne.0) stop 'allocation error in recalc_vgridxyz'
allocate(bcoeff_vel(nr), stat=err)
   if(err.ne.0) stop 'allocation error in recalc_vgridxyz'
allocate(ccoeff_vel(nr), stat=err)
   if(err.ne.0) stop 'allocation error in recalc_vgridxyz'
allocate(dcoeff_vel(nr), stat=err)
   if(err.ne.0) stop 'allocation error in recalc_vgridxyz'
allocate(z_new(ndzmax), stat=err)
   if(err.ne.0) stop 'allocation error in recalc_vgridxyz'
!
!set radial stratification
j=1
do i=ndzmax-nr, ndymax-1
   r_old(j) = y(i)
   vel_old(j) = velz3d(ndxmax/2+1,ndymax/2+1,i)
   j=j+1
enddo
!
r_orig=r_old
vel_orig=vel_old
!
!spline-interpolation coefficients
call spline_coeff(nr, r_orig, vel_orig, acoeff_vel, bcoeff_vel, ccoeff_vel, dcoeff_vel)
!
!check if radial stratification has too large velocity-steps
lx = .false.
do i=2, nr
   delv_x = abs(vel_old(i)-vel_old(i-1))
   if(delv_x.gt.delv) then
      lx=.true.
      exit
   endif
enddo
!
if(.not.lx) then
!input x-stratification is okay
   z_new = z
else
!calculate new x-stratification
!
   delv_x = (vel_old(nr)-vel_old(1))/(nr-1)
   delv_x = max(delv,delv_x)
!
   iterationz: do iit=1, itmax
      r_new1 = r_old
      vel_new = vel_old
      lgood = .true.
      do i=2, nr
         delv_old = vel_old(i)-vel_old(i-1)
         if(abs(delv_old)/delv_x.gt.1.01) then
!calculate new coordinates at (i-1) and at i
            delx_old = r_old(i) - r_old(i-1)
            delx_new = delv_x*delx_old/delv_old
            if(i.eq.2) then
               r_new1(i) = r_old(i) - delx_new
            elseif(i.eq.nr) then
               r_new1(i-1) = r_old(i-1) + delx_new
            else
               r_new1(i-1) = r_old(i-1) + (delx_old-delx_new)/two
               r_new1(i) = r_old(i) - (delx_old-delx_new)/two
            endif
!
!scale up rest of the grid towards left by keeping delx-ratios from original grid
            r_new2=r_new1
            do j=i-2, 1, -1
               r_new2(j) = r_new2(j+1) - abs((r_new2(j+2)-r_new2(j+1)))*(r_old(j+1)-r_old(j))/(r_old(j+2)-r_old(j+1))
            enddo
!scale up rest of the grid towards right by keeping delx-ratios from original grid
            do j=i+1, nr
               r_new2(j) = r_new2(j-1) + abs((r_new2(j-1)-r_new2(j-2)))*(r_old(j)-r_old(j-1))/(r_old(j-1)-r_old(j-2))
            enddo
!
!scale the grid to the left to ensure that left boundary is correct
            r_new3=r_new2
            cc = (r_new2(i-1)-r_old(1))/(r_new2(i-1)-r_new2(1))
            do j=i-2, 1, -1
            r_new3(j) = r_new3(j+1)-cc*(r_new2(j+1)-r_new2(j))
            enddo
!scale the grid to the left to ensure that left boundary is correct
            cc = (r_old(nr)-r_new2(i))/(r_new2(nr)-r_new2(i))
            do j=i+1, nr
               r_new3(j) = r_new3(j-1)+cc*(r_new2(j)-r_new2(j-1))
            enddo
!
!finally: interpolate velocities
            do j=1, nr
               call find_index(r_new3(j), r_orig, nr, iim2, iim1, ii, iip1)
!               vel_new(j) = vel_old(iim1) + (r_new3(j)-r_old(iim1))*(vel_old(ii)-vel_old(iim1))/(r_old(ii)-r_old(iim1))
                vel_new(j)=interpol_yp_spline(acoeff_vel(ii), bcoeff_vel(ii), ccoeff_vel(ii), dcoeff_vel(ii), r_orig(ii), r_new3(j))
            enddo
!
!            write(*,*) i
!            do j=1, nr
!               write(*,'(10es20.8)') r_orig(j), r_old(j), r_new1(j), r_new2(j), r_new3(j), vel_orig(j), vel_new(j)
!            enddo
!            write(*,*)
!
!overwrite original array and begin from start
            r_old=r_new3
            vel_old=vel_new
            lgood=.false.
!            if(abs(r_old(1)-one).gt.1.d-8) stop
            exit
         endif
!grid is okay and stop iteration
      enddo
!      write(*,*) iit, lgood
      if(lgood) exit iterationz
   enddo iterationz
!
!check if everything has worked
   if(.not.lgood) then
      write(*,*) 'error in recalc_vgridxyz: no solution found for nice z velocity-grid properties'
      write(*,*) '   either: increase number of iterations, or error in algorithm (e.g. for non-monotonic velocity fields?)'
      stop
   endif
   if(abs(r_old(1)-r_orig(1)).gt.1.d-10) stop 'error in recalc_vgridxyz: z minimum radius not conserved'
   if(abs(r_old(nr)-r_orig(nr)).gt.1.d-10) stop 'error in recalc_vgridxyz: z maximum radius not conserved'
   do i=2, nr
      if(r_old(i).le.r_old(i-1)) then
         write(*,*) 'error in recalc_vgridxyz: new radial z grid is not monotonic'
         write(*,*) '   either: trust original grid, and comment this procedure in main program'
         write(*,*) '   or find bug in algorithm (e.g. for non-monotonic velocity fields?'
         stop
      endif
   enddo
!
!if everything has worked: reset minimum and maximum radius to exact value, and reset x-coordinates
   r_old(1) = r_orig(1)
   r_old(nr) =r_orig(nr)
   z_new=z
   do i=2, nr+1
      z_new(i) = -r_old(nr-i+2)
   enddo
   j=1
   do i=ndzmax-nr, ndzmax-1
      z_new(i) = r_old(j)
      j=j+1
   enddo
!
endif
!
deallocate(r_orig)
deallocate(r_old)
deallocate(r_new1)
deallocate(r_new2)
deallocate(r_new3)
deallocate(vel_orig)
deallocate(vel_old)
deallocate(vel_new)
deallocate(acoeff_vel)
deallocate(bcoeff_vel)
deallocate(ccoeff_vel)
deallocate(dcoeff_vel)
!
!z-grid done
!
!-----------------------------------------------------------------------
!
!overwrite original x-variables
x=x_new
y=y_new
z=z_new
!
!
!
end subroutine recalc_vgridxyz
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine recalc_vgridxyz2
!
!------------set equidistant velocity steps on coordinate axes----------
!----------------------from stellar surface to rmax,--------------------
!----------------given a previous velocity stratification---------------
!
use prog_type
use fund_const
use mod_directories, only: output_dir
use dime3d, only: ndx, ndxmax, ndymax, ndzmax, x, y, z, velx3d, vely3d, velz3d
use dime1d, only: delv, n1d_dum, r1d_dum
use inf_reg, only: rmin
use params_input, only: vth_fiducial, vmax, vmin, beta
use mod_interp1d, only: find_index, cube_mono_coeff, interpol_yp_spline
use mod_grid, only: calc_pdf
!
implicit none
!
! ... local scalars
integer(i4b), parameter :: itmax=100000
integer(i4b) :: i, j, k, iit, err, nx_final
integer(i4b) :: iim2, iim1, ii, iip1, ir, nr
real(dp) :: delv_x, cc, delx_new, delx_old, delv_old, delv_new, b, xmin, xmax
!
! ... local logicals
logical :: lx, lgood
!
! ... local arrays
real(dp), dimension(:), allocatable :: r_orig, r_old, r_new1, r_new2, r_new3, &
                                       vel_orig, vel_old, vel_new, &
                                       acoeff_vel, bcoeff_vel, ccoeff_vel, dcoeff_vel, &
                                       x_new, y_new, z_new
real(dp), dimension(:), allocatable :: xfinal, xfinal_mid, pdf_xfinal, p_xfinal
!
! ... local functions
real(dp) :: bvel
!
!-------------------------calculate radial grid-------------------------
!
call gridr
!
allocate(r_orig(n1d_dum), stat=err)
   if(err.ne.0) stop 'allocation error in recalc_vgridxyz2'
allocate(vel_orig(n1d_dum), stat=err)
   if(err.ne.0) stop 'allocation error in recalc_vgridxyz2'
!
allocate(acoeff_vel(n1d_dum), stat=err)
   if(err.ne.0) stop 'allocation error in recalc_vgridxyz2'
allocate(bcoeff_vel(n1d_dum), stat=err)
   if(err.ne.0) stop 'allocation error in recalc_vgridxyz2'
allocate(ccoeff_vel(n1d_dum), stat=err)
   if(err.ne.0) stop 'allocation error in recalc_vgridxyz2'
allocate(dcoeff_vel(n1d_dum), stat=err)
   if(err.ne.0) stop 'allocation error in recalc_vgridxyz2'
allocate(x_new(ndxmax), stat=err)
   if(err.ne.0) stop 'allocation error in recalc_vgridxyz2'
!
r_orig=r1d_dum
!
!calculate corresponding velocities
b=one-(vmin/vmax)**(one/beta)
do i=1, n1d_dum
   vel_orig(i) = bvel(r_orig(i), vmax*1.d5, b, beta)/vth_fiducial
enddo
!
!spline-interpolation coefficients
!call spline_coeff(n1d_dum, r_orig, vel_orig, acoeff_vel, bcoeff_vel, ccoeff_vel, dcoeff_vel)
call cube_mono_coeff(n1d_dum, r_orig, vel_orig, acoeff_vel, bcoeff_vel, ccoeff_vel, dcoeff_vel)
!call lin_coeff(n1d_dum, r_orig, vel_orig, ccoeff_vel, dcoeff_vel)
!acoeff_vel=zero
!bcoeff_vel=zero
!
!
!-------------------------------x-grid----------------------------------
!
ir=0
do i=1, ndxmax
!   if(abs(x(i)-rmin).lt.1.d-14) ir = i
   if(x(i).ge.rmin) then
      ir=i
      exit
   endif
enddo
!
if(ir.eq.0) stop 'error in recalc_vgridxyz: x(i)=rmin not found'
!
nr=ndxmax-ir
!
allocate(r_old(nr), stat=err)
   if(err.ne.0) stop 'allocation error in recalc_vgridxyz2'
allocate(r_new1(nr), stat=err)
   if(err.ne.0) stop 'allocation error in recalc_vgridxyz2'
allocate(r_new2(nr), stat=err)
   if(err.ne.0) stop 'allocation error in recalc_vgridxyz2'
allocate(r_new3(nr), stat=err)
   if(err.ne.0) stop 'allocation error in recalc_vgridxyz2'
allocate(vel_old(nr), stat=err)
   if(err.ne.0) stop 'allocation error in recalc_vgridxyz2'
allocate(vel_new(nr), stat=err)
   if(err.ne.0) stop 'allocation error in recalc_vgridxyz2'
!
!set radial stratification
j=1
do i=ndxmax-nr, ndxmax-1
   r_old(j) = x(i)
   vel_old(j) = velx3d(i,ndymax/2+1,ndzmax/2+1)
   j=j+1
enddo
!
xmin=r_old(1)
xmax=r_old(nr)
!
!check if radial stratification has too large velocity-steps
lx = .false.
do i=2, nr
   delv_x = abs(vel_old(i)-vel_old(i-1))
   if(delv_x.gt.delv) then
      lx=.true.
      exit
   endif
enddo
!
if(.not.lx) then
!input x-stratification is okay
   x_new = x
else
!calculate new x-stratification
!
   delv_x = abs(vel_old(nr)-vel_old(1))/(nr-2) !(nr-1)
   delv_x = max(delv,delv_x)
!
   iterationx: do iit=1, itmax
      r_new1 = r_old
      vel_new = vel_old
      lgood = .true.
      do i=2, nr
         delv_old = vel_old(i)-vel_old(i-1)
         if(abs(delv_old)/delv_x.gt.1.01d0) then
!calculate new coordinates at (i-1) and at i
            delx_old = r_old(i) - r_old(i-1)
            delx_new = delv_x*delx_old/delv_old
            if(i.eq.2) then
               r_new1(i) = r_old(i) - delx_new
            elseif(i.eq.nr) then
               r_new1(i-1) = r_old(i-1) + delx_new
            else
               r_new1(i-1) = r_old(i-1) + (delx_old-delx_new)/two
               r_new1(i) = r_old(i) - (delx_old-delx_new)/two
            endif
!
!scale up rest of the grid towards left by keeping delx-ratios from original grid
            r_new2=r_new1
            do j=i-2, 1, -1
               r_new2(j) = r_new2(j+1) - abs((r_new2(j+2)-r_new2(j+1)))*(r_old(j+1)-r_old(j))/(r_old(j+2)-r_old(j+1))
            enddo
!scale up rest of the grid towards right by keeping delx-ratios from original grid
            do j=i+1, nr
               r_new2(j) = r_new2(j-1) + abs((r_new2(j-1)-r_new2(j-2)))*(r_old(j)-r_old(j-1))/(r_old(j-1)-r_old(j-2))
            enddo
!
!scale the grid to the left to ensure that left boundary is correct
            r_new3=r_new2
            do j=i-2, 1, -1
               cc = (r_new2(i-1)-r_old(1))/(r_new2(i-1)-r_new2(1))
               r_new3(j) = r_new3(j+1)-cc*(r_new2(j+1)-r_new2(j))
            enddo
!scale the grid to the left to ensure that left boundary is correct
            cc = (r_old(nr)-r_new2(i))/(r_new2(nr)-r_new2(i))
            do j=i+1, nr
               r_new3(j) = r_new3(j-1)+cc*(r_new2(j)-r_new2(j-1))
            enddo
!
!finally: interpolate velocities
            do j=1, nr
               call find_index(r_new3(j), r_orig, n1d_dum, iim2, iim1, ii, iip1)
!               vel_new(j) = vel_old(iim1) + (r_new3(j)-r_old(iim1))*(vel_old(ii)-vel_old(iim1))/(r_old(ii)-r_old(iim1))
                vel_new(j)=interpol_yp_spline(acoeff_vel(ii), bcoeff_vel(ii), ccoeff_vel(ii), dcoeff_vel(ii), r_orig(ii), r_new3(j))
            enddo
!
!overwrite original array and begin from start
            r_old=r_new3
            vel_old=vel_new
            lgood=.false.
!            if(abs(r_old(1)-one).gt.1.d-8) stop
!            write(*,*) delv_x, abs(delv_old), abs(delv_old)/delv_x
            exit
         endif
!grid is okay and stop iteration
      enddo
!      write(*,*) iit, lgood, delv_x, delv_old
!      write(*,*) r_new3
      if(lgood) exit iterationx
   enddo iterationx
!
!check if everything has worked
   if(.not.lgood) then
      write(*,*) 'error in recalc_vgridxyz2: no solution found for nice x velocity-grid properties'
      write(*,*) '   either: increase number of iterations, or error in algorithm (e.g. for non-monotonic velocity fields?)'
      stop
   endif
   if(abs(r_old(1)-xmin).gt.1.d-10) stop 'error in recalc_vgridxyz2: x minimum radius not conserved'
   if(abs(r_old(nr)-xmax).gt.1.d-10) stop 'error in recalc_vgridxyz2: x maximum radius not conserved'
   do i=2, nr
      if(r_old(i).le.r_old(i-1)) then
         write(*,*) 'error in recalc_vgridxyz2: new radial x grid is not monotonic'
         write(*,*) '   either: trust original grid, and comment this procedure in main program'
         write(*,*) '   or find bug in algorithm (e.g. for non-monotonic velocity fields?'
         stop
      endif
   enddo
!
!if everything has worked: reset minimum and maximum radius to exact value, and reset x-coordinates
   r_old(1) = xmin  !r_orig(1)
   r_old(nr) = xmax !r_orig(n1d_dum)
   x_new=x
   do i=2, nr+1
      x_new(i) = -r_old(nr-i+2)
   enddo
   j=1
   do i=ndxmax-nr, ndxmax-1
      x_new(i) = r_old(j)
      j=j+1
   enddo
!
!
endif
!
!
deallocate(r_orig)
deallocate(r_old)
deallocate(r_new1)
deallocate(r_new2)
deallocate(r_new3)
deallocate(vel_orig)
deallocate(vel_old)
deallocate(vel_new)
deallocate(acoeff_vel)
deallocate(bcoeff_vel)
deallocate(ccoeff_vel)
deallocate(dcoeff_vel)
!
!x-coord done
!
!overwrite original x,y,z grid
x=x_new
y=x_new
z=x_new
!
!------------------------calculate corresponding pdf--------------------
!
nx_final = ndx-1
!
allocate(xfinal(nx_final), stat=err)
allocate(xfinal_mid(nx_final-1), stat=err)
allocate(pdf_xfinal(nx_final-1), stat=err)
allocate(p_xfinal(nx_final-1), stat=err)
!
j=1
do i=ndxmax/2+1, ndxmax-1
   xfinal(j) = x(i)
   j=j+1
enddo
!
call calc_pdf(nx_final, xfinal, xfinal_mid, pdf_xfinal, p_xfinal)
!
open(1, file=output_dir//'/pdf_xgrid3.dat')
   write(1,*) 'actually used distribution at third step: recalc_vgridxyz2'
   write(1,'(3a20)') 'x_mid', 'pdf(x)', 'p(x)'
   do i=1, nx_final-1 
      write(1,'(3es20.8)') xfinal_mid(i), pdf_xfinal(i), p_xfinal(i)
   enddo
close(1)
!
!
!
end subroutine recalc_vgridxyz2
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine set_xyz_phantom
!
!----calculate phantom point on outer boundary from thomson-opacity-----
!---------with the ansatz that delta_tau shall be constant--------------
!--------store the phantom-point coordinates on x,y,z-grid--------------
!
use prog_type
use fund_const
use dime3d, only: ndxmax, ndymax, ndzmax, x, y, z, opac3d
!
implicit none
!
! ... local scalars
integer(i4b) :: i
real(dp) :: grad, x_phantom, y_phantom, z_phantom, dum
!
!-------------------phantom point along x-direction---------------------
!
!calculate phantom point at outer left boundary
if(opac3d(2,ndymax/2+1,ndzmax/2+1).le.zero) then
   write(*,*) 'error in set_xyz_phantom: opac3d(2,ndymax/2+1,ndzmax/2+1) <= 0.'
   write(*,*) '-> increase kcont'   
   write(*,*) '-> increase grid resolution'
   write(*,*) '-> increase rmax'   
   stop
endif

if(opac3d(3,ndymax/2+1,ndzmax/2+1).le.zero) then
   write(*,*) 'error in set_xyz_phantom: opac3d(2,ndymax/2+1,ndzmax/2+1) <= 0.'
   write(*,*) '-> increase kcont'      
   write(*,*) '-> increase grid resolution'
   write(*,*) '-> increase rmax'      
   stop
endif

grad = (log10(opac3d(2,ndymax/2+1,ndzmax/2+1)) - log10(opac3d(3,ndymax/2+1,ndzmax/2+1))) / &
        (log10(x(2)/x(3)))
grad = grad + one
dum = two*(-one*x(2))**grad - (-one*x(3))**grad
!
if(abs(grad).gt.one.and.dum.lt.zero) then
   !make equidistant step
   x(1)=x(2)-(x(3)-x(2))
else
   x_phantom = -one * (two*(-one*x(2))**grad - (-one*x(3))**grad)**(one/grad)
   x(1) = x_phantom
   if(x_phantom.ne.x_phantom) then
   !make equidistant grid step outside if phantom point can not be calculated (x(2)>>x(3))
      x(1)=x(2)-(x(3)-x(2))
   else
      x(1)=x_phantom
   endif
endif
!
!calculate phantom point at outer right boundary
grad = (log10(opac3d(ndxmax-1,ndymax/2+1,ndzmax/2+1))-log10(opac3d(ndxmax-2,ndymax/2+1,ndzmax/2+1))) / &
       (log10(x(ndxmax-1) / x(ndxmax-2)))
grad = grad + one
dum = two*(x(ndxmax-1))**grad - (x(ndxmax-2))**grad
if(abs(grad).gt.one.and.dum.lt.zero) then
   !make equidistant step
   x(ndxmax) = x(ndxmax-1) + (x(ndxmax-1)-x(ndxmax-2))
else
   x_phantom = (two*(x(ndxmax-1))**grad - (x(ndxmax-2))**grad)**(one/grad)
   if(x_phantom.ne.x_phantom) then
!make equidistant grid step outside if phantom point can not be calculated (x(2)>>x(3))
      x(ndxmax) = x(ndxmax-1) + (x(ndxmax-1)-x(ndxmax-2))
   else
      x(ndxmax)=x_phantom
   endif
endif
!
!-------------------phantom point along y-direction---------------------
!
!calculate phantom point at outer left boundary
grad = (log10(opac3d(ndxmax/2+1,2,ndzmax/2+1)) - log10(opac3d(ndxmax/2+1,3,ndzmax/2+1))) / &
        (log10(y(2)/y(3)))
grad = grad + one
y_phantom = -one * (two*(-one*y(2))**grad - (-one*y(3))**grad)**(one/grad)
!
if(y_phantom.ne.y_phantom) then
!make equidistant grid step outside if phantom point can not be calculated (y(2)>>y(3))
   y(1)=y(2)-(y(3)-y(2))
else
   y(1)=y_phantom
endif
!
!calculate phantom point at outer right boundary
grad = (log10(opac3d(ndxmax/2+1,ndymax-1,ndzmax/2+1))-log10(opac3d(ndxmax/2+1,ndymax-2,ndzmax/2+1))) / &
       (log10(y(ndymax-1) / y(ndymax-2)))
grad = grad + one
y_phantom = (two*(y(ndymax-1))**grad - (y(ndymax-2))**grad)**(one/grad)
!
if(y_phantom.ne.y_phantom) then
!make equidistant grid step outside if phantom point can not be calculated (y(ndymax-2)>>y(ndymax-3))
   y(ndymax) = y(ndymax-1) + (y(ndymax-1)-y(ndymax-2))
else
   y(ndymax)=y_phantom
endif
!
!-------------------phantom point along z-direction---------------------
!
!calculate phantom point at outer left boundary
grad = (log10(opac3d(ndxmax/2+1,ndymax/2+1,2)) - log10(opac3d(ndxmax/2+1,ndymax/2+1,3))) / &
        (log10(z(2)/z(3)))
grad = grad + one
z_phantom = -one * (two*(-one*z(2))**grad - (-one*z(3))**grad)**(one/grad)
!
if(z_phantom.ne.z_phantom) then
!make equidistant grid step outside if phantom point can not be calculated (z(2)>>z(3))
   z(1)=z(2)-(z(3)-z(2))
else
   z(1)=z_phantom
endif
!
!calculate phantom point at outer right boundary
grad = (log10(opac3d(ndxmax/2+1,ndymax/2+1,ndzmax-1))-log10(opac3d(ndxmax/2+1,ndymax/2+1,ndzmax-2))) / &
       (log10(z(ndzmax-1) / z(ndzmax-2)))
grad = grad + one
z_phantom = (two*(z(ndzmax-1))**grad - (z(ndzmax-2))**grad)**(one/grad)
!
if(z_phantom.ne.z_phantom) then
!make equidistant grid step outside if phantom point can not be calculated (z(ndzmax-2)>>z(ndzmax-3))
   z(ndzmax) = z(ndzmax-1) + (z(ndzmax-1)-z(ndzmax-2))
else
   z(ndzmax)=z_phantom
endif
!
!
!
end subroutine set_xyz_phantom
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine check_gridxyz
!
use prog_type
use fund_const
use params_input, only: rmax
use inf_reg, only: rmin
use dime3d, only: ndxmax, ndymax, ndzmax, x, y, z
!
implicit none
!
! ... local scalars
integer(i4b) :: i
real(dp) :: delta
logical :: lzer, lmax, lmin
!
!-----check for symmetry and if x-coordinates 0., rmax are included-----
!
!z=(/-12.7481211305995d0, -12.0000000000000d0, -11.3385746448511d0, -10.6771492897023d0, &
!  -10.0157239345534d0, -9.35429857940452d0, -8.69287322425565d0, -8.03144786910678d0, &
!  -7.37002251395791d0, -6.70859715880904d0, -6.04717180366017d0, -5.38574644851130d0, &
!  -4.72485661367735d0, -4.11932410576098d0, -3.58660457924528d0, -3.13555974600510d0, &
!  -2.75424862726760d0, -2.43475717569696d0, -2.16867407954251d0, -1.94337103598786d0, &
!  -1.75353812979766d0, -1.59120538618391d0, -1.45020844447888d0, -1.32941714503406d0, &
!  -1.22533542489622d0, -1.13720852325380d0, -1.07008508311000d0, -1.02573928057794d0, &
!  -1.00004941366441d0, -1.00000000000000d0, -0.962659950616609d0, -0.920644121675200d0, &
!  -0.874894379985053d0, -0.825789824258372d0, -0.773688422240171d0, -0.718386118816189d0, &
!  -0.660590848407663d0, -0.600119920277130d0, -0.537190310755168d0, -0.472108474400704d0, &
!  -0.408091686433751d0, -0.343469350308003d0, -0.276804144494936d0, -0.208209831161752d0, &
!  -0.136963645091191d0, -6.273065308998678d-2, zero, 6.273065308998678d-2, 0.136963645091191d0, &
!  0.208209831161752d0, 0.276804144494936d0, 0.343469350308003d0, 0.408091686433751d0, &
!  0.472108474400704d0, 0.537190310755168d0, 0.600119920277130d0, 0.660590848407663d0, &
!  0.718386118816189d0, 0.773688422240171d0, 0.825789824258372d0, 0.874894379985053d0, &
!  0.920644121675200d0, 0.962659950616609d0, 1.00000000000000d0, 1.00004941366441d0, &
!  1.02573928057794d0, 1.07008508311000d0, 1.13720852325380d0, 1.22533542489622d0, &
!  1.32941714503406d0, 1.45020844447888d0, 1.59120538618391d0,1.75353812979766d0, &
!  1.94337103598786d0, 2.16867407954251d0, 2.43475717569696d0, 2.75424862726760d0, &
!  3.13555974600510d0, 3.58660457924528d0, 4.11932410576098d0, 4.72485661367735d0, &
!  5.38574644851130d0, 6.04717180366017d0, 6.70859715880904d0,7.37002251395791d0, &
!  8.03144786910678d0, 8.69287322425565d0, 9.35429857940452d0, 10.0157239345534d0, &
!  10.6771492897023d0, 11.3385746448511d0, 12.0000000000000d0, 12.7481211305995d0 /)
!x=z
!y=z

!open(1, file='outputFILES_DEBUG/xtrash.dat')
!   do i=1, ndxmax
!      read(1,*) x(i)
!!      write(1,*) x(i)
!   enddo
!close(1)
!y=x
!z=x
!!
!do i=1, ndxmax
!   write(*,*) x(i)
!enddo
!stop 'go on in check_gridxyz'
!
lzer=.false.
lmax=.false.
lmin=.false.
!
do i=1, ndxmax
   if(x(i).eq.zero) lzer=.true.
   if(abs(abs(x(i))-rmax).lt.1.d-13) lmax=.true.
   if(abs(x(i)-rmin).lt.1.d-13) lmin=.true.
   if(abs(x(i)+x(ndxmax+1-i)).gt.1.d-13) stop 'error in check_gridxyz: x-grid not symmetric'
enddo
!
if(.not.lzer) stop 'error in check_gridxyz: zero coordinate of x-grid not included'
if(.not.lmax) stop 'error in check_gridxyz: rmax coordinate of x-grid not included'
if(.not.lmin) then
   write(*,*) 'warning in check_gridxyz: rmin coordinate of x-grid not included'
   write(*,*) '   -> play with input n1d, ncx, delx_max'
   write(*,*) '   -> or ignore this warning'
!   stop
endif
!
do i=2, ndxmax
   delta = x(i)-x(i-1)
   if(delta.lt.zero) stop 'error in check_gridxyz: x-grid not monotonic'
enddo
!
!------------------------------same for y-grid--------------------------
!
lzer=.false.
lmax=.false.
lmin=.false.
!
do i=1, ndymax
   if(y(i).eq.zero) lzer=.true.
   if(abs(abs(y(i))-rmax).lt.1.d-13) lmax=.true.
   if(abs(y(i)-rmin).lt.1.d-13) lmin=.true.
   if(abs(y(i)+y(ndymax+1-i)).gt.1.d-13) stop 'error in check_gridxyz: y-grid not symmetric'
enddo
!
if(.not.lzer) stop 'error in check_gridxyz: zero coordinate of y-grid not included'
if(.not.lmax) stop 'error in check_gridxyz: rmax coordinate of y-grid not included'
if(.not.lmin) then
   write(*,*) 'warning in check_gridxyz: rmin coordinate of y-grid not included'
   write(*,*) '   -> play with input n1d, ncy, dely_max'
   write(*,*) '   -> or ignore this warning'
!   stop
endif

do i=2, ndymax
   delta = y(i)-y(i-1)
   if(delta.lt.zero) stop 'error in check_gridxyz: y-grid not monotonic'
enddo
!
!------------------------------same for z-grid--------------------------
!
lzer=.false.
lmax=.false.
lmin=.false.
!
do i=1, ndzmax
   if(z(i).eq.zero) lzer=.true.
   if(abs(abs(z(i))-rmax).lt.1.d-13) lmax=.true.
   if(abs(z(i)-rmin).lt.1.d-13) lmin=.true.
   if(abs(z(i)+z(ndzmax+1-i)).gt.1.d-13) stop 'error in check_gridxyz: z-grid not symmetric'
enddo
!
if(.not.lzer) stop 'error in check_gridxyz: zero coordinate of z-grid not included'
if(.not.lmax) stop 'error in check_gridxyz: rmax coordinate of z-grid not included'
if(.not.lmin) then
   write(*,*) 'warningr in check_gridxyz: rmin coordinate of z-grid not included'
   write(*,*) '   -> play with input n1d, ncz, delz_max'
   write(*,*) '   -> or ignore this warning'
!   stop
endif

do i=2, ndzmax
   delta = z(i)-z(i-1)
   if(delta.lt.zero) stop 'error in check_gridxyz: z-grid not monotonic'
enddo
!
!
!
end subroutine check_gridxyz
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine gridr
!
use prog_type
use options, only: spatial_grid1d
!
implicit none
!
select case(spatial_grid1d)
   case(0)
      call grid1d_r_equi
   case(1)
      call grid1d_vel_equi
   case(2)
      call grid1d_tau_equi
   case(3)
      call grid1d_tau_log
   case(4)
      call grid1d_final
   case(5)
      call grid1d_final_2
   case(6)
      call grid1d_r_log
   case default
      stop 'error in gridr: unvalid spatial_grid1d'
end select
 
end subroutine gridr
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine grid1d_vel_equi
!
!-----------------------------------------------------------------------
!-------------------sets up 1-dimensional radial grid:------------------
!---------------equidistant in velocity space in steps of delv----------
!-----------------------------------------------------------------------
!
use prog_type
use fund_const
use dime1d, only: n1d, n1d_dum, n1d_t, n1d_r, r1d, r1d_dum, delv
use params_input, only: teff, tmin, vmax, vmin, beta, xmloss, yhe, hei, rmax, &
                        vth_fiducial
use inf_reg, only: rmin
use mod_sort, only: bubblesort
use mod_grid, only: recalc_grid1d
!
implicit none
!
! ... local scalars
integer(i4b) :: i
integer(i4b) :: err, indx
real(dp) :: b
real(dp) :: xmlosscgs
real(dp) :: c1, c2
real(dp) :: taumin, rmin_dum
real(dp) :: delvcgs
real(dp) :: vrmax
real(dp) :: max_delr_old, max_deltau_old, max_delr, max_deltau
real(dp) :: val1, val2, delr, delr1, maxr
!
! ... local arrays
real(dp), dimension(:), allocatable :: vel1d_dum
!
! ... local functions
!
!-------------------------calculate array dimensions--------------------
!number of 1d-grid points is (v(rmax)-v(rmin))/delv(fiducial)
!
allocate(r1d(n1d), stat=err)
   if(err.ne.0) stop 'error in grid1d_vel_equi: allocation'
!
write(*,*) 'v_min: ', vmin*1.d5
write(*,*) 'v_max: ', vmax*1.d5
write(*,*) 'v_thermal (fiducial) :', vth_fiducial
write(*,*) 'setting up velocity grid in units of', delv, 'v_thermal (fiducial)'
!
delvcgs=delv*vth_fiducial
!
!calculate b-factor for beta-velocity-law
b=one-(vmin/vmax)**(one/beta)
!calculate radius (vmin)
rmin_dum=b/(one-(vmin/vmax)**(one/beta))
!calculate velocity at rmax (beta-velocity-law)
vrmax=vmax*(one - b* rmin_dum/rmax)**beta
!
!calculate number of 1-d data points
n1d_dum=(vrmax-vmin)*1.d5/delvcgs + 1
!
!ensure that dummy velocity grid is larger than input grid (otherwise: no averaging in recalc_grid1d)
do i=1, 10
   if(n1d_dum.lt.n1d) then
      delvcgs=0.75*delvcgs
      n1d_dum=(vrmax-vmin)*1.d5/delvcgs + 1
   else
      exit
   endif
enddo
!
!allocate dummy arrays
allocate(vel1d_dum(n1d_dum), stat=err)
   if(err.ne.0) stop 'allocation error in grid1d_vel_equi: vel1d_dum'
allocate(r1d_dum(n1d_dum), stat=err)
   if(err.ne.0) stop 'allocation error in grid1d_vel_equi: r1d_dum'
!
!-------------------------setting up velocity grid----------------------
!
vel1d_dum(1)=vmin*1.d5
!
do i=2, n1d_dum
   vel1d_dum(i)=vel1d_dum(i-1)+delvcgs
enddo
!
!-------------------------setting up radius grid------------------------
!
b=one-(vmin/vmax)**(one/beta)
!
do i=1, n1d_dum
   r1d_dum(i)=b/(one-(vel1d_dum(i)/(vmax*1.d5))**(one/beta))
enddo
!
!set outermost point to maximum radius
!
r1d_dum(n1d_dum)=rmax
!
!------------setting up actual used radial grid from dummy grid---------
!
call recalc_grid1d(r1d_dum, n1d_dum, n1d, 2, r1d)
!
r1d(n1d-1)=rmin
r1d(n1d)=rmax
!
!sort r1d-grid
call bubblesort(r1d,n1d)
!
!--------------------check if grid-values occur twice-------------------
!
do i=2, n1d-1
   val1=r1d(i-1)
   val2=r1d(i)
   if(abs(val2-val1).lt.1.d-14) then
!insert value on half range (to upper grid point)
      r1d(i)=(r1d(i+1)+r1d(i))/two
   endif
enddo
!
!same on outer boundary
val1=r1d(n1d)
val2=r1d(n1d-1)
if(abs(val2-val1).lt.1.d-14) then
!insert value on half range (to lower grid point)
   r1d(n1d-1)=(r1d(n1d-2)+r1d(n1d-1))/two
endif
!
!-----make outer part of grid equidistant in order that delta(r1d)------
!---------------------is not decreasing---------------------------------
!
do i=2, n1d
   delr1=r1d(i)-r1d(i-1)
   max_delr=delr1*(n1d-i)
   maxr=r1d(i)+max_delr
   if(maxr.lt.rmax) then
      indx=i
   endif
enddo

if(n1d.ne.indx) then
   delr=(rmax-r1d(indx))/(n1d-indx)
   do i=indx+1, n1d
      r1d(i)=r1d(i-1)+delr
   enddo
endif
!
!
end subroutine grid1d_vel_equi
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine grid1d_tau_log
!
!-----------------------------------------------------------------------
!-------------------sets up 1-dimensional radial grid:------------------
!-----equidistant in log(tau_thomson) (estimated from beta-vel-law)-----
!-----------------------------------------------------------------------
!
use prog_type
use fund_const
use dime1d, only: n1d, n1d_dum, n1d_t, n1d_r, r1d, r1d_dum
use params_input, only: teff, vmax, vmin, beta, xmloss, yhe, hei, rmax
use params_stellar, only: sr
use inf_reg, only: rmin
use mod_sort, only: bubblesort
use mod_grid, only: recalc_grid1d
!
implicit none
!
! ... local scalars
integer(i4b) :: i, j, k
integer(i4b) :: gp_used, gp_free            !used and free grid points
integer(i4b) :: ind_max, del_ind, indx, indx_1, indx_2, err
real(dp) :: b, velr
real(dp) :: xmlosscgs
real(dp) :: c1, c2
real(dp) :: taumin, taumax, deltau, deltau_log, delr
real(dp) :: max_delr, max_delr_old, max_deltau, max_deltau_old, delr1, delr2, maxr
real(dp) :: delvcgs
real(dp) :: sigem
real(dp) :: taur_0
real(dp) :: val1, val2
!
! ... local characters
!
n1d_dum = n1d_t + n1d_r
!
!------------------calculate/define constants needed--------------------
!
!calculate b-factor for beta-velocity-law
b=one-(vmin/vmax)**(one/beta)
!
c1=(one+four*yhe)*cgs_mp
c2=(one+hei*yhe)/c1
!
!calculate mass loss rate in cgs-units
xmlosscgs=xmloss*xmsu/(365.25d0*24.d0*60.d0*60.d0)
!
!calculate tau_min, set tau_max
taumin=(xmlosscgs*sigmae*c2/(four*pi)) / (vmax*1.d5 * rmax*sr)
taumax=20
!
!allocate radial dummy grid
allocate(r1d_dum(n1d_dum), stat=err)
   if(err.ne.0) stop 'allocation error grid1d_tau_log: r1d_dum'
!
!-----------------------------------------------------------------------
!
!calculate n1d_t radial grid points that are equidistant in tau
!
r1d_dum=zero
!
r1d_dum(1)=rmin
!
outer: do j=1, 1000
!outer cycle adapts the optical depth-steps deltau(controlled by taumax)
!            until maximum radial grid point is less than rmax
!
   deltau_log=(log10(taumax)-log10(taumin))/(n1d_t-1)
!
!tau at grid point 1 is taumax
   taur_0=taumax
!
   do i=1, n1d_t-1
!calculate velocity at given radius
      velr=vmax*1.d5*(one - b*rmin/r1d_dum(i))**(beta)
!calculate opacity at given radius
      sigem=sigmae*c2*xmlosscgs / (four*pi*r1d_dum(i)*r1d_dum(i)*velr*sr)
!calculate tau steps at given radius
      deltau=(ten**(-deltau_log) - one) * taur_0
!calculate tau at given radius (+1 grid point)
      taur_0=taur_0-deltau
!calculate radial grid step at given radius
      delr= -deltau/sigem
!calculate radius for next grid point
      r1d_dum(i+1)=r1d_dum(i)+delr
!go to outer cycle if the calculated radius exceed rmax
      if(r1d_dum(i+1).gt.rmax) then
         taumax=taumax*99.d0/100.d0
         exit
      endif
   enddo
!
!stop outer cycle if the maximum calculated radius is less than rmax
   if(maxval(r1d_dum).lt.rmax) then
!      write(*,*) j, 'max', maxval(r1d_dum)/sr
      exit outer
   endif
!
enddo outer
!
if(maxval(r1d_dum).gt.rmax) then
   write(*,*) 'error in calculating grid: adapt cycle'
   stop
endif
!
!-----------------------------------------------------------------------
!
!insert n1d_r grid points at the center of maximum distance 
!                   between neighbouring radial grid points
!
do i=2, n1d_t
!calculate delta_r between neighbouring grid points
   delr1=r1d_dum(i)-r1d_dum(i-1)
!calculate maximum radius if grid points would be equidistant in r
   max_delr=delr1*(n1d_dum-i)
   maxr=r1d_dum(i)+max_delr
!save index where maximum radius (if equidistant in r) is lt rmax
   if(maxr.lt.rmax) then
      indx=i
   endif
enddo
!
!calculate new delta_r which is adopted from r(indx)
delr=(rmax-r1d_dum(indx))/(n1d_dum-indx)
!
do i=indx+1, n1d_dum
   r1d_dum(i)=r1d_dum(i-1)+delr
enddo
!
!------------setting up actual used radial grid from dummy grid---------
!
allocate(r1d(n1d), stat=err)
   if(err.ne.0) stop 'error in grid1d_tau_log: allocation'
call recalc_grid1d(r1d_dum, n1d_dum, n1d, 2, r1d)
!
r1d(n1d-1)=rmin
r1d(n1d)=rmax
!
!sort r1d-grid
call bubblesort(r1d,n1d)
!
!--------------------check if grid-values occur twice-------------------
!
do i=2, n1d-1
   val1=r1d(i-1)
   val2=r1d(i)
   if(abs(val2-val1).lt.1.d-14) then
!insert value on half range (to upper grid point)
      r1d(i)=(r1d(i+1)+r1d(i))/two
   endif
enddo
!
!same on outer boundary
val1=r1d(n1d)
val2=r1d(n1d-1)
if(abs(val2-val1).lt.1.d-14) then
!insert value on half range (to lower grid point)
   r1d(n1d-1)=(r1d(n1d-2)+r1d(n1d-1))/two
endif
!
!-----make outer part of grid equidistant in order that delta(r1d)------
!---------------------is not decreasing---------------------------------
!
do i=2, n1d
   delr1=r1d(i)-r1d(i-1)
   max_delr=delr1*(n1d-i)
   maxr=r1d(i)+max_delr
   if(maxr.lt.rmax*sr) then
      indx=i
   endif
enddo

if(n1d.ne.indx) then
   delr=(rmax*sr-r1d(indx))/(n1d-indx)
   do i=indx+1, n1d
      r1d(i)=r1d(i-1)+delr
   enddo
endif
!
!
!
end subroutine grid1d_tau_log
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine grid1d_final
!
!-----------------------------------------------------------------------
!-------------------sets up 1-dimensional radial grid:------------------
!-----------------------------------------------------------------------
!   n1d_t grid points equidistant in log(tau_thomson)
!         (estimated from beta-vel-law)
!         where increments of delr are adapted if delr(equidistant velocity) lower
!
!   n1d_r grid points equidistant in radius
!-----------------------------------------------------------------------
!
use prog_type
use dime1d, only: n1d, n1d_dum, n1d_t, n1d_r, r1d, r1d_dum, delv
use params_input, only: teff, tmin, vmax, vmin, beta, xmloss, yhe, hei, rmax, &
                        vth_fiducial
use fund_const
use params_stellar, only: sr
use inf_reg, only: rmin
use mod_sort, only: bubblesort
use mod_grid, only: recalc_grid1d
!
implicit none
!
! ... local scalars
integer(i4b) :: i, j, k
integer(i4b) :: gp_used, gp_free            !used and free grid points
integer(i4b) :: ind_max, del_ind, indx, indx_1, indx_2, err
integer(i4b) :: n1d_test
real(dp) :: b, velr
real(dp) :: xmlosscgs
real(dp) :: c1, c2
real(dp) :: taumin, taumax, deltau, deltau_log, delr, delr_v, delr_t
real(dp) :: v_0, v_1, v_rmax
real(dp) :: max_delr, max_delr_old, max_deltau, max_deltau_old, delr0, delr1, delr2, maxr
real(dp) :: delvcgs
real(dp) :: sigem
real(dp) :: taur_0
real(dp) :: val1, val2
!
! ... local characters
character(len=50) :: enter
!
! ... local functions
real(dp) :: bvel
!
n1d_dum = n1d_t + n1d_r
!
!------------------calculate/define constants needed--------------------
!
!calculate b-factor for beta-velocity-law
b=one-(vmin/vmax)**(one/beta)
!
c1=(one+four*yhe)*cgs_mp
c2=(one+hei*yhe)/c1
!
!calculate mass loss rate in cgs-units
xmlosscgs=xmloss*xmsu/(365.25d0*24.d0*60.d0*60.d0)
!
!calculate tau_min, set tau_max
taumin=(xmlosscgs*sigmae*c2/(four*pi)) / (vmax*1.d5 * rmax*sr)
taumax=20
!
!
!allocate radial dummy grid
allocate(r1d_dum(n1d_dum), stat=err)
   if(err.ne.0) stop 'allocation error in grid1d_final: r1d_dum'
!
!---------------check if number of data points needed ------------------
!---for equidistant velocity grid is on the order of pre-defined n1d----
!
v_rmax=bvel(rmax, vmax*1.d5, b, beta)/1.d5
delvcgs=delv*vth_fiducial
!
n1d_test=(v_rmax-vmin)*1.d5/delvcgs + 1
!
if(n1d_test.gt.n1d_dum) then
   write(*,*) 'minimum n1d_dum needed (for equi-vel-grid):', n1d_test
   stop
endif
!
!-----------------------------------------------------------------------
!
!calculate n1d_t radial grid points that are equidistant in tau
!
r1d_dum=zero
!
r1d_dum(1)=rmin
!
delvcgs=delv*vth_fiducial
!
!
outer: do j=1, 1000
!outer cycle adapts the optical depth-steps deltau(controlled by taumax)
!            until maximum radial grid point is less than rmax
!
   deltau_log=(log10(taumax)-log10(taumin))/(n1d_t-1)
!
!tau at grid point 1 is taumax
   taur_0=taumax
!
   do i=1, n1d_t-1
!calculate velocity at given radius
      velr=vmax*1.d5*(one - b*rmin/r1d_dum(i))**(beta)
!calculate opacity at given radius
      sigem=sigmae*c2*xmlosscgs / (four*pi*r1d_dum(i)*r1d_dum(i)*velr*sr)
!calculate tau steps at given radius
      deltau=(ten**(-deltau_log) - one) * taur_0
!calculate tau at given radius (+1 grid point)
      taur_0=taur_0-deltau
!calculate radial grid step at given radius if equidistant tau is used
      delr_t= -deltau/sigem
!calculate radial grid step at given radius if equidistant velocity is used
      v_0=bvel(r1d_dum(i), vmax*1.d5, b, beta)
      delr_v=(delvcgs*r1d_dum(i)*r1d_dum(i))/(b*beta*v_0) * (one-b/r1d_dum(i))
!calculate radius for next grid point where increments of delr are taken from
!     the smaller values of delr(equidistant in tau), delr(equidistant in velocity)
      if(delr_v.lt.delr_t) then
         r1d_dum(i+1)=r1d_dum(i)+delr_v
!         write(*,fmt='(i4, 3(e20.8))') i, r1d_dum(i-1)/sr, r1d_dum(i)/sr, r1d_dum(i+1)/sr
      else
         r1d_dum(i+1)=r1d_dum(i)+delr_t
      endif
!go to outer cycle if the calculated radius exceed rmax
      if(r1d_dum(i+1).gt.rmax) then
         taumax=taumax*99.d0/100.d0
         exit
      endif
   enddo
!
!stop outer cycle if the maximum calculated radius is less than rmax
   if(maxval(r1d_dum).lt.rmax) then
      write(*,*) j, 'max', maxval(r1d_dum)/sr
      exit outer
   endif
!
enddo outer
!
if(maxval(r1d_dum).gt.rmax) then
   write(*,*) 'error in calculating grid: adapt cycle'
   stop
endif
!
if(maxval(r1d_dum).lt.0.8d0*rmax) then
   write(*,*) 'max(r1d_dum) too low: enlarge n1d_dum'
   stop
endif
!
gp_used=n1d_t
gp_free=n1d_dum-n1d_t
!
!
!-----------------------------------------------------------------------
!
!insert n1d_r grid points at the center of maximum distance 
!                   between neighbouring radial grid points

do i=2, n1d_t
!calculate delta_r between neighbouring grid points
   delr1=r1d_dum(i)-r1d_dum(i-1)
!calculate maximum radius if grid points would be equidistant in r
   max_delr=delr1*(n1d_dum-i)
   maxr=r1d_dum(i)+max_delr
!save index where maximum radius (if equidistant in r) is lt rmax
   if(maxr.lt.rmax) then
      indx=i
   endif
enddo
!
!
!calculate new delta_r which is adopted from r(indx)
if(n1d_dum-indx.eq.0) then
   delr=zero
else
   delr=(rmax-r1d_dum(indx))/(n1d_dum-indx)
endif
!
do i=indx+1, n1d_dum
   r1d_dum(i)=r1d_dum(i-1)+delr
enddo
!
!------------setting up actual used radial grid from dummy grid---------
!
allocate(r1d(n1d), stat=err)
   if(err.ne.0) stop 'error in grid1d_final: allocation'
!
call recalc_grid1d(r1d_dum, n1d_dum, n1d, 2, r1d)
!
r1d(n1d-1)=rmin
r1d(n1d)=rmax
!
!sort r1d-grid
call bubblesort(r1d,n1d)
!
!--------------------check if grid-values occur twice-------------------
!
do i=2, n1d-1
   val1=r1d(i-1)
   val2=r1d(i)
   if(abs(val2-val1).lt.1.d-14) then
!insert value on half range (to upper grid point)
      r1d(i)=(r1d(i+1)+r1d(i))/two
   endif
enddo
!
!same on outer boundary
val1=r1d(n1d)
val2=r1d(n1d-1)
if(abs(val2-val1).lt.1.d-14) then
!insert value on half range (to lower grid point)
   r1d(n1d-1)=(r1d(n1d-2)+r1d(n1d-1))/two
endif
!
!-----make outer part of grid equidistant in order that delta(r1d)------
!---------------------is not decreasing---------------------------------
!
do i=2, n1d
   delr1=r1d(i)-r1d(i-1)
   max_delr=delr1*(n1d-i)
   maxr=r1d(i)+max_delr
   if(maxr.lt.rmax) then
      indx=i
   endif
enddo

if(n1d.ne.indx) then
   delr=(rmax-r1d(indx))/(n1d-indx)
   do i=indx+1, n1d
      r1d(i)=r1d(i-1)+delr
   enddo
endif
!
!
end subroutine grid1d_final
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine grid1d_final_2
!
!-----------------------------------------------------------------------
!-------------------sets up 1-dimensional radial grid:------------------
!-----------------------------------------------------------------------
!   n1d_t grid points equidistant in log(tau_thomson)
!         (estimated from beta-vel-law)
!         where increments of delr are adapted if delr(equidistant velocity) lower
!
!   n1d_r grid points equidistant in radius
!-----------------------------------------------------------------------
!
use prog_type
use dime1d, only: n1d, n1d_dum, n1d_t, n1d_r, r1d, r1d_dum, delv
use params_input, only: teff, tmin, vmax, vmin, beta, xmloss, yhe, hei, rmax, vth_fiducial
use fund_const
use params_stellar, only: sr
use inf_reg, only: rmin
use mod_sort, only: bubblesort
use mod_grid, only: recalc_grid1d
!
!
!
implicit none
!
! ... local scalars
integer(i4b) :: i, j, k
integer(i4b) :: gp_used, gp_free            !used and free grid points
integer(i4b) :: ind_max, del_ind, indx, indx_1, indx_2, err
real(dp) :: b, velr
real(dp) :: xmlosscgs
real(dp) :: c1, c2
real(dp) :: taumin, taumax, deltau, deltau_log, delr, delr_v, delr_t
real(dp) :: v_0, v_1
real(dp) :: max_delr, max_delr_old, max_deltau, max_deltau_old, delr0, delr1, delr2, maxr
real(dp) :: delvcgs
real(dp) :: sigem
real(dp) :: taur_0
real(dp) :: val1, val2
!
! ... local characters
character(len=50) :: enter
!
! ... local functions
real(dp) :: bvel
!
n1d_dum = n1d_t + n1d_r
!
!------------------calculate/define constants needed--------------------
!
!calculate b-factor for beta-velocity-law
b=one-(vmin/vmax)**(one/beta)
!
c1=(one+four*yhe)*cgs_mp
c2=(one+hei*yhe)/c1
!
!calculate mass loss rate in cgs-units
xmlosscgs=xmloss*xmsu/(365.25d0*24.d0*60.d0*60.d0)
!
!calculate tau_min, set tau_max
taumin=(xmlosscgs*sigmae*c2/(four*pi)) / (vmax*1.d5 * rmax*sr)
taumax=20
!
!
!allocate radial dummy grid
allocate(r1d_dum(n1d_dum), stat=err)
   if(err.ne.0) stop 'allocation error grid1d_final_2: r1d_dum'
!
!-----------------------------------------------------------------------
!
!calculate n1d_t radial grid points that are equidistant in tau
!
r1d_dum=zero
!
r1d_dum(1)=rmin
!
!!use a smaller velocity-width for radial grid,
!since 3d grid will stretch radial grid anyways (delv->delv/3.d0, or delv-> delv/two, or delv-> delv/1.5d0)
delvcgs=delv*vth_fiducial/3.
!
!
outer: do j=1, 1000
!outer cycle adapts the optical depth-steps deltau(controlled by taumax)
!            until maximum radial grid point is less than rmax
!
   deltau=(taumax-taumin)/(n1d_t-1)
!
   do i=1, n1d_t-1
!calculate velocity at given radius
      velr=vmax*1.d5*(one - b*rmin/r1d_dum(i))**(beta)
!calculate opacity at given radius
      sigem=sigmae*c2*xmlosscgs / (four*pi*r1d_dum(i)*r1d_dum(i)*velr*sr)
!calculate radial grid step at given radius if equidistant tau is used
      delr_t= -deltau/sigem
!calculate radial grid step at given radius if equidistant velocity is used
      v_0=bvel(r1d_dum(i), vmax*1.d5, b, beta)
      delr_v=(delvcgs*r1d_dum(i)*r1d_dum(i))/(b*beta*v_0) * (one-b/r1d_dum(i))
!calculate radius for next grid point where increments of delr are taken from
!     the smaller values of delr(equidistant in tau), delr(equidistant in velocity)
      if(abs(delr_v).lt.abs(delr_t)) then
         r1d_dum(i+1)=r1d_dum(i)+delr_v
!         write(*,fmt='(i4, 3(e20.8))') i, r1d_dum(i-1)/sr, r1d_dum(i)/sr, r1d_dum(i+1)/sr
      else
         r1d_dum(i+1)=r1d_dum(i)-delr_t
      endif
!go to outer cycle if the calculated radius exceed rmax
      if(r1d_dum(i+1).gt.rmax) then
         taumax=taumax*99.d0/100.d0
         exit
      endif
   enddo
!
!stop outer cycle if the maximum calculated radius is less than rmax
   if(maxval(r1d_dum).lt.rmax) then
!      write(*,*) j, 'max', maxval(r1d_dum)/sr
      exit outer
   endif
!
enddo outer
!
if(maxval(r1d_dum).gt.rmax) then
   write(*,*) 'error in calculating grid: adapt cycle'
   stop
endif
!
gp_used=n1d_t
gp_free=n1d_dum-n1d_t
!
!
!-----------------------------------------------------------------------
!
!insert n1d_r grid points at the center of maximum distance 
!                   between neighbouring radial grid points

do i=2, n1d_t
!calculate delta_r between neighbouring grid points
   delr1=r1d_dum(i)-r1d_dum(i-1)
!calculate maximum radius if grid points would be equidistant in r
   max_delr=delr1*(n1d_dum-i)
   maxr=r1d_dum(i)+max_delr
!save index where maximum radius (if equidistant in r) is lt rmax
   if(maxr.lt.rmax) then
      indx=i
   endif
enddo
!
!
!calculate new delta_r which is adopted from r(indx)
if(n1d_dum.ne.indx) then
   delr=(rmax-r1d_dum(indx))/(n1d_dum-indx)
!
   do i=indx+1, n1d_dum
      r1d_dum(i)=r1d_dum(i-1)+delr
   enddo
endif
!
!------------setting up actual used radial grid from dummy grid---------
!
allocate(r1d(n1d), stat=err)
   if(err.ne.0) stop 'error in grid1d_final2: allocation'
!
call recalc_grid1d(r1d_dum, n1d_dum, n1d, 2, r1d)
!
r1d(n1d-1)=rmin
r1d(n1d)=rmax
!
!sort r1d-grid
call bubblesort(r1d,n1d)
!
!--------------------check if grid-values occur twice-------------------
!
do i=2, n1d-1
   val1=r1d(i-1)
   val2=r1d(i)
   if(abs(val2-val1).lt.1.d-14) then
!insert value on half range (to upper grid point)
      r1d(i)=(r1d(i+1)+r1d(i))/two
   endif
enddo
!
!same on outer boundary
val1=r1d(n1d)
val2=r1d(n1d-1)
if(abs(val2-val1).lt.1.d-14) then
!insert value on half range (to lower grid point)
   r1d(n1d-1)=(r1d(n1d-2)+r1d(n1d-1))/two
endif
!
!-----make outer part of grid equidistant in order that delta(r1d)------
!---------------------is not decreasing---------------------------------
!
do i=2, n1d
   delr1=r1d(i)-r1d(i-1)
   max_delr=delr1*(n1d-i)
   maxr=r1d(i)+max_delr
   if(maxr.lt.rmax*sr) then
      indx=i
   endif
enddo

if(n1d.ne.indx) then
   delr=(rmax-r1d(indx))/(n1d-indx)
   do i=indx+1, n1d
      r1d(i)=r1d(i-1)+delr
   enddo
endif

!open(1, file='TRASH/grid.dat')
!   write(1,*), n1d
!   write(1,*) r1d
!   write(1,*)
!   write(1,*) n1d_dum
!   write(1,*) r1d_dum
!close(1)
!stop
!
end subroutine grid1d_final_2
!
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine grid1d_tau_equi
!
!-----------------------------------------------------------------------
!-------------------sets up 1-dimensional radial grid:------------------
!---------equidistant in tau_thomson (estimated from beta-vel-law)------
!-----------------------------------------------------------------------
!
use prog_type
use dime1d, only: n1d, n1d_dum, n1d_t, n1d_r, r1d, r1d_dum
use params_input, only: teff, vmax, vmin, beta, xmloss, yhe, hei, rmax
use params_stellar, only: sr
use fund_const
use inf_reg, only: rmin
use mod_sort, only: bubblesort
use mod_grid, only: recalc_grid1d
!
implicit none
!
! ... local scalars
integer(i4b) :: i, j, k
integer(i4b) :: gp_used, gp_free            !used and free grid points
integer(i4b) :: ind_max, del_ind, indx, err
real(dp) :: b, velr
real(dp) :: xmlosscgs
real(dp) :: c1, c2
real(dp) :: taumin, taumax, deltau, delr
real(dp) :: max_delr, max_delr_old, max_deltau, max_deltau_old, delr1, delr2, maxr
real(dp) :: delvcgs
real(dp) :: sigem
real(dp) :: val1, val2
!
! ... local characters
character(len=50) :: enter
!
n1d_dum = n1d_t + n1d_r
!
!------------------calculate/define constants needed--------------------
!
!calculate b-factor for beta-velocity-law
b=one-(vmin/vmax)**(one/beta)
!
c1=(one+four*yhe)*cgs_mp
c2=(one+hei*yhe)/c1
!
!calculate mass loss rate in cgs-units
xmlosscgs=xmloss*xmsu/(365.25d0*24.d0*60.d0*60.d0)
!
!calculate tau_min, set tau_max
taumin=(xmlosscgs*sigmae*c2/(four*pi)) / (vmax*1.d5 * rmax*sr)
taumax=20
!
!calculate equidistant tau-steps
deltau=(taumax-taumin)/(n1d_dum-1)
!
!allocate radial dummy grid
allocate(r1d_dum(n1d_dum), stat=err)
   if(err.ne.0) stop 'allocation error grid1d_tau_equi: r1d_dum'
!
!-----------------------------------------------------------------------
!
!calculate n1d_t radial grid points that are equidistant in tau
r1d_dum=zero
r1d_dum(1)=rmin
!
outer: do j=1, 1000
!outer cycle adapts the optical depth-steps deltau(controlled by taumax)
!            until maximum radial grid point is less than rmax
   deltau=(taumax-taumin)/(n1d_t-1)
!
   do i=1, n1d_t-1
!calculate velocity at given radius
      velr=vmax*1.d5*(one - b*rmin/r1d_dum(i))**(beta)
!calculate opacity at given radius
      sigem=sigmae*c2*xmlosscgs / (four*pi*r1d_dum(i)*r1d_dum(i)*velr*sr)
!calculate radial grid step at given radius
      delr= -deltau/sigem
!calculate radius for next grid point
      r1d_dum(i+1)=r1d_dum(i)-delr
!go to outer cycle if the calculated radius exceed rmax
      if(r1d_dum(i+1).gt.rmax) then
         taumax=taumax*99.d0/100.d0
         exit
      endif
   enddo
!
!stop outer cycle if the maximum calculated radius is less than rmax
   if(maxval(r1d_dum).lt.rmax) then
!      write(*,*) j, 'max', maxval(r1d_dum)
      exit outer
   endif
!
enddo outer
!
if(maxval(r1d_dum).gt.rmax) then
   stop 'error in grid1d_equi: adapt cycle'
   stop
endif
!
!-----------------------------------------------------------------------
!
!insert n1d_r grid points at the center of maximum distance 
!                   between neighbouring radial grid points
do i=2, n1d_t
!calculate delta_r between neighbouring grid points
   delr1=r1d_dum(i)-r1d_dum(i-1)
!calculate maximum radius if grid points would be equidistant in r
   max_delr=delr1*(n1d_dum-i)
   maxr=r1d_dum(i)+max_delr
!save index where maximum radius (if equidistant in r) is lt rmax
   if(maxr.lt.rmax) then
      indx=i
   endif
enddo
!
!calculate new delta_r which is adopted from r(indx)
delr=(rmax-r1d_dum(indx))/(n1d_dum-indx)
!
do i=indx+1, n1d_dum
   r1d_dum(i)=r1d_dum(i-1)+delr
enddo
!
!------------setting up actual used radial grid from dummy grid---------
!
allocate(r1d(n1d), stat=err)
   if(err.ne.0) stop 'error in grid1d_tau_equi: allocation'
!
call recalc_grid1d(r1d_dum, n1d_dum, n1d, 2, r1d)
!
r1d(n1d-1)=rmin
r1d(n1d)=rmax
!
!sort r1d-grid
call bubblesort(r1d,n1d)
!
!--------------------check if grid-values occur twice-------------------
!
do i=2, n1d-1
   val1=r1d(i-1)
   val2=r1d(i)
   if(abs(val2-val1).lt.1.d-14) then
!insert value on half range (to upper grid point)
      r1d(i)=(r1d(i+1)+r1d(i))/two
   endif
enddo
!
!same on outer boundary
val1=r1d(n1d)
val2=r1d(n1d-1)
if(abs(val2-val1).lt.1.d-14) then
!insert value on half range (to lower grid point)
   r1d(n1d-1)=(r1d(n1d-2)+r1d(n1d-1))/two
endif
!
!-----make outer part of grid equidistant in order that delta(r1d)------
!---------------------is not decreasing---------------------------------
!
do i=2, n1d
   delr1=r1d(i)-r1d(i-1)
   max_delr=delr1*(n1d-i)
   maxr=r1d(i)+max_delr
   if(maxr.lt.rmax*sr) then
      indx=i
   endif
enddo
!
if(n1d.ne.indx) then
   delr=(rmax-r1d(indx))/(n1d-indx)
   do i=indx+1, n1d
      r1d(i)=r1d(i-1)+delr
   enddo
endif
!
!
!
end subroutine grid1d_tau_equi
!!
!!-----------------------------------------------------------------------
!!-----------------------------------------------------------------------
!!-----------------------------------------------------------------------
!!
subroutine grid1d_r_equi
!
!-----------------------------------------------------------------------
!-------------------sets up 1-dimensional radial grid:------------------
!------------------------equidistant in radius--------------------------
!-----------------------------------------------------------------------
!
use prog_type
use dime1d, only: n1d, n1d_dum, n1d_t, n1d_r, r1d, r1d_dum
use params_input, only: rmax
use inf_reg, only: rmin
use fund_const

!
implicit none
!
! ... arguments
!
! ... local scalars
integer(i4b) :: i, err
real(dp) :: delr
!
n1d_dum = n1d_t + n1d_r
!
allocate(r1d(n1d), stat=err)
   if(err.ne.0) stop 'error in grid1d_r_equi: allocation'
allocate(r1d_dum(n1d_dum), stat=err)
   if(err.ne.0) stop 'error in grid1d_r_equi: allocation'
!
do i=1, n1d_dum
   r1d_dum(i) = rmin + (i-1)*(rmax-rmin)/float(n1d_dum-1)
enddo
!
delr=(rmax-rmin)/(n1d-1)
r1d(1)=rmin
do i=2, n1d
   r1d(i)=r1d(i-1)+delr
enddo
!
!
end subroutine grid1d_r_equi
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine grid1d_r_log
!
!-----------------------------------------------------------------------
!-------------------sets up 1-dimensional radial grid:------------------
!------------------------equidistant in log-space-----------------------
!
use prog_type
use fund_const
use dime1d, only: n1d_dum, n1d_t, n1d_r, r1d_dum, n1d, r1d
use params_input, only: rmax
use inf_reg, only: rmin
!
implicit none
!
! ... arguments
!
! ... local scalars
integer(i4b) :: i, err
real(dp) :: delr
!
n1d_dum = n1d_t + n1d_r
!
!------------------calculate radial dummy grid--------------------------
!
!allocate radial dummy grid
allocate(r1d_dum(n1d_dum), stat=err)
   if(err.ne.0) stop 'allocation error grid1d_r_equi: r1d_dum'
!
delr=log10(rmax/rmin)/(n1d_dum-1)
r1d_dum(1)=rmin
do i=2, n1d_dum
   r1d_dum(i)=r1d_dum(i-1)*ten**delr
enddo
!
!
!--------------------calculate actual radial grid-----------------------
!
allocate(r1d(n1d), stat=err)
   if(err.ne.0) stop 'error in grid1d_vel_equi: allocation'
!
delr=log10(rmax/rmin)/(n1d-1)
r1d(1)=rmin
do i=2, n1d
   r1d(i) = r1d(i-1)*ten**delr
enddo
!
!
!
!
end subroutine grid1d_r_log
!!
!!-----------------------------------------------------------------------
!!-----------------------------------------------------------------------
!!-----------------------------------------------------------------------
!!
!subroutine allocate_global1d
!!
!!--------------------allocates global 1d-grids--------------------------
!!
!use prog_type
!use dime1d, only: n1d, gradv1d, opalbar1d, opath1d, r1d, &
!                  t1d, tauth1d, vel1d, vth1d
!!
!implicit none
!!
!if(n1d.eq.0) then
!   write(*,*) 'need to specify n1d to allocate arrays'
!   stop
!endif
!!
!!allocate arrays from module dime1d, with dimension n1d
!!
!allocate(gradv1d(n1d))
!allocate(opath1d(n1d))
!allocate(opalbar1d(n1d))
!allocate(r1d(n1d))
!allocate(t1d(n1d))
!allocate(tauth1d(n1d))
!allocate(vel1d(n1d))
!allocate(vth1d(n1d))
!!
!!allocate arrays from module laop
!!   (not needed anymore)
!!
!!allocate arrays from module testlamb
!!   (not needed anymore)
!!
!!
!end subroutine allocate_global1d
!!
!!-----------------------------------------------------------------------
!!-----------------------------------------------------------------------
!!-----------------------------------------------------------------------
!!
!subroutine allocate_global1d_cr
!!
!!--------------------allocates global 1d-grids--------------------------
!!
!use prog_type
!use dime1d, only: n1d_cr, alocont_diag1d_cr, opalbar1d_cr, opath1d_cr, &
!                  r1d_cr, t1d_cr, vel1d_cr, betasob1d_cr, mint1d_cr, &
!                  mintbar1d_cr, scont1d_cr, sline1d_cr, ssobo1d_cr, &
!                  ssoboc1d_cr, normalization1d_cr, gradv1d_cr, &
!                  aloline_diag1d_cr, vth1d_cr
!use iter, only: eps1d_cr, itmaxc, itmaxl, epshistory1d_cr
!use dimz, only: ndzmax, z
!!
!implicit none
!!
!! ... local scalars
!integer(i4b) :: i
!integer(i4b) :: indx, itmax
!!
!!------------------------calculate n1d_cr-------------------------------
!!
!n1d_cr=0
!do i=ndzmax/2+1, ndzmax
!   if(z(i).ge.one) then
!      n1d_cr=n1d_cr+1
!   endif
!enddo
!!
!if(n1d_cr.eq.0) stop 'error in allocate_global1d_cr: n1d_cr = 0'
!!
!!-----------------------------------------------------------------------
!!
!!allocate arrays from module dime1d, with dimension n1d_cr
!!
!allocate(alocont_diag1d_cr(n1d_cr))
!allocate(aloline_diag1d_cr(n1d_cr))
!allocate(opath1d_cr(n1d_cr))
!allocate(opalbar1d_cr(n1d_cr))
!allocate(r1d_cr(n1d_cr))
!allocate(t1d_cr(n1d_cr))
!allocate(vel1d_cr(n1d_cr))
!allocate(gradv1d_cr(n1d_cr))
!allocate(vth1d_cr(n1d_cr))
!!
!!allocate arrays from module laop
!!   (not needed anymore)
!!
!!allocate arrays from module testlamb
!!   (not needed anymore)
!!
!allocate(betasob1d_cr(n1d_cr))
!allocate(mint1d_cr(n1d_cr))
!allocate(mintbar1d_cr(n1d_cr))
!allocate(scont1d_cr(n1d_cr))
!allocate(sline1d_cr(n1d_cr))
!allocate(ssobo1d_cr(n1d_cr))
!allocate(ssoboc1d_cr(n1d_cr))
!allocate(normalization1d_cr(n1d_cr))
!!
!!allocate arrays from module iter
!!
!itmax=max(itmaxc,itmaxl)
!allocate(eps1d_cr(n1d_cr))
!allocate(epshistory1d_cr(n1d_cr,itmax))
!epshistory1d_cr=zero
!
!!
!end subroutine allocate_global1d_cr
!!
!!-----------------------------------------------------------------------
!!-----------------------------------------------------------------------
!!-----------------------------------------------------------------------
!!
!subroutine grid1d_cr
!!
!use prog_type
!use dimz, only: z, ndzmax
!use dimx, only: ndxmax
!use dimy, only: ndymax
!use dime1d, only: n1d_cr, r1d_cr, vel1d_cr, gradv1d_cr, opath1d_cr, opalbar1d_cr, t1d_cr, vth1d_cr, &
!                  n1d, r1d, vel1d, gradv1d, opath1d, opalbar1d, t1d, vth1d
!use dime3d, only: opalbar3d, vel3d
!use params_stellar, only: sr
!use params_input, only: vmax, vmin, beta
!use diffapprox, only: opathboundary
!!
!implicit none
!!
!! ... local scalars
!integer(i4b) :: i
!real(dp) :: b
!!
!!calculate radial grid
!do i=1, n1d_cr
!   r1d_cr(i)=z(ndzmax-n1d_cr+i)*sr
!enddo
!!
!call interpol_lin(n1d_cr, n1d, log10(r1d_cr), log10(r1d), vel1d_cr, log10(vel1d))
!call interpol_lin(n1d_cr, n1d, log10(r1d_cr), log10(r1d), gradv1d_cr, log10(gradv1d))
!call interpol_lin(n1d_cr, n1d, log10(r1d_cr), log10(r1d), opath1d_cr, log10(opath1d))
!call interpol_lin(n1d_cr, n1d, log10(r1d_cr), log10(r1d), opalbar1d_cr, log10(opalbar1d))
!call interpol_lin(n1d_cr, n1d, log10(r1d_cr), log10(r1d), t1d_cr, log10(t1d))
!call interpol_lin(n1d_cr, n1d, log10(r1d_cr), log10(r1d), vth1d_cr, log10(vth1d))
!!
!!calculate real values since interpolation was performed in log-space
!vel1d_cr=ten**vel1d_cr
!gradv1d_cr=ten**gradv1d_cr
!opath1d_cr=ten**opath1d_cr
!opalbar1d_cr=ten**opalbar1d_cr
!t1d_cr=ten**t1d_cr
!vth1d_cr=ten**vth1d_cr
!!
!opathboundary=opath1d_cr(1)
!!
!!*************************debug open************************************
!!
!!
!!b=one-(vmin/vmax)**(one/beta)
!!
!!do i=1, n1d_cr
!!   vel1d_cr(i)=vel3d(ndxmax/2+1, ndymax/2+1, ndzmax-n1d_cr+i)
!!   opalbar1d_cr(i) = opalbar3d(ndxmax/2+1, ndymax/2+1, ndzmax-n1d_cr+i)
!!   gradv1d_cr(i)=(vmax*1.d5*beta*(one-b*(sr/r1d_cr(i)))**(beta-1))*b*sr/(r1d_cr(i)*r1d_cr(i))
!!enddo
!!
!!**************************debug close**********************************
!!
!!
!!
!end subroutine grid1d_cr
!!
!!-----------------------------------------------------------------------
!!-----------------------------------------------------------------------
!!-----------------------------------------------------------------------
!!
!subroutine prgrid
!!
!!-----------------------------------------------------------------------
!!--------------------prints out 1d model grid---------------------------
!!-----------------------------------------------------------------------
!!
!use prog_type
!use dime1d, only: n1d, vel1d, r1d, tauth1d, opath1d, opalbar1d, t1d
!use params_input, only: vmax, teff, na
!use params_stellar, only: sr
!!
!implicit none
!!
!integer(i4b) :: i
!!
!!
!!
!write(*,19) 'vel [vinf]', 'r [r_star]', 'tau_th', 'log(chi_th)', 'opalbar', 't [t_eff]'
!!
!do i=1, n1d
!   write(*,20) vel1d(i)/(vmax*1.d5), r1d(i)/sr, tauth1d(i), &
!               log10(opath1d(i)), opalbar1d(i), t1d(i)/teff
!enddo
!
!19 format (6a20)
!20 format (6es20.8)
!!
!!
!!
!end subroutine prgrid
!!
!!-----------------------------------------------------------------------
!!-----------------------------------------------------------------------
!!-----------------------------------------------------------------------
!!
!subroutine check_grid1d
!!
!use prog_type
!use fund_const
!use dime1d, only: tauth1d, n1d, r1d, vel1d, delv
!use warnings, only: warn_grid1, warn_grid2, warn_grid3
!use params_input, only: teff, tmin, na
!use params_stellar, only: vmicro, sr
!!
!implicit none
!!
!! ... local scalars
!integer(i4b) :: i
!integer(i4b) :: indx_1, indx_2
!real(dp) :: max_delr, max_delr_old, max_deltau, max_deltau_old, &
!            delr1, delr2, v_0, vth_min, delvcgs
!!
!! ... local functions
!real(dp) :: vthermal
!!
!!-----------find maxmum steps in delta-tau and delta-r------------------
!!
!max_delr_old=zero
!max_deltau_old=zero
!!
!do i=2, n1d
!   max_delr=r1d(i)-r1d(i-1)
!   max_deltau=tauth1d(i-1)-tauth1d(i)
!!
!   if(max_delr.gt.max_delr_old) then
!      max_delr_old=max_delr
!      indx_1=i
!   endif
!!
!   if(max_deltau.gt.max_deltau_old) then
!      max_deltau_old=max_deltau
!      indx_2=i
!   endif
!enddo
!!
!write(*,*) 'maximum delta tau:', max_deltau_old, indx_2
!write(*,*) 'maximum delta r  :', max_delr_old/sr, indx_1
!write(*,*)
!!
!if(max_deltau_old.gt.one/3.d0) then
!   warn_grid3=.true.
!endif
!!
!!
!!---------check if grid is spatially increasing (or constant)-----------
!!-------------- if velocity steps are okay------------------------------
!!
!vth_min=vthermal(vmicro, tmin, na)
!delvcgs=delv*vth_min
!!
!do i=2, n1d-1
!   delr1=r1d(i)-r1d(i-1)
!   delr2=r1d(i+1)-r1d(i)
!   v_0=vel1d(i)-vel1d(i-1)
!   if(delr2.lt.delr1-delr1*1.d-5) then
!      warn_grid1=.true.
!!      write(*,*) 'spatial grid is not increasing/constant: increase grid points n1d_t'
!!      write(*,*) delr2/sr, delr1/sr, (delr2-delr1)/sr
!!      stop
!   endif
!   if(v_0.gt.delvcgs+0.1d0*delvcgs) then
!      warn_grid2=.true.
!!      write(*,*) i, v_0
!!      write(*,*) 'steps in velocity space are too large: increase n1d_t or delv'
!!      stop
!   endif
!enddo
!!
!end subroutine check_grid1d
!!
!!-----------------------------------------------------------------------
!!-----------------------------------------------------------------------
!!-----------------------------------------------------------------------
!!
!subroutine interp_grid1d
!!
!!
!!-----------------------------------------------------------------------
!!---------------interpolating grids from external models----------------
!!-----------------------------------------------------------------------
!!
!!--------interpolation of model-grids onto radial grid to obtain--------
!!   velocity
!!   density
!!   electron-number-density
!!   hydrogen-number-density
!!   thomson opacity
!!   thomson optical depth
!!   temperature grid
!!-----------------------------------------------------------------------
!!
!use prog_type
!use fund_const
!use mod_ext_spherical
!use dime1d, only: n1d, gradv1d, opath1d, opalbar1d, r1d, t1d, tauth1d, vel1d, vth1d
!use params_input, only: teff, vmax, vmin, beta, xmloss, yhe, hei, rmax
!use params_stellar, only: sr
!use options, only: opt_interpol_1d
!
!!
!implicit none
!!
!! ... local scalars
!integer(i4b) :: i
!real(dp) :: b, c1, c2, xmlosscgs, taumin
!real(dp) :: delr
!real(dp) :: max_delr_old, max_deltau_old, max_delr, max_deltau
!!
!!-----------------------check allocation--------------------------------
!!
!if(.not.allocated(r1d_mod_ext)) stop 'ext_mod previously deallocated'
!if(.not.allocated(vel1d_mod_ext)) stop 'ext_mod previously deallocated'
!if(.not.allocated(gradv1d_mod_ext)) stop 'ext_mod previously deallocated'
!if(.not.allocated(opath1d_mod_ext)) stop 'ext_mod previously deallocated'
!if(.not.allocated(opalbar1d_mod_ext)) stop 'ext_mod previously deallocated'
!if(.not.allocated(t1d_mod_ext)) stop 'ext_mod previously deallocated'
!if(.not.allocated(vth1d_mod_ext)) stop 'ext_mod previously deallocated'
!!
!if(.not.allocated(r1d)) stop 'error interp_grid1d: r1d not allocated'
!if(.not.allocated(vel1d)) stop 'error interp_grid1d: vel1d not allocated'
!if(.not.allocated(t1d)) stop 'error interp_grid1d: t1d not allocated'
!if(.not.allocated(gradv1d)) stop 'error interp_grid1d: gradv1d not allocated'
!if(.not.allocated(opath1d)) stop 'error interp_grid1d: opath1d not allocated'
!if(.not.allocated(opalbar1d)) stop 'error interp_grid1d: opath1d not allocated'
!if(.not.allocated(tauth1d)) stop 'error interp_grid1d: tauth1d not allocated'
!if(.not.allocated(vth1d)) stop 'error interp_grid1d: vth1d not allocated'
!!
!!--------------------------interpolation--------------------------------
!!
!select case(opt_interpol_1d)
!   case(0)
!!
!!-------------perform interpolation linear in logspace------------------
!!
!      call interpol_lin(n1d, n1d_mod_ext, log10(r1d), log10(r1d_mod_ext), vel1d, log10(vel1d_mod_ext))
!      call interpol_lin(n1d, n1d_mod_ext, log10(r1d), log10(r1d_mod_ext), gradv1d, log10(gradv1d_mod_ext))
!      call interpol_lin(n1d, n1d_mod_ext, log10(r1d), log10(r1d_mod_ext), opath1d, log10(opath1d_mod_ext))
!      call interpol_lin(n1d, n1d_mod_ext, log10(r1d), log10(r1d_mod_ext), opalbar1d, log10(opalbar1d_mod_ext))
!      call interpol_lin(n1d, n1d_mod_ext, log10(r1d), log10(r1d_mod_ext), t1d, log10(t1d_mod_ext))
!      call interpol_lin(n1d, n1d_mod_ext, log10(r1d), log10(r1d_mod_ext), vth1d, log10(vth1d_mod_ext))
!      !calculate real values since interpolation was performed in log-space
!      vel1d=ten**vel1d
!      gradv1d=ten**gradv1d
!      opath1d=ten**opath1d
!      opalbar1d=ten**opalbar1d
!      t1d=ten**t1d
!      vth1d=ten**vth1d
!!
!   case(1)
!!
!!----------------------perform cubic interpolation----------------------
!!
!      call interpol_cube(n1d, n1d_mod_ext, r1d, r1d_mod_ext, vel1d, vel1d_mod_ext)
!      call interpol_cube(n1d, n1d_mod_ext, r1d, r1d_mod_ext, gradv1d, gradv1d_mod_ext)
!      call interpol_cube(n1d, n1d_mod_ext, r1d, r1d_mod_ext, opath1d, opath1d_mod_ext)
!      call interpol_cube(n1d, n1d_mod_ext, r1d, r1d_mod_ext, opalbar1d, opalbar1d_mod_ext)
!      call interpol_cube(n1d, n1d_mod_ext, r1d, r1d_mod_ext, t1d, t1d_mod_ext)
!      call interpol_cube(n1d, n1d_mod_ext, r1d, r1d_mod_ext, vth1d, vth1d_mod_ext)
!!
!   case(2)
!!
!!--------------------perform cubic spline interpolation-----------------
!!
!      call interpol_spline(n1d, n1d_mod_ext, r1d, r1d_mod_ext, vel1d, vel1d_mod_ext)
!      call interpol_spline(n1d, n1d_mod_ext, r1d, r1d_mod_ext, gradv1d, gradv1d_mod_ext)
!      call interpol_spline(n1d, n1d_mod_ext, r1d, r1d_mod_ext, opath1d, opath1d_mod_ext)
!      call interpol_spline(n1d, n1d_mod_ext, r1d, r1d_mod_ext, opalbar1d, opalbar1d_mod_ext)
!      call interpol_spline(n1d, n1d_mod_ext, r1d, r1d_mod_ext, t1d, t1d_mod_ext)
!      call interpol_spline(n1d, n1d_mod_ext, r1d, r1d_mod_ext, vth1d, vth1d_mod_ext)
!!
!   case(3)
!!
!!------------------perform monotonic cubic interpolation----------------
!!
!      call interpol_cube_mono(n1d, n1d_mod_ext, r1d, r1d_mod_ext, vel1d, vel1d_mod_ext)
!      call interpol_cube_mono(n1d, n1d_mod_ext, r1d, r1d_mod_ext, gradv1d, gradv1d_mod_ext)
!      call interpol_cube_mono(n1d, n1d_mod_ext, r1d, r1d_mod_ext, opath1d, opath1d_mod_ext)
!      call interpol_cube_mono(n1d, n1d_mod_ext, r1d, r1d_mod_ext, opalbar1d, opalbar1d_mod_ext)
!      call interpol_cube_mono(n1d, n1d_mod_ext, r1d, r1d_mod_ext, t1d, t1d_mod_ext)
!      call interpol_cube_mono(n1d, n1d_mod_ext, r1d, r1d_mod_ext, vth1d, vth1d_mod_ext)
!!
!   case default
!      stop 'set option opt_interpol_1d'
!end select
!!
!!----------------setting up thomson optical depth-----------------------
!!
!!calculate constants needed
!xmlosscgs=xmloss*xmsu/(365.25d0*24.d0*60.d0*60.d0)
!c1=(one+four*yhe)*cgs_mp
!c2=(one+hei*yhe)/c1
!!
!taumin=(xmlosscgs*sigmae*c2/(four*pi)) / (vmax*1.d5 * r1d(n1d))
!!
!tauth1d(n1d)=taumin
!!
!do i=n1d-1, 1, -1
!   tauth1d(i)=tauth1d(i+1)+(opath1d(i)+opath1d(i+1))/two * (r1d(i+1)-r1d(i))
!enddo
!!
!!-------------------------deallocate model arrays-----------------------
!!
!!from now on: can deallocate model-arrays (but for r1d_mod_ext, t1d_mod_ext, opath1d_mod_ext: diffus-approx)
!deallocate(vel1d_mod_ext)
!deallocate(gradv1d_mod_ext)
!deallocate(opalbar1d_mod_ext)
!deallocate(vth1d_mod_ext)
!
!end subroutine interp_grid1d
!
!!
!!-----------------------------------------------------------------------
!!-----------------------------------------------------------------------
!!-----------------------------------------------------------------------
!!
!subroutine test_interpol_1d(dirout)
!!
!!
!!-----------------------------------------------------------------------
!!---------------interpolating grids from external models----------------
!!-----------------------------------------------------------------------
!!
!!--------interpolation of model-grids onto radial grid to obtain--------
!!   velocity
!!   density
!!   frequency integrated line opacity
!!   thomson opacity
!!   thomson optical depth
!!   temperature grid
!!-----------------------------------------------------------------------
!!
!use prog_type
!use fund_const
!use mod_ext_spherical
!use dime1d, only: n1d, r1d
!use params_input, only: teff, vmax, vmin, beta, xmloss, yhe, hei, rmax
!use params_stellar, only: sr
!use options, only: opt_interpol_1d
!
!!
!implicit none
!!
!! ... arguments
!character(len=100), intent(in) :: dirout
!!
!! ... local scalars
!integer(i4b) :: i
!!
!! ... local arrays
!real(dp), dimension(:), allocatable :: vel1d_lin, opath1d_lin, opalbar1d_lin, gradv1d_lin, t1d_lin
!real(dp), dimension(:), allocatable :: vel1d_cube, opath1d_cube, opalbar1d_cube, gradv1d_cube, t1d_cube
!real(dp), dimension(:), allocatable :: vel1d_cube_mono, opath1d_cube_mono, opalbar1d_cube_mono, gradv1d_cube_mono, t1d_cube_mono
!real(dp), dimension(:), allocatable :: vel1d_spline, opath1d_spline, opalbar1d_spline, gradv1d_spline, t1d_spline
!!
!!-----------------------check allocation--------------------------------
!!
!if(.not.allocated(r1d_mod_ext)) stop 'ext_mod previously deallocated'
!if(.not.allocated(vel1d_mod_ext)) stop 'ext_mod previously deallocated'
!if(.not.allocated(opath1d_mod_ext)) stop 'ext_mod previously deallocated'
!if(.not.allocated(opalbar1d_mod_ext)) stop 'ext_mod previously deallocated'
!if(.not.allocated(gradv1d_mod_ext)) stop 'ext_mod previously deallocated'
!if(.not.allocated(t1d_mod_ext)) stop 'ext_mod previously deallocated'
!!
!if(.not.allocated(r1d)) stop 'error interp_grid1d: r1d not allocated'
!!
!!-----------------------allocate local arrays---------------------------
!!
!allocate(vel1d_lin(n1d))
!allocate(opath1d_lin(n1d))
!allocate(opalbar1d_lin(n1d))
!allocate(gradv1d_lin(n1d))
!allocate(t1d_lin(n1d))
!!
!allocate(vel1d_cube(n1d))
!allocate(opath1d_cube(n1d))
!allocate(opalbar1d_cube(n1d))
!allocate(gradv1d_cube(n1d))
!allocate(t1d_cube(n1d))
!!
!allocate(vel1d_cube_mono(n1d))
!allocate(opath1d_cube_mono(n1d))
!allocate(opalbar1d_cube_mono(n1d))
!allocate(gradv1d_cube_mono(n1d))
!allocate(t1d_cube_mono(n1d))
!!
!allocate(vel1d_spline(n1d))
!allocate(opath1d_spline(n1d))
!allocate(opalbar1d_spline(n1d))
!allocate(gradv1d_spline(n1d))
!allocate(t1d_spline(n1d))
!!
!!--------------------------interpolation--------------------------------
!!
!!-------------perform interpolation linear in logspace------------------
!!
!call interpol_lin(n1d, n1d_mod_ext, log10(r1d), log10(r1d_mod_ext), vel1d_lin, log10(vel1d_mod_ext))
!call interpol_lin(n1d, n1d_mod_ext, log10(r1d), log10(r1d_mod_ext), opath1d_lin, log10(opath1d_mod_ext))
!call interpol_lin(n1d, n1d_mod_ext, log10(r1d), log10(r1d_mod_ext), opalbar1d_lin, log10(opalbar1d_mod_ext))
!call interpol_lin(n1d, n1d_mod_ext, log10(r1d), log10(r1d_mod_ext), gradv1d_lin, log10(gradv1d_mod_ext))
!call interpol_lin(n1d, n1d_mod_ext, log10(r1d), log10(r1d_mod_ext), t1d_lin, log10(t1d_mod_ext))
!!calculate real values since interpolation was performed in log-space
!vel1d_lin=ten**vel1d_lin
!opath1d_lin=ten**opath1d_lin
!opalbar1d_lin=ten**opalbar1d_lin
!gradv1d_lin=ten**gradv1d_lin
!t1d_lin=ten**t1d_lin
!!
!!----------------------perform cubic interpolation----------------------
!!
!call interpol_cube(n1d, n1d_mod_ext, r1d, r1d_mod_ext, vel1d_cube, vel1d_mod_ext)
!call interpol_cube(n1d, n1d_mod_ext, r1d, r1d_mod_ext, opath1d_cube, opath1d_mod_ext)
!call interpol_cube(n1d, n1d_mod_ext, r1d, r1d_mod_ext, opalbar1d_cube, opalbar1d_mod_ext)
!call interpol_cube(n1d, n1d_mod_ext, r1d, r1d_mod_ext, gradv1d_cube, gradv1d_mod_ext)
!call interpol_cube(n1d, n1d_mod_ext, r1d, r1d_mod_ext, t1d_cube, t1d_mod_ext)
!!
!!----------------------perform cubic spline interpolation---------------
!!
!call interpol_spline(n1d, n1d_mod_ext, r1d, r1d_mod_ext, vel1d_spline, vel1d_mod_ext)
!call interpol_spline(n1d, n1d_mod_ext, r1d, r1d_mod_ext, opath1d_spline, opath1d_mod_ext)
!call interpol_spline(n1d, n1d_mod_ext, r1d, r1d_mod_ext, opalbar1d_spline, opalbar1d_mod_ext)
!call interpol_spline(n1d, n1d_mod_ext, r1d, r1d_mod_ext, gradv1d_spline, gradv1d_mod_ext)
!call interpol_spline(n1d, n1d_mod_ext, r1d, r1d_mod_ext, t1d_spline, t1d_mod_ext)
!!
!!
!!----------------------perform cubic interpolation----------------------
!!
!call interpol_cube_mono(n1d, n1d_mod_ext, r1d, r1d_mod_ext, vel1d_cube_mono, vel1d_mod_ext)
!call interpol_cube_mono(n1d, n1d_mod_ext, r1d, r1d_mod_ext, opath1d_cube_mono, opath1d_mod_ext)
!call interpol_cube_mono(n1d, n1d_mod_ext, r1d, r1d_mod_ext, opalbar1d_cube_mono, opalbar1d_mod_ext)
!call interpol_cube_mono(n1d, n1d_mod_ext, r1d, r1d_mod_ext, gradv1d_cube_mono, gradv1d_mod_ext)
!call interpol_cube_mono(n1d, n1d_mod_ext, r1d, r1d_mod_ext, t1d_cube_mono, t1d_mod_ext)
!!
!!---------------------------output to files-----------------------------
!!
!open(1, file=trim(dirout)//'/test_interp1d_model.dat', form='formatted')
!   do i=1, n1d_mod_ext
!      write(1,'(6e20.8)') r1d_mod_ext(i)/sr, vel1d_mod_ext(i), opath1d_mod_ext(i), opalbar1d_mod_ext(i), gradv1d_mod_ext(i), t1d_mod_ext(i)
!   enddo
!close(1)
!!
!open(1, file=trim(dirout)//'/test_interp1d.dat', form='formatted')
!   do i=1, n1d
!      write(1,'(21e20.8)') r1d(i)/sr, vel1d_lin(i), vel1d_cube(i), vel1d_cube_mono(i), vel1d_spline(i), &
!                                     opath1d_lin(i), opath1d_cube(i), opath1d_cube_mono(i), opath1d_spline(i), &
!                                     opalbar1d_lin(i), opalbar1d_cube(i), opalbar1d_cube_mono(i), opalbar1d_spline(i), &
!                                     gradv1d_lin(i), gradv1d_cube(i), gradv1d_cube_mono(i), gradv1d_spline(i), &
!                                     t1d_lin(i), t1d_cube(i), t1d_cube_mono(i), t1d_spline(i)
!   enddo
!close(1)
!
!write(*,*) 'done'
!
!end subroutine test_interpol_1d
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
