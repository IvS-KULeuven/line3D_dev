!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine rotate_x (angle, vector, dim)
!
!-------------rotates vector around x by specified angle----------------
!input: angle, vector, dim
!output: vector
!
use prog_type
use fund_const
!
implicit none
!
! ... arguments
integer(i4b) :: dim
real(dp) :: angle
real(dp), dimension(dim) :: vector
!
! ... local arrays
real(dp), dimension(3,3) :: r_x
!
if(dim.ne.3) stop 'wrong dimension in rotate_x'
!
r_x= reshape((/ one,   zero,        zero,     &
                zero, cos(angle), -sin(angle), &
                zero, sin(angle),  cos(angle) /), shape(r_x))
r_x=transpose(r_x)
!
vector=matmul(r_x,vector)
!
end subroutine rotate_x
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine rotate_y (angle, vector, dim)
!
!-------------rotates vector around y by specified angle----------------
!input: angle, vector, dim
!output: vector
!
use prog_type
use fund_const
!
implicit none
!
! ... arguments
integer(i4b) :: dim
real(dp) :: angle
real(dp), dimension(dim) :: vector
!
! ... local arrays
real(dp), dimension(3,3) :: r_y
!
if(dim.ne.3) stop 'wrong dimension in rotate_y'
!
r_y= reshape((/ cos(angle), zero, -sin(angle), &
                   zero,    one,    zero    , &
                sin(angle), zero,  cos(angle) /), shape(r_y))
r_y=transpose(r_y)
!
vector=matmul(r_y,vector)
!
end subroutine rotate_y
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine rotate_z (angle, vector, dim)
!
!-------------rotates vector around z by specified angle----------------
!input: angle, vector, dim
!output: vector
!
use prog_type
use fund_const
!
implicit none
!
! ... arguments
integer(i4b) :: dim
real(dp) :: angle
real(dp), dimension(dim) :: vector
!
! ... local arrays
real(dp), dimension(3,3) :: r_z
!
if(dim.ne.3) stop 'wrong dimension in rotate_z'
!
r_z= reshape((/ cos(angle), -sin(angle), zero, &
                sin(angle),  cos(angle), zero, &
                  zero,        zero,     one /), shape(r_z))
r_z=transpose(r_z)
!
vector=matmul(r_z,vector)
!
end subroutine rotate_z
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine imtql2 ( n, d, e, z, ierr )
!
!*****************************************************************************80
!
!! imtql2 computes all eigenvalues/vectors of a symmetric tridiagonal matrix.
!
!  discussion:
!
!    this subroutine finds the eigenvalues and eigenvectors
!    of a symmetric tridiagonal matrix by the implicit ql method.
!    the eigenvectors of a full symmetric matrix can also
!    be found if tred2 has been used to reduce this
!    full matrix to tridiagonal form.
!
!  licensing:
!
!    this code is distributed under the gnu lgpl license.
!
!  modified:
!
!    18 october 2009
!
!  author:
!
!    original fortran77 version by smith, boyle, dongarra, garbow, ikebe,
!    klema, moler.
!    fortran90 version by john burkardt.
!
!  reference:
!
!    james wilkinson, christian reinsch,
!    handbook for automatic computation,
!    volume ii, linear algebra, part 2,
!    springer, 1971,
!    isbn: 0387054146,
!    lc: qa251.w67.
!
!    brian smith, james boyle, jack dongarra, burton garbow,
!    yasuhiko ikebe, virginia klema, cleve moler,
!    matrix eigensystem routines, eispack guide,
!    lecture notes in computer science, volume 6,
!    springer verlag, 1976,
!    isbn13: 978-3540075462,
!    lc: qa193.m37.
!
!  parameters:
!
!    input, integer ( kind = 4 ) n, the order of the matrix.
!
!    input/output, real ( kind = 8 ) d(n).  on input, the diagonal elements of
!    the input matrix.  on output, the eigenvalues in ascending order.  if an
!    error exit is made, the eigenvalues are correct but
!    unordered for indices 1,2,...,ierr-1.
!
!    input/output, real ( kind = 8 ) e(n).  on input, the subdiagonal elements
!    of the input matrix in e(2:n).  e(1) is arbitrary.  on output, e is
!    overwritten.
!
!    input/output, real ( kind = 8 ) z(n,n).  on input, the transformation
!    matrix produced in the reduction by tred2, if performed.  if the
!    eigenvectors of the tridiagonal matrix are desired, z must contain the
!    identity matrix.  on output, z contains orthonormal eigenvectors of the
!    symmetric tridiagonal (or full) matrix.  if an error exit is made, z
!    contains the eigenvectors associated with the stored eigenvalues.
!
!    output, integer ( kind = 4 ) ierr, error flag.
!    0, for normal return,
!    j, if the j-th eigenvalue has not been determined after 30 iterations.
!
  use prog_type
  use fund_const

  implicit none

  integer ( kind = 4 ) :: n

!  real ( kind = 8 ):: b
!  real ( kind = 8 ):: c
!  real ( kind = 8 ):: d(n)
!  real ( kind = 8 ):: e(n)
!  real ( kind = 8 ):: f
!  real ( kind = 8 ):: g
!  integer ( kind = 4 ):: i
!  integer ( kind = 4 ):: ierr
!  integer ( kind = 4 ):: ii
!  integer ( kind = 4 ):: j
!  integer ( kind = 4 ):: k
!  integer ( kind = 4 ):: l
!  integer ( kind = 4 ):: m
!  integer ( kind = 4 ):: mml
!  real ( kind = 8 ) :: p
!  real ( kind = 8 ) :: pythag
!  real ( kind = 8 ) :: r
!  real ( kind = 8 ) :: s
!  real ( kind = 8 ) :: t(n)
!  real ( kind = 8 ) :: tst1
!  real ( kind = 8 ) :: tst2
!  real ( kind = 8 ) :: z(n,n)

  real ( dp ):: b
  real ( dp ):: c
  real ( dp ):: d(n)
  real ( dp ):: e(n)
  real ( dp ):: f
  real ( dp ):: g
  integer ( i4b ):: i
  integer ( i4b ):: ierr
  integer ( i4b ):: ii
  integer ( i4b ):: j
  integer ( i4b ):: k
  integer ( i4b ):: l
  integer ( i4b ):: m
  integer ( i4b ):: mml
  real ( dp ) :: p
  real ( dp ) :: pythag
  real ( dp ) :: r
  real ( dp ) :: s
  real ( dp ) :: t(n)
  real ( dp ) :: tst1
  real ( dp ) :: tst2
  real ( dp ) :: z(n,n)


  ierr = 0

  if ( n == 1 ) then
    return
  end if

  do i = 2, n
    e(i-1) = e(i)
  end do
  e(n) = zero

  do l = 1, n

    j = 0
!
!  look for a small sub-diagonal element.
!
105 continue

      m = l

      do m = l, n - 1

        tst1 = abs ( d(m) ) + abs ( d(m+1) )
        tst2 = tst1 + abs ( e(m) )

        if ( tst2.eq.tst1 ) then
          exit
        end if

      end do

      p = d(l)

      if ( m == l ) then
        cycle
      end if

      if ( 30 <= j ) then
        ierr = l
        return
      end if

      j = j + 1
!
!  form shift.
!
      g = ( d(l+1) - p ) / ( two * e(l) )
      r = pythag ( g, zero )
      g = d(m) - p + e(l) / ( g + sign ( r, g ) )
      s = one
      c = one
      p = zero
      mml = m - l

      do ii = 1, mml

        i = m - ii
        f = s * e(i)
        b = c * e(i)
        r = pythag ( f, g )
        e(i+1) = r
!
!  recover from underflow.
!
        if (r.eq.zero) then
          d(i+1) = d(i+1) - p
          e(m) = 0.0d+00
          go to 105
        end if

        s = f / r
        c = g / r
        g = d(i+1) - p
        r = ( d(i) - g ) * s + 2.0d+00 * c * b
        p = s * r
        d(i+1) = g + p
        g = c * r - b
!
!  form vector.
!
        do k = 1, n
          f = z(k,i+1)
          z(k,i+1) = s * z(k,i) + c * f
          z(k,i) = c * z(k,i) - s * f
        end do

      end do

      d(l) = d(l) - p
      e(l) = g
      e(m) = 0.0d+00

    go to 105

  end do
!
!  order eigenvalues and eigenvectors.
!
  do ii = 2, n

    i = ii - 1
    k = i
    p = d(i)

    do j = ii, n
      if ( d(j) < p ) then
        k = j
        p = d(j)
      end if
    end do

    if ( k /= i ) then

      d(k) = d(i)
      d(i) = p

      t(1:n)   = z(1:n,i)
      z(1:n,i) = z(1:n,k)
      z(1:n,k) = t(1:n)

    end if

  end do

  return
end subroutine imtql2


function pythag ( a, b )

!*****************************************************************************80
!
!! pythag computes sqrt ( a * a + b * b ) carefully.
!
!  discussion:
!
!    the formula
!
!      pythag = sqrt ( a * a + b * b )
!
!    is reasonably accurate, but can fail if, for example, a^2 is larger
!    than the machine overflow.  the formula can lose most of its accuracy
!    if the sum of the squares is very large or very small.
!
!  licensing:
!
!    this code is distributed under the gnu lgpl license.
!
!  modified:
!
!    18 october 2009
!
!  author:
!
!    original fortran77 version by smith, boyle, dongarra, garbow, ikebe,
!    klema, moler.
!    fortran90 version by john burkardt.
!
!  reference:
!
!    james wilkinson, christian reinsch,
!    handbook for automatic computation,
!    volume ii, linear algebra, part 2,
!    springer, 1971,
!    isbn: 0387054146,
!    lc: qa251.w67.
!
!    brian smith, james boyle, jack dongarra, burton garbow,
!    yasuhiko ikebe, virginia klema, cleve moler,
!    matrix eigensystem routines, eispack guide,
!    lecture notes in computer science, volume 6,
!    springer verlag, 1976,
!    isbn13: 978-3540075462,
!    lc: qa193.m37.
!
!  modified:
!
!    04 february 2003
!
!  parameters:
!
!    input, real ( kind = 8 ) a, b, the two legs of a right triangle.
!
!    output, real ( kind = 8 ) pythag, the length of the hypotenuse.
!
  use prog_type
  
  implicit none

  real ( dp ) a
  real ( dp ) b
  real ( dp ) p
  real ( dp ) pythag
  real ( dp ) r
  real ( dp ) s
  real ( dp ) t
  real ( dp ) u

  p = max ( abs ( a ), abs ( b ) )

  if ( p /= 0.0d+00 ) then

    r = ( min ( abs ( a ), abs ( b ) ) / p )**2

    do

      t = 4.0d+00 + r

      if ( t == 4.0d+00 ) then
        exit
      end if

      s = r / t
      u = 1.0d+00 + 2.0d+00 * s
      p = u * p
      r = ( s / u )**2 * r

    end do

  end if

  pythag = p

  return
end function pythag


!





!
subroutine nonlin(a,b,c, sol)
!
!-----------------------------------------------------------------------
!------------finds nodes of an expression of the form:------------------
!-----------------------x=a+b*exp(c*x)----------------------------------
!------------------------by iteration-----------------------------------
!-----------------------------------------------------------------------
!
use prog_type
!
implicit none
!
! ... arguments
real(dp) :: a,b,c, sol
!
! ... local scalars
integer(i4b) :: i
real(dp) :: xold, xnew
!
!   initial guess
xold=1.d0

do i=1, 1000
! 
   xnew=log(sqrt(((xold-a)/b)**2))/c
!
   if (abs(xold-xnew).lt.1.d-15) then
!      write(*,*) "found value after iteration: ", i
      exit
   endif

   xold=xnew

!
end do

sol=xold

if (i.eq.1001) then
   write(*,*) "error in convergence of iteration scheme"
   write(*,*) "calling regula falsi method"
   call regulafalsi(a,b,c,sol)
!   stop
end if


end subroutine nonlin
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine regulafalsi(a,b,c,sol)
!
!-----------------------------------------------------------------------
!------------finds nodes of an expression of the form:------------------
!-----------------------0=a+b*exp(c*x)-x--------------------------------
!-----------------------------------------------------------------------
!
use prog_type
!
implicit none
!
! ... arguments
real(dp) :: a,b,c, sol
!
! ... local scalars
integer(i4b) :: i
real(dp) :: x0, x1, x2
real(dp) :: y0, y1, y2
!
!   initial guess
x0=0.2d0
x1=0.5d0

do i=1, 1000
!
   y0=a+b*exp(c*x0)-x0
   y1=a+b*exp(c*x1)-x1
!
   if(y0*y1.gt.0.d0) then
      x0=x0-1.d0
      x1=x1+1.d0
   else
      exit
   end if
!
end do
!
y0=a+b*exp(c*x0)-x0
y1=a+b*exp(c*x1)-x1
!
if(y0*y1.gt.0.d0) then
   write(*,*) "error in regula falsi: bad nodes"
   stop
end if
!
do i=1, 1000
!
   y0=a+b*exp(c*x0)-x0
   y1=a+b*exp(c*x1)-x1
!
   x2=(x0*y1-x1*y0)/(y1-y0)
   y2=a+b*exp(c*x2)-x2

!
   if(abs(y2).lt.1.d-14) then
      sol=x2
!      write(*,*) "node found after calculation: ", i
      exit
   end if
!
   if(y2*y0.lt.0) then
      x1=x2
      y1=y2
   else
      x0=x2
      y0=y2
   end if
!
end do

sol=x2

if(i.eq.1001) then
   write(*,*) "error in regula falsi: too many iterations"
   write(*,*) "calling iteration method"
   call nonlin(a,b,c,sol)
end if
!
!write(*,*) "solution: ", sol
!
end subroutine regulafalsi
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine eigenvalues(dim,matrix,evr,evi)
!
use prog_type
!
!
! ... arguments
integer(i4b) :: dim
real(dp), dimension(dim,dim) :: matrix
real(dp), dimension(dim) :: evr, evi
!
! ... local scalars
integer(i4b) :: dummyerr
!
! ... local arrays
real(dp), dimension(dim,dim) :: dummymat, matev
real(dp), dimension(dim) :: dummyvec1, dummyvec2
!
!because matrix is destroyed in procedure rg
matev=matrix
!
call rg(dim,dim, matev, evr, evi, 0, dummymat, dummyvec1, dummyvec2, dummyerr)

!write(*,*) 'real ev', 'imag ev'
!do i=1, dim
!   write(*,*) evr(i), evi(i)
!enddo
!write(*,*)
!stop

!
end subroutine eigenvalues
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine conv_indx_2d_to_1d(indx_x, indx_z, nd_x, indx_1d)
!
!--------calculation of 1-d index corresponding to 2-d indices----------
!------------with: 1.: 1-d index increases by 1 along x-----------------
!------------------2.: 1-d index increases by nd_x along z--------------
!
use prog_type
!
implicit none
!
! ... arguments
integer(i4b), intent(in) :: indx_x, indx_z, nd_x
integer(i4b) :: indx_1d
!
indx_1d = indx_x + nd_x * (indx_z - 1)
!
end subroutine conv_indx_2d_to_1d
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine conv_indx_1d_to_2d (indx_1d, nd_x, indx_x, indx_z)
!
!--------calculation of 2-d indices corresponding to 1-d index----------
!------------with: 1.: 1-d index increases by 1 along x-----------------
!------------------2.: 1-d index increases by nd_x along z--------------
!
use prog_type
!
implicit none
!
! ... arguments
integer(i4b), intent(in) :: indx_1d, nd_x
integer(i4b) :: indx_x, indx_z
!
! ... local scalars
!
indx_x = mod(indx_1d, nd_x)
if(indx_x.eq.0) indx_x=nd_x
indx_z = (indx_1d - indx_x)/nd_x + 1
!
!
end subroutine conv_indx_1d_to_2d
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine conv_indx_3d_to_1d (indx_x, indx_y, indx_z, nd_x, nd_y, indx_1d)
!
!--------calculation of 1-d index corresponding to 3-d indices----------
!------------with: 1.: 1-d index increases by 1 along x-----------------
!------------------2.: 1-d index increases by nd_x along y--------------
!------------------3.: 1-d index increases by nd_x*nd_y along z---------
!
use prog_type
!
implicit none
!
! ... arguments
integer(i4b), intent(in) :: indx_x, indx_y, indx_z, nd_x, nd_y
integer(i4b) :: indx_1d
!
indx_1d = indx_x + nd_x * (indx_y - 1) + nd_x * nd_y * (indx_z - 1)
!
end subroutine conv_indx_3d_to_1d
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine conv_indx_1d_to_3d (indx_1d, nd_x, nd_y, indx_x, indx_y, indx_z)
!
!--------calculation of 3-d indices corresponding to 1-d index----------
!------------with: 1.: 1-d index increases by 1 along x-----------------
!------------------2.: 1-d index increases by nd_x along y--------------
!------------------3.: 1-d index increases by nd_x*nd_y along z---------
!
use prog_type
!
implicit none
!
! ... arguments
integer(i4b), intent(in) :: indx_1d, nd_x, nd_y
integer(i4b) :: indx_x, indx_y, indx_z
!
! ... local scalars
!
indx_x = mod(indx_1d, nd_x)
if(indx_x.eq.0) indx_x=nd_x
indx_y = mod((indx_1d - indx_x)/nd_x , nd_y) + 1
indx_z = ((indx_1d - indx_x)/nd_x - indx_y + 1) / nd_y + 1
!
!
end subroutine conv_indx_1d_to_3d

!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine calc_dev(x_old, x_new, nd, dev_max)
!
!-----------------------------------------------------------------------
!---------calculates max deviation between two iterates/vectors---------
!
!                   dev_max=max((x_old-x_new)/x_new)
!
!   input: x_old, x_new: old/new iterates (vectors)
!          nd: dimension of vectors x_new, x_old
!   output: dev_max: maximum deviation between both vectors
!           x_old is overwritten with the deviation
!
!-----------------------------------------------------------------------
!
use prog_type
!
implicit none
!
! ... arguments
integer(i4b), intent(in) :: nd
real(dp) :: dev_max
real(dp), dimension(nd), intent(in) :: x_new
real(dp), dimension(nd) :: x_old
!
! ... local scalars
integer(i4b) :: err
integer(i4b) :: indx_dev_max
!
! ... local arrays
real(dp), dimension(:), allocatable :: dev
!
!-----------------------------------------------------------------------
!
allocate(dev(nd), stat=err)
   if(err.ne.0) stop 'allocation error: calc_dev'
!
!-----------------------------------------------------------------------
!
!dev=(x_old-x_new)/x_new
!
dev=0.d0
!
where(x_new.ne.0.d0) dev=(x_old-x_new)/x_new
!
dev_max=maxval(abs(dev))
indx_dev_max=maxloc(abs(dev),1)
!
x_old=dev
!
write(*,'(a30, i10, es18.8)') 'max(dev) at 1d-grid-point:', indx_dev_max, dev_max
!
end subroutine calc_dev
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine calc_dev1d(x_old, x_new, imask, nd, dev_max, ix)
!
!
!-----------------------------------------------------------------------
!-----calculates max deviation between two iterates/vectors in 1-d------
!
!                   dev=(x_old-x_new)/x_new
!
!   input: x_old, x_new: old/new iterates (vectors)
!          mask: mask where vectors shall be compared
!          nd_x: dimension of vectors x_new, x_old
!
!   output: dev_max: maximum deviation (scalar)
!           x_old is overwritten with the deviation (1d array)
!           ix: indices where maximum deviation occurred
!           
!
!-----------------------------------------------------------------------
!
use prog_type
!
implicit none
!
! ... arguments
integer(i4b), intent(in) :: nd
integer(i4b), intent(out) :: ix
real(dp) :: dev_max
real(dp), dimension(nd), intent(in) :: x_new
real(dp), dimension(nd) :: x_old
integer, dimension(nd), intent(in) :: imask
!
! ... local scalars
integer(i4b) :: err
real(dp) :: dev_max1, dev_max2
!
! ... local arrays
integer(i4b), dimension(1) :: indx_max
real(dp), dimension(:), allocatable :: dev_1d
logical, dimension(:), allocatable :: lmask
!
!-----------------------allocate local arrays---------------------------
!
allocate(dev_1d(nd), stat=err)
   if(err.ne.0) stop 'allocation error: calc_dev_1d: dev_1d'
allocate(lmask(nd), stat=err)
   if(err.ne.0) stop 'allocation error: calc_dev_1d: lmask'

lmask=.false.
where(imask.eq.1.or.imask.eq.3) lmask=.true.
!
!-----------------------------------------------------------------------
!
dev_1d=0.d0
!
where(x_new.ne.0.d0) dev_1d = (x_new-x_old)/x_new
!
!------calculating maximum percentage-deviation of mean intensities-----
!
dev_max1=maxval(dev_1d,lmask)
dev_max2=minval(dev_1d,lmask)
!
if(abs(dev_max1).gt.abs(dev_max2)) then
   dev_max=dev_max1
   indx_max=maxloc(dev_1d,lmask)
else
   dev_max=dev_max2
   indx_max=minloc(dev_1d,lmask)
endif
!
x_old=dev_1d
!
ix=indx_max(1)
!
end subroutine calc_dev1d
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine calc_dev2d(x_old, x_new, imask, nd_x, nd_z, dev_max, ix, iz)
!
!
!-----------------------------------------------------------------------
!-----calculates max deviation between two iterates/vectors in 2-d------
!
!                   dev=(x_old-x_new)/x_new
!
!   input: x_old, x_new: old/new iterates (vectors)
!          mask: mask where vectors shall be compared
!          nd_x, nd_z: dimension of vectors x_new, x_old in x,z-direction
!
!   output: dev_max: maximum deviation (scalar)
!           x_old is overwritten with the deviation (3d array)
!           ix, iz: indices where maximum deviation occurred
!           
!
!-----------------------------------------------------------------------
!
use prog_type
!
implicit none
!
! ... arguments
integer(i4b), intent(in) :: nd_x, nd_z
integer(i4b), intent(out) :: ix, iz
real(dp) :: dev_max
real(dp), dimension(nd_x,nd_z), intent(in) :: x_new
real(dp), dimension(nd_x,nd_z) :: x_old
integer, dimension(nd_x,nd_z), intent(in) :: imask
!
! ... local scalars
integer(i4b) :: err
real(dp) :: dev_max1, dev_max2
!
! ... local arrays
integer(i4b), dimension(2) :: indx_max
real(dp), dimension(:,:), allocatable :: dev_2d
logical, dimension(:,:), allocatable :: lmask
!
!-----------------------allocate local arrays---------------------------
!
allocate(dev_2d(nd_x,nd_z), stat=err)
   if(err.ne.0) stop 'allocation error: calc_dev2d: dev_2d'
allocate(lmask(nd_x,nd_z), stat=err)
   if(err.ne.0) stop 'allocation error: calc_dev2d: lmask'

lmask=.true.
where(imask.eq.0) lmask=.false.
!
!-----------------------------------------------------------------------
!
dev_2d=0.d0
!
where(x_new.ne.0.d0) dev_2d = (x_new-x_old)/x_new
!
!------calculating maximum percentage-deviation of mean intensities-----
!
dev_max1=maxval(dev_2d,lmask)
dev_max2=minval(dev_2d,lmask)
!
if(abs(dev_max1).gt.abs(dev_max2)) then
   dev_max=dev_max1
   indx_max=maxloc(dev_2d,lmask)
else
   dev_max=dev_max2
   indx_max=minloc(dev_2d,lmask)
endif
!
x_old=dev_2d
!
ix=indx_max(1)
iz=indx_max(2)
!
end subroutine calc_dev2d
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine calc_dev3d(x_old, x_new, imask, nd_x, nd_y, nd_z, dev_max, ix, iy, iz)
!
!
!-----------------------------------------------------------------------
!-----calculates max deviation between two iterates/vectors in 3-d------
!
!                   dev=(x_old-x_new)/x_new
!
!   input: x_old, x_new: old/new iterates (vectors)
!          mask: mask where vectors shall be compared
!          nd_x, nd_y, nd_z: dimension of vectors x_new, x_old in x,y,z-direction
!
!   output: dev_max: maximum deviation (scalar)
!           x_old is overwritten with the deviation (3d array)
!           ix, iy, iz: indices where maximum deviation occurred
!           
!
!-----------------------------------------------------------------------
!
use prog_type
!
implicit none
!
! ... arguments
integer(i4b), intent(in) :: nd_x, nd_y, nd_z
integer(i4b), intent(out) :: ix, iy, iz
real(dp) :: dev_max
real(dp), dimension(nd_x,nd_y,nd_z), intent(in) :: x_new
real(dp), dimension(nd_x,nd_y,nd_z) :: x_old
integer(i1b), dimension(nd_x,nd_y,nd_z), intent(in) :: imask
!
! ... local scalars
integer(i4b) :: err
real(dp) :: dev_max1, dev_max2
!
! ... local arrays
integer(i4b), dimension(3) :: indx_max
real(dp), dimension(:,:,:), allocatable :: dev_3d
logical, dimension(:,:,:), allocatable :: lmask
!
!-----------------------allocate local arrays---------------------------
!
allocate(dev_3d(nd_x,nd_y,nd_z), stat=err)
   if(err.ne.0) stop 'allocation error: calc_dev_3d: dev_3d'
allocate(lmask(nd_x,nd_y,nd_z), stat=err)
   if(err.ne.0) stop 'allocation error: calc_dev3d: lmask'
!
lmask=.true.
where(imask.eq.0) lmask=.false.
!
!-----------------------------------------------------------------------
!
dev_3d=0.d0
!
where(x_new.ne.0.d0) dev_3d = (x_new-x_old)/x_new
!
!------calculating maximum percentage-deviation of mean intensities-----
!
dev_max1=maxval(dev_3d,lmask)
dev_max2=minval(dev_3d,lmask)
!
if(abs(dev_max1).gt.abs(dev_max2)) then
   dev_max=dev_max1
   indx_max=maxloc(dev_3d,lmask)
else
   dev_max=dev_max2
   indx_max=minloc(dev_3d,lmask)
endif
!
x_old=dev_3d
!
ix=indx_max(1)
iy=indx_max(2)
iz=indx_max(3)
!
end subroutine calc_dev3d
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine calc_mean3d(nd_x, nd_y, nd_z, imask3d, x3d, meanx)
!
!
!-----------------------------------------------------------------------
!-----calculates mean value of a 3d array x3d within masked region------
!
!   input: x3d: array of values for which mean shall be calculated
!          imask: mask where values shall be taken
!          nd_x, nd_y, nd_z: dimension of x3d, imask3d
!
!   output: meanx: mean value of x3d at imask3d
!           
!-----------------------------------------------------------------------
!
use prog_type
use fund_const, only: zero
!
implicit none
!
! ... arguments
integer(i4b), intent(in) :: nd_x, nd_y, nd_z
integer(i1b), dimension(nd_x,nd_y,nd_z), intent(in) :: imask3d
real(dp), dimension(nd_x,nd_y,nd_z), intent(in) :: x3d
real(dp), intent(out) :: meanx
!
! ... local scalars
integer(i4b) :: i, j, k, kk
!
! ... local arrays
!
!-----------------------allocate local arrays---------------------------
!
!-----------------------------------------------------------------------
!
meanx = zero
kk = 0
!
do i=1, nd_x
   do j=1, nd_y
      do k=1, nd_z
         if(imask3d(i,j,k).ne.0) then
            meanx = meanx + x3d(i,j,k)
            kk = kk + 1
         endif
      enddo
   enddo
enddo
!
meanx = meanx/float(kk)
!
!
end subroutine calc_mean3d
!
!-----------------------------------------------------------------------
!------------function to calculate r (xcoord, ycoord, zcoord)-----------
!-----------------------------------------------------------------------
!
function getr(xcoord, ycoord, zcoord)
!
use prog_type
!
implicit none
!
! ... arguments
real(dp), intent(in) :: xcoord, ycoord, zcoord
!
! ... local scalars
real(dp) :: getr
!
getr = sqrt(xcoord*xcoord + ycoord*ycoord + zcoord*zcoord)
!
return
end function getr
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
function dist_sphere(nx,ny,nz,xp,yp,zp)
!
!calculate distance of the intersection point
! of a line with direction n_x, n_y, n_z through point xp, yp, zp
!   with a sphere
!input: nx, ny, nz, xp, yp, zp
!output: dist_sphere:    calculated distance
!
!
use prog_type
use fund_const
!
implicit none
!
! ... arguments
real(dp), intent(in) :: nx, ny, nz, xp, yp, zp
real(dp) :: dist_sphere
!
! ... local scalars
real(dp) :: disc, dist1, dist2
!
!calculate discriminant
disc=(xp*nx+yp*ny+zp*nz)**2 - (xp**2+yp**2+zp**2-1.d0)
if(disc.ge.zero) then
!   dist1=-xp*nx-yp*ny-zp*nz+sqrt(disc)
!   dist2=-xp*nx-yp*ny-zp*nz-sqrt(disc)
   dist1=xp*nx+yp*ny+zp*nz+sqrt(disc)
   dist2=xp*nx+yp*ny+zp*nz-sqrt(disc)
   dist_sphere=min(dist1,dist2)
!if dist_sphere has negative values, set it to an arbitrary large value (dist_sphere not required)
   if(dist_sphere.lt.zero) dist_sphere=1.d10
!   dist_sphere=sign(1.d0,dist1)*min(abs(dist1),abs(dist2))
else
   dist_sphere=1.d10
endif
!
end function dist_sphere
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
function dist_ellipsoid(nx,ny,nz,xp,yp,zp,sma,smb,smc)
!
!calculate distance of the intersection point
! of a line with direction n_x, n_y, n_z through point xp, yp, zp
!   with an ellipsoid specified by semi-major-axes sma, smb, smc
!input: nx, ny, nz, xp, yp, zp, sma, smb, smc
!output: dist_sphere:    calculated distance
!
!
use prog_type
use fund_const  
!
implicit none
!
! ... arguments
real(dp), intent(in) :: nx, ny, nz, xp, yp, zp, sma, smb, smc
real(dp) :: dist_ellipsoid
!
! ... local scalars
real(dp) :: disc, dist1, dist2
real(dp) :: acoeff, bcoeff, ccoeff
!
!calculate discriminant
acoeff = (nx/sma)**2 + (ny/smb)**2 + (nz/smc)**2
bcoeff = xp*nx/sma**2 + yp*ny/smb**2 + zp*nz/smc**2
ccoeff = (xp/sma)**2 + (yp/smb)**2 + (zp/smc)**2 - 1.d0
!
disc=bcoeff**2 - acoeff*ccoeff
if(disc.ge.zero) then
   dist1 = (bcoeff + sqrt(disc))/acoeff
   dist2 = (bcoeff - sqrt(disc))/acoeff
   dist_ellipsoid=min(dist1,dist2)
!if dist_sphere has negative values, set it to an arbitrary large value (dist_sphere not required)
   if(dist_ellipsoid.lt.zero) dist_ellipsoid=1.d10
else
   dist_ellipsoid=1.d10
endif
!
end function dist_ellipsoid
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine cross_product(a, b, c)
!
!------------calculates the cross-product of two vectors----------------
!--------------------------a x b = c------------------------------------
!
use prog_type
!
implicit none
!
! ... arguments
real(dp), dimension(3), intent(in) :: a, b
real(dp), dimension(3) :: c
!
c(1)=a(2)*b(3)-a(3)*b(2)
c(2)=a(3)*b(1)-a(1)*b(3)
c(3)=a(1)*b(2)-a(2)*b(1)
!
end subroutine cross_product
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine invtri_lev(nd, a_vec, b_vec, c_vec, x_vec, y_vec)
!
!    calculates solution of linear tridiagonal system
!             a*x_vec = y_vec
!
!    with: sub-diag(a) = a_vec
!          super-diag(a) = c_vec
!          diag(a) = b_vec
!
!    note: a_vec, b_vec, c_vec are overwritten
!
!    using the thomas-algorithm
!
use prog_type
!
implicit none
!
! ... arguments
integer(i4b), intent(in) :: nd
real(dp), dimension(nd) :: a_vec, b_vec, c_vec
real(dp), dimension(nd) :: x_vec, y_vec
!
! ... local scalars
integer(i4b) :: i
real(dp) :: dum
!
!
do i=2, nd
   dum = a_vec(i)/b_vec(i-1)
   b_vec(i) = b_vec(i) - dum*c_vec(i-1)
   y_vec(i) = y_vec(i) - dum*y_vec(i-1)
enddo
!
!
!
x_vec(nd) = y_vec(nd)/b_vec(nd)
do i=nd-1, 1, -1
   x_vec(i) = (y_vec(i) - c_vec(i)*x_vec(i+1)) / b_vec(i)
enddo
!
!
!
end subroutine invtri_lev
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine get_angles_spc(x, y, z, theta, phi)
!
!--------------calculates angles in spherical coordinates---------------
!---------------for given carthesian coordinates x, y, z----------------
!----------------------for right-handed system--------------------------
!
!on input: coordinates x, y, z
!
!on output: angles theta, phi
!
!(from jon's formal solver)
!
!-----------------------------------------------------------------------
!
use prog_type
use fund_const, only: pi, zero, two, three
!
implicit none
!
! ... arguments
real(dp), intent(in) :: x, y, z
real(dp) :: theta, phi
!
! ... local scalars
real(dp) :: rad, x_test, y_test, z_test
!
! ... local functions
!
!-----------------------calculation of theta----------------------------
!
rad = sqrt(x*x+y*y+z*z)
!
if(rad.ne.zero) then
   theta = acos(z/rad)
else
   theta = zero
endif
!
!------------------------calculation of phi-----------------------------
!
if(z.eq.rad) then
   phi=zero
   return
endif
!
if(x.eq.zero) then
   if(y.eq.zero) then
      phi = zero
   else if(y.gt.zero) then
      phi = pi/two
   else
      phi = three*pi/two
   endif
else
   if(y.gt.zero.and.x.gt.zero) then
!first quadrant
      phi = atan(y/x)
   else if (x.gt.zero) then
!fourth quadrant
      phi = two*pi + atan(y/x)
   else
!second and third quadrant
      phi = pi + atan(y/x)
   endif
endif
!
!-------------------------test if everything worked---------------------
!
x_test=rad*sin(theta)*cos(phi)
y_test=rad*sin(theta)*sin(phi)
z_test=rad*cos(theta)
!
if(abs(x_test-x).gt.1.d-6) then
   write(*,*) 'error in get_angles_spc: x', x_test, x, x_test-x
   stop
endif
!
if(abs(y_test-y).gt.1.d-6) then
   write(*,*) 'error in get_angles_spc: y', y_test, y, y_test-y
   stop
endif
!
if(abs(z_test-z).gt.1.d-6) then
   write(*,*) 'error in get_angles_spc: z', z_test, z, z_test-z
   stop
endif
!
!
!
end subroutine get_angles_spc
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine transform_sigr_sigm(xr, yr, zr, xm, ym, zm, obliquity)
!
!   transforms carthesian coordinates (xr, yr, zr), in which z-axis is
!                   aligned with rotational axis
!                              to
!              carthesian coordinates (xm, ym, zm), in which z-axis is
!                   aligned with magnetic pole-axis
!
!   this is only a rotation around y-axis using the obliquity
!
use prog_type
!
implicit none
!
! ... arguments
real(dp), intent(in) :: xr, yr, zr, obliquity
real(dp) :: xm, ym, zm
!
! ... local scalars
real(dp) :: sinb, cosb
!
sinb=sin(obliquity)
cosb=cos(obliquity)
!
xm = cosb*xr - sinb*zr
ym = yr
zm = sinb*xr + cosb*zr
!
!
!
end subroutine transform_sigr_sigm
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine transform_sigm_sigr(xm, ym, zm, xr, yr, zr, obliquity)
!
!   transforms carthesian coordinates (xm, ym, zm), in which z-axis is
!                   aligned with magnetic pole-axis
!                              to
!              carthesian coordinates (xr, yr, zr), in which z-axis is
!                   aligned with rotational axis
!
!   this is only a rotation around y-axis using the obliquity
!
use prog_type
!
implicit none
!
! ... arguments
real(dp), intent(in) :: xm, ym, zm, obliquity
real(dp) :: xr, yr, zr
!
! ... local scalars
real(dp) :: sinb, cosb
!
sinb=sin(obliquity)
cosb=cos(obliquity)
!
xr = cosb*xm + sinb*zm
yr = ym
zr = -sinb*xm + cosb*zm
!
!
!
end subroutine transform_sigm_sigr
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine calc_phinorm(vproj, vth_local, vth_fiducial, xobs, phinorm)
!
!-------------calculates normalized profile function--------------------
!
! INPUT:   vproj:  projected line of sight velocity in vth_fiducial
!          vth_local: local thermal velocity (including micro-turbulent velocity) in cm/s
!          xobs:   frequency shift from line center in fiducial doppler units
!
! OUTPUT:  phinorm: normalized profile
!
use prog_type
use fund_const
!
implicit none
!
! ... arguments
real(dp), intent(in) :: vproj, vth_local, vth_fiducial, xobs
real(dp) :: phinorm
!
! ... local scalars
real(dp) :: delta, xcmf
!
!
delta=vth_local/vth_fiducial
!
xcmf = (xobs - vproj)/delta
!
!if(abs(xcmf).gt.5.d0) then
!   phinorm=0.d0
!else
   phinorm=exp(-xcmf*xcmf)/delta/spi
!endif
!
!
end subroutine calc_phinorm
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine calc_dtau_crit(dtu, dtd_a, dtd_b, dtd_c)
!
!calculate critical downwind delta-tau steps for a given upwind delta-tau step
!
use prog_type
!
implicit none
!
! ... arguments
real(dp), intent(in) :: dtu
real(dp), intent(out) :: dtd_a, dtd_b, dtd_c
!
! ... local scalars
integer(i4b), parameter :: nmax_a=100000, nmax_b=1000, nmax_c=100000
real(dp), parameter :: tol_a=1.d-4, tol_b=1.d-4, tol_c=1.d-4
integer(i4b) :: i
real(dp) :: eps, dtd_new, fim, grad, fi1, fi2, dtd1, dtd2, f0, a, b, c
!
!-----------------------------------------------------------------------
!
!find critical dtau_d with Newton-Method (1/2*(b_plus + b_minus) > 1, solution b)
!
if(dtu.lt.1.d-3) then
   a=0.5d0*dtu
   b=exp(-dtu)-dtu-1.d0
   c=2.d0*exp(-dtu)-2.d0+dtu+dtu*exp(-dtu)
   dtd_b=(-b-sqrt(b**2 - 4*a*c))/2.d0/a
else
!initial guess
   dtd_b=0.d0
   fim=(dtu+dtd_b+2.d0)*(exp(-dtu)+exp(-dtd_b))-4.d0

   do i=1, nmax_b
!calculate zero by newton method
      grad = exp(-dtu)+exp(-dtd_b) - (2.d0+dtu+dtd_b)*exp(-dtd_b)
      dtd_b = dtd_b - fim /grad
      fim = (dtu+dtd_b+2.d0)*(exp(-dtu)+exp(-dtd_b))-4.d0
      if(abs(fim).lt.tol_b) exit
   enddo
!
   if(i.gt.nmax_b) stop 'error in calc_dtau_crit: no convergence in case b'
endif
!
!-----------------------------------------------------------------------
!
!find critical dtau_d by iteration (1/2*(b_plus + b_minus) > 1, solution a)
!
!initial guess
dtd_a=dtu
if(dtu.lt.1.d-3) then
   dtd_a=1.d-2
elseif(dtu.lt.1.d-2) then
   dtd_a=1.d-1
elseif(dtu.lt.1.d-1) then
   dtd_a=1.d0
else
   dtd_a=dtu
endif

if(dtu.gt.230) then
!solution becomes infinite
   dtd_a=1.d100
else
!solution by iteration
   do i=1, nmax_a
      dtd_new = 4.d0/(exp(-dtu)+exp(-dtd_a))-dtu-2.d0
      eps=abs(dtd_a-dtd_new)/abs(dtd_a)
      dtd_a=dtd_new
      if(eps.lt.tol_a) exit
!      write(*,*) i, dtu, dtd_a, eps
   enddo
   if(i.gt.nmax_a) stop 'error in calc_dtau_crit: no convergence in case a'
endif
!
!-----------------------------------------------------------------------
!
!find critical dtau_d with bisection (a_plus+a_minus+c_plus+c_minus < 0)
dtd1=1.d-8
fi1=4.d0-(dtu+dtd1+dtu*dtd1+2.d0)*(exp(-dtu) + exp(-dtd1))
!
do i=1, 1000
   dtd2=-8.d0+(i-1)*(10.d2)/999.d0
   dtd2=10.d0**dtd2
   fi2=4.d0-(dtu+dtd2+dtu*dtd2+2.d0)*(exp(-dtu) + exp(-dtd2))
   if(fi1*fi2.lt.0.d0) exit
enddo
!if(dtu.gt.1.d2) dtd2=4.d0*dtd_b
!
if(fi1*fi2.gt.0.d0) stop 'error in calc_dtau_crit: initial guess is bad case c'
!
if(dtu.gt.1.d5) then
!solution by iteration
   dtd_c=15.d0
else
!
   do i=1, nmax_c
!
!calculate zero by bisection method
      dtd_c = dtd1 - fi1*(dtd2-dtd1)/(fi2-fi1)
      f0 = 4.d0-(dtu+dtd_c+dtu*dtd_c+2.d0)*(exp(-dtu) + exp(-dtd_c))
      if(abs(f0).lt.tol_c) exit
      if(f0*fi1.gt.0.d0) then 
         dtd1=dtd_c
         fi1=f0
      else
         dtd2=dtd_c
         fi2=f0
      endif
   enddo
!
   if(i.gt.nmax_c) stop 'error in calc_dtau_crit: no convergence in case c'
!
endif
!
end subroutine calc_dtau_crit
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine smooth_gauss(nd, x, y, sigma, ngauss, imask)
!
!------smoothes a function y(x) over sigma with gaussian in x space-----
!--------------NOTE: works only for equidistant x array-----------------
!
! input:   x:   x-coordinates of funcion y(x)
!          y:   y(x)
!          sigma:  width of smoothing function (dimension nd)
!          ngauss: number of datapoints to be smoothed over (dimension nd)
!          imask: mask where smoothing shall be performed (dimension nd, imsask(i)=1 if smoothing performed)  
!          nd:  number of datapoints
!
! output:  y:   smoothed function (input y is overwritten)
!
use prog_type
use fund_const
!
implicit none
!
! ... arguments
integer(i4b), intent(in) :: nd
integer(i4b), dimension(nd), intent(in) :: ngauss, imask
real(dp), dimension(nd), intent(in) :: x, sigma
real(dp), dimension(nd), intent(inout) :: y
!
! ... local scalars
integer(i4b) :: i, j, jstart, jend
real(dp) :: sum, norm, gauss, dx, dx_ref
!
! ... local arrays
real(dp), dimension(nd) :: ysmooth
!
!
!perform smoothing
do i=1, nd
   if(imask(i).eq.1) then
      if(modulo(ngauss(i),2).ne.1) stop 'error in smooth_gauss: ngauss has to be odd'
      if(ngauss(i).lt.3) stop  'error in smooth_gauss: ngauss has to be larger equals 3'
      
      jstart = maxval((/ 1, i-(ngauss(i)-1)/2 /))
      jend = minval((/ nd, i+(ngauss(i)-1)/2 /))

      sum=zero
      norm=zero
      do j=jstart,jend
         gauss = exp(-((x(i)-x(j))/sigma(i))**2)
         sum = sum + y(j)*gauss
         norm = norm + gauss
      enddo
      ysmooth(i) = sum/norm
   else
      ysmooth(i) = y(i)
   endif
enddo
!
y=ysmooth

end subroutine smooth_gauss
