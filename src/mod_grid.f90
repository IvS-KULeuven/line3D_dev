!
!-----------------------------------------------------------------------
!------module to define grids, probability density functions, etc.------
!-----------------------------------------------------------------------
!
module mod_grid
!
use prog_type
use fund_const
use mod_sort, only: quicksort  
!
implicit none
!
!
!
type grid1d
   !
   integer(i4b) :: nd
   real(dp), dimension(:), allocatable :: coord, coord_mid, coord_pdf, coord_p
   !
   !nd: number of grid points
   !coord: grid points
   !coord_mid: middle of grid points (dimension nd-1)
   !coord_pdf:   probability density function
   !coord_p:     probability for each interval (only for debugging reasons)
   !
end type grid1d
!
!
!
type grid
   !
   integer(i4b) :: ndim
   type(grid1d) :: xgrid, ygrid, zgrid
   !
end type grid
!
!
!
contains
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
  subroutine grid_equi(xmin, xmax, nd, x)
!
! ... arguments
    integer(i4b), intent(in) :: nd
    real(dp), intent(in) :: xmin, xmax
    real(dp), dimension(nd) :: x
    !
    ! ... local scalars
    integer(i4b) :: i
    !
    do i=1, nd
       x(i) = xmin + (i-1)*(xmax-xmin)/float(nd-1)
    enddo
    !
  end subroutine grid_equi
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
  subroutine grid_log(xmin, xmax, nd, x)
    !   
    ! ... arguments
    integer(i4b), intent(in) :: nd
    real(dp), intent(in) :: xmin, xmax
    real(dp), dimension(nd) :: x
    !
    ! ... local scalars
    integer(i4b) :: i
    real(dp) :: del
    !
    if(xmax/xmin.le.0.d0) stop 'error in grid_log: xmax/xmin le 0.'
    !
    del=log10(xmax/xmin)/float(nd-1)
    x(1)=xmin
    do i=2, nd
       x(i) = x(i-1)*10.d0**del
    enddo
    !
  end subroutine grid_log
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
  subroutine grid_loglog(x1, x2, xmax, nd, x)
!
    ! ... arguments
    integer(i4b), intent(in) :: nd
    real(dp), intent(in) :: x1, x2, xmax
    real(dp), dimension(nd) :: x
    !
    ! ... local scalars
    integer(i4b) :: i
    real(dp) :: del
    !
    if(x2.le.0.d0) stop 'error in grid_log: x2 le 0.'
    if(xmax.le.0.d0) stop 'error in grid_log: xmax le 0.'
    if(log10(xmax)/log10(x2).le.0.d0) stop 'error in grid_log: ration le 0.'
    !
    del=log10(log10(xmax)/log10(x2))/float(nd-2)
    x(1)=x1
    x(2)=x2
    do i=3, nd
       x(i) = 10.d0**(log10(x(i-1))*10.d0**del)
    enddo
    !
  end subroutine grid_loglog
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
  subroutine recalc_grid1d(coord_in, nin, nout, n_free, coord_out)
    !
    !--------recalculates 1d grid from a set of coordinates coord_in--------
    !----------to form a set of coordinates coord_out by averaging----------
    !
    !input:  coord_in: input coordinates, need to be unique, will be sorted
    !        nin:      number of input coordinates
    !        nout:     number of output coordinates
    !        n_free:   number of grid points that shall be kept free
    !                  (special grid points set to zero, or rmax or 1 later on)
    !output: coord_out: output coordinates
    !
    !
    ! ... arguments
    integer(i4b), intent(in) :: nin, nout, n_free
    real(dp), dimension(nin) :: coord_in
    real(dp), dimension(nout) :: coord_out
    !
    ! ... local scalars
    integer(i4b) :: i, j, k
    integer(i4b) :: nav, nrest_in, nrest_out
    real(dp) :: val1, val2, xav
    !
    !
    !sort coord_in
    !call bubblesort(coord_in, nin)
    call quicksort(coord_in, nin, 1, nin)
    !
    !set coord_out to -1 for debug reasons
    coord_out=-1.d0
    !
    !---------------check if coord_in is in fact unique---------------------
    !
    do i=2, nin
       val1=coord_in(i-1)
       val2=coord_in(i)
       if(abs(val2-val1).lt.1.d-14) stop 'error in recalc_grid1d: coord_in not uniq'
    enddo
    !
    !number of used x-coordinates for averaging
    nav=nin/(nout-n_free)
    if(nav.eq.0) stop 'error in recalc_grid1d: nout-n_free gt nin => need larger input grid, or smaller output grid'
    !
    nrest_in=nin
    nrest_out=nout-n_free
    !
    k=1
    do i=1, nout-n_free
       xav=0.d0
       do j=1, nav
          xav=xav+coord_in(k)
          k=k+1
       enddo
       xav=xav/float(nav)
       coord_out(i)=xav
       !
       !adapt nav in order that outer points are all resolved
       nrest_in=nrest_in-nav
       nrest_out=nrest_out-1
       if(nrest_out.ne.0) then
          nav=nrest_in/nrest_out
       endif
       !
    enddo
    !
    !
  end subroutine recalc_grid1d
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
  function frecalc_grid1d(grid_in, nd_out)
    !
    !----------recalculate the grid using nx_out grid points------------
    !----------------------following the pdf----------------------------
    !
    ! ... arguments
    integer(i4b) :: nd_out
    type(grid1d), intent(in) :: grid_in
    type(grid1d) :: grid_out, frecalc_grid1d
    !
    ! ... local scalars
    integer(i4b) :: i, j, err
    real(dp) :: pc, pc_i, xi, delp
    !
    grid_out%nd = nd_out
    !
    !probability that has to be matched for each interval
    pc = one/float(grid_out%nd-1)
    !
    if(allocated(grid_out%coord)) deallocate(grid_out%coord)
    allocate(grid_out%coord(grid_out%nd), stat=err)
    !
    !inner boundary
    grid_out%coord(1) = grid_in%coord(1)
    !
    pc_i = zero
    !
    do i=2, grid_out%nd-1
       !probability that has to be matched when integrating from grid point 1 to grid point i
       pc_i = pc_i + pc
       !
       delp = zero
       do j=1, grid_in%nd - 1
          !calculate a guess for the grid point assuming the pdf of input grid
          xi =  (pc_i-delp)/grid_in%coord_pdf(j)  + grid_in%coord(j)
          grid_out%coord(i) = xi
          if(xi.le.grid_in%coord(j+1)) then !correct position found: go to next x-coordinate
             exit
          else
             !go to next interval and add up previously found probabilities
             delp = delp + grid_in%coord_pdf(j)*(grid_in%coord(j+1)-grid_in%coord(j))
          endif
       enddo
    enddo

    !outer boundary
    grid_out%coord(grid_out%nd) = grid_in%coord(grid_in%nd)
    !
    call calc_pdf2(grid_out)
    !
    frecalc_grid1d = grid_out
    !
  end function frecalc_grid1d
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
  subroutine calc_pdf(nx, x, x_mid, pdf_x, p_x)
    !
    !calculates the probability density function for an input grid, where
    !
    !pdf_x at (x(i)+x(i-1)/2. is 1.d0/(x(i)-x(i-1))/(nx-1)
    !
    !on input: nx   number of coordinates of grid
    !          x    coordinates
    !
    !on output: x_mid     mid-points of coordinates (where pdf lives)
    !           pdf_x     probability density function
    !           p_x       probability of finding x in each interval
    !
    !
    ! ... arguments
    integer(i4b), intent(in) :: nx
    real(dp), dimension(nx), intent(in) :: x
    real(dp), dimension(nx-1), intent(out) :: x_mid, pdf_x, p_x
    !
    ! ... local scalars
    integer(i4b) :: i
    !
    do i=1, nx-1
       x_mid(i) = (x(i+1)+x(i))/2.d0
       pdf_x(i) = 1.d0/(x(i+1)-x(i))/(nx-1)
       p_x(i) = pdf_x(i)*(x(i+1)-x(i))
    enddo
    !
    !normalize pdf (numerical reasons)
    pdf_x=pdf_x/sum(p_x)
    !
  end subroutine calc_pdf
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!    
  subroutine calc_pdf2(grid)
    !
    !calculates the probability density function for an input grid, where
    !same as subroutine calc_pdf, however with the input of type grid
    !
    !pdf_x at (x(i)+x(i-1)/2. is 1.d0/(x(i)-x(i-1))/(nx-1)
    !
    !on input: gridx%nx  number of coordinates of grid
    !          grid%x    coordinates
    !
    !on output: grid%x_mid     mid-points of coordinates (where pdf lives)
    !           grid%pdf_x     probability density function
    !           grid%p_x       probability of finding x in each interval
    !
    ! ... arguments
    type(grid1d), intent(inout) :: grid
    !
    ! ... local scalars
    integer(i4b) :: i, err
    !
    if(allocated(grid%coord_mid)) deallocate(grid%coord_mid)
    if(allocated(grid%coord_pdf)) deallocate(grid%coord_pdf)
    if(allocated(grid%coord_p)) deallocate(grid%coord_p)
    !
    allocate(grid%coord_mid(grid%nd-1), stat=err)
    allocate(grid%coord_pdf(grid%nd-1), stat=err)
    allocate(grid%coord_p(grid%nd-1), stat=err)
    !
    do i=1, grid%nd-1
       grid%coord_mid(i) = (grid%coord(i+1)+grid%coord(i))/two
       grid%coord_pdf(i) = one/(grid%coord(i+1)-grid%coord(i))/float(grid%nd-1)
       grid%coord_p(i) = grid%coord_pdf(i)*(grid%coord(i+1)-grid%coord(i))
    enddo
    !
    !normalize pdf (numerical reasons)
    grid%coord_pdf=grid%coord_pdf/sum(grid%coord_p)
    !
  end subroutine calc_pdf2
!
!-----------------------------------------------------------------------
!
  subroutine print_grid1d(grid)
    !
    ! ... arguments
    type(grid1d) :: grid
    !
    ! ... local scalars
    integer(i4b) :: i
    !
    write(*,*) '-------------GRID-------------'
    write(*,*)
    do i=1, grid%nd-1
       write(*,'(es20.8)') grid%coord(i)
       write(*,'(a20,2es20.8)') '', grid%coord_pdf(i), grid%coord_p(i)
    enddo
    write(*,'(es20.8)') grid%coord(grid%nd)    
    write(*,*)
  end subroutine print_grid1d

end module mod_grid
!
