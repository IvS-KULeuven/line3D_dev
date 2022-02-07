!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
module mod_sort

use prog_type
implicit none
!
contains
!
!-----------------------------------------------------------------------
!
  subroutine bubblesort(array, n)
    !
    ! ... arguments
    integer(i4b) :: n
    real(dp), dimension(n) :: array
    !
    ! ... local scalars
    integer(i4b) :: i, j
    real(dp) :: temp1, temp2
    real(dp), dimension(n) :: array1
    !
    do j=1, n
       do i=1, n-1
          if(array(i+1).gt.array(i)) then
             temp1=array(i)
             temp2=array(i+1)
             array(i)=temp2
             array(i+1)=temp1
          end if
       end do
    end do
    !
    !------------rearranging such that array-values are increasing----------
    !
    do i=1, n
       array1(i)=array(n+1-i)
    end do
    !
    array=array1
    !
    !
  end subroutine bubblesort
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
  !
  recursive subroutine quicksort(arr, nd, istart_subarr, iend_subarr)
    !
    !from
    !
    ! quicksort.f -*-f90-*-
    ! author: t-nissie
    ! license: gplv3
    ! gist: https://gist.github.com/t-nissie/479f0f16966925fa29ea
    !
    !but variable names somewhat changed for more understanding :)
    !
    ! ... arguments
    integer(i4b), intent(in) :: nd
    real(dp), dimension(nd) :: arr
    integer(i4b) :: istart_subarr, iend_subarr
    !
    !local scalars
    integer(i4b) :: ipivot, ihi, ilo
    real(dp) :: pivot, temp
    !
    !
    !take middle element as pivot
    ipivot = (istart_subarr + iend_subarr)/2
    pivot= arr(ipivot)
    !
    !initialize index-variables
    ihi = istart_subarr
    ilo = iend_subarr
    !
    !infinite loop, until counter-variables cross
    do
       !search from left for elements which are higher than the pivot
       do while (arr(ihi).lt.pivot)
          ihi=ihi+1
       enddo
       !search from right for elements which are less than the pivot
       do while (arr(ilo).gt.pivot)
          ilo=ilo-1
       enddo
       !
       !exit outer loop if indices cross
       if(ihi.ge.ilo) exit
       !
       !swap elements if this is not the case
       temp = arr(ihi)
       arr(ihi) = arr(ilo)
       arr(ilo) = temp
       !and go one step on
       ihi=ihi+1
       ilo=ilo-1
       !
    enddo
    !
    if (istart_subarr.lt.ihi-1) call quicksort(arr, nd, istart_subarr, ihi-1)
    if (ilo+1.lt.iend_subarr) call quicksort(arr, nd, ilo+1, iend_subarr)
    !
    !
    !
  end subroutine quicksort
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
  subroutine uniq_elements(nd, nd_uniq, array, array_uniq)
    !
    !   finds unique elements of an array and stores them into
    !           array_uniq
    !note: input array has to be sorted
    !
    !
    ! ... arguments
    integer(i4b), intent(in) :: nd
    integer(i4b), intent(out) :: nd_uniq
    real(dp), dimension(nd) :: array
    real(dp), dimension(:), allocatable :: array_uniq
    !
    ! ... local scalars
    integer(i4b) :: i, err, indx
    real(dp) :: val1, val2
    !
    ! ... local arrays
    real(dp), dimension(nd) :: array_tmp
    !
    !
    indx=1
    val1=array(1)
    array_tmp(indx)=val1
    !
    do i=2, nd
       val2=array(i)
       if(val1.ne.val2) then
          indx=indx+1
          array_tmp(indx)=val2
       endif
       val1=val2
    enddo
    !
    !
    nd_uniq=indx
    !
    if(allocated(array_uniq)) deallocate(array_uniq)
    allocate(array_uniq(nd_uniq), stat=err)
    !
    array_uniq = array_tmp(1:nd_uniq)
    !
  end subroutine uniq_elements
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
  subroutine uniq_nvec(n_tot, dim_n_tot, dim_unique)
    !
    !--------------finds unique n-vectors for array n_tot ------------------
    !-----------(all rotated n-vectros are stored in n_tot)-----------------
    !
    !
    ! ... arguments
    integer(i4b) :: dim_n_tot, dim_unique
    real(dp), dimension(3,dim_n_tot) :: n_tot
    !
    ! ... local scalars
    integer(i4b) :: i, j
    !
    ! ... local functions
    !
    ! ... local arrays
    real(dp), dimension(3) :: n_test1, n_test2
    real(dp), dimension(3), parameter :: zero_vec=(/ 0.d0, 0.d0, 0.d0 /)
    !
    !
    do i=1, dim_n_tot
       n_test1=n_tot(:,i)
       do j=i+1, dim_n_tot
          n_test2=n_tot(:,j)
          if(arr_equal(n_test1,n_test2,3)) then
             n_tot(:,j)=0.d0
          endif
       enddo
    enddo
    !
    dim_unique=0
    !
    do i=1, dim_n_tot
       if(.not.arr_equal(n_tot(:,i),zero_vec,3)) then
          dim_unique=dim_unique+1
       endif
    enddo
    !
    !
  end subroutine uniq_nvec
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
  subroutine uniq_arr(array, dim_arr, dim_arr_unique)
    !
    !---------------finds unique elements of a sorrted array----------------
    !
    !
    ! ... arguments
    integer(i4b), intent(in) :: dim_arr
    integer(i4b), intent(out) :: dim_arr_unique
    real(dp), dimension(dim_arr), intent(in) :: array
    !
    ! ... local scalars
    integer(i4b) :: i
    real(dp) :: val1, val2
    !
    ! ... local functions
    !
    ! ... local arrays
    !
    !
    dim_arr_unique=1
    !
    val1=array(1)
    !
    do i=2, dim_arr
       val2=array(i)
       if(val1.ne.val2) then
          dim_arr_unique=dim_arr_unique+1
       endif
       val1=val2
    enddo
    !
    !
  end subroutine uniq_arr
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
  function arr_equal(array1, array2, dim)
    !
    !----------------check if array1 and array2 are equal-------------------
    !
    !
    ! ... arguments
    logical :: arr_equal
    integer(i4b) :: dim
    real(dp), dimension(dim) :: array1, array2
    !
    ! ... local scalars
    integer(i4b) :: i
    !
    !default
    arr_equal=.false.
    !
    !check if all elements are equal
    do i=1, size(array1)
       if(array1(i).eq.array2(i)) then
          arr_equal=.true.
       else
          arr_equal=.false.
          exit
       endif
    enddo
    !
    !
  end function arr_equal
!
!  
!
end module mod_sort  
  
