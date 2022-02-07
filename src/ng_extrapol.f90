!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine store_ng2d(iit, viit, ndxmax, ndzmax, source2d)
!
!--------stores subsequent iterates for ng/aitkens-extrapolation--------
!
! input: iit:   iteration step
!        ndxmax, ndzmax:   dimensions of 2d array
!        source2d:   2d source function
!
! output: viit: 3d source function at current iteration step iit as 1d array
!
use prog_type
!
implicit none
!
! ... arguments
integer(i4b), intent(in) :: iit, ndxmax, ndzmax
real(dp), dimension(ndxmax, ndzmax), intent(in) :: source2d
real(dp), dimension(4,ndxmax*ndzmax) :: viit
!
! ... local scalars
integer(i4b) :: i,k
integer(i4b) :: indx_1d
!
!-----------------------------------------------------------------------
!
select case(iit)
   case(1) 
      write(*,*) '----------------------storing source fct. at step n-3--------------------------'
      write(*,*)
   case(2) 
      write(*,*) '----------------------storing source fct. at step n-2--------------------------'
      write(*,*)
   case(3) 
      write(*,*) '----------------------storing source fct. at step n-1--------------------------'
      write(*,*)
   case(4) 
      write(*,*) '----------------------storing source fct. at step n----------------------------'
      write(*,*)
   case default
end select
!
!
!
do i=1, ndxmax
   do k=1, ndzmax
      call conv_indx_2d_to_1d (i, k, ndxmax, indx_1d)
      viit(iit,indx_1d) = source2d(i,k)
   enddo
enddo
!
end subroutine store_ng2d
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine ng_expol2d(viit, ndxmax, ndzmax, source2d)
!
use prog_type
!
implicit none
!
! ... arguments
integer(i4b), intent(in) :: ndxmax, ndzmax
real(dp), dimension(ndxmax,ndzmax) :: source2d
real(dp), dimension(4,ndxmax*ndzmax) :: viit
!
! ... local scalars
integer(i4b) :: i
integer(i4b) :: indx_x, indx_z
real(dp) :: a1, a2, b1, b2, c1, c2
real(dp) :: ap, bp, cp
!
! ... local arrays
!
! ... local logicals
logical :: negative
!
!-----------------------------------------------------------------------
!
write(*,*) '-------------------------calculating new iterate (ng)--------------------------'
write(*,*)
!
!
a1=0.d0
a2=0.d0
b1=0.d0
b2=0.d0
c1=0.d0
c2=0.d0
!
do i=1, ndxmax*ndzmax
   a1 = a1+(viit(4,i)-2.d0*viit(3,i)+viit(2,i))**2.d0
   b1 = b1+(viit(4,i)-viit(3,i)-viit(2,i)+viit(1,i))*(viit(4,i)-2.d0*viit(3,i)+viit(2,i))
   c1 = c1+(viit(4,i)-viit(3,i))*(viit(4,i)-2.d0*viit(3,i)+viit(2,i))
   a2 = a2+(viit(4,i)-viit(3,i)-viit(2,i)+viit(1,i))*(viit(4,i)-2.d0*viit(3,i)+viit(2,i))
   b2 = b2+(viit(4,i)-viit(3,i)-viit(2,i)+viit(1,i))**2.d0
   c2 = c2+(viit(4,i)-viit(3,i))*(viit(4,i)-viit(3,i)-viit(2,i)+viit(1,i))
end do
!
ap=(b2*c1-b1*c2)/(a1*b2-a2*b1)
bp=(a1*c2-a2*c1)/(a1*b2-a2*b1)

!extrapolation can yield to negative source-functions
!        if so, then skip extrapolation
negative=.false.

do i=1, ndxmax*ndzmax
   cp=(1.d0-ap-bp)*viit(4,i) + ap*viit(3,i)+bp*viit(2,i)
   if(cp.lt.0.d0) then
      write(*,*) 'ng-extrapolation yields negative source-functions => skip'
      negative=.true.
      exit
    endif
enddo
!
if(.not.negative) then
   viit(1,:)=(1.d0-ap-bp)*viit(4,:) + ap*viit(3,:)+bp*viit(2,:)
else
   viit(1,:) = viit(4,:)
endif
!
!--------------------back conversion on 2-d-array-----------------------
!
!source2d=0.d0
do i=1, ndxmax*ndzmax
   call conv_indx_1d_to_2d(i, ndxmax, indx_x, indx_z)
   source2d(indx_x,indx_z)=viit(1,i)
enddo
!
!
end subroutine ng_expol2d
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine ait_expol2d(viit,ndxmax, ndzmax, source2d)
!
use prog_type
!
implicit none
!
! ... arguments
integer(i4b), intent(in) :: ndxmax, ndzmax
real(dp), dimension(ndxmax,ndzmax) :: source2d
real(dp), dimension(4,ndxmax*ndzmax) :: viit
!
! ... local scalars
integer(i4b) :: i
integer(i4b) :: indx_x, indx_z
real(dp) :: ap
!
! ... debug
!
! ... local arrays
!
!
!
write(*,*) '-------------------------calculating new iterate (ait)-------------------------'
write(*,*)
!
!-----------------------------------------------------------------------
!
do i=1, ndxmax*ndzmax
!skip aitken extrapolation if ap = 0
   ap=2.d0*viit(3,i) - viit(4,i) - viit(2,i)
   if(abs(ap).lt.1.d-14) then
      viit(1,i)=viit(4,i)
   else
      viit(1,i)=viit(3,i) + (viit(3,i)-viit(2,i))*(viit(3,i)-viit(2,i))/ap
   endif
enddo
!
!--------------------back conversion on 3-d-array-----------------------
!
!source2d=0.d0
do i=1, ndxmax*ndzmax
   call conv_indx_1d_to_2d (i, ndxmax, indx_x, indx_z)
   if(viit(1,i).ge.0.d0) then
      source2d(indx_x,indx_z)=viit(1,i)
   endif
enddo
!
end subroutine ait_expol2d
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine store_ng3d(iit, viit, ndxmax, ndymax, ndzmax, source3d)
!
!--------stores subsequent iterates for ng/aitkens-extrapolation--------
!
! input: iit:   iteration step
!        ndxmax, ndymax, ndzmax:   dimensions of 3d array
!        source3d:   3d source function
!
! output: viit: 3d source function at current iteration step iit as 1d array
!
use prog_type
!
implicit none
!
! ... arguments
integer(i4b), intent(in) :: iit, ndxmax, ndymax, ndzmax
real(dp), dimension(ndxmax, ndymax, ndzmax), intent(in) :: source3d
real(dp), dimension(4,ndxmax*ndymax*ndzmax) :: viit
!
! ... local scalars
integer(i4b) :: i,j,k
integer(i4b) :: indx_1d
!
!-----------------------------------------------------------------------
!
select case(iit)
   case(1) 
      write(*,*) '----------------------storing source fct. at step n-3--------------------------'
      write(*,*)
   case(2) 
      write(*,*) '----------------------storing source fct. at step n-2--------------------------'
      write(*,*)
   case(3) 
      write(*,*) '----------------------storing source fct. at step n-1--------------------------'
      write(*,*)
   case(4) 
      write(*,*) '----------------------storing source fct. at step n----------------------------'
      write(*,*)
   case default
end select
!
!
!
do i=1, ndxmax
   do j=1, ndymax
      do k=1, ndzmax
         call conv_indx_3d_to_1d (i, j, k, ndxmax, ndymax, indx_1d)
         viit(iit,indx_1d) = source3d(i,j,k)
      enddo
   enddo
enddo
!
end subroutine store_ng3d
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine ng_expol3d(viit, ndxmax, ndymax, ndzmax, source3d)
!
use prog_type
!
implicit none
!
! ... arguments
integer(i4b), intent(in) :: ndxmax, ndymax, ndzmax
real(dp), dimension(ndxmax,ndymax,ndzmax) :: source3d
real(dp), dimension(4,ndxmax*ndymax*ndzmax) :: viit
!
! ... local scalars
integer(i4b) :: i
integer(i4b) :: indx_x, indx_y, indx_z
real(dp) :: a1, a2, b1, b2, c1, c2
real(dp) :: ap, bp, cp
!
! ... local arrays
!
! ... local logicals
logical :: negative
!
!-----------------------------------------------------------------------
!
write(*,*) '-------------------------calculating new iterate (ng)--------------------------'
write(*,*)
!
!
a1=0.d0
a2=0.d0
b1=0.d0
b2=0.d0
c1=0.d0
c2=0.d0
!
do i=1, ndxmax*ndymax*ndzmax
   a1 = a1+(viit(4,i)-2.d0*viit(3,i)+viit(2,i))**2.d0
   b1 = b1+(viit(4,i)-viit(3,i)-viit(2,i)+viit(1,i))*(viit(4,i)-2.d0*viit(3,i)+viit(2,i))
   c1 = c1+(viit(4,i)-viit(3,i))*(viit(4,i)-2.d0*viit(3,i)+viit(2,i))
   a2 = a2+(viit(4,i)-viit(3,i)-viit(2,i)+viit(1,i))*(viit(4,i)-2.d0*viit(3,i)+viit(2,i))
   b2 = b2+(viit(4,i)-viit(3,i)-viit(2,i)+viit(1,i))**2.d0
   c2 = c2+(viit(4,i)-viit(3,i))*(viit(4,i)-viit(3,i)-viit(2,i)+viit(1,i))
end do
!
ap=(b2*c1-b1*c2)/(a1*b2-a2*b1)
bp=(a1*c2-a2*c1)/(a1*b2-a2*b1)

!extrapolation can yield to negative source-functions
!        if so, then skip extrapolation
negative=.false.

do i=1, ndxmax*ndymax*ndzmax
   cp=(1.d0-ap-bp)*viit(4,i) + ap*viit(3,i)+bp*viit(2,i)
   if(cp.lt.0.d0) then
      write(*,*) 'ng-extrapolation yields negative source-functions => skip'
      negative=.true.
      exit
    endif
enddo
!
if(.not.negative) then
   viit(1,:)=(1.d0-ap-bp)*viit(4,:) + ap*viit(3,:)+bp*viit(2,:)
else
   viit(1,:) = viit(4,:)
endif
!
!--------------------back conversion on 3-d-array-----------------------
!
source3d=0.d0
do i=1, ndxmax*ndymax*ndzmax
   call conv_indx_1d_to_3d (i, ndxmax, ndymax, indx_x, indx_y, indx_z)
   source3d(indx_x,indx_y,indx_z)=viit(1,i)
enddo
!
!
end subroutine ng_expol3d
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine ait_expol3d(viit,ndxmax, ndymax, ndzmax, source3d)
!
use prog_type
!
implicit none
!
! ... arguments
integer(i4b), intent(in) :: ndxmax, ndymax, ndzmax
real(dp), dimension(ndxmax,ndymax,ndzmax) :: source3d
real(dp), dimension(4,ndxmax*ndymax*ndzmax) :: viit
!
! ... local scalars
integer(i4b) :: i
integer(i4b) :: indx_x, indx_y, indx_z
real(dp) :: ap
!
! ... debug
!
! ... local arrays
!
!
!
write(*,*) '-------------------------calculating new iterate (ait)-------------------------'
write(*,*)
!
!-----------------------------------------------------------------------
!
do i=1, ndxmax*ndymax*ndzmax
!skip aitken extrapolation if ap = 0
   ap=2.d0*viit(3,i) - viit(4,i) - viit(2,i)
   if(abs(ap).lt.1.d-14) then
      viit(1,i)=viit(4,i)
   else
      viit(1,i)=viit(3,i) + (viit(3,i)-viit(2,i))*(viit(3,i)-viit(2,i))/ap
   endif
enddo
!
!--------------------back conversion on 3-d-array-----------------------
!
source3d=0.d0
do i=1, ndxmax*ndymax*ndzmax
   call conv_indx_1d_to_3d (i, ndxmax, ndymax, indx_x, indx_y, indx_z)
   source3d(indx_x,indx_y,indx_z)=viit(1,i)
enddo
!
end subroutine ait_expol3d
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine ng_expol1d(viit, nd)
!
!-------same as ng_expol_3d, but: input and output in vector form-------
!
use prog_type
!
implicit none
!
! ... arguments
integer(i4b), intent(in) :: nd
real(dp), dimension(4,nd) :: viit
!
! ... local scalars
integer(i4b) :: i
real(dp) :: a1, a2, b1, b2, c1, c2
real(dp) :: ap, bp, cp
!
! ... local arrays
!
! ... local logicals
logical :: negative
!
!-----------------------------------------------------------------------
!
write(*,*) '-------------------------calculating new iterate (ng)--------------------------'
write(*,*)
!
!
a1=0.d0
a2=0.d0
b1=0.d0
b2=0.d0
c1=0.d0
c2=0.d0
!
do i=1, nd
   a1 = a1+(viit(4,i)-2.d0*viit(3,i)+viit(2,i))**2.d0
   b1 = b1+(viit(4,i)-viit(3,i)-viit(2,i)+viit(1,i))*(viit(4,i)-2.d0*viit(3,i)+viit(2,i))
   c1 = c1+(viit(4,i)-viit(3,i))*(viit(4,i)-2.d0*viit(3,i)+viit(2,i))
   a2 = a2+(viit(4,i)-viit(3,i)-viit(2,i)+viit(1,i))*(viit(4,i)-2.d0*viit(3,i)+viit(2,i))
   b2 = b2+(viit(4,i)-viit(3,i)-viit(2,i)+viit(1,i))**2.d0
   c2 = c2+(viit(4,i)-viit(3,i))*(viit(4,i)-viit(3,i)-viit(2,i)+viit(1,i))
end do
!
ap=(b2*c1-b1*c2)/(a1*b2-a2*b1)
bp=(a1*c2-a2*c1)/(a1*b2-a2*b1)

!extrapolation can yield to negative source-functions
!        if so, then skip extrapolation
negative=.false.

do i=1, nd
   cp=(1.d0-ap-bp)*viit(4,i) + ap*viit(3,i)+bp*viit(2,i)
   if(cp.lt.0.d0) then
      write(*,*) 'ng-extrapolation yields negative source-functions => skip'
      negative=.true.
      exit
    endif
enddo
!
if(.not.negative) then
   viit(1,:)=(1.d0-ap-bp)*viit(4,:) + ap*viit(3,:)+bp*viit(2,:)
endif
!
!
!
end subroutine ng_expol1d
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine ait_expol1d(viit,nd)
!
use prog_type
!
implicit none
!
! ... arguments
integer(i4b), intent(in) :: nd
real(dp), dimension(4,nd) :: viit
!
! ... local scalars
integer(i4b) :: i
real(dp), parameter :: cmax=10.d0, cmin=1.d-1
real(dp) :: qfac, del_si0, del_si1, corr, corrfac, snew
!
! ... local arrays
!
! ... local logicals
logical :: skip
!
!
!
write(*,*) '------------------------calculating new iterate (ait)--------------------------'
write(*,*)
!
snew=0.d0
!
do i=1, nd
   del_si0=viit(4,i)-viit(3,i)
   del_si1=viit(3,i)-viit(2,i)
!
   if(abs(del_si1).gt.1.d-14) then
!
      qfac=del_si0/del_si1
!
!skip extrapolation if qfac ge 1.d0
      if(qfac.ge.1.d0) then
         skip=.true.
      else
         skip=.false.
      endif
!
      if(.not.skip) then
         corr=del_si1/(1.d0-qfac)
         snew=viit(2,i)+corr
         corrfac=snew/viit(2,i)
!skip extrapolation if correction is too large or to low
         if(corrfac.gt.cmax) skip=.true.
         if(corrfac.lt.cmin) skip=.true.
      endif

   else
!skip extrapolation if del_si1 eq. 0.d0
      skip=.true.
   endif
!
   if(skip) then
      viit(1,i) = viit(4,i)
   else
      viit(1,i)=snew
   endif
!
enddo
!
!
end subroutine ait_expol1d
