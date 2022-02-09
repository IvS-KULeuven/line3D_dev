!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine grid_triangles
!
use prog_type
use fund_const
use mod_triangles, only: npoints, ntriangles, &
                         points_coords, points_indx, points_xcoord, points_ycoord, points_weight, &
                         triangles_ip, triangles_if
use dime_spec, only: cs1_np, cs1_nzeta, cs2_np, cs2_nzeta
use mod_spectrum, only: cs1_p, cs2_p, cs1_zeta, cs2_zeta, &
     transmat, transmat_inv, translvec1, translvec2, unit_length
use params_spec, only: rstar1, rstar2
!
implicit none
!
! ... local scalars
integer(i4b) :: i, j, k, ip1, ip2, ip3, err
real(dp) :: x1, x2, x3, y1, y2, y3, weight
!
! ... local arrays
real(dp), dimension(3) :: vec_cac, vec_cyc
!
! ... local characters
!
! ... local logicals
!
! ... local functions
!
!
!-----------------------------------------------------------------------
!
npoints=cs1_np*cs1_nzeta + cs2_np*cs2_nzeta
!
allocate(points_xcoord(npoints))
allocate(points_ycoord(npoints))
!
k=1
!
!---------------------------------star 1---------------------------------
!
do j=1, cs1_nzeta
   do i=1, cs1_np
!point in cylindrical coordinates of star 1
      vec_cyc(1)=cs1_p(i)*cos(cs1_zeta(j))
      vec_cyc(2)=cs1_p(i)*sin(cs1_zeta(j))
      vec_cyc(3)=zero
!
!point in cartesian coordinates of star 1
      vec_cac = matmul(transmat,vec_cyc)
!
!point in global cartesian coordinate system (account for possibly different units
      vec_cac = vec_cac*rstar1/unit_length + translvec1
!
!point in global cylindrical coordinates
      vec_cyc = matmul(transmat_inv,vec_cac)
!
      points_xcoord(k) = vec_cyc(1)
      points_ycoord(k) = vec_cyc(2)
      k=k+1
   enddo
enddo
!
!---------------------------------star 2--------------------------------
!
do j=1, cs2_nzeta
   do i=1, cs2_np
!point in cylindrical coordinates of star 2
      vec_cyc(1)=cs2_p(i)*cos(cs2_zeta(j))
      vec_cyc(2)=cs2_p(i)*sin(cs2_zeta(j))
      vec_cyc(3)=zero
!
!point in cartesian coordinates of star 1
      vec_cac = matmul(transmat,vec_cyc)
!
!point in global cartesian coordinate system (account for possibly different units)
      vec_cac = vec_cac*rstar2/unit_length + translvec2
!
!point in global cylindrical coordinates
      vec_cyc = matmul(transmat_inv,vec_cac)
!
      points_xcoord(k) = vec_cyc(1)
      points_ycoord(k) = vec_cyc(2)
!      write(*,*) cs2_p(i), points_xcoord(k), points_ycoord(k), vec_cac(1), vec_cac(2), vec_cac(3)
      k=k+1
   enddo
enddo
!
!-----------------------------------------------------------------------
!
call uniq_points
!
allocate(points_coords(2,npoints))
allocate(points_indx(npoints))
!
do i=1, npoints
   points_coords(1,i)=points_xcoord(i)
   points_coords(2,i)=points_ycoord(i)
   points_indx(i)=i
enddo


ntriangles=2*npoints-3-2
allocate(triangles_ip(3,ntriangles))
allocate(triangles_if(3,ntriangles))
!
call dtris2 (npoints, points_coords, points_indx, ntriangles, triangles_ip, triangles_if, err)
!
!------------------------create integration weights-----------------------
!
deallocate(points_xcoord)
deallocate(points_ycoord)

allocate(points_xcoord(npoints))
allocate(points_ycoord(npoints))
allocate(points_weight(npoints))
!
points_weight=zero
!
do i=1, ntriangles
   ip1 = triangles_ip(1,i)
   ip2 = triangles_ip(2,i)
   ip3 = triangles_ip(3,i)

   x1=points_coords(1,ip1)
   y1=points_coords(2,ip1)

   x2=points_coords(1,ip2)
   y2=points_coords(2,ip2)

   x3=points_coords(1,ip3)
   y3=points_coords(2,ip3)
   
!integration over this triangle
   weight = ((x2-x1)*(y3-y1) + (x3-x1)*(y1-y2))/6.d0
   points_weight(ip1)=points_weight(ip1)+weight
   points_weight(ip2)=points_weight(ip2)+weight
   points_weight(ip3)=points_weight(ip3)+weight

   points_xcoord(ip1) = points_coords(1,ip1)
   points_ycoord(ip1) = points_coords(2,ip1)
   points_xcoord(ip2) = points_coords(1,ip2)
   points_ycoord(ip2) = points_coords(2,ip2)
   points_xcoord(ip3) = points_coords(1,ip3)
   points_ycoord(ip3) = points_coords(2,ip3)   
!
enddo
!
!do i=1, npoints
!   write(*,*) i, points_xcoord(i), points_ycoord(i), sqrt(points_xcoord(i)**2+points_ycoord(i)**2)
!   if(sqrt(points_xcoord(i)**2+points_ycoord(i)**2).gt.2.&
!        .and.sqrt(points_xcoord(i)**2+points_ycoord(i)**2).lt.3.) stop
!enddo
!stop
!
!
end subroutine grid_triangles
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine uniq_points
!
!finds unique points
!
use prog_type
use fund_const
use mod_triangles
!
implicit none
!
! ... arguments
!
! ... local scalars
integer(i4b) :: i, j, npoints_uniq
real(dp) :: dij
!
! ... local arrays
real(dp), dimension(:), allocatable :: point_xcoord_uniq, point_ycoord_uniq
logical, dimension(:), allocatable :: point_ldelete
!
!-------------------find uniq values------------------------------------
!
allocate(point_ldelete(npoints))
point_ldelete=.false.
!
do i=1, npoints
   do j=i+1, npoints
!calculate distance
      dij = (points_xcoord(i)-points_xcoord(j))**2 + (points_ycoord(i)-points_ycoord(j))**2
      if(dij.lt.small_number) then
         point_ldelete(i)=.true.
      endif
   enddo
enddo
!
npoints_uniq=npoints
do i=1, npoints
   if(point_ldelete(i)) npoints_uniq=npoints_uniq-1
enddo
!
allocate(point_xcoord_uniq(npoints_uniq))
allocate(point_ycoord_uniq(npoints_uniq))
!
j=1
do i=1, npoints
   if(.not.point_ldelete(i)) then
      point_xcoord_uniq(j)=points_xcoord(i)
      point_ycoord_uniq(j)=points_ycoord(i)
      j=j+1
   endif
enddo
!
deallocate(points_xcoord)
deallocate(points_ycoord)
!
npoints=npoints_uniq
allocate(points_xcoord(npoints))
allocate(points_ycoord(npoints))
!
points_xcoord=point_xcoord_uniq
points_ycoord=point_ycoord_uniq
!
end subroutine uniq_points
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine deallocate_triangles
!
use prog_type
use fund_const
use mod_triangles, only: npoints, ntriangles, &
                         points_coords, points_indx, points_xcoord, points_ycoord, points_weight, &
                         triangles_ip, triangles_if
!
implicit none
!
deallocate(points_coords)
deallocate(points_indx)
deallocate(points_xcoord)
deallocate(points_ycoord)
deallocate(points_weight)
deallocate(triangles_ip)
deallocate(triangles_if)
!
end subroutine deallocate_triangles
