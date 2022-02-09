!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine calc_transmat
!
!---------------calculation of transformation matrix--------------------
!--------for system (ex, ey, ez) to system (eex, eey, eez)--------------
!
use prog_type
use fund_const
use params_spec, only: x01, y01, z01, x02, y02, z02, rstar1, rstar2, rmax0, rmax1, rmax2
use mod_spectrum, only: alpha, gamma, nhat, transmat, transmat_inv, unit_length, translvec1, translvec2
!
implicit none
!
! ... local scalars
real(dp) :: check1, check2, check3
!
! ... local arrays
real(dp), dimension(3) :: ex, ey, ez, eex, eey, eez
!
write(*,*) '----------calculate transformation matrix--------------'
write(*,*) '(alpha, gamma): ', alpha, gamma
write(*,*)
!
!calculate nhat (aligned with z-axis of cylindrical coordinate system)
nhat(1)=sin(alpha)*cos(gamma)
nhat(2)=sin(alpha)*sin(gamma)
nhat(3)=cos(alpha)
if(abs(nhat(1)).lt.1.d-14) nhat(1)=zero
if(abs(nhat(2)).lt.1.d-14) nhat(2)=zero
if(abs(nhat(3)).lt.1.d-14) nhat(3)=zero
!
write(*,*) 'nhat'
write(*,'(3es20.8)') nhat
write(*,*)
!calculate carthesian coordinate system that is aligned with 
!cylindrical coordinates
!
!ez predefined: direction to observer
ez=nhat
!
!calculate orthogonal ex from dot_product(ex,ez)=0
if(ez(1).eq.zero) then
   if(ez(2).eq.zero) then
      ex = (/ one, zero, zero /)
   else if(ez(3).eq.zero) then
      ex = (/ one, zero, zero /)
   else
      ex(1) = one
      ex(2) = one
      ex(3) = -ex(2)*ez(2)/ez(3)
   endif
else if (ez(2).eq.zero) then
   if(ez(1).eq.zero) then
      ex = (/ one, zero, zero /)
   else if(ez(3).eq.zero) then
      ex = (/ zero, one, zero /)
   else
      ex(3) = one
      ex(2) = one
      ex(1) = -ex(3)*ez(3)/ez(1)
   endif
else if (ez(3).eq.zero) then
   if(ez(1).eq.zero) then
      ex = (/ one, zero, zero /)
   else if(ez(2).eq.zero) then
      ex = (/ zero, one, zero /)
   else
      ex(3) = one
      ex(2) = one
      ex(1) = -ex(2)*ez(2)/ez(1)
   endif
else
   ex(1) = one
   ex(2) = one
   ex(3) = (-ex(1)*ez(1)-ex(2)*ez(2))/ez(3)
endif
!
!calculate orthogonal ey from cross-product
call cross_product(ez, ex, ey)
!
!normalize unit vectors
ex = ex/sqrt(dot_product(ex,ex))
ey = ey/sqrt(dot_product(ey,ey))
ez = ez/sqrt(dot_product(ez,ez))
!
!check for orthogonality
check1=dot_product(ex,ey)
check2=dot_product(ex,ez)
check3=dot_product(ey,ez)
!
if(abs(check1).gt.1.d-14) stop 'error in calc_transmat: ex, ey not orthogonal'
if(abs(check2).gt.1.d-14) stop 'error in calc_transmat: ex, ez not orthogonal'
if(abs(check3).gt.1.d-14) stop 'error in calc_transmat: ey, ez not orthogonal'
!
!----transformation matrix from system (ex,ey,ez) to (eex, eey, eez)----
!
eex = (/ one, zero, zero /)
eey = (/ zero, one, zero /)
eez = (/ zero, zero, one /)
!
transmat = reshape((/ dot_product(eex,ex), dot_product(eex,ey), dot_product(eex,ez), &
                      dot_product(eey,ex), dot_product(eey,ey), dot_product(eey,ez), &
                      dot_product(eez,ex), dot_product(eez,ey), dot_product(eez,ez) /), shape(transmat))
transmat = transpose(transmat)
!
transmat_inv = transpose(transmat)
!
write(*,*) 'basis of coordinate system'
write(*,'(a5, 3es20.8)') 'e_x', ex
write(*,'(a5, 3es20.8)') 'e_y', ey
write(*,'(a5, 3es20.8)') 'e_z', ez
write(*,*)
!
write(*,*) 'transformation matrix and its inverse'
write(*,'(3es20.8)') transmat
write(*,*)
write(*,'(3es20.8)') transmat_inv
write(*,*)
!
!---------------------define translation vectors-------------------------
!
!unit length defined in input reading procedure (when transforming opacities to correct units)
!
!translation vector to origin of star 1 (input is already in unit_length)
translvec1 = (/ x01, y01, z01 /)
!translvec1 = matmul(transmat,translvec1)
!
!translation vector to origin of star 2 (input is already in unit_length)
translvec2 = (/ x02, y02, z02 /)
!translvec2 = matmul(transmat,translvec2)
!
write(*,*) 'translation vector to star 1 and 2'
write(*,'(3es20.8)') translvec1
write(*,*)
write(*,'(3es20.8)') translvec2
write(*,*)
!
!calculate maximum allowed radius in global coordinate system (note, x01, y01, x02, y02 already measured in unit_length)
rmax0 = sqrt((x02-x01)**2 + (y02-y01)**2) + rmax1*rstar1/unit_length + rmax2*rstar2/unit_length
!
!
end subroutine calc_transmat
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine calc_rotmat
!
!---------------calculation of rotation matrix allowing-----------------
!----------------for tilts of the individual systems--------------------
!
!rotmat1 transforms from sigma_1 -> sigma_0
!rotmat2 transforms from sigma_2 -> sigma_0
!
use prog_type
use fund_const
use params_spec, only: ex01, ey01, ez01, ex02, ey02, ez02
use mod_spectrum, only: rotmat1, rotmat2, rotmat1_inv, rotmat2_inv
!
implicit none
!
! ... local scalars
real(dp) :: check1, check2, check3
!
! ... local arrays
real(dp), dimension(3) :: eex, eey, eez
!
write(*,*) '------------calculate rotation matrices----------------'
write(*,*)
!
!normalize unit vectors star 1
ex01 = ex01/sqrt(dot_product(ex01,ex01))
ey01 = ey01/sqrt(dot_product(ey01,ey01))
ez01 = ez01/sqrt(dot_product(ez01,ez01))
!
!check for orthogonality
check1=dot_product(ex01,ey01)
check2=dot_product(ex01,ez01)
check3=dot_product(ey01,ez01)
!
if(abs(check1).gt.1.d-14) stop 'error in calc_rotmat: ex01, ey01 not orthogonal'
if(abs(check2).gt.1.d-14) stop 'error in calc_rotmat: ex01, ez01 not orthogonal'
if(abs(check3).gt.1.d-14) stop 'error in calc_rotmat: ey01, ez01 not orthogonal'
!
!
!normalize unit vectors
ex02 = ex02/sqrt(dot_product(ex02,ex02))
ey02 = ey02/sqrt(dot_product(ey02,ey02))
ez02 = ez02/sqrt(dot_product(ez02,ez02))
!
!check for orthogonality
check1=dot_product(ex02,ey02)
check2=dot_product(ex02,ez02)
check3=dot_product(ey02,ez02)
!
if(abs(check1).gt.1.d-14) stop 'error in calc_rotmat: ex02, ey02 not orthogonal'
if(abs(check2).gt.1.d-14) stop 'error in calc_rotmat: ex02, ez02 not orthogonal'
if(abs(check3).gt.1.d-14) stop 'error in calc_rotmat: ey02, ez02 not orthogonal'
!
!----transformation matrix from system (ex,ey,ez) to (eex, eey, eez)----
!
eex = (/ one, zero, zero /)
eey = (/ zero, one, zero /)
eez = (/ zero, zero, one /)
!
!star 1
rotmat1 = reshape((/ dot_product(eex,ex01), dot_product(eex,ey01), dot_product(eex,ez01), &
                    dot_product(eey,ex01), dot_product(eey,ey01), dot_product(eey,ez01), &
                    dot_product(eez,ex01), dot_product(eez,ey01), dot_product(eez,ez01) /), shape(rotmat1))
rotmat1 = transpose(rotmat1)
!
rotmat1_inv = transpose(rotmat1)
!
!
!star 2
rotmat2 = reshape((/ dot_product(eex,ex02), dot_product(eex,ey02), dot_product(eex,ez02), &
                     dot_product(eey,ex02), dot_product(eey,ey02), dot_product(eey,ez02), &
                     dot_product(eez,ex02), dot_product(eez,ey02), dot_product(eez,ez02) /), shape(rotmat2))
rotmat2 = transpose(rotmat2)
!
rotmat2_inv = transpose(rotmat2)

!
write(*,*) 'basis of global coordinate system'
write(*,'(a5, 3es20.8)') 'e_x', eex
write(*,'(a5, 3es20.8)') 'e_y', eey
write(*,'(a5, 3es20.8)') 'e_z', eez
write(*,*)
!
write(*,*) 'basis of coordinate system 1 in global system'
write(*,'(a5, 3es20.8)') 'e_x', ex01
write(*,'(a5, 3es20.8)') 'e_y', ey01
write(*,'(a5, 3es20.8)') 'e_z', ez01
write(*,*)
!
write(*,*) 'basis of coordinate system 2 in global system'
write(*,'(a5, 3es20.8)') 'e_x', ex02
write(*,'(a5, 3es20.8)') 'e_y', ey02
write(*,'(a5, 3es20.8)') 'e_z', ez02
write(*,*)
!
write(*,*) 'rotation matrix 1 and its inverse'
write(*,'(3es20.8)') rotmat1
write(*,*)
write(*,'(3es20.8)') rotmat1_inv
write(*,*)
!
write(*,*) 'rotation matrix 2 and its inverse'
write(*,'(3es20.8)') rotmat2
write(*,*)
write(*,'(3es20.8)') rotmat2_inv
write(*,*)
!
!stop 'go on in calc_rotmat'
!
end subroutine calc_rotmat
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine transform_cs1_cs0(cs1_x, cs1_y, cs1_z, unit_length1, unit_length0, translvec1, transmat, transmat_inv, cs0_x, cs0_y, cs0_z)
!
!----------transform given coordinates cs1_x, cs1_y, cs1_z--------------
!----------------------in cylindrical system 1--------------------------
!-----------------------------to----------------------------------------
!------------------coordinates cs0_x, cs0_y, cs0_z----------------------
!-----------------------in cylindrical system 0-------------------------  
!
!input:
!   cs1_x, cs1_y, cs1_z:   coordinates in cylindrical system 1
!   unit_length1:          length scale of coordinate system 1
!   unit_length0:          length scale of coordinate system 0
!   translvec1:            translation-vector pointing to origin of system 1
!   transmat:              transformation matrix for cylindrical->cartesian
!   transmat_inv:          transformation matrix for cartesian->cylindrical
  
use prog_type
use fund_const

!
implicit none
!
! ... arguments
real(dp), intent(in) :: cs1_x, cs1_y, cs1_z, unit_length1, unit_length0
real(dp), dimension(3),intent(in) :: translvec1
real(dp), dimension(3,3),intent(in) :: transmat, transmat_inv
real(dp), intent(out) :: cs0_x, cs0_y, cs0_z
!
! ... local arrays
real(dp), dimension(3) :: vec_cyc, vec_cac
!
!
!coordinates in cylindrical system 1
vec_cyc(1) = cs1_x
vec_cyc(2) = cs1_y
vec_cyc(3) = cs1_z
!
!coordinates in cartesian system of system 1
vec_cac = matmul(transmat,vec_cyc)

!coordinates in global cartesian system
vec_cac = vec_cac*unit_length1/unit_length0 + translvec1
!
!coordinates in global cylindrical system  
vec_cyc = matmul(transmat_inv,vec_cac)
!
cs0_x = vec_cyc(1)
cs0_y = vec_cyc(2)
cs0_z = vec_cyc(3)
!
end subroutine transform_cs1_cs0
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine transform_cs0_cs1(cs0_x, cs0_y, cs0_z, unit_length0, unit_length1, translvec1, transmat, transmat_inv, cs1_x, cs1_y, cs1_z)
!
!----------transform given coordinates cs0_x, cs0_y, cs0_z--------------
!----------------------in cylindrical system 0--------------------------
!-----------------------------to----------------------------------------
!------------------coordinates cs1_x, cs1_y, cs1_z----------------------
!-----------------------in cylindrica system 1--------------------------
!
!input:
!   cs0_x, cs0_y, cs0_z:   coordinates in cylindrical system 0
!   unit_length0:          length scale of coordinate system 0
!   unit_length1:          length scale of coordinate system 1
!   translvec1:            translation-vector pointing to origin of system 1
!   transmat:              transformation matrix for cylindrical->cartesian
!   transmat_inv:          transformation matrix for cartesian->cylindrical
  
use prog_type
use fund_const

!
implicit none
!
! ... arguments
real(dp), intent(in) :: cs0_x, cs0_y, cs0_z, unit_length1, unit_length0
real(dp), dimension(3),intent(in) :: translvec1
real(dp), dimension(3,3),intent(in) :: transmat, transmat_inv
real(dp), intent(out) :: cs1_x, cs1_y, cs1_z
!
! ... local arrays
real(dp), dimension(3) :: vec_cyc, vec_cac
!
!
!coordinates in cylindrical system 0
vec_cyc(1) = cs0_x
vec_cyc(2) = cs0_y
vec_cyc(3) = cs0_z
!
!coordinates in cartesian system of system 0
vec_cac = matmul(transmat,vec_cyc)

!coordinates in cartesian system of 1 (accounting for different length scale)
vec_cac = (vec_cac - translvec1)*unit_length0/unit_length1
!
!coordinates in cylindrical system 1
vec_cyc = matmul(transmat_inv,vec_cac)
!
cs1_x = vec_cyc(1)
cs1_y = vec_cyc(2)
cs1_z = vec_cyc(3)
!
end subroutine transform_cs0_cs1
