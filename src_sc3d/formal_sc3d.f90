!
!***********************************************************************
!***********************************************************************
!
!        SOME ROUTINES USED FOR 3D CONTINUUM AND LINE TRANSPORT
!
!***********************************************************************
!***********************************************************************
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine set_boundary3d(n_x, n_y, n_z)
!
   use prog_type
   use dime3d, only: ndxmax, ndymax, ndzmax, x, y, z, imask3d, imask_innreg3d, int3d
   use bcondition, only: xic1
!
   implicit none
!
! ... arguments
   real(dp), intent(in) :: n_x, n_y, n_z
!
! ... local scalars
   integer(i4b) :: i, j, k, iim2, iim1, iip1, jjm2, jjm1, jjp1, kkm2, kkm1, kkp1, alpha, beta, gamma
   real(dp) :: mueff, rad
!
! ... local functions
   real(dp) :: calc_icore_gdark
!
   int3d=0.d0
!
   do k=3, ndzmax-2
      do j=3, ndymax-2
         do i=3, ndxmax-2
            if(imask3d(i,j,k).eq.4) then
!inside the star, set intensities correctly
!            write(*,'(3i5,10es20.8)') i, j, k, rad, x(i), y(j), z(k)
               rad=sqrt(x(i)**2+y(j)**2+z(k)**2)
               mueff=(n_x*x(i)+n_y*y(j)+n_z*z(k))/rad
               if(mueff.ge.0.d0) then
                  int3d(i,j,k) = calc_icore_gdark(z(k), rad)
!               int3d(i,j,k) = xic1
               else
                  int3d(i,j,k) = 0.d0
               endif
            endif
         enddo
      enddo
   enddo
!
!
!
end subroutine set_boundary3d
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
function calc_icore_gdark(zp, rad)
!
!this function calculates the intensity emerging from the core
!for a given xic1(theta), e.g. from gravity darkening, and for a given z-coordinate
!
   use prog_type
   use fund_const, only: pi
   use bcondition, only: ntheta_gdark, theta_gdark, xic1_gdark
   use mod_interp1d, only: interpol_yp
!
   implicit none
!
! ... arguments
   real(dp), intent(in) :: zp, rad
   real(dp) :: calc_icore_gdark
!
! ... local scalars
   integer(i4b) :: indx_gdark_iim1, indx_gdark_ii
   real(dp) :: theta, indx
!
! ... local functions
!
!calculate local theta
   theta = acos(abs(zp)/rad)
!
   if(abs(theta).lt.1.d-14) then
!polar value
      calc_icore_gdark = xic1_gdark(1)
   elseif(abs(theta-pi/2.d0).lt.1.d-14) then
!equatorial value
      calc_icore_gdark = xic1_gdark(ntheta_gdark)
   else
!interpolation
      indx = log10(1.d0-1.8d0*theta/pi)*(1.d0-ntheta_gdark)+1.d0
      indx_gdark_iim1=floor(indx)
      indx_gdark_ii=ceiling(indx)
      calc_icore_gdark = interpol_yp(theta_gdark(indx_gdark_iim1), theta_gdark(indx_gdark_ii), &
         xic1_gdark(indx_gdark_iim1), xic1_gdark(indx_gdark_ii), theta)
   endif
!
end function calc_icore_gdark
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine coeffcr1d_mbez(fc, f_im1, f_i, at, bt, ct, at2, bt2, ct2)
!
!calculates new coefficients for 2d bezier interpolation to ensure monotonicity
!                 when interpolating in right interval
!
!on input:
!   fc           - control point
!   f_im1, f_i   - function values that limit the interval
!   at, bt, ct   - coefficients that have been used for calculation of control point
!
!on output:
!   at2, bt2, ct2 - interpolation coefficients for control point sucht that monotonic interpolation can be used
!                   (e.g. at2=0.,bt2=0.,ct2=1. if f_c > f_i > f_im1)
!
!
   use prog_type
!
   implicit none
!
! ... arguments
   real(dp), intent(in) :: fc, f_im1, f_i, at, bt, ct
   real(dp), intent(out) :: at2, bt2, ct2
!
!--------------------------old version----------------------------------
!
!if(f_i.ge.f_im1) then
!   if(fc.gt.f_i) then
!      at2=0.
!      bt2=0.
!      ct2=1.
!   elseif(fc.lt.f_im1) then
!      at2=0.
!      bt2=1.
!      ct2=0.
!   else
!      at2=at
!      bt2=bt
!      ct2=ct
!   endif
!elseif(f_i.le.f_im1) then
!   if(fc.lt.f_i) then
!      at2=0.
!      bt2=0.
!      ct2=1.
!   elseif(fc.gt.f_im1) then
!      at2=0.
!      bt2=1.
!      ct2=0.
!   else
!      at2=at
!      bt2=bt
!      ct2=ct
!   endif
!endif
!
!--------------------------new version----------------------------------
!
   if((f_i-f_im1)*(f_i-fc).le.0.) then
      at2=0.
      bt2=0.
      ct2=1.
   elseif((f_i-f_im1)*(fc-f_im1).le.0.) then
      at2=0.
      bt2=1.
      ct2=0.
   else
      at2=at
      bt2=bt
      ct2=ct
   endif
!
end subroutine coeffcr1d_mbez
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine pointcr1d_mbez(fc, f_im1, f_i)
!
!calculates new control point for bezier interpolation to ensure monotonicity
!                 when interpolating in right interval
!
!on input:
!   fc           - current control point
!   f_im1, f_i   - function values that limit the interval
!
!on output:
!   fc - new control point such that monotonicity is ensured
!               (e.g. fc_new = f_i, if f_c > f_i > f_im1)
!
!
   use prog_type
!
   implicit none
!
! ... arguments
   real(dp), intent(in) :: f_im1, f_i
   real(dp), intent(inout) :: fc
!
!
!
   if((f_i-f_im1)*(f_i-fc).le.0.) then
      fc=f_i
   elseif((f_i-f_im1)*(fc-f_im1).le.0.) then
      fc=f_im1
   endif
!
end subroutine pointcr1d_mbez
!
!***********************************************************************
!***********************************************************************
!
!                 CONTINUUM TRASNFER ROUTINES
!
!***********************************************************************
!***********************************************************************
!
subroutine fsc_cont3d(oindx)
!
!-----------------------------------------------------------------------
!------short characteristics for continuum radiative transfer in 3d-----
!-----------calculating intensties for given mu,phi specified-----------
!---------------------------by input oindx------------------------------
!-----------------------------------------------------------------------
!
   use prog_type

   use fund_const
   use dime3d, only: ndxmax, ndzmax, ndymax, x, y, z
   use dime3d, only: int3d, opac3d, scont3d, imask3d, imask_totreg3d, &
      alocont_o_nn3d, alocont_nn3d_tmp, mint3d_tmp, normalization3d_tmp, &
      fcontx3d_tmp, fconty3d_tmp, fcontz3d_tmp, &
      kcontxx3d_tmp, kcontyy3d_tmp, kcontzz3d_tmp, kcontxy3d_tmp, kcontxz3d_tmp, kcontyz3d_tmp
   use angles, only: n_x, n_y, n_z, weight_omega, q_alo
   use bcondition, only: xic1
   use mod_debug, only: indxx, indxy, indxz
   use params_stellar, only: smajorax_a, smajorax_b, smajorax_c
!
   implicit none
!
! ... arguments
   integer(i4b), intent(in) :: oindx
!
! ... local scalars
   integer(i4b) :: i, j, k
   integer(i4b) :: alpha, beta, gamma
   integer(i4b) :: startx, starty, startz, endx, endy, endz
   integer(i4b) :: iim2, iim1, ii, iip1, jjm2, jjm1, jj, jjp1, kkm2, kkm1, kk, kkp1
   real(dp) :: nn_x, nn_y, nn_z, mueff, wall
   real(dp) :: dels_xyu, dels_xzu, dels_yzu, dels_u, dels_r
   real(dp) :: dels_xyd, dels_xzd, dels_yzd, dels_d
   real(dp) :: opac_p, scont_p
   real(dp) :: x_u, y_u, z_u, int_u, opac_u, scont_u
   real(dp) :: x_d, y_d, z_d, opac_d, scont_d
   real(dp) :: abs_sc, int_sc, contr_sc
   real(dp) :: alo_u, alo_p, alo_d
   integer :: q1, q2, q3, q4, q5, q6, q7, q8, q9, &
      q10, q11, q12, q13, q14, q15, q16, q17, q18, &
      q19, q20, q21, q22, q23, q24, q25, q26, q27
   real(dp) :: c01_scontu, c02_scontu, c03_scontu, c04_scontu, c05_scontu, c06_scontu, c07_scontu, &
      c08_scontu, c09_scontu, c10_scontu, c11_scontu, c12_scontu, c13_scontu, c14_scontu, &
      c15_scontu, c16_scontu, c17_scontu, c18_scontu, c19_scontu, c20_scontu, c21_scontu, &
      c22_scontu, c23_scontu, c24_scontu, c25_scontu, c26_scontu, c27_scontu, &
      c01_opacu, c02_opacu, c03_opacu, c04_opacu, c05_opacu, c06_opacu, c07_opacu, &
      c08_opacu, c09_opacu, c10_opacu, c11_opacu, c12_opacu, c13_opacu, c14_opacu, &
      c15_opacu, c16_opacu, c17_opacu, c18_opacu, c19_opacu, c20_opacu, c21_opacu, &
      c22_opacu, c23_opacu, c24_opacu, c25_opacu, c26_opacu, c27_opacu, &
      c02_intu, c04_intu, c05_intu, c06_intu, c08_intu, c10_intu, c11_intu, &
      c12_intu, c13_intu, c14_intu, c15_intu, c16_intu, c17_intu, c18_intu, &
      c20_intu, c22_intu, c23_intu, c24_intu, c26_intu, &
      c03_scontd, c06_scontd, c07_scontd, c08_scontd, c09_scontd, c12_scontd, c15_scontd, &
      c16_scontd, c17_scontd, c18_scontd, c19_scontd, c20_scontd, c21_scontd, c22_scontd, &
      c23_scontd, c24_scontd, c25_scontd, c26_scontd, c27_scontd, &
      c03_opacd, c06_opacd, c07_opacd, c08_opacd, c09_opacd, c12_opacd, c15_opacd, &
      c16_opacd, c17_opacd, c18_opacd, c19_opacd, c20_opacd, c21_opacd, c22_opacd, &
      c23_opacd, c24_opacd, c25_opacd, c26_opacd, c27_opacd
!
!for debugging
   real(dp) :: int_u2, scont_u2, opac_u2, int_u3, scont_u3, opac_u3, scont_d2, scont_d3, opac_d2, opac_d3, interpol2d_9p_quad, interpol2d_9p_bez
!
! ... local arrays
!
! ... local functions
   real(dp) :: dist_ellipsoid, calc_icore_gdark
!
! ... local characters
!character(len=50) :: enter
!
! ... local logicals
!
!directions
   nn_x=n_x(oindx)
   nn_y=n_y(oindx)
   nn_z=n_z(oindx)
!
!angulare integration weight
   wall=weight_omega(oindx)
!
!indices for nearest neighbour alo
   q1=q_alo(oindx,1)
   q2=q_alo(oindx,2)
   q3=q_alo(oindx,3)
   q4=q_alo(oindx,4)
   q5=q_alo(oindx,5)
   q6=q_alo(oindx,6)
   q7=q_alo(oindx,7)
   q8=q_alo(oindx,8)
   q9=q_alo(oindx,9)
   q10=q_alo(oindx,10)
   q11=q_alo(oindx,11)
   q12=q_alo(oindx,12)
   q13=q_alo(oindx,13)
   q14=q_alo(oindx,14)
   q15=q_alo(oindx,15)
   q16=q_alo(oindx,16)
   q17=q_alo(oindx,17)
   q18=q_alo(oindx,18)
   q19=q_alo(oindx,19)
   q20=q_alo(oindx,20)
   q21=q_alo(oindx,21)
   q22=q_alo(oindx,22)
   q23=q_alo(oindx,23)
   q24=q_alo(oindx,24)
   q25=q_alo(oindx,25)
   q26=q_alo(oindx,26)
   q27=q_alo(oindx,27)
!
!-----------------------------------------------------------------------
!
!set directional index-parameter (phantom points are excluded from calculation)
!
!index parameter:
!         if n_x >= 0                 if n_x < 0
!                startx = 2                  startx = ndxmax-1
!                endx = ndxmax-1             endx = 2
!                alpha=  1                   alpha=-1
!
!         if n_y >= 0                 if n_y < 0
!                starty = 2                  starty = ndymax-1
!                endy = ndymax-1             endy = 2
!                beta =  1                   beta =-1
!
!         if n_z >= 0                 if n_z < 0
!                startz = 2                  startz = ndzmax-1
!                endz = ndzmax-1             endz = 2
!                gamma = 1                   gamma = -1
!
   if(nn_x.gt.0.d0) then
      startx = 3
      endx = ndxmax-2
      alpha=  1
   elseif(nn_x.lt.0.d0) then
      startx = ndxmax-2
      endx = 3
      alpha=-1
   else
      stop 'error in fsc_cont3d: n_x = 0 not allowed'
   endif
!
   if(nn_y.gt.0.d0) then
      starty = 3
      endy = ndymax-2
      beta =  1
   elseif(nn_y.lt.0.d0) then
      starty = ndymax-2
      endy = 3
      beta =-1
   else
      stop 'error in fsc_cont3d: n_y = 0 not allowed'
   endif
!
   if(nn_z.gt.0.d0) then
      startz = 3
      endz = ndzmax-2
      gamma=  1
   elseif(nn_z.lt.0.d0) then
      startz = ndzmax-2
      endz = 3
      gamma=-1
   else
      stop 'error in fsc_cont3d: n_z = 0 not allowed'
   endif
!
!--------------------reset the intensities and alo----------------------
!
   call set_boundary3d(nn_x, nn_y, nn_z)
!
!for alo-tests
!scont3d=0.d0
!scont3d(indxx,indxy,indxz)=1.d0
!int3d=0.d0
!
   alocont_o_nn3d=0.d0
!
!-----------------------------------------------------------------------
!
   do k=startz, endz, gamma
      do j=starty, endy, beta
         do i=startx, endx, alpha
!
!         c02_scontu = 0.d0
!         c04_scontu = 0.d0
!         c05_scontu = 0.d0
!         c06_scontu = 0.d0
!         c08_scontu = 0.d0
!         c10_scontu = 0.d0
!         c11_scontu = 0.d0
!         c12_scontu = 0.d0
!         c13_scontu = 0.d0
!         c14_scontu = 0.d0
!         c15_scontu = 0.d0
!         c16_scontu = 0.d0
!         c17_scontu = 0.d0
!         c18_scontu = 0.d0
!         c20_scontu = 0.d0
!         c22_scontu = 0.d0
!         c23_scontu = 0.d0
!         c24_scontu = 0.d0
!         c26_scontu = 0.d0
!         c27_scontu = 0.d0
!         c02_intu = 0.d0
!         c04_intu = 0.d0
!         c05_intu = 0.d0
!         c06_intu = 0.d0
!         c08_intu = 0.d0
!         c10_intu = 0.d0
!         c11_intu = 0.d0
!         c12_intu = 0.d0
!         c13_intu = 0.d0
!         c14_intu = 0.d0
!         c15_intu = 0.d0
!         c16_intu = 0.d0
!         c17_intu = 0.d0
!         c18_intu = 0.d0
!         c20_intu = 0.d0
!         c22_intu = 0.d0
!         c23_intu = 0.d0
!         c24_intu = 0.d0
!         c26_intu = 0.d0
!         c03_scontd = 0.d0
!         c06_scontd = 0.d0
!         c07_scontd = 0.d0
!         c08_scontd = 0.d0
!         c09_scontd = 0.d0
!         c12_scontd = 0.d0
!         c15_scontd = 0.d0
!         c16_scontd = 0.d0
!         c17_scontd = 0.d0
!         c18_scontd = 0.d0
!         c19_scontd = 0.d0
!         c20_scontd = 0.d0
!         c21_scontd = 0.d0
!         c22_scontd = 0.d0
!         c23_scontd = 0.d0
!         c24_scontd = 0.d0
!         c25_scontd = 0.d0
!         c26_scontd = 0.d0
!         c27_scontd = 0.d0
!
            select case(imask3d(i,j,k))
!
!************************standard radiative transfer********************
!
             case(1)
!
!               call cpu_time(ts_case1)

               iim1=i-alpha
               jjm1=j-beta
               kkm1=k-gamma
               iip1=i+alpha
               jjp1=j+beta
               kkp1=k+gamma

!calculate distance of upwind point to (previous) xy-plane, xz-plane and yz-plane
               dels_xyu=(z(k)-z(kkm1))/nn_z
               dels_xzu=(y(j)-y(jjm1))/nn_y
               dels_yzu=(x(i)-x(iim1))/nn_x
               dels_u=min(dels_xyu,dels_xzu,dels_yzu)
!
!calculate distance of downwind point to (next) xy-plane, xz-plane and yz-plane
               dels_xyd=(z(kkp1)-z(k))/nn_z
               dels_xzd=(y(jjp1)-y(j))/nn_y
               dels_yzd=(x(iip1)-x(i))/nn_x
               dels_d=min(dels_xyd,dels_xzd,dels_yzd)
!
!----------------------------local point--------------------------------
!
               scont_p=scont3d(i,j,k)
               opac_p=opac3d(i,j,k)
!
!----------------------------upwind point-------------------------------
!
               if(dels_xyu.eq.dels_u) then
!intersection with x-y plane on level k-gamma
                  x_u = x(i) - dels_u*nn_x
                  y_u = y(j) - dels_u*nn_y
                  z_u = z(kkm1)
!
                  iim2=i-2*alpha
                  jjm2=j-2*beta
!
!                  call cpu_time(ts_interpu)
                  call coeff2d_contu(opac3d(iim2,jjm2,kkm1), opac3d(iim1,jjm2,kkm1), opac3d(i,jjm2,kkm1), &
                     opac3d(iim2,jjm1,kkm1), opac3d(iim1,jjm1,kkm1), opac3d(i,jjm1,kkm1), &
                     opac3d(iim2,j,kkm1), opac3d(iim1,j,kkm1), opac3d(i,j,kkm1), &
                     scont3d(iim2,jjm2,kkm1), scont3d(iim1,jjm2,kkm1), scont3d(i,jjm2,kkm1), &
                     scont3d(iim2,jjm1,kkm1), scont3d(iim1,jjm1,kkm1), scont3d(i,jjm1,kkm1), &
                     scont3d(iim2,j,kkm1), scont3d(iim1,j,kkm1), scont3d(i,j,kkm1), &
                     int3d(iim2,jjm2,kkm1), int3d(iim1,jjm2,kkm1), int3d(i,jjm2,kkm1), &
                     int3d(iim2,jjm1,kkm1), int3d(iim1,jjm1,kkm1), int3d(i,jjm1,kkm1), &
                     int3d(iim2,j,kkm1), int3d(iim1,j,kkm1), int3d(i,j,kkm1), &
                     x(iim2), x(iim1), x(i), y(jjm2), y(jjm1), y(j), x_u, y_u, &
                     c10_scontu, c11_scontu, c12_scontu, c13_scontu, c14_scontu, &
                     c15_scontu, c16_scontu, c17_scontu, c18_scontu, &
                     c10_intu, c11_intu, c12_intu, c13_intu, c14_intu, &
                     c15_intu, c16_intu, c17_intu, c18_intu, opac_u, scont_u, int_u)
!                  call cpu_time(te_interpu)
!                  tt_interpu=tt_interpu+te_interpu-ts_interpu
!
!set interpolation coefficients that are not used to zero
                  c23_scontu = 0.d0
                  c24_scontu = 0.d0
                  c26_scontu = 0.d0
                  c27_scontu = 0.d0
!
                  c02_intu = 0.d0
                  c04_intu = 0.d0
                  c05_intu = 0.d0
                  c06_intu = 0.d0
                  c08_intu = 0.d0
                  c20_intu = 0.d0
                  c22_intu = 0.d0
                  c23_intu = 0.d0
                  c24_intu = 0.d0
                  c26_intu = 0.d0
!
               elseif(dels_xzu.eq.dels_u) then
!intersection with x-z plane on level j-beta
                  x_u = x(i) - dels_u*nn_x
                  y_u = y(jjm1)
                  z_u = z(k) - dels_u*nn_z
!
                  iim2=i-2*alpha
                  kkm2=k-2*gamma
!
!                  call cpu_time(ts_interpu)
                  call coeff2d_contu(opac3d(iim2,jjm1,kkm2), opac3d(iim1,jjm1,kkm2), opac3d(i,jjm1,kkm2), &
                     opac3d(iim2,jjm1,kkm1), opac3d(iim1,jjm1,kkm1), opac3d(i,jjm1,kkm1), &
                     opac3d(iim2,jjm1,k), opac3d(iim1,jjm1,k), opac3d(i,jjm1,k), &
                     scont3d(iim2,jjm1,kkm2), scont3d(iim1,jjm1,kkm2), scont3d(i,jjm1,kkm2), &
                     scont3d(iim2,jjm1,kkm1), scont3d(iim1,jjm1,kkm1), scont3d(i,jjm1,kkm1), &
                     scont3d(iim2,jjm1,k), scont3d(iim1,jjm1,k), scont3d(i,jjm1,k), &
                     int3d(iim2,jjm1,kkm2), int3d(iim1,jjm1,kkm2), int3d(i,jjm1,kkm2), &
                     int3d(iim2,jjm1,kkm1), int3d(iim1,jjm1,kkm1), int3d(i,jjm1,kkm1), &
                     int3d(iim2,jjm1,k), int3d(iim1,jjm1,k), int3d(i,jjm1,k), &
                     x(iim2), x(iim1), x(i), z(kkm2), z(kkm1), z(k), x_u, z_u, &
                     c04_scontu, c05_scontu, c06_scontu, c13_scontu, c14_scontu, &
                     c15_scontu, c22_scontu, c23_scontu, c24_scontu, &
                     c04_intu, c05_intu, c06_intu, c13_intu, c14_intu, &
                     c15_intu, c22_intu, c23_intu, c24_intu, opac_u, scont_u, int_u)
!                  call cpu_time(te_interpu)
!                  tt_interpu=tt_interpu+te_interpu-ts_interpu
!
!set interpolation coefficients that are not used to zero
                  c17_scontu = 0.d0
                  c18_scontu = 0.d0
                  c26_scontu = 0.d0
                  c27_scontu = 0.d0
!
                  c02_intu = 0.d0
                  c08_intu = 0.d0
                  c10_intu = 0.d0
                  c11_intu = 0.d0
                  c12_intu = 0.d0
                  c16_intu = 0.d0
                  c17_intu = 0.d0
                  c18_intu = 0.d0
                  c20_intu = 0.d0
                  c26_intu = 0.d0
!
               elseif(dels_yzu.eq.dels_u) then
!intersection with y-z plane on level i-alpha
                  x_u = x(iim1)
                  y_u = y(j) - dels_u*nn_y
                  z_u = z(k) - dels_u*nn_z
!
                  jjm2=j-2*beta
                  kkm2=k-2*gamma
!
!                  call cpu_time(ts_interpu)
                  call coeff2d_contu(opac3d(iim1,jjm2,kkm2), opac3d(iim1,jjm1,kkm2), opac3d(iim1,j,kkm2), &
                     opac3d(iim1,jjm2,kkm1), opac3d(iim1,jjm1,kkm1), opac3d(iim1,j,kkm1), &
                     opac3d(iim1,jjm2,k), opac3d(iim1,jjm1,k), opac3d(iim1,j,k), &
                     scont3d(iim1,jjm2,kkm2), scont3d(iim1,jjm1,kkm2), scont3d(iim1,j,kkm2), &
                     scont3d(iim1,jjm2,kkm1), scont3d(iim1,jjm1,kkm1), scont3d(iim1,j,kkm1), &
                     scont3d(iim1,jjm2,k), scont3d(iim1,jjm1,k), scont3d(iim1,j,k), &
                     int3d(iim1,jjm2,kkm2), int3d(iim1,jjm1,kkm2), int3d(iim1,j,kkm2), &
                     int3d(iim1,jjm2,kkm1), int3d(iim1,jjm1,kkm1), int3d(iim1,j,kkm1), &
                     int3d(iim1,jjm2,k), int3d(iim1,jjm1,k), int3d(iim1,j,k), &
                     y(jjm2), y(jjm1), y(j), z(kkm2), z(kkm1), z(k), y_u, z_u, &
                     c02_scontu, c05_scontu, c08_scontu, c11_scontu, c14_scontu, &
                     c17_scontu, c20_scontu, c23_scontu, c26_scontu, &
                     c02_intu, c05_intu, c08_intu, c11_intu, c14_intu, &
                     c17_intu, c20_intu, c23_intu, c26_intu, opac_u, scont_u, int_u)
!                  call cpu_time(te_interpu)
!                  tt_interpu=tt_interpu+te_interpu-ts_interpu
!
!set interpolation coefficients that are not used to zero
                  c15_scontu = 0.d0
                  c18_scontu = 0.d0
                  c24_scontu = 0.d0
                  c27_scontu = 0.d0
!
                  c04_intu = 0.d0
                  c06_intu = 0.d0
                  c10_intu = 0.d0
                  c12_intu = 0.d0
                  c13_intu = 0.d0
                  c15_intu = 0.d0
                  c16_intu = 0.d0
                  c18_intu = 0.d0
                  c22_intu = 0.d0
                  c24_intu = 0.d0
!
               else
                  write(*,'(4es20.8,2l4)') dels_u, dels_xzu, dels_xyu, dels_yzu, dels_u.eq.dels_xzu, dels_xyu.eq.dels_u
                  stop 'error in fsc_cont3d: invalid dels_u'
               endif
!
!---------------------------downwind point------------------------------
!
               if(dels_xyd.eq.dels_d) then
!intersection with x-y plane on level k+gamma
                  x_d = x(i) + dels_d*nn_x
                  y_d = y(j) + dels_d*nn_y
                  z_d = z(kkp1)
!
!                  call cpu_time(ts_interpd)
                  call coeff2d_contd(opac3d(iim1,jjm1,kkp1), opac3d(i,jjm1,kkp1), opac3d(iip1,jjm1,kkp1), &
                     opac3d(iim1,j,kkp1), opac3d(i,j,kkp1), opac3d(iip1,j,kkp1), &
                     opac3d(iim1,jjp1,kkp1), opac3d(i,jjp1,kkp1), opac3d(iip1,jjp1,kkp1), &
                     scont3d(iim1,jjm1,kkp1), scont3d(i,jjm1,kkp1), scont3d(iip1,jjm1,kkp1), &
                     scont3d(iim1,j,kkp1), scont3d(i,j,kkp1), scont3d(iip1,j,kkp1), &
                     scont3d(iim1,jjp1,kkp1), scont3d(i,jjp1,kkp1), scont3d(iip1,jjp1,kkp1), &
                     x(iim1), x(i), x(iip1), y(jjm1), y(j), y(jjp1), x_d, y_d, &
                     c19_scontd, c20_scontd, c21_scontd, c22_scontd, c23_scontd, c24_scontd, &
                     c25_scontd, c26_scontd, c27_scontd, opac_d, scont_d)
!                  call cpu_time(te_interpd)
!                  tt_interpd=tt_interpd+te_interpd-ts_interpd
!
!set interpolation coefficients that are not used to zero
                  c03_scontd = 0.d0
                  c06_scontd = 0.d0
                  c07_scontd = 0.d0
                  c08_scontd = 0.d0
                  c09_scontd = 0.d0
                  c12_scontd = 0.d0
                  c15_scontd = 0.d0
                  c16_scontd = 0.d0
                  c17_scontd = 0.d0
                  c18_scontd = 0.d0
!
               elseif(dels_xzd.eq.dels_d) then
!intersection with x-z plane on level j+beta
                  x_d = x(i) + dels_d*nn_x
                  y_d = y(jjp1)
                  z_d = z(k) + dels_d*nn_z
!
!                  call cpu_time(ts_interpd)
                  call coeff2d_contd(opac3d(iim1,jjp1,kkm1), opac3d(i,jjp1,kkm1), opac3d(iip1,jjp1,kkm1), &
                     opac3d(iim1,jjp1,k), opac3d(i,jjp1,k), opac3d(iip1,jjp1,k), &
                     opac3d(iim1,jjp1,kkp1), opac3d(i,jjp1,kkp1), opac3d(iip1,jjp1,kkp1), &
                     scont3d(iim1,jjp1,kkm1), scont3d(i,jjp1,kkm1), scont3d(iip1,jjp1,kkm1), &
                     scont3d(iim1,jjp1,k), scont3d(i,jjp1,k), scont3d(iip1,jjp1,k), &
                     scont3d(iim1,jjp1,kkp1), scont3d(i,jjp1,kkp1), scont3d(iip1,jjp1,kkp1), &
                     x(iim1), x(i), x(iip1), z(kkm1), z(k), z(kkp1), x_d, z_d, &
                     c07_scontd, c08_scontd, c09_scontd, c16_scontd, c17_scontd, c18_scontd, &
                     c25_scontd, c26_scontd, c27_scontd, opac_d, scont_d)
!                  call cpu_time(te_interpd)
!                  tt_interpd=tt_interpd+te_interpd-ts_interpd
!
!set interpolation coefficients that are not used to zero
                  c03_scontd = 0.d0
                  c06_scontd = 0.d0
                  c12_scontd = 0.d0
                  c15_scontd = 0.d0
                  c19_scontd = 0.d0
                  c20_scontd = 0.d0
                  c21_scontd = 0.d0
                  c22_scontd = 0.d0
                  c23_scontd = 0.d0
                  c24_scontd = 0.d0
!
               elseif(dels_yzd.eq.dels_d) then
!intersection with y-z plane on level i+alpha
                  x_d = x(iip1)
                  y_d = y(j) + dels_d*nn_y
                  z_d = z(k) + dels_d*nn_z
!
!                  call cpu_time(ts_interpd)
                  call coeff2d_contd(opac3d(iip1,jjm1,kkm1), opac3d(iip1,j,kkm1), opac3d(iip1,jjp1,kkm1), &
                     opac3d(iip1,jjm1,k), opac3d(iip1,j,k), opac3d(iip1,jjp1,k), &
                     opac3d(iip1,jjm1,kkp1), opac3d(iip1,j,kkp1), opac3d(iip1,jjp1,kkp1), &
                     scont3d(iip1,jjm1,kkm1), scont3d(iip1,j,kkm1), scont3d(iip1,jjp1,kkm1), &
                     scont3d(iip1,jjm1,k), scont3d(iip1,j,k), scont3d(iip1,jjp1,k), &
                     scont3d(iip1,jjm1,kkp1), scont3d(iip1,j,kkp1), scont3d(iip1,jjp1,kkp1), &
                     y(jjm1), y(j), y(jjp1), z(kkm1), z(k), z(kkp1), y_d, z_d, &
                     c03_scontd, c06_scontd, c09_scontd, c12_scontd, c15_scontd, c18_scontd, &
                     c21_scontd, c24_scontd, c27_scontd, opac_d, scont_d)
!                  call cpu_time(te_interpd)
!                  tt_interpd=tt_interpd+te_interpd-ts_interpd
!
!set interpolation coefficients that are not used to zero
                  c07_scontd = 0.d0
                  c08_scontd = 0.d0
                  c16_scontd = 0.d0
                  c17_scontd = 0.d0
                  c19_scontd = 0.d0
                  c20_scontd = 0.d0
                  c22_scontd = 0.d0
                  c23_scontd = 0.d0
                  c25_scontd = 0.d0
                  c26_scontd = 0.d0
!
               else
                  write(*,'(4es20.8,l4)') dels_d, dels_xzd, dels_xyd, dels_yzd, dels_d.eq.dels_xzd
                  stop 'error in fsc_cont3d: invalid dels_d'
               endif
!
!--------------------------------radiative transfer---------------------
!
!               call cpu_time(ts_fs1d)
               call fsc_cont(int_u, opac_u, opac_p, opac_d, scont_u, scont_p, scont_d, &
                  dels_u, dels_d, abs_sc, contr_sc, int_sc, alo_u, alo_p, alo_d)
               int3d(i,j,k) = abs_sc*int_u + alo_u*scont_u + alo_p*scont_p + alo_d*scont_d
!               call cpu_time(te_fs1d)
!               tt_fs1d=tt_fs1d+te_fs1d-ts_fs1d
!
!               call cpu_time(ts_aloo)
               alocont_o_nn3d(iip1,jjp1,kkp1,q1) = imask_totreg3d(iip1,jjp1,kkp1)*&
                  alo_d*c27_scontd
               alocont_o_nn3d(i,jjp1,kkp1,q2) = imask_totreg3d(i,jjp1,kkp1) * &
                  (alo_d*c26_scontd + abs_sc*c26_intu*alocont_o_nn3d(i,jjp1,kkp1,q1))
               alocont_o_nn3d(iim1,jjp1,kkp1,q3) = imask_totreg3d(iim1,jjp1,kkp1) * &
                  (alo_d*c25_scontd + abs_sc*c26_intu*alocont_o_nn3d(iim1,jjp1,kkp1,q2))
               alocont_o_nn3d(iip1,j,kkp1,q4) = imask_totreg3d(iip1,j,kkp1) * &
                  (alo_d*c24_scontd + abs_sc*c24_intu*alocont_o_nn3d(iip1,j,kkp1,q1))
               alocont_o_nn3d(i,j,kkp1,q5) = imask_totreg3d(i,j,kkp1) * &
                  (alo_d*c23_scontd + abs_sc*(c23_intu*alocont_o_nn3d(i,j,kkp1,q1) + &
                  c24_intu*alocont_o_nn3d(i,j,kkp1,q2) + &
                  c26_intu*alocont_o_nn3d(i,j,kkp1,q4)))
               alocont_o_nn3d(iim1,j,kkp1,q6) = imask_totreg3d(iim1,j,kkp1) * &
                  (alo_d*c22_scontd + abs_sc*(c22_intu*alocont_o_nn3d(iim1,j,kkp1,q1) + &
                  c23_intu*alocont_o_nn3d(iim1,j,kkp1,q2) + &
                  c24_intu*alocont_o_nn3d(iim1,j,kkp1,q3) + &
                  c26_intu*alocont_o_nn3d(iim1,j,kkp1,q5)))
               alocont_o_nn3d(iip1,jjm1,kkp1,q7) = imask_totreg3d(iip1,jjm1,kkp1) * &
                  (alo_d*c21_scontd + abs_sc*c24_intu*alocont_o_nn3d(iip1,jjm1,kkp1,q4))
               alocont_o_nn3d(i,jjm1,kkp1,q8) = imask_totreg3d(i,jjm1,kkp1) * &
                  (alo_d*c20_scontd + abs_sc*(c23_intu*alocont_o_nn3d(i,jjm1,kkp1,q4) + &
                  c24_intu*alocont_o_nn3d(i,jjm1,kkp1,q5) + &
                  c20_intu*alocont_o_nn3d(i,jjm1,kkp1,q1) + &
                  c26_intu*alocont_o_nn3d(i,jjm1,kkp1,q7)))
               alocont_o_nn3d(iim1,jjm1,kkp1,q9) = imask_totreg3d(iim1,jjm1,kkp1) * &
                  (alo_d*c19_scontd + abs_sc*(c22_intu*alocont_o_nn3d(iim1,jjm1,kkp1,q4) + &
                  c23_intu*alocont_o_nn3d(iim1,jjm1,kkp1,q5) + &
                  c24_intu*alocont_o_nn3d(iim1,jjm1,kkp1,q6) + &
                  c20_intu*alocont_o_nn3d(iim1,jjm1,kkp1,q2) + &
                  c26_intu*alocont_o_nn3d(iim1,jjm1,kkp1,q8)))
               alocont_o_nn3d(iip1,jjp1,k,q10) = imask_totreg3d(iip1,jjp1,k) * &
                  (alo_d*c18_scontd + abs_sc*c18_intu*alocont_o_nn3d(iip1,jjp1,k,q1))
               alocont_o_nn3d(i,jjp1,k,q11) = imask_totreg3d(i,jjp1,k) * &
                  (alo_d*c17_scontd + abs_sc*(c17_intu*alocont_o_nn3d(i,jjp1,k,q1) + &
                  c18_intu*alocont_o_nn3d(i,jjp1,k,q2) + &
                  c26_intu*alocont_o_nn3d(i,jjp1,k,q10)))
               alocont_o_nn3d(iim1,jjp1,k,q12) = imask_totreg3d(iim1,jjp1,k) * &
                  (alo_d*c16_scontd + abs_sc*(c16_intu*alocont_o_nn3d(iim1,jjp1,k,q1) + &
                  c17_intu*alocont_o_nn3d(iim1,jjp1,k,q2) + &
                  c18_intu*alocont_o_nn3d(iim1,jjp1,k,q3) + &
                  c26_intu*alocont_o_nn3d(iim1,jjp1,k,q11)))
               alocont_o_nn3d(iip1,j,k,q13) = imask_totreg3d(iip1,j,k) * &
                  (alo_d*c15_scontd + abs_sc*(c15_intu*alocont_o_nn3d(iip1,j,k,q1) + &
                  c18_intu*alocont_o_nn3d(iip1,j,k,q4) + &
                  c24_intu*alocont_o_nn3d(iip1,j,k,q10)))
               alocont_o_nn3d(i,j,k,q14) = imask_totreg3d(i,j,k) * &
                  (alo_p + abs_sc*(c14_intu*alocont_o_nn3d(i,j,k,q1) + &
                  c15_intu*alocont_o_nn3d(i,j,k,q2) + &
                  c17_intu*alocont_o_nn3d(i,j,k,q4) + &
                  c26_intu*alocont_o_nn3d(i,j,k,q13) + &
                  c18_intu*alocont_o_nn3d(i,j,k,q5) + &
                  c23_intu*alocont_o_nn3d(i,j,k,q10) + &
                  c24_intu*alocont_o_nn3d(i,j,k,q11)))
               alocont_o_nn3d(iim1,j,k,q15) = imask_totreg3d(iim1,j,k) * &
                  (alo_u*c26_scontu + abs_sc*(c13_intu*alocont_o_nn3d(iim1,j,k,q1) + &
                  c14_intu*alocont_o_nn3d(iim1,j,k,q2) + &
                  c15_intu*alocont_o_nn3d(iim1,j,k,q3) + &
                  c16_intu*alocont_o_nn3d(iim1,j,k,q4) + &
                  c17_intu*alocont_o_nn3d(iim1,j,k,q5) + &
                  c18_intu*alocont_o_nn3d(iim1,j,k,q6) + &
                  c22_intu*alocont_o_nn3d(iim1,j,k,q10) + &
                  c23_intu*alocont_o_nn3d(iim1,j,k,q11) + &
                  c24_intu*alocont_o_nn3d(iim1,j,k,q12) + &
                  c26_intu*alocont_o_nn3d(iim1,j,k,q14)))
               alocont_o_nn3d(iip1,jjm1,k,q16) = imask_totreg3d(iip1,jjm1,k) * &
                  (alo_d*c12_scontd + abs_sc*(c12_intu*alocont_o_nn3d(iip1,jjm1,k,q1) + &
                  c15_intu*alocont_o_nn3d(iip1,jjm1,k,q4) + &
                  c18_intu*alocont_o_nn3d(iip1,jjm1,k,q7) + &
                  c24_intu*alocont_o_nn3d(iip1,jjm1,k,q13)))
               alocont_o_nn3d(i,jjm1,k,q17) = imask_totreg3d(i,jjm1,k) * &
                  (alo_u*c24_scontu + abs_sc*(c11_intu*alocont_o_nn3d(i,jjm1,k,q1) + &
                  c12_intu*alocont_o_nn3d(i,jjm1,k,q2) + &
                  c14_intu*alocont_o_nn3d(i,jjm1,k,q4) + &
                  c15_intu*alocont_o_nn3d(i,jjm1,k,q5) + &
                  c17_intu*alocont_o_nn3d(i,jjm1,k,q7) + &
                  c18_intu*alocont_o_nn3d(i,jjm1,k,q8) + &
                  c23_intu*alocont_o_nn3d(i,jjm1,k,q13) + &
                  c24_intu*alocont_o_nn3d(i,jjm1,k,q14) + &
                  c20_intu*alocont_o_nn3d(i,jjm1,k,q10) + &
                  c26_intu*alocont_o_nn3d(i,jjm1,k,q16)))
               alocont_o_nn3d(iim1,jjm1,k,q18) = imask_totreg3d(iim1,jjm1,k) * &
                  (alo_u*c23_scontu + abs_sc*(c10_intu*alocont_o_nn3d(iim1,jjm1,k,q1) + &
                  c11_intu*alocont_o_nn3d(iim1,jjm1,k,q2) + &
                  c12_intu*alocont_o_nn3d(iim1,jjm1,k,q3) + &
                  c13_intu*alocont_o_nn3d(iim1,jjm1,k,q4) + &
                  c14_intu*alocont_o_nn3d(iim1,jjm1,k,q5) + &
                  c15_intu*alocont_o_nn3d(iim1,jjm1,k,q6) + &
                  c16_intu*alocont_o_nn3d(iim1,jjm1,k,q7) + &
                  c17_intu*alocont_o_nn3d(iim1,jjm1,k,q8) + &
                  c18_intu*alocont_o_nn3d(iim1,jjm1,k,q9) + &
                  c22_intu*alocont_o_nn3d(iim1,jjm1,k,q13) + &
                  c23_intu*alocont_o_nn3d(iim1,jjm1,k,q14) + &
                  c24_intu*alocont_o_nn3d(iim1,jjm1,k,q15) + &
                  c20_intu*alocont_o_nn3d(iim1,jjm1,k,q11) + &
                  c26_intu*alocont_o_nn3d(iim1,jjm1,k,q17)))
               alocont_o_nn3d(iip1,jjp1,kkm1,q19) = imask_totreg3d(iip1,jjp1,kkm1) * &
                  (alo_d*c09_scontd + abs_sc*c18_intu*alocont_o_nn3d(iip1,jjp1,kkm1,q10))
               alocont_o_nn3d(i,jjp1,kkm1,q20) = imask_totreg3d(i,jjp1,kkm1) * &
                  (alo_d*c08_scontd + abs_sc*(c17_intu*alocont_o_nn3d(i,jjp1,kkm1,q10) + &
                  c18_intu*alocont_o_nn3d(i,jjp1,kkm1,q11) + &
                  c08_intu*alocont_o_nn3d(i,jjp1,kkm1,q1) + &
                  c26_intu*alocont_o_nn3d(i,jjp1,kkm1,q19)))
               alocont_o_nn3d(iim1,jjp1,kkm1,q21) = imask_totreg3d(iim1,jjp1,kkm1) * &
                  (alo_d*c07_scontd + abs_sc*(c16_intu*alocont_o_nn3d(iim1,jjp1,kkm1,q10) + &
                  c17_intu*alocont_o_nn3d(iim1,jjp1,kkm1,q11) + &
                  c18_intu*alocont_o_nn3d(iim1,jjp1,kkm1,q12) + &
                  c08_intu*alocont_o_nn3d(iim1,jjp1,kkm1,q2) + &
                  c26_intu*alocont_o_nn3d(iim1,jjp1,kkm1,q20)))
               alocont_o_nn3d(iip1,j,kkm1,q22) = imask_totreg3d(iip1,j,kkm1) * &
                  (alo_d*c06_scontd + abs_sc*(c15_intu*alocont_o_nn3d(iip1,j,kkm1,q10) + &
                  c18_intu*alocont_o_nn3d(iip1,j,kkm1,q13) + &
                  c06_intu*alocont_o_nn3d(iip1,j,kkm1,q1) + &
                  c24_intu*alocont_o_nn3d(iip1,j,kkm1,q19)))
               alocont_o_nn3d(i,j,kkm1,q23) = imask_totreg3d(i,j,kkm1) * &
                  (alo_u*c18_scontu + abs_sc*(c14_intu*alocont_o_nn3d(i,j,kkm1,q10) + &
                  c15_intu*alocont_o_nn3d(i,j,kkm1,q11) + &
                  c17_intu*alocont_o_nn3d(i,j,kkm1,q13) + &
                  c18_intu*alocont_o_nn3d(i,j,kkm1,q14) + &
                  c05_intu*alocont_o_nn3d(i,j,kkm1,q1) + &
                  c06_intu*alocont_o_nn3d(i,j,kkm1,q2) + &
                  c23_intu*alocont_o_nn3d(i,j,kkm1,q19) + &
                  c24_intu*alocont_o_nn3d(i,j,kkm1,q20) + &
                  c08_intu*alocont_o_nn3d(i,j,kkm1,q4) + &
                  c26_intu*alocont_o_nn3d(i,j,kkm1,q22)))
               alocont_o_nn3d(iim1,j,kkm1,q24) = imask_totreg3d(iim1,j,kkm1) * &
                  (alo_u*c17_scontu + abs_sc*(c13_intu*alocont_o_nn3d(iim1,j,kkm1,q10) + &
                  c14_intu*alocont_o_nn3d(iim1,j,kkm1,q11) + &
                  c15_intu*alocont_o_nn3d(iim1,j,kkm1,q12) + &
                  c16_intu*alocont_o_nn3d(iim1,j,kkm1,q13) + &
                  c17_intu*alocont_o_nn3d(iim1,j,kkm1,q14) + &
                  c18_intu*alocont_o_nn3d(iim1,j,kkm1,q15) + &
                  c04_intu*alocont_o_nn3d(iim1,j,kkm1,q1) + &
                  c05_intu*alocont_o_nn3d(iim1,j,kkm1,q2) + &
                  c06_intu*alocont_o_nn3d(iim1,j,kkm1,q3) + &
                  c22_intu*alocont_o_nn3d(iim1,j,kkm1,q19) + &
                  c23_intu*alocont_o_nn3d(iim1,j,kkm1,q20) + &
                  c24_intu*alocont_o_nn3d(iim1,j,kkm1,q21) + &
                  c08_intu*alocont_o_nn3d(iim1,j,kkm1,q5) + &
                  c26_intu*alocont_o_nn3d(iim1,j,kkm1,q23)))
               alocont_o_nn3d(iip1,jjm1,kkm1,q25) = imask_totreg3d(iip1,jjm1,kkm1) * &
                  (alo_d*c03_scontd + abs_sc*(c12_intu*alocont_o_nn3d(iip1,jjm1,kkm1,q10) + &
                  c15_intu*alocont_o_nn3d(iip1,jjm1,kkm1,q13) + &
                  c18_intu*alocont_o_nn3d(iip1,jjm1,kkm1,q16) + &
                  c06_intu*alocont_o_nn3d(iip1,jjm1,kkm1,q4) + &
                  c24_intu*alocont_o_nn3d(iip1,jjm1,kkm1,q22)))
               alocont_o_nn3d(i,jjm1,kkm1,q26) = imask_totreg3d(i,jjm1,kkm1) * &
                  (alo_u*c15_scontu + abs_sc*(c11_intu*alocont_o_nn3d(i,jjm1,kkm1,q10) + &
                  c12_intu*alocont_o_nn3d(i,jjm1,kkm1,q11) + &
                  c14_intu*alocont_o_nn3d(i,jjm1,kkm1,q13) + &
                  c15_intu*alocont_o_nn3d(i,jjm1,kkm1,q14) + &
                  c17_intu*alocont_o_nn3d(i,jjm1,kkm1,q16) + &
                  c18_intu*alocont_o_nn3d(i,jjm1,kkm1,q17) + &
                  c05_intu*alocont_o_nn3d(i,jjm1,kkm1,q4) + &
                  c06_intu*alocont_o_nn3d(i,jjm1,kkm1,q5) + &
                  c23_intu*alocont_o_nn3d(i,jjm1,kkm1,q22) + &
                  c24_intu*alocont_o_nn3d(i,jjm1,kkm1,q23) + &
                  c02_intu*alocont_o_nn3d(i,jjm1,kkm1,q1) + &
                  c08_intu*alocont_o_nn3d(i,jjm1,kkm1,q7) + &
                  c20_intu*alocont_o_nn3d(i,jjm1,kkm1,q19) + &
                  c26_intu*alocont_o_nn3d(i,jjm1,kkm1,q25)))
               alocont_o_nn3d(iim1,jjm1,kkm1,q27) = imask_totreg3d(iim1,jjm1,kkm1) * &
                  (alo_u*c14_scontu + abs_sc*(c10_intu*alocont_o_nn3d(iim1,jjm1,kkm1,q10) + &
                  c11_intu*alocont_o_nn3d(iim1,jjm1,kkm1,q11) + &
                  c12_intu*alocont_o_nn3d(iim1,jjm1,kkm1,q12) + &
                  c13_intu*alocont_o_nn3d(iim1,jjm1,kkm1,q13) + &
                  c14_intu*alocont_o_nn3d(iim1,jjm1,kkm1,q14) + &
                  c15_intu*alocont_o_nn3d(iim1,jjm1,kkm1,q15) + &
                  c16_intu*alocont_o_nn3d(iim1,jjm1,kkm1,q16) + &
                  c17_intu*alocont_o_nn3d(iim1,jjm1,kkm1,q17) + &
                  c18_intu*alocont_o_nn3d(iim1,jjm1,kkm1,q18) + &
                  c04_intu*alocont_o_nn3d(iim1,jjm1,kkm1,q4) + &
                  c05_intu*alocont_o_nn3d(iim1,jjm1,kkm1,q5) + &
                  c06_intu*alocont_o_nn3d(iim1,jjm1,kkm1,q6) + &
                  c22_intu*alocont_o_nn3d(iim1,jjm1,kkm1,q22) + &
                  c23_intu*alocont_o_nn3d(iim1,jjm1,kkm1,q23) + &
                  c24_intu*alocont_o_nn3d(iim1,jjm1,kkm1,q24) + &
                  c02_intu*alocont_o_nn3d(iim1,jjm1,kkm1,q2) + &
                  c08_intu*alocont_o_nn3d(iim1,jjm1,kkm1,q8) + &
                  c20_intu*alocont_o_nn3d(iim1,jjm1,kkm1,q20) + &
                  c26_intu*alocont_o_nn3d(iim1,jjm1,kkm1,q26)))
!               call cpu_time(te_aloo)
!               tt_aloo=tt_aloo+te_aloo-ts_aloo

!               if(i.eq.16.and.j.eq.18.and.k.eq.16) then
!                  write(*,*) c01_scontu+c02_scontu+c03_scontu+c04_scontu+c05_scontu+c06_scontu+c07_scontu+c08_scontu+c09_scontu+ &
!                             c10_scontu+c11_scontu+c12_scontu+c13_scontu+c14_scontu+c15_scontu+c16_scontu+c17_scontu+c18_scontu+ &
!                             c19_scontu+c20_scontu+c21_scontu+c22_scontu+c23_scontu+c24_scontu+c25_scontu+c26_scontu+c27_scontu
!                  write(*,*)
!                  write(*,*) c01_scontu,c02_scontu,c03_scontu,c04_scontu,c05_scontu,c06_scontu,c07_scontu,c08_scontu,c09_scontu, &
!                             c10_scontu,c11_scontu,c12_scontu,c13_scontu,c14_scontu,c15_scontu,c16_scontu,c17_scontu,c18_scontu, &
!                             c19_scontu,c20_scontu,c21_scontu,c22_scontu,c23_scontu,c24_scontu,c25_scontu,c26_scontu,c27_scontu
!                  write(*,*)
!                  write(*,*) alocont_o_nn3d(iip1,jjp1,kkp1,q1), alocont_o_nn3d(i,jjp1,kkp1,q2), &
!                  alocont_o_nn3d(iim1,jjp1,kkp1,q3), alocont_o_nn3d(iip1,j,kkp1,q4), &
!                  alocont_o_nn3d(i,j,kkp1,q5), alocont_o_nn3d(iim1,j,kkp1,q6), &
!                  alocont_o_nn3d(iip1,jjm1,kkp1,q7), alocont_o_nn3d(i,jjm1,kkp1,q8), &
!                  alocont_o_nn3d(iim1,jjm1,kkp1,q9), alocont_o_nn3d(iip1,jjp1,k,q10) , &
!                  alocont_o_nn3d(i,jjp1,k,q11), alocont_o_nn3d(iim1,jjp1,k,q12), &
!                  alocont_o_nn3d(iip1,j,k,q13), alocont_o_nn3d(i,j,k,q14), &
!                  alocont_o_nn3d(iim1,j,k,q15), alocont_o_nn3d(iip1,jjm1,k,q16), &
!                  alocont_o_nn3d(i,jjm1,k,q17), alocont_o_nn3d(iim1,jjm1,k,q18), &
!                  alocont_o_nn3d(iip1,jjp1,kkm1,q19), alocont_o_nn3d(i,jjp1,kkm1,q20), &
!                  alocont_o_nn3d(iim1,jjp1,kkm1,q21), alocont_o_nn3d(iip1,j,kkm1,q22), &
!                  alocont_o_nn3d(i,j,kkm1,q23), alocont_o_nn3d(iim1,j,kkm1,q24), &
!                  alocont_o_nn3d(iip1,jjm1,kkm1,q25), alocont_o_nn3d(i,jjm1,kkm1,q26), &
!                  alocont_o_nn3d(iim1,jjm1,kkm1,q27)
!               endif
!
!perform angular integration
!*********debug start
!               dels_r=dist_sphere(nn_x,nn_y,nn_z,x(i),y(j),z(k))
!               if(dels_r.ne.1.d10) then
!                  int3d(i,j,k)=xic1
!               else
!                  int3d(i,j,k)=0.d0
!               endif
!*********debug end
!
!               call cpu_time(ts_integ)
               mint3d_tmp(i,j,k) = mint3d_tmp(i,j,k) + int3d(i,j,k)*wall
               fcontx3d_tmp(i,j,k) = fcontx3d_tmp(i,j,k) + int3d(i,j,k)*nn_x*wall
               fconty3d_tmp(i,j,k) = fconty3d_tmp(i,j,k) + int3d(i,j,k)*nn_y*wall
               fcontz3d_tmp(i,j,k) = fcontz3d_tmp(i,j,k) + int3d(i,j,k)*nn_z*wall
               kcontxx3d_tmp(i,j,k) = kcontxx3d_tmp(i,j,k) + int3d(i,j,k)*nn_x*nn_x*wall
               kcontyy3d_tmp(i,j,k) = kcontyy3d_tmp(i,j,k) + int3d(i,j,k)*nn_y*nn_y*wall
               kcontzz3d_tmp(i,j,k) = kcontzz3d_tmp(i,j,k) + int3d(i,j,k)*nn_z*nn_z*wall
               kcontxy3d_tmp(i,j,k) = kcontxy3d_tmp(i,j,k) + int3d(i,j,k)*nn_x*nn_y*wall
               kcontxz3d_tmp(i,j,k) = kcontxz3d_tmp(i,j,k) + int3d(i,j,k)*nn_x*nn_z*wall
               kcontyz3d_tmp(i,j,k) = kcontyz3d_tmp(i,j,k) + int3d(i,j,k)*nn_y*nn_z*wall
!
               normalization3d_tmp(i,j,k) = normalization3d_tmp(i,j,k) + wall
!
               alocont_nn3d_tmp(iip1,jjp1,kkp1,q1) = alocont_nn3d_tmp(iip1,jjp1,kkp1,q1) + wall*alocont_o_nn3d(iip1,jjp1,kkp1,q1)
               alocont_nn3d_tmp(i,jjp1,kkp1,q2) = alocont_nn3d_tmp(i,jjp1,kkp1,q2) + wall*alocont_o_nn3d(i,jjp1,kkp1,q2)
               alocont_nn3d_tmp(iim1,jjp1,kkp1,q3) = alocont_nn3d_tmp(iim1,jjp1,kkp1,q3) + wall*alocont_o_nn3d(iim1,jjp1,kkp1,q3)
               alocont_nn3d_tmp(iip1,j,kkp1,q4) = alocont_nn3d_tmp(iip1,j,kkp1,q4) + wall*alocont_o_nn3d(iip1,j,kkp1,q4)
               alocont_nn3d_tmp(i,j,kkp1,q5) = alocont_nn3d_tmp(i,j,kkp1,q5) + wall*alocont_o_nn3d(i,j,kkp1,q5)
               alocont_nn3d_tmp(iim1,j,kkp1,q6) = alocont_nn3d_tmp(iim1,j,kkp1,q6) + wall*alocont_o_nn3d(iim1,j,kkp1,q6)
               alocont_nn3d_tmp(iip1,jjm1,kkp1,q7) = alocont_nn3d_tmp(iip1,jjm1,kkp1,q7) + wall*alocont_o_nn3d(iip1,jjm1,kkp1,q7)
               alocont_nn3d_tmp(i,jjm1,kkp1,q8) = alocont_nn3d_tmp(i,jjm1,kkp1,q8) + wall*alocont_o_nn3d(i,jjm1,kkp1,q8)
               alocont_nn3d_tmp(iim1,jjm1,kkp1,q9) = alocont_nn3d_tmp(iim1,jjm1,kkp1,q9) + wall*alocont_o_nn3d(iim1,jjm1,kkp1,q9)
               alocont_nn3d_tmp(iip1,jjp1,k,q10) = alocont_nn3d_tmp(iip1,jjp1,k,q10) + wall*alocont_o_nn3d(iip1,jjp1,k,q10)
               alocont_nn3d_tmp(i,jjp1,k,q11) = alocont_nn3d_tmp(i,jjp1,k,q11) + wall*alocont_o_nn3d(i,jjp1,k,q11)
               alocont_nn3d_tmp(iim1,jjp1,k,q12) = alocont_nn3d_tmp(iim1,jjp1,k,q12) + wall*alocont_o_nn3d(iim1,jjp1,k,q12)
               alocont_nn3d_tmp(iip1,j,k,q13) = alocont_nn3d_tmp(iip1,j,k,q13) + wall*alocont_o_nn3d(iip1,j,k,q13)
               alocont_nn3d_tmp(i,j,k,q14) = alocont_nn3d_tmp(i,j,k,q14) + wall*alocont_o_nn3d(i,j,k,q14)
               alocont_nn3d_tmp(iim1,j,k,q15) = alocont_nn3d_tmp(iim1,j,k,q15) + wall*alocont_o_nn3d(iim1,j,k,q15)
               alocont_nn3d_tmp(iip1,jjm1,k,q16) = alocont_nn3d_tmp(iip1,jjm1,k,q16) + wall*alocont_o_nn3d(iip1,jjm1,k,q16)
               alocont_nn3d_tmp(i,jjm1,k,q17) = alocont_nn3d_tmp(i,jjm1,k,q17) + wall*alocont_o_nn3d(i,jjm1,k,q17)
               alocont_nn3d_tmp(iim1,jjm1,k,q18) = alocont_nn3d_tmp(iim1,jjm1,k,q18) + wall*alocont_o_nn3d(iim1,jjm1,k,q18)
               alocont_nn3d_tmp(iip1,jjp1,kkm1,q19) = alocont_nn3d_tmp(iip1,jjp1,kkm1,q19) + wall*alocont_o_nn3d(iip1,jjp1,kkm1,q19)
               alocont_nn3d_tmp(i,jjp1,kkm1,q20) = alocont_nn3d_tmp(i,jjp1,kkm1,q20) + wall*alocont_o_nn3d(i,jjp1,kkm1,q20)
               alocont_nn3d_tmp(iim1,jjp1,kkm1,q21) = alocont_nn3d_tmp(iim1,jjp1,kkm1,q21) + wall*alocont_o_nn3d(iim1,jjp1,kkm1,q21)
               alocont_nn3d_tmp(iip1,j,kkm1,q22) = alocont_nn3d_tmp(iip1,j,kkm1,q22) + wall*alocont_o_nn3d(iip1,j,kkm1,q22)
               alocont_nn3d_tmp(i,j,kkm1,q23) = alocont_nn3d_tmp(i,j,kkm1,q23) + wall*alocont_o_nn3d(i,j,kkm1,q23)
               alocont_nn3d_tmp(iim1,j,kkm1,q24) = alocont_nn3d_tmp(iim1,j,kkm1,q24) + wall*alocont_o_nn3d(iim1,j,kkm1,q24)
               alocont_nn3d_tmp(iip1,jjm1,kkm1,q25) = alocont_nn3d_tmp(iip1,jjm1,kkm1,q25) + wall*alocont_o_nn3d(iip1,jjm1,kkm1,q25)
               alocont_nn3d_tmp(i,jjm1,kkm1,q26) = alocont_nn3d_tmp(i,jjm1,kkm1,q26) + wall*alocont_o_nn3d(i,jjm1,kkm1,q26)
               alocont_nn3d_tmp(iim1,jjm1,kkm1,q27) = alocont_nn3d_tmp(iim1,jjm1,kkm1,q27) + wall*alocont_o_nn3d(iim1,jjm1,kkm1,q27)
!               call cpu_time(te_integ)
!               tt_integ=tt_integ+te_integ-ts_integ
!
!               call cpu_time(te_case1)
!               tt_case1=tt_case1+te_case1-ts_case1
!
!***********special treatment: point is in vicinity of boundary*********
!
             case(2)
!
!               call cpu_time(ts_case2)
!
               iim1=i-alpha
               jjm1=j-beta
               kkm1=k-gamma
               iip1=i+alpha
               jjp1=j+beta
               kkp1=k+gamma
!
!calculate distance of upwind point to (previous) xy-plane, xz-plane and yz-plane
               dels_xyu=(z(k)-z(kkm1))/nn_z
               dels_xzu=(y(j)-y(jjm1))/nn_y
               dels_yzu=(x(i)-x(iim1))/nn_x
!
!calculate distance of downwind point to (next) xy-plane, xz-plane and yz-plane
               dels_xyd=(z(kkp1)-z(k))/nn_z
               dels_xzd=(y(jjp1)-y(j))/nn_y
               dels_yzd=(x(iip1)-x(i))/nn_x
!
!calculate distance to the boundary (here: sphere)
!               dels_r=dist_sphere(nn_x,nn_y,nn_z,x(i),y(j),z(k))
               dels_r=dist_ellipsoid(nn_x,nn_y,nn_z,x(i),y(j),z(k),smajorax_a,smajorax_b,smajorax_c)
!
               dels_u=min(dels_xyu,dels_xzu,dels_yzu,dels_r)
               dels_d=min(dels_xyd,dels_xzd,dels_yzd)
!
!----------------------------local point--------------------------------
!
               scont_p=scont3d(i,j,k)
               opac_p=opac3d(i,j,k)
!
!----------------------------upwind point-------------------------------
!
               if(dels_r.eq.dels_u) then
                  x_u = x(i) - dels_u*nn_x
                  y_u = y(j) - dels_u*nn_y
                  z_u = z(k) - dels_u*nn_z
!
                  iim2=i-2*alpha
                  jjm2=j-2*beta
                  kkm2=k-2*gamma
!
!for alo
!                  int_u=0.d0
!                  int_u=xic1
                  int_u=calc_icore_gdark(z_u, sqrt(x_u**2+y_u**2+z_u**2))
!
!                  call cpu_time(ts_interpu)
                  call coeff3d_contu(opac3d(iim2,jjm2,kkm2), opac3d(iim1,jjm2,kkm2), opac3d(i,jjm2,kkm2), &
                     opac3d(iim2,jjm1,kkm2), opac3d(iim1,jjm1,kkm2), opac3d(i,jjm1,kkm2), &
                     opac3d(iim2,j,kkm2),    opac3d(iim1,j,kkm2),    opac3d(i,j,kkm2), &
                     opac3d(iim2,jjm2,kkm1), opac3d(iim1,jjm2,kkm1), opac3d(i,jjm2,kkm1), &
                     opac3d(iim2,jjm1,kkm1), opac3d(iim1,jjm1,kkm1), opac3d(i,jjm1,kkm1), &
                     opac3d(iim2,j,kkm1),    opac3d(iim1,j,kkm1),    opac3d(i,j,kkm1), &
                     opac3d(iim2,jjm2,k),    opac3d(iim1,jjm2,k),    opac3d(i,jjm2,k), &
                     opac3d(iim2,jjm1,k),    opac3d(iim1,jjm1,k),    opac3d(i,jjm1,k), &
                     opac3d(iim2,j,k),       opac3d(iim1,j,k),       opac3d(i,j,k), &
                     scont3d(iim2,jjm2,kkm2), scont3d(iim1,jjm2,kkm2), scont3d(i,jjm2,kkm2), &
                     scont3d(iim2,jjm1,kkm2), scont3d(iim1,jjm1,kkm2), scont3d(i,jjm1,kkm2), &
                     scont3d(iim2,j,kkm2),    scont3d(iim1,j,kkm2),    scont3d(i,j,kkm2), &
                     scont3d(iim2,jjm2,kkm1), scont3d(iim1,jjm2,kkm1), scont3d(i,jjm2,kkm1), &
                     scont3d(iim2,jjm1,kkm1), scont3d(iim1,jjm1,kkm1), scont3d(i,jjm1,kkm1), &
                     scont3d(iim2,j,kkm1),    scont3d(iim1,j,kkm1),    scont3d(i,j,kkm1), &
                     scont3d(iim2,jjm2,k),    scont3d(iim1,jjm2,k),    scont3d(i,jjm2,k), &
                     scont3d(iim2,jjm1,k),    scont3d(iim1,jjm1,k),    scont3d(i,jjm1,k), &
                     scont3d(iim2,j,k),       scont3d(iim1,j,k),       scont3d(i,j,k), &
                     x(iim2), x(iim1), x(i), y(jjm2), y(jjm1), y(j), z(kkm2), z(kkm1), z(k), x_u, y_u, z_u, &
                     c01_scontu, c02_scontu, c03_scontu, c04_scontu, c05_scontu, c06_scontu, c07_scontu, c08_scontu, c09_scontu, &
                     c10_scontu, c11_scontu, c12_scontu, c13_scontu, c14_scontu, c15_scontu, c16_scontu, c17_scontu, c18_scontu, &
                     c19_scontu, c20_scontu, c21_scontu, c22_scontu, c23_scontu, c24_scontu, c25_scontu, c26_scontu, c27_scontu, &
                     opac_u, scont_u)
!                  call cpu_time(te_interpu)
!                  tt_interpu=tt_interpu+te_interpu-ts_interpu
!
!set interpolation coefficients that are not used to zero
                  c02_intu = 0.d0
                  c04_intu = 0.d0
                  c05_intu = 0.d0
                  c06_intu = 0.d0
                  c08_intu = 0.d0
                  c10_intu = 0.d0
                  c11_intu = 0.d0
                  c12_intu = 0.d0
                  c13_intu = 0.d0
                  c14_intu = 0.d0
                  c15_intu = 0.d0
                  c16_intu = 0.d0
                  c17_intu = 0.d0
                  c18_intu = 0.d0
                  c20_intu = 0.d0
                  c22_intu = 0.d0
                  c23_intu = 0.d0
                  c24_intu = 0.d0
                  c26_intu = 0.d0
!
               elseif(dels_xyu.eq.dels_u) then
!intersection with x-y plane on level k-gamma
                  x_u = x(i) - dels_u*nn_x
                  y_u = y(j) - dels_u*nn_y
                  z_u = z(kkm1)
!
                  iim2=i-2*alpha
                  jjm2=j-2*beta
!
!                  call cpu_time(ts_interpu)
                  call coeff2d_contu(opac3d(iim2,jjm2,kkm1), opac3d(iim1,jjm2,kkm1), opac3d(i,jjm2,kkm1), &
                     opac3d(iim2,jjm1,kkm1), opac3d(iim1,jjm1,kkm1), opac3d(i,jjm1,kkm1), &
                     opac3d(iim2,j,kkm1), opac3d(iim1,j,kkm1), opac3d(i,j,kkm1), &
                     scont3d(iim2,jjm2,kkm1), scont3d(iim1,jjm2,kkm1), scont3d(i,jjm2,kkm1), &
                     scont3d(iim2,jjm1,kkm1), scont3d(iim1,jjm1,kkm1), scont3d(i,jjm1,kkm1), &
                     scont3d(iim2,j,kkm1), scont3d(iim1,j,kkm1), scont3d(i,j,kkm1), &
                     int3d(iim2,jjm2,kkm1), int3d(iim1,jjm2,kkm1), int3d(i,jjm2,kkm1), &
                     int3d(iim2,jjm1,kkm1), int3d(iim1,jjm1,kkm1), int3d(i,jjm1,kkm1), &
                     int3d(iim2,j,kkm1), int3d(iim1,j,kkm1), int3d(i,j,kkm1), &
                     x(iim2), x(iim1), x(i), y(jjm2), y(jjm1), y(j), x_u, y_u, &
                     c10_scontu, c11_scontu, c12_scontu, c13_scontu, c14_scontu, &
                     c15_scontu, c16_scontu, c17_scontu, c18_scontu, &
                     c10_intu, c11_intu, c12_intu, c13_intu, c14_intu, &
                     c15_intu, c16_intu, c17_intu, c18_intu, opac_u, scont_u, int_u)
!                  call cpu_time(te_interpu)
!                  tt_interpu=tt_interpu+te_interpu-ts_interpu
!
!set interpolation coefficients that are not used to zero
                  c23_scontu = 0.d0
                  c24_scontu = 0.d0
                  c26_scontu = 0.d0
                  c27_scontu = 0.d0
!
                  c02_intu = 0.d0
                  c04_intu = 0.d0
                  c05_intu = 0.d0
                  c06_intu = 0.d0
                  c08_intu = 0.d0
                  c20_intu = 0.d0
                  c22_intu = 0.d0
                  c23_intu = 0.d0
                  c24_intu = 0.d0
                  c26_intu = 0.d0
!
               elseif(dels_xzu.eq.dels_u) then
!intersection with x-z plane on level j-beta
                  x_u = x(i) - dels_u*nn_x
                  y_u = y(jjm1)
                  z_u = z(k) - dels_u*nn_z
!
                  iim2=i-2*alpha
                  kkm2=k-2*gamma
!
!                  call cpu_time(ts_interpu)
                  call coeff2d_contu(opac3d(iim2,jjm1,kkm2), opac3d(iim1,jjm1,kkm2), opac3d(i,jjm1,kkm2), &
                     opac3d(iim2,jjm1,kkm1), opac3d(iim1,jjm1,kkm1), opac3d(i,jjm1,kkm1), &
                     opac3d(iim2,jjm1,k), opac3d(iim1,jjm1,k), opac3d(i,jjm1,k), &
                     scont3d(iim2,jjm1,kkm2), scont3d(iim1,jjm1,kkm2), scont3d(i,jjm1,kkm2), &
                     scont3d(iim2,jjm1,kkm1), scont3d(iim1,jjm1,kkm1), scont3d(i,jjm1,kkm1), &
                     scont3d(iim2,jjm1,k), scont3d(iim1,jjm1,k), scont3d(i,jjm1,k), &
                     int3d(iim2,jjm1,kkm2), int3d(iim1,jjm1,kkm2), int3d(i,jjm1,kkm2), &
                     int3d(iim2,jjm1,kkm1), int3d(iim1,jjm1,kkm1), int3d(i,jjm1,kkm1), &
                     int3d(iim2,jjm1,k), int3d(iim1,jjm1,k), int3d(i,jjm1,k), &
                     x(iim2), x(iim1), x(i), z(kkm2), z(kkm1), z(k), x_u, z_u, &
                     c04_scontu, c05_scontu, c06_scontu, c13_scontu, c14_scontu, &
                     c15_scontu, c22_scontu, c23_scontu, c24_scontu, &
                     c04_intu, c05_intu, c06_intu, c13_intu, c14_intu, &
                     c15_intu, c22_intu, c23_intu, c24_intu, opac_u, scont_u, int_u)
!                  call cpu_time(te_interpu)
!                  tt_interpu=tt_interpu+te_interpu-ts_interpu
!
!set interpolation coefficients that are not used to zero
                  c17_scontu = 0.d0
                  c18_scontu = 0.d0
                  c26_scontu = 0.d0
                  c27_scontu = 0.d0
!
                  c02_intu = 0.d0
                  c08_intu = 0.d0
                  c10_intu = 0.d0
                  c11_intu = 0.d0
                  c12_intu = 0.d0
                  c16_intu = 0.d0
                  c17_intu = 0.d0
                  c18_intu = 0.d0
                  c20_intu = 0.d0
                  c26_intu = 0.d0
!
               elseif(dels_yzu.eq.dels_u) then
!intersection with y-z plane on level iim1
                  x_u = x(iim1)
                  y_u = y(j) - dels_u*nn_y
                  z_u = z(k) - dels_u*nn_z
!
                  jjm2=j-2*beta
                  kkm2=k-2*gamma
!
!                  call cpu_time(ts_interpu)
                  call coeff2d_contu(opac3d(iim1,jjm2,kkm2), opac3d(iim1,jjm1,kkm2), opac3d(iim1,j,kkm2), &
                     opac3d(iim1,jjm2,kkm1), opac3d(iim1,jjm1,kkm1), opac3d(iim1,j,kkm1), &
                     opac3d(iim1,jjm2,k), opac3d(iim1,jjm1,k), opac3d(iim1,j,k), &
                     scont3d(iim1,jjm2,kkm2), scont3d(iim1,jjm1,kkm2), scont3d(iim1,j,kkm2), &
                     scont3d(iim1,jjm2,kkm1), scont3d(iim1,jjm1,kkm1), scont3d(iim1,j,kkm1), &
                     scont3d(iim1,jjm2,k), scont3d(iim1,jjm1,k), scont3d(iim1,j,k), &
                     int3d(iim1,jjm2,kkm2), int3d(iim1,jjm1,kkm2), int3d(iim1,j,kkm2), &
                     int3d(iim1,jjm2,kkm1), int3d(iim1,jjm1,kkm1), int3d(iim1,j,kkm1), &
                     int3d(iim1,jjm2,k), int3d(iim1,jjm1,k), int3d(iim1,j,k), &
                     y(jjm2), y(jjm1), y(j), z(kkm2), z(kkm1), z(k), y_u, z_u, &
                     c02_scontu, c05_scontu, c08_scontu, c11_scontu, c14_scontu, &
                     c17_scontu, c20_scontu, c23_scontu, c26_scontu, &
                     c02_intu, c05_intu, c08_intu, c11_intu, c14_intu, &
                     c17_intu, c20_intu, c23_intu, c26_intu, opac_u, scont_u, int_u)
!                  call cpu_time(te_interpu)
!                  tt_interpu=tt_interpu+te_interpu-ts_interpu
!
!set interpolation coefficients that are not used to zero
                  c15_scontu = 0.d0
                  c18_scontu = 0.d0
                  c24_scontu = 0.d0
                  c27_scontu = 0.d0
!
                  c04_intu = 0.d0
                  c06_intu = 0.d0
                  c10_intu = 0.d0
                  c12_intu = 0.d0
                  c13_intu = 0.d0
                  c15_intu = 0.d0
                  c16_intu = 0.d0
                  c18_intu = 0.d0
                  c22_intu = 0.d0
                  c24_intu = 0.d0
!
               else
                  write(*,'(4es20.8,2l4)') dels_u, dels_xzu, dels_xyu, dels_yzu, dels_u.eq.dels_xzu, dels_xyu.eq.dels_u
                  stop 'error in fsc_cont3d: invalid dels_u'
               endif
!
!---------------------------downwind point------------------------------
!
               if(dels_xyd.eq.dels_d) then
!intersection with x-y plane on level k+gamma
                  x_d = x(i) + dels_d*nn_x
                  y_d = y(j) + dels_d*nn_y
                  z_d = z(kkp1)
!
!                  call cpu_time(ts_interpd)
                  call coeff2d_contd(opac3d(iim1,jjm1,kkp1), opac3d(i,jjm1,kkp1), opac3d(iip1,jjm1,kkp1), &
                     opac3d(iim1,j,kkp1), opac3d(i,j,kkp1), opac3d(iip1,j,kkp1), &
                     opac3d(iim1,jjp1,kkp1), opac3d(i,jjp1,kkp1), opac3d(iip1,jjp1,kkp1), &
                     scont3d(iim1,jjm1,kkp1), scont3d(i,jjm1,kkp1), scont3d(iip1,jjm1,kkp1), &
                     scont3d(iim1,j,kkp1), scont3d(i,j,kkp1), scont3d(iip1,j,kkp1), &
                     scont3d(iim1,jjp1,kkp1), scont3d(i,jjp1,kkp1), scont3d(iip1,jjp1,kkp1), &
                     x(iim1), x(i), x(iip1), y(jjm1), y(j), y(jjp1), x_d, y_d, &
                     c19_scontd, c20_scontd, c21_scontd, c22_scontd, c23_scontd, c24_scontd, &
                     c25_scontd, c26_scontd, c27_scontd, opac_d, scont_d)
!                  call cpu_time(te_interpd)
!                  tt_interpd=tt_interpd+te_interpd-ts_interpd
!
!set interpolation coefficients that are not used to zero
                  c03_scontd = 0.d0
                  c06_scontd = 0.d0
                  c07_scontd = 0.d0
                  c08_scontd = 0.d0
                  c09_scontd = 0.d0
                  c12_scontd = 0.d0
                  c15_scontd = 0.d0
                  c16_scontd = 0.d0
                  c17_scontd = 0.d0
                  c18_scontd = 0.d0
!
               elseif(dels_xzd.eq.dels_d) then
!intersection with x-z plane on level j+beta
                  x_d = x(i) + dels_d*nn_x
                  y_d = y(jjp1)
                  z_d = z(k) + dels_d*nn_z
!
!                  call cpu_time(ts_interpd)
                  call coeff2d_contd(opac3d(iim1,jjp1,kkm1), opac3d(i,jjp1,kkm1), opac3d(iip1,jjp1,kkm1), &
                     opac3d(iim1,jjp1,k), opac3d(i,jjp1,k), opac3d(iip1,jjp1,k), &
                     opac3d(iim1,jjp1,kkp1), opac3d(i,jjp1,kkp1), opac3d(iip1,jjp1,kkp1), &
                     scont3d(iim1,jjp1,kkm1), scont3d(i,jjp1,kkm1), scont3d(iip1,jjp1,kkm1), &
                     scont3d(iim1,jjp1,k), scont3d(i,jjp1,k), scont3d(iip1,jjp1,k), &
                     scont3d(iim1,jjp1,kkp1), scont3d(i,jjp1,kkp1), scont3d(iip1,jjp1,kkp1), &
                     x(iim1), x(i), x(iip1), z(kkm1), z(k), z(kkp1), x_d, z_d, &
                     c07_scontd, c08_scontd, c09_scontd, c16_scontd, c17_scontd, c18_scontd, &
                     c25_scontd, c26_scontd, c27_scontd, opac_d, scont_d)
!                  call cpu_time(te_interpd)
!                  tt_interpd=tt_interpd+te_interpd-ts_interpd
!
!set interpolation coefficients that are not used to zero
                  c03_scontd = 0.d0
                  c06_scontd = 0.d0
                  c12_scontd = 0.d0
                  c15_scontd = 0.d0
                  c19_scontd = 0.d0
                  c20_scontd = 0.d0
                  c21_scontd = 0.d0
                  c22_scontd = 0.d0
                  c23_scontd = 0.d0
                  c24_scontd = 0.d0
!
               elseif(dels_yzd.eq.dels_d) then
!intersection with y-z plane on level i+alpha
                  x_d = x(iip1)
                  y_d = y(j) + dels_d*nn_y
                  z_d = z(k) + dels_d*nn_z
!
!                  call cpu_time(ts_interpd)
                  call coeff2d_contd(opac3d(iip1,jjm1,kkm1), opac3d(iip1,j,kkm1), opac3d(iip1,jjp1,kkm1), &
                     opac3d(iip1,jjm1,k), opac3d(iip1,j,k), opac3d(iip1,jjp1,k), &
                     opac3d(iip1,jjm1,kkp1), opac3d(iip1,j,kkp1), opac3d(iip1,jjp1,kkp1), &
                     scont3d(iip1,jjm1,kkm1), scont3d(iip1,j,kkm1), scont3d(iip1,jjp1,kkm1), &
                     scont3d(iip1,jjm1,k), scont3d(iip1,j,k), scont3d(iip1,jjp1,k), &
                     scont3d(iip1,jjm1,kkp1), scont3d(iip1,j,kkp1), scont3d(iip1,jjp1,kkp1), &
                     y(jjm1), y(j), y(jjp1), z(kkm1), z(k), z(kkp1), y_d, z_d, &
                     c03_scontd, c06_scontd, c09_scontd, c12_scontd, c15_scontd, c18_scontd, &
                     c21_scontd, c24_scontd, c27_scontd, opac_d, scont_d)
!                  call cpu_time(te_interpd)
!                  tt_interpd=tt_interpd+te_interpd-ts_interpd
!
!set interpolation coefficients that are not used to zero
                  c07_scontd = 0.d0
                  c08_scontd = 0.d0
                  c16_scontd = 0.d0
                  c17_scontd = 0.d0
                  c19_scontd = 0.d0
                  c20_scontd = 0.d0
                  c22_scontd = 0.d0
                  c23_scontd = 0.d0
                  c25_scontd = 0.d0
                  c26_scontd = 0.d0
!
               else
                  write(*,'(4es20.8,l4)') dels_d, dels_xzd, dels_xyd, dels_yzd, dels_d.eq.dels_xzd
                  stop 'error in fsc_cont3d: invalid dels_d'
               endif
!
!--------------------------------radiative transfer---------------------
!
!               call cpu_time(ts_fs1d)
               call fsc_cont(int_u, opac_u, opac_p, opac_d, scont_u, scont_p, scont_d, &
                  dels_u, dels_d, abs_sc, contr_sc, int_sc, alo_u, alo_p, alo_d)
               int3d(i,j,k) = abs_sc*int_u + alo_u*scont_u + alo_p*scont_p + alo_d*scont_d!int_sc
!               call cpu_time(te_fs1d)
!               tt_fs1d=tt_fs1d+te_fs1d-ts_fs1d
!
!               call cpu_time(ts_aloo)
               alocont_o_nn3d(iip1,jjp1,kkp1,q1) = imask_totreg3d(iip1,jjp1,kkp1)*alo_d*c27_scontd
               alocont_o_nn3d(i,jjp1,kkp1,q2) = imask_totreg3d(i,jjp1,kkp1) * &
                  (alo_d*c26_scontd + abs_sc*c26_intu*alocont_o_nn3d(i,jjp1,kkp1,q1))
               alocont_o_nn3d(iim1,jjp1,kkp1,q3) = imask_totreg3d(iim1,jjp1,kkp1) * &
                  (alo_d*c25_scontd + abs_sc*c26_intu*alocont_o_nn3d(iim1,jjp1,kkp1,q2))
               alocont_o_nn3d(iip1,j,kkp1,q4) = imask_totreg3d(iip1,j,kkp1) * &
                  (alo_d*c24_scontd + abs_sc*c24_intu*alocont_o_nn3d(iip1,j,kkp1,q1))
               alocont_o_nn3d(i,j,kkp1,q5) = imask_totreg3d(i,j,kkp1) * &
                  (alo_d*c23_scontd + abs_sc*(c23_intu*alocont_o_nn3d(i,j,kkp1,q1) + &
                  c24_intu*alocont_o_nn3d(i,j,kkp1,q2) + &
                  c26_intu*alocont_o_nn3d(i,j,kkp1,q4)))
               alocont_o_nn3d(iim1,j,kkp1,q6) = imask_totreg3d(iim1,j,kkp1) * &
                  (alo_d*c22_scontd + abs_sc*(c22_intu*alocont_o_nn3d(iim1,j,kkp1,q1) + &
                  c23_intu*alocont_o_nn3d(iim1,j,kkp1,q2) + &
                  c24_intu*alocont_o_nn3d(iim1,j,kkp1,q3) + &
                  c26_intu*alocont_o_nn3d(iim1,j,kkp1,q5)))
               alocont_o_nn3d(iip1,jjm1,kkp1,q7) = imask_totreg3d(iip1,jjm1,kkp1) * &
                  (alo_d*c21_scontd + abs_sc*c24_intu*alocont_o_nn3d(iip1,jjm1,kkp1,q4))
               alocont_o_nn3d(i,jjm1,kkp1,q8) = imask_totreg3d(i,jjm1,kkp1) * &
                  (alo_d*c20_scontd + abs_sc*(c23_intu*alocont_o_nn3d(i,jjm1,kkp1,q4) + &
                  c24_intu*alocont_o_nn3d(i,jjm1,kkp1,q5) + &
                  c20_intu*alocont_o_nn3d(i,jjm1,kkp1,q1) + &
                  c26_intu*alocont_o_nn3d(i,jjm1,kkp1,q7)))
               alocont_o_nn3d(iim1,jjm1,kkp1,q9) = imask_totreg3d(iim1,jjm1,kkp1) * &
                  (alo_d*c19_scontd + abs_sc*(c22_intu*alocont_o_nn3d(iim1,jjm1,kkp1,q4) + &
                  c23_intu*alocont_o_nn3d(iim1,jjm1,kkp1,q5) + &
                  c24_intu*alocont_o_nn3d(iim1,jjm1,kkp1,q6) + &
                  c20_intu*alocont_o_nn3d(iim1,jjm1,kkp1,q2) + &
                  c26_intu*alocont_o_nn3d(iim1,jjm1,kkp1,q8)))
               alocont_o_nn3d(iip1,jjp1,k,q10) = imask_totreg3d(iip1,jjp1,k) * &
                  (alo_d*c18_scontd + abs_sc*c18_intu*alocont_o_nn3d(iip1,jjp1,k,q1))
               alocont_o_nn3d(i,jjp1,k,q11) = imask_totreg3d(i,jjp1,k) * &
                  (alo_d*c17_scontd + abs_sc*(c17_intu*alocont_o_nn3d(i,jjp1,k,q1) + &
                  c18_intu*alocont_o_nn3d(i,jjp1,k,q2) + &
                  c26_intu*alocont_o_nn3d(i,jjp1,k,q10)))
               alocont_o_nn3d(iim1,jjp1,k,q12) = imask_totreg3d(iim1,jjp1,k) * &
                  (alo_d*c16_scontd + abs_sc*(c16_intu*alocont_o_nn3d(iim1,jjp1,k,q1) + &
                  c17_intu*alocont_o_nn3d(iim1,jjp1,k,q2) + &
                  c18_intu*alocont_o_nn3d(iim1,jjp1,k,q3) + &
                  c26_intu*alocont_o_nn3d(iim1,jjp1,k,q11)))
               alocont_o_nn3d(iip1,j,k,q13) = imask_totreg3d(iip1,j,k) * &
                  (alo_d*c15_scontd + abs_sc*(c15_intu*alocont_o_nn3d(iip1,j,k,q1) + &
                  c18_intu*alocont_o_nn3d(iip1,j,k,q4) + &
                  c24_intu*alocont_o_nn3d(iip1,j,k,q10)))
               alocont_o_nn3d(i,j,k,q14) = imask_totreg3d(i,j,k) * &
                  (alo_p + alo_u*c27_scontu + abs_sc*(c14_intu*alocont_o_nn3d(i,j,k,q1) + &
                  c15_intu*alocont_o_nn3d(i,j,k,q2) + &
                  c17_intu*alocont_o_nn3d(i,j,k,q4) + &
                  c26_intu*alocont_o_nn3d(i,j,k,q13) + &
                  c18_intu*alocont_o_nn3d(i,j,k,q5) + &
                  c23_intu*alocont_o_nn3d(i,j,k,q10) + &
                  c24_intu*alocont_o_nn3d(i,j,k,q11)))
               alocont_o_nn3d(iim1,j,k,q15) = imask_totreg3d(iim1,j,k) * &
                  (alo_u*c26_scontu + abs_sc*(c13_intu*alocont_o_nn3d(iim1,j,k,q1) + &
                  c14_intu*alocont_o_nn3d(iim1,j,k,q2) + &
                  c15_intu*alocont_o_nn3d(iim1,j,k,q3) + &
                  c16_intu*alocont_o_nn3d(iim1,j,k,q4) + &
                  c17_intu*alocont_o_nn3d(iim1,j,k,q5) + &
                  c18_intu*alocont_o_nn3d(iim1,j,k,q6) + &
                  c22_intu*alocont_o_nn3d(iim1,j,k,q10) + &
                  c23_intu*alocont_o_nn3d(iim1,j,k,q11) + &
                  c24_intu*alocont_o_nn3d(iim1,j,k,q12) + &
                  c26_intu*alocont_o_nn3d(iim1,j,k,q14)))
               alocont_o_nn3d(iip1,jjm1,k,q16) = imask_totreg3d(iip1,jjm1,k) * &
                  (alo_d*c12_scontd + abs_sc*(c12_intu*alocont_o_nn3d(iip1,jjm1,k,q1) + &
                  c15_intu*alocont_o_nn3d(iip1,jjm1,k,q4) + &
                  c18_intu*alocont_o_nn3d(iip1,jjm1,k,q7) + &
                  c24_intu*alocont_o_nn3d(iip1,jjm1,k,q13)))
               alocont_o_nn3d(i,jjm1,k,q17) = imask_totreg3d(i,jjm1,k) * &
                  (alo_u*c24_scontu + abs_sc*(c11_intu*alocont_o_nn3d(i,jjm1,k,q1) + &
                  c12_intu*alocont_o_nn3d(i,jjm1,k,q2) + &
                  c14_intu*alocont_o_nn3d(i,jjm1,k,q4) + &
                  c15_intu*alocont_o_nn3d(i,jjm1,k,q5) + &
                  c17_intu*alocont_o_nn3d(i,jjm1,k,q7) + &
                  c18_intu*alocont_o_nn3d(i,jjm1,k,q8) + &
                  c23_intu*alocont_o_nn3d(i,jjm1,k,q13) + &
                  c24_intu*alocont_o_nn3d(i,jjm1,k,q14) + &
                  c20_intu*alocont_o_nn3d(i,jjm1,k,q10) + &
                  c26_intu*alocont_o_nn3d(i,jjm1,k,q16)))
               alocont_o_nn3d(iim1,jjm1,k,q18) = imask_totreg3d(iim1,jjm1,k) * &
                  (alo_u*c23_scontu + abs_sc*(c10_intu*alocont_o_nn3d(iim1,jjm1,k,q1) + &
                  c11_intu*alocont_o_nn3d(iim1,jjm1,k,q2) + &
                  c12_intu*alocont_o_nn3d(iim1,jjm1,k,q3) + &
                  c13_intu*alocont_o_nn3d(iim1,jjm1,k,q4) + &
                  c14_intu*alocont_o_nn3d(iim1,jjm1,k,q5) + &
                  c15_intu*alocont_o_nn3d(iim1,jjm1,k,q6) + &
                  c16_intu*alocont_o_nn3d(iim1,jjm1,k,q7) + &
                  c17_intu*alocont_o_nn3d(iim1,jjm1,k,q8) + &
                  c18_intu*alocont_o_nn3d(iim1,jjm1,k,q9) + &
                  c22_intu*alocont_o_nn3d(iim1,jjm1,k,q13) + &
                  c23_intu*alocont_o_nn3d(iim1,jjm1,k,q14) + &
                  c24_intu*alocont_o_nn3d(iim1,jjm1,k,q15) + &
                  c20_intu*alocont_o_nn3d(iim1,jjm1,k,q11) + &
                  c26_intu*alocont_o_nn3d(iim1,jjm1,k,q17)))
               alocont_o_nn3d(iip1,jjp1,kkm1,q19) = imask_totreg3d(iip1,jjp1,kkm1) * &
                  (alo_d*c09_scontd + abs_sc*c18_intu*alocont_o_nn3d(iip1,jjp1,kkm1,q10))
               alocont_o_nn3d(i,jjp1,kkm1,q20) = imask_totreg3d(i,jjp1,kkm1) * &
                  (alo_d*c08_scontd + abs_sc*(c17_intu*alocont_o_nn3d(i,jjp1,kkm1,q10) + &
                  c18_intu*alocont_o_nn3d(i,jjp1,kkm1,q11) + &
                  c08_intu*alocont_o_nn3d(i,jjp1,kkm1,q1) + &
                  c26_intu*alocont_o_nn3d(i,jjp1,kkm1,q19)))
               alocont_o_nn3d(iim1,jjp1,kkm1,q21) = imask_totreg3d(iim1,jjp1,kkm1) * &
                  (alo_d*c07_scontd + abs_sc*(c16_intu*alocont_o_nn3d(iim1,jjp1,kkm1,q10) + &
                  c17_intu*alocont_o_nn3d(iim1,jjp1,kkm1,q11) + &
                  c18_intu*alocont_o_nn3d(iim1,jjp1,kkm1,q12) + &
                  c08_intu*alocont_o_nn3d(iim1,jjp1,kkm1,q2) + &
                  c26_intu*alocont_o_nn3d(iim1,jjp1,kkm1,q20)))
               alocont_o_nn3d(iip1,j,kkm1,q22) = imask_totreg3d(iip1,j,kkm1) * &
                  (alo_d*c06_scontd + abs_sc*(c15_intu*alocont_o_nn3d(iip1,j,kkm1,q10) + &
                  c18_intu*alocont_o_nn3d(iip1,j,kkm1,q13) + &
                  c06_intu*alocont_o_nn3d(iip1,j,kkm1,q1) + &
                  c24_intu*alocont_o_nn3d(iip1,j,kkm1,q19)))
               alocont_o_nn3d(i,j,kkm1,q23) = imask_totreg3d(i,j,kkm1) * &
                  (alo_u*c18_scontu + abs_sc*(c14_intu*alocont_o_nn3d(i,j,kkm1,q10) + &
                  c15_intu*alocont_o_nn3d(i,j,kkm1,q11) + &
                  c17_intu*alocont_o_nn3d(i,j,kkm1,q13) + &
                  c18_intu*alocont_o_nn3d(i,j,kkm1,q14) + &
                  c05_intu*alocont_o_nn3d(i,j,kkm1,q1) + &
                  c06_intu*alocont_o_nn3d(i,j,kkm1,q2) + &
                  c23_intu*alocont_o_nn3d(i,j,kkm1,q19) + &
                  c24_intu*alocont_o_nn3d(i,j,kkm1,q20) + &
                  c08_intu*alocont_o_nn3d(i,j,kkm1,q4) + &
                  c26_intu*alocont_o_nn3d(i,j,kkm1,q22)))
               alocont_o_nn3d(iim1,j,kkm1,q24) = imask_totreg3d(iim1,j,kkm1) * &
                  (alo_u*c17_scontu + abs_sc*(c13_intu*alocont_o_nn3d(iim1,j,kkm1,q10) + &
                  c14_intu*alocont_o_nn3d(iim1,j,kkm1,q11) + &
                  c15_intu*alocont_o_nn3d(iim1,j,kkm1,q12) + &
                  c16_intu*alocont_o_nn3d(iim1,j,kkm1,q13) + &
                  c17_intu*alocont_o_nn3d(iim1,j,kkm1,q14) + &
                  c18_intu*alocont_o_nn3d(iim1,j,kkm1,q15) + &
                  c04_intu*alocont_o_nn3d(iim1,j,kkm1,q1) + &
                  c05_intu*alocont_o_nn3d(iim1,j,kkm1,q2) + &
                  c06_intu*alocont_o_nn3d(iim1,j,kkm1,q3) + &
                  c22_intu*alocont_o_nn3d(iim1,j,kkm1,q19) + &
                  c23_intu*alocont_o_nn3d(iim1,j,kkm1,q20) + &
                  c24_intu*alocont_o_nn3d(iim1,j,kkm1,q21) + &
                  c08_intu*alocont_o_nn3d(iim1,j,kkm1,q5) + &
                  c26_intu*alocont_o_nn3d(iim1,j,kkm1,q23)))
               alocont_o_nn3d(iip1,jjm1,kkm1,q25) = imask_totreg3d(iip1,jjm1,kkm1) * &
                  (alo_d*c03_scontd + abs_sc*(c12_intu*alocont_o_nn3d(iip1,jjm1,kkm1,q10) + &
                  c15_intu*alocont_o_nn3d(iip1,jjm1,kkm1,q13) + &
                  c18_intu*alocont_o_nn3d(iip1,jjm1,kkm1,q16) + &
                  c06_intu*alocont_o_nn3d(iip1,jjm1,kkm1,q4) + &
                  c24_intu*alocont_o_nn3d(iip1,jjm1,kkm1,q22)))
               alocont_o_nn3d(i,jjm1,kkm1,q26) = imask_totreg3d(i,jjm1,kkm1) * &
                  (alo_u*c15_scontu + abs_sc*(c11_intu*alocont_o_nn3d(i,jjm1,kkm1,q10) + &
                  c12_intu*alocont_o_nn3d(i,jjm1,kkm1,q11) + &
                  c14_intu*alocont_o_nn3d(i,jjm1,kkm1,q13) + &
                  c15_intu*alocont_o_nn3d(i,jjm1,kkm1,q14) + &
                  c17_intu*alocont_o_nn3d(i,jjm1,kkm1,q16) + &
                  c18_intu*alocont_o_nn3d(i,jjm1,kkm1,q17) + &
                  c05_intu*alocont_o_nn3d(i,jjm1,kkm1,q4) + &
                  c06_intu*alocont_o_nn3d(i,jjm1,kkm1,q5) + &
                  c23_intu*alocont_o_nn3d(i,jjm1,kkm1,q22) + &
                  c24_intu*alocont_o_nn3d(i,jjm1,kkm1,q23) + &
                  c02_intu*alocont_o_nn3d(i,jjm1,kkm1,q1) + &
                  c08_intu*alocont_o_nn3d(i,jjm1,kkm1,q7) + &
                  c20_intu*alocont_o_nn3d(i,jjm1,kkm1,q19) + &
                  c26_intu*alocont_o_nn3d(i,jjm1,kkm1,q25)))
               alocont_o_nn3d(iim1,jjm1,kkm1,q27) = imask_totreg3d(iim1,jjm1,kkm1) * &
                  (alo_u*c14_scontu + abs_sc*(c10_intu*alocont_o_nn3d(iim1,jjm1,kkm1,q10) + &
                  c11_intu*alocont_o_nn3d(iim1,jjm1,kkm1,q11) + &
                  c12_intu*alocont_o_nn3d(iim1,jjm1,kkm1,q12) + &
                  c13_intu*alocont_o_nn3d(iim1,jjm1,kkm1,q13) + &
                  c14_intu*alocont_o_nn3d(iim1,jjm1,kkm1,q14) + &
                  c15_intu*alocont_o_nn3d(iim1,jjm1,kkm1,q15) + &
                  c16_intu*alocont_o_nn3d(iim1,jjm1,kkm1,q16) + &
                  c17_intu*alocont_o_nn3d(iim1,jjm1,kkm1,q17) + &
                  c18_intu*alocont_o_nn3d(iim1,jjm1,kkm1,q18) + &
                  c04_intu*alocont_o_nn3d(iim1,jjm1,kkm1,q4) + &
                  c05_intu*alocont_o_nn3d(iim1,jjm1,kkm1,q5) + &
                  c06_intu*alocont_o_nn3d(iim1,jjm1,kkm1,q6) + &
                  c22_intu*alocont_o_nn3d(iim1,jjm1,kkm1,q22) + &
                  c23_intu*alocont_o_nn3d(iim1,jjm1,kkm1,q23) + &
                  c24_intu*alocont_o_nn3d(iim1,jjm1,kkm1,q24) + &
                  c02_intu*alocont_o_nn3d(iim1,jjm1,kkm1,q2) + &
                  c08_intu*alocont_o_nn3d(iim1,jjm1,kkm1,q8) + &
                  c20_intu*alocont_o_nn3d(iim1,jjm1,kkm1,q20) + &
                  c26_intu*alocont_o_nn3d(iim1,jjm1,kkm1,q26)))
!               call cpu_time(te_aloo)
!               tt_aloo=tt_aloo+te_aloo-ts_aloo
!
!perform angular integration
!*********debug start
!               dels_r=dist_sphere(nn_x,nn_y,nn_z,x(i),y(j),z(k))
!               if(dels_r.ne.1.d10) then
!                  int3d(i,j,k)=xic1
!               else
!                  int3d(i,j,k)=0.d0
!               endif
!*********debug end
!
!               call cpu_time(ts_integ)
               mint3d_tmp(i,j,k) = mint3d_tmp(i,j,k) + int3d(i,j,k)*wall
               fcontx3d_tmp(i,j,k) = fcontx3d_tmp(i,j,k) + int3d(i,j,k)*nn_x*wall
               fconty3d_tmp(i,j,k) = fconty3d_tmp(i,j,k) + int3d(i,j,k)*nn_y*wall
               fcontz3d_tmp(i,j,k) = fcontz3d_tmp(i,j,k) + int3d(i,j,k)*nn_z*wall
               kcontxx3d_tmp(i,j,k) = kcontxx3d_tmp(i,j,k) + int3d(i,j,k)*nn_x*nn_x*wall
               kcontyy3d_tmp(i,j,k) = kcontyy3d_tmp(i,j,k) + int3d(i,j,k)*nn_y*nn_y*wall
               kcontzz3d_tmp(i,j,k) = kcontzz3d_tmp(i,j,k) + int3d(i,j,k)*nn_z*nn_z*wall
               kcontxy3d_tmp(i,j,k) = kcontxy3d_tmp(i,j,k) + int3d(i,j,k)*nn_x*nn_y*wall
               kcontxz3d_tmp(i,j,k) = kcontxz3d_tmp(i,j,k) + int3d(i,j,k)*nn_x*nn_z*wall
               kcontyz3d_tmp(i,j,k) = kcontyz3d_tmp(i,j,k) + int3d(i,j,k)*nn_y*nn_z*wall
               normalization3d_tmp(i,j,k) = normalization3d_tmp(i,j,k) + wall
!
               alocont_nn3d_tmp(iip1,jjp1,kkp1,q1) = alocont_nn3d_tmp(iip1,jjp1,kkp1,q1) + wall*alocont_o_nn3d(iip1,jjp1,kkp1,q1)
               alocont_nn3d_tmp(i,jjp1,kkp1,q2) = alocont_nn3d_tmp(i,jjp1,kkp1,q2) + wall*alocont_o_nn3d(i,jjp1,kkp1,q2)
               alocont_nn3d_tmp(iim1,jjp1,kkp1,q3) = alocont_nn3d_tmp(iim1,jjp1,kkp1,q3) + wall*alocont_o_nn3d(iim1,jjp1,kkp1,q3)
               alocont_nn3d_tmp(iip1,j,kkp1,q4) = alocont_nn3d_tmp(iip1,j,kkp1,q4) + wall*alocont_o_nn3d(iip1,j,kkp1,q4)
               alocont_nn3d_tmp(i,j,kkp1,q5) = alocont_nn3d_tmp(i,j,kkp1,q5) + wall*alocont_o_nn3d(i,j,kkp1,q5)
               alocont_nn3d_tmp(iim1,j,kkp1,q6) = alocont_nn3d_tmp(iim1,j,kkp1,q6) + wall*alocont_o_nn3d(iim1,j,kkp1,q6)
               alocont_nn3d_tmp(iip1,jjm1,kkp1,q7) = alocont_nn3d_tmp(iip1,jjm1,kkp1,q7) + wall*alocont_o_nn3d(iip1,jjm1,kkp1,q7)
               alocont_nn3d_tmp(i,jjm1,kkp1,q8) = alocont_nn3d_tmp(i,jjm1,kkp1,q8) + wall*alocont_o_nn3d(i,jjm1,kkp1,q8)
               alocont_nn3d_tmp(iim1,jjm1,kkp1,q9) = alocont_nn3d_tmp(iim1,jjm1,kkp1,q9) + wall*alocont_o_nn3d(iim1,jjm1,kkp1,q9)
               alocont_nn3d_tmp(iip1,jjp1,k,q10) = alocont_nn3d_tmp(iip1,jjp1,k,q10) + wall*alocont_o_nn3d(iip1,jjp1,k,q10)
               alocont_nn3d_tmp(i,jjp1,k,q11) = alocont_nn3d_tmp(i,jjp1,k,q11) + wall*alocont_o_nn3d(i,jjp1,k,q11)
               alocont_nn3d_tmp(iim1,jjp1,k,q12) = alocont_nn3d_tmp(iim1,jjp1,k,q12) + wall*alocont_o_nn3d(iim1,jjp1,k,q12)
               alocont_nn3d_tmp(iip1,j,k,q13) = alocont_nn3d_tmp(iip1,j,k,q13) + wall*alocont_o_nn3d(iip1,j,k,q13)
               alocont_nn3d_tmp(i,j,k,q14) = alocont_nn3d_tmp(i,j,k,q14) + wall*alocont_o_nn3d(i,j,k,q14)
               alocont_nn3d_tmp(iim1,j,k,q15) = alocont_nn3d_tmp(iim1,j,k,q15) + wall*alocont_o_nn3d(iim1,j,k,q15)
               alocont_nn3d_tmp(iip1,jjm1,k,q16) = alocont_nn3d_tmp(iip1,jjm1,k,q16) + wall*alocont_o_nn3d(iip1,jjm1,k,q16)
               alocont_nn3d_tmp(i,jjm1,k,q17) = alocont_nn3d_tmp(i,jjm1,k,q17) + wall*alocont_o_nn3d(i,jjm1,k,q17)
               alocont_nn3d_tmp(iim1,jjm1,k,q18) = alocont_nn3d_tmp(iim1,jjm1,k,q18) + wall*alocont_o_nn3d(iim1,jjm1,k,q18)
               alocont_nn3d_tmp(iip1,jjp1,kkm1,q19) = alocont_nn3d_tmp(iip1,jjp1,kkm1,q19) + wall*alocont_o_nn3d(iip1,jjp1,kkm1,q19)
               alocont_nn3d_tmp(i,jjp1,kkm1,q20) = alocont_nn3d_tmp(i,jjp1,kkm1,q20) + wall*alocont_o_nn3d(i,jjp1,kkm1,q20)
               alocont_nn3d_tmp(iim1,jjp1,kkm1,q21) = alocont_nn3d_tmp(iim1,jjp1,kkm1,q21) + wall*alocont_o_nn3d(iim1,jjp1,kkm1,q21)
               alocont_nn3d_tmp(iip1,j,kkm1,q22) = alocont_nn3d_tmp(iip1,j,kkm1,q22) + wall*alocont_o_nn3d(iip1,j,kkm1,q22)
               alocont_nn3d_tmp(i,j,kkm1,q23) = alocont_nn3d_tmp(i,j,kkm1,q23) + wall*alocont_o_nn3d(i,j,kkm1,q23)
               alocont_nn3d_tmp(iim1,j,kkm1,q24) = alocont_nn3d_tmp(iim1,j,kkm1,q24) + wall*alocont_o_nn3d(iim1,j,kkm1,q24)
               alocont_nn3d_tmp(iip1,jjm1,kkm1,q25) = alocont_nn3d_tmp(iip1,jjm1,kkm1,q25) + wall*alocont_o_nn3d(iip1,jjm1,kkm1,q25)
               alocont_nn3d_tmp(i,jjm1,kkm1,q26) = alocont_nn3d_tmp(i,jjm1,kkm1,q26) + wall*alocont_o_nn3d(i,jjm1,kkm1,q26)
               alocont_nn3d_tmp(iim1,jjm1,kkm1,q27) = alocont_nn3d_tmp(iim1,jjm1,kkm1,q27) + wall*alocont_o_nn3d(iim1,jjm1,kkm1,q27)
!               call cpu_time(te_integ)
!               tt_integ=tt_integ+te_integ-ts_integ
!
!               call cpu_time(te_case2)
!               tt_case2=tt_case2+te_case2-ts_case2
!
!***********special treatment: point lies directly on boundary**********
!
             case(3)
!
!               call cpu_time(ts_case3)
!
               mueff=nn_x*x(i)+nn_y*y(j)+nn_z*z(k)
               if(mueff.ge.0.d0) then
!no interpolations required
!for alo
!                 int3d(i,j,k) = 0.d0
!                  int3d(i,j,k) = xic1
                  int3d(i,j,k) = calc_icore_gdark(z(k), sqrt(x(i)**2+y(j)**2+z(k)**2))
                  alocont_o_nn3d(i,j,k,14) = 0.d0
!
!perform angular integration (alo not required, since boundary specified)
                  mint3d_tmp(i,j,k) = mint3d_tmp(i,j,k) + int3d(i,j,k)*wall
                  fcontx3d_tmp(i,j,k) = fcontx3d_tmp(i,j,k) + int3d(i,j,k)*nn_x*wall
                  fconty3d_tmp(i,j,k) = fconty3d_tmp(i,j,k) + int3d(i,j,k)*nn_y*wall
                  fcontz3d_tmp(i,j,k) = fcontz3d_tmp(i,j,k) + int3d(i,j,k)*nn_z*wall
                  kcontxx3d_tmp(i,j,k) = kcontxx3d_tmp(i,j,k) + int3d(i,j,k)*nn_x*nn_x*wall
                  kcontyy3d_tmp(i,j,k) = kcontyy3d_tmp(i,j,k) + int3d(i,j,k)*nn_y*nn_y*wall
                  kcontzz3d_tmp(i,j,k) = kcontzz3d_tmp(i,j,k) + int3d(i,j,k)*nn_z*nn_z*wall
                  kcontxy3d_tmp(i,j,k) = kcontxy3d_tmp(i,j,k) + int3d(i,j,k)*nn_x*nn_y*wall
                  kcontxz3d_tmp(i,j,k) = kcontxz3d_tmp(i,j,k) + int3d(i,j,k)*nn_x*nn_z*wall
                  kcontyz3d_tmp(i,j,k) = kcontyz3d_tmp(i,j,k) + int3d(i,j,k)*nn_y*nn_z*wall
                  normalization3d_tmp(i,j,k) = normalization3d_tmp(i,j,k) + wall
               else
!same interpolation as in (1)
!
                  iim1=i-alpha
                  jjm1=j-beta
                  kkm1=k-gamma
                  iip1=i+alpha
                  jjp1=j+beta
                  kkp1=k+gamma
!
!calculate distance of upwind point to (previous) xy-plane, xz-plane and yz-plane
                  dels_xyu=(z(k)-z(kkm1))/nn_z
                  dels_xzu=(y(j)-y(jjm1))/nn_y
                  dels_yzu=(x(i)-x(iim1))/nn_x
                  dels_u=min(dels_xyu,dels_xzu,dels_yzu)
!
!calculate distance of downwind point to (next) xy-plane, xz-plane and yz-plane
                  dels_xyd=(z(kkp1)-z(k))/nn_z
                  dels_xzd=(y(jjp1)-y(j))/nn_y
                  dels_yzd=(x(iip1)-x(i))/nn_x
                  dels_d=min(dels_xyd,dels_xzd,dels_yzd)
!
!----------------------------local point--------------------------------
!
                  scont_p=scont3d(i,j,k)
                  opac_p=opac3d(i,j,k)
!
!----------------------------upwind point-------------------------------
!
                  if(dels_xyu.eq.dels_u) then
!intersection with x-y plane on level kkm1
                     x_u = x(i) - dels_u*nn_x
                     y_u = y(j) - dels_u*nn_y
                     z_u = z(kkm1)
!
                     iim2=i-2*alpha
                     jjm2=j-2*beta
!
!                     call cpu_time(ts_interpu)
                     call coeff2d_contu(opac3d(iim2,jjm2,kkm1), opac3d(iim1,jjm2,kkm1), opac3d(i,jjm2,kkm1), &
                        opac3d(iim2,jjm1,kkm1), opac3d(iim1,jjm1,kkm1), opac3d(i,jjm1,kkm1), &
                        opac3d(iim2,j,kkm1), opac3d(iim1,j,kkm1), opac3d(i,j,kkm1), &
                        scont3d(iim2,jjm2,kkm1), scont3d(iim1,jjm2,kkm1), scont3d(i,jjm2,kkm1), &
                        scont3d(iim2,jjm1,kkm1), scont3d(iim1,jjm1,kkm1), scont3d(i,jjm1,kkm1), &
                        scont3d(iim2,j,kkm1), scont3d(iim1,j,kkm1), scont3d(i,j,kkm1), &
                        int3d(iim2,jjm2,kkm1), int3d(iim1,jjm2,kkm1), int3d(i,jjm2,kkm1), &
                        int3d(iim2,jjm1,kkm1), int3d(iim1,jjm1,kkm1), int3d(i,jjm1,kkm1), &
                        int3d(iim2,j,kkm1), int3d(iim1,j,kkm1), int3d(i,j,kkm1), &
                        x(iim2), x(iim1), x(i), y(jjm2), y(jjm1), y(j), x_u, y_u, &
                        c10_scontu, c11_scontu, c12_scontu, c13_scontu, c14_scontu, &
                        c15_scontu, c16_scontu, c17_scontu, c18_scontu, &
                        c10_intu, c11_intu, c12_intu, c13_intu, c14_intu, &
                        c15_intu, c16_intu, c17_intu, c18_intu, opac_u, scont_u, int_u)
!                     call cpu_time(te_interpu)
!                     tt_interpu=tt_interpu+te_interpu-ts_interpu
!
!set interpolation coefficients that are not used to zero
                     c23_scontu = 0.d0
                     c24_scontu = 0.d0
                     c26_scontu = 0.d0
                     c27_scontu = 0.d0
!
                     c02_intu = 0.d0
                     c04_intu = 0.d0
                     c05_intu = 0.d0
                     c06_intu = 0.d0
                     c08_intu = 0.d0
                     c20_intu = 0.d0
                     c22_intu = 0.d0
                     c23_intu = 0.d0
                     c24_intu = 0.d0
                     c26_intu = 0.d0
!
                  elseif(dels_xzu.eq.dels_u) then
!intersection with x-z plane on level jjm1
                     x_u = x(i) - dels_u*nn_x
                     y_u = y(jjm1)
                     z_u = z(k) - dels_u*nn_z
!
                     iim2=i-2*alpha
                     kkm2=k-2*gamma
!
!                     call cpu_time(ts_interpu)
                     call coeff2d_contu(opac3d(iim2,jjm1,kkm2), opac3d(iim1,jjm1,kkm2), opac3d(i,jjm1,kkm2), &
                        opac3d(iim2,jjm1,kkm1), opac3d(iim1,jjm1,kkm1), opac3d(i,jjm1,kkm1), &
                        opac3d(iim2,jjm1,k), opac3d(iim1,jjm1,k), opac3d(i,jjm1,k), &
                        scont3d(iim2,jjm1,kkm2), scont3d(iim1,jjm1,kkm2), scont3d(i,jjm1,kkm2), &
                        scont3d(iim2,jjm1,kkm1), scont3d(iim1,jjm1,kkm1), scont3d(i,jjm1,kkm1), &
                        scont3d(iim2,jjm1,k), scont3d(iim1,jjm1,k), scont3d(i,jjm1,k), &
                        int3d(iim2,jjm1,kkm2), int3d(iim1,jjm1,kkm2), int3d(i,jjm1,kkm2), &
                        int3d(iim2,jjm1,kkm1), int3d(iim1,jjm1,kkm1), int3d(i,jjm1,kkm1), &
                        int3d(iim2,jjm1,k), int3d(iim1,jjm1,k), int3d(i,jjm1,k), &
                        x(iim2), x(iim1), x(i), z(kkm2), z(kkm1), z(k), x_u, z_u, &
                        c04_scontu, c05_scontu, c06_scontu, c13_scontu, c14_scontu, &
                        c15_scontu, c22_scontu, c23_scontu, c24_scontu, &
                        c04_intu, c05_intu, c06_intu, c13_intu, c14_intu, &
                        c15_intu, c22_intu, c23_intu, c24_intu, opac_u, scont_u, int_u)
!                     call cpu_time(te_interpu)
!                     tt_interpu=tt_interpu+te_interpu-ts_interpu
!
!set interpolation coefficients that are not used to zero
                     c17_scontu = 0.d0
                     c18_scontu = 0.d0
                     c26_scontu = 0.d0
                     c27_scontu = 0.d0
!
                     c02_intu = 0.d0
                     c08_intu = 0.d0
                     c10_intu = 0.d0
                     c11_intu = 0.d0
                     c12_intu = 0.d0
                     c16_intu = 0.d0
                     c17_intu = 0.d0
                     c18_intu = 0.d0
                     c20_intu = 0.d0
                     c26_intu = 0.d0
!
                  elseif(dels_yzu.eq.dels_u) then
!intersection with y-z plane on level iim1
                     x_u = x(iim1)
                     y_u = y(j) - dels_u*nn_y
                     z_u = z(k) - dels_u*nn_z
!
                     jjm2=j-2*beta
                     kkm2=k-2*gamma
!
!                     call cpu_time(ts_interpu)
                     call coeff2d_contu(opac3d(iim1,jjm2,kkm2), opac3d(iim1,jjm1,kkm2), opac3d(iim1,j,kkm2), &
                        opac3d(iim1,jjm2,kkm1), opac3d(iim1,jjm1,kkm1), opac3d(iim1,j,kkm1), &
                        opac3d(iim1,jjm2,k), opac3d(iim1,jjm1,k), opac3d(iim1,j,k), &
                        scont3d(iim1,jjm2,kkm2), scont3d(iim1,jjm1,kkm2), scont3d(iim1,j,kkm2), &
                        scont3d(iim1,jjm2,kkm1), scont3d(iim1,jjm1,kkm1), scont3d(iim1,j,kkm1), &
                        scont3d(iim1,jjm2,k), scont3d(iim1,jjm1,k), scont3d(iim1,j,k), &
                        int3d(iim1,jjm2,kkm2), int3d(iim1,jjm1,kkm2), int3d(iim1,j,kkm2), &
                        int3d(iim1,jjm2,kkm1), int3d(iim1,jjm1,kkm1), int3d(iim1,j,kkm1), &
                        int3d(iim1,jjm2,k), int3d(iim1,jjm1,k), int3d(iim1,j,k), &
                        y(jjm2), y(jjm1), y(j), z(kkm2), z(kkm1), z(k), y_u, z_u, &
                        c02_scontu, c05_scontu, c08_scontu, c11_scontu, c14_scontu, &
                        c17_scontu, c20_scontu, c23_scontu, c26_scontu, &
                        c02_intu, c05_intu, c08_intu, c11_intu, c14_intu, &
                        c17_intu, c20_intu, c23_intu, c26_intu, opac_u, scont_u, int_u)
!                     call cpu_time(te_interpu)
!                     tt_interpu=tt_interpu+te_interpu-ts_interpu
!
!set interpolation coefficients that are not used to zero
                     c15_scontu = 0.d0
                     c18_scontu = 0.d0
                     c24_scontu = 0.d0
                     c27_scontu = 0.d0
!
                     c04_intu = 0.d0
                     c06_intu = 0.d0
                     c10_intu = 0.d0
                     c12_intu = 0.d0
                     c13_intu = 0.d0
                     c15_intu = 0.d0
                     c16_intu = 0.d0
                     c18_intu = 0.d0
                     c22_intu = 0.d0
                     c24_intu = 0.d0
!
                  else
                     write(*,'(4es20.8,2l4)') dels_u, dels_xzu, dels_xyu, dels_yzu, dels_u.eq.dels_xzu, dels_xyu.eq.dels_u
                     stop 'error in fsc_cont3d: invalid dels_u'
                  endif
!
!---------------------------downwind point------------------------------
!
                  if(dels_xyd.eq.dels_d) then
!intersection with x-y plane on level kkp1
                     x_d = x(i) + dels_d*nn_x
                     y_d = y(j) + dels_d*nn_y
                     z_d = z(kkp1)
!
!                     call cpu_time(ts_interpd)
                     call coeff2d_contd(opac3d(iim1,jjm1,kkp1), opac3d(i,jjm1,kkp1), opac3d(iip1,jjm1,kkp1), &
                        opac3d(iim1,j,kkp1), opac3d(i,j,kkp1), opac3d(iip1,j,kkp1), &
                        opac3d(iim1,jjp1,kkp1), opac3d(i,jjp1,kkp1), opac3d(iip1,jjp1,kkp1), &
                        scont3d(iim1,jjm1,kkp1), scont3d(i,jjm1,kkp1), scont3d(iip1,jjm1,kkp1), &
                        scont3d(iim1,j,kkp1), scont3d(i,j,kkp1), scont3d(iip1,j,kkp1), &
                        scont3d(iim1,jjp1,kkp1), scont3d(i,jjp1,kkp1), scont3d(iip1,jjp1,kkp1), &
                        x(iim1), x(i), x(iip1), y(jjm1), y(j), y(jjp1), x_d, y_d, &
                        c19_scontd, c20_scontd, c21_scontd, c22_scontd, c23_scontd, c24_scontd, &
                        c25_scontd, c26_scontd, c27_scontd, opac_d, scont_d)
!                     call cpu_time(te_interpd)
!                     tt_interpd=tt_interpd+te_interpd-ts_interpd
!
!set interpolation coefficients that are not used to zero
                     c03_scontd = 0.d0
                     c06_scontd = 0.d0
                     c07_scontd = 0.d0
                     c08_scontd = 0.d0
                     c09_scontd = 0.d0
                     c12_scontd = 0.d0
                     c15_scontd = 0.d0
                     c16_scontd = 0.d0
                     c17_scontd = 0.d0
                     c18_scontd = 0.d0
!
                  elseif(dels_xzd.eq.dels_d) then
!intersection with x-z plane on level jjp1
                     x_d = x(i) + dels_d*nn_x
                     y_d = y(jjp1)
                     z_d = z(k) + dels_d*nn_z
!
!                     call cpu_time(ts_interpd)
                     call coeff2d_contd(opac3d(iim1,jjp1,kkm1), opac3d(i,jjp1,kkm1), opac3d(iip1,jjp1,kkm1), &
                        opac3d(iim1,jjp1,k), opac3d(i,jjp1,k), opac3d(iip1,jjp1,k), &
                        opac3d(iim1,jjp1,kkp1), opac3d(i,jjp1,kkp1), opac3d(iip1,jjp1,kkp1), &
                        scont3d(iim1,jjp1,kkm1), scont3d(i,jjp1,kkm1), scont3d(iip1,jjp1,kkm1), &
                        scont3d(iim1,jjp1,k), scont3d(i,jjp1,k), scont3d(iip1,jjp1,k), &
                        scont3d(iim1,jjp1,kkp1), scont3d(i,jjp1,kkp1), scont3d(iip1,jjp1,kkp1), &
                        x(iim1), x(i), x(iip1), z(kkm1), z(k), z(kkp1), x_d, z_d, &
                        c07_scontd, c08_scontd, c09_scontd, c16_scontd, c17_scontd, c18_scontd, &
                        c25_scontd, c26_scontd, c27_scontd, opac_d, scont_d)
!                     call cpu_time(te_interpd)
!                     tt_interpd=tt_interpd+te_interpd-ts_interpd
!
!set interpolation coefficients that are not used to zero
                     c03_scontd = 0.d0
                     c06_scontd = 0.d0
                     c12_scontd = 0.d0
                     c15_scontd = 0.d0
                     c19_scontd = 0.d0
                     c20_scontd = 0.d0
                     c21_scontd = 0.d0
                     c22_scontd = 0.d0
                     c23_scontd = 0.d0
                     c24_scontd = 0.d0
!
                  elseif(dels_yzd.eq.dels_d) then
!intersection with y-z plane on level iip1
                     x_d = x(iip1)
                     y_d = y(j) + dels_d*nn_y
                     z_d = z(k) + dels_d*nn_z
!
!                     call cpu_time(ts_interpd)
                     call coeff2d_contd(opac3d(iip1,jjm1,kkm1), opac3d(iip1,j,kkm1), opac3d(iip1,jjp1,kkm1), &
                        opac3d(iip1,jjm1,k), opac3d(iip1,j,k), opac3d(iip1,jjp1,k), &
                        opac3d(iip1,jjm1,kkp1), opac3d(iip1,j,kkp1), opac3d(iip1,jjp1,kkp1), &
                        scont3d(iip1,jjm1,kkm1), scont3d(iip1,j,kkm1), scont3d(iip1,jjp1,kkm1), &
                        scont3d(iip1,jjm1,k), scont3d(iip1,j,k), scont3d(iip1,jjp1,k), &
                        scont3d(iip1,jjm1,kkp1), scont3d(iip1,j,kkp1), scont3d(iip1,jjp1,kkp1), &
                        y(jjm1), y(j), y(jjp1), z(kkm1), z(k), z(kkp1), y_d, z_d, &
                        c03_scontd, c06_scontd, c09_scontd, c12_scontd, c15_scontd, c18_scontd, &
                        c21_scontd, c24_scontd, c27_scontd, opac_d, scont_d)
!                     call cpu_time(te_interpd)
!                     tt_interpd=tt_interpd+te_interpd-ts_interpd
!
!set interpolation coefficients that are not used to zero
                     c07_scontd = 0.d0
                     c08_scontd = 0.d0
                     c16_scontd = 0.d0
                     c17_scontd = 0.d0
                     c19_scontd = 0.d0
                     c20_scontd = 0.d0
                     c22_scontd = 0.d0
                     c23_scontd = 0.d0
                     c25_scontd = 0.d0
                     c26_scontd = 0.d0
!
                  else
                     write(*,'(4es20.8,l4)') dels_d, dels_xzd, dels_xyd, dels_yzd, dels_d.eq.dels_xzd
                     stop 'error in fsc_cont3d: invalid dels_d'
                  endif
!
!--------------------------------radiative transfer---------------------
!
!                  call cpu_time(ts_fs1d)
                  call fsc_cont(int_u, opac_u, opac_p, opac_d, scont_u, scont_p, scont_d, &
                     dels_u, dels_d, abs_sc, contr_sc, int_sc, alo_u, alo_p, alo_d)
                  int3d(i,j,k) = abs_sc*int_u + alo_u*scont_u + alo_p*scont_p + alo_d*scont_d
!                  call cpu_time(te_fs1d)
!                  tt_fs1d=tt_fs1d+te_fs1d-ts_fs1d
!
!                  call cpu_time(ts_aloo)
                  alocont_o_nn3d(iip1,jjp1,kkp1,q1) = imask_totreg3d(iip1,jjp1,kkp1)*alo_d*c27_scontd
                  alocont_o_nn3d(i,jjp1,kkp1,q2) = imask_totreg3d(i,jjp1,kkp1) * &
                     (alo_d*c26_scontd + abs_sc*c26_intu*alocont_o_nn3d(i,jjp1,kkp1,q1))
                  alocont_o_nn3d(iim1,jjp1,kkp1,q3) = imask_totreg3d(iim1,jjp1,kkp1) * &
                     (alo_d*c25_scontd + abs_sc*c26_intu*alocont_o_nn3d(iim1,jjp1,kkp1,q2))
                  alocont_o_nn3d(iip1,j,kkp1,q4) = imask_totreg3d(iip1,j,kkp1) * &
                     (alo_d*c24_scontd + abs_sc*c24_intu*alocont_o_nn3d(iip1,j,kkp1,q1))
                  alocont_o_nn3d(i,j,kkp1,q5) = imask_totreg3d(i,j,kkp1) * &
                     (alo_d*c23_scontd + abs_sc*(c23_intu*alocont_o_nn3d(i,j,kkp1,q1) + &
                     c24_intu*alocont_o_nn3d(i,j,kkp1,q2) + &
                     c26_intu*alocont_o_nn3d(i,j,kkp1,q4)))
                  alocont_o_nn3d(iim1,j,kkp1,q6) = imask_totreg3d(iim1,j,kkp1) * &
                     (alo_d*c22_scontd + abs_sc*(c22_intu*alocont_o_nn3d(iim1,j,kkp1,q1) + &
                     c23_intu*alocont_o_nn3d(iim1,j,kkp1,q2) + &
                     c24_intu*alocont_o_nn3d(iim1,j,kkp1,q3) + &
                     c26_intu*alocont_o_nn3d(iim1,j,kkp1,q5)))
                  alocont_o_nn3d(iip1,jjm1,kkp1,q7) = imask_totreg3d(iip1,jjm1,kkp1) * &
                     (alo_d*c21_scontd + abs_sc*c24_intu*alocont_o_nn3d(iip1,jjm1,kkp1,q4))
                  alocont_o_nn3d(i,jjm1,kkp1,q8) = imask_totreg3d(i,jjm1,kkp1) * &
                     (alo_d*c20_scontd + abs_sc*(c23_intu*alocont_o_nn3d(i,jjm1,kkp1,q4) + &
                     c24_intu*alocont_o_nn3d(i,jjm1,kkp1,q5) + &
                     c20_intu*alocont_o_nn3d(i,jjm1,kkp1,q1) + &
                     c26_intu*alocont_o_nn3d(i,jjm1,kkp1,q7)))
                  alocont_o_nn3d(iim1,jjm1,kkp1,q9) = imask_totreg3d(iim1,jjm1,kkp1) * &
                     (alo_d*c19_scontd + abs_sc*(c22_intu*alocont_o_nn3d(iim1,jjm1,kkp1,q4) + &
                     c23_intu*alocont_o_nn3d(iim1,jjm1,kkp1,q5) + &
                     c24_intu*alocont_o_nn3d(iim1,jjm1,kkp1,q6) + &
                     c20_intu*alocont_o_nn3d(iim1,jjm1,kkp1,q2) + &
                     c26_intu*alocont_o_nn3d(iim1,jjm1,kkp1,q8)))
                  alocont_o_nn3d(iip1,jjp1,k,q10) = imask_totreg3d(iip1,jjp1,k) * &
                     (alo_d*c18_scontd + abs_sc*c18_intu*alocont_o_nn3d(iip1,jjp1,k,q1))
                  alocont_o_nn3d(i,jjp1,k,q11) = imask_totreg3d(i,jjp1,k) * &
                     (alo_d*c17_scontd + abs_sc*(c17_intu*alocont_o_nn3d(i,jjp1,k,q1) + &
                     c18_intu*alocont_o_nn3d(i,jjp1,k,q2) + &
                     c26_intu*alocont_o_nn3d(i,jjp1,k,q10)))
                  alocont_o_nn3d(iim1,jjp1,k,q12) = imask_totreg3d(iim1,jjp1,k) * &
                     (alo_d*c16_scontd + abs_sc*(c16_intu*alocont_o_nn3d(iim1,jjp1,k,q1) + &
                     c17_intu*alocont_o_nn3d(iim1,jjp1,k,q2) + &
                     c18_intu*alocont_o_nn3d(iim1,jjp1,k,q3) + &
                     c26_intu*alocont_o_nn3d(iim1,jjp1,k,q11)))
                  alocont_o_nn3d(iip1,j,k,q13) = imask_totreg3d(iip1,j,k) * &
                     (alo_d*c15_scontd + abs_sc*(c15_intu*alocont_o_nn3d(iip1,j,k,q1) + &
                     c18_intu*alocont_o_nn3d(iip1,j,k,q4) + &
                     c24_intu*alocont_o_nn3d(iip1,j,k,q10)))
                  alocont_o_nn3d(i,j,k,q14) = imask_totreg3d(i,j,k) * &
                     (alo_p + alo_u*c27_scontu + abs_sc*(c14_intu*alocont_o_nn3d(i,j,k,q1) + &
                     c15_intu*alocont_o_nn3d(i,j,k,q2) + &
                     c17_intu*alocont_o_nn3d(i,j,k,q4) + &
                     c26_intu*alocont_o_nn3d(i,j,k,q13) + &
                     c18_intu*alocont_o_nn3d(i,j,k,q5) + &
                     c23_intu*alocont_o_nn3d(i,j,k,q10) + &
                     c24_intu*alocont_o_nn3d(i,j,k,q11)))
                  alocont_o_nn3d(iim1,j,k,q15) = imask_totreg3d(iim1,j,k) * &
                     (alo_u*c26_scontu + abs_sc*(c13_intu*alocont_o_nn3d(iim1,j,k,q1) + &
                     c14_intu*alocont_o_nn3d(iim1,j,k,q2) + &
                     c15_intu*alocont_o_nn3d(iim1,j,k,q3) + &
                     c16_intu*alocont_o_nn3d(iim1,j,k,q4) + &
                     c17_intu*alocont_o_nn3d(iim1,j,k,q5) + &
                     c18_intu*alocont_o_nn3d(iim1,j,k,q6) + &
                     c22_intu*alocont_o_nn3d(iim1,j,k,q10) + &
                     c23_intu*alocont_o_nn3d(iim1,j,k,q11) + &
                     c24_intu*alocont_o_nn3d(iim1,j,k,q12) + &
                     c26_intu*alocont_o_nn3d(iim1,j,k,q14)))
                  alocont_o_nn3d(iip1,jjm1,k,q16) = imask_totreg3d(iip1,jjm1,k) * &
                     (alo_d*c12_scontd + abs_sc*(c12_intu*alocont_o_nn3d(iip1,jjm1,k,q1) + &
                     c15_intu*alocont_o_nn3d(iip1,jjm1,k,q4) + &
                     c18_intu*alocont_o_nn3d(iip1,jjm1,k,q7) + &
                     c24_intu*alocont_o_nn3d(iip1,jjm1,k,q13)))
                  alocont_o_nn3d(i,jjm1,k,q17) = imask_totreg3d(i,jjm1,k) * &
                     (alo_u*c24_scontu + abs_sc*(c11_intu*alocont_o_nn3d(i,jjm1,k,q1) + &
                     c12_intu*alocont_o_nn3d(i,jjm1,k,q2) + &
                     c14_intu*alocont_o_nn3d(i,jjm1,k,q4) + &
                     c15_intu*alocont_o_nn3d(i,jjm1,k,q5) + &
                     c17_intu*alocont_o_nn3d(i,jjm1,k,q7) + &
                     c18_intu*alocont_o_nn3d(i,jjm1,k,q8) + &
                     c23_intu*alocont_o_nn3d(i,jjm1,k,q13) + &
                     c24_intu*alocont_o_nn3d(i,jjm1,k,q14) + &
                     c20_intu*alocont_o_nn3d(i,jjm1,k,q10) + &
                     c26_intu*alocont_o_nn3d(i,jjm1,k,q16)))
                  alocont_o_nn3d(iim1,jjm1,k,q18) = imask_totreg3d(iim1,jjm1,k) * &
                     (alo_u*c23_scontu + abs_sc*(c10_intu*alocont_o_nn3d(iim1,jjm1,k,q1) + &
                     c11_intu*alocont_o_nn3d(iim1,jjm1,k,q2) + &
                     c12_intu*alocont_o_nn3d(iim1,jjm1,k,q3) + &
                     c13_intu*alocont_o_nn3d(iim1,jjm1,k,q4) + &
                     c14_intu*alocont_o_nn3d(iim1,jjm1,k,q5) + &
                     c15_intu*alocont_o_nn3d(iim1,jjm1,k,q6) + &
                     c16_intu*alocont_o_nn3d(iim1,jjm1,k,q7) + &
                     c17_intu*alocont_o_nn3d(iim1,jjm1,k,q8) + &
                     c18_intu*alocont_o_nn3d(iim1,jjm1,k,q9) + &
                     c22_intu*alocont_o_nn3d(iim1,jjm1,k,q13) + &
                     c23_intu*alocont_o_nn3d(iim1,jjm1,k,q14) + &
                     c24_intu*alocont_o_nn3d(iim1,jjm1,k,q15) + &
                     c20_intu*alocont_o_nn3d(iim1,jjm1,k,q11) + &
                     c26_intu*alocont_o_nn3d(iim1,jjm1,k,q17)))
                  alocont_o_nn3d(iip1,jjp1,kkm1,q19) = imask_totreg3d(iip1,jjp1,kkm1) * &
                     (alo_d*c09_scontd + abs_sc*c18_intu*alocont_o_nn3d(iip1,jjp1,kkm1,q10))
                  alocont_o_nn3d(i,jjp1,kkm1,q20) = imask_totreg3d(i,jjp1,kkm1) * &
                     (alo_d*c08_scontd + abs_sc*(c17_intu*alocont_o_nn3d(i,jjp1,kkm1,q10) + &
                     c18_intu*alocont_o_nn3d(i,jjp1,kkm1,q11) + &
                     c08_intu*alocont_o_nn3d(i,jjp1,kkm1,q1) + &
                     c26_intu*alocont_o_nn3d(i,jjp1,kkm1,q19)))
                  alocont_o_nn3d(iim1,jjp1,kkm1,q21) = imask_totreg3d(iim1,jjp1,kkm1) * &
                     (alo_d*c07_scontd + abs_sc*(c16_intu*alocont_o_nn3d(iim1,jjp1,kkm1,q10) + &
                     c17_intu*alocont_o_nn3d(iim1,jjp1,kkm1,q11) + &
                     c18_intu*alocont_o_nn3d(iim1,jjp1,kkm1,q12) + &
                     c08_intu*alocont_o_nn3d(iim1,jjp1,kkm1,q2) + &
                     c26_intu*alocont_o_nn3d(iim1,jjp1,kkm1,q20)))
                  alocont_o_nn3d(iip1,j,kkm1,q22) = imask_totreg3d(iip1,j,kkm1) * &
                     (alo_d*c06_scontd + abs_sc*(c15_intu*alocont_o_nn3d(iip1,j,kkm1,q10) + &
                     c18_intu*alocont_o_nn3d(iip1,j,kkm1,q13) + &
                     c06_intu*alocont_o_nn3d(iip1,j,kkm1,q1) + &
                     c24_intu*alocont_o_nn3d(iip1,j,kkm1,q19)))
                  alocont_o_nn3d(i,j,kkm1,q23) = imask_totreg3d(i,j,kkm1) * &
                     (alo_u*c18_scontu + abs_sc*(c14_intu*alocont_o_nn3d(i,j,kkm1,q10) + &
                     c15_intu*alocont_o_nn3d(i,j,kkm1,q11) + &
                     c17_intu*alocont_o_nn3d(i,j,kkm1,q13) + &
                     c18_intu*alocont_o_nn3d(i,j,kkm1,q14) + &
                     c05_intu*alocont_o_nn3d(i,j,kkm1,q1) + &
                     c06_intu*alocont_o_nn3d(i,j,kkm1,q2) + &
                     c23_intu*alocont_o_nn3d(i,j,kkm1,q19) + &
                     c24_intu*alocont_o_nn3d(i,j,kkm1,q20) + &
                     c08_intu*alocont_o_nn3d(i,j,kkm1,q4) + &
                     c26_intu*alocont_o_nn3d(i,j,kkm1,q22)))
                  alocont_o_nn3d(iim1,j,kkm1,q24) = imask_totreg3d(iim1,j,kkm1) * &
                     (alo_u*c17_scontu + abs_sc*(c13_intu*alocont_o_nn3d(iim1,j,kkm1,q10) + &
                     c14_intu*alocont_o_nn3d(iim1,j,kkm1,q11) + &
                     c15_intu*alocont_o_nn3d(iim1,j,kkm1,q12) + &
                     c16_intu*alocont_o_nn3d(iim1,j,kkm1,q13) + &
                     c17_intu*alocont_o_nn3d(iim1,j,kkm1,q14) + &
                     c18_intu*alocont_o_nn3d(iim1,j,kkm1,q15) + &
                     c04_intu*alocont_o_nn3d(iim1,j,kkm1,q1) + &
                     c05_intu*alocont_o_nn3d(iim1,j,kkm1,q2) + &
                     c06_intu*alocont_o_nn3d(iim1,j,kkm1,q3) + &
                     c22_intu*alocont_o_nn3d(iim1,j,kkm1,q19) + &
                     c23_intu*alocont_o_nn3d(iim1,j,kkm1,q20) + &
                     c24_intu*alocont_o_nn3d(iim1,j,kkm1,q21) + &
                     c08_intu*alocont_o_nn3d(iim1,j,kkm1,q5) + &
                     c26_intu*alocont_o_nn3d(iim1,j,kkm1,q23)))
                  alocont_o_nn3d(iip1,jjm1,kkm1,q25) = imask_totreg3d(iip1,jjm1,kkm1) * &
                     (alo_d*c03_scontd + abs_sc*(c12_intu*alocont_o_nn3d(iip1,jjm1,kkm1,q10) + &
                     c15_intu*alocont_o_nn3d(iip1,jjm1,kkm1,q13) + &
                     c18_intu*alocont_o_nn3d(iip1,jjm1,kkm1,q16) + &
                     c06_intu*alocont_o_nn3d(iip1,jjm1,kkm1,q4) + &
                     c24_intu*alocont_o_nn3d(iip1,jjm1,kkm1,q22)))
                  alocont_o_nn3d(i,jjm1,kkm1,q26) = imask_totreg3d(i,jjm1,kkm1) * &
                     (alo_u*c15_scontu + abs_sc*(c11_intu*alocont_o_nn3d(i,jjm1,kkm1,q10) + &
                     c12_intu*alocont_o_nn3d(i,jjm1,kkm1,q11) + &
                     c14_intu*alocont_o_nn3d(i,jjm1,kkm1,q13) + &
                     c15_intu*alocont_o_nn3d(i,jjm1,kkm1,q14) + &
                     c17_intu*alocont_o_nn3d(i,jjm1,kkm1,q16) + &
                     c18_intu*alocont_o_nn3d(i,jjm1,kkm1,q17) + &
                     c05_intu*alocont_o_nn3d(i,jjm1,kkm1,q4) + &
                     c06_intu*alocont_o_nn3d(i,jjm1,kkm1,q5) + &
                     c23_intu*alocont_o_nn3d(i,jjm1,kkm1,q22) + &
                     c24_intu*alocont_o_nn3d(i,jjm1,kkm1,q23) + &
                     c02_intu*alocont_o_nn3d(i,jjm1,kkm1,q1) + &
                     c08_intu*alocont_o_nn3d(i,jjm1,kkm1,q7) + &
                     c20_intu*alocont_o_nn3d(i,jjm1,kkm1,q19) + &
                     c26_intu*alocont_o_nn3d(i,jjm1,kkm1,q25)))
                  alocont_o_nn3d(iim1,jjm1,kkm1,q27) = imask_totreg3d(iim1,jjm1,kkm1) * &
                     (alo_u*c14_scontu + abs_sc*(c10_intu*alocont_o_nn3d(iim1,jjm1,kkm1,q10) + &
                     c11_intu*alocont_o_nn3d(iim1,jjm1,kkm1,q11) + &
                     c12_intu*alocont_o_nn3d(iim1,jjm1,kkm1,q12) + &
                     c13_intu*alocont_o_nn3d(iim1,jjm1,kkm1,q13) + &
                     c14_intu*alocont_o_nn3d(iim1,jjm1,kkm1,q14) + &
                     c15_intu*alocont_o_nn3d(iim1,jjm1,kkm1,q15) + &
                     c16_intu*alocont_o_nn3d(iim1,jjm1,kkm1,q16) + &
                     c17_intu*alocont_o_nn3d(iim1,jjm1,kkm1,q17) + &
                     c18_intu*alocont_o_nn3d(iim1,jjm1,kkm1,q18) + &
                     c04_intu*alocont_o_nn3d(iim1,jjm1,kkm1,q4) + &
                     c05_intu*alocont_o_nn3d(iim1,jjm1,kkm1,q5) + &
                     c06_intu*alocont_o_nn3d(iim1,jjm1,kkm1,q6) + &
                     c22_intu*alocont_o_nn3d(iim1,jjm1,kkm1,q22) + &
                     c23_intu*alocont_o_nn3d(iim1,jjm1,kkm1,q23) + &
                     c24_intu*alocont_o_nn3d(iim1,jjm1,kkm1,q24) + &
                     c02_intu*alocont_o_nn3d(iim1,jjm1,kkm1,q2) + &
                     c08_intu*alocont_o_nn3d(iim1,jjm1,kkm1,q8) + &
                     c20_intu*alocont_o_nn3d(iim1,jjm1,kkm1,q20) + &
                     c26_intu*alocont_o_nn3d(iim1,jjm1,kkm1,q26)))
!                  call cpu_time(te_aloo)
!                  tt_aloo=tt_aloo+te_aloo-ts_aloo

!
!perform angular integration
!*********debug start
!               dels_r=dist_sphere(nn_x,nn_y,nn_z,x(i),y(j),z(k))
!               if(dels_r.ne.1.d10) then
!                  int3d(i,j,k)=xic1
!               else
!                  int3d(i,j,k)=0.d0
!               endif
!*********debug end
!                  call cpu_time(ts_integ)
                  mint3d_tmp(i,j,k) = mint3d_tmp(i,j,k) + int3d(i,j,k)*wall
                  fcontx3d_tmp(i,j,k) = fcontx3d_tmp(i,j,k) + int3d(i,j,k)*nn_x*wall
                  fconty3d_tmp(i,j,k) = fconty3d_tmp(i,j,k) + int3d(i,j,k)*nn_y*wall
                  fcontz3d_tmp(i,j,k) = fcontz3d_tmp(i,j,k) + int3d(i,j,k)*nn_z*wall
                  kcontxx3d_tmp(i,j,k) = kcontxx3d_tmp(i,j,k) + int3d(i,j,k)*nn_x*nn_x*wall
                  kcontyy3d_tmp(i,j,k) = kcontyy3d_tmp(i,j,k) + int3d(i,j,k)*nn_y*nn_y*wall
                  kcontzz3d_tmp(i,j,k) = kcontzz3d_tmp(i,j,k) + int3d(i,j,k)*nn_z*nn_z*wall
                  kcontxy3d_tmp(i,j,k) = kcontxy3d_tmp(i,j,k) + int3d(i,j,k)*nn_x*nn_y*wall
                  kcontxz3d_tmp(i,j,k) = kcontxz3d_tmp(i,j,k) + int3d(i,j,k)*nn_x*nn_z*wall
                  kcontyz3d_tmp(i,j,k) = kcontyz3d_tmp(i,j,k) + int3d(i,j,k)*nn_y*nn_z*wall
                  normalization3d_tmp(i,j,k) = normalization3d_tmp(i,j,k) + wall
!
                  alocont_nn3d_tmp(iip1,jjp1,kkp1,q1) = alocont_nn3d_tmp(iip1,jjp1,kkp1,q1) + wall*alocont_o_nn3d(iip1,jjp1,kkp1,q1)
                  alocont_nn3d_tmp(i,jjp1,kkp1,q2) = alocont_nn3d_tmp(i,jjp1,kkp1,q2) + wall*alocont_o_nn3d(i,jjp1,kkp1,q2)
                  alocont_nn3d_tmp(iim1,jjp1,kkp1,q3) = alocont_nn3d_tmp(iim1,jjp1,kkp1,q3) + wall*alocont_o_nn3d(iim1,jjp1,kkp1,q3)
                  alocont_nn3d_tmp(iip1,j,kkp1,q4) = alocont_nn3d_tmp(iip1,j,kkp1,q4) + wall*alocont_o_nn3d(iip1,j,kkp1,q4)
                  alocont_nn3d_tmp(i,j,kkp1,q5) = alocont_nn3d_tmp(i,j,kkp1,q5) + wall*alocont_o_nn3d(i,j,kkp1,q5)
                  alocont_nn3d_tmp(iim1,j,kkp1,q6) = alocont_nn3d_tmp(iim1,j,kkp1,q6) + wall*alocont_o_nn3d(iim1,j,kkp1,q6)
                  alocont_nn3d_tmp(iip1,jjm1,kkp1,q7) = alocont_nn3d_tmp(iip1,jjm1,kkp1,q7) + wall*alocont_o_nn3d(iip1,jjm1,kkp1,q7)
                  alocont_nn3d_tmp(i,jjm1,kkp1,q8) = alocont_nn3d_tmp(i,jjm1,kkp1,q8) + wall*alocont_o_nn3d(i,jjm1,kkp1,q8)
                  alocont_nn3d_tmp(iim1,jjm1,kkp1,q9) = alocont_nn3d_tmp(iim1,jjm1,kkp1,q9) + wall*alocont_o_nn3d(iim1,jjm1,kkp1,q9)
                  alocont_nn3d_tmp(iip1,jjp1,k,q10) = alocont_nn3d_tmp(iip1,jjp1,k,q10) + wall*alocont_o_nn3d(iip1,jjp1,k,q10)
                  alocont_nn3d_tmp(i,jjp1,k,q11) = alocont_nn3d_tmp(i,jjp1,k,q11) + wall*alocont_o_nn3d(i,jjp1,k,q11)
                  alocont_nn3d_tmp(iim1,jjp1,k,q12) = alocont_nn3d_tmp(iim1,jjp1,k,q12) + wall*alocont_o_nn3d(iim1,jjp1,k,q12)
                  alocont_nn3d_tmp(iip1,j,k,q13) = alocont_nn3d_tmp(iip1,j,k,q13) + wall*alocont_o_nn3d(iip1,j,k,q13)
                  alocont_nn3d_tmp(i,j,k,q14) = alocont_nn3d_tmp(i,j,k,q14) + wall*alocont_o_nn3d(i,j,k,q14)
                  alocont_nn3d_tmp(iim1,j,k,q15) = alocont_nn3d_tmp(iim1,j,k,q15) + wall*alocont_o_nn3d(iim1,j,k,q15)
                  alocont_nn3d_tmp(iip1,jjm1,k,q16) = alocont_nn3d_tmp(iip1,jjm1,k,q16) + wall*alocont_o_nn3d(iip1,jjm1,k,q16)
                  alocont_nn3d_tmp(i,jjm1,k,q17) = alocont_nn3d_tmp(i,jjm1,k,q17) + wall*alocont_o_nn3d(i,jjm1,k,q17)
                  alocont_nn3d_tmp(iim1,jjm1,k,q18) = alocont_nn3d_tmp(iim1,jjm1,k,q18) + wall*alocont_o_nn3d(iim1,jjm1,k,q18)
                  alocont_nn3d_tmp(iip1,jjp1,kkm1,q19) = alocont_nn3d_tmp(iip1,jjp1,kkm1,q19) + wall*alocont_o_nn3d(iip1,jjp1,kkm1,q19)
                  alocont_nn3d_tmp(i,jjp1,kkm1,q20) = alocont_nn3d_tmp(i,jjp1,kkm1,q20) + wall*alocont_o_nn3d(i,jjp1,kkm1,q20)
                  alocont_nn3d_tmp(iim1,jjp1,kkm1,q21) = alocont_nn3d_tmp(iim1,jjp1,kkm1,q21) + wall*alocont_o_nn3d(iim1,jjp1,kkm1,q21)
                  alocont_nn3d_tmp(iip1,j,kkm1,q22) = alocont_nn3d_tmp(iip1,j,kkm1,q22) + wall*alocont_o_nn3d(iip1,j,kkm1,q22)
                  alocont_nn3d_tmp(i,j,kkm1,q23) = alocont_nn3d_tmp(i,j,kkm1,q23) + wall*alocont_o_nn3d(i,j,kkm1,q23)
                  alocont_nn3d_tmp(iim1,j,kkm1,q24) = alocont_nn3d_tmp(iim1,j,kkm1,q24) + wall*alocont_o_nn3d(iim1,j,kkm1,q24)
                  alocont_nn3d_tmp(iip1,jjm1,kkm1,q25) = alocont_nn3d_tmp(iip1,jjm1,kkm1,q25) + wall*alocont_o_nn3d(iip1,jjm1,kkm1,q25)
                  alocont_nn3d_tmp(i,jjm1,kkm1,q26) = alocont_nn3d_tmp(i,jjm1,kkm1,q26) + wall*alocont_o_nn3d(i,jjm1,kkm1,q26)
                  alocont_nn3d_tmp(iim1,jjm1,kkm1,q27) = alocont_nn3d_tmp(iim1,jjm1,kkm1,q27) + wall*alocont_o_nn3d(iim1,jjm1,kkm1,q27)
!                  call cpu_time(te_integ)
!                  tt_integ=tt_integ+te_integ-ts_integ
               endif
!
!               call cpu_time(te_case3)
!               tt_case3=tt_case3+te_case3-ts_case3
!
             case default
!
            end select
!
!
         enddo
      enddo
   enddo
!
!
!
end subroutine fsc_cont3d
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine fsc_cont3d_lin(oindx)
!
!-----------------------------------------------------------------------
!------short characteristics for continuum radiative transfer in 3d-----
!-----------calculating intensties for given mu,phi specified-----------
!---------------------------by input oindx------------------------------
!------------------only linear interpolations are used------------------
!-----------------------------------------------------------------------
!
   use prog_type

   use fund_const
   use dime3d, only: ndxmax, ndzmax, ndymax, x, y, z
   use dime3d, only: int3d, opac3d, scont3d, imask3d, imask_totreg3d, &
      alocont_o_nn3d, alocont_nn3d_tmp, mint3d_tmp, normalization3d_tmp, &
      fcontx3d_tmp, fconty3d_tmp, fcontz3d_tmp, &
      kcontxx3d_tmp, kcontyy3d_tmp, kcontzz3d_tmp, kcontxy3d_tmp, kcontxz3d_tmp, kcontyz3d_tmp
   use angles, only: n_x, n_y, n_z, weight_omega, q_alo
   use bcondition, only: xic1
   use mod_debug, only: indxx, indxy, indxz
   use params_stellar, only: smajorax_a, smajorax_b, smajorax_c
!
   implicit none
!
! ... arguments
   integer(i4b), intent(in) :: oindx
!
! ... local scalars
   integer(i4b) :: i, j, k
   integer(i4b) :: alpha, beta, gamma
   integer(i4b) :: startx, starty, startz, endx, endy, endz
   integer(i4b) :: iim2, iim1, ii, iip1, jjm2, jjm1, jj, jjp1, kkm2, kkm1, kk, kkp1
   real(dp) :: nn_x, nn_y, nn_z, mueff, wall
   real(dp) :: dels_xyu, dels_xzu, dels_yzu, dels_u, dels_r
   real(dp) :: dels_xyd, dels_xzd, dels_yzd, dels_d
   real(dp) :: opac_p, scont_p
   real(dp) :: x_u, y_u, z_u, int_u, opac_u, scont_u
   real(dp) :: x_d, y_d, z_d, opac_d, scont_d
   real(dp) :: abs_sc, int_sc, contr_sc
   real(dp) :: alo_u, alo_p, alo_d
   integer :: q14, q15, q17, q18, q23, q24, q26, q27
   real(dp) :: c14_scontu, c15_scontu, c17_scontu, c18_scontu, &
      c23_scontu, c24_scontu, c26_scontu, c27_scontu, &
      c14_intu, c15_intu, c17_intu, c18_intu, &
      c23_intu, c24_intu, c26_intu, &
      c15_scontd, c17_scontd, c18_scontd, &
      c23_scontd, c24_scontd, c26_scontd, c27_scontd
!
!for debugging
   real(dp) :: int_u2, scont_u2, opac_u2, int_u3, scont_u3, opac_u3, scont_d2, scont_d3, opac_d2, opac_d3, interpol2d_9p_quad, interpol2d_9p_bez
!
! ... local arrays
!
! ... local functions
   real(dp) :: dist_sphere, dist_ellipsoid, calc_icore_gdark
!
! ... local characters
!character(len=50) :: enter
!
! ... local logicals
!
!directions
   nn_x=n_x(oindx)
   nn_y=n_y(oindx)
   nn_z=n_z(oindx)
!
!angulare integration weight
   wall=weight_omega(oindx)
!
!indices for nearest neighbour alo
   q14=q_alo(oindx,14)
   q15=q_alo(oindx,15)
   q17=q_alo(oindx,17)
   q18=q_alo(oindx,18)
   q23=q_alo(oindx,23)
   q24=q_alo(oindx,24)
   q26=q_alo(oindx,26)
   q27=q_alo(oindx,27)
!
!-----------------------------------------------------------------------
!
!set directional index-parameter (phantom points are excluded from calculation)
!
!index parameter:
!         if n_x >= 0                 if n_x < 0
!                startx = 2                  startx = ndxmax-1
!                endx = ndxmax-1             endx = 2
!                alpha=  1                   alpha=-1
!
!         if n_y >= 0                 if n_y < 0
!                starty = 2                  starty = ndymax-1
!                endy = ndymax-1             endy = 2
!                beta =  1                   beta =-1
!
!         if n_z >= 0                 if n_z < 0
!                startz = 2                  startz = ndzmax-1
!                endz = ndzmax-1             endz = 2
!                gamma = 1                   gamma = -1
!
   if(nn_x.gt.0.d0) then
      startx = 3
      endx = ndxmax-2
      alpha=  1
   elseif(nn_x.lt.0.d0) then
      startx = ndxmax-2
      endx = 3
      alpha=-1
   else
      stop 'error in fsc_cont3d: n_x = 0 not allowed'
   endif
!
   if(nn_y.gt.0.d0) then
      starty = 3
      endy = ndymax-2
      beta =  1
   elseif(nn_y.lt.0.d0) then
      starty = ndymax-2
      endy = 3
      beta =-1
   else
      stop 'error in fsc_cont3d: n_y = 0 not allowed'
   endif
!
   if(nn_z.gt.0.d0) then
      startz = 3
      endz = ndzmax-2
      gamma=  1
   elseif(nn_z.lt.0.d0) then
      startz = ndzmax-2
      endz = 3
      gamma=-1
   else
      stop 'error in fsc_cont3d: n_z = 0 not allowed'
   endif
!
!--------------------reset the intensities and alo----------------------
!
   call set_boundary3d(nn_x, nn_y, nn_z)
!
!for alo-tests
!scont3d=0.d0
!scont3d(indxx,indxy,indxz)=1.d0
!int3d=0.d0
!
   alocont_o_nn3d=0.d0
!
!-----------------------------------------------------------------------
!
   do k=startz, endz, gamma
      do j=starty, endy, beta
         do i=startx, endx, alpha
!
            select case(imask3d(i,j,k))
!
!************************standard radiative transfer********************
!
             case(1)
!
               iip1=i+alpha
               jjp1=j+beta
               kkp1=k+gamma
               iim1=i-alpha
               jjm1=j-beta
               kkm1=k-gamma
!
!calculate distance of upwind point to (previous) xy-plane, xz-plane and yz-plane
               dels_xyu=(z(k)-z(kkm1))/nn_z
               dels_xzu=(y(j)-y(jjm1))/nn_y
               dels_yzu=(x(i)-x(iim1))/nn_x
               dels_u=min(dels_xyu,dels_xzu,dels_yzu)
!
!calculate distance of downwind point to (next) xy-plane, xz-plane and yz-plane
               dels_xyd=(z(kkp1)-z(k))/nn_z
               dels_xzd=(y(jjp1)-y(j))/nn_y
               dels_yzd=(x(iip1)-x(i))/nn_x
               dels_d=min(dels_xyd,dels_xzd,dels_yzd)
!
!----------------------------local point--------------------------------
!
               scont_p=scont3d(i,j,k)
               opac_p=opac3d(i,j,k)
!
!----------------------------upwind point-------------------------------
!
               if(dels_xyu.eq.dels_u) then
!intersection with x-y plane on level k-gamma
                  x_u = x(i) - dels_u*nn_x
                  y_u = y(j) - dels_u*nn_y
                  z_u = z(kkm1)
!
                  call coeff2d_contu_lin(opac3d(iim1,jjm1,kkm1), opac3d(i,jjm1,kkm1), opac3d(iim1,j,kkm1), opac3d(i,j,kkm1), &
                     scont3d(iim1,jjm1,kkm1), scont3d(i,jjm1,kkm1), scont3d(iim1,j,kkm1), scont3d(i,j,kkm1), &
                     int3d(iim1,jjm1,kkm1), int3d(i,jjm1,kkm1), int3d(iim1,j,kkm1), int3d(i,j,kkm1), &
                     x(iim1), x(i), y(jjm1), y(j), x_u, y_u, &
                     c14_scontu, c15_scontu, c17_scontu, c18_scontu, &
                     c14_intu, c15_intu, c17_intu, c18_intu, opac_u, scont_u, int_u)
                  c23_scontu=0.d0
                  c24_scontu=0.d0
                  c26_scontu=0.d0
                  c23_intu=0.d0
                  c24_intu=0.d0
                  c26_intu=0.d0
!
               elseif(dels_xzu.eq.dels_u) then
!intersection with x-z plane on level j-beta
                  x_u = x(i) - dels_u*nn_x
                  y_u = y(jjm1)
                  z_u = z(k) - dels_u*nn_z
!
                  call coeff2d_contu_lin(opac3d(iim1,jjm1,kkm1), opac3d(i,jjm1,kkm1), opac3d(iim1,jjm1,k), opac3d(i,jjm1,k), &
                     scont3d(iim1,jjm1,kkm1), scont3d(i,jjm1,kkm1), scont3d(iim1,jjm1,k), scont3d(i,jjm1,k), &
                     int3d(iim1,jjm1,kkm1), int3d(i,jjm1,kkm1), int3d(iim1,jjm1,k), int3d(i,jjm1,k), &
                     x(iim1), x(i), z(kkm1), z(k), x_u, z_u, &
                     c14_scontu, c15_scontu, c23_scontu, c24_scontu, &
                     c14_intu, c15_intu, c23_intu, c24_intu, opac_u, scont_u, int_u)
                  c17_scontu=0.d0
                  c18_scontu=0.d0
                  c26_scontu=0.d0
                  c17_intu=0.d0
                  c18_intu=0.d0
                  c26_intu=0.d0
!
               elseif(dels_yzu.eq.dels_u) then
!intersection with y-z plane on level i-alpha
                  x_u = x(iim1)
                  y_u = y(j) - dels_u*nn_y
                  z_u = z(k) - dels_u*nn_z
!
                  call coeff2d_contu_lin(opac3d(iim1,jjm1,kkm1), opac3d(iim1,j,kkm1), opac3d(iim1,jjm1,k), opac3d(iim1,j,k), &
                     scont3d(iim1,jjm1,kkm1), scont3d(iim1,j,kkm1), scont3d(iim1,jjm1,k), scont3d(iim1,j,k), &
                     int3d(iim1,jjm1,kkm1), int3d(iim1,j,kkm1), int3d(iim1,jjm1,k), int3d(iim1,j,k), &
                     y(jjm1), y(j), z(kkm1), z(k), y_u, z_u, &
                     c14_scontu, c17_scontu, c23_scontu, c26_scontu, &
                     c14_intu, c17_intu, c23_intu, c26_intu, opac_u, scont_u, int_u)
                  c15_scontu=0.d0
                  c18_scontu=0.d0
                  c24_scontu=0.d0
                  c15_intu=0.d0
                  c18_intu=0.d0
                  c24_intu=0.d0
!
               else
                  write(*,'(4es20.8,2l4)') dels_u, dels_xzu, dels_xyu, dels_yzu, dels_u.eq.dels_xzu, dels_xyu.eq.dels_u
                  stop 'error in fsc_cont3d: invalid dels_u'
               endif
!
!---------------------------downwind point------------------------------
!
               if(dels_xyd.eq.dels_d) then
!intersection with x-y plane on level k+gamma
                  x_d = x(i) + dels_d*nn_x
                  y_d = y(j) + dels_d*nn_y
                  z_d = z(kkp1)
!
                  call coeff2d_contd_lin(opac3d(i,j,kkp1), opac3d(iip1,j,kkp1), opac3d(i,jjp1,kkp1), opac3d(iip1,jjp1,kkp1), &
                     scont3d(i,j,kkp1), scont3d(iip1,j,kkp1), scont3d(i,jjp1,kkp1), scont3d(iip1,jjp1,kkp1), &
                     x(i), x(iip1), y(j), y(jjp1), x_d, y_d, &
                     c23_scontd, c24_scontd, c26_scontd, c27_scontd, opac_d, scont_d)
!
               elseif(dels_xzd.eq.dels_d) then
!intersection with x-z plane on level j+beta
                  x_d = x(i) + dels_d*nn_x
                  y_d = y(jjp1)
                  z_d = z(k) + dels_d*nn_z
!
                  call coeff2d_contd_lin(opac3d(i,jjp1,k), opac3d(iip1,jjp1,k), opac3d(i,jjp1,kkp1), opac3d(iip1,jjp1,kkp1), &
                     scont3d(i,jjp1,k), scont3d(iip1,jjp1,k), scont3d(i,jjp1,kkp1), scont3d(iip1,jjp1,kkp1), &
                     x(i), x(iip1), z(k), z(kkp1), x_d, z_d, &
                     c17_scontd, c18_scontd, c26_scontd, c27_scontd, opac_d, scont_d)
!
               elseif(dels_yzd.eq.dels_d) then
!intersection with y-z plane on level i+alpha
                  x_d = x(iip1)
                  y_d = y(j) + dels_d*nn_y
                  z_d = z(k) + dels_d*nn_z
!
                  call coeff2d_contd_lin(opac3d(iip1,j,k), opac3d(iip1,jjp1,k), opac3d(iip1,j,kkp1), opac3d(iip1,jjp1,kkp1), &
                     scont3d(iip1,j,k), scont3d(iip1,jjp1,k), scont3d(iip1,j,kkp1), scont3d(iip1,jjp1,kkp1), &
                     y(j), y(jjp1), z(k), z(kkp1), y_d, z_d, &
                     c15_scontd, c18_scontd, c24_scontd, c27_scontd, opac_d, scont_d)
!
               else
                  write(*,'(4es20.8,l4)') dels_d, dels_xzd, dels_xyd, dels_yzd, dels_d.eq.dels_xzd
                  stop 'error in fsc_cont3d: invalid dels_d'
               endif
!
!--------------------------------radiative transfer---------------------
!
               call fsc_cont_lin(int_u, opac_u, opac_p, opac_d, scont_u, scont_p, scont_d, &
                  dels_u, dels_d, abs_sc, contr_sc, int_sc, alo_u, alo_p)
               int3d(i,j,k) = abs_sc*int_u + alo_u*scont_u + alo_p*scont_p
!
               alocont_o_nn3d(i,j,k,q14) = imask_totreg3d(i,j,k) * alo_p
               alocont_o_nn3d(iim1,j,k,q15) = imask_totreg3d(iim1,j,k) * &
                  (alo_u*c26_scontu + abs_sc*c26_intu*alocont_o_nn3d(iim1,j,k,q14))

               alocont_o_nn3d(i,jjm1,k,q17) = imask_totreg3d(i,jjm1,k) * &
                  (alo_u*c24_scontu + abs_sc*c24_intu*alocont_o_nn3d(i,jjm1,k,q14))
               alocont_o_nn3d(iim1,jjm1,k,q18) = imask_totreg3d(iim1,jjm1,k) * &
                  (alo_u*c23_scontu + abs_sc*(c23_intu*alocont_o_nn3d(iim1,jjm1,k,q14) + &
                  c26_intu*alocont_o_nn3d(iim1,jjm1,k,q17)))
               alocont_o_nn3d(i,j,kkm1,q23) = imask_totreg3d(i,j,kkm1) * &
                  (alo_u*c18_scontu + abs_sc*c18_intu*alocont_o_nn3d(i,j,kkm1,q14))
               alocont_o_nn3d(iim1,j,kkm1,q24) = imask_totreg3d(iim1,j,kkm1) * &
                  (alo_u*c17_scontu + abs_sc*(c17_intu*alocont_o_nn3d(iim1,j,kkm1,q14) + &
                  c18_intu*alocont_o_nn3d(iim1,j,kkm1,q15) + &
                  c26_intu*alocont_o_nn3d(iim1,j,kkm1,q23)))
               alocont_o_nn3d(i,jjm1,kkm1,q26) = imask_totreg3d(i,jjm1,kkm1) * &
                  (alo_u*c15_scontu + abs_sc*(c15_intu*alocont_o_nn3d(i,jjm1,kkm1,q14) + &
                  c18_intu*alocont_o_nn3d(i,jjm1,kkm1,q17) + &
                  c24_intu*alocont_o_nn3d(i,jjm1,kkm1,q23)))
               alocont_o_nn3d(iim1,jjm1,kkm1,q27) = imask_totreg3d(iim1,jjm1,kkm1) * &
                  (alo_u*c14_scontu + abs_sc*(c14_intu*alocont_o_nn3d(iim1,jjm1,kkm1,q14) + &
                  c15_intu*alocont_o_nn3d(iim1,jjm1,kkm1,q15) + &
                  c17_intu*alocont_o_nn3d(iim1,jjm1,kkm1,q17) + &
                  c18_intu*alocont_o_nn3d(iim1,jjm1,kkm1,q18) + &
                  c23_intu*alocont_o_nn3d(iim1,jjm1,kkm1,q23) + &
                  c24_intu*alocont_o_nn3d(iim1,jjm1,kkm1,q24) + &
                  c26_intu*alocont_o_nn3d(iim1,jjm1,kkm1,q26)))
!
!perform angular integration
!*********debug start
!               dels_r=dist_sphere(nn_x,nn_y,nn_z,x(i),y(j),z(k))
!               if(dels_r.ne.1.d10) then
!                  int3d(i,j,k)=xic1
!               else
!                  int3d(i,j,k)=0.d0
!               endif
!*********debug end
!
               mint3d_tmp(i,j,k) = mint3d_tmp(i,j,k) + int3d(i,j,k)*wall
               fcontx3d_tmp(i,j,k) = fcontx3d_tmp(i,j,k) + int3d(i,j,k)*nn_x*wall
               fconty3d_tmp(i,j,k) = fconty3d_tmp(i,j,k) + int3d(i,j,k)*nn_y*wall
               fcontz3d_tmp(i,j,k) = fcontz3d_tmp(i,j,k) + int3d(i,j,k)*nn_z*wall
               kcontxx3d_tmp(i,j,k) = kcontxx3d_tmp(i,j,k) + int3d(i,j,k)*nn_x*nn_x*wall
               kcontyy3d_tmp(i,j,k) = kcontyy3d_tmp(i,j,k) + int3d(i,j,k)*nn_y*nn_y*wall
               kcontzz3d_tmp(i,j,k) = kcontzz3d_tmp(i,j,k) + int3d(i,j,k)*nn_z*nn_z*wall
               kcontxy3d_tmp(i,j,k) = kcontxy3d_tmp(i,j,k) + int3d(i,j,k)*nn_x*nn_y*wall
               kcontxz3d_tmp(i,j,k) = kcontxz3d_tmp(i,j,k) + int3d(i,j,k)*nn_x*nn_z*wall
               kcontyz3d_tmp(i,j,k) = kcontyz3d_tmp(i,j,k) + int3d(i,j,k)*nn_y*nn_z*wall
               normalization3d_tmp(i,j,k) = normalization3d_tmp(i,j,k) + wall
!
               alocont_nn3d_tmp(i,j,k,q14) = alocont_nn3d_tmp(i,j,k,q14) + wall*alocont_o_nn3d(i,j,k,q14)
               alocont_nn3d_tmp(iim1,j,k,q15) = alocont_nn3d_tmp(iim1,j,k,q15) + wall*alocont_o_nn3d(iim1,j,k,q15)
               alocont_nn3d_tmp(i,jjm1,k,q17) = alocont_nn3d_tmp(i,jjm1,k,q17) + wall*alocont_o_nn3d(i,jjm1,k,q17)
               alocont_nn3d_tmp(iim1,jjm1,k,q18) = alocont_nn3d_tmp(iim1,jjm1,k,q18) + wall*alocont_o_nn3d(iim1,jjm1,k,q18)
               alocont_nn3d_tmp(i,j,kkm1,q23) = alocont_nn3d_tmp(i,j,kkm1,q23) + wall*alocont_o_nn3d(i,j,kkm1,q23)
               alocont_nn3d_tmp(iim1,j,kkm1,q24) = alocont_nn3d_tmp(iim1,j,kkm1,q24) + wall*alocont_o_nn3d(iim1,j,kkm1,q24)
               alocont_nn3d_tmp(i,jjm1,kkm1,q26) = alocont_nn3d_tmp(i,jjm1,kkm1,q26) + wall*alocont_o_nn3d(i,jjm1,kkm1,q26)
               alocont_nn3d_tmp(iim1,jjm1,kkm1,q27) = alocont_nn3d_tmp(iim1,jjm1,kkm1,q27) + wall*alocont_o_nn3d(iim1,jjm1,kkm1,q27)
!
!***********special treatment: point is in vicinity of boundary*********
!
             case(2)
!
               iip1=i+alpha
               jjp1=j+beta
               kkp1=k+gamma
               iim1=i-alpha
               jjm1=j-beta
               kkm1=k-gamma
!
!calculate distance of upwind point to (previous) xy-plane, xz-plane and yz-plane
               dels_xyu=(z(k)-z(kkm1))/nn_z
               dels_xzu=(y(j)-y(jjm1))/nn_y
               dels_yzu=(x(i)-x(iim1))/nn_x
!
!calculate distance of downwind point to (next) xy-plane, xz-plane and yz-plane
               dels_xyd=(z(kkp1)-z(k))/nn_z
               dels_xzd=(y(jjp1)-y(j))/nn_y
               dels_yzd=(x(iip1)-x(i))/nn_x
!
!calculate distance to the boundary (here: sphere)
!               dels_r=dist_sphere(nn_x,nn_y,nn_z,x(i),y(j),z(k))
               dels_r=dist_ellipsoid(nn_x,nn_y,nn_z,x(i),y(j),z(k),smajorax_a,smajorax_b,smajorax_c)
!
               dels_u=min(dels_xyu,dels_xzu,dels_yzu,dels_r)
               dels_d=min(dels_xyd,dels_xzd,dels_yzd)
!
!----------------------------local point--------------------------------
!
               scont_p=scont3d(i,j,k)
               opac_p=opac3d(i,j,k)
!
!----------------------------upwind point-------------------------------
!
               if(dels_r.eq.dels_u) then
                  x_u = x(i) - dels_u*nn_x
                  y_u = y(j) - dels_u*nn_y
                  z_u = z(k) - dels_u*nn_z
!
!for alo
!                  int_u=0.d0
!                  int_u=xic1
                  int_u=calc_icore_gdark(z_u, sqrt(x_u**2+y_u**2+z_u**2))
!
                  call coeff3d_contu_lin(opac3d(iim1,jjm1,kkm1), opac3d(i,jjm1,kkm1), opac3d(iim1,j,kkm1), opac3d(i,j,kkm1), &
                     opac3d(iim1,jjm1,k),    opac3d(i,jjm1,k),    opac3d(iim1,j,k),    opac3d(i,j,k), &
                     scont3d(iim1,jjm1,kkm1), scont3d(i,jjm1,kkm1), scont3d(iim1,j,kkm1), scont3d(i,j,kkm1), &
                     scont3d(iim1,jjm1,k),    scont3d(i,jjm1,k),    scont3d(iim1,j,k),    scont3d(i,j,k), &
                     x(iim1), x(i), y(jjm1), y(j), z(kkm1), z(k), x_u, y_u, z_u, &
                     c14_scontu, c15_scontu, c17_scontu, c18_scontu, &
                     c23_scontu, c24_scontu, c26_scontu, c27_scontu, &
                     opac_u, scont_u)
!
               elseif(dels_xyu.eq.dels_u) then
!intersection with x-y plane on level k-gamma
                  x_u = x(i) - dels_u*nn_x
                  y_u = y(j) - dels_u*nn_y
                  z_u = z(kkm1)
!
                  call coeff2d_contu_lin(opac3d(iim1,jjm1,kkm1), opac3d(i,jjm1,kkm1), opac3d(iim1,j,kkm1), opac3d(i,j,kkm1), &
                     scont3d(iim1,jjm1,kkm1), scont3d(i,jjm1,kkm1), scont3d(iim1,j,kkm1), scont3d(i,j,kkm1), &
                     int3d(iim1,jjm1,kkm1), int3d(i,jjm1,kkm1), int3d(iim1,j,kkm1), int3d(i,j,kkm1), &
                     x(iim1), x(i), y(jjm1), y(j), x_u, y_u, &
                     c14_scontu, c15_scontu, c17_scontu, c18_scontu, &
                     c14_intu, c15_intu, c17_intu, c18_intu, opac_u, scont_u, int_u)
                  c23_scontu=0.d0
                  c24_scontu=0.d0
                  c26_scontu=0.d0
                  c27_scontu=0.d0
                  c23_intu=0.d0
                  c24_intu=0.d0
                  c26_intu=0.d0
!
               elseif(dels_xzu.eq.dels_u) then
!intersection with x-z plane on level j-beta
                  x_u = x(i) - dels_u*nn_x
                  y_u = y(jjm1)
                  z_u = z(k) - dels_u*nn_z
!
                  call coeff2d_contu_lin(opac3d(iim1,jjm1,kkm1), opac3d(i,jjm1,kkm1), opac3d(iim1,jjm1,k), opac3d(i,jjm1,k), &
                     scont3d(iim1,jjm1,kkm1), scont3d(i,jjm1,kkm1), scont3d(iim1,jjm1,k), scont3d(i,jjm1,k), &
                     int3d(iim1,jjm1,kkm1), int3d(i,jjm1,kkm1), int3d(iim1,jjm1,k), int3d(i,jjm1,k), &
                     x(iim1), x(i), z(kkm1), z(k), x_u, z_u, &
                     c14_scontu, c15_scontu, c23_scontu, c24_scontu, &
                     c14_intu, c15_intu, c23_intu, c24_intu, opac_u, scont_u, int_u)
                  c17_scontu=0.d0
                  c18_scontu=0.d0
                  c26_scontu=0.d0
                  c27_scontu=0.d0
                  c17_intu=0.d0
                  c18_intu=0.d0
                  c26_intu=0.d0
!
               elseif(dels_yzu.eq.dels_u) then
!intersection with y-z plane on level i-alpha
                  x_u = x(iim1)
                  y_u = y(j) - dels_u*nn_y
                  z_u = z(k) - dels_u*nn_z
!
                  call coeff2d_contu_lin(opac3d(iim1,jjm1,kkm1), opac3d(iim1,j,kkm1), opac3d(iim1,jjm1,k), opac3d(iim1,j,k), &
                     scont3d(iim1,jjm1,kkm1), scont3d(iim1,j,kkm1), scont3d(iim1,jjm1,k), scont3d(iim1,j,k), &
                     int3d(iim1,jjm1,kkm1), int3d(iim1,j,kkm1), int3d(iim1,jjm1,k), int3d(iim1,j,k), &
                     y(jjm1), y(j), z(kkm1), z(k), y_u, z_u, &
                     c14_scontu, c17_scontu, c23_scontu, c26_scontu, &
                     c14_intu, c17_intu, c23_intu, c26_intu, opac_u, scont_u, int_u)
                  c15_scontu=0.d0
                  c18_scontu=0.d0
                  c24_scontu=0.d0
                  c27_scontu=0.d0
                  c15_intu=0.d0
                  c18_intu=0.d0
                  c24_intu=0.d0
!
               else
                  write(*,'(4es20.8,2l4)') dels_u, dels_xzu, dels_xyu, dels_yzu, dels_u.eq.dels_xzu, dels_xyu.eq.dels_u
                  stop 'error in fsc_cont3d: invalid dels_u'
               endif
!
!---------------------------downwind point------------------------------
!
               if(dels_xyd.eq.dels_d) then
!intersection with x-y plane on level k+gamma
                  x_d = x(i) + dels_d*nn_x
                  y_d = y(j) + dels_d*nn_y
                  z_d = z(kkp1)
!
                  call coeff2d_contd_lin(opac3d(i,j,kkp1), opac3d(iip1,j,kkp1), opac3d(i,jjp1,kkp1), opac3d(iip1,jjp1,kkp1), &
                     scont3d(i,j,kkp1), scont3d(iip1,j,kkp1), scont3d(i,jjp1,kkp1), scont3d(iip1,jjp1,kkp1), &
                     x(i), x(iip1), y(j), y(jjp1), x_d, y_d, &
                     c23_scontd, c24_scontd, c26_scontd, c27_scontd, opac_d, scont_d)
!
               elseif(dels_xzd.eq.dels_d) then
!intersection with x-z plane on level j+beta
                  x_d = x(i) + dels_d*nn_x
                  y_d = y(jjp1)
                  z_d = z(k) + dels_d*nn_z
!
                  call coeff2d_contd_lin(opac3d(i,jjp1,k), opac3d(iip1,jjp1,k), opac3d(i,jjp1,kkp1), opac3d(iip1,jjp1,kkp1), &
                     scont3d(i,jjp1,k), scont3d(iip1,jjp1,k), scont3d(i,jjp1,kkp1), scont3d(iip1,jjp1,kkp1), &
                     x(i), x(iip1), z(k), z(kkp1), x_d, z_d, &
                     c17_scontd, c18_scontd, c26_scontd, c27_scontd, opac_d, scont_d)
!
               elseif(dels_yzd.eq.dels_d) then
!intersection with y-z plane on level i+alpha
                  x_d = x(iip1)
                  y_d = y(j) + dels_d*nn_y
                  z_d = z(k) + dels_d*nn_z
!
                  call coeff2d_contd_lin(opac3d(iip1,j,k), opac3d(iip1,jjp1,k), opac3d(iip1,j,kkp1), opac3d(iip1,jjp1,kkp1), &
                     scont3d(iip1,j,k), scont3d(iip1,jjp1,k), scont3d(iip1,j,kkp1), scont3d(iip1,jjp1,kkp1), &
                     y(j), y(jjp1), z(k), z(kkp1), y_d, z_d, &
                     c15_scontd, c18_scontd, c24_scontd, c27_scontd, opac_d, scont_d)
!
               else
                  write(*,'(4es20.8,l4)') dels_d, dels_xzd, dels_xyd, dels_yzd, dels_d.eq.dels_xzd
                  stop 'error in fsc_cont3d: invalid dels_d'
               endif
!--------------------------------radiative transfer---------------------
!

               call fsc_cont_lin(int_u, opac_u, opac_p, opac_d, scont_u, scont_p, scont_d, &
                  dels_u, dels_d, abs_sc, contr_sc, int_sc, alo_u, alo_p)
               int3d(i,j,k) = abs_sc*int_u + alo_u*scont_u + alo_p*scont_p
!
               alocont_o_nn3d(i,j,k,q14) = imask_totreg3d(i,j,k) * (alo_p + alo_u*c27_scontu)
               alocont_o_nn3d(iim1,j,k,q15) = imask_totreg3d(iim1,j,k) * &
                  (alo_u*c26_scontu + abs_sc*c26_intu*alocont_o_nn3d(iim1,j,k,q14))

               alocont_o_nn3d(i,jjm1,k,q17) = imask_totreg3d(i,jjm1,k) * &
                  (alo_u*c24_scontu + abs_sc*c24_intu*alocont_o_nn3d(i,jjm1,k,q14))
               alocont_o_nn3d(iim1,jjm1,k,q18) = imask_totreg3d(iim1,jjm1,k) * &
                  (alo_u*c23_scontu + abs_sc*(c23_intu*alocont_o_nn3d(iim1,jjm1,k,q14) + &
                  c26_intu*alocont_o_nn3d(iim1,jjm1,k,q17)))
               alocont_o_nn3d(i,j,kkm1,q23) = imask_totreg3d(i,j,kkm1) * &
                  (alo_u*c18_scontu + abs_sc*c18_intu*alocont_o_nn3d(i,j,kkm1,q14))
               alocont_o_nn3d(iim1,j,kkm1,q24) = imask_totreg3d(iim1,j,kkm1) * &
                  (alo_u*c17_scontu + abs_sc*(c17_intu*alocont_o_nn3d(iim1,j,kkm1,q14) + &
                  c18_intu*alocont_o_nn3d(iim1,j,kkm1,q15) + &
                  c26_intu*alocont_o_nn3d(iim1,j,kkm1,q23)))
               alocont_o_nn3d(i,jjm1,kkm1,q26) = imask_totreg3d(i,jjm1,kkm1) * &
                  (alo_u*c15_scontu + abs_sc*(c15_intu*alocont_o_nn3d(i,jjm1,kkm1,q14) + &
                  c18_intu*alocont_o_nn3d(i,jjm1,kkm1,q17) + &
                  c24_intu*alocont_o_nn3d(i,jjm1,kkm1,q23)))
               alocont_o_nn3d(iim1,jjm1,kkm1,q27) = imask_totreg3d(iim1,jjm1,kkm1) * &
                  (alo_u*c14_scontu + abs_sc*(c14_intu*alocont_o_nn3d(iim1,jjm1,kkm1,q14) + &
                  c15_intu*alocont_o_nn3d(iim1,jjm1,kkm1,q15) + &
                  c17_intu*alocont_o_nn3d(iim1,jjm1,kkm1,q17) + &
                  c18_intu*alocont_o_nn3d(iim1,jjm1,kkm1,q18) + &
                  c23_intu*alocont_o_nn3d(iim1,jjm1,kkm1,q23) + &
                  c24_intu*alocont_o_nn3d(iim1,jjm1,kkm1,q24) + &
                  c26_intu*alocont_o_nn3d(iim1,jjm1,kkm1,q26)))
!
!perform angular integration
!*********debug start
!               dels_r=dist_sphere(nn_x,nn_y,nn_z,x(i),y(j),z(k))
!               if(dels_r.ne.1.d10) then
!                  int3d(i,j,k)=xic1
!               else
!                  int3d(i,j,k)=0.d0
!               endif
!*********debug end
!
               mint3d_tmp(i,j,k) = mint3d_tmp(i,j,k) + int3d(i,j,k)*wall
               fcontx3d_tmp(i,j,k) = fcontx3d_tmp(i,j,k) + int3d(i,j,k)*nn_x*wall
               fconty3d_tmp(i,j,k) = fconty3d_tmp(i,j,k) + int3d(i,j,k)*nn_y*wall
               fcontz3d_tmp(i,j,k) = fcontz3d_tmp(i,j,k) + int3d(i,j,k)*nn_z*wall
               kcontxx3d_tmp(i,j,k) = kcontxx3d_tmp(i,j,k) + int3d(i,j,k)*nn_x*nn_x*wall
               kcontyy3d_tmp(i,j,k) = kcontyy3d_tmp(i,j,k) + int3d(i,j,k)*nn_y*nn_y*wall
               kcontzz3d_tmp(i,j,k) = kcontzz3d_tmp(i,j,k) + int3d(i,j,k)*nn_z*nn_z*wall
               kcontxy3d_tmp(i,j,k) = kcontxy3d_tmp(i,j,k) + int3d(i,j,k)*nn_x*nn_y*wall
               kcontxz3d_tmp(i,j,k) = kcontxz3d_tmp(i,j,k) + int3d(i,j,k)*nn_x*nn_z*wall
               kcontyz3d_tmp(i,j,k) = kcontyz3d_tmp(i,j,k) + int3d(i,j,k)*nn_y*nn_z*wall
               normalization3d_tmp(i,j,k) = normalization3d_tmp(i,j,k) + wall
!
               alocont_nn3d_tmp(i,j,k,q14) = alocont_nn3d_tmp(i,j,k,q14) + wall*alocont_o_nn3d(i,j,k,q14)
               alocont_nn3d_tmp(iim1,j,k,q15) = alocont_nn3d_tmp(iim1,j,k,q15) + wall*alocont_o_nn3d(iim1,j,k,q15)
               alocont_nn3d_tmp(i,jjm1,k,q17) = alocont_nn3d_tmp(i,jjm1,k,q17) + wall*alocont_o_nn3d(i,jjm1,k,q17)
               alocont_nn3d_tmp(iim1,jjm1,k,q18) = alocont_nn3d_tmp(iim1,jjm1,k,q18) + wall*alocont_o_nn3d(iim1,jjm1,k,q18)
               alocont_nn3d_tmp(i,j,kkm1,q23) = alocont_nn3d_tmp(i,j,kkm1,q23) + wall*alocont_o_nn3d(i,j,kkm1,q23)
               alocont_nn3d_tmp(iim1,j,kkm1,q24) = alocont_nn3d_tmp(iim1,j,kkm1,q24) + wall*alocont_o_nn3d(iim1,j,kkm1,q24)
               alocont_nn3d_tmp(i,jjm1,kkm1,q26) = alocont_nn3d_tmp(i,jjm1,kkm1,q26) + wall*alocont_o_nn3d(i,jjm1,kkm1,q26)
               alocont_nn3d_tmp(iim1,jjm1,kkm1,q27) = alocont_nn3d_tmp(iim1,jjm1,kkm1,q27) + wall*alocont_o_nn3d(iim1,jjm1,kkm1,q27)
!
!
!***********special treatment: point lies directly on boundary**********
!
             case(3)
               mueff=nn_x*x(i)+nn_y*y(j)+nn_z*z(k)
               if(mueff.ge.0.d0) then
!no interpolations required
!for alo
!                 int3d(i,j,k) = 0.d0
!                  int3d(i,j,k) = xic1
                  int3d(i,j,k)=calc_icore_gdark(z(k), sqrt(x(i)**2+y(j)**2+z(k)**2))
                  alocont_o_nn3d(i,j,k,14) = 0.d0
!
!perform angular integration (alo not required, since boundary specified)
                  mint3d_tmp(i,j,k) = mint3d_tmp(i,j,k) + int3d(i,j,k)*wall
                  fcontx3d_tmp(i,j,k) = fcontx3d_tmp(i,j,k) + int3d(i,j,k)*nn_x*wall
                  fconty3d_tmp(i,j,k) = fconty3d_tmp(i,j,k) + int3d(i,j,k)*nn_y*wall
                  fcontz3d_tmp(i,j,k) = fcontz3d_tmp(i,j,k) + int3d(i,j,k)*nn_z*wall
                  kcontxx3d_tmp(i,j,k) = kcontxx3d_tmp(i,j,k) + int3d(i,j,k)*nn_x*nn_x*wall
                  kcontyy3d_tmp(i,j,k) = kcontyy3d_tmp(i,j,k) + int3d(i,j,k)*nn_y*nn_y*wall
                  kcontzz3d_tmp(i,j,k) = kcontzz3d_tmp(i,j,k) + int3d(i,j,k)*nn_z*nn_z*wall
                  kcontxy3d_tmp(i,j,k) = kcontxy3d_tmp(i,j,k) + int3d(i,j,k)*nn_x*nn_y*wall
                  kcontxz3d_tmp(i,j,k) = kcontxz3d_tmp(i,j,k) + int3d(i,j,k)*nn_x*nn_z*wall
                  kcontyz3d_tmp(i,j,k) = kcontyz3d_tmp(i,j,k) + int3d(i,j,k)*nn_y*nn_z*wall
                  normalization3d_tmp(i,j,k) = normalization3d_tmp(i,j,k) + wall
               else
!same interpolation as in (1)
!
                  iip1=i+alpha
                  jjp1=j+beta
                  kkp1=k+gamma
                  iim1=i-alpha
                  jjm1=j-beta
                  kkm1=k-gamma
!
!calculate distance of upwind point to (previous) xy-plane, xz-plane and yz-plane
                  dels_xyu=(z(k)-z(kkm1))/nn_z
                  dels_xzu=(y(j)-y(jjm1))/nn_y
                  dels_yzu=(x(i)-x(iim1))/nn_x
                  dels_u=min(dels_xyu,dels_xzu,dels_yzu)
!
!calculate distance of downwind point to (next) xy-plane, xz-plane and yz-plane
                  dels_xyd=(z(kkp1)-z(k))/nn_z
                  dels_xzd=(y(jjp1)-y(j))/nn_y
                  dels_yzd=(x(iip1)-x(i))/nn_x
                  dels_d=min(dels_xyd,dels_xzd,dels_yzd)
!
!----------------------------local point--------------------------------
!
                  scont_p=scont3d(i,j,k)
                  opac_p=opac3d(i,j,k)
!
!----------------------------upwind point-------------------------------
!
                  if(dels_xyu.eq.dels_u) then
!intersection with x-y plane on level k-gamma
                     x_u = x(i) - dels_u*nn_x
                     y_u = y(j) - dels_u*nn_y
                     z_u = z(kkm1)
!
                     call coeff2d_contu_lin(opac3d(iim1,jjm1,kkm1), opac3d(i,jjm1,kkm1), opac3d(iim1,j,kkm1), opac3d(i,j,kkm1), &
                        scont3d(iim1,jjm1,kkm1), scont3d(i,jjm1,kkm1), scont3d(iim1,j,kkm1), scont3d(i,j,kkm1), &
                        int3d(iim1,jjm1,kkm1), int3d(i,jjm1,kkm1), int3d(iim1,j,kkm1), int3d(i,j,kkm1), &
                        x(iim1), x(i), y(jjm1), y(j), x_u, y_u, &
                        c14_scontu, c15_scontu, c17_scontu, c18_scontu, &
                        c14_intu, c15_intu, c17_intu, c18_intu, opac_u, scont_u, int_u)
                     c23_scontu=0.d0
                     c24_scontu=0.d0
                     c26_scontu=0.d0
                     c23_intu=0.d0
                     c24_intu=0.d0
                     c26_intu=0.d0
!
                  elseif(dels_xzu.eq.dels_u) then
!intersection with x-z plane on level j-beta
                     x_u = x(i) - dels_u*nn_x
                     y_u = y(jjm1)
                     z_u = z(k) - dels_u*nn_z
!
                     call coeff2d_contu_lin(opac3d(iim1,jjm1,kkm1), opac3d(i,jjm1,kkm1), opac3d(iim1,jjm1,k), opac3d(i,jjm1,k), &
                        scont3d(iim1,jjm1,kkm1), scont3d(i,jjm1,kkm1), scont3d(iim1,jjm1,k), scont3d(i,jjm1,k), &
                        int3d(iim1,jjm1,kkm1), int3d(i,jjm1,kkm1), int3d(iim1,jjm1,k), int3d(i,jjm1,k), &
                        x(iim1), x(i), z(kkm1), z(k), x_u, z_u, &
                        c14_scontu, c15_scontu, c23_scontu, c24_scontu, &
                        c14_intu, c15_intu, c23_intu, c24_intu, opac_u, scont_u, int_u)
                     c17_scontu=0.d0
                     c18_scontu=0.d0
                     c26_scontu=0.d0
                     c17_intu=0.d0
                     c18_intu=0.d0
                     c26_intu=0.d0
!
                  elseif(dels_yzu.eq.dels_u) then
!intersection with y-z plane on level i-alpha
                     x_u = x(iim1)
                     y_u = y(j) - dels_u*nn_y
                     z_u = z(k) - dels_u*nn_z
!
                     call coeff2d_contu_lin(opac3d(iim1,jjm1,kkm1), opac3d(iim1,j,kkm1), opac3d(iim1,jjm1,k), opac3d(iim1,j,k), &
                        scont3d(iim1,jjm1,kkm1), scont3d(iim1,j,kkm1), scont3d(iim1,jjm1,k), scont3d(iim1,j,k), &
                        int3d(iim1,jjm1,kkm1), int3d(iim1,j,kkm1), int3d(iim1,jjm1,k), int3d(iim1,j,k), &
                        y(jjm1), y(j), z(kkm1), z(k), y_u, z_u, &
                        c14_scontu, c17_scontu, c23_scontu, c26_scontu, &
                        c14_intu, c17_intu, c23_intu, c26_intu, opac_u, scont_u, int_u)
                     c15_scontu=0.d0
                     c18_scontu=0.d0
                     c24_scontu=0.d0
                     c15_intu=0.d0
                     c18_intu=0.d0
                     c24_intu=0.d0
!
                  else
                     write(*,'(4es20.8,2l4)') dels_u, dels_xzu, dels_xyu, dels_yzu, dels_u.eq.dels_xzu, dels_xyu.eq.dels_u
                     stop 'error in fsc_cont3d: invalid dels_u'
                  endif
!
!---------------------------downwind point------------------------------
!
                  if(dels_xyd.eq.dels_d) then
!intersection with x-y plane on level k+gamma
                     x_d = x(i) + dels_d*nn_x
                     y_d = y(j) + dels_d*nn_y
                     z_d = z(kkp1)
!
                     call coeff2d_contd_lin(opac3d(i,j,kkp1), opac3d(iip1,j,kkp1), opac3d(i,jjp1,kkp1), opac3d(iip1,jjp1,kkp1), &
                        scont3d(i,j,kkp1), scont3d(iip1,j,kkp1), scont3d(i,jjp1,kkp1), scont3d(iip1,jjp1,kkp1), &
                        x(i), x(iip1), y(j), y(jjp1), x_d, y_d, &
                        c23_scontd, c24_scontd, c26_scontd, c27_scontd, opac_d, scont_d)
!
                  elseif(dels_xzd.eq.dels_d) then
!intersection with x-z plane on level j+beta
                     x_d = x(i) + dels_d*nn_x
                     y_d = y(jjp1)
                     z_d = z(k) + dels_d*nn_z
!
                     call coeff2d_contd_lin(opac3d(i,jjp1,k), opac3d(iip1,jjp1,k), opac3d(i,jjp1,kkp1), opac3d(iip1,jjp1,kkp1), &
                        scont3d(i,jjp1,k), scont3d(iip1,jjp1,k), scont3d(i,jjp1,kkp1), scont3d(iip1,jjp1,kkp1), &
                        x(i), x(iip1), z(k), z(kkp1), x_d, z_d, &
                        c17_scontd, c18_scontd, c26_scontd, c27_scontd, opac_d, scont_d)
!
                  elseif(dels_yzd.eq.dels_d) then
!intersection with y-z plane on level i+alpha
                     x_d = x(iip1)
                     y_d = y(j) + dels_d*nn_y
                     z_d = z(k) + dels_d*nn_z
!
                     call coeff2d_contd_lin(opac3d(iip1,j,k), opac3d(iip1,jjp1,k), opac3d(iip1,j,kkp1), opac3d(iip1,jjp1,kkp1), &
                        scont3d(iip1,j,k), scont3d(iip1,jjp1,k), scont3d(iip1,j,kkp1), scont3d(iip1,jjp1,kkp1), &
                        y(j), y(jjp1), z(k), z(kkp1), y_d, z_d, &
                        c15_scontd, c18_scontd, c24_scontd, c27_scontd, opac_d, scont_d)
!
                  else
                     write(*,'(4es20.8,l4)') dels_d, dels_xzd, dels_xyd, dels_yzd, dels_d.eq.dels_xzd
                     stop 'error in fsc_cont3d: invalid dels_d'
                  endif
!
!--------------------------------radiative transfer---------------------
!
                  call fsc_cont_lin(int_u, opac_u, opac_p, opac_d, scont_u, scont_p, scont_d, &
                     dels_u, dels_d, abs_sc, contr_sc, int_sc, alo_u, alo_p)
                  int3d(i,j,k) = abs_sc*int_u + alo_u*scont_u + alo_p*scont_p
!
                  alocont_o_nn3d(i,j,k,q14) = imask_totreg3d(i,j,k) * alo_p
                  alocont_o_nn3d(iim1,j,k,q15) = imask_totreg3d(iim1,j,k) * &
                     (alo_u*c26_scontu + abs_sc*c26_intu*alocont_o_nn3d(iim1,j,k,q14))

                  alocont_o_nn3d(i,jjm1,k,q17) = imask_totreg3d(i,jjm1,k) * &
                     (alo_u*c24_scontu + abs_sc*c24_intu*alocont_o_nn3d(i,jjm1,k,q14))
                  alocont_o_nn3d(iim1,jjm1,k,q18) = imask_totreg3d(iim1,jjm1,k) * &
                     (alo_u*c23_scontu + abs_sc*(c23_intu*alocont_o_nn3d(iim1,jjm1,k,q14) + &
                     c26_intu*alocont_o_nn3d(iim1,jjm1,k,q17)))
                  alocont_o_nn3d(i,j,kkm1,q23) = imask_totreg3d(i,j,kkm1) * &
                     (alo_u*c18_scontu + abs_sc*c18_intu*alocont_o_nn3d(i,j,kkm1,q14))
                  alocont_o_nn3d(iim1,j,kkm1,q24) = imask_totreg3d(iim1,j,kkm1) * &
                     (alo_u*c17_scontu + abs_sc*(c17_intu*alocont_o_nn3d(iim1,j,kkm1,q14) + &
                     c18_intu*alocont_o_nn3d(iim1,j,kkm1,q15) + &
                     c26_intu*alocont_o_nn3d(iim1,j,kkm1,q23)))
                  alocont_o_nn3d(i,jjm1,kkm1,q26) = imask_totreg3d(i,jjm1,kkm1) * &
                     (alo_u*c15_scontu + abs_sc*(c15_intu*alocont_o_nn3d(i,jjm1,kkm1,q14) + &
                     c18_intu*alocont_o_nn3d(i,jjm1,kkm1,q17) + &
                     c24_intu*alocont_o_nn3d(i,jjm1,kkm1,q23)))
                  alocont_o_nn3d(iim1,jjm1,kkm1,q27) = imask_totreg3d(iim1,jjm1,kkm1) * &
                     (alo_u*c14_scontu + abs_sc*(c14_intu*alocont_o_nn3d(iim1,jjm1,kkm1,q14) + &
                     c15_intu*alocont_o_nn3d(iim1,jjm1,kkm1,q15) + &
                     c17_intu*alocont_o_nn3d(iim1,jjm1,kkm1,q17) + &
                     c18_intu*alocont_o_nn3d(iim1,jjm1,kkm1,q18) + &
                     c23_intu*alocont_o_nn3d(iim1,jjm1,kkm1,q23) + &
                     c24_intu*alocont_o_nn3d(iim1,jjm1,kkm1,q24) + &
                     c26_intu*alocont_o_nn3d(iim1,jjm1,kkm1,q26)))
!

!
!perform angular integration
!*********debug start
!               dels_r=dist_sphere(nn_x,nn_y,nn_z,x(i),y(j),z(k))
!               if(dels_r.ne.1.d10) then
!                  int3d(i,j,k)=xic1
!               else
!                  int3d(i,j,k)=0.d0
!               endif
!*********debug end
                  mint3d_tmp(i,j,k) = mint3d_tmp(i,j,k) + int3d(i,j,k)*wall
                  fcontx3d_tmp(i,j,k) = fcontx3d_tmp(i,j,k) + int3d(i,j,k)*nn_x*wall
                  fconty3d_tmp(i,j,k) = fconty3d_tmp(i,j,k) + int3d(i,j,k)*nn_y*wall
                  fcontz3d_tmp(i,j,k) = fcontz3d_tmp(i,j,k) + int3d(i,j,k)*nn_z*wall
                  kcontxx3d_tmp(i,j,k) = kcontxx3d_tmp(i,j,k) + int3d(i,j,k)*nn_x*nn_x*wall
                  kcontyy3d_tmp(i,j,k) = kcontyy3d_tmp(i,j,k) + int3d(i,j,k)*nn_y*nn_y*wall
                  kcontzz3d_tmp(i,j,k) = kcontzz3d_tmp(i,j,k) + int3d(i,j,k)*nn_z*nn_z*wall
                  kcontxy3d_tmp(i,j,k) = kcontxy3d_tmp(i,j,k) + int3d(i,j,k)*nn_x*nn_y*wall
                  kcontxz3d_tmp(i,j,k) = kcontxz3d_tmp(i,j,k) + int3d(i,j,k)*nn_x*nn_z*wall
                  kcontyz3d_tmp(i,j,k) = kcontyz3d_tmp(i,j,k) + int3d(i,j,k)*nn_y*nn_z*wall
                  normalization3d_tmp(i,j,k) = normalization3d_tmp(i,j,k) + wall
!
                  alocont_nn3d_tmp(i,j,k,q14) = alocont_nn3d_tmp(i,j,k,q14) + wall*alocont_o_nn3d(i,j,k,q14)
                  alocont_nn3d_tmp(iim1,j,k,q15) = alocont_nn3d_tmp(iim1,j,k,q15) + wall*alocont_o_nn3d(iim1,j,k,q15)
                  alocont_nn3d_tmp(i,jjm1,k,q17) = alocont_nn3d_tmp(i,jjm1,k,q17) + wall*alocont_o_nn3d(i,jjm1,k,q17)
                  alocont_nn3d_tmp(iim1,jjm1,k,q18) = alocont_nn3d_tmp(iim1,jjm1,k,q18) + wall*alocont_o_nn3d(iim1,jjm1,k,q18)
                  alocont_nn3d_tmp(i,j,kkm1,q23) = alocont_nn3d_tmp(i,j,kkm1,q23) + wall*alocont_o_nn3d(i,j,kkm1,q23)
                  alocont_nn3d_tmp(iim1,j,kkm1,q24) = alocont_nn3d_tmp(iim1,j,kkm1,q24) + wall*alocont_o_nn3d(iim1,j,kkm1,q24)
                  alocont_nn3d_tmp(i,jjm1,kkm1,q26) = alocont_nn3d_tmp(i,jjm1,kkm1,q26) + wall*alocont_o_nn3d(i,jjm1,kkm1,q26)
                  alocont_nn3d_tmp(iim1,jjm1,kkm1,q27) = alocont_nn3d_tmp(iim1,jjm1,kkm1,q27) + wall*alocont_o_nn3d(iim1,jjm1,kkm1,q27)
!
               endif
!
             case default
!
            end select
!
!
         enddo
      enddo
   enddo
!
!
!
end subroutine fsc_cont3d_lin
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine coeff2d_contu(opac_im2jm2, opac_im1jm2, opac_ijm2, &
   opac_im2jm1, opac_im1jm1, opac_ijm1, &
   opac_im2j,   opac_im1j,   opac_ij, &
   scont_im2jm2, scont_im1jm2, scont_ijm2, &
   scont_im2jm1, scont_im1jm1, scont_ijm1, &
   scont_im2j,   scont_im1j,   scont_ij, &
   int_im2jm2, int_im1jm2, int_ijm2, &
   int_im2jm1, int_im1jm1, int_ijm1, &
   int_im2j,   int_im1j,   int_ij, &
   x_im2, x_im1, x_i, y_jm2, y_jm1, y_j, x_p, y_p, &
   a_scont, b_scont, c_scont, d_scont, e_scont, &
   f_scont, g_scont, h_scont, i_scont, &
   a_inten, b_inten, c_inten, d_inten, e_inten, &
   f_inten, g_inten, h_inten, i_inten, opac_p, scont_p, int_p)
!
!         interpolates opacity, continuum source function and intensity
!               values given on a 2d grid onto point x_p, y_p
!
!on input (f_* stands for opac_*, scont_* and int_*, respectivly):
!
! y_j      f_im2j----------f_im1j--------------f_ij
!  |          |              |                  |
!  |          |              |                  |
!  |          |              |         x        |
!  |          |              |     (x_p,y_p)    |
!  |          |              |                  |
!y_jm1    f_im2jm1-------f_im1jm1------------f_ijm1
!  |          |              |                  |
!  |          |              |                  |
!  |          |              |                  |
!  |          |              |                  |
!  |          |              |                  |
!y_jm2    f_im2jm2-------f_im1jm2------------f_ijm2
!  |
!  --------x_im2-----------x_im1---------------x_i
!
!        x_p, y_p: coordinates of point onto which shall be interpolated
!
!on output:
!   1. interpolation coefficients for source function and intensity
!      (required for ALO calculations):
!         a_scont, b_scont, c_scont, d_scont, e_scont
!         f_scont, g_scont, h_scont, i_scont
!         a_inten, b_inten, c_inten, d_inten, e_inten
!         f_inten, g_inten, h_inten, i_inten
!
!      such that:
!         f_p = a*f_im2jm2 + b*f_im1jm2 + c*f_ijm2 +
!               d*f_im2jm1 + e*f_im1jm1 + f*f_ijm1 +
!               g*f_im2j   + h*f_im1j   + i*f_ij
!
!   2. interpolated values at point p: opac_p, scont_p, int_p
!
   use prog_type
   use mod_interp2d, only: wp_interp2d
!
   implicit none
!
! ... arguments
   real(dp), intent(in) :: opac_im2jm2, opac_im1jm2, opac_ijm2, &
      opac_im2jm1, opac_im1jm1, opac_ijm1, &
      opac_im2j,   opac_im1j,   opac_ij, &
      scont_im2jm2, scont_im1jm2, scont_ijm2, &
      scont_im2jm1, scont_im1jm1, scont_ijm1, &
      scont_im2j,   scont_im1j,   scont_ij, &
      int_im2jm2, int_im1jm2, int_ijm2, &
      int_im2jm1, int_im1jm1, int_ijm1, &
      int_im2j,   int_im1j,   int_ij, &
      x_im2, x_im1, x_i, y_jm2, y_jm1, y_j, x_p, y_p
   real(dp), intent(out) :: a_scont, b_scont, c_scont, d_scont, e_scont, &
      f_scont, g_scont, h_scont, i_scont, &
      a_inten, b_inten, c_inten, d_inten, e_inten, &
      f_inten, g_inten, h_inten, i_inten, &
      opac_p, int_p, scont_p
!
! ... local scalars
   real(dp) :: dxim1, dxi, dx, tx, ax, bx, cx, a, b, c, &
      dyjm1, dyj, dy, ty, ay, by, cy, &
      axt, bxt, cxt, axt2, bxt2, cxt2, &
      ayt, byt, cyt, ayt2, byt2, cyt2, &
      axj_opac, bxj_opac, cxj_opac, &
      axjm1_opac, bxjm1_opac, cxjm1_opac, &
      axjm2_opac, bxjm2_opac, cxjm2_opac, &
      axtj_opac, bxtj_opac, cxtj_opac, &
      axtjm1_opac, bxtjm1_opac, cxtjm1_opac, &
      axtjm2_opac, bxtjm2_opac, cxtjm2_opac, &
      opacc_jm2, opacc_jm1, opacc_j, &
      opac_jm2, opac_jm1, opac_j, opac_c, &
      ayt_opac, byt_opac, cyt_opac, &
      ay_opac, by_opac, cy_opac, &
      axj_scont, bxj_scont, cxj_scont, &
      axjm1_scont, bxjm1_scont, cxjm1_scont, &
      axjm2_scont, bxjm2_scont, cxjm2_scont, &
      axtj_scont, bxtj_scont, cxtj_scont, &
      axtjm1_scont, bxtjm1_scont, cxtjm1_scont, &
      axtjm2_scont, bxtjm2_scont, cxtjm2_scont, &
      scontc_jm2, scontc_jm1, scontc_j, &
      scont_jm2, scont_jm1, scont_j, scont_c, &
      ayt_scont, byt_scont, cyt_scont, &
      ay_scont, by_scont, cy_scont, &
      axj_int, bxj_int, cxj_int, &
      axjm1_int, bxjm1_int, cxjm1_int, &
      axjm2_int, bxjm2_int, cxjm2_int, &
      axtj_int, bxtj_int, cxtj_int, &
      axtjm1_int, bxtjm1_int, cxtjm1_int, &
      axtjm2_int, bxtjm2_int, cxtjm2_int, &
      intc_jm2, intc_jm1, intc_j, &
      int_jm2, int_jm1, int_j, int_c, &
      ayt_int, byt_int, cyt_int, &
      ay_int, by_int, cy_int, &
      rdxdy
   real(dp) :: d_im2jm2, d_im1jm2, d_ijm2, d_im2jm1, d_im1jm1, d_ijm1, d_im2j, d_im1j, d_ij, &
      w_im2jm2, w_im1jm2, w_ijm2, w_im2jm1, w_im1jm1, w_ijm1, w_im2j, w_im1j, w_ij, &
      norm
   real(dp) :: dx1, dx2, dx3, dx4, dx1s, dx2s, dx3s, dx4s, &
      dy1, dy2, dy3, dy4, dy1s, dy2s, dy3s, dy4s, &
      fac, fac2, norm01, norm10, norm12, norm21, norm11, &
      dxis, dyjs, d1, d2, d3, d4, d5
   real(dp) :: dxp, dxp2, dxp3, dxi2, dxim12, &
      dyp, dyp2, dyp3, dyj2, dyjm12
   real(dp) :: wim1, wi, wjm1, wj, mm, mp
!
!
!define deltax, deltay
   dxim1 = x_im1-x_im2
   dxi = x_i-x_im1
   dx = dxim1+dxi
   dyjm1 = y_jm1-y_jm2
   dyj = y_j-y_jm1
   dy = dyjm1+dyj
!
!define deltax, deltay-ratios
   tx = (x_p-x_im1)/dxi
   ty = (y_p-y_jm1)/dyj
!
!return
!10:  bilinear interpolation
!20:  inverse distance weighting
!30:  bezier inteprolation (with control points from derivatives)
!40:  bezier interpolation (with control points from derivatives, non-monotonic)
!50:  bezier interpolation (with control point from idw)
!60:  cubic interpolation with free central points
!70:  biquadratic interpolation
!80:  bezier interpolation (control point from assigning weight to each point)
!90:  bezier interpolation (control point from assigning weight to each point and approximated curvature)
!100: bicubic interpolation (with zero derivatives at all points)
!110: bezier interpolation (with conrol points from derivatives with minimum weight for left derivative, monotonic)
!
   goto 110
!
!
!-------------------------bilinear interpolation------------------------
!
10 continue
   rdxdy=tx*ty
!
   a_scont=0.d0
   b_scont=0.d0
   c_scont=0.d0
   d_scont=0.d0
   e_scont=1.d0-tx-ty+rdxdy
   f_scont=tx-rdxdy
   g_scont=0.d0
   h_scont=ty-rdxdy
   i_scont=rdxdy
!
   a_inten=a_scont
   b_inten=b_scont
   c_inten=c_scont
   d_inten=d_scont
   e_inten=e_scont
   f_inten=f_scont
   g_inten=g_scont
   h_inten=h_scont
   i_inten=i_scont
!
!opac_p = e_scont*opac_im1jm1 + f_scont*opac_ijm1 + h_scont*opac_im1j + i_scont*opac_ij
   scont_p = e_scont*scont_im1jm1 + f_scont*scont_ijm1 + h_scont*scont_im1j + i_scont*scont_ij
   int_p = e_scont*int_im1jm1 + f_scont*int_ijm1 + h_scont*int_im1j + i_scont*int_ij
   return
!
!-------------------------inverse distance weighting--------------------
!
20 continue
!
!define weighting factor
   fac=2.d0
!
!calculate distance to each point
   d_im2jm2 = (x_p-x_im2)**2+(y_p-y_jm2)**2
   d_im1jm2 = (x_p-x_im1)**2+(y_p-y_jm2)**2
   d_ijm2 = (x_p-x_i)**2+(y_p-y_jm2)**2
   d_im2jm1 = (x_p-x_im2)**2+(y_p-y_jm1)**2
   d_im1jm1 = (x_p-x_im1)**2+(y_p-y_jm1)**2
   d_ijm1 = (x_p-x_i)**2+(y_p-y_jm1)**2
   d_im2j = (x_p-x_im2)**2+(y_p-y_j)**2
   d_im1j = (x_p-x_im1)**2+(y_p-y_j)**2
   d_ij = (x_p-x_i)**2+(y_p-y_j)**2
!
!avoid division by zero if d_ij=0, or d_im1j=0, ...
   if(d_im2jm2.eq.0.) then
      w_im2jm2=1.d0
      w_im1jm2=0.d0
      w_ijm2=0.d0
      w_im2jm1=0.d0
      w_im1jm1=0.d0
      w_ijm1=0.d0
      w_im2j=0.d0
      w_im1j=0.d0
      w_ij=0.d0
   elseif(d_im1jm2.eq.0.) then
      w_im2jm2=0.d0
      w_im1jm2=1.d0
      w_ijm2=0.d0
      w_im2jm1=0.d0
      w_im1jm1=0.d0
      w_ijm1=0.d0
      w_im2j=0.d0
      w_im1j=0.d0
      w_ij=0.d0
   elseif(d_ijm2.eq.0.) then
      w_im2jm2=0.d0
      w_im1jm2=0.d0
      w_ijm2=1.d0
      w_im2jm1=0.d0
      w_im1jm1=0.d0
      w_ijm1=0.d0
      w_im2j=0.d0
      w_im1j=0.d0
      w_ij=0.d0
   elseif(d_im2jm1.eq.0.) then
      w_im2jm2=0.d0
      w_im1jm2=0.d0
      w_ijm2=0.d0
      w_im2jm1=1.d0
      w_im1jm1=0.d0
      w_ijm1=0.d0
      w_im2j=0.d0
      w_im1j=0.d0
      w_ij=0.d0
   elseif(d_im1jm1.eq.0.) then
      w_im2jm2=0.d0
      w_im1jm2=0.d0
      w_ijm2=0.d0
      w_im2jm1=0.d0
      w_im1jm1=1.d0
      w_ijm1=0.d0
      w_im2j=0.d0
      w_im1j=0.d0
      w_ij=0.d0
   elseif(d_ijm1.eq.0.) then
      w_im2jm2=0.d0
      w_im1jm2=0.d0
      w_ijm2=0.d0
      w_im2jm1=0.d0
      w_im1jm1=0.d0
      w_ijm1=1.d0
      w_im2j=0.d0
      w_im1j=0.d0
      w_ij=0.d0
   elseif(d_im2j.eq.0.) then
      w_im2jm2=0.d0
      w_im1jm2=0.d0
      w_ijm2=0.d0
      w_im2jm1=0.d0
      w_im1jm1=0.d0
      w_ijm1=0.d0
      w_im2j=1.d0
      w_im1j=0.d0
      w_ij=0.d0
   elseif(d_im1j.eq.0.) then
      w_im2jm2=0.d0
      w_im1jm2=0.d0
      w_ijm2=0.d0
      w_im2jm1=0.d0
      w_im1jm1=0.d0
      w_ijm1=0.d0
      w_im2j=0.d0
      w_im1j=1.d0
      w_ij=0.d0
   elseif(d_ij.eq.0.) then
      w_im2jm2=0.d0
      w_im1jm2=0.d0
      w_ijm2=0.d0
      w_im2jm1=0.d0
      w_im1jm1=0.d0
      w_ijm1=0.d0
      w_im2j=0.d0
      w_im1j=0.d0
      w_ij=1.d0
   else
!
      w_im2jm2=1.d0/d_im2jm2**(fac/2.d0)
      w_im1jm2=1.d0/d_im1jm2**(fac/2.d0)
      w_ijm2=1.d0/d_ijm2**(fac/2.d0)
      w_im2jm1=1.d0/d_im2jm1**(fac/2.d0)
      w_im1jm1=1.d0/d_im1jm1**(fac/2.d0)
      w_ijm1=1.d0/d_ijm1**(fac/2.d0)
      w_im2j=1.d0/d_im2j**(fac/2.d0)
      w_im1j=1.d0/d_im1j**(fac/2.d0)
      w_ij=1.d0/d_ij**(fac/2.d0)
!
   endif
!
   norm = w_im2jm2 + w_im1jm2 + w_ijm2 + w_im2jm1 + w_im1jm1 + w_ijm1 + w_im2j + w_im1j + w_ij
!
   a_scont = w_im2jm2/norm
   b_scont = w_im1jm2/norm
   c_scont = w_ijm2/norm
   d_scont = w_im2jm1/norm
   e_scont = w_im1jm1/norm
   f_scont = w_ijm1/norm
   g_scont = w_im2j/norm
   h_scont = w_im1j/norm
   i_scont = w_ij/norm
!
   a_inten = a_scont
   b_inten = b_scont
   c_inten = c_scont
   d_inten = d_scont
   e_inten = e_scont
   f_inten = f_scont
   g_inten = g_scont
   h_inten = h_scont
   i_inten = i_scont
!
!opac_p = a_scont*opac_im2jm2 + b_scont*opac_im1jm2 + c_scont*opac_ijm2 + &
!         d_scont*opac_im2jm1 + e_scont*opac_im1jm1 + f_scont*opac_ijm1 + &
!         g_scont*opac_im2j   + h_scont*opac_im1j   + i_scont*opac_ij
   scont_p = a_scont*scont_im2jm2 + b_scont*scont_im1jm2 + c_scont*scont_ijm2 + &
      d_scont*scont_im2jm1 + e_scont*scont_im1jm1 + f_scont*scont_ijm1 + &
      g_scont*scont_im2j   + h_scont*scont_im1j   + i_scont*scont_ij
   int_p = a_inten*int_im2jm2 + b_inten*int_im1jm2 + c_inten*int_ijm2 + &
      d_inten*int_im2jm1 + e_inten*int_im1jm1 + f_inten*int_ijm1 + &
      g_inten*int_im2j   + h_inten*int_im1j   + i_inten*int_ij

   return
!
!---------------monotonic quadratic bezier interpolation----------------
!--------(with derivatives for control point from weighted mean)--------
!
30 continue
!
   ax = (1.d0-tx)**2
   bx = 2.d0*tx*(1.d0-tx)
   cx = tx**2
!
!calculate control points on each j-level:
!   scontc_jm2, scontc_jm1, scontc_j
!   intc_jm2, intc_jm1, intc_j
!   opacc_jm2, opacc_jm1, opacc_j
   axt = -dxi**2/2.d0/dxim1/dx
   bxt = dx/2.d0/dxim1
   cxt = dxim1/2.d0/dx
!
!opacc_jm2 = opac_im2jm2*axt + opac_im1jm2*bxt + opac_ijm2*cxt
!opacc_jm1 = opac_im2jm1*axt + opac_im1jm1*bxt + opac_ijm1*cxt
!opacc_j   = opac_im2j*axt   + opac_im1j*bxt   + opac_ij*cxt
!
   scontc_jm2 = scont_im2jm2*axt + scont_im1jm2*bxt + scont_ijm2*cxt
   scontc_jm1 = scont_im2jm1*axt + scont_im1jm1*bxt + scont_ijm1*cxt
   scontc_j   = scont_im2j*axt   + scont_im1j*bxt   + scont_ij*cxt
!
   intc_jm2 = int_im2jm2*axt + int_im1jm2*bxt + int_ijm2*cxt
   intc_jm1 = int_im2jm1*axt + int_im1jm1*bxt + int_ijm1*cxt
   intc_j   = int_im2j*axt   + int_im1j*bxt   + int_ij*cxt
!
!ensure monotonicity on level j
!call coeffcr1d_mbez(opacc_j, opac_im1j, opac_ij, axt, bxt, cxt, axtj_opac, bxtj_opac, cxtj_opac)
   call coeffcr1d_mbez(scontc_j, scont_im1j, scont_ij, axt, bxt, cxt, axtj_scont, bxtj_scont, cxtj_scont)
   call coeffcr1d_mbez(intc_j, int_im1j, int_ij, axt, bxt, cxt, axtj_int, bxtj_int, cxtj_int)
!axj_opac = axtj_opac*bx
!bxj_opac = bxtj_opac*bx + ax
!cxj_opac = cxtj_opac*bx + cx
!opac_j = axj_opac*opac_im2j + bxj_opac*opac_im1j + cxj_opac*opac_ij
   axj_scont = axtj_scont*bx
   bxj_scont = bxtj_scont*bx + ax
   cxj_scont = cxtj_scont*bx + cx
   scont_j = axj_scont*scont_im2j + bxj_scont*scont_im1j + cxj_scont*scont_ij
   axj_int = axtj_int*bx
   bxj_int = bxtj_int*bx + ax
   cxj_int = cxtj_int*bx + cx
   int_j = axj_int*int_im2j + bxj_int*int_im1j + cxj_int*int_ij
!
!ensure monotonicity on level j-1
!call coeffcr1d_mbez(opacc_jm1, opac_im1jm1, opac_ijm1, axt, bxt, cxt, axtjm1_opac, bxtjm1_opac, cxtjm1_opac)
   call coeffcr1d_mbez(scontc_jm1, scont_im1jm1, scont_ijm1, axt, bxt, cxt, axtjm1_scont, bxtjm1_scont, cxtjm1_scont)
   call coeffcr1d_mbez(intc_jm1, int_im1jm1, int_ijm1, axt, bxt, cxt, axtjm1_int, bxtjm1_int, cxtjm1_int)
!axjm1_opac = axtjm1_opac*bx
!bxjm1_opac = bxtjm1_opac*bx + ax
!cxjm1_opac = cxtjm1_opac*bx + cx
!opac_jm1 = axjm1_opac*opac_im2jm1 + bxjm1_opac*opac_im1jm1 + cxjm1_opac*opac_ijm1
   axjm1_scont = axtjm1_scont*bx
   bxjm1_scont = bxtjm1_scont*bx + ax
   cxjm1_scont = cxtjm1_scont*bx + cx
   scont_jm1 = axjm1_scont*scont_im2jm1 + bxjm1_scont*scont_im1jm1 + cxjm1_scont*scont_ijm1
   axjm1_int = axtjm1_int*bx
   bxjm1_int = bxtjm1_int*bx + ax
   cxjm1_int = cxtjm1_int*bx + cx
   int_jm1 = axjm1_int*int_im2jm1 + bxjm1_int*int_im1jm1 + cxjm1_int*int_ijm1
!
!ensure monotonicity on level j-2
!call coeffcr1d_mbez(opacc_jm2, opac_im1jm2, opac_ijm2, axt, bxt, cxt, axtjm2_opac, bxtjm2_opac, cxtjm2_opac)
   call coeffcr1d_mbez(scontc_jm2, scont_im1jm2, scont_ijm2, axt, bxt, cxt, axtjm2_scont, bxtjm2_scont, cxtjm2_scont)
   call coeffcr1d_mbez(intc_jm2, int_im1jm2, int_ijm2, axt, bxt, cxt, axtjm2_int, bxtjm2_int, cxtjm2_int)
!axjm2_opac = axtjm2_opac*bx
!bxjm2_opac = bxtjm2_opac*bx + ax
!cxjm2_opac = cxtjm2_opac*bx + cx
!opac_jm2 = axjm2_opac*opac_im2jm2 + bxjm2_opac*opac_im1jm2 + cxjm2_opac*opac_ijm2
   axjm2_scont = axtjm2_scont*bx
   bxjm2_scont = bxtjm2_scont*bx + ax
   cxjm2_scont = cxtjm2_scont*bx + cx
   scont_jm2 = axjm2_scont*scont_im2jm2 + bxjm2_scont*scont_im1jm2 + cxjm2_scont*scont_ijm2
   axjm2_int = axtjm2_int*bx
   bxjm2_int = bxtjm2_int*bx + ax
   cxjm2_int = cxtjm2_int*bx + cx
   int_jm2 = axjm2_int*int_im2jm2 + bxjm2_int*int_im1jm2 + cxjm2_int*int_ijm2
!
!
!
!calculate control point for interpolation along y
   ayt = -dyj**2/2.d0/dyjm1/dy
   byt = dy/2.d0/dyjm1
   cyt = dyjm1/2.d0/dy
!opac_c = opac_jm2*ayt + opac_jm1*byt + opac_j*cyt
   scont_c = scont_jm2*ayt + scont_jm1*byt + scont_j*cyt
   int_c = int_jm2*ayt + int_jm1*byt + int_j*cyt
!
!ensure monotonicity
!call coeffcr1d_mbez(opac_c, opac_jm1, opac_j, ayt, byt, cyt, ayt_opac, byt_opac, cyt_opac)
   call coeffcr1d_mbez(scont_c, scont_jm1, scont_j, ayt, byt, cyt, ayt_scont, byt_scont, cyt_scont)
   call coeffcr1d_mbez(int_c, int_jm1, int_j, ayt, byt, cyt, ayt_int, byt_int, cyt_int)
!
   ay = (1.d0-ty)**2
   by = 2.d0*ty*(1.d0-ty)
   cy = ty**2
!
!ay_opac = ayt_opac*by
!by_opac = byt_opac*by + ay
!cy_opac = cyt_opac*by + cy
!opac_p = ay_opac*opac_jm2 + by_opac*opac_jm1 + cy_opac*opac_j
   ay_scont = ayt_scont*by
   by_scont = byt_scont*by + ay
   cy_scont = cyt_scont*by + cy
   scont_p = ay_scont*scont_jm2 + by_scont*scont_jm1 + cy_scont*scont_j
   ay_int = ayt_int*by
   by_int = byt_int*by + ay
   cy_int = cyt_int*by + cy
   int_p = ay_int*int_jm2 + by_int*int_jm1 + cy_int*int_j
!
   a_scont = ay_scont*axjm2_scont
   b_scont = ay_scont*bxjm2_scont
   c_scont = ay_scont*cxjm2_scont
   d_scont = by_scont*axjm1_scont
   e_scont = by_scont*bxjm1_scont
   f_scont = by_scont*cxjm1_scont
   g_scont = cy_scont*axj_scont
   h_scont = cy_scont*bxj_scont
   i_scont = cy_scont*cxj_scont
!
   a_inten = ay_int*axjm2_int
   b_inten = ay_int*bxjm2_int
   c_inten = ay_int*cxjm2_int
   d_inten = by_int*axjm1_int
   e_inten = by_int*bxjm1_int
   f_inten = by_int*cxjm1_int
   g_inten = cy_int*axj_int
   h_inten = cy_int*bxj_int
   i_inten = cy_int*cxj_int
!
!write(*,'(9es20.8)') a_scont, b_scont, c_scont, d_scont, e_scont, f_scont, g_scont, h_scont, i_scont
!
   return
!
!-------------quadratic bezier interpolation (non-monotonic)------------
!--------(with derivatives for control point from weighted mean)--------
!
40 continue
!
   ax = (1.d0-tx)**2
   bx = 2.d0*tx*(1.d0-tx)
   cx = tx**2
!
!calculate control points on each j-level:
!   scontc_jm2, scontc_jm1, scontc_j
!   intc_jm2, intc_jm1, intc_j
!   opacc_jm2, opacc_jm1, opacc_j
   axt = -dxi**2/2.d0/dxim1/dx
   bxt = dx/2.d0/dxim1
   cxt = dxim1/2.d0/dx
!
   axj_scont = axt*bx
   bxj_scont = ax+bx*bxt
   cxj_scont = cx+bx*cxt
!
   ay = (1.d0-ty)**2
   by = 2.d0*ty*(1.d0-ty)
   cy = ty**2
!
   ayt = -dyj**2/2.d0/dyjm1/dy
   byt = dy/2.d0/dyjm1
   cyt = dyjm1/2.d0/dy
!
   ay_scont = ayt*by
   by_scont = ay+by*byt
   cy_scont = cy+by*cyt
!
   a_scont = ay_scont*axj_scont
   b_scont = ay_scont*bxj_scont
   c_scont = ay_scont*cxj_scont
   d_scont = by_scont*axj_scont
   e_scont = by_scont*bxj_scont
   f_scont = by_scont*cxj_scont
   g_scont = cy_scont*axj_scont
   h_scont = cy_scont*bxj_scont
   i_scont = cy_scont*cxj_scont
!
   a_inten = a_scont
   b_inten = b_scont
   c_inten = c_scont
   d_inten = d_scont
   e_inten = e_scont
   f_inten = f_scont
   g_inten = g_scont
   h_inten = h_scont
   i_inten = i_scont
!
!
!opac_p = a_scont*opac_im2jm2 + b_scont*opac_im1jm2 + c_scont*opac_ijm2 + &
!         d_scont*opac_im2jm1 + e_scont*opac_im1jm1 + f_scont*opac_ijm1 + &
!         g_scont*opac_im2j   + h_scont*opac_im1j   + i_scont*opac_ij
   scont_p = a_scont*scont_im2jm2 + b_scont*scont_im1jm2 + c_scont*scont_ijm2 + &
      d_scont*scont_im2jm1 + e_scont*scont_im1jm1 + f_scont*scont_ijm1 + &
      g_scont*scont_im2j   + h_scont*scont_im1j   + i_scont*scont_ij
   int_p = a_inten*int_im2jm2 + b_inten*int_im1jm2 + c_inten*int_ijm2 + &
      d_inten*int_im2jm1 + e_inten*int_im1jm1 + f_inten*int_ijm1 + &
      g_inten*int_im2j   + h_inten*int_im1j   + i_inten*int_ij
!
!write(*,'(9es20.8)') a_scont, b_scont, c_scont, d_scont, e_scont, f_scont, g_scont, h_scont, i_scont
!
   return
!
!---------------monotonic quadratic bezier interpolation----------------
!----------(with control points from inverse distance weighing)---------
!
50 continue
!
   dxi=(x_i-x_im1)
   dyj=(y_j-y_jm1)
!
   tx = (x_p-x_im1)/dxi
   ty = (y_p-y_jm1)/dyj
!
   ax=(1.d0-tx)**2
   bx=2.d0*tx*(1.d0-tx)
   cx=tx**2
!
   ay=(1.d0-ty)**2
   by=2.d0*ty*(1.d0-ty)
   cy=ty**2
!
!calculate all control points via inverse distance weighting (weight = distance**fac)
   fac = -4.d0
!
!using only 4 local points
   dxis=dxi**2
   dyjs=dyj**2
   d1=(dxis/4.d0)**fac
   d2=(dyjs/4.d0)**fac
   d3=(dxis+dyjs/4.d0)**fac
   d4=(dxis/4.d0+dyjs)**fac
   d5=((dxis+dyjs)/4.d0)**fac
   norm01 = 2.d0*(d2+d3)
   norm10 = 2.d0*(d1+d4)
   norm12 = norm10
   norm21 = norm01
   norm11 = 4.d0*d5

   a_scont = 0.d0
   b_scont = 0.d0
   c_scont = 0.d0
   d_scont = 0.d0
   e_scont = ax*ay + ax*by*d2/norm01 + bx*ay*d1/norm10 + bx*by*d5/norm11 + bx*cy*d4/norm12 + cx*by*d3/norm21
   f_scont = cx*ay + ax*by*d3/norm01 + bx*ay*d1/norm10 + bx*by*d5/norm11 + bx*cy*d4/norm12 + cx*by*d2/norm21
   g_scont = 0.d0
   h_scont = ax*cy + ax*by*d2/norm01 + bx*ay*d4/norm10 + bx*by*d5/norm11 + bx*cy*d1/norm12 + cx*by*d3/norm21
   i_scont = cx*cy + ax*by*d3/norm01 + bx*ay*d4/norm10 + bx*by*d5/norm11 + bx*cy*d1/norm12 + cx*by*d2/norm21
!
!
!
!using all 9 points
!fac=-4.d0
!dx1 = x_i-x_im1
!dx2 = x_i-x_im2
!dx3 = (x_i+x_im1-2.d0*x_im2)/2.d0
!dx4 = x_im1-x_im2
!dx1s = dx1**2
!dx2s = dx2**2
!dx3s = dx3**2
!dx4s = dx4**2
!dy1 = y_j-y_jm1
!dy2 = y_j-y_jm2
!dy3 = (y_j+y_jm1-2.d0*y_jm2)/2.d0
!dy4 = y_jm1-y_jm2
!dy1s = dy1**2
!dy2s = dy2**2
!dy3s = dy3**2
!dy4s = dy4**2
!norm01 = (dx4s+dy3s)**fac + dy3s**fac + (dx1s+dy3s)**fac + (dx4s+dy1s/4.)**fac + &
!         (dy1s/4.)**fac + (dx1s+dy1s/4.)**fac + (dx4s+dy1s/4.)**fac + (dy1s/4.)**fac + (dx1s+dy1s/4.)**fac
!norm10 = (dx3s+dy4s)**fac + (dx1s/4.+dy4s)**fac + (dx1s/4.+dy4s)**fac + dx3s**fac + (dx1s/4.)**fac + &
!         (dx1s/4.)**fac + (dx3s+dy1s)**fac + (dx1s/4.+dy1s)**fac + (dx1s/4.+dy1s)**fac
!norm12 = (dx3s+dy2s)**fac + (dx1s/4.+dy2s)**fac + (dx1s/4.+dy2s)**fac + (dx3s+dy1s)**fac + &
!         (dx1s/4.+dy1s)**fac + (dx1s/4.+dy1s)**fac + dx3s**fac + (dx1s/4.)**fac + (dx1s/4.)**fac
!norm21 = (dx2s+dx3s)**fac + (dx1s+dy3s)**fac + dy3s**fac + (dx2s+dy1s/4.)**fac + (dx1s+dy1s/4.)**fac + &
!         (dy1s/4.)**fac + (dx2s+dy1s/4.)**fac + (dx1s+dy1s/4.)**fac + (dy1s/4.)**fac
!norm11 = (dx3s+dy3s)**fac + 2.*(dx1s/4.+dy3s)**fac + 2.*(dx3s+dy1s/4.)**fac + 4.*(dx1s/4.+dy1s/4.)**fac
!!
!a_scont = ax*by/norm01*(dx4s+dy3s)**fac + bx*ay/norm10*(dx3s+dy4s)**fac + bx*by/norm11*(dx3s+dy3s)**fac + &
!          bx*cy/norm12*(dx3s+dy2s)**fac + cx*by/norm21*(dx2s+dy3s)**fac
!b_scont = ax*by/norm01*dy3s**fac + bx*ay/norm10*(dx1s/4.+dy4s)**fac + bx*by/norm11*(dx1s/4.+dy3s)**fac + &
!          bx*cy/norm12*(dx1s/4.+dy2s)**fac + cx*by/norm21*(dx1s+dy3s)**fac
!c_scont = ax*by/norm01*(dx1s+dy3s)**fac + bx*ay/norm10*(dx1s/4.+dy4s)**fac + bx*by/norm11*(dx1s/4.+dy3s)**fac + &
!          bx*cy/norm12*(dx1s/4.+dy2s)**fac + cx*by/norm21*dy3s**fac
!d_scont = ax*by/norm01*(dx4s+dy1s/4.)**fac + bx*ay/norm10*dx3s**fac + bx*by/norm11*(dx3s+dy1s/4.)**fac + &
!          bx*cy/norm12*(dx3s+dy1s)**fac + cx*by/norm21*(dx2s+dy1s/4.)**fac
!e_scont = ax*ay + ax*by/norm01*(dy1s/4.)**fac + bx*ay/norm10*(dx1s/4.)**fac + bx*by/norm11*(dx1s/4.+dy1s/4.)**fac + &
!          bx*cy/norm12*(dx1s/4.+dy1s)**fac + cx*by/norm21*(dx1s+dy1s/4.)**fac
!f_scont = cx*ay + ax*by/norm01*(dx1s+dy1s/4.)**fac +bx*ay/norm10*(dx1s/4.)**fac + bx*by/norm11*(dx1s/4.+dy1s/4.)**fac + &
!          bx*cy/norm12*(dx1s/4.+dy1s)**fac + cx*by/norm21*(dy1s/4.)**fac
!g_scont = ax*by/norm01*(dx4s+dy1s/4.)**fac + bx*ay/norm10*(dx3s+dy1s)**fac + bx*by/norm11*(dx3s+dy1s/4.)**fac + &
!          bx*cy/norm12*dx3s**fac + cx*by/norm21*(dx2s+dy1s/4.)**fac
!h_scont = ax*cy + ax*by/norm01*(dy1s/4.)**fac + bx*ay/norm10*(dx1s/4.+dy1s)**fac + bx*by/norm11*(dx1s/4.+dy1s/4.)**fac + &
!          bx*cy/norm12*(dx1s/4.)**fac + cx*by/norm21*(dx1s+dy1s/4.)**fac
!i_scont = cx*cy + ax*by/norm01*(dx1s+dy1s/4.)**fac + bx*ay/norm10*(dx1s/4.+dy1s)**fac + bx*by/norm11*(dx1s/4.+dy1s/4.)**fac + &
!          bx*cy/norm12*(dx1s/4.)**fac + cx*by/norm21*(dy1s/4.)**fac
!
   a_inten = a_scont
   b_inten = b_scont
   c_inten = c_scont
   d_inten = d_scont
   e_inten = e_scont
   f_inten = f_scont
   g_inten = g_scont
   h_inten = h_scont
   i_inten = i_scont
!
!opac_p = a_scont*opac_im2jm2 + b_scont*opac_im1jm2 + c_scont*opac_ijm2 + &
!         d_scont*opac_im2jm1 + e_scont*opac_im1jm1 + f_scont*opac_ijm1 + &
!         g_scont*opac_im2j   + h_scont*opac_im1j   + i_scont*opac_ij
   scont_p = a_scont*scont_im2jm2 + b_scont*scont_im1jm2 + c_scont*scont_ijm2 + &
      d_scont*scont_im2jm1 + e_scont*scont_im1jm1 + f_scont*scont_ijm1 + &
      g_scont*scont_im2j   + h_scont*scont_im1j   + i_scont*scont_ij
   int_p = a_inten*int_im2jm2 + b_inten*int_im1jm2 + c_inten*int_ijm2 + &
      d_inten*int_im2jm1 + e_inten*int_im1jm1 + f_inten*int_ijm1 + &
      g_inten*int_im2j   + h_inten*int_im1j   + i_inten*int_ij

!write(*,'(9es20.8)') a_scont, b_scont, c_scont, d_scont, e_scont, f_scont, g_scont, h_scont, i_scont
!if(a_scont.lt.0.) stop
!if(b_scont.lt.0.) stop
!if(c_scont.lt.0.) stop
!if(d_scont.lt.0.) stop
!if(e_scont.lt.0.) stop
!if(f_scont.lt.0.) stop
!if(g_scont.lt.0.) stop
!if(h_scont.lt.0.) stop
!if(i_scont.lt.0.) stop
!write(*,*)

   return
!
!---------------cubic interpolation with free central points------------
!-----------------warning: might cause negative intensities-------------
!
60 continue
!
   dxp = x_p-x_im1
   dxp2 = dxp**2
   dxp3 = dxp**3
   dxi = x_i-x_im1
   dxi2 = dxi**2
   dxim1 = x_im1-x_im2
   dxim12 = dxim1**2
   dx = dxi+dxim1
   dx2 = dx**2
   dx3 = dx**3
!
   dyp = y_p-y_jm1
   dyp2 = dyp**2
   dyp3 = dyp**3
   dyj = y_j-y_jm1
   dyj2 = dyj**2
   dyjm1 = y_jm1-y_jm2
   dyjm12 = dyjm1**2
   dy = dyj+dyjm1
   dy2 = dy**2
   dy3 = dy**3

   ax = ((dxim1-dxi)*dxp3 + 2.d0*(dxi2+dxim12-dxi*dxim1)*dxp2+dxi*(dxi*dxim1-4.d0*dxim12-dxi2)*dxp+2.d0*dxi2*dxim12)/dxim1/dx3
   bx = ((dxi-dxim1)*dxp3 + 2.d0*(dxi*dxim1-dxi2-dxim12)*dxp2)/dxi/dxim1/dx2 + &
      (dxi2-dxim12)*(dxi2+dxim12-dxi*dxim1)*dxp/dxi/dxim1/dx3 + (dxi2+dxim12)/dx2
   cx = ((dxim1-dxi)*dxp3 + 2.d0*(dxi2+dxim12-dxi*dxim1)*dxp2 + dxim1*(4.d0*dxi2+dxim12-dxi*dxim1)*dxp + 2.d0*dxim12*dxi2)/dxi/dx3

   ay = ((dyjm1-dyj)*dyp3 + 2.d0*(dyj2+dyjm12-dyj*dyjm1)*dyp2+dyj*(dyj*dyjm1-4.d0*dyjm12-dyj2)*dyp+2.d0*dyj2*dyjm12)/dyjm1/dy3
   by = ((dyj-dyjm1)*dyp3 + 2.d0*(dyj*dyjm1-dyj2-dyjm12)*dyp2)/dyj/dyjm1/dy2 + &
      (dyj2-dyjm12)*(dyj2+dyjm12-dyj*dyjm1)*dyp/dyj/dyjm1/dy3 + (dyj2+dyjm12)/dy2
   cy = ((dyjm1-dyj)*dyp3 + 2.d0*(dyj2+dyjm12-dyj*dyjm1)*dyp2 + dyjm1*(4.d0*dyj2+dyjm12-dyj*dyjm1)*dyp + 2.d0*dyjm12*dyj2)/dyj/dy3

   a_scont = ay*ax
   b_scont = ay*bx
   c_scont = ay*cx
   d_scont = by*ax
   e_scont = by*bx
   f_scont = by*cx
   g_scont = cy*ax
   h_scont = cy*bx
   i_scont = cy*cx

   a_inten = a_scont
   b_inten = b_scont
   c_inten = c_scont
   d_inten = d_scont
   e_inten = e_scont
   f_inten = f_scont
   g_inten = g_scont
   h_inten = h_scont
   i_inten = i_scont
!
!opac_p = a_scont*opac_im2jm2 + b_scont*opac_im1jm2 + c_scont*opac_ijm2 + &
!         d_scont*opac_im2jm1 + e_scont*opac_im1jm1 + f_scont*opac_ijm1 + &
!         g_scont*opac_im2j   + h_scont*opac_im1j   + i_scont*opac_ij
   scont_p = a_scont*scont_im2jm2 + b_scont*scont_im1jm2 + c_scont*scont_ijm2 + &
      d_scont*scont_im2jm1 + e_scont*scont_im1jm1 + f_scont*scont_ijm1 + &
      g_scont*scont_im2j   + h_scont*scont_im1j   + i_scont*scont_ij
   int_p = a_inten*int_im2jm2 + b_inten*int_im1jm2 + c_inten*int_ijm2 + &
      d_inten*int_im2jm1 + e_inten*int_im1jm1 + f_inten*int_ijm1 + &
      g_inten*int_im2j   + h_inten*int_im1j   + i_inten*int_ij

   if(int_p.lt.0.) stop 'error in coeff2d_contu: negative intensities'

   return
!
!-----------------------biquadratic interpolation-----------------------
!-----------------warning: might cause negative intensities-------------
!
70 continue
!
   dxim1 = x_im1-x_im2
   dxi = x_i-x_im1
   dx = dxim1+dxi
!
   dyjm1 = y_jm1-y_jm2
   dyj = y_j-y_jm1
   dy = dyjm1+dyj
!
   tx = x_p-x_im1
   ty = y_p-y_jm1
!
   a = -dxi/dxim1/dx
   b = (dxi**2-dxim1**2)/dxi/dxim1/dx
   c = dxim1/dxi/dx
   ax = tx**2/dxim1/dx + a*tx
   bx = 1.d0-tx**2/dxi/dxim1 + b*tx
   cx = tx**2/dxi/dx + c*tx
!
   a = -dyj/dyjm1/dy
   b = (dyj**2-dyjm1**2)/dyj/dyjm1/dy
   c = dyjm1/dyj/dy
   ay = ty**2/dyjm1/dy + a*ty
   by = 1.d0-ty**2/dyj/dyjm1 + b*ty
   cy = ty**2/dyj/dy + c*ty
!
   a_scont=ax*ay
   b_scont=bx*ay
   c_scont=cx*ay
   d_scont=ax*by
   e_scont=bx*by
   f_scont=cx*by
   g_scont=ax*cy
   h_scont=bx*cy
   i_scont=cx*cy
!
   a_inten=a_scont
   b_inten=b_scont
   c_inten=c_scont
   d_inten=d_scont
   e_inten=e_scont
   f_inten=f_scont
   g_inten=g_scont
   h_inten=h_scont
   i_inten=i_scont

!a= a_scont+b_scont+c_scont+d_scont+e_scont+f_scont+g_scont+h_scont+i_scont
!write(*,*) a
!if(abs(a-1.d0).gt.1.d-14) stop 'error'
!
!opac_p = a_scont*opac_im2jm2 + b_scont*opac_im1jm2 + c_scont*opac_ijm2 + &
!         d_scont*opac_im2jm1 + e_scont*opac_im1jm1 + f_scont*opac_ijm1 + &
!         g_scont*opac_im2j   + h_scont*opac_im1j   + i_scont*opac_ij
   scont_p = a_scont*scont_im2jm2 + b_scont*scont_im1jm2 + c_scont*scont_ijm2 + &
      d_scont*scont_im2jm1 + e_scont*scont_im1jm1 + f_scont*scont_ijm1 + &
      g_scont*scont_im2j   + h_scont*scont_im1j   + i_scont*scont_ij
   int_p = a_inten*int_im2jm2 + b_inten*int_im1jm2 + c_inten*int_ijm2 + &
      d_inten*int_im2jm1 + e_inten*int_im1jm1 + f_inten*int_ijm1 + &
      g_inten*int_im2j   + h_inten*int_im1j   + i_inten*int_ij
   if(int_p.lt.0.) stop 'error in coeff2d_contu: negative intensities'
!
   return
!
!-------------------4-point quadratic bezier interpolation--------------
!-----------(with predefined control point by assigning weights)--------
!
80 continue
!
   wim1=0.7d0
   wi=1.d0-wim1
   wjm1=0.7d0
   wj=1.d0-wjm1
!
   tx = (x_p-x_im1)/(x_i-x_im1)
   ty = (y_p-y_jm1)/(y_j-y_jm1)
!
   ax = (1.d0-tx)**2 + 2.d0*tx*(1.d0-tx)*wim1
   bx = tx**2 + 2.d0*tx*(1.d0-tx)*wi
!
   ay = (1.d0-ty)**2 + 2.d0*ty*(1.d0-ty)*wjm1
   by = ty**2 + 2.d0*ty*(1.d0-ty)*wj
!
   a_scont=0.d0
   b_scont=0.d0
   c_scont=0.d0
   d_scont=0.d0
   e_scont=ax*ay
   f_scont=bx*ay
   g_scont=0.d0
   h_scont=ax*by
   i_scont=bx*by
!
   a_inten=a_scont
   b_inten=b_scont
   c_inten=c_scont
   d_inten=d_scont
   e_inten=e_scont
   f_inten=f_scont
   g_inten=g_scont
   h_inten=h_scont
   i_inten=i_scont
!
!opac_p = a_scont*opac_im2jm2 + b_scont*opac_im1jm2 + c_scont*opac_ijm2 + &
!         d_scont*opac_im2jm1 + e_scont*opac_im1jm1 + f_scont*opac_ijm1 + &
!         g_scont*opac_im2j   + h_scont*opac_im1j   + i_scont*opac_ij
   scont_p = e_scont*scont_im1jm1 + f_scont*scont_ijm1 + h_scont*scont_im1j + i_scont*scont_ij
   int_p = e_inten*int_im1jm1 + f_inten*int_ijm1 + h_inten*int_im1j + i_inten*int_ij
!
   return
!
!-------------------4-point quadratic bezier interpolation--------------
!------------with predefined control point by assigning weights---------
!------------and approximated curvature from negihboring points---------
!
90 continue
!
!wim1=fac, wi=1.d0-fac on each j level, or other way round
!(depending on curvature)
   fac=0.6d0
!
   dxim1=x_im1-x_im2
   dxi=x_i-x_im1
!
   dyjm1=y_jm1-y_jm2
   dyj=y_j-y_jm1
!
   tx = (x_p-x_im1)/dxi
   ty = (y_p-y_jm1)/dyj
!
!level j
   mm=abs((scont_im1j-scont_im2j)/dxim1)
   mp=abs((scont_ij-scont_im1j)/dxi)
   if(mm.gt.mp) then
      wim1=1.d0-fac
      wi=fac
   elseif(mm.lt.mp) then
      wim1=fac
      wi=1.d0-fac
   else
      wim1=0.5d0
      wi=0.5d0
   endif
   axj_scont = (1.d0-tx)**2 + 2.d0*tx*(1.d0-tx)*wim1
   bxj_scont = tx**2 + 2.d0*tx*(1.d0-tx)*wi
!
   mm=abs((int_im1j-int_im2j)/dxim1)
   mp=abs((int_ij-int_im1j)/dxi)
   if(mm.gt.mp) then
      wim1=1.d0-fac
      wi=fac
   elseif(mm.lt.mp) then
      wim1=fac
      wi=1.d0-fac
   else
      wim1=0.5d0
      wi=0.5d0
   endif
   axj_int = (1.d0-tx)**2 + 2.d0*tx*(1.d0-tx)*wim1
   bxj_int = tx**2 + 2.d0*tx*(1.d0-tx)*wi
!
!level j-1
   mm=abs((scont_im1jm1-scont_im2jm1)/dxim1)
   mp=abs((scont_ijm1-scont_im1jm1)/dxi)
   if(mm.gt.mp) then
      wim1=1.d0-fac
      wi=fac
   elseif(mm.lt.mp) then
      wim1=fac
      wi=1.d0-fac
   else
      wim1=0.5d0
      wi=0.5d0
   endif
   axjm1_scont = (1.d0-tx)**2 + 2.d0*tx*(1.d0-tx)*wim1
   bxjm1_scont = tx**2 + 2.d0*tx*(1.d0-tx)*wi
!
   mm=abs((int_im1jm1-int_im2jm1)/dxim1)
   mp=abs((int_ijm1-int_im1jm1)/dxi)
   if(mm.gt.mp) then
      wim1=1.d0-fac
      wi=fac
   elseif(mm.lt.mp) then
      wim1=fac
      wi=1.d0-fac
   else
      wim1=0.5d0
      wi=0.5d0
   endif
   axjm1_int = (1.d0-tx)**2 + 2.d0*tx*(1.d0-tx)*wim1
   bxjm1_int = tx**2 + 2.d0*tx*(1.d0-tx)*wi
!
!level j-2
   mm=abs((scont_im1jm2-scont_im2jm2)/dxim1)
   mp=abs((scont_ijm2-scont_im1jm2)/dxi)
   if(mm.gt.mp) then
      wim1=1.d0-fac
      wi=fac
   elseif(mm.lt.mp) then
      wim1=fac
      wi=1.d0-fac
   else
      wim1=0.5d0
      wi=0.5d0
   endif
   axjm2_scont = (1.d0-tx)**2 + 2.d0*tx*(1.d0-tx)*wim1
   bxjm2_scont = tx**2 + 2.d0*tx*(1.d0-tx)*wi
!
   mm=abs((int_im1jm2-int_im2jm2)/dxim1)
   mp=abs((int_ijm2-int_im1jm2)/dxi)
   if(mm.gt.mp) then
      wim1=1.d0-fac
      wi=fac
   elseif(mm.lt.mp) then
      wim1=fac
      wi=1.d0-fac
   else
      wim1=0.5d0
      wi=0.5d0
   endif
   axjm2_int = (1.d0-tx)**2 + 2.d0*tx*(1.d0-tx)*wim1
   bxjm2_int = tx**2 + 2.d0*tx*(1.d0-tx)*wi
!
!calculate value on each j-level
   scont_jm2 = axjm2_scont*scont_im1jm2 + bxjm2_scont*scont_ijm2
   scont_jm1 = axjm1_scont*scont_im1jm1 + bxjm1_scont*scont_ijm1
   scont_j   = axj_scont*scont_im1j     + bxj_scont*scont_ij
   int_jm2 = axjm2_int*int_im1jm2 + bxjm2_int*int_ijm2
   int_jm1 = axjm1_int*int_im1jm1 + bxjm1_int*int_ijm1
   int_j   = axj_int*int_im1j     + bxj_int*int_ij
!
!interpolation along y
   ay = (1.d0-ty)**2 + 2.d0*ty*(1.d0-ty)*wjm1
   by = ty**2 + 2.d0*ty*(1.d0-ty)*wj
!
   mm=abs((scont_jm1-scont_jm2)/dyjm1)
   mp=abs((scont_j-scont_jm1)/dyj)
   if(mm.gt.mp) then
      wim1=1.d0-fac
      wi=fac
   elseif(mm.lt.mp) then
      wim1=fac
      wi=1.d0-fac
   else
      wim1=0.5d0
      wi=0.5d0
   endif
   ay_scont = (1.d0-ty)**2 + 2.d0*ty*(1.d0-ty)*wim1
   by_scont = ty**2 + 2.d0*ty*(1.d0-ty)*wi
!
   mm=abs((int_jm1-int_jm2)/dyjm1)
   mp=abs((int_j-int_jm1)/dyj)
   if(mm.gt.mp) then
      wim1=1.d0-fac
      wi=fac
   elseif(mm.lt.mp) then
      wim1=fac
      wi=1.d0-fac
   else
      wim1=0.5d0
      wi=0.5d0
   endif
   ay_int = (1.d0-ty)**2 + 2.d0*ty*(1.d0-ty)*wim1
   by_int = ty**2 + 2.d0*ty*(1.d0-ty)*wi
!
   a_scont=0.d0
   b_scont=0.d0
   c_scont=0.d0
   d_scont=0.d0
   e_scont=ay_scont*axjm1_scont
   f_scont=ay_scont*bxjm1_scont
   g_scont=0.d0
   h_scont=by_scont*axj_scont
   i_scont=by_scont*bxj_scont
!
   a_inten=0.d0
   b_inten=0.d0
   c_inten=0.d0
   d_inten=0.d0
   e_inten=ay_int*axjm1_int
   f_inten=ay_int*bxjm1_int
   g_inten=0.d0
   h_inten=by_int*axj_int
   i_inten=by_int*bxj_int
!
   scont_p = ay_scont*scont_jm1 + by_scont*scont_j
   int_p   = ay_int*int_jm1 + by_int*int_j
!
!scont_p2 = e_scont*scont_im1jm1 + f_scont*scont_ijm1 + h_scont*scont_im1j + i_scont*scont_ij
!int_p2 = e_inten*int_im1jm1 + f_inten*int_ijm1 + h_inten*int_im1j + i_inten*int_ij
!
!write(*,'(4es20.8)') scont_p, scont_p2, int_p, int_p2
!if(abs(scont_p-scont_p2).gt.1.d-14) stop 'error in coeff2d_contu'
!if(abs(int_p-int_p2).gt.1.d-14) stop 'error in coeff2d_contu'
!
   return
!
!-------------------4-point bicubic interpolation-----------------------
!--------------with zero derivatives at all grid points-----------------
!
100 continue
!
   dxi=x_i-x_im1
   dyj=y_j-y_jm1
!
   tx = (x_p-x_im1)/dxi
   ty = (y_p-y_jm1)/dyj
!
   ax = 1.d0 - 3.d0*tx**2 + 2.d0*tx**3
   bx = 1.d0-ax
!
   ay = 1.d0 - 3.d0*ty**2 + 2.d0*ty**3
   by = 1.d0-ay

!
   a_scont=0.d0
   b_scont=0.d0
   c_scont=0.d0
   d_scont=0.d0
   e_scont=ay*ax
   f_scont=ay*bx
   g_scont=0.d0
   h_scont=by*ax
   i_scont=by*bx
!
   a_inten=0.d0
   b_inten=0.d0
   c_inten=0.d0
   d_inten=0.d0
   e_inten=e_scont
   f_inten=f_scont
   g_inten=0.d0
   h_inten=h_scont
   i_inten=i_scont
!
   scont_p = e_scont*scont_im1jm1 + f_scont*scont_ijm1 + h_scont*scont_im1j + i_scont*scont_ij
   int_p = e_inten*int_im1jm1 + f_inten*int_ijm1 + h_inten*int_im1j + i_inten*int_ij
!
   return
!
!---------------monotonic quadratic bezier interpolation----------------
!--------(with derivatives for control point from weighted mean)--------
!------------------where weights are assigned---------------------------
!
110 continue
!
   ax = (1.d0-tx)**2
   bx = 2.d0*tx*(1.d0-tx)
   cx = tx**2
!
!calculate control points on each j-level:
!   scontc_jm2, scontc_jm1, scontc_j
!   intc_jm2, intc_jm1, intc_j
!   opacc_jm2, opacc_jm1, opacc_j
!
   fac=max(wp_interp2d,dxim1/dx)
   fac2=dxim1/dx
   axt = (fac-1.d0)*dxi/2.d0/dxim1
   bxt = ((2.d0-fac)*dxim1 + (1.d0-fac)*dxi)/2.d0/dxim1
   cxt = fac/2.d0
   axt2 = (fac2-1.d0)*dxi/2.d0/dxim1
   bxt2 = ((2.d0-fac2)*dxim1 + (1.d0-fac2)*dxi)/2.d0/dxim1
   cxt2 = fac2/2.d0
!
   opacc_jm2 = opac_im2jm2*axt2 + opac_im1jm2*bxt2 + opac_ijm2*cxt2
   opacc_jm1 = opac_im2jm1*axt2 + opac_im1jm1*bxt2 + opac_ijm1*cxt2
   opacc_j   = opac_im2j*axt2   + opac_im1j*bxt2   + opac_ij*cxt2
!
   scontc_jm2 = scont_im2jm2*axt + scont_im1jm2*bxt + scont_ijm2*cxt
   scontc_jm1 = scont_im2jm1*axt + scont_im1jm1*bxt + scont_ijm1*cxt
   scontc_j   = scont_im2j*axt   + scont_im1j*bxt   + scont_ij*cxt
!
   intc_jm2 = int_im2jm2*axt + int_im1jm2*bxt + int_ijm2*cxt
   intc_jm1 = int_im2jm1*axt + int_im1jm1*bxt + int_ijm1*cxt
   intc_j   = int_im2j*axt   + int_im1j*bxt   + int_ij*cxt
!
!ensure monotonicity on level j
   call coeffcr1d_mbez(opacc_j, opac_im1j, opac_ij, axt2, bxt2, cxt2, axtj_opac, bxtj_opac, cxtj_opac)
   call coeffcr1d_mbez(scontc_j, scont_im1j, scont_ij, axt, bxt, cxt, axtj_scont, bxtj_scont, cxtj_scont)
   call coeffcr1d_mbez(intc_j, int_im1j, int_ij, axt, bxt, cxt, axtj_int, bxtj_int, cxtj_int)
   axj_opac = axtj_opac*bx
   bxj_opac = bxtj_opac*bx + ax
   cxj_opac = cxtj_opac*bx + cx
   opac_j = axj_opac*opac_im2j + bxj_opac*opac_im1j + cxj_opac*opac_ij
   axj_scont = axtj_scont*bx
   bxj_scont = bxtj_scont*bx + ax
   cxj_scont = cxtj_scont*bx + cx
   scont_j = axj_scont*scont_im2j + bxj_scont*scont_im1j + cxj_scont*scont_ij
   axj_int = axtj_int*bx
   bxj_int = bxtj_int*bx + ax
   cxj_int = cxtj_int*bx + cx
   int_j = axj_int*int_im2j + bxj_int*int_im1j + cxj_int*int_ij
!
!ensure monotonicity on level j-1
   call coeffcr1d_mbez(opacc_jm1, opac_im1jm1, opac_ijm1, axt2, bxt2, cxt2, axtjm1_opac, bxtjm1_opac, cxtjm1_opac)
   call coeffcr1d_mbez(scontc_jm1, scont_im1jm1, scont_ijm1, axt, bxt, cxt, axtjm1_scont, bxtjm1_scont, cxtjm1_scont)
   call coeffcr1d_mbez(intc_jm1, int_im1jm1, int_ijm1, axt, bxt, cxt, axtjm1_int, bxtjm1_int, cxtjm1_int)
   axjm1_opac = axtjm1_opac*bx
   bxjm1_opac = bxtjm1_opac*bx + ax
   cxjm1_opac = cxtjm1_opac*bx + cx
   opac_jm1 = axjm1_opac*opac_im2jm1 + bxjm1_opac*opac_im1jm1 + cxjm1_opac*opac_ijm1
   axjm1_scont = axtjm1_scont*bx
   bxjm1_scont = bxtjm1_scont*bx + ax
   cxjm1_scont = cxtjm1_scont*bx + cx
   scont_jm1 = axjm1_scont*scont_im2jm1 + bxjm1_scont*scont_im1jm1 + cxjm1_scont*scont_ijm1
   axjm1_int = axtjm1_int*bx
   bxjm1_int = bxtjm1_int*bx + ax
   cxjm1_int = cxtjm1_int*bx + cx
   int_jm1 = axjm1_int*int_im2jm1 + bxjm1_int*int_im1jm1 + cxjm1_int*int_ijm1
!
!ensure monotonicity on level j-2
   call coeffcr1d_mbez(opacc_jm2, opac_im1jm2, opac_ijm2, axt2, bxt2, cxt2, axtjm2_opac, bxtjm2_opac, cxtjm2_opac)
   call coeffcr1d_mbez(scontc_jm2, scont_im1jm2, scont_ijm2, axt, bxt, cxt, axtjm2_scont, bxtjm2_scont, cxtjm2_scont)
   call coeffcr1d_mbez(intc_jm2, int_im1jm2, int_ijm2, axt, bxt, cxt, axtjm2_int, bxtjm2_int, cxtjm2_int)
   axjm2_opac = axtjm2_opac*bx
   bxjm2_opac = bxtjm2_opac*bx + ax
   cxjm2_opac = cxtjm2_opac*bx + cx
   opac_jm2 = axjm2_opac*opac_im2jm2 + bxjm2_opac*opac_im1jm2 + cxjm2_opac*opac_ijm2
   axjm2_scont = axtjm2_scont*bx
   bxjm2_scont = bxtjm2_scont*bx + ax
   cxjm2_scont = cxtjm2_scont*bx + cx
   scont_jm2 = axjm2_scont*scont_im2jm2 + bxjm2_scont*scont_im1jm2 + cxjm2_scont*scont_ijm2
   axjm2_int = axtjm2_int*bx
   bxjm2_int = bxtjm2_int*bx + ax
   cxjm2_int = cxtjm2_int*bx + cx
   int_jm2 = axjm2_int*int_im2jm2 + bxjm2_int*int_im1jm2 + cxjm2_int*int_ijm2
!
!
!
!calculate control point for interpolation along y
   fac=max(wp_interp2d,dyjm1/dy)
   fac2=dyjm1/dy
   ayt = (fac-1.d0)*dyj/2.d0/dyjm1
   byt = ((2.d0-fac)*dyjm1 + (1.d0-fac)*dyj)/2.d0/dyjm1
   cyt = fac/2.d0
   ayt2 = (fac2-1.d0)*dyj/2.d0/dyjm1
   byt2 = ((2.d0-fac2)*dyjm1 + (1.d0-fac2)*dyj)/2.d0/dyjm1
   cyt2 = fac2/2.d0
   opac_c = opac_jm2*ayt2 + opac_jm1*byt2 + opac_j*cyt2
   scont_c = scont_jm2*ayt + scont_jm1*byt + scont_j*cyt
   int_c = int_jm2*ayt + int_jm1*byt + int_j*cyt
!
!ensure monotonicity
   call coeffcr1d_mbez(opac_c, opac_jm1, opac_j, ayt2, byt2, cyt2, ayt_opac, byt_opac, cyt_opac)
   call coeffcr1d_mbez(scont_c, scont_jm1, scont_j, ayt, byt, cyt, ayt_scont, byt_scont, cyt_scont)
   call coeffcr1d_mbez(int_c, int_jm1, int_j, ayt, byt, cyt, ayt_int, byt_int, cyt_int)
!
   ay = (1.d0-ty)**2
   by = 2.d0*ty*(1.d0-ty)
   cy = ty**2
!
   ay_opac = ayt_opac*by
   by_opac = byt_opac*by + ay
   cy_opac = cyt_opac*by + cy
   opac_p = ay_opac*opac_jm2 + by_opac*opac_jm1 + cy_opac*opac_j
   ay_scont = ayt_scont*by
   by_scont = byt_scont*by + ay
   cy_scont = cyt_scont*by + cy
   scont_p = ay_scont*scont_jm2 + by_scont*scont_jm1 + cy_scont*scont_j
   ay_int = ayt_int*by
   by_int = byt_int*by + ay
   cy_int = cyt_int*by + cy
   int_p = ay_int*int_jm2 + by_int*int_jm1 + cy_int*int_j
!
   a_scont = ay_scont*axjm2_scont
   b_scont = ay_scont*bxjm2_scont
   c_scont = ay_scont*cxjm2_scont
   d_scont = by_scont*axjm1_scont
   e_scont = by_scont*bxjm1_scont
   f_scont = by_scont*cxjm1_scont
   g_scont = cy_scont*axj_scont
   h_scont = cy_scont*bxj_scont
   i_scont = cy_scont*cxj_scont
!
   a_inten = ay_int*axjm2_int
   b_inten = ay_int*bxjm2_int
   c_inten = ay_int*cxjm2_int
   d_inten = by_int*axjm1_int
   e_inten = by_int*bxjm1_int
   f_inten = by_int*cxjm1_int
   g_inten = cy_int*axj_int
   h_inten = cy_int*bxj_int
   i_inten = cy_int*cxj_int
!
!write(*,'(9es20.8)') a_scont, b_scont, c_scont, d_scont, e_scont, f_scont, g_scont, h_scont, i_scont
!
   return
!
end subroutine coeff2d_contu
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine coeff2d_contd(opac_im2jm2, opac_im1jm2, opac_ijm2, &
   opac_im2jm1, opac_im1jm1, opac_ijm1, &
   opac_im2j,   opac_im1j,   opac_ij, &
   scont_im2jm2, scont_im1jm2, scont_ijm2, &
   scont_im2jm1, scont_im1jm1, scont_ijm1, &
   scont_im2j,   scont_im1j,   scont_ij, &
   x_im2, x_im1, x_i, y_jm2, y_jm1, y_j, x_p, y_p, &
   a_scont, b_scont, c_scont, d_scont, e_scont, &
   f_scont, g_scont, h_scont, i_scont, &
   opac_p, scont_p)
!
!            interpolates opacity and continuum source function
!               values given on a 2d grid onto point x_p, y_p
!
!on input (f_* stands for opac_* and scont_* respectivly):
!
! y_j      f_im2j----------f_im1j--------------f_ij
!  |          |              |                  |
!  |          |              |                  |
!  |          |              |         x        |
!  |          |              |     (x_p,y_p)    |
!  |          |              |                  |
!y_jm1    f_im2jm1-------f_im1jm1------------f_ijm1
!  |          |              |                  |
!  |          |              |                  |
!  |          |              |                  |
!  |          |              |                  |
!  |          |              |                  |
!y_jm2    f_im2jm2-------f_im1jm2------------f_ijm2
!  |
!  --------x_im2-----------x_im1---------------x_i
!
!        x_p, y_p: coordinates of point onto which shall be interpolated
!
!on output:
!   1. interpolation coefficients for source function
!      (required for ALO calculations):
!         a_scont, b_scont, c_scont, d_scont, e_scont
!         f_scont, g_scont, h_scont, i_scont
!
!      such that:
!         f_p = a*f_im2jm2 + b*f_im1jm2 + c*f_ijm2 +
!               d*f_im2jm1 + e*f_im1jm1 + f*f_ijm1 +
!               g*f_im2j   + h*f_im1j   + i*f_ij
!
!   2. interpolated values at point p: opac_p, scont_p
!
   use prog_type
   use mod_interp2d, only: wp_interp2d
!
   implicit none
!
! ... arguments
   real(dp), intent(in) :: opac_im2jm2, opac_im1jm2, opac_ijm2, &
      opac_im2jm1, opac_im1jm1, opac_ijm1, &
      opac_im2j,   opac_im1j,   opac_ij, &
      scont_im2jm2, scont_im1jm2, scont_ijm2, &
      scont_im2jm1, scont_im1jm1, scont_ijm1, &
      scont_im2j,   scont_im1j,   scont_ij, &
      x_im2, x_im1, x_i, y_jm2, y_jm1, y_j, x_p, y_p
   real(dp), intent(out) :: a_scont, b_scont, c_scont, d_scont, e_scont, &
      f_scont, g_scont, h_scont, i_scont, &
      opac_p, scont_p
!
! ... local scalars
   real(dp) :: dxim1, dxi, dx, tx, ax, bx, cx, &
      dyjm1, dyj, dy, ty, ay, by, cy, &
      axt, bxt, cxt, axt2, bxt2, cxt2, &
      ayt, byt, cyt, ayt2, byt2, cyt2, &
      axj_opac, bxj_opac, cxj_opac, &
      axjm1_opac, bxjm1_opac, cxjm1_opac, &
      axjm2_opac, bxjm2_opac, cxjm2_opac, &
      axtj_opac, bxtj_opac, cxtj_opac, &
      axtjm1_opac, bxtjm1_opac, cxtjm1_opac, &
      axtjm2_opac, bxtjm2_opac, cxtjm2_opac, &
      opacc_jm2, opacc_jm1, opacc_j, &
      opac_jm2, opac_jm1, opac_j, opac_c, &
      ayt_opac, byt_opac, cyt_opac, &
      ay_opac, by_opac, cy_opac, &
      axj_scont, bxj_scont, cxj_scont, &
      axjm1_scont, bxjm1_scont, cxjm1_scont, &
      axjm2_scont, bxjm2_scont, cxjm2_scont, &
      axtj_scont, bxtj_scont, cxtj_scont, &
      axtjm1_scont, bxtjm1_scont, cxtjm1_scont, &
      axtjm2_scont, bxtjm2_scont, cxtjm2_scont, &
      scontc_jm2, scontc_jm1, scontc_j, &
      scont_jm2, scont_jm1, scont_j, scont_c, &
      ayt_scont, byt_scont, cyt_scont, &
      ay_scont, by_scont, cy_scont, &
      rdxdy
   real(dp) :: d_im2jm2, d_im1jm2, d_ijm2, d_im2jm1, d_im1jm1, d_ijm1, d_im2j, d_im1j, d_ij, &
      w_im2jm2, w_im1jm2, w_ijm2, w_im2jm1, w_im1jm1, w_ijm1, w_im2j, w_im1j, w_ij, &
      norm, fac, fac2
   real(dp) :: dx1, dx2, dx3, dx4, dx1s, dx2s, dx3s, dx4s, &
      dy1, dy2, dy3, dy4, dy1s, dy2s, dy3s, dy4s, &
      norm01, norm10, norm12, norm21, norm11, &
      dxis, dyjs, d1, d2, d3, d4, d5
   real(dp) :: dxp, dxp2, dxp3, dxi2, dxim12, &
      dyp, dyp2, dyp3, dyj2, dyjm12, &
      wi, wim1, wj, wjm1, mm, mp, a, b, c
!


!define deltax, deltay
   dxim1 = x_im1-x_im2
   dxi = x_i-x_im1
   dx = dxim1+dxi
   dyjm1 = y_jm1-y_jm2
   dyj = y_j-y_jm1
   dy = dyjm1+dyj
!
!define deltax, deltay-ratios
   tx = (x_p-x_im1)/dxi
   ty = (y_p-y_jm1)/dyj
!
!
!10:  bilinear interpolation
!20:  inverse distance weighting
!30:  bezier inteprolation (with control points from derivatives)
!40:  bezier interpolation (with control points from derivatives, non-monotonic)
!50:  bezier interpolation (with control point from idw)
!60:  cubic interpolation with free central points
!70:  biquadratic interpolation
!80:  bezier interpolation (control point from assigning weight to each point)
!90:  bezier interpolation (control point from assigning weight to each point and approximated curvature)
!100: bicubic interpolation (with zero derivatives at all points)
!110: bezier interpolation (with conrol points from derivatives with minimum weight for left derivative, monotonic)
   goto 110
!
!-------------------------bilinear interpolation------------------------
!
10 continue
!
   rdxdy=tx*ty
!
   a_scont=0.d0
   b_scont=0.d0
   c_scont=0.d0
   d_scont=0.d0
   e_scont=1.d0-tx-ty+rdxdy
   f_scont=tx-rdxdy
   g_scont=0.d0
   h_scont=ty-rdxdy
   i_scont=rdxdy
!
!opac_p = e_scont*opac_im1jm1 + f_scont*opac_ijm1 + h_scont*opac_im1j + i_scont*opac_ij
   scont_p = e_scont*scont_im1jm1 + f_scont*scont_ijm1 + h_scont*scont_im1j + i_scont*scont_ij

   return
!
!-------------------------inverse distance weighting--------------------
!
20 continue
!
!define weighting factor
   fac=2.d0
!
!calculate distance to each point
   d_im2jm2 = (x_p-x_im2)**2+(y_p-y_jm2)**2
   d_im1jm2 = (x_p-x_im1)**2+(y_p-y_jm2)**2
   d_ijm2 = (x_p-x_i)**2+(y_p-y_jm2)**2
   d_im2jm1 = (x_p-x_im2)**2+(y_p-y_jm1)**2
   d_im1jm1 = (x_p-x_im1)**2+(y_p-y_jm1)**2
   d_ijm1 = (x_p-x_i)**2+(y_p-y_jm1)**2
   d_im2j = (x_p-x_im2)**2+(y_p-y_j)**2
   d_im1j = (x_p-x_im1)**2+(y_p-y_j)**2
   d_ij = (x_p-x_i)**2+(y_p-y_j)**2
!
!avoid division by zero if d_ij=0, or d_im1j=0, ...
   if(d_im2jm2.eq.0.) then
      w_im2jm2=1.d0
      w_im1jm2=0.d0
      w_ijm2=0.d0
      w_im2jm1=0.d0
      w_im1jm1=0.d0
      w_ijm1=0.d0
      w_im2j=0.d0
      w_im1j=0.d0
      w_ij=0.d0
   elseif(d_im1jm2.eq.0.) then
      w_im2jm2=0.d0
      w_im1jm2=1.d0
      w_ijm2=0.d0
      w_im2jm1=0.d0
      w_im1jm1=0.d0
      w_ijm1=0.d0
      w_im2j=0.d0
      w_im1j=0.d0
      w_ij=0.d0
   elseif(d_ijm2.eq.0.) then
      w_im2jm2=0.d0
      w_im1jm2=0.d0
      w_ijm2=1.d0
      w_im2jm1=0.d0
      w_im1jm1=0.d0
      w_ijm1=0.d0
      w_im2j=0.d0
      w_im1j=0.d0
      w_ij=0.d0
   elseif(d_im2jm1.eq.0.) then
      w_im2jm2=0.d0
      w_im1jm2=0.d0
      w_ijm2=0.d0
      w_im2jm1=1.d0
      w_im1jm1=0.d0
      w_ijm1=0.d0
      w_im2j=0.d0
      w_im1j=0.d0
      w_ij=0.d0
   elseif(d_im1jm1.eq.0.) then
      w_im2jm2=0.d0
      w_im1jm2=0.d0
      w_ijm2=0.d0
      w_im2jm1=0.d0
      w_im1jm1=1.d0
      w_ijm1=0.d0
      w_im2j=0.d0
      w_im1j=0.d0
      w_ij=0.d0
   elseif(d_ijm1.eq.0.) then
      w_im2jm2=0.d0
      w_im1jm2=0.d0
      w_ijm2=0.d0
      w_im2jm1=0.d0
      w_im1jm1=0.d0
      w_ijm1=1.d0
      w_im2j=0.d0
      w_im1j=0.d0
      w_ij=0.d0
   elseif(d_im2j.eq.0.) then
      w_im2jm2=0.d0
      w_im1jm2=0.d0
      w_ijm2=0.d0
      w_im2jm1=0.d0
      w_im1jm1=0.d0
      w_ijm1=0.d0
      w_im2j=1.d0
      w_im1j=0.d0
      w_ij=0.d0
   elseif(d_im1j.eq.0.) then
      w_im2jm2=0.d0
      w_im1jm2=0.d0
      w_ijm2=0.d0
      w_im2jm1=0.d0
      w_im1jm1=0.d0
      w_ijm1=0.d0
      w_im2j=0.d0
      w_im1j=1.d0
      w_ij=0.d0
   elseif(d_ij.eq.0.) then
      w_im2jm2=0.d0
      w_im1jm2=0.d0
      w_ijm2=0.d0
      w_im2jm1=0.d0
      w_im1jm1=0.d0
      w_ijm1=0.d0
      w_im2j=0.d0
      w_im1j=0.d0
      w_ij=1.d0
   else
!
      w_im2jm2=1.d0/d_im2jm2**(fac/2.d0)
      w_im1jm2=1.d0/d_im1jm2**(fac/2.d0)
      w_ijm2=1.d0/d_ijm2**(fac/2.d0)
      w_im2jm1=1.d0/d_im2jm1**(fac/2.d0)
      w_im1jm1=1.d0/d_im1jm1**(fac/2.d0)
      w_ijm1=1.d0/d_ijm1**(fac/2.d0)
      w_im2j=1.d0/d_im2j**(fac/2.d0)
      w_im1j=1.d0/d_im1j**(fac/2.d0)
      w_ij=1.d0/d_ij**(fac/2.d0)
!
   endif
!
   norm = w_im2jm2 + w_im1jm2 + w_ijm2 + w_im2jm1 + w_im1jm1 + w_ijm1 + w_im2j + w_im1j + w_ij
!
   a_scont = w_im2jm2/norm
   b_scont = w_im1jm2/norm
   c_scont = w_ijm2/norm
   d_scont = w_im2jm1/norm
   e_scont = w_im1jm1/norm
   f_scont = w_ijm1/norm
   g_scont = w_im2j/norm
   h_scont = w_im1j/norm
   i_scont = w_ij/norm
!

!opac_p = a_scont*opac_im2jm2 + b_scont*opac_im1jm2 + c_scont*opac_ijm2 + &
!         d_scont*opac_im2jm1 + e_scont*opac_im1jm1 + f_scont*opac_ijm1 + &
!         g_scont*opac_im2j   + h_scont*opac_im1j   + i_scont*opac_ij
   scont_p = a_scont*scont_im2jm2 + b_scont*scont_im1jm2 + c_scont*scont_ijm2 + &
      d_scont*scont_im2jm1 + e_scont*scont_im1jm1 + f_scont*scont_ijm1 + &
      g_scont*scont_im2j   + h_scont*scont_im1j   + i_scont*scont_ij
   return
!
!---------------monotonic quadratic bezier interpolation----------------
!--------(with derivatives for control point from weighted mean)--------
!
30 continue
!
   ax = (1.d0-tx)**2
   bx = 2.d0*tx*(1.d0-tx)
   cx = tx**2
!
!calculate control points on each j-level:
!   scontc_jm2, scontc_jm1, scontc_j
!   opacc_jm2, opacc_jm1, opacc_j
   axt = -dxi**2/2.d0/dxim1/dx
   bxt = dx/2.d0/dxim1
   cxt = dxim1/2.d0/dx
!
!opacc_jm2 = opac_im2jm2*axt + opac_im1jm2*bxt + opac_ijm2*cxt
!opacc_jm1 = opac_im2jm1*axt + opac_im1jm1*bxt + opac_ijm1*cxt
!opacc_j   = opac_im2j*axt   + opac_im1j*bxt   + opac_ij*cxt
!
   scontc_jm2 = scont_im2jm2*axt + scont_im1jm2*bxt + scont_ijm2*cxt
   scontc_jm1 = scont_im2jm1*axt + scont_im1jm1*bxt + scont_ijm1*cxt
   scontc_j   = scont_im2j*axt   + scont_im1j*bxt   + scont_ij*cxt
!
!ensure monotonicity on level j
!call coeffcr1d_mbez(opacc_j, opac_im1j, opac_ij, axt, bxt, cxt, axtj_opac, bxtj_opac, cxtj_opac)
   call coeffcr1d_mbez(scontc_j, scont_im1j, scont_ij, axt, bxt, cxt, axtj_scont, bxtj_scont, cxtj_scont)
!axj_opac = axtj_opac*bx
!bxj_opac = bxtj_opac*bx + ax
!cxj_opac = cxtj_opac*bx + cx
!opac_j = axj_opac*opac_im2j + bxj_opac*opac_im1j + cxj_opac*opac_ij
   axj_scont = axtj_scont*bx
   bxj_scont = bxtj_scont*bx + ax
   cxj_scont = cxtj_scont*bx + cx
   scont_j = axj_scont*scont_im2j + bxj_scont*scont_im1j + cxj_scont*scont_ij
!
!ensure monotonicity on level j-1
!call coeffcr1d_mbez(opacc_jm1, opac_im1jm1, opac_ijm1, axt, bxt, cxt, axtjm1_opac, bxtjm1_opac, cxtjm1_opac)
   call coeffcr1d_mbez(scontc_jm1, scont_im1jm1, scont_ijm1, axt, bxt, cxt, axtjm1_scont, bxtjm1_scont, cxtjm1_scont)
!axjm1_opac = axtjm1_opac*bx
!bxjm1_opac = bxtjm1_opac*bx + ax
!cxjm1_opac = cxtjm1_opac*bx + cx
!opac_jm1 = axjm1_opac*opac_im2jm1 + bxjm1_opac*opac_im1jm1 + cxjm1_opac*opac_ijm1
   axjm1_scont = axtjm1_scont*bx
   bxjm1_scont = bxtjm1_scont*bx + ax
   cxjm1_scont = cxtjm1_scont*bx + cx
   scont_jm1 = axjm1_scont*scont_im2jm1 + bxjm1_scont*scont_im1jm1 + cxjm1_scont*scont_ijm1
!
!ensure monotonicity on level j-2
!call coeffcr1d_mbez(opacc_jm2, opac_im1jm2, opac_ijm2, axt, bxt, cxt, axtjm2_opac, bxtjm2_opac, cxtjm2_opac)
   call coeffcr1d_mbez(scontc_jm2, scont_im1jm2, scont_ijm2, axt, bxt, cxt, axtjm2_scont, bxtjm2_scont, cxtjm2_scont)
!axjm2_opac = axtjm2_opac*bx
!bxjm2_opac = bxtjm2_opac*bx + ax
!cxjm2_opac = cxtjm2_opac*bx + cx
!opac_jm2 = axjm2_opac*opac_im2jm2 + bxjm2_opac*opac_im1jm2 + cxjm2_opac*opac_ijm2
   axjm2_scont = axtjm2_scont*bx
   bxjm2_scont = bxtjm2_scont*bx + ax
   cxjm2_scont = cxtjm2_scont*bx + cx
   scont_jm2 = axjm2_scont*scont_im2jm2 + bxjm2_scont*scont_im1jm2 + cxjm2_scont*scont_ijm2
!
!
!
!calculate control point for interpolation along y
   ayt = -dyj**2/2.d0/dyjm1/dy
   byt = dy/2.d0/dyjm1
   cyt = dyjm1/2.d0/dy
!opac_c = opac_jm2*ayt + opac_jm1*byt + opac_j*cyt
   scont_c = scont_jm2*ayt + scont_jm1*byt + scont_j*cyt
!
!ensure monotonicity
!call coeffcr1d_mbez(opac_c, opac_jm1, opac_j, ayt, byt, cyt, ayt_opac, byt_opac, cyt_opac)
   call coeffcr1d_mbez(scont_c, scont_jm1, scont_j, ayt, byt, cyt, ayt_scont, byt_scont, cyt_scont)
!
   ay = (1.d0-ty)**2
   by = 2.d0*ty*(1.d0-ty)
   cy = ty**2
!
!ay_opac = ayt_opac*by
!by_opac = byt_opac*by + ay
!cy_opac = cyt_opac*by + cy
!opac_p = ay_opac*opac_jm2 + by_opac*opac_jm1 + cy_opac*opac_j
   ay_scont = ayt_scont*by
   by_scont = byt_scont*by + ay
   cy_scont = cyt_scont*by + cy
   scont_p = ay_scont*scont_jm2 + by_scont*scont_jm1 + cy_scont*scont_j
!
   a_scont = ay_scont*axjm2_scont
   b_scont = ay_scont*bxjm2_scont
   c_scont = ay_scont*cxjm2_scont
   d_scont = by_scont*axjm1_scont
   e_scont = by_scont*bxjm1_scont
   f_scont = by_scont*cxjm1_scont
   g_scont = cy_scont*axj_scont
   h_scont = cy_scont*bxj_scont
   i_scont = cy_scont*cxj_scont
!
   return
!
!-------------quadratic bezier interpolation (non-monotonic)------------
!--------(with derivatives for control point from weighted mean)--------
!
40 continue
!
   ax = (1.d0-tx)**2
   bx = 2.d0*tx*(1.d0-tx)
   cx = tx**2
!
!calculate control points on each j-level:
!   scontc_jm2, scontc_jm1, scontc_j
!   intc_jm2, intc_jm1, intc_j
!   opacc_jm2, opacc_jm1, opacc_j
   axt = -dxi**2/2.d0/dxim1/dx
   bxt = dx/2.d0/dxim1
   cxt = dxim1/2.d0/dx
!
   axj_scont = axt*bx
   bxj_scont = ax+bx*bxt
   cxj_scont = cx+bx*cxt
!
   ay = (1.d0-ty)**2
   by = 2.d0*ty*(1.d0-ty)
   cy = ty**2
!
   ayt = -dyj**2/2.d0/dyjm1/dy
   byt = dy/2.d0/dyjm1
   cyt = dyjm1/2.d0/dy
!
   ay_scont = ayt*by
   by_scont = ay+by*byt
   cy_scont = cy+by*cyt
!
   a_scont = ay_scont*axj_scont
   b_scont = ay_scont*bxj_scont
   c_scont = ay_scont*cxj_scont
   d_scont = by_scont*axj_scont
   e_scont = by_scont*bxj_scont
   f_scont = by_scont*cxj_scont
   g_scont = cy_scont*axj_scont
   h_scont = cy_scont*bxj_scont
   i_scont = cy_scont*cxj_scont
!
!
!
!opac_p = a_scont*opac_im2jm2 + b_scont*opac_im1jm2 + c_scont*opac_ijm2 + &
!         d_scont*opac_im2jm1 + e_scont*opac_im1jm1 + f_scont*opac_ijm1 + &
!         g_scont*opac_im2j   + h_scont*opac_im1j   + i_scont*opac_ij
   scont_p = a_scont*scont_im2jm2 + b_scont*scont_im1jm2 + c_scont*scont_ijm2 + &
      d_scont*scont_im2jm1 + e_scont*scont_im1jm1 + f_scont*scont_ijm1 + &
      g_scont*scont_im2j   + h_scont*scont_im1j   + i_scont*scont_ij
!
!write(*,'(9es20.8)') a_scont, b_scont, c_scont, d_scont, e_scont, f_scont, g_scont, h_scont, i_scont
!
   return
!
!---------------monotonic quadratic bezier interpolation----------------
!----------(with control points from inverse distance weighing)---------
!
50 continue
!
   dxi=(x_i-x_im1)
   dyj=(y_j-y_jm1)
!
   tx = (x_p-x_im1)/dxi
   ty = (y_p-y_jm1)/dyj
!
   ax=(1.d0-tx)**2
   bx=2.d0*tx*(1.d0-tx)
   cx=tx**2
!
   ay=(1.d0-ty)**2
   by=2.d0*ty*(1.d0-ty)
   cy=ty**2
!
!calculate all control points via inverse distance weighting (weight = distance**fac)
   fac = -4.d0
!
!using only 4 local points
   dxis=dxi**2
   dyjs=dyj**2
   d1=(dxis/4.d0)**fac
   d2=(dyjs/4.d0)**fac
   d3=(dxis+dyjs/4.d0)**fac
   d4=(dxis/4.d0+dyjs)**fac
   d5=((dxis+dyjs)/4.d0)**fac
   norm01 = 2.d0*(d2+d3)
   norm10 = 2.d0*(d1+d4)
   norm12 = norm10
   norm21 = norm01
   norm11 = 4.d0*d5

   a_scont = 0.d0
   b_scont = 0.d0
   c_scont = 0.d0
   d_scont = 0.d0
   e_scont = ax*ay + ax*by*d2/norm01 + bx*ay*d1/norm10 + bx*by*d5/norm11 + bx*cy*d4/norm12 + cx*by*d3/norm21
   f_scont = cx*ay + ax*by*d3/norm01 + bx*ay*d1/norm10 + bx*by*d5/norm11 + bx*cy*d4/norm12 + cx*by*d2/norm21
   g_scont = 0.d0
   h_scont = ax*cy + ax*by*d2/norm01 + bx*ay*d4/norm10 + bx*by*d5/norm11 + bx*cy*d1/norm12 + cx*by*d3/norm21
   i_scont = cx*cy + ax*by*d3/norm01 + bx*ay*d4/norm10 + bx*by*d5/norm11 + bx*cy*d1/norm12 + cx*by*d2/norm21
!
!
!
!using all 9 points
!dx1 = x_i-x_im1
!dx2 = x_i-x_im2
!dx3 = (x_i+x_im1-2.d0*x_im2)/2.d0
!dx4 = x_im1-x_im2
!dx1s = dx1**2
!dx2s = dx2**2
!dx3s = dx3**2
!dx4s = dx4**2
!!
!dy1 = y_j-y_jm1
!dy2 = y_j-y_jm2
!dy3 = (y_j+y_jm1-2.d0*y_jm2)/2.d0
!dy4 = y_jm1-y_jm2
!dy1s = dy1**2
!dy2s = dy2**2
!dy3s = dy3**2
!dy4s = dy4**2
!!
!fac = -4.d0
!norm01 = (dx4s+dy3s)**fac + dy3s**fac + (dx1s+dy3s)**fac + (dx4s+dy1s/4.)**fac + &
!         (dy1s/4.)**fac + (dx1s+dy1s/4.)**fac + (dx4s+dy1s/4.)**fac + (dy1s/4.)**fac + (dx1s+dy1s/4.)**fac
!norm10 = (dx3s+dy4s)**fac + (dx1s/4.+dy4s)**fac + (dx1s/4.+dy4s)**fac + dx3s**fac + (dx1s/4.)**fac + &
!         (dx1s/4.)**fac + (dx3s+dy1s)**fac + (dx1s/4.+dy1s)**fac + (dx1s/4.+dy1s)**fac
!norm12 = (dx3s+dy2s)**fac + (dx1s/4.+dy2s)**fac + (dx1s/4.+dy2s)**fac + (dx3s+dy1s)**fac + &
!         (dx1s/4.+dy1s)**fac + (dx1s/4.+dy1s)**fac + dx3s**fac + (dx1s/4.)**fac + (dx1s/4.)**fac
!norm21 = (dx2s+dx3s)**fac + (dx1s+dy3s)**fac + dy3s**fac + (dx2s+dy1s/4.)**fac + (dx1s+dy1s/4.)**fac + &
!         (dy1s/4.)**fac + (dx2s+dy1s/4.)**fac + (dx1s+dy1s/4.)**fac + (dy1s/4.)**fac
!norm11 = (dx3s+dy3s)**fac + 2.*(dx1s/4.+dy3s)**fac + 2.*(dx3s+dy1s/4.)**fac + 4.*(dx1s/4.+dy1s/4.)**fac
!!
!a_scont = ax*by/norm01*(dx4s+dy3s)**fac + bx*ay/norm10*(dx3s+dy4s)**fac + bx*by/norm11*(dx3s+dy3s)**fac + &
!          bx*cy/norm12*(dx3s+dy2s)**fac + cx*by/norm21*(dx2s+dy3s)**fac
!b_scont = ax*by/norm01*dy3s**fac + bx*ay/norm10*(dx1s/4.+dy4s)**fac + bx*by/norm11*(dx1s/4.+dy3s)**fac + &
!          bx*cy/norm12*(dx1s/4.+dy2s)**fac + cx*by/norm21*(dx1s+dy3s)**fac
!c_scont = ax*by/norm01*(dx1s+dy3s)**fac + bx*ay/norm10*(dx1s/4.+dy4s)**fac + bx*by/norm11*(dx1s/4.+dy3s)**fac + &
!          bx*cy/norm12*(dx1s/4.+dy2s)**fac + cx*by/norm21*dy3s**fac
!d_scont = ax*by/norm01*(dx4s+dy1s/4.)**fac + bx*ay/norm10*dx3s**fac + bx*by/norm11*(dx3s+dy1s/4.)**fac + &
!          bx*cy/norm12*(dx3s+dy1s)**fac + cx*by/norm21*(dx2s+dy1s/4.)**fac
!e_scont = ax*ay + ax*by/norm01*(dy1s/4.)**fac + bx*ay/norm10*(dx1s/4.)**fac + bx*by/norm11*(dx1s/4.+dy1s/4.)**fac + &
!          bx*cy/norm12*(dx1s/4.+dy1s)**fac + cx*by/norm21*(dx1s+dy1s/4.)**fac
!f_scont = cx*ay + ax*by/norm01*(dx1s+dy1s/4.)**fac +bx*ay/norm10*(dx1s/4.)**fac + bx*by/norm11*(dx1s/4.+dy1s/4.)**fac + &
!          bx*cy/norm12*(dx1s/4.+dy1s)**fac + cx*by/norm21*(dy1s/4.)**fac
!g_scont = ax*by/norm01*(dx4s+dy1s/4.)**fac + bx*ay/norm10*(dx3s+dy1s)**fac + bx*by/norm11*(dx3s+dy1s/4.)**fac + &
!          bx*cy/norm12*dx3s**fac + cx*by/norm21*(dx2s+dy1s/4.)**fac
!h_scont = ax*cy + ax*by/norm01*(dy1s/4.)**fac + bx*ay/norm10*(dx1s/4.+dy1s)**fac + bx*by/norm11*(dx1s/4.+dy1s/4.)**fac + &
!          bx*cy/norm12*(dx1s/4.)**fac + cx*by/norm21*(dx1s+dy1s/4.)**fac
!i_scont = cx*cy + ax*by/norm01*(dx1s+dy1s/4.)**fac + bx*ay/norm10*(dx1s/4.+dy1s)**fac + bx*by/norm11*(dx1s/4.+dy1s/4.)**fac + &
!          bx*cy/norm12*(dx1s/4.)**fac + cx*by/norm21*(dy1s/4.)**fac
!
!opac_p = a_scont*opac_im2jm2 + b_scont*opac_im1jm2 + c_scont*opac_ijm2 + &
!         d_scont*opac_im2jm1 + e_scont*opac_im1jm1 + f_scont*opac_ijm1 + &
!         g_scont*opac_im2j   + h_scont*opac_im1j   + i_scont*opac_ij
   scont_p = a_scont*scont_im2jm2 + b_scont*scont_im1jm2 + c_scont*scont_ijm2 + &
      d_scont*scont_im2jm1 + e_scont*scont_im1jm1 + f_scont*scont_ijm1 + &
      g_scont*scont_im2j   + h_scont*scont_im1j   + i_scont*scont_ij
!
   return
!
!---------------cubic interpolation with free central points------------
!-----------------warning: might cause negative intensities-------------
!
60 continue
!
   dxp = x_p-x_im1
   dxp2 = dxp**2
   dxp3 = dxp**3
   dxi = x_i-x_im1
   dxi2 = dxi**2
   dxim1 = x_im1-x_im2
   dxim12 = dxim1**2
   dx = dxi+dxim1
   dx2 = dx**2
   dx3 = dx**3
!
   dyp = y_p-y_jm1
   dyp2 = dyp**2
   dyp3 = dyp**3
   dyj = y_j-y_jm1
   dyj2 = dyj**2
   dyjm1 = y_jm1-y_jm2
   dyjm12 = dyjm1**2
   dy = dyj+dyjm1
   dy2 = dy**2
   dy3 = dy**3

   ax = ((dxim1-dxi)*dxp3 + 2.d0*(dxi2+dxim12-dxi*dxim1)*dxp2+dxi*(dxi*dxim1-4.d0*dxim12-dxi2)*dxp+2.d0*dxi2*dxim12)/dxim1/dx3
   bx = ((dxi-dxim1)*dxp3 + 2.d0*(dxi*dxim1-dxi2-dxim12)*dxp2)/dxi/dxim1/dx2 + &
      (dxi2-dxim12)*(dxi2+dxim12-dxi*dxim1)*dxp/dxi/dxim1/dx3 + (dxi2+dxim12)/dx2
   cx = ((dxim1-dxi)*dxp3 + 2.d0*(dxi2+dxim12-dxi*dxim1)*dxp2 + dxim1*(4.d0*dxi2+dxim12-dxi*dxim1)*dxp + 2.d0*dxim12*dxi2)/dxi/dx3

   ay = ((dyjm1-dyj)*dyp3 + 2.d0*(dyj2+dyjm12-dyj*dyjm1)*dyp2+dyj*(dyj*dyjm1-4.d0*dyjm12-dyj2)*dyp+2.d0*dyj2*dyjm12)/dyjm1/dy3
   by = ((dyj-dyjm1)*dyp3 + 2.d0*(dyj*dyjm1-dyj2-dyjm12)*dyp2)/dyj/dyjm1/dy2 + &
      (dyj2-dyjm12)*(dyj2+dyjm12-dyj*dyjm1)*dyp/dyj/dyjm1/dy3 + (dyj2+dyjm12)/dy2
   cy = ((dyjm1-dyj)*dyp3 + 2.d0*(dyj2+dyjm12-dyj*dyjm1)*dyp2 + dyjm1*(4.d0*dyj2+dyjm12-dyj*dyjm1)*dyp + 2.d0*dyjm12*dyj2)/dyj/dy3

   a_scont = ay*ax
   b_scont = ay*bx
   c_scont = ay*cx
   d_scont = by*ax
   e_scont = by*bx
   f_scont = by*cx
   g_scont = cy*ax
   h_scont = cy*bx
   i_scont = cy*cx
!
!opac_p = a_scont*opac_im2jm2 + b_scont*opac_im1jm2 + c_scont*opac_ijm2 + &
!         d_scont*opac_im2jm1 + e_scont*opac_im1jm1 + f_scont*opac_ijm1 + &
!         g_scont*opac_im2j   + h_scont*opac_im1j   + i_scont*opac_ij
   scont_p = a_scont*scont_im2jm2 + b_scont*scont_im1jm2 + c_scont*scont_ijm2 + &
      d_scont*scont_im2jm1 + e_scont*scont_im1jm1 + f_scont*scont_ijm1 + &
      g_scont*scont_im2j   + h_scont*scont_im1j   + i_scont*scont_ij

   return
!
!-----------------------biquadratic interpolation-----------------------
!-----------------warning: might cause negative intensities-------------
!
70 continue
!
   dxim1 = x_im1-x_im2
   dxi = x_i-x_im1
   dx = dxim1+dxi
!
   dyjm1 = y_jm1-y_jm2
   dyj = y_j-y_jm1
   dy = dyjm1+dyj
!
   tx = x_p-x_im1
   ty = y_p-y_jm1
!
   a = -dxi/dxim1/dx
   b = (dxi**2-dxim1**2)/dxi/dxim1/dx
   c = dxim1/dxi/dx
   ax = tx**2/dxim1/dx + a*tx
   bx = 1.d0-tx**2/dxi/dxim1 + b*tx
   cx = tx**2/dxi/dx + c*tx
!
   a = -dyj/dyjm1/dy
   b = (dyj**2-dyjm1**2)/dyj/dyjm1/dy
   c = dyjm1/dyj/dy
   ay = ty**2/dyjm1/dy + a*ty
   by = 1.d0-ty**2/dyj/dyjm1 + b*ty
   cy = ty**2/dyj/dy + c*ty
!
   a_scont=ax*ay
   b_scont=bx*ay
   c_scont=cx*ay
   d_scont=ax*by
   e_scont=bx*by
   f_scont=cx*by
   g_scont=ax*cy
   h_scont=bx*cy
   i_scont=cx*cy
!
!opac_p = a_scont*opac_im2jm2 + b_scont*opac_im1jm2 + c_scont*opac_ijm2 + &
!         d_scont*opac_im2jm1 + e_scont*opac_im1jm1 + f_scont*opac_ijm1 + &
!         g_scont*opac_im2j   + h_scont*opac_im1j   + i_scont*opac_ij
   scont_p = a_scont*scont_im2jm2 + b_scont*scont_im1jm2 + c_scont*scont_ijm2 + &
      d_scont*scont_im2jm1 + e_scont*scont_im1jm1 + f_scont*scont_ijm1 + &
      g_scont*scont_im2j   + h_scont*scont_im1j   + i_scont*scont_ij
!
   return
!
!-------------------4-point quadratic bezier interpolation--------------
!-----------(with predefined control point by assigning weights)--------
!
80 continue
!
   wim1=0.7d0
   wi=1.d0-wim1
   wjm1=0.7d0
   wj=1.d0-wjm1
!
   tx = (x_p-x_im1)/(x_i-x_im1)
   ty = (y_p-y_jm1)/(y_j-y_jm1)
!
   ax = (1.d0-tx)**2 + 2.d0*tx*(1.d0-tx)*wim1
   bx = tx**2 + 2.d0*tx*(1.d0-tx)*wi
!
   ay = (1.d0-ty)**2 + 2.d0*ty*(1.d0-ty)*wjm1
   by = ty**2 + 2.d0*ty*(1.d0-ty)*wj
!
   a_scont=0.d0
   b_scont=0.d0
   c_scont=0.d0
   d_scont=0.d0
   e_scont=ax*ay
   f_scont=bx*ay
   g_scont=0.d0
   h_scont=ax*by
   i_scont=bx*by
!
!opac_p = e_scont*opac_im1jm1 + f_scont*opac_ijm1 + h_scont*opac_im1j   + i_scont*opac_ij
   scont_p = e_scont*scont_im1jm1 + f_scont*scont_ijm1 + h_scont*scont_im1j + i_scont*scont_ij
!
   return
!
!-------------------4-point quadratic bezier interpolation--------------
!------------with predefined control point by assigning weights---------
!------------and approximated curvature from negihboring points---------
!
90 continue
!
!wim1=fac, wi=1.d0-fac on each j level, or other way round
!(depending on curvature)
   fac=0.6d0
!
   dxim1=x_im1-x_im2
   dxi=x_i-x_im1
!
   dyjm1=y_jm1-y_jm2
   dyj=y_j-y_jm1
!
   tx = (x_p-x_im1)/dxi
   ty = (y_p-y_jm1)/dyj
!
!level j
   mm=abs((scont_im1j-scont_im2j)/dxim1)
   mp=abs((scont_ij-scont_im1j)/dxi)
   if(mm.gt.mp) then
      wim1=1.d0-fac
      wi=fac
   elseif(mm.lt.mp) then
      wim1=fac
      wi=1.d0-fac
   else
      wim1=0.5d0
      wi=0.5d0
   endif
   axj_scont = (1.d0-tx)**2 + 2.d0*tx*(1.d0-tx)*wim1
   bxj_scont = tx**2 + 2.d0*tx*(1.d0-tx)*wi
!
!level j-1
   mm=abs((scont_im1jm1-scont_im2jm1)/dxim1)
   mp=abs((scont_ijm1-scont_im1jm1)/dxi)
   if(mm.gt.mp) then
      wim1=1.d0-fac
      wi=fac
   elseif(mm.lt.mp) then
      wim1=fac
      wi=1.d0-fac
   else
      wim1=0.5d0
      wi=0.5d0
   endif
   axjm1_scont = (1.d0-tx)**2 + 2.d0*tx*(1.d0-tx)*wim1
   bxjm1_scont = tx**2 + 2.d0*tx*(1.d0-tx)*wi
!
!level j-2
   mm=abs((scont_im1jm2-scont_im2jm2)/dxim1)
   mp=abs((scont_ijm2-scont_im1jm2)/dxi)
   if(mm.gt.mp) then
      wim1=1.d0-fac
      wi=fac
   elseif(mm.lt.mp) then
      wim1=fac
      wi=1.d0-fac
   else
      wim1=0.5d0
      wi=0.5d0
   endif
   axjm2_scont = (1.d0-tx)**2 + 2.d0*tx*(1.d0-tx)*wim1
   bxjm2_scont = tx**2 + 2.d0*tx*(1.d0-tx)*wi
!
!calculate value on each j-level
   scont_jm2 = axjm2_scont*scont_im1jm2 + bxjm2_scont*scont_ijm2
   scont_jm1 = axjm1_scont*scont_im1jm1 + bxjm1_scont*scont_ijm1
   scont_j   = axj_scont*scont_im1j     + bxj_scont*scont_ij
!
!interpolation along y
   ay = (1.d0-ty)**2 + 2.d0*ty*(1.d0-ty)*wjm1
   by = ty**2 + 2.d0*ty*(1.d0-ty)*wj
!
   mm=abs((scont_jm1-scont_jm2)/dyjm1)
   mp=abs((scont_j-scont_jm1)/dyj)
   if(mm.gt.mp) then
      wim1=1.d0-fac
      wi=fac
   elseif(mm.lt.mp) then
      wim1=fac
      wi=1.d0-fac
   else
      wim1=0.5d0
      wi=0.5d0
   endif
   ay_scont = (1.d0-ty)**2 + 2.d0*ty*(1.d0-ty)*wim1
   by_scont = ty**2 + 2.d0*ty*(1.d0-ty)*wi
!
   a_scont=0.d0
   b_scont=0.d0
   c_scont=0.d0
   d_scont=0.d0
   e_scont=ay_scont*axjm1_scont
   f_scont=ay_scont*bxjm1_scont
   g_scont=0.d0
   h_scont=by_scont*axj_scont
   i_scont=by_scont*bxj_scont
!
!
   scont_p = ay_scont*scont_jm1 + by_scont*scont_j
!
   return
!
!-------------------4-point bicubic interpolation-----------------------
!--------------with zero derivatives at all grid points-----------------
!
100 continue
!
   dxi=x_i-x_im1
   dyj=y_j-y_jm1
!
   tx = (x_p-x_im1)/dxi
   ty = (y_p-y_jm1)/dyj
!
   ax = 1.d0 - 3.d0*tx**2 + 2.d0*tx**3
   bx = 1.d0-ax
!
   ay = 1.d0 - 3.d0*ty**2 + 2.d0*ty**3
   by = 1.d0-ay

!
   a_scont=0.d0
   b_scont=0.d0
   c_scont=0.d0
   d_scont=0.d0
   e_scont=ay*ax
   f_scont=ay*bx
   g_scont=0.d0
   h_scont=by*ax
   i_scont=by*bx
!
   scont_p = e_scont*scont_im1jm1 + f_scont*scont_ijm1 + h_scont*scont_im1j + i_scont*scont_ij
!
   return
!---------------monotonic quadratic bezier interpolation----------------
!--------(with derivatives for control point from weighted mean)--------
!----------------and weights for derivative assigned--------------------
!
110 continue
!
   ax = (1.d0-tx)**2
   bx = 2.d0*tx*(1.d0-tx)
   cx = tx**2
!
!calculate control points on each j-level:
!   scontc_jm2, scontc_jm1, scontc_j
!   opacc_jm2, opacc_jm1, opacc_j
   fac=max(wp_interp2d,dxim1/dx)
   fac2=dxim1/dx
   axt = (fac-1.d0)*dxi/2.d0/dxim1
   bxt = ((2.d0-fac)*dxim1 + (1.d0-fac)*dxi)/2.d0/dxim1
   cxt = fac/2.d0
   axt2 = (fac2-1.d0)*dxi/2.d0/dxim1
   bxt2 = ((2.d0-fac2)*dxim1 + (1.d0-fac2)*dxi)/2.d0/dxim1
   cxt2 = fac2/2.d0
!
   opacc_jm2 = opac_im2jm2*axt2 + opac_im1jm2*bxt2 + opac_ijm2*cxt2
   opacc_jm1 = opac_im2jm1*axt2 + opac_im1jm1*bxt2 + opac_ijm1*cxt2
   opacc_j   = opac_im2j*axt2   + opac_im1j*bxt2   + opac_ij*cxt2
!
   scontc_jm2 = scont_im2jm2*axt + scont_im1jm2*bxt + scont_ijm2*cxt
   scontc_jm1 = scont_im2jm1*axt + scont_im1jm1*bxt + scont_ijm1*cxt
   scontc_j   = scont_im2j*axt   + scont_im1j*bxt   + scont_ij*cxt
!
!ensure monotonicity on level j
   call coeffcr1d_mbez(opacc_j, opac_im1j, opac_ij, axt2, bxt2, cxt2, axtj_opac, bxtj_opac, cxtj_opac)
   call coeffcr1d_mbez(scontc_j, scont_im1j, scont_ij, axt, bxt, cxt, axtj_scont, bxtj_scont, cxtj_scont)
   axj_opac = axtj_opac*bx
   bxj_opac = bxtj_opac*bx + ax
   cxj_opac = cxtj_opac*bx + cx
   opac_j = axj_opac*opac_im2j + bxj_opac*opac_im1j + cxj_opac*opac_ij
   axj_scont = axtj_scont*bx
   bxj_scont = bxtj_scont*bx + ax
   cxj_scont = cxtj_scont*bx + cx
   scont_j = axj_scont*scont_im2j + bxj_scont*scont_im1j + cxj_scont*scont_ij
!
!ensure monotonicity on level j-1
   call coeffcr1d_mbez(opacc_jm1, opac_im1jm1, opac_ijm1, axt2, bxt2, cxt2, axtjm1_opac, bxtjm1_opac, cxtjm1_opac)
   call coeffcr1d_mbez(scontc_jm1, scont_im1jm1, scont_ijm1, axt, bxt, cxt, axtjm1_scont, bxtjm1_scont, cxtjm1_scont)
   axjm1_opac = axtjm1_opac*bx
   bxjm1_opac = bxtjm1_opac*bx + ax
   cxjm1_opac = cxtjm1_opac*bx + cx
   opac_jm1 = axjm1_opac*opac_im2jm1 + bxjm1_opac*opac_im1jm1 + cxjm1_opac*opac_ijm1
   axjm1_scont = axtjm1_scont*bx
   bxjm1_scont = bxtjm1_scont*bx + ax
   cxjm1_scont = cxtjm1_scont*bx + cx
   scont_jm1 = axjm1_scont*scont_im2jm1 + bxjm1_scont*scont_im1jm1 + cxjm1_scont*scont_ijm1
!
!ensure monotonicity on level j-2
   call coeffcr1d_mbez(opacc_jm2, opac_im1jm2, opac_ijm2, axt2, bxt2, cxt2, axtjm2_opac, bxtjm2_opac, cxtjm2_opac)
   call coeffcr1d_mbez(scontc_jm2, scont_im1jm2, scont_ijm2, axt, bxt, cxt, axtjm2_scont, bxtjm2_scont, cxtjm2_scont)
   axjm2_opac = axtjm2_opac*bx
   bxjm2_opac = bxtjm2_opac*bx + ax
   cxjm2_opac = cxtjm2_opac*bx + cx
   opac_jm2 = axjm2_opac*opac_im2jm2 + bxjm2_opac*opac_im1jm2 + cxjm2_opac*opac_ijm2
   axjm2_scont = axtjm2_scont*bx
   bxjm2_scont = bxtjm2_scont*bx + ax
   cxjm2_scont = cxtjm2_scont*bx + cx
   scont_jm2 = axjm2_scont*scont_im2jm2 + bxjm2_scont*scont_im1jm2 + cxjm2_scont*scont_ijm2
!
!
!
!calculate control point for interpolation along y
   fac=max(wp_interp2d,dyjm1/dy)
   fac2=dyjm1/dy
   ayt = (fac-1.d0)*dyj/2.d0/dyjm1
   byt = ((2.d0-fac)*dyjm1 + (1.d0-fac)*dyj)/2.d0/dyjm1
   cyt = fac/2.d0
   ayt2 = (fac2-1.d0)*dyj/2.d0/dyjm1
   byt2 = ((2.d0-fac2)*dyjm1 + (1.d0-fac2)*dyj)/2.d0/dyjm1
   cyt2 = fac2/2.d0
!
   opac_c = opac_jm2*ayt2 + opac_jm1*byt2 + opac_j*cyt2
   scont_c = scont_jm2*ayt + scont_jm1*byt + scont_j*cyt
!
!ensure monotonicity
   call coeffcr1d_mbez(opac_c, opac_jm1, opac_j, ayt2, byt2, cyt2, ayt_opac, byt_opac, cyt_opac)
   call coeffcr1d_mbez(scont_c, scont_jm1, scont_j, ayt, byt, cyt, ayt_scont, byt_scont, cyt_scont)
!
   ay = (1.d0-ty)**2
   by = 2.d0*ty*(1.d0-ty)
   cy = ty**2
!
   ay_opac = ayt_opac*by
   by_opac = byt_opac*by + ay
   cy_opac = cyt_opac*by + cy
   opac_p = ay_opac*opac_jm2 + by_opac*opac_jm1 + cy_opac*opac_j
   ay_scont = ayt_scont*by
   by_scont = byt_scont*by + ay
   cy_scont = cyt_scont*by + cy
   scont_p = ay_scont*scont_jm2 + by_scont*scont_jm1 + cy_scont*scont_j
!
   a_scont = ay_scont*axjm2_scont
   b_scont = ay_scont*bxjm2_scont
   c_scont = ay_scont*cxjm2_scont
   d_scont = by_scont*axjm1_scont
   e_scont = by_scont*bxjm1_scont
   f_scont = by_scont*cxjm1_scont
   g_scont = cy_scont*axj_scont
   h_scont = cy_scont*bxj_scont
   i_scont = cy_scont*cxj_scont
!
   return
!
!
end subroutine coeff2d_contd
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine coeff2d_contu_lin(opac_im1jm1, opac_ijm1, opac_im1j, opac_ij, &
   scont_im1jm1, scont_ijm1, scont_im1j, scont_ij, &
   int_im1jm1, int_ijm1, int_im1j, int_ij, &
   x_im1, x_i, y_jm1, y_j, x_p, y_p, &
   a_scont, b_scont, c_scont, d_scont, &
   a_inten, b_inten, c_inten, d_inten, &
   opac_p, scont_p, int_p)
!
!         interpolates opacity, continuum source function and intensity
!               values given on a 2d grid onto point x_p, y_p
!
!on input (f_* stands for opac_*, scont_* and int_*, respectivly):
!
! y_j      f_im1j--------------f_ij
!  |          |                  |
!  |          |                  |
!  |          |         x        |
!  |          |     (x_p,y_p)    |
!  |          |                  |
!y_jm1    f_im1jm1------------f_ijm1
!  |
!  --------x_im2-----------x_im1---------------x_i
!
!        x_p, y_p: coordinates of point onto which shall be interpolated
!
!on output:
!   1. interpolation coefficients for source function and intensity
!      (required for ALO calculations):
!         a_scont, b_scont, c_scont, d_scont
!         a_inten, b_inten, c_inten, d_inten
!
!      such that:
!         f_p = a*f_im1jm1 + b*f_ijm1 + c*f_im1j + d*f_ij
!
!   2. interpolated values at point p: opac_p, scont_p, int_p
!
   use prog_type
!
   implicit none
!
! ... arguments
   real(dp), intent(in) :: opac_im1jm1, opac_ijm1, opac_im1j,   opac_ij, &
      scont_im1jm1, scont_ijm1, scont_im1j, scont_ij, &
      int_im1jm1, int_ijm1, int_im1j, int_ij, &
      x_im1, x_i, y_jm1, y_j, x_p, y_p
   real(dp), intent(out) :: a_scont, b_scont, c_scont, d_scont, &
      a_inten, b_inten, c_inten, d_inten, &
      opac_p, int_p, scont_p
!
! ... local scalars
   real(dp) :: dxi, tx, dyj, ty, rdxdy
!
!
!define deltax, deltay
   dxi = x_i-x_im1
   dyj = y_j-y_jm1
!
!define deltax, deltay-ratios
   tx = (x_p-x_im1)/dxi
   ty = (y_p-y_jm1)/dyj
!
!-------------------------bilinear interpolation------------------------
!
   rdxdy=tx*ty
!
   a_scont=1.d0-tx-ty+rdxdy
   b_scont=tx-rdxdy
   c_scont=ty-rdxdy
   d_scont=rdxdy
!
   a_inten=a_scont
   b_inten=b_scont
   c_inten=c_scont
   d_inten=d_scont
!
   opac_p = a_scont*opac_im1jm1 + b_scont*opac_ijm1 + c_scont*opac_im1j + d_scont*opac_ij
   scont_p = a_scont*scont_im1jm1 + b_scont*scont_ijm1 + c_scont*scont_im1j + d_scont*scont_ij
   int_p = a_scont*int_im1jm1 + b_scont*int_ijm1 + c_scont*int_im1j + d_scont*int_ij
   return
!
end subroutine coeff2d_contu_lin
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine coeff2d_contd_lin(opac_im1jm1, opac_ijm1, opac_im1j,   opac_ij, &
   scont_im1jm1, scont_ijm1, scont_im1j,   scont_ij, &
   x_im1, x_i, y_jm1, y_j, x_p, y_p, &
   a_scont, b_scont, c_scont, d_scont, &
   opac_p, scont_p)
!
!            interpolates opacity and continuum source function
!               values given on a 2d grid onto point x_p, y_p
!
!on input (f_* stands for opac_* and scont_* respectivly):
!
! y_j      f_im1j--------------f_ij
!  |       |                  |
!  |       |                  |
!  |       |         x        |
!  |       |     (x_p,y_p)    |
!  |       |                  |
!y_jm1  f_im1jm1------------f_ijm1
!  |
!  ------x_im1---------------x_i
!
!        x_p, y_p: coordinates of point onto which shall be interpolated
!
!on output:
!   1. interpolation coefficients for source function
!      (required for ALO calculations):
!         a_scont, b_scont, c_scont, d_scont
!
!      such that:
!         f_p = a*f_im1jm1 + b*f_ijm1 + c*f_im1j + d*f_ij
!
!   2. interpolated values at point p: opac_p, scont_p
!
   use prog_type
!
   implicit none
!
! ... arguments
   real(dp), intent(in) :: opac_im1jm1, opac_ijm1, opac_im1j, opac_ij, &
      scont_im1jm1, scont_ijm1, scont_im1j, scont_ij, &
      x_im1, x_i, y_jm1, y_j, x_p, y_p
   real(dp), intent(out) :: a_scont, b_scont, c_scont, d_scont, opac_p, scont_p
!
! ... local scalars
   real(dp) :: dxi, tx, dyj, ty, rdxdy
!
!define deltax, deltay
   dxi = x_i-x_im1
   dyj = y_j-y_jm1
!
!define deltax, deltay-ratios
   tx = (x_p-x_im1)/dxi
   ty = (y_p-y_jm1)/dyj
!
!-------------------------bilinear interpolation------------------------
!
   rdxdy=tx*ty
!
   a_scont=1.d0-tx-ty+rdxdy
   b_scont=tx-rdxdy
   c_scont=ty-rdxdy
   d_scont=rdxdy
!
   opac_p = a_scont*opac_im1jm1 + b_scont*opac_ijm1 + c_scont*opac_im1j + d_scont*opac_ij
   scont_p = a_scont*scont_im1jm1 + b_scont*scont_ijm1 + c_scont*scont_im1j + d_scont*scont_ij
!
!
end subroutine coeff2d_contd_lin
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine coeff3d_contu(opac_im2jm2km2, opac_im1jm2km2, opac_ijm2km2, &
   opac_im2jm1km2, opac_im1jm1km2, opac_ijm1km2, &
   opac_im2jkm2,   opac_im1jkm2,   opac_ijkm2, &
   opac_im2jm2km1, opac_im1jm2km1, opac_ijm2km1, &
   opac_im2jm1km1, opac_im1jm1km1, opac_ijm1km1, &
   opac_im2jkm1,   opac_im1jkm1,   opac_ijkm1, &
   opac_im2jm2k, opac_im1jm2k,     opac_ijm2k, &
   opac_im2jm1k, opac_im1jm1k,     opac_ijm1k, &
   opac_im2jk,   opac_im1jk,       opac_ijk, &
   scont_im2jm2km2, scont_im1jm2km2, scont_ijm2km2, &
   scont_im2jm1km2, scont_im1jm1km2, scont_ijm1km2, &
   scont_im2jkm2,   scont_im1jkm2,   scont_ijkm2, &
   scont_im2jm2km1, scont_im1jm2km1, scont_ijm2km1, &
   scont_im2jm1km1, scont_im1jm1km1, scont_ijm1km1, &
   scont_im2jkm1,   scont_im1jkm1,   scont_ijkm1, &
   scont_im2jm2k, scont_im1jm2k,     scont_ijm2k, &
   scont_im2jm1k, scont_im1jm1k,     scont_ijm1k, &
   scont_im2jk,   scont_im1jk,       scont_ijk, &
   x_im2, x_im1, x_i, y_jm2, y_jm1, y_j, z_km2, z_km1, z_k, x_p, y_p, z_p, &
   c01, c02, c03, c04, c05, c06, c07, c08, c09, &
   c10, c11, c12, c13, c14, c15, c16, c17, c18, &
   c19, c20, c21, c22, c23, c24, c25, c26, c27, &
   opac_p, scont_p)
!
!         interpolates opacity and continuum source function
!           given on a 3d grid onto point x_p, y_p, z_p
!
!on input:
!
!                 f_im2jk---------------f_im1jk--------------f_ijk
!                    /|                    /|                 /|
!                   / |                   / |                / |
!                  /  |                  /  |               /  |
!                 /   |                 /   |              /   |
!                /    |                /    |             /    |
!               / f_im2jkm1-----------/-f_im1jkm1--------/--f_ijkm1
!              /     /|              /     /|           /     /|
!       f_im2jm1k--------------f_im1jm1k-----------f_ijm1k   / |
!            /|    /  |            /|    /  |         /|    /  |
!           / |   /   |           / |   /   |        / |   /   |
!          /  |  /    |          /  |  /    |       /  |  /    |
!         /   | / f_im2jkm2-----/---|-/-f_im1jkm2--/---|-/--f_ijkm2
!        /    |/     /         /    |/     /      /    |/     /
!       /f_im2jm1km1----------/-f_im1jm1km1------/--f_ijm1km1/
!      /     /|    /         /     /|    /      /     /|    /
!f_im2jm2k--------------f_im1jm2k------------f_ijm2k / |   /
!     |    /  |  /          |    /  |  /       |    /  |  /
!     |   /   | /           |   /   | /        |   /   | /
!     |  /    |/            |  /    |/         |  /    |/
!     | / f_im2jm1km2-------|-/-f_im1jm1km2----|-/--f_ijm1km2
!     |/     /              |/     /           |/     /
!f_im2jm2km1------------f_im1jm2km1---------f_ijm2km1/
!     |    /                |    /             |    /
!     |   /                 |   /              |   /
!     |  /                  |  /               |  /
!     | /                   | /                | /
!     |/                    |/                 |/
!f_im2jm2km2----------f_im1jm2km2---------f_ijm2km2
!
!        x_p, y_p, z_p: coordinates of point onto which shall be interpolated
!
!on output:
!   1. interpolation coefficients for source function and intensity
!      (required for ALO calculations):
!         c01, c02, c03, c04, c05, c06, c07, c08, c09,
!         c10, c11, c12, c13, c14, c15, c16, c17, c18,
!         c19, c20, c21, c22, c23, c24, c25, c26, c27
!
!   such that
!      f_p = c01*f_im2jm2km2 + c02*f_im1jm2km2 + c03*f_ijm2km2 +
!            c04*f_im2jm1km2 + c05*f_im1jm1km2 + c06*f_ijm1km2 +
!            c07*f_im2jkm2   + c08*f_im1jkm2   + c09*f_ijkm2 +
!
!            c10*f_im2jm2km1 + c11*f_im1jm2km1 + c12*f_ijm2km1 +
!            c13*f_im2jm1km1 + c14*f_im1jm1km1 + c15*f_ijm1km1 +
!            c16*f_im2jkm1   + c17*f_im1jkm1   + c18*f_ijkm1 +
!
!            c19*f_im2jm2k + c20*f_im1jm2k + c21*f_ijm2k +
!            c22*f_im2jm1k + c23*f_im1jm1k + c24*f_ijm1k +
!            c25*f_im2jk   + c26*f_im1jk   + c27*f_ijk

!         f_p = a*f_im2jm2 + b*f_im1jm2 + c*f_ijm2 +
!               d*f_im2jm1 + e*f_im1jm1 + f*f_ijm1 +
!               g*f_im2j   + h*f_im1j   + i*f_ij
!
!   2. interpolated values at point p: opac_p, scont_p
!
   use prog_type
!
   implicit none
!
! ... argments
   real(dp), intent(in) :: opac_im2jm2km2, opac_im1jm2km2, opac_ijm2km2, &
      opac_im2jm1km2, opac_im1jm1km2, opac_ijm1km2, &
      opac_im2jkm2,   opac_im1jkm2,   opac_ijkm2, &
      opac_im2jm2km1, opac_im1jm2km1, opac_ijm2km1, &
      opac_im2jm1km1, opac_im1jm1km1, opac_ijm1km1, &
      opac_im2jkm1,   opac_im1jkm1,   opac_ijkm1, &
      opac_im2jm2k, opac_im1jm2k,     opac_ijm2k, &
      opac_im2jm1k, opac_im1jm1k,     opac_ijm1k, &
      opac_im2jk,   opac_im1jk,       opac_ijk, &
      scont_im2jm2km2, scont_im1jm2km2, scont_ijm2km2, &
      scont_im2jm1km2, scont_im1jm1km2, scont_ijm1km2, &
      scont_im2jkm2,   scont_im1jkm2,   scont_ijkm2, &
      scont_im2jm2km1, scont_im1jm2km1, scont_ijm2km1, &
      scont_im2jm1km1, scont_im1jm1km1, scont_ijm1km1, &
      scont_im2jkm1,   scont_im1jkm1,   scont_ijkm1, &
      scont_im2jm2k, scont_im1jm2k,     scont_ijm2k, &
      scont_im2jm1k, scont_im1jm1k,     scont_ijm1k, &
      scont_im2jk,   scont_im1jk,       scont_ijk, &
      x_im2, x_im1, x_i, y_jm2, y_jm1, y_j, z_km2, z_km1, z_k, x_p, y_p, z_p
   real(dp), intent(out) :: c01, c02, c03, c04, c05, c06, c07, c08, c09, &
      c10, c11, c12, c13, c14, c15, c16, c17, c18, &
      c19, c20, c21, c22, c23, c24, c25, c26, c27, &
      opac_p, scont_p
!
! ... local scalars
   real(dp) :: dxim1, dxi, dx, dyjm1, dyj, dy, dzkm1, dzk, dz, &
      tx, ty, tz, ax, bx, cx, ay, by, cy, az, bz, cz, &
      axt, bxt, cxt, ayt, byt, cyt, azt, bzt, czt
   real(dp) :: scontc_jm2km2, scontc_jm1km2, scontc_jkm2, &
      scontc_jm2km1, scontc_jm1km1, scontc_jkm1, &
      scontc_jm2k,   scontc_jm1k,   scontc_jk, &
      scontc_km2, scontc_km1, scontc_k, scontc, &
      scont_jm2km2, scont_jm1km2, scont_jkm2, &
      scont_jm2km1, scont_jm1km1, scont_jkm1, &
      scont_jm2k,   scont_jm1k,   scont_jk, &
      scontc_jm2, scontc_jm1, scontc_j, &
      scont_km2, scont_km1, scont_k, &
      axjm2km2_scont, axjm1km2_scont, axjkm2_scont, &
      axjm2km1_scont, axjm1km1_scont, axjkm1_scont, &
      axjm2k_scont, axjm1k_scont, axjk_scont, &
      bxjm2km2_scont, bxjm1km2_scont, bxjkm2_scont, &
      bxjm2km1_scont, bxjm1km1_scont, bxjkm1_scont, &
      bxjm2k_scont, bxjm1k_scont, bxjk_scont, &
      cxjm2km2_scont, cxjm1km2_scont, cxjkm2_scont, &
      cxjm2km1_scont, cxjm1km1_scont, cxjkm1_scont, &
      cxjm2k_scont, cxjm1k_scont, cxjk_scont, &
      aykm2_scont, aykm1_scont, ayk_scont, &
      bykm2_scont, bykm1_scont, byk_scont, &
      cykm2_scont, cykm1_scont, cyk_scont, &
      az_scont, bz_scont, cz_scont
   real(dp) :: opacc_jm2km2, opacc_jm1km2, opacc_jkm2, &
      opacc_jm2km1, opacc_jm1km1, opacc_jkm1, &
      opacc_jm2k,   opacc_jm1k,   opacc_jk, &
      opacc_km2, opacc_km1, opacc_k, opacc, &
      opac_jm2km2, opac_jm1km2, opac_jkm2, &
      opac_jm2km1, opac_jm1km1, opac_jkm1, &
      opac_jm2k,   opac_jm1k,   opac_jk, &
      opacc_jm2, opacc_jm1, opacc_j, &
      opac_km2, opac_km1, opac_k, &
      axjm2km2_opac, axjm1km2_opac, axjkm2_opac, &
      axjm2km1_opac, axjm1km1_opac, axjkm1_opac, &
      axjm2k_opac, axjm1k_opac, axjk_opac, &
      bxjm2km2_opac, bxjm1km2_opac, bxjkm2_opac, &
      bxjm2km1_opac, bxjm1km1_opac, bxjkm1_opac, &
      bxjm2k_opac, bxjm1k_opac, bxjk_opac, &
      cxjm2km2_opac, cxjm1km2_opac, cxjkm2_opac, &
      cxjm2km1_opac, cxjm1km1_opac, cxjkm1_opac, &
      cxjm2k_opac, cxjm1k_opac, cxjk_opac, &
      aykm2_opac, aykm1_opac, ayk_opac, &
      bykm2_opac, bykm1_opac, byk_opac, &
      cykm2_opac, cykm1_opac, cyk_opac, &
      az_opac, bz_opac, cz_opac
!
!define deltax, deltay, deltaz
   dxim1 = x_im1-x_im2
   dxi = x_i-x_im1
   dx = dxim1+dxi
   tx = (x_p-x_im1)/dxi
!
   dyjm1 = y_jm1-y_jm2
   dyj = y_j-y_jm1
   dy = dyjm1+dyj
   ty = (y_p-y_jm1)/dyj
!
   dzkm1 = z_km1-z_km2
   dzk = z_k-z_km1
   dz = dzkm1+dzk
   tz = (z_p-z_km1)/dzk
!
!--------------------------trilinear interpolation----------------------
!
   ax=1.d0-tx
   ay=1.d0-ty
   az=1.d0-tz

   c01 = 0.d0
   c02 = 0.d0
   c03 = 0.d0
   c04 = 0.d0
   c05 = 0.d0
   c06 = 0.d0
   c07 = 0.d0
   c08 = 0.d0
   c09 = 0.d0
   c10 = 0.d0
   c11 = 0.d0
   c12 = 0.d0
   c13 = 0.d0
   c14 = ax*ay*az
   c15 = tx*ay*az
   c16 = 0.d0
   c17 = ax*ty*az
   c18 = tx*ty*az
   c19 = 0.d0
   c20 = 0.d0
   c21 = 0.d0
   c22 = 0.d0
   c23 = ay*ax*tz
   c24 = tx*ay*tz
   c25 = 0.d0
   c26 = ax*ty*tz
   c27 = tx*ty*tz

   opac_p = c14*opac_im1jm1km1 + c15*opac_ijm1km1 + c17*opac_im1jkm1   + c18*opac_ijkm1 + &
      c23*opac_im1jm1k + c24*opac_ijm1k + c26*opac_im1jk   + c27*opac_ijk

   scont_p = c14*scont_im1jm1km1 + c15*scont_ijm1km1 + c17*scont_im1jkm1   + c18*scont_ijkm1 + &
      c23*scont_im1jm1k + c24*scont_ijm1k + c26*scont_im1jk   + c27*scont_ijk
!
!return
!
   return

!NOTE: ONLY LINEAR INTERPOLATION IS ALLOWED BECAUSE OF ALO CALCULATIONS!!!!!

!
!***********************************************************************
!---------------monotonic quadratic bezier interpolation----------------
!--------(with derivatives for control point from weighted mean)--------
!***********************************************************************
!
   ax = (1.d0-tx)**2
   bx = 2.d0*tx*(1.d0-tx)
   cx = tx**2
!
   ay = (1.d0-ty)**2
   by = 2.d0*ty*(1.d0-ty)
   cy = ty**2
!
   az = (1.d0-tz)**2
   bz = 2.d0*tz*(1.d0-tz)
   cz = tz**2
!
   axt = -dxi**2/2.d0/dxim1/dx
   bxt = dx/2.d0/dxim1
   cxt = dxim1/2.d0/dx
!
   ayt = -dyj**2/2.d0/dyjm1/dy
   byt = dy/2.d0/dyjm1
   cyt = dyjm1/2.d0/dy
!
   azt = -dzk**2/2.d0/dzkm1/dz
   bzt = dz/2.d0/dzkm1
   czt = dzkm1/2.d0/dz
!
!-----------------------interpolation along x---------------------------
!----------------------------level k-2----------------------------------
!
!calculate control points on each j-level
!opacc_jm2km2 = opac_im2jm2km2*axt + opac_im1jm2km2*bxt + opac_ijm2km2*cxt
!opacc_jm1km2 = opac_im2jm1km2*axt + opac_im1jm1km2*bxt + opac_ijm1km2*cxt
!opacc_jkm2   = opac_im2jkm2*axt   + opac_im1jkm2*bxt   + opac_ijkm2*cxt
   scontc_jm2km2 = scont_im2jm2km2*axt + scont_im1jm2km2*bxt + scont_ijm2km2*cxt
   scontc_jm1km2 = scont_im2jm1km2*axt + scont_im1jm1km2*bxt + scont_ijm1km2*cxt
   scontc_jkm2   = scont_im2jkm2*axt   + scont_im1jkm2*bxt   + scont_ijkm2*cxt
!
!ensure monotonicity on level j
!call coeffcr1d_mbez(opacc_j, opac_im1jkm2, opac_ijkm2, axt, bxt, cxt, axjkm2_opac, bxjkm2_opac, cxjkm2_opac)
!axjkm2_opac = axjkm2_opac*bx
!bxjkm2_opac = bxjkm2_opac*bx + ax
!cxjkm2_opac = cxjkm2_opac*bx + cx
!opac_jkm2 = axjkm2_opac*opac_im2jkm2 + bxjkm2_opac*opac_im1jkm2 + cxjkm2_opac*opac_ijkm2
   call coeffcr1d_mbez(scontc_j, scont_im1jkm2, scont_ijkm2, axt, bxt, cxt, axjkm2_scont, bxjkm2_scont, cxjkm2_scont)
   axjkm2_scont = axjkm2_scont*bx
   bxjkm2_scont = bxjkm2_scont*bx + ax
   cxjkm2_scont = cxjkm2_scont*bx + cx
   scont_jkm2 = axjkm2_scont*scont_im2jkm2 + bxjkm2_scont*scont_im1jkm2 + cxjkm2_scont*scont_ijkm2
!
!ensure monotonicity on level j-1
!call coeffcr1d_mbez(opacc_jm1km2, opac_im1jm1km2, opac_ijm1km2, axt, bxt, cxt, axjm1km2_opac, bxjm1km2_opac, cxjm1km2_opac)
!axjm1km2_opac = axjm1km2_opac*bx
!bxjm1km2_opac = bxjm1km2_opac*bx + ax
!cxjm1km2_opac = cxjm1km2_opac*bx + cx
!opac_jm1km2 = axjm1km2_opac*opac_im2jm1km2 + bxjm1km2_opac*opac_im1jm1km2 + cxjm1km2_opac*opac_ijm1km2
   call coeffcr1d_mbez(scontc_jm1km2, scont_im1jm1km2, scont_ijm1km2, axt, bxt, cxt, axjm1km2_scont, bxjm1km2_scont, cxjm1km2_scont)
   axjm1km2_scont = axjm1km2_scont*bx
   bxjm1km2_scont = bxjm1km2_scont*bx + ax
   cxjm1km2_scont = cxjm1km2_scont*bx + cx
   scont_jm1km2 = axjm1km2_scont*scont_im2jm1km2 + bxjm1km2_scont*scont_im1jm1km2 + cxjm1km2_scont*scont_ijm1km2
!
!ensure monotonicity on level j-2
!call coeffcr1d_mbez(opacc_jm2km2, opac_im1jm2km2, opac_ijm2km2, axt, bxt, cxt, axjm2km2_opac, bxjm2km2_opac, cxjm2km2_opac)
!axjm2km2_opac = axjm2km2_opac*bx
!bxjm2km2_opac = bxjm2km2_opac*bx + ax
!cxjm2km2_opac = cxjm2km2_opac*bx + cx
!opac_jm2km2 = axjm2km2_opac*opac_im2jm2km2 + bxjm2km2_opac*opac_im1jm2km2 + cxjm2km2_opac*opac_ijm2km2
   call coeffcr1d_mbez(scontc_jm2km2, scont_im1jm2km2, scont_ijm2km2, axt, bxt, cxt, axjm2km2_scont, bxjm2km2_scont, cxjm2km2_scont)
   axjm2km2_scont = axjm2km2_scont*bx
   bxjm2km2_scont = bxjm2km2_scont*bx + ax
   cxjm2km2_scont = cxjm2km2_scont*bx + cx
   scont_jm2km2 = axjm2km2_scont*scont_im2jm2km2 + bxjm2km2_scont*scont_im1jm2km2 + cxjm2km2_scont*scont_ijm2km2
!
!---------------------------level k-1-----------------------------------
!
!calculate control points on each j-level
!opacc_jm2km1 = opac_im2jm2km1*axt + opac_im1jm2km1*bxt + opac_ijm2km1*cxt
!opacc_jm1km1 = opac_im2jm1km1*axt + opac_im1jm1km1*bxt + opac_ijm1km1*cxt
!opacc_jkm1   = opac_im2jkm1*axt   + opac_im1jkm1*bxt   + opac_ijkm1*cxt
   scontc_jm2km1 = scont_im2jm2km1*axt + scont_im1jm2km1*bxt + scont_ijm2km1*cxt
   scontc_jm1km1 = scont_im2jm1km1*axt + scont_im1jm1km1*bxt + scont_ijm1km1*cxt
   scontc_jkm1   = scont_im2jkm1*axt   + scont_im1jkm1*bxt   + scont_ijkm1*cxt
!
!ensure monotonicity on level j
!call coeffcr1d_mbez(opacc_j, opac_im1jkm1, opac_ijkm1, axt, bxt, cxt, axjkm1_opac, bxjkm1_opac, cxjkm1_opac)
!axjkm1_opac = axjkm1_opac*bx
!bxjkm1_opac = bxjkm1_opac*bx + ax
!cxjkm1_opac = cxjkm1_opac*bx + cx
!opac_jkm1 = axjkm1_opac*opac_im2jkm1 + bxjkm1_opac*opac_im1jkm1 + cxjkm1_opac*opac_ijkm1
   call coeffcr1d_mbez(scontc_j, scont_im1jkm1, scont_ijkm1, axt, bxt, cxt, axjkm1_scont, bxjkm1_scont, cxjkm1_scont)
   axjkm1_scont = axjkm1_scont*bx
   bxjkm1_scont = bxjkm1_scont*bx + ax
   cxjkm1_scont = cxjkm1_scont*bx + cx
   scont_jkm1 = axjkm1_scont*scont_im2jkm1 + bxjkm1_scont*scont_im1jkm1 + cxjkm1_scont*scont_ijkm1
!
!ensure monotonicity on level j-1
!call coeffcr1d_mbez(opacc_jm1km1, opac_im1jm1km1, opac_ijm1km1, axt, bxt, cxt, axjm1km1_opac, bxjm1km1_opac, cxjm1km1_opac)
!axjm1km1_opac = axjm1km1_opac*bx
!bxjm1km1_opac = bxjm1km1_opac*bx + ax
!cxjm1km1_opac = cxjm1km1_opac*bx + cx
!opac_jm1km1 = axjm1km1_opac*opac_im2jm1km1 + bxjm1km1_opac*opac_im1jm1km1 + cxjm1km1_opac*opac_ijm1km1
   call coeffcr1d_mbez(scontc_jm1km1, scont_im1jm1km1, scont_ijm1km1, axt, bxt, cxt, axjm1km1_scont, bxjm1km1_scont, cxjm1km1_scont)
   axjm1km1_scont = axjm1km1_scont*bx
   bxjm1km1_scont = bxjm1km1_scont*bx + ax
   cxjm1km1_scont = cxjm1km1_scont*bx + cx
   scont_jm1km1 = axjm1km1_scont*scont_im2jm1km1 + bxjm1km1_scont*scont_im1jm1km1 + cxjm1km1_scont*scont_ijm1km1
!
!ensure monotonicity on level j-2
!call coeffcr1d_mbez(opacc_jm2km1, opac_im1jm2km1, opac_ijm2km1, axt, bxt, cxt, axjm2km1_opac, bxjm2km1_opac, cxjm2km1_opac)
!axjm2km1_opac = axjm2km1_opac*bx
!bxjm2km1_opac = bxjm2km1_opac*bx + ax
!cxjm2km1_opac = cxjm2km1_opac*bx + cx
!opac_jm2km1 = axjm2km1_opac*opac_im2jm2km1 + bxjm2km1_opac*opac_im1jm2km1 + cxjm2km1_opac*opac_ijm2km1
   call coeffcr1d_mbez(scontc_jm2km1, scont_im1jm2km1, scont_ijm2km1, axt, bxt, cxt, axjm2km1_scont, bxjm2km1_scont, cxjm2km1_scont)
   axjm2km1_scont = axjm2km1_scont*bx
   bxjm2km1_scont = bxjm2km1_scont*bx + ax
   cxjm2km1_scont = cxjm2km1_scont*bx + cx
   scont_jm2km1 = axjm2km1_scont*scont_im2jm2km1 + bxjm2km1_scont*scont_im1jm2km1 + cxjm2km1_scont*scont_ijm2km1
!
!----------------------------level k------------------------------------
!
!calculate control points on each j-level
!opacc_jm2k = opac_im2jm2k*axt + opac_im1jm2k*bxt + opac_ijm2k*cxt
!opacc_jm1k = opac_im2jm1k*axt + opac_im1jm1k*bxt + opac_ijm1k*cxt
!opacc_jk   = opac_im2jk*axt   + opac_im1jk*bxt   + opac_ijk*cxt
   scontc_jm2k = scont_im2jm2k*axt + scont_im1jm2k*bxt + scont_ijm2k*cxt
   scontc_jm1k = scont_im2jm1k*axt + scont_im1jm1k*bxt + scont_ijm1k*cxt
   scontc_jk   = scont_im2jk*axt   + scont_im1jk*bxt   + scont_ijk*cxt
!
!ensure monotonicity on level j
!call coeffcr1d_mbez(opacc_j, opac_im1jk, opac_ijk, axt, bxt, cxt, axjk_opac, bxjk_opac, cxjk_opac)
!axjk_opac = axjk_opac*bx
!bxjk_opac = bxjk_opac*bx + ax
!cxjk_opac = cxjk_opac*bx + cx
!opac_jk = axjk_opac*opac_im2jk + bxjk_opac*opac_im1jk + cxjk_opac*opac_ijk
   call coeffcr1d_mbez(scontc_j, scont_im1jk, scont_ijk, axt, bxt, cxt, axjk_scont, bxjk_scont, cxjk_scont)
   axjk_scont = axjk_scont*bx
   bxjk_scont = bxjk_scont*bx + ax
   cxjk_scont = cxjk_scont*bx + cx
   scont_jk = axjk_scont*scont_im2jk + bxjk_scont*scont_im1jk + cxjk_scont*scont_ijk
!
!ensure monotonicity on level j-1
!call coeffcr1d_mbez(opacc_jm1k, opac_im1jm1k, opac_ijm1k, axt, bxt, cxt, axjm1k_opac, bxjm1k_opac, cxjm1k_opac)
!axjm1k_opac = axjm1k_opac*bx
!bxjm1k_opac = bxjm1k_opac*bx + ax
!cxjm1k_opac = cxjm1k_opac*bx + cx
!opac_jm1k = axjm1k_opac*opac_im2jm1k + bxjm1k_opac*opac_im1jm1k + cxjm1k_opac*opac_ijm1k
   call coeffcr1d_mbez(scontc_jm1k, scont_im1jm1k, scont_ijm1k, axt, bxt, cxt, axjm1k_scont, bxjm1k_scont, cxjm1k_scont)
   axjm1k_scont = axjm1k_scont*bx
   bxjm1k_scont = bxjm1k_scont*bx + ax
   cxjm1k_scont = cxjm1k_scont*bx + cx
   scont_jm1k = axjm1k_scont*scont_im2jm1k + bxjm1k_scont*scont_im1jm1k + cxjm1k_scont*scont_ijm1k
!
!ensure monotonicity on level j-2
!call coeffcr1d_mbez(opacc_jm2k, opac_im1jm2k, opac_ijm2k, axt, bxt, cxt, axjm2k_opac, bxjm2k_opac, cxjm2k_opac)
!axjm2k_opac = axjm2k_opac*bx
!bxjm2k_opac = bxjm2k_opac*bx + ax
!cxjm2k_opac = cxjm2k_opac*bx + cx
!opac_jm2k = axjm2k_opac*opac_im2jm2k + bxjm2k_opac*opac_im1jm2k + cxjm2k_opac*opac_ijm2k
   call coeffcr1d_mbez(scontc_jm2k, scont_im1jm2k, scont_ijm2k, axt, bxt, cxt, axjm2k_scont, bxjm2k_scont, cxjm2k_scont)
   axjm2k_scont = axjm2k_scont*bx
   bxjm2k_scont = bxjm2k_scont*bx + ax
   cxjm2k_scont = cxjm2k_scont*bx + cx
   scont_jm2k = axjm2k_scont*scont_im2jm2k + bxjm2k_scont*scont_im1jm2k + cxjm2k_scont*scont_ijm2k
!
!---------------------------interpolation along y-----------------------
!
!calculate control points on each k-level
!opacc_km2 = opac_jm2km2*ayt + opac_jm1km2*byt + opac_jkm2*cyt
!opacc_km1 = opac_jm2km1*ayt + opac_jm1km1*byt + opac_jkm1*cyt
!opacc_k   = opac_jm2k*ayt   + opac_jm1k*byt   + opac_jk*cyt
   scontc_km2 = scont_jm2km2*ayt + scont_jm1km2*byt + scont_jkm2*cyt
   scontc_km1 = scont_jm2km1*ayt + scont_jm1km1*byt + scont_jkm1*cyt
   scontc_k   = scont_jm2k*ayt   + scont_jm1k*byt   + scont_jk*cyt
!
!ensure monotonicity on level k
!call coeffcr1d_mbez(opacc_k, opac_jm1k, opac_jk, ayt, byt, cyt, ayk_opac, byk_opac, cyk_opac)
!ayk_opac = ayk_opac*by
!byk_opac = byk_opac*by + ay
!cyk_opac = cyk_opac*by + cy
!opac_k = ayk_opac*opac_jm2k + byk_opac*opac_jm1k + cyk_opac*opac_jk
   call coeffcr1d_mbez(scontc_k, scont_jm1k, scont_jk, ayt, byt, cyt, ayk_scont, byk_scont, cyk_scont)
   ayk_scont = ayk_scont*by
   byk_scont = byk_scont*by + ay
   cyk_scont = cyk_scont*by + cy
   scont_k = ayk_scont*scont_jm2k + byk_scont*scont_jm1k + cyk_scont*scont_jk
!
!ensure monotonicity on level k-1
!call coeffcr1d_mbez(opacc_km1, opac_jm1km1, opac_jkm1, ayt, byt, cyt, aykm1_opac, bykm1_opac, cykm1_opac)
!aykm1_opac = aykm1_opac*by
!bykm1_opac = bykm1_opac*by + ay
!cykm1_opac = cykm1_opac*by + cy
!opac_km1 = aykm1_opac*opac_jm2km1 + bykm1_opac*opac_jm1km1 + cykm1_opac*opac_jkm1
   call coeffcr1d_mbez(scontc_km1, scont_jm1km1, scont_jkm1, ayt, byt, cyt, aykm1_scont, bykm1_scont, cykm1_scont)
   aykm1_scont = aykm1_scont*by
   bykm1_scont = bykm1_scont*by + ay
   cykm1_scont = cykm1_scont*by + cy
   scont_km1 = aykm1_scont*scont_jm2km1 + bykm1_scont*scont_jm1km1 + cykm1_scont*scont_jkm1
!
!ensure monotonicity on level k-2
!call coeffcr1d_mbez(opacc_km2, opac_jm1km2, opac_jkm2, ayt, byt, cyt, aykm2_opac, bykm2_opac, cykm2_opac)
!aykm2_opac = aykm2_opac*by
!bykm2_opac = bykm2_opac*by + ay
!cykm2_opac = cykm2_opac*by + cy
!opac_km2 = aykm2_opac*opac_jm2km2 + bykm2_opac*opac_jm1km2 + cykm2_opac*opac_jkm2
   call coeffcr1d_mbez(scontc_km2, scont_jm1km2, scont_jkm2, ayt, byt, cyt, aykm2_scont, bykm2_scont, cykm2_scont)
   aykm2_scont = aykm2_scont*by
   bykm2_scont = bykm2_scont*by + ay
   cykm2_scont = cykm2_scont*by + cy
   scont_km2 = aykm2_scont*scont_jm2km2 + bykm2_scont*scont_jm1km2 + cykm2_scont*scont_jkm2
!
!-------------------------interpolation along z-------------------------
!
!calculate control point
!opacc = opac_km2*azt + opac_km1*bzt + opac_k*czt
   scontc = scont_km2*azt + scont_km1*bzt + scont_k*czt
!
!ensure monotonicity
!call coeffcr1d_mbez(opacc, opac_km1, opac_k, azt, bzt, czt, az_opac, bz_opac, cz_opac)
!az_opac = az_opac*bz
!bz_opac = bz_opac*bz + az
!cz_opac = cz_opac*bz + cz
!opac_p = az_opac*opac_km2 + bz_opac*opac_km1 + cz_opac*opac_k
   call coeffcr1d_mbez(scontc, scont_km1, scont_k, azt, bzt, czt, az_scont, bz_scont, cz_scont)
   az_scont = az_scont*bz
   bz_scont = bz_scont*bz + az
   cz_scont = cz_scont*bz + cz
   scont_p = az_scont*scont_km2 + bz_scont*scont_km1 + cz_scont*scont_k
!
   c01 = az_scont*aykm2_scont*axjm2km2_scont
   c02 = az_scont*aykm2_scont*bxjm2km2_scont
   c03 = az_scont*aykm2_scont*cxjm2km2_scont
   c04 = az_scont*bykm2_scont*axjm1km2_scont
   c05 = az_scont*bykm2_scont*bxjm1km2_scont
   c06 = az_scont*bykm2_scont*cxjm1km2_scont
   c07 = az_scont*cykm2_scont*axjkm2_scont
   c08 = az_scont*cykm2_scont*bxjkm2_scont
   c09 = az_scont*cykm2_scont*cxjkm2_scont
   c10 = bz_scont*aykm1_scont*axjm2km1_scont
   c11 = bz_scont*aykm1_scont*bxjm2km1_scont
   c12 = bz_scont*aykm1_scont*cxjm2km1_scont
   c13 = bz_scont*bykm1_scont*axjm1km1_scont
   c14 = bz_scont*bykm1_scont*bxjm1km1_scont
   c15 = bz_scont*bykm1_scont*cxjm1km1_scont
   c16 = bz_scont*cykm1_scont*axjkm1_scont
   c17 = bz_scont*cykm1_scont*bxjkm1_scont
   c18 = bz_scont*cykm1_scont*cxjkm1_scont
   c19 = cz_scont*ayk_scont*axjm2k_scont
   c20 = cz_scont*ayk_scont*bxjm2k_scont
   c21 = cz_scont*ayk_scont*cxjm2k_scont
   c22 = cz_scont*byk_scont*axjm1k_scont
   c23 = cz_scont*byk_scont*bxjm1k_scont
   c24 = cz_scont*byk_scont*cxjm1k_scont
   c25 = cz_scont*cyk_scont*axjk_scont
   c26 = cz_scont*cyk_scont*bxjk_scont
   c27 = cz_scont*cyk_scont*cxjk_scont
!
!
end subroutine coeff3d_contu
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine coeff3d_contu_lin(opac_im1jm1km1, opac_ijm1km1, opac_im1jkm1, opac_ijkm1, &
   opac_im1jm1k,   opac_ijm1k,   opac_im1jk,   opac_ijk, &
   scont_im1jm1km1, scont_ijm1km1, scont_im1jkm1, scont_ijkm1, &
   scont_im1jm1k,   scont_ijm1k,   scont_im1jk,   scont_ijk, &
   x_im1, x_i, y_jm1, y_j, z_km1, z_k, x_p, y_p, z_p, &
   c14, c15, c17, c18, c23, c24, c26, c27, &
   opac_p, scont_p)
!
!         interpolates opacity and continuum source function
!           given on a 3d grid onto point x_p, y_p, z_p
!
!on input:
!
!               f_im1jk--------------f_ijk
!                 /|                  / |
!                / |                 /  |
!               /  |                /   |
!              /   |               /    |
!             /    |              /     |
!            / f_im1jkm1---------/---f_ijkm1
!           /      /            /      /
!       f_im1jm1k------------f_ijm1k  /
!           |    /              |    /
!           |   /               |   /
!           |  /                |  /
!           | /                 | /
!           |/                  |/
!       f_im1jm1km1---------f_ijm1km1
!
!        x_p, y_p, z_p: coordinates of point onto which shall be interpolated
!
!on output:
!   1. interpolation coefficients for source function and intensity
!      (required for ALO calculations):
!         c14, c15, c17, c18,
!         c23, c24, c25, c26, c27
!
!   such that
!      f_p = c14*f_im1jm1km1 + c15*f_ijm1km1 +
!            c17*f_im1jkm1   + c18*f_ijkm1 +
!            c23*f_im1jm1k + c24*f_ijm1k +
!            c26*f_im1jk   + c27*f_ijk
!
!   2. interpolated values at point p: opac_p, scont_p
!
   use prog_type
!
   implicit none
!
! ... argments
   real(dp), intent(in) :: opac_im1jm1km1, opac_ijm1km1, opac_im1jkm1, opac_ijkm1, &
      opac_im1jm1k,   opac_ijm1k,   opac_im1jk,   opac_ijk, &
      scont_im1jm1km1, scont_ijm1km1, scont_im1jkm1, scont_ijkm1, &
      scont_im1jm1k,   scont_ijm1k,   scont_im1jk,   scont_ijk, &
      x_im1, x_i, y_jm1, y_j, z_km1, z_k, x_p, y_p, z_p
   real(dp), intent(out) :: c14, c15, c17, c18, c23, c24, c26, c27, &
      opac_p, scont_p
!
! ... local scalars
   real(dp) :: tx, ty, tz, ax, ay, az
!
!
!define deltax, deltay, deltaz
   tx = (x_p-x_im1)/(x_i-x_im1)
   ty = (y_p-y_jm1)/(y_j-y_jm1)
   tz = (z_p-z_km1)/(z_k-z_km1)
!
!--------------------------trilinear interpolation----------------------
!
   ax=1.d0-tx
   ay=1.d0-ty
   az=1.d0-tz
!
   c14 = ax*ay*az
   c15 = tx*ay*az
   c17 = ax*ty*az
   c18 = tx*ty*az
   c23 = ay*ax*tz
   c24 = tx*ay*tz
   c26 = ax*ty*tz
   c27 = tx*ty*tz

   opac_p = c14*opac_im1jm1km1 + c15*opac_ijm1km1 + c17*opac_im1jkm1   + c18*opac_ijkm1 + &
      c23*opac_im1jm1k + c24*opac_ijm1k + c26*opac_im1jk   + c27*opac_ijk

   scont_p = c14*scont_im1jm1km1 + c15*scont_ijm1km1 + c17*scont_im1jkm1   + c18*scont_ijkm1 + &
      c23*scont_im1jm1k + c24*scont_ijm1k + c26*scont_im1jk   + c27*scont_ijk
!
!return
!
   return
!
!
end subroutine coeff3d_contu_lin
!
!***********************************************************************
!***********************************************************************
!
!                    LINE TRANSFER ROUTINES
!
!***********************************************************************
!***********************************************************************
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine fsc_line3d(oindx,xobsindx)
!
!-----------------------------------------------------------------------
!--------short characteristics for line radiative transfer in 3d--------
!-----------calculating intensties for given mu,phi specified-----------
!-----------by input oindx, and for xobs specified by xobsindx----------
!---------------------using bezier interpolations-----------------------
!-----------------------------------------------------------------------
!
   use prog_type
   use fund_const
   use dime3d, only: ndxmax, ndzmax, ndymax, x, y, z
   use dime3d, only: int3d, opalbar3d, sline3d, velx3d, vely3d, velz3d, vth3d, imask3d, imask_totreg3d, &
      aloline_on_nn3d, aloline_nn3d_tmp, mintbar3d_tmp, normalization3d_tmp
   use angles, only: n_x, n_y, n_z, weight_omega, q_alo
   use freq, only: nodes_xobs, weight_xobs, nxobs, deltax, xcmf_min, xcmf_max
   use params_input, only: vth_fiducial
   use bcondition, only: xic1
   use mod_debug, only: indxx, indxy, indxz
   use params_stellar, only: smajorax_a, smajorax_b, smajorax_c
!
   implicit none
!
! ... arguments
   integer(i4b), intent(in) :: oindx, xobsindx
!
! ... local scalars
   integer(i4b) :: i, j, k
   integer(i4b) :: alpha, beta, gamma
   integer(i4b) :: startx, starty, startz, endx, endy, endz
   integer(i4b) :: iim2, iim1, ii, iip1, jjm2, jjm1, jj, jjp1, kkm2, kkm1, kk, kkp1
   real(dp) :: nn_x, nn_y, nn_z, mueff, xobs, wall, wall_global, phinorm, fdum
   real(dp) :: dels_xyu, dels_xzu, dels_yzu, dels_u, dels_r, &
      dels_xyd, dels_xzd, dels_yzd, dels_d
   real(dp) :: opalbar_p, velx_p, vely_p, velz_p, vel_p, vth_p, sline_p, rad_p
   real(dp) :: x_u, y_u, z_u, opalbar_u, velx_u, vely_u, velz_u, vel_u, vth_u, sline_u, int_u, rad_u
   real(dp) :: x_d, y_d, z_d, opalbar_d, velx_d, vely_d, velz_d, vel_d, vth_d, sline_d, rad_d
   real(dp) :: abs_sc, int_sc, contr_sc
   real(dp) :: alo_u, alo_p, alo_d
   integer :: q1,  q2,  q3,  q4,  q5,  q6,  q7,  q8,  q9, &
      q10, q11, q12, q13, q14, q15, q16, q17, q18, &
      q19, q20, q21, q22, q23, q24, q25, q26, q27
   real(dp) :: c02_slineu, c04_slineu, c05_slineu, c06_slineu, c08_slineu, &
      c10_slineu, c11_slineu, c12_slineu, c13_slineu, c14_slineu, c15_slineu, c16_slineu, c17_slineu, c18_slineu, &
      c20_slineu, c22_slineu, c23_slineu, c24_slineu, c26_slineu, c27_slineu, &
      c02_intu, c04_intu, c05_intu, c06_intu, c08_intu, &
      c10_intu, c11_intu, c12_intu, c13_intu, c14_intu, c15_intu, c16_intu, c17_intu, c18_intu, &
      c20_intu, c22_intu, c23_intu, c24_intu, c26_intu, &
      c03_slined, c06_slined, c07_slined, c08_slined, c09_slined, c12_slined, &
      c15_slined, c16_slined, c17_slined, c18_slined, c19_slined, c20_slined, &
      c21_slined, c22_slined, c23_slined, c24_slined, c25_slined, c26_slined, c27_slined
!
! ... local arrays
!
! ... local functions
   real(dp) :: dist_sphere, dist_ellipsoid, calc_icore_gdark
!
! ... local characters
!character(len=50) :: enter
!
! ... local logicals
!
!directions
   nn_x=n_x(oindx)
   nn_y=n_y(oindx)
   nn_z=n_z(oindx)
!
!frequency
   xobs=nodes_xobs(xobsindx)
!
!angular and frequency  integration weights
   wall_global=weight_omega(oindx)*weight_xobs(xobsindx)
!
!indices for nearest neighbour alo
   q1=q_alo(oindx,1)
   q2=q_alo(oindx,2)
   q3=q_alo(oindx,3)
   q4=q_alo(oindx,4)
   q5=q_alo(oindx,5)
   q6=q_alo(oindx,6)
   q7=q_alo(oindx,7)
   q8=q_alo(oindx,8)
   q9=q_alo(oindx,9)
   q10=q_alo(oindx,10)
   q11=q_alo(oindx,11)
   q12=q_alo(oindx,12)
   q13=q_alo(oindx,13)
   q14=q_alo(oindx,14)
   q15=q_alo(oindx,15)
   q16=q_alo(oindx,16)
   q17=q_alo(oindx,17)
   q18=q_alo(oindx,18)
   q19=q_alo(oindx,19)
   q20=q_alo(oindx,20)
   q21=q_alo(oindx,21)
   q22=q_alo(oindx,22)
   q23=q_alo(oindx,23)
   q24=q_alo(oindx,24)
   q25=q_alo(oindx,25)
   q26=q_alo(oindx,26)
   q27=q_alo(oindx,27)
!
!-----------------------------------------------------------------------
!
!set directional index-parameter (phantom points are excluded from calculation)
!
!index parameter:
!         if n_x >= 0                 if n_x < 0
!                startx = 2                  startx = ndxmax-1
!                endx = ndxmax-1             endx = 2
!                alpha=  1                   alpha=-1
!
!         if n_y >= 0                 if n_y < 0
!                starty = 2                  starty = ndymax-1
!                endy = ndymax-1             endy = 2
!                beta =  1                   beta =-1
!
!         if n_z >= 0                 if n_z < 0
!                startz = 2                  startz = ndzmax-1
!                endz = ndzmax-1             endz = 2
!                gamma = 1                   gamma = -1
!
   if(nn_x.gt.0.d0) then
      startx = 3
      endx = ndxmax-2
      alpha=  1
   elseif(nn_x.lt.0.d0) then
      startx = ndxmax-2
      endx = 3
      alpha=-1
   else
      stop 'error in fsc_line3d: n_x = 0 not allowed'
   endif
!
   if(nn_y.gt.0.d0) then
      starty = 3
      endy = ndymax-2
      beta =  1
   elseif(nn_y.lt.0.d0) then
      starty = ndymax-2
      endy = 3
      beta =-1
   else
      stop 'error in fsc_line3d: n_y = 0 not allowed'
   endif
!
   if(nn_z.gt.0.d0) then
      startz = 3
      endz = ndzmax-2
      gamma=  1
   elseif(nn_z.lt.0.d0) then
      startz = ndzmax-2
      endz = 3
      gamma=-1
   else
      stop 'error in fsc_line3d: n_z = 0 not allowed'
   endif
!
!--------------------reset the intensities and alo----------------------
!
   call set_boundary3d(nn_x, nn_y, nn_z)
!
   aloline_on_nn3d=0.d0
!
!-----------------------------------------------------------------------
!
!***debug start
!alo_d=0.d0
!c03_slined=0.d0
!c06_slined=0.d0
!c07_slined=0.d0
!c08_slined=0.d0
!c09_slined=0.d0
!c12_slined=0.d0
!c15_slined=0.d0
!c16_slined=0.d0
!c17_slined=0.d0
!c18_slined=0.d0
!c19_slined=0.d0
!c20_slined=0.d0
!c21_slined=0.d0
!c22_slined=0.d0
!c23_slined=0.d0
!c24_slined=0.d0
!c25_slined=0.d0
!c26_slined=0.d0
!c27_slined=0.d0
!open(1, file='TRASH/test_interp2d.dat', form='formatted', position='append')
!***debug end
!
   do k=startz, endz, gamma
      do j=starty, endy, beta
         do i=startx, endx, alpha
!
            select case(imask3d(i,j,k))
!
!************************standard radiative transfer********************
!
             case(1)
!
               iim1=i-alpha
               jjm1=j-beta
               kkm1=k-gamma
               iip1=i+alpha
               jjp1=j+beta
               kkp1=k+gamma
!
!calculate distance of upwind point to (previous) xy-plane, xz-plane and yz-plane
               dels_xyu=(z(k)-z(kkm1))/nn_z
               dels_xzu=(y(j)-y(jjm1))/nn_y
               dels_yzu=(x(i)-x(iim1))/nn_x
               dels_u=min(dels_xyu,dels_xzu,dels_yzu)
!
!calculate distance of downwind point to (next) xy-plane, xz-plane and yz-plane
               dels_xyd=(z(kkp1)-z(k))/nn_z
               dels_xzd=(y(jjp1)-y(j))/nn_y
               dels_yzd=(x(iip1)-x(i))/nn_x
               dels_d=min(dels_xyd,dels_xzd,dels_yzd)
!
!----------------------------local point--------------------------------
!
               sline_p=sline3d(i,j,k)
               opalbar_p=opalbar3d(i,j,k)
               velx_p=velx3d(i,j,k)
               vely_p=vely3d(i,j,k)
               velz_p=velz3d(i,j,k)
               vel_p = nn_x*velx_p + nn_y*vely_p + nn_z*velz_p
               vth_p = vth3d(i,j,k)
!
!----------------------------upwind point-------------------------------
!
               if(dels_xyu.eq.dels_u) then
!intersection with x-y plane on level k-gamma
                  x_u = x(i) - dels_u*nn_x
                  y_u = y(j) - dels_u*nn_y
                  z_u = z(kkm1)
!
                  iim2=i-2*alpha
                  jjm2=j-2*beta
!
                  call coeff2d_lineu(opalbar3d(iim2,jjm2,kkm1), opalbar3d(iim1,jjm2,kkm1), opalbar3d(i,jjm2,kkm1), &
                     opalbar3d(iim2,jjm1,kkm1), opalbar3d(iim1,jjm1,kkm1), opalbar3d(i,jjm1,kkm1), &
                     opalbar3d(iim2,j,kkm1),    opalbar3d(iim1,j,kkm1),    opalbar3d(i,j,kkm1), &
                     velx3d(iim2,jjm2,kkm1), velx3d(iim1,jjm2,kkm1), velx3d(i,jjm2,kkm1), &
                     velx3d(iim2,jjm1,kkm1), velx3d(iim1,jjm1,kkm1), velx3d(i,jjm1,kkm1), &
                     velx3d(iim2,j,kkm1),    velx3d(iim1,j,kkm1),    velx3d(i,j,kkm1), &
                     vely3d(iim2,jjm2,kkm1), vely3d(iim1,jjm2,kkm1), vely3d(i,jjm2,kkm1), &
                     vely3d(iim2,jjm1,kkm1), vely3d(iim1,jjm1,kkm1), vely3d(i,jjm1,kkm1), &
                     vely3d(iim2,j,kkm1),    vely3d(iim1,j,kkm1),    vely3d(i,j,kkm1), &
                     velz3d(iim2,jjm2,kkm1), velz3d(iim1,jjm2,kkm1), velz3d(i,jjm2,kkm1), &
                     velz3d(iim2,jjm1,kkm1), velz3d(iim1,jjm1,kkm1), velz3d(i,jjm1,kkm1), &
                     velz3d(iim2,j,kkm1),    velz3d(iim1,j,kkm1),    velz3d(i,j,kkm1), &
                     vth3d(iim2,jjm2,kkm1), vth3d(iim1,jjm2,kkm1), vth3d(i,jjm2,kkm1), &
                     vth3d(iim2,jjm1,kkm1), vth3d(iim1,jjm1,kkm1), vth3d(i,jjm1,kkm1), &
                     vth3d(iim2,j,kkm1),    vth3d(iim1,j,kkm1),    vth3d(i,j,kkm1), &
                     sline3d(iim2,jjm2,kkm1), sline3d(iim1,jjm2,kkm1), sline3d(i,jjm2,kkm1), &
                     sline3d(iim2,jjm1,kkm1), sline3d(iim1,jjm1,kkm1), sline3d(i,jjm1,kkm1), &
                     sline3d(iim2,j,kkm1),    sline3d(iim1,j,kkm1),    sline3d(i,j,kkm1), &
                     int3d(iim2,jjm2,kkm1), int3d(iim1,jjm2,kkm1), int3d(i,jjm2,kkm1), &
                     int3d(iim2,jjm1,kkm1), int3d(iim1,jjm1,kkm1), int3d(i,jjm1,kkm1), &
                     int3d(iim2,j,kkm1),    int3d(iim1,j,kkm1),    int3d(i,j,kkm1), &
                     x(iim2), x(iim1), x(i), y(jjm2), y(jjm1), y(j), x_u, y_u, &
                     c10_slineu, c11_slineu, c12_slineu, c13_slineu, c14_slineu, &
                     c15_slineu, c16_slineu, c17_slineu, c18_slineu, &
                     c10_intu, c11_intu, c12_intu, c13_intu, c14_intu, &
                     c15_intu, c16_intu, c17_intu, c18_intu, &
                     opalbar_u, velx_u, vely_u, velz_u, vth_u, sline_u, int_u)

!                  if(abs(c10_slineu*sline3d(iim2,jjm2,kkm1)+c11_slineu*sline3d(iim1,jjm2,kkm1)+c12_slineu*sline3d(i,jjm2,kkm1)+&
!                         c13_slineu*sline3d(iim2,jjm1,kkm1)+c14_slineu*sline3d(iim1,jjm1,kkm1)+c15_slineu*sline3d(i,jjm1,kkm1)+&
!                         c16_slineu*sline3d(iim2,j   ,kkm1)+c17_slineu*sline3d(iim1,j   ,kkm1)+c18_slineu*sline3d(i,j   ,kkm1)-sline_u).gt.1.d-14) then
!                     write(*,*) sline_u, (c10_slineu*sline3d(iim2,jjm2,kkm1)+c11_slineu*sline3d(iim1,jjm2,kkm1)+c12_slineu*sline3d(i,jjm2,kkm1)+&
!                                          c13_slineu*sline3d(iim2,jjm1,kkm1)+c14_slineu*sline3d(iim1,jjm1,kkm1)+c15_slineu*sline3d(i,jjm1,kkm1)+&
!                                          c16_slineu*sline3d(iim2,j   ,kkm1)+c17_slineu*sline3d(iim1,j   ,kkm1)+c18_slineu*sline3d(i,j   ,kkm1))
!                     stop
!                  endif
!                  if(abs(c10_intu*int3d(iim2,jjm2,kkm1)+c11_intu*int3d(iim1,jjm2,kkm1)+c12_intu*int3d(i,jjm2,kkm1)+&
!                         c13_intu*int3d(iim2,jjm1,kkm1)+c14_intu*int3d(iim1,jjm1,kkm1)+c15_intu*int3d(i,jjm1,kkm1)+&
!                         c16_intu*int3d(iim2,j   ,kkm1)+c17_intu*int3d(iim1,j   ,kkm1)+c18_intu*int3d(i,j   ,kkm1)-int_u).gt.1.d-14) then
!                     write(*,*) int_u, (c10_intu*int3d(iim2,jjm2,kkm1)+c11_intu*int3d(iim1,jjm2,kkm1)+c12_intu*int3d(i,jjm2,kkm1)+&
!                                        c13_intu*int3d(iim2,jjm1,kkm1)+c14_intu*int3d(iim1,jjm1,kkm1)+c15_intu*int3d(i,jjm1,kkm1)+&
!                                        c16_intu*int3d(iim2,j   ,kkm1)+c17_intu*int3d(iim1,j   ,kkm1)+c18_intu*int3d(i,j   ,kkm1))
!                     stop
!                  endif


!set interpolation coefficients that are not used to zero
                  c23_slineu=0.d0
                  c24_slineu=0.d0
                  c26_slineu=0.d0
                  c27_slineu=0.d0
!
                  c02_intu = 0.d0
                  c04_intu = 0.d0
                  c05_intu = 0.d0
                  c06_intu = 0.d0
                  c08_intu = 0.d0
                  c20_intu = 0.d0
                  c22_intu = 0.d0
                  c23_intu = 0.d0
                  c24_intu = 0.d0
                  c26_intu = 0.d0
!
               elseif(dels_xzu.eq.dels_u) then
!intersection with x-z plane on level j-beta
                  x_u = x(i) - dels_u*nn_x
                  y_u = y(jjm1)
                  z_u = z(k) - dels_u*nn_z
!
                  iim2=i-2*alpha
                  kkm2=k-2*gamma
!
                  call coeff2d_lineu(opalbar3d(iim2,jjm1,kkm2), opalbar3d(iim1,jjm1,kkm2), opalbar3d(i,jjm1,kkm2), &
                     opalbar3d(iim2,jjm1,kkm1), opalbar3d(iim1,jjm1,kkm1), opalbar3d(i,jjm1,kkm1), &
                     opalbar3d(iim2,jjm1,k),    opalbar3d(iim1,jjm1,k),    opalbar3d(i,jjm1,k), &
                     velx3d(iim2,jjm1,kkm2),    velx3d(iim1,jjm1,kkm2),    velx3d(i,jjm1,kkm2), &
                     velx3d(iim2,jjm1,kkm1),    velx3d(iim1,jjm1,kkm1),    velx3d(i,jjm1,kkm1), &
                     velx3d(iim2,jjm1,k),       velx3d(iim1,jjm1,k),       velx3d(i,jjm1,k), &
                     vely3d(iim2,jjm1,kkm2),    vely3d(iim1,jjm1,kkm2),    vely3d(i,jjm1,kkm2), &
                     vely3d(iim2,jjm1,kkm1),    vely3d(iim1,jjm1,kkm1),    vely3d(i,jjm1,kkm1), &
                     vely3d(iim2,jjm1,k),       vely3d(iim1,jjm1,k),       vely3d(i,jjm1,k), &
                     velz3d(iim2,jjm1,kkm2),    velz3d(iim1,jjm1,kkm2),    velz3d(i,jjm1,kkm2), &
                     velz3d(iim2,jjm1,kkm1),    velz3d(iim1,jjm1,kkm1),    velz3d(i,jjm1,kkm1), &
                     velz3d(iim2,jjm1,k),       velz3d(iim1,jjm1,k),       velz3d(i,jjm1,k), &
                     vth3d(iim2,jjm1,kkm2),     vth3d(iim1,jjm1,kkm2),     vth3d(i,jjm1,kkm2), &
                     vth3d(iim2,jjm1,kkm1),     vth3d(iim1,jjm1,kkm1),     vth3d(i,jjm1,kkm1), &
                     vth3d(iim2,jjm1,k),        vth3d(iim1,jjm1,k),        vth3d(i,jjm1,k), &
                     sline3d(iim2,jjm1,kkm2),   sline3d(iim1,jjm1,kkm2),   sline3d(i,jjm1,kkm2), &
                     sline3d(iim2,jjm1,kkm1),   sline3d(iim1,jjm1,kkm1),   sline3d(i,jjm1,kkm1), &
                     sline3d(iim2,jjm1,k),      sline3d(iim1,jjm1,k),      sline3d(i,jjm1,k), &
                     int3d(iim2,jjm1,kkm2),     int3d(iim1,jjm1,kkm2),     int3d(i,jjm1,kkm2), &
                     int3d(iim2,jjm1,kkm1),     int3d(iim1,jjm1,kkm1),     int3d(i,jjm1,kkm1), &
                     int3d(iim2,jjm1,k),        int3d(iim1,jjm1,k),        int3d(i,jjm1,k), &
                     x(iim2), x(iim1), x(i), z(kkm2), z(kkm1), z(k), x_u, z_u, &
                     c04_slineu, c05_slineu, c06_slineu, c13_slineu, c14_slineu, &
                     c15_slineu, c22_slineu, c23_slineu, c24_slineu, &
                     c04_intu, c05_intu, c06_intu, c13_intu, c14_intu, &
                     c15_intu, c22_intu, c23_intu, c24_intu, &
                     opalbar_u, velx_u, vely_u, velz_u, vth_u, sline_u, int_u)

!                  if(abs(c04_slineu*sline3d(iim2,jjm1,kkm2)+c05_slineu*sline3d(iim1,jjm1,kkm2)+c06_slineu*sline3d(i,jjm1,kkm2)+&
!                         c13_slineu*sline3d(iim2,jjm1,kkm1)+c14_slineu*sline3d(iim1,jjm1,kkm1)+c15_slineu*sline3d(i,jjm1,kkm1)+&
!                         c22_slineu*sline3d(iim2,jjm1,k   )+c23_slineu*sline3d(iim1,jjm1,k   )+c24_slineu*sline3d(i,jjm1,k   )-sline_u).gt.1.d-14) stop 'error coeff'
!                  if(abs(c04_intu*int3d(iim2,jjm1,kkm2)+c05_intu*int3d(iim1,jjm1,kkm2)+c06_intu*int3d(i,jjm1,kkm2)+&
!                         c13_intu*int3d(iim2,jjm1,kkm1)+c14_intu*int3d(iim1,jjm1,kkm1)+c15_intu*int3d(i,jjm1,kkm1)+&
!                         c22_intu*int3d(iim2,jjm1,k   )+c23_intu*int3d(iim1,jjm1,k   )+c24_intu*int3d(i,jjm1,k   )-int_u).gt.1.d-14) stop 'error coeff'
!
                  c17_slineu=0.d0
                  c18_slineu=0.d0
                  c26_slineu=0.d0
                  c27_slineu=0.d0
                  c02_intu = 0.d0
                  c08_intu = 0.d0
                  c10_intu = 0.d0
                  c11_intu = 0.d0
                  c12_intu = 0.d0
                  c16_intu = 0.d0
                  c17_intu = 0.d0
                  c18_intu = 0.d0
                  c20_intu = 0.d0
                  c26_intu = 0.d0
!
               elseif(dels_yzu.eq.dels_u) then
!intersection with y-z plane on level i-alpha
                  x_u = x(iim1)
                  y_u = y(j) - dels_u*nn_y
                  z_u = z(k) - dels_u*nn_z
!
                  jjm2=j-2*beta
                  kkm2=k-2*gamma
!
                  call coeff2d_lineu(opalbar3d(iim1,jjm2,kkm2), opalbar3d(iim1,jjm1,kkm2), opalbar3d(iim1,j,kkm2), &
                     opalbar3d(iim1,jjm2,kkm1), opalbar3d(iim1,jjm1,kkm1), opalbar3d(iim1,j,kkm1), &
                     opalbar3d(iim1,jjm2,k),    opalbar3d(iim1,jjm1,k),    opalbar3d(iim1,j,k), &
                     velx3d(iim1,jjm2,kkm2),    velx3d(iim1,jjm1,kkm2),    velx3d(iim1,j,kkm2), &
                     velx3d(iim1,jjm2,kkm1),    velx3d(iim1,jjm1,kkm1),    velx3d(iim1,j,kkm1), &
                     velx3d(iim1,jjm2,k),       velx3d(iim1,jjm1,k),       velx3d(iim1,j,k), &
                     vely3d(iim1,jjm2,kkm2),    vely3d(iim1,jjm1,kkm2),    vely3d(iim1,j,kkm2), &
                     vely3d(iim1,jjm2,kkm1),    vely3d(iim1,jjm1,kkm1),    vely3d(iim1,j,kkm1), &
                     vely3d(iim1,jjm2,k),       vely3d(iim1,jjm1,k),       vely3d(iim1,j,k), &
                     velz3d(iim1,jjm2,kkm2),    velz3d(iim1,jjm1,kkm2),    velz3d(iim1,j,kkm2), &
                     velz3d(iim1,jjm2,kkm1),    velz3d(iim1,jjm1,kkm1),    velz3d(iim1,j,kkm1), &
                     velz3d(iim1,jjm2,k),       velz3d(iim1,jjm1,k),       velz3d(iim1,j,k), &
                     vth3d(iim1,jjm2,kkm2),     vth3d(iim1,jjm1,kkm2),     vth3d(iim1,j,kkm2), &
                     vth3d(iim1,jjm2,kkm1),     vth3d(iim1,jjm1,kkm1),     vth3d(iim1,j,kkm1), &
                     vth3d(iim1,jjm2,k),        vth3d(iim1,jjm1,k),        vth3d(iim1,j,k), &
                     sline3d(iim1,jjm2,kkm2),   sline3d(iim1,jjm1,kkm2),   sline3d(iim1,j,kkm2), &
                     sline3d(iim1,jjm2,kkm1),   sline3d(iim1,jjm1,kkm1),   sline3d(iim1,j,kkm1), &
                     sline3d(iim1,jjm2,k),      sline3d(iim1,jjm1,k),      sline3d(iim1,j,k), &
                     int3d(iim1,jjm2,kkm2),     int3d(iim1,jjm1,kkm2),     int3d(iim1,j,kkm2), &
                     int3d(iim1,jjm2,kkm1),     int3d(iim1,jjm1,kkm1),     int3d(iim1,j,kkm1), &
                     int3d(iim1,jjm2,k),        int3d(iim1,jjm1,k),        int3d(iim1,j,k), &
                     y(jjm2), y(jjm1), y(j), z(kkm2), z(kkm1), z(k), y_u, z_u, &
                     c02_slineu, c05_slineu, c08_slineu, c11_slineu, c14_slineu, &
                     c17_slineu, c20_slineu, c23_slineu, c26_slineu, &
                     c02_intu, c05_intu, c08_intu, c11_intu, c14_intu, &
                     c17_intu, c20_intu, c23_intu, c26_intu, &
                     opalbar_u, velx_u, vely_u, velz_u, vth_u, sline_u, int_u)

!                  if(abs(c02_slineu*sline3d(iim1,jjm2,kkm2)+c05_slineu*sline3d(iim1,jjm1,kkm2)+c08_slineu*sline3d(iim1,j,kkm2)+&
!                         c11_slineu*sline3d(iim1,jjm2,kkm1)+c14_slineu*sline3d(iim1,jjm1,kkm1)+c17_slineu*sline3d(iim1,j,kkm1)+&
!                         c20_slineu*sline3d(iim1,jjm2,k   )+c23_slineu*sline3d(iim1,jjm1,k   )+c26_slineu*sline3d(iim1,j,k   )-sline_u).gt.1.d-14) stop 'error coeff'
!                  if(abs(c02_intu*int3d(iim1,jjm2,kkm2)+c05_intu*int3d(iim1,jjm1,kkm2)+c08_intu*int3d(iim1,j,kkm2)+&
!                         c11_intu*int3d(iim1,jjm2,kkm1)+c14_intu*int3d(iim1,jjm1,kkm1)+c17_intu*int3d(iim1,j,kkm1)+&
!                         c20_intu*int3d(iim1,jjm2,k   )+c23_intu*int3d(iim1,jjm1,k   )+c26_intu*int3d(iim1,j,k   )-int_u).gt.1.d-14) stop 'error coeff'

                  c15_slineu=0.d0
                  c18_slineu=0.d0
                  c24_slineu=0.d0
                  c27_slineu=0.d0
!
                  c04_intu = 0.d0
                  c06_intu = 0.d0
                  c10_intu = 0.d0
                  c12_intu = 0.d0
                  c13_intu = 0.d0
                  c15_intu = 0.d0
                  c16_intu = 0.d0
                  c18_intu = 0.d0
                  c22_intu = 0.d0
                  c24_intu = 0.d0
!
               else
                  write(*,'(4es20.8,2l4)') dels_u, dels_xzu, dels_xyu, dels_yzu, dels_u.eq.dels_xzu, dels_xyu.eq.dels_u
                  stop 'error in fsc_line3d: invalid dels_u'
               endif
!
!---------------------------downwind point------------------------------
!
               if(dels_xyd.eq.dels_d) then
!intersection with x-y plane on level k+gamma
                  x_d = x(i) + dels_d*nn_x
                  y_d = y(j) + dels_d*nn_y
                  z_d = z(kkp1)
!
                  call coeff2d_lined(opalbar3d(iim1,jjm1,kkp1), opalbar3d(i,jjm1,kkp1), opalbar3d(iip1,jjm1,kkp1), &
                     opalbar3d(iim1,j,kkp1),    opalbar3d(i,j,kkp1),    opalbar3d(iip1,j,kkp1), &
                     opalbar3d(iim1,jjp1,kkp1), opalbar3d(i,jjp1,kkp1), opalbar3d(iip1,jjp1,kkp1), &
                     velx3d(iim1,jjm1,kkp1),    velx3d(i,jjm1,kkp1),    velx3d(iip1,jjm1,kkp1), &
                     velx3d(iim1,j,kkp1),       velx3d(i,j,kkp1),       velx3d(iip1,j,kkp1), &
                     velx3d(iim1,jjp1,kkp1),    velx3d(i,jjp1,kkp1),    velx3d(iip1,jjp1,kkp1), &
                     vely3d(iim1,jjm1,kkp1),    vely3d(i,jjm1,kkp1),    vely3d(iip1,jjm1,kkp1), &
                     vely3d(iim1,j,kkp1),       vely3d(i,j,kkp1),       vely3d(iip1,j,kkp1), &
                     vely3d(iim1,jjp1,kkp1),    vely3d(i,jjp1,kkp1),    vely3d(iip1,jjp1,kkp1), &
                     velz3d(iim1,jjm1,kkp1),    velz3d(i,jjm1,kkp1),    velz3d(iip1,jjm1,kkp1), &
                     velz3d(iim1,j,kkp1),       velz3d(i,j,kkp1),       velz3d(iip1,j,kkp1), &
                     velz3d(iim1,jjp1,kkp1),    velz3d(i,jjp1,kkp1),    velz3d(iip1,jjp1,kkp1), &
                     vth3d(iim1,jjm1,kkp1),      vth3d(i,jjm1,kkp1),    vth3d(iip1,jjm1,kkp1), &
                     vth3d(iim1,j,kkp1),         vth3d(i,j,kkp1),       vth3d(iip1,j,kkp1), &
                     vth3d(iim1,jjp1,kkp1),      vth3d(i,jjp1,kkp1),    vth3d(iip1,jjp1,kkp1), &
                     sline3d(iim1,jjm1,kkp1),    sline3d(i,jjm1,kkp1),  sline3d(iip1,jjm1,kkp1), &
                     sline3d(iim1,j,kkp1),       sline3d(i,j,kkp1),     sline3d(iip1,j,kkp1), &
                     sline3d(iim1,jjp1,kkp1),    sline3d(i,jjp1,kkp1),  sline3d(iip1,jjp1,kkp1), &
                     x(iim1), x(i), x(iip1), y(jjm1), y(j), y(jjp1), x_d, y_d, &
                     c19_slined, c20_slined, c21_slined, c22_slined, c23_slined, c24_slined, &
                     c25_slined, c26_slined, c27_slined, &
                     opalbar_d, velx_d, vely_d, velz_d, vth_d, sline_d)

!                 if(abs(c19_slined*sline3d(iim1,jjm1,kkp1)+c20_slined*sline3d(i,jjm1,kkp1)+c21_slined*sline3d(iip1,jjm1,kkp1)+&
!                        c22_slined*sline3d(iim1,j   ,kkp1)+c23_slined*sline3d(i,j   ,kkp1)+c24_slined*sline3d(iip1,j   ,kkp1)+&
!                        c25_slined*sline3d(iim1,jjp1,kkp1)+c26_slined*sline3d(i,jjp1,kkp1)+c27_slined*sline3d(iip1,jjp1,kkp1)-sline_d).gt.1.d-14) stop 'error coeff'
!
!set interpolation coefficients that are not used to zero
                  c03_slined = 0.d0
                  c06_slined = 0.d0
                  c07_slined = 0.d0
                  c08_slined = 0.d0
                  c09_slined = 0.d0
                  c12_slined = 0.d0
                  c15_slined = 0.d0
                  c16_slined = 0.d0
                  c17_slined = 0.d0
                  c18_slined = 0.d0
!
               elseif(dels_xzd.eq.dels_d) then
!intersection with x-z plane on level j+beta
                  x_d = x(i) + dels_d*nn_x
                  y_d = y(jjp1)
                  z_d = z(k) + dels_d*nn_z
!
                  call coeff2d_lined(opalbar3d(iim1,jjp1,kkm1), opalbar3d(i,jjp1,kkm1), opalbar3d(iip1,jjp1,kkm1), &
                     opalbar3d(iim1,jjp1,k),    opalbar3d(i,jjp1,k),    opalbar3d(iip1,jjp1,k), &
                     opalbar3d(iim1,jjp1,kkp1), opalbar3d(i,jjp1,kkp1), opalbar3d(iip1,jjp1,kkp1), &
                     velx3d(iim1,jjp1,kkm1),    velx3d(i,jjp1,kkm1),    velx3d(iip1,jjp1,kkm1), &
                     velx3d(iim1,jjp1,k),       velx3d(i,jjp1,k),       velx3d(iip1,jjp1,k), &
                     velx3d(iim1,jjp1,kkp1),    velx3d(i,jjp1,kkp1),    velx3d(iip1,jjp1,kkp1), &
                     vely3d(iim1,jjp1,kkm1),    vely3d(i,jjp1,kkm1),    vely3d(iip1,jjp1,kkm1), &
                     vely3d(iim1,jjp1,k),       vely3d(i,jjp1,k),       vely3d(iip1,jjp1,k), &
                     vely3d(iim1,jjp1,kkp1),    vely3d(i,jjp1,kkp1),    vely3d(iip1,jjp1,kkp1), &
                     velz3d(iim1,jjp1,kkm1),    velz3d(i,jjp1,kkm1),    velz3d(iip1,jjp1,kkm1), &
                     velz3d(iim1,jjp1,k),       velz3d(i,jjp1,k),       velz3d(iip1,jjp1,k), &
                     velz3d(iim1,jjp1,kkp1),    velz3d(i,jjp1,kkp1),    velz3d(iip1,jjp1,kkp1), &
                     vth3d(iim1,jjp1,kkm1),     vth3d(i,jjp1,kkm1),     vth3d(iip1,jjp1,kkm1), &
                     vth3d(iim1,jjp1,k),        vth3d(i,jjp1,k),        vth3d(iip1,jjp1,k), &
                     vth3d(iim1,jjp1,kkp1),     vth3d(i,jjp1,kkp1),     vth3d(iip1,jjp1,kkp1), &
                     sline3d(iim1,jjp1,kkm1),   sline3d(i,jjp1,kkm1),   sline3d(iip1,jjp1,kkm1), &
                     sline3d(iim1,jjp1,k),      sline3d(i,jjp1,k),      sline3d(iip1,jjp1,k), &
                     sline3d(iim1,jjp1,kkp1),   sline3d(i,jjp1,kkp1),   sline3d(iip1,jjp1,kkp1), &
                     x(iim1), x(i), x(iip1), z(kkm1), z(k), z(kkp1), x_d, z_d, &
                     c07_slined, c08_slined, c09_slined, c16_slined, c17_slined, c18_slined, &
                     c25_slined, c26_slined, c27_slined, &
                     opalbar_d, velx_d, vely_d, velz_d, vth_d, sline_d)

!                 if(abs(c07_slined*sline3d(iim1,jjp1,kkm1)+c08_slined*sline3d(i,jjp1,kkm1)+c09_slined*sline3d(iip1,jjp1,kkm1)+&
!                        c16_slined*sline3d(iim1,jjp1,k   )+c17_slined*sline3d(i,jjp1,k   )+c18_slined*sline3d(iip1,jjp1,k   )+&
!                        c25_slined*sline3d(iim1,jjp1,kkp1)+c26_slined*sline3d(i,jjp1,kkp1)+c27_slined*sline3d(iip1,jjp1,kkp1)-sline_d).gt.1.d-14) stop 'error coeff'
!
!set interpolation coefficients that are not used to zero
                  c03_slined = 0.d0
                  c06_slined = 0.d0
                  c12_slined = 0.d0
                  c15_slined = 0.d0
                  c19_slined = 0.d0
                  c20_slined = 0.d0
                  c21_slined = 0.d0
                  c22_slined = 0.d0
                  c23_slined = 0.d0
                  c24_slined = 0.d0
!
               elseif(dels_yzd.eq.dels_d) then
!intersection with y-z plane on level i+alpha
                  x_d = x(iip1)
                  y_d = y(j) + dels_d*nn_y
                  z_d = z(k) + dels_d*nn_z
!
                  call coeff2d_lined(opalbar3d(iip1,jjm1,kkm1), opalbar3d(iip1,j,kkm1), opalbar3d(iip1,jjp1,kkm1), &
                     opalbar3d(iip1,jjm1,k),    opalbar3d(iip1,j,k),    opalbar3d(iip1,jjp1,k), &
                     opalbar3d(iip1,jjm1,kkp1), opalbar3d(iip1,j,kkp1), opalbar3d(iip1,jjp1,kkp1), &
                     velx3d(iip1,jjm1,kkm1),    velx3d(iip1,j,kkm1),    velx3d(iip1,jjp1,kkm1), &
                     velx3d(iip1,jjm1,k),       velx3d(iip1,j,k),       velx3d(iip1,jjp1,k), &
                     velx3d(iip1,jjm1,kkp1),    velx3d(iip1,j,kkp1),    velx3d(iip1,jjp1,kkp1), &
                     vely3d(iip1,jjm1,kkm1),    vely3d(iip1,j,kkm1),    vely3d(iip1,jjp1,kkm1), &
                     vely3d(iip1,jjm1,k),       vely3d(iip1,j,k),       vely3d(iip1,jjp1,k), &
                     vely3d(iip1,jjm1,kkp1),    vely3d(iip1,j,kkp1),    vely3d(iip1,jjp1,kkp1), &
                     velz3d(iip1,jjm1,kkm1),    velz3d(iip1,j,kkm1),    velz3d(iip1,jjp1,kkm1), &
                     velz3d(iip1,jjm1,k),       velz3d(iip1,j,k),       velz3d(iip1,jjp1,k), &
                     velz3d(iip1,jjm1,kkp1),    velz3d(iip1,j,kkp1),    velz3d(iip1,jjp1,kkp1), &
                     vth3d(iip1,jjm1,kkm1),     vth3d(iip1,j,kkm1),     vth3d(iip1,jjp1,kkm1), &
                     vth3d(iip1,jjm1,k),        vth3d(iip1,j,k),        vth3d(iip1,jjp1,k), &
                     vth3d(iip1,jjm1,kkp1),     vth3d(iip1,j,kkp1),     vth3d(iip1,jjp1,kkp1), &
                     sline3d(iip1,jjm1,kkm1),   sline3d(iip1,j,kkm1),   sline3d(iip1,jjp1,kkm1), &
                     sline3d(iip1,jjm1,k),      sline3d(iip1,j,k),      sline3d(iip1,jjp1,k), &
                     sline3d(iip1,jjm1,kkp1),   sline3d(iip1,j,kkp1),   sline3d(iip1,jjp1,kkp1), &
                     y(jjm1), y(j), y(jjp1), z(kkm1), z(k), z(kkp1), y_d, z_d, &
                     c03_slined, c06_slined, c09_slined, c12_slined, c15_slined, c18_slined, &
                     c21_slined, c24_slined, c27_slined, &
                     opalbar_d, velx_d, vely_d, velz_d, vth_d, sline_d)

!                 if(abs(c03_slined*sline3d(iip1,jjm1,kkm1)+c06_slined*sline3d(iip1,j   ,kkm1)+c09_slined*sline3d(iip1,jjp1,kkm1)+&
!                        c12_slined*sline3d(iip1,jjm1,k   )+c15_slined*sline3d(iip1,j   ,k   )+c18_slined*sline3d(iip1,jjp1,k   )+&
!                        c21_slined*sline3d(iip1,jjm1,kkp1)+c24_slined*sline3d(iip1,j   ,kkp1)+c27_slined*sline3d(iip1,jjp1,kkp1)-sline_d).gt.1.d-14) stop 'error coeff'
!
!set interpolation coefficients that are not used to zero
                  c07_slined = 0.d0
                  c08_slined = 0.d0
                  c16_slined = 0.d0
                  c17_slined = 0.d0
                  c19_slined = 0.d0
                  c20_slined = 0.d0
                  c22_slined = 0.d0
                  c23_slined = 0.d0
                  c25_slined = 0.d0
                  c26_slined = 0.d0
!
               else
                  write(*,'(4es20.8,l4)') dels_d, dels_xzd, dels_xyd, dels_yzd, dels_d.eq.dels_xzd
                  stop 'error in fsc_line3d: invalid dels_d'
               endif
!
!--------------------------------radiative transfer---------------------
!
!***debug start: analytical laws
!               call model_debug(x_u, y_u, z_u, velx_u, vely_u, velz_u, opalbar_u, fdum)
!               call model_debug(x_d, y_d, z_d, velx_d, vely_d, velz_d, opalbar_d, fdum)
!***debug end
               vel_u = velx_u*nn_x + vely_u*nn_y + velz_u*nn_z
               vel_d = velx_d*nn_x + vely_d*nn_y + velz_d*nn_z
               call fsc_line(xobs, vth_fiducial, deltax, xcmf_min, xcmf_max, int_u, opalbar_u, opalbar_p, opalbar_d, &
                  sline_u, sline_p, sline_d, vel_u, vel_p, vel_d, vth_u, vth_p, vth_d, &
                  dels_u, dels_d, abs_sc, contr_sc, int_sc, alo_u, alo_p, alo_d)
!               call fsc_line_lin(xobs, vth_fiducial, deltax, xcmf_min, xcmf_max, int_u, opalbar_u, opalbar_p, sline_u, sline_p, &
!                                   vel_u, vel_p, vth_u, vth_p, dels_u, abs_sc, contr_sc, int_sc, alo_u, alo_p)
!               alo_d=0.d0
               int3d(i,j,k) = abs_sc*int_u + alo_u*sline_u + alo_p*sline_p + alo_d*sline_d
!
               aloline_on_nn3d(iip1,jjp1,kkp1,q1) = imask_totreg3d(iip1,jjp1,kkp1)*&
                  alo_d*c27_slined
               aloline_on_nn3d(i,jjp1,kkp1,q2) = imask_totreg3d(i,jjp1,kkp1) * &
                  (alo_d*c26_slined + abs_sc*c26_intu*aloline_on_nn3d(i,jjp1,kkp1,q1))
               aloline_on_nn3d(iim1,jjp1,kkp1,q3) = imask_totreg3d(iim1,jjp1,kkp1) * &
                  (alo_d*c25_slined + abs_sc*c26_intu*aloline_on_nn3d(iim1,jjp1,kkp1,q2))
               aloline_on_nn3d(iip1,j,kkp1,q4) = imask_totreg3d(iip1,j,kkp1) * &
                  (alo_d*c24_slined + abs_sc*c24_intu*aloline_on_nn3d(iip1,j,kkp1,q1))
               aloline_on_nn3d(i,j,kkp1,q5) = imask_totreg3d(i,j,kkp1) * &
                  (alo_d*c23_slined + abs_sc*(c23_intu*aloline_on_nn3d(i,j,kkp1,q1) + &
                  c24_intu*aloline_on_nn3d(i,j,kkp1,q2) + &
                  c26_intu*aloline_on_nn3d(i,j,kkp1,q4)))
               aloline_on_nn3d(iim1,j,kkp1,q6) = imask_totreg3d(iim1,j,kkp1) * &
                  (alo_d*c22_slined + abs_sc*(c22_intu*aloline_on_nn3d(iim1,j,kkp1,q1) + &
                  c23_intu*aloline_on_nn3d(iim1,j,kkp1,q2) + &
                  c24_intu*aloline_on_nn3d(iim1,j,kkp1,q3) + &
                  c26_intu*aloline_on_nn3d(iim1,j,kkp1,q5)))
               aloline_on_nn3d(iip1,jjm1,kkp1,q7) = imask_totreg3d(iip1,jjm1,kkp1) * &
                  (alo_d*c21_slined + abs_sc*c24_intu*aloline_on_nn3d(iip1,jjm1,kkp1,q4))
               aloline_on_nn3d(i,jjm1,kkp1,q8) = imask_totreg3d(i,jjm1,kkp1) * &
                  (alo_d*c20_slined + abs_sc*(c23_intu*aloline_on_nn3d(i,jjm1,kkp1,q4) + &
                  c24_intu*aloline_on_nn3d(i,jjm1,kkp1,q5) + &
                  c20_intu*aloline_on_nn3d(i,jjm1,kkp1,q1) + &
                  c26_intu*aloline_on_nn3d(i,jjm1,kkp1,q7)))
               aloline_on_nn3d(iim1,jjm1,kkp1,q9) = imask_totreg3d(iim1,jjm1,kkp1) * &
                  (alo_d*c19_slined + abs_sc*(c22_intu*aloline_on_nn3d(iim1,jjm1,kkp1,q4) + &
                  c23_intu*aloline_on_nn3d(iim1,jjm1,kkp1,q5) + &
                  c24_intu*aloline_on_nn3d(iim1,jjm1,kkp1,q6) + &
                  c20_intu*aloline_on_nn3d(iim1,jjm1,kkp1,q2) + &
                  c26_intu*aloline_on_nn3d(iim1,jjm1,kkp1,q8)))
               aloline_on_nn3d(iip1,jjp1,k,q10) = imask_totreg3d(iip1,jjp1,k) * &
                  (alo_d*c18_slined + abs_sc*c18_intu*aloline_on_nn3d(iip1,jjp1,k,q1))
               aloline_on_nn3d(i,jjp1,k,q11) = imask_totreg3d(i,jjp1,k) * &
                  (alo_d*c17_slined + abs_sc*(c17_intu*aloline_on_nn3d(i,jjp1,k,q1) + &
                  c18_intu*aloline_on_nn3d(i,jjp1,k,q2) + &
                  c26_intu*aloline_on_nn3d(i,jjp1,k,q10)))
               aloline_on_nn3d(iim1,jjp1,k,q12) = imask_totreg3d(iim1,jjp1,k) * &
                  (alo_d*c16_slined + abs_sc*(c16_intu*aloline_on_nn3d(iim1,jjp1,k,q1) + &
                  c17_intu*aloline_on_nn3d(iim1,jjp1,k,q2) + &
                  c18_intu*aloline_on_nn3d(iim1,jjp1,k,q3) + &
                  c26_intu*aloline_on_nn3d(iim1,jjp1,k,q11)))
               aloline_on_nn3d(iip1,j,k,q13) = imask_totreg3d(iip1,j,k) * &
                  (alo_d*c15_slined + abs_sc*(c15_intu*aloline_on_nn3d(iip1,j,k,q1) + &
                  c18_intu*aloline_on_nn3d(iip1,j,k,q4) + &
                  c24_intu*aloline_on_nn3d(iip1,j,k,q10)))
               aloline_on_nn3d(i,j,k,q14) = imask_totreg3d(i,j,k) * &
                  (alo_p + abs_sc*(c14_intu*aloline_on_nn3d(i,j,k,q1) + &
                  c15_intu*aloline_on_nn3d(i,j,k,q2) + &
                  c17_intu*aloline_on_nn3d(i,j,k,q4) + &
                  c26_intu*aloline_on_nn3d(i,j,k,q13) + &
                  c18_intu*aloline_on_nn3d(i,j,k,q5) + &
                  c23_intu*aloline_on_nn3d(i,j,k,q10) + &
                  c24_intu*aloline_on_nn3d(i,j,k,q11)))
               aloline_on_nn3d(iim1,j,k,q15) = imask_totreg3d(iim1,j,k) * &
                  (alo_u*c26_slineu + abs_sc*(c13_intu*aloline_on_nn3d(iim1,j,k,q1) + &
                  c14_intu*aloline_on_nn3d(iim1,j,k,q2) + &
                  c15_intu*aloline_on_nn3d(iim1,j,k,q3) + &
                  c16_intu*aloline_on_nn3d(iim1,j,k,q4) + &
                  c17_intu*aloline_on_nn3d(iim1,j,k,q5) + &
                  c18_intu*aloline_on_nn3d(iim1,j,k,q6) + &
                  c22_intu*aloline_on_nn3d(iim1,j,k,q10) + &
                  c23_intu*aloline_on_nn3d(iim1,j,k,q11) + &
                  c24_intu*aloline_on_nn3d(iim1,j,k,q12) + &
                  c26_intu*aloline_on_nn3d(iim1,j,k,q14)))
               aloline_on_nn3d(iip1,jjm1,k,q16) = imask_totreg3d(iip1,jjm1,k) * &
                  (alo_d*c12_slined + abs_sc*(c12_intu*aloline_on_nn3d(iip1,jjm1,k,q1) + &
                  c15_intu*aloline_on_nn3d(iip1,jjm1,k,q4) + &
                  c18_intu*aloline_on_nn3d(iip1,jjm1,k,q7) + &
                  c24_intu*aloline_on_nn3d(iip1,jjm1,k,q13)))
               aloline_on_nn3d(i,jjm1,k,q17) = imask_totreg3d(i,jjm1,k) * &
                  (alo_u*c24_slineu + abs_sc*(c11_intu*aloline_on_nn3d(i,jjm1,k,q1) + &
                  c12_intu*aloline_on_nn3d(i,jjm1,k,q2) + &
                  c14_intu*aloline_on_nn3d(i,jjm1,k,q4) + &
                  c15_intu*aloline_on_nn3d(i,jjm1,k,q5) + &
                  c17_intu*aloline_on_nn3d(i,jjm1,k,q7) + &
                  c18_intu*aloline_on_nn3d(i,jjm1,k,q8) + &
                  c23_intu*aloline_on_nn3d(i,jjm1,k,q13) + &
                  c24_intu*aloline_on_nn3d(i,jjm1,k,q14) + &
                  c20_intu*aloline_on_nn3d(i,jjm1,k,q10) + &
                  c26_intu*aloline_on_nn3d(i,jjm1,k,q16)))
               aloline_on_nn3d(iim1,jjm1,k,q18) = imask_totreg3d(iim1,jjm1,k) * &
                  (alo_u*c23_slineu + abs_sc*(c10_intu*aloline_on_nn3d(iim1,jjm1,k,q1) + &
                  c11_intu*aloline_on_nn3d(iim1,jjm1,k,q2) + &
                  c12_intu*aloline_on_nn3d(iim1,jjm1,k,q3) + &
                  c13_intu*aloline_on_nn3d(iim1,jjm1,k,q4) + &
                  c14_intu*aloline_on_nn3d(iim1,jjm1,k,q5) + &
                  c15_intu*aloline_on_nn3d(iim1,jjm1,k,q6) + &
                  c16_intu*aloline_on_nn3d(iim1,jjm1,k,q7) + &
                  c17_intu*aloline_on_nn3d(iim1,jjm1,k,q8) + &
                  c18_intu*aloline_on_nn3d(iim1,jjm1,k,q9) + &
                  c22_intu*aloline_on_nn3d(iim1,jjm1,k,q13) + &
                  c23_intu*aloline_on_nn3d(iim1,jjm1,k,q14) + &
                  c24_intu*aloline_on_nn3d(iim1,jjm1,k,q15) + &
                  c20_intu*aloline_on_nn3d(iim1,jjm1,k,q11) + &
                  c26_intu*aloline_on_nn3d(iim1,jjm1,k,q17)))
               aloline_on_nn3d(iip1,jjp1,kkm1,q19) = imask_totreg3d(iip1,jjp1,kkm1) * &
                  (alo_d*c09_slined + abs_sc*c18_intu*aloline_on_nn3d(iip1,jjp1,kkm1,q10))
               aloline_on_nn3d(i,jjp1,kkm1,q20) = imask_totreg3d(i,jjp1,kkm1) * &
                  (alo_d*c08_slined + abs_sc*(c17_intu*aloline_on_nn3d(i,jjp1,kkm1,q10) + &
                  c18_intu*aloline_on_nn3d(i,jjp1,kkm1,q11) + &
                  c08_intu*aloline_on_nn3d(i,jjp1,kkm1,q1) + &
                  c26_intu*aloline_on_nn3d(i,jjp1,kkm1,q19)))
               aloline_on_nn3d(iim1,jjp1,kkm1,q21) = imask_totreg3d(iim1,jjp1,kkm1) * &
                  (alo_d*c07_slined + abs_sc*(c16_intu*aloline_on_nn3d(iim1,jjp1,kkm1,q10) + &
                  c17_intu*aloline_on_nn3d(iim1,jjp1,kkm1,q11) + &
                  c18_intu*aloline_on_nn3d(iim1,jjp1,kkm1,q12) + &
                  c08_intu*aloline_on_nn3d(iim1,jjp1,kkm1,q2) + &
                  c26_intu*aloline_on_nn3d(iim1,jjp1,kkm1,q20)))
               aloline_on_nn3d(iip1,j,kkm1,q22) = imask_totreg3d(iip1,j,kkm1) * &
                  (alo_d*c06_slined + abs_sc*(c15_intu*aloline_on_nn3d(iip1,j,kkm1,q10) + &
                  c18_intu*aloline_on_nn3d(iip1,j,kkm1,q13) + &
                  c06_intu*aloline_on_nn3d(iip1,j,kkm1,q1) + &
                  c24_intu*aloline_on_nn3d(iip1,j,kkm1,q19)))
               aloline_on_nn3d(i,j,kkm1,q23) = imask_totreg3d(i,j,kkm1) * &
                  (alo_u*c18_slineu + abs_sc*(c14_intu*aloline_on_nn3d(i,j,kkm1,q10) + &
                  c15_intu*aloline_on_nn3d(i,j,kkm1,q11) + &
                  c17_intu*aloline_on_nn3d(i,j,kkm1,q13) + &
                  c18_intu*aloline_on_nn3d(i,j,kkm1,q14) + &
                  c05_intu*aloline_on_nn3d(i,j,kkm1,q1) + &
                  c06_intu*aloline_on_nn3d(i,j,kkm1,q2) + &
                  c23_intu*aloline_on_nn3d(i,j,kkm1,q19) + &
                  c24_intu*aloline_on_nn3d(i,j,kkm1,q20) + &
                  c08_intu*aloline_on_nn3d(i,j,kkm1,q4) + &
                  c26_intu*aloline_on_nn3d(i,j,kkm1,q22)))
               aloline_on_nn3d(iim1,j,kkm1,q24) = imask_totreg3d(iim1,j,kkm1) * &
                  (alo_u*c17_slineu + abs_sc*(c13_intu*aloline_on_nn3d(iim1,j,kkm1,q10) + &
                  c14_intu*aloline_on_nn3d(iim1,j,kkm1,q11) + &
                  c15_intu*aloline_on_nn3d(iim1,j,kkm1,q12) + &
                  c16_intu*aloline_on_nn3d(iim1,j,kkm1,q13) + &
                  c17_intu*aloline_on_nn3d(iim1,j,kkm1,q14) + &
                  c18_intu*aloline_on_nn3d(iim1,j,kkm1,q15) + &
                  c04_intu*aloline_on_nn3d(iim1,j,kkm1,q1) + &
                  c05_intu*aloline_on_nn3d(iim1,j,kkm1,q2) + &
                  c06_intu*aloline_on_nn3d(iim1,j,kkm1,q3) + &
                  c22_intu*aloline_on_nn3d(iim1,j,kkm1,q19) + &
                  c23_intu*aloline_on_nn3d(iim1,j,kkm1,q20) + &
                  c24_intu*aloline_on_nn3d(iim1,j,kkm1,q21) + &
                  c08_intu*aloline_on_nn3d(iim1,j,kkm1,q5) + &
                  c26_intu*aloline_on_nn3d(iim1,j,kkm1,q23)))
               aloline_on_nn3d(iip1,jjm1,kkm1,q25) = imask_totreg3d(iip1,jjm1,kkm1) * &
                  (alo_d*c03_slined + abs_sc*(c12_intu*aloline_on_nn3d(iip1,jjm1,kkm1,q10) + &
                  c15_intu*aloline_on_nn3d(iip1,jjm1,kkm1,q13) + &
                  c18_intu*aloline_on_nn3d(iip1,jjm1,kkm1,q16) + &
                  c06_intu*aloline_on_nn3d(iip1,jjm1,kkm1,q4) + &
                  c24_intu*aloline_on_nn3d(iip1,jjm1,kkm1,q22)))
               aloline_on_nn3d(i,jjm1,kkm1,q26) = imask_totreg3d(i,jjm1,kkm1) * &
                  (alo_u*c15_slineu + abs_sc*(c11_intu*aloline_on_nn3d(i,jjm1,kkm1,q10) + &
                  c12_intu*aloline_on_nn3d(i,jjm1,kkm1,q11) + &
                  c14_intu*aloline_on_nn3d(i,jjm1,kkm1,q13) + &
                  c15_intu*aloline_on_nn3d(i,jjm1,kkm1,q14) + &
                  c17_intu*aloline_on_nn3d(i,jjm1,kkm1,q16) + &
                  c18_intu*aloline_on_nn3d(i,jjm1,kkm1,q17) + &
                  c05_intu*aloline_on_nn3d(i,jjm1,kkm1,q4) + &
                  c06_intu*aloline_on_nn3d(i,jjm1,kkm1,q5) + &
                  c23_intu*aloline_on_nn3d(i,jjm1,kkm1,q22) + &
                  c24_intu*aloline_on_nn3d(i,jjm1,kkm1,q23) + &
                  c02_intu*aloline_on_nn3d(i,jjm1,kkm1,q1) + &
                  c08_intu*aloline_on_nn3d(i,jjm1,kkm1,q7) + &
                  c20_intu*aloline_on_nn3d(i,jjm1,kkm1,q19) + &
                  c26_intu*aloline_on_nn3d(i,jjm1,kkm1,q25)))
               aloline_on_nn3d(iim1,jjm1,kkm1,q27) = imask_totreg3d(iim1,jjm1,kkm1) * &
                  (alo_u*c14_slineu + abs_sc*(c10_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q10) + &
                  c11_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q11) + &
                  c12_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q12) + &
                  c13_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q13) + &
                  c14_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q14) + &
                  c15_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q15) + &
                  c16_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q16) + &
                  c17_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q17) + &
                  c18_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q18) + &
                  c04_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q4) + &
                  c05_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q5) + &
                  c06_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q6) + &
                  c22_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q22) + &
                  c23_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q23) + &
                  c24_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q24) + &
                  c02_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q2) + &
                  c08_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q8) + &
                  c20_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q20) + &
                  c26_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q26)))

!
!perform angular integration
               call calc_phinorm(vel_p, vth_p, vth_fiducial, xobs, phinorm)
               wall=wall_global*phinorm
               mintbar3d_tmp(i,j,k) = mintbar3d_tmp(i,j,k) + int3d(i,j,k)*wall
               normalization3d_tmp(i,j,k) = normalization3d_tmp(i,j,k) + wall
!
               aloline_nn3d_tmp(iip1,jjp1,kkp1,q1) = aloline_nn3d_tmp(iip1,jjp1,kkp1,q1) + wall*aloline_on_nn3d(iip1,jjp1,kkp1,q1)
               aloline_nn3d_tmp(i,jjp1,kkp1,q2) = aloline_nn3d_tmp(i,jjp1,kkp1,q2) + wall*aloline_on_nn3d(i,jjp1,kkp1,q2)
               aloline_nn3d_tmp(iim1,jjp1,kkp1,q3) = aloline_nn3d_tmp(iim1,jjp1,kkp1,q3) + wall*aloline_on_nn3d(iim1,jjp1,kkp1,q3)
               aloline_nn3d_tmp(iip1,j,kkp1,q4) = aloline_nn3d_tmp(iip1,j,kkp1,q4) + wall*aloline_on_nn3d(iip1,j,kkp1,q4)
               aloline_nn3d_tmp(i,j,kkp1,q5) = aloline_nn3d_tmp(i,j,kkp1,q5) + wall*aloline_on_nn3d(i,j,kkp1,q5)
               aloline_nn3d_tmp(iim1,j,kkp1,q6) = aloline_nn3d_tmp(iim1,j,kkp1,q6) + wall*aloline_on_nn3d(iim1,j,kkp1,q6)
               aloline_nn3d_tmp(iip1,jjm1,kkp1,q7) = aloline_nn3d_tmp(iip1,jjm1,kkp1,q7) + wall*aloline_on_nn3d(iip1,jjm1,kkp1,q7)
               aloline_nn3d_tmp(i,jjm1,kkp1,q8) = aloline_nn3d_tmp(i,jjm1,kkp1,q8) + wall*aloline_on_nn3d(i,jjm1,kkp1,q8)
               aloline_nn3d_tmp(iim1,jjm1,kkp1,q9) = aloline_nn3d_tmp(iim1,jjm1,kkp1,q9) + wall*aloline_on_nn3d(iim1,jjm1,kkp1,q9)
               aloline_nn3d_tmp(iip1,jjp1,k,q10) = aloline_nn3d_tmp(iip1,jjp1,k,q10) + wall*aloline_on_nn3d(iip1,jjp1,k,q10)
               aloline_nn3d_tmp(i,jjp1,k,q11) = aloline_nn3d_tmp(i,jjp1,k,q11) + wall*aloline_on_nn3d(i,jjp1,k,q11)
               aloline_nn3d_tmp(iim1,jjp1,k,q12) = aloline_nn3d_tmp(iim1,jjp1,k,q12) + wall*aloline_on_nn3d(iim1,jjp1,k,q12)
               aloline_nn3d_tmp(iip1,j,k,q13) = aloline_nn3d_tmp(iip1,j,k,q13) + wall*aloline_on_nn3d(iip1,j,k,q13)
               aloline_nn3d_tmp(i,j,k,q14) = aloline_nn3d_tmp(i,j,k,q14) + wall*aloline_on_nn3d(i,j,k,q14)
               aloline_nn3d_tmp(iim1,j,k,q15) = aloline_nn3d_tmp(iim1,j,k,q15) + wall*aloline_on_nn3d(iim1,j,k,q15)
               aloline_nn3d_tmp(iip1,jjm1,k,q16) = aloline_nn3d_tmp(iip1,jjm1,k,q16) + wall*aloline_on_nn3d(iip1,jjm1,k,q16) !new
               aloline_nn3d_tmp(i,jjm1,k,q17) = aloline_nn3d_tmp(i,jjm1,k,q17) + wall*aloline_on_nn3d(i,jjm1,k,q17)
               aloline_nn3d_tmp(iim1,jjm1,k,q18) = aloline_nn3d_tmp(iim1,jjm1,k,q18) + wall*aloline_on_nn3d(iim1,jjm1,k,q18)
               aloline_nn3d_tmp(iip1,jjp1,kkm1,q19) = aloline_nn3d_tmp(iip1,jjp1,kkm1,q19) + wall*aloline_on_nn3d(iip1,jjp1,kkm1,q19) !new
               aloline_nn3d_tmp(i,jjp1,kkm1,q20) = aloline_nn3d_tmp(i,jjp1,kkm1,q20) + wall*aloline_on_nn3d(i,jjp1,kkm1,q20) !new
               aloline_nn3d_tmp(iim1,jjp1,kkm1,q21) = aloline_nn3d_tmp(iim1,jjp1,kkm1,q21) + wall*aloline_on_nn3d(iim1,jjp1,kkm1,q21) !new
               aloline_nn3d_tmp(iip1,j,kkm1,q22) = aloline_nn3d_tmp(iip1,j,kkm1,q22) + wall*aloline_on_nn3d(iip1,j,kkm1,q22) !new
               aloline_nn3d_tmp(i,j,kkm1,q23) = aloline_nn3d_tmp(i,j,kkm1,q23) + wall*aloline_on_nn3d(i,j,kkm1,q23)
               aloline_nn3d_tmp(iim1,j,kkm1,q24) = aloline_nn3d_tmp(iim1,j,kkm1,q24) + wall*aloline_on_nn3d(iim1,j,kkm1,q24)
               aloline_nn3d_tmp(iip1,jjm1,kkm1,q25) = aloline_nn3d_tmp(iip1,jjm1,kkm1,q25) + wall*aloline_on_nn3d(iip1,jjm1,kkm1,q25) !new
               aloline_nn3d_tmp(i,jjm1,kkm1,q26) = aloline_nn3d_tmp(i,jjm1,kkm1,q26) + wall*aloline_on_nn3d(i,jjm1,kkm1,q26)
               aloline_nn3d_tmp(iim1,jjm1,kkm1,q27) = aloline_nn3d_tmp(iim1,jjm1,kkm1,q27) + wall*aloline_on_nn3d(iim1,jjm1,kkm1,q27)
!
!***********special treatment: point is in vicinity of boundary*********
!
             case(2)
!
               iim1=i-alpha
               jjm1=j-beta
               kkm1=k-gamma
               iip1=i+alpha
               jjp1=j+beta
               kkp1=k+gamma
!
!calculate distance of upwind point to (previous) xy-plane, xz-plane and yz-plane
               dels_xyu=(z(k)-z(kkm1))/nn_z
               dels_xzu=(y(j)-y(jjm1))/nn_y
               dels_yzu=(x(i)-x(iim1))/nn_x
!
!calculate distance of downwind point to (next) xy-plane, xz-plane and yz-plane
               dels_xyd=(z(kkp1)-z(k))/nn_z
               dels_xzd=(y(jjp1)-y(j))/nn_y
               dels_yzd=(x(iip1)-x(i))/nn_x
!
!calculate distance to the boundary (here: sphere)
!               dels_r=dist_sphere(nn_x,nn_y,nn_z,x(i),y(j),z(k))
               dels_r=dist_ellipsoid(nn_x,nn_y,nn_z,x(i),y(j),z(k),smajorax_a,smajorax_b,smajorax_c)
!
               dels_u=min(dels_xyu,dels_xzu,dels_yzu,dels_r)
               dels_d=min(dels_xyd,dels_xzd,dels_yzd)
!
!----------------------------local point--------------------------------
!
               sline_p=sline3d(i,j,k)
               opalbar_p=opalbar3d(i,j,k)
               velx_p=velx3d(i,j,k)
               vely_p=vely3d(i,j,k)
               velz_p=velz3d(i,j,k)
               vel_p = nn_x*velx_p + nn_y*vely_p + nn_z*velz_p
               vth_p = vth3d(i,j,k)
!
!----------------------------upwind point-------------------------------
!
               if(dels_r.eq.dels_u) then
                  x_u = x(i) - dels_u*nn_x
                  y_u = y(j) - dels_u*nn_y
                  z_u = z(k) - dels_u*nn_z
!
!for alo
!                  int_u=0.d0
!                  int_u=xic1
                  int_u=calc_icore_gdark(z_u, sqrt(x_u**2+y_u**2+z_u**2))
!
                  call coeff3d_lineu_lin(iim1, i, jjm1, j, kkm1, k, x_u, y_u, z_u, &
                     c14_slineu, c15_slineu, c17_slineu, c18_slineu, &
                     c23_slineu, c24_slineu, c26_slineu, c27_slineu, &
                     opalbar_u, velx_u, vely_u, velz_u, vth_u, sline_u)
                  c02_slineu=0.d0
                  c04_slineu=0.d0
                  c05_slineu=0.d0
                  c06_slineu=0.d0
                  c08_slineu=0.d0
                  c10_slineu=0.d0
                  c11_slineu=0.d0
                  c12_slineu=0.d0
                  c13_slineu=0.d0
                  c16_slineu=0.d0
                  c20_slineu=0.d0
                  c22_slineu=0.d0
!set interpolation coefficients that are not used to zero
                  c02_intu = 0.d0
                  c04_intu = 0.d0
                  c05_intu = 0.d0
                  c06_intu = 0.d0
                  c08_intu = 0.d0
                  c10_intu = 0.d0
                  c11_intu = 0.d0
                  c12_intu = 0.d0
                  c13_intu = 0.d0
                  c14_intu = 0.d0
                  c15_intu = 0.d0
                  c16_intu = 0.d0
                  c17_intu = 0.d0
                  c18_intu = 0.d0
                  c20_intu = 0.d0
                  c22_intu = 0.d0
                  c23_intu = 0.d0
                  c24_intu = 0.d0
                  c26_intu = 0.d0
!
               elseif(dels_xyu.eq.dels_u) then
!intersection with x-y plane on level k-gamma
                  x_u = x(i) - dels_u*nn_x
                  y_u = y(j) - dels_u*nn_y
                  z_u = z(kkm1)
!
                  iim2=i-2*alpha
                  jjm2=j-2*beta
!
                  call coeff2d_lineu(opalbar3d(iim2,jjm2,kkm1), opalbar3d(iim1,jjm2,kkm1), opalbar3d(i,jjm2,kkm1), &
                     opalbar3d(iim2,jjm1,kkm1), opalbar3d(iim1,jjm1,kkm1), opalbar3d(i,jjm1,kkm1), &
                     opalbar3d(iim2,j,kkm1),    opalbar3d(iim1,j,kkm1),    opalbar3d(i,j,kkm1), &
                     velx3d(iim2,jjm2,kkm1), velx3d(iim1,jjm2,kkm1), velx3d(i,jjm2,kkm1), &
                     velx3d(iim2,jjm1,kkm1), velx3d(iim1,jjm1,kkm1), velx3d(i,jjm1,kkm1), &
                     velx3d(iim2,j,kkm1),    velx3d(iim1,j,kkm1),    velx3d(i,j,kkm1), &
                     vely3d(iim2,jjm2,kkm1), vely3d(iim1,jjm2,kkm1), vely3d(i,jjm2,kkm1), &
                     vely3d(iim2,jjm1,kkm1), vely3d(iim1,jjm1,kkm1), vely3d(i,jjm1,kkm1), &
                     vely3d(iim2,j,kkm1),    vely3d(iim1,j,kkm1),    vely3d(i,j,kkm1), &
                     velz3d(iim2,jjm2,kkm1), velz3d(iim1,jjm2,kkm1), velz3d(i,jjm2,kkm1), &
                     velz3d(iim2,jjm1,kkm1), velz3d(iim1,jjm1,kkm1), velz3d(i,jjm1,kkm1), &
                     velz3d(iim2,j,kkm1),    velz3d(iim1,j,kkm1),    velz3d(i,j,kkm1), &
                     vth3d(iim2,jjm2,kkm1), vth3d(iim1,jjm2,kkm1), vth3d(i,jjm2,kkm1), &
                     vth3d(iim2,jjm1,kkm1), vth3d(iim1,jjm1,kkm1), vth3d(i,jjm1,kkm1), &
                     vth3d(iim2,j,kkm1),    vth3d(iim1,j,kkm1),    vth3d(i,j,kkm1), &
                     sline3d(iim2,jjm2,kkm1), sline3d(iim1,jjm2,kkm1), sline3d(i,jjm2,kkm1), &
                     sline3d(iim2,jjm1,kkm1), sline3d(iim1,jjm1,kkm1), sline3d(i,jjm1,kkm1), &
                     sline3d(iim2,j,kkm1),    sline3d(iim1,j,kkm1),    sline3d(i,j,kkm1), &
                     int3d(iim2,jjm2,kkm1), int3d(iim1,jjm2,kkm1), int3d(i,jjm2,kkm1), &
                     int3d(iim2,jjm1,kkm1), int3d(iim1,jjm1,kkm1), int3d(i,jjm1,kkm1), &
                     int3d(iim2,j,kkm1),    int3d(iim1,j,kkm1),    int3d(i,j,kkm1), &
                     x(iim2), x(iim1), x(i), y(jjm2), y(jjm1), y(j), x_u, y_u, &
                     c10_slineu, c11_slineu, c12_slineu, c13_slineu, c14_slineu, &
                     c15_slineu, c16_slineu, c17_slineu, c18_slineu, &
                     c10_intu, c11_intu, c12_intu, c13_intu, c14_intu, &
                     c15_intu, c16_intu, c17_intu, c18_intu, &
                     opalbar_u, velx_u, vely_u, velz_u, vth_u, sline_u, int_u)
!set interpolation coefficients that are not used to zero
                  c23_slineu=0.d0
                  c24_slineu=0.d0
                  c26_slineu=0.d0
                  c27_slineu=0.d0
!
                  c02_intu = 0.d0
                  c04_intu = 0.d0
                  c05_intu = 0.d0
                  c06_intu = 0.d0
                  c08_intu = 0.d0
                  c20_intu = 0.d0
                  c22_intu = 0.d0
                  c23_intu = 0.d0
                  c24_intu = 0.d0
                  c26_intu = 0.d0
!
               elseif(dels_xzu.eq.dels_u) then
!intersection with x-z plane on level j-beta
                  x_u = x(i) - dels_u*nn_x
                  y_u = y(jjm1)
                  z_u = z(k) - dels_u*nn_z
!
                  iim2=i-2*alpha
                  kkm2=k-2*gamma
!
                  call coeff2d_lineu(opalbar3d(iim2,jjm1,kkm2), opalbar3d(iim1,jjm1,kkm2), opalbar3d(i,jjm1,kkm2), &
                     opalbar3d(iim2,jjm1,kkm1), opalbar3d(iim1,jjm1,kkm1), opalbar3d(i,jjm1,kkm1), &
                     opalbar3d(iim2,jjm1,k), opalbar3d(iim1,jjm1,k), opalbar3d(i,jjm1,k), &
                     velx3d(iim2,jjm1,kkm2), velx3d(iim1,jjm1,kkm2), velx3d(i,jjm1,kkm2), &
                     velx3d(iim2,jjm1,kkm1), velx3d(iim1,jjm1,kkm1), velx3d(i,jjm1,kkm1), &
                     velx3d(iim2,jjm1,k), velx3d(iim1,jjm1,k), velx3d(i,jjm1,k), &
                     vely3d(iim2,jjm1,kkm2), vely3d(iim1,jjm1,kkm2), vely3d(i,jjm1,kkm2), &
                     vely3d(iim2,jjm1,kkm1), vely3d(iim1,jjm1,kkm1), vely3d(i,jjm1,kkm1), &
                     vely3d(iim2,jjm1,k), vely3d(iim1,jjm1,k), vely3d(i,jjm1,k), &
                     velz3d(iim2,jjm1,kkm2), velz3d(iim1,jjm1,kkm2), velz3d(i,jjm1,kkm2), &
                     velz3d(iim2,jjm1,kkm1), velz3d(iim1,jjm1,kkm1), velz3d(i,jjm1,kkm1), &
                     velz3d(iim2,jjm1,k), velz3d(iim1,jjm1,k), velz3d(i,jjm1,k), &
                     vth3d(iim2,jjm1,kkm2), vth3d(iim1,jjm1,kkm2), vth3d(i,jjm1,kkm2), &
                     vth3d(iim2,jjm1,kkm1), vth3d(iim1,jjm1,kkm1), vth3d(i,jjm1,kkm1), &
                     vth3d(iim2,jjm1,k), vth3d(iim1,jjm1,k), vth3d(i,jjm1,k), &
                     sline3d(iim2,jjm1,kkm2), sline3d(iim1,jjm1,kkm2), sline3d(i,jjm1,kkm2), &
                     sline3d(iim2,jjm1,kkm1), sline3d(iim1,jjm1,kkm1), sline3d(i,jjm1,kkm1), &
                     sline3d(iim2,jjm1,k), sline3d(iim1,jjm1,k), sline3d(i,jjm1,k), &
                     int3d(iim2,jjm1,kkm2), int3d(iim1,jjm1,kkm2), int3d(i,jjm1,kkm2), &
                     int3d(iim2,jjm1,kkm1), int3d(iim1,jjm1,kkm1), int3d(i,jjm1,kkm1), &
                     int3d(iim2,jjm1,k), int3d(iim1,jjm1,k), int3d(i,jjm1,k), &
                     x(iim2), x(iim1), x(i), z(kkm2), z(kkm1), z(k), x_u, z_u, &
                     c04_slineu, c05_slineu, c06_slineu, c13_slineu, c14_slineu, &
                     c15_slineu, c22_slineu, c23_slineu, c24_slineu, &
                     c04_intu, c05_intu, c06_intu, c13_intu, c14_intu, &
                     c15_intu, c22_intu, c23_intu, c24_intu, &
                     opalbar_u, velx_u, vely_u, velz_u, vth_u, sline_u, int_u)
!set interpolation coefficients that are not used to zero
                  c17_slineu=0.d0
                  c18_slineu=0.d0
                  c26_slineu=0.d0
                  c27_slineu=0.d0
                  c02_intu = 0.d0
                  c08_intu = 0.d0
                  c10_intu = 0.d0
                  c11_intu = 0.d0
                  c12_intu = 0.d0
                  c16_intu = 0.d0
                  c17_intu = 0.d0
                  c18_intu = 0.d0
                  c20_intu = 0.d0
                  c26_intu = 0.d0
!
               elseif(dels_yzu.eq.dels_u) then
!intersection with y-z plane on level i-alpha
                  x_u = x(iim1)
                  y_u = y(j) - dels_u*nn_y
                  z_u = z(k) - dels_u*nn_z
!
                  jjm2=j-2*beta
                  kkm2=k-2*gamma
!
                  call coeff2d_lineu(opalbar3d(iim1,jjm2,kkm2), opalbar3d(iim1,jjm1,kkm2), opalbar3d(iim1,j,kkm2), &
                     opalbar3d(iim1,jjm2,kkm1), opalbar3d(iim1,jjm1,kkm1), opalbar3d(iim1,j,kkm1), &
                     opalbar3d(iim1,jjm2,k), opalbar3d(iim1,jjm1,k), opalbar3d(iim1,j,k), &
                     velx3d(iim1,jjm2,kkm2), velx3d(iim1,jjm1,kkm2), velx3d(iim1,j,kkm2), &
                     velx3d(iim1,jjm2,kkm1), velx3d(iim1,jjm1,kkm1), velx3d(iim1,j,kkm1), &
                     velx3d(iim1,jjm2,k), velx3d(iim1,jjm1,k), velx3d(iim1,j,k), &
                     vely3d(iim1,jjm2,kkm2), vely3d(iim1,jjm1,kkm2), vely3d(iim1,j,kkm2), &
                     vely3d(iim1,jjm2,kkm1), vely3d(iim1,jjm1,kkm1), vely3d(iim1,j,kkm1), &
                     vely3d(iim1,jjm2,k), vely3d(iim1,jjm1,k), vely3d(iim1,j,k), &
                     velz3d(iim1,jjm2,kkm2), velz3d(iim1,jjm1,kkm2), velz3d(iim1,j,kkm2), &
                     velz3d(iim1,jjm2,kkm1), velz3d(iim1,jjm1,kkm1), velz3d(iim1,j,kkm1), &
                     velz3d(iim1,jjm2,k), velz3d(iim1,jjm1,k), velz3d(iim1,j,k), &
                     vth3d(iim1,jjm2,kkm2), vth3d(iim1,jjm1,kkm2), vth3d(iim1,j,kkm2), &
                     vth3d(iim1,jjm2,kkm1), vth3d(iim1,jjm1,kkm1), vth3d(iim1,j,kkm1), &
                     vth3d(iim1,jjm2,k), vth3d(iim1,jjm1,k), vth3d(iim1,j,k), &
                     sline3d(iim1,jjm2,kkm2), sline3d(iim1,jjm1,kkm2), sline3d(iim1,j,kkm2), &
                     sline3d(iim1,jjm2,kkm1), sline3d(iim1,jjm1,kkm1), sline3d(iim1,j,kkm1), &
                     sline3d(iim1,jjm2,k), sline3d(iim1,jjm1,k), sline3d(iim1,j,k), &
                     int3d(iim1,jjm2,kkm2), int3d(iim1,jjm1,kkm2), int3d(iim1,j,kkm2), &
                     int3d(iim1,jjm2,kkm1), int3d(iim1,jjm1,kkm1), int3d(iim1,j,kkm1), &
                     int3d(iim1,jjm2,k), int3d(iim1,jjm1,k), int3d(iim1,j,k), &
                     y(jjm2), y(jjm1), y(j), z(kkm2), z(kkm1), z(k), y_u, z_u, &
                     c02_slineu, c05_slineu, c08_slineu, c11_slineu, c14_slineu, &
                     c17_slineu, c20_slineu, c23_slineu, c26_slineu, &
                     c02_intu, c05_intu, c08_intu, c11_intu, c14_intu, &
                     c17_intu, c20_intu, c23_intu, c26_intu, &
                     opalbar_u, velx_u, vely_u, velz_u, vth_u, sline_u, int_u)
!set interpolation coefficients that are not used to zero
                  c15_slineu=0.d0
                  c18_slineu=0.d0
                  c24_slineu=0.d0
                  c27_slineu=0.d0
!
                  c04_intu = 0.d0
                  c06_intu = 0.d0
                  c10_intu = 0.d0
                  c12_intu = 0.d0
                  c13_intu = 0.d0
                  c15_intu = 0.d0
                  c16_intu = 0.d0
                  c18_intu = 0.d0
                  c22_intu = 0.d0
                  c24_intu = 0.d0
!
               else
                  write(*,'(4es20.8,2l4)') dels_u, dels_xzu, dels_xyu, dels_yzu, dels_u.eq.dels_xzu, dels_xyu.eq.dels_u
                  stop 'error in fsc_line3d: invalid dels_u'
               endif
!
!---------------------------downwind point------------------------------
!
               if(dels_xyd.eq.dels_d) then
!intersection with x-y plane on level k+gamma
                  x_d = x(i) + dels_d*nn_x
                  y_d = y(j) + dels_d*nn_y
                  z_d = z(kkp1)
!
                  call coeff2d_lined(opalbar3d(iim1,jjm1,kkp1), opalbar3d(i,jjm1,kkp1), opalbar3d(iip1,jjm1,kkp1), &
                     opalbar3d(iim1,j,kkp1), opalbar3d(i,j,kkp1), opalbar3d(iip1,j,kkp1), &
                     opalbar3d(iim1,jjp1,kkp1), opalbar3d(i,jjp1,kkp1), opalbar3d(iip1,jjp1,kkp1), &
                     velx3d(iim1,jjm1,kkp1), velx3d(i,jjm1,kkp1), velx3d(iip1,jjm1,kkp1), &
                     velx3d(iim1,j,kkp1), velx3d(i,j,kkp1), velx3d(iip1,j,kkp1), &
                     velx3d(iim1,jjp1,kkp1), velx3d(i,jjp1,kkp1), velx3d(iip1,jjp1,kkp1), &
                     vely3d(iim1,jjm1,kkp1), vely3d(i,jjm1,kkp1), vely3d(iip1,jjm1,kkp1), &
                     vely3d(iim1,j,kkp1), vely3d(i,j,kkp1), vely3d(iip1,j,kkp1), &
                     vely3d(iim1,jjp1,kkp1), vely3d(i,jjp1,kkp1), vely3d(iip1,jjp1,kkp1), &
                     velz3d(iim1,jjm1,kkp1), velz3d(i,jjm1,kkp1), velz3d(iip1,jjm1,kkp1), &
                     velz3d(iim1,j,kkp1), velz3d(i,j,kkp1), velz3d(iip1,j,kkp1), &
                     velz3d(iim1,jjp1,kkp1), velz3d(i,jjp1,kkp1), velz3d(iip1,jjp1,kkp1), &
                     vth3d(iim1,jjm1,kkp1), vth3d(i,jjm1,kkp1), vth3d(iip1,jjm1,kkp1), &
                     vth3d(iim1,j,kkp1), vth3d(i,j,kkp1), vth3d(iip1,j,kkp1), &
                     vth3d(iim1,jjp1,kkp1), vth3d(i,jjp1,kkp1), vth3d(iip1,jjp1,kkp1), &
                     sline3d(iim1,jjm1,kkp1), sline3d(i,jjm1,kkp1), sline3d(iip1,jjm1,kkp1), &
                     sline3d(iim1,j,kkp1), sline3d(i,j,kkp1), sline3d(iip1,j,kkp1), &
                     sline3d(iim1,jjp1,kkp1), sline3d(i,jjp1,kkp1), sline3d(iip1,jjp1,kkp1), &
                     x(iim1), x(i), x(iip1), y(jjm1), y(j), y(jjp1), x_d, y_d, &
                     c19_slined, c20_slined, c21_slined, c22_slined, c23_slined, c24_slined, &
                     c25_slined, c26_slined, c27_slined, &
                     opalbar_d, velx_d, vely_d, velz_d, vth_d, sline_d)
!
!set interpolation coefficients that are not used to zero
                  c03_slined = 0.d0
                  c06_slined = 0.d0
                  c07_slined = 0.d0
                  c08_slined = 0.d0
                  c09_slined = 0.d0
                  c12_slined = 0.d0
                  c15_slined = 0.d0
                  c16_slined = 0.d0
                  c17_slined = 0.d0
                  c18_slined = 0.d0
!
               elseif(dels_xzd.eq.dels_d) then
!intersection with x-z plane on level j+beta
                  x_d = x(i) + dels_d*nn_x
                  y_d = y(jjp1)
                  z_d = z(k) + dels_d*nn_z
!
                  call coeff2d_lined(opalbar3d(iim1,jjp1,kkm1), opalbar3d(i,jjp1,kkm1), opalbar3d(iip1,jjp1,kkm1), &
                     opalbar3d(iim1,jjp1,k), opalbar3d(i,jjp1,k), opalbar3d(iip1,jjp1,k), &
                     opalbar3d(iim1,jjp1,kkp1), opalbar3d(i,jjp1,kkp1), opalbar3d(iip1,jjp1,kkp1), &
                     velx3d(iim1,jjp1,kkm1), velx3d(i,jjp1,kkm1), velx3d(iip1,jjp1,kkm1), &
                     velx3d(iim1,jjp1,k), velx3d(i,jjp1,k), velx3d(iip1,jjp1,k), &
                     velx3d(iim1,jjp1,kkp1), velx3d(i,jjp1,kkp1), velx3d(iip1,jjp1,kkp1), &
                     vely3d(iim1,jjp1,kkm1), vely3d(i,jjp1,kkm1), vely3d(iip1,jjp1,kkm1), &
                     vely3d(iim1,jjp1,k), vely3d(i,jjp1,k), vely3d(iip1,jjp1,k), &
                     vely3d(iim1,jjp1,kkp1), vely3d(i,jjp1,kkp1), vely3d(iip1,jjp1,kkp1), &
                     velz3d(iim1,jjp1,kkm1), velz3d(i,jjp1,kkm1), velz3d(iip1,jjp1,kkm1), &
                     velz3d(iim1,jjp1,k), velz3d(i,jjp1,k), velz3d(iip1,jjp1,k), &
                     velz3d(iim1,jjp1,kkp1), velz3d(i,jjp1,kkp1), velz3d(iip1,jjp1,kkp1), &
                     vth3d(iim1,jjp1,kkm1), vth3d(i,jjp1,kkm1), vth3d(iip1,jjp1,kkm1), &
                     vth3d(iim1,jjp1,k), vth3d(i,jjp1,k), vth3d(iip1,jjp1,k), &
                     vth3d(iim1,jjp1,kkp1), vth3d(i,jjp1,kkp1), vth3d(iip1,jjp1,kkp1), &
                     sline3d(iim1,jjp1,kkm1), sline3d(i,jjp1,kkm1), sline3d(iip1,jjp1,kkm1), &
                     sline3d(iim1,jjp1,k), sline3d(i,jjp1,k), sline3d(iip1,jjp1,k), &
                     sline3d(iim1,jjp1,kkp1), sline3d(i,jjp1,kkp1), sline3d(iip1,jjp1,kkp1), &
                     x(iim1), x(i), x(iip1), z(kkm1), z(k), z(kkp1), x_d, z_d, &
                     c07_slined, c08_slined, c09_slined, c16_slined, c17_slined, c18_slined, &
                     c25_slined, c26_slined, c27_slined, &
                     opalbar_d, velx_d, vely_d, velz_d, vth_d, sline_d)
!
!set interpolation coefficients that are not used to zero
                  c03_slined = 0.d0
                  c06_slined = 0.d0
                  c12_slined = 0.d0
                  c15_slined = 0.d0
                  c19_slined = 0.d0
                  c20_slined = 0.d0
                  c21_slined = 0.d0
                  c22_slined = 0.d0
                  c23_slined = 0.d0
                  c24_slined = 0.d0
!
               elseif(dels_yzd.eq.dels_d) then
!intersection with y-z plane on level i+alpha
                  x_d = x(iip1)
                  y_d = y(j) + dels_d*nn_y
                  z_d = z(k) + dels_d*nn_z
!
                  call coeff2d_lined(opalbar3d(iip1,jjm1,kkm1), opalbar3d(iip1,j,kkm1), opalbar3d(iip1,jjp1,kkm1), &
                     opalbar3d(iip1,jjm1,k), opalbar3d(iip1,j,k), opalbar3d(iip1,jjp1,k), &
                     opalbar3d(iip1,jjm1,kkp1), opalbar3d(iip1,j,kkp1), opalbar3d(iip1,jjp1,kkp1), &
                     velx3d(iip1,jjm1,kkm1), velx3d(iip1,j,kkm1), velx3d(iip1,jjp1,kkm1), &
                     velx3d(iip1,jjm1,k), velx3d(iip1,j,k), velx3d(iip1,jjp1,k), &
                     velx3d(iip1,jjm1,kkp1), velx3d(iip1,j,kkp1), velx3d(iip1,jjp1,kkp1), &
                     vely3d(iip1,jjm1,kkm1), vely3d(iip1,j,kkm1), vely3d(iip1,jjp1,kkm1), &
                     vely3d(iip1,jjm1,k), vely3d(iip1,j,k), vely3d(iip1,jjp1,k), &
                     vely3d(iip1,jjm1,kkp1), vely3d(iip1,j,kkp1), vely3d(iip1,jjp1,kkp1), &
                     velz3d(iip1,jjm1,kkm1), velz3d(iip1,j,kkm1), velz3d(iip1,jjp1,kkm1), &
                     velz3d(iip1,jjm1,k), velz3d(iip1,j,k), velz3d(iip1,jjp1,k), &
                     velz3d(iip1,jjm1,kkp1), velz3d(iip1,j,kkp1), velz3d(iip1,jjp1,kkp1), &
                     vth3d(iip1,jjm1,kkm1), vth3d(iip1,j,kkm1), vth3d(iip1,jjp1,kkm1), &
                     vth3d(iip1,jjm1,k), vth3d(iip1,j,k), vth3d(iip1,jjp1,k), &
                     vth3d(iip1,jjm1,kkp1), vth3d(iip1,j,kkp1), vth3d(iip1,jjp1,kkp1), &
                     sline3d(iip1,jjm1,kkm1), sline3d(iip1,j,kkm1), sline3d(iip1,jjp1,kkm1), &
                     sline3d(iip1,jjm1,k), sline3d(iip1,j,k), sline3d(iip1,jjp1,k), &
                     sline3d(iip1,jjm1,kkp1), sline3d(iip1,j,kkp1), sline3d(iip1,jjp1,kkp1), &
                     y(jjm1), y(j), y(jjp1), z(kkm1), z(k), z(kkp1), y_d, z_d, &
                     c03_slined, c06_slined, c09_slined, c12_slined, c15_slined, c18_slined, &
                     c21_slined, c24_slined, c27_slined, &
                     opalbar_d, velx_d, vely_d, velz_d, vth_d, sline_d)
!
!set interpolation coefficients that are not used to zero
                  c07_slined = 0.d0
                  c08_slined = 0.d0
                  c16_slined = 0.d0
                  c17_slined = 0.d0
                  c19_slined = 0.d0
                  c20_slined = 0.d0
                  c22_slined = 0.d0
                  c23_slined = 0.d0
                  c25_slined = 0.d0
                  c26_slined = 0.d0
!
               else
                  write(*,'(4es20.8,l4)') dels_d, dels_xzd, dels_xyd, dels_yzd, dels_d.eq.dels_xzd
                  stop 'error in fsc_line3d: invalid dels_d'
               endif
!
!--------------------------------radiative transfer---------------------
!
!***debug start: analytical laws
!               call model_debug(x_u, y_u, z_u, velx_u, vely_u, velz_u, opalbar_u, fdum)
!               call model_debug(x_d, y_d, z_d, velx_d, vely_d, velz_d, opalbar_d, fdum)
!***debug end
               vel_u = velx_u*nn_x + vely_u*nn_y + velz_u*nn_z
               vel_d = velx_d*nn_x + vely_d*nn_y + velz_d*nn_z
               call fsc_line(xobs, vth_fiducial, deltax, xcmf_min, xcmf_max, int_u, opalbar_u, opalbar_p, opalbar_d, &
                  sline_u, sline_p, sline_d, vel_u, vel_p, vel_d, vth_u, vth_p, vth_d, &
                  dels_u, dels_d, abs_sc, contr_sc, int_sc, alo_u, alo_p, alo_d)
!               call fsc_line_lin(xobs, vth_fiducial, deltax, xcmf_min, xcmf_max, int_u, opalbar_u, opalbar_p, sline_u, sline_p, &
!                                   vel_u, vel_p, vth_u, vth_p, dels_u, abs_sc, contr_sc, int_sc, alo_u, alo_p)
!               alo_d=0.d0
               int3d(i,j,k) = abs_sc*int_u + alo_u*sline_u + alo_p*sline_p + alo_d*sline_d
!
               aloline_on_nn3d(iip1,jjp1,kkp1,q1) = imask_totreg3d(iip1,jjp1,kkp1)*&
                  alo_d*c27_slined
               aloline_on_nn3d(i,jjp1,kkp1,q2) = imask_totreg3d(i,jjp1,kkp1) * &
                  (alo_d*c26_slined + abs_sc*c26_intu*aloline_on_nn3d(i,jjp1,kkp1,q1))
               aloline_on_nn3d(iim1,jjp1,kkp1,q3) = imask_totreg3d(iim1,jjp1,kkp1) * &
                  (alo_d*c25_slined + abs_sc*c26_intu*aloline_on_nn3d(iim1,jjp1,kkp1,q2))
               aloline_on_nn3d(iip1,j,kkp1,q4) = imask_totreg3d(iip1,j,kkp1) * &
                  (alo_d*c24_slined + abs_sc*c24_intu*aloline_on_nn3d(iip1,j,kkp1,q1))
               aloline_on_nn3d(i,j,kkp1,q5) = imask_totreg3d(i,j,kkp1) * &
                  (alo_d*c23_slined + abs_sc*(c23_intu*aloline_on_nn3d(i,j,kkp1,q1) + &
                  c24_intu*aloline_on_nn3d(i,j,kkp1,q2) + &
                  c26_intu*aloline_on_nn3d(i,j,kkp1,q4)))
               aloline_on_nn3d(iim1,j,kkp1,q6) = imask_totreg3d(iim1,j,kkp1) * &
                  (alo_d*c22_slined + abs_sc*(c22_intu*aloline_on_nn3d(iim1,j,kkp1,q1) + &
                  c23_intu*aloline_on_nn3d(iim1,j,kkp1,q2) + &
                  c24_intu*aloline_on_nn3d(iim1,j,kkp1,q3) + &
                  c26_intu*aloline_on_nn3d(iim1,j,kkp1,q5)))
               aloline_on_nn3d(iip1,jjm1,kkp1,q7) = imask_totreg3d(iip1,jjm1,kkp1) * &
                  (alo_d*c21_slined + abs_sc*c24_intu*aloline_on_nn3d(iip1,jjm1,kkp1,q4))
               aloline_on_nn3d(i,jjm1,kkp1,q8) = imask_totreg3d(i,jjm1,kkp1) * &
                  (alo_d*c20_slined + abs_sc*(c23_intu*aloline_on_nn3d(i,jjm1,kkp1,q4) + &
                  c24_intu*aloline_on_nn3d(i,jjm1,kkp1,q5) + &
                  c20_intu*aloline_on_nn3d(i,jjm1,kkp1,q1) + &
                  c26_intu*aloline_on_nn3d(i,jjm1,kkp1,q7)))
               aloline_on_nn3d(iim1,jjm1,kkp1,q9) = imask_totreg3d(iim1,jjm1,kkp1) * &
                  (alo_d*c19_slined + abs_sc*(c22_intu*aloline_on_nn3d(iim1,jjm1,kkp1,q4) + &
                  c23_intu*aloline_on_nn3d(iim1,jjm1,kkp1,q5) + &
                  c24_intu*aloline_on_nn3d(iim1,jjm1,kkp1,q6) + &
                  c20_intu*aloline_on_nn3d(iim1,jjm1,kkp1,q2) + &
                  c26_intu*aloline_on_nn3d(iim1,jjm1,kkp1,q8)))
               aloline_on_nn3d(iip1,jjp1,k,q10) = imask_totreg3d(iip1,jjp1,k) * &
                  (alo_d*c18_slined + abs_sc*c18_intu*aloline_on_nn3d(iip1,jjp1,k,q1))
               aloline_on_nn3d(i,jjp1,k,q11) = imask_totreg3d(i,jjp1,k) * &
                  (alo_d*c17_slined + abs_sc*(c17_intu*aloline_on_nn3d(i,jjp1,k,q1) + &
                  c18_intu*aloline_on_nn3d(i,jjp1,k,q2) + &
                  c26_intu*aloline_on_nn3d(i,jjp1,k,q10)))
               aloline_on_nn3d(iim1,jjp1,k,q12) = imask_totreg3d(iim1,jjp1,k) * &
                  (alo_d*c16_slined + abs_sc*(c16_intu*aloline_on_nn3d(iim1,jjp1,k,q1) + &
                  c17_intu*aloline_on_nn3d(iim1,jjp1,k,q2) + &
                  c18_intu*aloline_on_nn3d(iim1,jjp1,k,q3) + &
                  c26_intu*aloline_on_nn3d(iim1,jjp1,k,q11)))
               aloline_on_nn3d(iip1,j,k,q13) = imask_totreg3d(iip1,j,k) * &
                  (alo_d*c15_slined + abs_sc*(c15_intu*aloline_on_nn3d(iip1,j,k,q1) + &
                  c18_intu*aloline_on_nn3d(iip1,j,k,q4) + &
                  c24_intu*aloline_on_nn3d(iip1,j,k,q10)))
               aloline_on_nn3d(i,j,k,q14) = imask_totreg3d(i,j,k) * &
                  (alo_p + abs_sc*(c14_intu*aloline_on_nn3d(i,j,k,q1) + &
                  c15_intu*aloline_on_nn3d(i,j,k,q2) + &
                  c17_intu*aloline_on_nn3d(i,j,k,q4) + &
                  c26_intu*aloline_on_nn3d(i,j,k,q13) + &
                  c18_intu*aloline_on_nn3d(i,j,k,q5) + &
                  c23_intu*aloline_on_nn3d(i,j,k,q10) + &
                  c24_intu*aloline_on_nn3d(i,j,k,q11)))
               aloline_on_nn3d(iim1,j,k,q15) = imask_totreg3d(iim1,j,k) * &
                  (alo_u*c26_slineu + abs_sc*(c13_intu*aloline_on_nn3d(iim1,j,k,q1) + &
                  c14_intu*aloline_on_nn3d(iim1,j,k,q2) + &
                  c15_intu*aloline_on_nn3d(iim1,j,k,q3) + &
                  c16_intu*aloline_on_nn3d(iim1,j,k,q4) + &
                  c17_intu*aloline_on_nn3d(iim1,j,k,q5) + &
                  c18_intu*aloline_on_nn3d(iim1,j,k,q6) + &
                  c22_intu*aloline_on_nn3d(iim1,j,k,q10) + &
                  c23_intu*aloline_on_nn3d(iim1,j,k,q11) + &
                  c24_intu*aloline_on_nn3d(iim1,j,k,q12) + &
                  c26_intu*aloline_on_nn3d(iim1,j,k,q14)))
               aloline_on_nn3d(iip1,jjm1,k,q16) = imask_totreg3d(iip1,jjm1,k) * &
                  (alo_d*c12_slined + abs_sc*(c12_intu*aloline_on_nn3d(iip1,jjm1,k,q1) + &
                  c15_intu*aloline_on_nn3d(iip1,jjm1,k,q4) + &
                  c18_intu*aloline_on_nn3d(iip1,jjm1,k,q7) + &
                  c24_intu*aloline_on_nn3d(iip1,jjm1,k,q13)))
               aloline_on_nn3d(i,jjm1,k,q17) = imask_totreg3d(i,jjm1,k) * &
                  (alo_u*c24_slineu + abs_sc*(c11_intu*aloline_on_nn3d(i,jjm1,k,q1) + &
                  c12_intu*aloline_on_nn3d(i,jjm1,k,q2) + &
                  c14_intu*aloline_on_nn3d(i,jjm1,k,q4) + &
                  c15_intu*aloline_on_nn3d(i,jjm1,k,q5) + &
                  c17_intu*aloline_on_nn3d(i,jjm1,k,q7) + &
                  c18_intu*aloline_on_nn3d(i,jjm1,k,q8) + &
                  c23_intu*aloline_on_nn3d(i,jjm1,k,q13) + &
                  c24_intu*aloline_on_nn3d(i,jjm1,k,q14) + &
                  c20_intu*aloline_on_nn3d(i,jjm1,k,q10) + &
                  c26_intu*aloline_on_nn3d(i,jjm1,k,q16)))
               aloline_on_nn3d(iim1,jjm1,k,q18) = imask_totreg3d(iim1,jjm1,k) * &
                  (alo_u*c23_slineu + abs_sc*(c10_intu*aloline_on_nn3d(iim1,jjm1,k,q1) + &
                  c11_intu*aloline_on_nn3d(iim1,jjm1,k,q2) + &
                  c12_intu*aloline_on_nn3d(iim1,jjm1,k,q3) + &
                  c13_intu*aloline_on_nn3d(iim1,jjm1,k,q4) + &
                  c14_intu*aloline_on_nn3d(iim1,jjm1,k,q5) + &
                  c15_intu*aloline_on_nn3d(iim1,jjm1,k,q6) + &
                  c16_intu*aloline_on_nn3d(iim1,jjm1,k,q7) + &
                  c17_intu*aloline_on_nn3d(iim1,jjm1,k,q8) + &
                  c18_intu*aloline_on_nn3d(iim1,jjm1,k,q9) + &
                  c22_intu*aloline_on_nn3d(iim1,jjm1,k,q13) + &
                  c23_intu*aloline_on_nn3d(iim1,jjm1,k,q14) + &
                  c24_intu*aloline_on_nn3d(iim1,jjm1,k,q15) + &
                  c20_intu*aloline_on_nn3d(iim1,jjm1,k,q11) + &
                  c26_intu*aloline_on_nn3d(iim1,jjm1,k,q17)))
               aloline_on_nn3d(iip1,jjp1,kkm1,q19) = imask_totreg3d(iip1,jjp1,kkm1) * &
                  (alo_d*c09_slined + abs_sc*c18_intu*aloline_on_nn3d(iip1,jjp1,kkm1,q10))
               aloline_on_nn3d(i,jjp1,kkm1,q20) = imask_totreg3d(i,jjp1,kkm1) * &
                  (alo_d*c08_slined + abs_sc*(c17_intu*aloline_on_nn3d(i,jjp1,kkm1,q10) + &
                  c18_intu*aloline_on_nn3d(i,jjp1,kkm1,q11) + &
                  c08_intu*aloline_on_nn3d(i,jjp1,kkm1,q1) + &
                  c26_intu*aloline_on_nn3d(i,jjp1,kkm1,q19)))
               aloline_on_nn3d(iim1,jjp1,kkm1,q21) = imask_totreg3d(iim1,jjp1,kkm1) * &
                  (alo_d*c07_slined + abs_sc*(c16_intu*aloline_on_nn3d(iim1,jjp1,kkm1,q10) + &
                  c17_intu*aloline_on_nn3d(iim1,jjp1,kkm1,q11) + &
                  c18_intu*aloline_on_nn3d(iim1,jjp1,kkm1,q12) + &
                  c08_intu*aloline_on_nn3d(iim1,jjp1,kkm1,q2) + &
                  c26_intu*aloline_on_nn3d(iim1,jjp1,kkm1,q20)))
               aloline_on_nn3d(iip1,j,kkm1,q22) = imask_totreg3d(iip1,j,kkm1) * &
                  (alo_d*c06_slined + abs_sc*(c15_intu*aloline_on_nn3d(iip1,j,kkm1,q10) + &
                  c18_intu*aloline_on_nn3d(iip1,j,kkm1,q13) + &
                  c06_intu*aloline_on_nn3d(iip1,j,kkm1,q1) + &
                  c24_intu*aloline_on_nn3d(iip1,j,kkm1,q19)))
               aloline_on_nn3d(i,j,kkm1,q23) = imask_totreg3d(i,j,kkm1) * &
                  (alo_u*c18_slineu + abs_sc*(c14_intu*aloline_on_nn3d(i,j,kkm1,q10) + &
                  c15_intu*aloline_on_nn3d(i,j,kkm1,q11) + &
                  c17_intu*aloline_on_nn3d(i,j,kkm1,q13) + &
                  c18_intu*aloline_on_nn3d(i,j,kkm1,q14) + &
                  c05_intu*aloline_on_nn3d(i,j,kkm1,q1) + &
                  c06_intu*aloline_on_nn3d(i,j,kkm1,q2) + &
                  c23_intu*aloline_on_nn3d(i,j,kkm1,q19) + &
                  c24_intu*aloline_on_nn3d(i,j,kkm1,q20) + &
                  c08_intu*aloline_on_nn3d(i,j,kkm1,q4) + &
                  c26_intu*aloline_on_nn3d(i,j,kkm1,q22)))
               aloline_on_nn3d(iim1,j,kkm1,q24) = imask_totreg3d(iim1,j,kkm1) * &
                  (alo_u*c17_slineu + abs_sc*(c13_intu*aloline_on_nn3d(iim1,j,kkm1,q10) + &
                  c14_intu*aloline_on_nn3d(iim1,j,kkm1,q11) + &
                  c15_intu*aloline_on_nn3d(iim1,j,kkm1,q12) + &
                  c16_intu*aloline_on_nn3d(iim1,j,kkm1,q13) + &
                  c17_intu*aloline_on_nn3d(iim1,j,kkm1,q14) + &
                  c18_intu*aloline_on_nn3d(iim1,j,kkm1,q15) + &
                  c04_intu*aloline_on_nn3d(iim1,j,kkm1,q1) + &
                  c05_intu*aloline_on_nn3d(iim1,j,kkm1,q2) + &
                  c06_intu*aloline_on_nn3d(iim1,j,kkm1,q3) + &
                  c22_intu*aloline_on_nn3d(iim1,j,kkm1,q19) + &
                  c23_intu*aloline_on_nn3d(iim1,j,kkm1,q20) + &
                  c24_intu*aloline_on_nn3d(iim1,j,kkm1,q21) + &
                  c08_intu*aloline_on_nn3d(iim1,j,kkm1,q5) + &
                  c26_intu*aloline_on_nn3d(iim1,j,kkm1,q23)))
               aloline_on_nn3d(iip1,jjm1,kkm1,q25) = imask_totreg3d(iip1,jjm1,kkm1) * &
                  (alo_d*c03_slined + abs_sc*(c12_intu*aloline_on_nn3d(iip1,jjm1,kkm1,q10) + &
                  c15_intu*aloline_on_nn3d(iip1,jjm1,kkm1,q13) + &
                  c18_intu*aloline_on_nn3d(iip1,jjm1,kkm1,q16) + &
                  c06_intu*aloline_on_nn3d(iip1,jjm1,kkm1,q4) + &
                  c24_intu*aloline_on_nn3d(iip1,jjm1,kkm1,q22)))
               aloline_on_nn3d(i,jjm1,kkm1,q26) = imask_totreg3d(i,jjm1,kkm1) * &
                  (alo_u*c15_slineu + abs_sc*(c11_intu*aloline_on_nn3d(i,jjm1,kkm1,q10) + &
                  c12_intu*aloline_on_nn3d(i,jjm1,kkm1,q11) + &
                  c14_intu*aloline_on_nn3d(i,jjm1,kkm1,q13) + &
                  c15_intu*aloline_on_nn3d(i,jjm1,kkm1,q14) + &
                  c17_intu*aloline_on_nn3d(i,jjm1,kkm1,q16) + &
                  c18_intu*aloline_on_nn3d(i,jjm1,kkm1,q17) + &
                  c05_intu*aloline_on_nn3d(i,jjm1,kkm1,q4) + &
                  c06_intu*aloline_on_nn3d(i,jjm1,kkm1,q5) + &
                  c23_intu*aloline_on_nn3d(i,jjm1,kkm1,q22) + &
                  c24_intu*aloline_on_nn3d(i,jjm1,kkm1,q23) + &
                  c02_intu*aloline_on_nn3d(i,jjm1,kkm1,q1) + &
                  c08_intu*aloline_on_nn3d(i,jjm1,kkm1,q7) + &
                  c20_intu*aloline_on_nn3d(i,jjm1,kkm1,q19) + &
                  c26_intu*aloline_on_nn3d(i,jjm1,kkm1,q25)))
               aloline_on_nn3d(iim1,jjm1,kkm1,q27) = imask_totreg3d(iim1,jjm1,kkm1) * &
                  (alo_u*c14_slineu + abs_sc*(c10_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q10) + &
                  c11_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q11) + &
                  c12_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q12) + &
                  c13_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q13) + &
                  c14_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q14) + &
                  c15_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q15) + &
                  c16_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q16) + &
                  c17_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q17) + &
                  c18_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q18) + &
                  c04_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q4) + &
                  c05_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q5) + &
                  c06_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q6) + &
                  c22_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q22) + &
                  c23_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q23) + &
                  c24_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q24) + &
                  c02_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q2) + &
                  c08_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q8) + &
                  c20_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q20) + &
                  c26_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q26)))

!
!perform angular integration
               call calc_phinorm(vel_p, vth_p, vth_fiducial, xobs, phinorm)
               wall=wall_global*phinorm
               mintbar3d_tmp(i,j,k) = mintbar3d_tmp(i,j,k) + int3d(i,j,k)*wall
               normalization3d_tmp(i,j,k) = normalization3d_tmp(i,j,k) + wall
!
               aloline_nn3d_tmp(iip1,jjp1,kkp1,q1) = aloline_nn3d_tmp(iip1,jjp1,kkp1,q1) + wall*aloline_on_nn3d(iip1,jjp1,kkp1,q1)
               aloline_nn3d_tmp(i,jjp1,kkp1,q2) = aloline_nn3d_tmp(i,jjp1,kkp1,q2) + wall*aloline_on_nn3d(i,jjp1,kkp1,q2)
               aloline_nn3d_tmp(iim1,jjp1,kkp1,q3) = aloline_nn3d_tmp(iim1,jjp1,kkp1,q3) + wall*aloline_on_nn3d(iim1,jjp1,kkp1,q3)
               aloline_nn3d_tmp(iip1,j,kkp1,q4) = aloline_nn3d_tmp(iip1,j,kkp1,q4) + wall*aloline_on_nn3d(iip1,j,kkp1,q4)
               aloline_nn3d_tmp(i,j,kkp1,q5) = aloline_nn3d_tmp(i,j,kkp1,q5) + wall*aloline_on_nn3d(i,j,kkp1,q5)
               aloline_nn3d_tmp(iim1,j,kkp1,q6) = aloline_nn3d_tmp(iim1,j,kkp1,q6) + wall*aloline_on_nn3d(iim1,j,kkp1,q6)
               aloline_nn3d_tmp(iip1,jjm1,kkp1,q7) = aloline_nn3d_tmp(iip1,jjm1,kkp1,q7) + wall*aloline_on_nn3d(iip1,jjm1,kkp1,q7)
               aloline_nn3d_tmp(i,jjm1,kkp1,q8) = aloline_nn3d_tmp(i,jjm1,kkp1,q8) + wall*aloline_on_nn3d(i,jjm1,kkp1,q8)
               aloline_nn3d_tmp(iim1,jjm1,kkp1,q9) = aloline_nn3d_tmp(iim1,jjm1,kkp1,q9) + wall*aloline_on_nn3d(iim1,jjm1,kkp1,q9)
               aloline_nn3d_tmp(iip1,jjp1,k,q10) = aloline_nn3d_tmp(iip1,jjp1,k,q10) + wall*aloline_on_nn3d(iip1,jjp1,k,q10)
               aloline_nn3d_tmp(i,jjp1,k,q11) = aloline_nn3d_tmp(i,jjp1,k,q11) + wall*aloline_on_nn3d(i,jjp1,k,q11)
               aloline_nn3d_tmp(iim1,jjp1,k,q12) = aloline_nn3d_tmp(iim1,jjp1,k,q12) + wall*aloline_on_nn3d(iim1,jjp1,k,q12)
               aloline_nn3d_tmp(iip1,j,k,q13) = aloline_nn3d_tmp(iip1,j,k,q13) + wall*aloline_on_nn3d(iip1,j,k,q13)
               aloline_nn3d_tmp(i,j,k,q14) = aloline_nn3d_tmp(i,j,k,q14) + wall*aloline_on_nn3d(i,j,k,q14)
               aloline_nn3d_tmp(iim1,j,k,q15) = aloline_nn3d_tmp(iim1,j,k,q15) + wall*aloline_on_nn3d(iim1,j,k,q15)
               aloline_nn3d_tmp(iip1,jjm1,k,q16) = aloline_nn3d_tmp(iip1,jjm1,k,q16) + wall*aloline_on_nn3d(iip1,jjm1,k,q16) !new
               aloline_nn3d_tmp(i,jjm1,k,q17) = aloline_nn3d_tmp(i,jjm1,k,q17) + wall*aloline_on_nn3d(i,jjm1,k,q17)
               aloline_nn3d_tmp(iim1,jjm1,k,q18) = aloline_nn3d_tmp(iim1,jjm1,k,q18) + wall*aloline_on_nn3d(iim1,jjm1,k,q18)
               aloline_nn3d_tmp(iip1,jjp1,kkm1,q19) = aloline_nn3d_tmp(iip1,jjp1,kkm1,q19) + wall*aloline_on_nn3d(iip1,jjp1,kkm1,q19) !new
               aloline_nn3d_tmp(i,jjp1,kkm1,q20) = aloline_nn3d_tmp(i,jjp1,kkm1,q20) + wall*aloline_on_nn3d(i,jjp1,kkm1,q20) !new
               aloline_nn3d_tmp(iim1,jjp1,kkm1,q21) = aloline_nn3d_tmp(iim1,jjp1,kkm1,q21) + wall*aloline_on_nn3d(iim1,jjp1,kkm1,q21) !new
               aloline_nn3d_tmp(iip1,j,kkm1,q22) = aloline_nn3d_tmp(iip1,j,kkm1,q22) + wall*aloline_on_nn3d(iip1,j,kkm1,q22) !new
               aloline_nn3d_tmp(i,j,kkm1,q23) = aloline_nn3d_tmp(i,j,kkm1,q23) + wall*aloline_on_nn3d(i,j,kkm1,q23)
               aloline_nn3d_tmp(iim1,j,kkm1,q24) = aloline_nn3d_tmp(iim1,j,kkm1,q24) + wall*aloline_on_nn3d(iim1,j,kkm1,q24)
               aloline_nn3d_tmp(iip1,jjm1,kkm1,q25) = aloline_nn3d_tmp(iip1,jjm1,kkm1,q25) + wall*aloline_on_nn3d(iip1,jjm1,kkm1,q25) !new
               aloline_nn3d_tmp(i,jjm1,kkm1,q26) = aloline_nn3d_tmp(i,jjm1,kkm1,q26) + wall*aloline_on_nn3d(i,jjm1,kkm1,q26)
               aloline_nn3d_tmp(iim1,jjm1,kkm1,q27) = aloline_nn3d_tmp(iim1,jjm1,kkm1,q27) + wall*aloline_on_nn3d(iim1,jjm1,kkm1,q27)
!
!
!***********special treatment: point lies directly on boundary**********
!
             case(3)
               mueff=nn_x*x(i)+nn_y*y(j)+nn_z*z(k)
               if(mueff.ge.0.d0) then
!no interpolations required
!for alo
!                 int3d(i,j,k) = 0.d0
!                  int3d(i,j,k) = xic1
                  int3d(i,j,k)=calc_icore_gdark(z(k), sqrt(x(i)**2+y(j)**2+z(k)**2))
                  aloline_on_nn3d(i,j,k,14) = 0.d0
!
!perform angular integration (alo not required, since boundary specified)
                  vel_p = nn_x*velx3d(i,j,k) + nn_y*vely3d(i,j,k) + nn_z*velz3d(i,j,k)
                  vth_p = vth3d(i,j,k)
                  call calc_phinorm(vel_p, vth_p, vth_fiducial, xobs, phinorm)
                  wall=wall_global*phinorm
                  mintbar3d_tmp(i,j,k) = mintbar3d_tmp(i,j,k) + int3d(i,j,k)*wall
                  normalization3d_tmp(i,j,k) = normalization3d_tmp(i,j,k) + wall
               else
!same interpolation as in (1)
!
                  iim1=i-alpha
                  jjm1=j-beta
                  kkm1=k-gamma
                  iip1=i+alpha
                  jjp1=j+beta
                  kkp1=k+gamma
!
!calculate distance of upwind point to (previous) xy-plane, xz-plane and yz-plane
                  dels_xyu=(z(k)-z(kkm1))/nn_z
                  dels_xzu=(y(j)-y(jjm1))/nn_y
                  dels_yzu=(x(i)-x(iim1))/nn_x
                  dels_u=min(dels_xyu,dels_xzu,dels_yzu)
!
!calculate distance of downwind point to (next) xy-plane, xz-plane and yz-plane
                  dels_xyd=(z(kkp1)-z(k))/nn_z
                  dels_xzd=(y(jjp1)-y(j))/nn_y
                  dels_yzd=(x(iip1)-x(i))/nn_x
                  dels_d=min(dels_xyd,dels_xzd,dels_yzd)
!
!----------------------------local point--------------------------------
!
                  sline_p=sline3d(i,j,k)
                  opalbar_p=opalbar3d(i,j,k)
                  velx_p=velx3d(i,j,k)
                  vely_p=vely3d(i,j,k)
                  velz_p=velz3d(i,j,k)
                  vel_p = nn_x*velx_p + nn_y*vely_p + nn_z*velz_p
                  vth_p = vth3d(i,j,k)
!
!----------------------------upwind point-------------------------------
!
                  if(dels_xyu.eq.dels_u) then
!intersection with x-y plane on level k-gamma
                     x_u = x(i) - dels_u*nn_x
                     y_u = y(j) - dels_u*nn_y
                     z_u = z(kkm1)
!
                     iim2=i-2*alpha
                     jjm2=j-2*beta
!
                     call coeff2d_lineu(opalbar3d(iim2,jjm2,kkm1), opalbar3d(iim1,jjm2,kkm1), opalbar3d(i,jjm2,kkm1), &
                        opalbar3d(iim2,jjm1,kkm1), opalbar3d(iim1,jjm1,kkm1), opalbar3d(i,jjm1,kkm1), &
                        opalbar3d(iim2,j,kkm1),    opalbar3d(iim1,j,kkm1),    opalbar3d(i,j,kkm1), &
                        velx3d(iim2,jjm2,kkm1), velx3d(iim1,jjm2,kkm1), velx3d(i,jjm2,kkm1), &
                        velx3d(iim2,jjm1,kkm1), velx3d(iim1,jjm1,kkm1), velx3d(i,jjm1,kkm1), &
                        velx3d(iim2,j,kkm1),    velx3d(iim1,j,kkm1),    velx3d(i,j,kkm1), &
                        vely3d(iim2,jjm2,kkm1), vely3d(iim1,jjm2,kkm1), vely3d(i,jjm2,kkm1), &
                        vely3d(iim2,jjm1,kkm1), vely3d(iim1,jjm1,kkm1), vely3d(i,jjm1,kkm1), &
                        vely3d(iim2,j,kkm1),    vely3d(iim1,j,kkm1),    vely3d(i,j,kkm1), &
                        velz3d(iim2,jjm2,kkm1), velz3d(iim1,jjm2,kkm1), velz3d(i,jjm2,kkm1), &
                        velz3d(iim2,jjm1,kkm1), velz3d(iim1,jjm1,kkm1), velz3d(i,jjm1,kkm1), &
                        velz3d(iim2,j,kkm1),    velz3d(iim1,j,kkm1),    velz3d(i,j,kkm1), &
                        vth3d(iim2,jjm2,kkm1), vth3d(iim1,jjm2,kkm1), vth3d(i,jjm2,kkm1), &
                        vth3d(iim2,jjm1,kkm1), vth3d(iim1,jjm1,kkm1), vth3d(i,jjm1,kkm1), &
                        vth3d(iim2,j,kkm1),    vth3d(iim1,j,kkm1),    vth3d(i,j,kkm1), &
                        sline3d(iim2,jjm2,kkm1), sline3d(iim1,jjm2,kkm1), sline3d(i,jjm2,kkm1), &
                        sline3d(iim2,jjm1,kkm1), sline3d(iim1,jjm1,kkm1), sline3d(i,jjm1,kkm1), &
                        sline3d(iim2,j,kkm1),    sline3d(iim1,j,kkm1),    sline3d(i,j,kkm1), &
                        int3d(iim2,jjm2,kkm1), int3d(iim1,jjm2,kkm1), int3d(i,jjm2,kkm1), &
                        int3d(iim2,jjm1,kkm1), int3d(iim1,jjm1,kkm1), int3d(i,jjm1,kkm1), &
                        int3d(iim2,j,kkm1),    int3d(iim1,j,kkm1),    int3d(i,j,kkm1), &
                        x(iim2), x(iim1), x(i), y(jjm2), y(jjm1), y(j), x_u, y_u, &
                        c10_slineu, c11_slineu, c12_slineu, c13_slineu, c14_slineu, &
                        c15_slineu, c16_slineu, c17_slineu, c18_slineu, &
                        c10_intu, c11_intu, c12_intu, c13_intu, c14_intu, &
                        c15_intu, c16_intu, c17_intu, c18_intu, &
                        opalbar_u, velx_u, vely_u, velz_u, vth_u, sline_u, int_u)
!set interpolation coefficients that are not used to zero
                     c23_slineu=0.d0
                     c24_slineu=0.d0
                     c26_slineu=0.d0
                     c27_slineu=0.d0
!
                     c02_intu = 0.d0
                     c04_intu = 0.d0
                     c05_intu = 0.d0
                     c06_intu = 0.d0
                     c08_intu = 0.d0
                     c20_intu = 0.d0
                     c22_intu = 0.d0
                     c23_intu = 0.d0
                     c24_intu = 0.d0
                     c26_intu = 0.d0
!
                  elseif(dels_xzu.eq.dels_u) then
!intersection with x-z plane on level j-beta
                     x_u = x(i) - dels_u*nn_x
                     y_u = y(jjm1)
                     z_u = z(k) - dels_u*nn_z
!
                     iim2=i-2*alpha
                     kkm2=k-2*gamma

                     call coeff2d_lineu(opalbar3d(iim2,jjm1,kkm2), opalbar3d(iim1,jjm1,kkm2), opalbar3d(i,jjm1,kkm2), &
                        opalbar3d(iim2,jjm1,kkm1), opalbar3d(iim1,jjm1,kkm1), opalbar3d(i,jjm1,kkm1), &
                        opalbar3d(iim2,jjm1,k), opalbar3d(iim1,jjm1,k), opalbar3d(i,jjm1,k), &
                        velx3d(iim2,jjm1,kkm2), velx3d(iim1,jjm1,kkm2), velx3d(i,jjm1,kkm2), &
                        velx3d(iim2,jjm1,kkm1), velx3d(iim1,jjm1,kkm1), velx3d(i,jjm1,kkm1), &
                        velx3d(iim2,jjm1,k), velx3d(iim1,jjm1,k), velx3d(i,jjm1,k), &
                        vely3d(iim2,jjm1,kkm2), vely3d(iim1,jjm1,kkm2), vely3d(i,jjm1,kkm2), &
                        vely3d(iim2,jjm1,kkm1), vely3d(iim1,jjm1,kkm1), vely3d(i,jjm1,kkm1), &
                        vely3d(iim2,jjm1,k), vely3d(iim1,jjm1,k), vely3d(i,jjm1,k), &
                        velz3d(iim2,jjm1,kkm2), velz3d(iim1,jjm1,kkm2), velz3d(i,jjm1,kkm2), &
                        velz3d(iim2,jjm1,kkm1), velz3d(iim1,jjm1,kkm1), velz3d(i,jjm1,kkm1), &
                        velz3d(iim2,jjm1,k), velz3d(iim1,jjm1,k), velz3d(i,jjm1,k), &
                        vth3d(iim2,jjm1,kkm2), vth3d(iim1,jjm1,kkm2), vth3d(i,jjm1,kkm2), &
                        vth3d(iim2,jjm1,kkm1), vth3d(iim1,jjm1,kkm1), vth3d(i,jjm1,kkm1), &
                        vth3d(iim2,jjm1,k), vth3d(iim1,jjm1,k), vth3d(i,jjm1,k), &
                        sline3d(iim2,jjm1,kkm2), sline3d(iim1,jjm1,kkm2), sline3d(i,jjm1,kkm2), &
                        sline3d(iim2,jjm1,kkm1), sline3d(iim1,jjm1,kkm1), sline3d(i,jjm1,kkm1), &
                        sline3d(iim2,jjm1,k), sline3d(iim1,jjm1,k), sline3d(i,jjm1,k), &
                        int3d(iim2,jjm1,kkm2), int3d(iim1,jjm1,kkm2), int3d(i,jjm1,kkm2), &
                        int3d(iim2,jjm1,kkm1), int3d(iim1,jjm1,kkm1), int3d(i,jjm1,kkm1), &
                        int3d(iim2,jjm1,k), int3d(iim1,jjm1,k), int3d(i,jjm1,k), &
                        x(iim2), x(iim1), x(i), z(kkm2), z(kkm1), z(k), x_u, z_u, &
                        c04_slineu, c05_slineu, c06_slineu, c13_slineu, c14_slineu, &
                        c15_slineu, c22_slineu, c23_slineu, c24_slineu, &
                        c04_intu, c05_intu, c06_intu, c13_intu, c14_intu, &
                        c15_intu, c22_intu, c23_intu, c24_intu, &
                        opalbar_u, velx_u, vely_u, velz_u, vth_u, sline_u, int_u)
!set interpolation coefficients that are not used to zero
                     c17_slineu=0.d0
                     c18_slineu=0.d0
                     c26_slineu=0.d0
                     c27_slineu=0.d0
                     c02_intu = 0.d0
                     c08_intu = 0.d0
                     c10_intu = 0.d0
                     c11_intu = 0.d0
                     c12_intu = 0.d0
                     c16_intu = 0.d0
                     c17_intu = 0.d0
                     c18_intu = 0.d0
                     c20_intu = 0.d0
                     c26_intu = 0.d0
!
                  elseif(dels_yzu.eq.dels_u) then
!intersection with y-z plane on level i-alpha
                     x_u = x(iim1)
                     y_u = y(j) - dels_u*nn_y
                     z_u = z(k) - dels_u*nn_z
!
                     jjm2=j-2*beta
                     kkm2=k-2*gamma
!
                     call coeff2d_lineu(opalbar3d(iim1,jjm2,kkm2), opalbar3d(iim1,jjm1,kkm2), opalbar3d(iim1,j,kkm2), &
                        opalbar3d(iim1,jjm2,kkm1), opalbar3d(iim1,jjm1,kkm1), opalbar3d(iim1,j,kkm1), &
                        opalbar3d(iim1,jjm2,k), opalbar3d(iim1,jjm1,k), opalbar3d(iim1,j,k), &
                        velx3d(iim1,jjm2,kkm2), velx3d(iim1,jjm1,kkm2), velx3d(iim1,j,kkm2), &
                        velx3d(iim1,jjm2,kkm1), velx3d(iim1,jjm1,kkm1), velx3d(iim1,j,kkm1), &
                        velx3d(iim1,jjm2,k), velx3d(iim1,jjm1,k), velx3d(iim1,j,k), &
                        vely3d(iim1,jjm2,kkm2), vely3d(iim1,jjm1,kkm2), vely3d(iim1,j,kkm2), &
                        vely3d(iim1,jjm2,kkm1), vely3d(iim1,jjm1,kkm1), vely3d(iim1,j,kkm1), &
                        vely3d(iim1,jjm2,k), vely3d(iim1,jjm1,k), vely3d(iim1,j,k), &
                        velz3d(iim1,jjm2,kkm2), velz3d(iim1,jjm1,kkm2), velz3d(iim1,j,kkm2), &
                        velz3d(iim1,jjm2,kkm1), velz3d(iim1,jjm1,kkm1), velz3d(iim1,j,kkm1), &
                        velz3d(iim1,jjm2,k), velz3d(iim1,jjm1,k), velz3d(iim1,j,k), &
                        vth3d(iim1,jjm2,kkm2), vth3d(iim1,jjm1,kkm2), vth3d(iim1,j,kkm2), &
                        vth3d(iim1,jjm2,kkm1), vth3d(iim1,jjm1,kkm1), vth3d(iim1,j,kkm1), &
                        vth3d(iim1,jjm2,k), vth3d(iim1,jjm1,k), vth3d(iim1,j,k), &
                        sline3d(iim1,jjm2,kkm2), sline3d(iim1,jjm1,kkm2), sline3d(iim1,j,kkm2), &
                        sline3d(iim1,jjm2,kkm1), sline3d(iim1,jjm1,kkm1), sline3d(iim1,j,kkm1), &
                        sline3d(iim1,jjm2,k), sline3d(iim1,jjm1,k), sline3d(iim1,j,k), &
                        int3d(iim1,jjm2,kkm2), int3d(iim1,jjm1,kkm2), int3d(iim1,j,kkm2), &
                        int3d(iim1,jjm2,kkm1), int3d(iim1,jjm1,kkm1), int3d(iim1,j,kkm1), &
                        int3d(iim1,jjm2,k), int3d(iim1,jjm1,k), int3d(iim1,j,k), &
                        y(jjm2), y(jjm1), y(j), z(kkm2), z(kkm1), z(k), y_u, z_u, &
                        c02_slineu, c05_slineu, c08_slineu, c11_slineu, c14_slineu, &
                        c17_slineu, c20_slineu, c23_slineu, c26_slineu, &
                        c02_intu, c05_intu, c08_intu, c11_intu, c14_intu, &
                        c17_intu, c20_intu, c23_intu, c26_intu, &
                        opalbar_u, velx_u, vely_u, velz_u, vth_u, sline_u, int_u)
!set interpolation coefficients that are not used to zero
                     c15_slineu=0.d0
                     c18_slineu=0.d0
                     c24_slineu=0.d0
                     c27_slineu=0.d0
                     c04_intu = 0.d0
                     c06_intu = 0.d0
                     c10_intu = 0.d0
                     c12_intu = 0.d0
                     c13_intu = 0.d0
                     c15_intu = 0.d0
                     c16_intu = 0.d0
                     c18_intu = 0.d0
                     c22_intu = 0.d0
                     c24_intu = 0.d0
!
                  else
                     write(*,'(4es20.8,2l4)') dels_u, dels_xzu, dels_xyu, dels_yzu, dels_u.eq.dels_xzu, dels_xyu.eq.dels_u
                     stop 'error in fsc_line3d: invalid dels_u'
                  endif
!
!---------------------------downwind point------------------------------
!
                  if(dels_xyd.eq.dels_d) then
!intersection with x-y plane on level kkp1
                     x_d = x(i) + dels_d*nn_x
                     y_d = y(j) + dels_d*nn_y
                     z_d = z(kkp1)
!
                     call coeff2d_lined(opalbar3d(iim1,jjm1,kkp1), opalbar3d(i,jjm1,kkp1), opalbar3d(iip1,jjm1,kkp1), &
                        opalbar3d(iim1,j,kkp1), opalbar3d(i,j,kkp1), opalbar3d(iip1,j,kkp1), &
                        opalbar3d(iim1,jjp1,kkp1), opalbar3d(i,jjp1,kkp1), opalbar3d(iip1,jjp1,kkp1), &
                        velx3d(iim1,jjm1,kkp1), velx3d(i,jjm1,kkp1), velx3d(iip1,jjm1,kkp1), &
                        velx3d(iim1,j,kkp1), velx3d(i,j,kkp1), velx3d(iip1,j,kkp1), &
                        velx3d(iim1,jjp1,kkp1), velx3d(i,jjp1,kkp1), velx3d(iip1,jjp1,kkp1), &
                        vely3d(iim1,jjm1,kkp1), vely3d(i,jjm1,kkp1), vely3d(iip1,jjm1,kkp1), &
                        vely3d(iim1,j,kkp1), vely3d(i,j,kkp1), vely3d(iip1,j,kkp1), &
                        vely3d(iim1,jjp1,kkp1), vely3d(i,jjp1,kkp1), vely3d(iip1,jjp1,kkp1), &
                        velz3d(iim1,jjm1,kkp1), velz3d(i,jjm1,kkp1), velz3d(iip1,jjm1,kkp1), &
                        velz3d(iim1,j,kkp1), velz3d(i,j,kkp1), velz3d(iip1,j,kkp1), &
                        velz3d(iim1,jjp1,kkp1), velz3d(i,jjp1,kkp1), velz3d(iip1,jjp1,kkp1), &
                        vth3d(iim1,jjm1,kkp1), vth3d(i,jjm1,kkp1), vth3d(iip1,jjm1,kkp1), &
                        vth3d(iim1,j,kkp1), vth3d(i,j,kkp1), vth3d(iip1,j,kkp1), &
                        vth3d(iim1,jjp1,kkp1), vth3d(i,jjp1,kkp1), vth3d(iip1,jjp1,kkp1), &
                        sline3d(iim1,jjm1,kkp1), sline3d(i,jjm1,kkp1), sline3d(iip1,jjm1,kkp1), &
                        sline3d(iim1,j,kkp1), sline3d(i,j,kkp1), sline3d(iip1,j,kkp1), &
                        sline3d(iim1,jjp1,kkp1), sline3d(i,jjp1,kkp1), sline3d(iip1,jjp1,kkp1), &
                        x(iim1), x(i), x(iip1), y(jjm1), y(j), y(jjp1), x_d, y_d, &
                        c19_slined, c20_slined, c21_slined, c22_slined, c23_slined, c24_slined, &
                        c25_slined, c26_slined, c27_slined, &
                        opalbar_d, velx_d, vely_d, velz_d, vth_d, sline_d)
!
!set interpolation coefficients that are not used to zero
                     c03_slined = 0.d0
                     c06_slined = 0.d0
                     c07_slined = 0.d0
                     c08_slined = 0.d0
                     c09_slined = 0.d0
                     c12_slined = 0.d0
                     c15_slined = 0.d0
                     c16_slined = 0.d0
                     c17_slined = 0.d0
                     c18_slined = 0.d0
!
                  elseif(dels_xzd.eq.dels_d) then
!intersection with x-z plane on level jjp1
                     x_d = x(i) + dels_d*nn_x
                     y_d = y(jjp1)
                     z_d = z(k) + dels_d*nn_z
!
                     call coeff2d_lined(opalbar3d(iim1,jjp1,kkm1), opalbar3d(i,jjp1,kkm1), opalbar3d(iip1,jjp1,kkm1), &
                        opalbar3d(iim1,jjp1,k), opalbar3d(i,jjp1,k), opalbar3d(iip1,jjp1,k), &
                        opalbar3d(iim1,jjp1,kkp1), opalbar3d(i,jjp1,kkp1), opalbar3d(iip1,jjp1,kkp1), &
                        velx3d(iim1,jjp1,kkm1), velx3d(i,jjp1,kkm1), velx3d(iip1,jjp1,kkm1), &
                        velx3d(iim1,jjp1,k), velx3d(i,jjp1,k), velx3d(iip1,jjp1,k), &
                        velx3d(iim1,jjp1,kkp1), velx3d(i,jjp1,kkp1), velx3d(iip1,jjp1,kkp1), &
                        vely3d(iim1,jjp1,kkm1), vely3d(i,jjp1,kkm1), vely3d(iip1,jjp1,kkm1), &
                        vely3d(iim1,jjp1,k), vely3d(i,jjp1,k), vely3d(iip1,jjp1,k), &
                        vely3d(iim1,jjp1,kkp1), vely3d(i,jjp1,kkp1), vely3d(iip1,jjp1,kkp1), &
                        velz3d(iim1,jjp1,kkm1), velz3d(i,jjp1,kkm1), velz3d(iip1,jjp1,kkm1), &
                        velz3d(iim1,jjp1,k), velz3d(i,jjp1,k), velz3d(iip1,jjp1,k), &
                        velz3d(iim1,jjp1,kkp1), velz3d(i,jjp1,kkp1), velz3d(iip1,jjp1,kkp1), &
                        vth3d(iim1,jjp1,kkm1), vth3d(i,jjp1,kkm1), vth3d(iip1,jjp1,kkm1), &
                        vth3d(iim1,jjp1,k), vth3d(i,jjp1,k), vth3d(iip1,jjp1,k), &
                        vth3d(iim1,jjp1,kkp1), vth3d(i,jjp1,kkp1), vth3d(iip1,jjp1,kkp1), &
                        sline3d(iim1,jjp1,kkm1), sline3d(i,jjp1,kkm1), sline3d(iip1,jjp1,kkm1), &
                        sline3d(iim1,jjp1,k), sline3d(i,jjp1,k), sline3d(iip1,jjp1,k), &
                        sline3d(iim1,jjp1,kkp1), sline3d(i,jjp1,kkp1), sline3d(iip1,jjp1,kkp1), &
                        x(iim1), x(i), x(iip1), z(kkm1), z(k), z(kkp1), x_d, z_d, &
                        c07_slined, c08_slined, c09_slined, c16_slined, c17_slined, c18_slined, &
                        c25_slined, c26_slined, c27_slined, &
                        opalbar_d, velx_d, vely_d, velz_d, vth_d, sline_d)
!
!set interpolation coefficients that are not used to zero
                     c03_slined = 0.d0
                     c06_slined = 0.d0
                     c12_slined = 0.d0
                     c15_slined = 0.d0
                     c19_slined = 0.d0
                     c20_slined = 0.d0
                     c21_slined = 0.d0
                     c22_slined = 0.d0
                     c23_slined = 0.d0
                     c24_slined = 0.d0
!
                  elseif(dels_yzd.eq.dels_d) then
!intersection with y-z plane on level iip1
                     x_d = x(iip1)
                     y_d = y(j) + dels_d*nn_y
                     z_d = z(k) + dels_d*nn_z
!
                     call coeff2d_lined(opalbar3d(iip1,jjm1,kkm1), opalbar3d(iip1,j,kkm1), opalbar3d(iip1,jjp1,kkm1), &
                        opalbar3d(iip1,jjm1,k), opalbar3d(iip1,j,k), opalbar3d(iip1,jjp1,k), &
                        opalbar3d(iip1,jjm1,kkp1), opalbar3d(iip1,j,kkp1), opalbar3d(iip1,jjp1,kkp1), &
                        velx3d(iip1,jjm1,kkm1), velx3d(iip1,j,kkm1), velx3d(iip1,jjp1,kkm1), &
                        velx3d(iip1,jjm1,k), velx3d(iip1,j,k), velx3d(iip1,jjp1,k), &
                        velx3d(iip1,jjm1,kkp1), velx3d(iip1,j,kkp1), velx3d(iip1,jjp1,kkp1), &
                        vely3d(iip1,jjm1,kkm1), vely3d(iip1,j,kkm1), vely3d(iip1,jjp1,kkm1), &
                        vely3d(iip1,jjm1,k), vely3d(iip1,j,k), vely3d(iip1,jjp1,k), &
                        vely3d(iip1,jjm1,kkp1), vely3d(iip1,j,kkp1), vely3d(iip1,jjp1,kkp1), &
                        velz3d(iip1,jjm1,kkm1), velz3d(iip1,j,kkm1), velz3d(iip1,jjp1,kkm1), &
                        velz3d(iip1,jjm1,k), velz3d(iip1,j,k), velz3d(iip1,jjp1,k), &
                        velz3d(iip1,jjm1,kkp1), velz3d(iip1,j,kkp1), velz3d(iip1,jjp1,kkp1), &
                        vth3d(iip1,jjm1,kkm1), vth3d(iip1,j,kkm1), vth3d(iip1,jjp1,kkm1), &
                        vth3d(iip1,jjm1,k), vth3d(iip1,j,k), vth3d(iip1,jjp1,k), &
                        vth3d(iip1,jjm1,kkp1), vth3d(iip1,j,kkp1), vth3d(iip1,jjp1,kkp1), &
                        sline3d(iip1,jjm1,kkm1), sline3d(iip1,j,kkm1), sline3d(iip1,jjp1,kkm1), &
                        sline3d(iip1,jjm1,k), sline3d(iip1,j,k), sline3d(iip1,jjp1,k), &
                        sline3d(iip1,jjm1,kkp1), sline3d(iip1,j,kkp1), sline3d(iip1,jjp1,kkp1), &
                        y(jjm1), y(j), y(jjp1), z(kkm1), z(k), z(kkp1), y_d, z_d, &
                        c03_slined, c06_slined, c09_slined, c12_slined, c15_slined, c18_slined, &
                        c21_slined, c24_slined, c27_slined, &
                        opalbar_d, velx_d, vely_d, velz_d, vth_d, sline_d)
!
!set interpolation coefficients that are not used to zero
                     c07_slined = 0.d0
                     c08_slined = 0.d0
                     c16_slined = 0.d0
                     c17_slined = 0.d0
                     c19_slined = 0.d0
                     c20_slined = 0.d0
                     c22_slined = 0.d0
                     c23_slined = 0.d0
                     c25_slined = 0.d0
                     c26_slined = 0.d0
!
                  else
                     write(*,'(4es20.8,l4)') dels_d, dels_xzd, dels_xyd, dels_yzd, dels_d.eq.dels_xzd
                     stop 'error in fsc_line3d: invalid dels_d'
                  endif
!
!--------------------------------radiative transfer---------------------
!
!***debug start: analytical laws
!                  write(*,'(10es20.8)') velx_u, vely_u, velz_u, velx_d, vely_d, velz_d, opalbar_u, opalbar_d
!                  call model_debug(x_u, y_u, z_u, velx_u, vely_u, velz_u, opalbar_u, fdum)
!                  call model_debug(x_d, y_d, z_d, velx_d, vely_d, velz_d, opalbar_d, fdum)
!                  write(*,'(10es20.8)') velx_u, vely_u, velz_u, velx_d, vely_d, velz_d, opalbar_u, opalbar_d
!                  write(*,*)
!***debug end
                  vel_u = velx_u*nn_x + vely_u*nn_y + velz_u*nn_z
                  vel_d = velx_d*nn_x + vely_d*nn_y + velz_d*nn_z
                  call fsc_line(xobs, vth_fiducial, deltax, xcmf_min, xcmf_max, int_u, opalbar_u, opalbar_p, opalbar_d, &
                     sline_u, sline_p, sline_d, vel_u, vel_p, vel_d, vth_u, vth_p, vth_d, &
                     dels_u, dels_d, abs_sc, contr_sc, int_sc, alo_u, alo_p, alo_d)
!                 call fsc_line_lin(xobs, vth_fiducial, deltax, xcmf_min, xcmf_max, int_u, opalbar_u, opalbar_p, sline_u, sline_p, &
!                                     vel_u, vel_p, vth_u, vth_p, dels_u, abs_sc, contr_sc, int_sc, alo_u, alo_p)
!                 alo_d=0.d0
                  int3d(i,j,k) = abs_sc*int_u + alo_u*sline_u + alo_p*sline_p + alo_d*sline_d
!
                  aloline_on_nn3d(iip1,jjp1,kkp1,q1) = imask_totreg3d(iip1,jjp1,kkp1)*&
                     alo_d*c27_slined
                  aloline_on_nn3d(i,jjp1,kkp1,q2) = imask_totreg3d(i,jjp1,kkp1) * &
                     (alo_d*c26_slined + abs_sc*c26_intu*aloline_on_nn3d(i,jjp1,kkp1,q1))
                  aloline_on_nn3d(iim1,jjp1,kkp1,q3) = imask_totreg3d(iim1,jjp1,kkp1) * &
                     (alo_d*c25_slined + abs_sc*c26_intu*aloline_on_nn3d(iim1,jjp1,kkp1,q2))
                  aloline_on_nn3d(iip1,j,kkp1,q4) = imask_totreg3d(iip1,j,kkp1) * &
                     (alo_d*c24_slined + abs_sc*c24_intu*aloline_on_nn3d(iip1,j,kkp1,q1))
                  aloline_on_nn3d(i,j,kkp1,q5) = imask_totreg3d(i,j,kkp1) * &
                     (alo_d*c23_slined + abs_sc*(c23_intu*aloline_on_nn3d(i,j,kkp1,q1) + &
                     c24_intu*aloline_on_nn3d(i,j,kkp1,q2) + &
                     c26_intu*aloline_on_nn3d(i,j,kkp1,q4)))
                  aloline_on_nn3d(iim1,j,kkp1,q6) = imask_totreg3d(iim1,j,kkp1) * &
                     (alo_d*c22_slined + abs_sc*(c22_intu*aloline_on_nn3d(iim1,j,kkp1,q1) + &
                     c23_intu*aloline_on_nn3d(iim1,j,kkp1,q2) + &
                     c24_intu*aloline_on_nn3d(iim1,j,kkp1,q3) + &
                     c26_intu*aloline_on_nn3d(iim1,j,kkp1,q5)))
                  aloline_on_nn3d(iip1,jjm1,kkp1,q7) = imask_totreg3d(iip1,jjm1,kkp1) * &
                     (alo_d*c21_slined + abs_sc*c24_intu*aloline_on_nn3d(iip1,jjm1,kkp1,q4))
                  aloline_on_nn3d(i,jjm1,kkp1,q8) = imask_totreg3d(i,jjm1,kkp1) * &
                     (alo_d*c20_slined + abs_sc*(c23_intu*aloline_on_nn3d(i,jjm1,kkp1,q4) + &
                     c24_intu*aloline_on_nn3d(i,jjm1,kkp1,q5) + &
                     c20_intu*aloline_on_nn3d(i,jjm1,kkp1,q1) + &
                     c26_intu*aloline_on_nn3d(i,jjm1,kkp1,q7)))
                  aloline_on_nn3d(iim1,jjm1,kkp1,q9) = imask_totreg3d(iim1,jjm1,kkp1) * &
                     (alo_d*c19_slined + abs_sc*(c22_intu*aloline_on_nn3d(iim1,jjm1,kkp1,q4) + &
                     c23_intu*aloline_on_nn3d(iim1,jjm1,kkp1,q5) + &
                     c24_intu*aloline_on_nn3d(iim1,jjm1,kkp1,q6) + &
                     c20_intu*aloline_on_nn3d(iim1,jjm1,kkp1,q2) + &
                     c26_intu*aloline_on_nn3d(iim1,jjm1,kkp1,q8)))
                  aloline_on_nn3d(iip1,jjp1,k,q10) = imask_totreg3d(iip1,jjp1,k) * &
                     (alo_d*c18_slined + abs_sc*c18_intu*aloline_on_nn3d(iip1,jjp1,k,q1))
                  aloline_on_nn3d(i,jjp1,k,q11) = imask_totreg3d(i,jjp1,k) * &
                     (alo_d*c17_slined + abs_sc*(c17_intu*aloline_on_nn3d(i,jjp1,k,q1) + &
                     c18_intu*aloline_on_nn3d(i,jjp1,k,q2) + &
                     c26_intu*aloline_on_nn3d(i,jjp1,k,q10)))
                  aloline_on_nn3d(iim1,jjp1,k,q12) = imask_totreg3d(iim1,jjp1,k) * &
                     (alo_d*c16_slined + abs_sc*(c16_intu*aloline_on_nn3d(iim1,jjp1,k,q1) + &
                     c17_intu*aloline_on_nn3d(iim1,jjp1,k,q2) + &
                     c18_intu*aloline_on_nn3d(iim1,jjp1,k,q3) + &
                     c26_intu*aloline_on_nn3d(iim1,jjp1,k,q11)))
                  aloline_on_nn3d(iip1,j,k,q13) = imask_totreg3d(iip1,j,k) * &
                     (alo_d*c15_slined + abs_sc*(c15_intu*aloline_on_nn3d(iip1,j,k,q1) + &
                     c18_intu*aloline_on_nn3d(iip1,j,k,q4) + &
                     c24_intu*aloline_on_nn3d(iip1,j,k,q10)))
                  aloline_on_nn3d(i,j,k,q14) = imask_totreg3d(i,j,k) * &
                     (alo_p + abs_sc*(c14_intu*aloline_on_nn3d(i,j,k,q1) + &
                     c15_intu*aloline_on_nn3d(i,j,k,q2) + &
                     c17_intu*aloline_on_nn3d(i,j,k,q4) + &
                     c26_intu*aloline_on_nn3d(i,j,k,q13) + &
                     c18_intu*aloline_on_nn3d(i,j,k,q5) + &
                     c23_intu*aloline_on_nn3d(i,j,k,q10) + &
                     c24_intu*aloline_on_nn3d(i,j,k,q11)))
                  aloline_on_nn3d(iim1,j,k,q15) = imask_totreg3d(iim1,j,k) * &
                     (alo_u*c26_slineu + abs_sc*(c13_intu*aloline_on_nn3d(iim1,j,k,q1) + &
                     c14_intu*aloline_on_nn3d(iim1,j,k,q2) + &
                     c15_intu*aloline_on_nn3d(iim1,j,k,q3) + &
                     c16_intu*aloline_on_nn3d(iim1,j,k,q4) + &
                     c17_intu*aloline_on_nn3d(iim1,j,k,q5) + &
                     c18_intu*aloline_on_nn3d(iim1,j,k,q6) + &
                     c22_intu*aloline_on_nn3d(iim1,j,k,q10) + &
                     c23_intu*aloline_on_nn3d(iim1,j,k,q11) + &
                     c24_intu*aloline_on_nn3d(iim1,j,k,q12) + &
                     c26_intu*aloline_on_nn3d(iim1,j,k,q14)))
                  aloline_on_nn3d(iip1,jjm1,k,q16) = imask_totreg3d(iip1,jjm1,k) * &
                     (alo_d*c12_slined + abs_sc*(c12_intu*aloline_on_nn3d(iip1,jjm1,k,q1) + &
                     c15_intu*aloline_on_nn3d(iip1,jjm1,k,q4) + &
                     c18_intu*aloline_on_nn3d(iip1,jjm1,k,q7) + &
                     c24_intu*aloline_on_nn3d(iip1,jjm1,k,q13)))
                  aloline_on_nn3d(i,jjm1,k,q17) = imask_totreg3d(i,jjm1,k) * &
                     (alo_u*c24_slineu + abs_sc*(c11_intu*aloline_on_nn3d(i,jjm1,k,q1) + &
                     c12_intu*aloline_on_nn3d(i,jjm1,k,q2) + &
                     c14_intu*aloline_on_nn3d(i,jjm1,k,q4) + &
                     c15_intu*aloline_on_nn3d(i,jjm1,k,q5) + &
                     c17_intu*aloline_on_nn3d(i,jjm1,k,q7) + &
                     c18_intu*aloline_on_nn3d(i,jjm1,k,q8) + &
                     c23_intu*aloline_on_nn3d(i,jjm1,k,q13) + &
                     c24_intu*aloline_on_nn3d(i,jjm1,k,q14) + &
                     c20_intu*aloline_on_nn3d(i,jjm1,k,q10) + &
                     c26_intu*aloline_on_nn3d(i,jjm1,k,q16)))
                  aloline_on_nn3d(iim1,jjm1,k,q18) = imask_totreg3d(iim1,jjm1,k) * &
                     (alo_u*c23_slineu + abs_sc*(c10_intu*aloline_on_nn3d(iim1,jjm1,k,q1) + &
                     c11_intu*aloline_on_nn3d(iim1,jjm1,k,q2) + &
                     c12_intu*aloline_on_nn3d(iim1,jjm1,k,q3) + &
                     c13_intu*aloline_on_nn3d(iim1,jjm1,k,q4) + &
                     c14_intu*aloline_on_nn3d(iim1,jjm1,k,q5) + &
                     c15_intu*aloline_on_nn3d(iim1,jjm1,k,q6) + &
                     c16_intu*aloline_on_nn3d(iim1,jjm1,k,q7) + &
                     c17_intu*aloline_on_nn3d(iim1,jjm1,k,q8) + &
                     c18_intu*aloline_on_nn3d(iim1,jjm1,k,q9) + &
                     c22_intu*aloline_on_nn3d(iim1,jjm1,k,q13) + &
                     c23_intu*aloline_on_nn3d(iim1,jjm1,k,q14) + &
                     c24_intu*aloline_on_nn3d(iim1,jjm1,k,q15) + &
                     c20_intu*aloline_on_nn3d(iim1,jjm1,k,q11) + &
                     c26_intu*aloline_on_nn3d(iim1,jjm1,k,q17)))
                  aloline_on_nn3d(iip1,jjp1,kkm1,q19) = imask_totreg3d(iip1,jjp1,kkm1) * &
                     (alo_d*c09_slined + abs_sc*c18_intu*aloline_on_nn3d(iip1,jjp1,kkm1,q10))
                  aloline_on_nn3d(i,jjp1,kkm1,q20) = imask_totreg3d(i,jjp1,kkm1) * &
                     (alo_d*c08_slined + abs_sc*(c17_intu*aloline_on_nn3d(i,jjp1,kkm1,q10) + &
                     c18_intu*aloline_on_nn3d(i,jjp1,kkm1,q11) + &
                     c08_intu*aloline_on_nn3d(i,jjp1,kkm1,q1) + &
                     c26_intu*aloline_on_nn3d(i,jjp1,kkm1,q19)))
                  aloline_on_nn3d(iim1,jjp1,kkm1,q21) = imask_totreg3d(iim1,jjp1,kkm1) * &
                     (alo_d*c07_slined + abs_sc*(c16_intu*aloline_on_nn3d(iim1,jjp1,kkm1,q10) + &
                     c17_intu*aloline_on_nn3d(iim1,jjp1,kkm1,q11) + &
                     c18_intu*aloline_on_nn3d(iim1,jjp1,kkm1,q12) + &
                     c08_intu*aloline_on_nn3d(iim1,jjp1,kkm1,q2) + &
                     c26_intu*aloline_on_nn3d(iim1,jjp1,kkm1,q20)))
                  aloline_on_nn3d(iip1,j,kkm1,q22) = imask_totreg3d(iip1,j,kkm1) * &
                     (alo_d*c06_slined + abs_sc*(c15_intu*aloline_on_nn3d(iip1,j,kkm1,q10) + &
                     c18_intu*aloline_on_nn3d(iip1,j,kkm1,q13) + &
                     c06_intu*aloline_on_nn3d(iip1,j,kkm1,q1) + &
                     c24_intu*aloline_on_nn3d(iip1,j,kkm1,q19)))
                  aloline_on_nn3d(i,j,kkm1,q23) = imask_totreg3d(i,j,kkm1) * &
                     (alo_u*c18_slineu + abs_sc*(c14_intu*aloline_on_nn3d(i,j,kkm1,q10) + &
                     c15_intu*aloline_on_nn3d(i,j,kkm1,q11) + &
                     c17_intu*aloline_on_nn3d(i,j,kkm1,q13) + &
                     c18_intu*aloline_on_nn3d(i,j,kkm1,q14) + &
                     c05_intu*aloline_on_nn3d(i,j,kkm1,q1) + &
                     c06_intu*aloline_on_nn3d(i,j,kkm1,q2) + &
                     c23_intu*aloline_on_nn3d(i,j,kkm1,q19) + &
                     c24_intu*aloline_on_nn3d(i,j,kkm1,q20) + &
                     c08_intu*aloline_on_nn3d(i,j,kkm1,q4) + &
                     c26_intu*aloline_on_nn3d(i,j,kkm1,q22)))
                  aloline_on_nn3d(iim1,j,kkm1,q24) = imask_totreg3d(iim1,j,kkm1) * &
                     (alo_u*c17_slineu + abs_sc*(c13_intu*aloline_on_nn3d(iim1,j,kkm1,q10) + &
                     c14_intu*aloline_on_nn3d(iim1,j,kkm1,q11) + &
                     c15_intu*aloline_on_nn3d(iim1,j,kkm1,q12) + &
                     c16_intu*aloline_on_nn3d(iim1,j,kkm1,q13) + &
                     c17_intu*aloline_on_nn3d(iim1,j,kkm1,q14) + &
                     c18_intu*aloline_on_nn3d(iim1,j,kkm1,q15) + &
                     c04_intu*aloline_on_nn3d(iim1,j,kkm1,q1) + &
                     c05_intu*aloline_on_nn3d(iim1,j,kkm1,q2) + &
                     c06_intu*aloline_on_nn3d(iim1,j,kkm1,q3) + &
                     c22_intu*aloline_on_nn3d(iim1,j,kkm1,q19) + &
                     c23_intu*aloline_on_nn3d(iim1,j,kkm1,q20) + &
                     c24_intu*aloline_on_nn3d(iim1,j,kkm1,q21) + &
                     c08_intu*aloline_on_nn3d(iim1,j,kkm1,q5) + &
                     c26_intu*aloline_on_nn3d(iim1,j,kkm1,q23)))
                  aloline_on_nn3d(iip1,jjm1,kkm1,q25) = imask_totreg3d(iip1,jjm1,kkm1) * &
                     (alo_d*c03_slined + abs_sc*(c12_intu*aloline_on_nn3d(iip1,jjm1,kkm1,q10) + &
                     c15_intu*aloline_on_nn3d(iip1,jjm1,kkm1,q13) + &
                     c18_intu*aloline_on_nn3d(iip1,jjm1,kkm1,q16) + &
                     c06_intu*aloline_on_nn3d(iip1,jjm1,kkm1,q4) + &
                     c24_intu*aloline_on_nn3d(iip1,jjm1,kkm1,q22)))
                  aloline_on_nn3d(i,jjm1,kkm1,q26) = imask_totreg3d(i,jjm1,kkm1) * &
                     (alo_u*c15_slineu + abs_sc*(c11_intu*aloline_on_nn3d(i,jjm1,kkm1,q10) + &
                     c12_intu*aloline_on_nn3d(i,jjm1,kkm1,q11) + &
                     c14_intu*aloline_on_nn3d(i,jjm1,kkm1,q13) + &
                     c15_intu*aloline_on_nn3d(i,jjm1,kkm1,q14) + &
                     c17_intu*aloline_on_nn3d(i,jjm1,kkm1,q16) + &
                     c18_intu*aloline_on_nn3d(i,jjm1,kkm1,q17) + &
                     c05_intu*aloline_on_nn3d(i,jjm1,kkm1,q4) + &
                     c06_intu*aloline_on_nn3d(i,jjm1,kkm1,q5) + &
                     c23_intu*aloline_on_nn3d(i,jjm1,kkm1,q22) + &
                     c24_intu*aloline_on_nn3d(i,jjm1,kkm1,q23) + &
                     c02_intu*aloline_on_nn3d(i,jjm1,kkm1,q1) + &
                     c08_intu*aloline_on_nn3d(i,jjm1,kkm1,q7) + &
                     c20_intu*aloline_on_nn3d(i,jjm1,kkm1,q19) + &
                     c26_intu*aloline_on_nn3d(i,jjm1,kkm1,q25)))
                  aloline_on_nn3d(iim1,jjm1,kkm1,q27) = imask_totreg3d(iim1,jjm1,kkm1) * &
                     (alo_u*c14_slineu + abs_sc*(c10_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q10) + &
                     c11_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q11) + &
                     c12_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q12) + &
                     c13_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q13) + &
                     c14_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q14) + &
                     c15_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q15) + &
                     c16_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q16) + &
                     c17_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q17) + &
                     c18_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q18) + &
                     c04_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q4) + &
                     c05_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q5) + &
                     c06_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q6) + &
                     c22_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q22) + &
                     c23_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q23) + &
                     c24_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q24) + &
                     c02_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q2) + &
                     c08_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q8) + &
                     c20_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q20) + &
                     c26_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q26)))
!
!perform angular integration
                  call calc_phinorm(vel_p, vth_p, vth_fiducial, xobs, phinorm)
                  wall=wall_global*phinorm
                  mintbar3d_tmp(i,j,k) = mintbar3d_tmp(i,j,k) + int3d(i,j,k)*wall
                  normalization3d_tmp(i,j,k) = normalization3d_tmp(i,j,k) + wall
!
                  aloline_nn3d_tmp(iip1,jjp1,kkp1,q1) = aloline_nn3d_tmp(iip1,jjp1,kkp1,q1) + wall*aloline_on_nn3d(iip1,jjp1,kkp1,q1)
                  aloline_nn3d_tmp(i,jjp1,kkp1,q2) = aloline_nn3d_tmp(i,jjp1,kkp1,q2) + wall*aloline_on_nn3d(i,jjp1,kkp1,q2)
                  aloline_nn3d_tmp(iim1,jjp1,kkp1,q3) = aloline_nn3d_tmp(iim1,jjp1,kkp1,q3) + wall*aloline_on_nn3d(iim1,jjp1,kkp1,q3)
                  aloline_nn3d_tmp(iip1,j,kkp1,q4) = aloline_nn3d_tmp(iip1,j,kkp1,q4) + wall*aloline_on_nn3d(iip1,j,kkp1,q4)
                  aloline_nn3d_tmp(i,j,kkp1,q5) = aloline_nn3d_tmp(i,j,kkp1,q5) + wall*aloline_on_nn3d(i,j,kkp1,q5)
                  aloline_nn3d_tmp(iim1,j,kkp1,q6) = aloline_nn3d_tmp(iim1,j,kkp1,q6) + wall*aloline_on_nn3d(iim1,j,kkp1,q6)
                  aloline_nn3d_tmp(iip1,jjm1,kkp1,q7) = aloline_nn3d_tmp(iip1,jjm1,kkp1,q7) + wall*aloline_on_nn3d(iip1,jjm1,kkp1,q7)
                  aloline_nn3d_tmp(i,jjm1,kkp1,q8) = aloline_nn3d_tmp(i,jjm1,kkp1,q8) + wall*aloline_on_nn3d(i,jjm1,kkp1,q8)
                  aloline_nn3d_tmp(iim1,jjm1,kkp1,q9) = aloline_nn3d_tmp(iim1,jjm1,kkp1,q9) + wall*aloline_on_nn3d(iim1,jjm1,kkp1,q9)
                  aloline_nn3d_tmp(iip1,jjp1,k,q10) = aloline_nn3d_tmp(iip1,jjp1,k,q10) + wall*aloline_on_nn3d(iip1,jjp1,k,q10)
                  aloline_nn3d_tmp(i,jjp1,k,q11) = aloline_nn3d_tmp(i,jjp1,k,q11) + wall*aloline_on_nn3d(i,jjp1,k,q11)
                  aloline_nn3d_tmp(iim1,jjp1,k,q12) = aloline_nn3d_tmp(iim1,jjp1,k,q12) + wall*aloline_on_nn3d(iim1,jjp1,k,q12)
                  aloline_nn3d_tmp(iip1,j,k,q13) = aloline_nn3d_tmp(iip1,j,k,q13) + wall*aloline_on_nn3d(iip1,j,k,q13)
                  aloline_nn3d_tmp(i,j,k,q14) = aloline_nn3d_tmp(i,j,k,q14) + wall*aloline_on_nn3d(i,j,k,q14)
                  aloline_nn3d_tmp(iim1,j,k,q15) = aloline_nn3d_tmp(iim1,j,k,q15) + wall*aloline_on_nn3d(iim1,j,k,q15)
                  aloline_nn3d_tmp(iip1,jjm1,k,q16) = aloline_nn3d_tmp(iip1,jjm1,k,q16) + wall*aloline_on_nn3d(iip1,jjm1,k,q16) !new
                  aloline_nn3d_tmp(i,jjm1,k,q17) = aloline_nn3d_tmp(i,jjm1,k,q17) + wall*aloline_on_nn3d(i,jjm1,k,q17)
                  aloline_nn3d_tmp(iim1,jjm1,k,q18) = aloline_nn3d_tmp(iim1,jjm1,k,q18) + wall*aloline_on_nn3d(iim1,jjm1,k,q18)
                  aloline_nn3d_tmp(iip1,jjp1,kkm1,q19) = aloline_nn3d_tmp(iip1,jjp1,kkm1,q19) + wall*aloline_on_nn3d(iip1,jjp1,kkm1,q19) !new
                  aloline_nn3d_tmp(i,jjp1,kkm1,q20) = aloline_nn3d_tmp(i,jjp1,kkm1,q20) + wall*aloline_on_nn3d(i,jjp1,kkm1,q20) !new
                  aloline_nn3d_tmp(iim1,jjp1,kkm1,q21) = aloline_nn3d_tmp(iim1,jjp1,kkm1,q21) + wall*aloline_on_nn3d(iim1,jjp1,kkm1,q21) !new
                  aloline_nn3d_tmp(iip1,j,kkm1,q22) = aloline_nn3d_tmp(iip1,j,kkm1,q22) + wall*aloline_on_nn3d(iip1,j,kkm1,q22) !new
                  aloline_nn3d_tmp(i,j,kkm1,q23) = aloline_nn3d_tmp(i,j,kkm1,q23) + wall*aloline_on_nn3d(i,j,kkm1,q23)
                  aloline_nn3d_tmp(iim1,j,kkm1,q24) = aloline_nn3d_tmp(iim1,j,kkm1,q24) + wall*aloline_on_nn3d(iim1,j,kkm1,q24)
                  aloline_nn3d_tmp(iip1,jjm1,kkm1,q25) = aloline_nn3d_tmp(iip1,jjm1,kkm1,q25) + wall*aloline_on_nn3d(iip1,jjm1,kkm1,q25) !new
                  aloline_nn3d_tmp(i,jjm1,kkm1,q26) = aloline_nn3d_tmp(i,jjm1,kkm1,q26) + wall*aloline_on_nn3d(i,jjm1,kkm1,q26)
                  aloline_nn3d_tmp(iim1,jjm1,kkm1,q27) = aloline_nn3d_tmp(iim1,jjm1,kkm1,q27) + wall*aloline_on_nn3d(iim1,jjm1,kkm1,q27)
!
               endif
!
             case default
!
            end select
!
         enddo
      enddo
   enddo
!
!***debug start
!close(1)
!***debug end
!
!
end subroutine fsc_line3d
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine fsc_line3d_lin(oindx,xobsindx)
!
!-----------------------------------------------------------------------
!--------short characteristics for line radiative transfer in 3d--------
!-----------calculating intensties for given mu,phi specified-----------
!-----------by input oindx, and for xobs specified by xobsindx----------
!------------------only linear interpolations are used------------------
!-----------------------------------------------------------------------
!
   use prog_type
   use fund_const
   use dime3d, only: ndxmax, ndzmax, ndymax, x, y, z
   use dime3d, only: int3d, opalbar3d, sline3d, velx3d, vely3d, velz3d, vth3d, imask3d, imask_totreg3d, &
      aloline_on_nn3d, aloline_nn3d_tmp, mintbar3d_tmp, normalization3d_tmp
   use angles, only: n_x, n_y, n_z, weight_omega, q_alo
   use freq, only: nodes_xobs, weight_xobs, nxobs, deltax, xcmf_min, xcmf_max
   use params_input, only: vth_fiducial
   use bcondition, only: xic1
   use params_stellar, only: smajorax_a, smajorax_b, smajorax_c
!use mod_debug
!
   implicit none
!
! ... arguments
   integer(i4b), intent(in) :: oindx, xobsindx
!
! ... local scalars
   integer(i4b) :: i, j, k
   integer(i4b) :: alpha, beta, gamma
   integer(i4b) :: startx, starty, startz, endx, endy, endz
   integer(i4b) :: iim1, ii, jjm1, jj, kkm1, kk
   real(dp) :: nn_x, nn_y, nn_z, mueff, xobs, wall, wall_global, phinorm, fdum
   real(dp) :: dels_xyu, dels_xzu, dels_yzu, dels_u, dels_r
   real(dp) :: opalbar_p, velx_p, vely_p, velz_p, vel_p, vth_p, sline_p
   real(dp) :: x_u, y_u, z_u, int_u, opalbar_u, velx_u, vely_u, velz_u, vel_u, vth_u, sline_u
   real(dp) :: abs_sc, int_sc, contr_sc
   real(dp) :: alo_u, alo_p
   integer :: q14, q15, q17, q18, q23, q24, q26, q27
   real(dp) :: c14_slineu, c15_slineu, c17_slineu, c18_slineu, &
      c23_slineu, c24_slineu, c26_slineu, c27_slineu, &
      c14_intu, c15_intu, c17_intu, c18_intu, &
      c23_intu, c24_intu, c26_intu
!
! ... local arrays
!
! ... local functions
   real(dp) :: dist_sphere, dist_ellipsoid, calc_icore_gdark
!
! ... local characters
!character(len=50) :: enter
!
! ... local logicals
!
!
!directions
   nn_x=n_x(oindx)
   nn_y=n_y(oindx)
   nn_z=n_z(oindx)
!
!frequency
   xobs=nodes_xobs(xobsindx)
!
!angular and frequency  integration weights
   wall_global=weight_omega(oindx)*weight_xobs(xobsindx)
!
!indices for nearest neighbour alo
   q14=q_alo(oindx,14)
   q15=q_alo(oindx,15)
   q17=q_alo(oindx,17)
   q18=q_alo(oindx,18)
   q23=q_alo(oindx,23)
   q24=q_alo(oindx,24)
   q26=q_alo(oindx,26)
   q27=q_alo(oindx,27)
!
!-----------------------------------------------------------------------
!
!set directional index-parameter (phantom points are excluded from calculation)
!
!index parameter:
!         if n_x >= 0                 if n_x < 0
!                startx = 2                  startx = ndxmax-1
!                endx = ndxmax-1             endx = 2
!                alpha=  1                   alpha=-1
!
!         if n_y >= 0                 if n_y < 0
!                starty = 2                  starty = ndymax-1
!                endy = ndymax-1             endy = 2
!                beta =  1                   beta =-1
!
!         if n_z >= 0                 if n_z < 0
!                startz = 2                  startz = ndzmax-1
!                endz = ndzmax-1             endz = 2
!                gamma = 1                   gamma = -1
!
   if(nn_x.gt.0.d0) then
      startx = 3
      endx = ndxmax-2
      alpha=  1
   elseif(nn_x.lt.0.d0) then
      startx = ndxmax-2
      endx = 3
      alpha=-1
   else
      stop 'error in fsc_line3d: n_x = 0 not allowed'
   endif
!
   if(nn_y.gt.0.d0) then
      starty = 3
      endy = ndymax-2
      beta =  1
   elseif(nn_y.lt.0.d0) then
      starty = ndymax-2
      endy = 3
      beta =-1
   else
      stop 'error in fsc_line3d: n_y = 0 not allowed'
   endif
!
   if(nn_z.gt.0.d0) then
      startz = 3
      endz = ndzmax-2
      gamma=  1
   elseif(nn_z.lt.0.d0) then
      startz = ndzmax-2
      endz = 3
      gamma=-1
   else
      stop 'error in fsc_line3d: n_z = 0 not allowed'
   endif
!
!--------------------reset the intensities and alo----------------------
!
   call set_boundary3d(nn_x, nn_y, nn_z)
!
   aloline_on_nn3d=0.d0
!
!-----------------------------------------------------------------------
!
   do k=startz, endz, gamma
      do j=starty, endy, beta
         do i=startx, endx, alpha
!
            select case(imask3d(i,j,k))
!
!************************standard radiative transfer********************
!
             case(1)
!
               iim1=i-alpha
               jjm1=j-beta
               kkm1=k-gamma
!
!calculate distance of upwind point to (previous) xy-plane, xz-plane and yz-plane
               dels_xyu=(z(k)-z(kkm1))/nn_z
               dels_xzu=(y(j)-y(jjm1))/nn_y
               dels_yzu=(x(i)-x(iim1))/nn_x
               dels_u=min(dels_xyu,dels_xzu,dels_yzu)
!
!----------------------------local point--------------------------------
!
               sline_p=sline3d(i,j,k)
               opalbar_p=opalbar3d(i,j,k)
               velx_p=velx3d(i,j,k)
               vely_p=vely3d(i,j,k)
               velz_p=velz3d(i,j,k)
               vel_p = nn_x*velx_p + nn_y*vely_p + nn_z*velz_p
               vth_p = vth3d(i,j,k)
!
!----------------------------upwind point-------------------------------
!
               if(dels_xyu.eq.dels_u) then
!intersection with x-y plane on level k-gamma
                  x_u = x(i) - dels_u*nn_x
                  y_u = y(j) - dels_u*nn_y
                  z_u = z(kkm1)
!
                  call coeff2d_lineu_lin(opalbar3d(iim1,jjm1,kkm1), opalbar3d(i,jjm1,kkm1), opalbar3d(iim1,j,kkm1), opalbar3d(i,j,kkm1), &
                     velx3d(iim1,jjm1,kkm1), velx3d(i,jjm1,kkm1), velx3d(iim1,j,kkm1), velx3d(i,j,kkm1), &
                     vely3d(iim1,jjm1,kkm1), vely3d(i,jjm1,kkm1), vely3d(iim1,j,kkm1), vely3d(i,j,kkm1), &
                     velz3d(iim1,jjm1,kkm1), velz3d(i,jjm1,kkm1), velz3d(iim1,j,kkm1), velz3d(i,j,kkm1), &
                     vth3d(iim1,jjm1,kkm1), vth3d(i,jjm1,kkm1), vth3d(iim1,j,kkm1), vth3d(i,j,kkm1), &
                     sline3d(iim1,jjm1,kkm1), sline3d(i,jjm1,kkm1), sline3d(iim1,j,kkm1), sline3d(i,j,kkm1), &
                     int3d(iim1,jjm1,kkm1), int3d(i,jjm1,kkm1), int3d(iim1,j,kkm1), int3d(i,j,kkm1), &
                     x(iim1), x(i), y(jjm1), y(j), x_u, y_u, &
                     c14_slineu, c15_slineu, c17_slineu, c18_slineu, &
                     c14_intu, c15_intu, c17_intu, c18_intu, opalbar_u, velx_u, vely_u, velz_u, vth_u, sline_u, int_u)
                  c23_slineu=0.d0
                  c24_slineu=0.d0
                  c26_slineu=0.d0
                  c23_intu=0.d0
                  c24_intu=0.d0
                  c26_intu=0.d0
!
               elseif(dels_xzu.eq.dels_u) then
!intersection with x-z plane on level j-beta
                  x_u = x(i) - dels_u*nn_x
                  y_u = y(jjm1)
                  z_u = z(k) - dels_u*nn_z
!
                  call coeff2d_lineu_lin(opalbar3d(iim1,jjm1,kkm1), opalbar3d(i,jjm1,kkm1), opalbar3d(iim1,jjm1,k), opalbar3d(i,jjm1,k), &
                     velx3d(iim1,jjm1,kkm1), velx3d(i,jjm1,kkm1), velx3d(iim1,jjm1,k), velx3d(i,jjm1,k), &
                     vely3d(iim1,jjm1,kkm1), vely3d(i,jjm1,kkm1), vely3d(iim1,jjm1,k), vely3d(i,jjm1,k), &
                     velz3d(iim1,jjm1,kkm1), velz3d(i,jjm1,kkm1), velz3d(iim1,jjm1,k), velz3d(i,jjm1,k), &
                     vth3d(iim1,jjm1,kkm1), vth3d(i,jjm1,kkm1), vth3d(iim1,jjm1,k), vth3d(i,jjm1,k), &
                     sline3d(iim1,jjm1,kkm1), sline3d(i,jjm1,kkm1), sline3d(iim1,jjm1,k), sline3d(i,jjm1,k), &
                     int3d(iim1,jjm1,kkm1), int3d(i,jjm1,kkm1), int3d(iim1,jjm1,k), int3d(i,jjm1,k), &
                     x(iim1), x(i), z(kkm1), z(k), x_u, z_u, &
                     c14_slineu, c15_slineu, c23_slineu, c24_slineu, &
                     c14_intu, c15_intu, c23_intu, c24_intu, opalbar_u, velx_u, vely_u, velz_u, vth_u, sline_u, int_u)
                  c17_slineu=0.d0
                  c18_slineu=0.d0
                  c26_slineu=0.d0
                  c17_intu=0.d0
                  c18_intu=0.d0
                  c26_intu=0.d0
!
               elseif(dels_yzu.eq.dels_u) then
!intersection with y-z plane on level i-alpha
                  x_u = x(iim1)
                  y_u = y(j) - dels_u*nn_y
                  z_u = z(k) - dels_u*nn_z
!
                  call coeff2d_lineu_lin(opalbar3d(iim1,jjm1,kkm1), opalbar3d(iim1,j,kkm1), opalbar3d(iim1,jjm1,k), opalbar3d(iim1,j,k), &
                     velx3d(iim1,jjm1,kkm1), velx3d(iim1,j,kkm1), velx3d(iim1,jjm1,k), velx3d(iim1,j,k), &
                     vely3d(iim1,jjm1,kkm1), vely3d(iim1,j,kkm1), vely3d(iim1,jjm1,k), vely3d(iim1,j,k), &
                     velz3d(iim1,jjm1,kkm1), velz3d(iim1,j,kkm1), velz3d(iim1,jjm1,k), velz3d(iim1,j,k), &
                     vth3d(iim1,jjm1,kkm1), vth3d(iim1,j,kkm1), vth3d(iim1,jjm1,k), vth3d(iim1,j,k), &
                     sline3d(iim1,jjm1,kkm1), sline3d(iim1,j,kkm1), sline3d(iim1,jjm1,k), sline3d(iim1,j,k), &
                     int3d(iim1,jjm1,kkm1), int3d(iim1,j,kkm1), int3d(iim1,jjm1,k), int3d(iim1,j,k), &
                     y(jjm1), y(j), z(kkm1), z(k), y_u, z_u, &
                     c14_slineu, c17_slineu, c23_slineu, c26_slineu, &
                     c14_intu, c17_intu, c23_intu, c26_intu, opalbar_u, velx_u, vely_u, velz_u, vth_u, sline_u, int_u)
                  c15_slineu=0.d0
                  c18_slineu=0.d0
                  c24_slineu=0.d0
                  c15_intu=0.d0
                  c18_intu=0.d0
                  c24_intu=0.d0
!
               else
                  write(*,'(4es20.8,2l4)') dels_u, dels_xzu, dels_xyu, dels_yzu, dels_u.eq.dels_xzu, dels_xyu.eq.dels_u
                  stop 'error in fsc_line3d: invalid dels_u'
               endif
!
!--------------------------------radiative transfer---------------------
!
!***debug start: analytic law
!               write(*,'(4es20.8)') velx_u, vely_u, velz_u, opalbar_u
!               call model_debug(x_u, y_u, z_u, velx_u, vely_u, velz_u, opalbar_u, fdum)
!               write(*,'(4es20.8)') velx_u, vely_u, velz_u, opalbar_u
!               write(*,*)
!                xu_debug=x_u
!                yu_debug=y_u
!                zu_debug=z_u
!                nnx_debug=nn_x
!                nny_debug=nn_y
!                nnz_debug=nn_z
!***debug end
               vel_u = velx_u*nn_x + vely_u*nn_y + velz_u*nn_z
               call fsc_line_lin(xobs, vth_fiducial, deltax, xcmf_min, xcmf_max, int_u, opalbar_u, opalbar_p, sline_u, sline_p, &
                  vel_u, vel_p, vth_u, vth_p, dels_u, abs_sc, contr_sc, int_sc, alo_u, alo_p)
               int3d(i,j,k) = abs_sc*int_u + alo_u*sline_u + alo_p*sline_p
!***debug start
!               if(int3d(i,j,k).lt.0.) then
!                  write(*,*) sline_u, sline_p
!                  write(*,*) opalbar_u, opalbar_p
!                  write(*,*) vel_u, vel_p
!                  write(*,*) vth_u, vth_p
!                  write(*,*) dels_u
!                  write(*,*) alo_u, alo_p, abs_sc
!                  write(*,*) int_u, int3d(i,j,k)
!                  stop 'error'
!               endif
!***debug end
!
               aloline_on_nn3d(i,j,k,q14) = imask_totreg3d(i,j,k) * alo_p
               aloline_on_nn3d(iim1,j,k,q15) = imask_totreg3d(iim1,j,k) * &
                  (alo_u*c26_slineu + abs_sc*c26_intu*aloline_on_nn3d(iim1,j,k,q14))

               aloline_on_nn3d(i,jjm1,k,q17) = imask_totreg3d(i,jjm1,k) * &
                  (alo_u*c24_slineu + abs_sc*c24_intu*aloline_on_nn3d(i,jjm1,k,q14))
               aloline_on_nn3d(iim1,jjm1,k,q18) = imask_totreg3d(iim1,jjm1,k) * &
                  (alo_u*c23_slineu + abs_sc*(c23_intu*aloline_on_nn3d(iim1,jjm1,k,q14) + &
                  c26_intu*aloline_on_nn3d(iim1,jjm1,k,q17)))
               aloline_on_nn3d(i,j,kkm1,q23) = imask_totreg3d(i,j,kkm1) * &
                  (alo_u*c18_slineu + abs_sc*c18_intu*aloline_on_nn3d(i,j,kkm1,q14))
               aloline_on_nn3d(iim1,j,kkm1,q24) = imask_totreg3d(iim1,j,kkm1) * &
                  (alo_u*c17_slineu + abs_sc*(c17_intu*aloline_on_nn3d(iim1,j,kkm1,q14) + &
                  c18_intu*aloline_on_nn3d(iim1,j,kkm1,q15) + &
                  c26_intu*aloline_on_nn3d(iim1,j,kkm1,q23)))
               aloline_on_nn3d(i,jjm1,kkm1,q26) = imask_totreg3d(i,jjm1,kkm1) * &
                  (alo_u*c15_slineu + abs_sc*(c15_intu*aloline_on_nn3d(i,jjm1,kkm1,q14) + &
                  c18_intu*aloline_on_nn3d(i,jjm1,kkm1,q17) + &
                  c24_intu*aloline_on_nn3d(i,jjm1,kkm1,q23)))
               aloline_on_nn3d(iim1,jjm1,kkm1,q27) = imask_totreg3d(iim1,jjm1,kkm1) * &
                  (alo_u*c14_slineu + abs_sc*(c14_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q14) + &
                  c15_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q15) + &
                  c17_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q17) + &
                  c18_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q18) + &
                  c23_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q23) + &
                  c24_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q24) + &
                  c26_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q26)))
!
!perform angular integration
               call calc_phinorm(vel_p, vth_p, vth_fiducial, xobs, phinorm)
               wall=wall_global*phinorm
               mintbar3d_tmp(i,j,k) = mintbar3d_tmp(i,j,k) + int3d(i,j,k)*wall
               normalization3d_tmp(i,j,k) = normalization3d_tmp(i,j,k) + wall
!
               aloline_nn3d_tmp(i,j,k,q14) = aloline_nn3d_tmp(i,j,k,q14) + wall*aloline_on_nn3d(i,j,k,q14)
               aloline_nn3d_tmp(iim1,j,k,q15) = aloline_nn3d_tmp(iim1,j,k,q15) + wall*aloline_on_nn3d(iim1,j,k,q15)
               aloline_nn3d_tmp(i,jjm1,k,q17) = aloline_nn3d_tmp(i,jjm1,k,q17) + wall*aloline_on_nn3d(i,jjm1,k,q17)
               aloline_nn3d_tmp(iim1,jjm1,k,q18) = aloline_nn3d_tmp(iim1,jjm1,k,q18) + wall*aloline_on_nn3d(iim1,jjm1,k,q18)
               aloline_nn3d_tmp(i,j,kkm1,q23) = aloline_nn3d_tmp(i,j,kkm1,q23) + wall*aloline_on_nn3d(i,j,kkm1,q23)
               aloline_nn3d_tmp(iim1,j,kkm1,q24) = aloline_nn3d_tmp(iim1,j,kkm1,q24) + wall*aloline_on_nn3d(iim1,j,kkm1,q24)
               aloline_nn3d_tmp(i,jjm1,kkm1,q26) = aloline_nn3d_tmp(i,jjm1,kkm1,q26) + wall*aloline_on_nn3d(i,jjm1,kkm1,q26)
               aloline_nn3d_tmp(iim1,jjm1,kkm1,q27) = aloline_nn3d_tmp(iim1,jjm1,kkm1,q27) + wall*aloline_on_nn3d(iim1,jjm1,kkm1,q27)
!
!***********special treatment: point is in vicinity of boundary*********
!
             case(2)
!
               iim1=i-alpha
               jjm1=j-beta
               kkm1=k-gamma
!
!calculate distance of upwind point to (previous) xy-plane, xz-plane and yz-plane
               dels_xyu=(z(k)-z(kkm1))/nn_z
               dels_xzu=(y(j)-y(jjm1))/nn_y
               dels_yzu=(x(i)-x(iim1))/nn_x
!
!calculate distance to the boundary (here: sphere)
!               dels_r=dist_sphere(nn_x,nn_y,nn_z,x(i),y(j),z(k))
               dels_r=dist_ellipsoid(nn_x,nn_y,nn_z,x(i),y(j),z(k),smajorax_a,smajorax_b,smajorax_c)
!
               dels_u=min(dels_xyu,dels_xzu,dels_yzu,dels_r)
!
!----------------------------local point--------------------------------
!
               sline_p=sline3d(i,j,k)
               opalbar_p=opalbar3d(i,j,k)
               velx_p=velx3d(i,j,k)
               vely_p=vely3d(i,j,k)
               velz_p=velz3d(i,j,k)
               vel_p = nn_x*velx_p + nn_y*vely_p + nn_z*velz_p
               vth_p = vth3d(i,j,k)
!
!----------------------------upwind point-------------------------------
!
               if(dels_r.eq.dels_u) then
                  x_u = x(i) - dels_u*nn_x
                  y_u = y(j) - dels_u*nn_y
                  z_u = z(k) - dels_u*nn_z
!
!for alo
!                  int_u=0.d0
!                  int_u=xic1
                  int_u=calc_icore_gdark(z_u, sqrt(x_u**2+y_u**2+z_u**2))
!
                  call coeff3d_lineu_lin(iim1, i, jjm1, j, kkm1, k, x_u, y_u, z_u, &
                     c14_slineu, c15_slineu, c17_slineu, c18_slineu, &
                     c23_slineu, c24_slineu, c26_slineu, c27_slineu, &
                     opalbar_u, velx_u, vely_u, velz_u, vth_u, sline_u)
!
               elseif(dels_xyu.eq.dels_u) then
!intersection with x-y plane on level k-gamma
                  x_u = x(i) - dels_u*nn_x
                  y_u = y(j) - dels_u*nn_y
                  z_u = z(kkm1)
!
                  call coeff2d_lineu_lin(opalbar3d(iim1,jjm1,kkm1), opalbar3d(i,jjm1,kkm1), opalbar3d(iim1,j,kkm1), opalbar3d(i,j,kkm1), &
                     velx3d(iim1,jjm1,kkm1), velx3d(i,jjm1,kkm1), velx3d(iim1,j,kkm1), velx3d(i,j,kkm1), &
                     vely3d(iim1,jjm1,kkm1), vely3d(i,jjm1,kkm1), vely3d(iim1,j,kkm1), vely3d(i,j,kkm1), &
                     velz3d(iim1,jjm1,kkm1), velz3d(i,jjm1,kkm1), velz3d(iim1,j,kkm1), velz3d(i,j,kkm1), &
                     vth3d(iim1,jjm1,kkm1), vth3d(i,jjm1,kkm1), vth3d(iim1,j,kkm1), vth3d(i,j,kkm1), &
                     sline3d(iim1,jjm1,kkm1), sline3d(i,jjm1,kkm1), sline3d(iim1,j,kkm1), sline3d(i,j,kkm1), &
                     int3d(iim1,jjm1,kkm1), int3d(i,jjm1,kkm1), int3d(iim1,j,kkm1), int3d(i,j,kkm1), &
                     x(iim1), x(i), y(jjm1), y(j), x_u, y_u, &
                     c14_slineu, c15_slineu, c17_slineu, c18_slineu, &
                     c14_intu, c15_intu, c17_intu, c18_intu, opalbar_u, velx_u, vely_u, velz_u, vth_u, sline_u, int_u)
                  c23_slineu=0.d0
                  c24_slineu=0.d0
                  c26_slineu=0.d0
                  c27_slineu=0.d0
                  c23_intu=0.d0
                  c24_intu=0.d0
                  c26_intu=0.d0
!
               elseif(dels_xzu.eq.dels_u) then
!intersection with x-z plane on level j-beta
                  x_u = x(i) - dels_u*nn_x
                  y_u = y(jjm1)
                  z_u = z(k) - dels_u*nn_z
!
                  call coeff2d_lineu_lin(opalbar3d(iim1,jjm1,kkm1), opalbar3d(i,jjm1,kkm1), opalbar3d(iim1,jjm1,k), opalbar3d(i,jjm1,k), &
                     velx3d(iim1,jjm1,kkm1), velx3d(i,jjm1,kkm1), velx3d(iim1,jjm1,k), velx3d(i,jjm1,k), &
                     vely3d(iim1,jjm1,kkm1), vely3d(i,jjm1,kkm1), vely3d(iim1,jjm1,k), vely3d(i,jjm1,k), &
                     velz3d(iim1,jjm1,kkm1), velz3d(i,jjm1,kkm1), velz3d(iim1,jjm1,k), velz3d(i,jjm1,k), &
                     vth3d(iim1,jjm1,kkm1), vth3d(i,jjm1,kkm1), vth3d(iim1,jjm1,k), vth3d(i,jjm1,k), &
                     sline3d(iim1,jjm1,kkm1), sline3d(i,jjm1,kkm1), sline3d(iim1,jjm1,k), sline3d(i,jjm1,k), &
                     int3d(iim1,jjm1,kkm1), int3d(i,jjm1,kkm1), int3d(iim1,jjm1,k), int3d(i,jjm1,k), &
                     x(iim1), x(i), z(kkm1), z(k), x_u, z_u, &
                     c14_slineu, c15_slineu, c23_slineu, c24_slineu, &
                     c14_intu, c15_intu, c23_intu, c24_intu, opalbar_u, velx_u, vely_u, velz_u, vth_u, sline_u, int_u)
                  c17_slineu=0.d0
                  c18_slineu=0.d0
                  c26_slineu=0.d0
                  c27_slineu=0.d0
                  c17_intu=0.d0
                  c18_intu=0.d0
                  c26_intu=0.d0
!
               elseif(dels_yzu.eq.dels_u) then
!intersection with y-z plane on level i-alpha
                  x_u = x(iim1)
                  y_u = y(j) - dels_u*nn_y
                  z_u = z(k) - dels_u*nn_z
!
                  call coeff2d_lineu_lin(opalbar3d(iim1,jjm1,kkm1), opalbar3d(iim1,j,kkm1), opalbar3d(iim1,jjm1,k), opalbar3d(iim1,j,k), &
                     velx3d(iim1,jjm1,kkm1), velx3d(iim1,j,kkm1), velx3d(iim1,jjm1,k), velx3d(iim1,j,k), &
                     vely3d(iim1,jjm1,kkm1), vely3d(iim1,j,kkm1), vely3d(iim1,jjm1,k), vely3d(iim1,j,k), &
                     velz3d(iim1,jjm1,kkm1), velz3d(iim1,j,kkm1), velz3d(iim1,jjm1,k), velz3d(iim1,j,k), &
                     vth3d(iim1,jjm1,kkm1), vth3d(iim1,j,kkm1), vth3d(iim1,jjm1,k), vth3d(iim1,j,k), &
                     sline3d(iim1,jjm1,kkm1), sline3d(iim1,j,kkm1), sline3d(iim1,jjm1,k), sline3d(iim1,j,k), &
                     int3d(iim1,jjm1,kkm1), int3d(iim1,j,kkm1), int3d(iim1,jjm1,k), int3d(iim1,j,k), &
                     y(jjm1), y(j), z(kkm1), z(k), y_u, z_u, &
                     c14_slineu, c17_slineu, c23_slineu, c26_slineu, &
                     c14_intu, c17_intu, c23_intu, c26_intu, opalbar_u, velx_u, vely_u, velz_u, vth_u, sline_u, int_u)
                  c15_slineu=0.d0
                  c18_slineu=0.d0
                  c24_slineu=0.d0
                  c27_slineu=0.d0
                  c15_intu=0.d0
                  c18_intu=0.d0
                  c24_intu=0.d0
!
               else
                  write(*,'(4es20.8,2l4)') dels_u, dels_xzu, dels_xyu, dels_yzu, dels_u.eq.dels_xzu, dels_xyu.eq.dels_u
                  stop 'error in fsc_line3d: invalid dels_u'
               endif
!
!--------------------------------radiative transfer---------------------
!
!***debug start: analytic laws
!               call model_debug(x_u, y_u, z_u, velx_u, vely_u, velz_u, opalbar_u, fdum)
!                xu_debug=x_u
!                yu_debug=y_u
!                zu_debug=z_u
!                nnx_debug=nn_x
!                nny_debug=nn_y
!                nnz_debug=nn_z
!***debug end
               vel_u = velx_u*nn_x + vely_u*nn_y + velz_u*nn_z
               call fsc_line_lin(xobs, vth_fiducial, deltax, xcmf_min, xcmf_max, int_u, opalbar_u, opalbar_p, sline_u, sline_p, &
                  vel_u, vel_p, vth_u, vth_p, dels_u, abs_sc, contr_sc, int_sc, alo_u, alo_p)
               int3d(i,j,k) = abs_sc*int_u + alo_u*sline_u + alo_p*sline_p
!***debug start
!               if(int3d(i,j,k).lt.0.) then
!                  write(*,*) sline_u, sline_p
!                  write(*,*) opalbar_u, opalbar_p
!                  write(*,*) vel_u, vel_p
!                  write(*,*) vth_u, vth_p
!                  write(*,*) dels_u
!                  write(*,*) alo_u, alo_p, abs_sc
!                  write(*,*) int_u, int3d(i,j,k)
!                  stop 'error'
!               endif
!***debug end
!
               aloline_on_nn3d(i,j,k,q14) = imask_totreg3d(i,j,k) * (alo_p + alo_u*c27_slineu)
               aloline_on_nn3d(iim1,j,k,q15) = imask_totreg3d(iim1,j,k) * &
                  (alo_u*c26_slineu + abs_sc*c26_intu*aloline_on_nn3d(iim1,j,k,q14))

               aloline_on_nn3d(i,jjm1,k,q17) = imask_totreg3d(i,jjm1,k) * &
                  (alo_u*c24_slineu + abs_sc*c24_intu*aloline_on_nn3d(i,jjm1,k,q14))
               aloline_on_nn3d(iim1,jjm1,k,q18) = imask_totreg3d(iim1,jjm1,k) * &
                  (alo_u*c23_slineu + abs_sc*(c23_intu*aloline_on_nn3d(iim1,jjm1,k,q14) + &
                  c26_intu*aloline_on_nn3d(iim1,jjm1,k,q17)))
               aloline_on_nn3d(i,j,kkm1,q23) = imask_totreg3d(i,j,kkm1) * &
                  (alo_u*c18_slineu + abs_sc*c18_intu*aloline_on_nn3d(i,j,kkm1,q14))
               aloline_on_nn3d(iim1,j,kkm1,q24) = imask_totreg3d(iim1,j,kkm1) * &
                  (alo_u*c17_slineu + abs_sc*(c17_intu*aloline_on_nn3d(iim1,j,kkm1,q14) + &
                  c18_intu*aloline_on_nn3d(iim1,j,kkm1,q15) + &
                  c26_intu*aloline_on_nn3d(iim1,j,kkm1,q23)))
               aloline_on_nn3d(i,jjm1,kkm1,q26) = imask_totreg3d(i,jjm1,kkm1) * &
                  (alo_u*c15_slineu + abs_sc*(c15_intu*aloline_on_nn3d(i,jjm1,kkm1,q14) + &
                  c18_intu*aloline_on_nn3d(i,jjm1,kkm1,q17) + &
                  c24_intu*aloline_on_nn3d(i,jjm1,kkm1,q23)))
               aloline_on_nn3d(iim1,jjm1,kkm1,q27) = imask_totreg3d(iim1,jjm1,kkm1) * &
                  (alo_u*c14_slineu + abs_sc*(c14_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q14) + &
                  c15_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q15) + &
                  c17_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q17) + &
                  c18_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q18) + &
                  c23_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q23) + &
                  c24_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q24) + &
                  c26_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q26)))
!
!perform integration
               call calc_phinorm(vel_p, vth_p, vth_fiducial, xobs, phinorm)
               wall=wall_global*phinorm
               mintbar3d_tmp(i,j,k) = mintbar3d_tmp(i,j,k) + int3d(i,j,k)*wall
               normalization3d_tmp(i,j,k) = normalization3d_tmp(i,j,k) + wall
!
               aloline_nn3d_tmp(i,j,k,q14) = aloline_nn3d_tmp(i,j,k,q14) + wall*aloline_on_nn3d(i,j,k,q14)
               aloline_nn3d_tmp(iim1,j,k,q15) = aloline_nn3d_tmp(iim1,j,k,q15) + wall*aloline_on_nn3d(iim1,j,k,q15)
               aloline_nn3d_tmp(i,jjm1,k,q17) = aloline_nn3d_tmp(i,jjm1,k,q17) + wall*aloline_on_nn3d(i,jjm1,k,q17)
               aloline_nn3d_tmp(iim1,jjm1,k,q18) = aloline_nn3d_tmp(iim1,jjm1,k,q18) + wall*aloline_on_nn3d(iim1,jjm1,k,q18)
               aloline_nn3d_tmp(i,j,kkm1,q23) = aloline_nn3d_tmp(i,j,kkm1,q23) + wall*aloline_on_nn3d(i,j,kkm1,q23)
               aloline_nn3d_tmp(iim1,j,kkm1,q24) = aloline_nn3d_tmp(iim1,j,kkm1,q24) + wall*aloline_on_nn3d(iim1,j,kkm1,q24)
               aloline_nn3d_tmp(i,jjm1,kkm1,q26) = aloline_nn3d_tmp(i,jjm1,kkm1,q26) + wall*aloline_on_nn3d(i,jjm1,kkm1,q26)
               aloline_nn3d_tmp(iim1,jjm1,kkm1,q27) = aloline_nn3d_tmp(iim1,jjm1,kkm1,q27) + wall*aloline_on_nn3d(iim1,jjm1,kkm1,q27)
!
!
!***********special treatment: point lies directly on boundary**********
!
             case(3)
               mueff=nn_x*x(i)+nn_y*y(j)+nn_z*z(k)
               if(mueff.ge.0.d0) then
!no interpolations required
!for alo
!                 int3d(i,j,k) = 0.d0
!                  int3d(i,j,k) = xic1
                  int3d(i,j,k)=calc_icore_gdark(z(k), sqrt(x(i)**2+y(j)**2+z(k)**2))
                  aloline_on_nn3d(i,j,k,14) = 0.d0
!
!perform angular integration (alo not required, since boundary specified)
                  vel_p = nn_x*velx3d(i,j,k) + nn_y*vely3d(i,j,k) + nn_z*velz3d(i,j,k)
                  vth_p = vth3d(i,j,k)
                  call calc_phinorm(vel_p, vth_p, vth_fiducial, xobs, phinorm)
                  wall=wall_global*phinorm
                  mintbar3d_tmp(i,j,k) = mintbar3d_tmp(i,j,k) + int3d(i,j,k)*wall
                  normalization3d_tmp(i,j,k) = normalization3d_tmp(i,j,k) + wall
               else
!same interpolation as in (1)
!
                  iim1=i-alpha
                  jjm1=j-beta
                  kkm1=k-gamma
!
!calculate distance of upwind point to (previous) xy-plane, xz-plane and yz-plane
                  dels_xyu=(z(k)-z(kkm1))/nn_z
                  dels_xzu=(y(j)-y(jjm1))/nn_y
                  dels_yzu=(x(i)-x(iim1))/nn_x
                  dels_u=min(dels_xyu,dels_xzu,dels_yzu)
!
!----------------------------local point--------------------------------
!
                  sline_p=sline3d(i,j,k)
                  opalbar_p=opalbar3d(i,j,k)
                  velx_p=velx3d(i,j,k)
                  vely_p=vely3d(i,j,k)
                  velz_p=velz3d(i,j,k)
                  vel_p = nn_x*velx_p + nn_y*vely_p + nn_z*velz_p
                  vth_p = vth3d(i,j,k)
!
!----------------------------upwind point-------------------------------
!
                  if(dels_xyu.eq.dels_u) then
!intersection with x-y plane on level k-gamma
                     x_u = x(i) - dels_u*nn_x
                     y_u = y(j) - dels_u*nn_y
                     z_u = z(kkm1)
!
                     call coeff2d_lineu_lin(opalbar3d(iim1,jjm1,kkm1), opalbar3d(i,jjm1,kkm1), opalbar3d(iim1,j,kkm1), opalbar3d(i,j,kkm1), &
                        velx3d(iim1,jjm1,kkm1), velx3d(i,jjm1,kkm1), velx3d(iim1,j,kkm1), velx3d(i,j,kkm1), &
                        vely3d(iim1,jjm1,kkm1), vely3d(i,jjm1,kkm1), vely3d(iim1,j,kkm1), vely3d(i,j,kkm1), &
                        velz3d(iim1,jjm1,kkm1), velz3d(i,jjm1,kkm1), velz3d(iim1,j,kkm1), velz3d(i,j,kkm1), &
                        vth3d(iim1,jjm1,kkm1), vth3d(i,jjm1,kkm1), vth3d(iim1,j,kkm1), vth3d(i,j,kkm1), &
                        sline3d(iim1,jjm1,kkm1), sline3d(i,jjm1,kkm1), sline3d(iim1,j,kkm1), sline3d(i,j,kkm1), &
                        int3d(iim1,jjm1,kkm1), int3d(i,jjm1,kkm1), int3d(iim1,j,kkm1), int3d(i,j,kkm1), &
                        x(iim1), x(i), y(jjm1), y(j), x_u, y_u, &
                        c14_slineu, c15_slineu, c17_slineu, c18_slineu, &
                        c14_intu, c15_intu, c17_intu, c18_intu, opalbar_u, velx_u, vely_u, velz_u, vth_u, sline_u, int_u)
                     c23_slineu=0.d0
                     c24_slineu=0.d0
                     c26_slineu=0.d0
                     c23_intu=0.d0
                     c24_intu=0.d0
                     c26_intu=0.d0
!
                  elseif(dels_xzu.eq.dels_u) then
!intersection with x-z plane on level j-beta
                     x_u = x(i) - dels_u*nn_x
                     y_u = y(jjm1)
                     z_u = z(k) - dels_u*nn_z
!
                     call coeff2d_lineu_lin(opalbar3d(iim1,jjm1,kkm1), opalbar3d(i,jjm1,kkm1), opalbar3d(iim1,jjm1,k), opalbar3d(i,jjm1,k), &
                        velx3d(iim1,jjm1,kkm1), velx3d(i,jjm1,kkm1), velx3d(iim1,jjm1,k), velx3d(i,jjm1,k), &
                        vely3d(iim1,jjm1,kkm1), vely3d(i,jjm1,kkm1), vely3d(iim1,jjm1,k), vely3d(i,jjm1,k), &
                        velz3d(iim1,jjm1,kkm1), velz3d(i,jjm1,kkm1), velz3d(iim1,jjm1,k), velz3d(i,jjm1,k), &
                        vth3d(iim1,jjm1,kkm1), vth3d(i,jjm1,kkm1), vth3d(iim1,jjm1,k), vth3d(i,jjm1,k), &
                        sline3d(iim1,jjm1,kkm1), sline3d(i,jjm1,kkm1), sline3d(iim1,jjm1,k), sline3d(i,jjm1,k), &
                        int3d(iim1,jjm1,kkm1), int3d(i,jjm1,kkm1), int3d(iim1,jjm1,k), int3d(i,jjm1,k), &
                        x(iim1), x(i), z(kkm1), z(k), x_u, z_u, &
                        c14_slineu, c15_slineu, c23_slineu, c24_slineu, &
                        c14_intu, c15_intu, c23_intu, c24_intu, opalbar_u, velx_u, vely_u, velz_u, vth_u, sline_u, int_u)
                     c17_slineu=0.d0
                     c18_slineu=0.d0
                     c26_slineu=0.d0
                     c17_intu=0.d0
                     c18_intu=0.d0
                     c26_intu=0.d0
!
                  elseif(dels_yzu.eq.dels_u) then
!intersection with y-z plane on level i-alpha
                     x_u = x(iim1)
                     y_u = y(j) - dels_u*nn_y
                     z_u = z(k) - dels_u*nn_z
!
                     call coeff2d_lineu_lin(opalbar3d(iim1,jjm1,kkm1), opalbar3d(iim1,j,kkm1), opalbar3d(iim1,jjm1,k), opalbar3d(iim1,j,k), &
                        velx3d(iim1,jjm1,kkm1), velx3d(iim1,j,kkm1), velx3d(iim1,jjm1,k), velx3d(iim1,j,k), &
                        vely3d(iim1,jjm1,kkm1), vely3d(iim1,j,kkm1), vely3d(iim1,jjm1,k), vely3d(iim1,j,k), &
                        velz3d(iim1,jjm1,kkm1), velz3d(iim1,j,kkm1), velz3d(iim1,jjm1,k), velz3d(iim1,j,k), &
                        vth3d(iim1,jjm1,kkm1), vth3d(iim1,j,kkm1), vth3d(iim1,jjm1,k), vth3d(iim1,j,k), &
                        sline3d(iim1,jjm1,kkm1), sline3d(iim1,j,kkm1), sline3d(iim1,jjm1,k), sline3d(iim1,j,k), &
                        int3d(iim1,jjm1,kkm1), int3d(iim1,j,kkm1), int3d(iim1,jjm1,k), int3d(iim1,j,k), &
                        y(jjm1), y(j), z(kkm1), z(k), y_u, z_u, &
                        c14_slineu, c17_slineu, c23_slineu, c26_slineu, &
                        c14_intu, c17_intu, c23_intu, c26_intu, opalbar_u, velx_u, vely_u, velz_u, vth_u, sline_u, int_u)
                     c15_slineu=0.d0
                     c18_slineu=0.d0
                     c24_slineu=0.d0
                     c15_intu=0.d0
                     c18_intu=0.d0
                     c24_intu=0.d0
!
                  else
                     write(*,'(4es20.8,2l4)') dels_u, dels_xzu, dels_xyu, dels_yzu, dels_u.eq.dels_xzu, dels_xyu.eq.dels_u
                     stop 'error in fsc_line3d: invalid dels_u'
                  endif
!
!--------------------------------radiative transfer---------------------
!
!***debug start: analytic laws
!               call model_debug(x_u, y_u, z_u, velx_u, vely_u, velz_u, opalbar_u, fdum)
!                xu_debug=x_u
!                yu_debug=y_u
!                zu_debug=z_u
!                nnx_debug=nn_x
!                nny_debug=nn_y
!                nnz_debug=nn_z
!***debug end
                  vel_u = velx_u*nn_x + vely_u*nn_y + velz_u*nn_z
                  call fsc_line_lin(xobs, vth_fiducial, deltax, xcmf_min, xcmf_max, int_u, opalbar_u, opalbar_p, sline_u, sline_p, &
                     vel_u, vel_p, vth_u, vth_p, dels_u, abs_sc, contr_sc, int_sc, alo_u, alo_p)
                  int3d(i,j,k) = abs_sc*int_u + alo_u*sline_u + alo_p*sline_p
!***debug start
!               if(int3d(i,j,k).lt.0.) then
!                  write(*,*) sline_u, sline_p
!                  write(*,*) opalbar_u, opalbar_p
!                  write(*,*) vel_u, vel_p
!                  write(*,*) vth_u, vth_p
!                  write(*,*) dels_u
!                  write(*,*) alo_u, alo_p, abs_sc
!                  write(*,*) int_u, int3d(i,j,k)
!                  stop 'error'
!               endif
!***debug end
!
                  aloline_on_nn3d(i,j,k,q14) = imask_totreg3d(i,j,k) * alo_p
                  aloline_on_nn3d(iim1,j,k,q15) = imask_totreg3d(iim1,j,k) * &
                     (alo_u*c26_slineu + abs_sc*c26_intu*aloline_on_nn3d(iim1,j,k,q14))

                  aloline_on_nn3d(i,jjm1,k,q17) = imask_totreg3d(i,jjm1,k) * &
                     (alo_u*c24_slineu + abs_sc*c24_intu*aloline_on_nn3d(i,jjm1,k,q14))
                  aloline_on_nn3d(iim1,jjm1,k,q18) = imask_totreg3d(iim1,jjm1,k) * &
                     (alo_u*c23_slineu + abs_sc*(c23_intu*aloline_on_nn3d(iim1,jjm1,k,q14) + &
                     c26_intu*aloline_on_nn3d(iim1,jjm1,k,q17)))
                  aloline_on_nn3d(i,j,kkm1,q23) = imask_totreg3d(i,j,kkm1) * &
                     (alo_u*c18_slineu + abs_sc*c18_intu*aloline_on_nn3d(i,j,kkm1,q14))
                  aloline_on_nn3d(iim1,j,kkm1,q24) = imask_totreg3d(iim1,j,kkm1) * &
                     (alo_u*c17_slineu + abs_sc*(c17_intu*aloline_on_nn3d(iim1,j,kkm1,q14) + &
                     c18_intu*aloline_on_nn3d(iim1,j,kkm1,q15) + &
                     c26_intu*aloline_on_nn3d(iim1,j,kkm1,q23)))
                  aloline_on_nn3d(i,jjm1,kkm1,q26) = imask_totreg3d(i,jjm1,kkm1) * &
                     (alo_u*c15_slineu + abs_sc*(c15_intu*aloline_on_nn3d(i,jjm1,kkm1,q14) + &
                     c18_intu*aloline_on_nn3d(i,jjm1,kkm1,q17) + &
                     c24_intu*aloline_on_nn3d(i,jjm1,kkm1,q23)))
                  aloline_on_nn3d(iim1,jjm1,kkm1,q27) = imask_totreg3d(iim1,jjm1,kkm1) * &
                     (alo_u*c14_slineu + abs_sc*(c14_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q14) + &
                     c15_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q15) + &
                     c17_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q17) + &
                     c18_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q18) + &
                     c23_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q23) + &
                     c24_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q24) + &
                     c26_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q26)))
!
!integration
                  call calc_phinorm(vel_p, vth_p, vth_fiducial, xobs, phinorm)
                  wall=wall_global*phinorm
!
                  mintbar3d_tmp(i,j,k) = mintbar3d_tmp(i,j,k) + int3d(i,j,k)*wall
                  normalization3d_tmp(i,j,k) = normalization3d_tmp(i,j,k) + wall
!
                  aloline_nn3d_tmp(i,j,k,q14) = aloline_nn3d_tmp(i,j,k,q14) + wall*aloline_on_nn3d(i,j,k,q14)
                  aloline_nn3d_tmp(iim1,j,k,q15) = aloline_nn3d_tmp(iim1,j,k,q15) + wall*aloline_on_nn3d(iim1,j,k,q15)
                  aloline_nn3d_tmp(i,jjm1,k,q17) = aloline_nn3d_tmp(i,jjm1,k,q17) + wall*aloline_on_nn3d(i,jjm1,k,q17)
                  aloline_nn3d_tmp(iim1,jjm1,k,q18) = aloline_nn3d_tmp(iim1,jjm1,k,q18) + wall*aloline_on_nn3d(iim1,jjm1,k,q18)
                  aloline_nn3d_tmp(i,j,kkm1,q23) = aloline_nn3d_tmp(i,j,kkm1,q23) + wall*aloline_on_nn3d(i,j,kkm1,q23)
                  aloline_nn3d_tmp(iim1,j,kkm1,q24) = aloline_nn3d_tmp(iim1,j,kkm1,q24) + wall*aloline_on_nn3d(iim1,j,kkm1,q24)
                  aloline_nn3d_tmp(i,jjm1,kkm1,q26) = aloline_nn3d_tmp(i,jjm1,kkm1,q26) + wall*aloline_on_nn3d(i,jjm1,kkm1,q26)
                  aloline_nn3d_tmp(iim1,jjm1,kkm1,q27) = aloline_nn3d_tmp(iim1,jjm1,kkm1,q27) + wall*aloline_on_nn3d(iim1,jjm1,kkm1,q27)
!
               endif
!
             case default
!
            end select
!
!
         enddo
      enddo
   enddo
!
!
!
end subroutine fsc_line3d_lin
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine fsc_line3d_linb(oindx,xobsindx)
!
!-----------------------------------------------------------------------
!--------short characteristics for line radiative transfer in 3d--------
!-----------calculating intensties for given mu,phi specified-----------
!-----------by input oindx, and for xobs specified by xobsindx----------
!---------using linear interpolations for upwind/downwind point---------
!----------------and bezier integration along ray segment---------------
!-----------------------------------------------------------------------
!
   use prog_type
   use fund_const
   use dime3d, only: ndxmax, ndzmax, ndymax, x, y, z
   use dime3d, only: int3d, opalbar3d, sline3d, velx3d, vely3d, velz3d, vth3d, imask3d, imask_totreg3d, &
      aloline_on_nn3d, aloline_nn3d_tmp, mintbar3d_tmp, normalization3d_tmp
   use angles, only: n_x, n_y, n_z, weight_omega, q_alo
   use freq, only: nodes_xobs, weight_xobs, nxobs, deltax, xcmf_min, xcmf_max
   use params_input, only: vth_fiducial
   use bcondition, only: xic1
   use mod_debug, only: indxx, indxy, indxz
   use params_stellar, only: smajorax_a, smajorax_b, smajorax_c
!
   implicit none
!
! ... arguments
   integer(i4b), intent(in) :: oindx, xobsindx
!
! ... local scalars
   integer(i4b) :: i, j, k
   integer(i4b) :: alpha, beta, gamma
   integer(i4b) :: startx, starty, startz, endx, endy, endz
   integer(i4b) :: iim2, iim1, ii, iip1, jjm2, jjm1, jj, jjp1, kkm2, kkm1, kk, kkp1
   real(dp) :: nn_x, nn_y, nn_z, mueff, xobs, wall, wall_global, phinorm
   real(dp) :: dels_xyu, dels_xzu, dels_yzu, dels_u, dels_r, &
      dels_xyd, dels_xzd, dels_yzd, dels_d
   real(dp) :: opalbar_p, velx_p, vely_p, velz_p, vel_p, vth_p, sline_p, rad_p
   real(dp) :: x_u, y_u, z_u, opalbar_u, velx_u, vely_u, velz_u, vel_u, vth_u, sline_u, int_u, rad_u
   real(dp) :: x_d, y_d, z_d, opalbar_d, velx_d, vely_d, velz_d, vel_d, vth_d, sline_d, rad_d
   real(dp) :: abs_sc, int_sc, contr_sc
   real(dp) :: alo_u, alo_p, alo_d
   integer :: q1,  q2,  q3,  q4,  q5,  q6,  q7,  q8,  q9, &
      q10, q11, q12, q13, q14, q15, q16, q17, q18, &
      q19, q20, q21, q22, q23, q24, q25, q26, q27
   real(dp) :: c14_slineu, c15_slineu, c17_slineu, c18_slineu, &
      c23_slineu, c24_slineu, c26_slineu, c27_slineu, &
      c14_intu, c15_intu, c17_intu, c18_intu, &
      c23_intu, c24_intu, c26_intu, &
      c15_slined, c17_slined, c18_slined, c23_slined, &
      c24_slined, c26_slined, c27_slined
!
!not required
!real(dp) :: c25_slined, c22_slined, c21_slined, c20_slined, c19_slined, c16_slined, c13_slined, c12_slined, &
!            c09_slined, c08_slined, c07_slined, c06_slined, c03_slined, &
!            c22_intu, c20_intu, c16_intu, c13_intu, c12_intu, c11_intu, c10_intu, &
!            c08_intu, c06_intu, c05_intu, c04_intu, c03_intu, c02_intu

! ... local arrays
!
! ... local functions
   real(dp) :: dist_sphere, dist_ellipsoid, calc_icore_gdark
!
! ... local characters
!character(len=50) :: enter
!
! ... local logicals
!
!directions
   nn_x=n_x(oindx)
   nn_y=n_y(oindx)
   nn_z=n_z(oindx)
!
!frequency
   xobs=nodes_xobs(xobsindx)
!
!angular and frequency  integration weights
   wall_global=weight_omega(oindx)*weight_xobs(xobsindx)
!
!indices for nearest neighbour alo
   q1=q_alo(oindx,1)
   q2=q_alo(oindx,2)
   q3=q_alo(oindx,3)
   q4=q_alo(oindx,4)
   q5=q_alo(oindx,5)
   q6=q_alo(oindx,6)
   q7=q_alo(oindx,7)
   q8=q_alo(oindx,8)
   q9=q_alo(oindx,9)
   q10=q_alo(oindx,10)
   q11=q_alo(oindx,11)
   q12=q_alo(oindx,12)
   q13=q_alo(oindx,13)
   q14=q_alo(oindx,14)
   q15=q_alo(oindx,15)
   q16=q_alo(oindx,16)
   q17=q_alo(oindx,17)
   q18=q_alo(oindx,18)
   q19=q_alo(oindx,19)
   q20=q_alo(oindx,20)
   q21=q_alo(oindx,21)
   q22=q_alo(oindx,22)
   q23=q_alo(oindx,23)
   q24=q_alo(oindx,24)
   q25=q_alo(oindx,25)
   q26=q_alo(oindx,26)
   q27=q_alo(oindx,27)
!
!-----------------------------------------------------------------------
!
!set directional index-parameter (phantom points are excluded from calculation)
!
!index parameter:
!         if n_x >= 0                 if n_x < 0
!                startx = 2                  startx = ndxmax-1
!                endx = ndxmax-1             endx = 2
!                alpha=  1                   alpha=-1
!
!         if n_y >= 0                 if n_y < 0
!                starty = 2                  starty = ndymax-1
!                endy = ndymax-1             endy = 2
!                beta =  1                   beta =-1
!
!         if n_z >= 0                 if n_z < 0
!                startz = 2                  startz = ndzmax-1
!                endz = ndzmax-1             endz = 2
!                gamma = 1                   gamma = -1
!
   if(nn_x.gt.0.d0) then
      startx = 3
      endx = ndxmax-2
      alpha=  1
   elseif(nn_x.lt.0.d0) then
      startx = ndxmax-2
      endx = 3
      alpha=-1
   else
      stop 'error in fsc_line3d: n_x = 0 not allowed'
   endif
!
   if(nn_y.gt.0.d0) then
      starty = 3
      endy = ndymax-2
      beta =  1
   elseif(nn_y.lt.0.d0) then
      starty = ndymax-2
      endy = 3
      beta =-1
   else
      stop 'error in fsc_line3d: n_y = 0 not allowed'
   endif
!
   if(nn_z.gt.0.d0) then
      startz = 3
      endz = ndzmax-2
      gamma=  1
   elseif(nn_z.lt.0.d0) then
      startz = ndzmax-2
      endz = 3
      gamma=-1
   else
      stop 'error in fsc_line3d: n_z = 0 not allowed'
   endif
!
!--------------------reset the intensities and alo----------------------
!
   call set_boundary3d(nn_x, nn_y, nn_z)
!
   aloline_on_nn3d=0.d0
!
!-----------------------------------------------------------------------
!
   do k=startz, endz, gamma
      do j=starty, endy, beta
         do i=startx, endx, alpha
!
            select case(imask3d(i,j,k))
!
!************************standard radiative transfer********************
!
             case(1)
!
               iim1=i-alpha
               jjm1=j-beta
               kkm1=k-gamma
               iip1=i+alpha
               jjp1=j+beta
               kkp1=k+gamma
!
!calculate distance of upwind point to (previous) xy-plane, xz-plane and yz-plane
               dels_xyu=(z(k)-z(kkm1))/nn_z
               dels_xzu=(y(j)-y(jjm1))/nn_y
               dels_yzu=(x(i)-x(iim1))/nn_x
               dels_u=min(dels_xyu,dels_xzu,dels_yzu)
!
!calculate distance of downwind point to (next) xy-plane, xz-plane and yz-plane
               dels_xyd=(z(kkp1)-z(k))/nn_z
               dels_xzd=(y(jjp1)-y(j))/nn_y
               dels_yzd=(x(iip1)-x(i))/nn_x
               dels_d=min(dels_xyd,dels_xzd,dels_yzd)
!
!----------------------------local point--------------------------------
!
               sline_p=sline3d(i,j,k)
               opalbar_p=opalbar3d(i,j,k)
               velx_p=velx3d(i,j,k)
               vely_p=vely3d(i,j,k)
               velz_p=velz3d(i,j,k)
               vel_p = nn_x*velx_p + nn_y*vely_p + nn_z*velz_p
               vth_p = vth3d(i,j,k)
!
!----------------------------upwind point-------------------------------
!
               if(dels_xyu.eq.dels_u) then
!intersection with x-y plane on level k-gamma
                  x_u = x(i) - dels_u*nn_x
                  y_u = y(j) - dels_u*nn_y
                  z_u = z(kkm1)
!
                  call coeff2d_lineu_lin(opalbar3d(iim1,jjm1,kkm1), opalbar3d(i,jjm1,kkm1), opalbar3d(iim1,j,kkm1), opalbar3d(i,j,kkm1), &
                     velx3d(iim1,jjm1,kkm1), velx3d(i,jjm1,kkm1), velx3d(iim1,j,kkm1), velx3d(i,j,kkm1), &
                     vely3d(iim1,jjm1,kkm1), vely3d(i,jjm1,kkm1), vely3d(iim1,j,kkm1), vely3d(i,j,kkm1), &
                     velz3d(iim1,jjm1,kkm1), velz3d(i,jjm1,kkm1), velz3d(iim1,j,kkm1), velz3d(i,j,kkm1), &
                     vth3d(iim1,jjm1,kkm1), vth3d(i,jjm1,kkm1), vth3d(iim1,j,kkm1), vth3d(i,j,kkm1), &
                     sline3d(iim1,jjm1,kkm1), sline3d(i,jjm1,kkm1), sline3d(iim1,j,kkm1), sline3d(i,j,kkm1), &
                     int3d(iim1,jjm1,kkm1), int3d(i,jjm1,kkm1), int3d(iim1,j,kkm1), int3d(i,j,kkm1), &
                     x(iim1), x(i), y(jjm1), y(j), x_u, y_u, &
                     c14_slineu, c15_slineu, c17_slineu, c18_slineu, &
                     c14_intu, c15_intu, c17_intu, c18_intu, opalbar_u, velx_u, vely_u, velz_u, vth_u, sline_u, int_u)
                  c23_slineu=0.d0
                  c24_slineu=0.d0
                  c26_slineu=0.d0
                  c23_intu=0.d0
                  c24_intu=0.d0
                  c26_intu=0.d0
!
               elseif(dels_xzu.eq.dels_u) then
!intersection with x-z plane on level j-beta
                  x_u = x(i) - dels_u*nn_x
                  y_u = y(jjm1)
                  z_u = z(k) - dels_u*nn_z
!
                  call coeff2d_lineu_lin(opalbar3d(iim1,jjm1,kkm1), opalbar3d(i,jjm1,kkm1), opalbar3d(iim1,jjm1,k), opalbar3d(i,jjm1,k), &
                     velx3d(iim1,jjm1,kkm1), velx3d(i,jjm1,kkm1), velx3d(iim1,jjm1,k), velx3d(i,jjm1,k), &
                     vely3d(iim1,jjm1,kkm1), vely3d(i,jjm1,kkm1), vely3d(iim1,jjm1,k), vely3d(i,jjm1,k), &
                     velz3d(iim1,jjm1,kkm1), velz3d(i,jjm1,kkm1), velz3d(iim1,jjm1,k), velz3d(i,jjm1,k), &
                     vth3d(iim1,jjm1,kkm1), vth3d(i,jjm1,kkm1), vth3d(iim1,jjm1,k), vth3d(i,jjm1,k), &
                     sline3d(iim1,jjm1,kkm1), sline3d(i,jjm1,kkm1), sline3d(iim1,jjm1,k), sline3d(i,jjm1,k), &
                     int3d(iim1,jjm1,kkm1), int3d(i,jjm1,kkm1), int3d(iim1,jjm1,k), int3d(i,jjm1,k), &
                     x(iim1), x(i), z(kkm1), z(k), x_u, z_u, &
                     c14_slineu, c15_slineu, c23_slineu, c24_slineu, &
                     c14_intu, c15_intu, c23_intu, c24_intu, opalbar_u, velx_u, vely_u, velz_u, vth_u, sline_u, int_u)
                  c17_slineu=0.d0
                  c18_slineu=0.d0
                  c26_slineu=0.d0
                  c17_intu=0.d0
                  c18_intu=0.d0
                  c26_intu=0.d0
!
               elseif(dels_yzu.eq.dels_u) then
!intersection with y-z plane on level i-alpha
                  x_u = x(iim1)
                  y_u = y(j) - dels_u*nn_y
                  z_u = z(k) - dels_u*nn_z
!
                  call coeff2d_lineu_lin(opalbar3d(iim1,jjm1,kkm1), opalbar3d(iim1,j,kkm1), opalbar3d(iim1,jjm1,k), opalbar3d(iim1,j,k), &
                     velx3d(iim1,jjm1,kkm1), velx3d(iim1,j,kkm1), velx3d(iim1,jjm1,k), velx3d(iim1,j,k), &
                     vely3d(iim1,jjm1,kkm1), vely3d(iim1,j,kkm1), vely3d(iim1,jjm1,k), vely3d(iim1,j,k), &
                     velz3d(iim1,jjm1,kkm1), velz3d(iim1,j,kkm1), velz3d(iim1,jjm1,k), velz3d(iim1,j,k), &
                     vth3d(iim1,jjm1,kkm1), vth3d(iim1,j,kkm1), vth3d(iim1,jjm1,k), vth3d(iim1,j,k), &
                     sline3d(iim1,jjm1,kkm1), sline3d(iim1,j,kkm1), sline3d(iim1,jjm1,k), sline3d(iim1,j,k), &
                     int3d(iim1,jjm1,kkm1), int3d(iim1,j,kkm1), int3d(iim1,jjm1,k), int3d(iim1,j,k), &
                     y(jjm1), y(j), z(kkm1), z(k), y_u, z_u, &
                     c14_slineu, c17_slineu, c23_slineu, c26_slineu, &
                     c14_intu, c17_intu, c23_intu, c26_intu, opalbar_u, velx_u, vely_u, velz_u, vth_u, sline_u, int_u)
                  c15_slineu=0.d0
                  c18_slineu=0.d0
                  c24_slineu=0.d0
                  c15_intu=0.d0
                  c18_intu=0.d0
                  c24_intu=0.d0
!
               else
                  write(*,'(4es20.8,2l4)') dels_u, dels_xzu, dels_xyu, dels_yzu, dels_u.eq.dels_xzu, dels_xyu.eq.dels_u
                  stop 'error in fsc_line3d: invalid dels_u'
               endif
!
!---------------------------downwind point------------------------------
!
               if(dels_xyd.eq.dels_d) then
!intersection with x-y plane on level k+gamma
                  x_d = x(i) + dels_d*nn_x
                  y_d = y(j) + dels_d*nn_y
                  z_d = z(kkp1)
!
                  call coeff2d_lined_lin(opalbar3d(i,j,kkp1), opalbar3d(iip1,j,kkp1), opalbar3d(i,jjp1,kkp1), opalbar3d(iip1,jjp1,kkp1), &
                     velx3d(i,j,kkp1), velx3d(iip1,j,kkp1), velx3d(i,jjp1,kkp1), velx3d(iip1,jjp1,kkp1), &
                     vely3d(i,j,kkp1), vely3d(iip1,j,kkp1), vely3d(i,jjp1,kkp1), vely3d(iip1,jjp1,kkp1), &
                     velz3d(i,j,kkp1), velz3d(iip1,j,kkp1), velz3d(i,jjp1,kkp1), velz3d(iip1,jjp1,kkp1), &
                     vth3d(i,j,kkp1), vth3d(iip1,j,kkp1), vth3d(i,jjp1,kkp1), vth3d(iip1,jjp1,kkp1), &
                     sline3d(i,j,kkp1), sline3d(iip1,j,kkp1), sline3d(i,jjp1,kkp1), sline3d(iip1,jjp1,kkp1), &
                     x(i), x(iip1), y(j), y(jjp1), x_d, y_d, &
                     c23_slined, c24_slined, c26_slined, c27_slined, &
                     opalbar_d, velx_d, vely_d, velz_d, vth_d, sline_d)
!set interpolation coefficients that are not used to zero
                  c15_slined = 0.d0
                  c17_slined = 0.d0
                  c18_slined = 0.d0
!
               elseif(dels_xzd.eq.dels_d) then
!intersection with x-z plane on level j+beta
                  x_d = x(i) + dels_d*nn_x
                  y_d = y(jjp1)
                  z_d = z(k) + dels_d*nn_z
!
                  call coeff2d_lined_lin(opalbar3d(i,jjp1,k), opalbar3d(iip1,jjp1,k), opalbar3d(i,jjp1,kkp1), opalbar3d(iip1,jjp1,kkp1), &
                     velx3d(i,jjp1,k), velx3d(iip1,jjp1,k), velx3d(i,jjp1,kkp1), velx3d(iip1,jjp1,kkp1), &
                     vely3d(i,jjp1,k), vely3d(iip1,jjp1,k), vely3d(i,jjp1,kkp1), vely3d(iip1,jjp1,kkp1), &
                     velz3d(i,jjp1,k), velz3d(iip1,jjp1,k), velz3d(i,jjp1,kkp1), velz3d(iip1,jjp1,kkp1), &
                     vth3d(i,jjp1,k), vth3d(iip1,jjp1,k), vth3d(i,jjp1,kkp1), vth3d(iip1,jjp1,kkp1), &
                     sline3d(i,jjp1,k), sline3d(iip1,jjp1,k), sline3d(i,jjp1,kkp1), sline3d(iip1,jjp1,kkp1), &
                     x(i), x(iip1), z(k), z(kkp1), x_d, z_d, &
                     c17_slined, c18_slined, c26_slined, c27_slined, &
                     opalbar_d, velx_d, vely_d, velz_d, vth_d, sline_d)
!set interpolation coefficients that are not used to zero
                  c15_slined = 0.d0
                  c23_slined = 0.d0
                  c24_slined = 0.d0
!
               elseif(dels_yzd.eq.dels_d) then
!intersection with y-z plane on level i+alpha
                  x_d = x(iip1)
                  y_d = y(j) + dels_d*nn_y
                  z_d = z(k) + dels_d*nn_z
!
                  call coeff2d_lined_lin(opalbar3d(iip1,j,k), opalbar3d(iip1,jjp1,k), opalbar3d(iip1,j,kkp1), opalbar3d(iip1,jjp1,kkp1), &
                     velx3d(iip1,j,k), velx3d(iip1,jjp1,k), velx3d(iip1,j,kkp1), velx3d(iip1,jjp1,kkp1), &
                     vely3d(iip1,j,k), vely3d(iip1,jjp1,k), vely3d(iip1,j,kkp1), vely3d(iip1,jjp1,kkp1), &
                     velz3d(iip1,j,k), velz3d(iip1,jjp1,k), velz3d(iip1,j,kkp1), velz3d(iip1,jjp1,kkp1), &
                     vth3d(iip1,j,k), vth3d(iip1,jjp1,k), vth3d(iip1,j,kkp1), vth3d(iip1,jjp1,kkp1), &
                     sline3d(iip1,j,k), sline3d(iip1,jjp1,k), sline3d(iip1,j,kkp1), sline3d(iip1,jjp1,kkp1), &
                     y(j), y(jjp1), z(k), z(kkp1), y_d, z_d, &
                     c15_slined, c18_slined, c24_slined, c27_slined, &
                     opalbar_d, velx_d, vely_d, velz_d, vth_d, sline_d)
!set interpolation coefficients that are not used to zero
                  c17_slined = 0.d0
                  c23_slined = 0.d0
                  c26_slined = 0.d0
!
               else
                  write(*,'(4es20.8,l4)') dels_d, dels_xzd, dels_xyd, dels_yzd, dels_d.eq.dels_xzd
                  stop 'error in fsc_line3d: invalid dels_d'
               endif
!
!--------------------------------radiative transfer---------------------
!
               vel_u = velx_u*nn_x + vely_u*nn_y + velz_u*nn_z
               vel_d = velx_d*nn_x + vely_d*nn_y + velz_d*nn_z
               call fsc_line(xobs, vth_fiducial, deltax, xcmf_min, xcmf_max, int_u, opalbar_u, opalbar_p, opalbar_d, &
                  sline_u, sline_p, sline_d, vel_u, vel_p, vel_d, vth_u, vth_p, vth_d, &
                  dels_u, dels_d, abs_sc, contr_sc, int_sc, alo_u, alo_p, alo_d)
!                  call fsc_line_lin(xobs, vth_fiducial, deltax, xcmf_min, xcmf_max, int_u, opalbar_u, opalbar_p, sline_u, sline_p, &
!                                      vel_u, vel_p, vth_u, vth_p, dels_u, abs_sc, contr_sc, int_sc, alo_u, alo_p)
!                 alo_d=0.d0
               int3d(i,j,k) = abs_sc*int_u + alo_u*sline_u + alo_p*sline_p + alo_d*sline_d
!
               aloline_on_nn3d(iip1,jjp1,kkp1,q1) = imask_totreg3d(iip1,jjp1,kkp1)*&
                  alo_d*c27_slined
               aloline_on_nn3d(i,jjp1,kkp1,q2) = imask_totreg3d(i,jjp1,kkp1) * &
                  (alo_d*c26_slined + abs_sc*c26_intu*aloline_on_nn3d(i,jjp1,kkp1,q1))
               aloline_on_nn3d(iim1,jjp1,kkp1,q3) = imask_totreg3d(iim1,jjp1,kkp1) * &
                  (abs_sc*c26_intu*aloline_on_nn3d(iim1,jjp1,kkp1,q2))
               aloline_on_nn3d(iip1,j,kkp1,q4) = imask_totreg3d(iip1,j,kkp1) * &
                  (alo_d*c24_slined + abs_sc*c24_intu*aloline_on_nn3d(iip1,j,kkp1,q1))
               aloline_on_nn3d(i,j,kkp1,q5) = imask_totreg3d(i,j,kkp1) * &
                  (alo_d*c23_slined + abs_sc*(c23_intu*aloline_on_nn3d(i,j,kkp1,q1) + &
                  c24_intu*aloline_on_nn3d(i,j,kkp1,q2) + &
                  c26_intu*aloline_on_nn3d(i,j,kkp1,q4)))
               aloline_on_nn3d(iim1,j,kkp1,q6) = imask_totreg3d(iim1,j,kkp1) * &
                  (abs_sc*(c23_intu*aloline_on_nn3d(iim1,j,kkp1,q2) + &
                  c24_intu*aloline_on_nn3d(iim1,j,kkp1,q3) + &
                  c26_intu*aloline_on_nn3d(iim1,j,kkp1,q5)))
               aloline_on_nn3d(iip1,jjm1,kkp1,q7) = imask_totreg3d(iip1,jjm1,kkp1) * &
                  (abs_sc*c24_intu*aloline_on_nn3d(iip1,jjm1,kkp1,q4))
               aloline_on_nn3d(i,jjm1,kkp1,q8) = imask_totreg3d(i,jjm1,kkp1) * &
                  (abs_sc*(c23_intu*aloline_on_nn3d(i,jjm1,kkp1,q4) + &
                  c24_intu*aloline_on_nn3d(i,jjm1,kkp1,q5) + &
                  c26_intu*aloline_on_nn3d(i,jjm1,kkp1,q7)))
               aloline_on_nn3d(iim1,jjm1,kkp1,q9) = imask_totreg3d(iim1,jjm1,kkp1) * &
                  (abs_sc*(c23_intu*aloline_on_nn3d(iim1,jjm1,kkp1,q5) + &
                  c24_intu*aloline_on_nn3d(iim1,jjm1,kkp1,q6) + &
                  c26_intu*aloline_on_nn3d(iim1,jjm1,kkp1,q8)))
               aloline_on_nn3d(iip1,jjp1,k,q10) = imask_totreg3d(iip1,jjp1,k) * &
                  (alo_d*c18_slined + abs_sc*c18_intu*aloline_on_nn3d(iip1,jjp1,k,q1))
               aloline_on_nn3d(i,jjp1,k,q11) = imask_totreg3d(i,jjp1,k) * &
                  (alo_d*c17_slined + abs_sc*(c17_intu*aloline_on_nn3d(i,jjp1,k,q1) + &
                  c18_intu*aloline_on_nn3d(i,jjp1,k,q2) + &
                  c26_intu*aloline_on_nn3d(i,jjp1,k,q10)))
               aloline_on_nn3d(iim1,jjp1,k,q12) = imask_totreg3d(iim1,jjp1,k) * &
                  (abs_sc*(c17_intu*aloline_on_nn3d(iim1,jjp1,k,q2) + &
                  c18_intu*aloline_on_nn3d(iim1,jjp1,k,q3) + &
                  c26_intu*aloline_on_nn3d(iim1,jjp1,k,q11)))
               aloline_on_nn3d(iip1,j,k,q13) = imask_totreg3d(iip1,j,k) * &
                  (alo_d*c15_slined + abs_sc*(c15_intu*aloline_on_nn3d(iip1,j,k,q1) + &
                  c18_intu*aloline_on_nn3d(iip1,j,k,q4) + &
                  c24_intu*aloline_on_nn3d(iip1,j,k,q10)))
               aloline_on_nn3d(i,j,k,q14) = imask_totreg3d(i,j,k) * &
                  (alo_p + abs_sc*(c14_intu*aloline_on_nn3d(i,j,k,q1) + &
                  c15_intu*aloline_on_nn3d(i,j,k,q2) + &
                  c17_intu*aloline_on_nn3d(i,j,k,q4) + &
                  c26_intu*aloline_on_nn3d(i,j,k,q13) + &
                  c18_intu*aloline_on_nn3d(i,j,k,q5) + &
                  c23_intu*aloline_on_nn3d(i,j,k,q10) + &
                  c24_intu*aloline_on_nn3d(i,j,k,q11)))
               aloline_on_nn3d(iim1,j,k,q15) = imask_totreg3d(iim1,j,k) * &
                  (alo_u*c26_slineu + abs_sc*(c14_intu*aloline_on_nn3d(iim1,j,k,q2) + &
                  c15_intu*aloline_on_nn3d(iim1,j,k,q3) + &
                  c17_intu*aloline_on_nn3d(iim1,j,k,q5) + &
                  c18_intu*aloline_on_nn3d(iim1,j,k,q6) + &
                  c23_intu*aloline_on_nn3d(iim1,j,k,q11) + &
                  c24_intu*aloline_on_nn3d(iim1,j,k,q12) + &
                  c26_intu*aloline_on_nn3d(iim1,j,k,q14)))
               aloline_on_nn3d(iip1,jjm1,k,q16) = imask_totreg3d(iip1,jjm1,k) * &
                  (abs_sc*(c15_intu*aloline_on_nn3d(iip1,jjm1,k,q4) + &
                  c18_intu*aloline_on_nn3d(iip1,jjm1,k,q7) + &
                  c24_intu*aloline_on_nn3d(iip1,jjm1,k,q13)))
               aloline_on_nn3d(i,jjm1,k,q17) = imask_totreg3d(i,jjm1,k) * &
                  (alo_u*c24_slineu + abs_sc*(c14_intu*aloline_on_nn3d(i,jjm1,k,q4) + &
                  c15_intu*aloline_on_nn3d(i,jjm1,k,q5) + &
                  c17_intu*aloline_on_nn3d(i,jjm1,k,q7) + &
                  c18_intu*aloline_on_nn3d(i,jjm1,k,q8) + &
                  c23_intu*aloline_on_nn3d(i,jjm1,k,q13) + &
                  c24_intu*aloline_on_nn3d(i,jjm1,k,q14) + &
                  c26_intu*aloline_on_nn3d(i,jjm1,k,q16)))
               aloline_on_nn3d(iim1,jjm1,k,q18) = imask_totreg3d(iim1,jjm1,k) * &
                  (alo_u*c23_slineu + abs_sc*(c14_intu*aloline_on_nn3d(iim1,jjm1,k,q5) + &
                  c15_intu*aloline_on_nn3d(iim1,jjm1,k,q6) + &
                  c17_intu*aloline_on_nn3d(iim1,jjm1,k,q8) + &
                  c18_intu*aloline_on_nn3d(iim1,jjm1,k,q9) + &
                  c23_intu*aloline_on_nn3d(iim1,jjm1,k,q14) + &
                  c24_intu*aloline_on_nn3d(iim1,jjm1,k,q15) + &
                  c26_intu*aloline_on_nn3d(iim1,jjm1,k,q17)))
               aloline_on_nn3d(iip1,jjp1,kkm1,q19) = imask_totreg3d(iip1,jjp1,kkm1) * &
                  (abs_sc*c18_intu*aloline_on_nn3d(iip1,jjp1,kkm1,q10))
               aloline_on_nn3d(i,jjp1,kkm1,q20) = imask_totreg3d(i,jjp1,kkm1) * &
                  (abs_sc*(c17_intu*aloline_on_nn3d(i,jjp1,kkm1,q10) + &
                  c18_intu*aloline_on_nn3d(i,jjp1,kkm1,q11) + &
                  c26_intu*aloline_on_nn3d(i,jjp1,kkm1,q19)))
               aloline_on_nn3d(iim1,jjp1,kkm1,q21) = imask_totreg3d(iim1,jjp1,kkm1) * &
                  (abs_sc*(c17_intu*aloline_on_nn3d(iim1,jjp1,kkm1,q11) + &
                  c18_intu*aloline_on_nn3d(iim1,jjp1,kkm1,q12) + &
                  c26_intu*aloline_on_nn3d(iim1,jjp1,kkm1,q20)))
               aloline_on_nn3d(iip1,j,kkm1,q22) = imask_totreg3d(iip1,j,kkm1) * &
                  (abs_sc*(c15_intu*aloline_on_nn3d(iip1,j,kkm1,q10) + &
                  c18_intu*aloline_on_nn3d(iip1,j,kkm1,q13) + &
                  c24_intu*aloline_on_nn3d(iip1,j,kkm1,q19)))
               aloline_on_nn3d(i,j,kkm1,q23) = imask_totreg3d(i,j,kkm1) * &
                  (alo_u*c18_slineu + abs_sc*(c14_intu*aloline_on_nn3d(i,j,kkm1,q10) + &
                  c15_intu*aloline_on_nn3d(i,j,kkm1,q11) + &
                  c17_intu*aloline_on_nn3d(i,j,kkm1,q13) + &
                  c18_intu*aloline_on_nn3d(i,j,kkm1,q14) + &
                  c23_intu*aloline_on_nn3d(i,j,kkm1,q19) + &
                  c24_intu*aloline_on_nn3d(i,j,kkm1,q20) + &
                  c26_intu*aloline_on_nn3d(i,j,kkm1,q22)))
               aloline_on_nn3d(iim1,j,kkm1,q24) = imask_totreg3d(iim1,j,kkm1) * &
                  (alo_u*c17_slineu + abs_sc*(c14_intu*aloline_on_nn3d(iim1,j,kkm1,q11) + &
                  c15_intu*aloline_on_nn3d(iim1,j,kkm1,q12) + &
                  c17_intu*aloline_on_nn3d(iim1,j,kkm1,q14) + &
                  c18_intu*aloline_on_nn3d(iim1,j,kkm1,q15) + &
                  c23_intu*aloline_on_nn3d(iim1,j,kkm1,q20) + &
                  c24_intu*aloline_on_nn3d(iim1,j,kkm1,q21) + &
                  c26_intu*aloline_on_nn3d(iim1,j,kkm1,q23)))
               aloline_on_nn3d(iip1,jjm1,kkm1,q25) = imask_totreg3d(iip1,jjm1,kkm1) * &
                  (abs_sc*(c15_intu*aloline_on_nn3d(iip1,jjm1,kkm1,q13) + &
                  c18_intu*aloline_on_nn3d(iip1,jjm1,kkm1,q16) + &
                  c24_intu*aloline_on_nn3d(iip1,jjm1,kkm1,q22)))
               aloline_on_nn3d(i,jjm1,kkm1,q26) = imask_totreg3d(i,jjm1,kkm1) * &
                  (alo_u*c15_slineu + abs_sc*(c14_intu*aloline_on_nn3d(i,jjm1,kkm1,q13) + &
                  c15_intu*aloline_on_nn3d(i,jjm1,kkm1,q14) + &
                  c17_intu*aloline_on_nn3d(i,jjm1,kkm1,q16) + &
                  c18_intu*aloline_on_nn3d(i,jjm1,kkm1,q17) + &
                  c23_intu*aloline_on_nn3d(i,jjm1,kkm1,q22) + &
                  c24_intu*aloline_on_nn3d(i,jjm1,kkm1,q23) + &
                  c26_intu*aloline_on_nn3d(i,jjm1,kkm1,q25)))
               aloline_on_nn3d(iim1,jjm1,kkm1,q27) = imask_totreg3d(iim1,jjm1,kkm1) * &
                  (alo_u*c14_slineu + abs_sc*(c14_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q14) + &
                  c15_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q15) + &
                  c17_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q17) + &
                  c18_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q18) + &
                  c23_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q23) + &
                  c24_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q24) + &
                  c26_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q26)))

!
!perform angular integration
               call calc_phinorm(vel_p, vth_p, vth_fiducial, xobs, phinorm)
               wall=wall_global*phinorm
               mintbar3d_tmp(i,j,k) = mintbar3d_tmp(i,j,k) + int3d(i,j,k)*wall
               normalization3d_tmp(i,j,k) = normalization3d_tmp(i,j,k) + wall
!
               aloline_nn3d_tmp(iip1,jjp1,kkp1,q1) = aloline_nn3d_tmp(iip1,jjp1,kkp1,q1) + wall*aloline_on_nn3d(iip1,jjp1,kkp1,q1)
               aloline_nn3d_tmp(i,jjp1,kkp1,q2) = aloline_nn3d_tmp(i,jjp1,kkp1,q2) + wall*aloline_on_nn3d(i,jjp1,kkp1,q2)
               aloline_nn3d_tmp(iim1,jjp1,kkp1,q3) = aloline_nn3d_tmp(iim1,jjp1,kkp1,q3) + wall*aloline_on_nn3d(iim1,jjp1,kkp1,q3)
               aloline_nn3d_tmp(iip1,j,kkp1,q4) = aloline_nn3d_tmp(iip1,j,kkp1,q4) + wall*aloline_on_nn3d(iip1,j,kkp1,q4)
               aloline_nn3d_tmp(i,j,kkp1,q5) = aloline_nn3d_tmp(i,j,kkp1,q5) + wall*aloline_on_nn3d(i,j,kkp1,q5)
               aloline_nn3d_tmp(iim1,j,kkp1,q6) = aloline_nn3d_tmp(iim1,j,kkp1,q6) + wall*aloline_on_nn3d(iim1,j,kkp1,q6)
               aloline_nn3d_tmp(iip1,jjm1,kkp1,q7) = aloline_nn3d_tmp(iip1,jjm1,kkp1,q7) + wall*aloline_on_nn3d(iip1,jjm1,kkp1,q7)
               aloline_nn3d_tmp(i,jjm1,kkp1,q8) = aloline_nn3d_tmp(i,jjm1,kkp1,q8) + wall*aloline_on_nn3d(i,jjm1,kkp1,q8)
               aloline_nn3d_tmp(iim1,jjm1,kkp1,q9) = aloline_nn3d_tmp(iim1,jjm1,kkp1,q9) + wall*aloline_on_nn3d(iim1,jjm1,kkp1,q9)
               aloline_nn3d_tmp(iip1,jjp1,k,q10) = aloline_nn3d_tmp(iip1,jjp1,k,q10) + wall*aloline_on_nn3d(iip1,jjp1,k,q10)
               aloline_nn3d_tmp(i,jjp1,k,q11) = aloline_nn3d_tmp(i,jjp1,k,q11) + wall*aloline_on_nn3d(i,jjp1,k,q11)
               aloline_nn3d_tmp(iim1,jjp1,k,q12) = aloline_nn3d_tmp(iim1,jjp1,k,q12) + wall*aloline_on_nn3d(iim1,jjp1,k,q12)
               aloline_nn3d_tmp(iip1,j,k,q13) = aloline_nn3d_tmp(iip1,j,k,q13) + wall*aloline_on_nn3d(iip1,j,k,q13)
               aloline_nn3d_tmp(i,j,k,q14) = aloline_nn3d_tmp(i,j,k,q14) + wall*aloline_on_nn3d(i,j,k,q14)
               aloline_nn3d_tmp(iim1,j,k,q15) = aloline_nn3d_tmp(iim1,j,k,q15) + wall*aloline_on_nn3d(iim1,j,k,q15)
               aloline_nn3d_tmp(iip1,jjm1,k,q16) = aloline_nn3d_tmp(iip1,jjm1,k,q16) + wall*aloline_on_nn3d(iip1,jjm1,k,q16) !new
               aloline_nn3d_tmp(i,jjm1,k,q17) = aloline_nn3d_tmp(i,jjm1,k,q17) + wall*aloline_on_nn3d(i,jjm1,k,q17)
               aloline_nn3d_tmp(iim1,jjm1,k,q18) = aloline_nn3d_tmp(iim1,jjm1,k,q18) + wall*aloline_on_nn3d(iim1,jjm1,k,q18)
               aloline_nn3d_tmp(iip1,jjp1,kkm1,q19) = aloline_nn3d_tmp(iip1,jjp1,kkm1,q19) + wall*aloline_on_nn3d(iip1,jjp1,kkm1,q19) !new
               aloline_nn3d_tmp(i,jjp1,kkm1,q20) = aloline_nn3d_tmp(i,jjp1,kkm1,q20) + wall*aloline_on_nn3d(i,jjp1,kkm1,q20) !new
               aloline_nn3d_tmp(iim1,jjp1,kkm1,q21) = aloline_nn3d_tmp(iim1,jjp1,kkm1,q21) + wall*aloline_on_nn3d(iim1,jjp1,kkm1,q21) !new
               aloline_nn3d_tmp(iip1,j,kkm1,q22) = aloline_nn3d_tmp(iip1,j,kkm1,q22) + wall*aloline_on_nn3d(iip1,j,kkm1,q22) !new
               aloline_nn3d_tmp(i,j,kkm1,q23) = aloline_nn3d_tmp(i,j,kkm1,q23) + wall*aloline_on_nn3d(i,j,kkm1,q23)
               aloline_nn3d_tmp(iim1,j,kkm1,q24) = aloline_nn3d_tmp(iim1,j,kkm1,q24) + wall*aloline_on_nn3d(iim1,j,kkm1,q24)
               aloline_nn3d_tmp(iip1,jjm1,kkm1,q25) = aloline_nn3d_tmp(iip1,jjm1,kkm1,q25) + wall*aloline_on_nn3d(iip1,jjm1,kkm1,q25) !new
               aloline_nn3d_tmp(i,jjm1,kkm1,q26) = aloline_nn3d_tmp(i,jjm1,kkm1,q26) + wall*aloline_on_nn3d(i,jjm1,kkm1,q26)
               aloline_nn3d_tmp(iim1,jjm1,kkm1,q27) = aloline_nn3d_tmp(iim1,jjm1,kkm1,q27) + wall*aloline_on_nn3d(iim1,jjm1,kkm1,q27)
!
!***********special treatment: point is in vicinity of boundary*********
!
             case(2)
!
               iim1=i-alpha
               jjm1=j-beta
               kkm1=k-gamma
               iip1=i+alpha
               jjp1=j+beta
               kkp1=k+gamma
!
!calculate distance of upwind point to (previous) xy-plane, xz-plane and yz-plane
               dels_xyu=(z(k)-z(kkm1))/nn_z
               dels_xzu=(y(j)-y(jjm1))/nn_y
               dels_yzu=(x(i)-x(iim1))/nn_x
!
!calculate distance of downwind point to (next) xy-plane, xz-plane and yz-plane
               dels_xyd=(z(kkp1)-z(k))/nn_z
               dels_xzd=(y(jjp1)-y(j))/nn_y
               dels_yzd=(x(iip1)-x(i))/nn_x
!
!calculate distance to the boundary (here: sphere)
!               dels_r=dist_sphere(nn_x,nn_y,nn_z,x(i),y(j),z(k))
               dels_r=dist_ellipsoid(nn_x,nn_y,nn_z,x(i),y(j),z(k),smajorax_a,smajorax_b,smajorax_c)
!
               dels_u=min(dels_xyu,dels_xzu,dels_yzu,dels_r)
               dels_d=min(dels_xyd,dels_xzd,dels_yzd)
!
!----------------------------local point--------------------------------
!
               sline_p=sline3d(i,j,k)
               opalbar_p=opalbar3d(i,j,k)
               velx_p=velx3d(i,j,k)
               vely_p=vely3d(i,j,k)
               velz_p=velz3d(i,j,k)
               vel_p = nn_x*velx_p + nn_y*vely_p + nn_z*velz_p
               vth_p = vth3d(i,j,k)
!
!----------------------------upwind point-------------------------------
!
               if(dels_r.eq.dels_u) then
                  x_u = x(i) - dels_u*nn_x
                  y_u = y(j) - dels_u*nn_y
                  z_u = z(k) - dels_u*nn_z
!
!for alo
!                  int_u=0.d0
!                  int_u=xic1
                  int_u=calc_icore_gdark(z_u, sqrt(x_u**2+y_u**2+z_u**2))
!
                  call coeff3d_lineu_lin(iim1, i, jjm1, j, kkm1, k, x_u, y_u, z_u, &
                     c14_slineu, c15_slineu, c17_slineu, c18_slineu, &
                     c23_slineu, c24_slineu, c26_slineu, c27_slineu, &
                     opalbar_u, velx_u, vely_u, velz_u, vth_u, sline_u)
!
               elseif(dels_xyu.eq.dels_u) then
!intersection with x-y plane on level k-gamma
                  x_u = x(i) - dels_u*nn_x
                  y_u = y(j) - dels_u*nn_y
                  z_u = z(kkm1)
!
                  call coeff2d_lineu_lin(opalbar3d(iim1,jjm1,kkm1), opalbar3d(i,jjm1,kkm1), opalbar3d(iim1,j,kkm1), opalbar3d(i,j,kkm1), &
                     velx3d(iim1,jjm1,kkm1), velx3d(i,jjm1,kkm1), velx3d(iim1,j,kkm1), velx3d(i,j,kkm1), &
                     vely3d(iim1,jjm1,kkm1), vely3d(i,jjm1,kkm1), vely3d(iim1,j,kkm1), vely3d(i,j,kkm1), &
                     velz3d(iim1,jjm1,kkm1), velz3d(i,jjm1,kkm1), velz3d(iim1,j,kkm1), velz3d(i,j,kkm1), &
                     vth3d(iim1,jjm1,kkm1), vth3d(i,jjm1,kkm1), vth3d(iim1,j,kkm1), vth3d(i,j,kkm1), &
                     sline3d(iim1,jjm1,kkm1), sline3d(i,jjm1,kkm1), sline3d(iim1,j,kkm1), sline3d(i,j,kkm1), &
                     int3d(iim1,jjm1,kkm1), int3d(i,jjm1,kkm1), int3d(iim1,j,kkm1), int3d(i,j,kkm1), &
                     x(iim1), x(i), y(jjm1), y(j), x_u, y_u, &
                     c14_slineu, c15_slineu, c17_slineu, c18_slineu, &
                     c14_intu, c15_intu, c17_intu, c18_intu, opalbar_u, velx_u, vely_u, velz_u, vth_u, sline_u, int_u)
                  c23_slineu=0.d0
                  c24_slineu=0.d0
                  c26_slineu=0.d0
                  c27_slineu=0.d0
                  c23_intu=0.d0
                  c24_intu=0.d0
                  c26_intu=0.d0
!
               elseif(dels_xzu.eq.dels_u) then
!intersection with x-z plane on level j-beta
                  x_u = x(i) - dels_u*nn_x
                  y_u = y(jjm1)
                  z_u = z(k) - dels_u*nn_z
!
                  call coeff2d_lineu_lin(opalbar3d(iim1,jjm1,kkm1), opalbar3d(i,jjm1,kkm1), opalbar3d(iim1,jjm1,k), opalbar3d(i,jjm1,k), &
                     velx3d(iim1,jjm1,kkm1), velx3d(i,jjm1,kkm1), velx3d(iim1,jjm1,k), velx3d(i,jjm1,k), &
                     vely3d(iim1,jjm1,kkm1), vely3d(i,jjm1,kkm1), vely3d(iim1,jjm1,k), vely3d(i,jjm1,k), &
                     velz3d(iim1,jjm1,kkm1), velz3d(i,jjm1,kkm1), velz3d(iim1,jjm1,k), velz3d(i,jjm1,k), &
                     vth3d(iim1,jjm1,kkm1), vth3d(i,jjm1,kkm1), vth3d(iim1,jjm1,k), vth3d(i,jjm1,k), &
                     sline3d(iim1,jjm1,kkm1), sline3d(i,jjm1,kkm1), sline3d(iim1,jjm1,k), sline3d(i,jjm1,k), &
                     int3d(iim1,jjm1,kkm1), int3d(i,jjm1,kkm1), int3d(iim1,jjm1,k), int3d(i,jjm1,k), &
                     x(iim1), x(i), z(kkm1), z(k), x_u, z_u, &
                     c14_slineu, c15_slineu, c23_slineu, c24_slineu, &
                     c14_intu, c15_intu, c23_intu, c24_intu, opalbar_u, velx_u, vely_u, velz_u, vth_u, sline_u, int_u)
                  c17_slineu=0.d0
                  c18_slineu=0.d0
                  c26_slineu=0.d0
                  c27_slineu=0.d0
                  c17_intu=0.d0
                  c18_intu=0.d0
                  c26_intu=0.d0
!
               elseif(dels_yzu.eq.dels_u) then
!intersection with y-z plane on level i-alpha
                  x_u = x(iim1)
                  y_u = y(j) - dels_u*nn_y
                  z_u = z(k) - dels_u*nn_z
!
                  call coeff2d_lineu_lin(opalbar3d(iim1,jjm1,kkm1), opalbar3d(iim1,j,kkm1), opalbar3d(iim1,jjm1,k), opalbar3d(iim1,j,k), &
                     velx3d(iim1,jjm1,kkm1), velx3d(iim1,j,kkm1), velx3d(iim1,jjm1,k), velx3d(iim1,j,k), &
                     vely3d(iim1,jjm1,kkm1), vely3d(iim1,j,kkm1), vely3d(iim1,jjm1,k), vely3d(iim1,j,k), &
                     velz3d(iim1,jjm1,kkm1), velz3d(iim1,j,kkm1), velz3d(iim1,jjm1,k), velz3d(iim1,j,k), &
                     vth3d(iim1,jjm1,kkm1), vth3d(iim1,j,kkm1), vth3d(iim1,jjm1,k), vth3d(iim1,j,k), &
                     sline3d(iim1,jjm1,kkm1), sline3d(iim1,j,kkm1), sline3d(iim1,jjm1,k), sline3d(iim1,j,k), &
                     int3d(iim1,jjm1,kkm1), int3d(iim1,j,kkm1), int3d(iim1,jjm1,k), int3d(iim1,j,k), &
                     y(jjm1), y(j), z(kkm1), z(k), y_u, z_u, &
                     c14_slineu, c17_slineu, c23_slineu, c26_slineu, &
                     c14_intu, c17_intu, c23_intu, c26_intu, opalbar_u, velx_u, vely_u, velz_u, vth_u, sline_u, int_u)
                  c15_slineu=0.d0
                  c18_slineu=0.d0
                  c24_slineu=0.d0
                  c27_slineu=0.d0
                  c15_intu=0.d0
                  c18_intu=0.d0
                  c24_intu=0.d0
!
               else
                  write(*,'(4es20.8,2l4)') dels_u, dels_xzu, dels_xyu, dels_yzu, dels_u.eq.dels_xzu, dels_xyu.eq.dels_u
                  stop 'error in fsc_line3d: invalid dels_u'
               endif
!
!---------------------------downwind point------------------------------
!
               if(dels_xyd.eq.dels_d) then
!intersection with x-y plane on level k+gamma
                  x_d = x(i) + dels_d*nn_x
                  y_d = y(j) + dels_d*nn_y
                  z_d = z(kkp1)
!
                  call coeff2d_lined_lin(opalbar3d(i,j,kkp1), opalbar3d(iip1,j,kkp1), opalbar3d(i,jjp1,kkp1), opalbar3d(iip1,jjp1,kkp1), &
                     velx3d(i,j,kkp1), velx3d(iip1,j,kkp1), velx3d(i,jjp1,kkp1), velx3d(iip1,jjp1,kkp1), &
                     vely3d(i,j,kkp1), vely3d(iip1,j,kkp1), vely3d(i,jjp1,kkp1), vely3d(iip1,jjp1,kkp1), &
                     velz3d(i,j,kkp1), velz3d(iip1,j,kkp1), velz3d(i,jjp1,kkp1), velz3d(iip1,jjp1,kkp1), &
                     vth3d(i,j,kkp1), vth3d(iip1,j,kkp1), vth3d(i,jjp1,kkp1), vth3d(iip1,jjp1,kkp1), &
                     sline3d(i,j,kkp1), sline3d(iip1,j,kkp1), sline3d(i,jjp1,kkp1), sline3d(iip1,jjp1,kkp1), &
                     x(i), x(iip1), y(j), y(jjp1), x_d, y_d, &
                     c23_slined, c24_slined, c26_slined, c27_slined, &
                     opalbar_d, velx_d, vely_d, velz_d, vth_d, sline_d)
!set interpolation coefficients that are not used to zero
                  c15_slined = 0.d0
                  c17_slined = 0.d0
                  c18_slined = 0.d0
!
               elseif(dels_xzd.eq.dels_d) then
!intersection with x-z plane on level j+beta
                  x_d = x(i) + dels_d*nn_x
                  y_d = y(jjp1)
                  z_d = z(k) + dels_d*nn_z
!
                  call coeff2d_lined_lin(opalbar3d(i,jjp1,k), opalbar3d(iip1,jjp1,k), opalbar3d(i,jjp1,kkp1), opalbar3d(iip1,jjp1,kkp1), &
                     velx3d(i,jjp1,k), velx3d(iip1,jjp1,k), velx3d(i,jjp1,kkp1), velx3d(iip1,jjp1,kkp1), &
                     vely3d(i,jjp1,k), vely3d(iip1,jjp1,k), vely3d(i,jjp1,kkp1), vely3d(iip1,jjp1,kkp1), &
                     velz3d(i,jjp1,k), velz3d(iip1,jjp1,k), velz3d(i,jjp1,kkp1), velz3d(iip1,jjp1,kkp1), &
                     vth3d(i,jjp1,k), vth3d(iip1,jjp1,k), vth3d(i,jjp1,kkp1), vth3d(iip1,jjp1,kkp1), &
                     sline3d(i,jjp1,k), sline3d(iip1,jjp1,k), sline3d(i,jjp1,kkp1), sline3d(iip1,jjp1,kkp1), &
                     x(i), x(iip1), z(k), z(kkp1), x_d, z_d, &
                     c17_slined, c18_slined, c26_slined, c27_slined, &
                     opalbar_d, velx_d, vely_d, velz_d, vth_d, sline_d)
!set interpolation coefficients that are not used to zero
                  c15_slined = 0.d0
                  c23_slined = 0.d0
                  c24_slined = 0.d0
!
               elseif(dels_yzd.eq.dels_d) then
!intersection with y-z plane on level i+alpha
                  x_d = x(iip1)
                  y_d = y(j) + dels_d*nn_y
                  z_d = z(k) + dels_d*nn_z
!
                  call coeff2d_lined_lin(opalbar3d(iip1,j,k), opalbar3d(iip1,jjp1,k), opalbar3d(iip1,j,kkp1), opalbar3d(iip1,jjp1,kkp1), &
                     velx3d(iip1,j,k), velx3d(iip1,jjp1,k), velx3d(iip1,j,kkp1), velx3d(iip1,jjp1,kkp1), &
                     vely3d(iip1,j,k), vely3d(iip1,jjp1,k), vely3d(iip1,j,kkp1), vely3d(iip1,jjp1,kkp1), &
                     velz3d(iip1,j,k), velz3d(iip1,jjp1,k), velz3d(iip1,j,kkp1), velz3d(iip1,jjp1,kkp1), &
                     vth3d(iip1,j,k), vth3d(iip1,jjp1,k), vth3d(iip1,j,kkp1), vth3d(iip1,jjp1,kkp1), &
                     sline3d(iip1,j,k), sline3d(iip1,jjp1,k), sline3d(iip1,j,kkp1), sline3d(iip1,jjp1,kkp1), &
                     y(j), y(jjp1), z(k), z(kkp1), y_d, z_d, &
                     c15_slined, c18_slined, c24_slined, c27_slined, &
                     opalbar_d, velx_d, vely_d, velz_d, vth_d, sline_d)
!set interpolation coefficients that are not used to zero
                  c17_slined = 0.d0
                  c23_slined = 0.d0
                  c26_slined = 0.d0
!
               else
                  write(*,'(4es20.8,l4)') dels_d, dels_xzd, dels_xyd, dels_yzd, dels_d.eq.dels_xzd
                  stop 'error in fsc_line3d: invalid dels_d'
               endif
!
!--------------------------------radiative transfer---------------------
!
               vel_u = velx_u*nn_x + vely_u*nn_y + velz_u*nn_z
               vel_d = velx_d*nn_x + vely_d*nn_y + velz_d*nn_z
               call fsc_line(xobs, vth_fiducial, deltax, xcmf_min, xcmf_max, int_u, opalbar_u, opalbar_p, opalbar_d, &
                  sline_u, sline_p, sline_d, vel_u, vel_p, vel_d, vth_u, vth_p, vth_d, &
                  dels_u, dels_d, abs_sc, contr_sc, int_sc, alo_u, alo_p, alo_d)
!                  call fsc_line_lin(xobs, vth_fiducial, deltax, xcmf_min, xcmf_max, int_u, opalbar_u, opalbar_p, sline_u, sline_p, &
!                                      vel_u, vel_p, vth_u, vth_p, dels_u, abs_sc, contr_sc, int_sc, alo_u, alo_p)
!                 alo_d=0.d0
               int3d(i,j,k) = abs_sc*int_u + alo_u*sline_u + alo_p*sline_p + alo_d*sline_d
!
               aloline_on_nn3d(iip1,jjp1,kkp1,q1) = imask_totreg3d(iip1,jjp1,kkp1)*&
                  alo_d*c27_slined
               aloline_on_nn3d(i,jjp1,kkp1,q2) = imask_totreg3d(i,jjp1,kkp1) * &
                  (alo_d*c26_slined + abs_sc*c26_intu*aloline_on_nn3d(i,jjp1,kkp1,q1))
               aloline_on_nn3d(iim1,jjp1,kkp1,q3) = imask_totreg3d(iim1,jjp1,kkp1) * &
                  (abs_sc*c26_intu*aloline_on_nn3d(iim1,jjp1,kkp1,q2))
               aloline_on_nn3d(iip1,j,kkp1,q4) = imask_totreg3d(iip1,j,kkp1) * &
                  (alo_d*c24_slined + abs_sc*c24_intu*aloline_on_nn3d(iip1,j,kkp1,q1))
               aloline_on_nn3d(i,j,kkp1,q5) = imask_totreg3d(i,j,kkp1) * &
                  (alo_d*c23_slined + abs_sc*(c23_intu*aloline_on_nn3d(i,j,kkp1,q1) + &
                  c24_intu*aloline_on_nn3d(i,j,kkp1,q2) + &
                  c26_intu*aloline_on_nn3d(i,j,kkp1,q4)))
               aloline_on_nn3d(iim1,j,kkp1,q6) = imask_totreg3d(iim1,j,kkp1) * &
                  (abs_sc*(c23_intu*aloline_on_nn3d(iim1,j,kkp1,q2) + &
                  c24_intu*aloline_on_nn3d(iim1,j,kkp1,q3) + &
                  c26_intu*aloline_on_nn3d(iim1,j,kkp1,q5)))
               aloline_on_nn3d(iip1,jjm1,kkp1,q7) = imask_totreg3d(iip1,jjm1,kkp1) * &
                  (abs_sc*c24_intu*aloline_on_nn3d(iip1,jjm1,kkp1,q4))
               aloline_on_nn3d(i,jjm1,kkp1,q8) = imask_totreg3d(i,jjm1,kkp1) * &
                  (abs_sc*(c23_intu*aloline_on_nn3d(i,jjm1,kkp1,q4) + &
                  c24_intu*aloline_on_nn3d(i,jjm1,kkp1,q5) + &
                  c26_intu*aloline_on_nn3d(i,jjm1,kkp1,q7)))
               aloline_on_nn3d(iim1,jjm1,kkp1,q9) = imask_totreg3d(iim1,jjm1,kkp1) * &
                  (abs_sc*(c23_intu*aloline_on_nn3d(iim1,jjm1,kkp1,q5) + &
                  c24_intu*aloline_on_nn3d(iim1,jjm1,kkp1,q6) + &
                  c26_intu*aloline_on_nn3d(iim1,jjm1,kkp1,q8)))
               aloline_on_nn3d(iip1,jjp1,k,q10) = imask_totreg3d(iip1,jjp1,k) * &
                  (alo_d*c18_slined + abs_sc*c18_intu*aloline_on_nn3d(iip1,jjp1,k,q1))
               aloline_on_nn3d(i,jjp1,k,q11) = imask_totreg3d(i,jjp1,k) * &
                  (alo_d*c17_slined + abs_sc*(c17_intu*aloline_on_nn3d(i,jjp1,k,q1) + &
                  c18_intu*aloline_on_nn3d(i,jjp1,k,q2) + &
                  c26_intu*aloline_on_nn3d(i,jjp1,k,q10)))
               aloline_on_nn3d(iim1,jjp1,k,q12) = imask_totreg3d(iim1,jjp1,k) * &
                  (abs_sc*(c17_intu*aloline_on_nn3d(iim1,jjp1,k,q2) + &
                  c18_intu*aloline_on_nn3d(iim1,jjp1,k,q3) + &
                  c26_intu*aloline_on_nn3d(iim1,jjp1,k,q11)))
               aloline_on_nn3d(iip1,j,k,q13) = imask_totreg3d(iip1,j,k) * &
                  (alo_d*c15_slined + abs_sc*(c15_intu*aloline_on_nn3d(iip1,j,k,q1) + &
                  c18_intu*aloline_on_nn3d(iip1,j,k,q4) + &
                  c24_intu*aloline_on_nn3d(iip1,j,k,q10)))
               aloline_on_nn3d(i,j,k,q14) = imask_totreg3d(i,j,k) * &
                  (alo_p + abs_sc*(c14_intu*aloline_on_nn3d(i,j,k,q1) + &
                  c15_intu*aloline_on_nn3d(i,j,k,q2) + &
                  c17_intu*aloline_on_nn3d(i,j,k,q4) + &
                  c26_intu*aloline_on_nn3d(i,j,k,q13) + &
                  c18_intu*aloline_on_nn3d(i,j,k,q5) + &
                  c23_intu*aloline_on_nn3d(i,j,k,q10) + &
                  c24_intu*aloline_on_nn3d(i,j,k,q11)))
               aloline_on_nn3d(iim1,j,k,q15) = imask_totreg3d(iim1,j,k) * &
                  (alo_u*c26_slineu + abs_sc*(c14_intu*aloline_on_nn3d(iim1,j,k,q2) + &
                  c15_intu*aloline_on_nn3d(iim1,j,k,q3) + &
                  c17_intu*aloline_on_nn3d(iim1,j,k,q5) + &
                  c18_intu*aloline_on_nn3d(iim1,j,k,q6) + &
                  c23_intu*aloline_on_nn3d(iim1,j,k,q11) + &
                  c24_intu*aloline_on_nn3d(iim1,j,k,q12) + &
                  c26_intu*aloline_on_nn3d(iim1,j,k,q14)))
               aloline_on_nn3d(iip1,jjm1,k,q16) = imask_totreg3d(iip1,jjm1,k) * &
                  (abs_sc*(c15_intu*aloline_on_nn3d(iip1,jjm1,k,q4) + &
                  c18_intu*aloline_on_nn3d(iip1,jjm1,k,q7) + &
                  c24_intu*aloline_on_nn3d(iip1,jjm1,k,q13)))
               aloline_on_nn3d(i,jjm1,k,q17) = imask_totreg3d(i,jjm1,k) * &
                  (alo_u*c24_slineu + abs_sc*(c14_intu*aloline_on_nn3d(i,jjm1,k,q4) + &
                  c15_intu*aloline_on_nn3d(i,jjm1,k,q5) + &
                  c17_intu*aloline_on_nn3d(i,jjm1,k,q7) + &
                  c18_intu*aloline_on_nn3d(i,jjm1,k,q8) + &
                  c23_intu*aloline_on_nn3d(i,jjm1,k,q13) + &
                  c24_intu*aloline_on_nn3d(i,jjm1,k,q14) + &
                  c26_intu*aloline_on_nn3d(i,jjm1,k,q16)))
               aloline_on_nn3d(iim1,jjm1,k,q18) = imask_totreg3d(iim1,jjm1,k) * &
                  (alo_u*c23_slineu + abs_sc*(c14_intu*aloline_on_nn3d(iim1,jjm1,k,q5) + &
                  c15_intu*aloline_on_nn3d(iim1,jjm1,k,q6) + &
                  c17_intu*aloline_on_nn3d(iim1,jjm1,k,q8) + &
                  c18_intu*aloline_on_nn3d(iim1,jjm1,k,q9) + &
                  c23_intu*aloline_on_nn3d(iim1,jjm1,k,q14) + &
                  c24_intu*aloline_on_nn3d(iim1,jjm1,k,q15) + &
                  c26_intu*aloline_on_nn3d(iim1,jjm1,k,q17)))
               aloline_on_nn3d(iip1,jjp1,kkm1,q19) = imask_totreg3d(iip1,jjp1,kkm1) * &
                  (abs_sc*c18_intu*aloline_on_nn3d(iip1,jjp1,kkm1,q10))
               aloline_on_nn3d(i,jjp1,kkm1,q20) = imask_totreg3d(i,jjp1,kkm1) * &
                  (abs_sc*(c17_intu*aloline_on_nn3d(i,jjp1,kkm1,q10) + &
                  c18_intu*aloline_on_nn3d(i,jjp1,kkm1,q11) + &
                  c26_intu*aloline_on_nn3d(i,jjp1,kkm1,q19)))
               aloline_on_nn3d(iim1,jjp1,kkm1,q21) = imask_totreg3d(iim1,jjp1,kkm1) * &
                  (abs_sc*(c17_intu*aloline_on_nn3d(iim1,jjp1,kkm1,q11) + &
                  c18_intu*aloline_on_nn3d(iim1,jjp1,kkm1,q12) + &
                  c26_intu*aloline_on_nn3d(iim1,jjp1,kkm1,q20)))
               aloline_on_nn3d(iip1,j,kkm1,q22) = imask_totreg3d(iip1,j,kkm1) * &
                  (abs_sc*(c15_intu*aloline_on_nn3d(iip1,j,kkm1,q10) + &
                  c18_intu*aloline_on_nn3d(iip1,j,kkm1,q13) + &
                  c24_intu*aloline_on_nn3d(iip1,j,kkm1,q19)))
               aloline_on_nn3d(i,j,kkm1,q23) = imask_totreg3d(i,j,kkm1) * &
                  (alo_u*c18_slineu + abs_sc*(c14_intu*aloline_on_nn3d(i,j,kkm1,q10) + &
                  c15_intu*aloline_on_nn3d(i,j,kkm1,q11) + &
                  c17_intu*aloline_on_nn3d(i,j,kkm1,q13) + &
                  c18_intu*aloline_on_nn3d(i,j,kkm1,q14) + &
                  c23_intu*aloline_on_nn3d(i,j,kkm1,q19) + &
                  c24_intu*aloline_on_nn3d(i,j,kkm1,q20) + &
                  c26_intu*aloline_on_nn3d(i,j,kkm1,q22)))
               aloline_on_nn3d(iim1,j,kkm1,q24) = imask_totreg3d(iim1,j,kkm1) * &
                  (alo_u*c17_slineu + abs_sc*(c14_intu*aloline_on_nn3d(iim1,j,kkm1,q11) + &
                  c15_intu*aloline_on_nn3d(iim1,j,kkm1,q12) + &
                  c17_intu*aloline_on_nn3d(iim1,j,kkm1,q14) + &
                  c18_intu*aloline_on_nn3d(iim1,j,kkm1,q15) + &
                  c23_intu*aloline_on_nn3d(iim1,j,kkm1,q20) + &
                  c24_intu*aloline_on_nn3d(iim1,j,kkm1,q21) + &
                  c26_intu*aloline_on_nn3d(iim1,j,kkm1,q23)))
               aloline_on_nn3d(iip1,jjm1,kkm1,q25) = imask_totreg3d(iip1,jjm1,kkm1) * &
                  (abs_sc*(c15_intu*aloline_on_nn3d(iip1,jjm1,kkm1,q13) + &
                  c18_intu*aloline_on_nn3d(iip1,jjm1,kkm1,q16) + &
                  c24_intu*aloline_on_nn3d(iip1,jjm1,kkm1,q22)))
               aloline_on_nn3d(i,jjm1,kkm1,q26) = imask_totreg3d(i,jjm1,kkm1) * &
                  (alo_u*c15_slineu + abs_sc*(c14_intu*aloline_on_nn3d(i,jjm1,kkm1,q13) + &
                  c15_intu*aloline_on_nn3d(i,jjm1,kkm1,q14) + &
                  c17_intu*aloline_on_nn3d(i,jjm1,kkm1,q16) + &
                  c18_intu*aloline_on_nn3d(i,jjm1,kkm1,q17) + &
                  c23_intu*aloline_on_nn3d(i,jjm1,kkm1,q22) + &
                  c24_intu*aloline_on_nn3d(i,jjm1,kkm1,q23) + &
                  c26_intu*aloline_on_nn3d(i,jjm1,kkm1,q25)))
               aloline_on_nn3d(iim1,jjm1,kkm1,q27) = imask_totreg3d(iim1,jjm1,kkm1) * &
                  (alo_u*c14_slineu + abs_sc*(c14_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q14) + &
                  c15_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q15) + &
                  c17_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q17) + &
                  c18_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q18) + &
                  c23_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q23) + &
                  c24_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q24) + &
                  c26_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q26)))
!
!perform angular integration
               call calc_phinorm(vel_p, vth_p, vth_fiducial, xobs, phinorm)
               wall=wall_global*phinorm
               mintbar3d_tmp(i,j,k) = mintbar3d_tmp(i,j,k) + int3d(i,j,k)*wall
               normalization3d_tmp(i,j,k) = normalization3d_tmp(i,j,k) + wall
!
               aloline_nn3d_tmp(iip1,jjp1,kkp1,q1) = aloline_nn3d_tmp(iip1,jjp1,kkp1,q1) + wall*aloline_on_nn3d(iip1,jjp1,kkp1,q1)
               aloline_nn3d_tmp(i,jjp1,kkp1,q2) = aloline_nn3d_tmp(i,jjp1,kkp1,q2) + wall*aloline_on_nn3d(i,jjp1,kkp1,q2)
               aloline_nn3d_tmp(iim1,jjp1,kkp1,q3) = aloline_nn3d_tmp(iim1,jjp1,kkp1,q3) + wall*aloline_on_nn3d(iim1,jjp1,kkp1,q3)
               aloline_nn3d_tmp(iip1,j,kkp1,q4) = aloline_nn3d_tmp(iip1,j,kkp1,q4) + wall*aloline_on_nn3d(iip1,j,kkp1,q4)
               aloline_nn3d_tmp(i,j,kkp1,q5) = aloline_nn3d_tmp(i,j,kkp1,q5) + wall*aloline_on_nn3d(i,j,kkp1,q5)
               aloline_nn3d_tmp(iim1,j,kkp1,q6) = aloline_nn3d_tmp(iim1,j,kkp1,q6) + wall*aloline_on_nn3d(iim1,j,kkp1,q6)
               aloline_nn3d_tmp(iip1,jjm1,kkp1,q7) = aloline_nn3d_tmp(iip1,jjm1,kkp1,q7) + wall*aloline_on_nn3d(iip1,jjm1,kkp1,q7)
               aloline_nn3d_tmp(i,jjm1,kkp1,q8) = aloline_nn3d_tmp(i,jjm1,kkp1,q8) + wall*aloline_on_nn3d(i,jjm1,kkp1,q8)
               aloline_nn3d_tmp(iim1,jjm1,kkp1,q9) = aloline_nn3d_tmp(iim1,jjm1,kkp1,q9) + wall*aloline_on_nn3d(iim1,jjm1,kkp1,q9)
               aloline_nn3d_tmp(iip1,jjp1,k,q10) = aloline_nn3d_tmp(iip1,jjp1,k,q10) + wall*aloline_on_nn3d(iip1,jjp1,k,q10)
               aloline_nn3d_tmp(i,jjp1,k,q11) = aloline_nn3d_tmp(i,jjp1,k,q11) + wall*aloline_on_nn3d(i,jjp1,k,q11)
               aloline_nn3d_tmp(iim1,jjp1,k,q12) = aloline_nn3d_tmp(iim1,jjp1,k,q12) + wall*aloline_on_nn3d(iim1,jjp1,k,q12)
               aloline_nn3d_tmp(iip1,j,k,q13) = aloline_nn3d_tmp(iip1,j,k,q13) + wall*aloline_on_nn3d(iip1,j,k,q13)
               aloline_nn3d_tmp(i,j,k,q14) = aloline_nn3d_tmp(i,j,k,q14) + wall*aloline_on_nn3d(i,j,k,q14)
               aloline_nn3d_tmp(iim1,j,k,q15) = aloline_nn3d_tmp(iim1,j,k,q15) + wall*aloline_on_nn3d(iim1,j,k,q15)
               aloline_nn3d_tmp(iip1,jjm1,k,q16) = aloline_nn3d_tmp(iip1,jjm1,k,q16) + wall*aloline_on_nn3d(iip1,jjm1,k,q16) !new
               aloline_nn3d_tmp(i,jjm1,k,q17) = aloline_nn3d_tmp(i,jjm1,k,q17) + wall*aloline_on_nn3d(i,jjm1,k,q17)
               aloline_nn3d_tmp(iim1,jjm1,k,q18) = aloline_nn3d_tmp(iim1,jjm1,k,q18) + wall*aloline_on_nn3d(iim1,jjm1,k,q18)
               aloline_nn3d_tmp(iip1,jjp1,kkm1,q19) = aloline_nn3d_tmp(iip1,jjp1,kkm1,q19) + wall*aloline_on_nn3d(iip1,jjp1,kkm1,q19) !new
               aloline_nn3d_tmp(i,jjp1,kkm1,q20) = aloline_nn3d_tmp(i,jjp1,kkm1,q20) + wall*aloline_on_nn3d(i,jjp1,kkm1,q20) !new
               aloline_nn3d_tmp(iim1,jjp1,kkm1,q21) = aloline_nn3d_tmp(iim1,jjp1,kkm1,q21) + wall*aloline_on_nn3d(iim1,jjp1,kkm1,q21) !new
               aloline_nn3d_tmp(iip1,j,kkm1,q22) = aloline_nn3d_tmp(iip1,j,kkm1,q22) + wall*aloline_on_nn3d(iip1,j,kkm1,q22) !new
               aloline_nn3d_tmp(i,j,kkm1,q23) = aloline_nn3d_tmp(i,j,kkm1,q23) + wall*aloline_on_nn3d(i,j,kkm1,q23)
               aloline_nn3d_tmp(iim1,j,kkm1,q24) = aloline_nn3d_tmp(iim1,j,kkm1,q24) + wall*aloline_on_nn3d(iim1,j,kkm1,q24)
               aloline_nn3d_tmp(iip1,jjm1,kkm1,q25) = aloline_nn3d_tmp(iip1,jjm1,kkm1,q25) + wall*aloline_on_nn3d(iip1,jjm1,kkm1,q25) !new
               aloline_nn3d_tmp(i,jjm1,kkm1,q26) = aloline_nn3d_tmp(i,jjm1,kkm1,q26) + wall*aloline_on_nn3d(i,jjm1,kkm1,q26)
               aloline_nn3d_tmp(iim1,jjm1,kkm1,q27) = aloline_nn3d_tmp(iim1,jjm1,kkm1,q27) + wall*aloline_on_nn3d(iim1,jjm1,kkm1,q27)
!
!
!***********special treatment: point lies directly on boundary**********
!
             case(3)
               mueff=nn_x*x(i)+nn_y*y(j)+nn_z*z(k)
               if(mueff.ge.0.d0) then
!no interpolations required
!for alo
!                 int3d(i,j,k) = 0.d0
!                  int3d(i,j,k) = xic1
                  int3d(i,j,k)=calc_icore_gdark(z(k), sqrt(x(i)**2+y(j)**2+z(k)**2))
                  aloline_on_nn3d(i,j,k,14) = 0.d0
!
!perform angular integration (alo not required, since boundary specified)
                  vel_p = nn_x*velx3d(i,j,k) + nn_y*vely3d(i,j,k) + nn_z*velz3d(i,j,k)
                  vth_p = vth3d(i,j,k)
                  call calc_phinorm(vel_p, vth_p, vth_fiducial, xobs, phinorm)
                  wall=wall_global*phinorm
                  mintbar3d_tmp(i,j,k) = mintbar3d_tmp(i,j,k) + int3d(i,j,k)*wall
                  normalization3d_tmp(i,j,k) = normalization3d_tmp(i,j,k) + wall
               else
!same interpolation as in (1)
!
                  iim1=i-alpha
                  jjm1=j-beta
                  kkm1=k-gamma
                  iip1=i+alpha
                  jjp1=j+beta
                  kkp1=k+gamma
!
!calculate distance of upwind point to (previous) xy-plane, xz-plane and yz-plane
                  dels_xyu=(z(k)-z(kkm1))/nn_z
                  dels_xzu=(y(j)-y(jjm1))/nn_y
                  dels_yzu=(x(i)-x(iim1))/nn_x
                  dels_u=min(dels_xyu,dels_xzu,dels_yzu)
!
!calculate distance of downwind point to (next) xy-plane, xz-plane and yz-plane
                  dels_xyd=(z(kkp1)-z(k))/nn_z
                  dels_xzd=(y(jjp1)-y(j))/nn_y
                  dels_yzd=(x(iip1)-x(i))/nn_x
                  dels_d=min(dels_xyd,dels_xzd,dels_yzd)
!
!----------------------------local point--------------------------------
!
                  sline_p=sline3d(i,j,k)
                  opalbar_p=opalbar3d(i,j,k)
                  velx_p=velx3d(i,j,k)
                  vely_p=vely3d(i,j,k)
                  velz_p=velz3d(i,j,k)
                  vel_p = nn_x*velx_p + nn_y*vely_p + nn_z*velz_p
                  vth_p = vth3d(i,j,k)
!
!----------------------------upwind point-------------------------------
!
                  if(dels_xyu.eq.dels_u) then
!intersection with x-y plane on level k-gamma
                     x_u = x(i) - dels_u*nn_x
                     y_u = y(j) - dels_u*nn_y
                     z_u = z(kkm1)
!
                     call coeff2d_lineu_lin(opalbar3d(iim1,jjm1,kkm1), opalbar3d(i,jjm1,kkm1), opalbar3d(iim1,j,kkm1), opalbar3d(i,j,kkm1), &
                        velx3d(iim1,jjm1,kkm1), velx3d(i,jjm1,kkm1), velx3d(iim1,j,kkm1), velx3d(i,j,kkm1), &
                        vely3d(iim1,jjm1,kkm1), vely3d(i,jjm1,kkm1), vely3d(iim1,j,kkm1), vely3d(i,j,kkm1), &
                        velz3d(iim1,jjm1,kkm1), velz3d(i,jjm1,kkm1), velz3d(iim1,j,kkm1), velz3d(i,j,kkm1), &
                        vth3d(iim1,jjm1,kkm1), vth3d(i,jjm1,kkm1), vth3d(iim1,j,kkm1), vth3d(i,j,kkm1), &
                        sline3d(iim1,jjm1,kkm1), sline3d(i,jjm1,kkm1), sline3d(iim1,j,kkm1), sline3d(i,j,kkm1), &
                        int3d(iim1,jjm1,kkm1), int3d(i,jjm1,kkm1), int3d(iim1,j,kkm1), int3d(i,j,kkm1), &
                        x(iim1), x(i), y(jjm1), y(j), x_u, y_u, &
                        c14_slineu, c15_slineu, c17_slineu, c18_slineu, &
                        c14_intu, c15_intu, c17_intu, c18_intu, opalbar_u, velx_u, vely_u, velz_u, vth_u, sline_u, int_u)
                     c23_slineu=0.d0
                     c24_slineu=0.d0
                     c26_slineu=0.d0
                     c23_intu=0.d0
                     c24_intu=0.d0
                     c26_intu=0.d0
!
                  elseif(dels_xzu.eq.dels_u) then
!intersection with x-z plane on level j-beta
                     x_u = x(i) - dels_u*nn_x
                     y_u = y(jjm1)
                     z_u = z(k) - dels_u*nn_z
!
                     call coeff2d_lineu_lin(opalbar3d(iim1,jjm1,kkm1), opalbar3d(i,jjm1,kkm1), opalbar3d(iim1,jjm1,k), opalbar3d(i,jjm1,k), &
                        velx3d(iim1,jjm1,kkm1), velx3d(i,jjm1,kkm1), velx3d(iim1,jjm1,k), velx3d(i,jjm1,k), &
                        vely3d(iim1,jjm1,kkm1), vely3d(i,jjm1,kkm1), vely3d(iim1,jjm1,k), vely3d(i,jjm1,k), &
                        velz3d(iim1,jjm1,kkm1), velz3d(i,jjm1,kkm1), velz3d(iim1,jjm1,k), velz3d(i,jjm1,k), &
                        vth3d(iim1,jjm1,kkm1), vth3d(i,jjm1,kkm1), vth3d(iim1,jjm1,k), vth3d(i,jjm1,k), &
                        sline3d(iim1,jjm1,kkm1), sline3d(i,jjm1,kkm1), sline3d(iim1,jjm1,k), sline3d(i,jjm1,k), &
                        int3d(iim1,jjm1,kkm1), int3d(i,jjm1,kkm1), int3d(iim1,jjm1,k), int3d(i,jjm1,k), &
                        x(iim1), x(i), z(kkm1), z(k), x_u, z_u, &
                        c14_slineu, c15_slineu, c23_slineu, c24_slineu, &
                        c14_intu, c15_intu, c23_intu, c24_intu, opalbar_u, velx_u, vely_u, velz_u, vth_u, sline_u, int_u)
                     c17_slineu=0.d0
                     c18_slineu=0.d0
                     c26_slineu=0.d0
                     c17_intu=0.d0
                     c18_intu=0.d0
                     c26_intu=0.d0
!
                  elseif(dels_yzu.eq.dels_u) then
!intersection with y-z plane on level i-alpha
                     x_u = x(iim1)
                     y_u = y(j) - dels_u*nn_y
                     z_u = z(k) - dels_u*nn_z
!
                     call coeff2d_lineu_lin(opalbar3d(iim1,jjm1,kkm1), opalbar3d(iim1,j,kkm1), opalbar3d(iim1,jjm1,k), opalbar3d(iim1,j,k), &
                        velx3d(iim1,jjm1,kkm1), velx3d(iim1,j,kkm1), velx3d(iim1,jjm1,k), velx3d(iim1,j,k), &
                        vely3d(iim1,jjm1,kkm1), vely3d(iim1,j,kkm1), vely3d(iim1,jjm1,k), vely3d(iim1,j,k), &
                        velz3d(iim1,jjm1,kkm1), velz3d(iim1,j,kkm1), velz3d(iim1,jjm1,k), velz3d(iim1,j,k), &
                        vth3d(iim1,jjm1,kkm1), vth3d(iim1,j,kkm1), vth3d(iim1,jjm1,k), vth3d(iim1,j,k), &
                        sline3d(iim1,jjm1,kkm1), sline3d(iim1,j,kkm1), sline3d(iim1,jjm1,k), sline3d(iim1,j,k), &
                        int3d(iim1,jjm1,kkm1), int3d(iim1,j,kkm1), int3d(iim1,jjm1,k), int3d(iim1,j,k), &
                        y(jjm1), y(j), z(kkm1), z(k), y_u, z_u, &
                        c14_slineu, c17_slineu, c23_slineu, c26_slineu, &
                        c14_intu, c17_intu, c23_intu, c26_intu, opalbar_u, velx_u, vely_u, velz_u, vth_u, sline_u, int_u)
                     c15_slineu=0.d0
                     c18_slineu=0.d0
                     c24_slineu=0.d0
                     c15_intu=0.d0
                     c18_intu=0.d0
                     c24_intu=0.d0
!
                  else
                     write(*,'(4es20.8,2l4)') dels_u, dels_xzu, dels_xyu, dels_yzu, dels_u.eq.dels_xzu, dels_xyu.eq.dels_u
                     stop 'error in fsc_line3d: invalid dels_u'
                  endif
!
!---------------------------downwind point------------------------------
!
                  if(dels_xyd.eq.dels_d) then
!intersection with x-y plane on level k+gamma
                     x_d = x(i) + dels_d*nn_x
                     y_d = y(j) + dels_d*nn_y
                     z_d = z(kkp1)
!
                     call coeff2d_lined_lin(opalbar3d(i,j,kkp1), opalbar3d(iip1,j,kkp1), opalbar3d(i,jjp1,kkp1), opalbar3d(iip1,jjp1,kkp1), &
                        velx3d(i,j,kkp1), velx3d(iip1,j,kkp1), velx3d(i,jjp1,kkp1), velx3d(iip1,jjp1,kkp1), &
                        vely3d(i,j,kkp1), vely3d(iip1,j,kkp1), vely3d(i,jjp1,kkp1), vely3d(iip1,jjp1,kkp1), &
                        velz3d(i,j,kkp1), velz3d(iip1,j,kkp1), velz3d(i,jjp1,kkp1), velz3d(iip1,jjp1,kkp1), &
                        vth3d(i,j,kkp1), vth3d(iip1,j,kkp1), vth3d(i,jjp1,kkp1), vth3d(iip1,jjp1,kkp1), &
                        sline3d(i,j,kkp1), sline3d(iip1,j,kkp1), sline3d(i,jjp1,kkp1), sline3d(iip1,jjp1,kkp1), &
                        x(i), x(iip1), y(j), y(jjp1), x_d, y_d, &
                        c23_slined, c24_slined, c26_slined, c27_slined, &
                        opalbar_d, velx_d, vely_d, velz_d, vth_d, sline_d)
!set interpolation coefficients that are not used to zero
                     c15_slined = 0.d0
                     c17_slined = 0.d0
                     c18_slined = 0.d0
!
                  elseif(dels_xzd.eq.dels_d) then
!intersection with x-z plane on level j+beta
                     x_d = x(i) + dels_d*nn_x
                     y_d = y(jjp1)
                     z_d = z(k) + dels_d*nn_z
!
                     call coeff2d_lined_lin(opalbar3d(i,jjp1,k), opalbar3d(iip1,jjp1,k), opalbar3d(i,jjp1,kkp1), opalbar3d(iip1,jjp1,kkp1), &
                        velx3d(i,jjp1,k), velx3d(iip1,jjp1,k), velx3d(i,jjp1,kkp1), velx3d(iip1,jjp1,kkp1), &
                        vely3d(i,jjp1,k), vely3d(iip1,jjp1,k), vely3d(i,jjp1,kkp1), vely3d(iip1,jjp1,kkp1), &
                        velz3d(i,jjp1,k), velz3d(iip1,jjp1,k), velz3d(i,jjp1,kkp1), velz3d(iip1,jjp1,kkp1), &
                        vth3d(i,jjp1,k), vth3d(iip1,jjp1,k), vth3d(i,jjp1,kkp1), vth3d(iip1,jjp1,kkp1), &
                        sline3d(i,jjp1,k), sline3d(iip1,jjp1,k), sline3d(i,jjp1,kkp1), sline3d(iip1,jjp1,kkp1), &
                        x(i), x(iip1), z(k), z(kkp1), x_d, z_d, &
                        c17_slined, c18_slined, c26_slined, c27_slined, &
                        opalbar_d, velx_d, vely_d, velz_d, vth_d, sline_d)
!set interpolation coefficients that are not used to zero
                     c15_slined = 0.d0
                     c23_slined = 0.d0
                     c24_slined = 0.d0
!
                  elseif(dels_yzd.eq.dels_d) then
!intersection with y-z plane on level i+alpha
                     x_d = x(iip1)
                     y_d = y(j) + dels_d*nn_y
                     z_d = z(k) + dels_d*nn_z
!
                     call coeff2d_lined_lin(opalbar3d(iip1,j,k), opalbar3d(iip1,jjp1,k), opalbar3d(iip1,j,kkp1), opalbar3d(iip1,jjp1,kkp1), &
                        velx3d(iip1,j,k), velx3d(iip1,jjp1,k), velx3d(iip1,j,kkp1), velx3d(iip1,jjp1,kkp1), &
                        vely3d(iip1,j,k), vely3d(iip1,jjp1,k), vely3d(iip1,j,kkp1), vely3d(iip1,jjp1,kkp1), &
                        velz3d(iip1,j,k), velz3d(iip1,jjp1,k), velz3d(iip1,j,kkp1), velz3d(iip1,jjp1,kkp1), &
                        vth3d(iip1,j,k), vth3d(iip1,jjp1,k), vth3d(iip1,j,kkp1), vth3d(iip1,jjp1,kkp1), &
                        sline3d(iip1,j,k), sline3d(iip1,jjp1,k), sline3d(iip1,j,kkp1), sline3d(iip1,jjp1,kkp1), &
                        y(j), y(jjp1), z(k), z(kkp1), y_d, z_d, &
                        c15_slined, c18_slined, c24_slined, c27_slined, &
                        opalbar_d, velx_d, vely_d, velz_d, vth_d, sline_d)
!set interpolation coefficients that are not used to zero
                     c17_slined = 0.d0
                     c23_slined = 0.d0
                     c26_slined = 0.d0
!
                  else
                     write(*,'(4es20.8,l4)') dels_d, dels_xzd, dels_xyd, dels_yzd, dels_d.eq.dels_xzd
                     stop 'error in fsc_line3d: invalid dels_d'
                  endif
!
!--------------------------------radiative transfer---------------------
!
                  vel_u = velx_u*nn_x + vely_u*nn_y + velz_u*nn_z
                  vel_d = velx_d*nn_x + vely_d*nn_y + velz_d*nn_z
                  call fsc_line(xobs, vth_fiducial, deltax, xcmf_min, xcmf_max, int_u, opalbar_u, opalbar_p, opalbar_d, &
                     sline_u, sline_p, sline_d, vel_u, vel_p, vel_d, vth_u, vth_p, vth_d, &
                     dels_u, dels_d, abs_sc, contr_sc, int_sc, alo_u, alo_p, alo_d)
!                  call fsc_line_lin(xobs, vth_fiducial, deltax, xcmf_min, xcmf_max, int_u, opalbar_u, opalbar_p, sline_u, sline_p, &
!                                      vel_u, vel_p, vth_u, vth_p, dels_u, abs_sc, contr_sc, int_sc, alo_u, alo_p)
!                 alo_d=0.d0
                  int3d(i,j,k) = abs_sc*int_u + alo_u*sline_u + alo_p*sline_p + alo_d*sline_d
!
                  aloline_on_nn3d(iip1,jjp1,kkp1,q1) = imask_totreg3d(iip1,jjp1,kkp1)*&
                     alo_d*c27_slined
                  aloline_on_nn3d(i,jjp1,kkp1,q2) = imask_totreg3d(i,jjp1,kkp1) * &
                     (alo_d*c26_slined + abs_sc*c26_intu*aloline_on_nn3d(i,jjp1,kkp1,q1))
                  aloline_on_nn3d(iim1,jjp1,kkp1,q3) = imask_totreg3d(iim1,jjp1,kkp1) * &
                     (abs_sc*c26_intu*aloline_on_nn3d(iim1,jjp1,kkp1,q2))
                  aloline_on_nn3d(iip1,j,kkp1,q4) = imask_totreg3d(iip1,j,kkp1) * &
                     (alo_d*c24_slined + abs_sc*c24_intu*aloline_on_nn3d(iip1,j,kkp1,q1))
                  aloline_on_nn3d(i,j,kkp1,q5) = imask_totreg3d(i,j,kkp1) * &
                     (alo_d*c23_slined + abs_sc*(c23_intu*aloline_on_nn3d(i,j,kkp1,q1) + &
                     c24_intu*aloline_on_nn3d(i,j,kkp1,q2) + &
                     c26_intu*aloline_on_nn3d(i,j,kkp1,q4)))
                  aloline_on_nn3d(iim1,j,kkp1,q6) = imask_totreg3d(iim1,j,kkp1) * &
                     (abs_sc*(c23_intu*aloline_on_nn3d(iim1,j,kkp1,q2) + &
                     c24_intu*aloline_on_nn3d(iim1,j,kkp1,q3) + &
                     c26_intu*aloline_on_nn3d(iim1,j,kkp1,q5)))
                  aloline_on_nn3d(iip1,jjm1,kkp1,q7) = imask_totreg3d(iip1,jjm1,kkp1) * &
                     (abs_sc*c24_intu*aloline_on_nn3d(iip1,jjm1,kkp1,q4))
                  aloline_on_nn3d(i,jjm1,kkp1,q8) = imask_totreg3d(i,jjm1,kkp1) * &
                     (abs_sc*(c23_intu*aloline_on_nn3d(i,jjm1,kkp1,q4) + &
                     c24_intu*aloline_on_nn3d(i,jjm1,kkp1,q5) + &
                     c26_intu*aloline_on_nn3d(i,jjm1,kkp1,q7)))
                  aloline_on_nn3d(iim1,jjm1,kkp1,q9) = imask_totreg3d(iim1,jjm1,kkp1) * &
                     (abs_sc*(c23_intu*aloline_on_nn3d(iim1,jjm1,kkp1,q5) + &
                     c24_intu*aloline_on_nn3d(iim1,jjm1,kkp1,q6) + &
                     c26_intu*aloline_on_nn3d(iim1,jjm1,kkp1,q8)))
                  aloline_on_nn3d(iip1,jjp1,k,q10) = imask_totreg3d(iip1,jjp1,k) * &
                     (alo_d*c18_slined + abs_sc*c18_intu*aloline_on_nn3d(iip1,jjp1,k,q1))
                  aloline_on_nn3d(i,jjp1,k,q11) = imask_totreg3d(i,jjp1,k) * &
                     (alo_d*c17_slined + abs_sc*(c17_intu*aloline_on_nn3d(i,jjp1,k,q1) + &
                     c18_intu*aloline_on_nn3d(i,jjp1,k,q2) + &
                     c26_intu*aloline_on_nn3d(i,jjp1,k,q10)))
                  aloline_on_nn3d(iim1,jjp1,k,q12) = imask_totreg3d(iim1,jjp1,k) * &
                     (abs_sc*(c17_intu*aloline_on_nn3d(iim1,jjp1,k,q2) + &
                     c18_intu*aloline_on_nn3d(iim1,jjp1,k,q3) + &
                     c26_intu*aloline_on_nn3d(iim1,jjp1,k,q11)))
                  aloline_on_nn3d(iip1,j,k,q13) = imask_totreg3d(iip1,j,k) * &
                     (alo_d*c15_slined + abs_sc*(c15_intu*aloline_on_nn3d(iip1,j,k,q1) + &
                     c18_intu*aloline_on_nn3d(iip1,j,k,q4) + &
                     c24_intu*aloline_on_nn3d(iip1,j,k,q10)))
                  aloline_on_nn3d(i,j,k,q14) = imask_totreg3d(i,j,k) * &
                     (alo_p + abs_sc*(c14_intu*aloline_on_nn3d(i,j,k,q1) + &
                     c15_intu*aloline_on_nn3d(i,j,k,q2) + &
                     c17_intu*aloline_on_nn3d(i,j,k,q4) + &
                     c26_intu*aloline_on_nn3d(i,j,k,q13) + &
                     c18_intu*aloline_on_nn3d(i,j,k,q5) + &
                     c23_intu*aloline_on_nn3d(i,j,k,q10) + &
                     c24_intu*aloline_on_nn3d(i,j,k,q11)))
                  aloline_on_nn3d(iim1,j,k,q15) = imask_totreg3d(iim1,j,k) * &
                     (alo_u*c26_slineu + abs_sc*(c14_intu*aloline_on_nn3d(iim1,j,k,q2) + &
                     c15_intu*aloline_on_nn3d(iim1,j,k,q3) + &
                     c17_intu*aloline_on_nn3d(iim1,j,k,q5) + &
                     c18_intu*aloline_on_nn3d(iim1,j,k,q6) + &
                     c23_intu*aloline_on_nn3d(iim1,j,k,q11) + &
                     c24_intu*aloline_on_nn3d(iim1,j,k,q12) + &
                     c26_intu*aloline_on_nn3d(iim1,j,k,q14)))
                  aloline_on_nn3d(iip1,jjm1,k,q16) = imask_totreg3d(iip1,jjm1,k) * &
                     (abs_sc*(c15_intu*aloline_on_nn3d(iip1,jjm1,k,q4) + &
                     c18_intu*aloline_on_nn3d(iip1,jjm1,k,q7) + &
                     c24_intu*aloline_on_nn3d(iip1,jjm1,k,q13)))
                  aloline_on_nn3d(i,jjm1,k,q17) = imask_totreg3d(i,jjm1,k) * &
                     (alo_u*c24_slineu + abs_sc*(c14_intu*aloline_on_nn3d(i,jjm1,k,q4) + &
                     c15_intu*aloline_on_nn3d(i,jjm1,k,q5) + &
                     c17_intu*aloline_on_nn3d(i,jjm1,k,q7) + &
                     c18_intu*aloline_on_nn3d(i,jjm1,k,q8) + &
                     c23_intu*aloline_on_nn3d(i,jjm1,k,q13) + &
                     c24_intu*aloline_on_nn3d(i,jjm1,k,q14) + &
                     c26_intu*aloline_on_nn3d(i,jjm1,k,q16)))
                  aloline_on_nn3d(iim1,jjm1,k,q18) = imask_totreg3d(iim1,jjm1,k) * &
                     (alo_u*c23_slineu + abs_sc*(c14_intu*aloline_on_nn3d(iim1,jjm1,k,q5) + &
                     c15_intu*aloline_on_nn3d(iim1,jjm1,k,q6) + &
                     c17_intu*aloline_on_nn3d(iim1,jjm1,k,q8) + &
                     c18_intu*aloline_on_nn3d(iim1,jjm1,k,q9) + &
                     c23_intu*aloline_on_nn3d(iim1,jjm1,k,q14) + &
                     c24_intu*aloline_on_nn3d(iim1,jjm1,k,q15) + &
                     c26_intu*aloline_on_nn3d(iim1,jjm1,k,q17)))
                  aloline_on_nn3d(iip1,jjp1,kkm1,q19) = imask_totreg3d(iip1,jjp1,kkm1) * &
                     (abs_sc*c18_intu*aloline_on_nn3d(iip1,jjp1,kkm1,q10))
                  aloline_on_nn3d(i,jjp1,kkm1,q20) = imask_totreg3d(i,jjp1,kkm1) * &
                     (abs_sc*(c17_intu*aloline_on_nn3d(i,jjp1,kkm1,q10) + &
                     c18_intu*aloline_on_nn3d(i,jjp1,kkm1,q11) + &
                     c26_intu*aloline_on_nn3d(i,jjp1,kkm1,q19)))
                  aloline_on_nn3d(iim1,jjp1,kkm1,q21) = imask_totreg3d(iim1,jjp1,kkm1) * &
                     (abs_sc*(c17_intu*aloline_on_nn3d(iim1,jjp1,kkm1,q11) + &
                     c18_intu*aloline_on_nn3d(iim1,jjp1,kkm1,q12) + &
                     c26_intu*aloline_on_nn3d(iim1,jjp1,kkm1,q20)))
                  aloline_on_nn3d(iip1,j,kkm1,q22) = imask_totreg3d(iip1,j,kkm1) * &
                     (abs_sc*(c15_intu*aloline_on_nn3d(iip1,j,kkm1,q10) + &
                     c18_intu*aloline_on_nn3d(iip1,j,kkm1,q13) + &
                     c24_intu*aloline_on_nn3d(iip1,j,kkm1,q19)))
                  aloline_on_nn3d(i,j,kkm1,q23) = imask_totreg3d(i,j,kkm1) * &
                     (alo_u*c18_slineu + abs_sc*(c14_intu*aloline_on_nn3d(i,j,kkm1,q10) + &
                     c15_intu*aloline_on_nn3d(i,j,kkm1,q11) + &
                     c17_intu*aloline_on_nn3d(i,j,kkm1,q13) + &
                     c18_intu*aloline_on_nn3d(i,j,kkm1,q14) + &
                     c23_intu*aloline_on_nn3d(i,j,kkm1,q19) + &
                     c24_intu*aloline_on_nn3d(i,j,kkm1,q20) + &
                     c26_intu*aloline_on_nn3d(i,j,kkm1,q22)))
                  aloline_on_nn3d(iim1,j,kkm1,q24) = imask_totreg3d(iim1,j,kkm1) * &
                     (alo_u*c17_slineu + abs_sc*(c14_intu*aloline_on_nn3d(iim1,j,kkm1,q11) + &
                     c15_intu*aloline_on_nn3d(iim1,j,kkm1,q12) + &
                     c17_intu*aloline_on_nn3d(iim1,j,kkm1,q14) + &
                     c18_intu*aloline_on_nn3d(iim1,j,kkm1,q15) + &
                     c23_intu*aloline_on_nn3d(iim1,j,kkm1,q20) + &
                     c24_intu*aloline_on_nn3d(iim1,j,kkm1,q21) + &
                     c26_intu*aloline_on_nn3d(iim1,j,kkm1,q23)))
                  aloline_on_nn3d(iip1,jjm1,kkm1,q25) = imask_totreg3d(iip1,jjm1,kkm1) * &
                     (abs_sc*(c15_intu*aloline_on_nn3d(iip1,jjm1,kkm1,q13) + &
                     c18_intu*aloline_on_nn3d(iip1,jjm1,kkm1,q16) + &
                     c24_intu*aloline_on_nn3d(iip1,jjm1,kkm1,q22)))
                  aloline_on_nn3d(i,jjm1,kkm1,q26) = imask_totreg3d(i,jjm1,kkm1) * &
                     (alo_u*c15_slineu + abs_sc*(c14_intu*aloline_on_nn3d(i,jjm1,kkm1,q13) + &
                     c15_intu*aloline_on_nn3d(i,jjm1,kkm1,q14) + &
                     c17_intu*aloline_on_nn3d(i,jjm1,kkm1,q16) + &
                     c18_intu*aloline_on_nn3d(i,jjm1,kkm1,q17) + &
                     c23_intu*aloline_on_nn3d(i,jjm1,kkm1,q22) + &
                     c24_intu*aloline_on_nn3d(i,jjm1,kkm1,q23) + &
                     c26_intu*aloline_on_nn3d(i,jjm1,kkm1,q25)))
                  aloline_on_nn3d(iim1,jjm1,kkm1,q27) = imask_totreg3d(iim1,jjm1,kkm1) * &
                     (alo_u*c14_slineu + abs_sc*(c14_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q14) + &
                     c15_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q15) + &
                     c17_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q17) + &
                     c18_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q18) + &
                     c23_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q23) + &
                     c24_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q24) + &
                     c26_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q26)))
!
!perform angular integration
                  call calc_phinorm(vel_p, vth_p, vth_fiducial, xobs, phinorm)
                  wall=wall_global*phinorm
                  mintbar3d_tmp(i,j,k) = mintbar3d_tmp(i,j,k) + int3d(i,j,k)*wall
                  normalization3d_tmp(i,j,k) = normalization3d_tmp(i,j,k) + wall
!
                  aloline_nn3d_tmp(iip1,jjp1,kkp1,q1) = aloline_nn3d_tmp(iip1,jjp1,kkp1,q1) + wall*aloline_on_nn3d(iip1,jjp1,kkp1,q1)
                  aloline_nn3d_tmp(i,jjp1,kkp1,q2) = aloline_nn3d_tmp(i,jjp1,kkp1,q2) + wall*aloline_on_nn3d(i,jjp1,kkp1,q2)
                  aloline_nn3d_tmp(iim1,jjp1,kkp1,q3) = aloline_nn3d_tmp(iim1,jjp1,kkp1,q3) + wall*aloline_on_nn3d(iim1,jjp1,kkp1,q3)
                  aloline_nn3d_tmp(iip1,j,kkp1,q4) = aloline_nn3d_tmp(iip1,j,kkp1,q4) + wall*aloline_on_nn3d(iip1,j,kkp1,q4)
                  aloline_nn3d_tmp(i,j,kkp1,q5) = aloline_nn3d_tmp(i,j,kkp1,q5) + wall*aloline_on_nn3d(i,j,kkp1,q5)
                  aloline_nn3d_tmp(iim1,j,kkp1,q6) = aloline_nn3d_tmp(iim1,j,kkp1,q6) + wall*aloline_on_nn3d(iim1,j,kkp1,q6)
                  aloline_nn3d_tmp(iip1,jjm1,kkp1,q7) = aloline_nn3d_tmp(iip1,jjm1,kkp1,q7) + wall*aloline_on_nn3d(iip1,jjm1,kkp1,q7)
                  aloline_nn3d_tmp(i,jjm1,kkp1,q8) = aloline_nn3d_tmp(i,jjm1,kkp1,q8) + wall*aloline_on_nn3d(i,jjm1,kkp1,q8)
                  aloline_nn3d_tmp(iim1,jjm1,kkp1,q9) = aloline_nn3d_tmp(iim1,jjm1,kkp1,q9) + wall*aloline_on_nn3d(iim1,jjm1,kkp1,q9)
                  aloline_nn3d_tmp(iip1,jjp1,k,q10) = aloline_nn3d_tmp(iip1,jjp1,k,q10) + wall*aloline_on_nn3d(iip1,jjp1,k,q10)
                  aloline_nn3d_tmp(i,jjp1,k,q11) = aloline_nn3d_tmp(i,jjp1,k,q11) + wall*aloline_on_nn3d(i,jjp1,k,q11)
                  aloline_nn3d_tmp(iim1,jjp1,k,q12) = aloline_nn3d_tmp(iim1,jjp1,k,q12) + wall*aloline_on_nn3d(iim1,jjp1,k,q12)
                  aloline_nn3d_tmp(iip1,j,k,q13) = aloline_nn3d_tmp(iip1,j,k,q13) + wall*aloline_on_nn3d(iip1,j,k,q13)
                  aloline_nn3d_tmp(i,j,k,q14) = aloline_nn3d_tmp(i,j,k,q14) + wall*aloline_on_nn3d(i,j,k,q14)
                  aloline_nn3d_tmp(iim1,j,k,q15) = aloline_nn3d_tmp(iim1,j,k,q15) + wall*aloline_on_nn3d(iim1,j,k,q15)
                  aloline_nn3d_tmp(iip1,jjm1,k,q16) = aloline_nn3d_tmp(iip1,jjm1,k,q16) + wall*aloline_on_nn3d(iip1,jjm1,k,q16) !new
                  aloline_nn3d_tmp(i,jjm1,k,q17) = aloline_nn3d_tmp(i,jjm1,k,q17) + wall*aloline_on_nn3d(i,jjm1,k,q17)
                  aloline_nn3d_tmp(iim1,jjm1,k,q18) = aloline_nn3d_tmp(iim1,jjm1,k,q18) + wall*aloline_on_nn3d(iim1,jjm1,k,q18)
                  aloline_nn3d_tmp(iip1,jjp1,kkm1,q19) = aloline_nn3d_tmp(iip1,jjp1,kkm1,q19) + wall*aloline_on_nn3d(iip1,jjp1,kkm1,q19) !new
                  aloline_nn3d_tmp(i,jjp1,kkm1,q20) = aloline_nn3d_tmp(i,jjp1,kkm1,q20) + wall*aloline_on_nn3d(i,jjp1,kkm1,q20) !new
                  aloline_nn3d_tmp(iim1,jjp1,kkm1,q21) = aloline_nn3d_tmp(iim1,jjp1,kkm1,q21) + wall*aloline_on_nn3d(iim1,jjp1,kkm1,q21) !new
                  aloline_nn3d_tmp(iip1,j,kkm1,q22) = aloline_nn3d_tmp(iip1,j,kkm1,q22) + wall*aloline_on_nn3d(iip1,j,kkm1,q22) !new
                  aloline_nn3d_tmp(i,j,kkm1,q23) = aloline_nn3d_tmp(i,j,kkm1,q23) + wall*aloline_on_nn3d(i,j,kkm1,q23)
                  aloline_nn3d_tmp(iim1,j,kkm1,q24) = aloline_nn3d_tmp(iim1,j,kkm1,q24) + wall*aloline_on_nn3d(iim1,j,kkm1,q24)
                  aloline_nn3d_tmp(iip1,jjm1,kkm1,q25) = aloline_nn3d_tmp(iip1,jjm1,kkm1,q25) + wall*aloline_on_nn3d(iip1,jjm1,kkm1,q25) !new
                  aloline_nn3d_tmp(i,jjm1,kkm1,q26) = aloline_nn3d_tmp(i,jjm1,kkm1,q26) + wall*aloline_on_nn3d(i,jjm1,kkm1,q26)
                  aloline_nn3d_tmp(iim1,jjm1,kkm1,q27) = aloline_nn3d_tmp(iim1,jjm1,kkm1,q27) + wall*aloline_on_nn3d(iim1,jjm1,kkm1,q27)
!
               endif
!
             case default
!
            end select
!
         enddo
      enddo
   enddo
!
!***debug start
!close(1)
!***debug end
!
!
end subroutine fsc_line3d_linb
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine coeff2d_lineu(opalbar_im2jm2, opalbar_im1jm2, opalbar_ijm2, &
   opalbar_im2jm1, opalbar_im1jm1, opalbar_ijm1, &
   opalbar_im2j,   opalbar_im1j,   opalbar_ij, &
   velx_im2jm2, velx_im1jm2, velx_ijm2, &
   velx_im2jm1, velx_im1jm1, velx_ijm1, &
   velx_im2j,   velx_im1j,   velx_ij, &
   vely_im2jm2, vely_im1jm2, vely_ijm2, &
   vely_im2jm1, vely_im1jm1, vely_ijm1, &
   vely_im2j,   vely_im1j,   vely_ij, &
   velz_im2jm2, velz_im1jm2, velz_ijm2, &
   velz_im2jm1, velz_im1jm1, velz_ijm1, &
   velz_im2j,   velz_im1j,   velz_ij, &
   vth_im2jm2, vth_im1jm2, vth_ijm2, &
   vth_im2jm1, vth_im1jm1, vth_ijm1, &
   vth_im2j,   vth_im1j,   vth_ij, &
   sline_im2jm2, sline_im1jm2, sline_ijm2, &
   sline_im2jm1, sline_im1jm1, sline_ijm1, &
   sline_im2j,   sline_im1j,   sline_ij, &
   int_im2jm2, int_im1jm2, int_ijm2, &
   int_im2jm1, int_im1jm1, int_ijm1, &
   int_im2j,   int_im1j,   int_ij, &
   x_im2, x_im1, x_i, y_jm2, y_jm1, y_j, x_p, y_p, &
   a_sline, b_sline, c_sline, d_sline, e_sline, &
   f_sline, g_sline, h_sline, i_sline, &
   a_inten, b_inten, c_inten, d_inten, e_inten, &
   f_inten, g_inten, h_inten, i_inten, opalbar_p, velx_p, vely_p, velz_p, vth_p, sline_p, int_p)
!
!         interpolates opaciy, velocity components, thermal velocity,
!                   line source function and intensity
!               values given on a 2d grid onto point x_p, y_p
!
!on input (f_* stands for opalbar_*, sline_* velx_*, vely_*, velz_*, vth_*, and int_*, respectivly):
!
! y_j      f_im2j----------f_im1j--------------f_ij
!  |          |              |                  |
!  |          |              |                  |
!  |          |              |         x        |
!  |          |              |     (x_p,y_p)    |
!  |          |              |                  |
!y_jm1    f_im2jm1-------f_im1jm1------------f_ijm1
!  |          |              |                  |
!  |          |              |                  |
!  |          |              |                  |
!  |          |              |                  |
!  |          |              |                  |
!y_jm2    f_im2jm2-------f_im1jm2------------f_ijm2
!  |
!  --------x_im2-----------x_im1---------------x_i
!
!        x_p, y_p: coordinates of point onto which shall be interpolated
!
!on output:
!   1. interpolation coefficients for source function and intensity
!      (required for ALO calculations):
!         a_sline, b_sline, c_sline, d_sline, e_sline,
!         f_sline, g_sline, h_sline, i_sline
!         a_inten, b_inten, c_inten, d_inten, e_inten
!         f_inten, g_inten, h_inten, i_inten
!
!      such that:
!         f_p = a*f_im2jm2 + b*f_im1jm2 + c*f_ijm2 +
!               d*f_im2jm1 + e*f_im1jm1 + f*f_ijm1 +
!               g*f_im2j   + h*f_im1j   + i*f_ij
!
!   2. interpolated values at point p: opalbar_p, sline_p, velx_p, vely_p, velz_p,
!                                      vth_p, int_p
!
   use prog_type
   use mod_interp2d, only: wp_interp2d
!
   implicit none
!
! ... arguments
   real(dp), intent(in) :: opalbar_im2jm2, opalbar_im1jm2, opalbar_ijm2, &
      opalbar_im2jm1, opalbar_im1jm1, opalbar_ijm1, &
      opalbar_im2j,   opalbar_im1j,   opalbar_ij, &
      velx_im2jm2, velx_im1jm2, velx_ijm2, &
      velx_im2jm1, velx_im1jm1, velx_ijm1, &
      velx_im2j,   velx_im1j,   velx_ij, &
      vely_im2jm2, vely_im1jm2, vely_ijm2, &
      vely_im2jm1, vely_im1jm1, vely_ijm1, &
      vely_im2j,   vely_im1j,   vely_ij, &
      velz_im2jm2, velz_im1jm2, velz_ijm2, &
      velz_im2jm1, velz_im1jm1, velz_ijm1, &
      velz_im2j,   velz_im1j,   velz_ij, &
      vth_im2jm2, vth_im1jm2, vth_ijm2, &
      vth_im2jm1, vth_im1jm1, vth_ijm1, &
      vth_im2j,   vth_im1j,   vth_ij, &
      sline_im2jm2, sline_im1jm2, sline_ijm2, &
      sline_im2jm1, sline_im1jm1, sline_ijm1, &
      sline_im2j,   sline_im1j,   sline_ij, &
      int_im2jm2, int_im1jm2, int_ijm2, &
      int_im2jm1, int_im1jm1, int_ijm1, &
      int_im2j,   int_im1j,   int_ij, &
      x_im2, x_im1, x_i, y_jm2, y_jm1, y_j, x_p, y_p
   real(dp), intent(out) :: a_sline, b_sline, c_sline, d_sline, e_sline, &
      f_sline, g_sline, h_sline, i_sline, &
      a_inten, b_inten, c_inten, d_inten, e_inten, &
      f_inten, g_inten, h_inten, i_inten, &
      opalbar_p, velx_p, vely_p, velz_p, vth_p, int_p, sline_p
!
! ... local scalars
   real(dp) :: dxim1, dxi, dx, tx, ax, bx, cx, &
      dyjm1, dyj, dy, ty, ay, by, cy, &
      axt, bxt, cxt, axt2, bxt2, cxt2, &
      ayt, byt, cyt, ayt2, byt2, cyt2, &
      axtjm2_sline, bxtjm2_sline, cxtjm2_sline, &
      axtjm1_sline, bxtjm1_sline, cxtjm1_sline, &
      axtj_sline, bxtj_sline, cxtj_sline, &
      axjm2_sline, bxjm2_sline, cxjm2_sline, &
      axjm1_sline, bxjm1_sline, cxjm1_sline, &
      axj_sline, bxj_sline, cxj_sline, &
      ayt_sline, byt_sline, cyt_sline, &
      ay_sline, by_sline, cy_sline, &
      axtjm2_int, bxtjm2_int, cxtjm2_int, &
      axtjm1_int, bxtjm1_int, cxtjm1_int, &
      axtj_int, bxtj_int, cxtj_int, &
      axjm2_int, bxjm2_int, cxjm2_int, &
      axjm1_int, bxjm1_int, cxjm1_int, &
      axj_int, bxj_int, cxj_int, &
      ay_int, by_int, cy_int, &
      ayt_int, byt_int, cyt_int, &
      opalbar_jm2, opalbar_jm1, opalbar_j, opalbarc_jm2, opalbarc_jm1, opalbarc_j, opalbar_c, &
      velx_jm2, velx_jm1, velx_j, velxc_jm2, velxc_jm1, velxc_j, velx_c, &
      vely_jm2, vely_jm1, vely_j, velyc_jm2, velyc_jm1, velyc_j, vely_c, &
      velz_jm2, velz_jm1, velz_j, velzc_jm2, velzc_jm1, velzc_j, velz_c, &
      vth_jm2, vth_jm1, vth_j, vthc_jm2, vthc_jm1, vthc_j, vth_c, &
      sline_jm2, sline_jm1, sline_j, slinec_jm2, slinec_jm1, slinec_j, sline_c, &
      int_jm2, int_jm1, int_j, intc_jm2, intc_jm1, intc_j, int_c!, &

   real(dp) :: fac, fac2
!
!
!for linear approach:
!call coeff2d_lineu_lin(opalbar_im1jm1, opalbar_ijm1, opalbar_im1j, opalbar_ij, &
!                       velx_im1jm1, velx_ijm1, velx_im1j, velx_ij, &
!                       vely_im1jm1, vely_ijm1, vely_im1j, vely_ij, &
!                       velz_im1jm1, velz_ijm1, velz_im1j, velz_ij, &
!                       vth_im1jm1, vth_ijm1, vth_im1j, vth_ij, &
!                       sline_im1jm1, sline_ijm1, sline_im1j, sline_ij, &
!                       int_im1jm1, int_ijm1, int_im1j, int_ij, &
!                       x_im1, x_i, y_jm1, y_j, x_p, y_p, &
!                       e_sline, f_sline, h_sline, i_sline, &
!                       e_inten, f_inten, h_inten, i_inten, &
!                       opalbar_p, velx_p, vely_p, velz_p, vth_p, sline_p, int_p)
!!
!a_sline=0.d0
!b_sline=0.d0
!c_sline=0.d0
!d_sline=0.d0
!g_sline=0.d0
!!
!a_inten=0.d0
!b_inten=0.d0
!c_inten=0.d0
!d_inten=0.d0
!g_inten=0.d0
!!
!return
!
!---------------monotonic quadratic bezier interpolation----------------
!--------(with derivatives for control point from weighted mean)--------
!--------------------where weights are assigned-------------------------
!
!define deltax, deltay
   dxim1 = x_im1-x_im2
   dxi = x_i-x_im1
   dx = dxim1+dxi
   dyjm1 = y_jm1-y_jm2
   dyj = y_j-y_jm1
   dy = dyjm1+dyj
!
!define deltax, deltay-ratios
   tx = (x_p-x_im1)/dxi
   ty = (y_p-y_jm1)/dyj
!
   ax = (1.d0-tx)**2
   bx = 2.d0*tx*(1.d0-tx)
   cx = tx**2
!
!----------------calculate control points on each j-level---------------
!   slinec_jm2, slinec_jm1, slinec_j
!   intc_jm2, intc_jm1, intc_j
!   opacc_jm2, opacc_jm1, opacc_j
!
!derivative weights for source function and intensity
   fac=max(wp_interp2d,dxim1/dx)
!derivative weights for velocity components and opacity
   fac2=dxim1/dx
!
   axt = (fac-1.d0)*dxi/2.d0/dxim1
   bxt = ((2.d0-fac)*dxim1 + (1.d0-fac)*dxi)/2.d0/dxim1
   cxt = fac/2.d0
   axt2 = (fac2-1.d0)*dxi/2.d0/dxim1
   bxt2 = ((2.d0-fac2)*dxim1 + (1.d0-fac2)*dxi)/2.d0/dxim1
   cxt2 = fac2/2.d0
!
   opalbarc_jm2 = opalbar_im2jm2*axt2 + opalbar_im1jm2*bxt2 + opalbar_ijm2*cxt2
   opalbarc_jm1 = opalbar_im2jm1*axt2 + opalbar_im1jm1*bxt2 + opalbar_ijm1*cxt2
   opalbarc_j   = opalbar_im2j*axt2   + opalbar_im1j*bxt2   + opalbar_ij*cxt2
!
   velxc_jm2 = velx_im2jm2*axt2 + velx_im1jm2*bxt2 + velx_ijm2*cxt2
   velxc_jm1 = velx_im2jm1*axt2 + velx_im1jm1*bxt2 + velx_ijm1*cxt2
   velxc_j   = velx_im2j*axt2   + velx_im1j*bxt2   + velx_ij*cxt2
!
   velyc_jm2 = vely_im2jm2*axt2 + vely_im1jm2*bxt2 + vely_ijm2*cxt2
   velyc_jm1 = vely_im2jm1*axt2 + vely_im1jm1*bxt2 + vely_ijm1*cxt2
   velyc_j   = vely_im2j*axt2   + vely_im1j*bxt2   + vely_ij*cxt2
!
   velzc_jm2 = velz_im2jm2*axt2 + velz_im1jm2*bxt2 + velz_ijm2*cxt2
   velzc_jm1 = velz_im2jm1*axt2 + velz_im1jm1*bxt2 + velz_ijm1*cxt2
   velzc_j   = velz_im2j*axt2   + velz_im1j*bxt2   + velz_ij*cxt2
!
   vthc_jm2 = vth_im2jm2*axt2 + vth_im1jm2*bxt2 + vth_ijm2*cxt2
   vthc_jm1 = vth_im2jm1*axt2 + vth_im1jm1*bxt2 + vth_ijm1*cxt2
   vthc_j   = vth_im2j*axt2   + vth_im1j*bxt2   + vth_ij*cxt2
!
   slinec_jm2 = sline_im2jm2*axt + sline_im1jm2*bxt + sline_ijm2*cxt
   slinec_jm1 = sline_im2jm1*axt + sline_im1jm1*bxt + sline_ijm1*cxt
   slinec_j   = sline_im2j*axt   + sline_im1j*bxt   + sline_ij*cxt
!
   intc_jm2 = int_im2jm2*axt + int_im1jm2*bxt + int_ijm2*cxt
   intc_jm1 = int_im2jm1*axt + int_im1jm1*bxt + int_ijm1*cxt
   intc_j   = int_im2j*axt   + int_im1j*bxt   + int_ij*cxt
!
!------------------ensure monotonicity on level j-----------------------
!
!velocity components (no monotonicity required)
   velx_j = ax*velx_im1j + bx*velxc_j + cx*velx_ij
   vely_j = ax*vely_im1j + bx*velyc_j + cx*vely_ij
   velz_j = ax*velz_im1j + bx*velzc_j + cx*velz_ij
!
!opacity and thermal velocity (no interpolation coefficients required)
   call pointcr1d_mbez(vthc_j, vth_im1j, vth_ij)
   call pointcr1d_mbez(opalbarc_j, opalbar_im1j, opalbar_ij)
   vth_j = ax*vth_im1j + bx*vthc_j + cx*vth_ij
   opalbar_j = ax*opalbar_im1j + bx*opalbarc_j + cx*opalbar_ij
!
!line source function and intensity
   call coeffcr1d_mbez(slinec_j, sline_im1j, sline_ij, axt, bxt, cxt, axtj_sline, bxtj_sline, cxtj_sline)
   call coeffcr1d_mbez(intc_j, int_im1j, int_ij, axt, bxt, cxt, axtj_int, bxtj_int, cxtj_int)
   axj_sline = axtj_sline*bx
   bxj_sline = bxtj_sline*bx + ax
   cxj_sline = cxtj_sline*bx + cx
   sline_j = axj_sline*sline_im2j + bxj_sline*sline_im1j + cxj_sline*sline_ij
   axj_int = axtj_int*bx
   bxj_int = bxtj_int*bx + ax
   cxj_int = cxtj_int*bx + cx
   int_j = axj_int*int_im2j + bxj_int*int_im1j + cxj_int*int_ij
!
!---------------ensure monotonicity on level j-1------------------------
!
!velocity components (no monotonicity required)
   velx_jm1 = ax*velx_im1jm1 + bx*velxc_jm1 + cx*velx_ijm1
   vely_jm1 = ax*vely_im1jm1 + bx*velyc_jm1 + cx*vely_ijm1
   velz_jm1 = ax*velz_im1jm1 + bx*velzc_jm1 + cx*velz_ijm1
!
!opacity and thermal velocity (no interpolation coefficients required)
   call pointcr1d_mbez(vthc_jm1, vth_im1jm1, vth_ijm1)
   call pointcr1d_mbez(opalbarc_jm1, opalbar_im1jm1, opalbar_ijm1)
   vth_jm1 = ax*vth_im1jm1 + bx*vthc_jm1 + cx*vth_ijm1
   opalbar_jm1 = ax*opalbar_im1jm1 + bx*opalbarc_jm1 + cx*opalbar_ijm1
!
!line source function and intensity
   call coeffcr1d_mbez(slinec_jm1, sline_im1jm1, sline_ijm1, axt, bxt, cxt, axtjm1_sline, bxtjm1_sline, cxtjm1_sline)
   call coeffcr1d_mbez(intc_jm1, int_im1jm1, int_ijm1, axt, bxt, cxt, axtjm1_int, bxtjm1_int, cxtjm1_int)
   axjm1_sline = axtjm1_sline*bx
   bxjm1_sline = bxtjm1_sline*bx + ax
   cxjm1_sline = cxtjm1_sline*bx + cx
   sline_jm1 = axjm1_sline*sline_im2jm1 + bxjm1_sline*sline_im1jm1 + cxjm1_sline*sline_ijm1
   axjm1_int = axtjm1_int*bx
   bxjm1_int = bxtjm1_int*bx + ax
   cxjm1_int = cxtjm1_int*bx + cx
   int_jm1 = axjm1_int*int_im2jm1 + bxjm1_int*int_im1jm1 + cxjm1_int*int_ijm1
!
!---------------ensure monotonicity on level j-2------------------------
!
!velocity components (no monotonicity required)
   velx_jm2 = ax*velx_im1jm2 + bx*velxc_jm2 + cx*velx_ijm2
   vely_jm2 = ax*vely_im1jm2 + bx*velyc_jm2 + cx*vely_ijm2
   velz_jm2 = ax*velz_im1jm2 + bx*velzc_jm2 + cx*velz_ijm2
!
!opacity and thermal velocity (no interpolation coefficients required)
   call pointcr1d_mbez(vthc_jm2, vth_im1jm2, vth_ijm2)
   call pointcr1d_mbez(opalbarc_jm2, opalbar_im1jm2, opalbar_ijm2)
   vth_jm2 = ax*vth_im1jm2 + bx*vthc_jm2 + cx*vth_ijm2
   opalbar_jm2 = ax*opalbar_im1jm2 + bx*opalbarc_jm2 + cx*opalbar_ijm2
!
!line source function and intensity
   call coeffcr1d_mbez(slinec_jm2, sline_im1jm2, sline_ijm2, axt, bxt, cxt, axtjm2_sline, bxtjm2_sline, cxtjm2_sline)
   call coeffcr1d_mbez(intc_jm2, int_im1jm2, int_ijm2, axt, bxt, cxt, axtjm2_int, bxtjm2_int, cxtjm2_int)
   axjm2_sline = axtjm2_sline*bx
   bxjm2_sline = bxtjm2_sline*bx + ax
   cxjm2_sline = cxtjm2_sline*bx + cx
   sline_jm2 = axjm2_sline*sline_im2jm2 + bxjm2_sline*sline_im1jm2 + cxjm2_sline*sline_ijm2
   axjm2_int = axtjm2_int*bx
   bxjm2_int = bxtjm2_int*bx + ax
   cxjm2_int = cxtjm2_int*bx + cx
   int_jm2 = axjm2_int*int_im2jm2 + bxjm2_int*int_im1jm2 + cxjm2_int*int_ijm2
!
!------------calculate control point for interpolation along y----------
!
!derivative weights for source function and intensity
   fac=max(wp_interp2d,dyjm1/dy)
!derivative weights for velocity components and opacity
   fac2=dyjm1/dy
!
   ayt = (fac-1.d0)*dyj/2.d0/dyjm1
   byt = ((2.d0-fac)*dyjm1 + (1.d0-fac)*dyj)/2.d0/dyjm1
   cyt = fac/2.d0
   ayt2 = (fac2-1.d0)*dyj/2.d0/dyjm1
   byt2 = ((2.d0-fac2)*dyjm1 + (1.d0-fac2)*dyj)/2.d0/dyjm1
   cyt2 = fac2/2.d0
!
   ay = (1.d0-ty)**2
   by = 2.d0*ty*(1.d0-ty)
   cy = ty**2
!
   velx_c = velx_jm2*ayt2 + velx_jm1*byt2 + velx_j*cyt2
   vely_c = vely_jm2*ayt2 + vely_jm1*byt2 + vely_j*cyt2
   velz_c = velz_jm2*ayt2 + velz_jm1*byt2 + velz_j*cyt2
   vth_c = vth_jm2*ayt2 + vth_jm1*byt2 + vth_j*cyt2
   opalbar_c = opalbar_jm2*ayt2 + opalbar_jm1*byt2 + opalbar_j*cyt2
!
   sline_c = sline_jm2*ayt + sline_jm1*byt + sline_j*cyt
   int_c = int_jm2*ayt + int_jm1*byt + int_j*cyt
!
!------------------------ensure monotonicity----------------------------
!
!velocity components (no monotonicity required)
   velx_p = ay*velx_jm1 + by*velx_c + cy*velx_j
   vely_p = ay*vely_jm1 + by*vely_c + cy*vely_j
   velz_p = ay*velz_jm1 + by*velz_c + cy*velz_j
!
!opacity and thermal velocity (no interpolation coefficients required)
   call pointcr1d_mbez(vth_c, vth_jm1, vth_j)
   call pointcr1d_mbez(opalbar_c, opalbar_jm1, opalbar_j)
   vth_p = ay*vth_jm1 + by*vth_c + cy*vth_j
   opalbar_p = ay*opalbar_jm1 + by*opalbar_c + cy*opalbar_j
!
!line source function and intensity
   call coeffcr1d_mbez(sline_c, sline_jm1, sline_j, ayt, byt, cyt, ayt_sline, byt_sline, cyt_sline)
   call coeffcr1d_mbez(int_c, int_jm1, int_j, ayt, byt, cyt, ayt_int, byt_int, cyt_int)
   ay_sline = ayt_sline*by
   by_sline = byt_sline*by + ay
   cy_sline = cyt_sline*by + cy
   sline_p = ay_sline*sline_jm2 + by_sline*sline_jm1 + cy_sline*sline_j
   ay_int = ayt_int*by
   by_int = byt_int*by + ay
   cy_int = cyt_int*by + cy
   int_p = ay_int*int_jm2 + by_int*int_jm1 + cy_int*int_j
!
   a_sline = ay_sline*axjm2_sline
   b_sline = ay_sline*bxjm2_sline
   c_sline = ay_sline*cxjm2_sline
   d_sline = by_sline*axjm1_sline
   e_sline = by_sline*bxjm1_sline
   f_sline = by_sline*cxjm1_sline
   g_sline = cy_sline*axj_sline
   h_sline = cy_sline*bxj_sline
   i_sline = cy_sline*cxj_sline
!
   a_inten = ay_int*axjm2_int
   b_inten = ay_int*bxjm2_int
   c_inten = ay_int*cxjm2_int
   d_inten = by_int*axjm1_int
   e_inten = by_int*bxjm1_int
   f_inten = by_int*cxjm1_int
   g_inten = cy_int*axj_int
   h_inten = cy_int*bxj_int
   i_inten = cy_int*cxj_int
!
!write(*,'(11es20.8)') x_p, y_p, z_p, opalbar_p, opalbar_p2, velx_p, velx_p2, vely_p, vely_p2, velz_p, velz_p2
!if(int_p.ne.int_p) then
!   write(*,*) 'error: int_p=', int_p
!   write(*,*)
!   write(*,*) int_im2jm2, int_im1jm2, int_ijm2, &
!              int_im2jm1, int_im1jm1, int_ijm1, &
!              int_im2j,   int_im1j, int_ij
!   stop
!endif
!
   return
!
end subroutine coeff2d_lineu

!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine coeff2d_lined(opalbar_im2jm2, opalbar_im1jm2, opalbar_ijm2, &
   opalbar_im2jm1, opalbar_im1jm1, opalbar_ijm1, &
   opalbar_im2j,   opalbar_im1j,   opalbar_ij, &
   velx_im2jm2, velx_im1jm2, velx_ijm2, &
   velx_im2jm1, velx_im1jm1, velx_ijm1, &
   velx_im2j,   velx_im1j,   velx_ij, &
   vely_im2jm2,  vely_im1jm2,  vely_ijm2, &
   vely_im2jm1,  vely_im1jm1, vely_ijm1, &
   vely_im2j,   vely_im1j,   vely_ij, &
   velz_im2jm2, velz_im1jm2, velz_ijm2, &
   velz_im2jm1, velz_im1jm1, velz_ijm1, &
   velz_im2j,   velz_im1j,   velz_ij, &
   vth_im2jm2, vth_im1jm2, vth_ijm2, &
   vth_im2jm1, vth_im1jm1, vth_ijm1, &
   vth_im2j,   vth_im1j,   vth_ij, &
   sline_im2jm2, sline_im1jm2, sline_ijm2, &
   sline_im2jm1, sline_im1jm1, sline_ijm1, &
   sline_im2j,   sline_im1j,   sline_ij, &
   x_im2, x_im1, x_i, y_jm2, y_jm1, y_j, x_p, y_p, &
   a_sline, b_sline, c_sline, d_sline, e_sline, &
   f_sline, g_sline, h_sline, i_sline, &
   opalbar_p, velx_p, vely_p, velz_p, vth_p, sline_p)
!
!         interpolates opaciy, velocity components, thermal velocity,
!                   line source function and intensity
!               values given on a 2d grid onto point x_p, y_p
!
!on input (f_* stands for opalbar_*, sline_* velx_*, vely_*, velz_*, vth_*, and int_*, respectivly):
!
! y_j      f_im2j----------f_im1j--------------f_ij
!  |          |              |                  |
!  |          |              |                  |
!  |          |              |         x        |
!  |          |              |     (x_p,y_p)    |
!  |          |              |                  |
!y_jm1    f_im2jm1-------f_im1jm1------------f_ijm1
!  |          |              |                  |
!  |          |              |                  |
!  |          |              |                  |
!  |          |              |                  |
!  |          |              |                  |
!y_jm2    f_im2jm2-------f_im1jm2------------f_ijm2
!  |
!  --------x_im2-----------x_im1---------------x_i
!
!        x_p, y_p: coordinates of point onto which shall be interpolated
!
!on output:
!   1. interpolation coefficients for source function
!      (required for ALO calculations):
!         a_sline, b_sline, c_sline, d_sline, e_sline
!         f_sline, g_sline, h_sline, i_sline
!
!      such that:
!         f_p = a*f_im2jm2 + b*f_im1jm2 + c*f_ijm2 +
!               d*f_im2jm1 + e*f_im1jm1 + f*f_ijm1 +
!               g*f_im2j   + h*f_im1j   + i*f_ij
!
!   2. interpolated values at point p: opalbar_p, sline_p, velx_p, vely_p, velz_p,
!                                      vth_p
!
   use prog_type
   use mod_interp2d, only: wp_interp2d
!
   implicit none
!
! ... arguments
   real(dp), intent(in) :: opalbar_im2jm2, opalbar_im1jm2, opalbar_ijm2, &
      opalbar_im2jm1, opalbar_im1jm1, opalbar_ijm1, &
      opalbar_im2j,   opalbar_im1j,   opalbar_ij, &
      velx_im2jm2, velx_im1jm2, velx_ijm2, &
      velx_im2jm1, velx_im1jm1, velx_ijm1, &
      velx_im2j,   velx_im1j,   velx_ij, &
      vely_im2jm2, vely_im1jm2, vely_ijm2, &
      vely_im2jm1, vely_im1jm1, vely_ijm1, &
      vely_im2j,   vely_im1j,   vely_ij, &
      velz_im2jm2, velz_im1jm2, velz_ijm2, &
      velz_im2jm1, velz_im1jm1, velz_ijm1, &
      velz_im2j,   velz_im1j,   velz_ij, &
      vth_im2jm2, vth_im1jm2, vth_ijm2, &
      vth_im2jm1, vth_im1jm1, vth_ijm1, &
      vth_im2j,   vth_im1j,   vth_ij, &
      sline_im2jm2, sline_im1jm2, sline_ijm2, &
      sline_im2jm1, sline_im1jm1, sline_ijm1, &
      sline_im2j,   sline_im1j,   sline_ij, &
      x_im2, x_im1, x_i, y_jm2, y_jm1, y_j, x_p, y_p
   real(dp), intent(out) :: a_sline, b_sline, c_sline, d_sline, e_sline, &
      f_sline, g_sline, h_sline, i_sline, &
      opalbar_p, velx_p, vely_p, velz_p, vth_p, sline_p
!
! ... local scalars
   real(dp) :: dxim1, dxi, dx, tx, ax, bx, cx, &
      dyjm1, dyj, dy, ty, ay, by, cy, &
      axt, bxt, cxt, axt2, bxt2, cxt2, &
      ayt, byt, cyt, ayt2, byt2, cyt2, &
      axtjm2_sline, bxtjm2_sline, cxtjm2_sline, &
      axtjm1_sline, bxtjm1_sline, cxtjm1_sline, &
      axtj_sline, bxtj_sline, cxtj_sline, &
      axjm2_sline, bxjm2_sline, cxjm2_sline, &
      axjm1_sline, bxjm1_sline, cxjm1_sline, &
      axj_sline, bxj_sline, cxj_sline, &
      ayt_sline, byt_sline, cyt_sline, &
      ay_sline, by_sline, cy_sline, &
      opalbar_jm2, opalbar_jm1, opalbar_j, opalbarc_jm2, opalbarc_jm1, opalbarc_j, opalbar_c, &
      velx_jm2, velx_jm1, velx_j, velxc_jm2, velxc_jm1, velxc_j, velx_c, &
      vely_jm2, vely_jm1, vely_j, velyc_jm2, velyc_jm1, velyc_j, vely_c, &
      velz_jm2, velz_jm1, velz_j, velzc_jm2, velzc_jm1, velzc_j, velz_c, &
      vth_jm2, vth_jm1, vth_j, vthc_jm2, vthc_jm1, vthc_j, vth_c, &
      sline_jm2, sline_jm1, sline_j, slinec_jm2, slinec_jm1, slinec_j, sline_c
!
   real(dp) :: fac, fac2
!


!for linear approach:
!call coeff2d_lined_lin(opalbar_im1jm1, opalbar_ijm1, opalbar_im1j, opalbar_ij, &
!                       velx_im1jm1, velx_ijm1, velx_im1j, velx_ij, &
!                       vely_im1jm1, vely_ijm1, vely_im1j, vely_ij, &
!                       velz_im1jm1, velz_ijm1, velz_im1j, velz_ij, &
!                       vth_im1jm1, vth_ijm1, vth_im1j, vth_ij, &
!                       sline_im1jm1, sline_ijm1, sline_im1j, sline_ij, &
!                       x_im1, x_i, y_jm1, y_j, x_p, y_p, &
!                       e_sline, f_sline, h_sline, i_sline, &
!                       opalbar_p, velx_p, vely_p, velz_p, vth_p, sline_p)
!a_sline=0.d0
!b_sline=0.d0
!c_sline=0.d0
!d_sline=0.d0
!g_sline=0.d0
!
!return
!
!---------------monotonic quadratic bezier interpolation----------------
!--------(with derivatives for control point from weighted mean)--------
!--------------------where weights are assigned-------------------------
!
!define deltax, deltay
   dxim1 = x_im1-x_im2
   dxi = x_i-x_im1
   dx = dxim1+dxi
   dyjm1 = y_jm1-y_jm2
   dyj = y_j-y_jm1
   dy = dyjm1+dyj
!
!define deltax, deltay-ratios
   tx = (x_p-x_im1)/dxi
   ty = (y_p-y_jm1)/dyj
!
   ax = (1.d0-tx)**2
   bx = 2.d0*tx*(1.d0-tx)
   cx = tx**2
!
!----------------calculate control points on each j-level---------------
!   slinec_jm2, slinec_jm1, slinec_j
!   opalbarc_jm2, opalbarc_jm1, opalbarc_j
!
!derivative weights for source function and intensity
   fac=max(wp_interp2d,dxim1/dx)
!derivative weights for velocity components and opacity
   fac2=dxim1/dx
!
   axt = (fac-1.d0)*dxi/2.d0/dxim1
   bxt = ((2.d0-fac)*dxim1 + (1.d0-fac)*dxi)/2.d0/dxim1
   cxt = fac/2.d0
   axt2 = (fac2-1.d0)*dxi/2.d0/dxim1
   bxt2 = ((2.d0-fac2)*dxim1 + (1.d0-fac2)*dxi)/2.d0/dxim1
   cxt2 = fac2/2.d0
!
   opalbarc_jm2 = opalbar_im2jm2*axt2 + opalbar_im1jm2*bxt2 + opalbar_ijm2*cxt2
   opalbarc_jm1 = opalbar_im2jm1*axt2 + opalbar_im1jm1*bxt2 + opalbar_ijm1*cxt2
   opalbarc_j   = opalbar_im2j*axt2   + opalbar_im1j*bxt2   + opalbar_ij*cxt2
!
   velxc_jm2 = velx_im2jm2*axt2 + velx_im1jm2*bxt2 + velx_ijm2*cxt2
   velxc_jm1 = velx_im2jm1*axt2 + velx_im1jm1*bxt2 + velx_ijm1*cxt2
   velxc_j   = velx_im2j*axt2   + velx_im1j*bxt2   + velx_ij*cxt2
!
   velyc_jm2 = vely_im2jm2*axt2 + vely_im1jm2*bxt2 + vely_ijm2*cxt2
   velyc_jm1 = vely_im2jm1*axt2 + vely_im1jm1*bxt2 + vely_ijm1*cxt2
   velyc_j   = vely_im2j*axt2   + vely_im1j*bxt2   + vely_ij*cxt2
!
   velzc_jm2 = velz_im2jm2*axt2 + velz_im1jm2*bxt2 + velz_ijm2*cxt2
   velzc_jm1 = velz_im2jm1*axt2 + velz_im1jm1*bxt2 + velz_ijm1*cxt2
   velzc_j   = velz_im2j*axt2   + velz_im1j*bxt2   + velz_ij*cxt2
!
   vthc_jm2 = vth_im2jm2*axt2 + vth_im1jm2*bxt2 + vth_ijm2*cxt2
   vthc_jm1 = vth_im2jm1*axt2 + vth_im1jm1*bxt2 + vth_ijm1*cxt2
   vthc_j   = vth_im2j*axt2   + vth_im1j*bxt2   + vth_ij*cxt2
!
   slinec_jm2 = sline_im2jm2*axt + sline_im1jm2*bxt + sline_ijm2*cxt
   slinec_jm1 = sline_im2jm1*axt + sline_im1jm1*bxt + sline_ijm1*cxt
   slinec_j   = sline_im2j*axt   + sline_im1j*bxt   + sline_ij*cxt
!
!------------------ensure monotonicity on level j-----------------------
!
!velocity components (no monotonicity required)
   velx_j = ax*velx_im1j + bx*velxc_j + cx*velx_ij
   vely_j = ax*vely_im1j + bx*velyc_j + cx*vely_ij
   velz_j = ax*velz_im1j + bx*velzc_j + cx*velz_ij
!
!opacity and thermal velocity (no interpolation coefficients required)
   call pointcr1d_mbez(vthc_j, vth_im1j, vth_ij)
   call pointcr1d_mbez(opalbarc_j, opalbar_im1j, opalbar_ij)
   vth_j = ax*vth_im1j + bx*vthc_j + cx*vth_ij
   opalbar_j = ax*opalbar_im1j + bx*opalbarc_j + cx*opalbar_ij
!
!line source function
   call coeffcr1d_mbez(slinec_j, sline_im1j, sline_ij, axt, bxt, cxt, axtj_sline, bxtj_sline, cxtj_sline)
   axj_sline = axtj_sline*bx
   bxj_sline = bxtj_sline*bx + ax
   cxj_sline = cxtj_sline*bx + cx
   sline_j = axj_sline*sline_im2j + bxj_sline*sline_im1j + cxj_sline*sline_ij
!
!---------------ensure monotonicity on level j-1------------------------
!
!velocity components (no monotonicity required)
   velx_jm1 = ax*velx_im1jm1 + bx*velxc_jm1 + cx*velx_ijm1
   vely_jm1 = ax*vely_im1jm1 + bx*velyc_jm1 + cx*vely_ijm1
   velz_jm1 = ax*velz_im1jm1 + bx*velzc_jm1 + cx*velz_ijm1
!
!opacity and thermal velocity (no interpolation coefficients required)
   call pointcr1d_mbez(vthc_jm1, vth_im1jm1, vth_ijm1)
   call pointcr1d_mbez(opalbarc_jm1, opalbar_im1jm1, opalbar_ijm1)
   vth_jm1 = ax*vth_im1jm1 + bx*vthc_jm1 + cx*vth_ijm1
   opalbar_jm1 = ax*opalbar_im1jm1 + bx*opalbarc_jm1 + cx*opalbar_ijm1
!
!line source function
   call coeffcr1d_mbez(slinec_jm1, sline_im1jm1, sline_ijm1, axt, bxt, cxt, axtjm1_sline, bxtjm1_sline, cxtjm1_sline)
   axjm1_sline = axtjm1_sline*bx
   bxjm1_sline = bxtjm1_sline*bx + ax
   cxjm1_sline = cxtjm1_sline*bx + cx
   sline_jm1 = axjm1_sline*sline_im2jm1 + bxjm1_sline*sline_im1jm1 + cxjm1_sline*sline_ijm1
!
!---------------ensure monotonicity on level j-2------------------------
!
!velocity components (no monotonicity required)
   velx_jm2 = ax*velx_im1jm2 + bx*velxc_jm2 + cx*velx_ijm2
   vely_jm2 = ax*vely_im1jm2 + bx*velyc_jm2 + cx*vely_ijm2
   velz_jm2 = ax*velz_im1jm2 + bx*velzc_jm2 + cx*velz_ijm2
!
!opacity and thermal velocity (no interpolation coefficients required)
   call pointcr1d_mbez(vthc_jm2, vth_im1jm2, vth_ijm2)
   call pointcr1d_mbez(opalbarc_jm2, opalbar_im1jm2, opalbar_ijm2)
   vth_jm2 = ax*vth_im1jm2 + bx*vthc_jm2 + cx*vth_ijm2
   opalbar_jm2 = ax*opalbar_im1jm2 + bx*opalbarc_jm2 + cx*opalbar_ijm2
!
!line source function and intensity
   call coeffcr1d_mbez(slinec_jm2, sline_im1jm2, sline_ijm2, axt, bxt, cxt, axtjm2_sline, bxtjm2_sline, cxtjm2_sline)
   axjm2_sline = axtjm2_sline*bx
   bxjm2_sline = bxtjm2_sline*bx + ax
   cxjm2_sline = cxtjm2_sline*bx + cx
   sline_jm2 = axjm2_sline*sline_im2jm2 + bxjm2_sline*sline_im1jm2 + cxjm2_sline*sline_ijm2
!
!------------calculate control point for interpolation along y----------
!
!derivative weights for source function
   fac=max(wp_interp2d,dyjm1/dy)
!derivative weights for velocity components and opacity
   fac2=dyjm1/dy
!
   ayt = (fac-1.d0)*dyj/2.d0/dyjm1
   byt = ((2.d0-fac)*dyjm1 + (1.d0-fac)*dyj)/2.d0/dyjm1
   cyt = fac/2.d0
   ayt2 = (fac2-1.d0)*dyj/2.d0/dyjm1
   byt2 = ((2.d0-fac2)*dyjm1 + (1.d0-fac2)*dyj)/2.d0/dyjm1
   cyt2 = fac2/2.d0
!
   ay = (1.d0-ty)**2
   by = 2.d0*ty*(1.d0-ty)
   cy = ty**2
!
   velx_c = velx_jm2*ayt2 + velx_jm1*byt2 + velx_j*cyt2
   vely_c = vely_jm2*ayt2 + vely_jm1*byt2 + vely_j*cyt2
   velz_c = velz_jm2*ayt2 + velz_jm1*byt2 + velz_j*cyt2
   vth_c = vth_jm2*ayt2 + vth_jm1*byt2 + vth_j*cyt2
   opalbar_c = opalbar_jm2*ayt2 + opalbar_jm1*byt2 + opalbar_j*cyt2
!
   sline_c = sline_jm2*ayt + sline_jm1*byt + sline_j*cyt
!
!------------------------ensure monotonicity----------------------------
!
!velocity components (no monotonicity required)
   velx_p = ay*velx_jm1 + by*velx_c + cy*velx_j
   vely_p = ay*vely_jm1 + by*vely_c + cy*vely_j
   velz_p = ay*velz_jm1 + by*velz_c + cy*velz_j
!
!opacity and thermal velocity (no interpolation coefficients required)
   call pointcr1d_mbez(vth_c, vth_jm1, vth_j)
   call pointcr1d_mbez(opalbar_c, opalbar_jm1, opalbar_j)
   vth_p = ay*vth_jm1 + by*vth_c + cy*vth_j
   opalbar_p = ay*opalbar_jm1 + by*opalbar_c + cy*opalbar_j
!
!line source function and intensity
   call coeffcr1d_mbez(sline_c, sline_jm1, sline_j, ayt, byt, cyt, ayt_sline, byt_sline, cyt_sline)
   ay_sline = ayt_sline*by
   by_sline = byt_sline*by + ay
   cy_sline = cyt_sline*by + cy
   sline_p = ay_sline*sline_jm2 + by_sline*sline_jm1 + cy_sline*sline_j
!
   a_sline = ay_sline*axjm2_sline
   b_sline = ay_sline*bxjm2_sline
   c_sline = ay_sline*cxjm2_sline
   d_sline = by_sline*axjm1_sline
   e_sline = by_sline*bxjm1_sline
   f_sline = by_sline*cxjm1_sline
   g_sline = cy_sline*axj_sline
   h_sline = cy_sline*bxj_sline
   i_sline = cy_sline*cxj_sline
!
!write(*,'(11es20.8)') x_p, y_p, z_p, opalbar_p, opalbar_p2, velx_p, velx_p2, vely_p, vely_p2, velz_p, velz_p2
!
   return


!
!
end subroutine coeff2d_lined
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine coeff2d_lineu_lin(opalbar_im1jm1, opalbar_ijm1, opalbar_im1j, opalbar_ij, &
   velx_im1jm1, velx_ijm1, velx_im1j, velx_ij, &
   vely_im1jm1, vely_ijm1, vely_im1j, vely_ij, &
   velz_im1jm1, velz_ijm1, velz_im1j, velz_ij, &
   vth_im1jm1, vth_ijm1, vth_im1j, vth_ij, &
   sline_im1jm1, sline_ijm1, sline_im1j, sline_ij, &
   int_im1jm1, int_ijm1, int_im1j, int_ij, &
   x_im1, x_i, y_jm1, y_j, x_p, y_p, &
   a_sline, b_sline, c_sline, d_sline, &
   a_inten, b_inten, c_inten, d_inten, &
   opalbar_p, velx_p, vely_p, velz_p, vth_p, sline_p, int_p)
!
!         interpolates opaciy, velocity components, thermal velocity,
!                   line source function and intensity
!               values given on a 2d grid onto point x_p, y_p
!
!on input (f_* stands for opalbar_*, sline_* and int_*, respectivly):
!
! y_j      f_im1j--------------f_ij
!  |          |                  |
!  |          |                  |
!  |          |         x        |
!  |          |     (x_p,y_p)    |
!  |          |                  |
!y_jm1    f_im1jm1------------f_ijm1
!  |
!  --------x_im2-----------x_im1---------------x_i
!
!        x_p, y_p: coordinates of point onto which shall be interpolated
!
!on output:
!   1. interpolation coefficients for source function and intensity
!      (required for ALO calculations):
!         a_sline, b_sline, c_sline, d_sline
!         a_inten, b_inten, c_inten, d_inten
!
!      such that:
!         f_p = a*f_im1jm1 + b*f_ijm1 + c*f_im1j + d*f_ij
!
!   2. interpolated values at point p: opalbar_p, sline_p, velx_p, vely_p, velz_p,
!                                      vth_p, int_p
!
   use prog_type
!
   implicit none
!
! ... arguments
   real(dp), intent(in) :: opalbar_im1jm1, opalbar_ijm1, opalbar_im1j, opalbar_ij, &
      velx_im1jm1, velx_ijm1, velx_im1j, velx_ij, &
      vely_im1jm1, vely_ijm1, vely_im1j, vely_ij, &
      velz_im1jm1, velz_ijm1, velz_im1j, velz_ij, &
      vth_im1jm1, vth_ijm1, vth_im1j, vth_ij, &
      sline_im1jm1, sline_ijm1, sline_im1j, sline_ij, &
      int_im1jm1, int_ijm1, int_im1j, int_ij, &
      x_im1, x_i, y_jm1, y_j, x_p, y_p
   real(dp), intent(out) :: a_sline, b_sline, c_sline, d_sline, &
      a_inten, b_inten, c_inten, d_inten, &
      opalbar_p, velx_p, vely_p, velz_p, vth_p, int_p, sline_p
!
! ... local scalars
   real(dp) :: tx, ty, rdxdy
!
!-------------------------bilinear interpolation------------------------
!
!define deltax, deltay-ratios
   tx = (x_p-x_im1)/(x_i-x_im1)
   ty = (y_p-y_jm1)/(y_j-y_jm1)
!
   rdxdy=tx*ty
!
   a_sline=1.d0-tx-ty+rdxdy
   b_sline=tx-rdxdy
   c_sline=ty-rdxdy
   d_sline=rdxdy
!
   a_inten=a_sline
   b_inten=b_sline
   c_inten=c_sline
   d_inten=d_sline
!
   opalbar_p = a_sline*opalbar_im1jm1 + b_sline*opalbar_ijm1 + c_sline*opalbar_im1j + d_sline*opalbar_ij
   velx_p = a_sline*velx_im1jm1 + b_sline*velx_ijm1 + c_sline*velx_im1j + d_sline*velx_ij
   vely_p = a_sline*vely_im1jm1 + b_sline*vely_ijm1 + c_sline*vely_im1j + d_sline*vely_ij
   velz_p = a_sline*velz_im1jm1 + b_sline*velz_ijm1 + c_sline*velz_im1j + d_sline*velz_ij
   vth_p = a_sline*vth_im1jm1 + b_sline*vth_ijm1 + c_sline*vth_im1j + d_sline*vth_ij
!
   sline_p = a_sline*sline_im1jm1 + b_sline*sline_ijm1 + c_sline*sline_im1j + d_sline*sline_ij
   int_p = a_sline*int_im1jm1 + b_sline*int_ijm1 + c_sline*int_im1j + d_sline*int_ij
!
   return
!
end subroutine coeff2d_lineu_lin
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine coeff2d_lined_lin(opalbar_im1jm1, opalbar_ijm1, opalbar_im1j, opalbar_ij, &
   velx_im1jm1, velx_ijm1, velx_im1j, velx_ij, &
   vely_im1jm1, vely_ijm1, vely_im1j, vely_ij, &
   velz_im1jm1, velz_ijm1, velz_im1j, velz_ij, &
   vth_im1jm1, vth_ijm1, vth_im1j, vth_ij, &
   sline_im1jm1, sline_ijm1, sline_im1j, sline_ij, &
   x_im1, x_i, y_jm1, y_j, x_p, y_p, &
   a_sline, b_sline, c_sline, d_sline, &
   opalbar_p, velx_p, vely_p, velz_p, vth_p, sline_p)
!
!         interpolates opaciy, velocity components, thermal velocity,
!                           line source function
!               values given on a 2d grid onto point x_p, y_p
!
!on input (f_* stands for opalbar_*, sline_* and respectivly):
!
! y_j      f_im1j--------------f_ij
!  |          |                  |
!  |          |                  |
!  |          |         x        |
!  |          |     (x_p,y_p)    |
!  |          |                  |
!y_jm1    f_im1jm1------------f_ijm1
!  |
!  --------x_im2-----------x_im1---------------x_i
!
!        x_p, y_p: coordinates of point onto which shall be interpolated
!
!on output:
!   1. interpolation coefficients for source function and intensity
!      (required for ALO calculations):
!         a_sline, b_sline, c_sline, d_sline
!
!      such that:
!         f_p = a*f_im1jm1 + b*f_ijm1 + c*f_im1j + d*f_ij
!
!   2. interpolated values at point p: opalbar_p, sline_p, velx_p, vely_p, velz_p,
!                                      vth_p
!
   use prog_type
!
   implicit none
!
! ... arguments
   real(dp), intent(in) :: opalbar_im1jm1, opalbar_ijm1, opalbar_im1j, opalbar_ij, &
      velx_im1jm1, velx_ijm1, velx_im1j, velx_ij, &
      vely_im1jm1, vely_ijm1, vely_im1j, vely_ij, &
      velz_im1jm1, velz_ijm1, velz_im1j, velz_ij, &
      vth_im1jm1, vth_ijm1, vth_im1j, vth_ij, &
      sline_im1jm1, sline_ijm1, sline_im1j, sline_ij, &
      x_im1, x_i, y_jm1, y_j, x_p, y_p
   real(dp), intent(out) :: a_sline, b_sline, c_sline, d_sline, &
      opalbar_p, velx_p, vely_p, velz_p, vth_p, sline_p
!
! ... local scalars
   real(dp) :: tx, ty, rdxdy
!
!-------------------------bilinear interpolation------------------------
!
!define deltax, deltay-ratios
   tx = (x_p-x_im1)/(x_i-x_im1)
   ty = (y_p-y_jm1)/(y_j-y_jm1)
!
   rdxdy=tx*ty
!
   a_sline=1.d0-tx-ty+rdxdy
   b_sline=tx-rdxdy
   c_sline=ty-rdxdy
   d_sline=rdxdy
!
   opalbar_p = a_sline*opalbar_im1jm1 + b_sline*opalbar_ijm1 + c_sline*opalbar_im1j + d_sline*opalbar_ij
   velx_p = a_sline*velx_im1jm1 + b_sline*velx_ijm1 + c_sline*velx_im1j + d_sline*velx_ij
   vely_p = a_sline*vely_im1jm1 + b_sline*vely_ijm1 + c_sline*vely_im1j + d_sline*vely_ij
   velz_p = a_sline*velz_im1jm1 + b_sline*velz_ijm1 + c_sline*velz_im1j + d_sline*velz_ij
   vth_p = a_sline*vth_im1jm1 + b_sline*vth_ijm1 + c_sline*vth_im1j + d_sline*vth_ij
!
   sline_p = a_sline*sline_im1jm1 + b_sline*sline_ijm1 + c_sline*sline_im1j + d_sline*sline_ij
!
   return
!
end subroutine coeff2d_lined_lin
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine coeff3d_lineu_lin(iim1, ii, jjm1, jj, kkm1, kk, x_p, y_p, z_p, &
   c14, c15, c17, c18, c23, c24, c26, c27, &
   opalbar_p, velx_p, vely_p, velz_p, vth_p, sline_p)
!
!         interpolates opacity, velocity components, thermal velocity,
!                     and line source function
!           given on a 3d grid onto point x_p, y_p, z_p
!
!on input:
!
!               f_im1jk--------------f_ijk
!                 /|                  / |
!                / |                 /  |
!               /  |                /   |
!              /   |               /    |
!             /    |              /     |
!            / f_im1jkm1---------/---f_ijkm1
!           /      /            /      /
!       f_im1jm1k------------f_ijm1k  /
!           |    /              |    /
!           |   /               |   /
!           |  /                |  /
!           | /                 | /
!           |/                  |/
!       f_im1jm1km1---------f_ijm1km1
!
!        x_p, y_p, z_p: coordinates of point onto which shall be interpolated
!
!on output:
!   1. interpolation coefficients for source function and intensity
!      (required for ALO calculations):
!         c14, c15, c17, c18,
!         c23, c24, c25, c26, c27
!
!   such that
!      f_p = c14*f_im1jm1km1 + c15*f_ijm1km1 +
!            c17*f_im1jkm1   + c18*f_ijkm1 +
!            c23*f_im1jm1k + c24*f_ijm1k +
!            c26*f_im1jk   + c27*f_ijk
!
!   2. interpolated values at point p: opalbar_p, sline_p
!
   use prog_type
   use dime3d, only: x, y, z, opalbar3d, velx3d, vely3d, velz3d, vth3d, sline3d
!***debug start
   use fund_const, only: xmsu, pi
   use params_input, only: vmin, vmax, yhe, hei, kcont, kappa0, alpha, xmloss, beta, vth_fiducial
   use params_stellar, only: sr
   use mod_opacities, only: opalbar_model_hamann, opalbar_model_kline
!***debug end
!
   implicit none
!
! ... argments
   integer(i4b), intent(in) :: iim1, ii, jjm1, jj, kkm1, kk
   real(dp), intent(in) :: x_p, y_p, z_p
   real(dp), intent(out) :: c14, c15, c17, c18, c23, c24, c26, c27, &
      opalbar_p, velx_p, vely_p, velz_p, vth_p, sline_p
!
! ... local scalars
   real(dp) :: tx, ty, tz, ax, ay, az

!
!***debug start
   real(dp) :: rad, vinf, xmloss_cgs, bconst, velr, rho, c1, c2, ne, opalbar_p2, sline_p2, int_p2, &
      velx_p2, vely_p2, velz_p2
!***debug end
!
!****************************debug start********************************
!----------------analytic expressions for velocity and opacity----------
!
!velocity
!vinf=vmax*1.d5
!bconst = 1.d0-(vmin/vmax)**(1.d0/beta)
!rad=sqrt(x_p**2+y_p**2+z_p**2)
!velr = vinf*(1.d0-bconst/rad)**beta
!velx_p=velr*x_p/rad / vth_fiducial
!vely_p=velr*y_p/rad / vth_fiducial
!velz_p=velr*z_p/rad / vth_fiducial
!
!xmloss_cgs=xmloss*xmsu/(365.25d0*24.d0*60.d0*60.d0)
!rho = xmloss_cgs/4.d0/pi/velr/(rad*sr)**2

!!opalbar_p2 = opalbar_model_kline(yhe, hei, rho, kline)*sr
!opalbar_p = opalbar_model_hamann(sr, vinf, xmloss_cgs, kappa0, alpha, vth_fiducial, rad*sr, rho) * sr
!
!****************************debug end**********************************
!
!--------------------------trilinear interpolation----------------------
!
!define deltax, deltay, deltaz
   tx = (x_p-x(iim1))/(x(ii)-x(iim1))
   ty = (y_p-y(jjm1))/(y(jj)-y(jjm1))
   tz = (z_p-z(kkm1))/(z(kk)-z(kkm1))
!
   ax=1.d0-tx
   ay=1.d0-ty
   az=1.d0-tz
!
   c14 = ax*ay*az
   c15 = tx*ay*az
   c17 = ax*ty*az
   c18 = tx*ty*az
   c23 = ax*ay*tz
   c24 = tx*ay*tz
   c26 = ax*ty*tz
   c27 = tx*ty*tz

   opalbar_p = c14*opalbar3d(iim1,jjm1,kkm1) + c15*opalbar3d(ii,jjm1,kkm1) + c17*opalbar3d(iim1,jj,kkm1) + c18*opalbar3d(ii,jj,kkm1) + &
      c23*opalbar3d(iim1,jjm1,kk)   + c24*opalbar3d(ii,jjm1,kk)   + c26*opalbar3d(iim1,jj,kk)   + c27*opalbar3d(ii,jj,kk)
   velx_p = c14*velx3d(iim1,jjm1,kkm1) + c15*velx3d(ii,jjm1,kkm1) + c17*velx3d(iim1,jj,kkm1) + c18*velx3d(ii,jj,kkm1) + &
      c23*velx3d(iim1,jjm1,kk)   + c24*velx3d(ii,jjm1,kk)   + c26*velx3d(iim1,jj,kk)   + c27*velx3d(ii,jj,kk)
   vely_p = c14*vely3d(iim1,jjm1,kkm1) + c15*vely3d(ii,jjm1,kkm1) + c17*vely3d(iim1,jj,kkm1) + c18*vely3d(ii,jj,kkm1) + &
      c23*vely3d(iim1,jjm1,kk)   + c24*vely3d(ii,jjm1,kk)   + c26*vely3d(iim1,jj,kk)   + c27*vely3d(ii,jj,kk)
   velz_p = c14*velz3d(iim1,jjm1,kkm1) + c15*velz3d(ii,jjm1,kkm1) + c17*velz3d(iim1,jj,kkm1) + c18*velz3d(ii,jj,kkm1) + &
      c23*velz3d(iim1,jjm1,kk)   + c24*velz3d(ii,jjm1,kk)   + c26*velz3d(iim1,jj,kk)   + c27*velz3d(ii,jj,kk)
   vth_p = c14*vth3d(iim1,jjm1,kkm1) + c15*vth3d(ii,jjm1,kkm1) + c17*vth3d(iim1,jj,kkm1) + c18*vth3d(ii,jj,kkm1) + &
      c23*vth3d(iim1,jjm1,kk)   + c24*vth3d(ii,jjm1,kk)   + c26*vth3d(iim1,jj,kk)   + c27*vth3d(ii,jj,kk)
   sline_p = c14*sline3d(iim1,jjm1,kkm1) + c15*sline3d(ii,jjm1,kkm1) + c17*sline3d(iim1,jj,kkm1) + c18*sline3d(ii,jj,kkm1) + &
      c23*sline3d(iim1,jjm1,kk)   + c24*sline3d(ii,jjm1,kk)   + c26*sline3d(iim1,jj,kk)   + c27*sline3d(ii,jj,kk)


!write(*,'(12es20.8)') x_p, y_p, z_p, rad, opalbar_p2, opalbar_p, velx_p2, velx_p, vely_p2, vely_p, velz_p2, velz_p
!write(*,*)
!
   return
!
!
end subroutine coeff3d_lineu_lin
!
!***********************************************************************
!***********************************************************************
!
!               LINE + CONTINUUM TRANSFER ROUTINES
!
!***********************************************************************
!***********************************************************************
!
subroutine fsc_linec3d_lin(oindx,xobsindx)
!
!-----------------------------------------------------------------------
!--short characteristics for line + continuum radiative transfer in 3d--
!-----------calculating intensties for given mu,phi specified-----------
!-----------by input oindx, and for xobs specified by xobsindx----------
!------------------only linear interpolations are used------------------
!-----------------------------------------------------------------------
!
   use prog_type
   use fund_const
   use dime3d, only: ndxmax, ndzmax, ndymax, x, y, z
   use dime3d, only: int3d, opac3d, opalbar3d, scont3d, sline3d, velx3d, vely3d, velz3d, vth3d, &
      imask3d, imask_totreg3d, aloline_on_nn3d, aloline_nn3d_tmp, mintbar3d_tmp, &
      normalization3d_tmp
   use angles, only: n_x, n_y, n_z, weight_omega, q_alo
   use freq, only: nodes_xobs, weight_xobs, nxobs, deltax, xcmf_min, xcmf_max
   use params_input, only: vth_fiducial
   use bcondition, only: xic1
   use mod_debug, only: indxx, indxy, indxz
   use params_stellar, only: smajorax_a, smajorax_b, smajorax_c
!
   implicit none
!
! ... arguments
   integer(i4b), intent(in) :: oindx, xobsindx
!
! ... local scalars
   integer(i4b) :: i, j, k
   integer(i4b) :: alpha, beta, gamma
   integer(i4b) :: startx, starty, startz, endx, endy, endz
   integer(i4b) :: iim1, ii, jjm1, jj, kkm1, kk
   real(dp) :: nn_x, nn_y, nn_z, mueff, xobs, wall, wall_global, phinorm
   real(dp) :: dels_xyu, dels_xzu, dels_yzu, dels_u, dels_r
   real(dp) :: opac_p, opalbar_p, velx_p, vely_p, velz_p, vel_p, vth_p, scont_p, sline_p
   real(dp) :: x_u, y_u, z_u, int_u, opac_u, opalbar_u, velx_u, vely_u, velz_u, vel_u, vth_u, scont_u, sline_u
   real(dp) :: abs_sc, int_sc, contr_sc
   real(dp) :: alo_u, alo_p
   integer :: q14, q15, q17, q18, q23, q24, q26, q27
   real(dp) :: c14_slineu, c15_slineu, c17_slineu, c18_slineu, &
      c23_slineu, c24_slineu, c26_slineu, c27_slineu, &
      c14_intu, c15_intu, c17_intu, c18_intu, &
      c23_intu, c24_intu, c26_intu
!
! *** for debugging
   real(dp) :: alo_u2, alo_p2, int_sc2, contr_sc2, abs_sc2
!
! ... local arrays
!
! ... local functions
   real(dp) :: dist_sphere, dist_ellipsoid, calc_icore_gdark
!
! ... local characters
!character(len=50) :: enter
!
! ... local logicals
!
!directions
   nn_x=n_x(oindx)
   nn_y=n_y(oindx)
   nn_z=n_z(oindx)
!
!frequency
   xobs=nodes_xobs(xobsindx)
!
!angular and frequency  integration weights
   wall_global=weight_omega(oindx)*weight_xobs(xobsindx)
!
!indices for nearest neighbour alo
   q14=q_alo(oindx,14)
   q15=q_alo(oindx,15)
   q17=q_alo(oindx,17)
   q18=q_alo(oindx,18)
   q23=q_alo(oindx,23)
   q24=q_alo(oindx,24)
   q26=q_alo(oindx,26)
   q27=q_alo(oindx,27)
!
!-----------------------------------------------------------------------
!
!set directional index-parameter (phantom points are excluded from calculation)
!
!index parameter:
!         if n_x >= 0                 if n_x < 0
!                startx = 2                  startx = ndxmax-1
!                endx = ndxmax-1             endx = 2
!                alpha=  1                   alpha=-1
!
!         if n_y >= 0                 if n_y < 0
!                starty = 2                  starty = ndymax-1
!                endy = ndymax-1             endy = 2
!                beta =  1                   beta =-1
!
!         if n_z >= 0                 if n_z < 0
!                startz = 2                  startz = ndzmax-1
!                endz = ndzmax-1             endz = 2
!                gamma = 1                   gamma = -1
!
   if(nn_x.gt.0.d0) then
      startx = 3
      endx = ndxmax-2
      alpha=  1
   elseif(nn_x.lt.0.d0) then
      startx = ndxmax-2
      endx = 3
      alpha=-1
   else
      stop 'error in fsc_linec3d_lin: n_x = 0 not allowed'
   endif
!
   if(nn_y.gt.0.d0) then
      starty = 3
      endy = ndymax-2
      beta =  1
   elseif(nn_y.lt.0.d0) then
      starty = ndymax-2
      endy = 3
      beta =-1
   else
      stop 'error in fsc_linec3d_lin: n_y = 0 not allowed'
   endif
!
   if(nn_z.gt.0.d0) then
      startz = 3
      endz = ndzmax-2
      gamma=  1
   elseif(nn_z.lt.0.d0) then
      startz = ndzmax-2
      endz = 3
      gamma=-1
   else
      stop 'error in fsc_linec3d_lin: n_z = 0 not allowed'
   endif
!
!--------------------reset the intensities and alo----------------------
!
   call set_boundary3d(nn_x, nn_y, nn_z)
!
   aloline_on_nn3d=0.d0
!
!-----------------------------------------------------------------------
!
   do k=startz, endz, gamma
      do j=starty, endy, beta
         do i=startx, endx, alpha
!
            select case(imask3d(i,j,k))
!
!************************standard radiative transfer********************
!
             case(1)
!
               iim1=i-alpha
               jjm1=j-beta
               kkm1=k-gamma
!
!calculate distance of upwind point to (previous) xy-plane, xz-plane and yz-plane
               dels_xyu=(z(k)-z(kkm1))/nn_z
               dels_xzu=(y(j)-y(jjm1))/nn_y
               dels_yzu=(x(i)-x(iim1))/nn_x
               dels_u=min(dels_xyu,dels_xzu,dels_yzu)
!
!----------------------------local point--------------------------------
!
               sline_p=sline3d(i,j,k)
               scont_p=scont3d(i,j,k)
               opalbar_p=opalbar3d(i,j,k)
               opac_p=opac3d(i,j,k)
               velx_p=velx3d(i,j,k)
               vely_p=vely3d(i,j,k)
               velz_p=velz3d(i,j,k)
               vel_p = nn_x*velx_p + nn_y*vely_p + nn_z*velz_p
               vth_p = vth3d(i,j,k)
!
!----------------------------upwind point-------------------------------
!
               if(dels_xyu.eq.dels_u) then
!intersection with x-y plane on level k-gamma
                  x_u = x(i) - dels_u*nn_x
                  y_u = y(j) - dels_u*nn_y
                  z_u = z(kkm1)
!
                  call coeff2d_linecu_lin(opac3d(iim1,jjm1,kkm1), opac3d(i,jjm1,kkm1), opac3d(iim1,j,kkm1), opac3d(i,j,kkm1), &
                     opalbar3d(iim1,jjm1,kkm1), opalbar3d(i,jjm1,kkm1), opalbar3d(iim1,j,kkm1), opalbar3d(i,j,kkm1), &
                     velx3d(iim1,jjm1,kkm1), velx3d(i,jjm1,kkm1), velx3d(iim1,j,kkm1), velx3d(i,j,kkm1), &
                     vely3d(iim1,jjm1,kkm1), vely3d(i,jjm1,kkm1), vely3d(iim1,j,kkm1), vely3d(i,j,kkm1), &
                     velz3d(iim1,jjm1,kkm1), velz3d(i,jjm1,kkm1), velz3d(iim1,j,kkm1), velz3d(i,j,kkm1), &
                     vth3d(iim1,jjm1,kkm1), vth3d(i,jjm1,kkm1), vth3d(iim1,j,kkm1), vth3d(i,j,kkm1), &
                     scont3d(iim1,jjm1,kkm1), scont3d(i,jjm1,kkm1), scont3d(iim1,j,kkm1), scont3d(i,j,kkm1), &
                     sline3d(iim1,jjm1,kkm1), sline3d(i,jjm1,kkm1), sline3d(iim1,j,kkm1), sline3d(i,j,kkm1), &
                     int3d(iim1,jjm1,kkm1), int3d(i,jjm1,kkm1), int3d(iim1,j,kkm1), int3d(i,j,kkm1), &
                     x(iim1), x(i), y(jjm1), y(j), x_u, y_u, &
                     c14_slineu, c15_slineu, c17_slineu, c18_slineu, &
                     c14_intu, c15_intu, c17_intu, c18_intu, opac_u, opalbar_u, velx_u, vely_u, velz_u, vth_u, scont_u, sline_u, int_u)
                  c23_slineu=0.d0
                  c24_slineu=0.d0
                  c26_slineu=0.d0
                  c23_intu=0.d0
                  c24_intu=0.d0
                  c26_intu=0.d0
!
               elseif(dels_xzu.eq.dels_u) then
!intersection with x-z plane on level j-beta
                  x_u = x(i) - dels_u*nn_x
                  y_u = y(jjm1)
                  z_u = z(k) - dels_u*nn_z
!
                  call coeff2d_linecu_lin(opac3d(iim1,jjm1,kkm1), opac3d(i,jjm1,kkm1), opac3d(iim1,jjm1,k), opac3d(i,jjm1,k), &
                     opalbar3d(iim1,jjm1,kkm1), opalbar3d(i,jjm1,kkm1), opalbar3d(iim1,jjm1,k), opalbar3d(i,jjm1,k), &
                     velx3d(iim1,jjm1,kkm1), velx3d(i,jjm1,kkm1), velx3d(iim1,jjm1,k), velx3d(i,jjm1,k), &
                     vely3d(iim1,jjm1,kkm1), vely3d(i,jjm1,kkm1), vely3d(iim1,jjm1,k), vely3d(i,jjm1,k), &
                     velz3d(iim1,jjm1,kkm1), velz3d(i,jjm1,kkm1), velz3d(iim1,jjm1,k), velz3d(i,jjm1,k), &
                     vth3d(iim1,jjm1,kkm1), vth3d(i,jjm1,kkm1), vth3d(iim1,jjm1,k), vth3d(i,jjm1,k), &
                     scont3d(iim1,jjm1,kkm1), scont3d(i,jjm1,kkm1), scont3d(iim1,jjm1,k), scont3d(i,jjm1,k), &
                     sline3d(iim1,jjm1,kkm1), sline3d(i,jjm1,kkm1), sline3d(iim1,jjm1,k), sline3d(i,jjm1,k), &
                     int3d(iim1,jjm1,kkm1), int3d(i,jjm1,kkm1), int3d(iim1,jjm1,k), int3d(i,jjm1,k), &
                     x(iim1), x(i), z(kkm1), z(k), x_u, z_u, &
                     c14_slineu, c15_slineu, c23_slineu, c24_slineu, &
                     c14_intu, c15_intu, c23_intu, c24_intu, &
                     opac_u, opalbar_u, velx_u, vely_u, velz_u, vth_u, scont_u, sline_u, int_u)
                  c17_slineu=0.d0
                  c18_slineu=0.d0
                  c26_slineu=0.d0
                  c17_intu=0.d0
                  c18_intu=0.d0
                  c26_intu=0.d0
!
               elseif(dels_yzu.eq.dels_u) then
!intersection with y-z plane on level i-alpha
                  x_u = x(iim1)
                  y_u = y(j) - dels_u*nn_y
                  z_u = z(k) - dels_u*nn_z
!
                  call coeff2d_linecu_lin(opac3d(iim1,jjm1,kkm1), opac3d(iim1,j,kkm1), opac3d(iim1,jjm1,k), opac3d(iim1,j,k), &
                     opalbar3d(iim1,jjm1,kkm1), opalbar3d(iim1,j,kkm1), opalbar3d(iim1,jjm1,k), opalbar3d(iim1,j,k), &
                     velx3d(iim1,jjm1,kkm1), velx3d(iim1,j,kkm1), velx3d(iim1,jjm1,k), velx3d(iim1,j,k), &
                     vely3d(iim1,jjm1,kkm1), vely3d(iim1,j,kkm1), vely3d(iim1,jjm1,k), vely3d(iim1,j,k), &
                     velz3d(iim1,jjm1,kkm1), velz3d(iim1,j,kkm1), velz3d(iim1,jjm1,k), velz3d(iim1,j,k), &
                     vth3d(iim1,jjm1,kkm1), vth3d(iim1,j,kkm1), vth3d(iim1,jjm1,k), vth3d(iim1,j,k), &
                     scont3d(iim1,jjm1,kkm1), scont3d(iim1,j,kkm1), scont3d(iim1,jjm1,k), scont3d(iim1,j,k), &
                     sline3d(iim1,jjm1,kkm1), sline3d(iim1,j,kkm1), sline3d(iim1,jjm1,k), sline3d(iim1,j,k), &
                     int3d(iim1,jjm1,kkm1), int3d(iim1,j,kkm1), int3d(iim1,jjm1,k), int3d(iim1,j,k), &
                     y(jjm1), y(j), z(kkm1), z(k), y_u, z_u, &
                     c14_slineu, c17_slineu, c23_slineu, c26_slineu, &
                     c14_intu, c17_intu, c23_intu, c26_intu, &
                     opac_u, opalbar_u, velx_u, vely_u, velz_u, vth_u, scont_u, sline_u, int_u)
                  c15_slineu=0.d0
                  c18_slineu=0.d0
                  c24_slineu=0.d0
                  c15_intu=0.d0
                  c18_intu=0.d0
                  c24_intu=0.d0
!
               else
                  write(*,'(4es20.8,2l4)') dels_u, dels_xzu, dels_xyu, dels_yzu, dels_u.eq.dels_xzu, dels_xyu.eq.dels_u
                  stop 'error in fsc_linec3d_lin: invalid dels_u'
               endif
!
!--------------------------------radiative transfer---------------------
!
               vel_u = velx_u*nn_x + vely_u*nn_y + velz_u*nn_z
               call fsc_linec_lin(xobs, vth_fiducial, deltax, xcmf_max, xcmf_min, int_u, opac_u, opac_p, opalbar_u, opalbar_p, scont_u, scont_p, &
                  sline_u, sline_p, vel_u, vel_p, vth_u, vth_p, dels_u, abs_sc, contr_sc, int_sc, alo_u, alo_p)
               int3d(i,j,k) = int_sc
!
               aloline_on_nn3d(i,j,k,q14) = imask_totreg3d(i,j,k) * alo_p
               aloline_on_nn3d(iim1,j,k,q15) = imask_totreg3d(iim1,j,k) * &
                  (alo_u*c26_slineu + abs_sc*c26_intu*aloline_on_nn3d(iim1,j,k,q14))

               aloline_on_nn3d(i,jjm1,k,q17) = imask_totreg3d(i,jjm1,k) * &
                  (alo_u*c24_slineu + abs_sc*c24_intu*aloline_on_nn3d(i,jjm1,k,q14))
               aloline_on_nn3d(iim1,jjm1,k,q18) = imask_totreg3d(iim1,jjm1,k) * &
                  (alo_u*c23_slineu + abs_sc*(c23_intu*aloline_on_nn3d(iim1,jjm1,k,q14) + &
                  c26_intu*aloline_on_nn3d(iim1,jjm1,k,q17)))
               aloline_on_nn3d(i,j,kkm1,q23) = imask_totreg3d(i,j,kkm1) * &
                  (alo_u*c18_slineu + abs_sc*c18_intu*aloline_on_nn3d(i,j,kkm1,q14))
               aloline_on_nn3d(iim1,j,kkm1,q24) = imask_totreg3d(iim1,j,kkm1) * &
                  (alo_u*c17_slineu + abs_sc*(c17_intu*aloline_on_nn3d(iim1,j,kkm1,q14) + &
                  c18_intu*aloline_on_nn3d(iim1,j,kkm1,q15) + &
                  c26_intu*aloline_on_nn3d(iim1,j,kkm1,q23)))
               aloline_on_nn3d(i,jjm1,kkm1,q26) = imask_totreg3d(i,jjm1,kkm1) * &
                  (alo_u*c15_slineu + abs_sc*(c15_intu*aloline_on_nn3d(i,jjm1,kkm1,q14) + &
                  c18_intu*aloline_on_nn3d(i,jjm1,kkm1,q17) + &
                  c24_intu*aloline_on_nn3d(i,jjm1,kkm1,q23)))
               aloline_on_nn3d(iim1,jjm1,kkm1,q27) = imask_totreg3d(iim1,jjm1,kkm1) * &
                  (alo_u*c14_slineu + abs_sc*(c14_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q14) + &
                  c15_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q15) + &
                  c17_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q17) + &
                  c18_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q18) + &
                  c23_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q23) + &
                  c24_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q24) + &
                  c26_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q26)))
!
!perform angular integration
               call calc_phinorm(vel_p, vth_p, vth_fiducial, xobs, phinorm)
               wall=wall_global*phinorm
               mintbar3d_tmp(i,j,k) = mintbar3d_tmp(i,j,k) + int3d(i,j,k)*wall
               normalization3d_tmp(i,j,k) = normalization3d_tmp(i,j,k) + wall
!
               aloline_nn3d_tmp(i,j,k,q14) = aloline_nn3d_tmp(i,j,k,q14) + wall*aloline_on_nn3d(i,j,k,q14)
               aloline_nn3d_tmp(iim1,j,k,q15) = aloline_nn3d_tmp(iim1,j,k,q15) + wall*aloline_on_nn3d(iim1,j,k,q15)
               aloline_nn3d_tmp(i,jjm1,k,q17) = aloline_nn3d_tmp(i,jjm1,k,q17) + wall*aloline_on_nn3d(i,jjm1,k,q17)
               aloline_nn3d_tmp(iim1,jjm1,k,q18) = aloline_nn3d_tmp(iim1,jjm1,k,q18) + wall*aloline_on_nn3d(iim1,jjm1,k,q18)
               aloline_nn3d_tmp(i,j,kkm1,q23) = aloline_nn3d_tmp(i,j,kkm1,q23) + wall*aloline_on_nn3d(i,j,kkm1,q23)
               aloline_nn3d_tmp(iim1,j,kkm1,q24) = aloline_nn3d_tmp(iim1,j,kkm1,q24) + wall*aloline_on_nn3d(iim1,j,kkm1,q24)
               aloline_nn3d_tmp(i,jjm1,kkm1,q26) = aloline_nn3d_tmp(i,jjm1,kkm1,q26) + wall*aloline_on_nn3d(i,jjm1,kkm1,q26)
               aloline_nn3d_tmp(iim1,jjm1,kkm1,q27) = aloline_nn3d_tmp(iim1,jjm1,kkm1,q27) + wall*aloline_on_nn3d(iim1,jjm1,kkm1,q27)
!
!***********special treatment: point is in vicinity of boundary*********
!
             case(2)
!
               iim1=i-alpha
               jjm1=j-beta
               kkm1=k-gamma
!
!calculate distance of upwind point to (previous) xy-plane, xz-plane and yz-plane
               dels_xyu=(z(k)-z(kkm1))/nn_z
               dels_xzu=(y(j)-y(jjm1))/nn_y
               dels_yzu=(x(i)-x(iim1))/nn_x
!
!calculate distance to the boundary (here: sphere)
!               dels_r=dist_sphere(nn_x,nn_y,nn_z,x(i),y(j),z(k))
               dels_r=dist_ellipsoid(nn_x,nn_y,nn_z,x(i),y(j),z(k),smajorax_a,smajorax_b,smajorax_c)
!
               dels_u=min(dels_xyu,dels_xzu,dels_yzu,dels_r)
!
!----------------------------local point--------------------------------
!
               sline_p=sline3d(i,j,k)
               scont_p=scont3d(i,j,k)
               opalbar_p=opalbar3d(i,j,k)
               opac_p=opac3d(i,j,k)
               velx_p=velx3d(i,j,k)
               vely_p=vely3d(i,j,k)
               velz_p=velz3d(i,j,k)
               vel_p = nn_x*velx_p + nn_y*vely_p + nn_z*velz_p
               vth_p = vth3d(i,j,k)
!
!----------------------------upwind point-------------------------------
!
               if(dels_r.eq.dels_u) then
                  x_u = x(i) - dels_u*nn_x
                  y_u = y(j) - dels_u*nn_y
                  z_u = z(k) - dels_u*nn_z
!
!for alo
!                  int_u=0.d0
!                  int_u=xic1
                  int_u=calc_icore_gdark(z_u, sqrt(x_u**2+y_u**2+z_u**2))
!
                  call coeff3d_linecu_lin(iim1, i, jjm1, j, kkm1, k, x_u, y_u, z_u, &
                     c14_slineu, c15_slineu, c17_slineu, c18_slineu, &
                     c23_slineu, c24_slineu, c26_slineu, c27_slineu, &
                     opac_u, opalbar_u, velx_u, vely_u, velz_u, vth_u, scont_u, sline_u)
!
               elseif(dels_xyu.eq.dels_u) then
!intersection with x-y plane on level k-gamma
                  x_u = x(i) - dels_u*nn_x
                  y_u = y(j) - dels_u*nn_y
                  z_u = z(kkm1)
!
                  call coeff2d_linecu_lin(opac3d(iim1,jjm1,kkm1), opac3d(i,jjm1,kkm1), opac3d(iim1,j,kkm1), opac3d(i,j,kkm1), &
                     opalbar3d(iim1,jjm1,kkm1), opalbar3d(i,jjm1,kkm1), opalbar3d(iim1,j,kkm1), opalbar3d(i,j,kkm1), &
                     velx3d(iim1,jjm1,kkm1), velx3d(i,jjm1,kkm1), velx3d(iim1,j,kkm1), velx3d(i,j,kkm1), &
                     vely3d(iim1,jjm1,kkm1), vely3d(i,jjm1,kkm1), vely3d(iim1,j,kkm1), vely3d(i,j,kkm1), &
                     velz3d(iim1,jjm1,kkm1), velz3d(i,jjm1,kkm1), velz3d(iim1,j,kkm1), velz3d(i,j,kkm1), &
                     vth3d(iim1,jjm1,kkm1), vth3d(i,jjm1,kkm1), vth3d(iim1,j,kkm1), vth3d(i,j,kkm1), &
                     scont3d(iim1,jjm1,kkm1), scont3d(i,jjm1,kkm1), scont3d(iim1,j,kkm1), scont3d(i,j,kkm1), &
                     sline3d(iim1,jjm1,kkm1), sline3d(i,jjm1,kkm1), sline3d(iim1,j,kkm1), sline3d(i,j,kkm1), &
                     int3d(iim1,jjm1,kkm1), int3d(i,jjm1,kkm1), int3d(iim1,j,kkm1), int3d(i,j,kkm1), &
                     x(iim1), x(i), y(jjm1), y(j), x_u, y_u, &
                     c14_slineu, c15_slineu, c17_slineu, c18_slineu, &
                     c14_intu, c15_intu, c17_intu, c18_intu, &
                     opac_u, opalbar_u, velx_u, vely_u, velz_u, vth_u, scont_u, sline_u, int_u)
                  c23_slineu=0.d0
                  c24_slineu=0.d0
                  c26_slineu=0.d0
                  c27_slineu=0.d0
                  c23_intu=0.d0
                  c24_intu=0.d0
                  c26_intu=0.d0
!
               elseif(dels_xzu.eq.dels_u) then
!intersection with x-z plane on level j-beta
                  x_u = x(i) - dels_u*nn_x
                  y_u = y(jjm1)
                  z_u = z(k) - dels_u*nn_z
!
                  call coeff2d_linecu_lin(opac3d(iim1,jjm1,kkm1), opac3d(i,jjm1,kkm1), opac3d(iim1,jjm1,k), opac3d(i,jjm1,k), &
                     opalbar3d(iim1,jjm1,kkm1), opalbar3d(i,jjm1,kkm1), opalbar3d(iim1,jjm1,k), opalbar3d(i,jjm1,k), &
                     velx3d(iim1,jjm1,kkm1), velx3d(i,jjm1,kkm1), velx3d(iim1,jjm1,k), velx3d(i,jjm1,k), &
                     vely3d(iim1,jjm1,kkm1), vely3d(i,jjm1,kkm1), vely3d(iim1,jjm1,k), vely3d(i,jjm1,k), &
                     velz3d(iim1,jjm1,kkm1), velz3d(i,jjm1,kkm1), velz3d(iim1,jjm1,k), velz3d(i,jjm1,k), &
                     vth3d(iim1,jjm1,kkm1), vth3d(i,jjm1,kkm1), vth3d(iim1,jjm1,k), vth3d(i,jjm1,k), &
                     scont3d(iim1,jjm1,kkm1), scont3d(i,jjm1,kkm1), scont3d(iim1,jjm1,k), scont3d(i,jjm1,k), &
                     sline3d(iim1,jjm1,kkm1), sline3d(i,jjm1,kkm1), sline3d(iim1,jjm1,k), sline3d(i,jjm1,k), &
                     int3d(iim1,jjm1,kkm1), int3d(i,jjm1,kkm1), int3d(iim1,jjm1,k), int3d(i,jjm1,k), &
                     x(iim1), x(i), z(kkm1), z(k), x_u, z_u, &
                     c14_slineu, c15_slineu, c23_slineu, c24_slineu, &
                     c14_intu, c15_intu, c23_intu, c24_intu, &
                     opac_u, opalbar_u, velx_u, vely_u, velz_u, vth_u, scont_u, sline_u, int_u)
                  c17_slineu=0.d0
                  c18_slineu=0.d0
                  c26_slineu=0.d0
                  c27_slineu=0.d0
                  c17_intu=0.d0
                  c18_intu=0.d0
                  c26_intu=0.d0
!
               elseif(dels_yzu.eq.dels_u) then
!intersection with y-z plane on level i-alpha
                  x_u = x(iim1)
                  y_u = y(j) - dels_u*nn_y
                  z_u = z(k) - dels_u*nn_z
!
                  call coeff2d_linecu_lin(opac3d(iim1,jjm1,kkm1), opac3d(iim1,j,kkm1), opac3d(iim1,jjm1,k), opac3d(iim1,j,k), &
                     opalbar3d(iim1,jjm1,kkm1), opalbar3d(iim1,j,kkm1), opalbar3d(iim1,jjm1,k), opalbar3d(iim1,j,k), &
                     velx3d(iim1,jjm1,kkm1), velx3d(iim1,j,kkm1), velx3d(iim1,jjm1,k), velx3d(iim1,j,k), &
                     vely3d(iim1,jjm1,kkm1), vely3d(iim1,j,kkm1), vely3d(iim1,jjm1,k), vely3d(iim1,j,k), &
                     velz3d(iim1,jjm1,kkm1), velz3d(iim1,j,kkm1), velz3d(iim1,jjm1,k), velz3d(iim1,j,k), &
                     vth3d(iim1,jjm1,kkm1), vth3d(iim1,j,kkm1), vth3d(iim1,jjm1,k), vth3d(iim1,j,k), &
                     scont3d(iim1,jjm1,kkm1), scont3d(iim1,j,kkm1), scont3d(iim1,jjm1,k), scont3d(iim1,j,k), &
                     sline3d(iim1,jjm1,kkm1), sline3d(iim1,j,kkm1), sline3d(iim1,jjm1,k), sline3d(iim1,j,k), &
                     int3d(iim1,jjm1,kkm1), int3d(iim1,j,kkm1), int3d(iim1,jjm1,k), int3d(iim1,j,k), &
                     y(jjm1), y(j), z(kkm1), z(k), y_u, z_u, &
                     c14_slineu, c17_slineu, c23_slineu, c26_slineu, &
                     c14_intu, c17_intu, c23_intu, c26_intu, &
                     opac_u, opalbar_u, velx_u, vely_u, velz_u, vth_u, scont_u, sline_u, int_u)
                  c15_slineu=0.d0
                  c18_slineu=0.d0
                  c24_slineu=0.d0
                  c27_slineu=0.d0
                  c15_intu=0.d0
                  c18_intu=0.d0
                  c24_intu=0.d0
!
               else
                  write(*,'(4es20.8,2l4)') dels_u, dels_xzu, dels_xyu, dels_yzu, dels_u.eq.dels_xzu, dels_xyu.eq.dels_u
                  stop 'error in fsc_linec3d: invalid dels_u'
               endif
!
!--------------------------------radiative transfer---------------------
!
               vel_u = velx_u*nn_x + vely_u*nn_y + velz_u*nn_z
               call fsc_linec_lin(xobs, vth_fiducial, deltax, xcmf_max, xcmf_min, int_u, opac_u, opac_p, opalbar_u, opalbar_p, scont_u, scont_p, &
                  sline_u, sline_p, vel_u, vel_p, vth_u, vth_p, dels_u, abs_sc, contr_sc, int_sc, alo_u, alo_p)
               int3d(i,j,k) = int_sc
!
               aloline_on_nn3d(i,j,k,q14) = imask_totreg3d(i,j,k) * (alo_p + alo_u*c27_slineu)
               aloline_on_nn3d(iim1,j,k,q15) = imask_totreg3d(iim1,j,k) * &
                  (alo_u*c26_slineu + abs_sc*c26_intu*aloline_on_nn3d(iim1,j,k,q14))

               aloline_on_nn3d(i,jjm1,k,q17) = imask_totreg3d(i,jjm1,k) * &
                  (alo_u*c24_slineu + abs_sc*c24_intu*aloline_on_nn3d(i,jjm1,k,q14))
               aloline_on_nn3d(iim1,jjm1,k,q18) = imask_totreg3d(iim1,jjm1,k) * &
                  (alo_u*c23_slineu + abs_sc*(c23_intu*aloline_on_nn3d(iim1,jjm1,k,q14) + &
                  c26_intu*aloline_on_nn3d(iim1,jjm1,k,q17)))
               aloline_on_nn3d(i,j,kkm1,q23) = imask_totreg3d(i,j,kkm1) * &
                  (alo_u*c18_slineu + abs_sc*c18_intu*aloline_on_nn3d(i,j,kkm1,q14))
               aloline_on_nn3d(iim1,j,kkm1,q24) = imask_totreg3d(iim1,j,kkm1) * &
                  (alo_u*c17_slineu + abs_sc*(c17_intu*aloline_on_nn3d(iim1,j,kkm1,q14) + &
                  c18_intu*aloline_on_nn3d(iim1,j,kkm1,q15) + &
                  c26_intu*aloline_on_nn3d(iim1,j,kkm1,q23)))
               aloline_on_nn3d(i,jjm1,kkm1,q26) = imask_totreg3d(i,jjm1,kkm1) * &
                  (alo_u*c15_slineu + abs_sc*(c15_intu*aloline_on_nn3d(i,jjm1,kkm1,q14) + &
                  c18_intu*aloline_on_nn3d(i,jjm1,kkm1,q17) + &
                  c24_intu*aloline_on_nn3d(i,jjm1,kkm1,q23)))
               aloline_on_nn3d(iim1,jjm1,kkm1,q27) = imask_totreg3d(iim1,jjm1,kkm1) * &
                  (alo_u*c14_slineu + abs_sc*(c14_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q14) + &
                  c15_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q15) + &
                  c17_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q17) + &
                  c18_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q18) + &
                  c23_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q23) + &
                  c24_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q24) + &
                  c26_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q26)))
!
!perform integration
               call calc_phinorm(vel_p, vth_p, vth_fiducial, xobs, phinorm)
               wall=wall_global*phinorm
               mintbar3d_tmp(i,j,k) = mintbar3d_tmp(i,j,k) + int3d(i,j,k)*wall
               normalization3d_tmp(i,j,k) = normalization3d_tmp(i,j,k) + wall
!
               aloline_nn3d_tmp(i,j,k,q14) = aloline_nn3d_tmp(i,j,k,q14) + wall*aloline_on_nn3d(i,j,k,q14)
               aloline_nn3d_tmp(iim1,j,k,q15) = aloline_nn3d_tmp(iim1,j,k,q15) + wall*aloline_on_nn3d(iim1,j,k,q15)
               aloline_nn3d_tmp(i,jjm1,k,q17) = aloline_nn3d_tmp(i,jjm1,k,q17) + wall*aloline_on_nn3d(i,jjm1,k,q17)
               aloline_nn3d_tmp(iim1,jjm1,k,q18) = aloline_nn3d_tmp(iim1,jjm1,k,q18) + wall*aloline_on_nn3d(iim1,jjm1,k,q18)
               aloline_nn3d_tmp(i,j,kkm1,q23) = aloline_nn3d_tmp(i,j,kkm1,q23) + wall*aloline_on_nn3d(i,j,kkm1,q23)
               aloline_nn3d_tmp(iim1,j,kkm1,q24) = aloline_nn3d_tmp(iim1,j,kkm1,q24) + wall*aloline_on_nn3d(iim1,j,kkm1,q24)
               aloline_nn3d_tmp(i,jjm1,kkm1,q26) = aloline_nn3d_tmp(i,jjm1,kkm1,q26) + wall*aloline_on_nn3d(i,jjm1,kkm1,q26)
               aloline_nn3d_tmp(iim1,jjm1,kkm1,q27) = aloline_nn3d_tmp(iim1,jjm1,kkm1,q27) + wall*aloline_on_nn3d(iim1,jjm1,kkm1,q27)
!
!
!***********special treatment: point lies directly on boundary**********
!
             case(3)
               mueff=nn_x*x(i)+nn_y*y(j)+nn_z*z(k)
               if(mueff.ge.0.d0) then
!no interpolations required
!for alo
!                 int3d(i,j,k) = 0.d0
!                  int3d(i,j,k) = xic1
                  int3d(i,j,k)=calc_icore_gdark(z(k), sqrt(x(i)**2+y(j)**2+z(k)**2))
                  aloline_on_nn3d(i,j,k,14) = 0.d0
!
!perform angular integration (alo not required, since boundary specified)
                  vel_p = nn_x*velx3d(i,j,k) + nn_y*vely3d(i,j,k) + nn_z*velz3d(i,j,k)
                  vth_p = vth3d(i,j,k)
                  call calc_phinorm(vel_p, vth_p, vth_fiducial, xobs, phinorm)
                  wall=wall_global*phinorm
                  mintbar3d_tmp(i,j,k) = mintbar3d_tmp(i,j,k) + int3d(i,j,k)*wall
                  normalization3d_tmp(i,j,k) = normalization3d_tmp(i,j,k) + wall
               else
!same interpolation as in (1)
!
                  iim1=i-alpha
                  jjm1=j-beta
                  kkm1=k-gamma
!
!calculate distance of upwind point to (previous) xy-plane, xz-plane and yz-plane
                  dels_xyu=(z(k)-z(kkm1))/nn_z
                  dels_xzu=(y(j)-y(jjm1))/nn_y
                  dels_yzu=(x(i)-x(iim1))/nn_x
                  dels_u=min(dels_xyu,dels_xzu,dels_yzu)
!
!----------------------------local point--------------------------------
!
                  sline_p=sline3d(i,j,k)
                  scont_p=scont3d(i,j,k)
                  opalbar_p=opalbar3d(i,j,k)
                  opac_p=opac3d(i,j,k)
                  velx_p=velx3d(i,j,k)
                  vely_p=vely3d(i,j,k)
                  velz_p=velz3d(i,j,k)
                  vel_p = nn_x*velx_p + nn_y*vely_p + nn_z*velz_p
                  vth_p = vth3d(i,j,k)
!
!----------------------------upwind point-------------------------------
!
                  if(dels_xyu.eq.dels_u) then
!intersection with x-y plane on level k-gamma
                     x_u = x(i) - dels_u*nn_x
                     y_u = y(j) - dels_u*nn_y
                     z_u = z(kkm1)
!
                     call coeff2d_linecu_lin(opac3d(iim1,jjm1,kkm1), opac3d(i,jjm1,kkm1), opac3d(iim1,j,kkm1), opac3d(i,j,kkm1), &
                        opalbar3d(iim1,jjm1,kkm1), opalbar3d(i,jjm1,kkm1), opalbar3d(iim1,j,kkm1), opalbar3d(i,j,kkm1), &
                        velx3d(iim1,jjm1,kkm1), velx3d(i,jjm1,kkm1), velx3d(iim1,j,kkm1), velx3d(i,j,kkm1), &
                        vely3d(iim1,jjm1,kkm1), vely3d(i,jjm1,kkm1), vely3d(iim1,j,kkm1), vely3d(i,j,kkm1), &
                        velz3d(iim1,jjm1,kkm1), velz3d(i,jjm1,kkm1), velz3d(iim1,j,kkm1), velz3d(i,j,kkm1), &
                        vth3d(iim1,jjm1,kkm1), vth3d(i,jjm1,kkm1), vth3d(iim1,j,kkm1), vth3d(i,j,kkm1), &
                        scont3d(iim1,jjm1,kkm1), scont3d(i,jjm1,kkm1), scont3d(iim1,j,kkm1), scont3d(i,j,kkm1), &
                        sline3d(iim1,jjm1,kkm1), sline3d(i,jjm1,kkm1), sline3d(iim1,j,kkm1), sline3d(i,j,kkm1), &
                        int3d(iim1,jjm1,kkm1), int3d(i,jjm1,kkm1), int3d(iim1,j,kkm1), int3d(i,j,kkm1), &
                        x(iim1), x(i), y(jjm1), y(j), x_u, y_u, &
                        c14_slineu, c15_slineu, c17_slineu, c18_slineu, &
                        c14_intu, c15_intu, c17_intu, c18_intu, &
                        opac_u, opalbar_u, velx_u, vely_u, velz_u, vth_u, scont_u, sline_u, int_u)
                     c23_slineu=0.d0
                     c24_slineu=0.d0
                     c26_slineu=0.d0
                     c23_intu=0.d0
                     c24_intu=0.d0
                     c26_intu=0.d0
!
                  elseif(dels_xzu.eq.dels_u) then
!intersection with x-z plane on level j-beta
                     x_u = x(i) - dels_u*nn_x
                     y_u = y(jjm1)
                     z_u = z(k) - dels_u*nn_z
!
                     call coeff2d_linecu_lin(opac3d(iim1,jjm1,kkm1), opac3d(i,jjm1,kkm1), opac3d(iim1,jjm1,k), opac3d(i,jjm1,k), &
                        opalbar3d(iim1,jjm1,kkm1), opalbar3d(i,jjm1,kkm1), opalbar3d(iim1,jjm1,k), opalbar3d(i,jjm1,k), &
                        velx3d(iim1,jjm1,kkm1), velx3d(i,jjm1,kkm1), velx3d(iim1,jjm1,k), velx3d(i,jjm1,k), &
                        vely3d(iim1,jjm1,kkm1), vely3d(i,jjm1,kkm1), vely3d(iim1,jjm1,k), vely3d(i,jjm1,k), &
                        velz3d(iim1,jjm1,kkm1), velz3d(i,jjm1,kkm1), velz3d(iim1,jjm1,k), velz3d(i,jjm1,k), &
                        vth3d(iim1,jjm1,kkm1), vth3d(i,jjm1,kkm1), vth3d(iim1,jjm1,k), vth3d(i,jjm1,k), &
                        scont3d(iim1,jjm1,kkm1), scont3d(i,jjm1,kkm1), scont3d(iim1,jjm1,k), scont3d(i,jjm1,k), &
                        sline3d(iim1,jjm1,kkm1), sline3d(i,jjm1,kkm1), sline3d(iim1,jjm1,k), sline3d(i,jjm1,k), &
                        int3d(iim1,jjm1,kkm1), int3d(i,jjm1,kkm1), int3d(iim1,jjm1,k), int3d(i,jjm1,k), &
                        x(iim1), x(i), z(kkm1), z(k), x_u, z_u, &
                        c14_slineu, c15_slineu, c23_slineu, c24_slineu, &
                        c14_intu, c15_intu, c23_intu, c24_intu, &
                        opac_u, opalbar_u, velx_u, vely_u, velz_u, vth_u, scont_u, sline_u, int_u)
                     c17_slineu=0.d0
                     c18_slineu=0.d0
                     c26_slineu=0.d0
                     c17_intu=0.d0
                     c18_intu=0.d0
                     c26_intu=0.d0
!
                  elseif(dels_yzu.eq.dels_u) then
!intersection with y-z plane on level i-alpha
                     x_u = x(iim1)
                     y_u = y(j) - dels_u*nn_y
                     z_u = z(k) - dels_u*nn_z
!
                     call coeff2d_linecu_lin(opac3d(iim1,jjm1,kkm1), opac3d(iim1,j,kkm1), opac3d(iim1,jjm1,k), opac3d(iim1,j,k), &
                        opalbar3d(iim1,jjm1,kkm1), opalbar3d(iim1,j,kkm1), opalbar3d(iim1,jjm1,k), opalbar3d(iim1,j,k), &
                        velx3d(iim1,jjm1,kkm1), velx3d(iim1,j,kkm1), velx3d(iim1,jjm1,k), velx3d(iim1,j,k), &
                        vely3d(iim1,jjm1,kkm1), vely3d(iim1,j,kkm1), vely3d(iim1,jjm1,k), vely3d(iim1,j,k), &
                        velz3d(iim1,jjm1,kkm1), velz3d(iim1,j,kkm1), velz3d(iim1,jjm1,k), velz3d(iim1,j,k), &
                        vth3d(iim1,jjm1,kkm1), vth3d(iim1,j,kkm1), vth3d(iim1,jjm1,k), vth3d(iim1,j,k), &
                        scont3d(iim1,jjm1,kkm1), scont3d(iim1,j,kkm1), scont3d(iim1,jjm1,k), scont3d(iim1,j,k), &
                        sline3d(iim1,jjm1,kkm1), sline3d(iim1,j,kkm1), sline3d(iim1,jjm1,k), sline3d(iim1,j,k), &
                        int3d(iim1,jjm1,kkm1), int3d(iim1,j,kkm1), int3d(iim1,jjm1,k), int3d(iim1,j,k), &
                        y(jjm1), y(j), z(kkm1), z(k), y_u, z_u, &
                        c14_slineu, c17_slineu, c23_slineu, c26_slineu, &
                        c14_intu, c17_intu, c23_intu, c26_intu, &
                        opac_u, opalbar_u, velx_u, vely_u, velz_u, vth_u, scont_u, sline_u, int_u)
                     c15_slineu=0.d0
                     c18_slineu=0.d0
                     c24_slineu=0.d0
                     c15_intu=0.d0
                     c18_intu=0.d0
                     c24_intu=0.d0
!
                  else
                     write(*,'(4es20.8,2l4)') dels_u, dels_xzu, dels_xyu, dels_yzu, dels_u.eq.dels_xzu, dels_xyu.eq.dels_u
                     stop 'error in fsc_linec3d: invalid dels_u'
                  endif
!
!--------------------------------radiative transfer---------------------
!
                  vel_u = velx_u*nn_x + vely_u*nn_y + velz_u*nn_z
                  call fsc_linec_lin(xobs, vth_fiducial, deltax, xcmf_max, xcmf_min, int_u, opac_u, opac_p, opalbar_u, opalbar_p, scont_u, scont_p, &
                     sline_u, sline_p, vel_u, vel_p, vth_u, vth_p, dels_u, abs_sc, contr_sc, int_sc, alo_u, alo_p)
                  int3d(i,j,k) = int_sc
!
                  aloline_on_nn3d(i,j,k,q14) = imask_totreg3d(i,j,k) * alo_p
                  aloline_on_nn3d(iim1,j,k,q15) = imask_totreg3d(iim1,j,k) * &
                     (alo_u*c26_slineu + abs_sc*c26_intu*aloline_on_nn3d(iim1,j,k,q14))

                  aloline_on_nn3d(i,jjm1,k,q17) = imask_totreg3d(i,jjm1,k) * &
                     (alo_u*c24_slineu + abs_sc*c24_intu*aloline_on_nn3d(i,jjm1,k,q14))
                  aloline_on_nn3d(iim1,jjm1,k,q18) = imask_totreg3d(iim1,jjm1,k) * &
                     (alo_u*c23_slineu + abs_sc*(c23_intu*aloline_on_nn3d(iim1,jjm1,k,q14) + &
                     c26_intu*aloline_on_nn3d(iim1,jjm1,k,q17)))
                  aloline_on_nn3d(i,j,kkm1,q23) = imask_totreg3d(i,j,kkm1) * &
                     (alo_u*c18_slineu + abs_sc*c18_intu*aloline_on_nn3d(i,j,kkm1,q14))
                  aloline_on_nn3d(iim1,j,kkm1,q24) = imask_totreg3d(iim1,j,kkm1) * &
                     (alo_u*c17_slineu + abs_sc*(c17_intu*aloline_on_nn3d(iim1,j,kkm1,q14) + &
                     c18_intu*aloline_on_nn3d(iim1,j,kkm1,q15) + &
                     c26_intu*aloline_on_nn3d(iim1,j,kkm1,q23)))
                  aloline_on_nn3d(i,jjm1,kkm1,q26) = imask_totreg3d(i,jjm1,kkm1) * &
                     (alo_u*c15_slineu + abs_sc*(c15_intu*aloline_on_nn3d(i,jjm1,kkm1,q14) + &
                     c18_intu*aloline_on_nn3d(i,jjm1,kkm1,q17) + &
                     c24_intu*aloline_on_nn3d(i,jjm1,kkm1,q23)))
                  aloline_on_nn3d(iim1,jjm1,kkm1,q27) = imask_totreg3d(iim1,jjm1,kkm1) * &
                     (alo_u*c14_slineu + abs_sc*(c14_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q14) + &
                     c15_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q15) + &
                     c17_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q17) + &
                     c18_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q18) + &
                     c23_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q23) + &
                     c24_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q24) + &
                     c26_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q26)))
!
!integration
                  call calc_phinorm(vel_p, vth_p, vth_fiducial, xobs, phinorm)
                  wall=wall_global*phinorm
!
                  mintbar3d_tmp(i,j,k) = mintbar3d_tmp(i,j,k) + int3d(i,j,k)*wall
                  normalization3d_tmp(i,j,k) = normalization3d_tmp(i,j,k) + wall
!
                  aloline_nn3d_tmp(i,j,k,q14) = aloline_nn3d_tmp(i,j,k,q14) + wall*aloline_on_nn3d(i,j,k,q14)
                  aloline_nn3d_tmp(iim1,j,k,q15) = aloline_nn3d_tmp(iim1,j,k,q15) + wall*aloline_on_nn3d(iim1,j,k,q15)
                  aloline_nn3d_tmp(i,jjm1,k,q17) = aloline_nn3d_tmp(i,jjm1,k,q17) + wall*aloline_on_nn3d(i,jjm1,k,q17)
                  aloline_nn3d_tmp(iim1,jjm1,k,q18) = aloline_nn3d_tmp(iim1,jjm1,k,q18) + wall*aloline_on_nn3d(iim1,jjm1,k,q18)
                  aloline_nn3d_tmp(i,j,kkm1,q23) = aloline_nn3d_tmp(i,j,kkm1,q23) + wall*aloline_on_nn3d(i,j,kkm1,q23)
                  aloline_nn3d_tmp(iim1,j,kkm1,q24) = aloline_nn3d_tmp(iim1,j,kkm1,q24) + wall*aloline_on_nn3d(iim1,j,kkm1,q24)
                  aloline_nn3d_tmp(i,jjm1,kkm1,q26) = aloline_nn3d_tmp(i,jjm1,kkm1,q26) + wall*aloline_on_nn3d(i,jjm1,kkm1,q26)
                  aloline_nn3d_tmp(iim1,jjm1,kkm1,q27) = aloline_nn3d_tmp(iim1,jjm1,kkm1,q27) + wall*aloline_on_nn3d(iim1,jjm1,kkm1,q27)
!
               endif
!
             case default
!
            end select
!
!         if(mintbar3d_tmp(i,j,k).ne.mintbar3d_tmp(i,j,k)) then
!            write(*,*) mintbar3d_tmp(i,j,k), int3d(i,j,k), imask3d(i,j,k)
!            stop
!         endif
!
         enddo
      enddo
   enddo
!
!
!
end subroutine fsc_linec3d_lin
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine fsc_linec3d(oindx,xobsindx)
!
!-----------------------------------------------------------------------
!--------short characteristics for line radiative transfer in 3d--------
!-----------calculating intensties for given mu,phi specified-----------
!-----------by input oindx, and for xobs specified by xobsindx----------
!---------------------using bezier interpolations-----------------------
!-----------------------------------------------------------------------
!
   use prog_type
   use fund_const
   use dime3d, only: ndxmax, ndzmax, ndymax, x, y, z
   use dime3d, only: int3d, opac3d, opalbar3d, scont3d, sline3d, velx3d, vely3d, velz3d, &
      vth3d, imask3d, imask_totreg3d, aloline_on_nn3d, aloline_nn3d_tmp, &
      mintbar3d_tmp, normalization3d_tmp
   use angles, only: n_x, n_y, n_z, weight_omega, q_alo
   use freq, only: nodes_xobs, weight_xobs, nxobs, deltax, xcmf_min, xcmf_max
   use params_input, only: vth_fiducial
   use bcondition, only: xic1
   use mod_debug, only: indxx, indxy, indxz
   use params_stellar, only: smajorax_a, smajorax_b, smajorax_c
!
   implicit none
!
! ... arguments
   integer(i4b), intent(in) :: oindx, xobsindx
!
! ... local scalars
   integer(i4b) :: i, j, k
   integer(i4b) :: alpha, beta, gamma
   integer(i4b) :: startx, starty, startz, endx, endy, endz
   integer(i4b) :: iim2, iim1, ii, iip1, jjm2, jjm1, jj, jjp1, kkm2, kkm1, kk, kkp1
   real(dp) :: nn_x, nn_y, nn_z, mueff, xobs, wall, wall_global, phinorm
   real(dp) :: dels_xyu, dels_xzu, dels_yzu, dels_u, dels_r, &
      dels_xyd, dels_xzd, dels_yzd, dels_d
   real(dp) :: opac_p, opalbar_p, velx_p, vely_p, velz_p, vel_p, vth_p, scont_p, sline_p, rad_p
   real(dp) :: x_u, y_u, z_u, opac_u, opalbar_u, velx_u, vely_u, velz_u, vel_u, vth_u, scont_u, sline_u, int_u, rad_u
   real(dp) :: x_d, y_d, z_d, opac_d, opalbar_d, velx_d, vely_d, velz_d, vel_d, vth_d, scont_d, sline_d, rad_d
   real(dp) :: abs_sc, int_sc, contr_sc
   real(dp) :: alo_u, alo_p, alo_d
   integer :: q1,  q2,  q3,  q4,  q5,  q6,  q7,  q8,  q9, &
      q10, q11, q12, q13, q14, q15, q16, q17, q18, &
      q19, q20, q21, q22, q23, q24, q25, q26, q27
   real(dp) :: c02_slineu, c04_slineu, c05_slineu, c06_slineu, c08_slineu, &
      c10_slineu, c11_slineu, c12_slineu, c13_slineu, c14_slineu, c15_slineu, c16_slineu, c17_slineu, c18_slineu, &
      c20_slineu, c22_slineu, c23_slineu, c24_slineu, c26_slineu, c27_slineu, &
      c02_intu, c04_intu, c05_intu, c06_intu, c08_intu, &
      c10_intu, c11_intu, c12_intu, c13_intu, c14_intu, c15_intu, c16_intu, c17_intu, c18_intu, &
      c20_intu, c22_intu, c23_intu, c24_intu, c26_intu, &
      c03_slined, c06_slined, c07_slined, c08_slined, c09_slined, c12_slined, &
      c15_slined, c16_slined, c17_slined, c18_slined, c19_slined, c20_slined, &
      c21_slined, c22_slined, c23_slined, c24_slined, c25_slined, c26_slined, c27_slined
!
! ... local arrays
!
! ... local functions
   real(dp) :: dist_sphere, dist_ellipsoid, calc_icore_gdark
!
! ... local characters
!character(len=50) :: enter
!
! ... local logicals
!
!directions
   nn_x=n_x(oindx)
   nn_y=n_y(oindx)
   nn_z=n_z(oindx)
!
!frequency
   xobs=nodes_xobs(xobsindx)
!
!angular and frequency  integration weights
   wall_global=weight_omega(oindx)*weight_xobs(xobsindx)
!
!indices for nearest neighbour alo
   q1=q_alo(oindx,1)
   q2=q_alo(oindx,2)
   q3=q_alo(oindx,3)
   q4=q_alo(oindx,4)
   q5=q_alo(oindx,5)
   q6=q_alo(oindx,6)
   q7=q_alo(oindx,7)
   q8=q_alo(oindx,8)
   q9=q_alo(oindx,9)
   q10=q_alo(oindx,10)
   q11=q_alo(oindx,11)
   q12=q_alo(oindx,12)
   q13=q_alo(oindx,13)
   q14=q_alo(oindx,14)
   q15=q_alo(oindx,15)
   q16=q_alo(oindx,16)
   q17=q_alo(oindx,17)
   q18=q_alo(oindx,18)
   q19=q_alo(oindx,19)
   q20=q_alo(oindx,20)
   q21=q_alo(oindx,21)
   q22=q_alo(oindx,22)
   q23=q_alo(oindx,23)
   q24=q_alo(oindx,24)
   q25=q_alo(oindx,25)
   q26=q_alo(oindx,26)
   q27=q_alo(oindx,27)
!
!-----------------------------------------------------------------------
!
!set directional index-parameter (phantom points are excluded from calculation)
!
!index parameter:
!         if n_x >= 0                 if n_x < 0
!                startx = 2                  startx = ndxmax-1
!                endx = ndxmax-1             endx = 2
!                alpha=  1                   alpha=-1
!
!         if n_y >= 0                 if n_y < 0
!                starty = 2                  starty = ndymax-1
!                endy = ndymax-1             endy = 2
!                beta =  1                   beta =-1
!
!         if n_z >= 0                 if n_z < 0
!                startz = 2                  startz = ndzmax-1
!                endz = ndzmax-1             endz = 2
!                gamma = 1                   gamma = -1
!
   if(nn_x.gt.0.d0) then
      startx = 3
      endx = ndxmax-2
      alpha=  1
   elseif(nn_x.lt.0.d0) then
      startx = ndxmax-2
      endx = 3
      alpha=-1
   else
      stop 'error in fsc_linec3d: n_x = 0 not allowed'
   endif
!
   if(nn_y.gt.0.d0) then
      starty = 3
      endy = ndymax-2
      beta =  1
   elseif(nn_y.lt.0.d0) then
      starty = ndymax-2
      endy = 3
      beta =-1
   else
      stop 'error in fsc_linec3d: n_y = 0 not allowed'
   endif
!
   if(nn_z.gt.0.d0) then
      startz = 3
      endz = ndzmax-2
      gamma=  1
   elseif(nn_z.lt.0.d0) then
      startz = ndzmax-2
      endz = 3
      gamma=-1
   else
      stop 'error in fsc_linec3d: n_z = 0 not allowed'
   endif
!
!--------------------reset the intensities and alo----------------------
!
   call set_boundary3d(nn_x, nn_y, nn_z)
!
   aloline_on_nn3d=0.d0
!
!-----------------------------------------------------------------------
!
   do k=startz, endz, gamma
      do j=starty, endy, beta
         do i=startx, endx, alpha
!
            select case(imask3d(i,j,k))
!
!************************standard radiative transfer********************
!
             case(1)
!
               iim1=i-alpha
               jjm1=j-beta
               kkm1=k-gamma
               iip1=i+alpha
               jjp1=j+beta
               kkp1=k+gamma
!
!calculate distance of upwind point to (previous) xy-plane, xz-plane and yz-plane
               dels_xyu=(z(k)-z(kkm1))/nn_z
               dels_xzu=(y(j)-y(jjm1))/nn_y
               dels_yzu=(x(i)-x(iim1))/nn_x
               dels_u=min(dels_xyu,dels_xzu,dels_yzu)
!
!calculate distance of downwind point to (next) xy-plane, xz-plane and yz-plane
               dels_xyd=(z(kkp1)-z(k))/nn_z
               dels_xzd=(y(jjp1)-y(j))/nn_y
               dels_yzd=(x(iip1)-x(i))/nn_x
               dels_d=min(dels_xyd,dels_xzd,dels_yzd)
!
!----------------------------local point--------------------------------
!
               sline_p=sline3d(i,j,k)
               scont_p=scont3d(i,j,k)
               opalbar_p=opalbar3d(i,j,k)
               opac_p=opac3d(i,j,k)
               velx_p=velx3d(i,j,k)
               vely_p=vely3d(i,j,k)
               velz_p=velz3d(i,j,k)
               vel_p = nn_x*velx_p + nn_y*vely_p + nn_z*velz_p
               vth_p = vth3d(i,j,k)
!
!----------------------------upwind point-------------------------------
!
               if(dels_xyu.eq.dels_u) then
!intersection with x-y plane on level k-gamma
                  x_u = x(i) - dels_u*nn_x
                  y_u = y(j) - dels_u*nn_y
                  z_u = z(kkm1)
!
                  iim2=i-2*alpha
                  jjm2=j-2*beta
!
                  call coeff2d_linecu(opac3d(iim2,jjm2,kkm1), opac3d(iim1,jjm2,kkm1), opac3d(i,jjm2,kkm1), &
                     opac3d(iim2,jjm1,kkm1), opac3d(iim1,jjm1,kkm1), opac3d(i,jjm1,kkm1), &
                     opac3d(iim2,j,kkm1),    opac3d(iim1,j,kkm1),    opac3d(i,j,kkm1), &
                     opalbar3d(iim2,jjm2,kkm1), opalbar3d(iim1,jjm2,kkm1), opalbar3d(i,jjm2,kkm1), &
                     opalbar3d(iim2,jjm1,kkm1), opalbar3d(iim1,jjm1,kkm1), opalbar3d(i,jjm1,kkm1), &
                     opalbar3d(iim2,j,kkm1),    opalbar3d(iim1,j,kkm1),    opalbar3d(i,j,kkm1), &
                     velx3d(iim2,jjm2,kkm1), velx3d(iim1,jjm2,kkm1), velx3d(i,jjm2,kkm1), &
                     velx3d(iim2,jjm1,kkm1), velx3d(iim1,jjm1,kkm1), velx3d(i,jjm1,kkm1), &
                     velx3d(iim2,j,kkm1),    velx3d(iim1,j,kkm1),    velx3d(i,j,kkm1), &
                     vely3d(iim2,jjm2,kkm1), vely3d(iim1,jjm2,kkm1), vely3d(i,jjm2,kkm1), &
                     vely3d(iim2,jjm1,kkm1), vely3d(iim1,jjm1,kkm1), vely3d(i,jjm1,kkm1), &
                     vely3d(iim2,j,kkm1),    vely3d(iim1,j,kkm1),    vely3d(i,j,kkm1), &
                     velz3d(iim2,jjm2,kkm1), velz3d(iim1,jjm2,kkm1), velz3d(i,jjm2,kkm1), &
                     velz3d(iim2,jjm1,kkm1), velz3d(iim1,jjm1,kkm1), velz3d(i,jjm1,kkm1), &
                     velz3d(iim2,j,kkm1),    velz3d(iim1,j,kkm1),    velz3d(i,j,kkm1), &
                     vth3d(iim2,jjm2,kkm1), vth3d(iim1,jjm2,kkm1), vth3d(i,jjm2,kkm1), &
                     vth3d(iim2,jjm1,kkm1), vth3d(iim1,jjm1,kkm1), vth3d(i,jjm1,kkm1), &
                     vth3d(iim2,j,kkm1),    vth3d(iim1,j,kkm1),    vth3d(i,j,kkm1), &
                     scont3d(iim2,jjm2,kkm1), scont3d(iim1,jjm2,kkm1), scont3d(i,jjm2,kkm1), &
                     scont3d(iim2,jjm1,kkm1), scont3d(iim1,jjm1,kkm1), scont3d(i,jjm1,kkm1), &
                     scont3d(iim2,j,kkm1),    scont3d(iim1,j,kkm1),    scont3d(i,j,kkm1), &
                     sline3d(iim2,jjm2,kkm1), sline3d(iim1,jjm2,kkm1), sline3d(i,jjm2,kkm1), &
                     sline3d(iim2,jjm1,kkm1), sline3d(iim1,jjm1,kkm1), sline3d(i,jjm1,kkm1), &
                     sline3d(iim2,j,kkm1),    sline3d(iim1,j,kkm1),    sline3d(i,j,kkm1), &
                     int3d(iim2,jjm2,kkm1), int3d(iim1,jjm2,kkm1), int3d(i,jjm2,kkm1), &
                     int3d(iim2,jjm1,kkm1), int3d(iim1,jjm1,kkm1), int3d(i,jjm1,kkm1), &
                     int3d(iim2,j,kkm1),    int3d(iim1,j,kkm1),    int3d(i,j,kkm1), &
                     x(iim2), x(iim1), x(i), y(jjm2), y(jjm1), y(j), x_u, y_u, &
                     c10_slineu, c11_slineu, c12_slineu, c13_slineu, c14_slineu, &
                     c15_slineu, c16_slineu, c17_slineu, c18_slineu, &
                     c10_intu, c11_intu, c12_intu, c13_intu, c14_intu, &
                     c15_intu, c16_intu, c17_intu, c18_intu, &
                     opac_u, opalbar_u, velx_u, vely_u, velz_u, vth_u, scont_u, sline_u, int_u)
!set interpolation coefficients that are not used to zero
                  c23_slineu=0.d0
                  c24_slineu=0.d0
                  c26_slineu=0.d0
                  c27_slineu=0.d0
!
                  c02_intu = 0.d0
                  c04_intu = 0.d0
                  c05_intu = 0.d0
                  c06_intu = 0.d0
                  c08_intu = 0.d0
                  c20_intu = 0.d0
                  c22_intu = 0.d0
                  c23_intu = 0.d0
                  c24_intu = 0.d0
                  c26_intu = 0.d0
!
               elseif(dels_xzu.eq.dels_u) then
!intersection with x-z plane on level j-beta
                  x_u = x(i) - dels_u*nn_x
                  y_u = y(jjm1)
                  z_u = z(k) - dels_u*nn_z
!
                  iim2=i-2*alpha
                  kkm2=k-2*gamma
!
                  call coeff2d_linecu(   opac3d(iim2,jjm1,kkm2),    opac3d(iim1,jjm1,kkm2),    opac3d(i,jjm1,kkm2), &
                     opac3d(iim2,jjm1,kkm1),    opac3d(iim1,jjm1,kkm1),    opac3d(i,jjm1,kkm1), &
                     opac3d(iim2,jjm1,k),       opac3d(iim1,jjm1,k),       opac3d(i,jjm1,k), &
                     opalbar3d(iim2,jjm1,kkm2), opalbar3d(iim1,jjm1,kkm2), opalbar3d(i,jjm1,kkm2), &
                     opalbar3d(iim2,jjm1,kkm1), opalbar3d(iim1,jjm1,kkm1), opalbar3d(i,jjm1,kkm1), &
                     opalbar3d(iim2,jjm1,k),    opalbar3d(iim1,jjm1,k),    opalbar3d(i,jjm1,k), &
                     velx3d(iim2,jjm1,kkm2),    velx3d(iim1,jjm1,kkm2),    velx3d(i,jjm1,kkm2), &
                     velx3d(iim2,jjm1,kkm1),    velx3d(iim1,jjm1,kkm1),    velx3d(i,jjm1,kkm1), &
                     velx3d(iim2,jjm1,k),       velx3d(iim1,jjm1,k),       velx3d(i,jjm1,k), &
                     vely3d(iim2,jjm1,kkm2),    vely3d(iim1,jjm1,kkm2),    vely3d(i,jjm1,kkm2), &
                     vely3d(iim2,jjm1,kkm1),    vely3d(iim1,jjm1,kkm1),    vely3d(i,jjm1,kkm1), &
                     vely3d(iim2,jjm1,k),       vely3d(iim1,jjm1,k),       vely3d(i,jjm1,k), &
                     velz3d(iim2,jjm1,kkm2),    velz3d(iim1,jjm1,kkm2),    velz3d(i,jjm1,kkm2), &
                     velz3d(iim2,jjm1,kkm1),    velz3d(iim1,jjm1,kkm1),    velz3d(i,jjm1,kkm1), &
                     velz3d(iim2,jjm1,k),       velz3d(iim1,jjm1,k),       velz3d(i,jjm1,k), &
                     vth3d(iim2,jjm1,kkm2),     vth3d(iim1,jjm1,kkm2),     vth3d(i,jjm1,kkm2), &
                     vth3d(iim2,jjm1,kkm1),     vth3d(iim1,jjm1,kkm1),     vth3d(i,jjm1,kkm1), &
                     vth3d(iim2,jjm1,k),        vth3d(iim1,jjm1,k),        vth3d(i,jjm1,k), &
                     scont3d(iim2,jjm1,kkm2),   scont3d(iim1,jjm1,kkm2),   scont3d(i,jjm1,kkm2), &
                     scont3d(iim2,jjm1,kkm1),   scont3d(iim1,jjm1,kkm1),   scont3d(i,jjm1,kkm1), &
                     scont3d(iim2,jjm1,k),      scont3d(iim1,jjm1,k),      scont3d(i,jjm1,k), &
                     sline3d(iim2,jjm1,kkm2),   sline3d(iim1,jjm1,kkm2),   sline3d(i,jjm1,kkm2), &
                     sline3d(iim2,jjm1,kkm1),   sline3d(iim1,jjm1,kkm1),   sline3d(i,jjm1,kkm1), &
                     sline3d(iim2,jjm1,k),      sline3d(iim1,jjm1,k),      sline3d(i,jjm1,k), &
                     int3d(iim2,jjm1,kkm2),     int3d(iim1,jjm1,kkm2),     int3d(i,jjm1,kkm2), &
                     int3d(iim2,jjm1,kkm1),     int3d(iim1,jjm1,kkm1),     int3d(i,jjm1,kkm1), &
                     int3d(iim2,jjm1,k),        int3d(iim1,jjm1,k),        int3d(i,jjm1,k), &
                     x(iim2), x(iim1), x(i), z(kkm2), z(kkm1), z(k), x_u, z_u, &
                     c04_slineu, c05_slineu, c06_slineu, c13_slineu, c14_slineu, &
                     c15_slineu, c22_slineu, c23_slineu, c24_slineu, &
                     c04_intu, c05_intu, c06_intu, c13_intu, c14_intu, &
                     c15_intu, c22_intu, c23_intu, c24_intu, &
                     opac_u, opalbar_u, velx_u, vely_u, velz_u, vth_u, scont_u, sline_u, int_u)
                  c17_slineu=0.d0
                  c18_slineu=0.d0
                  c26_slineu=0.d0
                  c27_slineu=0.d0
                  c02_intu = 0.d0
                  c08_intu = 0.d0
                  c10_intu = 0.d0
                  c11_intu = 0.d0
                  c12_intu = 0.d0
                  c16_intu = 0.d0
                  c17_intu = 0.d0
                  c18_intu = 0.d0
                  c20_intu = 0.d0
                  c26_intu = 0.d0
!
               elseif(dels_yzu.eq.dels_u) then
!intersection with y-z plane on level i-alpha
                  x_u = x(iim1)
                  y_u = y(j) - dels_u*nn_y
                  z_u = z(k) - dels_u*nn_z
!
                  jjm2=j-2*beta
                  kkm2=k-2*gamma
!
                  call coeff2d_linecu(   opac3d(iim1,jjm2,kkm2),    opac3d(iim1,jjm1,kkm2),    opac3d(iim1,j,kkm2), &
                     opac3d(iim1,jjm2,kkm1),    opac3d(iim1,jjm1,kkm1),    opac3d(iim1,j,kkm1), &
                     opac3d(iim1,jjm2,k),       opac3d(iim1,jjm1,k),       opac3d(iim1,j,k), &
                     opalbar3d(iim1,jjm2,kkm2), opalbar3d(iim1,jjm1,kkm2), opalbar3d(iim1,j,kkm2), &
                     opalbar3d(iim1,jjm2,kkm1), opalbar3d(iim1,jjm1,kkm1), opalbar3d(iim1,j,kkm1), &
                     opalbar3d(iim1,jjm2,k),    opalbar3d(iim1,jjm1,k),    opalbar3d(iim1,j,k), &
                     velx3d(iim1,jjm2,kkm2),    velx3d(iim1,jjm1,kkm2),    velx3d(iim1,j,kkm2), &
                     velx3d(iim1,jjm2,kkm1),    velx3d(iim1,jjm1,kkm1),    velx3d(iim1,j,kkm1), &
                     velx3d(iim1,jjm2,k),       velx3d(iim1,jjm1,k),       velx3d(iim1,j,k), &
                     vely3d(iim1,jjm2,kkm2),    vely3d(iim1,jjm1,kkm2),    vely3d(iim1,j,kkm2), &
                     vely3d(iim1,jjm2,kkm1),    vely3d(iim1,jjm1,kkm1),    vely3d(iim1,j,kkm1), &
                     vely3d(iim1,jjm2,k),       vely3d(iim1,jjm1,k),       vely3d(iim1,j,k), &
                     velz3d(iim1,jjm2,kkm2),    velz3d(iim1,jjm1,kkm2),    velz3d(iim1,j,kkm2), &
                     velz3d(iim1,jjm2,kkm1),    velz3d(iim1,jjm1,kkm1),    velz3d(iim1,j,kkm1), &
                     velz3d(iim1,jjm2,k),       velz3d(iim1,jjm1,k),       velz3d(iim1,j,k), &
                     vth3d(iim1,jjm2,kkm2),     vth3d(iim1,jjm1,kkm2),     vth3d(iim1,j,kkm2), &
                     vth3d(iim1,jjm2,kkm1),     vth3d(iim1,jjm1,kkm1),     vth3d(iim1,j,kkm1), &
                     vth3d(iim1,jjm2,k),        vth3d(iim1,jjm1,k),        vth3d(iim1,j,k), &
                     scont3d(iim1,jjm2,kkm2),   scont3d(iim1,jjm1,kkm2),   scont3d(iim1,j,kkm2), &
                     scont3d(iim1,jjm2,kkm1),   scont3d(iim1,jjm1,kkm1),   scont3d(iim1,j,kkm1), &
                     scont3d(iim1,jjm2,k),      scont3d(iim1,jjm1,k),      scont3d(iim1,j,k), &
                     sline3d(iim1,jjm2,kkm2),   sline3d(iim1,jjm1,kkm2),   sline3d(iim1,j,kkm2), &
                     sline3d(iim1,jjm2,kkm1),   sline3d(iim1,jjm1,kkm1),   sline3d(iim1,j,kkm1), &
                     sline3d(iim1,jjm2,k),      sline3d(iim1,jjm1,k),      sline3d(iim1,j,k), &
                     int3d(iim1,jjm2,kkm2),     int3d(iim1,jjm1,kkm2),     int3d(iim1,j,kkm2), &
                     int3d(iim1,jjm2,kkm1),     int3d(iim1,jjm1,kkm1),     int3d(iim1,j,kkm1), &
                     int3d(iim1,jjm2,k),        int3d(iim1,jjm1,k),        int3d(iim1,j,k), &
                     y(jjm2), y(jjm1), y(j), z(kkm2), z(kkm1), z(k), y_u, z_u, &
                     c02_slineu, c05_slineu, c08_slineu, c11_slineu, c14_slineu, &
                     c17_slineu, c20_slineu, c23_slineu, c26_slineu, &
                     c02_intu, c05_intu, c08_intu, c11_intu, c14_intu, &
                     c17_intu, c20_intu, c23_intu, c26_intu, &
                     opac_u, opalbar_u, velx_u, vely_u, velz_u, vth_u, scont_u, sline_u, int_u)
                  c15_slineu=0.d0
                  c18_slineu=0.d0
                  c24_slineu=0.d0
                  c27_slineu=0.d0
!
                  c04_intu = 0.d0
                  c06_intu = 0.d0
                  c10_intu = 0.d0
                  c12_intu = 0.d0
                  c13_intu = 0.d0
                  c15_intu = 0.d0
                  c16_intu = 0.d0
                  c18_intu = 0.d0
                  c22_intu = 0.d0
                  c24_intu = 0.d0
!
               else
                  write(*,'(4es20.8,2l4)') dels_u, dels_xzu, dels_xyu, dels_yzu, dels_u.eq.dels_xzu, dels_xyu.eq.dels_u
                  stop 'error in fsc_linec3d: invalid dels_u'
               endif
!
!---------------------------downwind point------------------------------
!
               if(dels_xyd.eq.dels_d) then
!intersection with x-y plane on level k+gamma
                  x_d = x(i) + dels_d*nn_x
                  y_d = y(j) + dels_d*nn_y
                  z_d = z(kkp1)
!
                  call coeff2d_linecd(   opac3d(iim1,jjm1,kkp1),    opac3d(i,jjm1,kkp1),    opac3d(iip1,jjm1,kkp1), &
                     opac3d(iim1,j,kkp1),       opac3d(i,j,kkp1),       opac3d(iip1,j,kkp1), &
                     opac3d(iim1,jjp1,kkp1),    opac3d(i,jjp1,kkp1),    opac3d(iip1,jjp1,kkp1), &
                     opalbar3d(iim1,jjm1,kkp1), opalbar3d(i,jjm1,kkp1), opalbar3d(iip1,jjm1,kkp1), &
                     opalbar3d(iim1,j,kkp1),    opalbar3d(i,j,kkp1),    opalbar3d(iip1,j,kkp1), &
                     opalbar3d(iim1,jjp1,kkp1), opalbar3d(i,jjp1,kkp1), opalbar3d(iip1,jjp1,kkp1), &
                     velx3d(iim1,jjm1,kkp1),    velx3d(i,jjm1,kkp1),    velx3d(iip1,jjm1,kkp1), &
                     velx3d(iim1,j,kkp1),       velx3d(i,j,kkp1),       velx3d(iip1,j,kkp1), &
                     velx3d(iim1,jjp1,kkp1),    velx3d(i,jjp1,kkp1),    velx3d(iip1,jjp1,kkp1), &
                     vely3d(iim1,jjm1,kkp1),    vely3d(i,jjm1,kkp1),    vely3d(iip1,jjm1,kkp1), &
                     vely3d(iim1,j,kkp1),       vely3d(i,j,kkp1),       vely3d(iip1,j,kkp1), &
                     vely3d(iim1,jjp1,kkp1),    vely3d(i,jjp1,kkp1),    vely3d(iip1,jjp1,kkp1), &
                     velz3d(iim1,jjm1,kkp1),    velz3d(i,jjm1,kkp1),    velz3d(iip1,jjm1,kkp1), &
                     velz3d(iim1,j,kkp1),       velz3d(i,j,kkp1),       velz3d(iip1,j,kkp1), &
                     velz3d(iim1,jjp1,kkp1),    velz3d(i,jjp1,kkp1),    velz3d(iip1,jjp1,kkp1), &
                     vth3d(iim1,jjm1,kkp1),      vth3d(i,jjm1,kkp1),    vth3d(iip1,jjm1,kkp1), &
                     vth3d(iim1,j,kkp1),         vth3d(i,j,kkp1),       vth3d(iip1,j,kkp1), &
                     vth3d(iim1,jjp1,kkp1),      vth3d(i,jjp1,kkp1),    vth3d(iip1,jjp1,kkp1), &
                     scont3d(iim1,jjm1,kkp1),    scont3d(i,jjm1,kkp1),  scont3d(iip1,jjm1,kkp1), &
                     scont3d(iim1,j,kkp1),       scont3d(i,j,kkp1),     scont3d(iip1,j,kkp1), &
                     scont3d(iim1,jjp1,kkp1),    scont3d(i,jjp1,kkp1),  scont3d(iip1,jjp1,kkp1), &
                     sline3d(iim1,jjm1,kkp1),    sline3d(i,jjm1,kkp1),  sline3d(iip1,jjm1,kkp1), &
                     sline3d(iim1,j,kkp1),       sline3d(i,j,kkp1),     sline3d(iip1,j,kkp1), &
                     sline3d(iim1,jjp1,kkp1),    sline3d(i,jjp1,kkp1),  sline3d(iip1,jjp1,kkp1), &
                     x(iim1), x(i), x(iip1), y(jjm1), y(j), y(jjp1), x_d, y_d, &
                     c19_slined, c20_slined, c21_slined, c22_slined, c23_slined, c24_slined, &
                     c25_slined, c26_slined, c27_slined, &
                     opac_d, opalbar_d, velx_d, vely_d, velz_d, vth_d, scont_d, sline_d)
!set interpolation coefficients that are not used to zero
                  c03_slined = 0.d0
                  c06_slined = 0.d0
                  c07_slined = 0.d0
                  c08_slined = 0.d0
                  c09_slined = 0.d0
                  c12_slined = 0.d0
                  c15_slined = 0.d0
                  c16_slined = 0.d0
                  c17_slined = 0.d0
                  c18_slined = 0.d0
!
               elseif(dels_xzd.eq.dels_d) then
!intersection with x-z plane on level j+beta
                  x_d = x(i) + dels_d*nn_x
                  y_d = y(jjp1)
                  z_d = z(k) + dels_d*nn_z
!
                  call coeff2d_linecd(   opac3d(iim1,jjp1,kkm1),    opac3d(i,jjp1,kkm1),    opac3d(iip1,jjp1,kkm1), &
                     opac3d(iim1,jjp1,k),       opac3d(i,jjp1,k),       opac3d(iip1,jjp1,k), &
                     opac3d(iim1,jjp1,kkp1),    opac3d(i,jjp1,kkp1),    opac3d(iip1,jjp1,kkp1), &
                     opalbar3d(iim1,jjp1,kkm1), opalbar3d(i,jjp1,kkm1), opalbar3d(iip1,jjp1,kkm1), &
                     opalbar3d(iim1,jjp1,k),    opalbar3d(i,jjp1,k),    opalbar3d(iip1,jjp1,k), &
                     opalbar3d(iim1,jjp1,kkp1), opalbar3d(i,jjp1,kkp1), opalbar3d(iip1,jjp1,kkp1), &
                     velx3d(iim1,jjp1,kkm1),    velx3d(i,jjp1,kkm1),    velx3d(iip1,jjp1,kkm1), &
                     velx3d(iim1,jjp1,k),       velx3d(i,jjp1,k),       velx3d(iip1,jjp1,k), &
                     velx3d(iim1,jjp1,kkp1),    velx3d(i,jjp1,kkp1),    velx3d(iip1,jjp1,kkp1), &
                     vely3d(iim1,jjp1,kkm1),    vely3d(i,jjp1,kkm1),    vely3d(iip1,jjp1,kkm1), &
                     vely3d(iim1,jjp1,k),       vely3d(i,jjp1,k),       vely3d(iip1,jjp1,k), &
                     vely3d(iim1,jjp1,kkp1),    vely3d(i,jjp1,kkp1),    vely3d(iip1,jjp1,kkp1), &
                     velz3d(iim1,jjp1,kkm1),    velz3d(i,jjp1,kkm1),    velz3d(iip1,jjp1,kkm1), &
                     velz3d(iim1,jjp1,k),       velz3d(i,jjp1,k),       velz3d(iip1,jjp1,k), &
                     velz3d(iim1,jjp1,kkp1),    velz3d(i,jjp1,kkp1),    velz3d(iip1,jjp1,kkp1), &
                     vth3d(iim1,jjp1,kkm1),     vth3d(i,jjp1,kkm1),     vth3d(iip1,jjp1,kkm1), &
                     vth3d(iim1,jjp1,k),        vth3d(i,jjp1,k),        vth3d(iip1,jjp1,k), &
                     vth3d(iim1,jjp1,kkp1),     vth3d(i,jjp1,kkp1),     vth3d(iip1,jjp1,kkp1), &
                     scont3d(iim1,jjp1,kkm1),   scont3d(i,jjp1,kkm1),   scont3d(iip1,jjp1,kkm1), &
                     scont3d(iim1,jjp1,k),      scont3d(i,jjp1,k),      scont3d(iip1,jjp1,k), &
                     scont3d(iim1,jjp1,kkp1),   scont3d(i,jjp1,kkp1),   scont3d(iip1,jjp1,kkp1), &
                     sline3d(iim1,jjp1,kkm1),   sline3d(i,jjp1,kkm1),   sline3d(iip1,jjp1,kkm1), &
                     sline3d(iim1,jjp1,k),      sline3d(i,jjp1,k),      sline3d(iip1,jjp1,k), &
                     sline3d(iim1,jjp1,kkp1),   sline3d(i,jjp1,kkp1),   sline3d(iip1,jjp1,kkp1), &
                     x(iim1), x(i), x(iip1), z(kkm1), z(k), z(kkp1), x_d, z_d, &
                     c07_slined, c08_slined, c09_slined, c16_slined, c17_slined, c18_slined, &
                     c25_slined, c26_slined, c27_slined, &
                     opac_d, opalbar_d, velx_d, vely_d, velz_d, vth_d, scont_d, sline_d)

!                 if(abs(c07_slined*sline3d(iim1,jjp1,kkm1)+c08_slined*sline3d(i,jjp1,kkm1)+c09_slined*sline3d(iip1,jjp1,kkm1)+&
!                        c16_slined*sline3d(iim1,jjp1,k   )+c17_slined*sline3d(i,jjp1,k   )+c18_slined*sline3d(iip1,jjp1,k   )+&
!                        c25_slined*sline3d(iim1,jjp1,kkp1)+c26_slined*sline3d(i,jjp1,kkp1)+c27_slined*sline3d(iip1,jjp1,kkp1)-sline_d).gt.1.d-14) stop 'error coeff'
!
!set interpolation coefficients that are not used to zero
                  c03_slined = 0.d0
                  c06_slined = 0.d0
                  c12_slined = 0.d0
                  c15_slined = 0.d0
                  c19_slined = 0.d0
                  c20_slined = 0.d0
                  c21_slined = 0.d0
                  c22_slined = 0.d0
                  c23_slined = 0.d0
                  c24_slined = 0.d0
!
               elseif(dels_yzd.eq.dels_d) then
!intersection with y-z plane on level i+alpha
                  x_d = x(iip1)
                  y_d = y(j) + dels_d*nn_y
                  z_d = z(k) + dels_d*nn_z
!
                  call coeff2d_linecd(   opac3d(iip1,jjm1,kkm1),    opac3d(iip1,j,kkm1),    opac3d(iip1,jjp1,kkm1), &
                     opac3d(iip1,jjm1,k),       opac3d(iip1,j,k),       opac3d(iip1,jjp1,k), &
                     opac3d(iip1,jjm1,kkp1),    opac3d(iip1,j,kkp1),    opac3d(iip1,jjp1,kkp1), &
                     opalbar3d(iip1,jjm1,kkm1), opalbar3d(iip1,j,kkm1), opalbar3d(iip1,jjp1,kkm1), &
                     opalbar3d(iip1,jjm1,k),    opalbar3d(iip1,j,k),    opalbar3d(iip1,jjp1,k), &
                     opalbar3d(iip1,jjm1,kkp1), opalbar3d(iip1,j,kkp1), opalbar3d(iip1,jjp1,kkp1), &
                     velx3d(iip1,jjm1,kkm1),    velx3d(iip1,j,kkm1),    velx3d(iip1,jjp1,kkm1), &
                     velx3d(iip1,jjm1,k),       velx3d(iip1,j,k),       velx3d(iip1,jjp1,k), &
                     velx3d(iip1,jjm1,kkp1),    velx3d(iip1,j,kkp1),    velx3d(iip1,jjp1,kkp1), &
                     vely3d(iip1,jjm1,kkm1),    vely3d(iip1,j,kkm1),    vely3d(iip1,jjp1,kkm1), &
                     vely3d(iip1,jjm1,k),       vely3d(iip1,j,k),       vely3d(iip1,jjp1,k), &
                     vely3d(iip1,jjm1,kkp1),    vely3d(iip1,j,kkp1),    vely3d(iip1,jjp1,kkp1), &
                     velz3d(iip1,jjm1,kkm1),    velz3d(iip1,j,kkm1),    velz3d(iip1,jjp1,kkm1), &
                     velz3d(iip1,jjm1,k),       velz3d(iip1,j,k),       velz3d(iip1,jjp1,k), &
                     velz3d(iip1,jjm1,kkp1),    velz3d(iip1,j,kkp1),    velz3d(iip1,jjp1,kkp1), &
                     vth3d(iip1,jjm1,kkm1),     vth3d(iip1,j,kkm1),     vth3d(iip1,jjp1,kkm1), &
                     vth3d(iip1,jjm1,k),        vth3d(iip1,j,k),        vth3d(iip1,jjp1,k), &
                     vth3d(iip1,jjm1,kkp1),     vth3d(iip1,j,kkp1),     vth3d(iip1,jjp1,kkp1), &
                     scont3d(iip1,jjm1,kkm1),   scont3d(iip1,j,kkm1),   scont3d(iip1,jjp1,kkm1), &
                     scont3d(iip1,jjm1,k),      scont3d(iip1,j,k),      scont3d(iip1,jjp1,k), &
                     scont3d(iip1,jjm1,kkp1),   scont3d(iip1,j,kkp1),   scont3d(iip1,jjp1,kkp1), &
                     sline3d(iip1,jjm1,kkm1),   sline3d(iip1,j,kkm1),   sline3d(iip1,jjp1,kkm1), &
                     sline3d(iip1,jjm1,k),      sline3d(iip1,j,k),      sline3d(iip1,jjp1,k), &
                     sline3d(iip1,jjm1,kkp1),   sline3d(iip1,j,kkp1),   sline3d(iip1,jjp1,kkp1), &
                     y(jjm1), y(j), y(jjp1), z(kkm1), z(k), z(kkp1), y_d, z_d, &
                     c03_slined, c06_slined, c09_slined, c12_slined, c15_slined, c18_slined, &
                     c21_slined, c24_slined, c27_slined, &
                     opac_d, opalbar_d, velx_d, vely_d, velz_d, vth_d, scont_d, sline_d)
!set interpolation coefficients that are not used to zero
                  c07_slined = 0.d0
                  c08_slined = 0.d0
                  c16_slined = 0.d0
                  c17_slined = 0.d0
                  c19_slined = 0.d0
                  c20_slined = 0.d0
                  c22_slined = 0.d0
                  c23_slined = 0.d0
                  c25_slined = 0.d0
                  c26_slined = 0.d0
!
               else
                  write(*,'(4es20.8,l4)') dels_d, dels_xzd, dels_xyd, dels_yzd, dels_d.eq.dels_xzd
                  stop 'error in fsc_linec3d: invalid dels_d'
               endif
!
!--------------------------------radiative transfer---------------------
!
               vel_u = velx_u*nn_x + vely_u*nn_y + velz_u*nn_z
               vel_d = velx_d*nn_x + vely_d*nn_y + velz_d*nn_z
               call fsc_linec(xobs, vth_fiducial, deltax, xcmf_max, xcmf_min, int_u, opac_u, opac_p, opac_d, &
                  opalbar_u, opalbar_p, opalbar_d, scont_u, scont_p, scont_d, &
                  sline_u, sline_p, sline_d, vel_u, vel_p, vel_d, vth_u, vth_p, vth_d, &
                  dels_u, dels_d, abs_sc, contr_sc, int_sc, alo_u, alo_p, alo_d)!
               int3d(i,j,k) = int_sc
!
               aloline_on_nn3d(iip1,jjp1,kkp1,q1) = imask_totreg3d(iip1,jjp1,kkp1)*&
                  alo_d*c27_slined
               aloline_on_nn3d(i,jjp1,kkp1,q2) = imask_totreg3d(i,jjp1,kkp1) * &
                  (alo_d*c26_slined + abs_sc*c26_intu*aloline_on_nn3d(i,jjp1,kkp1,q1))
               aloline_on_nn3d(iim1,jjp1,kkp1,q3) = imask_totreg3d(iim1,jjp1,kkp1) * &
                  (alo_d*c25_slined + abs_sc*c26_intu*aloline_on_nn3d(iim1,jjp1,kkp1,q2))
               aloline_on_nn3d(iip1,j,kkp1,q4) = imask_totreg3d(iip1,j,kkp1) * &
                  (alo_d*c24_slined + abs_sc*c24_intu*aloline_on_nn3d(iip1,j,kkp1,q1))
               aloline_on_nn3d(i,j,kkp1,q5) = imask_totreg3d(i,j,kkp1) * &
                  (alo_d*c23_slined + abs_sc*(c23_intu*aloline_on_nn3d(i,j,kkp1,q1) + &
                  c24_intu*aloline_on_nn3d(i,j,kkp1,q2) + &
                  c26_intu*aloline_on_nn3d(i,j,kkp1,q4)))
               aloline_on_nn3d(iim1,j,kkp1,q6) = imask_totreg3d(iim1,j,kkp1) * &
                  (alo_d*c22_slined + abs_sc*(c22_intu*aloline_on_nn3d(iim1,j,kkp1,q1) + &
                  c23_intu*aloline_on_nn3d(iim1,j,kkp1,q2) + &
                  c24_intu*aloline_on_nn3d(iim1,j,kkp1,q3) + &
                  c26_intu*aloline_on_nn3d(iim1,j,kkp1,q5)))
               aloline_on_nn3d(iip1,jjm1,kkp1,q7) = imask_totreg3d(iip1,jjm1,kkp1) * &
                  (alo_d*c21_slined + abs_sc*c24_intu*aloline_on_nn3d(iip1,jjm1,kkp1,q4))
               aloline_on_nn3d(i,jjm1,kkp1,q8) = imask_totreg3d(i,jjm1,kkp1) * &
                  (alo_d*c20_slined + abs_sc*(c23_intu*aloline_on_nn3d(i,jjm1,kkp1,q4) + &
                  c24_intu*aloline_on_nn3d(i,jjm1,kkp1,q5) + &
                  c20_intu*aloline_on_nn3d(i,jjm1,kkp1,q1) + &
                  c26_intu*aloline_on_nn3d(i,jjm1,kkp1,q7)))
               aloline_on_nn3d(iim1,jjm1,kkp1,q9) = imask_totreg3d(iim1,jjm1,kkp1) * &
                  (alo_d*c19_slined + abs_sc*(c22_intu*aloline_on_nn3d(iim1,jjm1,kkp1,q4) + &
                  c23_intu*aloline_on_nn3d(iim1,jjm1,kkp1,q5) + &
                  c24_intu*aloline_on_nn3d(iim1,jjm1,kkp1,q6) + &
                  c20_intu*aloline_on_nn3d(iim1,jjm1,kkp1,q2) + &
                  c26_intu*aloline_on_nn3d(iim1,jjm1,kkp1,q8)))
               aloline_on_nn3d(iip1,jjp1,k,q10) = imask_totreg3d(iip1,jjp1,k) * &
                  (alo_d*c18_slined + abs_sc*c18_intu*aloline_on_nn3d(iip1,jjp1,k,q1))
               aloline_on_nn3d(i,jjp1,k,q11) = imask_totreg3d(i,jjp1,k) * &
                  (alo_d*c17_slined + abs_sc*(c17_intu*aloline_on_nn3d(i,jjp1,k,q1) + &
                  c18_intu*aloline_on_nn3d(i,jjp1,k,q2) + &
                  c26_intu*aloline_on_nn3d(i,jjp1,k,q10)))
               aloline_on_nn3d(iim1,jjp1,k,q12) = imask_totreg3d(iim1,jjp1,k) * &
                  (alo_d*c16_slined + abs_sc*(c16_intu*aloline_on_nn3d(iim1,jjp1,k,q1) + &
                  c17_intu*aloline_on_nn3d(iim1,jjp1,k,q2) + &
                  c18_intu*aloline_on_nn3d(iim1,jjp1,k,q3) + &
                  c26_intu*aloline_on_nn3d(iim1,jjp1,k,q11)))
               aloline_on_nn3d(iip1,j,k,q13) = imask_totreg3d(iip1,j,k) * &
                  (alo_d*c15_slined + abs_sc*(c15_intu*aloline_on_nn3d(iip1,j,k,q1) + &
                  c18_intu*aloline_on_nn3d(iip1,j,k,q4) + &
                  c24_intu*aloline_on_nn3d(iip1,j,k,q10)))
               aloline_on_nn3d(i,j,k,q14) = imask_totreg3d(i,j,k) * &
                  (alo_p + abs_sc*(c14_intu*aloline_on_nn3d(i,j,k,q1) + &
                  c15_intu*aloline_on_nn3d(i,j,k,q2) + &
                  c17_intu*aloline_on_nn3d(i,j,k,q4) + &
                  c26_intu*aloline_on_nn3d(i,j,k,q13) + &
                  c18_intu*aloline_on_nn3d(i,j,k,q5) + &
                  c23_intu*aloline_on_nn3d(i,j,k,q10) + &
                  c24_intu*aloline_on_nn3d(i,j,k,q11)))
               aloline_on_nn3d(iim1,j,k,q15) = imask_totreg3d(iim1,j,k) * &
                  (alo_u*c26_slineu + abs_sc*(c13_intu*aloline_on_nn3d(iim1,j,k,q1) + &
                  c14_intu*aloline_on_nn3d(iim1,j,k,q2) + &
                  c15_intu*aloline_on_nn3d(iim1,j,k,q3) + &
                  c16_intu*aloline_on_nn3d(iim1,j,k,q4) + &
                  c17_intu*aloline_on_nn3d(iim1,j,k,q5) + &
                  c18_intu*aloline_on_nn3d(iim1,j,k,q6) + &
                  c22_intu*aloline_on_nn3d(iim1,j,k,q10) + &
                  c23_intu*aloline_on_nn3d(iim1,j,k,q11) + &
                  c24_intu*aloline_on_nn3d(iim1,j,k,q12) + &
                  c26_intu*aloline_on_nn3d(iim1,j,k,q14)))
               aloline_on_nn3d(iip1,jjm1,k,q16) = imask_totreg3d(iip1,jjm1,k) * &
                  (alo_d*c12_slined + abs_sc*(c12_intu*aloline_on_nn3d(iip1,jjm1,k,q1) + &
                  c15_intu*aloline_on_nn3d(iip1,jjm1,k,q4) + &
                  c18_intu*aloline_on_nn3d(iip1,jjm1,k,q7) + &
                  c24_intu*aloline_on_nn3d(iip1,jjm1,k,q13)))
               aloline_on_nn3d(i,jjm1,k,q17) = imask_totreg3d(i,jjm1,k) * &
                  (alo_u*c24_slineu + abs_sc*(c11_intu*aloline_on_nn3d(i,jjm1,k,q1) + &
                  c12_intu*aloline_on_nn3d(i,jjm1,k,q2) + &
                  c14_intu*aloline_on_nn3d(i,jjm1,k,q4) + &
                  c15_intu*aloline_on_nn3d(i,jjm1,k,q5) + &
                  c17_intu*aloline_on_nn3d(i,jjm1,k,q7) + &
                  c18_intu*aloline_on_nn3d(i,jjm1,k,q8) + &
                  c23_intu*aloline_on_nn3d(i,jjm1,k,q13) + &
                  c24_intu*aloline_on_nn3d(i,jjm1,k,q14) + &
                  c20_intu*aloline_on_nn3d(i,jjm1,k,q10) + &
                  c26_intu*aloline_on_nn3d(i,jjm1,k,q16)))
               aloline_on_nn3d(iim1,jjm1,k,q18) = imask_totreg3d(iim1,jjm1,k) * &
                  (alo_u*c23_slineu + abs_sc*(c10_intu*aloline_on_nn3d(iim1,jjm1,k,q1) + &
                  c11_intu*aloline_on_nn3d(iim1,jjm1,k,q2) + &
                  c12_intu*aloline_on_nn3d(iim1,jjm1,k,q3) + &
                  c13_intu*aloline_on_nn3d(iim1,jjm1,k,q4) + &
                  c14_intu*aloline_on_nn3d(iim1,jjm1,k,q5) + &
                  c15_intu*aloline_on_nn3d(iim1,jjm1,k,q6) + &
                  c16_intu*aloline_on_nn3d(iim1,jjm1,k,q7) + &
                  c17_intu*aloline_on_nn3d(iim1,jjm1,k,q8) + &
                  c18_intu*aloline_on_nn3d(iim1,jjm1,k,q9) + &
                  c22_intu*aloline_on_nn3d(iim1,jjm1,k,q13) + &
                  c23_intu*aloline_on_nn3d(iim1,jjm1,k,q14) + &
                  c24_intu*aloline_on_nn3d(iim1,jjm1,k,q15) + &
                  c20_intu*aloline_on_nn3d(iim1,jjm1,k,q11) + &
                  c26_intu*aloline_on_nn3d(iim1,jjm1,k,q17)))
               aloline_on_nn3d(iip1,jjp1,kkm1,q19) = imask_totreg3d(iip1,jjp1,kkm1) * &
                  (alo_d*c09_slined + abs_sc*c18_intu*aloline_on_nn3d(iip1,jjp1,kkm1,q10))
               aloline_on_nn3d(i,jjp1,kkm1,q20) = imask_totreg3d(i,jjp1,kkm1) * &
                  (alo_d*c08_slined + abs_sc*(c17_intu*aloline_on_nn3d(i,jjp1,kkm1,q10) + &
                  c18_intu*aloline_on_nn3d(i,jjp1,kkm1,q11) + &
                  c08_intu*aloline_on_nn3d(i,jjp1,kkm1,q1) + &
                  c26_intu*aloline_on_nn3d(i,jjp1,kkm1,q19)))
               aloline_on_nn3d(iim1,jjp1,kkm1,q21) = imask_totreg3d(iim1,jjp1,kkm1) * &
                  (alo_d*c07_slined + abs_sc*(c16_intu*aloline_on_nn3d(iim1,jjp1,kkm1,q10) + &
                  c17_intu*aloline_on_nn3d(iim1,jjp1,kkm1,q11) + &
                  c18_intu*aloline_on_nn3d(iim1,jjp1,kkm1,q12) + &
                  c08_intu*aloline_on_nn3d(iim1,jjp1,kkm1,q2) + &
                  c26_intu*aloline_on_nn3d(iim1,jjp1,kkm1,q20)))
               aloline_on_nn3d(iip1,j,kkm1,q22) = imask_totreg3d(iip1,j,kkm1) * &
                  (alo_d*c06_slined + abs_sc*(c15_intu*aloline_on_nn3d(iip1,j,kkm1,q10) + &
                  c18_intu*aloline_on_nn3d(iip1,j,kkm1,q13) + &
                  c06_intu*aloline_on_nn3d(iip1,j,kkm1,q1) + &
                  c24_intu*aloline_on_nn3d(iip1,j,kkm1,q19)))
               aloline_on_nn3d(i,j,kkm1,q23) = imask_totreg3d(i,j,kkm1) * &
                  (alo_u*c18_slineu + abs_sc*(c14_intu*aloline_on_nn3d(i,j,kkm1,q10) + &
                  c15_intu*aloline_on_nn3d(i,j,kkm1,q11) + &
                  c17_intu*aloline_on_nn3d(i,j,kkm1,q13) + &
                  c18_intu*aloline_on_nn3d(i,j,kkm1,q14) + &
                  c05_intu*aloline_on_nn3d(i,j,kkm1,q1) + &
                  c06_intu*aloline_on_nn3d(i,j,kkm1,q2) + &
                  c23_intu*aloline_on_nn3d(i,j,kkm1,q19) + &
                  c24_intu*aloline_on_nn3d(i,j,kkm1,q20) + &
                  c08_intu*aloline_on_nn3d(i,j,kkm1,q4) + &
                  c26_intu*aloline_on_nn3d(i,j,kkm1,q22)))
               aloline_on_nn3d(iim1,j,kkm1,q24) = imask_totreg3d(iim1,j,kkm1) * &
                  (alo_u*c17_slineu + abs_sc*(c13_intu*aloline_on_nn3d(iim1,j,kkm1,q10) + &
                  c14_intu*aloline_on_nn3d(iim1,j,kkm1,q11) + &
                  c15_intu*aloline_on_nn3d(iim1,j,kkm1,q12) + &
                  c16_intu*aloline_on_nn3d(iim1,j,kkm1,q13) + &
                  c17_intu*aloline_on_nn3d(iim1,j,kkm1,q14) + &
                  c18_intu*aloline_on_nn3d(iim1,j,kkm1,q15) + &
                  c04_intu*aloline_on_nn3d(iim1,j,kkm1,q1) + &
                  c05_intu*aloline_on_nn3d(iim1,j,kkm1,q2) + &
                  c06_intu*aloline_on_nn3d(iim1,j,kkm1,q3) + &
                  c22_intu*aloline_on_nn3d(iim1,j,kkm1,q19) + &
                  c23_intu*aloline_on_nn3d(iim1,j,kkm1,q20) + &
                  c24_intu*aloline_on_nn3d(iim1,j,kkm1,q21) + &
                  c08_intu*aloline_on_nn3d(iim1,j,kkm1,q5) + &
                  c26_intu*aloline_on_nn3d(iim1,j,kkm1,q23)))
               aloline_on_nn3d(iip1,jjm1,kkm1,q25) = imask_totreg3d(iip1,jjm1,kkm1) * &
                  (alo_d*c03_slined + abs_sc*(c12_intu*aloline_on_nn3d(iip1,jjm1,kkm1,q10) + &
                  c15_intu*aloline_on_nn3d(iip1,jjm1,kkm1,q13) + &
                  c18_intu*aloline_on_nn3d(iip1,jjm1,kkm1,q16) + &
                  c06_intu*aloline_on_nn3d(iip1,jjm1,kkm1,q4) + &
                  c24_intu*aloline_on_nn3d(iip1,jjm1,kkm1,q22)))
               aloline_on_nn3d(i,jjm1,kkm1,q26) = imask_totreg3d(i,jjm1,kkm1) * &
                  (alo_u*c15_slineu + abs_sc*(c11_intu*aloline_on_nn3d(i,jjm1,kkm1,q10) + &
                  c12_intu*aloline_on_nn3d(i,jjm1,kkm1,q11) + &
                  c14_intu*aloline_on_nn3d(i,jjm1,kkm1,q13) + &
                  c15_intu*aloline_on_nn3d(i,jjm1,kkm1,q14) + &
                  c17_intu*aloline_on_nn3d(i,jjm1,kkm1,q16) + &
                  c18_intu*aloline_on_nn3d(i,jjm1,kkm1,q17) + &
                  c05_intu*aloline_on_nn3d(i,jjm1,kkm1,q4) + &
                  c06_intu*aloline_on_nn3d(i,jjm1,kkm1,q5) + &
                  c23_intu*aloline_on_nn3d(i,jjm1,kkm1,q22) + &
                  c24_intu*aloline_on_nn3d(i,jjm1,kkm1,q23) + &
                  c02_intu*aloline_on_nn3d(i,jjm1,kkm1,q1) + &
                  c08_intu*aloline_on_nn3d(i,jjm1,kkm1,q7) + &
                  c20_intu*aloline_on_nn3d(i,jjm1,kkm1,q19) + &
                  c26_intu*aloline_on_nn3d(i,jjm1,kkm1,q25)))
               aloline_on_nn3d(iim1,jjm1,kkm1,q27) = imask_totreg3d(iim1,jjm1,kkm1) * &
                  (alo_u*c14_slineu + abs_sc*(c10_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q10) + &
                  c11_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q11) + &
                  c12_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q12) + &
                  c13_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q13) + &
                  c14_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q14) + &
                  c15_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q15) + &
                  c16_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q16) + &
                  c17_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q17) + &
                  c18_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q18) + &
                  c04_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q4) + &
                  c05_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q5) + &
                  c06_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q6) + &
                  c22_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q22) + &
                  c23_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q23) + &
                  c24_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q24) + &
                  c02_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q2) + &
                  c08_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q8) + &
                  c20_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q20) + &
                  c26_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q26)))

!
!perform angular integration
               call calc_phinorm(vel_p, vth_p, vth_fiducial, xobs, phinorm)
               wall=wall_global*phinorm
               mintbar3d_tmp(i,j,k) = mintbar3d_tmp(i,j,k) + int3d(i,j,k)*wall
               normalization3d_tmp(i,j,k) = normalization3d_tmp(i,j,k) + wall
!
               aloline_nn3d_tmp(iip1,jjp1,kkp1,q1) = aloline_nn3d_tmp(iip1,jjp1,kkp1,q1) + wall*aloline_on_nn3d(iip1,jjp1,kkp1,q1)
               aloline_nn3d_tmp(i,jjp1,kkp1,q2) = aloline_nn3d_tmp(i,jjp1,kkp1,q2) + wall*aloline_on_nn3d(i,jjp1,kkp1,q2)
               aloline_nn3d_tmp(iim1,jjp1,kkp1,q3) = aloline_nn3d_tmp(iim1,jjp1,kkp1,q3) + wall*aloline_on_nn3d(iim1,jjp1,kkp1,q3)
               aloline_nn3d_tmp(iip1,j,kkp1,q4) = aloline_nn3d_tmp(iip1,j,kkp1,q4) + wall*aloline_on_nn3d(iip1,j,kkp1,q4)
               aloline_nn3d_tmp(i,j,kkp1,q5) = aloline_nn3d_tmp(i,j,kkp1,q5) + wall*aloline_on_nn3d(i,j,kkp1,q5)
               aloline_nn3d_tmp(iim1,j,kkp1,q6) = aloline_nn3d_tmp(iim1,j,kkp1,q6) + wall*aloline_on_nn3d(iim1,j,kkp1,q6)
               aloline_nn3d_tmp(iip1,jjm1,kkp1,q7) = aloline_nn3d_tmp(iip1,jjm1,kkp1,q7) + wall*aloline_on_nn3d(iip1,jjm1,kkp1,q7)
               aloline_nn3d_tmp(i,jjm1,kkp1,q8) = aloline_nn3d_tmp(i,jjm1,kkp1,q8) + wall*aloline_on_nn3d(i,jjm1,kkp1,q8)
               aloline_nn3d_tmp(iim1,jjm1,kkp1,q9) = aloline_nn3d_tmp(iim1,jjm1,kkp1,q9) + wall*aloline_on_nn3d(iim1,jjm1,kkp1,q9)
               aloline_nn3d_tmp(iip1,jjp1,k,q10) = aloline_nn3d_tmp(iip1,jjp1,k,q10) + wall*aloline_on_nn3d(iip1,jjp1,k,q10)
               aloline_nn3d_tmp(i,jjp1,k,q11) = aloline_nn3d_tmp(i,jjp1,k,q11) + wall*aloline_on_nn3d(i,jjp1,k,q11)
               aloline_nn3d_tmp(iim1,jjp1,k,q12) = aloline_nn3d_tmp(iim1,jjp1,k,q12) + wall*aloline_on_nn3d(iim1,jjp1,k,q12)
               aloline_nn3d_tmp(iip1,j,k,q13) = aloline_nn3d_tmp(iip1,j,k,q13) + wall*aloline_on_nn3d(iip1,j,k,q13)
               aloline_nn3d_tmp(i,j,k,q14) = aloline_nn3d_tmp(i,j,k,q14) + wall*aloline_on_nn3d(i,j,k,q14)
               aloline_nn3d_tmp(iim1,j,k,q15) = aloline_nn3d_tmp(iim1,j,k,q15) + wall*aloline_on_nn3d(iim1,j,k,q15)
               aloline_nn3d_tmp(iip1,jjm1,k,q16) = aloline_nn3d_tmp(iip1,jjm1,k,q16) + wall*aloline_on_nn3d(iip1,jjm1,k,q16) !new
               aloline_nn3d_tmp(i,jjm1,k,q17) = aloline_nn3d_tmp(i,jjm1,k,q17) + wall*aloline_on_nn3d(i,jjm1,k,q17)
               aloline_nn3d_tmp(iim1,jjm1,k,q18) = aloline_nn3d_tmp(iim1,jjm1,k,q18) + wall*aloline_on_nn3d(iim1,jjm1,k,q18)
               aloline_nn3d_tmp(iip1,jjp1,kkm1,q19) = aloline_nn3d_tmp(iip1,jjp1,kkm1,q19) + wall*aloline_on_nn3d(iip1,jjp1,kkm1,q19) !new
               aloline_nn3d_tmp(i,jjp1,kkm1,q20) = aloline_nn3d_tmp(i,jjp1,kkm1,q20) + wall*aloline_on_nn3d(i,jjp1,kkm1,q20) !new
               aloline_nn3d_tmp(iim1,jjp1,kkm1,q21) = aloline_nn3d_tmp(iim1,jjp1,kkm1,q21) + wall*aloline_on_nn3d(iim1,jjp1,kkm1,q21) !new
               aloline_nn3d_tmp(iip1,j,kkm1,q22) = aloline_nn3d_tmp(iip1,j,kkm1,q22) + wall*aloline_on_nn3d(iip1,j,kkm1,q22) !new
               aloline_nn3d_tmp(i,j,kkm1,q23) = aloline_nn3d_tmp(i,j,kkm1,q23) + wall*aloline_on_nn3d(i,j,kkm1,q23)
               aloline_nn3d_tmp(iim1,j,kkm1,q24) = aloline_nn3d_tmp(iim1,j,kkm1,q24) + wall*aloline_on_nn3d(iim1,j,kkm1,q24)
               aloline_nn3d_tmp(iip1,jjm1,kkm1,q25) = aloline_nn3d_tmp(iip1,jjm1,kkm1,q25) + wall*aloline_on_nn3d(iip1,jjm1,kkm1,q25) !new
               aloline_nn3d_tmp(i,jjm1,kkm1,q26) = aloline_nn3d_tmp(i,jjm1,kkm1,q26) + wall*aloline_on_nn3d(i,jjm1,kkm1,q26)
               aloline_nn3d_tmp(iim1,jjm1,kkm1,q27) = aloline_nn3d_tmp(iim1,jjm1,kkm1,q27) + wall*aloline_on_nn3d(iim1,jjm1,kkm1,q27)
!
!***********special treatment: point is in vicinity of boundary*********
!
             case(2)
!
               iim1=i-alpha
               jjm1=j-beta
               kkm1=k-gamma
               iip1=i+alpha
               jjp1=j+beta
               kkp1=k+gamma
!
!calculate distance of upwind point to (previous) xy-plane, xz-plane and yz-plane
               dels_xyu=(z(k)-z(kkm1))/nn_z
               dels_xzu=(y(j)-y(jjm1))/nn_y
               dels_yzu=(x(i)-x(iim1))/nn_x
!
!calculate distance of downwind point to (next) xy-plane, xz-plane and yz-plane
               dels_xyd=(z(kkp1)-z(k))/nn_z
               dels_xzd=(y(jjp1)-y(j))/nn_y
               dels_yzd=(x(iip1)-x(i))/nn_x
!
!calculate distance to the boundary (here: sphere)
!               dels_r=dist_sphere(nn_x,nn_y,nn_z,x(i),y(j),z(k))
               dels_r=dist_ellipsoid(nn_x,nn_y,nn_z,x(i),y(j),z(k),smajorax_a,smajorax_b,smajorax_c)
!
               dels_u=min(dels_xyu,dels_xzu,dels_yzu,dels_r)
               dels_d=min(dels_xyd,dels_xzd,dels_yzd)
!
!----------------------------local point--------------------------------
!
               sline_p=sline3d(i,j,k)
               scont_p=scont3d(i,j,k)
               opac_p=opac3d(i,j,k)
               opalbar_p=opalbar3d(i,j,k)
               velx_p=velx3d(i,j,k)
               vely_p=vely3d(i,j,k)
               velz_p=velz3d(i,j,k)
               vel_p = nn_x*velx_p + nn_y*vely_p + nn_z*velz_p
               vth_p = vth3d(i,j,k)
!
!----------------------------upwind point-------------------------------
!
               if(dels_r.eq.dels_u) then
                  x_u = x(i) - dels_u*nn_x
                  y_u = y(j) - dels_u*nn_y
                  z_u = z(k) - dels_u*nn_z
!
!for alo
!                  int_u=0.d0
!                  int_u=xic1
                  int_u=calc_icore_gdark(z_u, sqrt(x_u**2+y_u**2+z_u**2))
!
                  call coeff3d_linecu_lin(iim1, i, jjm1, j, kkm1, k, x_u, y_u, z_u, &
                     c14_slineu, c15_slineu, c17_slineu, c18_slineu, &
                     c23_slineu, c24_slineu, c26_slineu, c27_slineu, &
                     opac_u, opalbar_u, velx_u, vely_u, velz_u, vth_u, scont_u, sline_u)
                  c02_slineu=0.d0
                  c04_slineu=0.d0
                  c05_slineu=0.d0
                  c06_slineu=0.d0
                  c08_slineu=0.d0
                  c10_slineu=0.d0
                  c11_slineu=0.d0
                  c12_slineu=0.d0
                  c13_slineu=0.d0
                  c16_slineu=0.d0
                  c20_slineu=0.d0
                  c22_slineu=0.d0
!set interpolation coefficients that are not used to zero
                  c02_intu = 0.d0
                  c04_intu = 0.d0
                  c05_intu = 0.d0
                  c06_intu = 0.d0
                  c08_intu = 0.d0
                  c10_intu = 0.d0
                  c11_intu = 0.d0
                  c12_intu = 0.d0
                  c13_intu = 0.d0
                  c14_intu = 0.d0
                  c15_intu = 0.d0
                  c16_intu = 0.d0
                  c17_intu = 0.d0
                  c18_intu = 0.d0
                  c20_intu = 0.d0
                  c22_intu = 0.d0
                  c23_intu = 0.d0
                  c24_intu = 0.d0
                  c26_intu = 0.d0
!
               elseif(dels_xyu.eq.dels_u) then
!intersection with x-y plane on level k-gamma
                  x_u = x(i) - dels_u*nn_x
                  y_u = y(j) - dels_u*nn_y
                  z_u = z(kkm1)
!
                  iim2=i-2*alpha
                  jjm2=j-2*beta
!
                  call coeff2d_linecu(opac3d(iim2,jjm2,kkm1), opac3d(iim1,jjm2,kkm1), opac3d(i,jjm2,kkm1), &
                     opac3d(iim2,jjm1,kkm1), opac3d(iim1,jjm1,kkm1), opac3d(i,jjm1,kkm1), &
                     opac3d(iim2,j,kkm1),    opac3d(iim1,j,kkm1),    opac3d(i,j,kkm1), &
                     opalbar3d(iim2,jjm2,kkm1), opalbar3d(iim1,jjm2,kkm1), opalbar3d(i,jjm2,kkm1), &
                     opalbar3d(iim2,jjm1,kkm1), opalbar3d(iim1,jjm1,kkm1), opalbar3d(i,jjm1,kkm1), &
                     opalbar3d(iim2,j,kkm1),    opalbar3d(iim1,j,kkm1),    opalbar3d(i,j,kkm1), &
                     velx3d(iim2,jjm2,kkm1), velx3d(iim1,jjm2,kkm1), velx3d(i,jjm2,kkm1), &
                     velx3d(iim2,jjm1,kkm1), velx3d(iim1,jjm1,kkm1), velx3d(i,jjm1,kkm1), &
                     velx3d(iim2,j,kkm1),    velx3d(iim1,j,kkm1),    velx3d(i,j,kkm1), &
                     vely3d(iim2,jjm2,kkm1), vely3d(iim1,jjm2,kkm1), vely3d(i,jjm2,kkm1), &
                     vely3d(iim2,jjm1,kkm1), vely3d(iim1,jjm1,kkm1), vely3d(i,jjm1,kkm1), &
                     vely3d(iim2,j,kkm1),    vely3d(iim1,j,kkm1),    vely3d(i,j,kkm1), &
                     velz3d(iim2,jjm2,kkm1), velz3d(iim1,jjm2,kkm1), velz3d(i,jjm2,kkm1), &
                     velz3d(iim2,jjm1,kkm1), velz3d(iim1,jjm1,kkm1), velz3d(i,jjm1,kkm1), &
                     velz3d(iim2,j,kkm1),    velz3d(iim1,j,kkm1),    velz3d(i,j,kkm1), &
                     vth3d(iim2,jjm2,kkm1), vth3d(iim1,jjm2,kkm1), vth3d(i,jjm2,kkm1), &
                     vth3d(iim2,jjm1,kkm1), vth3d(iim1,jjm1,kkm1), vth3d(i,jjm1,kkm1), &
                     vth3d(iim2,j,kkm1),    vth3d(iim1,j,kkm1),    vth3d(i,j,kkm1), &
                     scont3d(iim2,jjm2,kkm1), scont3d(iim1,jjm2,kkm1), scont3d(i,jjm2,kkm1), &
                     scont3d(iim2,jjm1,kkm1), scont3d(iim1,jjm1,kkm1), scont3d(i,jjm1,kkm1), &
                     scont3d(iim2,j,kkm1),    scont3d(iim1,j,kkm1),    scont3d(i,j,kkm1), &
                     sline3d(iim2,jjm2,kkm1), sline3d(iim1,jjm2,kkm1), sline3d(i,jjm2,kkm1), &
                     sline3d(iim2,jjm1,kkm1), sline3d(iim1,jjm1,kkm1), sline3d(i,jjm1,kkm1), &
                     sline3d(iim2,j,kkm1),    sline3d(iim1,j,kkm1),    sline3d(i,j,kkm1), &
                     int3d(iim2,jjm2,kkm1), int3d(iim1,jjm2,kkm1), int3d(i,jjm2,kkm1), &
                     int3d(iim2,jjm1,kkm1), int3d(iim1,jjm1,kkm1), int3d(i,jjm1,kkm1), &
                     int3d(iim2,j,kkm1),    int3d(iim1,j,kkm1),    int3d(i,j,kkm1), &
                     x(iim2), x(iim1), x(i), y(jjm2), y(jjm1), y(j), x_u, y_u, &
                     c10_slineu, c11_slineu, c12_slineu, c13_slineu, c14_slineu, &
                     c15_slineu, c16_slineu, c17_slineu, c18_slineu, &
                     c10_intu, c11_intu, c12_intu, c13_intu, c14_intu, &
                     c15_intu, c16_intu, c17_intu, c18_intu, &
                     opac_u, opalbar_u, velx_u, vely_u, velz_u, vth_u, scont_u, sline_u, int_u)
!set interpolation coefficients that are not used to zero
                  c23_slineu=0.d0
                  c24_slineu=0.d0
                  c26_slineu=0.d0
                  c27_slineu=0.d0
!
                  c02_intu = 0.d0
                  c04_intu = 0.d0
                  c05_intu = 0.d0
                  c06_intu = 0.d0
                  c08_intu = 0.d0
                  c20_intu = 0.d0
                  c22_intu = 0.d0
                  c23_intu = 0.d0
                  c24_intu = 0.d0
                  c26_intu = 0.d0
!
               elseif(dels_xzu.eq.dels_u) then
!intersection with x-z plane on level j-beta
                  x_u = x(i) - dels_u*nn_x
                  y_u = y(jjm1)
                  z_u = z(k) - dels_u*nn_z
!
                  iim2=i-2*alpha
                  kkm2=k-2*gamma
!
                  call coeff2d_linecu(opac3d(iim2,jjm1,kkm2), opac3d(iim1,jjm1,kkm2), opac3d(i,jjm1,kkm2), &
                     opac3d(iim2,jjm1,kkm1), opac3d(iim1,jjm1,kkm1), opac3d(i,jjm1,kkm1), &
                     opac3d(iim2,jjm1,k), opac3d(iim1,jjm1,k), opac3d(i,jjm1,k), &
                     opalbar3d(iim2,jjm1,kkm2), opalbar3d(iim1,jjm1,kkm2), opalbar3d(i,jjm1,kkm2), &
                     opalbar3d(iim2,jjm1,kkm1), opalbar3d(iim1,jjm1,kkm1), opalbar3d(i,jjm1,kkm1), &
                     opalbar3d(iim2,jjm1,k), opalbar3d(iim1,jjm1,k), opalbar3d(i,jjm1,k), &
                     velx3d(iim2,jjm1,kkm2), velx3d(iim1,jjm1,kkm2), velx3d(i,jjm1,kkm2), &
                     velx3d(iim2,jjm1,kkm1), velx3d(iim1,jjm1,kkm1), velx3d(i,jjm1,kkm1), &
                     velx3d(iim2,jjm1,k), velx3d(iim1,jjm1,k), velx3d(i,jjm1,k), &
                     vely3d(iim2,jjm1,kkm2), vely3d(iim1,jjm1,kkm2), vely3d(i,jjm1,kkm2), &
                     vely3d(iim2,jjm1,kkm1), vely3d(iim1,jjm1,kkm1), vely3d(i,jjm1,kkm1), &
                     vely3d(iim2,jjm1,k), vely3d(iim1,jjm1,k), vely3d(i,jjm1,k), &
                     velz3d(iim2,jjm1,kkm2), velz3d(iim1,jjm1,kkm2), velz3d(i,jjm1,kkm2), &
                     velz3d(iim2,jjm1,kkm1), velz3d(iim1,jjm1,kkm1), velz3d(i,jjm1,kkm1), &
                     velz3d(iim2,jjm1,k), velz3d(iim1,jjm1,k), velz3d(i,jjm1,k), &
                     vth3d(iim2,jjm1,kkm2), vth3d(iim1,jjm1,kkm2), vth3d(i,jjm1,kkm2), &
                     vth3d(iim2,jjm1,kkm1), vth3d(iim1,jjm1,kkm1), vth3d(i,jjm1,kkm1), &
                     vth3d(iim2,jjm1,k), vth3d(iim1,jjm1,k), vth3d(i,jjm1,k), &
                     scont3d(iim2,jjm1,kkm2), scont3d(iim1,jjm1,kkm2), scont3d(i,jjm1,kkm2), &
                     scont3d(iim2,jjm1,kkm1), scont3d(iim1,jjm1,kkm1), scont3d(i,jjm1,kkm1), &
                     scont3d(iim2,jjm1,k), scont3d(iim1,jjm1,k), scont3d(i,jjm1,k), &
                     sline3d(iim2,jjm1,kkm2), sline3d(iim1,jjm1,kkm2), sline3d(i,jjm1,kkm2), &
                     sline3d(iim2,jjm1,kkm1), sline3d(iim1,jjm1,kkm1), sline3d(i,jjm1,kkm1), &
                     sline3d(iim2,jjm1,k), sline3d(iim1,jjm1,k), sline3d(i,jjm1,k), &
                     int3d(iim2,jjm1,kkm2), int3d(iim1,jjm1,kkm2), int3d(i,jjm1,kkm2), &
                     int3d(iim2,jjm1,kkm1), int3d(iim1,jjm1,kkm1), int3d(i,jjm1,kkm1), &
                     int3d(iim2,jjm1,k), int3d(iim1,jjm1,k), int3d(i,jjm1,k), &
                     x(iim2), x(iim1), x(i), z(kkm2), z(kkm1), z(k), x_u, z_u, &
                     c04_slineu, c05_slineu, c06_slineu, c13_slineu, c14_slineu, &
                     c15_slineu, c22_slineu, c23_slineu, c24_slineu, &
                     c04_intu, c05_intu, c06_intu, c13_intu, c14_intu, &
                     c15_intu, c22_intu, c23_intu, c24_intu, &
                     opac_u, opalbar_u, velx_u, vely_u, velz_u, vth_u, scont_u, sline_u, int_u)
!set interpolation coefficients that are not used to zero
                  c17_slineu=0.d0
                  c18_slineu=0.d0
                  c26_slineu=0.d0
                  c27_slineu=0.d0
                  c02_intu = 0.d0
                  c08_intu = 0.d0
                  c10_intu = 0.d0
                  c11_intu = 0.d0
                  c12_intu = 0.d0
                  c16_intu = 0.d0
                  c17_intu = 0.d0
                  c18_intu = 0.d0
                  c20_intu = 0.d0
                  c26_intu = 0.d0
!
               elseif(dels_yzu.eq.dels_u) then
!intersection with y-z plane on level i-alpha
                  x_u = x(iim1)
                  y_u = y(j) - dels_u*nn_y
                  z_u = z(k) - dels_u*nn_z
!
                  jjm2=j-2*beta
                  kkm2=k-2*gamma
!
                  call coeff2d_linecu(opac3d(iim1,jjm2,kkm2), opac3d(iim1,jjm1,kkm2), opac3d(iim1,j,kkm2), &
                     opac3d(iim1,jjm2,kkm1), opac3d(iim1,jjm1,kkm1), opac3d(iim1,j,kkm1), &
                     opac3d(iim1,jjm2,k), opac3d(iim1,jjm1,k), opac3d(iim1,j,k), &
                     opalbar3d(iim1,jjm2,kkm2), opalbar3d(iim1,jjm1,kkm2), opalbar3d(iim1,j,kkm2), &
                     opalbar3d(iim1,jjm2,kkm1), opalbar3d(iim1,jjm1,kkm1), opalbar3d(iim1,j,kkm1), &
                     opalbar3d(iim1,jjm2,k), opalbar3d(iim1,jjm1,k), opalbar3d(iim1,j,k), &
                     velx3d(iim1,jjm2,kkm2), velx3d(iim1,jjm1,kkm2), velx3d(iim1,j,kkm2), &
                     velx3d(iim1,jjm2,kkm1), velx3d(iim1,jjm1,kkm1), velx3d(iim1,j,kkm1), &
                     velx3d(iim1,jjm2,k), velx3d(iim1,jjm1,k), velx3d(iim1,j,k), &
                     vely3d(iim1,jjm2,kkm2), vely3d(iim1,jjm1,kkm2), vely3d(iim1,j,kkm2), &
                     vely3d(iim1,jjm2,kkm1), vely3d(iim1,jjm1,kkm1), vely3d(iim1,j,kkm1), &
                     vely3d(iim1,jjm2,k), vely3d(iim1,jjm1,k), vely3d(iim1,j,k), &
                     velz3d(iim1,jjm2,kkm2), velz3d(iim1,jjm1,kkm2), velz3d(iim1,j,kkm2), &
                     velz3d(iim1,jjm2,kkm1), velz3d(iim1,jjm1,kkm1), velz3d(iim1,j,kkm1), &
                     velz3d(iim1,jjm2,k), velz3d(iim1,jjm1,k), velz3d(iim1,j,k), &
                     vth3d(iim1,jjm2,kkm2), vth3d(iim1,jjm1,kkm2), vth3d(iim1,j,kkm2), &
                     vth3d(iim1,jjm2,kkm1), vth3d(iim1,jjm1,kkm1), vth3d(iim1,j,kkm1), &
                     vth3d(iim1,jjm2,k), vth3d(iim1,jjm1,k), vth3d(iim1,j,k), &
                     scont3d(iim1,jjm2,kkm2), scont3d(iim1,jjm1,kkm2), scont3d(iim1,j,kkm2), &
                     scont3d(iim1,jjm2,kkm1), scont3d(iim1,jjm1,kkm1), scont3d(iim1,j,kkm1), &
                     scont3d(iim1,jjm2,k), scont3d(iim1,jjm1,k), scont3d(iim1,j,k), &
                     sline3d(iim1,jjm2,kkm2), sline3d(iim1,jjm1,kkm2), sline3d(iim1,j,kkm2), &
                     sline3d(iim1,jjm2,kkm1), sline3d(iim1,jjm1,kkm1), sline3d(iim1,j,kkm1), &
                     sline3d(iim1,jjm2,k), sline3d(iim1,jjm1,k), sline3d(iim1,j,k), &
                     int3d(iim1,jjm2,kkm2), int3d(iim1,jjm1,kkm2), int3d(iim1,j,kkm2), &
                     int3d(iim1,jjm2,kkm1), int3d(iim1,jjm1,kkm1), int3d(iim1,j,kkm1), &
                     int3d(iim1,jjm2,k), int3d(iim1,jjm1,k), int3d(iim1,j,k), &
                     y(jjm2), y(jjm1), y(j), z(kkm2), z(kkm1), z(k), y_u, z_u, &
                     c02_slineu, c05_slineu, c08_slineu, c11_slineu, c14_slineu, &
                     c17_slineu, c20_slineu, c23_slineu, c26_slineu, &
                     c02_intu, c05_intu, c08_intu, c11_intu, c14_intu, &
                     c17_intu, c20_intu, c23_intu, c26_intu, &
                     opac_u, opalbar_u, velx_u, vely_u, velz_u, vth_u, scont_u, sline_u, int_u)
!set interpolation coefficients that are not used to zero
                  c15_slineu=0.d0
                  c18_slineu=0.d0
                  c24_slineu=0.d0
                  c27_slineu=0.d0
!
                  c04_intu = 0.d0
                  c06_intu = 0.d0
                  c10_intu = 0.d0
                  c12_intu = 0.d0
                  c13_intu = 0.d0
                  c15_intu = 0.d0
                  c16_intu = 0.d0
                  c18_intu = 0.d0
                  c22_intu = 0.d0
                  c24_intu = 0.d0
!
               else
                  write(*,'(4es20.8,2l4)') dels_u, dels_xzu, dels_xyu, dels_yzu, dels_u.eq.dels_xzu, dels_xyu.eq.dels_u
                  stop 'error in fsc_linec3d: invalid dels_u'
               endif
!
!---------------------------downwind point------------------------------
!
               if(dels_xyd.eq.dels_d) then
!intersection with x-y plane on level k+gamma
                  x_d = x(i) + dels_d*nn_x
                  y_d = y(j) + dels_d*nn_y
                  z_d = z(kkp1)
!
                  call coeff2d_linecd(opac3d(iim1,jjm1,kkp1), opac3d(i,jjm1,kkp1), opac3d(iip1,jjm1,kkp1), &
                     opac3d(iim1,j,kkp1), opac3d(i,j,kkp1), opac3d(iip1,j,kkp1), &
                     opac3d(iim1,jjp1,kkp1), opac3d(i,jjp1,kkp1), opac3d(iip1,jjp1,kkp1), &
                     opalbar3d(iim1,jjm1,kkp1), opalbar3d(i,jjm1,kkp1), opalbar3d(iip1,jjm1,kkp1), &
                     opalbar3d(iim1,j,kkp1), opalbar3d(i,j,kkp1), opalbar3d(iip1,j,kkp1), &
                     opalbar3d(iim1,jjp1,kkp1), opalbar3d(i,jjp1,kkp1), opalbar3d(iip1,jjp1,kkp1), &
                     velx3d(iim1,jjm1,kkp1), velx3d(i,jjm1,kkp1), velx3d(iip1,jjm1,kkp1), &
                     velx3d(iim1,j,kkp1), velx3d(i,j,kkp1), velx3d(iip1,j,kkp1), &
                     velx3d(iim1,jjp1,kkp1), velx3d(i,jjp1,kkp1), velx3d(iip1,jjp1,kkp1), &
                     vely3d(iim1,jjm1,kkp1), vely3d(i,jjm1,kkp1), vely3d(iip1,jjm1,kkp1), &
                     vely3d(iim1,j,kkp1), vely3d(i,j,kkp1), vely3d(iip1,j,kkp1), &
                     vely3d(iim1,jjp1,kkp1), vely3d(i,jjp1,kkp1), vely3d(iip1,jjp1,kkp1), &
                     velz3d(iim1,jjm1,kkp1), velz3d(i,jjm1,kkp1), velz3d(iip1,jjm1,kkp1), &
                     velz3d(iim1,j,kkp1), velz3d(i,j,kkp1), velz3d(iip1,j,kkp1), &
                     velz3d(iim1,jjp1,kkp1), velz3d(i,jjp1,kkp1), velz3d(iip1,jjp1,kkp1), &
                     vth3d(iim1,jjm1,kkp1), vth3d(i,jjm1,kkp1), vth3d(iip1,jjm1,kkp1), &
                     vth3d(iim1,j,kkp1), vth3d(i,j,kkp1), vth3d(iip1,j,kkp1), &
                     vth3d(iim1,jjp1,kkp1), vth3d(i,jjp1,kkp1), vth3d(iip1,jjp1,kkp1), &
                     scont3d(iim1,jjm1,kkp1), scont3d(i,jjm1,kkp1), scont3d(iip1,jjm1,kkp1), &
                     scont3d(iim1,j,kkp1), scont3d(i,j,kkp1), scont3d(iip1,j,kkp1), &
                     scont3d(iim1,jjp1,kkp1), scont3d(i,jjp1,kkp1), scont3d(iip1,jjp1,kkp1), &
                     sline3d(iim1,jjm1,kkp1), sline3d(i,jjm1,kkp1), sline3d(iip1,jjm1,kkp1), &
                     sline3d(iim1,j,kkp1), sline3d(i,j,kkp1), sline3d(iip1,j,kkp1), &
                     sline3d(iim1,jjp1,kkp1), sline3d(i,jjp1,kkp1), sline3d(iip1,jjp1,kkp1), &
                     x(iim1), x(i), x(iip1), y(jjm1), y(j), y(jjp1), x_d, y_d, &
                     c19_slined, c20_slined, c21_slined, c22_slined, c23_slined, c24_slined, &
                     c25_slined, c26_slined, c27_slined, &
                     opac_d, opalbar_d, velx_d, vely_d, velz_d, vth_d, scont_d, sline_d)
!
!set interpolation coefficients that are not used to zero
                  c03_slined = 0.d0
                  c06_slined = 0.d0
                  c07_slined = 0.d0
                  c08_slined = 0.d0
                  c09_slined = 0.d0
                  c12_slined = 0.d0
                  c15_slined = 0.d0
                  c16_slined = 0.d0
                  c17_slined = 0.d0
                  c18_slined = 0.d0
!
               elseif(dels_xzd.eq.dels_d) then
!intersection with x-z plane on level j+beta
                  x_d = x(i) + dels_d*nn_x
                  y_d = y(jjp1)
                  z_d = z(k) + dels_d*nn_z
!
                  call coeff2d_linecd(opac3d(iim1,jjp1,kkm1), opac3d(i,jjp1,kkm1), opac3d(iip1,jjp1,kkm1), &
                     opac3d(iim1,jjp1,k), opac3d(i,jjp1,k), opac3d(iip1,jjp1,k), &
                     opac3d(iim1,jjp1,kkp1), opac3d(i,jjp1,kkp1), opac3d(iip1,jjp1,kkp1), &
                     opalbar3d(iim1,jjp1,kkm1), opalbar3d(i,jjp1,kkm1), opalbar3d(iip1,jjp1,kkm1), &
                     opalbar3d(iim1,jjp1,k), opalbar3d(i,jjp1,k), opalbar3d(iip1,jjp1,k), &
                     opalbar3d(iim1,jjp1,kkp1), opalbar3d(i,jjp1,kkp1), opalbar3d(iip1,jjp1,kkp1), &
                     velx3d(iim1,jjp1,kkm1), velx3d(i,jjp1,kkm1), velx3d(iip1,jjp1,kkm1), &
                     velx3d(iim1,jjp1,k), velx3d(i,jjp1,k), velx3d(iip1,jjp1,k), &
                     velx3d(iim1,jjp1,kkp1), velx3d(i,jjp1,kkp1), velx3d(iip1,jjp1,kkp1), &
                     vely3d(iim1,jjp1,kkm1), vely3d(i,jjp1,kkm1), vely3d(iip1,jjp1,kkm1), &
                     vely3d(iim1,jjp1,k), vely3d(i,jjp1,k), vely3d(iip1,jjp1,k), &
                     vely3d(iim1,jjp1,kkp1), vely3d(i,jjp1,kkp1), vely3d(iip1,jjp1,kkp1), &
                     velz3d(iim1,jjp1,kkm1), velz3d(i,jjp1,kkm1), velz3d(iip1,jjp1,kkm1), &
                     velz3d(iim1,jjp1,k), velz3d(i,jjp1,k), velz3d(iip1,jjp1,k), &
                     velz3d(iim1,jjp1,kkp1), velz3d(i,jjp1,kkp1), velz3d(iip1,jjp1,kkp1), &
                     vth3d(iim1,jjp1,kkm1), vth3d(i,jjp1,kkm1), vth3d(iip1,jjp1,kkm1), &
                     vth3d(iim1,jjp1,k), vth3d(i,jjp1,k), vth3d(iip1,jjp1,k), &
                     vth3d(iim1,jjp1,kkp1), vth3d(i,jjp1,kkp1), vth3d(iip1,jjp1,kkp1), &
                     scont3d(iim1,jjp1,kkm1), scont3d(i,jjp1,kkm1), scont3d(iip1,jjp1,kkm1), &
                     scont3d(iim1,jjp1,k), scont3d(i,jjp1,k), scont3d(iip1,jjp1,k), &
                     scont3d(iim1,jjp1,kkp1), scont3d(i,jjp1,kkp1), scont3d(iip1,jjp1,kkp1), &
                     sline3d(iim1,jjp1,kkm1), sline3d(i,jjp1,kkm1), sline3d(iip1,jjp1,kkm1), &
                     sline3d(iim1,jjp1,k), sline3d(i,jjp1,k), sline3d(iip1,jjp1,k), &
                     sline3d(iim1,jjp1,kkp1), sline3d(i,jjp1,kkp1), sline3d(iip1,jjp1,kkp1), &
                     x(iim1), x(i), x(iip1), z(kkm1), z(k), z(kkp1), x_d, z_d, &
                     c07_slined, c08_slined, c09_slined, c16_slined, c17_slined, c18_slined, &
                     c25_slined, c26_slined, c27_slined, &
                     opac_d, opalbar_d, velx_d, vely_d, velz_d, vth_d, scont_d, sline_d)
!
!set interpolation coefficients that are not used to zero
                  c03_slined = 0.d0
                  c06_slined = 0.d0
                  c12_slined = 0.d0
                  c15_slined = 0.d0
                  c19_slined = 0.d0
                  c20_slined = 0.d0
                  c21_slined = 0.d0
                  c22_slined = 0.d0
                  c23_slined = 0.d0
                  c24_slined = 0.d0
!
               elseif(dels_yzd.eq.dels_d) then
!intersection with y-z plane on level i+alpha
                  x_d = x(iip1)
                  y_d = y(j) + dels_d*nn_y
                  z_d = z(k) + dels_d*nn_z
!
                  call coeff2d_linecd(opac3d(iip1,jjm1,kkm1), opac3d(iip1,j,kkm1), opac3d(iip1,jjp1,kkm1), &
                     opac3d(iip1,jjm1,k), opac3d(iip1,j,k), opac3d(iip1,jjp1,k), &
                     opac3d(iip1,jjm1,kkp1), opac3d(iip1,j,kkp1), opac3d(iip1,jjp1,kkp1), &
                     opalbar3d(iip1,jjm1,kkm1), opalbar3d(iip1,j,kkm1), opalbar3d(iip1,jjp1,kkm1), &
                     opalbar3d(iip1,jjm1,k), opalbar3d(iip1,j,k), opalbar3d(iip1,jjp1,k), &
                     opalbar3d(iip1,jjm1,kkp1), opalbar3d(iip1,j,kkp1), opalbar3d(iip1,jjp1,kkp1), &
                     velx3d(iip1,jjm1,kkm1), velx3d(iip1,j,kkm1), velx3d(iip1,jjp1,kkm1), &
                     velx3d(iip1,jjm1,k), velx3d(iip1,j,k), velx3d(iip1,jjp1,k), &
                     velx3d(iip1,jjm1,kkp1), velx3d(iip1,j,kkp1), velx3d(iip1,jjp1,kkp1), &
                     vely3d(iip1,jjm1,kkm1), vely3d(iip1,j,kkm1), vely3d(iip1,jjp1,kkm1), &
                     vely3d(iip1,jjm1,k), vely3d(iip1,j,k), vely3d(iip1,jjp1,k), &
                     vely3d(iip1,jjm1,kkp1), vely3d(iip1,j,kkp1), vely3d(iip1,jjp1,kkp1), &
                     velz3d(iip1,jjm1,kkm1), velz3d(iip1,j,kkm1), velz3d(iip1,jjp1,kkm1), &
                     velz3d(iip1,jjm1,k), velz3d(iip1,j,k), velz3d(iip1,jjp1,k), &
                     velz3d(iip1,jjm1,kkp1), velz3d(iip1,j,kkp1), velz3d(iip1,jjp1,kkp1), &
                     vth3d(iip1,jjm1,kkm1), vth3d(iip1,j,kkm1), vth3d(iip1,jjp1,kkm1), &
                     vth3d(iip1,jjm1,k), vth3d(iip1,j,k), vth3d(iip1,jjp1,k), &
                     vth3d(iip1,jjm1,kkp1), vth3d(iip1,j,kkp1), vth3d(iip1,jjp1,kkp1), &
                     scont3d(iip1,jjm1,kkm1), scont3d(iip1,j,kkm1), scont3d(iip1,jjp1,kkm1), &
                     scont3d(iip1,jjm1,k), scont3d(iip1,j,k), scont3d(iip1,jjp1,k), &
                     scont3d(iip1,jjm1,kkp1), scont3d(iip1,j,kkp1), scont3d(iip1,jjp1,kkp1), &
                     sline3d(iip1,jjm1,kkm1), sline3d(iip1,j,kkm1), sline3d(iip1,jjp1,kkm1), &
                     sline3d(iip1,jjm1,k), sline3d(iip1,j,k), sline3d(iip1,jjp1,k), &
                     sline3d(iip1,jjm1,kkp1), sline3d(iip1,j,kkp1), sline3d(iip1,jjp1,kkp1), &
                     y(jjm1), y(j), y(jjp1), z(kkm1), z(k), z(kkp1), y_d, z_d, &
                     c03_slined, c06_slined, c09_slined, c12_slined, c15_slined, c18_slined, &
                     c21_slined, c24_slined, c27_slined, &
                     opac_d, opalbar_d, velx_d, vely_d, velz_d, vth_d, scont_d, sline_d)
!
!set interpolation coefficients that are not used to zero
                  c07_slined = 0.d0
                  c08_slined = 0.d0
                  c16_slined = 0.d0
                  c17_slined = 0.d0
                  c19_slined = 0.d0
                  c20_slined = 0.d0
                  c22_slined = 0.d0
                  c23_slined = 0.d0
                  c25_slined = 0.d0
                  c26_slined = 0.d0
!
               else
                  write(*,'(4es20.8,l4)') dels_d, dels_xzd, dels_xyd, dels_yzd, dels_d.eq.dels_xzd
                  stop 'error in fsc_linec3d: invalid dels_d'
               endif
!
!--------------------------------radiative transfer---------------------
!
               vel_u = velx_u*nn_x + vely_u*nn_y + velz_u*nn_z
               vel_d = velx_d*nn_x + vely_d*nn_y + velz_d*nn_z
               call fsc_linec(xobs, vth_fiducial, deltax, xcmf_max, xcmf_min, int_u, opac_u, opac_p, opac_d, &
                  opalbar_u, opalbar_p, opalbar_d, scont_u, scont_p, scont_d, &
                  sline_u, sline_p, sline_d, vel_u, vel_p, vel_d, vth_u, vth_p, vth_d, &
                  dels_u, dels_d, abs_sc, contr_sc, int_sc, alo_u, alo_p, alo_d)
               int3d(i,j,k) = int_sc
!
               aloline_on_nn3d(iip1,jjp1,kkp1,q1) = imask_totreg3d(iip1,jjp1,kkp1)*&
                  alo_d*c27_slined
               aloline_on_nn3d(i,jjp1,kkp1,q2) = imask_totreg3d(i,jjp1,kkp1) * &
                  (alo_d*c26_slined + abs_sc*c26_intu*aloline_on_nn3d(i,jjp1,kkp1,q1))
               aloline_on_nn3d(iim1,jjp1,kkp1,q3) = imask_totreg3d(iim1,jjp1,kkp1) * &
                  (alo_d*c25_slined + abs_sc*c26_intu*aloline_on_nn3d(iim1,jjp1,kkp1,q2))
               aloline_on_nn3d(iip1,j,kkp1,q4) = imask_totreg3d(iip1,j,kkp1) * &
                  (alo_d*c24_slined + abs_sc*c24_intu*aloline_on_nn3d(iip1,j,kkp1,q1))
               aloline_on_nn3d(i,j,kkp1,q5) = imask_totreg3d(i,j,kkp1) * &
                  (alo_d*c23_slined + abs_sc*(c23_intu*aloline_on_nn3d(i,j,kkp1,q1) + &
                  c24_intu*aloline_on_nn3d(i,j,kkp1,q2) + &
                  c26_intu*aloline_on_nn3d(i,j,kkp1,q4)))
               aloline_on_nn3d(iim1,j,kkp1,q6) = imask_totreg3d(iim1,j,kkp1) * &
                  (alo_d*c22_slined + abs_sc*(c22_intu*aloline_on_nn3d(iim1,j,kkp1,q1) + &
                  c23_intu*aloline_on_nn3d(iim1,j,kkp1,q2) + &
                  c24_intu*aloline_on_nn3d(iim1,j,kkp1,q3) + &
                  c26_intu*aloline_on_nn3d(iim1,j,kkp1,q5)))
               aloline_on_nn3d(iip1,jjm1,kkp1,q7) = imask_totreg3d(iip1,jjm1,kkp1) * &
                  (alo_d*c21_slined + abs_sc*c24_intu*aloline_on_nn3d(iip1,jjm1,kkp1,q4))
               aloline_on_nn3d(i,jjm1,kkp1,q8) = imask_totreg3d(i,jjm1,kkp1) * &
                  (alo_d*c20_slined + abs_sc*(c23_intu*aloline_on_nn3d(i,jjm1,kkp1,q4) + &
                  c24_intu*aloline_on_nn3d(i,jjm1,kkp1,q5) + &
                  c20_intu*aloline_on_nn3d(i,jjm1,kkp1,q1) + &
                  c26_intu*aloline_on_nn3d(i,jjm1,kkp1,q7)))
               aloline_on_nn3d(iim1,jjm1,kkp1,q9) = imask_totreg3d(iim1,jjm1,kkp1) * &
                  (alo_d*c19_slined + abs_sc*(c22_intu*aloline_on_nn3d(iim1,jjm1,kkp1,q4) + &
                  c23_intu*aloline_on_nn3d(iim1,jjm1,kkp1,q5) + &
                  c24_intu*aloline_on_nn3d(iim1,jjm1,kkp1,q6) + &
                  c20_intu*aloline_on_nn3d(iim1,jjm1,kkp1,q2) + &
                  c26_intu*aloline_on_nn3d(iim1,jjm1,kkp1,q8)))
               aloline_on_nn3d(iip1,jjp1,k,q10) = imask_totreg3d(iip1,jjp1,k) * &
                  (alo_d*c18_slined + abs_sc*c18_intu*aloline_on_nn3d(iip1,jjp1,k,q1))
               aloline_on_nn3d(i,jjp1,k,q11) = imask_totreg3d(i,jjp1,k) * &
                  (alo_d*c17_slined + abs_sc*(c17_intu*aloline_on_nn3d(i,jjp1,k,q1) + &
                  c18_intu*aloline_on_nn3d(i,jjp1,k,q2) + &
                  c26_intu*aloline_on_nn3d(i,jjp1,k,q10)))
               aloline_on_nn3d(iim1,jjp1,k,q12) = imask_totreg3d(iim1,jjp1,k) * &
                  (alo_d*c16_slined + abs_sc*(c16_intu*aloline_on_nn3d(iim1,jjp1,k,q1) + &
                  c17_intu*aloline_on_nn3d(iim1,jjp1,k,q2) + &
                  c18_intu*aloline_on_nn3d(iim1,jjp1,k,q3) + &
                  c26_intu*aloline_on_nn3d(iim1,jjp1,k,q11)))
               aloline_on_nn3d(iip1,j,k,q13) = imask_totreg3d(iip1,j,k) * &
                  (alo_d*c15_slined + abs_sc*(c15_intu*aloline_on_nn3d(iip1,j,k,q1) + &
                  c18_intu*aloline_on_nn3d(iip1,j,k,q4) + &
                  c24_intu*aloline_on_nn3d(iip1,j,k,q10)))
               aloline_on_nn3d(i,j,k,q14) = imask_totreg3d(i,j,k) * &
                  (alo_p + abs_sc*(c14_intu*aloline_on_nn3d(i,j,k,q1) + &
                  c15_intu*aloline_on_nn3d(i,j,k,q2) + &
                  c17_intu*aloline_on_nn3d(i,j,k,q4) + &
                  c26_intu*aloline_on_nn3d(i,j,k,q13) + &
                  c18_intu*aloline_on_nn3d(i,j,k,q5) + &
                  c23_intu*aloline_on_nn3d(i,j,k,q10) + &
                  c24_intu*aloline_on_nn3d(i,j,k,q11)))
               aloline_on_nn3d(iim1,j,k,q15) = imask_totreg3d(iim1,j,k) * &
                  (alo_u*c26_slineu + abs_sc*(c13_intu*aloline_on_nn3d(iim1,j,k,q1) + &
                  c14_intu*aloline_on_nn3d(iim1,j,k,q2) + &
                  c15_intu*aloline_on_nn3d(iim1,j,k,q3) + &
                  c16_intu*aloline_on_nn3d(iim1,j,k,q4) + &
                  c17_intu*aloline_on_nn3d(iim1,j,k,q5) + &
                  c18_intu*aloline_on_nn3d(iim1,j,k,q6) + &
                  c22_intu*aloline_on_nn3d(iim1,j,k,q10) + &
                  c23_intu*aloline_on_nn3d(iim1,j,k,q11) + &
                  c24_intu*aloline_on_nn3d(iim1,j,k,q12) + &
                  c26_intu*aloline_on_nn3d(iim1,j,k,q14)))
               aloline_on_nn3d(iip1,jjm1,k,q16) = imask_totreg3d(iip1,jjm1,k) * &
                  (alo_d*c12_slined + abs_sc*(c12_intu*aloline_on_nn3d(iip1,jjm1,k,q1) + &
                  c15_intu*aloline_on_nn3d(iip1,jjm1,k,q4) + &
                  c18_intu*aloline_on_nn3d(iip1,jjm1,k,q7) + &
                  c24_intu*aloline_on_nn3d(iip1,jjm1,k,q13)))
               aloline_on_nn3d(i,jjm1,k,q17) = imask_totreg3d(i,jjm1,k) * &
                  (alo_u*c24_slineu + abs_sc*(c11_intu*aloline_on_nn3d(i,jjm1,k,q1) + &
                  c12_intu*aloline_on_nn3d(i,jjm1,k,q2) + &
                  c14_intu*aloline_on_nn3d(i,jjm1,k,q4) + &
                  c15_intu*aloline_on_nn3d(i,jjm1,k,q5) + &
                  c17_intu*aloline_on_nn3d(i,jjm1,k,q7) + &
                  c18_intu*aloline_on_nn3d(i,jjm1,k,q8) + &
                  c23_intu*aloline_on_nn3d(i,jjm1,k,q13) + &
                  c24_intu*aloline_on_nn3d(i,jjm1,k,q14) + &
                  c20_intu*aloline_on_nn3d(i,jjm1,k,q10) + &
                  c26_intu*aloline_on_nn3d(i,jjm1,k,q16)))
               aloline_on_nn3d(iim1,jjm1,k,q18) = imask_totreg3d(iim1,jjm1,k) * &
                  (alo_u*c23_slineu + abs_sc*(c10_intu*aloline_on_nn3d(iim1,jjm1,k,q1) + &
                  c11_intu*aloline_on_nn3d(iim1,jjm1,k,q2) + &
                  c12_intu*aloline_on_nn3d(iim1,jjm1,k,q3) + &
                  c13_intu*aloline_on_nn3d(iim1,jjm1,k,q4) + &
                  c14_intu*aloline_on_nn3d(iim1,jjm1,k,q5) + &
                  c15_intu*aloline_on_nn3d(iim1,jjm1,k,q6) + &
                  c16_intu*aloline_on_nn3d(iim1,jjm1,k,q7) + &
                  c17_intu*aloline_on_nn3d(iim1,jjm1,k,q8) + &
                  c18_intu*aloline_on_nn3d(iim1,jjm1,k,q9) + &
                  c22_intu*aloline_on_nn3d(iim1,jjm1,k,q13) + &
                  c23_intu*aloline_on_nn3d(iim1,jjm1,k,q14) + &
                  c24_intu*aloline_on_nn3d(iim1,jjm1,k,q15) + &
                  c20_intu*aloline_on_nn3d(iim1,jjm1,k,q11) + &
                  c26_intu*aloline_on_nn3d(iim1,jjm1,k,q17)))
               aloline_on_nn3d(iip1,jjp1,kkm1,q19) = imask_totreg3d(iip1,jjp1,kkm1) * &
                  (alo_d*c09_slined + abs_sc*c18_intu*aloline_on_nn3d(iip1,jjp1,kkm1,q10))
               aloline_on_nn3d(i,jjp1,kkm1,q20) = imask_totreg3d(i,jjp1,kkm1) * &
                  (alo_d*c08_slined + abs_sc*(c17_intu*aloline_on_nn3d(i,jjp1,kkm1,q10) + &
                  c18_intu*aloline_on_nn3d(i,jjp1,kkm1,q11) + &
                  c08_intu*aloline_on_nn3d(i,jjp1,kkm1,q1) + &
                  c26_intu*aloline_on_nn3d(i,jjp1,kkm1,q19)))
               aloline_on_nn3d(iim1,jjp1,kkm1,q21) = imask_totreg3d(iim1,jjp1,kkm1) * &
                  (alo_d*c07_slined + abs_sc*(c16_intu*aloline_on_nn3d(iim1,jjp1,kkm1,q10) + &
                  c17_intu*aloline_on_nn3d(iim1,jjp1,kkm1,q11) + &
                  c18_intu*aloline_on_nn3d(iim1,jjp1,kkm1,q12) + &
                  c08_intu*aloline_on_nn3d(iim1,jjp1,kkm1,q2) + &
                  c26_intu*aloline_on_nn3d(iim1,jjp1,kkm1,q20)))
               aloline_on_nn3d(iip1,j,kkm1,q22) = imask_totreg3d(iip1,j,kkm1) * &
                  (alo_d*c06_slined + abs_sc*(c15_intu*aloline_on_nn3d(iip1,j,kkm1,q10) + &
                  c18_intu*aloline_on_nn3d(iip1,j,kkm1,q13) + &
                  c06_intu*aloline_on_nn3d(iip1,j,kkm1,q1) + &
                  c24_intu*aloline_on_nn3d(iip1,j,kkm1,q19)))
               aloline_on_nn3d(i,j,kkm1,q23) = imask_totreg3d(i,j,kkm1) * &
                  (alo_u*c18_slineu + abs_sc*(c14_intu*aloline_on_nn3d(i,j,kkm1,q10) + &
                  c15_intu*aloline_on_nn3d(i,j,kkm1,q11) + &
                  c17_intu*aloline_on_nn3d(i,j,kkm1,q13) + &
                  c18_intu*aloline_on_nn3d(i,j,kkm1,q14) + &
                  c05_intu*aloline_on_nn3d(i,j,kkm1,q1) + &
                  c06_intu*aloline_on_nn3d(i,j,kkm1,q2) + &
                  c23_intu*aloline_on_nn3d(i,j,kkm1,q19) + &
                  c24_intu*aloline_on_nn3d(i,j,kkm1,q20) + &
                  c08_intu*aloline_on_nn3d(i,j,kkm1,q4) + &
                  c26_intu*aloline_on_nn3d(i,j,kkm1,q22)))
               aloline_on_nn3d(iim1,j,kkm1,q24) = imask_totreg3d(iim1,j,kkm1) * &
                  (alo_u*c17_slineu + abs_sc*(c13_intu*aloline_on_nn3d(iim1,j,kkm1,q10) + &
                  c14_intu*aloline_on_nn3d(iim1,j,kkm1,q11) + &
                  c15_intu*aloline_on_nn3d(iim1,j,kkm1,q12) + &
                  c16_intu*aloline_on_nn3d(iim1,j,kkm1,q13) + &
                  c17_intu*aloline_on_nn3d(iim1,j,kkm1,q14) + &
                  c18_intu*aloline_on_nn3d(iim1,j,kkm1,q15) + &
                  c04_intu*aloline_on_nn3d(iim1,j,kkm1,q1) + &
                  c05_intu*aloline_on_nn3d(iim1,j,kkm1,q2) + &
                  c06_intu*aloline_on_nn3d(iim1,j,kkm1,q3) + &
                  c22_intu*aloline_on_nn3d(iim1,j,kkm1,q19) + &
                  c23_intu*aloline_on_nn3d(iim1,j,kkm1,q20) + &
                  c24_intu*aloline_on_nn3d(iim1,j,kkm1,q21) + &
                  c08_intu*aloline_on_nn3d(iim1,j,kkm1,q5) + &
                  c26_intu*aloline_on_nn3d(iim1,j,kkm1,q23)))
               aloline_on_nn3d(iip1,jjm1,kkm1,q25) = imask_totreg3d(iip1,jjm1,kkm1) * &
                  (alo_d*c03_slined + abs_sc*(c12_intu*aloline_on_nn3d(iip1,jjm1,kkm1,q10) + &
                  c15_intu*aloline_on_nn3d(iip1,jjm1,kkm1,q13) + &
                  c18_intu*aloline_on_nn3d(iip1,jjm1,kkm1,q16) + &
                  c06_intu*aloline_on_nn3d(iip1,jjm1,kkm1,q4) + &
                  c24_intu*aloline_on_nn3d(iip1,jjm1,kkm1,q22)))
               aloline_on_nn3d(i,jjm1,kkm1,q26) = imask_totreg3d(i,jjm1,kkm1) * &
                  (alo_u*c15_slineu + abs_sc*(c11_intu*aloline_on_nn3d(i,jjm1,kkm1,q10) + &
                  c12_intu*aloline_on_nn3d(i,jjm1,kkm1,q11) + &
                  c14_intu*aloline_on_nn3d(i,jjm1,kkm1,q13) + &
                  c15_intu*aloline_on_nn3d(i,jjm1,kkm1,q14) + &
                  c17_intu*aloline_on_nn3d(i,jjm1,kkm1,q16) + &
                  c18_intu*aloline_on_nn3d(i,jjm1,kkm1,q17) + &
                  c05_intu*aloline_on_nn3d(i,jjm1,kkm1,q4) + &
                  c06_intu*aloline_on_nn3d(i,jjm1,kkm1,q5) + &
                  c23_intu*aloline_on_nn3d(i,jjm1,kkm1,q22) + &
                  c24_intu*aloline_on_nn3d(i,jjm1,kkm1,q23) + &
                  c02_intu*aloline_on_nn3d(i,jjm1,kkm1,q1) + &
                  c08_intu*aloline_on_nn3d(i,jjm1,kkm1,q7) + &
                  c20_intu*aloline_on_nn3d(i,jjm1,kkm1,q19) + &
                  c26_intu*aloline_on_nn3d(i,jjm1,kkm1,q25)))
               aloline_on_nn3d(iim1,jjm1,kkm1,q27) = imask_totreg3d(iim1,jjm1,kkm1) * &
                  (alo_u*c14_slineu + abs_sc*(c10_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q10) + &
                  c11_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q11) + &
                  c12_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q12) + &
                  c13_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q13) + &
                  c14_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q14) + &
                  c15_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q15) + &
                  c16_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q16) + &
                  c17_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q17) + &
                  c18_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q18) + &
                  c04_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q4) + &
                  c05_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q5) + &
                  c06_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q6) + &
                  c22_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q22) + &
                  c23_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q23) + &
                  c24_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q24) + &
                  c02_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q2) + &
                  c08_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q8) + &
                  c20_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q20) + &
                  c26_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q26)))

!
!perform angular integration
               call calc_phinorm(vel_p, vth_p, vth_fiducial, xobs, phinorm)
               wall=wall_global*phinorm
               mintbar3d_tmp(i,j,k) = mintbar3d_tmp(i,j,k) + int3d(i,j,k)*wall
               normalization3d_tmp(i,j,k) = normalization3d_tmp(i,j,k) + wall
!
               aloline_nn3d_tmp(iip1,jjp1,kkp1,q1) = aloline_nn3d_tmp(iip1,jjp1,kkp1,q1) + wall*aloline_on_nn3d(iip1,jjp1,kkp1,q1)
               aloline_nn3d_tmp(i,jjp1,kkp1,q2) = aloline_nn3d_tmp(i,jjp1,kkp1,q2) + wall*aloline_on_nn3d(i,jjp1,kkp1,q2)
               aloline_nn3d_tmp(iim1,jjp1,kkp1,q3) = aloline_nn3d_tmp(iim1,jjp1,kkp1,q3) + wall*aloline_on_nn3d(iim1,jjp1,kkp1,q3)
               aloline_nn3d_tmp(iip1,j,kkp1,q4) = aloline_nn3d_tmp(iip1,j,kkp1,q4) + wall*aloline_on_nn3d(iip1,j,kkp1,q4)
               aloline_nn3d_tmp(i,j,kkp1,q5) = aloline_nn3d_tmp(i,j,kkp1,q5) + wall*aloline_on_nn3d(i,j,kkp1,q5)
               aloline_nn3d_tmp(iim1,j,kkp1,q6) = aloline_nn3d_tmp(iim1,j,kkp1,q6) + wall*aloline_on_nn3d(iim1,j,kkp1,q6)
               aloline_nn3d_tmp(iip1,jjm1,kkp1,q7) = aloline_nn3d_tmp(iip1,jjm1,kkp1,q7) + wall*aloline_on_nn3d(iip1,jjm1,kkp1,q7)
               aloline_nn3d_tmp(i,jjm1,kkp1,q8) = aloline_nn3d_tmp(i,jjm1,kkp1,q8) + wall*aloline_on_nn3d(i,jjm1,kkp1,q8)
               aloline_nn3d_tmp(iim1,jjm1,kkp1,q9) = aloline_nn3d_tmp(iim1,jjm1,kkp1,q9) + wall*aloline_on_nn3d(iim1,jjm1,kkp1,q9)
               aloline_nn3d_tmp(iip1,jjp1,k,q10) = aloline_nn3d_tmp(iip1,jjp1,k,q10) + wall*aloline_on_nn3d(iip1,jjp1,k,q10)
               aloline_nn3d_tmp(i,jjp1,k,q11) = aloline_nn3d_tmp(i,jjp1,k,q11) + wall*aloline_on_nn3d(i,jjp1,k,q11)
               aloline_nn3d_tmp(iim1,jjp1,k,q12) = aloline_nn3d_tmp(iim1,jjp1,k,q12) + wall*aloline_on_nn3d(iim1,jjp1,k,q12)
               aloline_nn3d_tmp(iip1,j,k,q13) = aloline_nn3d_tmp(iip1,j,k,q13) + wall*aloline_on_nn3d(iip1,j,k,q13)
               aloline_nn3d_tmp(i,j,k,q14) = aloline_nn3d_tmp(i,j,k,q14) + wall*aloline_on_nn3d(i,j,k,q14)
               aloline_nn3d_tmp(iim1,j,k,q15) = aloline_nn3d_tmp(iim1,j,k,q15) + wall*aloline_on_nn3d(iim1,j,k,q15)
               aloline_nn3d_tmp(iip1,jjm1,k,q16) = aloline_nn3d_tmp(iip1,jjm1,k,q16) + wall*aloline_on_nn3d(iip1,jjm1,k,q16) !new
               aloline_nn3d_tmp(i,jjm1,k,q17) = aloline_nn3d_tmp(i,jjm1,k,q17) + wall*aloline_on_nn3d(i,jjm1,k,q17)
               aloline_nn3d_tmp(iim1,jjm1,k,q18) = aloline_nn3d_tmp(iim1,jjm1,k,q18) + wall*aloline_on_nn3d(iim1,jjm1,k,q18)
               aloline_nn3d_tmp(iip1,jjp1,kkm1,q19) = aloline_nn3d_tmp(iip1,jjp1,kkm1,q19) + wall*aloline_on_nn3d(iip1,jjp1,kkm1,q19) !new
               aloline_nn3d_tmp(i,jjp1,kkm1,q20) = aloline_nn3d_tmp(i,jjp1,kkm1,q20) + wall*aloline_on_nn3d(i,jjp1,kkm1,q20) !new
               aloline_nn3d_tmp(iim1,jjp1,kkm1,q21) = aloline_nn3d_tmp(iim1,jjp1,kkm1,q21) + wall*aloline_on_nn3d(iim1,jjp1,kkm1,q21) !new
               aloline_nn3d_tmp(iip1,j,kkm1,q22) = aloline_nn3d_tmp(iip1,j,kkm1,q22) + wall*aloline_on_nn3d(iip1,j,kkm1,q22) !new
               aloline_nn3d_tmp(i,j,kkm1,q23) = aloline_nn3d_tmp(i,j,kkm1,q23) + wall*aloline_on_nn3d(i,j,kkm1,q23)
               aloline_nn3d_tmp(iim1,j,kkm1,q24) = aloline_nn3d_tmp(iim1,j,kkm1,q24) + wall*aloline_on_nn3d(iim1,j,kkm1,q24)
               aloline_nn3d_tmp(iip1,jjm1,kkm1,q25) = aloline_nn3d_tmp(iip1,jjm1,kkm1,q25) + wall*aloline_on_nn3d(iip1,jjm1,kkm1,q25) !new
               aloline_nn3d_tmp(i,jjm1,kkm1,q26) = aloline_nn3d_tmp(i,jjm1,kkm1,q26) + wall*aloline_on_nn3d(i,jjm1,kkm1,q26)
               aloline_nn3d_tmp(iim1,jjm1,kkm1,q27) = aloline_nn3d_tmp(iim1,jjm1,kkm1,q27) + wall*aloline_on_nn3d(iim1,jjm1,kkm1,q27)
!
!
!***********special treatment: point lies directly on boundary**********
!
             case(3)
               mueff=nn_x*x(i)+nn_y*y(j)+nn_z*z(k)
               if(mueff.ge.0.d0) then
!no interpolations required
!for alo
!                 int3d(i,j,k) = 0.d0
!                  int3d(i,j,k) = xic1
                  int3d(i,j,k)=calc_icore_gdark(z(k), sqrt(x(i)**2+y(j)**2+z(k)**2))
                  aloline_on_nn3d(i,j,k,14) = 0.d0
!
!perform angular integration (alo not required, since boundary specified)
                  vel_p = nn_x*velx3d(i,j,k) + nn_y*vely3d(i,j,k) + nn_z*velz3d(i,j,k)
                  vth_p = vth3d(i,j,k)
                  call calc_phinorm(vel_p, vth_p, vth_fiducial, xobs, phinorm)
                  wall=wall_global*phinorm
                  mintbar3d_tmp(i,j,k) = mintbar3d_tmp(i,j,k) + int3d(i,j,k)*wall
                  normalization3d_tmp(i,j,k) = normalization3d_tmp(i,j,k) + wall
               else
!same interpolation as in (1)
!
                  iim1=i-alpha
                  jjm1=j-beta
                  kkm1=k-gamma
                  iip1=i+alpha
                  jjp1=j+beta
                  kkp1=k+gamma
!
!calculate distance of upwind point to (previous) xy-plane, xz-plane and yz-plane
                  dels_xyu=(z(k)-z(kkm1))/nn_z
                  dels_xzu=(y(j)-y(jjm1))/nn_y
                  dels_yzu=(x(i)-x(iim1))/nn_x
                  dels_u=min(dels_xyu,dels_xzu,dels_yzu)
!
!calculate distance of downwind point to (next) xy-plane, xz-plane and yz-plane
                  dels_xyd=(z(kkp1)-z(k))/nn_z
                  dels_xzd=(y(jjp1)-y(j))/nn_y
                  dels_yzd=(x(iip1)-x(i))/nn_x
                  dels_d=min(dels_xyd,dels_xzd,dels_yzd)
!
!----------------------------local point--------------------------------
!
                  scont_p=scont3d(i,j,k)
                  sline_p=sline3d(i,j,k)
                  opac_p=opac3d(i,j,k)
                  opalbar_p=opalbar3d(i,j,k)
                  velx_p=velx3d(i,j,k)
                  vely_p=vely3d(i,j,k)
                  velz_p=velz3d(i,j,k)
                  vel_p = nn_x*velx_p + nn_y*vely_p + nn_z*velz_p
                  vth_p = vth3d(i,j,k)
!
!----------------------------upwind point-------------------------------
!
                  if(dels_xyu.eq.dels_u) then
!intersection with x-y plane on level k-gamma
                     x_u = x(i) - dels_u*nn_x
                     y_u = y(j) - dels_u*nn_y
                     z_u = z(kkm1)
!
                     iim2=i-2*alpha
                     jjm2=j-2*beta
!
                     call coeff2d_linecu(opac3d(iim2,jjm2,kkm1), opac3d(iim1,jjm2,kkm1), opac3d(i,jjm2,kkm1), &
                        opac3d(iim2,jjm1,kkm1), opac3d(iim1,jjm1,kkm1), opac3d(i,jjm1,kkm1), &
                        opac3d(iim2,j,kkm1),    opac3d(iim1,j,kkm1),    opac3d(i,j,kkm1), &
                        opalbar3d(iim2,jjm2,kkm1), opalbar3d(iim1,jjm2,kkm1), opalbar3d(i,jjm2,kkm1), &
                        opalbar3d(iim2,jjm1,kkm1), opalbar3d(iim1,jjm1,kkm1), opalbar3d(i,jjm1,kkm1), &
                        opalbar3d(iim2,j,kkm1),    opalbar3d(iim1,j,kkm1),    opalbar3d(i,j,kkm1), &
                        velx3d(iim2,jjm2,kkm1), velx3d(iim1,jjm2,kkm1), velx3d(i,jjm2,kkm1), &
                        velx3d(iim2,jjm1,kkm1), velx3d(iim1,jjm1,kkm1), velx3d(i,jjm1,kkm1), &
                        velx3d(iim2,j,kkm1),    velx3d(iim1,j,kkm1),    velx3d(i,j,kkm1), &
                        vely3d(iim2,jjm2,kkm1), vely3d(iim1,jjm2,kkm1), vely3d(i,jjm2,kkm1), &
                        vely3d(iim2,jjm1,kkm1), vely3d(iim1,jjm1,kkm1), vely3d(i,jjm1,kkm1), &
                        vely3d(iim2,j,kkm1),    vely3d(iim1,j,kkm1),    vely3d(i,j,kkm1), &
                        velz3d(iim2,jjm2,kkm1), velz3d(iim1,jjm2,kkm1), velz3d(i,jjm2,kkm1), &
                        velz3d(iim2,jjm1,kkm1), velz3d(iim1,jjm1,kkm1), velz3d(i,jjm1,kkm1), &
                        velz3d(iim2,j,kkm1),    velz3d(iim1,j,kkm1),    velz3d(i,j,kkm1), &
                        vth3d(iim2,jjm2,kkm1), vth3d(iim1,jjm2,kkm1), vth3d(i,jjm2,kkm1), &
                        vth3d(iim2,jjm1,kkm1), vth3d(iim1,jjm1,kkm1), vth3d(i,jjm1,kkm1), &
                        vth3d(iim2,j,kkm1),    vth3d(iim1,j,kkm1),    vth3d(i,j,kkm1), &
                        scont3d(iim2,jjm2,kkm1), scont3d(iim1,jjm2,kkm1), scont3d(i,jjm2,kkm1), &
                        scont3d(iim2,jjm1,kkm1), scont3d(iim1,jjm1,kkm1), scont3d(i,jjm1,kkm1), &
                        scont3d(iim2,j,kkm1),    scont3d(iim1,j,kkm1),    scont3d(i,j,kkm1), &
                        sline3d(iim2,jjm2,kkm1), sline3d(iim1,jjm2,kkm1), sline3d(i,jjm2,kkm1), &
                        sline3d(iim2,jjm1,kkm1), sline3d(iim1,jjm1,kkm1), sline3d(i,jjm1,kkm1), &
                        sline3d(iim2,j,kkm1),    sline3d(iim1,j,kkm1),    sline3d(i,j,kkm1), &
                        int3d(iim2,jjm2,kkm1), int3d(iim1,jjm2,kkm1), int3d(i,jjm2,kkm1), &
                        int3d(iim2,jjm1,kkm1), int3d(iim1,jjm1,kkm1), int3d(i,jjm1,kkm1), &
                        int3d(iim2,j,kkm1),    int3d(iim1,j,kkm1),    int3d(i,j,kkm1), &
                        x(iim2), x(iim1), x(i), y(jjm2), y(jjm1), y(j), x_u, y_u, &
                        c10_slineu, c11_slineu, c12_slineu, c13_slineu, c14_slineu, &
                        c15_slineu, c16_slineu, c17_slineu, c18_slineu, &
                        c10_intu, c11_intu, c12_intu, c13_intu, c14_intu, &
                        c15_intu, c16_intu, c17_intu, c18_intu, &
                        opac_u, opalbar_u, velx_u, vely_u, velz_u, vth_u, scont_u, sline_u, int_u)
!set interpolation coefficients that are not used to zero
                     c23_slineu=0.d0
                     c24_slineu=0.d0
                     c26_slineu=0.d0
                     c27_slineu=0.d0
!
                     c02_intu = 0.d0
                     c04_intu = 0.d0
                     c05_intu = 0.d0
                     c06_intu = 0.d0
                     c08_intu = 0.d0
                     c20_intu = 0.d0
                     c22_intu = 0.d0
                     c23_intu = 0.d0
                     c24_intu = 0.d0
                     c26_intu = 0.d0
!
                  elseif(dels_xzu.eq.dels_u) then
!intersection with x-z plane on level j-beta
                     x_u = x(i) - dels_u*nn_x
                     y_u = y(jjm1)
                     z_u = z(k) - dels_u*nn_z
!
                     iim2=i-2*alpha
                     kkm2=k-2*gamma

                     call coeff2d_linecu(opac3d(iim2,jjm1,kkm2), opac3d(iim1,jjm1,kkm2), opac3d(i,jjm1,kkm2), &
                        opac3d(iim2,jjm1,kkm1), opac3d(iim1,jjm1,kkm1), opac3d(i,jjm1,kkm1), &
                        opac3d(iim2,jjm1,k), opac3d(iim1,jjm1,k), opac3d(i,jjm1,k), &
                        opalbar3d(iim2,jjm1,kkm2), opalbar3d(iim1,jjm1,kkm2), opalbar3d(i,jjm1,kkm2), &
                        opalbar3d(iim2,jjm1,kkm1), opalbar3d(iim1,jjm1,kkm1), opalbar3d(i,jjm1,kkm1), &
                        opalbar3d(iim2,jjm1,k), opalbar3d(iim1,jjm1,k), opalbar3d(i,jjm1,k), &
                        velx3d(iim2,jjm1,kkm2), velx3d(iim1,jjm1,kkm2), velx3d(i,jjm1,kkm2), &
                        velx3d(iim2,jjm1,kkm1), velx3d(iim1,jjm1,kkm1), velx3d(i,jjm1,kkm1), &
                        velx3d(iim2,jjm1,k), velx3d(iim1,jjm1,k), velx3d(i,jjm1,k), &
                        vely3d(iim2,jjm1,kkm2), vely3d(iim1,jjm1,kkm2), vely3d(i,jjm1,kkm2), &
                        vely3d(iim2,jjm1,kkm1), vely3d(iim1,jjm1,kkm1), vely3d(i,jjm1,kkm1), &
                        vely3d(iim2,jjm1,k), vely3d(iim1,jjm1,k), vely3d(i,jjm1,k), &
                        velz3d(iim2,jjm1,kkm2), velz3d(iim1,jjm1,kkm2), velz3d(i,jjm1,kkm2), &
                        velz3d(iim2,jjm1,kkm1), velz3d(iim1,jjm1,kkm1), velz3d(i,jjm1,kkm1), &
                        velz3d(iim2,jjm1,k), velz3d(iim1,jjm1,k), velz3d(i,jjm1,k), &
                        vth3d(iim2,jjm1,kkm2), vth3d(iim1,jjm1,kkm2), vth3d(i,jjm1,kkm2), &
                        vth3d(iim2,jjm1,kkm1), vth3d(iim1,jjm1,kkm1), vth3d(i,jjm1,kkm1), &
                        vth3d(iim2,jjm1,k), vth3d(iim1,jjm1,k), vth3d(i,jjm1,k), &
                        scont3d(iim2,jjm1,kkm2), scont3d(iim1,jjm1,kkm2), scont3d(i,jjm1,kkm2), &
                        scont3d(iim2,jjm1,kkm1), scont3d(iim1,jjm1,kkm1), scont3d(i,jjm1,kkm1), &
                        scont3d(iim2,jjm1,k), scont3d(iim1,jjm1,k), scont3d(i,jjm1,k), &
                        sline3d(iim2,jjm1,kkm2), sline3d(iim1,jjm1,kkm2), sline3d(i,jjm1,kkm2), &
                        sline3d(iim2,jjm1,kkm1), sline3d(iim1,jjm1,kkm1), sline3d(i,jjm1,kkm1), &
                        sline3d(iim2,jjm1,k), sline3d(iim1,jjm1,k), sline3d(i,jjm1,k), &
                        int3d(iim2,jjm1,kkm2), int3d(iim1,jjm1,kkm2), int3d(i,jjm1,kkm2), &
                        int3d(iim2,jjm1,kkm1), int3d(iim1,jjm1,kkm1), int3d(i,jjm1,kkm1), &
                        int3d(iim2,jjm1,k), int3d(iim1,jjm1,k), int3d(i,jjm1,k), &
                        x(iim2), x(iim1), x(i), z(kkm2), z(kkm1), z(k), x_u, z_u, &
                        c04_slineu, c05_slineu, c06_slineu, c13_slineu, c14_slineu, &
                        c15_slineu, c22_slineu, c23_slineu, c24_slineu, &
                        c04_intu, c05_intu, c06_intu, c13_intu, c14_intu, &
                        c15_intu, c22_intu, c23_intu, c24_intu, &
                        opac_u, opalbar_u, velx_u, vely_u, velz_u, vth_u, scont_u, sline_u, int_u)
!set interpolation coefficients that are not used to zero
                     c17_slineu=0.d0
                     c18_slineu=0.d0
                     c26_slineu=0.d0
                     c27_slineu=0.d0
                     c02_intu = 0.d0
                     c08_intu = 0.d0
                     c10_intu = 0.d0
                     c11_intu = 0.d0
                     c12_intu = 0.d0
                     c16_intu = 0.d0
                     c17_intu = 0.d0
                     c18_intu = 0.d0
                     c20_intu = 0.d0
                     c26_intu = 0.d0
!
                  elseif(dels_yzu.eq.dels_u) then
!intersection with y-z plane on level i-alpha
                     x_u = x(iim1)
                     y_u = y(j) - dels_u*nn_y
                     z_u = z(k) - dels_u*nn_z
!
                     jjm2=j-2*beta
                     kkm2=k-2*gamma
!
                     call coeff2d_linecu(opac3d(iim1,jjm2,kkm2), opac3d(iim1,jjm1,kkm2), opac3d(iim1,j,kkm2), &
                        opac3d(iim1,jjm2,kkm1), opac3d(iim1,jjm1,kkm1), opac3d(iim1,j,kkm1), &
                        opac3d(iim1,jjm2,k), opac3d(iim1,jjm1,k), opac3d(iim1,j,k), &
                        opalbar3d(iim1,jjm2,kkm2), opalbar3d(iim1,jjm1,kkm2), opalbar3d(iim1,j,kkm2), &
                        opalbar3d(iim1,jjm2,kkm1), opalbar3d(iim1,jjm1,kkm1), opalbar3d(iim1,j,kkm1), &
                        opalbar3d(iim1,jjm2,k), opalbar3d(iim1,jjm1,k), opalbar3d(iim1,j,k), &
                        velx3d(iim1,jjm2,kkm2), velx3d(iim1,jjm1,kkm2), velx3d(iim1,j,kkm2), &
                        velx3d(iim1,jjm2,kkm1), velx3d(iim1,jjm1,kkm1), velx3d(iim1,j,kkm1), &
                        velx3d(iim1,jjm2,k), velx3d(iim1,jjm1,k), velx3d(iim1,j,k), &
                        vely3d(iim1,jjm2,kkm2), vely3d(iim1,jjm1,kkm2), vely3d(iim1,j,kkm2), &
                        vely3d(iim1,jjm2,kkm1), vely3d(iim1,jjm1,kkm1), vely3d(iim1,j,kkm1), &
                        vely3d(iim1,jjm2,k), vely3d(iim1,jjm1,k), vely3d(iim1,j,k), &
                        velz3d(iim1,jjm2,kkm2), velz3d(iim1,jjm1,kkm2), velz3d(iim1,j,kkm2), &
                        velz3d(iim1,jjm2,kkm1), velz3d(iim1,jjm1,kkm1), velz3d(iim1,j,kkm1), &
                        velz3d(iim1,jjm2,k), velz3d(iim1,jjm1,k), velz3d(iim1,j,k), &
                        vth3d(iim1,jjm2,kkm2), vth3d(iim1,jjm1,kkm2), vth3d(iim1,j,kkm2), &
                        vth3d(iim1,jjm2,kkm1), vth3d(iim1,jjm1,kkm1), vth3d(iim1,j,kkm1), &
                        vth3d(iim1,jjm2,k), vth3d(iim1,jjm1,k), vth3d(iim1,j,k), &
                        scont3d(iim1,jjm2,kkm2), scont3d(iim1,jjm1,kkm2), scont3d(iim1,j,kkm2), &
                        scont3d(iim1,jjm2,kkm1), scont3d(iim1,jjm1,kkm1), scont3d(iim1,j,kkm1), &
                        scont3d(iim1,jjm2,k), scont3d(iim1,jjm1,k), scont3d(iim1,j,k), &
                        sline3d(iim1,jjm2,kkm2), sline3d(iim1,jjm1,kkm2), sline3d(iim1,j,kkm2), &
                        sline3d(iim1,jjm2,kkm1), sline3d(iim1,jjm1,kkm1), sline3d(iim1,j,kkm1), &
                        sline3d(iim1,jjm2,k), sline3d(iim1,jjm1,k), sline3d(iim1,j,k), &
                        int3d(iim1,jjm2,kkm2), int3d(iim1,jjm1,kkm2), int3d(iim1,j,kkm2), &
                        int3d(iim1,jjm2,kkm1), int3d(iim1,jjm1,kkm1), int3d(iim1,j,kkm1), &
                        int3d(iim1,jjm2,k), int3d(iim1,jjm1,k), int3d(iim1,j,k), &
                        y(jjm2), y(jjm1), y(j), z(kkm2), z(kkm1), z(k), y_u, z_u, &
                        c02_slineu, c05_slineu, c08_slineu, c11_slineu, c14_slineu, &
                        c17_slineu, c20_slineu, c23_slineu, c26_slineu, &
                        c02_intu, c05_intu, c08_intu, c11_intu, c14_intu, &
                        c17_intu, c20_intu, c23_intu, c26_intu, &
                        opac_u, opalbar_u, velx_u, vely_u, velz_u, vth_u, scont_u, sline_u, int_u)
!set interpolation coefficients that are not used to zero
                     c15_slineu=0.d0
                     c18_slineu=0.d0
                     c24_slineu=0.d0
                     c27_slineu=0.d0
                     c04_intu = 0.d0
                     c06_intu = 0.d0
                     c10_intu = 0.d0
                     c12_intu = 0.d0
                     c13_intu = 0.d0
                     c15_intu = 0.d0
                     c16_intu = 0.d0
                     c18_intu = 0.d0
                     c22_intu = 0.d0
                     c24_intu = 0.d0
!
                  else
                     write(*,'(4es20.8,2l4)') dels_u, dels_xzu, dels_xyu, dels_yzu, dels_u.eq.dels_xzu, dels_xyu.eq.dels_u
                     stop 'error in fsc_linec3d: invalid dels_u'
                  endif
!
!---------------------------downwind point------------------------------
!
                  if(dels_xyd.eq.dels_d) then
!intersection with x-y plane on level kkp1
                     x_d = x(i) + dels_d*nn_x
                     y_d = y(j) + dels_d*nn_y
                     z_d = z(kkp1)
!
                     call coeff2d_linecd(opac3d(iim1,jjm1,kkp1), opac3d(i,jjm1,kkp1), opac3d(iip1,jjm1,kkp1), &
                        opac3d(iim1,j,kkp1), opac3d(i,j,kkp1), opac3d(iip1,j,kkp1), &
                        opac3d(iim1,jjp1,kkp1), opac3d(i,jjp1,kkp1), opac3d(iip1,jjp1,kkp1), &
                        opalbar3d(iim1,jjm1,kkp1), opalbar3d(i,jjm1,kkp1), opalbar3d(iip1,jjm1,kkp1), &
                        opalbar3d(iim1,j,kkp1), opalbar3d(i,j,kkp1), opalbar3d(iip1,j,kkp1), &
                        opalbar3d(iim1,jjp1,kkp1), opalbar3d(i,jjp1,kkp1), opalbar3d(iip1,jjp1,kkp1), &
                        velx3d(iim1,jjm1,kkp1), velx3d(i,jjm1,kkp1), velx3d(iip1,jjm1,kkp1), &
                        velx3d(iim1,j,kkp1), velx3d(i,j,kkp1), velx3d(iip1,j,kkp1), &
                        velx3d(iim1,jjp1,kkp1), velx3d(i,jjp1,kkp1), velx3d(iip1,jjp1,kkp1), &
                        vely3d(iim1,jjm1,kkp1), vely3d(i,jjm1,kkp1), vely3d(iip1,jjm1,kkp1), &
                        vely3d(iim1,j,kkp1), vely3d(i,j,kkp1), vely3d(iip1,j,kkp1), &
                        vely3d(iim1,jjp1,kkp1), vely3d(i,jjp1,kkp1), vely3d(iip1,jjp1,kkp1), &
                        velz3d(iim1,jjm1,kkp1), velz3d(i,jjm1,kkp1), velz3d(iip1,jjm1,kkp1), &
                        velz3d(iim1,j,kkp1), velz3d(i,j,kkp1), velz3d(iip1,j,kkp1), &
                        velz3d(iim1,jjp1,kkp1), velz3d(i,jjp1,kkp1), velz3d(iip1,jjp1,kkp1), &
                        vth3d(iim1,jjm1,kkp1), vth3d(i,jjm1,kkp1), vth3d(iip1,jjm1,kkp1), &
                        vth3d(iim1,j,kkp1), vth3d(i,j,kkp1), vth3d(iip1,j,kkp1), &
                        vth3d(iim1,jjp1,kkp1), vth3d(i,jjp1,kkp1), vth3d(iip1,jjp1,kkp1), &
                        scont3d(iim1,jjm1,kkp1), scont3d(i,jjm1,kkp1), scont3d(iip1,jjm1,kkp1), &
                        scont3d(iim1,j,kkp1), scont3d(i,j,kkp1), scont3d(iip1,j,kkp1), &
                        scont3d(iim1,jjp1,kkp1), scont3d(i,jjp1,kkp1), scont3d(iip1,jjp1,kkp1), &
                        sline3d(iim1,jjm1,kkp1), sline3d(i,jjm1,kkp1), sline3d(iip1,jjm1,kkp1), &
                        sline3d(iim1,j,kkp1), sline3d(i,j,kkp1), sline3d(iip1,j,kkp1), &
                        sline3d(iim1,jjp1,kkp1), sline3d(i,jjp1,kkp1), sline3d(iip1,jjp1,kkp1), &
                        x(iim1), x(i), x(iip1), y(jjm1), y(j), y(jjp1), x_d, y_d, &
                        c19_slined, c20_slined, c21_slined, c22_slined, c23_slined, c24_slined, &
                        c25_slined, c26_slined, c27_slined, &
                        opac_d, opalbar_d, velx_d, vely_d, velz_d, vth_d, scont_d, sline_d)
!
!set interpolation coefficients that are not used to zero
                     c03_slined = 0.d0
                     c06_slined = 0.d0
                     c07_slined = 0.d0
                     c08_slined = 0.d0
                     c09_slined = 0.d0
                     c12_slined = 0.d0
                     c15_slined = 0.d0
                     c16_slined = 0.d0
                     c17_slined = 0.d0
                     c18_slined = 0.d0
!
                  elseif(dels_xzd.eq.dels_d) then
!intersection with x-z plane on level jjp1
                     x_d = x(i) + dels_d*nn_x
                     y_d = y(jjp1)
                     z_d = z(k) + dels_d*nn_z
!
                     call coeff2d_linecd(opac3d(iim1,jjp1,kkm1), opac3d(i,jjp1,kkm1), opac3d(iip1,jjp1,kkm1), &
                        opac3d(iim1,jjp1,k), opac3d(i,jjp1,k), opac3d(iip1,jjp1,k), &
                        opac3d(iim1,jjp1,kkp1), opac3d(i,jjp1,kkp1), opac3d(iip1,jjp1,kkp1), &
                        opalbar3d(iim1,jjp1,kkm1), opalbar3d(i,jjp1,kkm1), opalbar3d(iip1,jjp1,kkm1), &
                        opalbar3d(iim1,jjp1,k), opalbar3d(i,jjp1,k), opalbar3d(iip1,jjp1,k), &
                        opalbar3d(iim1,jjp1,kkp1), opalbar3d(i,jjp1,kkp1), opalbar3d(iip1,jjp1,kkp1), &
                        velx3d(iim1,jjp1,kkm1), velx3d(i,jjp1,kkm1), velx3d(iip1,jjp1,kkm1), &
                        velx3d(iim1,jjp1,k), velx3d(i,jjp1,k), velx3d(iip1,jjp1,k), &
                        velx3d(iim1,jjp1,kkp1), velx3d(i,jjp1,kkp1), velx3d(iip1,jjp1,kkp1), &
                        vely3d(iim1,jjp1,kkm1), vely3d(i,jjp1,kkm1), vely3d(iip1,jjp1,kkm1), &
                        vely3d(iim1,jjp1,k), vely3d(i,jjp1,k), vely3d(iip1,jjp1,k), &
                        vely3d(iim1,jjp1,kkp1), vely3d(i,jjp1,kkp1), vely3d(iip1,jjp1,kkp1), &
                        velz3d(iim1,jjp1,kkm1), velz3d(i,jjp1,kkm1), velz3d(iip1,jjp1,kkm1), &
                        velz3d(iim1,jjp1,k), velz3d(i,jjp1,k), velz3d(iip1,jjp1,k), &
                        velz3d(iim1,jjp1,kkp1), velz3d(i,jjp1,kkp1), velz3d(iip1,jjp1,kkp1), &
                        vth3d(iim1,jjp1,kkm1), vth3d(i,jjp1,kkm1), vth3d(iip1,jjp1,kkm1), &
                        vth3d(iim1,jjp1,k), vth3d(i,jjp1,k), vth3d(iip1,jjp1,k), &
                        vth3d(iim1,jjp1,kkp1), vth3d(i,jjp1,kkp1), vth3d(iip1,jjp1,kkp1), &
                        scont3d(iim1,jjp1,kkm1), scont3d(i,jjp1,kkm1), scont3d(iip1,jjp1,kkm1), &
                        scont3d(iim1,jjp1,k), scont3d(i,jjp1,k), scont3d(iip1,jjp1,k), &
                        scont3d(iim1,jjp1,kkp1), scont3d(i,jjp1,kkp1), scont3d(iip1,jjp1,kkp1), &
                        sline3d(iim1,jjp1,kkm1), sline3d(i,jjp1,kkm1), sline3d(iip1,jjp1,kkm1), &
                        sline3d(iim1,jjp1,k), sline3d(i,jjp1,k), sline3d(iip1,jjp1,k), &
                        sline3d(iim1,jjp1,kkp1), sline3d(i,jjp1,kkp1), sline3d(iip1,jjp1,kkp1), &
                        x(iim1), x(i), x(iip1), z(kkm1), z(k), z(kkp1), x_d, z_d, &
                        c07_slined, c08_slined, c09_slined, c16_slined, c17_slined, c18_slined, &
                        c25_slined, c26_slined, c27_slined, &
                        opac_d, opalbar_d, velx_d, vely_d, velz_d, vth_d, scont_d, sline_d)
!
!set interpolation coefficients that are not used to zero
                     c03_slined = 0.d0
                     c06_slined = 0.d0
                     c12_slined = 0.d0
                     c15_slined = 0.d0
                     c19_slined = 0.d0
                     c20_slined = 0.d0
                     c21_slined = 0.d0
                     c22_slined = 0.d0
                     c23_slined = 0.d0
                     c24_slined = 0.d0
!
                  elseif(dels_yzd.eq.dels_d) then
!intersection with y-z plane on level iip1
                     x_d = x(iip1)
                     y_d = y(j) + dels_d*nn_y
                     z_d = z(k) + dels_d*nn_z
!
                     call coeff2d_linecd(opac3d(iip1,jjm1,kkm1), opac3d(iip1,j,kkm1), opac3d(iip1,jjp1,kkm1), &
                        opac3d(iip1,jjm1,k), opac3d(iip1,j,k), opac3d(iip1,jjp1,k), &
                        opac3d(iip1,jjm1,kkp1), opac3d(iip1,j,kkp1), opac3d(iip1,jjp1,kkp1), &
                        opalbar3d(iip1,jjm1,kkm1), opalbar3d(iip1,j,kkm1), opalbar3d(iip1,jjp1,kkm1), &
                        opalbar3d(iip1,jjm1,k), opalbar3d(iip1,j,k), opalbar3d(iip1,jjp1,k), &
                        opalbar3d(iip1,jjm1,kkp1), opalbar3d(iip1,j,kkp1), opalbar3d(iip1,jjp1,kkp1), &
                        velx3d(iip1,jjm1,kkm1), velx3d(iip1,j,kkm1), velx3d(iip1,jjp1,kkm1), &
                        velx3d(iip1,jjm1,k), velx3d(iip1,j,k), velx3d(iip1,jjp1,k), &
                        velx3d(iip1,jjm1,kkp1), velx3d(iip1,j,kkp1), velx3d(iip1,jjp1,kkp1), &
                        vely3d(iip1,jjm1,kkm1), vely3d(iip1,j,kkm1), vely3d(iip1,jjp1,kkm1), &
                        vely3d(iip1,jjm1,k), vely3d(iip1,j,k), vely3d(iip1,jjp1,k), &
                        vely3d(iip1,jjm1,kkp1), vely3d(iip1,j,kkp1), vely3d(iip1,jjp1,kkp1), &
                        velz3d(iip1,jjm1,kkm1), velz3d(iip1,j,kkm1), velz3d(iip1,jjp1,kkm1), &
                        velz3d(iip1,jjm1,k), velz3d(iip1,j,k), velz3d(iip1,jjp1,k), &
                        velz3d(iip1,jjm1,kkp1), velz3d(iip1,j,kkp1), velz3d(iip1,jjp1,kkp1), &
                        vth3d(iip1,jjm1,kkm1), vth3d(iip1,j,kkm1), vth3d(iip1,jjp1,kkm1), &
                        vth3d(iip1,jjm1,k), vth3d(iip1,j,k), vth3d(iip1,jjp1,k), &
                        vth3d(iip1,jjm1,kkp1), vth3d(iip1,j,kkp1), vth3d(iip1,jjp1,kkp1), &
                        scont3d(iip1,jjm1,kkm1), scont3d(iip1,j,kkm1), scont3d(iip1,jjp1,kkm1), &
                        scont3d(iip1,jjm1,k), scont3d(iip1,j,k), scont3d(iip1,jjp1,k), &
                        scont3d(iip1,jjm1,kkp1), scont3d(iip1,j,kkp1), scont3d(iip1,jjp1,kkp1), &
                        sline3d(iip1,jjm1,kkm1), sline3d(iip1,j,kkm1), sline3d(iip1,jjp1,kkm1), &
                        sline3d(iip1,jjm1,k), sline3d(iip1,j,k), sline3d(iip1,jjp1,k), &
                        sline3d(iip1,jjm1,kkp1), sline3d(iip1,j,kkp1), sline3d(iip1,jjp1,kkp1), &
                        y(jjm1), y(j), y(jjp1), z(kkm1), z(k), z(kkp1), y_d, z_d, &
                        c03_slined, c06_slined, c09_slined, c12_slined, c15_slined, c18_slined, &
                        c21_slined, c24_slined, c27_slined, &
                        opac_d, opalbar_d, velx_d, vely_d, velz_d, vth_d, scont_d, sline_d)
!
!set interpolation coefficients that are not used to zero
                     c07_slined = 0.d0
                     c08_slined = 0.d0
                     c16_slined = 0.d0
                     c17_slined = 0.d0
                     c19_slined = 0.d0
                     c20_slined = 0.d0
                     c22_slined = 0.d0
                     c23_slined = 0.d0
                     c25_slined = 0.d0
                     c26_slined = 0.d0
!
                  else
                     write(*,'(4es20.8,l4)') dels_d, dels_xzd, dels_xyd, dels_yzd, dels_d.eq.dels_xzd
                     stop 'error in fsc_linec3d: invalid dels_d'
                  endif
!
!--------------------------------radiative transfer---------------------
!
                  vel_u = velx_u*nn_x + vely_u*nn_y + velz_u*nn_z
                  vel_d = velx_d*nn_x + vely_d*nn_y + velz_d*nn_z
                  call fsc_linec(xobs, vth_fiducial, deltax, xcmf_max, xcmf_min, int_u, opac_u, opac_p, opac_d, &
                     opalbar_u, opalbar_p, opalbar_d, scont_u, scont_p, scont_d, &
                     sline_u, sline_p, sline_d, vel_u, vel_p, vel_d, vth_u, vth_p, vth_d, &
                     dels_u, dels_d, abs_sc, contr_sc, int_sc, alo_u, alo_p, alo_d)
                  int3d(i,j,k) = int_sc
!
                  aloline_on_nn3d(iip1,jjp1,kkp1,q1) = imask_totreg3d(iip1,jjp1,kkp1)*&
                     alo_d*c27_slined
                  aloline_on_nn3d(i,jjp1,kkp1,q2) = imask_totreg3d(i,jjp1,kkp1) * &
                     (alo_d*c26_slined + abs_sc*c26_intu*aloline_on_nn3d(i,jjp1,kkp1,q1))
                  aloline_on_nn3d(iim1,jjp1,kkp1,q3) = imask_totreg3d(iim1,jjp1,kkp1) * &
                     (alo_d*c25_slined + abs_sc*c26_intu*aloline_on_nn3d(iim1,jjp1,kkp1,q2))
                  aloline_on_nn3d(iip1,j,kkp1,q4) = imask_totreg3d(iip1,j,kkp1) * &
                     (alo_d*c24_slined + abs_sc*c24_intu*aloline_on_nn3d(iip1,j,kkp1,q1))
                  aloline_on_nn3d(i,j,kkp1,q5) = imask_totreg3d(i,j,kkp1) * &
                     (alo_d*c23_slined + abs_sc*(c23_intu*aloline_on_nn3d(i,j,kkp1,q1) + &
                     c24_intu*aloline_on_nn3d(i,j,kkp1,q2) + &
                     c26_intu*aloline_on_nn3d(i,j,kkp1,q4)))
                  aloline_on_nn3d(iim1,j,kkp1,q6) = imask_totreg3d(iim1,j,kkp1) * &
                     (alo_d*c22_slined + abs_sc*(c22_intu*aloline_on_nn3d(iim1,j,kkp1,q1) + &
                     c23_intu*aloline_on_nn3d(iim1,j,kkp1,q2) + &
                     c24_intu*aloline_on_nn3d(iim1,j,kkp1,q3) + &
                     c26_intu*aloline_on_nn3d(iim1,j,kkp1,q5)))
                  aloline_on_nn3d(iip1,jjm1,kkp1,q7) = imask_totreg3d(iip1,jjm1,kkp1) * &
                     (alo_d*c21_slined + abs_sc*c24_intu*aloline_on_nn3d(iip1,jjm1,kkp1,q4))
                  aloline_on_nn3d(i,jjm1,kkp1,q8) = imask_totreg3d(i,jjm1,kkp1) * &
                     (alo_d*c20_slined + abs_sc*(c23_intu*aloline_on_nn3d(i,jjm1,kkp1,q4) + &
                     c24_intu*aloline_on_nn3d(i,jjm1,kkp1,q5) + &
                     c20_intu*aloline_on_nn3d(i,jjm1,kkp1,q1) + &
                     c26_intu*aloline_on_nn3d(i,jjm1,kkp1,q7)))
                  aloline_on_nn3d(iim1,jjm1,kkp1,q9) = imask_totreg3d(iim1,jjm1,kkp1) * &
                     (alo_d*c19_slined + abs_sc*(c22_intu*aloline_on_nn3d(iim1,jjm1,kkp1,q4) + &
                     c23_intu*aloline_on_nn3d(iim1,jjm1,kkp1,q5) + &
                     c24_intu*aloline_on_nn3d(iim1,jjm1,kkp1,q6) + &
                     c20_intu*aloline_on_nn3d(iim1,jjm1,kkp1,q2) + &
                     c26_intu*aloline_on_nn3d(iim1,jjm1,kkp1,q8)))
                  aloline_on_nn3d(iip1,jjp1,k,q10) = imask_totreg3d(iip1,jjp1,k) * &
                     (alo_d*c18_slined + abs_sc*c18_intu*aloline_on_nn3d(iip1,jjp1,k,q1))
                  aloline_on_nn3d(i,jjp1,k,q11) = imask_totreg3d(i,jjp1,k) * &
                     (alo_d*c17_slined + abs_sc*(c17_intu*aloline_on_nn3d(i,jjp1,k,q1) + &
                     c18_intu*aloline_on_nn3d(i,jjp1,k,q2) + &
                     c26_intu*aloline_on_nn3d(i,jjp1,k,q10)))
                  aloline_on_nn3d(iim1,jjp1,k,q12) = imask_totreg3d(iim1,jjp1,k) * &
                     (alo_d*c16_slined + abs_sc*(c16_intu*aloline_on_nn3d(iim1,jjp1,k,q1) + &
                     c17_intu*aloline_on_nn3d(iim1,jjp1,k,q2) + &
                     c18_intu*aloline_on_nn3d(iim1,jjp1,k,q3) + &
                     c26_intu*aloline_on_nn3d(iim1,jjp1,k,q11)))
                  aloline_on_nn3d(iip1,j,k,q13) = imask_totreg3d(iip1,j,k) * &
                     (alo_d*c15_slined + abs_sc*(c15_intu*aloline_on_nn3d(iip1,j,k,q1) + &
                     c18_intu*aloline_on_nn3d(iip1,j,k,q4) + &
                     c24_intu*aloline_on_nn3d(iip1,j,k,q10)))
                  aloline_on_nn3d(i,j,k,q14) = imask_totreg3d(i,j,k) * &
                     (alo_p + abs_sc*(c14_intu*aloline_on_nn3d(i,j,k,q1) + &
                     c15_intu*aloline_on_nn3d(i,j,k,q2) + &
                     c17_intu*aloline_on_nn3d(i,j,k,q4) + &
                     c26_intu*aloline_on_nn3d(i,j,k,q13) + &
                     c18_intu*aloline_on_nn3d(i,j,k,q5) + &
                     c23_intu*aloline_on_nn3d(i,j,k,q10) + &
                     c24_intu*aloline_on_nn3d(i,j,k,q11)))
                  aloline_on_nn3d(iim1,j,k,q15) = imask_totreg3d(iim1,j,k) * &
                     (alo_u*c26_slineu + abs_sc*(c13_intu*aloline_on_nn3d(iim1,j,k,q1) + &
                     c14_intu*aloline_on_nn3d(iim1,j,k,q2) + &
                     c15_intu*aloline_on_nn3d(iim1,j,k,q3) + &
                     c16_intu*aloline_on_nn3d(iim1,j,k,q4) + &
                     c17_intu*aloline_on_nn3d(iim1,j,k,q5) + &
                     c18_intu*aloline_on_nn3d(iim1,j,k,q6) + &
                     c22_intu*aloline_on_nn3d(iim1,j,k,q10) + &
                     c23_intu*aloline_on_nn3d(iim1,j,k,q11) + &
                     c24_intu*aloline_on_nn3d(iim1,j,k,q12) + &
                     c26_intu*aloline_on_nn3d(iim1,j,k,q14)))
                  aloline_on_nn3d(iip1,jjm1,k,q16) = imask_totreg3d(iip1,jjm1,k) * &
                     (alo_d*c12_slined + abs_sc*(c12_intu*aloline_on_nn3d(iip1,jjm1,k,q1) + &
                     c15_intu*aloline_on_nn3d(iip1,jjm1,k,q4) + &
                     c18_intu*aloline_on_nn3d(iip1,jjm1,k,q7) + &
                     c24_intu*aloline_on_nn3d(iip1,jjm1,k,q13)))
                  aloline_on_nn3d(i,jjm1,k,q17) = imask_totreg3d(i,jjm1,k) * &
                     (alo_u*c24_slineu + abs_sc*(c11_intu*aloline_on_nn3d(i,jjm1,k,q1) + &
                     c12_intu*aloline_on_nn3d(i,jjm1,k,q2) + &
                     c14_intu*aloline_on_nn3d(i,jjm1,k,q4) + &
                     c15_intu*aloline_on_nn3d(i,jjm1,k,q5) + &
                     c17_intu*aloline_on_nn3d(i,jjm1,k,q7) + &
                     c18_intu*aloline_on_nn3d(i,jjm1,k,q8) + &
                     c23_intu*aloline_on_nn3d(i,jjm1,k,q13) + &
                     c24_intu*aloline_on_nn3d(i,jjm1,k,q14) + &
                     c20_intu*aloline_on_nn3d(i,jjm1,k,q10) + &
                     c26_intu*aloline_on_nn3d(i,jjm1,k,q16)))
                  aloline_on_nn3d(iim1,jjm1,k,q18) = imask_totreg3d(iim1,jjm1,k) * &
                     (alo_u*c23_slineu + abs_sc*(c10_intu*aloline_on_nn3d(iim1,jjm1,k,q1) + &
                     c11_intu*aloline_on_nn3d(iim1,jjm1,k,q2) + &
                     c12_intu*aloline_on_nn3d(iim1,jjm1,k,q3) + &
                     c13_intu*aloline_on_nn3d(iim1,jjm1,k,q4) + &
                     c14_intu*aloline_on_nn3d(iim1,jjm1,k,q5) + &
                     c15_intu*aloline_on_nn3d(iim1,jjm1,k,q6) + &
                     c16_intu*aloline_on_nn3d(iim1,jjm1,k,q7) + &
                     c17_intu*aloline_on_nn3d(iim1,jjm1,k,q8) + &
                     c18_intu*aloline_on_nn3d(iim1,jjm1,k,q9) + &
                     c22_intu*aloline_on_nn3d(iim1,jjm1,k,q13) + &
                     c23_intu*aloline_on_nn3d(iim1,jjm1,k,q14) + &
                     c24_intu*aloline_on_nn3d(iim1,jjm1,k,q15) + &
                     c20_intu*aloline_on_nn3d(iim1,jjm1,k,q11) + &
                     c26_intu*aloline_on_nn3d(iim1,jjm1,k,q17)))
                  aloline_on_nn3d(iip1,jjp1,kkm1,q19) = imask_totreg3d(iip1,jjp1,kkm1) * &
                     (alo_d*c09_slined + abs_sc*c18_intu*aloline_on_nn3d(iip1,jjp1,kkm1,q10))
                  aloline_on_nn3d(i,jjp1,kkm1,q20) = imask_totreg3d(i,jjp1,kkm1) * &
                     (alo_d*c08_slined + abs_sc*(c17_intu*aloline_on_nn3d(i,jjp1,kkm1,q10) + &
                     c18_intu*aloline_on_nn3d(i,jjp1,kkm1,q11) + &
                     c08_intu*aloline_on_nn3d(i,jjp1,kkm1,q1) + &
                     c26_intu*aloline_on_nn3d(i,jjp1,kkm1,q19)))
                  aloline_on_nn3d(iim1,jjp1,kkm1,q21) = imask_totreg3d(iim1,jjp1,kkm1) * &
                     (alo_d*c07_slined + abs_sc*(c16_intu*aloline_on_nn3d(iim1,jjp1,kkm1,q10) + &
                     c17_intu*aloline_on_nn3d(iim1,jjp1,kkm1,q11) + &
                     c18_intu*aloline_on_nn3d(iim1,jjp1,kkm1,q12) + &
                     c08_intu*aloline_on_nn3d(iim1,jjp1,kkm1,q2) + &
                     c26_intu*aloline_on_nn3d(iim1,jjp1,kkm1,q20)))
                  aloline_on_nn3d(iip1,j,kkm1,q22) = imask_totreg3d(iip1,j,kkm1) * &
                     (alo_d*c06_slined + abs_sc*(c15_intu*aloline_on_nn3d(iip1,j,kkm1,q10) + &
                     c18_intu*aloline_on_nn3d(iip1,j,kkm1,q13) + &
                     c06_intu*aloline_on_nn3d(iip1,j,kkm1,q1) + &
                     c24_intu*aloline_on_nn3d(iip1,j,kkm1,q19)))
                  aloline_on_nn3d(i,j,kkm1,q23) = imask_totreg3d(i,j,kkm1) * &
                     (alo_u*c18_slineu + abs_sc*(c14_intu*aloline_on_nn3d(i,j,kkm1,q10) + &
                     c15_intu*aloline_on_nn3d(i,j,kkm1,q11) + &
                     c17_intu*aloline_on_nn3d(i,j,kkm1,q13) + &
                     c18_intu*aloline_on_nn3d(i,j,kkm1,q14) + &
                     c05_intu*aloline_on_nn3d(i,j,kkm1,q1) + &
                     c06_intu*aloline_on_nn3d(i,j,kkm1,q2) + &
                     c23_intu*aloline_on_nn3d(i,j,kkm1,q19) + &
                     c24_intu*aloline_on_nn3d(i,j,kkm1,q20) + &
                     c08_intu*aloline_on_nn3d(i,j,kkm1,q4) + &
                     c26_intu*aloline_on_nn3d(i,j,kkm1,q22)))
                  aloline_on_nn3d(iim1,j,kkm1,q24) = imask_totreg3d(iim1,j,kkm1) * &
                     (alo_u*c17_slineu + abs_sc*(c13_intu*aloline_on_nn3d(iim1,j,kkm1,q10) + &
                     c14_intu*aloline_on_nn3d(iim1,j,kkm1,q11) + &
                     c15_intu*aloline_on_nn3d(iim1,j,kkm1,q12) + &
                     c16_intu*aloline_on_nn3d(iim1,j,kkm1,q13) + &
                     c17_intu*aloline_on_nn3d(iim1,j,kkm1,q14) + &
                     c18_intu*aloline_on_nn3d(iim1,j,kkm1,q15) + &
                     c04_intu*aloline_on_nn3d(iim1,j,kkm1,q1) + &
                     c05_intu*aloline_on_nn3d(iim1,j,kkm1,q2) + &
                     c06_intu*aloline_on_nn3d(iim1,j,kkm1,q3) + &
                     c22_intu*aloline_on_nn3d(iim1,j,kkm1,q19) + &
                     c23_intu*aloline_on_nn3d(iim1,j,kkm1,q20) + &
                     c24_intu*aloline_on_nn3d(iim1,j,kkm1,q21) + &
                     c08_intu*aloline_on_nn3d(iim1,j,kkm1,q5) + &
                     c26_intu*aloline_on_nn3d(iim1,j,kkm1,q23)))
                  aloline_on_nn3d(iip1,jjm1,kkm1,q25) = imask_totreg3d(iip1,jjm1,kkm1) * &
                     (alo_d*c03_slined + abs_sc*(c12_intu*aloline_on_nn3d(iip1,jjm1,kkm1,q10) + &
                     c15_intu*aloline_on_nn3d(iip1,jjm1,kkm1,q13) + &
                     c18_intu*aloline_on_nn3d(iip1,jjm1,kkm1,q16) + &
                     c06_intu*aloline_on_nn3d(iip1,jjm1,kkm1,q4) + &
                     c24_intu*aloline_on_nn3d(iip1,jjm1,kkm1,q22)))
                  aloline_on_nn3d(i,jjm1,kkm1,q26) = imask_totreg3d(i,jjm1,kkm1) * &
                     (alo_u*c15_slineu + abs_sc*(c11_intu*aloline_on_nn3d(i,jjm1,kkm1,q10) + &
                     c12_intu*aloline_on_nn3d(i,jjm1,kkm1,q11) + &
                     c14_intu*aloline_on_nn3d(i,jjm1,kkm1,q13) + &
                     c15_intu*aloline_on_nn3d(i,jjm1,kkm1,q14) + &
                     c17_intu*aloline_on_nn3d(i,jjm1,kkm1,q16) + &
                     c18_intu*aloline_on_nn3d(i,jjm1,kkm1,q17) + &
                     c05_intu*aloline_on_nn3d(i,jjm1,kkm1,q4) + &
                     c06_intu*aloline_on_nn3d(i,jjm1,kkm1,q5) + &
                     c23_intu*aloline_on_nn3d(i,jjm1,kkm1,q22) + &
                     c24_intu*aloline_on_nn3d(i,jjm1,kkm1,q23) + &
                     c02_intu*aloline_on_nn3d(i,jjm1,kkm1,q1) + &
                     c08_intu*aloline_on_nn3d(i,jjm1,kkm1,q7) + &
                     c20_intu*aloline_on_nn3d(i,jjm1,kkm1,q19) + &
                     c26_intu*aloline_on_nn3d(i,jjm1,kkm1,q25)))
                  aloline_on_nn3d(iim1,jjm1,kkm1,q27) = imask_totreg3d(iim1,jjm1,kkm1) * &
                     (alo_u*c14_slineu + abs_sc*(c10_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q10) + &
                     c11_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q11) + &
                     c12_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q12) + &
                     c13_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q13) + &
                     c14_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q14) + &
                     c15_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q15) + &
                     c16_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q16) + &
                     c17_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q17) + &
                     c18_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q18) + &
                     c04_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q4) + &
                     c05_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q5) + &
                     c06_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q6) + &
                     c22_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q22) + &
                     c23_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q23) + &
                     c24_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q24) + &
                     c02_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q2) + &
                     c08_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q8) + &
                     c20_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q20) + &
                     c26_intu*aloline_on_nn3d(iim1,jjm1,kkm1,q26)))
!
!perform angular integration
                  call calc_phinorm(vel_p, vth_p, vth_fiducial, xobs, phinorm)
                  wall=wall_global*phinorm
                  mintbar3d_tmp(i,j,k) = mintbar3d_tmp(i,j,k) + int3d(i,j,k)*wall
                  normalization3d_tmp(i,j,k) = normalization3d_tmp(i,j,k) + wall
!
                  aloline_nn3d_tmp(iip1,jjp1,kkp1,q1) = aloline_nn3d_tmp(iip1,jjp1,kkp1,q1) + wall*aloline_on_nn3d(iip1,jjp1,kkp1,q1)
                  aloline_nn3d_tmp(i,jjp1,kkp1,q2) = aloline_nn3d_tmp(i,jjp1,kkp1,q2) + wall*aloline_on_nn3d(i,jjp1,kkp1,q2)
                  aloline_nn3d_tmp(iim1,jjp1,kkp1,q3) = aloline_nn3d_tmp(iim1,jjp1,kkp1,q3) + wall*aloline_on_nn3d(iim1,jjp1,kkp1,q3)
                  aloline_nn3d_tmp(iip1,j,kkp1,q4) = aloline_nn3d_tmp(iip1,j,kkp1,q4) + wall*aloline_on_nn3d(iip1,j,kkp1,q4)
                  aloline_nn3d_tmp(i,j,kkp1,q5) = aloline_nn3d_tmp(i,j,kkp1,q5) + wall*aloline_on_nn3d(i,j,kkp1,q5)
                  aloline_nn3d_tmp(iim1,j,kkp1,q6) = aloline_nn3d_tmp(iim1,j,kkp1,q6) + wall*aloline_on_nn3d(iim1,j,kkp1,q6)
                  aloline_nn3d_tmp(iip1,jjm1,kkp1,q7) = aloline_nn3d_tmp(iip1,jjm1,kkp1,q7) + wall*aloline_on_nn3d(iip1,jjm1,kkp1,q7)
                  aloline_nn3d_tmp(i,jjm1,kkp1,q8) = aloline_nn3d_tmp(i,jjm1,kkp1,q8) + wall*aloline_on_nn3d(i,jjm1,kkp1,q8)
                  aloline_nn3d_tmp(iim1,jjm1,kkp1,q9) = aloline_nn3d_tmp(iim1,jjm1,kkp1,q9) + wall*aloline_on_nn3d(iim1,jjm1,kkp1,q9)
                  aloline_nn3d_tmp(iip1,jjp1,k,q10) = aloline_nn3d_tmp(iip1,jjp1,k,q10) + wall*aloline_on_nn3d(iip1,jjp1,k,q10)
                  aloline_nn3d_tmp(i,jjp1,k,q11) = aloline_nn3d_tmp(i,jjp1,k,q11) + wall*aloline_on_nn3d(i,jjp1,k,q11)
                  aloline_nn3d_tmp(iim1,jjp1,k,q12) = aloline_nn3d_tmp(iim1,jjp1,k,q12) + wall*aloline_on_nn3d(iim1,jjp1,k,q12)
                  aloline_nn3d_tmp(iip1,j,k,q13) = aloline_nn3d_tmp(iip1,j,k,q13) + wall*aloline_on_nn3d(iip1,j,k,q13)
                  aloline_nn3d_tmp(i,j,k,q14) = aloline_nn3d_tmp(i,j,k,q14) + wall*aloline_on_nn3d(i,j,k,q14)
                  aloline_nn3d_tmp(iim1,j,k,q15) = aloline_nn3d_tmp(iim1,j,k,q15) + wall*aloline_on_nn3d(iim1,j,k,q15)
                  aloline_nn3d_tmp(iip1,jjm1,k,q16) = aloline_nn3d_tmp(iip1,jjm1,k,q16) + wall*aloline_on_nn3d(iip1,jjm1,k,q16) !new
                  aloline_nn3d_tmp(i,jjm1,k,q17) = aloline_nn3d_tmp(i,jjm1,k,q17) + wall*aloline_on_nn3d(i,jjm1,k,q17)
                  aloline_nn3d_tmp(iim1,jjm1,k,q18) = aloline_nn3d_tmp(iim1,jjm1,k,q18) + wall*aloline_on_nn3d(iim1,jjm1,k,q18)
                  aloline_nn3d_tmp(iip1,jjp1,kkm1,q19) = aloline_nn3d_tmp(iip1,jjp1,kkm1,q19) + wall*aloline_on_nn3d(iip1,jjp1,kkm1,q19) !new
                  aloline_nn3d_tmp(i,jjp1,kkm1,q20) = aloline_nn3d_tmp(i,jjp1,kkm1,q20) + wall*aloline_on_nn3d(i,jjp1,kkm1,q20) !new
                  aloline_nn3d_tmp(iim1,jjp1,kkm1,q21) = aloline_nn3d_tmp(iim1,jjp1,kkm1,q21) + wall*aloline_on_nn3d(iim1,jjp1,kkm1,q21) !new
                  aloline_nn3d_tmp(iip1,j,kkm1,q22) = aloline_nn3d_tmp(iip1,j,kkm1,q22) + wall*aloline_on_nn3d(iip1,j,kkm1,q22) !new
                  aloline_nn3d_tmp(i,j,kkm1,q23) = aloline_nn3d_tmp(i,j,kkm1,q23) + wall*aloline_on_nn3d(i,j,kkm1,q23)
                  aloline_nn3d_tmp(iim1,j,kkm1,q24) = aloline_nn3d_tmp(iim1,j,kkm1,q24) + wall*aloline_on_nn3d(iim1,j,kkm1,q24)
                  aloline_nn3d_tmp(iip1,jjm1,kkm1,q25) = aloline_nn3d_tmp(iip1,jjm1,kkm1,q25) + wall*aloline_on_nn3d(iip1,jjm1,kkm1,q25) !new
                  aloline_nn3d_tmp(i,jjm1,kkm1,q26) = aloline_nn3d_tmp(i,jjm1,kkm1,q26) + wall*aloline_on_nn3d(i,jjm1,kkm1,q26)
                  aloline_nn3d_tmp(iim1,jjm1,kkm1,q27) = aloline_nn3d_tmp(iim1,jjm1,kkm1,q27) + wall*aloline_on_nn3d(iim1,jjm1,kkm1,q27)
!
               endif
!
             case default
!
            end select
!
         enddo
      enddo
   enddo
!
!***debug start
!close(1)
!***debug end
!
!
end subroutine fsc_linec3d
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine coeff2d_linecu_lin(opac_im1jm1, opac_ijm1, opac_im1j, opac_ij, &
   opalbar_im1jm1, opalbar_ijm1, opalbar_im1j, opalbar_ij, &
   velx_im1jm1, velx_ijm1, velx_im1j, velx_ij, &
   vely_im1jm1, vely_ijm1, vely_im1j, vely_ij, &
   velz_im1jm1, velz_ijm1, velz_im1j, velz_ij, &
   vth_im1jm1, vth_ijm1, vth_im1j, vth_ij, &
   scont_im1jm1, scont_ijm1, scont_im1j, scont_ij, &
   sline_im1jm1, sline_ijm1, sline_im1j, sline_ij, &
   int_im1jm1, int_ijm1, int_im1j, int_ij, &
   x_im1, x_i, y_jm1, y_j, x_p, y_p, &
   a_sline, b_sline, c_sline, d_sline, &
   a_inten, b_inten, c_inten, d_inten, &
   opac_p, opalbar_p, velx_p, vely_p, velz_p, vth_p, scont_p, sline_p, int_p)
!
!         interpolates opaciy, velocity components, thermal velocity,
!               line + continuum source function and intensity
!               values given on a 2d grid onto point x_p, y_p
!
!on input (f_* stands for opalbar_*, sline_* and int_*, respectivly):
!
! y_j      f_im1j--------------f_ij
!  |          |                  |
!  |          |                  |
!  |          |         x        |
!  |          |     (x_p,y_p)    |
!  |          |                  |
!y_jm1    f_im1jm1------------f_ijm1
!  |
!  --------x_im2-----------x_im1---------------x_i
!
!        x_p, y_p: coordinates of point onto which shall be interpolated
!
!on output:
!   1. interpolation coefficients for source function and intensity
!      (required for ALO calculations):
!         a_sline, b_sline, c_sline, d_sline
!         a_inten, b_inten, c_inten, d_inten
!
!      such that:
!         f_p = a*f_im1jm1 + b*f_ijm1 + c*f_im1j + d*f_ij
!
!   2. interpolated values at point p: opalbar_p, sline_p, velx_p, vely_p, velz_p,
!                                      vth_p, int_p
!
   use prog_type
!
   implicit none
!
! ... arguments
   real(dp), intent(in) :: opalbar_im1jm1, opalbar_ijm1, opalbar_im1j, opalbar_ij, &
      opac_im1jm1, opac_ijm1, opac_im1j, opac_ij, &
      velx_im1jm1, velx_ijm1, velx_im1j, velx_ij, &
      vely_im1jm1, vely_ijm1, vely_im1j, vely_ij, &
      velz_im1jm1, velz_ijm1, velz_im1j, velz_ij, &
      vth_im1jm1, vth_ijm1, vth_im1j, vth_ij, &
      scont_im1jm1, scont_ijm1, scont_im1j, scont_ij, &
      sline_im1jm1, sline_ijm1, sline_im1j, sline_ij, &
      int_im1jm1, int_ijm1, int_im1j, int_ij, &
      x_im1, x_i, y_jm1, y_j, x_p, y_p
   real(dp), intent(out) :: a_sline, b_sline, c_sline, d_sline, &
      a_inten, b_inten, c_inten, d_inten, &
      opac_p, opalbar_p, velx_p, vely_p, velz_p, vth_p, int_p, scont_p, sline_p
!
! ... local scalars
   real(dp) :: tx, ty, rdxdy
!
!-------------------------bilinear interpolation------------------------
!
!define deltax, deltay-ratios
   tx = (x_p-x_im1)/(x_i-x_im1)
   ty = (y_p-y_jm1)/(y_j-y_jm1)
!
   rdxdy=tx*ty
!
   a_sline=1.d0-tx-ty+rdxdy
   b_sline=tx-rdxdy
   c_sline=ty-rdxdy
   d_sline=rdxdy
!
   a_inten=a_sline
   b_inten=b_sline
   c_inten=c_sline
   d_inten=d_sline
!
   opac_p = a_sline*opac_im1jm1 + b_sline*opac_ijm1 + c_sline*opac_im1j + d_sline*opac_ij
   opac_p = max(opac_p,1.d-20) !avoid division by zero in radiative transport
   opalbar_p = a_sline*opalbar_im1jm1 + b_sline*opalbar_ijm1 + c_sline*opalbar_im1j + d_sline*opalbar_ij
   velx_p = a_sline*velx_im1jm1 + b_sline*velx_ijm1 + c_sline*velx_im1j + d_sline*velx_ij
   vely_p = a_sline*vely_im1jm1 + b_sline*vely_ijm1 + c_sline*vely_im1j + d_sline*vely_ij
   velz_p = a_sline*velz_im1jm1 + b_sline*velz_ijm1 + c_sline*velz_im1j + d_sline*velz_ij
   vth_p = a_sline*vth_im1jm1 + b_sline*vth_ijm1 + c_sline*vth_im1j + d_sline*vth_ij
   scont_p = a_sline*scont_im1jm1 + b_sline*scont_ijm1 + c_sline*scont_im1j + d_sline*scont_ij
!
   sline_p = a_sline*sline_im1jm1 + b_sline*sline_ijm1 + c_sline*sline_im1j + d_sline*sline_ij
   int_p = a_sline*int_im1jm1 + b_sline*int_ijm1 + c_sline*int_im1j + d_sline*int_ij
!
   return
!
end subroutine coeff2d_linecu_lin
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine coeff3d_linecu_lin(iim1, ii, jjm1, jj, kkm1, kk, x_p, y_p, z_p, &
   c14, c15, c17, c18, c23, c24, c26, c27, &
   opac_p, opalbar_p, velx_p, vely_p, velz_p, vth_p, scont_p, sline_p)
!
!         interpolates opacity, velocity components, thermal velocity,
!                     and line + continuum source function
!           given on a 3d grid onto point x_p, y_p, z_p
!
!on input:
!
!               f_im1jk--------------f_ijk
!                 /|                  / |
!                / |                 /  |
!               /  |                /   |
!              /   |               /    |
!             /    |              /     |
!            / f_im1jkm1---------/---f_ijkm1
!           /      /            /      /
!       f_im1jm1k------------f_ijm1k  /
!           |    /              |    /
!           |   /               |   /
!           |  /                |  /
!           | /                 | /
!           |/                  |/
!       f_im1jm1km1---------f_ijm1km1
!
!        x_p, y_p, z_p: coordinates of point onto which shall be interpolated
!
!on output:
!   1. interpolation coefficients for source function and intensity
!      (required for ALO calculations):
!         c14, c15, c17, c18,
!         c23, c24, c25, c26, c27
!
!   such that
!      f_p = c14*f_im1jm1km1 + c15*f_ijm1km1 +
!            c17*f_im1jkm1   + c18*f_ijkm1 +
!            c23*f_im1jm1k + c24*f_ijm1k +
!            c26*f_im1jk   + c27*f_ijk
!
!   2. interpolated values at point p: opalbar_p, sline_p
!
   use prog_type
   use dime3d, only: x, y, z, opac3d, opalbar3d, velx3d, vely3d, velz3d, vth3d, scont3d, sline3d
!***debug start
   use fund_const, only: xmsu, pi
   use params_input, only: vmin, vmax, yhe, hei, kcont, kappa0, alpha, xmloss, beta, vth_fiducial
   use params_stellar, only: sr
   use mod_opacities, only: opalbar_model_hamann, opalbar_model_kline
!***debug end
!
   implicit none
!
! ... argments
   integer(i4b), intent(in) :: iim1, ii, jjm1, jj, kkm1, kk
   real(dp), intent(in) :: x_p, y_p, z_p
   real(dp), intent(out) :: c14, c15, c17, c18, c23, c24, c26, c27, &
      opac_p, opalbar_p, velx_p, vely_p, velz_p, vth_p, sline_p, scont_p
!
! ... local scalars
   real(dp) :: tx, ty, tz, ax, ay, az
!
!***debug start
   real(dp) :: rad, vinf, xmloss_cgs, bconst, velr, rho, c1, c2, ne, opalbar_p2, sline_p2, int_p2, &
      velx_p2, vely_p2, velz_p2
!***debug end
!
!****************************debug start********************************
!----------------analytic expressions for velocity and opacity----------
!
!velocity
   vinf=vmax*1.d5
   bconst = 1.d0-(vmin/vmax)**(1.d0/beta)
   rad=sqrt(x_p**2+y_p**2+z_p**2)
   velr = vinf*(1.d0-bconst/rad)**beta
   velx_p=velr*x_p/rad / vth_fiducial
   vely_p=velr*y_p/rad / vth_fiducial
   velz_p=velr*z_p/rad / vth_fiducial
!
   xmloss_cgs=xmloss*xmsu/(365.25d0*24.d0*60.d0*60.d0)
   rho = xmloss_cgs/4.d0/pi/velr/(rad*sr)**2
!
!opalbar_p2 = opalbar_model_kline(yhe, hei, rho, kline)*sr
   opalbar_p = opalbar_model_hamann(sr, vinf, xmloss_cgs, kappa0, alpha, vth_fiducial, rad*sr, rho) * sr
!
!****************************debug end**********************************
!
!--------------------------trilinear interpolation----------------------
!
!define deltax, deltay, deltaz
   tx = (x_p-x(iim1))/(x(ii)-x(iim1))
   ty = (y_p-y(jjm1))/(y(jj)-y(jjm1))
   tz = (z_p-z(kkm1))/(z(kk)-z(kkm1))
!
   ax=1.d0-tx
   ay=1.d0-ty
   az=1.d0-tz
!
   c14 = ax*ay*az
   c15 = tx*ay*az
   c17 = ax*ty*az
   c18 = tx*ty*az
   c23 = ax*ay*tz
   c24 = tx*ay*tz
   c26 = ax*ty*tz
   c27 = tx*ty*tz

   opac_p = c14*opac3d(iim1,jjm1,kkm1) + c15*opac3d(ii,jjm1,kkm1) + c17*opac3d(iim1,jj,kkm1) + c18*opac3d(ii,jj,kkm1) + &
      c23*opac3d(iim1,jjm1,kk)   + c24*opac3d(ii,jjm1,kk)   + c26*opac3d(iim1,jj,kk)   + c27*opac3d(ii,jj,kk)
   opac_p = max(opac_p,1.d-20) !avoid division by zero in radiative transport
!opalbar_p = c14*opalbar3d(iim1,jjm1,kkm1) + c15*opalbar3d(ii,jjm1,kkm1) + c17*opalbar3d(iim1,jj,kkm1) + c18*opalbar3d(ii,jj,kkm1) + &
!            c23*opalbar3d(iim1,jjm1,kk)   + c24*opalbar3d(ii,jjm1,kk)   + c26*opalbar3d(iim1,jj,kk)   + c27*opalbar3d(ii,jj,kk)
!velx_p = c14*velx3d(iim1,jjm1,kkm1) + c15*velx3d(ii,jjm1,kkm1) + c17*velx3d(iim1,jj,kkm1) + c18*velx3d(ii,jj,kkm1) + &
!         c23*velx3d(iim1,jjm1,kk)   + c24*velx3d(ii,jjm1,kk)   + c26*velx3d(iim1,jj,kk)   + c27*velx3d(ii,jj,kk)
!vely_p = c14*vely3d(iim1,jjm1,kkm1) + c15*vely3d(ii,jjm1,kkm1) + c17*vely3d(iim1,jj,kkm1) + c18*vely3d(ii,jj,kkm1) + &
!         c23*vely3d(iim1,jjm1,kk)   + c24*vely3d(ii,jjm1,kk)   + c26*vely3d(iim1,jj,kk)   + c27*vely3d(ii,jj,kk)
!velz_p = c14*velz3d(iim1,jjm1,kkm1) + c15*velz3d(ii,jjm1,kkm1) + c17*velz3d(iim1,jj,kkm1) + c18*velz3d(ii,jj,kkm1) + &
!         c23*velz3d(iim1,jjm1,kk)   + c24*velz3d(ii,jjm1,kk)   + c26*velz3d(iim1,jj,kk)   + c27*velz3d(ii,jj,kk)
   vth_p = c14*vth3d(iim1,jjm1,kkm1) + c15*vth3d(ii,jjm1,kkm1) + c17*vth3d(iim1,jj,kkm1) + c18*vth3d(ii,jj,kkm1) + &
      c23*vth3d(iim1,jjm1,kk)   + c24*vth3d(ii,jjm1,kk)   + c26*vth3d(iim1,jj,kk)   + c27*vth3d(ii,jj,kk)
   scont_p = c14*scont3d(iim1,jjm1,kkm1) + c15*scont3d(ii,jjm1,kkm1) + c17*scont3d(iim1,jj,kkm1) + c18*scont3d(ii,jj,kkm1) + &
      c23*scont3d(iim1,jjm1,kk)   + c24*scont3d(ii,jjm1,kk)   + c26*scont3d(iim1,jj,kk)   + c27*scont3d(ii,jj,kk)
   sline_p = c14*sline3d(iim1,jjm1,kkm1) + c15*sline3d(ii,jjm1,kkm1) + c17*sline3d(iim1,jj,kkm1) + c18*sline3d(ii,jj,kkm1) + &
      c23*sline3d(iim1,jjm1,kk)   + c24*sline3d(ii,jjm1,kk)   + c26*sline3d(iim1,jj,kk)   + c27*sline3d(ii,jj,kk)


!write(*,'(12es20.8)') x_p, y_p, z_p, rad, opalbar_p2, opalbar_p, velx_p2, velx_p, vely_p2, vely_p, velz_p2, velz_p
!write(*,*)
!
   return
!
!
end subroutine coeff3d_linecu_lin
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine coeff2d_linecu(opac_im2jm2, opac_im1jm2, opac_ijm2, &
   opac_im2jm1, opac_im1jm1, opac_ijm1, &
   opac_im2j,   opac_im1j,   opac_ij, &
   opalbar_im2jm2, opalbar_im1jm2, opalbar_ijm2, &
   opalbar_im2jm1, opalbar_im1jm1, opalbar_ijm1, &
   opalbar_im2j,   opalbar_im1j,   opalbar_ij, &
   velx_im2jm2, velx_im1jm2, velx_ijm2, &
   velx_im2jm1, velx_im1jm1, velx_ijm1, &
   velx_im2j,   velx_im1j,   velx_ij, &
   vely_im2jm2, vely_im1jm2, vely_ijm2, &
   vely_im2jm1, vely_im1jm1, vely_ijm1, &
   vely_im2j,   vely_im1j,   vely_ij, &
   velz_im2jm2, velz_im1jm2, velz_ijm2, &
   velz_im2jm1, velz_im1jm1, velz_ijm1, &
   velz_im2j,   velz_im1j,   velz_ij, &
   vth_im2jm2, vth_im1jm2, vth_ijm2, &
   vth_im2jm1, vth_im1jm1, vth_ijm1, &
   vth_im2j,   vth_im1j,   vth_ij, &
   scont_im2jm2, scont_im1jm2, scont_ijm2, &
   scont_im2jm1, scont_im1jm1, scont_ijm1, &
   scont_im2j,   scont_im1j,   scont_ij, &
   sline_im2jm2, sline_im1jm2, sline_ijm2, &
   sline_im2jm1, sline_im1jm1, sline_ijm1, &
   sline_im2j,   sline_im1j,   sline_ij, &
   int_im2jm2, int_im1jm2, int_ijm2, &
   int_im2jm1, int_im1jm1, int_ijm1, &
   int_im2j,   int_im1j,   int_ij, &
   x_im2, x_im1, x_i, y_jm2, y_jm1, y_j, x_p, y_p, &
   a_sline, b_sline, c_sline, d_sline, e_sline, &
   f_sline, g_sline, h_sline, i_sline, &
   a_inten, b_inten, c_inten, d_inten, e_inten, &
   f_inten, g_inten, h_inten, i_inten, &
   opac_p, opalbar_p, velx_p, vely_p, velz_p, vth_p, scont_p, sline_p, int_p)
!
!         interpolates opaciy, velocity components, thermal velocity,
!               line + continuum source function and intensity
!               values given on a 2d grid onto point x_p, y_p
!
!on input (f_* stands for opalbar_*, sline_* velx_*, vely_*, velz_*, vth_*, and int_*, respectivly):
!
! y_j      f_im2j----------f_im1j--------------f_ij
!  |          |              |                  |
!  |          |              |                  |
!  |          |              |         x        |
!  |          |              |     (x_p,y_p)    |
!  |          |              |                  |
!y_jm1    f_im2jm1-------f_im1jm1------------f_ijm1
!  |          |              |                  |
!  |          |              |                  |
!  |          |              |                  |
!  |          |              |                  |
!  |          |              |                  |
!y_jm2    f_im2jm2-------f_im1jm2------------f_ijm2
!  |
!  --------x_im2-----------x_im1---------------x_i
!
!        x_p, y_p: coordinates of point onto which shall be interpolated
!
!on output:
!   1. interpolation coefficients for source function and intensity
!      (required for ALO calculations):
!         a_sline, b_sline, c_sline, d_sline, e_sline,
!         f_sline, g_sline, h_sline, i_sline
!         a_inten, b_inten, c_inten, d_inten, e_inten
!         f_inten, g_inten, h_inten, i_inten
!
!      such that:
!         f_p = a*f_im2jm2 + b*f_im1jm2 + c*f_ijm2 +
!               d*f_im2jm1 + e*f_im1jm1 + f*f_ijm1 +
!               g*f_im2j   + h*f_im1j   + i*f_ij
!
!   2. interpolated values at point p: opalbar_p, sline_p, velx_p, vely_p, velz_p,
!                                      vth_p, int_p, opac_p, scont_p
!
   use prog_type
   use mod_interp2d, only: wp_interp2d
!
   implicit none
!
! ... arguments
   real(dp), intent(in) :: opac_im2jm2, opac_im1jm2, opac_ijm2, &
      opac_im2jm1, opac_im1jm1, opac_ijm1, &
      opac_im2j,   opac_im1j,   opac_ij, &
      opalbar_im2jm2, opalbar_im1jm2, opalbar_ijm2, &
      opalbar_im2jm1, opalbar_im1jm1, opalbar_ijm1, &
      opalbar_im2j,   opalbar_im1j,   opalbar_ij, &
      velx_im2jm2, velx_im1jm2, velx_ijm2, &
      velx_im2jm1, velx_im1jm1, velx_ijm1, &
      velx_im2j,   velx_im1j,   velx_ij, &
      vely_im2jm2, vely_im1jm2, vely_ijm2, &
      vely_im2jm1, vely_im1jm1, vely_ijm1, &
      vely_im2j,   vely_im1j,   vely_ij, &
      velz_im2jm2, velz_im1jm2, velz_ijm2, &
      velz_im2jm1, velz_im1jm1, velz_ijm1, &
      velz_im2j,   velz_im1j,   velz_ij, &
      vth_im2jm2, vth_im1jm2, vth_ijm2, &
      vth_im2jm1, vth_im1jm1, vth_ijm1, &
      vth_im2j,   vth_im1j,   vth_ij, &
      scont_im2jm2, scont_im1jm2, scont_ijm2, &
      scont_im2jm1, scont_im1jm1, scont_ijm1, &
      scont_im2j,   scont_im1j,   scont_ij, &
      sline_im2jm2, sline_im1jm2, sline_ijm2, &
      sline_im2jm1, sline_im1jm1, sline_ijm1, &
      sline_im2j,   sline_im1j,   sline_ij, &
      int_im2jm2, int_im1jm2, int_ijm2, &
      int_im2jm1, int_im1jm1, int_ijm1, &
      int_im2j,   int_im1j,   int_ij, &
      x_im2, x_im1, x_i, y_jm2, y_jm1, y_j, x_p, y_p
   real(dp), intent(out) :: a_sline, b_sline, c_sline, d_sline, e_sline, &
      f_sline, g_sline, h_sline, i_sline, &
      a_inten, b_inten, c_inten, d_inten, e_inten, &
      f_inten, g_inten, h_inten, i_inten, &
      opac_p, opalbar_p, velx_p, vely_p, velz_p, vth_p, int_p, scont_p, sline_p
!
! ... local scalars
   real(dp) :: dxim1, dxi, dx, tx, ax, bx, cx, &
      dyjm1, dyj, dy, ty, ay, by, cy, &
      axt, bxt, cxt, axt2, bxt2, cxt2, &
      ayt, byt, cyt, ayt2, byt2, cyt2, &
      axtjm2_sline, bxtjm2_sline, cxtjm2_sline, &
      axtjm1_sline, bxtjm1_sline, cxtjm1_sline, &
      axtj_sline, bxtj_sline, cxtj_sline, &
      axjm2_sline, bxjm2_sline, cxjm2_sline, &
      axjm1_sline, bxjm1_sline, cxjm1_sline, &
      axj_sline, bxj_sline, cxj_sline, &
      ayt_sline, byt_sline, cyt_sline, &
      ay_sline, by_sline, cy_sline, &
      axtjm2_int, bxtjm2_int, cxtjm2_int, &
      axtjm1_int, bxtjm1_int, cxtjm1_int, &
      axtj_int, bxtj_int, cxtj_int, &
      axjm2_int, bxjm2_int, cxjm2_int, &
      axjm1_int, bxjm1_int, cxjm1_int, &
      axj_int, bxj_int, cxj_int, &
      ay_int, by_int, cy_int, &
      ayt_int, byt_int, cyt_int, &
      opac_jm2, opac_jm1, opac_j, opacc_jm2, opacc_jm1, opacc_j, opac_c, &
      opalbar_jm2, opalbar_jm1, opalbar_j, opalbarc_jm2, opalbarc_jm1, opalbarc_j, opalbar_c, &
      velx_jm2, velx_jm1, velx_j, velxc_jm2, velxc_jm1, velxc_j, velx_c, &
      vely_jm2, vely_jm1, vely_j, velyc_jm2, velyc_jm1, velyc_j, vely_c, &
      velz_jm2, velz_jm1, velz_j, velzc_jm2, velzc_jm1, velzc_j, velz_c, &
      vth_jm2, vth_jm1, vth_j, vthc_jm2, vthc_jm1, vthc_j, vth_c, &
      scont_jm2, scont_jm1, scont_j, scontc_jm2, scontc_jm1, scontc_j, scont_c, &
      sline_jm2, sline_jm1, sline_j, slinec_jm2, slinec_jm1, slinec_j, sline_c, &
      int_jm2, int_jm1, int_j, intc_jm2, intc_jm1, intc_j, int_c!, &

   real(dp) :: fac, fac2
!
!---------------monotonic quadratic bezier interpolation----------------
!--------(with derivatives for control point from weighted mean)--------
!--------------------where weights are assigned-------------------------
!
!define deltax, deltay
   dxim1 = x_im1-x_im2
   dxi = x_i-x_im1
   dx = dxim1+dxi
   dyjm1 = y_jm1-y_jm2
   dyj = y_j-y_jm1
   dy = dyjm1+dyj
!
!define deltax, deltay-ratios
   tx = (x_p-x_im1)/dxi
   ty = (y_p-y_jm1)/dyj
!
   ax = (1.d0-tx)**2
   bx = 2.d0*tx*(1.d0-tx)
   cx = tx**2
!
!----------------calculate control points on each j-level---------------
!   scontc_jm2, scontc_jm1, scontc_j
!   slinec_jm2, slinec_jm1, slinec_j
!   intc_jm2, intc_jm1, intc_j
!   opacc_jm2, opacc_jm1, opacc_j
!   opalbarc_jm2, opalbarc_jm1, opalbarc_j
!
!derivative weights for source function and intensity
   fac=max(wp_interp2d,dxim1/dx)
!derivative weights for velocity components and opacity
   fac2=dxim1/dx
!
   axt = (fac-1.d0)*dxi/2.d0/dxim1
   bxt = ((2.d0-fac)*dxim1 + (1.d0-fac)*dxi)/2.d0/dxim1
   cxt = fac/2.d0
   axt2 = (fac2-1.d0)*dxi/2.d0/dxim1
   bxt2 = ((2.d0-fac2)*dxim1 + (1.d0-fac2)*dxi)/2.d0/dxim1
   cxt2 = fac2/2.d0
!
   opacc_jm2 = opac_im2jm2*axt2 + opac_im1jm2*bxt2 + opac_ijm2*cxt2
   opacc_jm1 = opac_im2jm1*axt2 + opac_im1jm1*bxt2 + opac_ijm1*cxt2
   opacc_j   = opac_im2j*axt2   + opac_im1j*bxt2   + opac_ij*cxt2
!
   opalbarc_jm2 = opalbar_im2jm2*axt2 + opalbar_im1jm2*bxt2 + opalbar_ijm2*cxt2
   opalbarc_jm1 = opalbar_im2jm1*axt2 + opalbar_im1jm1*bxt2 + opalbar_ijm1*cxt2
   opalbarc_j   = opalbar_im2j*axt2   + opalbar_im1j*bxt2   + opalbar_ij*cxt2
!
   velxc_jm2 = velx_im2jm2*axt2 + velx_im1jm2*bxt2 + velx_ijm2*cxt2
   velxc_jm1 = velx_im2jm1*axt2 + velx_im1jm1*bxt2 + velx_ijm1*cxt2
   velxc_j   = velx_im2j*axt2   + velx_im1j*bxt2   + velx_ij*cxt2
!
   velyc_jm2 = vely_im2jm2*axt2 + vely_im1jm2*bxt2 + vely_ijm2*cxt2
   velyc_jm1 = vely_im2jm1*axt2 + vely_im1jm1*bxt2 + vely_ijm1*cxt2
   velyc_j   = vely_im2j*axt2   + vely_im1j*bxt2   + vely_ij*cxt2
!
   velzc_jm2 = velz_im2jm2*axt2 + velz_im1jm2*bxt2 + velz_ijm2*cxt2
   velzc_jm1 = velz_im2jm1*axt2 + velz_im1jm1*bxt2 + velz_ijm1*cxt2
   velzc_j   = velz_im2j*axt2   + velz_im1j*bxt2   + velz_ij*cxt2
!
   vthc_jm2 = vth_im2jm2*axt2 + vth_im1jm2*bxt2 + vth_ijm2*cxt2
   vthc_jm1 = vth_im2jm1*axt2 + vth_im1jm1*bxt2 + vth_ijm1*cxt2
   vthc_j   = vth_im2j*axt2   + vth_im1j*bxt2   + vth_ij*cxt2
!
   scontc_jm2 = scont_im2jm2*axt2 + scont_im1jm2*bxt2 + scont_ijm2*cxt2
   scontc_jm1 = scont_im2jm1*axt2 + scont_im1jm1*bxt2 + scont_ijm1*cxt2
   scontc_j   = scont_im2j*axt2   + scont_im1j*bxt2   + scont_ij*cxt2
!
   slinec_jm2 = sline_im2jm2*axt + sline_im1jm2*bxt + sline_ijm2*cxt
   slinec_jm1 = sline_im2jm1*axt + sline_im1jm1*bxt + sline_ijm1*cxt
   slinec_j   = sline_im2j*axt   + sline_im1j*bxt   + sline_ij*cxt
!
   intc_jm2 = int_im2jm2*axt + int_im1jm2*bxt + int_ijm2*cxt
   intc_jm1 = int_im2jm1*axt + int_im1jm1*bxt + int_ijm1*cxt
   intc_j   = int_im2j*axt   + int_im1j*bxt   + int_ij*cxt
!
!------------------ensure monotonicity on level j-----------------------
!
!velocity components (no monotonicity required)
   velx_j = ax*velx_im1j + bx*velxc_j + cx*velx_ij
   vely_j = ax*vely_im1j + bx*velyc_j + cx*vely_ij
   velz_j = ax*velz_im1j + bx*velzc_j + cx*velz_ij
!
!opacities, continuum source function and thermal velocity (no interpolation coefficients required)
   call pointcr1d_mbez(vthc_j, vth_im1j, vth_ij)
   call pointcr1d_mbez(opalbarc_j, opalbar_im1j, opalbar_ij)
   call pointcr1d_mbez(opacc_j, opac_im1j, opac_ij)
   call pointcr1d_mbez(scontc_j, scont_im1j, scont_ij)
   vth_j = ax*vth_im1j + bx*vthc_j + cx*vth_ij
   opac_j = ax*opac_im1j + bx*opacc_j + cx*opac_ij
   scont_j = ax*scont_im1j + bx*scontc_j + cx*scont_ij
   opalbar_j = ax*opalbar_im1j + bx*opalbarc_j + cx*opalbar_ij
!
!line source function and intensity
   call coeffcr1d_mbez(slinec_j, sline_im1j, sline_ij, axt, bxt, cxt, axtj_sline, bxtj_sline, cxtj_sline)
   call coeffcr1d_mbez(intc_j, int_im1j, int_ij, axt, bxt, cxt, axtj_int, bxtj_int, cxtj_int)
   axj_sline = axtj_sline*bx
   bxj_sline = bxtj_sline*bx + ax
   cxj_sline = cxtj_sline*bx + cx
   sline_j = axj_sline*sline_im2j + bxj_sline*sline_im1j + cxj_sline*sline_ij
   axj_int = axtj_int*bx
   bxj_int = bxtj_int*bx + ax
   cxj_int = cxtj_int*bx + cx
   int_j = axj_int*int_im2j + bxj_int*int_im1j + cxj_int*int_ij
!
!---------------ensure monotonicity on level j-1------------------------
!
!velocity components (no monotonicity required)
   velx_jm1 = ax*velx_im1jm1 + bx*velxc_jm1 + cx*velx_ijm1
   vely_jm1 = ax*vely_im1jm1 + bx*velyc_jm1 + cx*vely_ijm1
   velz_jm1 = ax*velz_im1jm1 + bx*velzc_jm1 + cx*velz_ijm1
!
!opacity and thermal velocity (no interpolation coefficients required)
   call pointcr1d_mbez(vthc_jm1, vth_im1jm1, vth_ijm1)
   call pointcr1d_mbez(opacc_jm1, opac_im1jm1, opac_ijm1)
   call pointcr1d_mbez(opalbarc_jm1, opalbar_im1jm1, opalbar_ijm1)
   call pointcr1d_mbez(scontc_jm1, scont_im1jm1, scont_ijm1)
   vth_jm1 = ax*vth_im1jm1 + bx*vthc_jm1 + cx*vth_ijm1
   opac_jm1 = ax*opac_im1jm1 + bx*opacc_jm1 + cx*opac_ijm1
   opalbar_jm1 = ax*opalbar_im1jm1 + bx*opalbarc_jm1 + cx*opalbar_ijm1
   scont_jm1 = ax*scont_im1jm1 + bx*scontc_jm1 + cx*scont_ijm1
!
!line source function and intensity
   call coeffcr1d_mbez(slinec_jm1, sline_im1jm1, sline_ijm1, axt, bxt, cxt, axtjm1_sline, bxtjm1_sline, cxtjm1_sline)
   call coeffcr1d_mbez(intc_jm1, int_im1jm1, int_ijm1, axt, bxt, cxt, axtjm1_int, bxtjm1_int, cxtjm1_int)
   axjm1_sline = axtjm1_sline*bx
   bxjm1_sline = bxtjm1_sline*bx + ax
   cxjm1_sline = cxtjm1_sline*bx + cx
   sline_jm1 = axjm1_sline*sline_im2jm1 + bxjm1_sline*sline_im1jm1 + cxjm1_sline*sline_ijm1
   axjm1_int = axtjm1_int*bx
   bxjm1_int = bxtjm1_int*bx + ax
   cxjm1_int = cxtjm1_int*bx + cx
   int_jm1 = axjm1_int*int_im2jm1 + bxjm1_int*int_im1jm1 + cxjm1_int*int_ijm1
!
!---------------ensure monotonicity on level j-2------------------------
!
!velocity components (no monotonicity required)
   velx_jm2 = ax*velx_im1jm2 + bx*velxc_jm2 + cx*velx_ijm2
   vely_jm2 = ax*vely_im1jm2 + bx*velyc_jm2 + cx*vely_ijm2
   velz_jm2 = ax*velz_im1jm2 + bx*velzc_jm2 + cx*velz_ijm2
!
!opacity and thermal velocity (no interpolation coefficients required)
   call pointcr1d_mbez(vthc_jm2, vth_im1jm2, vth_ijm2)
   call pointcr1d_mbez(opalbarc_jm2, opalbar_im1jm2, opalbar_ijm2)
   call pointcr1d_mbez(opacc_jm2, opac_im1jm2, opac_ijm2)
   call pointcr1d_mbez(scontc_jm2, scont_im1jm2, scont_ijm2)
   vth_jm2 = ax*vth_im1jm2 + bx*vthc_jm2 + cx*vth_ijm2
   opac_jm2 = ax*opac_im1jm2 + bx*opacc_jm2 + cx*opac_ijm2
   scont_jm2 = ax*scont_im1jm2 + bx*scontc_jm2 + cx*scont_ijm2
   opalbar_jm2 = ax*opalbar_im1jm2 + bx*opalbarc_jm2 + cx*opalbar_ijm2
!
!line source function and intensity
   call coeffcr1d_mbez(slinec_jm2, sline_im1jm2, sline_ijm2, axt, bxt, cxt, axtjm2_sline, bxtjm2_sline, cxtjm2_sline)
   call coeffcr1d_mbez(intc_jm2, int_im1jm2, int_ijm2, axt, bxt, cxt, axtjm2_int, bxtjm2_int, cxtjm2_int)
   axjm2_sline = axtjm2_sline*bx
   bxjm2_sline = bxtjm2_sline*bx + ax
   cxjm2_sline = cxtjm2_sline*bx + cx
   sline_jm2 = axjm2_sline*sline_im2jm2 + bxjm2_sline*sline_im1jm2 + cxjm2_sline*sline_ijm2
   axjm2_int = axtjm2_int*bx
   bxjm2_int = bxtjm2_int*bx + ax
   cxjm2_int = cxtjm2_int*bx + cx
   int_jm2 = axjm2_int*int_im2jm2 + bxjm2_int*int_im1jm2 + cxjm2_int*int_ijm2
!
!------------calculate control point for interpolation along y----------
!
!derivative weights for source function and intensity
   fac=max(wp_interp2d,dyjm1/dy)
!derivative weights for velocity components and opacity
   fac2=dyjm1/dy
!
   ayt = (fac-1.d0)*dyj/2.d0/dyjm1
   byt = ((2.d0-fac)*dyjm1 + (1.d0-fac)*dyj)/2.d0/dyjm1
   cyt = fac/2.d0
   ayt2 = (fac2-1.d0)*dyj/2.d0/dyjm1
   byt2 = ((2.d0-fac2)*dyjm1 + (1.d0-fac2)*dyj)/2.d0/dyjm1
   cyt2 = fac2/2.d0
!
   ay = (1.d0-ty)**2
   by = 2.d0*ty*(1.d0-ty)
   cy = ty**2
!
   velx_c = velx_jm2*ayt2 + velx_jm1*byt2 + velx_j*cyt2
   vely_c = vely_jm2*ayt2 + vely_jm1*byt2 + vely_j*cyt2
   velz_c = velz_jm2*ayt2 + velz_jm1*byt2 + velz_j*cyt2
   vth_c = vth_jm2*ayt2 + vth_jm1*byt2 + vth_j*cyt2
   opac_c = opac_jm2*ayt2 + opac_jm1*byt2 + opac_j*cyt2
   opalbar_c = opalbar_jm2*ayt2 + opalbar_jm1*byt2 + opalbar_j*cyt2
   scont_c = scont_jm2*ayt2 + scont_jm1*byt2 + scont_j*cyt2
!
   sline_c = sline_jm2*ayt + sline_jm1*byt + sline_j*cyt
   int_c = int_jm2*ayt + int_jm1*byt + int_j*cyt
!
!------------------------ensure monotonicity----------------------------
!
!velocity components (no monotonicity required)
   velx_p = ay*velx_jm1 + by*velx_c + cy*velx_j
   vely_p = ay*vely_jm1 + by*vely_c + cy*vely_j
   velz_p = ay*velz_jm1 + by*velz_c + cy*velz_j
!
!opacity and thermal velocity (no interpolation coefficients required)
   call pointcr1d_mbez(vth_c, vth_jm1, vth_j)
   call pointcr1d_mbez(opalbar_c, opalbar_jm1, opalbar_j)
   call pointcr1d_mbez(opac_c, opac_jm1, opac_j)
   call pointcr1d_mbez(scont_c, scont_jm1, scont_j)
   vth_p = ay*vth_jm1 + by*vth_c + cy*vth_j
   opac_p = ay*opac_jm1 + by*opac_c + cy*opac_j
   opac_p = max(opac_p,1.d-20) !avoid division by zero in radiative transport
   scont_p = ay*scont_jm1 + by*scont_c + cy*scont_j
   opalbar_p = ay*opalbar_jm1 + by*opalbar_c + cy*opalbar_j
!
!line source function and intensity
   call coeffcr1d_mbez(sline_c, sline_jm1, sline_j, ayt, byt, cyt, ayt_sline, byt_sline, cyt_sline)
   call coeffcr1d_mbez(int_c, int_jm1, int_j, ayt, byt, cyt, ayt_int, byt_int, cyt_int)
   ay_sline = ayt_sline*by
   by_sline = byt_sline*by + ay
   cy_sline = cyt_sline*by + cy
   sline_p = ay_sline*sline_jm2 + by_sline*sline_jm1 + cy_sline*sline_j
   ay_int = ayt_int*by
   by_int = byt_int*by + ay
   cy_int = cyt_int*by + cy
   int_p = ay_int*int_jm2 + by_int*int_jm1 + cy_int*int_j
!
   a_sline = ay_sline*axjm2_sline
   b_sline = ay_sline*bxjm2_sline
   c_sline = ay_sline*cxjm2_sline
   d_sline = by_sline*axjm1_sline
   e_sline = by_sline*bxjm1_sline
   f_sline = by_sline*cxjm1_sline
   g_sline = cy_sline*axj_sline
   h_sline = cy_sline*bxj_sline
   i_sline = cy_sline*cxj_sline
!
   a_inten = ay_int*axjm2_int
   b_inten = ay_int*bxjm2_int
   c_inten = ay_int*cxjm2_int
   d_inten = by_int*axjm1_int
   e_inten = by_int*bxjm1_int
   f_inten = by_int*cxjm1_int
   g_inten = cy_int*axj_int
   h_inten = cy_int*bxj_int
   i_inten = cy_int*cxj_int
!
!write(*,'(11es20.8)') x_p, y_p, z_p, opalbar_p, opalbar_p2, velx_p, velx_p2, vely_p, vely_p2, velz_p, velz_p2
!if(int_p.ne.int_p) then
!   write(*,*) 'error: int_p=', int_p
!   write(*,*)
!   write(*,*) int_im2jm2, int_im1jm2, int_ijm2, &
!              int_im2jm1, int_im1jm1, int_ijm1, &
!              int_im2j,   int_im1j, int_ij
!   stop
!endif
!
   return
!
end subroutine coeff2d_linecu
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine coeff2d_linecd(opac_im2jm2, opac_im1jm2, opac_ijm2, &
   opac_im2jm1, opac_im1jm1, opac_ijm1, &
   opac_im2j,   opac_im1j,   opac_ij, &
   opalbar_im2jm2, opalbar_im1jm2, opalbar_ijm2, &
   opalbar_im2jm1, opalbar_im1jm1, opalbar_ijm1, &
   opalbar_im2j,   opalbar_im1j,   opalbar_ij, &
   velx_im2jm2, velx_im1jm2, velx_ijm2, &
   velx_im2jm1, velx_im1jm1, velx_ijm1, &
   velx_im2j,   velx_im1j,   velx_ij, &
   vely_im2jm2,  vely_im1jm2,  vely_ijm2, &
   vely_im2jm1,  vely_im1jm1, vely_ijm1, &
   vely_im2j,   vely_im1j,   vely_ij, &
   velz_im2jm2, velz_im1jm2, velz_ijm2, &
   velz_im2jm1, velz_im1jm1, velz_ijm1, &
   velz_im2j,   velz_im1j,   velz_ij, &
   vth_im2jm2, vth_im1jm2, vth_ijm2, &
   vth_im2jm1, vth_im1jm1, vth_ijm1, &
   vth_im2j,   vth_im1j,   vth_ij, &
   scont_im2jm2, scont_im1jm2, scont_ijm2, &
   scont_im2jm1, scont_im1jm1, scont_ijm1, &
   scont_im2j,   scont_im1j,   scont_ij, &
   sline_im2jm2, sline_im1jm2, sline_ijm2, &
   sline_im2jm1, sline_im1jm1, sline_ijm1, &
   sline_im2j,   sline_im1j,   sline_ij, &
   x_im2, x_im1, x_i, y_jm2, y_jm1, y_j, x_p, y_p, &
   a_sline, b_sline, c_sline, d_sline, e_sline, &
   f_sline, g_sline, h_sline, i_sline, &
   opac_p, opalbar_p, velx_p, vely_p, velz_p, vth_p, scont_p, sline_p)
!
!         interpolates opaciy, velocity components, thermal velocity,
!                   line source function and intensity
!               values given on a 2d grid onto point x_p, y_p
!
!on input (f_* stands for opalbar_*, sline_* velx_*, vely_*, velz_*, vth_*, and int_*, respectivly):
!
! y_j      f_im2j----------f_im1j--------------f_ij
!  |          |              |                  |
!  |          |              |                  |
!  |          |              |         x        |
!  |          |              |     (x_p,y_p)    |
!  |          |              |                  |
!y_jm1    f_im2jm1-------f_im1jm1------------f_ijm1
!  |          |              |                  |
!  |          |              |                  |
!  |          |              |                  |
!  |          |              |                  |
!  |          |              |                  |
!y_jm2    f_im2jm2-------f_im1jm2------------f_ijm2
!  |
!  --------x_im2-----------x_im1---------------x_i
!
!        x_p, y_p: coordinates of point onto which shall be interpolated
!
!on output:
!   1. interpolation coefficients for source function
!      (required for ALO calculations):
!         a_sline, b_sline, c_sline, d_sline, e_sline
!         f_sline, g_sline, h_sline, i_sline
!
!      such that:
!         f_p = a*f_im2jm2 + b*f_im1jm2 + c*f_ijm2 +
!               d*f_im2jm1 + e*f_im1jm1 + f*f_ijm1 +
!               g*f_im2j   + h*f_im1j   + i*f_ij
!
!   2. interpolated values at point p: opalbar_p, sline_p, velx_p, vely_p, velz_p,
!                                      vth_p
!
   use prog_type
   use mod_interp2d, only: wp_interp2d
!
   implicit none
!
! ... arguments
   real(dp), intent(in) :: opac_im2jm2, opac_im1jm2, opac_ijm2, &
      opac_im2jm1, opac_im1jm1, opac_ijm1, &
      opac_im2j,   opac_im1j,   opac_ij, &
      opalbar_im2jm2, opalbar_im1jm2, opalbar_ijm2, &
      opalbar_im2jm1, opalbar_im1jm1, opalbar_ijm1, &
      opalbar_im2j,   opalbar_im1j,   opalbar_ij, &
      velx_im2jm2, velx_im1jm2, velx_ijm2, &
      velx_im2jm1, velx_im1jm1, velx_ijm1, &
      velx_im2j,   velx_im1j,   velx_ij, &
      vely_im2jm2, vely_im1jm2, vely_ijm2, &
      vely_im2jm1, vely_im1jm1, vely_ijm1, &
      vely_im2j,   vely_im1j,   vely_ij, &
      velz_im2jm2, velz_im1jm2, velz_ijm2, &
      velz_im2jm1, velz_im1jm1, velz_ijm1, &
      velz_im2j,   velz_im1j,   velz_ij, &
      vth_im2jm2, vth_im1jm2, vth_ijm2, &
      vth_im2jm1, vth_im1jm1, vth_ijm1, &
      vth_im2j,   vth_im1j,   vth_ij, &
      scont_im2jm2, scont_im1jm2, scont_ijm2, &
      scont_im2jm1, scont_im1jm1, scont_ijm1, &
      scont_im2j,   scont_im1j,   scont_ij, &
      sline_im2jm2, sline_im1jm2, sline_ijm2, &
      sline_im2jm1, sline_im1jm1, sline_ijm1, &
      sline_im2j,   sline_im1j,   sline_ij, &
      x_im2, x_im1, x_i, y_jm2, y_jm1, y_j, x_p, y_p
   real(dp), intent(out) :: a_sline, b_sline, c_sline, d_sline, e_sline, &
      f_sline, g_sline, h_sline, i_sline, &
      opalbar_p, velx_p, vely_p, velz_p, vth_p, sline_p, opac_p, scont_p
!
! ... local scalars
   real(dp) :: dxim1, dxi, dx, tx, ax, bx, cx, &
      dyjm1, dyj, dy, ty, ay, by, cy, &
      axt, bxt, cxt, axt2, bxt2, cxt2, &
      ayt, byt, cyt, ayt2, byt2, cyt2, &
      axtjm2_sline, bxtjm2_sline, cxtjm2_sline, &
      axtjm1_sline, bxtjm1_sline, cxtjm1_sline, &
      axtj_sline, bxtj_sline, cxtj_sline, &
      axjm2_sline, bxjm2_sline, cxjm2_sline, &
      axjm1_sline, bxjm1_sline, cxjm1_sline, &
      axj_sline, bxj_sline, cxj_sline, &
      ayt_sline, byt_sline, cyt_sline, &
      ay_sline, by_sline, cy_sline, &
      opac_jm2, opac_jm1, opac_j, opacc_jm2, opacc_jm1, opacc_j, opac_c, &
      opalbar_jm2, opalbar_jm1, opalbar_j, opalbarc_jm2, opalbarc_jm1, opalbarc_j, opalbar_c, &
      velx_jm2, velx_jm1, velx_j, velxc_jm2, velxc_jm1, velxc_j, velx_c, &
      vely_jm2, vely_jm1, vely_j, velyc_jm2, velyc_jm1, velyc_j, vely_c, &
      velz_jm2, velz_jm1, velz_j, velzc_jm2, velzc_jm1, velzc_j, velz_c, &
      vth_jm2, vth_jm1, vth_j, vthc_jm2, vthc_jm1, vthc_j, vth_c, &
      scont_jm2, scont_jm1, scont_j, scontc_jm2, scontc_jm1, scontc_j, scont_c, &
      sline_jm2, sline_jm1, sline_j, slinec_jm2, slinec_jm1, slinec_j, sline_c
!
   real(dp) :: fac, fac2
!
!---------------monotonic quadratic bezier interpolation----------------
!--------(with derivatives for control point from weighted mean)--------
!--------------------where weights are assigned-------------------------
!
!define deltax, deltay
   dxim1 = x_im1-x_im2
   dxi = x_i-x_im1
   dx = dxim1+dxi
   dyjm1 = y_jm1-y_jm2
   dyj = y_j-y_jm1
   dy = dyjm1+dyj
!
!define deltax, deltay-ratios
   tx = (x_p-x_im1)/dxi
   ty = (y_p-y_jm1)/dyj
!
   ax = (1.d0-tx)**2
   bx = 2.d0*tx*(1.d0-tx)
   cx = tx**2
!
!----------------calculate control points on each j-level---------------
!   slinec_jm2, slinec_jm1, slinec_j
!   opalbarc_jm2, opalbarc_jm1, opalbarc_j
!
!derivative weights for source function and intensity
   fac=max(wp_interp2d,dxim1/dx)
!derivative weights for velocity components and opacity
   fac2=dxim1/dx
!
   axt = (fac-1.d0)*dxi/2.d0/dxim1
   bxt = ((2.d0-fac)*dxim1 + (1.d0-fac)*dxi)/2.d0/dxim1
   cxt = fac/2.d0
   axt2 = (fac2-1.d0)*dxi/2.d0/dxim1
   bxt2 = ((2.d0-fac2)*dxim1 + (1.d0-fac2)*dxi)/2.d0/dxim1
   cxt2 = fac2/2.d0
!
   opalbarc_jm2 = opalbar_im2jm2*axt2 + opalbar_im1jm2*bxt2 + opalbar_ijm2*cxt2
   opalbarc_jm1 = opalbar_im2jm1*axt2 + opalbar_im1jm1*bxt2 + opalbar_ijm1*cxt2
   opalbarc_j   = opalbar_im2j*axt2   + opalbar_im1j*bxt2   + opalbar_ij*cxt2
!
   opacc_jm2 = opac_im2jm2*axt2 + opac_im1jm2*bxt2 + opac_ijm2*cxt2
   opacc_jm1 = opac_im2jm1*axt2 + opac_im1jm1*bxt2 + opac_ijm1*cxt2
   opacc_j   = opac_im2j*axt2   + opac_im1j*bxt2   + opac_ij*cxt2
!
   velxc_jm2 = velx_im2jm2*axt2 + velx_im1jm2*bxt2 + velx_ijm2*cxt2
   velxc_jm1 = velx_im2jm1*axt2 + velx_im1jm1*bxt2 + velx_ijm1*cxt2
   velxc_j   = velx_im2j*axt2   + velx_im1j*bxt2   + velx_ij*cxt2
!
   velyc_jm2 = vely_im2jm2*axt2 + vely_im1jm2*bxt2 + vely_ijm2*cxt2
   velyc_jm1 = vely_im2jm1*axt2 + vely_im1jm1*bxt2 + vely_ijm1*cxt2
   velyc_j   = vely_im2j*axt2   + vely_im1j*bxt2   + vely_ij*cxt2
!
   velzc_jm2 = velz_im2jm2*axt2 + velz_im1jm2*bxt2 + velz_ijm2*cxt2
   velzc_jm1 = velz_im2jm1*axt2 + velz_im1jm1*bxt2 + velz_ijm1*cxt2
   velzc_j   = velz_im2j*axt2   + velz_im1j*bxt2   + velz_ij*cxt2
!
   vthc_jm2 = vth_im2jm2*axt2 + vth_im1jm2*bxt2 + vth_ijm2*cxt2
   vthc_jm1 = vth_im2jm1*axt2 + vth_im1jm1*bxt2 + vth_ijm1*cxt2
   vthc_j   = vth_im2j*axt2   + vth_im1j*bxt2   + vth_ij*cxt2
!
   scontc_jm2 = scont_im2jm2*axt2 + scont_im1jm2*bxt2 + scont_ijm2*cxt2
   scontc_jm1 = scont_im2jm1*axt2 + scont_im1jm1*bxt2 + scont_ijm1*cxt2
   scontc_j   = scont_im2j*axt2   + scont_im1j*bxt2   + scont_ij*cxt2
!
   slinec_jm2 = sline_im2jm2*axt + sline_im1jm2*bxt + sline_ijm2*cxt
   slinec_jm1 = sline_im2jm1*axt + sline_im1jm1*bxt + sline_ijm1*cxt
   slinec_j   = sline_im2j*axt   + sline_im1j*bxt   + sline_ij*cxt
!
!------------------ensure monotonicity on level j-----------------------
!
!velocity components (no monotonicity required)
   velx_j = ax*velx_im1j + bx*velxc_j + cx*velx_ij
   vely_j = ax*vely_im1j + bx*velyc_j + cx*vely_ij
   velz_j = ax*velz_im1j + bx*velzc_j + cx*velz_ij
!
!opacity and thermal velocity (no interpolation coefficients required)
   call pointcr1d_mbez(vthc_j, vth_im1j, vth_ij)
   call pointcr1d_mbez(opacc_j, opac_im1j, opac_ij)
   call pointcr1d_mbez(scontc_j, scont_im1j, scont_ij)
   call pointcr1d_mbez(opalbarc_j, opalbar_im1j, opalbar_ij)
   vth_j = ax*vth_im1j + bx*vthc_j + cx*vth_ij
   opac_j = ax*opac_im1j + bx*opacc_j + cx*opac_ij
   scont_j = ax*scont_im1j + bx*scontc_j + cx*scont_ij
   opalbar_j = ax*opalbar_im1j + bx*opalbarc_j + cx*opalbar_ij
!
!line source function
   call coeffcr1d_mbez(slinec_j, sline_im1j, sline_ij, axt, bxt, cxt, axtj_sline, bxtj_sline, cxtj_sline)
   axj_sline = axtj_sline*bx
   bxj_sline = bxtj_sline*bx + ax
   cxj_sline = cxtj_sline*bx + cx
   sline_j = axj_sline*sline_im2j + bxj_sline*sline_im1j + cxj_sline*sline_ij
!
!---------------ensure monotonicity on level j-1------------------------
!
!velocity components (no monotonicity required)
   velx_jm1 = ax*velx_im1jm1 + bx*velxc_jm1 + cx*velx_ijm1
   vely_jm1 = ax*vely_im1jm1 + bx*velyc_jm1 + cx*vely_ijm1
   velz_jm1 = ax*velz_im1jm1 + bx*velzc_jm1 + cx*velz_ijm1
!
!opacity and thermal velocity (no interpolation coefficients required)
   call pointcr1d_mbez(vthc_jm1, vth_im1jm1, vth_ijm1)
   call pointcr1d_mbez(opalbarc_jm1, opalbar_im1jm1, opalbar_ijm1)
   call pointcr1d_mbez(opacc_jm1, opac_im1jm1, opac_ijm1)
   call pointcr1d_mbez(scontc_jm1, scont_im1jm1, scont_ijm1)
   vth_jm1 = ax*vth_im1jm1 + bx*vthc_jm1 + cx*vth_ijm1
   opac_jm1 = ax*opac_im1jm1 + bx*opacc_jm1 + cx*opac_ijm1
   scont_jm1 = ax*scont_im1jm1 + bx*scontc_jm1 + cx*scont_ijm1
   opalbar_jm1 = ax*opalbar_im1jm1 + bx*opalbarc_jm1 + cx*opalbar_ijm1
!
!line source function
   call coeffcr1d_mbez(slinec_jm1, sline_im1jm1, sline_ijm1, axt, bxt, cxt, axtjm1_sline, bxtjm1_sline, cxtjm1_sline)
   axjm1_sline = axtjm1_sline*bx
   bxjm1_sline = bxtjm1_sline*bx + ax
   cxjm1_sline = cxtjm1_sline*bx + cx
   sline_jm1 = axjm1_sline*sline_im2jm1 + bxjm1_sline*sline_im1jm1 + cxjm1_sline*sline_ijm1
!
!---------------ensure monotonicity on level j-2------------------------
!
!velocity components (no monotonicity required)
   velx_jm2 = ax*velx_im1jm2 + bx*velxc_jm2 + cx*velx_ijm2
   vely_jm2 = ax*vely_im1jm2 + bx*velyc_jm2 + cx*vely_ijm2
   velz_jm2 = ax*velz_im1jm2 + bx*velzc_jm2 + cx*velz_ijm2
!
!opacity and thermal velocity (no interpolation coefficients required)
   call pointcr1d_mbez(vthc_jm2, vth_im1jm2, vth_ijm2)
   call pointcr1d_mbez(opalbarc_jm2, opalbar_im1jm2, opalbar_ijm2)
   call pointcr1d_mbez(opacc_jm2, opac_im1jm2, opac_ijm2)
   call pointcr1d_mbez(scontc_jm2, scont_im1jm2, scont_ijm2)
   vth_jm2 = ax*vth_im1jm2 + bx*vthc_jm2 + cx*vth_ijm2
   opac_jm2 = ax*opac_im1jm2 + bx*opacc_jm2 + cx*opac_ijm2
   scont_jm2 = ax*scont_im1jm2 + bx*scontc_jm2 + cx*scont_ijm2
   opalbar_jm2 = ax*opalbar_im1jm2 + bx*opalbarc_jm2 + cx*opalbar_ijm2
!
!line source function and intensity
   call coeffcr1d_mbez(slinec_jm2, sline_im1jm2, sline_ijm2, axt, bxt, cxt, axtjm2_sline, bxtjm2_sline, cxtjm2_sline)
   axjm2_sline = axtjm2_sline*bx
   bxjm2_sline = bxtjm2_sline*bx + ax
   cxjm2_sline = cxtjm2_sline*bx + cx
   sline_jm2 = axjm2_sline*sline_im2jm2 + bxjm2_sline*sline_im1jm2 + cxjm2_sline*sline_ijm2
!
!------------calculate control point for interpolation along y----------
!
!derivative weights for source function
   fac=max(wp_interp2d,dyjm1/dy)
!derivative weights for velocity components and opacity
   fac2=dyjm1/dy
!
   ayt = (fac-1.d0)*dyj/2.d0/dyjm1
   byt = ((2.d0-fac)*dyjm1 + (1.d0-fac)*dyj)/2.d0/dyjm1
   cyt = fac/2.d0
   ayt2 = (fac2-1.d0)*dyj/2.d0/dyjm1
   byt2 = ((2.d0-fac2)*dyjm1 + (1.d0-fac2)*dyj)/2.d0/dyjm1
   cyt2 = fac2/2.d0
!
   ay = (1.d0-ty)**2
   by = 2.d0*ty*(1.d0-ty)
   cy = ty**2
!
   velx_c = velx_jm2*ayt2 + velx_jm1*byt2 + velx_j*cyt2
   vely_c = vely_jm2*ayt2 + vely_jm1*byt2 + vely_j*cyt2
   velz_c = velz_jm2*ayt2 + velz_jm1*byt2 + velz_j*cyt2
   vth_c = vth_jm2*ayt2 + vth_jm1*byt2 + vth_j*cyt2
   opalbar_c = opalbar_jm2*ayt2 + opalbar_jm1*byt2 + opalbar_j*cyt2
   opac_c = opac_jm2*ayt2 + opac_jm1*byt2 + opac_j*cyt2
   scont_c = scont_jm2*ayt2 + scont_jm1*byt2 + scont_j*cyt2
!
   sline_c = sline_jm2*ayt + sline_jm1*byt + sline_j*cyt
!
!------------------------ensure monotonicity----------------------------
!
!velocity components (no monotonicity required)
   velx_p = ay*velx_jm1 + by*velx_c + cy*velx_j
   vely_p = ay*vely_jm1 + by*vely_c + cy*vely_j
   velz_p = ay*velz_jm1 + by*velz_c + cy*velz_j
!
!opacity and thermal velocity (no interpolation coefficients required)
   call pointcr1d_mbez(vth_c, vth_jm1, vth_j)
   call pointcr1d_mbez(opalbar_c, opalbar_jm1, opalbar_j)
   call pointcr1d_mbez(opac_c, opac_jm1, opac_j)
   call pointcr1d_mbez(scont_c, scont_jm1, scont_j)
   vth_p = ay*vth_jm1 + by*vth_c + cy*vth_j
   opac_p = ay*opac_jm1 + by*opac_c + cy*opac_j
   opac_p = max(opac_p,1.d-20) !avoid division by zero in radiative transport
   scont_p = ay*scont_jm1 + by*scont_c + cy*scont_j
   opalbar_p = ay*opalbar_jm1 + by*opalbar_c + cy*opalbar_j
!
!line source function and intensity
   call coeffcr1d_mbez(sline_c, sline_jm1, sline_j, ayt, byt, cyt, ayt_sline, byt_sline, cyt_sline)
   ay_sline = ayt_sline*by
   by_sline = byt_sline*by + ay
   cy_sline = cyt_sline*by + cy
   sline_p = ay_sline*sline_jm2 + by_sline*sline_jm1 + cy_sline*sline_j
!
   a_sline = ay_sline*axjm2_sline
   b_sline = ay_sline*bxjm2_sline
   c_sline = ay_sline*cxjm2_sline
   d_sline = by_sline*axjm1_sline
   e_sline = by_sline*bxjm1_sline
   f_sline = by_sline*cxjm1_sline
   g_sline = cy_sline*axj_sline
   h_sline = cy_sline*bxj_sline
   i_sline = cy_sline*cxj_sline
!
!
   return


!
!
end subroutine coeff2d_linecd
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine model_debug(xp, yp, zp, velx, vely, velz, opalbar, opac)
!
!----------------analytic expressions for velocity and opacity----------
!
   use prog_type
   use fund_const, only: pi, sigmae, xmsu, cgs_mp
   use params_input, only: vmin, vmax, beta, xmloss, vth_fiducial, hei, yhe, kline, kcont, alpha, kappa0
   use params_stellar, only: sr
   use mod_opacities, only: opalbar_model_hamann, opalbar_model_kline
!
   implicit none
!
! ... arguments
   real(dp), intent(in) :: xp, yp, zp
   real(dp), intent(out) :: velx, vely, velz, opalbar, opac
!
! ... local scalars
   real(dp) :: rad, vinf, xmloss_cgs, bconst, velr, rho, c1, c2, ne
!
! ... local functions
!
!
!****************************debug start********************************
!
!velocity
   vinf=vmax*1.d5
   bconst = 1.d0-(vmin/vmax)**(1.d0/beta)
   rad=sqrt(xp**2+yp**2+zp**2)
!rad=max(rad,bconst)   !for downwind point: inside the star
   rad=max(rad,1.d0) !if downwind point inside star
   velr = vinf*(1.d0-bconst/rad)**beta
   velx=velr*xp/rad / vth_fiducial
   vely=velr*yp/rad / vth_fiducial
   velz=velr*zp/rad / vth_fiducial
!
   xmloss_cgs=xmloss*xmsu/(365.25d0*24.d0*60.d0*60.d0)
   rho = xmloss_cgs/4.d0/pi/velr/(rad*sr)**2
!
!opalbar = opalbar_model_kline(yhe, hei, rho, kline)*sr
   opalbar = opalbar_model_hamann(sr, vinf, xmloss_cgs, kappa0, alpha, vth_fiducial, rad*sr, rho) * sr
!
   c1=(1.d0+4.d0*yhe)*cgs_mp
   c2=(1.d0+hei*yhe)/c1
   ne=c2*rho
   opac=sigmae*ne*kcont*(velr/vinf)**alpha/sr
   opac=0.d0
!
end subroutine model_debug
