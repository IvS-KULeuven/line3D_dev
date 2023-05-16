!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine fsc_cont2d(muindx)
!
!-----------------------------------------------------------------------
!------short characteristics for continuum radiative transfer in 2d-----
!---------------calculating intensties for given mu---------------------
!-----------------------------------------------------------------------
!
   use prog_type

   use fund_const
   use dime3d, only: ndxmax, ndzmax, ndymax, x, y, z
   use dime3d, only: int3d, opac3d, scont3d, alocont_nn3d, imask3d, imask_innreg3d
   use angles, only: nodes_mu
   use bcondition, only: xic1
   use mod_debug, only: indx1, indx2
   use mod_interp1d, only: interpol_yp, interpol_ypl, interpol_typ_quad2b, &
      interpol_typ_quad2b, interpol_typ_quad3b
   use mod_interp2d, only: interpol2d_4p_lin
!
   implicit none
!
! ... arguments
   integer(i4b), intent(in) :: muindx
   real(dp) :: n_x, n_z
!
! ... local scalars
   integer(i4b) :: i, j, k
   integer(i4b) :: alpha, gamma
   integer(i4b) :: startx, startz, endx, endz
   integer(i4b) :: iim2, iim1, ii, iip1
   real(dp) ::  mueff
   real(dp) :: opac_u, scont_u, dels_u, dels_xu, dels_zu, int_u, &
      opac_d, scont_d, dels_d, dels_xd, dels_zd, &
      opac_u1, opac_u2, scont_u1, scont_u2, &
      opac_d1, opac_d2, scont_d1, scont_d2, &
      opac_p, scont_p, dels_r
   real(dp) :: x_u, x_u1, x_u2, z_u, z_u1, z_u2, t_u1, t_u2, t_u, &
      x_d, x_d1, x_d2, z_d, z_d1, z_d2, t_d1, t_d2, t_d
   real(dp) :: abs_sc, int_sc, contr_sc, alo_p, alo_u, alo_d
   real(dp) :: acoeff, bcoeff, ccoeff, dcoeff
   real(dp) :: ts1, te1, ts2, te2
!
! ... for debugging
!
! ... local arrays
!
! ... local functions
   real(dp) :: dist_sphere
!
! ... local characters
   character(len=50) :: enter
!
! ... local logicals
!
   n_x=sqrt(1.d0 - nodes_mu(muindx)*nodes_mu(muindx))
   n_z=nodes_mu(muindx)
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
!         if n_z >= 0                 if n_z < 0
!                startz = 2                  startz = ndzmax-1
!                endz = ndzmax-1             endz = 2
!                gamma = 1                   gamma = -1
!
!could make an own array for this!!!
   if(n_x.gt.0.d0) then
      startx = 3
      endx = ndxmax-1
      alpha=  1
   elseif(n_x.lt.0.d0) then
      startx = ndxmax-2
      endx = 2
      alpha=-1
   else
      stop 'error in fsc_cont2d: n_x = 0 not allowed'
   endif
!
   if(n_z.gt.0.d0) then
      startz = 3
      endz = ndzmax-1
      gamma=  1
   elseif(n_z.lt.0.d0) then
      startz = ndzmax-2
      endz = 2
      gamma=-1
   else
      stop 'error in fsc_cont2d: n_z = 0 not allowed'
   endif
!
!-----------------------reset the intensities---------------------------
!
!call cpu_time(ts1)
   call set_boundary2d(muindx)
!call cpu_time(te1)
!
!-----------------------------------------------------------------------
!
   alocont_nn3d=0.d0
!
!-----------------------------------------------------------------------
!
   j=ndymax/2+1

!call cpu_time(ts2)
!
   do k=startz, endz, gamma
      do i=startx, endx, alpha

         select case(imask3d(i,j,k))

!
!*************************standard rt procedure*************************
!
          case(1)
!
!calculate distance to intersection point with upwind x-grid-coordinate
            dels_xu=(x(i)-x(i-alpha))/n_x
!calculate distance to instersection point with upwind z-grid-coordinate
            dels_zu=(z(k)-z(k-gamma))/n_z
            dels_u=min(dels_xu,dels_zu)
!calculate distance to intersection point with downwind x-grid-coordinate
            dels_xd=(x(i+alpha)-x(i))/n_x
!calculate distance to instersection point with downwind z-grid-coordinate
            dels_zd=(z(k+gamma)-z(k))/n_z
            dels_d=min(dels_xd,dels_zd)
!
!----------------------------local point--------------------------------
!
            opac_p=opac3d(i,j,k)
            scont_p=scont3d(i,j,k)
!
!----------------------------upwind point-------------------------------
!
            if(dels_xu.eq.dels_u) then
!upwind point intersects on z-level
               x_u=x(i)-dels_u*n_x
               z_u=z(k)-dels_u*n_z
!interpolation with bezier spline or with log-log
               opac_u = interpol_ypl(z(k), z(k-gamma), opac3d(i-alpha,j,k), opac3d(i-alpha,j,k-gamma), z_u)
!               scont_u = interpol_typ_quad3b(scont3d(i-alpha,j,k-2*gamma), scont3d(i-alpha,j,k-gamma), scont3d(i-alpha,j,k), &
!                                             z(k-2*gamma), z(k-gamma), z(k), z_u)
               int_u = interpol_typ_quad2b(int3d(i-alpha,j,k-2*gamma), int3d(i-alpha,j,k-gamma), &
                  int3d(i-alpha,j,k), z(k-2*gamma), z(k-gamma), z(k), z_u)
               scont_u = interpol_yp(z(k), z(k-gamma), scont3d(i-alpha,j,k), scont3d(i-alpha,j,k-gamma), z_u)
!               int_u = interpol_yp(z(k), z(k-gamma), int3d(i-alpha,j,k), int3d(i-alpha,j,k-gamma), z_u)
            else
!upwind point intersects on x-level
               x_u=x(i)-dels_u*n_x
               z_u=z(k)-dels_u*n_z
!interpolation with bezier spline or with log-log
               opac_u = interpol_ypl(x(i), x(i-alpha), opac3d(i,j,k-gamma), opac3d(i-alpha,j,k-gamma), x_u)
!               scont_u = interpol_typ_quad3b(scont3d(i-2*alpha,j,k-gamma), scont3d(i-alpha,j,k-gamma), scont3d(i,j,k-gamma), &
!                                             x(i-2*alpha), x(i-alpha), x(i), x_u)
               int_u = interpol_typ_quad2b(int3d(i-2*alpha,j,k-gamma), int3d(i-alpha,j,k-gamma), &
                  int3d(i,j,k-gamma), x(i-2*alpha), x(i-alpha), x(i), x_u)
               scont_u = interpol_yp(x(i), x(i-alpha), scont3d(i,j,k-gamma), scont3d(i-alpha,j,k-gamma), x_u)
!               int_u = interpol_yp(x(i), x(i-alpha), int3d(i,j,k-gamma), int3d(i-alpha,j,k-gamma), x_u)
            endif
!
!--------------------------downwind point-------------------------------
!
            if(dels_xd.eq.dels_d) then
!downwind point intersects on z-level
               x_d=x(i)+dels_d*n_x
               z_d=z(k)+dels_d*n_z
               opac_d = interpol_ypl(z(k+gamma), z(k), opac3d(i+alpha,j,k+gamma), opac3d(i+alpha,j,k), z_d)
!               scont_d = interpol_typ_quad3b(scont3d(i+alpha,j,k-gamma), scont3d(i+alpha,j,k), scont3d(i+alpha,j,k+gamma), &
!                                             z(k-gamma), z(k), z(k+gamma), z_d)
               scont_d = interpol_yp(z(k+gamma), z(k), scont3d(i+alpha,j,k+gamma), scont3d(i+alpha,j,k), z_d)
            else
!downwind point intersects on x-level
               dels_d=dels_zd
               x_d=x(i)+dels_d*n_x
               z_d=z(k)+dels_d*n_z
               opac_d = interpol_ypl(x(i+alpha), x(i), opac3d(i+alpha,j,k+gamma), opac3d(i,j,k+gamma), x_d)
!               scont_d = interpol_typ_quad3b(scont3d(i-alpha,j,k+gamma), scont3d(i,j,k+gamma), scont3d(i+alpha,j,k+gamma), &
!                                             x(i-alpha), x(i), x(i+alpha), x_d)
               scont_d = interpol_ypl(x(i+alpha), x(i), scont3d(i+alpha,j,k+gamma), scont3d(i,j,k+gamma), x_d)
            endif
!
!-----------------------radiative transfer------------------------------
!
            call fsc_cont(int_u, opac_u, opac_p, opac_d, scont_u, scont_p, scont_d, &
               dels_u, dels_d, abs_sc, contr_sc, int_sc, alo_p, alo_u, alo_d)
!            call fsc_cont_lin(int_u, opac_u, opac_p, opac_d, scont_u, scont_p, scont_d, &
!                                dels_u, dels_d, abs_sc, contr_sc, int_sc, alo_u, alo_p)
            int3d(i,j,k) = int_sc
            alocont_nn3d(i,j,k,14) = alo_p
!
!***********special treatment: point is in vicinity of boundary*********
!
          case(2)
!
!calculate distance to intersection point with upwind x-grid-coordinate
            dels_xu=(x(i)-x(i-alpha))/n_x
!calculate distance to instersection point with upwind z-grid-coordinate
            dels_zu=(z(k)-z(k-gamma))/n_z
!calculate distance to intersection point with downwind x-grid-coordinate
            dels_xd=(x(i+alpha)-x(i))/n_x
!calculate distance to instersection point with downwind z-grid-coordinate
            dels_zd=(z(k+gamma)-z(k))/n_z
!calculate distance to the boundary (here: sphere)
            dels_r=dist_sphere(n_x,0.d0,n_z,x(i),0.d0,z(k))
!
            dels_u=min(dels_xu,dels_zu,dels_r)
            dels_d=min(dels_xd,dels_zd)
!
!----------------------------local point--------------------------------
!
            opac_p=opac3d(i,j,k)
            scont_p=scont3d(i,j,k)
!
!----------------------------upwind point-------------------------------
!
            if(dels_r.eq.dels_u) then
!
!upwind point intersects with sphere
!               indx1=indx1+1
               x_u=x(i)-dels_r*n_x
               z_u=z(k)-dels_r*n_z
               int_u=xic1
!               opac_u=opac3d(ndxmax/2+1,ndymax/2+1,64)
               scont_u=xic1
!               scont_u=interpol2d_4p_lin(scont3d(i-alpha,j,k-gamma), scont3d(i,j,k-gamma), &
!                                         scont3d(i-alpha,j,k), scont3d(i,j,k), &
!                                         x(i-alpha), x(i), z(k-gamma), z(k), x_u, z_u)
!               call coeff2d_4p_lin(x(i-alpha), x(i), z(k-gamma), z(k), x_u, z_u, &
!                                   acoeff, bcoeff, ccoeff, dcoeff)
!               scont_u = acoeff*scont3d(i-alpha,j,k-gamma) + bcoeff*scont3d(i,j,k-gamma) + &
!                         ccoeff*scont3d(i-alpha,j,k) + dcoeff*scont3d(i,j,k)
               opac_u = interpol2d_4p_lin(opac3d(i-alpha,j,k-gamma), opac3d(i,j,k-gamma), &
                  opac3d(i-alpha,j,k), opac3d(i,j,k), &
                  x(i-alpha), x(i), z(k-gamma), z(k), x_u, z_u)
            elseif(dels_xu.eq.dels_u) then
!               indx2=indx2+1
!upwind point intersects on z-level
               x_u=x(i)-dels_u*n_x
               z_u=z(k)-dels_u*n_z
!interpolation with bezier spline or with log-log
               opac_u = interpol_ypl(z(k), z(k-gamma), opac3d(i-alpha,j,k), opac3d(i-alpha,j,k-gamma), z_u)
!               scont_u = interpol_typ_quad3b(scont3d(i-alpha,j,k-2*gamma), scont3d(i-alpha,j,k-gamma), scont3d(i-alpha,j,k), &
!                                             z(k-2*gamma), z(k-gamma), z(k), z_u)
               int_u = interpol_typ_quad2b(int3d(i-alpha,j,k-2*gamma), int3d(i-alpha,j,k-gamma), &
                  int3d(i-alpha,j,k), z(k-2*gamma), z(k-gamma), z(k), z_u)
               scont_u = interpol_yp(z(k), z(k-gamma), scont3d(i-alpha,j,k), scont3d(i-alpha,j,k-gamma), z_u)
!               int_u = interpol_yp(z(k), z(k-gamma), int3d(i-alpha,j,k), int3d(i-alpha,j,k-gamma), z_u)
               dcoeff=0.d0
!***debug start
!               if(imask_innreg3d(i-alpha,j,k-gamma).eq.1) then
!                  x_u1=x(i-alpha)
!                  z_u1=sign(1.,z_u)*sqrt(1.d0-x_u1**2)
!                  mueff = n_x*x_u1 + n_z*z_u1
!                  if(mueff.ge.0.d0) then
!                     int_u = interpol_yp(z(k), z_u1, int3d(i-alpha,j,k), xic1, z_u)
!                  else
!                     int_u = interpol_yp(z(k), z_u1, int3d(i-alpha,j,k), 0.d0, z_u)
!!                     int_u = interpol_yp(z(k), z_u1, int3d(i-alpha,j,k), xic1, z_u)
!                  endif
!               elseif(imask_innreg3d(i-alpha,j,k).eq.1) then
!                  x_u1=x(i-alpha)
!                  z_u1=sign(1.,z_u)*sqrt(1.d0-x_u1**2)
!                  mueff = n_x*x_u1 + n_z*z_u1
!                  if(mueff.ge.0.d0) then
!                     int_u = interpol_yp(z(k-gamma), z_u1, int3d(i-alpha,j,k-gamma), xic1, z_u)
!                  else
!                     int_u = interpol_yp(z(k-gamma), z_u1, int3d(i-alpha,j,k-gamma), 0.d0, z_u)
!!                     int_u = interpol_yp(z(k-gamma), z_u1, int3d(i-alpha,j,k-gamma), xic1, z_u)
!                  endif
!               endif
!***debug end
            else
!              indx2=indx2+1
!upwind point intersects on x-level
               x_u=x(i)-dels_u*n_x
               z_u=z(k)-dels_u*n_z
!interpolation with bezier spline or with log-log
               opac_u = interpol_ypl(x(i), x(i-alpha), opac3d(i,j,k-gamma), opac3d(i-alpha,j,k-gamma), x_u)
!               scont_u = interpol_typ_quad3b(scont3d(i-2*alpha,j,k-gamma), scont3d(i-alpha,j,k-gamma), scont3d(i,j,k-gamma), &
!                                             x(i-2*alpha), x(i-alpha), x(i), x_u)
               int_u = interpol_typ_quad2b(int3d(i-2*alpha,j,k-gamma), int3d(i-alpha,j,k-gamma), &
                  int3d(i,j,k-gamma), x(i-2*alpha), x(i-alpha), x(i), x_u)
               scont_u = interpol_yp(x(i), x(i-alpha), scont3d(i,j,k-gamma), scont3d(i-alpha,j,k-gamma), x_u)
!               int_u = interpol_yp(x(i), x(i-alpha), int3d(i,j,k-gamma), int3d(i-alpha,j,k-gamma), x_u)
               dcoeff=0.d0
!***debug start
!               if(imask_innreg3d(i-alpha,j,k-gamma).eq.1) then
!                  z_u1=z(k-gamma)
!                  x_u1=sign(1.,x_u)*sqrt(1.d0-z_u1**2)
!                  mueff = n_x*x_u1 + n_z*z_u1
!                  if(mueff.ge.0.d0) then
!                     int_u = interpol_yp(x(i), x_u1, int3d(i,j,k-gamma), xic1, x_u)
!                  else
!                     int_u = interpol_yp(x(i), x_u1, int3d(i,j,k-gamma), 0.d0, x_u)
!!                     int_u = interpol_yp(x(i), x_u1, int3d(i,j,k-gamma), xic1, x_u)
!                  endif
!               elseif(imask_innreg3d(i,j,k-gamma).eq.1) then
!                  z_u1=z(k-gamma)
!                  x_u1=sign(1.,x_u)*sqrt(1.d0-z_u1**2)
!                  mueff = n_x*x_u1 + n_z*z_u1
!                  if(mueff.ge.0.d0) then
!                     int_u = interpol_yp(x(i-alpha), x_u1, int3d(i-alpha,j,k-gamma), xic1, x_u)
!                  else
!                     int_u = interpol_yp(x(i-alpha), x_u1, int3d(i-alpha,j,k-gamma), 0.d0, x_u)
!!                     int_u = interpol_yp(x(i-alpha), x_u1, int3d(i-alpha,j,k-gamma), xic1, x_u)
!                  endif
!               endif
!***debug end
            endif
!
!--------------------------downwind point-------------------------------
!
            if(dels_xd.eq.dels_d) then
!downwind point intersects on z-level
               x_d=x(i)+dels_d*n_x
               z_d=z(k)+dels_d*n_z
               opac_d = interpol_ypl(z(k+gamma), z(k), opac3d(i+alpha,j,k+gamma), opac3d(i+alpha,j,k), z_d)
!               scont_d = interpol_typ_quad3b(scont3d(i+alpha,j,k-gamma), scont3d(i+alpha,j,k), scont3d(i+alpha,j,k+gamma), &
!                                             z(k-gamma), z(k), z(k+gamma), z_d)
               scont_d = interpol_yp(z(k+gamma), z(k), scont3d(i+alpha,j,k+gamma), scont3d(i+alpha,j,k), z_d)
            else
!downwind point intersects on x-level
               x_d=x(i)+dels_d*n_x
               z_d=z(k)+dels_d*n_z
               opac_d = interpol_ypl(x(i+alpha), x(i), opac3d(i+alpha,j,k+gamma), opac3d(i,j,k+gamma), x_d)
!               scont_d = interpol_typ_quad3b(scont3d(i-alpha,j,k+gamma), scont3d(i,j,k+gamma), scont3d(i+alpha,j,k+gamma), &
!                                             x(i-alpha), x(i), x(i+alpha), x_d)
               scont_d = interpol_yp(x(i+alpha), x(i), scont3d(i+alpha,j,k+gamma), scont3d(i,j,k+gamma), x_d)
            endif
!
!-----------------------radiative transfer------------------------------
!
            call fsc_cont(int_u, opac_u, opac_p, opac_d, scont_u, scont_p, scont_d, &
               dels_u, dels_d, abs_sc, contr_sc, int_sc, alo_p, alo_u, alo_d)
!            call fsc_cont_lin(int_u, opac_u, opac_p, opac_d, scont_u, scont_p, scont_d, &
!                                dels_u, dels_d, abs_sc, contr_sc, int_sc, alo_u, alo_p)
            int3d(i,j,k) = int_sc
            alocont_nn3d(i,j,k,14) = alo_p! + dcoeff*alo_u
!
!----------------------point lies directly on boundary------------------
!***********special treatment: point lies directly on boundary**********
!
          case(3)
            mueff=n_x*x(i)+n_z*z(k)
            if(mueff.ge.0.d0) then
               int3d(i,j,k)=xic1
               alocont_nn3d(i,j,k,14) = 0.d0
            else
!
!calculate distance to intersection point with upwind x-grid-coordinate
               dels_xu=(x(i)-x(i-alpha))/n_x
!calculate distance to instersection point with upwind z-grid-coordinate
               dels_zu=(z(k)-z(k-gamma))/n_z
!calculate distance to intersection point with downwind x-grid-coordinate
               dels_xd=(x(i+alpha)-x(i))/n_x
!calculate distance to instersection point with downwind z-grid-coordinate
               dels_zd=(z(k+gamma)-z(k))/n_z
!
               dels_u=min(dels_xu,dels_zu)
               dels_d=min(dels_xd,dels_zd)
!
!----------------------------local point--------------------------------
!
               opac_p=opac3d(i,j,k)
               scont_p=scont3d(i,j,k)
!
!----------------------------upwind point-------------------------------
!
               if(dels_xu.eq.dels_u) then
!upwind point intersects on z-level
                  x_u=x(i)-dels_u*n_x
                  z_u=z(k)-dels_u*n_z
!interpolation with bezier spline or with log-log
                  opac_u = interpol_ypl(z(k), z(k-gamma), opac3d(i-alpha,j,k), opac3d(i-alpha,j,k-gamma), z_u)
!                  scont_u = interpol_typ_quad3b(scont3d(i-alpha,j,k-2*gamma), scont3d(i-alpha,j,k-gamma), scont3d(i-alpha,j,k), &
!                                                z(k-2*gamma), z(k-gamma), z(k), z_u)
                  int_u = interpol_typ_quad2b(int3d(i-alpha,j,k-2*gamma), int3d(i-alpha,j,k-gamma), &
                     int3d(i-alpha,j,k), z(k-2*gamma), z(k-gamma), z(k), z_u)
                  scont_u = interpol_yp(z(k), z(k-gamma), scont3d(i-alpha,j,k), scont3d(i-alpha,j,k-gamma), z_u)
!                  int_u = interpol_yp(z(k), z(k-gamma), int3d(i-alpha,j,k), int3d(i-alpha,j,k-gamma), z_u)
               else
!upwind point intersects on x-level
                  x_u=x(i)-dels_u*n_x
                  z_u=z(k)-dels_u*n_z
!interpolation with bezier spline or with log-log
                  opac_u = interpol_ypl(x(i), x(i-alpha), opac3d(i,j,k-gamma), opac3d(i-alpha,j,k-gamma), x_u)
!                  scont_u = interpol_typ_quad3b(scont3d(i-2*alpha,j,k-gamma), scont3d(i-alpha,j,k-gamma), scont3d(i,j,k-gamma), &
!                                                x(i-2*alpha), x(i-alpha), x(i), x_u)
                  int_u = interpol_typ_quad2b(int3d(i-2*alpha,j,k-gamma), int3d(i-alpha,j,k-gamma), &
                     int3d(i,j,k-gamma), x(i-2*alpha), x(i-alpha), x(i), x_u)
                  scont_u = interpol_yp(x(i), x(i-alpha), scont3d(i,j,k-gamma), scont3d(i-alpha,j,k-gamma), x_u)
!                  int_u = interpol_yp(x(i), x(i-alpha), int3d(i,j,k-gamma), int3d(i-alpha,j,k-gamma), x_u)
               endif
!
!--------------------------downwind point-------------------------------
!
               if(dels_xd.eq.dels_d) then
!downwind point intersects on z-level
                  x_d=x(i)+dels_d*n_x
                  z_d=z(k)+dels_d*n_z
                  opac_d = interpol_ypl(z(k+gamma), z(k), opac3d(i+alpha,j,k+gamma), opac3d(i+alpha,j,k), z_d)
!                  scont_d = interpol_typ_quad3b(scont3d(i+alpha,j,k-gamma), scont3d(i+alpha,j,k), scont3d(i+alpha,j,k+gamma), &
!                                             z(k-gamma), z(k), z(k+gamma), z_d)
                  scont_d = interpol_yp(z(k+gamma), z(k), scont3d(i+alpha,j,k+gamma), scont3d(i+alpha,j,k), z_d)
               else
!downwind point intersects on x-level
                  x_d=x(i)+dels_d*n_x
                  z_d=z(k)+dels_d*n_z
                  opac_d = interpol_ypl(x(i+alpha), x(i), opac3d(i+alpha,j,k+gamma), opac3d(i,j,k+gamma), x_d)
!                  scont_d = interpol_typ_quad3b(scont3d(i-alpha,j,k+gamma), scont3d(i,j,k+gamma), scont3d(i+alpha,j,k+gamma), &
!                                                x(i-alpha), x(i), x(i+alpha), x_d)
                  scont_d = interpol_yp(x(i+alpha), x(i), scont3d(i+alpha,j,k+gamma), scont3d(i,j,k+gamma), x_d)
               endif
!
!-----------------------radiative transfer------------------------------
!
               call fsc_cont(int_u, opac_u, opac_p, opac_d, scont_u, scont_p, scont_d, &
                  dels_u, dels_d, abs_sc, contr_sc, int_sc, alo_p, alo_u, alo_d)
!               call fsc_cont_lin(int_u, opac_u, opac_p, opac_d, scont_u, scont_p, scont_d, &
!                                   dels_u, dels_d, abs_sc, contr_sc, int_sc, alo_u, alo_p)
               int3d(i,j,k) = int_sc
               alocont_nn3d(i,j,k,14) = alo_p
            endif
!
          case default
!
         end select
      enddo
   enddo
!
!call cpu_time(te2)

!write(*,'(a10,3es20.8)') 'timing', te1-ts1, te2-ts2, (te2-ts2)/(te1-ts1)
!
end subroutine fsc_cont2d
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine fsc_cont2d_lin(muindx)
!
!-----------------------------------------------------------------------
!------short characteristics for continuum radiative transfer in 2d-----
!---------------calculating intensties for given mu---------------------
!-----------------------------------------------------------------------
!
   use prog_type

   use fund_const
   use dime3d, only: ndxmax, ndzmax, ndymax, x, y, z
   use dime3d, only: int3d, opac3d, scont3d, alocont_nn3d, imask3d, imask_innreg3d
   use angles, only: nodes_mu
   use bcondition, only: xic1
   use mod_debug, only: indx1, indx2
   use mod_interp1d, only: interpol_yp, interpol_ypl, interpol_typ_quad2b, interpol_typ_quad3b
   use mod_interp2d, only: interpol2d_4p_lin
!
   implicit none
!
! ... arguments
   integer(i4b), intent(in) :: muindx
   real(dp) :: n_x, n_z
!
! ... local scalars
   integer(i4b) :: i, j, k
   integer(i4b) :: alpha, gamma
   integer(i4b) :: startx, startz, endx, endz
   integer(i4b) :: iim2, iim1, ii, iip1
   real(dp) ::  mueff
   real(dp) :: opac_u, scont_u, dels_u, dels_xu, dels_zu, int_u, &
      opac_d, scont_d, dels_d, dels_xd, dels_zd, &
      opac_u1, opac_u2, scont_u1, scont_u2, &
      opac_d1, opac_d2, scont_d1, scont_d2, &
      opac_p, scont_p, dels_r
   real(dp) :: x_u, x_u1, x_u2, z_u, z_u1, z_u2, t_u1, t_u2, t_u, &
      x_d, x_d1, x_d2, z_d, z_d1, z_d2, t_d1, t_d2, t_d
   real(dp) :: abs_sc, int_sc, contr_sc, alo_p, alo_u, alo_d
   real(dp) :: acoeff, bcoeff, ccoeff, dcoeff
   real(dp) :: ts1, te1, ts2, te2
!
! ... for debugging
!
! ... local arrays
!
! ... local functions
   real(dp) :: dist_sphere
!
! ... local characters
   character(len=50) :: enter
!
! ... local logicals
!
   n_x=sqrt(1.d0 - nodes_mu(muindx)*nodes_mu(muindx))
   n_z=nodes_mu(muindx)
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
!         if n_z >= 0                 if n_z < 0
!                startz = 2                  startz = ndzmax-1
!                endz = ndzmax-1             endz = 2
!                gamma = 1                   gamma = -1
!
!could make an own array for this!!!
   if(n_x.gt.0.d0) then
      startx = 3
      endx = ndxmax-1
      alpha=  1
   elseif(n_x.lt.0.d0) then
      startx = ndxmax-2
      endx = 2
      alpha=-1
   else
      stop 'error in fsc_cont2d: n_x = 0 not allowed'
   endif
!
   if(n_z.gt.0.d0) then
      startz = 3
      endz = ndzmax-1
      gamma=  1
   elseif(n_z.lt.0.d0) then
      startz = ndzmax-2
      endz = 2
      gamma=-1
   else
      stop 'error in fsc_cont2d: n_z = 0 not allowed'
   endif
!
!-----------------------reset the intensities---------------------------
!
!call cpu_time(ts1)
   call set_boundary2d(muindx)
!call cpu_time(te1)
!
!-----------------------------------------------------------------------
!
   alocont_nn3d=0.d0
!
!-----------------------------------------------------------------------
!
   j=ndymax/2+1

!call cpu_time(ts2)
!
   do k=startz, endz, gamma
      do i=startx, endx, alpha

         select case(imask3d(i,j,k))

!
!*************************standard rt procedure*************************
!
          case(1)
!
!calculate distance to intersection point with upwind x-grid-coordinate
            dels_xu=(x(i)-x(i-alpha))/n_x
!calculate distance to instersection point with upwind z-grid-coordinate
            dels_zu=(z(k)-z(k-gamma))/n_z
            dels_u=min(dels_xu,dels_zu)
!calculate distance to intersection point with downwind x-grid-coordinate
            dels_xd=(x(i+alpha)-x(i))/n_x
!calculate distance to instersection point with downwind z-grid-coordinate
            dels_zd=(z(k+gamma)-z(k))/n_z
            dels_d=min(dels_xd,dels_zd)
!
!----------------------------local point--------------------------------
!
            opac_p=opac3d(i,j,k)
            scont_p=scont3d(i,j,k)
!
!----------------------------upwind point-------------------------------
!
            if(dels_xu.eq.dels_u) then
!upwind point intersects on z-level
               x_u=x(i)-dels_u*n_x
               z_u=z(k)-dels_u*n_z
!interpolation with bezier spline or with log-log
               opac_u = interpol_ypl(z(k), z(k-gamma), opac3d(i-alpha,j,k), opac3d(i-alpha,j,k-gamma), z_u)
!               scont_u = interpol_typ_quad3b(scont3d(i-alpha,j,k-2*gamma), scont3d(i-alpha,j,k-gamma), scont3d(i-alpha,j,k), &
!                                             z(k-2*gamma), z(k-gamma), z(k), z_u)
               int_u = interpol_typ_quad2b(int3d(i-alpha,j,k-2*gamma), int3d(i-alpha,j,k-gamma), &
                  int3d(i-alpha,j,k), z(k-2*gamma), z(k-gamma), z(k), z_u)
               scont_u = interpol_yp(z(k), z(k-gamma), scont3d(i-alpha,j,k), scont3d(i-alpha,j,k-gamma), z_u)
!               int_u = interpol_yp(z(k), z(k-gamma), int3d(i-alpha,j,k), int3d(i-alpha,j,k-gamma), z_u)
            else
!upwind point intersects on x-level
               x_u=x(i)-dels_u*n_x
               z_u=z(k)-dels_u*n_z
!interpolation with bezier spline or with log-log
               opac_u = interpol_ypl(x(i), x(i-alpha), opac3d(i,j,k-gamma), opac3d(i-alpha,j,k-gamma), x_u)
!               scont_u = interpol_typ_quad3b(scont3d(i-2*alpha,j,k-gamma), scont3d(i-alpha,j,k-gamma), scont3d(i,j,k-gamma), &
!                                             x(i-2*alpha), x(i-alpha), x(i), x_u)
               int_u = interpol_typ_quad2b(int3d(i-2*alpha,j,k-gamma), int3d(i-alpha,j,k-gamma), &
                  int3d(i,j,k-gamma), x(i-2*alpha), x(i-alpha), x(i), x_u)
               scont_u = interpol_yp(x(i), x(i-alpha), scont3d(i,j,k-gamma), scont3d(i-alpha,j,k-gamma), x_u)
!               int_u = interpol_yp(x(i), x(i-alpha), int3d(i,j,k-gamma), int3d(i-alpha,j,k-gamma), x_u)
            endif
!
!--------------------------downwind point-------------------------------
!
            if(dels_xd.eq.dels_d) then
!downwind point intersects on z-level
               x_d=x(i)+dels_d*n_x
               z_d=z(k)+dels_d*n_z
               opac_d = interpol_ypl(z(k+gamma), z(k), opac3d(i+alpha,j,k+gamma), opac3d(i+alpha,j,k), z_d)
!               scont_d = interpol_typ_quad3b(scont3d(i+alpha,j,k-gamma), scont3d(i+alpha,j,k), scont3d(i+alpha,j,k+gamma), &
!                                             z(k-gamma), z(k), z(k+gamma), z_d)
               scont_d = interpol_yp(z(k+gamma), z(k), scont3d(i+alpha,j,k+gamma), scont3d(i+alpha,j,k), z_d)
            else
!downwind point intersects on x-level
               dels_d=dels_zd
               x_d=x(i)+dels_d*n_x
               z_d=z(k)+dels_d*n_z
               opac_d = interpol_ypl(x(i+alpha), x(i), opac3d(i+alpha,j,k+gamma), opac3d(i,j,k+gamma), x_d)
!               scont_d = interpol_typ_quad3b(scont3d(i-alpha,j,k+gamma), scont3d(i,j,k+gamma), scont3d(i+alpha,j,k+gamma), &
!                                             x(i-alpha), x(i), x(i+alpha), x_d)
               scont_d = interpol_ypl(x(i+alpha), x(i), scont3d(i+alpha,j,k+gamma), scont3d(i,j,k+gamma), x_d)
            endif
!
!-----------------------radiative transfer------------------------------
!
!            call fsc_cont(int_u, opac_u, opac_p, opac_d, scont_u, scont_p, scont_d, &
!                            dels_u, dels_d, abs_sc, contr_sc, int_sc, alo_p, alo_u, alo_d)
            call fsc_cont_lin(int_u, opac_u, opac_p, opac_d, scont_u, scont_p, scont_d, &
               dels_u, dels_d, abs_sc, contr_sc, int_sc, alo_u, alo_p)
            int3d(i,j,k) = int_sc
            alocont_nn3d(i,j,k,14) = alo_p
!
!***********special treatment: point is in vicinity of boundary*********
!
          case(2)
!
!calculate distance to intersection point with upwind x-grid-coordinate
            dels_xu=(x(i)-x(i-alpha))/n_x
!calculate distance to instersection point with upwind z-grid-coordinate
            dels_zu=(z(k)-z(k-gamma))/n_z
!calculate distance to intersection point with downwind x-grid-coordinate
            dels_xd=(x(i+alpha)-x(i))/n_x
!calculate distance to instersection point with downwind z-grid-coordinate
            dels_zd=(z(k+gamma)-z(k))/n_z
!calculate distance to the boundary (here: sphere)
            dels_r=dist_sphere(n_x,0.d0,n_z,x(i),0.d0,z(k))
!
            dels_u=min(dels_xu,dels_zu,dels_r)
            dels_d=min(dels_xd,dels_zd)
!
!----------------------------local point--------------------------------
!
            opac_p=opac3d(i,j,k)
            scont_p=scont3d(i,j,k)
!
!----------------------------upwind point-------------------------------
!
            if(dels_r.eq.dels_u) then
!
!upwind point intersects with sphere
!               indx1=indx1+1
               x_u=x(i)-dels_r*n_x
               z_u=z(k)-dels_r*n_z
               int_u=xic1
!               opac_u=opac3d(ndxmax/2+1,ndymax/2+1,64)
               scont_u=xic1
!               scont_u=interpol2d_4p_lin(scont3d(i-alpha,j,k-gamma), scont3d(i,j,k-gamma), &
!                                         scont3d(i-alpha,j,k), scont3d(i,j,k), &
!                                         x(i-alpha), x(i), z(k-gamma), z(k), x_u, z_u)
!               call coeff2d_4p_lin(x(i-alpha), x(i), z(k-gamma), z(k), x_u, z_u, &
!                                   acoeff, bcoeff, ccoeff, dcoeff)
!               scont_u = acoeff*scont3d(i-alpha,j,k-gamma) + bcoeff*scont3d(i,j,k-gamma) + &
!                         ccoeff*scont3d(i-alpha,j,k) + dcoeff*scont3d(i,j,k)
               opac_u = interpol2d_4p_lin(opac3d(i-alpha,j,k-gamma), opac3d(i,j,k-gamma), &
                  opac3d(i-alpha,j,k), opac3d(i,j,k), &
                  x(i-alpha), x(i), z(k-gamma), z(k), x_u, z_u)
            elseif(dels_xu.eq.dels_u) then
!               indx2=indx2+1
!upwind point intersects on z-level
               x_u=x(i)-dels_u*n_x
               z_u=z(k)-dels_u*n_z
!interpolation with bezier spline or with log-log
               opac_u = interpol_ypl(z(k), z(k-gamma), opac3d(i-alpha,j,k), opac3d(i-alpha,j,k-gamma), z_u)
!               scont_u = interpol_typ_quad3b(scont3d(i-alpha,j,k-2*gamma), scont3d(i-alpha,j,k-gamma), scont3d(i-alpha,j,k), &
!                                             z(k-2*gamma), z(k-gamma), z(k), z_u)
               int_u = interpol_typ_quad2b(int3d(i-alpha,j,k-2*gamma), int3d(i-alpha,j,k-gamma), &
                  int3d(i-alpha,j,k), z(k-2*gamma), z(k-gamma), z(k), z_u)
               scont_u = interpol_yp(z(k), z(k-gamma), scont3d(i-alpha,j,k), scont3d(i-alpha,j,k-gamma), z_u)
!               int_u = interpol_yp(z(k), z(k-gamma), int3d(i-alpha,j,k), int3d(i-alpha,j,k-gamma), z_u)
               dcoeff=0.d0
!***debug start
!               if(imask_innreg3d(i-alpha,j,k-gamma).eq.1) then
!                  x_u1=x(i-alpha)
!                  z_u1=sign(1.,z_u)*sqrt(1.d0-x_u1**2)
!                  mueff = n_x*x_u1 + n_z*z_u1
!                  if(mueff.ge.0.d0) then
!                     int_u = interpol_yp(z(k), z_u1, int3d(i-alpha,j,k), xic1, z_u)
!                  else
!                     int_u = interpol_yp(z(k), z_u1, int3d(i-alpha,j,k), 0.d0, z_u)
!!                     int_u = interpol_yp(z(k), z_u1, int3d(i-alpha,j,k), xic1, z_u)
!                  endif
!               elseif(imask_innreg3d(i-alpha,j,k).eq.1) then
!                  x_u1=x(i-alpha)
!                  z_u1=sign(1.,z_u)*sqrt(1.d0-x_u1**2)
!                  mueff = n_x*x_u1 + n_z*z_u1
!                  if(mueff.ge.0.d0) then
!                     int_u = interpol_yp(z(k-gamma), z_u1, int3d(i-alpha,j,k-gamma), xic1, z_u)
!                  else
!                     int_u = interpol_yp(z(k-gamma), z_u1, int3d(i-alpha,j,k-gamma), 0.d0, z_u)
!!                     int_u = interpol_yp(z(k-gamma), z_u1, int3d(i-alpha,j,k-gamma), xic1, z_u)
!                  endif
!               endif
!***debug end
            else
!              indx2=indx2+1
!upwind point intersects on x-level
               x_u=x(i)-dels_u*n_x
               z_u=z(k)-dels_u*n_z
!interpolation with bezier spline or with log-log
               opac_u = interpol_ypl(x(i), x(i-alpha), opac3d(i,j,k-gamma), opac3d(i-alpha,j,k-gamma), x_u)
!               scont_u = interpol_typ_quad3b(scont3d(i-2*alpha,j,k-gamma), scont3d(i-alpha,j,k-gamma), scont3d(i,j,k-gamma), &
!                                             x(i-2*alpha), x(i-alpha), x(i), x_u)
               int_u = interpol_typ_quad2b(int3d(i-2*alpha,j,k-gamma), int3d(i-alpha,j,k-gamma), &
                  int3d(i,j,k-gamma), x(i-2*alpha), x(i-alpha), x(i), x_u)
               scont_u = interpol_yp(x(i), x(i-alpha), scont3d(i,j,k-gamma), scont3d(i-alpha,j,k-gamma), x_u)
!               int_u = interpol_yp(x(i), x(i-alpha), int3d(i,j,k-gamma), int3d(i-alpha,j,k-gamma), x_u)
               dcoeff=0.d0
!***debug start
!               if(imask_innreg3d(i-alpha,j,k-gamma).eq.1) then
!                  z_u1=z(k-gamma)
!                  x_u1=sign(1.,x_u)*sqrt(1.d0-z_u1**2)
!                  mueff = n_x*x_u1 + n_z*z_u1
!                  if(mueff.ge.0.d0) then
!                     int_u = interpol_yp(x(i), x_u1, int3d(i,j,k-gamma), xic1, x_u)
!                  else
!                     int_u = interpol_yp(x(i), x_u1, int3d(i,j,k-gamma), 0.d0, x_u)
!!                     int_u = interpol_yp(x(i), x_u1, int3d(i,j,k-gamma), xic1, x_u)
!                  endif
!               elseif(imask_innreg3d(i,j,k-gamma).eq.1) then
!                  z_u1=z(k-gamma)
!                  x_u1=sign(1.,x_u)*sqrt(1.d0-z_u1**2)
!                  mueff = n_x*x_u1 + n_z*z_u1
!                  if(mueff.ge.0.d0) then
!                     int_u = interpol_yp(x(i-alpha), x_u1, int3d(i-alpha,j,k-gamma), xic1, x_u)
!                  else
!                     int_u = interpol_yp(x(i-alpha), x_u1, int3d(i-alpha,j,k-gamma), 0.d0, x_u)
!!                     int_u = interpol_yp(x(i-alpha), x_u1, int3d(i-alpha,j,k-gamma), xic1, x_u)
!                  endif
!               endif
!***debug end
            endif
!
!--------------------------downwind point-------------------------------
!
            if(dels_xd.eq.dels_d) then
!downwind point intersects on z-level
               x_d=x(i)+dels_d*n_x
               z_d=z(k)+dels_d*n_z
               opac_d = interpol_ypl(z(k+gamma), z(k), opac3d(i+alpha,j,k+gamma), opac3d(i+alpha,j,k), z_d)
!               scont_d = interpol_typ_quad3b(scont3d(i+alpha,j,k-gamma), scont3d(i+alpha,j,k), scont3d(i+alpha,j,k+gamma), &
!                                             z(k-gamma), z(k), z(k+gamma), z_d)
               scont_d = interpol_yp(z(k+gamma), z(k), scont3d(i+alpha,j,k+gamma), scont3d(i+alpha,j,k), z_d)
            else
!downwind point intersects on x-level
               x_d=x(i)+dels_d*n_x
               z_d=z(k)+dels_d*n_z
               opac_d = interpol_ypl(x(i+alpha), x(i), opac3d(i+alpha,j,k+gamma), opac3d(i,j,k+gamma), x_d)
!               scont_d = interpol_typ_quad3b(scont3d(i-alpha,j,k+gamma), scont3d(i,j,k+gamma), scont3d(i+alpha,j,k+gamma), &
!                                             x(i-alpha), x(i), x(i+alpha), x_d)
               scont_d = interpol_yp(x(i+alpha), x(i), scont3d(i+alpha,j,k+gamma), scont3d(i,j,k+gamma), x_d)
            endif
!
!-----------------------radiative transfer------------------------------
!
!            call fsc_cont(int_u, opac_u, opac_p, opac_d, scont_u, scont_p, scont_d, &
!                            dels_u, dels_d, abs_sc, contr_sc, int_sc, alo_p, alo_u, alo_d)
            call fsc_cont_lin(int_u, opac_u, opac_p, opac_d, scont_u, scont_p, scont_d, &
               dels_u, dels_d, abs_sc, contr_sc, int_sc, alo_u, alo_p)
            int3d(i,j,k) = int_sc
            alocont_nn3d(i,j,k,14) = alo_p! + dcoeff*alo_u
!
!----------------------point lies directly on boundary------------------
!***********special treatment: point lies directly on boundary**********
!
          case(3)
            mueff=n_x*x(i)+n_z*z(k)
            if(mueff.ge.0.d0) then
               int3d(i,j,k)=xic1
               alocont_nn3d(i,j,k,14) = 0.d0
            else
!
!calculate distance to intersection point with upwind x-grid-coordinate
               dels_xu=(x(i)-x(i-alpha))/n_x
!calculate distance to instersection point with upwind z-grid-coordinate
               dels_zu=(z(k)-z(k-gamma))/n_z
!calculate distance to intersection point with downwind x-grid-coordinate
               dels_xd=(x(i+alpha)-x(i))/n_x
!calculate distance to instersection point with downwind z-grid-coordinate
               dels_zd=(z(k+gamma)-z(k))/n_z
!
               dels_u=min(dels_xu,dels_zu)
               dels_d=min(dels_xd,dels_zd)
!
!----------------------------local point--------------------------------
!
               opac_p=opac3d(i,j,k)
               scont_p=scont3d(i,j,k)
!
!----------------------------upwind point-------------------------------
!
               if(dels_xu.eq.dels_u) then
!upwind point intersects on z-level
                  x_u=x(i)-dels_u*n_x
                  z_u=z(k)-dels_u*n_z
!interpolation with bezier spline or with log-log
                  opac_u = interpol_ypl(z(k), z(k-gamma), opac3d(i-alpha,j,k), opac3d(i-alpha,j,k-gamma), z_u)
!                  scont_u = interpol_typ_quad3b(scont3d(i-alpha,j,k-2*gamma), scont3d(i-alpha,j,k-gamma), scont3d(i-alpha,j,k), &
!                                                z(k-2*gamma), z(k-gamma), z(k), z_u)
                  int_u = interpol_typ_quad2b(int3d(i-alpha,j,k-2*gamma), int3d(i-alpha,j,k-gamma), &
                     int3d(i-alpha,j,k), z(k-2*gamma), z(k-gamma), z(k), z_u)
                  scont_u = interpol_yp(z(k), z(k-gamma), scont3d(i-alpha,j,k), scont3d(i-alpha,j,k-gamma), z_u)
!                  int_u = interpol_yp(z(k), z(k-gamma), int3d(i-alpha,j,k), int3d(i-alpha,j,k-gamma), z_u)
               else
!upwind point intersects on x-level
                  x_u=x(i)-dels_u*n_x
                  z_u=z(k)-dels_u*n_z
!interpolation with bezier spline or with log-log
                  opac_u = interpol_ypl(x(i), x(i-alpha), opac3d(i,j,k-gamma), opac3d(i-alpha,j,k-gamma), x_u)
!                  scont_u = interpol_typ_quad3b(scont3d(i-2*alpha,j,k-gamma), scont3d(i-alpha,j,k-gamma), scont3d(i,j,k-gamma), &
!                                                x(i-2*alpha), x(i-alpha), x(i), x_u)
                  int_u = interpol_typ_quad2b(int3d(i-2*alpha,j,k-gamma), int3d(i-alpha,j,k-gamma), &
                     int3d(i,j,k-gamma), x(i-2*alpha), x(i-alpha), x(i), x_u)
                  scont_u = interpol_yp(x(i), x(i-alpha), scont3d(i,j,k-gamma), scont3d(i-alpha,j,k-gamma), x_u)
!                  int_u = interpol_yp(x(i), x(i-alpha), int3d(i,j,k-gamma), int3d(i-alpha,j,k-gamma), x_u)
               endif
!
!--------------------------downwind point-------------------------------
!
               if(dels_xd.eq.dels_d) then
!downwind point intersects on z-level
                  x_d=x(i)+dels_d*n_x
                  z_d=z(k)+dels_d*n_z
                  opac_d = interpol_ypl(z(k+gamma), z(k), opac3d(i+alpha,j,k+gamma), opac3d(i+alpha,j,k), z_d)
!                  scont_d = interpol_typ_quad3b(scont3d(i+alpha,j,k-gamma), scont3d(i+alpha,j,k), scont3d(i+alpha,j,k+gamma), &
!                                             z(k-gamma), z(k), z(k+gamma), z_d)
                  scont_d = interpol_yp(z(k+gamma), z(k), scont3d(i+alpha,j,k+gamma), scont3d(i+alpha,j,k), z_d)
               else
!downwind point intersects on x-level
                  x_d=x(i)+dels_d*n_x
                  z_d=z(k)+dels_d*n_z
                  opac_d = interpol_ypl(x(i+alpha), x(i), opac3d(i+alpha,j,k+gamma), opac3d(i,j,k+gamma), x_d)
!                  scont_d = interpol_typ_quad3b(scont3d(i-alpha,j,k+gamma), scont3d(i,j,k+gamma), scont3d(i+alpha,j,k+gamma), &
!                                                x(i-alpha), x(i), x(i+alpha), x_d)
                  scont_d = interpol_yp(x(i+alpha), x(i), scont3d(i+alpha,j,k+gamma), scont3d(i,j,k+gamma), x_d)
               endif
!
!-----------------------radiative transfer------------------------------
!
!               call fsc_cont(int_u, opac_u, opac_p, opac_d, scont_u, scont_p, scont_d, &
!                               dels_u, dels_d, abs_sc, contr_sc, int_sc, alo_p, alo_u, alo_d)
               call fsc_cont_lin(int_u, opac_u, opac_p, opac_d, scont_u, scont_p, scont_d, &
                  dels_u, dels_d, abs_sc, contr_sc, int_sc, alo_u, alo_p)
               int3d(i,j,k) = int_sc
               alocont_nn3d(i,j,k,14) = alo_p
            endif
!
          case default
!
         end select
      enddo
   enddo
!
!call cpu_time(te2)

!write(*,'(a10,3es20.8)') 'timing', te1-ts1, te2-ts2, (te2-ts2)/(te1-ts1)
!
end subroutine fsc_cont2d_lin
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine set_boundary2d(muindx)
!
   use prog_type
   use dime3d, only: ndxmax, ndymax, ndzmax, x, y, z, imask3d, imask_innreg3d, int3d
   use angles, only: nodes_mu
   use bcondition, only: xic1, int_inner, indx_xinner, indx_yinner, indx_zinner, n_inner
!
   implicit none
!
! ... arguments
   integer(i4b), intent(in) :: muindx
!
! ... local scalars
   integer(i4b) :: i, j, k
   integer(i4b) :: phiindx
   real(dp) :: mueff
   real(dp) :: n_x, n_z
!
! ... local arrays
!real(dp), dimension(:,:,:), allocatable :: int3d_test
!
! ... local characters
   character(len=32) :: fname
!
!
   int3d=0.d0
!
!-------------------method 1: calculate mueff directly------------------
!
   n_x=sqrt(1.d0 - nodes_mu(muindx)**2)
   n_z=nodes_mu(muindx)
!
   j=ndymax/2+1
!
   do k=3, ndzmax-2
      do i=3, ndxmax-2
         if(imask3d(i,j,k).eq.4) then
!inside the star, set intensities correctly
            mueff=(n_x*x(i)+n_z*z(k))/sqrt(x(i)**2+z(k)**2)
            if(mueff.ge.0.d0) then
               int3d(i,j,k) = xic1
            else
               int3d(i,j,k) = 0.d0
            endif
         endif
      enddo
   enddo
!
!--------------------method2: read in everything------------------------
!
!allocate(int3d_test(ndxmax,ndymax,ndzmax))
!int3d_test=0.d0
!
!write(*,*) '---------------reading in inner boundary condition: intensity------------------'
!
!phiindx=1
!write(fname,'(a19,i4.4,a5,i4.4)') 'data_bcondition/imu', muindx, '_iphi', phiindx
!open(1, file=fname, form='unformatted')
!   read(1) int_inner
!   read(1) indx_xinner
!   read(1) indx_yinner
!   read(1) indx_zinner
!close(1)
!!
!do i=1, n_inner
!   int3d_test(indx_xinner(i),indx_yinner(i),indx_zinner(i)) = int_inner(i)
!enddo!
!
!j=ndymax/2+1
!do i=1, ndxmax
!   do k=1, ndzmax
!      if(int3d(i,j,k).ne.int3d_test(i,j,k)) then
!         write(*,*) 'error: intensities different'
!         write(*,*) i, j, k, x(i), y(j), z(k)
!         write(*,*) int3d(i,j,k), int3d_test(i,j,k)
!         write(*,*) muindx
!         stop
!      endif
!   enddo
!enddo
!stop
!
!
!
end subroutine set_boundary2d
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine fsc_line2d(muindx,xobsindx)
!
!-----------------------------------------------------------------------
!-------short characteristics for line radiative transfer in 2d---------
!-----------calculating intensities for given mu and xobs---------------
!-----------------------------------------------------------------------
!
   use prog_type

   use fund_const
   use dime3d, only: ndxmax, ndzmax, ndymax, x, y, z
   use dime3d, only: int3d, mintbar3d, opalbar3d, sline3d, velx3d, vely3d, velz3d, vth3d, &
      normalization3d, imask3d, aloline_nn3d, aloline_on_nn3d
   use dimecr, only: n1d_cr, r1d_cr, sline1d_cr
   use angles, only: nodes_mu, weight_mu, dim_mu
   use freq, only: nodes_xobs, weight_xobs, nxobs, deltax, xcmf_min, xcmf_max
   use params_input, only: vth_fiducial
   use bcondition, only: xic1
   use mod_interp1d, only: interpol_yp, interpol_ypl, interpol_typ_quad2b, interpol_typ_quad3, interpol_typ_quad3b
   use mod_interp2d, only: coeff2d_4p_lin
!
   implicit none
!
! ... arguments
   integer(i4b), intent(in) :: muindx, xobsindx
!
! ... local scalars
   integer(i4b) :: i, j, k
   integer(i4b) :: alpha, gamma
   integer(i4b) :: startx, startz, endx, endz
   integer(i4b) :: iim2, iim1, iip1, kkm2, kkm1, kkp1
   real(dp) :: nn_x, nn_y, nn_z, xobs, mueff, wmu, wxobs, wall, phinorm
   real(dp) :: opalbar_u, sline_u, dels_u, dels_xu, dels_zu, int_u, &
      opalbar_d, sline_d, dels_d, dels_xd, dels_zd, &
      opalbar_p, sline_p, dels_r, &
      vel_u, vel_p, vel_d, velx_u, velx_p, velx_d, &
      vely_u, vely_p, vely_d, velz_u, velz_p, velz_d, &
      vth_u, vth_p, vth_d, &
      opalbar_u2, opac_p, opac_u2, velx_u2, vely_u2, velz_u2, opac_u
   real(dp) :: x_u, z_u, x_d, z_d
   real(dp) :: acoeff, bcoeff, ccoeff, dcoeff
   real(dp) :: abs_sc, contr_sc, int_sc, alo_u, alo_p, alo_d
!
! ... local arrays
!
! ... local functions
   real(dp) :: dist_sphere
!
! ... local characters
   character(len=50) :: enter
!
! ... local logicals
!
!
!stop 'there might be errors due to new version of function dist_sphere: check this!!!'
!
   nn_x=sqrt(1.d0 - nodes_mu(muindx)*nodes_mu(muindx))
   nn_y=0.d0
   nn_z=nodes_mu(muindx)
!
   xobs=nodes_xobs(xobsindx)
!
   wmu=weight_mu(muindx)
   wxobs=weight_xobs(xobsindx)
   wall=wmu*wxobs
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
!         if n_z >= 0                 if n_z < 0
!                startz = 2                  startz = ndzmax-1
!                endz = ndzmax-1             endz = 2
!                gamma = 1                   gamma = -1
!
!could make an own array for this!!!
   if(nn_x.gt.0.d0) then
      startx = 3
      endx = ndxmax-2
      alpha=  1
   elseif(nn_x.lt.0.d0) then
      startx = ndxmax-2
      endx = 3
      alpha=-1
   else
      stop 'error in fsc_line2d: n_x = 0 not allowed'
   endif
!
   if(nn_z.gt.0.d0) then
      startz = 2
      endz = ndzmax-2
      gamma=  1
   elseif(nn_z.lt.0.d0) then
      startz = ndzmax-2
      endz = 3
      gamma=-1
   else
      stop 'error in fsc_line2d: n_z = 0 not allowed'
   endif
!
!--------------------reset the intensities and alo----------------------
!
   call set_boundary2d(muindx)
!
   aloline_on_nn3d=0.d0
!
!-----------------------------------------------------------------------
!
   j=ndymax/2+1
!
!***debug
!open(1, file='TRASH/interp2.dat', form='formatted', position='append')
!
   do k=startz, endz, gamma
      do i=startx, endx, alpha
!
         select case(imask3d(i,j,k))
!
!************************standard radiative transfer********************
!
          case(1)
!
            iim1=i-alpha
            kkm1=k-gamma
            iip1=i+alpha
            kkp1=k+gamma
!
!calculate distance of upwind point (at previous x, z level)
            dels_xu=(x(i)-x(iim1))/nn_x
            dels_zu=(z(k)-z(kkm1))/nn_z
            dels_u=min(dels_xu,dels_zu)
!calculate distance of downwind point (at next x, z level)
            dels_xd=(x(iip1)-x(i))/nn_x
            dels_zd=(z(kkp1)-z(k))/nn_z
            dels_d=min(dels_xd,dels_zd)
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
            if(dels_xu.eq.dels_u) then
!upwind point intersects on x-level at i-alpha
               x_u=x(iim1)
               z_u=z(k) - dels_u*nn_z
!
               opalbar_u = interpol_yp(z(kkm1), z(k), opalbar3d(iim1,j,kkm1), opalbar3d(iim1,j,k), z_u)
               sline_u =   interpol_yp(z(kkm1), z(k),   sline3d(iim1,j,kkm1),   sline3d(iim1,j,k), z_u)
               velx_u =    interpol_yp(z(kkm1), z(k),    velx3d(iim1,j,kkm1),    velx3d(iim1,j,k), z_u)
               vely_u =    interpol_yp(z(kkm1), z(k),    vely3d(iim1,j,kkm1),    vely3d(iim1,j,k), z_u)
               velz_u =    interpol_yp(z(kkm1), z(k),    velz3d(iim1,j,kkm1),    velz3d(iim1,j,k), z_u)
               vth_u =     interpol_yp(z(kkm1), z(k),     vth3d(iim1,j,kkm1),     vth3d(iim1,j,k), z_u)
               int_u =     interpol_yp(z(kkm1), z(k),     int3d(iim1,j,kkm1),     int3d(iim1,j,k), z_u)

            elseif(dels_zu.eq.dels_u) then
!upwind point intersects on z-level at k-gamma
               x_u=x(i) - dels_u*nn_x
               z_u=z(kkm1)
!
               opalbar_u = interpol_yp(x(iim1), x(i), opalbar3d(iim1,j,kkm1), opalbar3d(i,j,kkm1), x_u)
               sline_u =   interpol_yp(x(iim1), x(i),   sline3d(iim1,j,kkm1),   sline3d(i,j,kkm1), x_u)
               velx_u =    interpol_yp(x(iim1), x(i),    velx3d(iim1,j,kkm1),    velx3d(i,j,kkm1), x_u)
               vely_u =    interpol_yp(x(iim1), x(i),    vely3d(iim1,j,kkm1),    vely3d(i,j,kkm1), x_u)
               velz_u =    interpol_yp(x(iim1), x(i),    velz3d(iim1,j,kkm1),    velz3d(i,j,kkm1), x_u)
               vth_u =     interpol_yp(x(iim1), x(i),     vth3d(iim1,j,kkm1),     vth3d(i,j,kkm1), x_u)
               int_u =     interpol_yp(x(iim1), x(i),     int3d(iim1,j,kkm1),     int3d(i,j,kkm1), x_u)
!
            else
               write(*,'(3es20.8,2l4)') dels_u, dels_xu, dels_zu, dels_u.eq.dels_xu, dels_u.eq.dels_zu
               stop 'error in fsc_line2d: invalid dels_u'
            endif
!
!----------------------------downwind point-----------------------------
!
!            if(dels_xd.eq.dels_d) then
!!upwind point intersects on x-level at i+alpha
!               x_d=x(iip1)
!               z_d=z(k) + dels_d*nn_z
!
!               opalbar_d = interpol_yp(z(k), z(kkp1), opalbar3d(iip1,j,k), opalbar3d(iip1,j,kkp1), z_d)
!               sline_d =   interpol_yp(z(k), z(kkp1),   sline3d(iip1,j,k),   sline3d(iip1,j,kkp1), z_d)
!               velx_d =    interpol_yp(z(k), z(kkp1),    velx3d(iip1,j,k),    velx3d(iip1,j,kkp1), z_d)
!               vely_d =    interpol_yp(z(k), z(kkp1),    vely3d(iip1,j,k),    vely3d(iip1,j,kkp1), z_d)
!               velz_d =    interpol_yp(z(k), z(kkp1),    velz3d(iip1,j,k),    velz3d(iip1,j,kkp1), z_d)
!               vth_d =     interpol_yp(z(k), z(kkp1),     vth3d(iip1,j,k),     vth3d(iip1,j,kkp1), z_d)
!
!            elseif(dels_zd.eq.dels_d) then
!!upwind point intersects on z-level at k+gamma
!               x_d=x(i) + dels_d*nn_x
!               z_d=z(kkp1)
!!
!               opalbar_d = interpol_yp(x(i), x(iip1), opalbar3d(i,j,kkp1), opalbar3d(iip1,j,kkp1), x_d)
!               sline_d =   interpol_yp(x(i), x(iip1),   sline3d(i,j,kkp1),   sline3d(iip1,j,kkp1), x_d)
!               velx_d =    interpol_yp(x(i), x(iip1),    velx3d(i,j,kkp1),    velx3d(iip1,j,kkp1), x_d)
!               vely_d =    interpol_yp(x(i), x(iip1),    vely3d(i,j,kkp1),    vely3d(iip1,j,kkp1), x_d)
!               velz_d =    interpol_yp(x(i), x(iip1),    velz3d(i,j,kkp1),    velz3d(iip1,j,kkp1), x_d)
!               vth_d =     interpol_yp(x(i), x(iip1),     vth3d(i,j,kkp1),     vth3d(iip1,j,kkp1), x_d)
!!
!            else
!               write(*,'(3es20.8,2l4)') dels_d, dels_xd, dels_zd, dels_d.eq.dels_xd, dels_d.eq.dels_zd
!               stop 'error in fsc_line2d: invalid dels_u'
!            endif
!
!-----------------------radiative transfer------------------------------
!
!***debug start
            call model_debug(x_u, 0.d0, z_u, velx_u2, vely_u2, velz_u2, opalbar_u2, opac_u)
!            write(*,'(8es20.8)') opalbar_p2/opalbar_p
!***debug end
            vel_u = nn_x*velx_u + nn_y*vely_u + nn_z*velz_u
!            vel_d = nn_x*velx_d + nn_y*vely_d + nn_z*velz_d
            call fsc_line_lin(xobs, vth_fiducial, deltax, xcmf_min, xcmf_max, int_u, opalbar_u, opalbar_p, &
               sline_u, sline_p, vel_u, vel_p, vth_u, vth_p, &
               dels_u, abs_sc, contr_sc, int_sc, alo_u, alo_p)
!            call fsc_line(xobs, vth_fiducial, deltax, xcmf_min, xcmf_max, int_u, &
!                            opalbar_u, opalbar_p, opalbar_d, &
!                            sline_u, sline_p, sline_d, &
!                            vel_u, vel_p, vel_d, &
!                            vth_u, vth_p, vth_d, &
!                            dels_u, dels_d, &
!                            abs_sc, contr_sc, int_sc, alo_u, alo_p, alo_d)
            int3d(i,j,k) = int_sc
            aloline_on_nn3d(i,j,k,14) = alo_p
!
!perform integration
            call calc_phinorm(vel_p, vth_p, vth_fiducial, xobs, phinorm)
            mintbar3d(i,j,k) = mintbar3d(i,j,k) + int3d(i,j,k)*wall*phinorm
            normalization3d(i,j,k) = normalization3d(i,j,k) + wall*phinorm
            aloline_nn3d(i,j,k,14) = aloline_nn3d(i,j,k,14) + aloline_on_nn3d(i,j,k,14)*wall*phinorm
!
!***********special treatment: point is in vicinity of boundary*********
!
          case(2)
!
            iim1=i-alpha
            kkm1=k-gamma
            iip1=i+alpha
            kkp1=k+gamma
!
!calculate distance of upwind point (at previous x, z level)
            dels_xu=(x(i)-x(iim1))/nn_x
            dels_zu=(z(k)-z(kkm1))/nn_z
!calculate distance of downwind point (at next x, z level)
            dels_xd=(x(iip1)-x(i))/nn_x
            dels_zd=(z(kkp1)-z(k))/nn_z
!calculate distance to the boundary (here: sphere)
            dels_r=dist_sphere(nn_x,0.d0,nn_z,x(i),0.d0,z(k))
!
            dels_u=min(dels_xu,dels_zu,dels_r)
            dels_d=min(dels_xd,dels_zd)
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
!
!upwind point intersects with sphere
               x_u=x(i)-dels_r*nn_x
               z_u=z(k)-dels_r*nn_z
               int_u=xic1
               call coeff2d_4p_lin(x(iim1), x(i), z(kkm1), z(k), x_u, z_u, &
                  acoeff, bcoeff, ccoeff, dcoeff)
               sline_u = acoeff*sline3d(iim1,j,kkm1) + bcoeff*sline3d(i,j,kkm1) + &
                  ccoeff*sline3d(iim1,j,k) + dcoeff*sline3d(i,j,k)
               opalbar_u = acoeff*opalbar3d(iim1,j,kkm1) + bcoeff*opalbar3d(i,j,kkm1) + &
                  ccoeff*opalbar3d(iim1,j,k) + dcoeff*opalbar3d(i,j,k)
               velx_u = acoeff*velx3d(iim1,j,kkm1) + bcoeff*velx3d(i,j,kkm1) + &
                  ccoeff*velx3d(iim1,j,k) + dcoeff*velx3d(i,j,k)
               vely_u = acoeff*vely3d(iim1,j,kkm1) + bcoeff*vely3d(i,j,kkm1) + &
                  ccoeff*vely3d(iim1,j,k) + dcoeff*vely3d(i,j,k)
               velz_u = acoeff*velz3d(iim1,j,kkm1) + bcoeff*velz3d(i,j,kkm1) + &
                  ccoeff*velz3d(iim1,j,k) + dcoeff*velz3d(i,j,k)
               vth_u = acoeff*vth3d(iim1,j,kkm1) + bcoeff*vth3d(i,j,kkm1) + &
                  ccoeff*vth3d(iim1,j,k) + dcoeff*vth3d(i,j,k)
!               sline_u = xic1
!               dcoeff = 0.d0
            elseif(dels_xu.eq.dels_u) then
!upwind point intersects on x-level at i-alpha
               x_u=x(iim1)
               z_u=z(k)-dels_u*nn_z
!
               opalbar_u = interpol_yp(z(kkm1), z(k), opalbar3d(iim1,j,kkm1), opalbar3d(iim1,j,k), z_u)
               sline_u =   interpol_yp(z(kkm1), z(k),   sline3d(iim1,j,kkm1),   sline3d(iim1,j,k), z_u)
               velx_u =    interpol_yp(z(kkm1), z(k),    velx3d(iim1,j,kkm1),    velx3d(iim1,j,k), z_u)
               vely_u =    interpol_yp(z(kkm1), z(k),    vely3d(iim1,j,kkm1),    vely3d(iim1,j,k), z_u)
               velz_u =    interpol_yp(z(kkm1), z(k),    velz3d(iim1,j,kkm1),    velz3d(iim1,j,k), z_u)
               vth_u =     interpol_yp(z(kkm1), z(k),     vth3d(iim1,j,kkm1),     vth3d(iim1,j,k), z_u)
               int_u =     interpol_yp(z(kkm1), z(k),     int3d(iim1,j,kkm1),     int3d(iim1,j,k), z_u)
               dcoeff=0.d0

            elseif(dels_zu.eq.dels_u) then
!upwind point intersects on z-level at k-gamma
               x_u=x(i) - dels_u*nn_x
               z_u=z(kkm1)
!
               opalbar_u = interpol_yp(x(iim1), x(i), opalbar3d(iim1,j,kkm1), opalbar3d(i,j,kkm1), x_u)
               sline_u =   interpol_yp(x(iim1), x(i),   sline3d(iim1,j,kkm1),   sline3d(i,j,kkm1), x_u)
               velx_u =    interpol_yp(x(iim1), x(i),    velx3d(iim1,j,kkm1),    velx3d(i,j,kkm1), x_u)
               vely_u =    interpol_yp(x(iim1), x(i),    vely3d(iim1,j,kkm1),    vely3d(i,j,kkm1), x_u)
               velz_u =    interpol_yp(x(iim1), x(i),    velz3d(iim1,j,kkm1),    velz3d(i,j,kkm1), x_u)
               vth_u =     interpol_yp(x(iim1), x(i),     vth3d(iim1,j,kkm1),     vth3d(i,j,kkm1), x_u)
               int_u =     interpol_yp(x(iim1), x(i),     int3d(iim1,j,kkm1),     int3d(i,j,kkm1), x_u)
               dcoeff=0.d0
!
            else
               write(*,'(4es20.8,2l4)') dels_u, dels_xu, dels_zu, dels_r, dels_u.eq.dels_xu, dels_u.eq.dels_zu, dels_u.eq.dels_r
               stop 'error in fsc_line2d: invalid dels_u'
            endif
!
!--------------------------downwind point-------------------------------
!
!            if(dels_xd.eq.dels_d) then
!!downwind point intersects on z-level
!               x_d=x(iip1)
!               z_d=z(k)+dels_d*nn_z
!               sline_d = interpol_yp(z(k), z(kkp1), sline3d(iip1,j,k), sline3d(iip1,j,kkp1), z_d)
!               opalbar_d = interpol_yp(z(k), z(kkp1), opalbar3d(iip1,j,k), opalbar3d(iip1,j,kkp1), z_d)
!               velx_d = interpol_yp(z(k), z(kkp1), velx3d(iip1,j,k), velx3d(iip1,j,kkp1), z_d)
!               vely_d = interpol_yp(z(k), z(kkp1), vely3d(iip1,j,k), vely3d(iip1,j,kkp1), z_d)
!               velz_d = interpol_yp(z(k), z(kkp1), velz3d(iip1,j,k), velz3d(iip1,j,kkp1), z_d)
!               vth_d = interpol_yp(z(k), z(kkp1), vth3d(iip1,j,k), vth3d(iip1,j,kkp1), z_d)
!            elseif(dels_zd.eq.dels_d) then
!!downwind point intersects on x-level
!               x_d=x(i)+dels_d*nn_x
!               z_d=z(kkp1)
!               sline_d = interpol_yp(x(i), x(iip1), sline3d(i,j,kkp1), sline3d(iip1,j,kkp1), x_d)
!               opalbar_d = interpol_yp(x(i), x(iip1), opalbar3d(i,j,kkp1), opalbar3d(iip1,j,kkp1), x_d)
!               velx_d = interpol_yp(x(i), x(iip1), velx3d(i,j,kkp1), velx3d(iip1,j,kkp1), x_d)
!               vely_d = interpol_yp(x(i), x(iip1), vely3d(i,j,kkp1), vely3d(iip1,j,kkp1), x_d)
!               velz_d = interpol_yp(x(i), x(iip1), velz3d(i,j,kkp1), velz3d(iip1,j,kkp1), x_d)
!               vth_d = interpol_yp(x(i), x(iip1), vth3d(i,j,kkp1), vth3d(iip1,j,kkp1), x_d)
!            else
!               write(*,'(3es20.8,2l4)') dels_d, dels_xd, dels_zd, dels_u.eq.dels_xu, dels_u.eq.dels_zu
!               stop 'error in fsc_line2d: invalid dels_u'
!            endif
!
!-----------------------radiative transfer------------------------------
!
!***debug start
            call model_debug(x_u, 0.d0, z_u, velx_u2, vely_u2, velz_u2, opalbar_u2, opac_u)
!            write(*,'(8es20.8)') opalbar_p2/opalbar_p
!***debug end
            vel_u = nn_x*velx_u + nn_y*vely_u + nn_z*velz_u
!            vel_d = nn_x*velx_d + nn_y*vely_d + nn_z*velz_d
            call fsc_line_lin(xobs, vth_fiducial, deltax, xcmf_min, xcmf_max, int_u, opalbar_u, opalbar_p, &
               sline_u, sline_p, vel_u, vel_p, vth_u, vth_p, &
               dels_u, abs_sc, contr_sc, int_sc, alo_u, alo_p)
!            call fsc_line(xobs, vth_fiducial, deltax, xcmf_min, xcmf_max, int_u, &
!                            opalbar_u, opalbar_p, opalbar_d, &
!                            sline_u, sline_p, sline_d, &
!                            vel_u, vel_p, vel_d, &
!                            vth_u, vth_p, vth_d, &
!                            dels_u, dels_d, &
!                            abs_sc, contr_sc, int_sc, alo_u, alo_p, alo_d)
            int3d(i,j,k) = int_sc
            aloline_on_nn3d(i,j,k,14) = alo_p +dcoeff*alo_u
!
!perform integration
            call calc_phinorm(vel_p, vth_p, vth_fiducial, xobs, phinorm)
            mintbar3d(i,j,k) = mintbar3d(i,j,k) + int3d(i,j,k)*wall*phinorm
            normalization3d(i,j,k) = normalization3d(i,j,k) + wall*phinorm
            aloline_nn3d(i,j,k,14) = aloline_nn3d(i,j,k,14) + aloline_on_nn3d(i,j,k,14)*wall*phinorm
!
!***********special treatment: point lies directly on boundary**********
!
          case(3)
            mueff=nn_x*x(i)+nn_z*z(k)
            if(mueff.ge.0.d0) then
               int3d(i,j,k)=xic1
               aloline_on_nn3d(i,j,k,14) = 0.d0
!
!perform integration (alo not required, since boundary specified)
               vel_p = nn_x*velx3d(i,j,k) + nn_y*vely3d(i,j,k) + nn_z*velz3d(i,j,k)
               vth_p = vth3d(i,j,k)
               call calc_phinorm(vel_p, vth_p, vth_fiducial, xobs, phinorm)
               mintbar3d(i,j,k) = mintbar3d(i,j,k) + int3d(i,j,k)*wall*phinorm
               normalization3d(i,j,k) = normalization3d(i,j,k) + wall*phinorm
!
            else
!
               iim1=i-alpha
               kkm1=k-gamma
               iip1=i+alpha
               kkp1=k+gamma
!
!calculate distance of upwind point (at previous x, z level)
               dels_xu=(x(i)-x(iim1))/nn_x
               dels_zu=(z(k)-z(kkm1))/nn_z
               dels_u=min(dels_xu,dels_zu)
!calculate distance of downwind point (at next x, z level)
               dels_xd=(x(iip1)-x(i))/nn_x
               dels_zd=(z(kkp1)-z(k))/nn_z
               dels_d=min(dels_xd,dels_zd)
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
               if(dels_xu.eq.dels_u) then
!upwind point intersects on x-level at i-alpha
                  x_u=x(iim1)
                  z_u=z(k) - dels_u*nn_z

                  opalbar_u = interpol_yp(z(kkm1), z(k), opalbar3d(iim1,j,kkm1), opalbar3d(iim1,j,k), z_u)
                  sline_u =   interpol_yp(z(kkm1), z(k),   sline3d(iim1,j,kkm1),   sline3d(iim1,j,k), z_u)
                  velx_u =    interpol_yp(z(kkm1), z(k),    velx3d(iim1,j,kkm1),    velx3d(iim1,j,k), z_u)
                  vely_u =    interpol_yp(z(kkm1), z(k),    vely3d(iim1,j,kkm1),    vely3d(iim1,j,k), z_u)
                  velz_u =    interpol_yp(z(kkm1), z(k),    velz3d(iim1,j,kkm1),    velz3d(iim1,j,k), z_u)
                  vth_u =     interpol_yp(z(kkm1), z(k),     vth3d(iim1,j,kkm1),     vth3d(iim1,j,k), z_u)
                  int_u =     interpol_yp(z(kkm1), z(k),     int3d(iim1,j,kkm1),     int3d(iim1,j,k), z_u)

               elseif(dels_zu.eq.dels_u) then
!upwind point intersects on z-level at k-gamma
                  x_u=x(i) - dels_u*nn_x
                  z_u=z(kkm1)
!
                  opalbar_u = interpol_yp(x(iim1), x(i), opalbar3d(iim1,j,kkm1), opalbar3d(i,j,kkm1), x_u)
                  sline_u =   interpol_yp(x(iim1), x(i),   sline3d(iim1,j,kkm1),   sline3d(i,j,kkm1), x_u)
                  velx_u =    interpol_yp(x(iim1), x(i),    velx3d(iim1,j,kkm1),    velx3d(i,j,kkm1), x_u)
                  vely_u =    interpol_yp(x(iim1), x(i),    vely3d(iim1,j,kkm1),    vely3d(i,j,kkm1), x_u)
                  velz_u =    interpol_yp(x(iim1), x(i),    velz3d(iim1,j,kkm1),    velz3d(i,j,kkm1), x_u)
                  vth_u =     interpol_yp(x(iim1), x(i),     vth3d(iim1,j,kkm1),     vth3d(i,j,kkm1), x_u)
                  int_u =     interpol_yp(x(iim1), x(i),     int3d(iim1,j,kkm1),     int3d(i,j,kkm1), x_u)
!
               else
                  write(*,'(3es20.8,2l4)') dels_u, dels_xu, dels_zu, dels_u.eq.dels_xu, dels_u.eq.dels_zu
                  stop 'error in fsc_line2d: invalid dels_u'
               endif
!
!----------------------------downwind point-----------------------------
!
!               if(dels_xd.eq.dels_d) then
!!upwind point intersects on x-level at i+alpha
!                  x_d=x(iip1)
!                  z_d=z(k) + dels_d*nn_z
!
!                  opalbar_d = interpol_yp(z(k), z(kkp1), opalbar3d(iip1,j,k), opalbar3d(iip1,j,kkp1), z_d)
!                  sline_d =   interpol_yp(z(k), z(kkp1),   sline3d(iip1,j,k),   sline3d(iip1,j,kkp1), z_d)
!                  velx_d =    interpol_yp(z(k), z(kkp1),    velx3d(iip1,j,k),    velx3d(iip1,j,kkp1), z_d)
!                  vely_d =    interpol_yp(z(k), z(kkp1),    vely3d(iip1,j,k),    vely3d(iip1,j,kkp1), z_d)
!                  velz_d =    interpol_yp(z(k), z(kkp1),    velz3d(iip1,j,k),    velz3d(iip1,j,kkp1), z_d)
!                  vth_d =     interpol_yp(z(k), z(kkp1),     vth3d(iip1,j,k),     vth3d(iip1,j,kkp1), z_d)
!
!               elseif(dels_zd.eq.dels_d) then
!!upwind point intersects on z-level at k+gamma
!                  x_d=x(i) + dels_d*nn_x
!                  z_d=z(kkp1)
!!
!                  opalbar_d = interpol_yp(x(i), x(iip1), opalbar3d(i,j,kkp1), opalbar3d(iip1,j,kkp1), x_d)
!                  sline_d =   interpol_yp(x(i), x(iip1),   sline3d(i,j,kkp1),   sline3d(iip1,j,kkp1), x_d)
!                  velx_d =    interpol_yp(x(i), x(iip1),    velx3d(i,j,kkp1),    velx3d(iip1,j,kkp1), x_d)
!                  vely_d =    interpol_yp(x(i), x(iip1),    vely3d(i,j,kkp1),    vely3d(iip1,j,kkp1), x_d)
!                  velz_d =    interpol_yp(x(i), x(iip1),    velz3d(i,j,kkp1),    velz3d(iip1,j,kkp1), x_d)
!                  vth_d =     interpol_yp(x(i), x(iip1),     vth3d(i,j,kkp1),     vth3d(iip1,j,kkp1), x_d)
!!
!               else
!                  write(*,'(3es20.8,2l4)') dels_d, dels_xd, dels_zd, dels_d.eq.dels_xd, dels_d.eq.dels_zd
!                  stop 'error in fsc_line2d: invalid dels_u'
!               endif
!
!-----------------------radiative transfer------------------------------
!
!***debug start
               call model_debug(x_u, 0.d0, z_u, velx_u2, vely_u2, velz_u2, opalbar_u2, opac_u)
!            write(*,'(8es20.8)') opalbar_p2/opalbar_p
!***debug end
               vel_u = nn_x*velx_u + nn_y*vely_u + nn_z*velz_u
!               vel_d = nn_x*velx_d + nn_y*vely_d + nn_z*velz_d
               call fsc_line_lin(xobs, vth_fiducial, deltax, xcmf_min, xcmf_max, int_u, opalbar_u, opalbar_p, &
                  sline_u, sline_p, vel_u, vel_p, vth_u, vth_p, &
                  dels_u, abs_sc, contr_sc, int_sc, alo_u, alo_p)
!               call fsc_line(xobs, vth_fiducial, deltax, xcmf_min, xcmf_max,, int_u, &
!                               opalbar_u, opalbar_p, opalbar_d, &
!                               sline_u, sline_p, sline_d, &
!                               vel_u, vel_p, vel_d, &
!                               vth_u, vth_p, vth_d, &
!                               dels_u, dels_d, &
!                               abs_sc, contr_sc, int_sc, alo_u, alo_p, alo_d)
               int3d(i,j,k) = int_sc
               aloline_on_nn3d(i,j,k,14) = alo_p
!
!perform integration
               call calc_phinorm(vel_p, vth_p, vth_fiducial, xobs, phinorm)
               mintbar3d(i,j,k) = mintbar3d(i,j,k) + int3d(i,j,k)*wall*phinorm
               normalization3d(i,j,k) = normalization3d(i,j,k) + wall*phinorm
               aloline_nn3d(i,j,k,14) = aloline_nn3d(i,j,k,14) + aloline_on_nn3d(i,j,k,14)*wall*phinorm

            endif
!
          case default
!
         end select
!

      enddo
   enddo
!
!***debug
!close(1)
!
end subroutine fsc_line2d
