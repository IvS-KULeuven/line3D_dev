subroutine grid_delxyz
!
!-----------------------------------------------------------------------
!------------same as grid_delxyz but for symmetric deltas---------------
!-----------------------------------------------------------------------
!
use prog_type
use dime3d, only: ndxmax, ndymax, ndzmax, x, y, z, delx_arr, dely_arr, delz_arr, &
                  imask_totreg3d, imask_innreg3d, imask_bpoint3d, opac3d
use options, only: opt_incl_sdist
!
implicit none
!
! ... local scalars
integer(i4b) :: i,j,k, err
integer(i4b) :: startx, starty, startz, endx, endy, endz, pos_indx
integer(i4b) :: alpha, beta, gamma, aalpha, bbeta, ggamma
real(dp) :: grad, x_phantom, y_phantom, z_phantom
real(dp) :: xn, yn, zn, delx, dely, delz
!
! ... loacl logicals
!
!
write(*,*) '--------------------calculating delx, dely, delz (array)-----------------------'
write(*,*)

if(opt_incl_sdist) stop 'error in grid_delxyz: surface distortion does not work with FVM yet'
!
allocate(delx_arr(ndxmax,ndymax,ndzmax,2), stat=err)
   if(err.ne.0) stop 'allocation error grid_delxyz'
allocate(dely_arr(ndxmax,ndymax,ndzmax,2), stat=err)
   if(err.ne.0) stop 'allocation error grid_delxyz'
allocate(delz_arr(ndxmax,ndymax,ndzmax,2), stat=err)
   if(err.ne.0) stop 'allocation error grid_delxyz'
!
!-------------------for positive n_x, n_y, n_z--------------------------
!
!need delx, dely, delz-arrays only on index-range [2,nd(xyz)max-1] 
!   because of phantom point
delx_arr=0.d0
dely_arr=0.d0
delz_arr=0.d0
!
alpha= 1
beta= 1
gamma= 1
startx= 2
starty= 2
startz= 2
endx=ndxmax-1
endy=ndymax-1
endz=ndzmax-1
pos_indx=1
!
delx_arr=0.d0
dely_arr=0.d0
delz_arr=0.d0
!
!
do i=startx, endx, alpha
   do j=starty, endy, beta
      do k=startz, endz, gamma
!
         if(imask_totreg3d(i,j,k).eq.1) then
!store delx, dely, delz only if information region
!
!----------------------------calculate delx-----------------------------
!
!calculate central (or backward) differences on complete grid
            delx= 0.5d0*(x(i+alpha)-x(i-alpha))
            !delx= (x(i)-x(i-alpha))
!
!check if inner boundary is hit, then overwrite central differences from above
            if(imask_innreg3d(i-alpha,j,k).eq.1) then
               !calculating intersection point with photosphere
               !note: only for spherical star
               xn=alpha*sqrt(1.d0-y(j)**2-z(k)**2)
               !delx=0.5d0*(x(i+alpha)-xn)
               delx=x(i)-xn
            end if
!check if inner boundary is hit, then overwrite central differences from above
            if(imask_innreg3d(i+alpha,j,k).eq.1) then
               !calculating intersection point with photosphere
               !note: only for spherical star
               !xn=-1.d0*alpha*sqrt(1.d0-y(j)*y(j)-z(k)*z(k))
               !delx=0.5d0*(xn-x(i-alpha))
               delx=x(i)-x(i-alpha)
            end if
!
!----------------------------calculate dely-----------------------------
!
!calculate central (or backward) differences
            dely= 0.5d0*(y(j+beta)-y(j- beta))
            !dely= (y(j)-y(j-beta))
!
!check if inner boundary is hit, then overwrite central differences from above
            if(imask_innreg3d(i,j-beta,k).eq.1) then
               !calculating intersection point with photosphere
               !note: only for spherical star
               yn=beta*sqrt(1.d0-x(i)*x(i)-z(k)*z(k))
               !dely=0.5d0*(y(j+beta)-yn)
               dely=y(j)-yn
            end if
!
!check if inner boundary is hit, then overwrite central differences from above
            if(imask_innreg3d(i,j+beta,k).eq.1) then
               !calculating intersection point with photosphere
               !note: only for spherical star
               !yn=-1.d0*beta*sqrt(1.d0-x(i)*x(i)-z(k)*z(k))
               !dely=0.5d0*(yn-y(j-beta))
               dely=y(j)-y(j-beta)
            end if
!
!----------------------------calculate delz-----------------------------
!
!calculate central (or backward) differences
            delz= 0.5d0*(z(k+gamma)-z(k- gamma))
            !delz= (z(k)-z(k-gamma))
!
!check if inner boundary is hit, then overwrite central differences from above
            if(imask_innreg3d(i,j,k-gamma).eq.1) then
               !calculating intersection point with photosphere
               !note: only for spherical star
               zn=gamma*sqrt(1.d0-x(i)*x(i)-y(j)*y(j))
               !delz=0.5d0*(z(k+gamma)-zn)
               delz=z(k)-zn
            end if
!
!check if inner boundary is hit, then overwrite central differences from above
            if(imask_innreg3d(i,j,k+gamma).eq.1) then
               !calculating intersection point with photosphere
               !note: only for spherical star
               zn=-1.d0*gamma*sqrt(1.d0-x(i)*x(i)-y(j)*y(j))
               !delz=0.5d0*(zn-z(k-gamma))
               delz=z(k)-z(k-gamma)
            end if
!
!----------------check if point lies directly on boundary---------------
!
            !directly on boundaries: backward differences
            if(imask_bpoint3d(i,j,k).eq.1) then
               aalpha=nint(sign(1.d0, x(i)))
               bbeta=nint(sign(1.d0, y(j)))
               ggamma=nint(sign(1.d0, z(k)))
!
               delx=aalpha*(x(i+aalpha)-x(i))
               dely=bbeta*(y(j+bbeta)-y(j))
               delz=ggamma*(z(k+ggamma)-z(k))
            end if
!
            delx_arr(i,j,k, pos_indx)=delx
            dely_arr(i,j,k, pos_indx)=dely
            delz_arr(i,j,k, pos_indx)=delz
!
         endif
!
      enddo
   enddo
enddo
!
!-------------------for negative n_x, n_y, n_z--------------------------
!
alpha= -1
beta= -1
gamma= -1
startx= ndxmax-1
starty= ndymax-1
startz= ndzmax-1
endx=2
endy=2
endz=2
pos_indx=2
!
!
!
do i=startx, endx, alpha
   do j=starty, endy, beta
      do k=startz, endz, gamma

         if(imask_totreg3d(i,j,k).eq.1) then
!store delx, dely, delz only if information region
!
!----------------------------calculate delx-----------------------------
!
!calculate central (or backward) differences
            delx= 0.5d0*(x(i+alpha)-x(i-alpha))
            !delx= (x(i)-x(i-alpha))
!
!check if inner boundary is hit, then overwrite central differences from above
            if(imask_innreg3d(i-alpha,j,k).eq.1) then
               !calculating intersection point with photosphere
               !note: only for spherical star
               xn=alpha*sqrt(1.d0-y(j)*y(j)-z(k)*z(k))
               !delx=0.5d0*(x(i+alpha)-xn)
               delx=x(i)-xn
            end if
!
!check if inner boundary is hit, then overwrite central differences from above
            if(imask_innreg3d(i+alpha,j,k).eq.1) then
               !calculating intersection point with photosphere
               !note: only for spherical star
               !xn=-1.d0*alpha*sqrt(1.d0-y(j)*y(j)-z(k)*z(k))
               !delx=0.5d0*(xn-x(i-alpha))
               delx=x(i)-x(i-alpha)
            end if
!
!----------------------------calculate dely-----------------------------
!
!calculate central (or backward) differences
            dely= 0.5d0*(y(j+beta)-y(j- beta))
            !dely= (y(j)-y(j-beta))
!
!check if inner boundary is hit, then overwrite central differences from above
            if(imask_innreg3d(i,j-beta,k).eq.1) then
               !calculating intersection point with photosphere
               !note: only for spherical star
               yn=beta*sqrt(1.d0-x(i)*x(i)-z(k)*z(k))
               !dely=0.5d0*(y(j+beta)-yn)
               dely=y(j)-yn
            end if
!
!check if inner boundary is hit, then overwrite central differences from above
            if(imask_innreg3d(i,j+beta,k).eq.1) then
               !calculating intersection point with photosphere
               !note: only for spherical star
               !yn=-1.d0*beta*sqrt(1.d0-x(i)*x(i)-z(k)*z(k))
               !dely=0.5d0*(yn-y(j-beta))
               dely=y(j)-y(j-beta)
            end if
!
!----------------------------calculate delz-----------------------------
!
!calculate central (or backward) differences
            delz= 0.5d0*(z(k+gamma)-z(k- gamma))
            !delz= (z(k)-z(k-gamma))
!
!check if inner boundary is hit, then overwrite central differences from above
            if(imask_innreg3d(i,j,k-gamma).eq.1) then
               !calculating intersection point with photosphere
               !note: only for spherical star
               zn=gamma*sqrt(1.d0-x(i)*x(i)-y(j)*y(j))
               !delz=0.5d0*(z(k+gamma)-zn)
               delz=z(k)-zn
            endif
!
!check if inner boundary is hit, then overwrite central differences from above
            if(imask_innreg3d(i,j,k+gamma).eq.1) then
               !calculating intersection point with photosphere
               !note: only for spherical star
               !zn=-1.d0*gamma*sqrt(1.d0-x(i)*x(i)-y(j)*y(j))
               !delz=0.5d0*(zn-z(k-gamma))
               delz=z(k)-z(k-gamma)
            end if
!
!----------------check if point lies directly on boundary---------------
!
            !directly on boundaries: backward differences
            if(imask_bpoint3d(i,j,k).eq.1) then
               aalpha=nint(sign(1.d0, x(i)))
               bbeta=nint(sign(1.d0, y(j)))
               ggamma=nint(sign(1.d0, z(k)))
!
               delx=aalpha*(x(i)-x(i+aalpha))
               dely=bbeta*(y(j)-y(j+bbeta))
               delz=ggamma*(z(k)-z(k+ggamma))
            end if
!
            delx_arr(i,j,k, pos_indx)=delx
            dely_arr(i,j,k, pos_indx)=dely
            delz_arr(i,j,k, pos_indx)=delz
!
         endif
!
      enddo
   enddo
enddo
!
!-----------------------------------------------------------------------
!
end subroutine grid_delxyz
