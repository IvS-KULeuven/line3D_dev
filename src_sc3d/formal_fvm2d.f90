subroutine ffvm_cont2d(muindx)
!
!-----------------------------------------------------------------------
!---------------------finite volume method in 2d------------------------
!-----------calculating intensties for given mu, phi=0------------------
!-----------------------------------------------------------------------
!
use prog_type

use fund_const
use dime3d, only: ndxmax, ndymax, ndzmax, x, y, z, delx_arr, dely_arr, delz_arr
use dime3d, only: int3d, opac3d, scont3d, imask_totreg3d, imask_innreg3d, imask_bpoint3d, alocont_nn3d
use angles, only: nodes_mu, dim_mu
use bcondition, only: xic1
!
implicit none
!
! ... arguments
integer(i4b) :: muindx
!
! ... local scalars
integer(i4b) :: i, j, k
integer(i4b) :: alpha, beta, gamma, aalpha, bbeta, ggamma
integer(i4b) :: startx, starty, startz, endx, endy, endz
integer(i4b) :: posx_indx, posy_indx, posz_indx
real(dp) :: mueff
real(dp) :: xn, yn, zn, intm, intn, into, alocontm, alocontn, aloconto, &
            phibm, phibn, phibo
real(dp) :: denominator, nn_x, nn_y, nn_z
real(dp) :: duma, dumb, dumc, dumd, consta, constb, constc, constd
!
! ... local arrays
!
! ... local functions
!
! ... local characters
character(len=50) :: enter
!
! ... local logicals
!
!----------calculate intensities via normal difference method-----------
!
nn_x=sqrt(1.d0 - nodes_mu(muindx)**2)
nn_y=0.d0
nn_z=nodes_mu(muindx)
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
!                beta = 1                    beta = -1          
!
!         if n_z >= 0                 if n_z < 0                
!                startz = 2                  startz = ndzmax-1  
!                endz = ndzmax-1             endz = 2           
!                gamma = 1                   gamma = -1
!
if (nn_x.ge.0.d0) then
   startx = 2
   endx = ndxmax-1
   alpha=  1
   aalpha= 15
   posx_indx=1
else 
   startx = ndxmax-1
   endx = 2
   alpha=-1
   aalpha= 13
   posx_indx=2
endif
!
if (nn_y.ge.0.d0) then
   starty = 2
   endy = ndymax-1
   beta=  1
   bbeta= 17
   posy_indx=1
else
   starty = ndymax-1
   endy = 2
   beta =-1
   bbeta= 11
   posy_indx=2
endif
!
if (nn_z.ge.0.d0) then
   startz = 2
   endz = ndzmax-1
   gamma=  1
   ggamma= 23
   posz_indx=1
else
   startz = ndzmax-1
   endz = 2
   gamma=-1
   ggamma= 5
   posz_indx=2
endif
!
int3d=0.d0
alocont_nn3d=0.d0
!
!-----------------------------------------------------------------------
!
j=ndymax/2+1
do k=startz, endz, gamma
      do i=startx, endx, alpha
!
      if(imask_totreg3d(i,j,k).eq.1) then
!
!-----calculating boundary intensities on intermesh-points if needed----
!
         if(imask_innreg3d(i-alpha,j,k).eq.1) then
            xn=alpha*sqrt(1.d0-y(j)**2-z(k)**2)
            mueff=xn*nn_x + y(j)*nn_y + z(k) * nn_z
            if(mueff.lt.0.d0) then
               intn = 0.d0
            else
               intn = xic1
            endif
         else
           intn=int3d(i-alpha,j,k)
         end if
!
         if(imask_innreg3d(i,j-beta,k).eq.1) then
            yn=beta*sqrt(1.d0-x(i)*x(i)-z(k)*z(k))
            mueff=x(i)*nn_x + yn*nn_y + z(k)*nn_z
            if(mueff.lt.0.d0) then
               intm=0.d0
            else
               intm = xic1
            endif
         else
            intm=int3d(i,j-beta, k)
         end if
!
         if(imask_innreg3d(i, j, k-gamma).eq.1) then
            zn=gamma*sqrt(1.d0-x(i)*x(i)-y(j)*y(j))
            mueff=x(i)*nn_x + y(j)*nn_y + zn*nn_z
            if(mueff.lt.0.d0) then
               into = 0.d0
            else
               into = xic1
            endif
         else
            into=int3d(i,j, k-gamma)
         end if
!
         duma=opac3d(i,j,k)
         dumb=nn_x/delx_arr(i,j,k,posx_indx)
         dumc=nn_y/dely_arr(i,j,k,posy_indx)
         dumd=nn_z/delz_arr(i,j,k,posz_indx)
!
         denominator= duma+dumb+dumc+dumd
         consta=duma/denominator
         constb=dumb/denominator
         constc=dumc/denominator
         constd=dumd/denominator
!
!-------------calculate intensities via difference method---------------
!
         int3d(i,j,k) = consta * scont3d(i,j,k) + constb * intn + &
                        constc * intm + constd * into
!
!-----------calculate alo (per direction) for nearest neighbours--------
!
         alocont_nn3d(i,j,k,14)=consta
!         alocont_nn3d(i-alpha,j,k,aalpha)=alocont_nn3d(i-alpha,j,k,14)*constb
!         alocont_nn3d(i,j-beta,k,bbeta)=alocont_nn3d(i,j-beta,k,14)*constc
!         alocont_nn3d(i,j,k-gamma,ggamma) = alocont_nn3d(i,j,k-gamma,14)*constd
!
!-----------calculating boundary intensities on grid-points-------------
!-------------overwriting intensities calculated above------------------
!
         if(imask_bpoint3d(i,j,k).eq.1) then
            !calculate mueff
            mueff=nn_x*x(i) + nn_y*y(j) + nn_z*z(k)
            if(mueff.ge.0.d0) then
               int3d(i,j,k)=  xic1
               alocont_nn3d(i,j,k,:)= 0.d0
!               alocont_nn3d(i-alpha,j,k,aalpha)=0.d0
!               alocont_nn3d(i,j-beta,k,bbeta)=0.d0
!               alocont_nn3d(i,j,k-gamma,ggamma) =0.d0
            end if
         end if
      end if
!
      enddo
!   enddo
enddo
!
!
end subroutine ffvm_cont2d
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine ffvm_cont2d_debug(muindx, phiindx)
!
!-----------------------------------------------------------------------
!---------------------finite volume method in 2d------------------------
!-----------calculating intensties for given mu, phi=0------------------
!------------------storing absorption and emission part explicitly
!-----------------------------------------------------------------------
!
use prog_type

use fund_const
use dime3d, only: ndxmax, ndymax, ndzmax, x, y, z, delx_arr, dely_arr, delz_arr
use dime3d, only: int3d, opac3d, scont3d, imask_totreg3d, imask_innreg3d, imask_bpoint3d
use angles, only: nodes_mu, dim_mu
use bcondition, only: xic1
use mod_benchmark, only: int2d_fvm, abs2d_fvm, contr2d_fvm
!
implicit none
!
! ... arguments
integer(i4b) :: muindx, phiindx
!
! ... local scalars
integer(i4b) :: i, j, k
integer(i4b) :: alpha, beta, gamma, aalpha, bbeta, ggamma
integer(i4b) :: startx, starty, startz, endx, endy, endz
integer(i4b) :: posx_indx, posy_indx, posz_indx
real(dp) :: mueff
real(dp) :: xn, yn, zn, intm, intn, into, alocontm, alocontn, aloconto, &
            phibm, phibn, phibo
real(dp) :: denominator, nn_x, nn_y, nn_z
real(dp) :: duma, dumb, dumc, dumd, consta, constb, constc, constd
!
! ... local arrays
!
! ... local functions
!
! ... local characters
character(len=50) :: enter
!
! ... local logicals
!
!----------calculate intensities via normal difference method-----------
!
!nn_x=n_x(muindx,phiindx)
!nn_y=n_y(muindx,phiindx)
!nn_z=n_z(muindx,phiindx)
nn_x=sqrt(1.d0 - nodes_mu(muindx)*nodes_mu(muindx))
nn_y=0.d0
nn_z=nodes_mu(muindx)
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
!                beta = 1                    beta = -1          
!
!         if n_z >= 0                 if n_z < 0                
!                startz = 2                  startz = ndzmax-1  
!                endz = ndzmax-1             endz = 2           
!                gamma = 1                   gamma = -1
!
if (nn_x.ge.0.d0) then
   startx = 2
   endx = ndxmax-1
   alpha=  1
   aalpha= 15
   posx_indx=1
else 
   startx = ndxmax-1
   endx = 2
   alpha=-1
   aalpha= 13
   posx_indx=2
endif
!
if (nn_y.ge.0.d0) then
   starty = 2
   endy = ndymax-1
   beta=  1
   bbeta= 17
   posy_indx=1
else
   starty = ndymax-1
   endy = 2
   beta =-1
   bbeta= 11
   posy_indx=2
endif
!
if (nn_z.ge.0.d0) then
   startz = 2
   endz = ndzmax-1
   gamma=  1
   ggamma= 23
   posz_indx=1
else
   startz = ndzmax-1
   endz = 2
   gamma=-1
   ggamma= 5
   posz_indx=2
endif
!
int3d=0.d0
!alocont_nn3d=0.d0
!
!-----------------------------------------------------------------------
!
do k=startz, endz, gamma
!   do j=starty, endy, beta
   j=ndymax/2+1
      do i=startx, endx, alpha
!
      if(imask_totreg3d(i,j,k).eq.1) then
!
!-----calculating boundary intensities on intermesh-points if needed----
!
         if(imask_innreg3d(i-alpha,j,k).eq.1) then
            intn = xic1
         else
           intn=int3d(i-alpha,j,k)
         end if
!
         if(imask_innreg3d(i,j-beta,k).eq.1) then
            intm = xic1
         else
            intm=int3d(i,j-beta, k)
         end if
!
         if(imask_innreg3d(i, j, k-gamma).eq.1) then
            into = xic1
         else
            into=int3d(i,j, k-gamma)
         end if
!
         duma=opac3d(i,j,k)
         dumb=nn_x/delx_arr(i,j,k,posx_indx)
         dumc=nn_y/dely_arr(i,j,k,posy_indx)
         dumd=nn_z/delz_arr(i,j,k,posz_indx)
!
         denominator= duma+dumb+dumc+dumd
         consta=duma/denominator
         constb=dumb/denominator
         constc=dumc/denominator
         constd=dumd/denominator
!
!-------------calculate intensities via difference method---------------
!
         int3d(i,j,k) = consta * scont3d(i,j,k) + constb * intn + &
                      constc * intm + constd * into
         int2d_fvm(i,k) = int3d(i,j,k)
         abs2d_fvm(i,k) = constb + constc + constd
         contr2d_fvm(i,k) = consta * scont3d(i,j,k)
!
!-----------calculate alo (per direction) for nearest neighbours--------
!
!         alocont_nn3d(i,j,k,14)=consta
!         alocont_nn3d(i-alpha,j,k,aalpha)=alocont_nn3d(i-alpha,j,k,14)*constb
!         alocont_nn3d(i,j-beta,k,bbeta)=alocont_nn3d(i,j-beta,k,14)*constc
!         alocont_nn3d(i,j,k-gamma,ggamma) = alocont_nn3d(i,j,k-gamma,14)*constd
!
!-----------calculating boundary intensities on grid-points-------------
!-------------overwriting intensities calculated above------------------
!
         if(imask_bpoint3d(i,j,k).eq.1) then
            !calculate mueff
            mueff=nn_x*x(i) + nn_y*y(j) + nn_z*z(k)
            if(mueff.ge.0.d0) then
               int3d(i,j,k)=  xic1
               int2d_fvm(i,k) = xic1
!               alocont_nn3d(i,j,k,:)= 0.d0
!               alocont_nn3d(i-alpha,j,k,aalpha)=0.d0
!               alocont_nn3d(i,j-beta,k,bbeta)=0.d0
!               alocont_nn3d(i,j,k-gamma,ggamma) =0.d0
            end if
         end if
      end if
!
      enddo
!   enddo
enddo
!
!
end subroutine ffvm_cont2d_debug
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine ffvm_line2d(muindx, xobsindx)
!
!-----------------------------------------------------------------------
!-----------------------difference method-------------------------------
!-----------calculating intensties for given mu, xobs-------------------
!-------------------------in x-z-plane----------------------------------
!-----------------------------------------------------------------------
!
use prog_type
use fund_const
use dime3d, only: ndxmax, ndymax, ndzmax, x, y, z, delx_arr, dely_arr, delz_arr, &
                  int3d, imask_totreg3d, imask_innreg3d, imask_bpoint3d, &
                  opac3d, opalbar3d, scont3d, sline3d, mintbar3d, velx3d, vely3d, velz3d, &
                  vth3d, normalization3d, aloline_nn3d, aloline_on_nn3d
use bcondition, only: xic1
use params_input, only: vth_fiducial
use angles, only: nodes_mu, weight_mu
use freq, only: nodes_xobs, weight_xobs
!
implicit none
!
! ... arguments
integer(i4b), intent(in) :: muindx, xobsindx
!
! ... local scalars
integer(i4b) :: i, j, k
integer(i4b) :: alpha, beta, gamma, aalpha, bbeta, ggamma
integer(i4b) :: startx, starty, startz, endx, endy, endz
integer(i4b) :: posx_indx, posy_indx, posz_indx
real(dp) :: mueff, val_profile, velz
real(dp) :: xn, yn, zn, intm, intn, into
real(dp) :: nn_x, nn_y, nn_z, xobs
real(dp) :: duma, dumb, dumc, dumd, dume
real(dp) :: consta, constb, constc, constd, conste, denominator
real(dp) :: wmu, wxobs
real(dp) :: ic
!
! ... local arrays
!
! ... local functions
!
! ... local characters
character(len=50) :: enter
!
! ... local logicals
!
!----------calculate intensities via normal difference method-----------
!
nn_x=sqrt(1.d0 - nodes_mu(muindx)*nodes_mu(muindx))
nn_y=0.d0
nn_z=nodes_mu(muindx)
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
!                beta = 1                    beta = -1          
!
!         if n_z >= 0                 if n_z < 0                
!                startz = 2                  startz = ndzmax-1  
!                endz = ndzmax-1             endz = 2           
!                gamma = 1                   gamma = -1
!
if (nn_x.ge.0.d0) then
   startx = 2
   endx = ndxmax-1
   alpha=  1
   aalpha= 15
   posx_indx=1
else 
   startx = ndxmax-1
   endx = 2
   alpha=-1
   aalpha= 13
   posx_indx=2
endif
!
if (nn_y.ge.0.d0) then
   starty = 2
   endy = ndymax-1
   beta=  1
   bbeta= 17
   posy_indx=1
else
   starty = ndymax-1
   endy = 2
   beta =-1
   bbeta= 11
   posy_indx=2
endif
!
if (nn_z.ge.0.d0) then
   startz = 2
   endz = ndzmax-1
   gamma=  1
   ggamma= 23
   posz_indx=1
else
   startz = ndzmax-1
   endz = 2
   gamma=-1
   ggamma= 5
   posz_indx=2
endif
!
int3d=0.d0
aloline_on_nn3d=0.d0
!
xobs=nodes_xobs(xobsindx)
!
wmu=weight_mu(muindx)
wxobs=weight_xobs(xobsindx)
!
!-----------------------------------------------------------------------
!
j=ndymax/2+1
!
do k=startz, endz, gamma
   do i=startx, endx, alpha
!
      if(imask_totreg3d(i,j,k).eq.1) then
!
!-----calculating boundary intensities on intermesh-points if needed----
!
         if(imask_innreg3d(i-alpha,j,k).eq.1) then
            intn = xic1
         else
           intn=int3d(i-alpha,j,k)
         end if
!
         if(imask_innreg3d(i,j-beta,k).eq.1) then
            intm = xic1
         else
            intm=int3d(i,j-beta, k)
         end if
!
         if(imask_innreg3d(i,j,k-gamma).eq.1) then
            into = xic1
         else
            into=int3d(i,j, k-gamma)
         end if
!
!-------------------calculate normalized profile function---------------
!
         velz = velx3d(i,j,k)*nn_x + vely3d(i,j,k)*nn_y + velz3d(i,j,k)*nn_z
         call calc_phinorm(velz, vth3d(i,j,k), vth_fiducial, xobs, val_profile)
!
!*****************************new version*******************************
!
         duma=opac3d(i,j,k)
         dumb=opalbar3d(i,j,k)*val_profile
         dumc=nn_x/delx_arr(i,j,k,posx_indx)
         dumd=nn_y/dely_arr(i,j,k,posy_indx)
         dume=nn_z/delz_arr(i,j,k,posz_indx)
!
         denominator=duma+dumb+dumc+dumd+dume
!
         consta=duma/denominator
         constb=dumb/denominator
         constc=dumc/denominator
         constd=dumd/denominator
         conste=dume/denominator
!
!-------------calculate intensities via difference method---------------
!
         int3d(i,j,k) = consta*scont3d(i,j,k) + &
                        constb*sline3d(i,j,k) + &
                        constc*intn + constd*intm + conste*into
!
!-------calculate nearest neighbour alo per direction and frequency-----
!
          aloline_on_nn3d(i,j,k,14)=constb
!          aloline_on_nn3d(i-alpha,j,k,aalpha)=aloline_on_nn3d(i-alpha,j,k,14)*constc
!          aloline_on_nn3d(i,j-beta,k,bbeta)=aloline_on_nn3d(i,j-beta,k,14)*constd
!          aloline_on_nn3d(i,j,k-gamma,ggamma) = aloline_on_nn3d(i,j,k-gamma,14)*conste
!
!-----------calculate boundary intensities on grid-points---------------
!---------------overwrite intensities calculated above------------------
!
         if(imask_bpoint3d(i,j,k).eq.1) then
            !calculate mueff
            mueff=nn_x*x(i) + nn_y*y(j) + nn_z*z(k)
            if(mueff.ge.0.d0) then
               int3d(i,j,k)=  xic1
                aloline_on_nn3d(i,j,k,:)= 0.d0
!               aloline_on_nn3d(i-alpha,j,k,aalpha)=0.d0
!               aloline_on_nn3d(i,j-beta,k,bbeta)=0.d0
!               aloline_on_nn3d(i,j,k-gamma,ggamma) =0.d0
            endif
         endif
!
!----------------perform angular and frequency integration--------------
!
         mintbar3d(i,j,k) = mintbar3d(i,j,k) + val_profile * int3d(i,j,k) * wmu * wxobs
!         if(i.eq.ndxmax/2+1.and.k.eq.20) then
!            write(*,'(9es20.8)') wmu, wxobs, val_profile, int3d(i,j,k), nn_x, nn_z, x(i), z(k), mintbar3d(i,j,k)
!         endif
!
         aloline_nn3d(i,j,k,14) = aloline_nn3d(i,j,k,14) + val_profile * aloline_on_nn3d(i,j,k,14) * &
                                         wmu * wxobs
!
         normalization3d(i,j,k) = normalization3d(i,j,k) + val_profile * wmu * wxobs
!!
      endif
!
   enddo
enddo
!
!
!
end subroutine ffvm_line2d
