subroutine ffvm_cont3d(oindx)
!
!-----------------------------------------------------------------------
!---------------------finite volume method in 3d------------------------
!-----------calculating intensties for given mu, phi--------------------
!---------------------specified by inpud oindx--------------------------
!-----------------------------------------------------------------------
!
use prog_type

use fund_const
use dime3d, only: ndxmax, ndymax, ndzmax, x, y, z, delx_arr, dely_arr, delz_arr
use dime3d, only: int3d, alocont_o_nn3d, opac3d, scont3d, imask_totreg3d, imask_innreg3d, imask_bpoint3d, &
                  mint3d_tmp, normalization3d_tmp, alocont_nn3d_tmp, &
                  fcontx3d_tmp, fconty3d_tmp, fcontz3d_tmp, &
                  kcontxx3d_tmp, kcontyy3d_tmp, kcontzz3d_tmp, kcontxy3d_tmp, kcontxz3d_tmp, kcontyz3d_tmp
use angles, only: n_x, n_y, n_z, weight_omega
use bcondition, only: xic1
!
implicit none
!
! ... arguments
integer(i4b), intent(in) :: oindx
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
real(dp) :: wall
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
nn_x=n_x(oindx)
nn_y=n_y(oindx)
nn_z=n_z(oindx)
!
!integration weight for solid angle
wall=weight_omega(oindx)
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
alocont_o_nn3d=0.d0
!
!-----------------------------------------------------------------------
!
do k=startz, endz, gamma
   do j=starty, endy, beta
      do i=startx, endx, alpha
!
         if(imask_totreg3d(i,j,k).eq.1) then
!
!-----calculating boundary intensities on intermesh-points if needed----
!
            if(imask_innreg3d(i-alpha,j,k).eq.1) then
               xn=alpha*sqrt(1.d0-y(j)*y(j)-z(k)*z(k))
               mueff=xn*nn_x + y(j)*nn_y + z(k) * nn_z
               if(mueff.lt.0.d0) then
                  intn=0.d0
               else
                  intn=xic1
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
                  intm=xic1
               endif
            else
               intm=int3d(i,j-beta, k)
            end if
!
            if(imask_innreg3d(i, j, k-gamma).eq.1) then
               zn=gamma*sqrt(1.d0-x(i)*x(i)-y(j)*y(j))
               mueff=x(i)*nn_x + y(j)*nn_y + zn*nn_z
               if(mueff.lt.0.d0) then
                  into=0.d0
               else
                  into=xic1
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
            alocont_o_nn3d(i,j,k,14)=consta
            alocont_o_nn3d(i-alpha,j,k,aalpha)=alocont_o_nn3d(i-alpha,j,k,14)*constb
            alocont_o_nn3d(i,j-beta,k,bbeta)=alocont_o_nn3d(i,j-beta,k,14)*constc
            alocont_o_nn3d(i,j,k-gamma,ggamma) = alocont_o_nn3d(i,j,k-gamma,14)*constd
!
!-----------calculating boundary intensities on grid-points-------------
!-------------overwriting intensities calculated above------------------
!
            if(imask_bpoint3d(i,j,k).eq.1) then
            !calculate mueff
               mueff=nn_x*x(i) + nn_y*y(j) + nn_z*z(k)
               if(mueff.ge.0.d0) then
                  int3d(i,j,k)=  xic1
                  alocont_o_nn3d(i,j,k,:)= 0.d0
                  alocont_o_nn3d(i-alpha,j,k,aalpha)=0.d0
                  alocont_o_nn3d(i,j-beta,k,bbeta)=0.d0
                  alocont_o_nn3d(i,j,k-gamma,ggamma) =0.d0
               end if
            end if
!
!------------------------perform angular integration--------------------
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
            alocont_nn3d_tmp(i,j,k,14) = alocont_nn3d_tmp(i,j,k,14) + alocont_o_nn3d(i,j,k,14)*wall
            alocont_nn3d_tmp(i-alpha,j,k,aalpha) = alocont_nn3d_tmp(i-alpha,j,k,aalpha) + alocont_o_nn3d(i-alpha,j,k,aalpha)*wall
            alocont_nn3d_tmp(i,j-beta,k,bbeta) = alocont_nn3d_tmp(i,j-beta,k,bbeta) + alocont_o_nn3d(i,j-beta,k,bbeta)*wall
            alocont_nn3d_tmp(i,j,k-gamma,ggamma) = alocont_nn3d_tmp(i,j,k-gamma,ggamma) + alocont_o_nn3d(i,j,k-gamma,ggamma)*wall
            normalization3d_tmp(i,j,k) = normalization3d_tmp(i,j,k) + wall
!
         end if
!
      enddo
   enddo
enddo
!
!
end subroutine ffvm_cont3d
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine ffvm_line3d(oindx, xobsindx)
!
!-----------------------------------------------------------------------
!---------------------finite volume method in 2d------------------------
!-------calculating intensties for given mu, phi and frequency----------
!--------------specified by input oindx and xobsindx--------------------
!-----------------------------------------------------------------------
!
use prog_type
use fund_const
use dime3d, only: ndxmax, ndymax, ndzmax, x, y, z, delx_arr, dely_arr, delz_arr, &
                  int3d, velx3d, vely3d, velz3d, vth3d, scont3d, sline3d, &
                  opac3d, opalbar3d, aloline_on_nn3d, &
                  imask_totreg3d, imask_innreg3d, imask_bpoint3d, &
                  mintbar3d_tmp, normalization3d_tmp, aloline_nn3d_tmp
!                  mintbar3d, normalization3d, aloline_nn3d
use angles, only: n_x, n_y, n_z, weight_omega
use freq, only: nodes_xobs, weight_xobs, nxobs
use params_input, only: vth_fiducial
use bcondition, only: xic1
use omp_lib
!
implicit none
!
! ... arguments
integer(i4b), intent(in) :: oindx, xobsindx
!
! ... local scalars
integer(i4b) :: i, j, k
integer(i4b) :: alpha, beta, gamma, aalpha, bbeta, ggamma
integer(i4b) :: startx, starty, startz, endx, endy, endz
integer(i4b) :: posx_indx, posy_indx, posz_indx
real(dp) :: nn_x, nn_y, nn_z, mueff, xobs, wall, wall_global, phinorm, velz
real(dp) :: xn, yn, zn, intm, intn, into
real(dp) :: duma, dumb, dumc, dumd, dume
real(dp) :: consta, constb, constc, constd, conste, denominator
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
nn_x=n_x(oindx)
nn_y=n_y(oindx)
nn_z=n_z(oindx)
!
xobs=nodes_xobs(xobsindx)
!
wall_global=weight_omega(oindx)*weight_xobs(xobsindx)
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
!-----------------------------------------------------------------------
!
do k=startz, endz, gamma
   do j=starty, endy, beta
      do i=startx, endx, alpha
!
      if(imask_totreg3d(i,j,k).eq.1) then
!
!-----calculating boundary intensities on intermesh-points if needed----
!
         if(imask_innreg3d(i-alpha,j,k).eq.1) then
            !calculating intersection point with photosphere (only for spherical star)
            xn=alpha*sqrt(1.d0-y(j)*y(j)-z(k)*z(k))
            mueff=xn*nn_x + y(j)*nn_y + z(k)*nn_z
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
            !calculating intersection point with photosphere (only for spherical star)
            yn=beta*sqrt(1.d0-x(i)*x(i)-z(k)*z(k))
            mueff=x(i)*nn_x + yn*nn_y + z(k)*nn_z
            if(mueff.lt.0.d0) then
               intm = 0.d0
            else
               intm = xic1
            endif
         else
            intm=int3d(i,j-beta, k)
         end if
!
         if(imask_innreg3d(i,j,k-gamma).eq.1) then
            !calculating intersection point with photosphere (only for spherical star)
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
!-------------------calculate normalized profile function---------------
!
         velz = velx3d(i,j,k)*nn_x + vely3d(i,j,k)*nn_y + velz3d(i,j,k)*nn_z
         call calc_phinorm(velz, vth3d(i,j,k), vth_fiducial, xobs, phinorm)
!
!no refinement required
         duma=opac3d(i,j,k)
         dumb=opalbar3d(i,j,k)*phinorm
         dumc=nn_x/delx_arr(i,j,k,posx_indx)
         dumd=nn_y/dely_arr(i,j,k,posy_indx)
         dume=nn_z/delz_arr(i,j,k,posz_indx)

         denominator=duma+dumb+dumc+dumd+dume

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
          aloline_on_nn3d(i-alpha,j,k,aalpha)=aloline_on_nn3d(i-alpha,j,k,14)*constc
          aloline_on_nn3d(i,j-beta,k,bbeta)=aloline_on_nn3d(i,j-beta,k,14)*constd
          aloline_on_nn3d(i,j,k-gamma,ggamma) = aloline_on_nn3d(i,j,k-gamma,14)*conste
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
               aloline_on_nn3d(i-alpha,j,k,aalpha)=0.d0
               aloline_on_nn3d(i,j-beta,k,bbeta)=0.d0
               aloline_on_nn3d(i,j,k-gamma,ggamma) =0.d0
            end if
         end if
!
!----------------perform angular and frequency integration--------------
!
         wall=wall_global*phinorm
         mintbar3d_tmp(i,j,k) = mintbar3d_tmp(i,j,k) + int3d(i,j,k)*wall
         aloline_nn3d_tmp(i,j,k,14) = aloline_nn3d_tmp(i,j,k,14) + aloline_on_nn3d(i,j,k,14) * wall
         aloline_nn3d_tmp(i-alpha,j,k,aalpha) = aloline_nn3d_tmp(i-alpha,j,k,aalpha) + aloline_on_nn3d(i-alpha,j,k,aalpha) * wall
         aloline_nn3d_tmp(i,j-beta,k,bbeta) = aloline_nn3d_tmp(i,j-beta,k,bbeta) + aloline_on_nn3d(i,j-beta,k,bbeta) * wall
         aloline_nn3d_tmp(i,j,k-gamma,ggamma) = aloline_nn3d_tmp(i,j,k-gamma,ggamma) + aloline_on_nn3d(i,j,k-gamma,ggamma) * wall
         normalization3d_tmp(i,j,k) = normalization3d_tmp(i,j,k) + wall

!         mintbar3d(i,j,k) = mintbar3d(i,j,k) + int3d(i,j,k)*wall
!         aloline_nn3d(i,j,k,14) = aloline_nn3d(i,j,k,14) + aloline_on_nn3d(i,j,k,14) * wall
!         aloline_nn3d(i-alpha,j,k,aalpha) = aloline_nn3d(i-alpha,j,k,aalpha) + aloline_on_nn3d(i-alpha,j,k,aalpha) * wall
!         aloline_nn3d(i,j-beta,k,bbeta) = aloline_nn3d(i,j-beta,k,bbeta) + aloline_on_nn3d(i,j-beta,k,bbeta) * wall
!         aloline_nn3d(i,j,k-gamma,ggamma) = aloline_nn3d(i,j,k-gamma,ggamma) + aloline_on_nn3d(i,j,k-gamma,ggamma) * wall
!         normalization3d(i,j,k) = normalization3d(i,j,k) + wall
!
      end if
!
      enddo
   enddo
enddo
!
!
!
end subroutine ffvm_line3d
