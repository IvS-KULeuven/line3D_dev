subroutine benchmark06_solution
!
!-----calculating line transport along a direction for a single cell----
!--------testing different velocity laws and interpolation schemes------
!--------testing refinement of the spatial grid along ray---------------
!
use prog_type
use mod_benchmark, only: im_vel, xobs1, xobs2, xobs3
use mod_benchmark, only: s1d_method01, tau1d_method01, vel1d_method01, vth1d_method01, &
                         xcmf1d_method01, profile1d_method01, opalbar1d_method01, &
                         opal1d_method01, sline1d_method01, int1d_method01, nd_method01, &
                         s1d_method02, tau1d_method02, vel1d_method02, vth1d_method02, &
                         xcmf1d_method02, profile1d_method02, opalbar1d_method02, &
                         opal1d_method02, sline1d_method02, int1d_method02, nd_method02, &
                         s1d_method03, tau1d_method03, vel1d_method03, vth1d_method03, &
                         xcmf1d_method03, profile1d_method03, opalbar1d_method03, &
                         opal1d_method03, sline1d_method03, int1d_method03, nd_method03, &
                         s1d_method04, tau1d_method04, vel1d_method04, vth1d_method04, &
                         xcmf1d_method04, profile1d_method04, opalbar1d_method04, &
                         opal1d_method04, sline1d_method04, int1d_method04, nd_method04, &
                         s1d_fine, tau1d_fine, vel1d_fine, vth1d_fine, &
                         xcmf1d_fine, profile1d_fine, opalbar1d_fine, &
                         opal1d_fine, sline1d_fine, int1d_fine, nd_fine, &
                         s1dc_method01, tau1dc_method01, vel1dc_method01, vth1dc_method01, &
                         xcmf1dc_method01, profile1dc_method01, opalbar1dc_method01, &
                         opal1dc_method01, sline1dc_method01, int1dc_method01, ndc_method01, &
                         opac1dc_method01, scont1dc_method01, stot1dc_method01, opatot1dc_method01, &
                         s1dc_method02, tau1dc_method02, vel1dc_method02, vth1dc_method02, &
                         xcmf1dc_method02, profile1dc_method02, opalbar1dc_method02, &
                         opal1dc_method02, sline1dc_method02, int1dc_method02, ndc_method02, &
                         opac1dc_method02, scont1dc_method02, stot1dc_method02, opatot1dc_method02, &
                         s1dc_method03, tau1dc_method03, vel1dc_method03, vth1dc_method03, &
                         xcmf1dc_method03, profile1dc_method03, opalbar1dc_method03, &
                         opal1dc_method03, sline1dc_method03, int1dc_method03, ndc_method03, &
                         opac1dc_method03, scont1dc_method03, stot1dc_method03, opatot1dc_method03, &
                         s1dc_method04, tau1dc_method04, vel1dc_method04, vth1dc_method04, &
                         xcmf1dc_method04, profile1dc_method04, opalbar1dc_method04, &
                         opal1dc_method04, sline1dc_method04, int1dc_method04, ndc_method04, &
                         opac1dc_method04, scont1dc_method04, stot1dc_method04, opatot1dc_method04, &
                         s1dc_fine, tau1dc_fine, vel1dc_fine, vth1dc_fine, &
                         xcmf1dc_fine, profile1dc_fine, opalbar1dc_fine, &
                         opal1dc_fine, sline1dc_fine, int1dc_fine, ndc_fine, &
                         opac1dc_fine, scont1dc_fine, stot1dc_fine, opatot1dc_fine
use mod_interp1d, only: interpol_yp, interpol_ypl, interpol_typ_quad3
!
implicit none
!
!interface required to allocate arrays in fsc_line1d_debug
interface
   subroutine fsc_line_debug(xobs, vth_fiducial, deltax, xcmf_min, xcmf_max, int_u, opalbar_u, opalbar_p, opalbar_d, &
                         sline_u, sline_p, sline_d, vel_u, vel_p, vel_d, &
                         vth_u, vth_p, vth_d, dels_u, dels_d, abs_sc, contr_sc, int_sc, alo_u, alo_p, alo_d, luse_refine, &
                         xcmf1d, vel1d, opalbar1d, opal1d, vth1d, sline1d, profile1d, int1d, s1d, tau1d, nd)
      use prog_type
      logical, intent(in) :: luse_refine
      real(dp), intent(in) :: xobs, deltax, xcmf_min, xcmf_max, vth_fiducial, int_u, opalbar_p, sline_p, vel_p, vth_p, &
                              opalbar_u, sline_u, vel_u, vth_u, dels_u, &
                              opalbar_d, sline_d, vel_d, vth_d, dels_d
      real(dp), intent(out) :: alo_u, alo_p, alo_d
      integer(i4b), intent(inout) :: nd
      real(dp), intent(out) :: abs_sc, contr_sc, int_sc
      real(dp), dimension(:), allocatable, intent(inout) :: xcmf1d, vel1d, s1d, tau1d, opalbar1d, sline1d, &
                                                            vth1d, opal1d, int1d, profile1d
   end subroutine fsc_line_debug
!
   subroutine fsc_line_lin_debug(xobs, vth_fiducial, deltax, xcmf_min, xcmf_max, int_u, &
                                   opalbar_u, opalbar_p, &
                                   sline_u, sline_p, &
                                   vel_u, vel_p, &
                                   vth_u, vth_p, &
                                   dels_u, &
                                   abs_sc, contr_sc, int_sc, alo_u, alo_p, luse_refine, &
                                   xcmf1d, vel1d, opalbar1d, opal1d, vth1d, sline1d, profile1d, int1d, s1d, tau1d, nd)
      use prog_type
      logical, intent(in) :: luse_refine
      real(dp), intent(in) :: xobs, deltax, xcmf_min, xcmf_max, vth_fiducial, int_u, opalbar_p, sline_p, vel_p, vth_p, &
                              opalbar_u, sline_u, vel_u, vth_u, dels_u
      real(dp), intent(out) :: abs_sc, contr_sc, int_sc, alo_u, alo_p
      integer(i4b), intent(out) :: nd
      real(dp), dimension(:), allocatable, intent(inout) :: xcmf1d, vel1d, s1d, tau1d, opalbar1d, sline1d, &
                                                            vth1d, opal1d, int1d, profile1d
   end subroutine fsc_line_lin_debug
!
   subroutine fsc_linec_lin_debug(xobs, vth_fiducial, deltax, xcmf_min, xcmf_max, int_u, &
                                    opac_u, opac_p, &
                                    opalbar_u, opalbar_p, &
                                    scont_u, scont_p, &
                                    sline_u, sline_p, &
                                    vel_u, vel_p, &
                                    vth_u, vth_p, &
                                    dels_u, &
                                    abs_sc, contr_sc, int_sc, alo_u, alo_p, luse_refine, &
                                    xcmf1d, vel1d, opac1d, opalbar1d, opal1d, opatot1d, &
                                    vth1d, scont1d, sline1d, stot1d, profile1d, int1d, s1d, tau1d, nd)
      use prog_type
      logical, intent(in) :: luse_refine
      integer, intent(out) :: nd
      real(dp), intent(in) :: xobs, vth_fiducial, deltax, xcmf_min, xcmf_max, int_u, dels_u, &
                              opac_u, opalbar_u, scont_u, sline_u, vel_u, vth_u, & 
                              opac_p, opalbar_p, scont_p, sline_p, vel_p, vth_p
      real(dp), intent(out) :: abs_sc, contr_sc, int_sc, alo_u, alo_p
      real(dp), dimension(:), allocatable, intent(inout) :: xcmf1d, vel1d, s1d, tau1d, opalbar1d, opac1d, opatot1d, opal1d, &
                                                            sline1d, scont1d, stot1d, vth1d, int1d, profile1d
   end subroutine fsc_linec_lin_debug
!
   subroutine fsc_linec_debug(xobs, vth_fiducial, deltax, xcmf_min, xcmf_max, int_u, &
                                opac_u, opac_p, opac_d, &
                                opalbar_u, opalbar_p, opalbar_d, &
                                scont_u, scont_p, scont_d, &
                                sline_u, sline_p, sline_d, &
                                vel_u, vel_p, vel_d, &
                                vth_u, vth_p, vth_d, &
                                dels_u, dels_d, &
                                abs_sc, contr_sc, int_sc, alo_u, alo_p, alo_d, luse_refine, &
                                xcmf1d, vel1d, opac1d, opalbar1d, opal1d, opatot1d, &
                                vth1d, scont1d, sline1d, stot1d, profile1d, int1d, s1d, tau1d, nd)
      use prog_type
      logical, intent(in) :: luse_refine
      integer, intent(out) :: nd
      real(dp), intent(in) :: xobs, vth_fiducial, deltax, xcmf_min, xcmf_max, int_u, dels_u, dels_d, &
                              opac_u, opalbar_u, scont_u, sline_u, vel_u, vth_u, & 
                              opac_p, opalbar_p, scont_p, sline_p, vel_p, vth_p, &
                              opac_d, opalbar_d, scont_d, sline_d, vel_d, vth_d
      real(dp), intent(out) :: abs_sc, contr_sc, int_sc, alo_u, alo_p, alo_d
      real(dp), dimension(:), allocatable, intent(inout) :: xcmf1d, vel1d, s1d, tau1d, opalbar1d, opac1d, opatot1d, opal1d, &
                                                            sline1d, scont1d, stot1d, vth1d, int1d, profile1d
   end subroutine fsc_linec_debug
!
end interface
!
! ... local scalars
integer(i4b) :: i, nd_dum
real(dp) :: vth_fiducial, deltax, xcmf_min, xcmf_max, vmin, vmax, opalbar_min, opalbar_max, sline_min, sline_max, &
            opac_min, opac_max, scont_min, scont_max, &
            s_u, s_p, s_d, vel_u, vel_p, vel_d, opac_u, opac_p, opac_d, opalbar_u, opalbar_p, opalbar_d, &
            scont_u, scont_p, scont_d, sline_u, sline_p, sline_d, vth_u, vth_p, vth_d, dels_u, dels_d, int_u
real(dp) :: abs_sc, contr_sc, int_sc, int_sc_sspace, dtau, alo_u, alo_p, alo_d
!
! ... local arrays
real(dp), dimension(:), allocatable :: s1d_dum, tau1d_dum, vel1d_dum, vth1d_dum, &
                                       xcmf1d_dum, profile1d_dum, opalbar1d_dum, &
                                       opal1d_dum, sline1d_dum, int1d_dum, &
                                       opac1d_dum, opatot1d_dum, scont1d_dum, stot1d_dum
!
! ... local logicals
logical :: luse_refine
!
! ... local functions
real(dp) :: integral3, vlaw, slinelaw, opalbarlaw, opaclaw, scontlaw
!
!---------------------------define model--------------------------------
!
deltax = 1.d0/3.d0
xcmf_min =-3.d0
xcmf_max = 3.d0
!
!
!
s_u=0.d0
s_p=1.d0
s_d=2.d0
dels_u=s_p-s_u
dels_d=s_d-s_p
!
!opacities
opalbar_min=0.d0
opalbar_max=20.d0
opalbar_u=opalbarlaw(opalbar_min,opalbar_max,s_u,s_p,s_u)
opalbar_p=opalbarlaw(opalbar_min,opalbar_max,s_u,s_p,s_p)
opalbar_d=opalbarlaw(opalbar_min,opalbar_max,s_u,s_p,s_d)
!
opac_min=0.1d0
opac_max=1.d0
opac_u=opaclaw(opac_min,opac_max,s_u,s_p,s_u)
opac_p=opaclaw(opac_min,opac_max,s_u,s_p,s_p)
opac_d=opaclaw(opac_min,opac_max,s_u,s_p,s_d)
!
!line and continuum source functions
sline_min=0.d0
sline_max=0.5d0
sline_u=slinelaw(sline_min,sline_max,s_u,s_p,s_u)
sline_p=slinelaw(sline_min,sline_max,s_u,s_p,s_p)
sline_d=slinelaw(sline_min,sline_max,s_u,s_p,s_d)
!
scont_min=0.d0
scont_max=0.5d0
scont_u=scontlaw(scont_min,scont_max,s_u,s_p,s_u)
scont_p=scontlaw(scont_min,scont_max,s_u,s_p,s_p)
scont_d=scontlaw(scont_min,scont_max,s_u,s_p,s_d)
!
!thermal velocities
vth_fiducial=1.d0
vth_u=vth_fiducial
vth_p=vth_fiducial
vth_d=vth_fiducial
!
vmin=0.d0
vmax=12.d0
vel_u=vlaw(vmin,vmax,s_u,s_p,s_u,im_vel)
vel_p=vlaw(vmin,vmax,s_u,s_p,s_p,im_vel)
vel_d=vlaw(vmin,vmax,s_u,s_p,s_d,im_vel)
!
!for resonance region right in the center of u and p, and exactly at u, p, respectively
xobs1=(vel_u+vel_p)/2.d0
!xobs1=vel_u
!xobs1=vel_p
!
!********************solutions without continuum************************
!
!--------------method00: fine grid and linear approach------------------
!----------------(gives roughly the exact solution)---------------------
!
nd_fine=1000
allocate(vel1d_fine(nd_fine))
allocate(s1d_fine(nd_fine))
allocate(opalbar1d_fine(nd_fine))
allocate(opal1d_fine(nd_fine))
allocate(vth1d_fine(nd_fine))
allocate(sline1d_fine(nd_fine))
allocate(int1d_fine(nd_fine))
allocate(profile1d_fine(nd_fine))
allocate(tau1d_fine(nd_fine))
allocate(xcmf1d_fine(nd_fine))
!
do i=1, nd_fine
   s1d_fine(i)=s_u+(i-1)*(s_p-s_u)/(nd_fine-1)
enddo
!
!velocity law defined in funcion vlaw
!linear thermal velocity law
!linear opacity law
!quadratic line source functions
do i=1, nd_fine
   vel1d_fine(i)=vlaw(vmin,vmax,s_u,s_p,s1d_fine(i),im_vel)
   vth1d_fine(i)=interpol_yp(s_p,s_u,vth_p,vth_u,s1d_fine(i))
   opalbar1d_fine(i)=opalbarlaw(opalbar_min,opalbar_max,s_u,s_p,s1d_fine(i))
   sline1d_fine(i)=slinelaw(sline_min,sline_max,s_u,s_p,s1d_fine(i))
enddo
!
int_u=1.d0
!
int1d_fine(1)=int_u
tau1d_fine(1)=0.d0
do i=2, nd_fine
   call fsc_line_lin_debug(xobs1, vth_fiducial, deltax, xcmf_min, xcmf_max, int1d_fine(i-1), &
                             opalbar1d_fine(i-1), opalbar1d_fine(i), &
                             sline1d_fine(i-1), sline1d_fine(i), &
                             vel1d_fine(i-1), vel1d_fine(i), &
                             vth1d_fine(i-1), vth1d_fine(i), &
                             abs(s1d_fine(i)-s1d_fine(i-1)), abs_sc, contr_sc, int_sc, alo_u, alo_p, .false., &
                             xcmf1d_dum, vel1d_dum, opalbar1d_dum, opal1d_dum, &
                             vth1d_dum, sline1d_dum, profile1d_dum, int1d_dum, &
                             s1d_dum, tau1d_dum, nd_dum)
   if(i.eq.2) then
      int1d_fine(i-1) = int1d_dum(1)
      xcmf1d_fine(i-1) = xcmf1d_dum(1)
      opal1d_fine(i-1) = opal1d_dum(1)
      profile1d_fine(i-1) = profile1d_dum(1)
   endif
   int1d_fine(i) = int1d_dum(2)
   xcmf1d_fine(i) = xcmf1d_dum(2)
   opal1d_fine(i) = opal1d_dum(2)
   profile1d_fine(i) = profile1d_dum(2)
   tau1d_fine(i) = tau1d_fine(i-1)+tau1d_dum(2)
enddo
!
!-----------------------------method01----------------------------------
!--------linear approach without continuum and without refinement-------
!
luse_refine=.false.
call fsc_line_lin_debug(xobs1, vth_fiducial, deltax, xcmf_min, xcmf_max, int_u, &
                          opalbar_u, opalbar_p, &
                          sline_u, sline_p, &
                          vel_u, vel_p, &
                          vth_u, vth_p, &
                          dels_u, abs_sc, contr_sc, int_sc, alo_u, alo_p, luse_refine, &
                          xcmf1d_method01, vel1d_method01, opalbar1d_method01, opal1d_method01, &
                          vth1d_method01, sline1d_method01, profile1d_method01, int1d_method01, &
                          s1d_method01, tau1d_method01, nd_method01)
!
!-----------------------------method02-----------------------------------
!--------linear approach without continuum and including refinement------
!
luse_refine=.true.
call fsc_line_lin_debug(xobs1, vth_fiducial, deltax, xcmf_min, xcmf_max, int_u, &
                          opalbar_u, opalbar_p, &
                          sline_u, sline_p, &
                          vel_u, vel_p, &
                          vth_u, vth_p, &
                          dels_u, abs_sc, contr_sc, int_sc, alo_u, alo_p, luse_refine, &
                          xcmf1d_method02, vel1d_method02, opalbar1d_method02, opal1d_method02, &
                          vth1d_method02, sline1d_method02, profile1d_method02, int1d_method02, &
                          s1d_method02, tau1d_method02, nd_method02)

call fsc_line_lin(xobs1, vth_fiducial, deltax, xcmf_min, xcmf_max, int_u, &
                    opalbar_u, opalbar_p, &
                    sline_u, sline_p, &
                    vel_u, vel_p, &
                    vth_u, vth_p, &
                    dels_u, &
                    abs_sc, contr_sc, int_sc, alo_u, alo_p)
!
write(*,'(a40,es20.8)') 'debugging routine lin (no refinement)', int1d_method01(2)
write(*,'(a40,es20.8)') 'debugging routine lin (refinement)', int1d_method02(nd_method02)
write(*,'(a40,es20.8)') 'normal routine lin', int_sc
write(*,*)
if(abs(int1d_method02(nd_method02)-int_sc).gt.1.d-14) stop 'error in benchmark06: fsc_line1d_lin_debug and fsc_line1d_lin give different solutions'
!
!-----------------------------method03-----------------------------------
!--------bezier approach without continuum and without refinement--------
!
luse_refine=.false.
call fsc_line_debug(xobs1, vth_fiducial, deltax, xcmf_min, xcmf_max, int_u, &
                      opalbar_u, opalbar_p, opalbar_d, &
                      sline_u, sline_p, sline_d, &
                      vel_u, vel_p, vel_d, &
                      vth_u, vth_p, vth_d, &
                      dels_u, dels_d, abs_sc, contr_sc, int_sc, alo_u, alo_p, alo_d, luse_refine, &
                      xcmf1d_method03, vel1d_method03, opalbar1d_method03, opal1d_method03, &
                      vth1d_method03, sline1d_method03, profile1d_method03, int1d_method03, &
                      s1d_method03, tau1d_method03, nd_method03)
!
!-----------------------------method04-----------------------------------
!--------bezier approach without continuum and including refinement------
!
luse_refine=.true.
call fsc_line_debug(xobs1, vth_fiducial, deltax, xcmf_min, xcmf_max, int_u, &
                      opalbar_u, opalbar_p, opalbar_d, &
                      sline_u, sline_p, sline_d, &
                      vel_u, vel_p, vel_d, &
                      vth_u, vth_p, vth_d, &
                      dels_u, dels_d, abs_sc, contr_sc, int_sc, alo_u, alo_p, alo_d, luse_refine, &
                      xcmf1d_method04, vel1d_method04, opalbar1d_method04, opal1d_method04, &
                      vth1d_method04, sline1d_method04, profile1d_method04, int1d_method04, &
                      s1d_method04, tau1d_method04, nd_method04)
!
call fsc_line(xobs1, vth_fiducial, deltax, xcmf_min, xcmf_max, int_u, &
                opalbar_u, opalbar_p, opalbar_d, &
                sline_u, sline_p, sline_d, &
                vel_u, vel_p, vel_d, &
                vth_u, vth_p, vth_d, &
                dels_u, dels_d, abs_sc, contr_sc, int_sc, alo_u, alo_p, alo_d)

!
write(*,'(a40,es20.8)') 'debugging routine bez (no refinement)', int1d_method03(2)
write(*,'(a40,es20.8)') 'debugging routine bez (refinement)', int1d_method04(nd_method04)
write(*,'(a40,es20.8)') 'normal routine bez ', int_sc
write(*,*)
if(abs(int1d_method04(nd_method04)-int_sc).gt.1.d-14) stop 'error in benchmark06: fsc_line1d_debug and fsc_line1d give different solutions'
!
!
!********************solutions including continuum**********************
!
!--------------method00c: fine grid and linear approach-----------------
!----------------(gives roughly the exact solution)---------------------
!
ndc_fine=1000
allocate(vel1dc_fine(ndc_fine))
allocate(s1dc_fine(ndc_fine))
allocate(opalbar1dc_fine(ndc_fine))
allocate(opal1dc_fine(ndc_fine))
allocate(opac1dc_fine(ndc_fine))
allocate(opatot1dc_fine(ndc_fine))
allocate(vth1dc_fine(ndc_fine))
allocate(sline1dc_fine(ndc_fine))
allocate(scont1dc_fine(ndc_fine))
allocate(stot1dc_fine(ndc_fine))
allocate(int1dc_fine(ndc_fine))
allocate(profile1dc_fine(ndc_fine))
allocate(tau1dc_fine(ndc_fine))
allocate(xcmf1dc_fine(ndc_fine))
!
do i=1, nd_fine
   s1dc_fine(i)=s_u+(i-1)*(s_p-s_u)/(nd_fine-1)
enddo
!
!velocity law defined in funcion vlaw
!linear thermal velocity law
!linear opacity law
!quadratic line source functions
do i=1, nd_fine
   vel1dc_fine(i)=vlaw(vmin,vmax,s_u,s_p,s1dc_fine(i),im_vel)
   vth1dc_fine(i)=interpol_yp(s_p,s_u,vth_p,vth_u,s1dc_fine(i))
   opalbar1dc_fine(i)=opalbarlaw(opalbar_min,opalbar_max,s_u,s_p,s1dc_fine(i))
   opac1dc_fine(i)=opaclaw(opac_min,opac_max,s_u,s_p,s1dc_fine(i))
   sline1dc_fine(i)=slinelaw(sline_min,sline_max,s_u,s_p,s1dc_fine(i))
   scont1dc_fine(i)=scontlaw(scont_min,scont_max,s_u,s_p,s1dc_fine(i))
enddo
!
int_u=1.d0
!
int1dc_fine(1)=int_u
tau1dc_fine(1)=0.d0
do i=2, ndc_fine
   call fsc_linec_lin_debug(xobs1, vth_fiducial, deltax, xcmf_min, xcmf_max, int1dc_fine(i-1), &
                              opac1dc_fine(i-1), opac1dc_fine(i), &
                              opalbar1dc_fine(i-1), opalbar1dc_fine(i), &
                              scont1dc_fine(i-1), scont1dc_fine(i), &
                              sline1dc_fine(i-1), sline1dc_fine(i), &
                              vel1dc_fine(i-1), vel1dc_fine(i), &
                              vth1dc_fine(i-1), vth1dc_fine(i), &
                              abs(s1dc_fine(i)-s1dc_fine(i-1)), abs_sc, contr_sc, int_sc, alo_u, alo_p, .false., &
                              xcmf1d_dum, vel1d_dum, opac1d_dum, opalbar1d_dum, opal1d_dum, opatot1d_dum, &
                              vth1d_dum, scont1d_dum, sline1d_dum, stot1d_dum, profile1d_dum, int1d_dum, &
                              s1d_dum, tau1d_dum, nd_dum)
   if(i.eq.2) then
      int1dc_fine(i-1) = int1d_dum(1)
      xcmf1dc_fine(i-1) = xcmf1d_dum(1)
      opal1dc_fine(i-1) = opal1d_dum(1)
      profile1dc_fine(i-1) = profile1d_dum(1)
      opatot1dc_fine(i-1) = opatot1d_dum(1)
      stot1dc_fine(i-1) = stot1d_dum(1)
   endif
   int1dc_fine(i) = int1d_dum(2)
   xcmf1dc_fine(i) = xcmf1d_dum(2)
   opal1dc_fine(i) = opal1d_dum(2)
   profile1dc_fine(i) = profile1d_dum(2)
   tau1dc_fine(i) = tau1dc_fine(i-1)+tau1d_dum(2)
   opatot1dc_fine(i) = opatot1d_dum(2)
   stot1dc_fine(i) = stot1d_dum(2)
enddo
!
!-----------------------------method01c---------------------------------
!--------linear approach including continuum and without refinement-----
!
luse_refine=.false.

call fsc_linec_lin_debug(xobs1, vth_fiducial, deltax, xcmf_min, xcmf_max, int_u, &
                           opac_u, opac_p, &
                           opalbar_u, opalbar_p, &
                           scont_u, scont_p, &
                           sline_u, sline_p, &
                           vel_u, vel_p, &
                           vth_u, vth_p, &
                           dels_u, abs_sc, contr_sc, int_sc, alo_u, alo_p, luse_refine, &
                           xcmf1dc_method01, vel1dc_method01, opac1dc_method01, opalbar1dc_method01, opal1dc_method01, opatot1dc_method01, &
                           vth1dc_method01, scont1dc_method01, sline1dc_method01, stot1dc_method01, profile1dc_method01, int1dc_method01, &
                           s1dc_method01, tau1dc_method01, ndc_method01)
!
!-----------------------------method02c----------------------------------
!--------linear approach including continuum and including refinement----
!
luse_refine=.true.

call fsc_linec_lin_debug(xobs1, vth_fiducial, deltax, xcmf_min, xcmf_max, int_u, &
                           opac_u, opac_p, &
                           opalbar_u, opalbar_p, &
                           scont_u, scont_p, &
                           sline_u, sline_p, &
                           vel_u, vel_p, &
                           vth_u, vth_p, &
                           dels_u, abs_sc, contr_sc, int_sc, alo_u, alo_p, luse_refine, &
                           xcmf1dc_method02, vel1dc_method02, opac1dc_method02, opalbar1dc_method02, opal1dc_method02, opatot1dc_method02, &
                           vth1dc_method02, scont1dc_method02, sline1dc_method02, stot1dc_method02, profile1dc_method02, int1dc_method02, &
                           s1dc_method02, tau1dc_method02, ndc_method02)

call fsc_linec_lin(xobs1, vth_fiducial, deltax, xcmf_min, xcmf_max, int_u, &
                     opac_u, opac_p, &
                     opalbar_u, opalbar_p, &
                     scont_u, scont_p, &
                     sline_u, sline_p, &
                     vel_u, vel_p, &
                     vth_u, vth_p, &
                     dels_u, abs_sc, contr_sc, int_sc, alo_u, alo_p)
!
write(*,'(a40,es20.8)') 'debugging routine lin (no refinement)', int1dc_method01(2)
write(*,'(a40,es20.8)') 'debugging routine lin (refinement)', int1dc_method02(ndc_method02)
write(*,'(a40,es20.8)') 'normal routine lin', int_sc
write(*,*)
if(abs(int1dc_method02(ndc_method02)-int_sc).gt.1.d-14) stop 'error in benchmark06: fsc_linec1d_lin_debug and fsc_linec1d_lin give different solutions'
!
!-----------------------------method03c----------------------------------
!--------bezier approach including continuum and without refinement------
!
luse_refine=.false.
call fsc_linec_debug(xobs1, vth_fiducial, deltax, xcmf_min, xcmf_max, int_u, &
                       opac_u, opac_p, opac_d, &
                       opalbar_u, opalbar_p, opalbar_d, &
                       scont_u, scont_p, scont_d, &
                       sline_u, sline_p, sline_d, &
                       vel_u, vel_p, vel_d, &
                       vth_u, vth_p, vth_d, &
                       dels_u, dels_d, abs_sc, contr_sc, int_sc, alo_u, alo_p, alo_d, luse_refine, &
                       xcmf1dc_method03, vel1dc_method03, opac1dc_method03, opalbar1dc_method03, opal1dc_method03, opatot1dc_method03, &
                       vth1dc_method03, scont1dc_method03, sline1dc_method03, stot1dc_method03, profile1dc_method03, int1dc_method03, &
                       s1dc_method03, tau1dc_method03, ndc_method03)
!
!-----------------------------method04c----------------------------------
!------bezier approach including continuum and including refinement------
!
luse_refine=.true.
call fsc_linec_debug(xobs1, vth_fiducial, deltax, xcmf_min, xcmf_max, int_u, &
                       opac_u, opac_p, opac_d, &
                       opalbar_u, opalbar_p, opalbar_d, &
                       scont_u, scont_p, scont_d, &
                       sline_u, sline_p, sline_d, &
                       vel_u, vel_p, vel_d, &
                       vth_u, vth_p, vth_d, &
                       dels_u, dels_d, abs_sc, contr_sc, int_sc, alo_u, alo_p, alo_d, luse_refine, &
                       xcmf1dc_method04, vel1dc_method04, opac1dc_method04, opalbar1dc_method04, opal1dc_method04, opatot1dc_method04, &
                       vth1dc_method04, scont1dc_method04, sline1dc_method04, stot1dc_method04, profile1dc_method04, int1dc_method04, &
                       s1dc_method04, tau1dc_method04, ndc_method04)
!
call fsc_linec(xobs1, vth_fiducial, deltax, xcmf_min, xcmf_max, int_u, &
                 opac_u, opac_p, opac_d, &
                 opalbar_u, opalbar_p, opalbar_d, &
                 scont_u, scont_p, scont_d, &
                 sline_u, sline_p, sline_d, &
                 vel_u, vel_p, vel_d, &
                 vth_u, vth_p, vth_d, &
                 dels_u, dels_d, &
                 abs_sc, contr_sc, int_sc, alo_u, alo_p, alo_d)
!
!
write(*,'(a40,es20.8)') 'debugging routine bez (no refinement)', int1dc_method03(2)
write(*,'(a40,es20.8)') 'debugging routine bez (refinement)', int1dc_method04(ndc_method04)
write(*,'(a40,es20.8)') 'normal routine bez ', int_sc
write(*,*)
if(abs(int1dc_method04(ndc_method04)-int_sc).gt.1.d-14) stop 'error in benchmark06: fsc_linec1d_debug and fsc_linec1d give different solutions'
!
!
write(*,'(3a30)') 'method', 'solution', 'error'
write(*,*) 'without continuum'
write(*,'(a30,es30.8)') 'fine grid (lin)',  int1d_fine(nd_fine)
write(*,'(a30,2es30.8)') 'without refinement (lin)', int1d_method01(nd_method01), abs(int1d_fine(nd_fine)-int1d_method01(nd_method01))/int1d_fine(nd_fine)
write(*,'(a30,2es30.8)') 'with    refinement (lin)', int1d_method02(nd_method02), abs(int1d_fine(nd_fine)-int1d_method02(nd_method02))/int1d_fine(nd_fine)
write(*,'(a30,2es30.8)') 'without refinement (bez)', int1d_method03(nd_method03), abs(int1d_fine(nd_fine)-int1d_method03(nd_method03))/int1d_fine(nd_fine)
write(*,'(a30,2es30.8)') 'with    refinement (bez)', int1d_method04(nd_method04), abs(int1d_fine(nd_fine)-int1d_method04(nd_method04))/int1d_fine(nd_fine)
write(*,*) 'with    continuum'
write(*,'(a30,es30.8)') 'fine grid (lin)',  int1dc_fine(ndc_fine)
write(*,'(a30,2es30.8)') 'without refinement (lin)', int1dc_method01(ndc_method01), abs(int1dc_fine(ndc_fine)-int1dc_method01(ndc_method01))/int1dc_fine(ndc_fine)
write(*,'(a30,2es30.8)') 'with    refinement (lin)', int1dc_method02(ndc_method02), abs(int1dc_fine(ndc_fine)-int1dc_method02(ndc_method02))/int1dc_fine(ndc_fine)
write(*,'(a30,2es30.8)') 'without refinement (bez)', int1dc_method03(ndc_method03), abs(int1dc_fine(ndc_fine)-int1dc_method03(ndc_method03))/int1dc_fine(ndc_fine)
write(*,'(a30,2es30.8)') 'with    refinement (bez)', int1dc_method04(ndc_method04), abs(int1dc_fine(ndc_fine)-int1dc_method04(ndc_method04))/int1dc_fine(ndc_fine)
!
!-----------------------------------------------------------------------
!
!int_u = 5.689863897292737d-3
!vth_fiducial = 10033425.4684825d0
!
!scont_u = 4.869925881736490d-3
!scont_p = 5.932179089584421d-3 
!scont_d = 1.379710699999998d-2
!
!sline_u = 5.231458618546425d-3
!sline_p = 6.867156543779988d-3
!sline_d = 1.379710699999998d-2
!
!opac_u = 4.145138673836960d-8
!opac_p = 1.648993046375141d-8
!opac_d = 1.243915534735381d-8
!
!opalbar_u = 3183362256.98250d0
!opalbar_p = 3163381.69546548d0
!opalbar_d = 751503.367231688d0
!
!vel_u = -1.81176631637354d0
!vel_p = -0.468625757253121d0
!vel_d = 0.d0
!
!vth_u = 10033425.4684825d0
!vth_p = 10033425.4684825d0
!vth_d = 10033425.4684825d0   
!
!dels_u = 0.148716481171059d0
!dels_d = 0.146808830968520d0   
!
!xobs1 = -2.97983350771677d0   
!
!call fsc_linec(xobs1, vth_fiducial, deltax, xcmf_min, xcmf_max, int_u, &
!                 opac_u, opac_p, opac_d, &
!                 opalbar_u, opalbar_p, opalbar_d, &
!                 scont_u, scont_p, scont_d, &
!                 sline_u, sline_p, sline_d, &
!                 vel_u, vel_p, vel_d, &
!                 vth_u, vth_p, vth_d, &
!                 dels_u, dels_d, &
!                 abs_sc, contr_sc, int_sc, alo_u, alo_p, alo_d)
!write(*,*) int_sc, abs_sc, contr_sc
!write(*,*) alo_u, alo_p, alo_d
!
!call fsc_line(xobs1, vth_fiducial, deltax, xcmf_min, xcmf_max, int_u, &
!                opalbar_u, opalbar_p, opalbar_d, &
!                sline_u, sline_p, sline_d, &
!                vel_u, vel_p, vel_d, &
!                vth_u, vth_p, vth_d, &
!                dels_u, dels_d, abs_sc, contr_sc, int_sc, alo_u, alo_p, alo_d)
!write(*,*) int_sc, abs_sc, contr_sc
!write(*,*) alo_u, alo_p, alo_d
!stop
!
!
end subroutine benchmark06_solution
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
function vlaw(vmin,vmax,smin,smax,p,indx_vel)
!
!v(s)=m*s+t
!
use prog_type
!
implicit none
!
! ... arguments
integer(i4b), intent(in) :: indx_vel
real(dp), intent(in) :: vmin, vmax, smin, smax, p
real(dp) :: vlaw
!
select case(indx_vel)
   case(0)
!linear velocity law
      vlaw=vmin+(vmax-vmin)/(smax-smin)*(p-smin)
   case(1)
!step function
      vlaw=vmax
      if(p.lt.(smin+smax)/2.d0) vlaw=vmin
   case default
      stop 'error in benchmark06_solution: im_vel not specified properly'
end select
!
end function vlaw
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
function slinelaw(sline_min, sline_max, smin, smax, p)
!
use prog_type
!
implicit none
!
! ... arguments
real(dp), intent(in) :: sline_min, sline_max, smin, smax, p
real(dp) :: slinelaw
!
!for constant source function
slinelaw = (sline_min+sline_max)/2.d0
!
!for linear law
slinelaw = sline_min + (sline_max-sline_min)/(smax-smin)*(p-smin)
!
!for quadratic law: sline=a*s^2+c
slinelaw = sline_min + (sline_max-sline_min)*(p/smax)**2
!
end function slinelaw
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
function scontlaw(scont_min, scont_max, smin, smax, p)
!
use prog_type
!
implicit none
!
! ... arguments
real(dp), intent(in) :: scont_min, scont_max, smin, smax, p
real(dp) :: scontlaw
!
!for constant source function
scontlaw = (scont_min+scont_max)/2.d0
!
!for linear law
scontlaw = scont_min + (scont_max-scont_min)/(smax-smin)*(p-smin)
!
!for quadratic law: scont=a*s^2+c
scontlaw = scont_min + (scont_max-scont_min)*(p/smax)**2
!
end function scontlaw
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
function opalbarlaw(opalbar_min, opalbar_max, smin, smax, p)
!
use prog_type
!
implicit none
!
! ... arguments
real(dp), intent(in) :: opalbar_min, opalbar_max, smin, smax, p
real(dp) :: opalbarlaw
!
!for constant values
opalbarlaw=(opalbar_min+opalbar_max)/2.d0
!
!for linear law
opalbarlaw = opalbar_min + (opalbar_max-opalbar_min)/(smax-smin)*(p-smin)
!
!for quadratic law: opalbar=a*s^2+c
opalbarlaw = opalbar_min + (opalbar_max-opalbar_min)*(p/smax)**2

end function opalbarlaw
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
function opaclaw(opac_min, opac_max, smin, smax, p)
!
use prog_type
!
implicit none
!
! ... arguments
real(dp), intent(in) :: opac_min, opac_max, smin, smax, p
real(dp) :: opaclaw
!
!for constant values
opaclaw=(opac_min+opac_max)/2.d0
!
!for linear law
opaclaw = opac_min + (opac_max-opac_min)/(smax-smin)*(p-smin)
!
!for quadratic law: opac=a*s^2+c
opaclaw = opac_min + (opac_max-opac_min)*(p/smax)**2

end function opaclaw
