!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine integration_test
!
!-------------------integrates weight functions of----------------------
!-----------------legendre and chebyshev integration--------------------
!----------to check if nodes/weights are calculated correctly-----------
!
   use prog_type
   use fund_const, only: pi
   use angles, only: nodes_mu, weight_mu, weight_omega, dim_omega, dim_mu
   use options, only: opt_angint_method, opt_sol2d
!
   implicit none
!
! ... local scalars
   integer(i4b) :: i, j, err
   real(dp) :: integral, integral_exact
!
! ... local arrays
   real(dp), dimension(:), allocatable :: yvalue
!
   write(*,*) '----------------------check integration nodes/weights--------------------------'
   write(*,*)
!
   if(opt_sol2d) then
      allocate(yvalue(dim_mu), stat=err)
      if(err.ne.0) stop 'allocation error in integration_test'
      yvalue=1.d0
      integral_exact=1.d0
      integral=sum(yvalue*weight_mu)
      write(*,fmt='(3(a20))') 'exact', 'integration', 'relative error'
      write(*,fmt='(3(es20.8))') integral_exact, integral, (integral_exact-integral)/integral_exact
      write(*,*)
   else
      allocate(yvalue(dim_omega), stat=err)
      if(err.ne.0) stop 'allocation error in integration_test'
      yvalue=1.d0
      integral_exact=1.d0
      integral=sum(yvalue*weight_omega)
      write(*,fmt='(3(a20))') 'exact', 'integration', 'relative error'
      write(*,fmt='(3(es20.8))') integral_exact, integral, (integral_exact-integral)/integral_exact
      write(*,*)
   endif
!
!
end subroutine integration_test
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine integration_test2
!
   use prog_type
   use mod_integ1d, only: precalc_weight_trapez, precalc_weight_simps, &
      precalc_weight_boole, precalc_weight_spline, precalc_oweight_trapez, &
      precalc_oweight_simps, precalc_oweight_boole
!
   implicit none
!
! ... local scalars
   integer(i4b), parameter :: nd=13
   integer(i4b) :: i
   real(dp) :: xmax, xmin, xmaxo, xmino
   real(dp) :: integ_theo, integ_trapez, integ_trapez2, integ_trapezo, &
      integ_simps, integ_simps2, integ_simpso, integ_simpso2, integ_spline
!
! ... local arrays
   real(dp), dimension(nd) :: x_val, x_valo, y_theo, y_theoo, w_trapez, w_trapezo, &
      w_simps, w_simps2, w_simpso, w_simpso2, w_spline
!
!----------------make x,y-grid on which integration is tested-----------
!
   xmax=3.d0
   xmin=0.d0
!
!for open boundaries
   xmaxo=2.9d0
   xmino=0.1d0
!
   do i=1, nd
      x_val(i) = xmin + (i-1)*(xmax-xmin)/(nd-1)
      x_valo(i) = xmino + (i-1)*(xmaxo-xmino)/(nd-1)
   enddo
!
!constant function
   y_theo=2.d0
   y_theoo=2.d0
   integ_theo=2.d0*(xmax - xmin)
!
!linear function
   y_theo=x_val
   y_theoo=x_valo
   integ_theo=0.5d0*(xmax*xmax - xmin*xmin)
!
!quadratic function
   y_theo=x_val**2
   y_theoo=x_valo**2
   integ_theo=(xmax**3 - xmin**3)/3.d0
!
!exponential function
   y_theo=exp(x_val)
   y_theoo=exp(x_valo)
   integ_theo=exp(xmax)-exp(xmin)
!
!sine function
!y_theo=sin(x_val)
!y_theoo=sin(x_valo)
!integ_theo=-cos(xmax)+cos(xmin)
!
!---------------------------calculate weights---------------------------
!
   call precalc_weight_trapez(x_val, nd, w_trapez)
   call precalc_weight_simps(x_val, nd, w_simps)
   call precalc_weight_boole(x_val, nd, w_simps2)
   call precalc_weight_spline(x_val, nd, w_spline, .false.)
!
!for open boundaries
   call precalc_oweight_trapez(x_valo, nd, xmin, xmax, w_trapezo)
   call precalc_oweight_simps(x_valo, nd, xmin, xmax, w_simpso)
   call precalc_oweight_boole(x_valo, nd, xmin, xmax, w_simpso2)
!
!--------------------calculate integral for test-function---------------
!
   integ_trapez=sum(y_theo*w_trapez)
   integ_simps=sum(y_theo*w_simps)
   integ_simps2=sum(y_theo*w_simps2)
   integ_spline=sum(y_theo*w_spline)
!
!for open boundaries
   integ_trapezo=sum(y_theo*w_trapezo)
   integ_simpso=sum(y_theo*w_simpso)
   integ_simpso2=sum(y_theo*w_simpso2)
!
!--------------------------output---------------------------------------
!
   write(*,*) '-------------several integrals----------------'
   write(*,'(a30, es20.8)') 'integ theo', integ_theo
   write(*,'(a30, 2es20.8)') 'integ trapez', integ_trapez, (integ_trapez-integ_theo)/integ_theo
   write(*,'(a30, 2es20.8)') 'integ simps', integ_simps, (integ_simps-integ_theo)/integ_theo
   write(*,'(a30, 2es20.8)') 'integ simps2', integ_simps2, (integ_simps2-integ_theo)/integ_theo
   write(*,'(a30, 2es20.8)') 'integ spline', integ_spline, (integ_spline-integ_theo)/integ_theo
   write(*,'(a30, 2es20.8)') 'integ trapez (open boundary)', integ_trapezo, (integ_trapezo-integ_theo)/integ_theo
   write(*,'(a30, 2es20.8)') 'integ simps (open boundary)', integ_simpso, (integ_simpso-integ_theo)/integ_theo
   write(*,'(a30, 2es20.8)') 'integ simps2 (open boundary)', integ_simpso2, (integ_simpso2-integ_theo)/integ_theo
   write(*,*)
!
!
end subroutine integration_test2
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine interpolation_test
!
!--------------interpolation of different test functions----------------
!------------to check which interpolation scheme works best-------------
!
   use prog_type
   use mod_directories, only: output_dir_test
   use hdf5
   use mod_interp1d, only: interpol_yp, interpol_typ_quad, interpol_typ_quad2, interpol_fyp_cube, &
      interpol_typ_cube, interpol_fyp_cube2, interpol_fyp_cube3, interpol_fyp_cube4
!
   implicit none
!
! ... local scalars
   integer(i4b) :: i
   integer(i4b), parameter :: n_data=4, n_out=100
   real(dp), parameter :: sim1=1.d0, si=2.d0
   real(dp) :: yim2, yim1, yi, yip1, sim2, sip1
   real(dp) :: dels
!
! ... local arrays
   real(dp), dimension(n_data) :: s_data, y_data1, y_data2, y_data3, y_data4, y_data5, &
      y_data6, y_data7, y_data8
   real(dp), dimension(n_out) :: s_interp, &
      y_interp11, y_interp12, y_interp13, y_interp14, y_interp15, y_interp16, y_interp17, y_interp18, &
      y_interp21, y_interp22, y_interp23, y_interp24, y_interp25, y_interp26, y_interp27, y_interp28, &
      y_interp31, y_interp32, y_interp33, y_interp34, y_interp35, y_interp36, y_interp37, y_interp38, &
      y_interp41, y_interp42, y_interp43, y_interp44, y_interp45, y_interp46, y_interp47, y_interp48, &
      y_interp51, y_interp52, y_interp53, y_interp54, y_interp55, y_interp56, y_interp57, y_interp58, &
      y_interp61, y_interp62, y_interp63, y_interp64, y_interp65, y_interp66, y_interp67, y_interp68, &
      y_interp71, y_interp72, y_interp73, y_interp74, y_interp75, y_interp76, y_interp77, y_interp78, &
      y_interp81, y_interp82, y_interp83, y_interp84, y_interp85, y_interp86, y_interp87, y_interp88
   real(dp), dimension(n_out) :: s_theo, y_theo1, y_theo2, y_theo3, y_theo4, y_theo5, &
      y_theo6, y_theo7, y_theo8
!
! ... local characters
   character(len=21) :: fname='test_interpolation.h5'
!
! ... local functions
!
! ... for output to hdf5
   integer(i4b) :: err
   integer(hid_t) :: file_id, dspace_id, dset_id, aspace_id, attr_id, group_id
   integer(hsize_t), dimension(1) :: dims_scalars = (/1/)
   integer(hsize_t), dimension(1) :: dims_ndata , dims_nout
!
!
!
   write(*,*) '--------------testing interpolation schemes for different functions------------'
   write(*,*)
!
   yim1=1.d0
   yi=7.d0
!
!prepare grid along direction and theoretical solution
   dels=si-sim1
   s_data=(/ sim1-dels, sim1, si, si+dels /)
   do i=1, n_out
      s_interp(i) = sim1 + (i-1)*(si-sim1)/float(n_out-1)
   enddo
   do i=1, n_out
      s_theo(i) = sim1-dels + (i-1)*(si-sim1+2.d0*dels)/float(n_out-1)
      call functions_test(s_theo(i), sim1, si, yim1, yi, y_theo1(i), y_theo2(i), &
         y_theo3(i), y_theo4(i), y_theo5(i), y_theo6(i), y_theo7(i), y_theo8(i))
   enddo
!
   s_interp=s_theo
!prepare data on the grid along direction
   do i=1, n_data
      call functions_test(s_data(i), sim1, si, yim1, yi, y_data1(i), y_data2(i), &
         y_data3(i), y_data4(i), y_data5(i), y_data6(i), y_data7(i), y_data8(i))
   enddo
!
!perform interpolations
   sim2=s_data(1)
   sip1=s_data(4)
!
   do i=1, n_out
!function 1
      yim2=y_data1(1)
      yim1=y_data1(2)
      yi=y_data1(3)
      yip1=y_data1(4)
      y_interp11(i)=interpol_yp(sim1, si, yim1, yi, s_interp(i))
      y_interp12(i)=interpol_typ_quad(yim1, yi, yip1, sim1, si, sip1, s_interp(i))
      y_interp13(i)=interpol_typ_quad2(yim1, yi, yip1, sim1, si, sip1, s_interp(i))
      y_interp14(i)=interpol_fyp_cube(yim2, yim1, yi, yip1, sim2, sim1, si, sip1, s_interp(i))
      y_interp15(i)=interpol_typ_cube(yim1, yi, yip1, sim1, si, sip1, s_interp(i))
      y_interp16(i)=interpol_fyp_cube2(yim2, yim1, yi, yip1, sim2, sim1, si, sip1, s_interp(i))
      y_interp17(i)=interpol_fyp_cube3(yim2, yim1, yi, yip1, sim2, sim1, si, sip1, s_interp(i), 3, 4)
      y_interp18(i)=interpol_fyp_cube4(yim2, yim1, yi, yip1, sim2, sim1, si, sip1, s_interp(i))
!function 2
      yim2=y_data2(1)
      yim1=y_data2(2)
      yi=y_data2(3)
      yip1=y_data2(4)
      y_interp21(i)=interpol_yp(sim1, si, yim1, yi, s_interp(i))
      y_interp22(i)=interpol_typ_quad(yim1, yi, yip1, sim1, si, sip1, s_interp(i))
      y_interp23(i)=interpol_typ_quad2(yim1, yi, yip1, sim1, si, sip1, s_interp(i))
      y_interp24(i)=interpol_fyp_cube(yim2, yim1, yi, yip1, sim2, sim1, si, sip1, s_interp(i))
      y_interp25(i)=interpol_typ_cube(yim1, yi, yip1, sim1, si, sip1, s_interp(i))
      y_interp26(i)=interpol_fyp_cube2(yim2, yim1, yi, yip1, sim2, sim1, si, sip1, s_interp(i))
      y_interp27(i)=interpol_fyp_cube3(yim2, yim1, yi, yip1, sim2, sim1, si, sip1, s_interp(i), 3, 4)
      y_interp28(i)=interpol_fyp_cube4(yim2, yim1, yi, yip1, sim2, sim1, si, sip1, s_interp(i))
!function 3
      yim2=y_data3(1)
      yim1=y_data3(2)
      yi=y_data3(3)
      yip1=y_data3(4)
      y_interp31(i)=interpol_yp(sim1, si, yim1, yi, s_interp(i))
      y_interp32(i)=interpol_typ_quad(yim1, yi, yip1, sim1, si, sip1, s_interp(i))
      y_interp33(i)=interpol_typ_quad2(yim1, yi, yip1, sim1, si, sip1, s_interp(i))
      y_interp34(i)=interpol_fyp_cube(yim2, yim1, yi, yip1, sim2, sim1, si, sip1, s_interp(i))
      y_interp35(i)=interpol_typ_cube(yim1, yi, yip1, sim1, si, sip1, s_interp(i))
      y_interp36(i)=interpol_fyp_cube2(yim2, yim1, yi, yip1, sim2, sim1, si, sip1, s_interp(i))
      y_interp37(i)=interpol_fyp_cube3(yim2, yim1, yi, yip1, sim2, sim1, si, sip1, s_interp(i), 3, 4)
      y_interp38(i)=interpol_fyp_cube4(yim2, yim1, yi, yip1, sim2, sim1, si, sip1, s_interp(i))
!function 4
      yim2=y_data4(1)
      yim1=y_data4(2)
      yi=y_data4(3)
      yip1=y_data4(4)
      y_interp41(i)=interpol_yp(sim1, si, yim1, yi, s_interp(i))
      y_interp42(i)=interpol_typ_quad(yim1, yi, yip1, sim1, si, sip1, s_interp(i))
      y_interp43(i)=interpol_typ_quad2(yim1, yi, yip1, sim1, si, sip1, s_interp(i))
      y_interp44(i)=interpol_fyp_cube(yim2, yim1, yi, yip1, sim2, sim1, si, sip1, s_interp(i))
      y_interp45(i)=interpol_typ_cube(yim1, yi, yip1, sim1, si, sip1, s_interp(i))
      y_interp46(i)=interpol_fyp_cube2(yim2, yim1, yi, yip1, sim2, sim1, si, sip1, s_interp(i))
      y_interp47(i)=interpol_fyp_cube3(yim2, yim1, yi, yip1, sim2, sim1, si, sip1, s_interp(i), 3, 4)
      y_interp48(i)=interpol_fyp_cube4(yim2, yim1, yi, yip1, sim2, sim1, si, sip1, s_interp(i))
!function 5
      yim2=y_data5(1)
      yim1=y_data5(2)
      yi=y_data5(3)
      yip1=y_data5(4)
      y_interp51(i)=interpol_yp(sim1, si, yim1, yi, s_interp(i))
      y_interp52(i)=interpol_typ_quad(yim1, yi, yip1, sim1, si, sip1, s_interp(i))
      y_interp53(i)=interpol_typ_quad2(yim1, yi, yip1, sim1, si, sip1, s_interp(i))
      y_interp54(i)=interpol_fyp_cube(yim2, yim1, yi, yip1, sim2, sim1, si, sip1, s_interp(i))
      y_interp55(i)=interpol_typ_cube(yim1, yi, yip1, sim1, si, sip1, s_interp(i))
      y_interp56(i)=interpol_fyp_cube2(yim2, yim1, yi, yip1, sim2, sim1, si, sip1, s_interp(i))
      y_interp57(i)=interpol_fyp_cube3(yim2, yim1, yi, yip1, sim2, sim1, si, sip1, s_interp(i), 3, 4)
      y_interp58(i)=interpol_fyp_cube4(yim2, yim1, yi, yip1, sim2, sim1, si, sip1, s_interp(i))
!function 6
      yim2=y_data6(1)
      yim1=y_data6(2)
      yi=y_data6(3)
      yip1=y_data6(4)
      y_interp61(i)=interpol_yp(sim1, si, yim1, yi, s_interp(i))
      y_interp62(i)=interpol_typ_quad(yim1, yi, yip1, sim1, si, sip1, s_interp(i))
      y_interp63(i)=interpol_typ_quad2(yim1, yi, yip1, sim1, si, sip1, s_interp(i))
      y_interp64(i)=interpol_fyp_cube(yim2, yim1, yi, yip1, sim2, sim1, si, sip1, s_interp(i))
      y_interp65(i)=interpol_typ_cube(yim1, yi, yip1, sim1, si, sip1, s_interp(i))
      y_interp66(i)=interpol_fyp_cube2(yim2, yim1, yi, yip1, sim2, sim1, si, sip1, s_interp(i))
      y_interp67(i)=interpol_fyp_cube3(yim2, yim1, yi, yip1, sim2, sim1, si, sip1, s_interp(i), 3, 4)
      y_interp68(i)=interpol_fyp_cube4(yim2, yim1, yi, yip1, sim2, sim1, si, sip1, s_interp(i))
!function 7
      yim2=y_data7(1)
      yim1=y_data7(2)
      yi=y_data7(3)
      yip1=y_data7(4)
      y_interp71(i)=interpol_yp(sim1, si, yim1, yi, s_interp(i))
      y_interp72(i)=interpol_typ_quad(yim1, yi, yip1, sim1, si, sip1, s_interp(i))
      y_interp73(i)=interpol_typ_quad2(yim1, yi, yip1, sim1, si, sip1, s_interp(i))
      y_interp74(i)=interpol_fyp_cube(yim2, yim1, yi, yip1, sim2, sim1, si, sip1, s_interp(i))
      y_interp75(i)=interpol_typ_cube(yim1, yi, yip1, sim1, si, sip1, s_interp(i))
      y_interp76(i)=interpol_fyp_cube2(yim2, yim1, yi, yip1, sim2, sim1, si, sip1, s_interp(i))
      y_interp77(i)=interpol_fyp_cube3(yim2, yim1, yi, yip1, sim2, sim1, si, sip1, s_interp(i), 3, 4)
      y_interp78(i)=interpol_fyp_cube4(yim2, yim1, yi, yip1, sim2, sim1, si, sip1, s_interp(i))
!function 8
      yim2=y_data8(1)
      yim1=y_data8(2)
      yi=y_data8(3)
      yip1=y_data8(4)
      y_interp81(i)=interpol_yp(sim1, si, yim1, yi, s_interp(i))
      y_interp82(i)=interpol_typ_quad(yim1, yi, yip1, sim1, si, sip1, s_interp(i))
      y_interp83(i)=interpol_typ_quad2(yim1, yi, yip1, sim1, si, sip1, s_interp(i))
      y_interp84(i)=interpol_fyp_cube(yim2, yim1, yi, yip1, sim2, sim1, si, sip1, s_interp(i))
      y_interp85(i)=interpol_typ_cube(yim1, yi, yip1, sim1, si, sip1, s_interp(i))
      y_interp86(i)=interpol_fyp_cube2(yim2, yim1, yi, yip1, sim2, sim1, si, sip1, s_interp(i))
      y_interp87(i)=interpol_fyp_cube3(yim2, yim1, yi, yip1, sim2, sim1, si, sip1, s_interp(i), 3, 4)
      y_interp88(i)=interpol_fyp_cube4(yim2, yim1, yi, yip1, sim2, sim1, si, sip1, s_interp(i))
   enddo
!
!linear interolation for all functions: interpol_yp
!quadratic bezier spline: interpol_typ_quad
!monotonic quadratic bezier spline: interpol_typ_quad2
!cubic catmull rom spline: interpol_fyp_cube
!cubic spline using only three points: interpol_typ_cube
!cubic spline using harmonic mean for derivatives: interpol_fyp_cube2 (monotonic from Ibgui 2013)
!cubic spline using harmonic mean for derivatives (if appropriate, monotonic spline from Steffen 1990):
!   interpol_fyp_cube3
!cubic spline using weighted mean for derivatives: interpol_fyp_cube4
!
   write(*,*) 'output to file ', output_dir_test//'/'//fname
   write(*,*)
!
!
   dims_ndata = (/ n_data /)
   dims_nout = (/ n_out /)
!
   call h5open_f (err)
   call h5fcreate_f(output_dir_test//'/'//fname, h5f_acc_trunc_f, file_id, err)
!
!----------------------------dimensions---------------------------------
!
   call h5gcreate_f(file_id, 'dimensions', group_id, err)
   call h5screate_simple_f(1, dims_scalars, aspace_id, err)
   call h5acreate_f(group_id, 'n_data', h5t_native_integer, aspace_id, &
      attr_id, err)
   call h5awrite_f(attr_id, h5t_native_integer, n_data, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5acreate_f(group_id, 'n_grid', h5t_native_integer, aspace_id, &
      attr_id, err)
   call h5awrite_f(attr_id, h5t_native_integer, n_out, dims_scalars, err)
   call h5aclose_f(attr_id, err)
   call h5sclose_f(aspace_id, err)
   call h5gclose_f(group_id, err)
!
!----------------------------coordinates--------------------------------
!
   call h5gcreate_f(file_id, 'coordinates', group_id, err)
   call h5screate_simple_f(1, dims_ndata, dspace_id, err)
   call h5dcreate_f(group_id, 'grid_data', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, s_data, dims_ndata, err)
   call h5dclose_f(dset_id, err)
   call h5sclose_f(dspace_id, err)
   call h5screate_simple_f(1, dims_nout, dspace_id, err)
   call h5dcreate_f(group_id, 'grid_interp', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, s_interp, dims_nout, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'grid_theo', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, s_theo, dims_nout, err)
   call h5dclose_f(dset_id, err)
   call h5sclose_f(dspace_id, err)
   call h5gclose_f(group_id, err)
!
!----------------------------solution-----------------------------------
!
   call h5gcreate_f(file_id, 'function1', group_id, err)
   call h5screate_simple_f(1, dims_nout, dspace_id, err)
   call h5dcreate_f(group_id, 'y_interp11', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, y_interp11, dims_nout, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'y_interp12', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, y_interp12, dims_nout, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'y_interp13', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, y_interp13, dims_nout, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'y_interp14', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, y_interp14, dims_nout, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'y_interp15', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, y_interp15, dims_nout, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'y_interp16', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, y_interp16, dims_nout, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'y_interp17', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, y_interp17, dims_nout, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'y_interp18', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, y_interp18, dims_nout, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'y_theo1', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, y_theo1, dims_nout, err)
   call h5dclose_f(dset_id, err)
   call h5sclose_f(dspace_id, err)
   call h5screate_simple_f(1, dims_ndata, dspace_id, err)
   call h5dcreate_f(group_id, 'y_data1', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, y_data1, dims_ndata, err)
   call h5dclose_f(dset_id, err)
   call h5sclose_f(dspace_id, err)
   call h5gclose_f(group_id, err)
!
   call h5gcreate_f(file_id, 'function2', group_id, err)
   call h5screate_simple_f(1, dims_nout, dspace_id, err)
   call h5dcreate_f(group_id, 'y_interp21', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, y_interp21, dims_nout, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'y_interp22', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, y_interp22, dims_nout, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'y_interp23', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, y_interp23, dims_nout, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'y_interp24', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, y_interp24, dims_nout, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'y_interp25', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, y_interp25, dims_nout, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'y_interp26', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, y_interp26, dims_nout, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'y_interp27', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, y_interp27, dims_nout, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'y_interp28', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, y_interp28, dims_nout, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'y_theo2', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, y_theo2, dims_nout, err)
   call h5dclose_f(dset_id, err)
   call h5sclose_f(dspace_id, err)
   call h5screate_simple_f(1, dims_ndata, dspace_id, err)
   call h5dcreate_f(group_id, 'y_data2', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, y_data2, dims_ndata, err)
   call h5dclose_f(dset_id, err)
   call h5sclose_f(dspace_id, err)
   call h5gclose_f(group_id, err)
!
   call h5gcreate_f(file_id, 'function3', group_id, err)
   call h5screate_simple_f(1, dims_nout, dspace_id, err)
   call h5dcreate_f(group_id, 'y_interp31', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, y_interp31, dims_nout, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'y_interp32', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, y_interp32, dims_nout, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'y_interp33', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, y_interp33, dims_nout, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'y_interp34', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, y_interp34, dims_nout, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'y_interp35', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, y_interp35, dims_nout, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'y_interp36', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, y_interp36, dims_nout, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'y_interp37', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, y_interp37, dims_nout, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'y_interp38', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, y_interp38, dims_nout, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'y_theo3', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, y_theo3, dims_nout, err)
   call h5dclose_f(dset_id, err)
   call h5sclose_f(dspace_id, err)
   call h5screate_simple_f(1, dims_ndata, dspace_id, err)
   call h5dcreate_f(group_id, 'y_data3', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, y_data3, dims_ndata, err)
   call h5dclose_f(dset_id, err)
   call h5sclose_f(dspace_id, err)
   call h5gclose_f(group_id, err)
!
   call h5gcreate_f(file_id, 'function4', group_id, err)
   call h5screate_simple_f(1, dims_nout, dspace_id, err)
   call h5dcreate_f(group_id, 'y_interp41', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, y_interp41, dims_nout, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'y_interp42', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, y_interp42, dims_nout, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'y_interp43', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, y_interp43, dims_nout, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'y_interp44', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, y_interp44, dims_nout, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'y_interp45', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, y_interp45, dims_nout, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'y_interp46', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, y_interp46, dims_nout, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'y_interp47', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, y_interp47, dims_nout, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'y_interp48', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, y_interp48, dims_nout, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'y_theo4', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, y_theo4, dims_nout, err)
   call h5dclose_f(dset_id, err)
   call h5sclose_f(dspace_id, err)
   call h5screate_simple_f(1, dims_ndata, dspace_id, err)
   call h5dcreate_f(group_id, 'y_data4', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, y_data4, dims_ndata, err)
   call h5dclose_f(dset_id, err)
   call h5sclose_f(dspace_id, err)
   call h5gclose_f(group_id, err)
!
   call h5gcreate_f(file_id, 'function5', group_id, err)
   call h5screate_simple_f(1, dims_nout, dspace_id, err)
   call h5dcreate_f(group_id, 'y_interp51', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, y_interp51, dims_nout, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'y_interp52', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, y_interp52, dims_nout, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'y_interp53', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, y_interp53, dims_nout, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'y_interp54', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, y_interp54, dims_nout, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'y_interp55', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, y_interp55, dims_nout, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'y_interp56', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, y_interp56, dims_nout, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'y_interp57', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, y_interp57, dims_nout, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'y_interp58', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, y_interp58, dims_nout, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'y_theo5', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, y_theo5, dims_nout, err)
   call h5dclose_f(dset_id, err)
   call h5sclose_f(dspace_id, err)
   call h5screate_simple_f(1, dims_ndata, dspace_id, err)
   call h5dcreate_f(group_id, 'y_data5', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, y_data5, dims_ndata, err)
   call h5dclose_f(dset_id, err)
   call h5sclose_f(dspace_id, err)
   call h5gclose_f(group_id, err)
!
   call h5gcreate_f(file_id, 'function6', group_id, err)
   call h5screate_simple_f(1, dims_nout, dspace_id, err)
   call h5dcreate_f(group_id, 'y_interp61', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, y_interp61, dims_nout, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'y_interp62', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, y_interp62, dims_nout, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'y_interp63', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, y_interp63, dims_nout, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'y_interp64', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, y_interp64, dims_nout, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'y_interp65', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, y_interp65, dims_nout, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'y_interp66', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, y_interp66, dims_nout, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'y_interp67', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, y_interp67, dims_nout, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'y_interp68', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, y_interp68, dims_nout, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'y_theo6', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, y_theo6, dims_nout, err)
   call h5dclose_f(dset_id, err)
   call h5sclose_f(dspace_id, err)
   call h5screate_simple_f(1, dims_ndata, dspace_id, err)
   call h5dcreate_f(group_id, 'y_data6', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, y_data6, dims_ndata, err)
   call h5dclose_f(dset_id, err)
   call h5sclose_f(dspace_id, err)
   call h5gclose_f(group_id, err)
!
   call h5gcreate_f(file_id, 'function7', group_id, err)
   call h5screate_simple_f(1, dims_nout, dspace_id, err)
   call h5dcreate_f(group_id, 'y_interp71', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, y_interp71, dims_nout, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'y_interp72', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, y_interp72, dims_nout, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'y_interp73', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, y_interp73, dims_nout, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'y_interp74', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, y_interp74, dims_nout, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'y_interp75', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, y_interp75, dims_nout, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'y_interp76', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, y_interp76, dims_nout, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'y_interp77', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, y_interp77, dims_nout, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'y_interp78', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, y_interp78, dims_nout, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'y_theo7', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, y_theo7, dims_nout, err)
   call h5dclose_f(dset_id, err)
   call h5sclose_f(dspace_id, err)
   call h5screate_simple_f(1, dims_ndata, dspace_id, err)
   call h5dcreate_f(group_id, 'y_data7', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, y_data7, dims_ndata, err)
   call h5dclose_f(dset_id, err)
   call h5sclose_f(dspace_id, err)
   call h5gclose_f(group_id, err)
!
   call h5gcreate_f(file_id, 'function8', group_id, err)
   call h5screate_simple_f(1, dims_nout, dspace_id, err)
   call h5dcreate_f(group_id, 'y_interp81', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, y_interp81, dims_nout, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'y_interp82', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, y_interp82, dims_nout, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'y_interp83', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, y_interp83, dims_nout, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'y_interp84', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, y_interp84, dims_nout, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'y_interp85', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, y_interp85, dims_nout, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'y_interp86', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, y_interp86, dims_nout, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'y_interp87', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, y_interp87, dims_nout, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'y_interp88', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, y_interp88, dims_nout, err)
   call h5dclose_f(dset_id, err)
   call h5dcreate_f(group_id, 'y_theo8', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, y_theo8, dims_nout, err)
   call h5dclose_f(dset_id, err)
   call h5sclose_f(dspace_id, err)
   call h5screate_simple_f(1, dims_ndata, dspace_id, err)
   call h5dcreate_f(group_id, 'y_data8', h5t_native_real, dspace_id, &
      dset_id, err)
   call h5dwrite_f(dset_id, h5t_native_double, y_data8, dims_ndata, err)
   call h5dclose_f(dset_id, err)
   call h5sclose_f(dspace_id, err)
   call h5gclose_f(group_id, err)
!
   call h5fclose_f(file_id, err)
   call h5close_f(err)
!
!
!
end subroutine interpolation_test
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine interpolation_test2
!
   use prog_type
   use mod_interp2d, only: interpol2d_9p_quad, interpol2d_9p_bez, coeff2d_9p_bez2
!
   implicit none
!
! ... local scalars
   integer(i4b) :: i,j
   integer(i4b) :: iim2, iim1, ii, jjm2, jjm1, jj
   integer(i4b), parameter :: nxp=3, nyp=3, nx=10, ny=10
   real(dp) :: xmin, xmax, ymin, ymax, rad, xdum, ydum
   real(dp) :: c01, c02, c03, c04, c05, c06, c07, c08, c09
!
! ... local arrays
   real(dp), dimension(nxp) :: xp
   real(dp), dimension(nyp) :: yp
   real(dp), dimension(nxp,nyp) :: fp2d
   real(dp), dimension(nx) :: x
   real(dp), dimension(ny) :: y
   real(dp), dimension(nx,ny) :: f2d_method1, f2d_method2
!
! ... local functions
!
!-----------------------------create grids and data---------------------
   xmin=1.d0
   xmax=2.d0
   ymin=1.d0
   ymax=2.d0
!
   do i=1, nxp
      xp(i) = xmin + (i-1)*(xmax-xmin)/(nxp-1)
   enddo
   do i=1, nyp
      yp(i) = ymin + (i-1)*(ymax-ymin)/(nyp-1)
   enddo
   do i=1, nxp
      do j=1, nyp
         rad=sqrt(xp(i)**2+yp(j)**2)
         if(rad.ge.2.2d0) then
            fp2d(i,j)=1.d0
         else
            fp2d(i,j)=0.d0
         endif
!     fp2d(i,j)=1.d0/rad
      enddo
   enddo
!
   open(1, file='TRASH/out_debug.dat', form='formatted')
   read(1,'(17es40.20)') fp2d(1,1), fp2d(2,1), fp2d(3,1), &
      fp2d(1,2), fp2d(2,2), fp2d(3,2), &
      fp2d(1,3), fp2d(2,3), fp2d(3,3), &
      xp(1), xp(2), xp(3), yp(1), yp(2), yp(3), xdum, ydum
   close(1)
!
   x=xdum
   y=ydum
!
!do i=1, nx
!   x(i) = xp(2) + (i-1)*(xp(3)-xp(2))/(nx-1)
!enddo
!do i=1, ny
!   y(i) = yp(2) + (i-1)*(yp(3)-yp(2))/(ny-1)
!enddo
!
!-------------------perform interpolations------------------------------
!
   iim2=1
   iim1=2
   ii=3
   jjm2=1
   jjm1=2
   jj=3
!
   do i=1, nx
      do j=1, ny

         call coeff2d_9p_bez2(fp2d(iim2,jjm2), fp2d(iim1,jjm2), fp2d(ii,jjm2), &
            fp2d(iim2,jjm1), fp2d(iim1,jjm1), fp2d(ii,jjm1), &
            fp2d(iim2,jj),   fp2d(iim1,jj),   fp2d(ii,jj), &
            xp(iim2), xp(iim1), xp(ii), yp(jjm2), yp(jjm1), yp(jj), &
            x(i), y(j), c01, c02, c03, c04, c05, c06, c07, c08, c09)
         f2d_method2(i,j) = c01*fp2d(iim2,jjm2) + c02*fp2d(iim1,jjm2) + c03*fp2d(ii,jjm2) + &
            c04*fp2d(iim2,jjm1) + c05*fp2d(iim1,jjm1) + c06*fp2d(ii,jjm1) + &
            c07*fp2d(iim2,jj)   + c08*fp2d(iim1,jj)   + c09*fp2d(ii,jj)

         f2d_method1(i,j) = interpol2d_9p_bez(fp2d(iim2,jjm2), fp2d(iim1,jjm2), fp2d(ii,jjm2), &
            fp2d(iim2,jjm1), fp2d(iim1,jjm1), fp2d(ii,jjm1), &
            fp2d(iim2,jj),   fp2d(iim1,jj),   fp2d(ii,jj), &
            xp(iim2), xp(iim1), xp(ii), yp(jjm2), yp(jjm1), yp(jj), &
            x(i), y(j))

         write(*,*) f2d_method1(i,j), f2d_method2(i,j)
         write(*,*) c01, c02, c03, c04, c05, c06, c07, c08, c09
         write(*,*)
      enddo
   enddo
   stop




end subroutine interpolation_test2
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine interpolation_test2d
!
!----------------2d interpolation tests on upwind points----------------
!------------------for a specified direction theta, phi-----------------
!
   use prog_type
   use fund_const, only: pi
   use dime3d, only: x, y, z, ndxmax, ndymax, ndzmax, imask3d
   use mod_directories, only: output_dir_test
   use mod_interp2d, only: interpol2d_4p_lin, interpol2d_9p_quad, interpol2d_9p_bez
   use mod_interp3d, only: interpol3d_8p_lin
!
   implicit none
!
! ... local scalars
   integer(i4b) :: i, j, k
   integer(i4b) :: startx, endx, alpha, starty, endy, beta, startz, endz, gamma
   integer(i4b) :: iim2, iim1, ii, iip1, jjm2, jjm1, jj, jjp1, kkm2, kkm1, kk, kkp1
   real(dp), parameter :: theta=54.7d0*pi/180.d0, phi=45.d0*pi/180.d0
   real(dp) :: n_x, n_y, n_z, rad
   real(dp) :: dels_xyu, dels_xzu, dels_yzu, dels_u, dels_r
   real(dp) :: dels_xyd, dels_xzd, dels_yzd, dels_d
   real(dp) :: x_u, y_u, z_u, f_u
   real(dp) :: x_d, y_d, z_d, f_d
   real(dp) :: f_test, mueff
!
! ... local arrays
   real(dp), dimension(:,:,:), allocatable :: f3d
!
! ... local characters
   character(len=25) :: fname1='test_interpolation2da.dat'
   character(len=25) :: fname2='test_interpolation2db.dat'
!
! ... local functions
   real(dp) :: dist_sphere
!
! ... for output to hdf5
!
!
   write(*,*) '------------------testing 2d interpolation scheme------------------------------'
   write(*,*)
!y
!given direction
   n_x=sin(theta)*cos(phi)
   n_y=sin(theta)*sin(phi)
   n_z=cos(theta)
   write(*,*) 'n_x', n_x
   write(*,*) 'n_y', n_y
   write(*,*) 'n_z', n_z
!
!---------------calculate function that shall be interpolated-----------
!
   allocate(f3d(ndxmax,ndymax,ndzmax))
   f3d=0.d0
   do i=1, ndxmax
      do j=1, ndymax
         do k=1, ndzmax
            rad=sqrt(x(i)**2+y(j)**2+z(k)**2)
            select case(imask3d(i,j,k))
             case(1)
               f3d(i,j,k)=exp(1.d0-rad)!1.d0/rad^3
             case(2)
               f3d(i,j,k)=exp(1.d0-rad)!1.d0/rad^3
             case(3)
               f3d(i,j,k)=1.d0
             case(4)
               f3d(i,j,k)=1.d0
             case default
            end select
         enddo
      enddo
   enddo
!
   write(*,*) 'output to file ', output_dir_test//'/'//fname1
   write(*,*) 'output to file ', output_dir_test//'/'//fname2
   write(*,*)
!
   open(1, file=output_dir_test//'/'//fname1, form='formatted')
   do i=1, ndxmax
      do j=1, ndymax
         do k=1, ndzmax
            rad=sqrt(x(i)**2+y(j)**2+z(k)**2)
            write(1,'(2es20.8)') rad, f3d(i,j,k)
         enddo
      enddo
   enddo
   close(1)

!
!------------------actual interpolation on upwind points----------------
!
   if(n_x.gt.0.d0) then
      startx = 3
      endx = ndxmax-2
      alpha=  1
   elseif(n_x.lt.0.d0) then
      startx = ndxmax-2
      endx = 3
      alpha=-1
   else
      stop 'error in interpolation_test2d: n_x = 0 not allowed'
   endif
!
   if(n_y.gt.0.d0) then
      starty = 3
      endy = ndymax-2
      beta =  1
   elseif(n_y.lt.0.d0) then
      starty = ndymax-2
      endy = 3
      beta =-1
   else
      stop 'error in interpolation_test2d: n_y = 0 not allowed'
   endif
!
   if(n_z.gt.0.d0) then
      startz = 3
      endz = ndzmax-2
      gamma=  1
   elseif(n_z.lt.0.d0) then
      startz = ndzmax-2
      endz = 3
      gamma=-1
   else
      stop 'error in interpolation_test2d: n_z = 0 not allowed'
   endif
!
!-----------------------------------------------------------------------
!
   open(1, file=output_dir_test//'/'//fname2, form='formatted')
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
!calculate distance of upwind point to (previous) xy-plane, xz-plane and yz-plane
               dels_xyu=(z(k)-z(k-gamma))/n_z
               dels_xzu=(y(j)-y(j-beta))/n_y
               dels_yzu=(x(i)-x(i-alpha))/n_x
               dels_u=min(dels_xyu,dels_xzu,dels_yzu)
!
!calculate distance of downwind point to (next) xy-plane, xz-plane and yz-plane
               dels_xyd=(z(k+gamma)-z(k))/n_z
               dels_xzd=(y(j+beta)-y(j))/n_y
               dels_yzd=(x(i+alpha)-x(i))/n_x
               dels_d=min(dels_xyd,dels_xzd,dels_yzd)
!
!----------------------------upwind point-------------------------------
!
               if(dels_xyu.eq.dels_u) then
!intersection with x-y plane on level k-gamma
                  x_u = x(i) - dels_u*n_x
                  y_u = y(j) - dels_u*n_y
                  z_u = z(k-gamma)
!
                  iim2=i-2*alpha
                  iim1=i-alpha
                  jjm2=j-2*beta
                  jjm1=j-beta
                  kkm1=k-gamma
!
!                  f_u = interpol2d_4p_lin(f3d(iim1,jjm1,kkm1), f3d(i,jjm1,kkm1), &
!                                          f3d(iim1,j,   kkm1), f3d(i,j,   kkm1), &
!                                          x(iim1), x(i), y(jjm1), y(j), x_u, y_u)
                  f_u = interpol2d_9p_quad(f3d(iim2,jjm2,kkm1), f3d(iim1,jjm2,kkm1), f3d(i,jjm2,kkm1), &
                     f3d(iim2,jjm1,kkm1), f3d(iim1,jjm1,kkm1), f3d(i,jjm1,kkm1), &
                     f3d(iim2,j,   kkm1), f3d(iim1,j,   kkm1), f3d(i,j,   kkm1), &
                     x(iim2), x(iim1), x(i), y(jjm2), y(jjm1), y(j), x_u, y_u)
                  rad=sqrt(x_u**2+y_u**2+z_u**2)
                  write(1,'(2e20.8,2a10,i10)') rad, f_u, 'u', 'xy', imask3d(i,j,k)
!
               elseif(dels_xzu.eq.dels_u) then
!intersection with x-z plane on level j-beta
                  x_u = x(i) - dels_u*n_x
                  y_u = y(j-beta)
                  z_u = z(k) - dels_u*n_z
!
                  iim2=i-2*alpha
                  iim1=i-alpha
                  jjm1=j-beta
                  kkm2=k-2*gamma
                  kkm1=k-gamma
!
!                  f_u = interpol2d_4p_lin(f3d(iim1,jjm1,kkm1), f3d(i,jjm1,kkm1), &
!                                          f3d(iim1,jjm1,k),    f3d(i,jjm1,k), &
!                                          x(iim1), x(i), z(kkm1), z(k), x_u, z_u)
                  f_u = interpol2d_9p_quad(f3d(iim2,jjm1,kkm2), f3d(iim1,jjm1,kkm2), f3d(i,jjm1,kkm2), &
                     f3d(iim2,jjm1,kkm1), f3d(iim1,jjm1,kkm1), f3d(i,jjm1,kkm1), &
                     f3d(iim2,jjm1,k),    f3d(iim1,jjm1,k),    f3d(i,jjm1,k), &
                     x(iim2), x(iim1), x(i), z(kkm2), z(kkm1), z(k), x_u, z_u)
                  rad=sqrt(x_u**2+y_u**2+z_u**2)
                  write(1,'(2e20.8,2a10,i10)') rad, f_u, 'u', 'xz', imask3d(i,j,k)
!
               elseif(dels_yzu.eq.dels_u) then
!intersection with y-z plane on level i-alpha
                  x_u = x(i-alpha)
                  y_u = y(j) - dels_u*n_y
                  z_u = z(k) - dels_u*n_z
!
                  iim1=i-alpha
                  jjm2=j-2*beta
                  jjm1=j-beta
                  kkm2=k-2*gamma
                  kkm1=k-gamma
!
!                  f_u = interpol2d_4p_lin(f3d(iim1,jjm1,kkm1), f3d(iim1,j,kkm1), &
!                                          f3d(iim1,jjm1,k),    f3d(iim1,j,k), &
!                                          y(jjm1), y(j), z(kkm1), z(k), y_u, z_u)
                  f_u = interpol2d_9p_quad(f3d(iim1,jjm2,kkm2), f3d(iim1,jjm1,kkm2), f3d(iim1,j,kkm2), &
                     f3d(iim1,jjm2,kkm1), f3d(iim1,jjm1,kkm1), f3d(iim1,j,kkm1), &
                     f3d(iim1,jjm2,k),    f3d(iim1,jjm1,k),    f3d(iim1,j,k), &
                     y(jjm2), y(jjm1), y(j), z(kkm2), z(kkm1), z(k), y_u, z_u)
                  rad=sqrt(x_u**2+y_u**2+z_u**2)
                  write(1,'(2e20.8,2a10,i10)') rad, f_u, 'u', 'yz', imask3d(i,j,k)
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
                  x_d = x(i) + dels_d*n_x
                  y_d = y(j) + dels_d*n_y
                  z_d = z(k+gamma)
!
                  iim1=i-alpha
                  iip1=i+alpha
                  jjm1=j-beta
                  jjp1=j+beta
                  kkp1=k+gamma
!
!                  f_d = interpol2d_4p_lin(f3d(i,j,  kkp1), f3d(iip1,j,  kkp1), &
!                                          f3d(i,jjp1,kkp1), f3d(iip1,jjp1,kkp1), &
!                                          x(i), x(iip1), y(j), y(jjp1), x_d, y_d)
                  f_d = interpol2d_9p_quad(f3d(iim1,jjm1,kkp1), f3d(i,jjm1,kkp1), f3d(iip1,jjm1,kkp1), &
                     f3d(iim1,j,   kkp1), f3d(i,j,   kkp1), f3d(iip1,j,   kkp1), &
                     f3d(iim1,jjp1,kkp1), f3d(i,jjp1,kkp1), f3d(iip1,jjp1,kkp1), &
                     x(iim1), x(i), x(iip1), y(jjm1), y(j), y(jjp1), x_d, y_d)
                  rad=sqrt(x_d**2+y_d**2+z_d**2)
                  write(1,'(2e20.8,2a10,i10)') rad, f_d, 'd', 'xy', imask3d(i,j,k)
!
               elseif(dels_xzd.eq.dels_d) then
!intersection with x-z plane on level j+beta
                  x_d = x(i) + dels_d*n_x
                  y_d = y(j+beta)
                  z_d = z(k) + dels_d*n_z
!
                  iim1=i-alpha
                  iip1=i+alpha
                  jjp1=j+beta
                  kkm1=k-gamma
                  kkp1=k+gamma
!
!                  f_d = interpol2d_4p_lin(f3d(i,jjp1,k),    f3d(iip1,jjp1,k), &
!                                          f3d(i,jjp1,kkp1), f3d(iip1,jjp1,kkp1), &
!                                          x(i), x(iip1), z(k), z(kkp1), x_d, z_d)
                  f_d = interpol2d_9p_quad(f3d(iim1,jjp1,kkm1), f3d(i,jjp1,kkm1), f3d(iip1,jjp1,kkm1), &
                     f3d(iim1,jjp1,k),    f3d(i,jjp1,k),    f3d(iip1,jjp1,k), &
                     f3d(iim1,jjp1,kkp1), f3d(i,jjp1,kkp1), f3d(iip1,jjp1,kkp1), &
                     x(iim1), x(i), x(iip1), z(kkm1), z(k), z(kkp1), x_d, z_d)
                  rad=sqrt(x_d**2+y_d**2+z_d**2)
                  write(1,'(2e20.8,2a10,i10)') rad, f_d, 'd', 'xz', imask3d(i,j,k)
!
               elseif(dels_yzd.eq.dels_d) then
!intersection with y-z plane on level i+alpha
                  x_d = x(i+alpha)
                  y_d = y(j) + dels_d*n_y
                  z_d = z(k) + dels_d*n_z
!
                  iip1=i+alpha
                  jjm1=j-beta
                  jjp1=j+beta
                  kkm1=k-gamma
                  kkp1=k+gamma
!
!                  f_d = interpol2d_4p_lin(f3d(iip1,j,k),    f3d(iip1,jjp1,k), &
!                                          f3d(iip1,j,kkp1), f3d(iip1,jjp1,kkp1), &
!                                          y(j), y(jjp1), z(k), z(kkp1), y_d, z_d)
                  f_d = interpol2d_9p_quad(f3d(iip1,jjm1,kkm1), f3d(iip1,j,kkm1), f3d(iip1,jjp1,kkm1), &
                     f3d(iip1,jjm1,k),    f3d(iip1,j,k),    f3d(iip1,jjp1,k), &
                     f3d(iip1,jjm1,kkp1), f3d(iip1,j,kkp1), f3d(iip1,jjp1,kkp1), &
                     y(jjm1), y(j), y(jjp1), z(kkm1), z(k), z(kkp1), y_d, z_d)
                  rad=sqrt(x_d**2+y_d**2+z_d**2)
                  write(1,'(2e20.8,2a10,i10)') rad, f_d, 'd', 'yz', imask3d(i,j,k)
!
               else
                  write(*,'(4es20.8,l4)') dels_d, dels_xzd, dels_xyd, dels_yzd, dels_d.eq.dels_xzd
                  stop 'error in fsc_cont3d: invalid dels_d'
               endif
!
!***********special treatment: point is in vicinity of boundary*********
!
             case(2)
!
!calculate distance of upwind point to (previous) xy-plane, xz-plane and yz-plane
               dels_xyu=(z(k)-z(k-gamma))/n_z
               dels_xzu=(y(j)-y(j-beta))/n_y
               dels_yzu=(x(i)-x(i-alpha))/n_x
!
!calculate distance of downwind point to (next) xy-plane, xz-plane and yz-plane
               dels_xyd=(z(k+gamma)-z(k))/n_z
               dels_xzd=(y(j+beta)-y(j))/n_y
               dels_yzd=(x(i+alpha)-x(i))/n_x
!
!calculate distance to the boundary (here: sphere)
               dels_r=dist_sphere(n_x,n_y,n_z,x(i),y(j),z(k))
!
               dels_u=min(dels_xyu,dels_xzu,dels_yzu,dels_r)
               dels_d=min(dels_xyd,dels_xzd,dels_yzd)
!
!----------------------------upwind point-------------------------------
!
               if(dels_r.eq.dels_u) then
                  x_u = x(i) - dels_u*n_x
                  y_u = y(j) - dels_u*n_y
                  z_u = z(k) - dels_u*n_z
!
                  iim1=i-alpha
                  jjm1=j-beta
                  kkm1=k-gamma
!
                  f_u = interpol3d_8p_lin(f3d(iim1,jjm1,kkm1), f3d(i,jjm1,kkm1), f3d(iim1,j,kkm1), f3d(i,j,kkm1), &
                     f3d(iim1,jjm1,k), f3d(i,jjm1,k), f3d(iim1,j,k), f3d(i,j,k), &
                     x(iim1), x(i), y(jjm1), y(j), z(kkm1), z(k), x_u, y_u, z_u)
!
                  rad=sqrt(x_u**2+y_u**2+z_u**2)
                  write(1,'(2e20.8,2a10,i10)') rad, f_u, 'u', 's', imask3d(i,j,k)
!
!                  if(rad.ne.1.d0) then
!                     write(*,*) 'error for dels_r: radius not correct'
!                     write(*,*) rad
!                     write(*,*) dels_r
!                     write(*,*) x_u, y_u, z_u
!                     stop
!                  endif

               elseif(dels_xyu.eq.dels_u) then
!intersection with x-y plane on level k-gamma
                  x_u = x(i) - dels_u*n_x
                  y_u = y(j) - dels_u*n_y
                  z_u = z(k-gamma)
!
                  iim2=i-2*alpha
                  iim1=i-alpha
                  jjm2=j-2*beta
                  jjm1=j-beta
                  kkm1=k-gamma
!
!                  f_u = interpol2d_4p_lin(f3d(iim1,jjm1,kkm1), f3d(i,jjm1,kkm1), &
!                                          f3d(iim1,j,   kkm1), f3d(i,j,   kkm1), &
!                                          x(iim1), x(i), y(jjm1), y(j), x_u, y_u)
                  f_u = interpol2d_9p_quad(f3d(iim2,jjm2,kkm1), f3d(iim1,jjm2,kkm1), f3d(i,jjm2,kkm1), &
                     f3d(iim2,jjm1,kkm1), f3d(iim1,jjm1,kkm1), f3d(i,jjm1,kkm1), &
                     f3d(iim2,j,   kkm1), f3d(iim1,j,   kkm1), f3d(i,j,   kkm1), &
                     x(iim2), x(iim1), x(i), y(jjm2), y(jjm1), y(j), x_u, y_u)
!
                  rad=sqrt(x_u**2+y_u**2+z_u**2)
                  write(1,'(2e20.8,2a10,i10)') rad, f_u, 'u', 'xy', imask3d(i,j,k)
!
               elseif(dels_xzu.eq.dels_u) then
!intersection with x-z plane on level j-beta
                  x_u = x(i) - dels_u*n_x
                  y_u = y(j-beta)
                  z_u = z(k) - dels_u*n_z
!
                  iim2=i-2*alpha
                  iim1=i-alpha
                  jjm1=j-beta
                  kkm2=k-2*gamma
                  kkm1=k-gamma
!
!                  f_u = interpol2d_4p_lin(f3d(iim1,jjm1,kkm1), f3d(i,jjm1,kkm1), &
!                                          f3d(iim1,jjm1,k),    f3d(i,jjm1,k), &
!                                          x(iim1), x(i), z(kkm1), z(k), x_u, z_u)
                  f_u = interpol2d_9p_quad(f3d(iim2,jjm1,kkm2), f3d(iim1,jjm1,kkm2), f3d(i,jjm1,kkm2), &
                     f3d(iim2,jjm1,kkm1), f3d(iim1,jjm1,kkm1), f3d(i,jjm1,kkm1), &
                     f3d(iim2,jjm1,k),    f3d(iim1,jjm1,k),    f3d(i,jjm1,k), &
                     x(iim2), x(iim1), x(i), z(kkm2), z(kkm1), z(k), x_u, z_u)
!
                  rad=sqrt(x_u**2+y_u**2+z_u**2)
                  write(1,'(2e20.8,2a10,i10)') rad, f_u, 'u', 'xz', imask3d(i,j,k)
!
               elseif(dels_yzu.eq.dels_u) then
!intersection with y-z plane on level i-alpha
                  x_u = x(i-alpha)
                  y_u = y(j) - dels_u*n_y
                  z_u = z(k) - dels_u*n_z
!
                  iim1=i-alpha
                  jjm2=j-2*beta
                  jjm1=j-beta
                  kkm2=k-2*gamma
                  kkm1=k-gamma
!
!                  f_u = interpol2d_4p_lin(f3d(iim1,jjm1,kkm1), f3d(iim1,j,kkm1), &
!                                          f3d(iim1,jjm1,k),    f3d(iim1,j,k), &
!                                          y(jjm1), y(j), z(kkm1), z(k), y_u, z_u)
                  f_u = interpol2d_9p_quad(f3d(iim1,jjm2,kkm2), f3d(iim1,jjm1,kkm2), f3d(iim1,j,kkm2), &
                     f3d(iim1,jjm2,kkm1), f3d(iim1,jjm1,kkm1), f3d(iim1,j,kkm1), &
                     f3d(iim1,jjm2,k),    f3d(iim1,jjm1,k),    f3d(iim1,j,k), &
                     y(jjm2), y(jjm1), y(j), z(kkm2), z(kkm1), z(k), y_u, z_u)
!
                  rad=sqrt(x_u**2+y_u**2+z_u**2)
                  write(1,'(2e20.8,2a10,i10)') rad, f_u, 'u', 'yz', imask3d(i,j,k)
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
                  x_d = x(i) + dels_d*n_x
                  y_d = y(j) + dels_d*n_y
                  z_d = z(k+gamma)
!
                  iim1=i-alpha
                  iip1=i+alpha
                  jjm1=j-beta
                  jjp1=j+beta
                  kkp1=k+gamma
!
!                  f_d = interpol2d_4p_lin(f3d(i,j,  kkp1), f3d(iip1,j,  kkp1), &
!                                          f3d(i,jjp1,kkp1), f3d(iip1,jjp1,kkp1), &
!                                          x(i), x(iip1), y(j), y(jjp1), x_d, y_d)
                  f_d = interpol2d_9p_quad(f3d(iim1,jjm1,kkp1), f3d(i,jjm1,kkp1), f3d(iip1,jjm1,kkp1), &
                     f3d(iim1,j,   kkp1), f3d(i,j,   kkp1), f3d(iip1,j,   kkp1), &
                     f3d(iim1,jjp1,kkp1), f3d(i,jjp1,kkp1), f3d(iip1,jjp1,kkp1), &
                     x(iim1), x(i), x(iip1), y(jjm1), y(j), y(jjp1), x_d, y_d)!
                  rad=sqrt(x_d**2+y_d**2+z_d**2)
                  write(1,'(2e20.8,2a10,i10)') rad, f_d, 'd', 'xy', imask3d(i,j,k)
!
               elseif(dels_xzd.eq.dels_d) then
!intersection with x-z plane on level j+beta
                  x_d = x(i) + dels_d*n_x
                  y_d = y(j+beta)
                  z_d = z(k) + dels_d*n_z
!
                  iim1=i-alpha
                  iip1=i+alpha
                  jjp1=j+beta
                  kkm1=k-gamma
                  kkp1=k+gamma
!
!                  f_d = interpol2d_4p_lin(f3d(i,jjp1,k),    f3d(iip1,jjp1,k), &
!                                          f3d(i,jjp1,kkp1), f3d(iip1,jjp1,kkp1), &
!                                          x(i), x(iip1), z(k), z(kkp1), x_d, z_d)
                  f_d = interpol2d_9p_quad(f3d(iim1,jjp1,kkm1), f3d(i,jjp1,kkm1), f3d(iip1,jjp1,kkm1), &
                     f3d(iim1,jjp1,k),    f3d(i,jjp1,k),    f3d(iip1,jjp1,k), &
                     f3d(iim1,jjp1,kkp1), f3d(i,jjp1,kkp1), f3d(iip1,jjp1,kkp1), &
                     x(iim1), x(i), x(iip1), z(kkm1), z(k), z(kkp1), x_d, z_d)
                  rad=sqrt(x_d**2+y_d**2+z_d**2)
                  write(1,'(2e20.8,2a10,i10)') rad, f_d, 'd', 'xz', imask3d(i,j,k)
!
               elseif(dels_yzd.eq.dels_d) then
!intersection with y-z plane on level i+alpha
                  x_d = x(i+alpha)
                  y_d = y(j) + dels_d*n_y
                  z_d = z(k) + dels_d*n_z
!
                  iip1=i+alpha
                  jjm1=j-beta
                  jjp1=j+beta
                  kkm1=k-gamma
                  kkp1=k+gamma
!
!                  f_d = interpol2d_4p_lin(f3d(iip1,j,k),    f3d(iip1,jjp1,k), &
!                                          f3d(iip1,j,kkp1), f3d(iip1,jjp1,kkp1), &
!                                          y(j), y(jjp1), z(k), z(kkp1), y_d, z_d)
                  f_d = interpol2d_9p_quad(f3d(iip1,jjm1,kkm1), f3d(iip1,j,kkm1), f3d(iip1,jjp1,kkm1), &
                     f3d(iip1,jjm1,k),    f3d(iip1,j,k),    f3d(iip1,jjp1,k), &
                     f3d(iip1,jjm1,kkp1), f3d(iip1,j,kkp1), f3d(iip1,jjp1,kkp1), &
                     y(jjm1), y(j), y(jjp1), z(kkm1), z(k), z(kkp1), y_d, z_d)
                  rad=sqrt(x_d**2+y_d**2+z_d**2)
                  write(1,'(2e20.8,2a10,i10)') rad, f_d, 'd', 'yz', imask3d(i,j,k)
!
               else
                  write(*,'(4es20.8,l4)') dels_d, dels_xzd, dels_xyd, dels_yzd, dels_d.eq.dels_xzd
                  stop 'error in fsc_cont3d: invalid dels_d'
               endif
!
!***********special treatment: point lies directly on boundary**********
!
             case(3)
               mueff=n_x*x(i)+n_y*y(j)+n_z*z(k)
               if(mueff.ge.0.d0) then
!no interpolations required
               else
!same interpolation as in (1)
!
!calculate distance of upwind point to (previous) xy-plane, xz-plane and yz-plane
                  dels_xyu=(z(k)-z(k-gamma))/n_z
                  dels_xzu=(y(j)-y(j-beta))/n_y
                  dels_yzu=(x(i)-x(i-alpha))/n_x
                  dels_u=min(dels_xyu,dels_xzu,dels_yzu)
!
!calculate distance of downwind point to (next) xy-plane, xz-plane and yz-plane
                  dels_xyd=(z(k+gamma)-z(k))/n_z
                  dels_xzd=(y(j+beta)-y(j))/n_y
                  dels_yzd=(x(i+alpha)-x(i))/n_x
                  dels_d=min(dels_xyd,dels_xzd,dels_yzd)
!
!----------------------------upwind point-------------------------------
!
                  if(dels_xyu.eq.dels_u) then
!intersection with x-y plane on level k-gamma
                     x_u = x(i) - dels_u*n_x
                     y_u = y(j) - dels_u*n_y
                     z_u = z(k-gamma)
!
                     iim2=i-2*alpha
                     iim1=i-alpha
                     jjm2=j-2*beta
                     jjm1=j-beta
                     kkm1=k-gamma
!
!                     f_u = interpol2d_4p_lin(f3d(iim1,jjm1,kkm1), f3d(i,jjm1,kkm1), &
!                                             f3d(iim1,j,   kkm1), f3d(i,j,   kkm1), &
!                                             x(iim1), x(i), y(jjm1), y(j), x_u, y_u)
                     f_u = interpol2d_9p_quad(f3d(iim2,jjm2,kkm1), f3d(iim1,jjm2,kkm1), f3d(i,jjm2,kkm1), &
                        f3d(iim2,jjm1,kkm1), f3d(iim1,jjm1,kkm1), f3d(i,jjm1,kkm1), &
                        f3d(iim2,j,   kkm1), f3d(iim1,j,   kkm1), f3d(i,j,   kkm1), &
                        x(iim2), x(iim1), x(i), y(jjm2), y(jjm1), y(j), x_u, y_u)
                     rad=sqrt(x_u**2+y_u**2+z_u**2)
                     write(1,'(2e20.8,2a10,i10)') rad, f_u, 'u', 'xy', imask3d(i,j,k)
!
                  elseif(dels_xzu.eq.dels_u) then
!intersection with x-z plane on level j-beta
                     x_u = x(i) - dels_u*n_x
                     y_u = y(j-beta)
                     z_u = z(k) - dels_u*n_z
!
                     iim2=i-2*alpha
                     iim1=i-alpha
                     jjm1=j-beta
                     kkm2=k-2*gamma
                     kkm1=k-gamma
!
!                     f_u = interpol2d_4p_lin(f3d(iim1,jjm1,kkm1), f3d(i,jjm1,kkm1), &
!                                             f3d(iim1,jjm1,k),    f3d(i,jjm1,k), &
!                                             x(iim1), x(i), z(kkm1), z(k), x_u, z_u)
                     f_u = interpol2d_9p_quad(f3d(iim2,jjm1,kkm2), f3d(iim1,jjm1,kkm2), f3d(i,jjm1,kkm2), &
                        f3d(iim2,jjm1,kkm1), f3d(iim1,jjm1,kkm1), f3d(i,jjm1,kkm1), &
                        f3d(iim2,jjm1,k),    f3d(iim1,jjm1,k),    f3d(i,jjm1,k), &
                        x(iim2), x(iim1), x(i), z(kkm2), z(kkm1), z(k), x_u, z_u)
                     rad=sqrt(x_u**2+y_u**2+z_u**2)
                     write(1,'(2e20.8,2a10,i10)') rad, f_u, 'u', 'xz', imask3d(i,j,k)
!
                  elseif(dels_yzu.eq.dels_u) then
!intersection with y-z plane on level i-alpha
                     x_u = x(i-alpha)
                     y_u = y(j) - dels_u*n_y
                     z_u = z(k) - dels_u*n_z
!
                     iim1=i-alpha
                     jjm2=j-2*beta
                     jjm1=j-beta
                     kkm2=k-2*gamma
                     kkm1=k-gamma
!
!                     f_u = interpol2d_4p_lin(f3d(iim1,jjm1,kkm1), f3d(iim1,j,kkm1), &
!                                             f3d(iim1,jjm1,k),    f3d(iim1,j,k), &
!                                             y(jjm1), y(j), z(kkm1), z(k), y_u, z_u)
                     f_u = interpol2d_9p_quad(f3d(iim1,jjm2,kkm2), f3d(iim1,jjm1,kkm2), f3d(iim1,j,kkm2), &
                        f3d(iim1,jjm2,kkm1), f3d(iim1,jjm1,kkm1), f3d(iim1,j,kkm1), &
                        f3d(iim1,jjm2,k),    f3d(iim1,jjm1,k),    f3d(iim1,j,k), &
                        y(jjm2), y(jjm1), y(j), z(kkm2), z(kkm1), z(k), y_u, z_u)
                     rad=sqrt(x_u**2+y_u**2+z_u**2)
                     write(1,'(2e20.8,2a10,i10)') rad, f_u, 'u', 'yz', imask3d(i,j,k)
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
                     x_d = x(i) + dels_d*n_x
                     y_d = y(j) + dels_d*n_y
                     z_d = z(k+gamma)
!
                     iim1=i-alpha
                     iip1=i+alpha
                     jjm1=j-beta
                     jjp1=j+beta
                     kkp1=k+gamma
!
!                     f_d = interpol2d_4p_lin(f3d(i,j,  kkp1), f3d(iip1,j,  kkp1), &
!                                             f3d(i,jjp1,kkp1), f3d(iip1,jjp1,kkp1), &
!                                             x(i), x(iip1), y(j), y(jjp1), x_d, y_d)
                     f_d = interpol2d_9p_quad(f3d(iim1,jjm1,kkp1), f3d(i,jjm1,kkp1), f3d(iip1,jjm1,kkp1), &
                        f3d(iim1,j,   kkp1), f3d(i,j,   kkp1), f3d(iip1,j,   kkp1), &
                        f3d(iim1,jjp1,kkp1), f3d(i,jjp1,kkp1), f3d(iip1,jjp1,kkp1), &
                        x(iim1), x(i), x(iip1), y(jjm1), y(j), y(jjp1), x_d, y_d)
                     rad=sqrt(x_d**2+y_d**2+z_d**2)
                     write(1,'(2e20.8,2a10,i10)') rad, f_d, 'd', 'xy', imask3d(i,j,k)
!
                  elseif(dels_xzd.eq.dels_d) then
!intersection with x-z plane on level j+beta
                     x_d = x(i) + dels_d*n_x
                     y_d = y(j+beta)
                     z_d = z(k) + dels_d*n_z
!
                     iim1=i-alpha
                     iip1=i+alpha
                     jjp1=j+beta
                     kkm1=k-gamma
                     kkp1=k+gamma
!
!                     f_d = interpol2d_4p_lin(f3d(i,jjp1,k),    f3d(iip1,jjp1,k), &
!                                             f3d(i,jjp1,kkp1), f3d(iip1,jjp1,kkp1), &
!                                             x(i), x(iip1), z(k), z(kkp1), x_d, z_d)
                     f_d = interpol2d_9p_quad(f3d(iim1,jjp1,kkm1), f3d(i,jjp1,kkm1), f3d(iip1,jjp1,kkm1), &
                        f3d(iim1,jjp1,k),    f3d(i,jjp1,k),    f3d(iip1,jjp1,k), &
                        f3d(iim1,jjp1,kkp1), f3d(i,jjp1,kkp1), f3d(iip1,jjp1,kkp1), &
                        x(iim1), x(i), x(iip1), z(kkm1), z(k), z(kkp1), x_d, z_d)
                     rad=sqrt(x_d**2+y_d**2+z_d**2)
                     write(1,'(2e20.8,2a10,i10)') rad, f_d, 'd', 'xz', imask3d(i,j,k)
!
                  elseif(dels_yzd.eq.dels_d) then
!intersection with y-z plane on level i+alpha
                     x_d = x(i+alpha)
                     y_d = y(j) + dels_d*n_y
                     z_d = z(k) + dels_d*n_z
!
                     iip1=i+alpha
                     jjm1=j-beta
                     jjp1=j+beta
                     kkm1=k-gamma
                     kkp1=k+gamma
!
!                     f_d = interpol2d_4p_lin(f3d(iip1,j,k),    f3d(iip1,jjp1,k), &
!                                             f3d(iip1,j,kkp1), f3d(iip1,jjp1,kkp1), &
!                                             y(j), y(jjp1), z(k), z(kkp1), y_d, z_d)
                     f_d = interpol2d_9p_quad(f3d(iip1,jjm1,kkm1), f3d(iip1,j,kkm1), f3d(iip1,jjp1,kkm1), &
                        f3d(iip1,jjm1,k),    f3d(iip1,j,k),    f3d(iip1,jjp1,k), &
                        f3d(iip1,jjm1,kkp1), f3d(iip1,j,kkp1), f3d(iip1,jjp1,kkp1), &
                        y(jjm1), y(j), y(jjp1), z(kkm1), z(k), z(kkp1), y_d, z_d)
                     rad=sqrt(x_d**2+y_d**2+z_d**2)
                     write(1,'(2e20.8,2a10,i10)') rad, f_d, 'd', 'yz', imask3d(i,j,k)
!
                  else
                     write(*,'(4es20.8,l4)') dels_d, dels_xzd, dels_xyd, dels_yzd, dels_d.eq.dels_xzd
                     stop 'error in fsc_cont3d: invalid dels_d'
                  endif
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


   close(1)
   stop
!
!
end subroutine interpolation_test2d
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine functions_test(x, xim1, xi, yim1, yi, f1, f2, f3, f4, f5, f6, f7, f8)
!
   use prog_type
!
   implicit none
!
! ... arguments
   real(dp), intent(in) :: x, xim1, xi, yim1, yi
   real(dp) :: f1, f2, f3, f4, f5, f6, f7, f8
!
! ... local scalars
   real(dp) :: m, t, dely, delx
!
!function 1: linear; f(x)=m*x+t
   f1 = yim1 + (yi-yim1)/(xi-xim1)*(x-xim1)
!
!function 2: eqponential increase; f(x)=exp(-m*x+t)
   m=log(yim1/yi)/(xi-xim1)
   t=log(yi)+m*xi
   f2 = exp(-m*x+t)
!
!function 3: step function at center
   dely=yi-yim1
   delx=xi-xim1
   if(x.le.xim1+delx/2.d0) f3=-dely/2.d0
   if(x.gt.xim1+delx/2.d0) f3=dely/2.d0
!
!function 4: step function near left boundary
   if(x.le.xim1+delx/4.d0) f4=-dely/2.d0
   if(x.gt.xim1+delx/4.d0) f4=dely/2.d0
!
!function 5: step function near right boundary
   if(x.le.xim1+delx*3.d0/4.d0) f5=-dely/2.d0
   if(x.gt.xim1+delx*3.d0/4.d0) f5=dely/2.d0
!
!function 6: quadratic function with monotonically increasing interval
   f6=(x-xim1)**2
!
!function 7: quadratic function with monotonically decreasing interval
   f7=(x-xi)**2
!
!function 8: quadratic function, non-monotonic in interval
   f8=(x-(xi-xi/4.d0+xim1/4.d0))**2

end subroutine functions_test
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine profile_test
!
!--------------------calculates doppler-profile-------------------------
!---------------------check for normalization---------------------------
!---------including output to files (grid needs to be small enough)-----
!
   use prog_type
   use dime3d, only: ndxmax, ndymax, ndzmax, velx3d, vely3d, velz3d, vth3d, imask_totreg3d
   use angles, only: dim_omega, n_x, n_y, n_z
   use freq, only: nxobs, nodes_xobs, weight_xobs
   use mod_directories, only: output_dir_test
   use params_input, only: vth_fiducial
!
   implicit none
!
! ... local scalars
   integer(i4b) :: i,j,k,l,m,n
   integer(i4b) :: indx_x, indx_y, indx_z, indx_omega, indx_xobs, indx_xt1, indx_yt1, indx_zt1
   integer(i4b) :: err
   real(dp) :: integ, val_profile, velz
!
! ... local arrays
   real(dp), dimension(:,:,:,:), allocatable :: prbar_angdep
   real(dp), dimension(:,:,:), allocatable :: pr_angdep
!
!
!*******************checking profile: output to terminal****************
!
   write(*,*) '---------------------checking profile (output to file)-------------------------'
   write(*,*)
!
!-------------------at a given point in the atmosphere------------------
!
   indx_x=ndxmax/2+1
   indx_y=ndymax/2+1
   indx_z=ndzmax-1
!
   write(*,'(a15, 3(i5))') 'at (i,j,k)', indx_x, indx_y, indx_z
   write(*,'(a15, e15.8)') 'vel_x:', velx3d(indx_x,indx_y,indx_z)
   write(*,'(a15, e15.8)') 'vel_y:', vely3d(indx_x,indx_y,indx_z)
   write(*,'(a15, e15.8)') 'vel_z:', velz3d(indx_x,indx_y,indx_z)
   write(*,*)
!
!------------------for a given angle pair mu, phi-----------------------
!
   indx_omega=dim_omega/2
!
   write(*,'(a20, 2(i5))') 'for angles', indx_omega
   write(*,'(a20, 3(e20.8))') 'n_x, n_y, n_z', n_x(indx_omega), n_y(indx_omega), n_z(indx_omega)
   write(*,*)
!
!-----------------------frequency integration---------------------------
!
   integ=0.d0
!
   write(*,'(a10, 4(a20))') 'freq-point', 'xobs', 'weight', 'profile', 'integral'
   do i=1, nxobs
      velz=n_x(indx_omega)*velx3d(indx_x,indx_y,indx_z) + &
         n_y(indx_omega)*vely3d(indx_x,indx_y,indx_z) + &
         n_z(indx_omega)*velz3d(indx_x,indx_y,indx_z)
      call calc_phinorm(velz, vth3d(indx_x,indx_y,indx_z), vth_fiducial, nodes_xobs(i), val_profile)
      integ=integ+val_profile*weight_xobs(i)
      write(*,'(i10, 4(e20.8))') i,  nodes_xobs(i), weight_xobs(i), val_profile, integ
   enddo
!
!******************checking profile (for idl routines)******************
!
!-------------test1: profile on a given axis of the atmosphere----------
!
   indx_yt1=ndymax/2+1
   indx_zt1=ndzmax/2+1
!
!
   allocate(pr_angdep(ndxmax, nxobs, dim_omega), stat=err)
   if(err.ne.0) stop 'allocation error profile_test: pr_angdep'
   pr_angdep=0.d0
!
   do indx_xt1=1, ndxmax
      if(imask_totreg3d(indx_xt1,indx_yt1,indx_zt1).ne.0) then
         do indx_xobs=1, nxobs
            do indx_omega=1, dim_omega
               velz=n_x(indx_omega)*velx3d(indx_xt1,indx_yt1,indx_zt1) + &
                  n_y(indx_omega)*vely3d(indx_xt1,indx_yt1,indx_zt1) + &
                  n_z(indx_omega)*velz3d(indx_xt1,indx_yt1,indx_zt1)
               call calc_phinorm(velz, vth3d(indx_xt1, indx_yt1, indx_zt1), vth_fiducial, nodes_xobs(indx_xobs), val_profile)
               pr_angdep(indx_xt1, indx_xobs, indx_omega)=val_profile
            enddo
         enddo
      endif
   enddo
!
!-------test2: frequency integrated as function of pos and angle--------
!
   allocate(prbar_angdep(ndxmax, ndymax, ndzmax, indx_omega), stat=err)
   if(err.ne.0) stop 'allocation error profile_test: prbar_angdep'
   prbar_angdep=0.d0
!
   do indx_z=1, ndzmax
      do indx_y=1, ndymax
         do indx_x=1, ndxmax
            if(imask_totreg3d(indx_x,indx_y,indx_z).ne.0) then
               do indx_omega=1, dim_omega
                  integ=0.d0
                  velz=n_x(indx_omega)*velx3d(indx_x,indx_y,indx_z) + &
                     n_y(indx_omega)*vely3d(indx_x,indx_y,indx_z) + &
                     n_z(indx_omega)*velz3d(indx_x,indx_y,indx_z)
                  do indx_xobs=1,nxobs
                     call calc_phinorm(velz, vth3d(indx_x,indx_y,indx_z), vth_fiducial, nodes_xobs(indx_xobs), val_profile)
                     integ=integ+val_profile*weight_xobs(indx_xobs)
                  enddo
                  prbar_angdep(indx_x, indx_y, indx_z, indx_omega)=integ
               enddo
            endif
         enddo
      enddo
   enddo
!
!--------------------------output to files------------------------------
!
!setting indx_(x,y,z)t1 to zero if test 1 is calculated along that axis
   indx_xt1=0
!!
   open(1, file=trim(output_dir_test)//'/tprof_in1.dat', form='formatted')
   write(1,'(3(a6))') 'indx_x', 'indx_y', 'indx_z'
   write(1,'(3(i6))') indx_xt1, indx_yt1, indx_zt1
   close(1)
!
   open(1, file=trim(output_dir_test)//'/tprof_t1.dat', form='unformatted')
   write(1) pr_angdep
   close(1)
!
   open(1, file=trim(output_dir_test)//'/tprof_t2.dat', form='unformatted')
   write(1) prbar_angdep
   close(1)
!
!
!
end subroutine profile_test
