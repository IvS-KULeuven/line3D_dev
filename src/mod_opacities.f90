module mod_opacities

  use prog_type
  use fund_const, only: pi, cgs_mp, sigmae
  use mod_interp1d
  use mod_interp2d, only: interpol2d_4p_lin
  use mod_lte, only: rho_lte, temp_lte, nlower_lte, nrho_lte, ntemp_lte
  use mod_iline, only: gf, gl

  implicit none

contains
  !
  !-----------------------------------------------------------------------  
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !
  function opac_thomson(yhe, hei, rho, kcont)

    !calculate thomson opacity, assuming fully ionized hydrogen
    !
    !  INPUT: yhe:     helium abundance by number
    !         hei:     number of free electrons per helium atom (=2 for full ionization)
    !         rho:     density in cgs
    !         kcont:   multiplication factor
    !  OUTPUT: thomson opacity in 1/cm
    
    
    ! ... arguments
    real(dp), intent(in) :: yhe, hei, rho, kcont
    real(dp) :: opac_thomson

    ! ... local scalars
    real(dp) :: c1, c2, ne

    c1=(1.d0+4.d0*yhe)*cgs_mp
    c2=(1.d0+hei*yhe)/c1
    ne=c2*rho
    opac_thomson=sigmae*ne*kcont

  end function opac_thomson
  !
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !
  function opac_thomson2(sr, vinf, mdot, yhe, hei, rho, kcont, alpha, rad)

    !same as opac_thomson, however, opacity multiplied by additional factor (v(r)/vinf)**alpha
    !
    !  INPUT: sr:      stellar radius in cgs
    !         vinf:    terminal velocity in cgs
    !         mdot:    mass-loss rate in cgs
    !         rad:     radial coordinate in cgs
    !         rho:     density in cgs
    !  OUTPUT: opacity in 1/cm

    ! ... arguments
    real(dp), intent(in) :: sr, vinf, mdot, yhe, hei, rho, kcont, alpha, rad
    real(dp) :: opac_thomson2

    ! ... local scalars
    real(dp) :: vr, c1, c2, ne

    vr = mdot/4.d0/pi/rad/rad/rho

    c1=(1.d0+4.d0*yhe)*cgs_mp
    c2=(1.d0+hei*yhe)/c1
    ne=c2*rho
    opac_thomson2=sigmae*ne*kcont*(vr/vinf)**alpha

    !open(1, file='outputFILES_DEBUG/data_test.dat', access='append')
    !   write(1,*), rad/sr, vr/vinf, opac_thomson2
    !close(1)

  end function opac_thomson2
  !
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !
  function opac_thomson3(yhe, hei, hi, rho, kcont)

    !calculate thomson opacity
    !
    !  INPUT: yhe:     helium abundance by number
    !         hei:     number of free electrons per helium atom (=2 for full ionization)
    !         hi:      number of free electrons per hydrogen atom (=1 for full ionization)
    !         rho:     density in cgs
    !         kcont:   multiplication factor
    !  OUTPUT: thomson opacity in 1/cm

    ! ... arguments
    real(dp), intent(in) :: yhe, hei, hi, rho, kcont
    real(dp) :: opac_thomson3
    
    ! ... local scalars
    real(dp) :: c1, c2, ne

    c1=(1.d0+4.d0*yhe)*cgs_mp
    c2=(hi+hei*yhe)/c1
    ne=c2*rho
    opac_thomson3=sigmae*ne*kcont

  end function opac_thomson3
  !
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !
  function opac_opal(kcont, yhe, hei, rho, temp, nrho, ntemp, rho_opal, temp_opal, kappa_opal)
    !
    !all temperatures and densities in log-space
    !kappa-table in log-space
    !
    ! ... arguments
    integer(i4b), intent(in) :: nrho, ntemp
    real(dp), intent(in) :: kcont, rho, temp, yhe, hei
    real(dp), dimension(nrho), intent(in) :: rho_opal
    real(dp), dimension(ntemp), intent(in) :: temp_opal
    real(dp), dimension(ntemp,nrho) :: kappa_opal
    real(dp) :: opac_opal
    !
    ! ... local scalars
    integer(i4b) :: i, iim2, iim1, ii, iip1, jjm2, jjm1, jj, jjp1
    real(dp) :: kappa
    real(dp) :: rho_local, t_local
    !
    ! ... local functions
    !
    !
    t_local=10.**temp
    rho_local = 10.d0**rho/(t_local*1.d-6)**3
    rho_local = log10(rho_local)
    t_local=temp
    !
    !
    !use thomson opacity if temperatures too large
    if(t_local.gt.maxval(temp_opal)) then
       opac_opal = opac_thomson(yhe, hei, 10.d0**rho, kcont)
       write(*,*) opac_opal/10.d0**rho
    else
       !
       !adapt rho-local values if outside of opacity table
       do i=1, 1000
          if(rho_local.lt.minval(rho_opal)) then
             rho_local = rho_local + 0.5
          elseif(rho_local.gt.maxval(rho_opal)) then
             rho_local = rho_local - 0.5
          else
             exit
          endif
       enddo
       !
       !adapt t_local values if outside of opacity table
       do i=1, 1000
          if(t_local.lt.minval(temp_opal)) then
             t_local = t_local + 0.05
             !   elseif(t_local.gt.maxval(temp_opal)) then
             !      t_local = t_local - 0.05
          else
             exit
          endif
       enddo
       !
       !
       !
       call find_index(temp, temp_opal, ntemp, iim2, iim1, ii, iip1)
       call find_index(rho_local, rho_opal, nrho, jjm2, jjm1, jj, jjp1)

       kappa = interpol2d_4p_lin(kappa_opal(iim1, jjm1), kappa_opal(ii,jjm1), &
                                 kappa_opal(iim1, jj),   kappa_opal(ii,jj), &
                                 temp_opal(iim1), temp_opal(ii), &
                                 rho_opal(jjm1), rho_opal(jj), temp, rho_local)
       !
       !
       !write(*,*) rho, temp, rho_local, rho_local_orig, t_local, kappa, 10.d0**kappa
       !
       opac_opal = kcont * 10.d0**kappa * 10.d0**rho
       !
    endif
    !
  end function opac_opal
  !
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !
  function opalbar_model_kline(yhe, hei, rho, kline)
    !
    !  calculates line opacity according to line-strength parameter
    !          (see definition in puls & springmann 2000)
    !
    !  input: yhe   helium abundance by number
    !         hei   number of free electrons per helium atom
    !         rho   density in cgs
    !         kline  line-strength parameter
    !  output: frequency integrated line opacity in frequency-space (in 1/cm)
    !
    ! ... arguments
    real(dp), intent(in) :: yhe, hei, rho, kline
    real(dp) :: opalbar_model_kline
    !
    ! ... local scalars
    real(dp) :: c1, c2, ne, chi_thomson
    !
    ! ... local arrays
    !
    !calculate thomson-opacity
    c1=(1.d0+4.d0*yhe)*cgs_mp
    c2=(1.d0+hei*yhe)/c1
    ne=c2*rho
    chi_thomson=sigmae*ne
    !
    !calculate line opacity
    opalbar_model_kline = kline*chi_thomson
    !
  end function opalbar_model_kline
  !
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !
  function opalbar_model_kline2(yhe, hei, rho, kline)
    !
    !  calculates line opacity according to line-strength parameter
    !               for density-squared opacities
    !                  (opalbar = kline*rho^2)
    !
    !  input: yhe   helium abundance by number
    !         hei   number of free electrons per helium atom
    !         rho   density in cgs
    !         kline  line-strength parameter
    !  output: frequency integrated line opacity in frequency-space (in 1/cm)
    !
    !
    ! ... arguments
    real(dp), intent(in) :: yhe, hei, rho, kline
    real(dp) :: opalbar_model_kline2
    !
    ! ... local scalars
    real(dp) :: c1, c2, ne, chi_thomson
    !
    ! ... local arrays
    !
    !calculate thomson-opacity
    !c1=(1.d0+4.d0*yhe)*cgs_mp
    !c2=(1.d0+hei*yhe)/c1
    !ne=c2*rho
    !chi_thomson=sigmae*ne
    !
    !calculate line opacity
    opalbar_model_kline2 = kline*rho**2 !chi_thomson
    !
  end function opalbar_model_kline2
  !
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !
  function opalbar_model_hamann(sr, vinf, mdot, kappa0, alpha, vth_fid, rad, rho)
    !
    !  calculates line opacity according to definition by hamann 1981
    !(see also Puls,Owocki,Fullerton 1993, Sundqvist Diss (Appendix))
    !
    !  INPUT: sr:      stellar radius in cgs
    !         vinf:    terminal velocity in cgs
    !         mdot:    mass-loss rate in cgs
    !         kappa0:  line-strength parameter (dimensionless) as defined by Hamann
    !         alpha:   see Hamann
    !         vth_fid: fiducial thermal velocity in cgs
    !         rad:     radial coordinate in cgs
    !         rho:     density in cgs
    !  OUTPUT: frequency integrated line opacity in frequency space (in 1/cm)
    !
    !
    ! ... arguments
    real(dp), intent(in) :: sr, vinf, mdot, kappa0, vth_fid, rho, rad, alpha
    real(dp) :: opalbar_model_hamann
    !
    ! ... local scalars
    !
    ! ... local arrays
    !
    !-----------------------------------------------------------------------
    !
    !old version: without alpha
    !in own units
    !opalbar_model_hamann = 4.d0*pi*sr*vinf*vinf*kappa0*rho/mdot/vth_fid
    !
    !new version: including alpha: opalbar=kappa0*v(r)^alpha / r^2 / v(r) (with v in vmax, r in rstar)
    !in own units
    opalbar_model_hamann = 4.d0*pi*sr*vinf*vinf*kappa0*rho/mdot/vth_fid * (mdot/4.d0/pi/rad/rad/rho/vinf)**alpha
    !
    !
  end function opalbar_model_hamann
  !
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !
  subroutine depcoeff_petrenz(velr, b2, b3)
    !
    !--calculate departure coefficients according to petrenz & puls 1995----
    !-----------------from puls et al 1995, eq 45---------------------------
    !
    !input:  v_r in units of v_inf
    !output: departure coefficients for h-alpha: b2, b3
    !
    !
    ! ... arguments
    real(dp), intent(in) :: velr
    real(dp) :: b2, b3
    !
    ! ... local scalars
    real(dp), parameter :: b2_in=1.5d0, b3_in=1.2d0, &
                           b2_min=1.2d0, b3_min=1.1d0, &
                           b2_inf=1.3d0, b3_inf=1.1d0
    !
    ! ... local functions
    !
    if(velr.lt.0.01d0) then
       b2 = 1.d0 + (b2_in - 1.d0)*velr/0.01d0
       b3 = 0.9d0 + (b3_in - .9d0)*velr/0.01d0
    else if (velr.lt.0.1d0) then
       b2 = b2_in + (b2_min-b2_in)*(velr-0.01d0)/0.09d0
       b3 = b3_in + (b3_min-b3_in)*(velr-0.01d0)/0.09d0
    else if (velr.le.1.d0) then
       b2 = b2_min + (b2_inf-b2_min)*(velr-0.1d0)/(0.9d0)
       b3 = b3_min + (b3_inf-b3_min)*(velr-0.1d0)/(0.9d0)
    else
       stop 'error in depcoeff_petrenz: wrong velocity range'
    endif
    !
    !
  end subroutine depcoeff_petrenz
  !
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !
  function opalbar_petrenz(sr, yhe, hei, temp, vth_fiducial, xnue0, b2, b3, rho)
    !
    !-----------------according to petrenz & puls 1995----------------------
    !
    ! ... arguments
    real(dp), intent(in) :: rho, temp, b2, b3, yhe, hei, sr, vth_fiducial, xnue0
    real(dp) :: opalbar_petrenz
    !
    ! ... local scalars
    real(dp), parameter :: c1=2.07d-16, gf=5.1256d0
    real(dp) ::deldop_fiducial, dum1
    !
    ! ... local arrays
    !
    dum1=(1.d0+yhe*hei)/(1.d0+4.d0*yhe)**2
    dum1=dum1*rho*rho*temp**(-1.5d0)
    dum1=dum1*c1*gf*pi*cgs_e*cgs_e/cgs_me/cgs_clight/cgs_mp/cgs_mp
    dum1=dum1*(b2*exp(3.954d4/temp) - b3*exp(1.753d4/temp))
    !
    !so far in cgs, now in own units
    !
    !finally need to divide this by deldop_fiducial
    deldop_fiducial=xnue0*vth_fiducial/cgs_clight
    opalbar_petrenz=dum1*sr/deldop_fiducial
    !
    !
  end function opalbar_petrenz
  !
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !
  function opalbar_halpha(sr, yhe, hei, temp, vth_fiducial, xnue0, b2, b3, rho)
    !
    !assuming complete H-ionization and recombination line
    !all input in cgs
    !
    !-----------------see personal notes on opacities-----------------------
    !
    !
    ! ... arguments
    real(dp), intent(in) :: rho, temp, b2, b3, yhe, hei, sr, vth_fiducial, xnue0
    real(dp) :: opalbar_halpha
    !
    ! ... local scalars
    real(dp), parameter :: energy_leveln=cgs_clight*cgs_planck*109678.77d0, &
                           energy_level2=cgs_clight*cgs_planck*82259.158d0, &
                           energy_level3=cgs_clight*cgs_planck*97492.305d0, &
                           c1 = half*(cgs_planck**2/two/pi/cgs_me/cgs_kb)**(three/two), &
                           c2 = pi*cgs_e**2 / cgs_me/cgs_clight /cgs_mp**2, &
                           gf=5.1286d0
    real(dp) ::deldop_fiducial, dum1
    !
    ! ... local arrays
    !
    if(temp.lt.6.d3) stop 'error in opalbar_halpha: t<6000 does not make sense for within this formulation'

    dum1=(one+yhe*hei)/(one+four*yhe)**2
    dum1=dum1*rho*rho*temp**(-three/two)
    dum1=dum1*c1*c2*gf
    dum1=dum1*(b2*exp((energy_leveln-energy_level2)/cgs_kb/temp) - &
         b3*exp((energy_leveln-energy_level3)/cgs_kb/temp))
    !
    !write(*,*) rho, temp, yhe, hei, dum1
    !so far in cgs, now in own units
    !
    !finally need to divide this by deldop_fiducial
    deldop_fiducial=xnue0*vth_fiducial/cgs_clight
    opalbar_halpha=dum1*sr/deldop_fiducial
    !
    !
  end function opalbar_halpha
  !
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !
  function opalbar_hbeta(sr, yhe, hei, temp, vth_fiducial, xnue0, b2, b4, rho)
    !
    !assuming complete H-ioniazation and recombination line
    !all input in cgs
    !
    !-----------------see personal notes on opacities-----------------------
    !
    !
    ! ... arguments
    real(dp), intent(in) :: rho, temp, b2, b4, yhe, hei, sr, vth_fiducial, xnue0
    real(dp) :: opalbar_hbeta
    !
    ! ... local scalars
    real(dp), parameter :: energy_leveln=cgs_clight*cgs_planck*109678.77d0, &
                           energy_level2=cgs_clight*cgs_planck*82259.158d0, &
                           energy_level4=cgs_clight*cgs_planck*102832.904d0, &
                           c1 = half*(cgs_planck**2/two/pi/cgs_me/cgs_kb)**(three/two), &
                           c2 = pi*cgs_e**2 / cgs_me/cgs_clight /cgs_mp**2, &
                           gf = 0.95508055
    real(dp) ::deldop_fiducial, dum1
    !
    ! ... local arrays
    !
    dum1=(one+yhe*hei)/(one+four*yhe)**2
    dum1=dum1*rho*rho*temp**(-three/two)
    dum1=dum1*c1*c2*gf
    dum1=dum1*(b2*exp((energy_leveln-energy_level2)/cgs_kb/temp) - &
         b4*exp((energy_leveln-energy_level4)/cgs_kb/temp))
    !
    !write(*,*) rho, temp, yhe, hei, dum1
    !so far in cgs, now in own units
    !
    !finally need to divide this by deldop_fiducial
    deldop_fiducial=xnue0*vth_fiducial/cgs_clight
    opalbar_hbeta=dum1*sr/deldop_fiducial
    !
    !
  end function opalbar_hbeta
  !
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !
  function opalbar_halpha2(sr, yhe, temp, vth_fiducial, xnue0, b2, b3, rho)
    !
    !assuming neutral H, and all elements in ground state
    !all input in cgs
    !
    !-----------------see personal notes on opacities-----------------------
    !
    !
    ! ... arguments
    real(dp), intent(in) :: rho, temp, b2, b3, yhe, sr, vth_fiducial, xnue0
    real(dp) :: opalbar_halpha2
    !
    ! ... local scalars
    real(dp), parameter :: energy_level2=cgs_clight*cgs_planck*82259.158d0, &
                           energy_level3=cgs_clight*cgs_planck*97492.305d0, &
                           f23 = 6.4108d-1, &
                           g1 = 2.d0, g2 = 8.d0, &
                           c1 = pi*cgs_e**2 / cgs_me/cgs_clight /cgs_mp
    real(dp) ::deldop_fiducial, dum1
    !
    ! ... local arrays
    !
    dum1 = c1 * f23 * g2/g1 * rho/(one+four*yhe)
    dum1=dum1*(b2*exp(-energy_level2/cgs_kb/temp) - &
         b3*exp(-energy_level3/cgs_kb/temp))
    !write(*,*) rho, temp, yhe, hei, dum1
    !so far in cgs, now in own units
    !
    !
    !finally need to divide this by deldop_fiducial
    deldop_fiducial=xnue0*vth_fiducial/cgs_clight
    opalbar_halpha2=dum1*sr/deldop_fiducial
    !
    !
  end function opalbar_halpha2
  !
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !
  function opalbar_halpha3(sr, yhe, hei, temp, vth_fiducial, xnue0, b2, b3, rho)
    !
    !calculate hydrogen ionization balance from LTE
    !calculate level populations with departure coefficients
    !all input in cgs
    !
    !-----------------see personal notes on opacities-----------------------
    !
    !
    ! ... arguments
    real(dp), intent(in) :: rho, temp, b2, b3, yhe, hei, sr, vth_fiducial, xnue0
    real(dp) :: opalbar_halpha3
    !
    ! ... atomic data (from NIST)
    integer(i4b), parameter :: nlist = 40
    real(dp), dimension(nlist), parameter :: gvalue= (/ 2.d0, 8.d0, 18.d0, 32.d0, 50.d0, 72.d0, &
         98.d0, 128.d0, 162.d0, 200.d0, 242.d0, 288.d0, 338.d0, 392.d0, 450.d0, 512.d0, 578.d0, &
         648.d0, 722.d0, 800.d0, 882.d0, 968.d0, 1058.d0, 1152.d0, 1250.d0, 1352.d0, 1458.d0, &
         1568.d0, 1682.d0, 1800.d0, 1922.d0, 2048.d0, 2178.d0, 2312.d0, 2450.d0, 2592.d0, 2738.d0, &
         2888.d0, 3042.d0, 3200d0 /)
    real(dp), dimension(nlist), parameter :: energy= (/ 0.d0, 82259.158d0, 97492.304d0, 102823.904d0, &
         105291.657d0, 106632.1681d0, 107440.4508d0, 107965.0568d0, 108324.7253d0, 108581.9945d0, &
         108772.3445d0, 108917.1209d0, 109029.7913d0, 109119.1917d0, 109191.3154d0, 109250.3433d0, &
         109299.2642d0, 109340.2604d0, 109374.9555d0, 109404.5776d0, 109430.0696d0, 109452.1650d0, &
         109471.4416d0, 109488.3592d0, 109503.2875d0, 109516.5267d0, 109528.3223d0, 109538.8768d0, &
         109548.3584d0, 109556.9077d0, 109564.6431d0, 109571.6647d0, 109578.0577d0, 109583.8949d0, &
         109589.2390d0, 109594.1439d0, 109598.6566d0, 109602.8177d0, 109606.6628d0, 109610.2232d0 /) * cgs_clight*cgs_planck
    real(dp), parameter :: energy_ion = 109678.77d0 * cgs_clight*cgs_planck
    real(dp), parameter :: f23 = 6.4108d-1, &
         c1 = half*(cgs_planck**2/two/pi/cgs_me/cgs_kb)**(three/two), &
         c2 = pi*cgs_e**2 / cgs_me/cgs_clight

    ! ... local scalars
    integer(i4b) :: i
    real(dp) :: energy_level2, energy_level3, g2, g3
    real(dp) :: deldop_fiducial, fdum1, fdum2, dum1
    real(dp) :: nh, nhi, nhii, n2, n3, nel, zhi, zhii, q, disc
    !
    ! ... local arrays
    !
    !hydrogen density
    nh = rho/(one+four*yhe)/cgs_mp
    !
    !partition functions
    zhi = zero
    zhii = one
    do i=1, nlist
       zhi = zhi + gvalue(i) * exp(-energy(i)/cgs_kb/temp)
    enddo
    !
    !calculate q value (see personal notes)
    q = temp**(3./2.)/c1 * zhii/zhi * exp(-energy_ion/cgs_kb/temp)
    !
    !calculate electron density
    fdum1 = q - hei*yhe*nh
    fdum2 = -q*nh*(one+hei*yhe)
    disc = fdum1**2 - four*fdum2
    if(disc.lt.zero) stop 'error in opalbar_halpha3: cannot calculate electron density with discriminant < zero'
    nel = (-fdum1 + sqrt(disc))/two

    !write(*,*) temp, nh, q, disc, nel, zhii, zhi
    !
    !calculate hydrogen ionization stages
    nhii = nel - hei*yhe*nh
    nhi = nh - nhii
    !
    !calculate level populations of level two and three
    energy_level2 = energy(2)
    energy_level3 = energy(3)
    g2 = gvalue(2)
    g3 = gvalue(3)
    n2 = b2 * g2 * nhi/zhi * exp(-energy_level2/cgs_kb/temp)
    n3 = b3 * g3 * nhi/zhi * exp(-energy_level3/cgs_kb/temp)

    dum1 = c2*f23*(n2 - g2/g3 * n3)


    !write(*,'(10es20.8)') temp, nel, nh, nhi, nhii, n2, n3, dum1
    !write(*,*) rho, temp, yhe, hei, dum1
    !so far in cgs, now in own units
    !
    !
    !finally need to divide this by deldop_fiducial
    deldop_fiducial=xnue0*vth_fiducial/cgs_clight
    opalbar_halpha3=dum1*sr/deldop_fiducial
    !
    if(opalbar_halpha3.lt.zero) stop 'error in opalbar_halpha3: opacity < zero not allowed'
    !
  end function opalbar_halpha3
  !
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !
  function opalbar_lte_table(sr, vth_fiducial, xnue0, bl, bu, rho, temp)
    !
    !all temperatures and densities in log-space
    !nlower-table in log-space
    !
    !bu, bl: NLTE departure coefficients
    !xnue0: transition frequency
    !gl: lower level degeneracy
    !gf: gf value (NOT in log-space)
    !vth_fiducial: fiducial termal velocity in cm/s
    !sr: length scale (to have opacities in 1/sr on output instead of 1/cm)
    !
    ! ... arguments
    real(dp), intent(in) :: rho, temp, bu, bl, xnue0, vth_fiducial, sr
    real(dp) :: opalbar_lte_table
    !
    ! ... local scalars
    integer(i4b) :: i, iim2, iim1, ii, iip1, jjm2, jjm1, jj, jjp1
    real(dp) :: nlower, deldop_fiducial
    real(dp) :: rho_local, t_local, fdum
    real(dp), parameter :: c1 = pi*cgs_e**2 / cgs_me/cgs_clight
    !
    ! ... local functions
    !
    rho_local = log10(rho)
    t_local = log10(temp)
    !
    !
    if(t_local.gt.maxval(temp_lte).or.t_local.lt.minval(temp_lte)) then
       stop 'error in opalbar_lte_table: temperature out of range'
    endif
    if(rho_local.gt.maxval(rho_lte).or.rho_local.lt.minval(rho_lte)) then
       stop 'error in opalbar_lte_table: density out of range'
    endif    
    !
    call find_index(t_local, temp_lte, ntemp_lte, iim2, iim1, ii, iip1)
    call find_index(rho_local, rho_lte, nrho_lte, jjm2, jjm1, jj, jjp1)

    nlower = interpol2d_4p_lin(nlower_lte(iim1, jjm1), nlower_lte(ii,jjm1), &
                              nlower_lte(iim1, jj),   nlower_lte(ii,jj), &
                              temp_lte(iim1), temp_lte(ii), &
                              rho_lte(jjm1), rho_lte(jj), t_local, rho_local)
    !
    t_local = ten**t_local
    fdum = one - bu/bl * exp(-cgs_planck*xnue0/cgs_kb/t_local)
    
    opalbar_lte_table = c1*gf * bl/gl * ten**nlower * fdum

    !finally need to divide this by deldop_fiducial and in units of 1/sr
    deldop_fiducial=xnue0*vth_fiducial/cgs_clight
    opalbar_lte_table=opalbar_lte_table*sr/deldop_fiducial
    !
  end function opalbar_lte_table
  !
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !
  function get_opalbar(iline, kline, sr, yhe, hei, temp, vth_fiducial, xnue0, bl, bu, rho)
    !
    !get opalbar for a given line transition described by iline
    !opacity law for balmer lines by assuming completely ionized hydrogen
    !
    !input:  sr   stellar radius in cgs
    !        yhe, hei:  helium abundance and number of free electrons per helium
    !        temp:  temperature in K
    !        vth_fiducial:  fiducial thermal velocity in cm/s
    !        xnue0:   transition frequency
    !        bl,bu:   NLTE-departure coefficients for lower and upper level
    !        rho:    density in cgs
    !        kline:   just a scaling factor to increase or decrease opacity (line strength)
    !
    ! ... arguments
    integer(i4b), intent(in) :: iline
    real(dp), intent(in) :: sr, yhe, hei, temp, vth_fiducial, xnue0, bl, bu, rho, kline
    real(dp) :: get_opalbar
    !
    ! ... local functions
    select case(iline)
       case(0)
          get_opalbar = opalbar_lte_table(sr, vth_fiducial, xnue0, bl, bu, rho, temp)*kline
       case(1)
          get_opalbar = opalbar_halpha(sr, yhe, hei, temp, vth_fiducial, xnue0, bl, bu, rho)*kline
       case(2)
          get_opalbar = opalbar_hbeta(sr, yhe, hei, temp, vth_fiducial, xnue0, bl, bu, rho)*kline
       case(10)
          get_opalbar = opalbar_model_kline(yhe, hei, rho, kline)*sr
       case(11)
          get_opalbar = opalbar_model_kline2(yhe, hei, rho, kline)*sr
       case default
          stop 'error in get_opalbar: iline not properly specified'
    end select
    !write(*,*) sr, yhe, hei, temp, vth_fiducial, xnue0, bl, bu, rho, kline
    !
  end function get_opalbar
  !
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !
  function get_opalbar2(iline, kline, sr, yhe, hei, temp, vth_fiducial, xnue0, bl, bu, rho)
    !
    !get opalbar for a given line transition described by iline
    !opacity law for balmer lines by assuming neutral hydrogen in ground state
    !
    !input:  sr   stellar radius in cgs
    !        yhe, hei:  helium abundance and number of free electrons per helium
    !        temp:  temperature in K
    !        vth_fiducial:  fiducial thermal velocity in cm/s
    !        xnue0:   transition frequency
    !        bl,bu:   NLTE-departure coefficients for lower and upper level
    !        rho:    density in cgs
    !        kline:   just a scaling factor to increase or decrease opacity (line strength)
    !
    ! ... arguments
    integer(i4b), intent(in) :: iline
    real(dp), intent(in) :: sr, yhe, hei, temp, vth_fiducial, xnue0, bl, bu, rho, kline
    real(dp) :: get_opalbar2
    !
    select case(iline)
       case(1)
          get_opalbar2 = opalbar_halpha2(sr, yhe, temp, vth_fiducial, xnue0, bl, bu, rho)*kline
   !   case(2)
   !      get_opalbar2 = opalbar_hbeta2(sr, yhe, temp, vth_fiducial, xnue0, bl, bu, rho)*kline
       case(10)
          get_opalbar2 = opalbar_model_kline(yhe, hei, rho, kline)*sr
       case default
          stop 'error in get_opalbar2: iline not properly specified'
    end select

    !write(*,*) sr, yhe, hei, temp, vth_fiducial, xnue0, bl, bu, rho, kline
    !
  end function get_opalbar2
  !
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !
  function get_opalbar3(iline, kline, sr, yhe, hei, temp, vth_fiducial, xnue0, bl, bu, rho)
    !
    !get opalbar for a given line transition described by iline
    !opacity law for balmer lines by assuming LTE ionization balance for hydrogen (given hei)
    !
    !input:  sr   stellar radius in cgs
    !        yhe, hei:  helium abundance and number of free electrons per helium
    !        temp:  temperature in K
    !        vth_fiducial:  fiducial thermal velocity in cm/s
    !        xnue0:   transition frequency
    !        bl,bu:   NLTE-departure coefficients for lower and upper level
    !        rho:    density in cgs
    !        kline:   just a scaling factor to increase or decrease opacity (line strength)
    !
    !
    ! ... arguments
    integer(i4b), intent(in) :: iline
    real(dp), intent(in) :: sr, yhe, hei, temp, vth_fiducial, xnue0, bl, bu, rho, kline
    real(dp) :: get_opalbar3
    !
    !
    !
    select case(iline)
       case(1)
          get_opalbar3 = opalbar_halpha3(sr, yhe, hei, temp, vth_fiducial, xnue0, bl, bu, rho)*kline
   !   case(2)
   !      get_opalbar3 = opalbar_hbeta3(sr, yhe, hei, temp, vth_fiducial, xnue0, bl, bu, rho)*kline
       case(10)
          get_opalbar3 = opalbar_model_kline(yhe, hei, rho, kline)*sr
       case default
          stop 'error in get_opalbar3: iline not properly specified'
    end select

    !write(*,*) sr, yhe, hei, temp, vth_fiducial, xnue0, bl, bu, rho, kline
    !
  end function get_opalbar3


  

end module mod_opacities
