!
!-----------------------------------------------------------------------
!------------------------TYPE-DEFINITION--------------------------------
!-----------------------------------------------------------------------
!
module prog_type  

integer, parameter :: i1b = selected_int_kind(1)
integer, parameter :: i4b = selected_int_kind(9)
integer, parameter :: i8b = selected_int_kind(10)
integer, parameter :: sp = kind(1.0)                        
integer, parameter :: dp = kind(1.d0)
integer, parameter :: qp = selected_real_kind(33,4931)
!integer, parameter :: dp = qp

!sp: single precision ( 7 significant digits, range [10^-38 , 10^38 ])
!dp: double precision (15 significant digits, range [10^-308, 10^308])
!qp: high   precision (33 significant digits, range [10^-4931, 10^4932]
!i1b: 1-byte integer (range [-128,127])
!i4b: 4-byte integer (range [10^-8 , 10^8])
!i8b: 8-byte integer (range [10^-18 , 10^18]))

! PL migrating the parameter definitions to accomodate the NLTE and Line modules of MForce
INTEGER, PARAMETER :: I2B = SELECTED_INT_KIND(2)

INTEGER, PARAMETER :: stdlen = 80
INTEGER, PARAMETER :: ibmlen = 132

INTEGER, PARAMETER :: stdin  = 5
INTEGER, PARAMETER :: stdout = 6
INTEGER, PARAMETER :: stderr = 0

end module prog_type
!
!-----------------------------------------------------------------------
!---------------------FUNDAMENTAL CONSTANTS-----------------------------
!-----------------------------------------------------------------------
!
module fund_const
!
use prog_type
!
implicit none
!
real(dp), parameter :: pi=4.d0*atan(1.d0) !pi = 3.141592653589793
real(dp), parameter :: spi=sqrt(pi)
real(dp), parameter :: almost_one=0.9999999999d0
real(dp), parameter :: small_number=1.d-14
real(dp), parameter :: zero=0.d0
real(dp), parameter :: half=0.5d0
real(dp), parameter :: one=1.d0
real(dp), parameter :: two=2.d0
real(dp), parameter :: three=3.d0
real(dp), parameter :: four=4.d0
real(dp), parameter :: five=5.d0
real(dp), parameter :: six=6.d0
real(dp), parameter :: seven=7.d0
real(dp), parameter :: eight=8.d0
real(dp), parameter :: ten=10.d0
real(dp), parameter :: twelve=12.d0
real(dp), parameter :: fourteen=14.d0
real(dp), parameter :: fifteen=15.d0
real(dp), parameter :: twentyfour=24.d0
real(dp), parameter :: twentyfive=25.d0
real(dp), parameter :: twentyseven=27.d0
real(dp), parameter :: fortyfive=45.d0
real(dp), parameter :: fortyeight=48.d0
real(dp), parameter :: sixtyfour=64.d0
real(dp), parameter :: euln=exp(1.d0)
real(dp), parameter :: cgs_clight=2.99792458d10
real(dp), parameter :: cgs_planck=6.6262d-27
real(dp), parameter :: cgs_kb=1.38062d-16
real(dp), parameter :: cgs_sb=5.670367d-5
real(dp), parameter :: cgs_mp=1.67265d-24
real(dp), parameter :: cgs_me=9.10938d-28
real(dp), parameter :: cgs_grav=6.6740800d-8
real(dp), parameter :: sigmae = 6.65d-25
real(dp), parameter :: cgs_au = 1.49598d13
real(dp), parameter :: rsu = 6.96d10, xmsu = 1.989d33, xlsu = 3.82d33
real(dp), parameter :: cgs_e=4.80320427d-10
real(dp), parameter :: yr=365.25d0*24.d0*3600.d0
!
!
!cgs_clight: speed of light
!cgs_planck: planck-constant
!cgs_kb: boltzmann-constant
!cgs_sb: stefan-boltzmann constant
!sigmae: cross-section thomson scattering (see radprocesses p. 66)
!cgs_mp: mass of hydrogen atom / proton
!cgs_me: mass of electron
!cgs_e: elementary charge
!rsu: radius sun
!xmsu: mass sun
!xlsu: luminosity sun

end module fund_const
