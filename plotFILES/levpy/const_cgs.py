import numpy as np
#
#define some constants in cgs units and store as system variables
#from nist (or wikipedia, if explicitly noted)
#
#clight
cgs_clight = 2.99792458e10
#
#planck-constant
cgs_planck = 6.62607e-27
#
#gravitation constant
cgs_grav = 6.67408e-8
#
#proton mass
cgs_mp = 1.672621898e-24
#
#electron mass
cgs_me = 9.10938356e-28
#
#elementary charge (from wiki)
cgs_e = 4.80320427e-10
#
#boltzmann-constant
cgs_kb = 1.38064852e-16
#
#stefan-boltzmann-constant
cgs_sb = 5.6704e-5
#
#thomson-cross section
cgs_sigmae = 6.65e-25
#
#
#solar parameter
cgs_rsu = 6.96e10
cgs_msu = 1.989e33
cgs_lsu = 3.82e33
#
#distance parameters
cgs_parsec = 3.0857e18 #;from wiki
cgs_ly = 9.4607e17
cgs_au = 1.49598e13
#
#time parameters
cgs_yr = 365.25e0*24.e0*60.e0*60.e0
#
#------------------output if help is needed-----------------------------
#
#if(keyword_set(print_help)) then begin
print('----------------name of constants-----------------')
print('')
print('{vstr:30s} {vval:}'.format(vstr='vacuum speed of light', vval='cgs_clight'))
print('{vstr:30s} {vval:}'.format(vstr='planck constant', vval='cgs_planck'))
print('{vstr:30s} {vval:}'.format(vstr='gravitational constant', vval='cgs_grav'))
print('{vstr:30s} {vval:}'.format(vstr='proton mass', vval='cgs_mp'))
print('{vstr:30s} {vval:}'.format(vstr='electron mass', vval='cgs_me'))
print('{vstr:30s} {vval:}'.format(vstr='elementary charge', vval='cgs_e'))
print('{vstr:30s} {vval:}'.format(vstr='boltzmann constant', vval='cgs_kb'))
print('{vstr:30s} {vval:}'.format(vstr='stefan-boltzmann constant', vval='cgs_sb'))
print('{vstr:30s} {vval:}'.format(vstr='thomson cross section', vval='cgs_sigmae'))
print('')
print('{vstr:30s} {vval:}'.format(vstr='solar radius', vval='cgs_rsu'))
print('{vstr:30s} {vval:}'.format(vstr='solar mass', vval='cgs_msu'))
print('{vstr:30s} {vval:}'.format(vstr='solar luminosity', vval='cgs_lsu'))
print('')
print('{vstr:30s} {vval:}'.format(vstr='parsec', vval='cgs_parsec'))
print('{vstr:30s} {vval:}'.format(vstr='lightyr', vval='cgs_ly'))
print('{vstr:30s} {vval:}'.format(vstr='astr. unit', vval='cgs_au'))
print('')
print('{vstr:30s} {vval:}'.format(vstr='year', vval='cgs_yr'))
