pro const_cgs, help=print_help
;
;define some constants in cgs units and store as system variables
;from nist (or wikipedia, if explicitly noted)
;
;clight
defsysv, '!cgs_clight', 2.99792458d10, 1
;
;planck-constant
defsysv, '!cgs_planck', 6.62607d-27, 1
;
;gravitation constant
defsysv, '!cgs_grav', 6.67408d-8, 1
;
;proton mass
defsysv, '!cgs_mp', 1.672621898d-24, 1
;
;electron mass
defsysv, '!cgs_me', 9.10938356d-28, 1
;
;elementary charge (from wiki)
defsysv, '!cgs_e', 4.80320427d-10, 1
;
;boltzmann-constant
defsysv, '!cgs_kb', 1.38064852d-16
;
;stefan-boltzmann-constant
defsysv, '!cgs_sb', 5.6704d-5
;
;thomson-cross section
defsysv, '!cgs_sigmae', 6.65d-25
;
;
;solar parameter
defsysv, '!rsu', 6.96d10
defsysv, '!msu', 1.989d33
defsysv, '!lsu', 3.82d33
;
;distance parameters
defsysv, '!parsec', 3.0857d18 ;from wiki
defsysv, '!ly', 9.4607d17
defsysv, '!au', 1.49598d13
;
;time parameters
defsysv, '!yr', 365.25d0*24.d0*60.d0*60.d0
;
;------------------output if help is needed-----------------------------
;
if(keyword_set(print_help)) then begin
   print, '----------------name of constants-----------------'
   print, ''
   print, format='(2A25)', 'vacuum speed of light', '!cgs_clight'
   print, format='(2A25)', 'planck constant', '!cgs_planck'
   print, format='(2A25)', 'gravitational constant', '!cgs_grav'
   print, format='(2A25)', 'proton mass', '!cgs_mp'
   print, format='(2A25)', 'electron mass', '!cgs_me'
   print, format='(2A25)', 'elementary charge', '!cgs_e'
   print, format='(2A25)', 'boltzmann constant', '!cgs_kb'
   print, format='(2A25)', 'stefan-boltzmann constant', '!cgs_sb'
   print, format='(2A25)', 'thomson cross section', '!cgs_sigmae'
   print, ''
   print, format='(2A25)', 'solar radius', '!rsu'
   print, format='(2A25)', 'solar mass', '!msu'
   print, format='(2A25)', 'solar luminosity', '!msu'
   print, ''
   print, format='(2A25)', 'parsec', '!parsec'
   print, format='(2A25)', 'lightyr', '!ly'
   print, format='(2A25)', 'astr. unit', '!au'
   print, ''
   print, format='(2A25)', 'year', '!yr'

endif

end
