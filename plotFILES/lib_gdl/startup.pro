HOME = './'

!PATH = '.:' + $
        HOME+'cmsvlib' + ':' + $        
        HOME+'proLEV' + ':' + $
        HOME+'proLEV/misc' + ':' + $
        HOME+'proLEV/models' + ':' + $                
        HOME+'proLEV/lte' + ':' + $
        HOME+'proLEV/grids' + ':' + $
        HOME+'proLEV/getall/3d' + ':' + $        
        HOME+'proLEV/getall/version_sc3d' + ':' + $
        HOME+'proLEV/getall/sc_cont2d' + ':' + $        
        HOME+'proLEV/getall/sc_cont3d' + ':' + $
        HOME+'proLEV/getall/forsol' + ':' + $
        HOME+'proLEV/fluxem' + ':' + $                
        HOME+'proLEV/contour_plots' + ':' + $
        HOME+'proLEV/math' + ':' + $
        HOME+'proLEV/hdf5_conversion' + ':' + $        
        HOME+'proLEV/math/interpolation' + ':' + $
        HOME+'proLEV/math/integration' + ':' + $        
        HOME+'proLEV/math/statistics' + ':' + $
        HOME+'proLEV/math/sparse' + ':' + $        
        HOME+'proLEV/colors' + ':' + $
        HOME+'proLEV/downsize' + ':' + $
        HOME+'proLEV/plots_small' + ':' + $                
        HOME+'proJO' + ':' + $
        HOME+'proGENERAL' + ':' + $
        HOME+'testsuite' + ':' + $
        HOME+'proLEV/textoidl' + ':' + $
        !PATH

;now define those programs that shall be compiled at startup
;load constants in cgs units
  const_cgs, /HELP
;load globally stored line transitions
  line_transitions, /HELP


