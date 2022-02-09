pro loadct_lev, ci_black, ci_blue, ci_red, ci_cyan, ci_green, ci_magenta, ci_yellow, ci_white, ci_llgrey=ci_llgrey, ci_dgreen=ci_dgreen, help=print_help
;+
; NAME:
;       loadct_lev
;
; PURPOSE:
;       This procedure loads some basic colors for the current graphics device
;       and returns the corresponding color-indices
;
; CALLING SEQUENCE:
;       loadct_lev, ci_black, ci_blue, ci_red, ci_cyan, ci_green, ci_magenta, ci_yellow, ci_white
;
; INPUTS:
;
; KEYWORD PARAMETERS:
;       help:   Set this keyword (flag) tho show the documentation of this function
;
; EXAMPLE:
;       loadct_lev, ci_black, ci_blue, ci_red, ci_cyan, ci_green, ci_magenta, ci_yellow, ci_white, help=print_help
;-
;
;--------------------if help is needed----------------------------------
;
if(keyword_set(print_help)) then begin
   doc_library, 'loadct_lev'
   return
endif
;
;-----------------------------------------------------------------------

red=[0,0,0,255,255,0,255,255,0,254]
green=[0,0,255,0,255,255,0,255,150,255]
blue=[0,255,0,0,0,255,255,255,0,255]
color_name=['black','blue','green','red','yellow','cyan','magenta','white','dgreen','llgrey']
ci_black=where(color_name eq 'black')
ci_blue=where(color_name eq 'blue')
ci_red=where(color_name eq 'red')
ci_cyan=where(color_name eq 'cyan')
ci_green=where(color_name eq 'green')
ci_white=where(color_name eq 'white')
ci_magenta=where(color_name eq 'magenta')
ci_yellow=where(color_name eq 'yellow')
ci_dgreen=where(color_name eq 'dgreen')
ci_llgrey=where(color_name eq 'llgrey')
;
tvlct, red, green, blue

end
