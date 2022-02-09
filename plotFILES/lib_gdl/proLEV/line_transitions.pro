pro line_transitions, help=print_help
;
defsysv, '!lam_halpha', 6562.8399d0, 1
defsysv, '!lam_hbeta', 4861.333d0, 1  
;
;-------------------output if help is needed----------------------------
;
if(keyword_set(print_help)) then begin
   print, '-----globally stored line transitions (in A)------'
   print, ''
   print, format='(2A25)', 'h-alpha (2 <-> 3) transition', '!lam_halpha'
   print, format='(2A25)', 'h-beta  (2 <-> 4) transition', '!lam_hbeta'   
endif
;
;
;
end
