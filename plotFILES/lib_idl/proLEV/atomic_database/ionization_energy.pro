pro ionization_energy, element, eion, help=print_help
;
;-------------------output if help is needed----------------------------
;
if(keyword_set(print_help)) then begin
   print, '----------stored ionization energies--------------'
   print, ''
   print, format='(a6,a4,a6)', 'ci', '->', ' cii'
   print, format='(a6,a4,a6)', 'cii', '->', ' ciii'
   print, format='(a6,a4,a6)', 'ciii', '->', ' civ'
   print, format='(a6,a4,a6)', 'civ', '->',  'cv'
   print, format='(a6,a4,a6)', 'cv', '->',  'cvi'
endif
;
;in erg
if(element eq 'ci') then eion = 1.8040509d-11
if(element eq 'cii') then eion = 3.9065872d-11
if(element eq 'ciii') then eion = 7.6608d-11
if(element eq 'civ') then eion = 1.031896d-10
if(element eq 'cv') then eion = 6.2805324d-10
;
;
end
