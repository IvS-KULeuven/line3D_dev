function sign, x
;
nd = n_elements(x)
;
if(nd eq 1) then begin
   if(x ge 0.) then begin
      signs=1
   endif else begin
      signs=-1
   endelse
endif else begin
   signs=indgen(nd)
   for i=0, nd-1 do begin
      if(x(i) ge 0.) then begin
         signs(i)=1
      endif else begin
         signs(i)=-1
      endelse
   endfor
endelse
;
return, signs
;
end
