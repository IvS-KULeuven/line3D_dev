PRO psfullps,type=type,name=name,xsize=xsize,ysize=ysize,col=col,reverse=reverse,brewer=brewer

if keyword_set(col) then begin 
   if (col eq 100) then cgloadct,0 else cgloadct,col
   IF KEYWORD_SET(reverse) THEN cgloadct,col,/reverse
endif else begin 
   colors
endelse
IF keyword_set(brewer) THEN cgloadct,col,/brewer
IF keyword_set(brewer) and keyword_set(reverse) THEN cgloadct,col,/brewer,/reverse

if not keyword_set(type) then type='l'
if not keyword_set(name) then name='idl.eps'


if not keyword_set(xsize) and type eq 'l' then xsize=26 
if not keyword_set(ysize) and type eq 'l' then ysize=20
if not keyword_set(xsize) and type eq 'p' then xsize=20 
if not keyword_set(ysize) and type eq 'p' then ysize=26

print,'Starting plot to: '+name

if type eq 'l' then begin
set_plot,'ps'
device,filename=name,xsize=xsize,ysize=ysize,xoffset=0.0,/landscape,/color,$
bits_per_pixel = 8,/ISOLATIN1,/encapsul
endif

if type eq 'p' then begin
set_plot,'ps'
device,filename=name,xsize=xsize,ysize=ysize,xoffset=0.0,/color,$
bits_per_pixel = 8,/ISOLATIN1,/encapsul
;device,filename=name,xsize=xsize,ysize=ysize,xoffset=0.0,zoffset=0.0,/color,$
;bits_per_pixel = 8,/ISOLATIN1
endif

END
