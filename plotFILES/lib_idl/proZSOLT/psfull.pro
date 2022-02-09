PRO psfull,type=type,name=name,xsize=xsize,ysize=ysize,col=col

if keyword_set(col) then begin 
   loadct,col
endif else begin 
   colors
endelse

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
/encapsulated,bits_per_pixel = 8,/ISOLATIN1
endif

if type eq 'p' then begin
set_plot,'ps'
device,filename=name,xsize=xsize,ysize=ysize,xoffset=0.0,zoffset=0.0,/color,$
/encapsulated,bits_per_pixel = 8,/ISOLATIN1
endif

END
