function interpol_ypl, x_im1, x_i, f_im1, f_i, x_p, help=print_help
;
;+
; name:
;	interpol_ypl
;
; purpose:
;	this function interpolates linearly in log-log between two points.
;       if log-log is not allowed, interpolation is linear
;       
;
; calling sequence:
;	interpol_ypl, x_im1, x_i, f_im1, f_i, x_p
;
; inputs:
;	x_im1, x_im2:	coordinates of given points
;	f_im1, f_im2:	function values at given points
;	x_p:     coordinate of unknown point
;	
; keyword parameters:
;	help:	set this keyword (flag) tho show the documentation of this procedure
;
; outputs:
;	f_p:     function value at unknown point
;
; example:
;	interpol_ypl, 1., 2., 1., 4., 1.5
;-
;
;-------------------------OUTPUT IF HELP NEEDED-------------------------
;
if(keyword_set(print_help)) then begin
   doc_library, 'interpol_ypl'
   return, 0
endif
;
;-----------------------------------------------------------------------
;
llog=1
;logarithmic interpolation not allowed when
if(x_im1 eq 0.) then begin
   llog=0
endif else begin
   if(x_p/x_im1 le 0.) then llog=0
   if(x_i/x_im1 le 0.) then llog=0
endelse
;
if(f_im1 eq 0.) then begin
   llog=0
endif else begin
   if(f_i/f_im1 le 0.) then llogf=0
endelse
;
if(llog eq 1) then begin
   grad = alog10(f_i/f_im1)/alog10(x_i/x_im1)
   return, f_im1*(x_p/x_im1)^grad
endif else begin
   return, f_im1 + (f_i-f_im1)/(x_i-x_im1) * (x_p-x_im1)
endelse
;
end
