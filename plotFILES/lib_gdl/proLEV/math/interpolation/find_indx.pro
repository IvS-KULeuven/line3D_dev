pro find_indx, xin, xarr, nd, iim1, ii, help=print_help
;
;-----------------------------------------------------------------------
;
;  finds index for interpolation of two points
;
;  input:  
;          
;
;

;+
; NAME:
;       find_indx
; PURPOSE:
;       finds index for interpolation of two points
;
;       for interpolation: xarr(iim1)----xin-----xarr(ii)
;       for extrapolation: xin----xarr(iim1)-----xarr(ii)
;                          xarr(iim)----xarr(ii)-----xin
;
; CALLING SEQUENCE:
;       find_indx, xin, xarr, nd, iim1, ii
;
; INPUTS:
;       coordinate of point:             xin
;       dimension of grid:               nd
;       grid:                            xarr
;
; OUTPUTS: 
;       indices:   iim1, ii
;
; KEYWORDS:
;
;-
;
;-----------------------------------------------------------------------
;
if(keyword_set(print_help)) then begin
   doc_library, 'find_indx'
   return
endif
;
;-----------------------------------------------------------------------
;
if(xin ge xarr(nd-1)) then begin
   iim1=nd-2
   ii=nd-1
   return
endif
;
iim1=0
ii=1
for i=1, nd-1 do begin
   if(xarr(i) ge xin) then begin
      ii=i
      iim1=i-1
      break
   endif
endfor
;
end
