FUNCTION DOT_PRODUCT, V1, V2, help=print_help
;
;+
; NAME:
;	DOT_PRODUCT
;
; PURPOSE:
;	This function calculates the dot-product for two given vectors.
;
; CALLING SEQUENCE:
;	Result = DOT_PRODUCT(v1, v2)
;
; INPUTS:
;	v1:	Vector 1 of given dimension
;	v2:	Vector 2 (with same dimension if v1)
;	
; KEYWORD PARAMETERS:
;	help:	Set this keyword (flag) tho show the documentation of this function
;
; EXAMPLE:
;	Result = DOT_PRODUCT([1.,1.,1.],[2.,3.,4])
;-
;
;--------------------if help is needed----------------------------------
;
if(keyword_set(print_help)) then begin
   doc_library, 'dot_product'
   return, 0
endif
;
;----------------------------check dimensions---------------------------
;
nd=n_elements(v1)
nd1=n_elements(v2)
;
if(nd ne nd1) then begin
   print, 'error dot_product: input vectors v1, v2 need to have same number of elements'
   return, 0
endif
;
;-----------------------------------------------------------------------
;
vtemp = v1*v2
;
sum=0.d0
for i=0, nd-1 do begin
   sum=sum+vtemp(i)
endfor
;
RETURN, sum
;
END
