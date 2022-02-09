;+****************************************************************************
;
;*NAME:
;
;     within
;
;*PURPOSE:
;
;     to determine whether the elements of a test array lie within the range
;     of another array.
;
;*CALLING SEQUENCE:
;
;     WITHIN,LIMARRAY,TESTARRAY,RESULT
;
;*PARAMETERS:
;
;     limarray    (req) (i) (1) (bidlf)
;                 Array whose minimum and maximum values will serve as the
;                 lower and upper limits for the test array.
;
;     testarray   (req) (i) (1) (bidlf)
;                 Array to be tested to see whether each element falls within
;                 the limits defined by the limarray.
;
;     result      (req) (o) (1) (i)
;                 Resulting array with the same dimensions as testarray.
;                 Elements less than the minimum limit are set equal to -1.
;                 Elements greater than the maximum limit are set equal to 1.
;                 Elements within the limits are set equal to zero.
;
;*SYSTEM VARIABLES USED:
;
;     none
;
;*FILES USED:
;
;     none
;
;*SUBROUTINES USED:
;
;     parcheck
;
;*PROCEDURE:
;
;     The limits are set equal to max(limarray) and min(limarray).  The
;     result vector is created which will have the same structure as the
;     given testarray, elements being set equal to zero.  For any elements
;     in the testarray less than the minimum limit, the corresponding elements
;     in the result array are set equal to -1.  For any elements in the 
;     testarray greater than the maximum limit, the corresponding elements in
;     the result array are set equal to 1.  
;
;*MODIFICATION HISTORY:
;
;     Written 8 February 1994 by LLT.
;
;-*****************************************************************************
pro within,limarray,testarray,result

npar=n_params(0)
if npar eq 0 then begin
   print,'WITHIN,LIMARRAY,TESTARRAY,RESULT'
   retall
endif ;npar eq 0

parcheck,npar,[3],'WITHIN'

result=fix(testarray)*0
upper=max(limarray)
lower=min(limarray)
temp=where(testarray lt lower,count)
if count gt 0 then result(temp)=-1
temp=where(testarray gt upper,count)
if count gt 0 then result(temp)=1

return
end ;within
