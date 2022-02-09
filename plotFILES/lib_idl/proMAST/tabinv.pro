;**********************************************************************
;+
;
;*NAME:
;
;    TABINV     JUNE 15, 1981
; 
;*CLASS:
;
;*CATEGORY:
;
;*PURPOSE:  
;
;    To find the effective index of a function value in the
;    domain of a vector array representing the function.
; 
;*CALLING SEQUENCE:
;
;    TABINV,XARR,X,IEFF
; 
;*PARAMETERS:
;
;    XARR   (REQ)  (I) (1)   (BILFD)
;           The vector array to be searched
;
;    X      (REQ)  (I) (0 1) (BILFD)
;           The function value(s) whose effective
;           index is sought (scalar or vector)
; 
;    IEFF   (REQ)  (0) (0 1) (FD)
;           The effective index or indices of x in xarr
;           Double precision if either XARR or X is double precision;
;           floating point otherwise.
; 
;*EXAMPLES:
;
;    Set all flux values of a spectrum (wave vs flux) to zero
;    for wavelengths less than 1150 angstroms:
;
;         TABINV,WAVE,1150.0,i
;         for j=0,fix(i) do flux(i)=0.
;
;*SUBROUTINES CALLED:
;
;*FILES USED:
;
;*SIDE EFFECTS:
;
;*RESTRICTIONS:
;
;*NOTES:
;
;    -  The binary search technique used in tabinv requires that
;       the input array xarr is monotonic (increasing or decreasing).
;       This means that input vectors with padded zeroes will
;       cause the routine to abort.
;    -  Note that iue wavelength arrays which include the 
;       vacuum-to-air correction may not be monotonic.
;    -  If x is outside the range of xarr, ieff is set to the 
;       indice of the closest endpoint (i.e. ieff = ieff >0.0 <npt).
;    -  Due to precision differences on different machines, this routine
;       may return slightly different results.  For example, IEFF on a Vax
;       might be 25.0000019, whereas on a Sun it might be 24.9999981---FIX 
;       will truncate the results to 25 and 24 respectively.  To use directly 
;       as an array index, it may be best to round IEFF by adding 0.5 
;       beforehand.
;    -  It was discovered that the 1994 version of tabinv could return
;       truncated integer indice values if the input parameters xarr & x 
;       were vectors of byte or integer data type. No problem occurred if
;       x was a scalar value. For example, if xarr = [0,1,3,5] and x=2, 
;       then ieff = 1.5. If x = [2,4], tabinv incorrectly returned ieff = 
;       [1.0,2.0] instead  of [1.5,2.5]. The current (i.e., 1997) version 
;       converts the X parameter to floating point format in these cases 
;       and thereby prevents the problem. Note the 1994 version could have 
;       affected programs such as linterp and quadterp which call tabinv. 
;       Most programs however used tabinv to interpolate wavelengths which 
;       were rarely (if ever) specified as integers.
;
;*PROCEDURE:
;
;    The input array XARR is first tested for monotonicity
;    by evaluating the array XARR - SHIFT(XARR,1). If all
;    but one value are positive the array is considered to be
;    monotonicly increasing; if all but one are negative
;    the array is monotonicly decreasing. (Because the shifts are
;    circular, the difference between the first point and the 
;    last point is ignored.) Any other result will cause
;    tabinv to abort.
;
;    A binary search is then used to find the values XARR(I)
;    and XARR(I+1) where XARR(I) LE X LT XARR(I+1).
;    IEFF is then computed using linear interpolation 
;    between I and I+1.
; 
;         IEFF = I + (X-XARR(I)) / (XARR(I+1)-XARR(I))
; 
;    Let N = number of elements in XARR
; 
;         if X < XARR(0) then IEFF is set to 0
;         if X > XARR(N-1) then IEFF is set to N-1
; 
;*SUBROUTINES CALLED:
; 
;*MODIFICATION HISTORY:
;
;    written by d. lindler         feb. 17 1980
;    21 apr 1982 by f.h. schiffer 3rd cr#039 - monotonic functions
;         copied into [177001] rwt 1-17-84
;    7 mar 1984 by rwt - test for monotonicity added (see smr #34)
;    9 mar 1984 by rwt - monotonicity test improved (now aborts when
;         2 or more pts have same value in array xarr)
;    15 mar 1984 by rwt - more efficient coding of monotonicity test
;    5 jan 1988 rwt modify for large arrays (i.e. use longword integers),
;         use idl shift command to test monotonicity, & add procedure 
;         call listing
;    21-sep-88:   converted to sun idl, john hoegy.
;    10 May 1991 PJL corrected prolog format
;    21 Jun 1991 GRA cleaned up; tested on SUN, DEC, VAX;
;                updated prolog;
;    14 Jan 1992 RWT use byte & halfword integers instead of longword
;           integers (for arrays with < 32,000 elements) to conserve memory
;    19 Sep 1993 LLT Change 32000 limit to 16000 limit.
;     3 Apr 1994 LLT add note to prolog
;    23 Oct 1997 RWT convert byte and integer XARR arrays to floating point
;                (otherwise ieff will also be byte or integer)
;-
;***********************************************************************
 pro tabinv,xarr,x,ieff
;
 npar = n_params()
 if npar eq 0 then begin
    print,' TABINV,XARR,X,IEFF'
    retall
 endif  ; npar 
;
; get size of the arrays
;
 s  = size(xarr)
 sx = size(x)
 npoints = s(1)
 if (s(1) lt 16000) then  npoints = fix(npoints) 
 npt = npoints - 1
;
; make x into a vector if it is a scalar
; make x real if its a byte or integer
;
 if sx(0) eq 0 then begin
    xsamp = fltarr(1) + x
    nx = 1   
 end else begin
    if (sx(sx(0)+1) le 3) then xsamp = float(x) $
                          else xsamp = x
    nx = sx(1)
 endelse  ; sx(0) eq 0
;
; initialize binary search area and compute number of divisions needed
;
 ileft = intarr(nx) 
 if (nx gt 16000) then ileft = ileft * 1L
 iright = ileft
 ndivisions = fix( alog10(npoints) / alog10(2.0) + 1.0 )
;
; test for monotonicity (clever)
;
 i = xarr - shift(xarr,1)
 a = where((i gt 0), tmp_pts)
 if (tmp_pts eq npt) then begin        ; increasing array
    iright = iright + npt
 endif else begin
    a = where((i lt 0), tmp_pts)       ; test for decreasing array
    if (tmp_pts eq npt) then ileft = ileft + npt else begin
       print, 'TMP_PTS = ', tmp_pts
       print, 'NPT = ', npt
       print,'Aborting TABINV: XARR vector not monotonic' 
       retall
    endelse  ; monotonicity
 endelse  ; monotonicity 
;
; perform binary search by dividing search interval in
; half ndivision times
;
 for i = 1,ndivisions do begin
    idiv = (ileft + iright) / 2        ;split interval in half
    xval = xarr(idiv)                  ;find function values at center
    greater = (xsamp gt xval)          ;determine which side xsamp is on
    less   = (xsamp le xval)  
    ileft = ileft*less + idiv*greater  ;compute new search area
    iright = iright*greater + idiv*less
 endfor  ; i loop
;
; interpolate between interval of width = 1
;
 xleft = xarr(ileft)                   ;value on left side
 xright = xarr(iright)                 ;value on right side
 ieff = (xright-xsamp)*ileft + (xsamp-xleft)*iright + ileft*(xright eq xleft)
 ieff = ieff / (xright - xleft + (xright eq xleft))  ;interpolate
 ieff = ieff > 0.0 < npt               ;don't allow extrapolation beyond ends
 if sx(0) eq 0 then ieff = ieff(0)     ;if needed change to scalar
;
 return  
 end  ; tabinv
