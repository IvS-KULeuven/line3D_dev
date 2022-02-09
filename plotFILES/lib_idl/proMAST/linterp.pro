;*************************************************************************
;+
;
;*NAME:   
;
;    LINTERP  (General IDL Library 01)      MARCH 15, 1981
;
;*CLASS:
;
;    Interpolation
;
;*CATEGORY:
;    
;*PURPOSE: 
;
;    To linearly interpolate tabulated data from one data grid to another.
;
;*CALLING SEQUENCE:
;
;    LINTERP,XTAB,YTAB,XINT,YINT
;
;*PARAMETERS: 
;
;    XTAB     (REQ) (I) (1) (I L F D)
;             Required input vector containing the current independent
;             variable grid.
;    YTAB     (REQ) (I) (1) (I L F D)
;             Required input vector containing the current dependent
;             variable values at the XTAB grid points.
;    XINT     (REQ) (I) (0 1) (I L F D)
;             Required input scalar or vector containing the new
;             independent variable grid points for which interpolated
;             value(s) of the dependent variable are sought.
;
;    YINT     (REQ) (O) (0 1) (F D) 
;             Required output scalar or vector with the interpolated
;             value(s) of the dependent variable at the XINT grid points.
; 
;*EXAMPLES:
;
;    To linearly interpolate from an IUE spectrum wavelength grid to
;    another grid defined as:
;    WGRID=[1540., 1541., 1542., 1543., 1544, 1545.]
;
;    LINTERP,WAVE,FLUX,WGRID,FGRID
;
;*SYSTEM VARIABLES CALLED:
;
;*INTERACTIVE INPUT:
; 
;*SUBROUTINES CALLED:
;
;    TABINV
;    PARCHECK
;
;*FILES USED: 
;
;*SIDE EFFECTS:
;
;*RESTRICTIONS:
;    - XTAB must be monotonically increasing or decreasing (a requirement
;      for subroutine TABINV). 
;
;*NOTES:
;    - Points outside the range of XTAB are set equal to the first or last
;      value of YTAB. This is different than the program GRIDTERP which 
;      sets points outside the range of the input vectors to zero.
;    - Spacing between data points in input vectors may vary.
;    - Input vectors are preserved (unlike GRIDTERP).
;
;*PROCEDURE: 
;
;    Uses TABINV to calculate the effective index of the values
;    in XINT in the table XTAB.  The resulting index is used
;    to calculate the interpolated values YINT from the values
;    in YTAB. 
;
;*MODIFICATION HISTORY:
;
;    Mar 15 1981  D. Lindler   initial program
;    Dec 22 1981  FHS3    GSFC to allow scalar XINT and to use TABINV in the
;                              interpolation
;    Oct 23 1985  JKF     GSFC DIDL compatible...replaced function REORDER
;                              with vector subscripting and indirect 
;                              compilations
;    Jun  5 1987  RWT     GSFC add PARCHECK
;    Mar 14 1988  CAG     GSFC add VAX RDAF-style prolog
;                              and printing of calling sequence when
;                              executed without parameters.
;    Jun 21 1991  GRA     CASA cleaned up; tested on SUN, DEC, VAX;
;                              lower case; updated prolog
;    Jun 30, 1997 RWT     added documentation
;-
;****************************************************************************
 pro linterp,xtab,ytab,xint,yint
;
 npar = n_params()
 if npar eq 0 then begin
    print,'LINTERP,XTAB,YTAB,XINT,YINT'
    retall
 endif  ; npar
 parcheck,npar,4,'LINTERP'
 pcheck,xtab,1,010,0111
 pcheck,ytab,2,010,0111
 pcheck,xint,3,110,0111
;
 x = xint
 s = size(x)
 if s(0) eq 0 then  x = fltarr(1) + x 
;
; determine index of data-points from which interpolation is made
;
 tabinv,xtab,x,r
 index = fix(r)
 r = r - index
;
; perform linear interpolation
;
 yint = ytab(index + 1)*r + ytab(index)*(1 - r)
;
; return scalar if input was scalar
;
 if s(0) eq 0 then yint = yint(0)
;
 return
 end  ; linterp
