;******************************************************************************
;+
;*NAME: 
;
;	NEW_SPLICE         (IUE Experimental Library)  July 22, 1988
;
;*CLASS:
;
;	IUE calibration
;
;*CATEGORY:
;
;	NEWSIPS
;
;*PURPOSE:
;
;	To calculate splice points for NEWSIPS orders based on the echelle 
;       ripple correction.
;
;*CALLING SEQUENCE:
;
;	NEW_SPLICE,H,M0,W0,M1,W1,WSPLICE
;
;*PARAMETERS:
;
;	H	(REQ) (IO) (1) (S)
;		FITS primary header as a string array.
;
;	M0	(REQ) (I) (0) (I)
;               last order number
;
;	W0	(REQ) (I) (1) (F)
;               last wavelength array (corresponding to M0)
;
;	M1	(REQ) (I) (0) (I)
;		New spectral order (must be adjacent to m0).
;
;	W1      (REQ) (I) (1) (F)
;               Required input wavelengths corresponding to order M1.
;
;	WSPLICE (REQ) (O) (0) (F)
;		Required output scalar giving the wavelength where
;               the difference in the ripple corrections are a
;               minimum.
;
;*EXAMPLES:
;
;
;*SYSTEM VARIABLES USED:
;
;	none
;
;*INTERACTIVE INPUT:
;
;	none
;
;*SUBROUTINES CALLED:
;
;       within
;       iuerip
;	PARCHECK
;       STPAR
;
;*FILES USED:
;
;	none
;
;*SIDE EFFECTS:
;
;	none
;
;*RESTRICTIONS:
;
;       - For NEWSIPS data only.
;	- Only SWP parameters currently available. (See iuerip.pro)
;    
;
;*NOTES:
;
;    - If no overlap is found in wavelengths, WSPLICE = 0.
;    - NEWSIPS reference not yet available (used code from EBRIPL.FOR)
;
;*PROCEDURE:
;
;	M and K values are calculated using IUERIP. Wavelength overlay region
;       is determined and the ripple correction is calculated for this region
;       only. WSPLICE is defined as the minimum of |cor0-cor1| where cor0
;       cor1 are the ripple corrections for the two overlapping orders.
;       If the minimum happens to be the first point in one order, WSPLICE
;       is redefined to be 10th element of wave1 as long as it is <
;       (n_elements(wave1)-1). This is to avoid target ring noise at the 
;       beginning (or end) of an extracted order.
;
;*MODIFICATION HISTORY:
;   21 Feb 97 RWT written
;   08 Jul 97 RWT remove heliocentric velocity correction from 
;                 wavelengths used in ripple correction calculation
;   16 Jan 98 RWT add delta lambda term for new LWR ripple correction
;-
;******************************************************************************
pro new_splice,h,m0,w0,m1,w1,wsplice
;
npar = n_params(0)
if (npar eq 0) then begin
    print,'NEW_SPLICE,H,M0,W0,M1,W1,WSPLICE'
    retall 
endif  ; print calling sequence
parcheck,npar,6,'NEW_SPLICE'
wsplice = 0.0
;
;  call IUERIP to get K & ALPHA parameters
;  (and dlam for new LWR correction)
;
order = [m0,m1]

iuerip,h,order,k,alpha,dlam
; 
;  uncorrect wavelengths for heliocentric velocity correction 
;  before ripple correction
;
stpar,h,'aperture',ap,err
ap = strlowcase(strtrim(ap,2))
ap = strmid(ap,0,1)
if (ap eq 'b') then begin
    ap = 'l'
    print,' aperture = "BOTH", setting ap = "l" by default'
endif
stpar,h,ap+'radvelo',vel,err
vcorr = 1.0D0 + (vel/2.99792501D+5)
;
;  save only portions which overlap
;
   wave0 = w0
   wave1 = w1
   within,wave0,wave1,res
   ind = where(res eq 0,nres1)
   if (nres1 eq 0) then return      ; no overlap then return
   wave1 = wave1(ind)
   within,wave1,wave0,res
   ind = where(res eq 0,nres0)
   wave0 = wave0(ind)      
;
;  calculate ripple correction for both orders (in overlap region)
;
pi = !dpi
wave0 = wave0 / vcorr
wave1 = wave1 / vcorr
for i=0,1 do begin
    wblaze = (k(i)/order(i)) + dlam(i)
    if (i eq 0) then $
      x = pi * order(i) * alpha(i) * ((wave0 - wblaze)/wave0) else $
      x = pi * order(i) * alpha(i) * ((wave1 - wblaze)/wave1) 
    absx = abs(x)
    x = absx + .01 * (absx lt .01)
    x = x < 2.61
    sincx = sin(x)/x
    cor = sincx * sincx
    if (i eq 0) then cor0 = cor else cor1 = cor
endfor
;
;  calculate minimum of |cor0-cor1|
;
 diff = abs(cor0 - cor1)
 ind = where(diff eq min(diff))
;
; discard first 10 points of order (if that many pts are in overlap)
;
 if (ind(0) eq 0) then ind(0) = 10 < (n_elements(wave1)-1)
 wsplice = wave1(ind(0)) * vcorr    ; add velocity correction for splice point
;
 return
 end  ; new_splice
