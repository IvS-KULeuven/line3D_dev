;******************************************************************************
;+
;*NAME: 
;
;	IUERIP         (IUE Experimental Library)  July 22, 1988
;
;*CLASS:
;
;	IUE calibration
;
;*CATEGORY:
;
;	IUESIPS,NEWSIPS
;
;*PURPOSE:
;
;	To calculate the echelle ripple parameters K and ALPHA for order M
;	for a given IUE spectrum.
;
;*CALLING SEQUENCE:
;
;	IUERIP,H,M,K,ALPHA,dlam,swpt=swpt
;
;*PARAMETERS:
;
;	H	(REQ) (IO) (1) (I)
;		Required input vector containing the scale factor information,
;		camera number, and ITF number for IUESIPS data, or the FITS
;               header for NEWSIPS data.  New K and Alpha values are
;		stored in H(97) and H(99) as done in QIUEHI3, or as FITS 
;               keywords.
;
;	M	(REQ) (I) (01) (I)
;		Required input parameter specifying the spectral order(s)
;		for which the echelle parameters are desired.
;
;	K	(REQ) (O) (01) (F)
;		Required output parameter giving the echelle blaze
;		parameter for order M, and iue camera number given by H(3).
;
;	ALPHA	(REQ) (O) (01) (F)
;		Required output parameter giving the echelle blaze
;		scale factor for camera given in H(3).
;
;       DLAM    (OPT) (O) (01) (D)
;               Output parameter containing the delta lambda term
;               required for deriving the blaze wavelength for the latest
;               LWR NEWSIPS riplpe correction (i.e., Wblaze= (K/M) + dlam).
;               Set to 0.0D0 for all other cameras.
;
;	SWPT    (KEY) (O) (0) (BILFD)
;               If specified, it will be used as the THDA value in the 
;               SWP ripple correction. (Normally the THDA at end of exposure
;               is ectracted from the header.) Keyword has no effect on other
;               cameras.
;
;*EXAMPLES:
;
;	To calculate the K value for order 100 for the SWP:
;
;		iuerip,h,100,k,alpha
;
;       To calculate K and alpha for orders 92-94:
;
;               iuerip,h,[92,93,94],k,alpha
;
;       To calculate K and alpha for SWP orders 70 using a THDA
;       value of 8.5:
;
;               iuerip,h,70,k,alpha,swpt=8.5
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
;	PARCHECK
;       STPAR
;       DATECONV
;       TIMELWR
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
;       - For IUESIPS data,
;         For wavelengths >2000 A, the ripple data are given for AIR wavelengths.
;         For wavelengths <2000 A, the ripple data are given for VACUUM
;         wavelengths.
;         Since the LWR and LWP cameras have very low sensitivity <2000 A, 
;         ripple data for these cameras are given only for AIR wavelengths. 
;       - For NEWSIPS data,
;         All wavelengths values are assumed to be in vacuum and in the 
;         satellite reference frame (i.e., no heliocentric velocity correction)
;
;*NOTES:
;
;       - References for the IUESIPS ripple parameterizations:
;		LWP  (old ITF/ air wavelengths:) Ake (1985)
;		     (ITF 2 / air wavelengths:) Ake (1985) 
;
;		LWR  (old ITF/ air wavelengths:) Ake (1982)
;
;		SWP  (old ITF/ vacuum wavelengths:) Grady and Garhart 1989,
;						    NASA IUE Newsletter 37.
;       - NEWSIPS reference not yet available (IUERIPR is based on the
;       code from EBRIPL.FOR, and the latest LWR ripple correction came 
;       from the write-up from Angelo Cassetella dated Sept. 16th, 1997.
;
;	- Note that LWP spectra processed 1989 and earlier are subject to
;	spectral mis-registration and extraction errors which cause loss of
;	light from the gross spectrum, and the introduction of spectral light
;	into the inter-order background at the short wavelength side of each
;	order.  The effect is most pronounced for those portions of a spectrum
;	where the gross spectrum is at most a few time the interorder
;	background level, and may result in systematic photometric errors of
;	10-20%.  The ripple correction in use here will not correct for this
;	effect, with the result that the net ripple corrected spectrum appears
;	"saw-toothed".  See Grady, et al. IUE NASA Newsletter 38 for details. 
;       
;       - The SWP NEWSIPS thda dependence uss the temperature at the end
;        of exposure not the start of exposure. The prolog of the NEWSIPS 
;        module EBRIPL incorrectly states that the THDA at start of exposure 
;        was used. Actually, according to Angelo Cassatella, the analysis 
;        was intended to be based on the THDA at time of read, but the
;        THDAs at time of read and end of exposure for the images he used 
;        were basically the same.
;
;       - It was discovered that a 2 degree difference in the THDA value can
;        change the SWP ripple correction by more than 7 percent! Increasing 
;        the THDA value shifts the blaze wavelength (k/m) to a higher 
;        wavelength. The effect on the ripple correction is largest at 
;        the ends of the orders and increases with order number.
;
;*EXAMPLES:
;
;       To calculate K and alpha for orders 92-94:
;               iuerip,h,[92,93,94],k,alpha
;
;*PROCEDURE:
;
;       IUERIP uses the H vector to distinguish IUESIPS and NEWSIPS data. 
;	M and K values are calculated using the polynominal parameterizations
;	as a function of camera, ITF, and order number given in the references
;	listed above. NEWSIPS SWP camera K values have a THDA correction, LWP 
;       and LWR NEWSIPS K values use a time-dependent correction. The LWR 
;       NEWSIPS alpha term is different for orders above and below order 100.
;
;*MODIFICATION HISTORY:
;
;    7-22-88 CAG GSFC RDAF  initial program, based on code in IUEHI, but 
;            generalized to accomodate multiple ITFs.
;    9-21-88 RWT store K and Alpha values in H vector and add defaults for 
;            new LWP
;    3-08-89 CAG GSFC RDAF  replace Grady and Fireman ripple for LWP with
;            earlier Ake (1985) ripple, since that function performs better
;            (for discussion, see Grady et al. 1989).
;    9 Jul 91 LLT convert to lowercase, clean up, update prolog, tested on VAX
;   23 Jul 91 PJL tested on SUN; updated prolog
;   08 Nov 93 PJL add H parameter check and IUESIPS restriction
;   13 May 96 RWT add newsips support (only SWP for now)
;   21 Feb 97 RWT updated SWP NEWSIPS parameters based on changes to EBRIPL.FOR
;             dated 9/17/96
;   24 Feb 97 RWT allow M (and therefore K & alpha) to be vectors for NEWSIPS
;             data.
;   30 Jun 97 RWT add LWP NEWSIPS parameters (based on EBRIPL.FOR)
;   12 Aug 97 RWT add LWR NEWSIPS parameters (based on EBRIPL.FOR)   
;   19 Sep 97 RWT make M and thda double precision
;   23 Sep 97 RWT extract THDA at END of exposure (not at START of exposure)
;   02 Oct 97 RWT add SWPT keyword 
;   09 Oct 97 RWT add new LWR NEWSIPS parameters & call TIMELWR to derive new
;                 delta lambda term.
;-
;******************************************************************************
pro iuerip,h,order,k,alpha,dlam,swpt=swpt
;
npar = n_params(0)
if (npar eq 0) then begin
    print,'IUERIP,H,ORDER,K,ALPHA,dlam,swpt=swpt'
    retall 
endif  ; print calling sequence
parcheck,npar,[4,5],'IUERIP'
;
;  check header characteristics
;  to distinguish iuesips from newsips data
;
temp = size(h)
if (temp(temp(0)+1) eq 7) then newsips = 1 else newsips = 0
;
if (newsips) then begin        ; newsips
;
; initialize parameters & see if order is scalar or vector
;
  temp = size(order)
  if (temp(0) eq 0) then begin
     k = dblarr(1)
     m = intarr(1)
     m(0) = double(order)
  endif else begin
     k = double(order) * 0.0D+00
     m = double(order)
  endelse
  alpha = k
;
  stpar,h,'camera',cam,err
  cam = strlowcase(strtrim(cam,2))
  stpar,h,'aperture',ap,err
  ap = strlowcase(strtrim(ap,2))
  ap = strmid(ap,0,1)
  if (ap eq 'b') then begin
    ap = 'l'
    print,' aperture = "BOTH", setting ap = "l" by default'
  endif
;
; read THDA at end of exposure (i.e., not start of exposure)
; (or use keyword value)
;  stpar,h,ap+'thdastr',thda,err
  if (keyword_set(swpt)) then thda = double(swpt) else begin
     stpar,h,ap+'thdaend',thda,err
     thda = double(thda)
  endelse
;
; LWP K term varies with observation date,
; as does the LWR dlam term
; use dateconv to convert DATEOBS keyword to a decimal year
;
stpar,h,ap+'dateobs',date,err         ; observation date
dateconv,date,'v',vd                    ; convert to vector
date = vd(0) + vd(1)/365.0D0            ; convert to decimal date
;
;  call timelwr to derive dlam 
;  
  if (cam eq 'lwr') then timelwr,date,m,dlam $
                    else dlam = m * 0.0D0
;
; assign ripple parameters as in NEWSIPS routine EBRIPL
;
  for i=0,n_elements(m)-1 do begin
    case cam of
      'lwp': begin
;
; values as of 5/7/97

             alpha(i) =  0.406835D0 + 0.01077191D0 * m(i) -  $
                         5.945406D-5 * m(i) * m(i)
             a1    =  230868.1770D0
             a2    =  56.433405D0
             t1    =  0.0D0
             t2    =  -0.0263910D0
             k(i) = a1 + (a2 + t1 + t2*date) * m(i) 
             end
;
      'lwr': begin
;
; values as of 8/11/97

             if (m(i) le 100) then alpha(i) = 3.757863D0 -     $
                0.0640201D0 * m(i) + 3.5664390D-4 * m(i) * m(i) $
             else  alpha(i) = 1.360633D0 - 4.252626D-3 * m(i)
;             a1 = 230538.5180D0
;             a2 = 90.768579D0
;             t1 = 0.0D0
;             t2 = -0.0425003D0
;             k(i) = a1 + (a2 + t1 + t2*date) * m(i)        
;             print,'old K =',k(i)
; new K value (as of 10/9/97)

             a1 = 0.281749635D+06
             a2 = -0.223565585D+04
             a3 = 0.365319482D+02
             a4 = -0.262477775D+00
             a5 = 0.701464055D-03
             k(i) = a1 + ( a2 + (a3 + a4*m(i) + a5*m(i)*m(i)) * m(i) ) * m(i)
;             print,'new K =',k(i)
             end
;
      'swp': begin
;
; old values
;
;               alpha = 0.9779439D+00 - 0.0013619907D+00 * m
;               a1 = 137574.506D+00
;               a2 = 1.29145D+00
;               t1 = 0.0D+00
;               t2 = 0.0584406D+00
;
; new values (as of 9/17/96)
;
             alpha(i) = 0.926208D+00 - 7.890132D-04 * m(i)
             a1 = 137508.316D+00
             a2 = 2.111841D+00
             t1 = 0.0D+00
             t2 = 0.0321729D+00
             k(i) = a1 + (a2 + t1 + t2*thda) * m(i)        
             end
      else: begin
              print,' Invalid camera name'
              print,' Returning from IUERIP'
              retall
            end
    endcase ; cam
  endfor

endif else begin      ; iuesips
  m = order
  cam = h(3)
  itf = h(580)
  k0_val = [231150.0,231150.0,137725.0,137725.0]
  k0 = k0_val(cam - 1)
  case cam of
     1: begin
           if (itf lt 2) then begin   ; Ake (1985) for ITF 1
              k = 230648. + 5.391 * m
              alpha= 0.896
           endif else begin   ; Ake (1985) for ITF 2/air wavelengths
              k = 230648. + 5.391 * m 
              alpha = 0.896
           endelse  ; itf
        end  ; case of LWP camera
     2: begin  ; Ake 1982 ripple correction parameters
           k = 230036. + 15.3456 * m - 0.050638 * m * m
           alpha = 0.896
        end  ; case of LWR camera
     3: begin  ; Grady and Garhart 1989
           k = 137730.0 - 3.0644019 * m + 0.0335206184 * m * m
           alpha = 0.856
        end  ; case of SWP camera
     else: begin   
              print,'Invalid camera id stored in H record.'
              print,'Action: Returning'
              retall 
           end  ; invalid camera entered
  endcase  ; cam
  h([97,98,99]) = [fix(k-k0),0,fix(alpha*1000)]  ;  store values in H vector
endelse
;
 return
 end  ; iuerip
