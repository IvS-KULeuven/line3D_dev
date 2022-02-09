;******************************************************************************
;+
;*NAME: 
;
;	TIMELWR
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
;	To calculate the delta lambda term needed to derive the blaze 
;       wavelength for NEWSIPS LWR high dispersion spectra.
;
;*CALLING SEQUENCE:
;
;	TIMELWR,TIME,M,DLAM
;
;*PARAMETERS:
;
;	TIME	(REQ) (I) (0) (FD)
;		Required input vector containing the fractional year
;               of a given observation (i.e., 1981.7).
;
;	M	(REQ) (I) (01) (I)
;		Required input parameter specifying the spectral order(s)
;		for which the echelle parameters are desired.
;
;       DLAM    (OPT) (O) (01) (D)
;               Output parameter containing the delta lambda term
;               required for deriving the blaze wavelength for the latest
;               LWR NEWSIPS riplpe correction (i.e., Wblaze= (K/M) + dlam).
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
;       linterp
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
;
;*NOTES:
;
;       - From the write-up by Angelo Cassetella dated Sept. 16th, 1997.
;
;*EXAMPLES:
;
;
;*PROCEDURE:
;
;      Define coefficients and use 3rd order polynomial to calculate
;      delta lambda (dlam) as a function of date for each reference order. 
;      Then use LINTERP to linear interpolate to the input order numbers. 
;      See paper by Angelo Cassetella dated 9/16/97 for more information.       
;
;*MODIFICATION HISTORY:
;
;   09 Oct 97 RWT written 
;   17 Oct 97 RWT modify for orders 80 - 82
;   22 Oct 97 RWT use original coding for orders 80-82 but pass the 
;             order number to linterp as a real number.
;-
;******************************************************************************
pro timelwr,date,m,dlam
;
npar = n_params(0)
if (npar eq 0) then begin
    print,'timelwr,date,m,dlam'
    retall 
endif  ; print calling sequence
;
; define coefficients and reference order numbers
;
order=[115,111,107,103,99,95,91,87,83,79]
a=[35736.699D0,62433.199D0,53287.133D0,42742.709D0,28040.843D0,10463.439D0, $
   6223.1919D0,4478.2512D0,112.32014D0,156.16074D0]
b=[-35.970790D0,-62.859885D0,-53.638642D0,-43.014583D0,-28.206111D0, $
   -10.501169D0,-6.2320255D0,-4.4662128D0,-0.0566704D0,-0.0787367D0]
c=[0.0090515542D0,0.0158223297D0,0.0134980459D0,0.0108219635D0, $
   0.0070930274D0,0.0026347077D0,0.0015601476D0,0.0011134034D0,0.0D0,0.0D0]
;
; calculate delta lambda for reference orders as a function of
; observation date
;
dlref = a + (b + c*date) * date
;
; use LINTERP to linear interpolate to input orders
; (values outside range are given values of end orders)
;
;
linterp,order,dlref,float(m),dlam
;
 return
 end  ; timelwr
