;******************************************************************************
;+
;*NAME:
;
;    FLAGNAN
;
;*CLASS:
;
;*CATEGORY:
;
;*PURPOSE:
;
;    To flag and/or replace IEEE NaN values in an IDL array.
;
;*CALLING SEQUENCE:
;
;    FLAGNAN,ARRAY,IND,COUNT,nvalue
;
;*PARAMETERS:
;
;    ARRAY (REQ) (IO) (012) (RD)
;       IDL array. 
;
;    IND (REQ) (O) (1) (I) 
;       indices of found NaN values
;
;    COUNT (REQ) (O) (0) (IL)
;       number of found NaN values. (0 if none found.)
;
;    nvalue (OPT) (I) (0) (BILF)
;       optional parameter specifying replacement value.
;       If not specified, user is prompted for value.
;
;*EXAMPLES:
;
;*SYSTEM VARIABLES USED:
;
;*INTERACTIVE INPUT:
;
;*SUBROUTINES CALLED:
;
;   PARCHECK
;
;*FILES USED:
;
;*SIDE EFFECTS:
;
;	If no NaN values are found, or if ARRAY is not of type float, double
;	precision, or complex, then -1 is returned, and COUNT is set to 0.
;
;*RESTRICTIONS:
;
;	ARRAY must be of type float, double-precision, or complex.
;
;*PROCEDURE:
;
;	The bit patterns of the numbers being tested are compared against the
;	IEEE NaN standard.
;
;*INF_1:
;
;*NOTES:
;
;	tested with IDL Version 2.2.0 (sunos sparc)   14 Nov 91
;       tested with IDL Version 2.2.0 (ultrix mipsel) 14 Nov 91
;       tested with IDL Version 2.2.0 (ultrix vax)    14 Nov 91
;       tested with IDL Version 2.3.2 (vax vms)       19 Oct 92
;
;*MODIFICATION HISTORY:
;
;	written by Bill Thompson 2/92
;        8-14-92 rwt add nvalue option, change from a function to a 
;               procedure & print out results
;       10-19-92 rwt make ind and count required output parameters
;-
;********************************************************************
 pro flagnan,array,ind,count,nvalue
;
 npar = n_params()
 if (npar eq 0) then begin
    print,' FLAGNAN,ARRAY,IND,COUNT,nvalue'
    retall
 endif  ; npar
 parcheck,npar,[3,4],'FLAGNAN'
;
; call wherenan
;
;
	ON_ERROR,2
;
;  Parse the input array based on the datatype.
;
	SZ = SIZE(ARRAY)
	CASE SZ(SZ(0)+1) OF
;
;  Single precision floating point.
;
		4:  BEGIN
			LARRAY = LONG(ARRAY,0,N_ELEMENTS(ARRAY))
                        BYTEORDER,LARRAY,/NTOHL
			E0 = '7F800000'X
			E = LARRAY AND E0
			F = LARRAY AND '7FFFFF'X
			RESULT = WHERE((E EQ E0) AND (F NE 0), COUNT)
			END
;
;  Double precision floating point.
;
		5:  BEGIN
			LARRAY = LONG(ARRAY,0,2,N_ELEMENTS(ARRAY))
                        BYTEORDER,LARRAY,/NTOHL
			E0 = '7FF00000'X
			E = LARRAY(0,*) AND E0
			F1 = LARRAY(0,*) AND 'FFFFF'X
			RESULT = WHERE((E EQ E0) AND ((F1 NE 0) OR	$
				(LARRAY(1,*) NE 0)), COUNT)
			END
;
;  Single precision complex floating point.
;
		6:  BEGIN
			LARRAY = LONG(ARRAY,0,2,N_ELEMENTS(ARRAY))
                        BYTEORDER,LARRAY,/NTOHL
			E0 = '7F800000'X
			E1 = LARRAY(0,*) AND E0
			E2 = LARRAY(1,*) AND E0
			F1 = LARRAY(0,*) AND '7FFFFF'X
			F2 = LARRAY(1,*) AND '7FFFFF'X
			RESULT = WHERE(((E1 EQ E0) AND (F1 NE 0)) OR	$
				((E2 EQ E0) AND (F2 NE 0)), COUNT)
			END
		ELSE:  BEGIN
			MESSAGE,'Data type must be floating point',/INFORMATIONAL
			RESULT = -1
			COUNT = 0
			END
	ENDCASE
ind = result
;
; if any found, reassign with user-specified value
;
 if (count gt 0) and (npar eq 4) then begin
     print,'  ',strtrim(string(count),2),' special character(s) found,', $
              ' reassigned value =',strtrim(string(nvalue),2)
     array(ind) = nvalue 
 endif
;
return
end     ;flagnan
