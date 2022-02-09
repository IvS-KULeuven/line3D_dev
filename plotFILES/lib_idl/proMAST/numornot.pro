;************************************************************************
;+
;*NAME:
;
;	NUMORNOT
;
;*CLASS:
;
;*CATEGORY:
;
;*PURPOSE:
;
;	Determine if a string consists of only numbers (0-9) or not.
;
;*CALLING SEQUENCE:
;
;	NUMORNOT,INSTR,ANSWER,fopt
;
;*PARAMETERS:
;
;	INSTR	(REQ) (I) (string)
;		The string to be checked.
;
;	ANSWER	(REQ) (O) (integer)
;		If equals 0, then the string is not a number.  If equals 1,
;		then the string is a number.
;
;	FOPT	(OPT) (I) (integer)
;		Allows the instr value to be a floating point.
;
;*EXAMPLES:
;
;	numornot,'1234',answer
;	answer = 1
;
;	numornot,'123t',answer
;	answer = 0
;
;	numornot,'4.5',answer
;	answer = 0
;
;	numornot,'4.5',answer,1
;	answer = 1
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
;
;*FILES USED:
;
;	none
;
;*SIDE EFFECTS:
;
;*RESTRICTIONS:
;
;*NOTES:
;
;*PROCEDURE:
;
;	The byte array equivalent of the string is used.  Each element is
;	checked to be sure that it is greater than or equal to 48 (0) and
;	less than or equal to 57 (9).
;
;	If the fopt is set and there is a decimal point in the string, then
;	the string is divided at the decimal point and each part is checked.
;
;*I_HELP  nn:
;
;*MODIFICATION HISTORY:
;
;	12 May 1992  PJL  wrote
;	15 May 1992  PJL  added floating point option
;	19 May 1992  PJL  added check for '-' and '+'
;
;-
;*************************************************************************
 pro numornot,teststr,answer,fopt
;
 npar = n_params(0)
 if (npar eq 0) then begin
    print,'NUMORNOT,INSTR,ANSWER,fopt'
    retall
 endif  ; npar eq 0
 parcheck,npar,[2,3],'NUMORNOT'
;
 instr = teststr
 linstr = strlen(instr)
;
;  see if the first character is a '-' or '+'.  If so, remove it.
;
 if ( (strmid(instr,0,1) eq '-') or (strmid(instr,0,1) eq '+') ) then begin
    instr = strtrim(strmid(instr,1,linstr),2)
    linstr = strlen(instr)
 endif  ; '-' or '+'
;
 sepa = strpos(instr,'.')
 if ( (sepa gt -1) and (npar eq 3) ) then begin
    part1 = strmid(instr,0,sepa)
    part2 = strmid(instr,sepa+1,linstr)
    if (part1 eq '') then count1 = 0 else   $
       type = where( (byte(part1) lt 48) or (byte(part1) gt 57),count1 )
    if (part2 eq '') then count2 = 0 else   $
       type = where( (byte(part2) lt 48) or (byte(part2) gt 57),count2 )
    if ( (count1 ne 0) or (count2 ne 0) ) then answer = 0 else answer = 1
    if ( (part1 eq '') and (part2 eq '') ) then answer = 0
 endif else begin
    type = where( (byte(instr) lt 48) or (byte(instr) gt 57),count )
    if (count ne 0) then answer = 0 else answer = 1
 endelse  ; npar eq 3
;
 return
 end  ; numornot
