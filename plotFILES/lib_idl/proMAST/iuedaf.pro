;******************************************************************************
;+
;*NAME:
;
;	IUEDAF
;
;*CLASS:
;
;*CATEGORY:
;
;	NEWSIPS
;
;*PURPOSE:
;
;	Determine if the IUEDAC keyword has been added to the main fits header.
;	If not, it attempts to add it.
;
;*CALLING SEQUENCE:
;
;	IUEDAF,MAIN,FLAG,/date,/silent,/override
;
;*PARAMETERS:
;
;	MAIN	(REQ) (I/O) (1) (S)
;		The main fits header (may include the vicar label and the
;		history portion of the fits header).
;
;	FLAG	(REQ) (O) (0) (I)
;		Set to 1 if the IUEDAC keyword exists or has been added.  Set
;		to 0 if the IUEDAC keyword does not exist and there is an
;		error when attempting to add it.
;
;	DATE	(KEY) (I) (0) (I)
;		If set, then the fits keyword EXTDATE will be added/updated
;		in the main fits header.
;
;	SILENT	(KEY) (I) (0) (I)
;               If set, program will not print out message if IUEDAC section
;               already exists.
;
;	OVERRIDE  (KEY) (I) (0) (I)
;		If set, procedure will not insist on placing the IUEDAC section
;		before the line 'COMMENT * CORE DATA ITEMS - COMMON SET'.
;
;*EXAMPLES:
;
;	iuefhrd,filename,params,main,/silent,/full
;	iuedaf,main,flag,/date
;
;*SYSTEM VARIABLES USED:
;
;	!stime
;
;*INTERACTIVE INPUT:
;
;	none
;
;*SUBROUTINES CALLED:
;
;	PARCHECK
;	STPAR
;	ADDPAR
;       DATECONV
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
;       If /override is not specified, then the line 
;	'COMMENT * CORE DATA ITEMS - COMMON SET' must exist or the
;       program will abort.
;
;*PROCEDURE:
;
;	STPAR is used to determine if the IUEDAC keyword is already in the
;	fits header.  If it does not exist, then the line
;	'COMMENT * CORE DATA ITEMS - COMMON SET' is searched for.  If the line
;	is not found and /override is not specified, the IUEDAC keyword is 
;       not added and the procedure returns (flag = 0). If the line is found
;       or /override is specified, then the IUEDAC keyword and four
;	comment lines are added two lines above that line (or at the end of
;       the header.  If the DATE keyword is set, then the fits keyword 
;       EXTDATE is added or updated.
;
;*I_HELP  nn:
;
;*MODIFICATION HISTORY:
;
;	 2 Jun 93  PJL  wrote
;	 3 Jun 93  PJL  add date keyword
;       10 Jun 93  RWT  add silent keyword
;	20 Jul 93  PJL  added override keyword
;	 6 Aug 93  PJL  expanded use of silent keyword
;        5 Sep 93  RWT  use DATECONV to correct error in calculating date
;                       and update prolog.
;       18 Oct 93  RWT  correct bug in calculating extdate when IUEDAC 
;                       section is found
;        1 Apr 98  RWT  correct bug in adding IUEDAC keyword 
;       04 Jan 99  RWT  use NEWFITS mode in DATECONV
;-
;******************************************************************************
 pro iuedaf,main,flag,date=date,silent=silent,override=override
;
 npar = n_params(0)
 if (npar eq 0) then begin
    print,'IUEDAF,MAIN,FLAG,/date,/silent,/override'
    retall
 endif  ; npar eq 0
 parcheck,npar,2,'IUEDAF'
;
 flag = 0
;
;  make sure that IUEDAC section does not exist
;
 stpar,main,'iuedac',junk,check
;
 if (check eq 0) then begin
;
;  section does exist
;
    if (not (keyword_set(silent))) then begin
       print,'The IUEDAC section has been added to the main fits header.'
       print,'It will not be added again.'
    endif  ; not (keyword_set(silent)
 endif else begin
;
;  add IUEDAC section
;
    findline = 'COMMENT * CORE DATA ITEMS - COMMON SET'
    linenum = where(strtrim(main,2) eq findline,ct)
;
;  can not find start line
;
    if (ct le 0) then begin
       if (keyword_set(override)) then begin
          findline = 'END'
          linenum = where(strtrim(main,2) eq findline,ct2)
          case (ct2) of
             0:  begin
                    print,' '
                    print,"Can not locate 'END' keyword."
                    print,'ACTION:  Returning'
                    return
                 end  ; ct eq 0
             1:  linenum = linenum(0) + 1
             else:  linenum = max(linenum) + 1
          endcase  ; ct2
       endif else begin
          if (not (keyword_set(silent))) then begin
             print,"Can not find the 'COMMENT * CORE DATA ITEMS - " +   $
                "COMMON SET' line."
             print,'Can not determine where to add the IUEDAC section.'
             print,'ACTION:  Returning'
          endif  ; not (keyword_set(silent))
          return
       endelse  ; not (keyword_set(override))
    endif  ; ct
;
    if (linenum(0) lt 2) then begin
       if (not (keyword_set(silent))) then begin
          print,'The line ' + findline + ' is oddly located.'
          print,'This variable may not be a fits header.'
          print,'You may want to look at the fits header.'
          print,'ACTION:  Returning'
       endif  ; not (keyword_set(silent))
       return
    endif  ; linenum(0) lt 2
;
;  make new main header array
;
    smain = size(main)
    newmain = strarr(smain(1) + 5)
    newmain(0) = main(0:linenum(0)-2)
    newmain(linenum(0)+4) = main(linenum(0)-1:smain(1)-1)
;
;  add IUEDAC section
;
    blank = string(replicate(32b,80))           ; blank line
;
    linec = blank
    strput,linec,'COMMENT'
    strput,linec,string(replicate(42b,70)),8
    newmain(linenum(0)-1) = linec
    newmain(linenum(0)+3) = linec
;
    linec = blank
    strput,linec,'COMMENT START IUEDAC PROCESSING INFORMATION'
    newmain(linenum(0)) = linec
;
    linec = blank
    strput,linec,"IUEDAC  = 'END KEYWORDS'"
    newmain(linenum(0)+1) = linec
;    addpar,temp,'IUEDAC','END KEYWORDS'
;    newmain(linenum(0)+1) = temp(0)
;
    linec = blank
    strput,linec,'COMMENT END IUEDAC PROCESSING INFORMATION'
    newmain(linenum(0)+2) = linec
;
    main = newmain
 endelse  ; check
 flag = 1
;
 if (keyword_set(date)) then begin        ; add/update extdate keyword
       dateconv,!stime,'newfits',out
       addpar,main,'extdate',out, $
       ' date data was extracted from FITS file','iuedac' 
 endif  ; keyword_set(date)
;
 return 
 end  ; iuedaf

