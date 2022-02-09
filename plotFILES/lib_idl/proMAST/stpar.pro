;******************************************************************************
;+
;*NAME:
;
;	STPAR
;
;*CLASS:
;
;*CATEGORY:
;
;	FITS
;
;*PURPOSE:
;
;	Obtain the value(s) for the specified keyword from a fits 
;       header.
;
;*CALLING SEQUENCE:
;
;	STPAR,HDR,NAME,RESULT,err,nfields=nfields,valid=valid,typecode=typecode
;
;*PARAMETERS:
;
;	HDR	(REQ) (I) (1) (S)
;		FITS header stored as an IDL string array.
;
;	NAME	(REQ) (I) (0) (S)
;		String name of the parameter to return. If NAME is of
;               form 'keyword*' then an array is returned containing the
;               values of KEYWORDn where n is an integer. The data type
;               of RESULT is determined by the first valid match.
;		
;	RESULT	(REQ) (O) (0,1) (I L R D S)
;		Scalar value of parameter.  If parameter is double precision,
;               floating, long or string, the result is of that type.  
;               Apostrophes are stripped from strings.  If the parameter 
;               is logical, 1 is returned for T, and 0 is returned for F.
;               If NAME is of form 'keyword*' then RESULT will be a vector
;               of size n where n is the largest number (from KEYWORDn
;               parameters located).
;
;	ERR	(opt) (O) (0) (I)
;		Status of errors. A value of 0 indicates paramter was found,
;               and -1 means it was not found. If multiple values are
;               returned (i.e., when NAME = keyword*), err specifies the 
;               number of entries found. 
;
;       nfields (key) (i) (0) (i)
;               For NAME of form 'keyword*' this is the number of elements in
;               the output array.  Note that ERR will still specify the number
;               entries actually found, not including the padded values.  The
;               MAXIMUM of nfields and the largest n (from KEYWORDn) value is
;               used and if n is greater than nfields, nfields is reset to n.
;
;       valid   (key) (o) (0,1) (i)
;               For vector results, this contains the indices of the valid
;               values (due to the possible presence of padded zeroes).  For
;               scalar results, it just equals ERR.
;
;      typecode (key) (i) (0) (i)
;               If no matches are found, RESULT is normally 0 (integer).  To
;               return RESULT as a data type other than integer, set typecode
;               to the type code as used by the SIZE function.  For example,
;               to return a null string instead of 0, typecode=7.  If matches
;               are found, the data type of RESULT is determined automatically
;               and this keyword is ignored.
;
;*EXAMPLES:
;
;	readri,'swp23456.rilo',main,image 
;	stpar,main,'TELESCOP',result,err  ; return value of TELESCOP
;         if (err lt 0) then print,' telescop keyword not found'
;       stpar,main,'naxis*',result,err    ; return all NAXISn keyword values
;       stpar,main,'comment',result       ; return all comments in header
;
;       iuefhrd,filename,p,h,eh                    ;Read extension header
;       stpar,eh,'tform*',tform,ntform             ;TFORM is required
;       stpar,eh,'tdim*',tdim,ntdim,nfields=ntform ;TDIM is not required and
;                                                  ;may or may not be present.
;                                                  ;Make sure the result is the
;                                                  ;same size as TFORM.
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
;       NUMORNOT
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
;	tested with IDL Version 2.1.0 (sunos sparc)	 7 Aug 91
;	tested with IDL Version 2.1.0 (vms vax)    	 7 Aug 91
;
;*PROCEDURE:
;
;	The first 8 characters of each element of HDR are searched for a match
;	to NAME.  The value from the next 20 characters is returned unless a
;       commentary keyword is specified in which case the next 72 characters
;       are returned.  An error occurs if there is no parameter with the 
;       given name.
;
;       The COMMENT, HISTORY, and blank commentary keywords are treated 
;       differently than other keywords. All entries found and all characters 
;       after byte 8 are returned as a string array.
;
;       If a numeric value has no decimal point it is returned as type long.
;       If it contains >8 characters or contains the letter 'D' then it is
;       returned as type DOUBLE. Otherwise it is type FLOAT.
;
;*I_HELP  nn:
;
;*MODIFICATION HISTORY:
;
;	DMS, May, 1983, Written.
;	D. Lindler Oct 86 Returns !err=-1 if name not found,
;			and allows longer string values
;	27 Mar 1991	PJL modified for unix/sun; added PARCHECK;
;			replaced !ERR with ERR
;	10 May 1991	PJL corrected prolog format
;	29 Sep 93  PJL  prolog
;       06 Sep 94  RWT  add updates as contained in UIT's SXPAR routine 
;                       including support for double precision keyword values, 
;                       commentary keywords, and use of wildcards, but 
;                       removed strtrim for commentary keywords, and 
;                       corrected test for blank keywords.
;       23 Jan 95  LLT  add nfields, typecode, and valid keywords.
;       30 Jan 95  RWT  assume exponential notation using E implies floating 
;                       point (not double) precision
;       19 Jul 95  RWT  remove /nozero keyword from make_array command to fix 
;                       unix bug when wild cards are used.
;-
;******************************************************************************
 pro stpar,hdr,name,output,err,nfields=nfields,valid=valid,typecode=typecode
;
 if n_params(0) eq 0 then begin
   print,'STPAR,HDR,NAME,OUTPUT,err,nfields=nfields,valid=valid,typecode=typecode'
    retall
 endif  ; n_params(0)
 parcheck,n_params(0),[3,4],'STPAR'
;
 on_error,2		; return to caller if error.
 if n_params(0) eq 0 then begin
   print,'STPAR,HDR,NAME,OUTPUT,err,nfields=nfields,valid=valid,typecode=typecode'
    retall
 endif  ; n_params(0)
 parcheck,n_params(0),[3,4],'STPAR'
;
; check for valid header
;
 s = size(hdr)		; check header for proper attributes.
 if (s(0) ne 1) or (s(2) ne 7) then begin
    print,'STPAR - Hdr array is not the correct type'
    output = 0
    return
 endif  ; strlen(hdr)
;
; determine if name is of form keyword*
;
 nam = strtrim(strupcase(name))   	; copy name, make upper case
 name_length = strlen(nam) - 1
 if ( strpos(nam,'*') eq (name_length > 1) ) then begin
    nam = strmid(nam,0,name_length)
    vector = 1
    num_length = 8 - name_length
    if num_length le 0 then $
       print,' Keyword length must be 8 characters or less'
 endif else begin
    while strlen(nam) lt 8 do nam = nam + ' ' ;make 8 chars long
    vector = 0
 endelse
;
; find number of occurances
;
 keyword = strmid(hdr,0,8)
 histnam = (nam eq 'HISTORY ') or (nam eq 'COMMENT ') or (nam eq '        ')
 if vector then begin
    nfound = where(strpos(keyword,nam) ge 0, matches)
    if (matches gt 0) then begin
       numst = strmid(hdr(nfound),name_length,num_length)
       number = intarr(matches) - 1
       for i=0,matches-1 do begin
           numornot,strtrim(numst(i),2),num,1
           if (num) then number(i) = double(strtrim(numst(i),2))
       endfor
       igood = where(number ge 0,matches)
       if (matches gt 0) then begin
          nfound = nfound(igood) 
          number = number(igood)
       endif
    endif else number=-1
    if not keyword_set(nfields) then nfields=1
    nfields=max(number)>nfields 
    valid=number-1
 endif else begin
    nfields=1
    nfound = where(keyword eq nam, matches)
    if not histnam then if matches gt 1 then print, $
       'Warning - Keyword '+nam+' located more than once in hdr'
 endelse
;
; Loop through lines of header
;
 if (matches gt 0) then begin    
    line = hdr(nfound)
    svalue = strtrim( strmid(line,9,71),2)   ; rest of line after keyword
    if histnam then value = strmid(line,8,72) else $
    for i = 0,matches-1 do begin
       if ( strmid(svalue(i),0,1) EQ "'" ) then begin   ;Is it a string?
          test = strmid( svalue(i),1,strlen( svalue(i) )-1)
	  next_char = 0
	  value = '' 
          NEXT_APOST:
	  endap = strpos(test, "'", next_char)      ;Ending apostrophe  
	  if endap LT 0 then $ 
	      print,'Value of '+name+' invalid in hdr'
	  value = value + strmid( test, next_char, endap-next_char )  
;
;  Test to see if the next character is also an apostrophe.  If so, then the
;  string isn't completed yet.  Apostrophes in the text string are signalled as
;  two apostrophes in a row.
;
          if strmid( test, endap+1, 1) EQ "'" then begin    
	     value = value + "'"
             next_char = endap+2         
             goto, NEXT_APOST
          endif      
;
; Process non-string value  
;
       endif else begin
	  value = strtrim( strmid(line(i), 10, 20), 2)   ;Extract value    
	  if ( value EQ 'T' ) then value = 1 else $
	  if ( value EQ 'F' ) then value = 0 else begin
;
;  Test to see if a complex number.  It's not a complex number if columns 31-50
;  are blank, or if a slash character occurs in columns 11-50.  Otherwise, try
;  to interpret columns 31-50 as a number.
;
	     value2 = strtrim( strmid(line(i),30,20), 2)  ; Imaginary part
             if (strpos( strmid( line(i),10,40), '/') GE 0) or $
             (value2 EQ '') then begin           ; not complex
               On_ioerror, GOT_VALUE
               if strpos(value,'.') GE 0 then begin      
                 if ( strpos(value,'D') GT 0 ) then value = double(value) else $
                 if ( strpos(value,'E') GT 0 ) then value = float(value) else $
                 if ( strlen(value) GE 8 ) then value = double(value) $
                 else value = float(value)
               endif else value = long(value)
             endif else begin                    ; complex
                  On_ioerror, NOT_COMPLEX
                  value2 = float(value2)
                  value = complex(value,value2)
             endelse
             goto, GOT_VALUE
;
;  not complex number(?)  Decide if it is a floating point, double precision,
;  or integer number (again).
;
             NOT_COMPLEX:
	     On_ioerror, GOT_VALUE
	     if strpos(value,'.') GE 0 then begin      
               if ( strpos(value,'D') GT 0 ) then value = double(value) else $
               if ( strpos(value,'E') GT 0 ) then value = float(value) else $
	       if ( strlen(value) GE 8 ) then value = double(value) $
               else value = float(value)
             endif else value = long(value)
;
             GOT_VALUE:
	     On_IOerror, NULL
	   endelse   ; value ne F
        endelse      ; if c eq apost
;
;  output value as scalar or vector 
;
	if vector then begin
	   if (i eq 0) then begin
               sz_value = size(value)
               output = make_array( nfields, type=sz_value(1))
	   endif 
           output( number(i)-1 ) =  value
	endif else output = value
    endfor
    if vector then err = matches else begin
       err = 0
       output = value
    endelse
 endif else  begin
    err = -1
    if not keyword_set(typecode) then typecode=2
    output=make_array(nfields,type=typecode)
    if not vector then output=output(0)
    valid=err
 endelse
 return
 end    ;stpar
