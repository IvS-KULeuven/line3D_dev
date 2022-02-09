;************************************************************************
;+
;*NAME:
;
;	ADDPAR
;
;*PURPOSE:
;
;	Add or modify a parameter in a FITS header array.
;
;*CALLING SEQUENCE:
;
;	ADDPAR, HEADER, NAME, VALUE, comment, location, format=format,
;               before = before, after = after
;
;*PARAMETERS:
;
;	HEADER    (REQ)  (IO)  (1)  (S)
;             String array containing FITS header. Max string length 
;             must be equal to 80. If not defined, then ADDPAR will 
;             create an empty FITS header array.
;
;	NAME      (REQ)  (I)  (0)  (S)
;             Name of parameter. If Name is already in the header the 
;             value and possibly comment fields are modified. Otherwise 
;             a new record is added to the header. Commentary keywords,
;             (i.e. 'HISTORY', 'COMMENT' or ' ') will be added to the 
;             header without replacement.  In these cases the comment 
;             parameter is ignored.
;
;	VALUE     (REQ)  (I)  (0)  (FIS)
;             Value for parameter.  The value expression must be of the 
;             correct type, e.g. integer, floating or string.  String 
;             values of 'T' or 'F' are considered logical values.
;
;	comment   (OPT)  (I)  (0)  (S)
;             String field.  The '/' is added by this routine.  
;             Added starting in position 31. If not supplied, or set 
;             equal to '', then the previous comment field is retained 
;             (when found). The comment field is ignored for commentary
;             keywords.
;
;	location  (OPT)  (I)  (0)  (S)
;              Keyword string name.  The parameter will be placed
;	       before the first location of this keyword.  For example, if 
;              Location = 'HISTORY' then the parameter will be placed 
;              before the first history location.   This applies only 
;              when adding a new keyword; keywords already in the header 
;              (except commentary keywords) are kept in the same position.
;              Note this parameter is identical to using the BEFORE keyword
;              described below.
;
;	format    (KEY)  (I)  (0)   (S)  
;              Keyword for specifying format for parameter.  A scalar
;              string should be used.
;
;	before    (KEY)  (I)  (0)   (S)  
;              The Keyword specified by NAME is inserted before the location
;              of this keyword. For example, before='extend' will place
;              the new keyword before the extend keyword.
;
;	after     (KEY)  (I)  (0)   (S)  
;              Same as above except keyword specified by NAME is inserted
;              after this keyword. (If optional parameter LOCATION is 
;              specified, the default is to place the new keyword BEFORE 
;              not AFTER.)
;
;*SUBROUTINES CALLED:
;
;	PARCHECK
;
;*COMMON BLOCKS:
;
;*SIDE EFFECTS:
;
;*RESTRICTIONS:
;
;	Warning -- Parameters and names are not checked
;		against valid FITS parameter names, values and types.
;
;*NOTES:
;
;       If location is specified, keywords BEFORE and AFTER are ignored,
;       and the keyword will be put before location. If BEFORE is specified, 
;       then AFTER would be ignored. If none of the above are specified, 
;       keyword is written at end of input header.
;
;       If LOCATION, BEFORE, or AFTER refer to a commentary keyword,
;       the first occurance will be used as the location for the 
;       new keyword.
;
;	tested with IDL version 2.1.2  (sunos sparc)    12 Nov 1991
;       tested with IDL version 2.1.2  (ultrix mipsel)  12 Nov 1991
;       tested with IDL version 2.2.0  (ultrix vax)     12 Nov 1991
;       tested with IDL version 2.1.2  (vms vax)        14 Dec 1993
;
;*PROCEDURE:
;
;	Straightforward. Note that the special commentary keywords (i.e., 
;       'HISTORY', 'COMMENT' and ' ') are treated slightly differently 
;       in that they will not overwrite an existing 
;       entry. This allows multiple commentary keywords to be included 
;       within a FITS header. Also, the '=' and '/' characters are not 
;       included for these keywords, the VALUE paramter is written as a
;       comment, and the COMMENT parameter is ignored.
;       
;
;*MODIFICATION HISTORY:
;
;	DMS, RSI, July, 1983.
;	D. Lindler Oct. 86  Added longer string value capability
;	Converted to NEWIDL  D. Lindler April 90
;       Added Format keyword, J. Isensee, July, 1990
;       12 Nov 91  GRA updated IUERDAF version of sxaddpar.pro, which
;                      fixed errors when updating existing fits header
;                      values; added PARCHECK; cleaned up; updated prolog.
;       13 Dec 93  RWT allow multiple commentary keywords to be added 
;                      anywhere, (without '=' or '/' characters), and add
;                      BEFORE and AFTER parameters as suggested by
;                      K. Venkatakrishna (May, 92).
;       22 Feb 99  RWT add format keyword in listed calling sequence
;-
;************************************************************************
 pro addpar, header, name, value, comment, location, format=format, $
             before = before, after = after
;
 on_error,2                                  ; return to caller
;
 npar = n_params()
 if npar eq 0 then begin                    
   print,' ADDPAR,HEADER,NAME,VALUE,comment,location,format=format,'
   print,' before=before,after=after'
   return
 endif  ; npar
 parcheck,npar,[3,4,5],'ADDPAR'
;
; define a blank line and the end line
;
 endline = 'END' + string(replicate(32b,77)) ; end line
 blank = string(replicate(32b,80))	     ; blank line
;
;  if location parameter not defined, set it equal to 'END     '
;
 if n_params() gt 4 then loc = strupcase(location) else $
 if keyword_set(before) then loc = strupcase(before) else $
 if keyword_set(after) then loc = strupcase(after) else $
 loc = 'END     '
 while strlen(loc) lt 8 do loc = loc + ' '
;
 if n_params() lt 4 then comment = ''        ; is comment field specified?
;
 n = n_elements(header)	                     ; # of lines in fits header
 if n eq 0 then begin	                     ; header defined?
    header = strarr(10)                      ; no, make it.
    header(0) = endline
    n = 10
 endif else begin
    s = size(header)                         ; check for string type
    if (s(0) ne 1) or (s(2) ne 7) then begin
       print,'FITS header (first parameter) must be a string array'
       return
    endif  ; s(0) ne 1 or s(2) ne 7
 endelse  ; n eq 0
;
;  make sure name is 8 characters long
;
 nn = string(replicate(32b,8))	             ; 8 char name
 strput,nn,strupcase(name)                   ; insert name
;
;  extract first 8 characters of each line of header, and locate END line
;
 keywrd = strmid(header,0,8)                 ; header keywords
 iend = where(keywrd eq 'END     ',nfound)
 if nfound eq 0 then header(0) = endline     ; no end, insert at beginning
 iend = iend(0) > 0                          ; make scalar
;
; find location i to insert keyword
; (but don't replace commentary keywords)
;
 if (nn eq 'HISTORY ') or (nn eq 'COMMENT ') or (nn eq '        ') then begin
    nfound = 0
    nohc = 0
 endif else begin
    ipos  = where(keywrd eq nn,nfound)
    nohc = 1
    if nfound gt 0 then begin
       i = ipos(0)
       if comment eq '' then comment=strmid(header(i),32,48)  ; save comment?
       goto,replace
    endif
 endelse
 if loc ne '' then begin
    iloc =  where(keywrd eq loc,nloc)
    if nloc gt 0 then begin
       i = iloc(0)
       if keyword_set(after) then i = i+1 < iend
       if i gt 0 then header=[header(0:i-1),blank,header(i:n-1)]  $
       else header=[blank,header(0:n-1)]
       goto,replace
    endif  ; nloc gt 0
 endif  ; loc ne ''
;
; at this point keyword and location parameters were not found, so a new
; line is added at the end of the fits header
;
 if iend lt (n-1) then begin	             ; not found, add more?
    header(iend+1) = endline	             ; no, already long enough.
    i = iend		                     ; position to add.
 endif else begin		             ; must lengthen.
    header = [header,replicate(blank,5)]     ; add an element on the end
    header(n) = endline		             ; save "END"
    i = n-1			             ; add to end
 endelse
;
; now put value into keyword at line i
;
 replace: 
    h=blank			             ; 80 blanks
    if (nohc) then strput,h,nn+'= ' else strput,h,nn   ; insert name (and '='?)
    apost = "'"	                             ; quote a quote
    type = size(value)	                     ; get type of value parameter
    if type(0) ne 0 then print,'keyword value (3rd parameter) must be scalar'
   ;
    case type(1) of		             ; which type?
       7: begin
             upval = strupcase(value)	     ; force upper case.
	     if (upval eq 'T') or (upval eq 'F') then strput,h,upval,29 $
             else begin                        ; not T or F
               if (nohc) then begin            ; not history or comment
                  if strlen(value) gt 18 then begin	; long string
  		    strput,h,apost+value+apost+' /'+comment,10
		    header(i) = h
		    return
                  endif  ; strlen(value) gt 18
                  strput,h,apost+strmid(value,0,18)+apost,10 ;insert string val
               endif else strput,h,value,10    ; history or comment keyword
             endelse                           ; not T or F
	  endcase  ; 7
    else: begin
	     if (n_elements(format) eq 1) then begin    ; use format keyword
	         v = string(value,'('+format+')' ) 
             endif else v = strtrim(value,2)      ; convert to string, default 
	     s = strlen(v)                   ; right justify
             strput,h,v,(30-s)>10            ; insert
	  endcase  ; else
    endcase  ; type(1)
   ;
    if (nohc) then strput,h,' /',30	     ; add ' /'
    if (nohc) then strput,h,comment,32	     ; add comment
    header(i)=h		                     ; save line
;
 return
 end
