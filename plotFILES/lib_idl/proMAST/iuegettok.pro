;************************************************************************
;+
;*NAME:
;
;	IUEGETTOK	(formerly GETTOK)
;
;*CLASS:
;
;*CATEGORY:
;
;*PURPOSE:
;
;	Procedure to retrieve the first part of the string
;	until the character char (input parameter) is encountered.
;
;*CALLING SEQUENCE:
;
;	IUEGETTOK,ST,CHAR,OUT
;
;*PARAMETERS:
;
;	CHAR	(I) (REQ)    character separating tokens
;
;	ST	(I/O) (REQ)  string to get token from.  On output token
;			     is removed (Leading blanks and tabs are also
;			     removed)
;
;	OUT	(O) (REQ)    taken value
;
;*EXAMPLES:
;
;	idl> st = 'abc=999' 
;	idl> iuegettok,st,'=',new_var 
;	idl> print,new_var
;		 abc
;	idl> print,st
;		 999
;
;*SYSTEM VARIABLES USED:
;
;*INTERACTIVE INPUT:
;
;*SUBROUTINES CALLED:
;
;	PARCHECK
;
;*FILES USED:
;
;*SIDE EFFECTS:
;
;*RESTRICTIONS:
;
;*NOTES:
;
;	tested with IDL Version 2.1.0 (sunos sparc)	21 Jun 91
;	tested with IDL Version 2.1.0 (ultrix mipsel)	N/A
;	tested with IDL Version 2.1.0 (vms vax)		21 Jun 91
;
;*PROCEDURE:
;
;*I_HELP  nn:
;
;*MODIFICATION HISTORY:
;
;	version 1  by D. Lindler APR,86
;	version 2  by D. Lindler Sept'88  - look for TABS
;	 2-21-91 PJL converted to sun/unix; renamed from GETTOK to 
;		 IUEGETTOK; RDAF prolog added; added PARCHECK
;	 5-23-91 PJL changed procedure from a function to a procedure
;        5-24-91 GRA changed tab removal section
;	 6-21-91 PJL tested on SUN and VAX; updated prolog
;
;-
;*************************************************************************
 pro iuegettok,st,char,out
;
; if char is a blank treat tabs as blanks
;
 if n_params(0) eq 0 then begin
    print,'IUEGETTOK,ST,CHAR,OUT'
    retall
 endif  ; n_params(0)
 parcheck,n_params(0),3,'IUEGETTOK'
;
; replace tabs with blanks
;
 bst = byte(st)
 tabs = where(bst eq 9b,count)
 if count ne 0 then bst(tabs) = 32b
 st = string(bst)
;
; find character in string
;
 pos = strpos(st,char)
 if pos eq -1 then begin	;char not found?
    token = st
    st = ''
    out = token
    return
 endif  ; pos
;
; extract token
;
 token = strmid(st,0,pos)
 len = strlen(st)
 if pos eq (len-1) then st = '' else st = strmid(st,pos+1,len-pos-1)
;
; if char is ' ' then remove multiple blanks
;
 if char eq ' ' then st = strtrim(st,1)
 out = token
 return
 end  ; iuegettok
