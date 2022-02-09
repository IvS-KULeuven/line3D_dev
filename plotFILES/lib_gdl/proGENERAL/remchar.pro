pro remchar,st,char	;Remove character
;+
; Project     : SOHO - CDS
;
; Name        : 
;	REMCHAR
; Purpose     : 
;	Remove all appearances of a character from a string.
; Explanation : 
;	Remove all appearances of character (char) from string (st).
; Use         : 
;	REMCHAR,ST,CHAR
; Inputs      : 
;	ST  - String from which character will be removed.  
;	CHAR- Character to be removed from string. 
; Opt. Inputs : 
;	None.
; Outputs     : 
;	ST	= The modified string is returned in ST.
; Opt. Outputs: 
;	None.
; Keywords    : 
;	None.
; Calls       : 
;	None.
; Common      : 
;	None.
; Restrictions: 
;	None.
; Side effects: 
;	None.
; Category    : 
;	Utilities, strings.
; Prev. Hist. : 
;	Written D. Lindler October 1986
;	Test if empty string needs to be returned   W. Landsman  Feb 1991
; Written     : 
;	Don Lindler, GSFC/HRS, October 1986.
; Modified    : 
;	Version 1, William Thompson, GSFC, 12 April 1993.
;		Incorporated into CDS library.
; Version     : 
;	Version 1, 12 April 1993.
;-
;
bst = byte(st)                                 ;Convert string to byte

bchar = byte(char) & bchar = bchar(0)          ;Convert character to byte

good = where(bst NE bchar,ngood)
if ngood GT 0 then st=string(bst(good)) else st = ''

return
end
