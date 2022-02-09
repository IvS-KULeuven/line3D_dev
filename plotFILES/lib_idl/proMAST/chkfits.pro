;******************************************************************************
;+
;*NAME:
;
;	CHKFITS
;
;*CLASS:
;
;*CATEGORY:
;
;*PURPOSE:
;
;	Reads the first 80 bytes of a file.  Based on the values of the
;       first six and 30th bytes, CHKFITS determines if the file is a FITS file.
;
;*CALLING SEQUENCE:
;
;	CHKFITS,FILENAME,RESULT,silent=silent
;
;*PARAMETERS:
;
;	FILENAME  (REQ) (I) (0) (S)
;		  The filename (including path if necessary) of the file to be
;		  inspected.
;
;	RESULT	(REQ) (O) (0) (I)
;		The result of the file inspection with the following 
;               definitions:
;                  0 - file is not a standard FITS file
;                  1 - file is a standard FITS file
;                  2 - file is not found 
;                  3 - file exists but can not be opened and/or read.
;
;	SILENT  (KEY) (I) (0) (I)
;		Keyword to prevent messages from being printed to the screen.
;		If not set or set equal to 0, then messages are printed to the
;		screen.  Otherwise, the messages are not printed.
;
;*EXAMPLES:
;
;	chkfits,'swp32525.silo',res
;
;*SYSTEM VARIABLES USED:
;
;	!err_string
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
;	filename in calling sequence
;
;*SIDE EFFECTS:
;
;	none
;
;*RESTRICTIONS:
;
;	none
;
;*NOTES:
;
;       - Although SIMPLE = 'F' can represent a valid FITS format, it is 
;       considered non-standard and RESULT is set to 0.
;       - RESULT = 3 can occur if the input file is already opened, or
;       the file protection prevents read access.
;
;*PROCEDURE:
;
;	The file is opened and the first 80 bytes are read.  The first six (6)
;	characters are compared to 'SIMPLE'.  If they are equal, then the 30th
;	character is compared to 'T'.  If both tests are passed, then the file
;	is a FITS file, and result is set to 1.  If either test is failed, or 
;       an unexpected end-of-file is encountered, then the file is NOT a 
;       (simple) FITS file, and result is set equal to 0.  
;       If the file can not be found, result is set equal to 2, and if any
;       other error occurs during the reading of a file (known to exist), 
;       then result is set to 3. If the keyword silent is set and not equal 
;       to 0, then the error messages are not printed. 
;
;*I_HELP  nn:
;
;*MODIFICATION HISTORY:
;
;	 3 Mar 93  PJL  wrote
;        6 Aug 93  LLT  added code so that program won't bomb for files whose
;                       formats aren't compatible with OPENR and ASSOC stmts.
;       23 Sep 93  RWT  set result = 2 if file is not found
;       13 Jan 94  PJL  replace use of assoc with readu  (assoc has problems
;			with VMS variable length records)
;       11 Apr 94  PJL  added findfile
;       28 Sep 94  PJL  correct silent bug when no files found;  free_lun after
;                       on_ioerror error
;       30 Sep 94  RWT  add result = 3 if file exists but is already open 
;                       and/or can't be read to determine file format
;
;-
;******************************************************************************
 pro chkfits,filename,result,silent=silent
;
 npar = n_params(0)
 if (npar eq 0) then begin
    print,'CHKFITS,FILENAME,RESULT,silent=silent'
    retall
 endif  ; npar eq 0
 parcheck,npar,2,'CHKFITS'
 if keyword_set(silent) then display = 0 else display = 1
;
;  see if file exists
;
 temp = findfile(filename,count=ct)
 if (ct le 0) then begin
    if (display) then print,filename + ' not found.'
    result = 2
    return
 endif  ; ct le 0
;
;  retrieve first 80 bytes
;
 on_ioerror,oops
 get_lun,lun
 openr,lun,filename
 line = bytarr(80)
 readu,lun,line 
 on_ioerror,null
 free_lun,lun
 line = string(line)
;
 if (strmid(line,0,6) eq 'SIMPLE') then begin
    if (strmid(line,29,1) eq 'T') then begin
       if (display) then begin
          print,' '
          print,filename + ' is a FITS file.'
          print,' '
       endif  ; display
       result = 1
    endif else begin
       if (display) then begin
          print,' '
          print,filename + ' is NOT a FITS file.'
          print,"Keyword SIMPLE must equal 'T' (in the 30th column)."
          print,' '
       endif  ; display
       result = 0
    endelse  ; strmid(line,29,1) eq 'T'
 endif else begin
    if (display) then begin
       print,' '
       print,filename + ' is NOT a FITS file.'
       print,"'SIMPLE' must be the first six characters of the file."
       print,' '
    endif  ; display
    result = 0
 endelse  ; strmid(line,0,6) eq 'SIMPLE'
;
 return
 
 oops:
    free_lun,lun
    if (strpos(!err_string,'End of file') gt 0) then begin
       if display then print,!err_string
       result = 0    ; unexpected end of file
    endif else $     
    if (strpos(!err_string,'already open') gt 0) then begin
       if display then print,!err_string
       result = 3    ; logical unit being used
    endif else begin 
       if display then print,!err_string
       result = 3    ; general error reading file
    endelse     
;
 return
 end  ; chkfits
