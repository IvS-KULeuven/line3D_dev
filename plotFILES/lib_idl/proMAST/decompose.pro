;******************************************************************************
;+
;*NAME:
;
;	DECOMPOSE     (General IDL Library 01) 19-JUL-81
;
;*CLASS:
;
;	PARSING
; 
;*CATEGORY:
;
;*PURPOSE:
;
;	Break file name into component parts, using VMS and RSX11 punctuation
;	rules.
; 
;*CALLING SEQUENCE:
;
;	DECOMPOSE,FILE,DISK,PATH,NAME,EXTN,VERSION
; 
;*PARAMETERS:
;
;	FILE	(REQ) (I) (0) (S)
;		Required input string variable giving the file name in RSX11
;		and VMS file description format, e.g.
;		diska:[myaccount]file.ext;vers
; 
;	DISK	(REQ) (O) (0) (S)
;		Required output string giving the name of the disk.  If not 
;		specified in FILE, the value returned will be a null string.
;
;	PATH	(REQ) (O) (0) (S)
;		Required output string giving the path name.  If not specified
;		in FILE, the value returned will be a null string.
;
;	NAME	(REQ) (O) (0) (S)
;		Required output string giving the name of the file. 
;
;	EXTN	(REQ) (O) (0) (S)
;		Required output string giving the extension, or file type, 
;		associated with FILE.
;
;	VERS	(REQ) (O) (0) (S)
;		Required output string giving the version number for FILE.
;
;*EXAMPLES:
;
;	in VMS to break up 'DISKA:[IUERDAF.IUELIB.DAT]EBCASC.FIT'
;
;	   decompose,'DISKA:[IUERDAF.IUELIB.DAT]EBCASC.FIT;1',d,path,name,e,v
;
;	where,
;		d = 'DISKA:'
;		path = '[IUERDAF.PRODUCTION.DAT]'
;		name = 'EBCASC'
;		e = '.FIT'
;		v = ';1'
;
;
;	in UNIX or ULTRIX to break up '/home/iuerdaf/production/dat/ebcasc.fit'
;
;	   decompose,'/home/iuerdaf/production/dat/ebcasc.fit',d,path,name,e,v
;
;	where,
;		d =''
;		path = '/home/iuerdaf/production/dat/'
;		name = 'ebcasc'
;		extn = '.fit'
;		vers = ''
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
;       PLATFORM
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
;	All parameters are strings; default parameters are null strings.
;	On UNIX systems, disk and version will always be null strings.
;
;*PROCEDURE: 
;
;	DECOMPOSE looks for the '/' and '.' marks used as delimeters in UNIX
;	file names, and through string manipulations, breaks down the filename
;	into its component parts.
;
;*MODIFICATION HISTORY:
;
;    Jul 19 1981  FHS 3rd  GSFC program written
;    Apr 15 1987  RWT      GSFC add PARCHECK and STRUPCASE
;    Mar  8 1988  CAG      GSFC add VAX RDAF-style prolog
;    ?/89   SUN IDL version written at Univ. of Michigan ?
;    8/31/89 RWT correctly calculates path name 
;    11/21/90 GRA merged vms and unix decompose programs using !version.os 
;     4/10/91 GRA added parcheck to unix section
;     6/21/91 PJL cleaned up; tested on SUN and VAX; updated prolog
;	16 Aug 93  PJL  updated example; changed !version.os if to case
;       14 Jun 94  PJL  replace !version with PLATFORM
;        5 Dec 94  LLT  fix error when null filename entered.
;-
;******************************************************************************
 pro decompose,file,disk,path,name,extn,version
;
 npar = n_params(0)
;
 if (npar eq 0) then begin
    print,'DECOMPOSE,FILE,DISK,PATH,NAME,EXTN,VERSION'
    retall
 endif  ; npar eq 0
 parcheck,npar,6,'DECOMPOSE'
;
;  opbatin system specific information
;
 platform,dummy,syntax=syntax
;
;   file = strupcase(file)
;
 len = strlen(file)
 na = file
;
 disk = ''                      ; default disk name
;;; sc = strpos(na,':')
 sc = strpos(na,syntax.disksep)
 if (sc gt 0) then begin        ; parse disk
    disk = strmid(na,0,sc+1)
    len = len - sc - 1
    na = strmid(na,sc+1,len)
 endif  ; parse disk name
;
 path = ''                      ; default path
;
;;; if (!version eq 'vms') then sc = strpos(na,']') else begin
;;;    a = 0
;;;    pos = 0
;;;    repeat begin
;;;       sc = a
;;;       a = strpos(na,'/',pos)
;;;       pos = a + 1
;;;    end until (a eq -1)
;;; endelse  ; !version.os
 if (syntax.endpath ne ']') then begin
    a = 0
    pos = 0
    repeat begin
       sc = a
       a = strpos(na,syntax.endpath,pos)
       pos = a + 1
    end until (a eq -1)
 endif else sc = strpos(na,syntax.endpath)
;
 if (sc gt 0) then begin        ; parse path
    path = strmid(na,0,sc+1)
    len = len - sc - 1
    na = strmid(na,sc+1,len)
 endif  ; parse path
;
 version = ''                       ; default version
;;; sc = strpos(na,';')
 sc = strpos(na,syntax.versep)
 if (sc gt 0) then begin            ; parse version number
    version = strmid(na,sc,len-sc)
    na = strmid(na,0,sc)
    len = sc
 endif  ; parse version number
;
 extn = ''                          ; default extension
;;; sc = strpos(na,'.')
 sc = strpos(na,syntax.extsep)
 if (sc ge 0) then begin            ; parse extension
    extn = strmid(na,sc,len-sc)
    na = strmid(na,0,sc)
    len = sc
 endif  ; parse extension
;
 name = na                ; remainder after the parse = name!
;
 return
 end  ; decompose
