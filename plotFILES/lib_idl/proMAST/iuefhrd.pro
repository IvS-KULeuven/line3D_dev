;******************************************************************************
;+
;*NAME:
;
;    IUEFHRD
;
;*CLASS:
;
;*CATEGORY:
;
;*PURPOSE:
;
;    To read the header(s) of a FITS disk file, extract the keywords
;    (and IUE label if present) as IDL string arrays, and determine the 
;    overall structure of the file. Up to 15 extension headers can be 
;    extracted, although the entire file is read.
;
;*CALLING SEQUENCE:
;
;    IUEFHRD,FNAME,PARMS,HEADER,LAB,EHD1,EHD2,...,EHD15,SILENT=silent,FULL=full
;
;*PARAMETERS:
;
;    FNAME (REQ) (I) (0) (S)
;       FITS disk file name. If no extension, .FIT is assumed.
;
;    PARMS   (REQ) (O) (1) (1)
;       Vector with up to 100 parameters describing the FITS file format,
;       defined in following order:
;          element              description
;             0       number of records in primary FITS header
;             1       number of data records in primary array 
;             2       number of header records in 1st extension 
;             3       number of data records in 1st extension
;             4       number of header records in 2nd extension 
;             5       number of data records in 2nd extension
;             .....
;
;    HEADER   (REQ) (O) (1) (S)
;       string array of FITS primary header. Up to 900 keywords are allowed.
;
;    LAB      (OPT) (O) (1) (S)
;       string array of IUE VICAR label (including processing history).
;       If no label is present, LAB = ''.
;
;    EHD1    (OPT) (O) (1) (S)
;       optional string array of FITS keywords from 1st extension header 
;
;    EHD2    (OPT) (O) (1) (S)
;       optional string array of FITS keywords from 2nd extension header 
;
;    EHD3    (OPT) (O) (1) (S)
;       optional string array of FITS keywords from 3rd extension header 
;
;    ....
;
;    EHD15    (OPT) (O) (1) (S)
;       optional string array of FITS keywords from 15th extension header 
;
;    SILENT   (KEY) (I) (0) (S)
;       If specified, file structure information will not be displayed
;
;    FULL     (KEY) (I) (0) (S)
;       If specified, label & history portions (if any) will not be removed
;       from output parameter HEADER. Label & history portions will still be
;       written to LAB, if LAB is included in procedure call.
;
;*EXAMPLES:
;    To read headers from an IUE MXLO file:
;       IUEFHRD,'SWP12345.MXLO',PAR,HEADER,LABEL,EXTHEADER
;    To read an MXLO file but keep VICAR label in primary header:
;       IUEFHRD,'SWP12345.MXLO',PAR,HEADER,LABEL,EXTHEADER,/FULL
;
;*SYSTEM VARIABLES USED:
;
;*INTERACTIVE INPUT:
;
;*SUBROUTINES CALLED:
;
;   PARCHECK
;   DECOMPOSE
;   STPAR
;
;*FILES USED:
;
;*SIDE EFFECTS:
;
;*PROCEDURE:
;      
;    FITS primary header file is read until 'END' is found. IUE VICAR label
;    is separated by searching for comment fields 'THE IUE_VICAR HEADER'
;    and 'IUE-LOG FINISHED' and is stored in LAB. Remaining keywords 
;    are stored in string array called HEADER. The required FITS keywords 
;    are read using STPAR. The extension header(s) are extracted and 
;    stored in EXHD1, EXHD2, ..., and EXHD15. The file is read until the
;    end-of-file is found.
;
;*INF_1:
;
;*RESTRICTIONS:
;       Program is currently limited to extracting up to 900 FITS 
;       keywords per header, and up to 15 extension headers.
;
;*NOTES:
;
;      The IUE VICAR label is found by searching for comment fields 
;      'THE IUE_VICAR HEADER' and 'IUE-LOG FINISHED'. If these fields are
;      not found, the program assumes the VICAR label was not included in
;      the FITS header. If only the first string is found, the label is 
;      assumed to continue till the end of the header.
;
;	tested with IDL Version 2.1.0 (sunos sparc)	 7 Aug 91
;	tested with IDL Version 2.3.2 (vms vax)    	15 Dec 92
;
;*MODIFICATION HISTORY:
;
;	written by R Thompson 2/6/91
;	27 Mar 1991 PJL modified for unix/sun
;       30 Jul 1991 RWT include version number in file name (for vms)
;       15 Nov 1991 RWT increase size of HEADER to allow 600 (instead
;              of 300) fits keywords.
;       19 Nov 1991 GRA correct typo in prolog.
;       03 Mar 1992 RWT remove message about no VICAR label found
;                   and correctly handle floating point arrays
;       02 Jul 1992 RWT increase max keyword number to 900
;       08 Sep 1992 RWT read up to 4 extensions
;       16 Sep 1992 RWT make sure end of vicar label (elab) is not greater 
;              than end of header.
;       29 Sep 1992 RWT check record # before searching for extensions
;       15 Dec 1992 RWT add SILENT keyword option
;       29 APR 1993 RWT add FULL keyword option & allow LAB to be optional
;       30 JUL 1997 RWT add test for SIMPLE keyword
;        2 DEC 1997 RWT allow parms to have 100 (rather than 10) entries
;                       & read up to 10 (rather than 4) extensions
;       12 Jan 1998 RWT modify stpar call to always return string values
;       24 Apr 1998 RWT add default values for gcount & pcount (even
;                       though they're required keywords!) & use a
;                       wildcard in the stpar call for reading NAXISn
;                       keywords
;       15 Jul 1998 RWT read up to 15 extensions
;       10 Nov 1998 RWT modify test for END keyword (read all 8 characters)
;       23 Sep 1999 RWT modify to read files with > 2 billion bytes in
;                       data records (e.g., MOSAIC data)
;       25 Oct 1999 RWT make skiprec a longword integer to fix a
;                       bug seen in ifitsrd when reading extensions 
;-
;********************************************************************
 pro iuefhrd,fname,parms,header,lab,ehd1,ehd2,ehd3,ehd4,ehd5,ehd6,ehd7, $
     ehd8,ehd9,ehd10,eh11,eh12,eh13,eh14,eh15,silent=silent,full=full
;
 npar = n_params(0)
   if (npar eq 0) then begin
   print, $
   ' IUEFHRD,FNAME,PARMS,HEADER,lab,ehd1,ehd2,...,ehd15,silent=silent,full=full'
    retall
 endif  ; npar
 parcheck,npar,indgen(17)+3,'IUEFHRD'
;
; process input parameters
;
 decompose,fname,disk,acc,fnam,ext,vers
 if ext eq '' then ext = '.fit'
 fil = disk+acc+fnam+ext+vers
;
; initialize internal parameters
;
 parms = lonarr(100)                     ; array of output parameters
 lab = ' '                               ; IUE VICAR label
 header = strarr(900)		 	 ; fits header array
 ehd1 = ' '                              ; fits extension header
 nhead = 0L      			 ; number of lines in header
 x = bytarr(2880)			 ; all fits records are 2880 bytes
 nr = -1L                                ; record number of fits file
 lastr = 0L                              ; last record flag
 slab = 0L                               ; starting line of IUE label
 elab = 0L                               ; end line of IUE label
 skiprec = 0L                            ; primary data records
 display = 1                             ; display file structure
 separate = 1                            ; separate vicar label & history
 if keyword_set(silent) then display = 0 ; don't display "  "
 if keyword_set(full) then separate = 0 ; don't separate header sections
;
; Open file and check record number
;
 openr,lun,fil,/get_lun,/block
 status = fstat(lun)
 totrec = (status.size) / 2880L           ; total # of records in file
 if display then begin
    print,' '
    print,'FITS File Structure:
    print,'Total number of 2880-byte logical records =',totrec
    print,' '
 endif
;
; Read header of fits file 
;
 rec = assoc(lun,bytarr(2880))	 	 ; define record type
 repeat begin                            ; read records until 'END' found
    nr = nr + 1L
    x = rec(nr)                          ; read next record
    if (nr eq 0) then begin              ; look for SIMPLE keyword
       if (string(x(0:5)) ne 'SIMPLE') then begin
          print,' Not a valid FITS file, returning'
          retall
       endif else if (string(x(29)) ne 'T') then $
             print,' Not a standard FITS file, but continuing'
    endif
    i = -1                               ; keyword number in record
    repeat begin                         ; process keywords til 'END' found
       i = i + 1L
       st = i * 80L
       h = x(st:st+79L)                   ; extract next entry
       sh = string(h)                    
       kw = strtrim(string(h(0:7)))      ; keyword name
       header(nhead) = sh                ; add to header string array
       if (kw eq 'END') then lastr = 1  ; check for end of header
       if lastr eq 0 then begin          ; find specific comments
          if (strpos(sh,'THE IUE VICAR HEADER',8) ne -1) then slab = nhead
          if (slab ne 0) then if (strpos(sh,'IUE-LOG FINISHED',8) ne -1)  $
             then elab = nhead
       endif  ; lastr
       nhead = nhead + 1L
    endrep until (lastr) or (i eq 35)        ; continue up to 36 times
 endrep until (lastr) or ((nr+1) ge totrec)  ; 'end' (or end of file) found 
 parms(0) = nr+1L                            ; number of header records
 header = header(0:nhead-1)
;
;  if requested, separate keywords and IUE VICAR label 
;
 if ((slab+elab) ne 0) and (npar gt 3) then begin
    if (elab eq 0) then elab = nhead-2
    lab = header(slab-1:elab) 
    newh = [header(0:slab-2),header(elab+1:*)]
    nhead = n_elements(newh)
    if (separate) then header = newh
 endif else if ((display) and (separate)) then print,' no vicar label found'
 if display then begin
   print,' PRIMARY ARRAY'
   snh = strtrim(string(nhead),2)
   print,'      header:',nr+1L,' records with ',snh,' keywords'
 endif
;
; get needed keywords from header
;
 stpar,header,'bitpix',bitpix
 bpp = abs(bitpix)
 stpar,header,'naxis',naxis
 stpar,header,'gcount',gcount
 stpar,header,'extend',extend,err
 if err lt 0 then gcount=1
 stpar,header,'psize',psize,err
 if err lt 0 then psize=0
 if naxis gt 0 then begin
    stpar,header,'naxis*',naxp
    npoints=1L
    for i=0L,naxis-1 do npoints = npoints * naxp(i)
 endif else npoints=0L
 skiprec = long((npoints/8)*bpp) 
 skiprec = long(( skiprec / 2880.0 + ((skiprec mod 2880) ne 0)))
 parms(1) = skiprec
 if display then print,'      data:  ',long(skiprec),' records'
 nr = nr + skiprec 
 pnum = 1
;
; if more records exist , read headers & skip data 
; store number of header & data records in parms
;
 while ((nr+1) lt totrec) do begin       ; search for extensions
    exhd = strarr(900)                   ; fits extension header
    tmp = nr
    nhead = 0L                           ; reinitialize 
    lastr = 0L
    tablename = ' '
    repeat begin                         ; read until 'END'
       nr = nr + 1L
       x = rec(nr)                       ; read next record
       i = -1L
       repeat begin                      ; process 36 header lines
          i = i + 1L
          st = i * 80L
          h = x(st:st+79)                ; extract next keyword
          sh = string(h)                 
          kw = strtrim(string(h(0:7)))   ; keyword name
          exhd(nhead) = sh               ; add to header array
          nhead = nhead + 1L
          if (kw eq 'END') then lastr = 1 ; check for end of header
       endrep until (lastr) or (i eq 35)    ; read keywords til 'end' or i=35
       if ((nr-tmp) eq 1) then stpar,exhd,'xtension',tablename,errx,typecode=7
    endrep until (lastr) or ((nr+1) ge totrec) or (tablename eq ' ')
                                               ; read records till 'end' or eof 
    stpar,exhd,'extname',extname,err
    if (err lt 0) then extname = '' 
    if (errx ge 0) and (display) then $
        print,' ',strtrim(tablename,2),' extension:  '+extname 
    if (errx lt 0) then begin
        tablename = ' '
        if (display) then print,' non-standard extension found'
    endif
    if (tablename ne ' ') then begin
       pnum = pnum + 1
       parms(pnum) = nr - tmp                  ; # of extension header records
       snh = strtrim(string(nhead),2)
       if display then print, $
          '      header:',long(nr-tmp),' records with ',snh,' keywords'
       if (pnum eq 2) then ehd1 = exhd(0:nhead-1)
       if (pnum eq 4) then ehd2 = exhd(0:nhead-1)
       if (pnum eq 6) then ehd3 = exhd(0:nhead-1)
       if (pnum eq 8) then ehd4 = exhd(0:nhead-1)
       if (pnum eq 10) then ehd5 = exhd(0:nhead-1)
       if (pnum eq 12) then ehd6 = exhd(0:nhead-1)
       if (pnum eq 14) then ehd7 = exhd(0:nhead-1)
       if (pnum eq 16) then ehd8 = exhd(0:nhead-1)
       if (pnum eq 18) then ehd9 = exhd(0:nhead-1)
       if (pnum eq 20) then ehd10 = exhd(0:nhead-1)
       if (pnum eq 22) then ehd11 = exhd(0:nhead-1)
       if (pnum eq 24) then ehd12 = exhd(0:nhead-1)
       if (pnum eq 26) then ehd13 = exhd(0:nhead-1)
       if (pnum eq 28) then ehd14 = exhd(0:nhead-1)
       if (pnum eq 30) then ehd15 = exhd(0:nhead-1)
;
;      read header to determine data records to skip
;
       stpar,exhd,'bitpix',bitpix
       bpp = abs(bitpix)
       stpar,exhd,'naxis',naxis
       stpar,exhd,'gcount',gcount,err
       if (err lt 0) then begin
          print,'      assume gcount = 1'
          gcount = 1
       endif
       stpar,exhd,'pcount',pcount,err
       if (err lt 0) then begin
          print,'      assume pcount = 0'
          pcount = 0
       endif
       stpar,exhd,'psize',psize,err
       if err lt 0 then psize=0
       if naxis gt 0 then begin
          stpar,exhd,'naxis*',naxp
          npoints=1L
          for i=0,naxis-1 do npoints = npoints * naxp(i)
       endif else npoints=0
       skiprec = long(((npoints/8)*bpp + pcount) * gcount) 
       skiprec = long( skiprec / 2880.0 + ((skiprec mod 2880) ne 0))
       if display then print,'      data:  ',long(skiprec),' records'
       pnum = pnum + 1
       parms(pnum) = skiprec
       nr = nr + skiprec 
    endif
 endwhile       ; end extension search
 parms = parms(0:pnum)     ; remove trailing zeroes
 free_lun,lun
;
 return
 end  ; iuefhrd
