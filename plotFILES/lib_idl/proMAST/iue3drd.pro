;******************************************************************************
;+
;*NAME:
;
;    IUE3DRD
;
;*CLASS:
;
;*CATEGORY:
;
;*PURPOSE:
;
;    To read up to 20 columns from the specified row of a FITS binary table. 
;
;*CALLING SEQUENCE:
;
;    IUE3DRD,FILENAME,ROW,SREC,B3DH,COL1,col2,...col20,silent=silent,efld=efld
;           ,npts=npts
;
;*PARAMETERS:
;
;    FILENAME (REQ) (I) (0) (S)
;       IUE FITS disk file name. If no extension, .fit is assumed.
;       If specified as a number rather than a string, it will interpret 
;       the value as a logical unit number for a FITS file already opened. 
;       (This can reduce execution time for repeated calls.)
;
;    ROW      (REQ) (I) (0) (BIL)
;       Row number to be read (e.g., first row = 1).
;   
;    SREC      (REQ) (I) (0) (BIL)
;      starting data record for 3D table (assumes 2880 byte records,
;      and first record = 0).
;
;    B3DH     (REQ) (I) (1) (S)
;      string array containing binary table extension header keywords.
;      (created by routine iuefhrd.pro).
;
;    COLn     (OPT) (O) (01) (BILRS)
;       Extracted fields (up to 20) from binary table. Data type and size 
;       are determined by table format. Sclars or vectors are allowed.
;       Columns are read sequentially starting with first field or, if
;       EFLD is specified as a vector, the column numbers listed in EFLD
;       are extracted and stored sequentially in col1, col2,...
;
;    SILENT  (KEY) (I) (0) (S)
;       Optional keyword for turning off printout.
;
;    EFLD    (KEY) (I) (01) (I)
;       Optional keyword for specifying the particular columns (i.e., fields)
;       to be extracted. Up to 20 columns can be selected with the results
;       written to parameters COLn. If not specified, the first n columns
;       will be extracted.  EFLD will NOT be sorted.
;
;    NPTS (KEY) (IO) (01) (I)
;       Optional input/output keyword vector containing the number of points 
;       stored in extracted fields which contain variable length arrays. On
;       input, if NPTS is assigned a value, iue3drd will pad variable length 
;       arrays with 0's up to the size specified by maxlen. On output,
;       NPTS will contain a vector with 1 element for each extracted field.
;       A value of 0 designates fields with fixed length entries and
;       non-zero values represent the actual number of points in the variable
;       length array. Note: NPTS must be assigned some value on input or
;       no padding will be performed. 
;       
;
;*EXAMPLES:
;
;    read first 5 fields from the first row of the binary table starting in 
;    the 7th record of 'file.fit' with the extension header stored in exhd:
;        iue3drd,'file.fit',1,6,exhd,p1,p2,p3,p4,p5
;
;    as above but read only 4th field of first row, and turn off printout:
;        iue3drd,'file.fit',1,6,exhd,p4,efld=4,/silent
;
;    as above but read the 4th, 7th, and 9th columns from the 3rd row
;    of table:
;        iue3drd,'file.fit',3,6,exhd,p4,p7,p9,efld=[4,7,9],/silent;
;
;    read 3 fields (some of which are variable length arrays) 
;    store actual number of points for var. len. arrays in npts:
;        npts = 1
;        iue3drd,'file.fit',3,6,exhd,p10,p11,p12,efld=[10,11,12],npts=npts
;
;        
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
;   FITSCON 
;
;*FILES USED:
;
;*SIDE EFFECTS:
;
;*RESTRICTIONS:
;
;*PROCEDURE:
;      
;    FITS extension header is read to determine format of 3D table.
;    Fields are read and stored in COLn parameters. If any fields are
;    of type 'P' (i.e., variable length array pointers) then the file is
;    reopened and the data heap following the binary table is read. Data
;    is converted to the current operating system using the routine 
;    fitscon.
;
;*INF_1:
;
;*NOTES:
;
;     - if efld values are not in increasing order they will NOT be sorted and 
;     the output parameters will refer to the original field numbers. For
;     example, if 
;        iue3drd,'file.fit',3,6,exhd,a,b,c,efld=[4,9,7]
;     then a = field 4, b = field 9 and c = field 7.
;     Note that earlier versions of this program sorted efld.
;
;     - If npts is defined, variable length array entries are padded with 
;     0s according to the value of maxlen as stored in the binary table header
;     for the TFORMn keyword. The output NPTS vector can be used to remove 
;     the 0s but they are initially added so IFITSRD can store variable 
;     length array entries from multiple rows in one parameter.
;
;     - If the NPTS keyword is not included in the procedure call, or it
;     is included but not previously defined, then the var. len. arrays 
;     are not padded with 0's. (This is actually preferable when only one 
;     row of data is needed.)
;
;	tested with IDL Version 2.2.0 (sunos sparc)   14 Nov 91
;       tested with IDL Version 2.2.0 (ultrix mipsel) 14 Nov 91
;       tested with IDL Version 2.2.0 (ultrix vax)    14 Nov 91
;       tested with IDL Version 3.1   (vax vms)       23 Nov 93
;
;*MODIFICATION HISTORY:
;
;	written by R Thompson 2/7/91
;	3-15-91	RWT replace ieeetov routine with fitstov
;	3-20-91	RWT make totbyte array longword
;	4-16-91	PJL modified for unix/sun; replaced fitstov with fitssun
;       7-25-91 RWT modified for mulit-cpu: replaced fitssun with fitscon.
;       8-05-91 RWT allow version numbers (for vms), trim blanks from
;                   keywords, and correct st and sdrec calculations
;      11/14/91 GRA changed !version.os to !version.arch for FITSCON
;                   tests; modified EXTFITS and IUE3DRD to write high
;                   dispersion orders to FITS 3d table rows padded to
;                   a multiple of 2880 bytes; tested on sun, vax, dec.
;       3/11/92 RWT remove last change which made rows a multiple of 2880
;                   bytes (wasn't compatible with binary table proposal).
;       7/06/92 RWT accept B,C,M,P data types, TFORM without r specified,
;                   and 12 instead of 9 columns.
;       7/29/92 RWT use /block in openr statements
;       8/02/92 RWT allow 0 element fields
;       8/14/92 RWT read TNULL's and test for NaN's using FLAGNAN
;       9/15/92 RWT allow up to 15 columns
;       9/18/92 RWT allow both BINTABLE & A3DTABLE extension names
;      10/19/92 RWT remove flagnan call (added to fitscon)
;       2/18/93 RWT fix bug when vectors are exactly 2880 bytes long
;       8/10/93 RWT add silent keyword
;       8/30/93 PJL add efld keyword
;       9/15/93 RWT correct error in calculating newtot parameter
;       9/22/93 RWT add support for multi-dimensional arrays
;      11/18/93 RWT add support for variable length arrays
;      11/29/93 RWT set default value for THEAP if keyword not specified
;      12/05/93 RWT correct total data heap record calculation (i.e., hrtot)
;       6/23/94 RWT allow up to 15 columns from anywhere in the table to be
;                   extracted and read only requested data from table & heap
;                   to reduce memory usage
;       1/11/95 RWT leave undefined parameters undefined (for ifitsrd)
;       1/17/95 RWT pad P type variables with 0's according to maxlen and store
;                   actual # of points in NPTS vector.
;     24 Jan 95 LLT Remove calls to stpar from FOR loop, replace strpos with
;                   byte comparision and remove from FOR loop.  These changes
;                   increase speed.  Use new version of STPAR.  Remove sorting
;                   of efld vector.
;     14 Feb 95 RWT pad P type variables only if NPTS is defined
;     21 Apr 95 RWT correct problem when var. length array ends on record
;                   boundary (was reading an extra record)
;     21 Aug 95 RWT allow up to 20 columns
;     07 May 99 RWT strip off blanks in extension name (shouldn't be
;                   necessary, but FUSE test files have blanks after BINTABLE)
;
;******************************************************************************
 pro iue3drd,filename,row,srec,b3dh,col1,col2,col3,col4,col5,col6,col7,col8, $
 col9,col10,col11,col12,col13,col14,col15,col16,col17,col18,col19,col20, $
 silent=silent,efld=efld,npts=npts
;
 npar = n_params(0)
 optpar = npar - 4         ; # of optional parameters
 if (npar eq 0) then begin
    print,'iue3drd,FILENAME,ROW,SREC,B3DH,COL1,col2,...col20,' 
    print,'  silent=silent,efld=efld,npts=npts'
    retall
 endif  ; npar
 parcheck,npar,indgen(20)+5,'IUE3DRD'
;
; process input parameters
;
 tmp = size(filename)
 tc = tmp(tmp(0)+1)                  ; type code for filename parameter
 if (tc eq 7) then begin             ; filename is a string
    decompose,filename,disk,acc,fname,ext,vers
    if ext eq '' then ext = '.fit'
    fil = disk+acc+fname+ext+vers
 end else lun = filename             ; if filename is a number assume lun
 if keyword_set(silent) then prnt=0 else prnt=1
 if keyword_set(npts) then pad = 1.0 else pad = 0.0
;
; use efld to store actual columns to be read
;
 stpar,b3dh,'tfields',temptf        ; number of fields per row
 if (keyword_set(efld)) then begin  ;efld specified
;     efld = efld(sort(efld))
     indf = where(efld le temptf,count)
     count = (count < 20) < temptf
     if (count lt 1) then begin
        print,' all requested field number(s) are outside allowable range'
        print,' fields requested =',efld,' # of fields found =',temptf
        print,' returning from IUE3DRD'
        retall
     endif
     efield=efld(indf)
;    efld = efld(ind(0:count-1))
 endif else begin                   ;efld not specified
     count = (npar - 4) < temptf
     indf = indgen(count) 
     efld=indf+1
     efield=efld
 endelse
 minefield=min(efield)
 maxefield=max(efield)
;
; number of fields in table (up to 20 or efld value)
;
 nfields = temptf < maxefield
 r = lonarr(nfields)                          ; repeat count
 ppp=intarr(nfields)                          ; first letter of tform P or not
 maxlen = r                                   ; max. len for var. len. arrays
 dt = bytarr(nfields)                         ; data type
 vt = dt                                      ; var. len. array data type
 td = dt                                      ; array dimensions (tdim)
; ts = fltarr(nfields)                         ; tscale scale factor
; t0 = ts                                      ; tzero scale factor
 totbyte = long(r)                            ; total bytes per field
; type = strarr(nfields)                       ; field description
; unit = type                                  ; field units
; null = type                                  ; null values
 bc = lonarr(nfields)                         ; bytes per data point
 vbc = bc                                     ; same as bc for var. len. arr.
 heap = 0                                     ; number of var. len. arrays
 neloff = lonarr(2,nfields)                   ; heap pointers
;
; Read input extension header to determine table format
;
 stpar,b3dh,'xtension',tempxt
 tempxt = strtrim(tempxt,2)     ; added for FUSE
 if (tempxt ne 'BINTABLE') and (tempxt ne 'A3DTABLE') then begin
    print,' standard binary 3D table header not recognized'
    retall
 endif  ; tempxt
 stpar,b3dh,'naxis1',nbytes                   ; bytes per row
 nbytes = long(nbytes)
 stpar,b3dh,'naxis2',nrows                    ; number of rows
 tabb = nbytes * nrows                        ; total bytes in table (exc heap)
;
; read required tform keywords (up to last requested field)
;
 stpar,b3dh,'tform*',form,ntform
 form = strcompress(form,/remove_all)  ; remove blanks
;
; read optional keywords
;
 stpar,b3dh,'tdim*',td,ntd,nfields=ntform,typecode=7     ; look for tdim keyword
 if ntd gt 0 then begin
    temptyp = strcompress(td,/remove_all) ; remove blanks
    for i=0,n_elements(td)-1 do $
       td(i) = strmid(temptyp(i),1,strlen(temptyp(i))-2) ; remove ( & )
;    td = strmid(temptyp,1,strlen(temptyp)-2)    ; remove ( and )
 endif ;ntd gt 0
 stpar,b3dh,'ttype*',type,ntemp,nfields=ntform,typecode=7
 stpar,b3dh,'tunit*',unit,nfields=ntform,typecode=7
 stpar,b3dh,'tscal*',ts,nfields=ntform
 stpar,b3dh,'tzero*',t0,nfields=ntform
 stpar,b3dh,'tnull*',null,nfields=ntform,typecode=7
;
; evaluate TFORM keyword for each field
;
  for n=0,nfields-1 do begin
    bform = byte(form(n))                       ; convert to byte array
    ind = where( (bform lt 48) or (bform gt 57))  ; find first letter
    pos = ind(0) > 0                            ; 1st char. position
    dt(n) = bform(pos)                          ; data type
    r(n) = (pos eq 0)+(pos gt 0)*long( strmid(form(n),0,pos) )
    ppp(n)=dt(n) eq 80b                         ; is the first letter P?
    if ppp(n) then begin                        ; variable length arrays
       vt(n) = bform(pos+1)                     ; data type
       maxlen(n) = fix(strmid(form(n),pos+3,strlen(form(n))-4-pos)) ; max. len.
    endif ;ppp(n)
 endfor ;n=0,nfields-1

 bbb=fix(dt eq 88b)                                    ; does first letter = X? 
 r = r*(bbb eq 0)+bbb*(long(r/8.0) + ((r mod 8) ne 0))
 bc=(dt eq 76b)+(dt eq 65b)+(dt eq 66b)+(dt eq 88b)+$  ;first letter is L,A,B,X
    (dt eq 73b)*2 + $                                  ;first letter is I
    ((dt eq 74b)+(dt eq 69b))*4 + $                    ;first letter is J,E
    ((dt eq 67b)+(dt eq 80b)+(dt eq 68b))*8 + $        ;first letter is C,P,D
    (dt eq 77b)*16                                     ;first letter is M
 bc=long(bc)
;
; check for variable length array information
; Calculate max. length for column (bytes/element) for variable length items
;
 heap = fix(total(ppp))
 vbc=ppp*$                                             ; 0 if first letter not P
     ((vt eq 76b)+(vt eq 65b)+(vt eq 66b)+(vt eq 88b)+$ + $  ; L,A,B,X
      (vt eq 73b)*2 +$                                       ; I
      ((vt eq 74b)+(vt eq 69b))*4 + $                        ; J,E
      ((vt eq 67b)+(vt eq 80b)+(vt eq 68b))*8 + $            ; C,P,D
      (vt eq 77b)*16)                                        ; M
 vbc=long(vbc)
 totbyte = r * long(bc)
 dt=string(transpose(dt))
 vt=string(transpose(vt))
;
; If variable length arrays are found, find additional keywords
;
 if (heap gt 0) then begin
     stpar,b3dh,'theap',hboffset,err          ; byte offset to heap
     if err eq -1 then hboffset = tabb        ; default offset
     sbh = hboffset mod 2880                  ; starting byte
     stpar,b3dh,'pcount',hbtot                ; total bytes gap+heap
     gap = hboffset - tabb                    ; bytes in gap
     temp = sbh + hbtot - gap                 ; tot. bytes left from rec bound
     eod = (tabb + hbtot) mod 2880            ; last heap byte in last record
     hrtot = long(temp/2880) + (eod ne 0)     ; total records in heap
     if (temp le 2880) then hrtot = 1         ; necessary?     
     frheap = srec + long(hboffset/2880)      ; first heap record number
 endif 
;
; Open disk file to read appropriate table records (if filename is a string)
;
 if (tc eq 7) then  openr,lun,fil,/get_lun,/block
 rec=assoc(lun,bytarr(2880))	              ; define record type
;
; determine proper table record(s) to read
; total bytes to read from efld(0) to efld(count-1)
;
 startf = 0
; if (efld(0) gt 1) then startf = total(totbyte(0:efld(0)-2)) 
 if minefield gt 1 then startf=total(totbyte(0:minefield-2))
 totb = nbytes * (row-1) + startf             ; total bytes of data to skip
 st =  (totb mod 2880)                        ; starting index
; if (count le 1) then lastf = totbyte(efld(0)-1) else $
 if count le 1 then lastf=totbyte(minefield-1) else $
    lastf=total(totbyte(minefield-1:maxefield-1))
; lastf = total(totbyte(efld(0)-1:efld(count-1)-1))    ; total bytes to read
 temp = st + lastf                            ; bytes from last record boundary
                                              ;  to end of requested fields
 eod =  temp mod 2880                         ; ending index
 nrec = (long(temp/2880) + (eod ne 0))        ; # of data records to read 
 if (temp le 2880) then nrec = 1              ; read one record (in this case)
 sdrec = long(totb/2880)                      ; # of data records to skip
 startr = srec + sdrec                        ; starting data record
 buf = bytarr(lastf)                          ; buffer for reading row of table
 diff = 0
;
;
; loop through FITS file records till one row is read
;
 nrm1 = nrec - 1
 for i=0,nrm1 do begin
    x = rec(startr + i)                       ; reading data records
    init = st * (i eq 0)                      ; offset to 1st byte
    fin = (eod + 2880*(eod eq 0)) * (i eq nrm1)  +  2880 * (i ne (nrm1))
                                              ; offset to last byte
    buf(diff) = x(init:fin-1)
    diff = diff + (fin-init)
 endfor  ; i
 if (tc eq 7) then free_lun,lun    ; free lun if filename was a string
;
; separate fields & convert data to appropriate format 
;
    nn = 0
    for n=0,count-1 do begin
      fno = efield(n)-1
;      if (n eq 0) then st = 0 else st = total(totbyte(efld(0)-1:fno-1))
      if efield(n) eq minefield then st=0  $
          else st=total(totbyte(minefield-1:fno-1))
      if (prnt) then print,'Column',byte(efield(n)),':   field = ',type(fno),  $
          ' data type = ',dt(fno),'  npoints = ',strtrim(string(r(fno)),2)
      if (totbyte(fno) gt 0) then begin
        newtot = (st + totbyte(fno) -1) < (nbytes - 1)
        tmp = buf(st:newtot) 
        if (td(fno) ne '') then fitscon,tmp,dt(fno),gencol,tdim=td(fno) $
           else fitscon,tmp,dt(fno),gencol
        if (dt(fno) eq 'P') then begin
            neloff(0,fno) = gencol   ; define w/o scaling
        endif
        if (ts(fno) ne 0.0) then gencol = gencol * ts(fno) + t0(fno)
;
            case indf(n) of
              0:  col1 = gencol
              1:  col2 = gencol
              2:  col3 = gencol
              3:  col4 = gencol
              4:  col5 = gencol
              5:  col6 = gencol
              6:  col7 = gencol
              7:  col8 = gencol
              8:  col9 = gencol
              9:  col10 = gencol
              10:  col11 = gencol
              11:  col12 = gencol
              12:  col13 = gencol
              13:  col14 = gencol
              14:  col15 = gencol
              15:  col16 = gencol
              16:  col17 = gencol
              17:  col18 = gencol
              18:  col19 = gencol
              19:  col20 = gencol
              else:  print,'This should not happen.'
           endcase  ; n
       endif  ; totbyte(fno) gt 0
;
       if (prnt) then begin
         if (ts(fno) ne 0.0) then $
            print,'  scale factors used: tscale=',ts(fno),' tzero=',t0(fno)
         if (null(fno) ne 0) then $
            print,'  null flag defined, value =',null(fno)
       endif  ; prnt
    endfor  ; n
    npts = neloff(0,*)
    npts = rotate(npts,3)   ; create a row vector
    npts = npts(efield-1)       ; save only requested fields
;
;
; If any variable length array fields exist, reopen file and read additional
; records in heap. Assume data heap follows binary table. Data does not
; necessarily start on record boudary or immediately follow table data.
; (i.e., a gap exists if THEAP does not equal NAXIS1*NAXIS2).
;
;
 if (heap gt 0) then begin
    if (prnt) then begin
       print,' '
       print,' Reading variable length array data'
       print,' '
    endif
;
; determine proper range of heap data to read
;
    ind = where(neloff(0,*) ne 0) ; requested fields using heap
    a = neloff(1,*)        ; vector of offsets
    startf = min(a(ind))   ; first data point to read (from start of heap)
    b=vbc*neloff(0,*)+neloff(1,*) ; bytes/point * # of data points + offset
    lastf = max(b)         ; last data point to read from heap
;    
; store desired data from heap in buf
; redefine parameters to reduce portion of heap to be read
;
    if (tc eq 7) then openr,lun,fil,/get_lun,/block
    rec=assoc(lun,bytarr(2880))	              
    diff = 0
    buf = bytarr(lastf - startf)
    frheap = srec + long((hboffset+startf)/2880) ; new first heap record
    lrheap = srec + long((hboffset+lastf)/2880)  ; new last heap record
    sbh = (hboffset + startf) mod 2880           ; new start offset
    eod = (hboffset + lastf) mod 2880            ; new last offset
    nrm1 = lrheap - frheap - (eod eq 0)          ; new records to read
;   nrm1 = lrheap - frheap                       ; new records to read
;   nrm1 = hrtot - 1
;   buf = bytarr(hbtot-gap)
;
; read heap
;
    for i=0,nrm1 do begin                        ; read heap records
       x = rec(frheap + i)                       ; ith heap record
       init = sbh * (i eq 0)                      ; offset to 1st byte
       fin = (eod + 2880*(eod eq 0)) * (i eq nrm1)  +  2880 * (i ne (nrm1))
                                                 ; offset to last byte
       buf(diff) = x(init:fin-1)                 ; store in buf
       diff = diff + (fin-init)
    endfor  ; i
    if (tc eq 7) then free_lun,lun
;
; separate fields and convert data to appropriate format (as above)
;
     nn = 0
     ind = where(maxlen,count)   ; find fields with var. len. arrays
     for i=0,count-1 do begin
         n = ind(i)
         fno = where((n+1) eq efld) 
         if (fno(0) ge 0) then begin
            nn = nn + 1
            if (prnt) then begin
               nvla = neloff(0,n)
               print,'Column',byte(n+1),':   field = ',type(n),   $
               ' data type = ',vt(n),'  npoints = ',strtrim(string(nvla),2)
            endif
            if (nn eq 1) then ind0 = n
            st = neloff(1,n) - neloff(1,ind0)
            if (neloff(0,n) gt 0) then begin
               newtot = ( st - 1 + neloff(0,n)*vbc(n) )
               buflen = maxlen(n) * vbc(n) * (pad) + (newtot-st+1) * (not pad)
               tmp = bytarr(buflen)      ; if (pad) then pad tmp to max length 
               tmp(0) = buf(st:newtot) 
               if (td(n) ne '') then fitscon,tmp,vt(n),gencol,tdim=td(n) $
                  else fitscon,tmp,vt(n),gencol
               if (ts(n) ne 0.0) then gencol = gencol * ts(n) + t0(n)
;
               case fno(0) of
                   0:  col1 = gencol
                   1:  col2 = gencol
                   2:  col3 = gencol
                   3:  col4 = gencol
                   4:  col5 = gencol
                   5:  col6 = gencol
                   6:  col7 = gencol
                   7:  col8 = gencol
                   8:  col9 = gencol
                   9:  col10 = gencol
                  10:  col11 = gencol
                  11:  col12 = gencol
                  12:  col13 = gencol
                  13:  col14 = gencol
                  14:  col15 = gencol
                  15:  col16 = gencol
                  16:  col17 = gencol
                  17:  col18 = gencol
                  18:  col19 = gencol
                  19:  col20 = gencol
                  else:  print,'This should not happen.'
              endcase  ; n
            endif  ; neloff(0,n) gt 0
;
          if (prnt) then begin
            if (ts(n) ne 0.0) then $
              print,'  scale factors used: tscale=',ts(n),' tzero=',t0(n)
            if (null(n) ne 0) then $
              print,'  null flag defined, value =',null(n)
          endif ; prnt
;
       endif    ; fno ge 0 (i.e., n eq efld)
    endfor      ; i
 endif          ; heap

;
 return  
 end  ; iue3drd
