;******************************************************************************
;+
;*NAME:
;
;	READMX
;
;*CLASS:
;
;*CATEGORY:
;
;	NEWSIPS
;
;*PURPOSE:
;
;	Read IUE merged extracted spectrum fits files and return the
;	main FITS header (including vicar label) wavelengths and absolute 
;       flux. NEWSIPS flags, sigmas/noise vector,  background flux, net flux,
;       and ripple-corrected net fluxes may also be optionally returned.
;
;*CALLING SEQUENCE:
;
;	READMX,FILENAME,MAIN_HEADER,WAVE,ABS_FLUX,FLAGS,SIGMA,bkgrd,net,  $
;	   ripple,reqaper=reqaper,noninter=noninter,orange=orange, $
;          wrange=wrange,uncalib=uncalib,noisecal=noisecal
;
;*PARAMETERS:
;
;	FILENAME  (REQ) (I) (0) (S)
;		The filename - including extension.
;
;	MAIN_HEADER	(REQ) (O) (1) (S)
;		The primary fits keywords, the vicar label, and the history
;		portion of the fits header.  If there is an error, it will be
;		set equal to strarr(1).
;
;	WAVE	(REQ) (O) (1) (D)
;		The wavelength vector.  If there is an error, it will be set to
;		dblarr(1).
;
;	ABS_FLUX	(REQ) (O) (1) (R)
;		The absolute flux vector.  If there is an error, it will be
;		set to fltarr(1).
;
;	FLAGS	(opt) (O) (1) (I)
;		The NEWSIPS error flags vector.  If there is an error, it will
;		be set to intarr(1).
;
;	SIGMA	(opt) (O) (1) (R)
;               For low dispersion files the NEWSIPS sigma vector is returned. 
;               For high dispersion, the (uncalibrated) "noise" vector
;               is returned. (See NOISECAL keyword to calibrate noise vector.)
;               If an error occurs, SIGMA will be set to fltarr(1). 
;
;	BKGRD	(OPT) (O) (1) (R)
;		The background flux vector.  If there is an error, it will be
;		set to fltarr(1).
;
;	NET	(OPT) (O) (1) (R)
;		The net flux vector.  If there is an error, it will be set to
;		fltarr(1).
;
;       ripple  (opt) (o) (1) (f)
;               Ripple corrected data for high dispersion.  If not 
;               present (e.g., low dispersion), it is set to fltarr(1).
;
;	REQAPER	(KEY) (I) (0) (S)
;		Keyword for handling double aperture images.  May be equal
;		to 'LARGE' or 'SMALL'.  Not needed if the image is single
;		aperture or high dispersion. If it is not given for a double 
;               aperture image and the NONINTER keyword is not set, the user 
;               will be prompted. If it is not given for a double aperture 
;               image and the NONINTER keyword is set, the LARGE aperture 
;               will be assumed.  Note that this keyword is intended for use 
;               with low dispersion data.
;
;	NONINTER  (KEY) (I) (0) (I)
;		When set the user will not be prompted for information. For a
;		single aperture image, if the incorrect aperture is set with
;		the REQAPER keyword, the other aperture's data will not be
;		obtained (default 'n').  For a double aperture image, if the
;		REQAPER is not set, then the large aperture data is obtained.
;
;       wrange  (key) (i) (1) (f)
;               Wavelength range to be returned.  For high dispersion this is
;               used to determine which orders to extract---at present, each
;               order's wavelength range is determined and compared to WRANGE
;               and flagged for extraction as necessary.  This is slower than
;               using the ORANGE keyword.  For both dispersions, the data are
;               trimmed to match WRANGE.
;
;       orange   (key) (i) (1) (i)
;               Order range to extract.  If the MXHI file does not have all
;               of the expected orders (based on the number of orders present
;               and the values of the first and last order---e.g., if REDUCEMX
;               was used to discard orders), then it is slightly slower to
;               determine which orders to extract.
;               Note if both orange and wrange are specified, readmx uses
;               orange to determine which rows of the binary table to read
;               and then uses wrange to truncate the resulting vectors.
;               Generally, the keywords should not be specified together.
;
;       uncalib (key) (i) (1) (i)
;               By default, READMX removes points which are not calibrated
;               (i.e., points are removed where the nu flag value = -2).
;               Specifying /uncalib causes all extracted points to be output
;               (although splicing is still performed when multiple orders
;               are requested).
;
;      noisecal (key) (i) (1) (i)
;               The high dispersion noise vector stored in the MXHI file is 
;               uncalibrated. Specifying this keyword will cause the noise
;               vector to be scaled by the ratio of the absolute flux to 
;               the net flux; giving it the same units as the absolute
;               flux vector.
;
;*EXAMPLES:
;
;	For a single aperture low dispersion image:
;
;	   readmx,'swp32525.mxlo',main,wave,flux,flags,sigma,bkgrd,net
;
;       To read the region of a high dispersion spectrum containing Mg II:
;
;          readmx,'lwp12345.mxhi',h,w,f,q,s,wrange=[2790,2810]
;  
;       To output all extracted points for high dispersion order 100 and
;          calibrate the noise vector. Display the results using the 
;          calibrated noise vector as an error bar:
;
;          readmx,'lwp12345.mxhi',h,w,f,q,sc,orange=100,/uncalib,/noisecal
;          plot,w,f
;          errbar,w,f,w*0,sc,psy=1
;       
;*SYSTEM VARIABLES USED:
;
;	none
;
;*INTERACTIVE INPUT:
;
;	If the image is a single aperture image, and the keyword REQAPER is 
;	set but does not equal the aperture present, then the user is asked
;       about extracting the data present (low dispersion).
;
;	If the image is a low dispersion double aperture image, and the 
;	keywords REQAPER and  NONINTER are not set, the user will be prompted
;	with the options:
;
;			1 - large
;			2 - small
;			3 - exit
;
;       For high dispersion, if the NONINTER, WRANGE, and ORANGE keywords are 
;       not set, the user will be asked to enter an order range to extract.
;
;*SUBROUTINES CALLED:
;
;	PARCHECK
;	IUEFHRD
;	STPAR
;	CHKFITS
;	IUEDAF
;       within
;       new_splice
;
;*FILES USED:
;
;	filename given in calling sequence
;
;*SIDE EFFECTS:
;
;*RESTRICTIONS:
;
;      The input file must be a binary 3-D table FITS file with the following
;      fields (actual field location is unimportant):
;
;      ORDER  (high dispersion, required)
;      APERTURE (low dispersion, required)
;      NPOINTS (both dispersions, required)
;      WAVELENGTH (both dispersions, required)
;      DELTAW (both dispersions, required)
;      NET (both dispersions, optional)
;      BACKGROUND (both dispersions, optional)
;      NOISE (high dispersion, optional, will be returned in SIGMA vector)
;      SIGMA (low dispersion, optional)
;      QUALITY (both dispersions, optional)
;      FLUX (low dispersion, required)
;      ABS_CAL (high dispersion, required)
;      RIPPLE (high dispersion, optional)
;      STARTPIX (high dispersion, required)
;
;      For SWP high dispersion data, no absolutely calibrated fluxes above
;      about 1980 A are available, and that vector contains zeroes beyond
;      that point.  If the minimum wavelength is below 1980 A, and the
;      maximum is above that value, all extracted vectors will be truncated
;      to 1980 A.  To examine those vectors, the user should ensure that the
;      minimum wavelength is above 1980 A.
;
;*NOTES:
;
;      - The NOISECAL keyword has no effect on low dispersion data.
;      - Users are cautioned that the calibrated noise vector may be an
;        inaccurate estimate of the error for high dispersion fluxes.
;
;*PROCEDURE:
;
;       IUEFHRD is used to obtain the primary and extension fits keywords.
;	IUEDAF is used to add an IUEDAC section (if necessary) and add/update
;	the fits keyword EXTDATE.  If the fits keyword FILENAME is not in the
;	main fits header, it is added.  The fits keyword FILENAME is check 
;       to make sure 'MX' are the first two characters of the
;	extension.  If not, then the file is assumed not to be for a merged
;	extracted spectrum image.  The parameters of the file - obtained via
;	IUEFHRD - are checked.  If they do not match the expected values for a
;	merged extracted spectrum image file, the procedure returns.
;
;       This program looks at the TTYPE keywords in the extension header and
;       finds the location of each item it needs.  If the APERTURE/ORDER,
;       WAVELENGTH, DELTAW, and NPOINTS keywords are missing, the procedure
;       returns.  
;
;       Low Dispersion:
;	For a single aperture image, if the keyword REQAPER is set and does
;	not equal the aperture present, the user is asked about obtaining the
;	data present.  
;	For a double aperture image, the keyword REQAPER needs to be set.  If
;	it is not set and the keyword NONINTER is not set, the user is prompted
;	with the options of '1 - large', '2 - small', and '3 - exit'.  When
;	'3 - exit' is selected the procedure returns.  When '1 - large' or
;	'2 - small' is selected, IF REQAPER is not set and NONINTER is set,
;	then the large aperture is used.  IUE3DRD is used to obtain the
;	aperture and determine which row the requested data is in.   
;
;       High Dispersion:
;       The range of orders, and the number of rows present, are determined.
;       If the user set the WRANGE keyword, the wavelength range for every 
;       order in the table is determined; those that contain data within the
;       limits are flagged for extraction.  This is slow.  If the user set the
;       ORANGE keyword, the program can assume it knows the location of the
;       desired orders IF the number of rows present matches the difference
;       between the first and last orders.  Otherwise, it will have to look
;       at each row and flag those orders that are in the desired range.
;       The orders are extracted and the data trimmed to their wavelength
;       ranges using the STARTPIX keyword.  They are spliced at the 
;       point where the difference in the ripple correction used in adjacent
;       orders is a minimum. If the user did not set either the WRANGE or 
;       the ORANGE keywords, the program prompts for an order range.
;
;	ADDPAR is used to add/update the EXTFILE fits keyword.  The uncalibrated
;	data is removed from the wavelength, absolute flux, nu flags, sigma,
;	and net flux vectors unless /uncalib is specified.  If WRANGE is 
;       set the data are trimmed.
;
;*I_HELP  nn:
;
;*MODIFICATION HISTORY:
;
;	 2 Dec 92  PJL  wrote
;        4 Dec 92  PJL  split off merged extracted spectrum image specific
;			code from IUEREAD
;	 7 Dec 92  PJL  changed the wavelength arrays from single to double
;			precision
;	15 Jan 93  PJL  made vicar label optional; changed parameter order
;	18 Jan 93  PJL  added check for reqaper when only one aperture's data
;			available
;	 4 Mar 93  PJL  added CHKFITS and telescope check
;	17 Mar 93  PJL  removed both option; if returning because of an error,
;			then main = strarr(1) and flux = fltarr(1)
;	 3 Jun 93  PJL  add IUEDAF, CURESTR, findfile, and fits keyword EXTFILE;
;			remove vicar as optional parameter
;	 4 Jun 93  PJL  trim uncalibrated data points; update prolog
;	14 Jun 93  PJL  replaced IFITSREAD with DECIPHMX
;	15 Jun 93  PJL  do not trim vectors that do not exist; FILENAME fits
;			keyword should be in IUEDAC section
;	16 Jun 93  PJL  added message about trimming uncalibrated data
;	25 Jun 93  PJL  corrected calculation of points removed for message
;       16 Sep 94  RWT  add blanks to make EXTAPER keyword value 8 characters
;       14 Dec 94  RWT  optimize to reduce execution time by calling iue3drd
;                       directly when possible and calling FDMX instead of 
;                       DECIPHMX
;       23 Feb 95  RWT  fix call to FDMX to extract proper row for dbl aper data
;        2 Oct 95  LLT  High dispersion support, consolidated low dispersion 
;                       code and removed FDMX, flags and sigma optional, orange,
;                       wrange keywords, ripple parameter.
;        1 Oct 96  RWT  correct nu flag vector length when orange and wrange
;                       both specified. Also use dindgen command with wave.
;       18 Nov 96  jrc  added comment in prolog about SWP high dispersion
;                       abs. fluxes and vectors above 1980 A.
;       24 Dec 96  RWT  add /uncalib keyword option
;       24 Feb 97  RWT  make splice point where difference in ripple correction
;                       is a minimum
;       28 Jul 97  RWT  add /noisecal keyword option
;       05 Aug 97  RWT  show flags & sigma vectors as optional in print command
;       18 Sep 97  RWT  if wrange is outside extracted wavelength region,
;                       set all vectors = 0 and return
;       05 Nov 97  RWT  add header comments about order and wavelength range 
;       17 Nov 97  RWT  fix bug when only 1 order is input 
;       02 Feb 98  RWT  return all orders when no range is specified and
;                       noninter is set.
;-
;******************************************************************************
 pro readmx,filename,main,wave,flux,flags,sigma,bkgrd,net,ripple, $
        reqaper=reqaper,noninter=noninter,wrange=wrange,orange=orange, $
        uncalib=uncalib,noisecal=noisecal
;
 npar = n_params(0)
 if (npar eq 0) then begin
    print,'READMX,FILENAME,MAIN_HEADER,WAVE,ABS_FLUX,flags,sigma,bkgrd,net,  $'
    print,'  ripple,reqaper=reqaper,/noninter,wrange=wrange,orange=orange, $'
    print,'  uncalib=uncalib,noisecal=noisecal'
    retall
 endif  ; npar eq 0
 parcheck,npar,[4,5,6,7,8,9],'READMX'
;
;  set defaults
;
 main = strarr(1)
 wave = dblarr(1)
 flux = fltarr(1)
 flags = intarr(1)
 sigma = flux
 bkgrd = flux
 net = flux
 ripple=flux
;
;  check that the file in question is an existing fits file
;
 chkfits,filename,checked,/silent
 case checked of
 0: begin
       print,filename + 'is not a fits file.'
       print,'ACTION: Returning.'
       return
    end
 1: 
 2: begin
      print,'Unable to find ' + filename
      print,'ACTION:  Returning
      retall
    end
 3: begin
      print,'File ',filename,' found but not readable'
      print,'ACTION:  Returning
      retall
    end
 endcase   ; checked

;
;  retrieve main fits header, vicar label, and fits history portion - and the
;  extension header
;
 iuefhrd,filename,params,main,vicar,extn,/silent,/full
 srec = total(params(0:2))                   ; starting data record

;
;  determine if there is an IUEDAC section
;
 iuedaf,main,daflag,/date
 if (not(daflag)) then begin
    print,' '
    print,'There is no IUEDAC section in the fits header.'
    print,'ACTION:  Continuing.'
 endif  ; daflag

;
;  check telescope keyword
;
 stpar,main,'telescop',tel,err
 if (strtrim(tel,2) ne 'IUE') then begin
    print,'Error in extracting telescope name from fits header.'
    print,'ACTION:  Returning'
    return
 endif  ; strtrim(tel,2) ne 'IUE'
;
;  extract filename from fits header
;  if the FILENAME keyword is not already in the main header, find it
;  in the extension header and add it to the main header
;
 stpar,main,'filename',name,merr
 if (merr ne 0) then begin
    stpar,extn,'filename',name,err
    if (err ne 0) then begin
       print,'Error in extracting filename from fits header.'
       print,'ACTION:  Returning'
       return
    endif  ; err ne 0
    addpar,main,'filename',name,' Filename(camera)(number).MX(disp)','IUEDAC'
 endif  ; merr ne 0
;
;  extract file type from filename
;  check that type indicates a merged extracted spectral file
;
 pdpos = strpos(name,'.')
 type = strlowcase(strmid(name,pdpos+1,4))
 if (strmid(type,0,2) ne 'mx') then begin
    print,filename + ' is not a merged extracted spectral file.'
    print,'ACTION: Returning.'
    return
 endif  ; strmid(type,0,2) ne 'mx'

 stpar,extn,'naxis2',naxis2,err
 if err ne 0 then begin
    print,'Unknown number of apertures/orders in FITS table.  Returning.'
    return
 endif ;naxis2 error

 stpar,main,'disptype',disp,derr
 if derr ne 0 then stpar,main,'dispersn',disp,derr
 if derr ne 0 then begin
    print,'Unknown dispersion.  Returning.'
    return
 endif else disp=strtrim(disp,2)

;Extract all the TTYPE keywords; find out if the fields we need exist.

 stpar,extn,'ttype*',ttype,err

 if disp eq 'HIGH' then $
    field=['ORDER   ','NPOINTS ','WAVELENGTH','DELTAW  ','ABS_CAL ',$
            'BACKGROUND','NOISE   ','QUALITY ','NET     ','RIPPLE  ',$
            'STARTPIX'] $
 else field=['APERTURE','NPOINTS ','WAVELENGTH','DELTAW  ','FLUX    ',$
              'BACKGROUND','SIGMA   ','QUALITY ','NET     ']

;nf will contain the location of the desired field in the file.  This is so
;that files written with WRITEMX (which may not include all fields) can be
;read.  Fields not found are assigned position -1.  In these cases, the first
;field (should be a scalar field) will be read in their places.  This will also
;be done for optional parameters not specified in the calling sequence, EXCEPT
;for the quality flags, which are needed to trim the uncalibrated data.

;NOTE:  This program expects certain items in the FIELD vector to be in the
;same position for high and low dispersion (e.g., QUALITY).  Also, subscripts
;referring to various items occur when vectors are sorted etc.  Thus, if new
;items are to be added, it would be easiest to add them to the end of FIELD.
;Be careful to keep the correspondence between FIELD and the calling sequence
;of IUE3DRD.

 n=n_elements(field)-1
 nf=replicate(1,n+1)
 if npar lt 8 then nf(where(field eq 'NET     '))=-1     ;optional parameters
 if npar lt 7 then nf(where(field eq 'BACKGROUND'))=-1
 if disp eq 'HIGH' then begin
    if npar lt 9 then nf(where(field eq 'RIPPLE  '))=-1
    if npar lt 6 then nf(where(field eq 'NOISE   '))=-1 $
    else begin       ; need net flux to calibrate noise vector
       if (keyword_set(noisecal)) then nf(where(field eq 'NET     '))=1
    endelse 
 endif else if npar lt 6 then nf(where(field eq 'SIGMA   '))=-1 

 for i=0,n do if nf(i) gt 0 then begin
     temp=where(ttype eq field(i))
     if temp(0) eq -1 then print,field(i)+' not found in FITS table.'
     nf(i)=temp(0)+nf(i)
;print,field(i),nf(i),' ',ttype(temp(0))
 endif ;
 quality=nf(where(field eq 'QUALITY '))
 if total(nf(0:4) gt 0) lt 5 then begin
    print,'WAVELENGTH, NPOINTS, DELTAW, FLUX/ABS_CAL, and APERTURE/ORDER'
    print,'must ALL be present for this program to work.  Returning.'
    return
 endif ;

;Read the first field (aperture for low disp, order for high).  Also (if more
;than one row is present) read the first field in the last row.  Here it is
;assumed that the data are arranged so that the orders are monotonic (either
;increasing or decreasing).

 iue3drd,filename,1,srec,extn,ordap,efld=nf(0),/silent
 if naxis2 gt 1 then begin
    iue3drd,filename,naxis2,srec,extn,tempao,efld=nf(0),/silent
    ordap=[ordap,tempao]
 endif ;naxis2 gt 1
 row=-1

 if disp eq 'HIGH' then begin
    maxord=fix(max(ordap)) 
    minord=fix(min(ordap))

;Since we don't have K values for determining the blaze wavelengths for each
;order, currently we're reading all the wavelengths fromall the rows, 
;calculating the min and max per order, and determining what rows to extract
;using that information.  It is expected that this will be time consuming and
;should be removed once a better way presents itself.

    if keyword_set(wrange) and not keyword_set(orange) then begin
       if total(nf(1:3) lt 0) gt 0 then begin
          print,'Unable to determine wavelength range without NPOINTS,'
          print,'WAVELENGTH, and DELTAW.  Returning.'
          return
       endif ;
       print,'Determining orders, try to be patient.....'
       w0=fltarr(naxis2)
       wf=fltarr(naxis2)
       i=1
       for i=1,naxis2 do begin
          iue3drd,filename,i,srec,extn,n,tempw,dw,efld=nf(1:3),/silent
          w0(i-1)=tempw          ;Min wave for order
          wf(i-1)=tempw+dw*n     ;Max wave for order
       endfor

;Find out what orders, if any, overlap with the desired range.

       within,wrange,w0,res0                            
       within,wrange,wf,resf                            
       temp=where((res0 eq 0) or (resf eq 0) or ((resf+res0) eq 0),ntemp)     
       if ntemp eq 0 then begin
          print,'No data within specified wavelength range.  Returning.'
          return
       endif else row=temp+1

    endif else begin

;If neither a wavelength range nor an order range was given, prompt for it
;(or return if the user set NONINTER).

       if not keyword_set(orange) then begin
          print,'There are '+strtrim(naxis2,2)+' orders available numbered'+$
                ' from '+ strtrim(minord,2)+' to '+strtrim(maxord,2)+'.'
          order=intarr(2)
          if not keyword_set(noninter) then begin
             stl = 'Enter desired min and max order number (or a single order'
             print,stl+' and 0),'
             read,' separated by a blank or comma: ',order
             temp=where(order ne 0,ntemp)
             if ntemp eq 1 then order=order(temp)
          endif else begin
             print,'No orders selected.  Read all.'
             order=[minord,maxord]
          endelse ;noninter
       endif else order=orange ;

;Test to see whether all possible orders in the range from minord to maxord
;are actually likely to exist (as they should in an original MXHI file).  If
;not, then each order will have to be read so that the desired rows can be
;located.

       if (maxord-minord+1) ne naxis2 then begin
          ord=intarr(naxis2)
          ord(0)=ordap(0) 
          ord(naxis2-1)=ordap(1)
          print,'Need to read orders.  Stand by....'
          for i=2,naxis2-1 do begin
              iue3drd,filename,i,srec,extn,temp,/silent
              ord(i-1)=temp
          endfor
          within,order,ord,res
          row=where(res eq 0,n)
          if n gt 0 then row=row+1 else begin
             print,'None of the desired orders are available.  Returning.'
             return
          endelse  ;
      endif else begin

;Test to see if the orders desired are in range.  Return if all are out of
;range.  Otherwise, discard those that are out of range and keep the rest.

          within,[maxord,minord],order,res
          if total(res ne 0) ne 0 then begin
             print,'Out of range order(s): ',order(where(res ne 0))
             print,'Available orders: '+strtrim(fix(minord),2)+' - '+$
                 strtrim(fix(maxord),2)
             if total(res ne 0) eq n_elements(res) then begin
                print,'Returning.'
                return
             endif ;
          endif ;
          order=order(where(res eq 0))
          row=indgen(max(order)-min(order)+1)+(maxord-max(order))+1
       endelse ;are all orders included in file or not
    endelse ; wrange or orange

;Read the first desired row, initializing the output vectors.

    nrows=n_elements(row)-1
    iue3drd,filename,row(0),srec,extn,ord,n,wtemp,dw,flux,bkgrd,sigma,$
         flags,net,ripple,spix,efld=abs(nf),/silent
    wave=wtemp+dindgen(n)*dw
    print,'Order: ',strtrim(fix(ord),2)+$
          '   '+strtrim(wtemp,2)+' - '+strtrim(wave(n-1),2)
    firstm = fix(ord)
    lastm = firstm
    indo=indgen(n)+spix-1
    flux=flux(indo)
    if nf(8) gt 0 then net=net(indo)                   ;Don't bother with
    if nf(5) gt 0 then bkgrd=bkgrd(indo)               ;vectors not extracted
    if nf(6) gt 0 then begin
       sigma=sigma(indo)
       if (keyword_set(noisecal)) then sigma = sigma * (flux/net)
    endif
    if nf(7) gt 0 then flags=flags(indo)
    if nf(9) gt 0 then ripple=ripple(indo)
    wtemp = wave                               ; save first wavelengths

;Read each successive desired row, comparing the wavelength ranges.  Orders
;will be spliced at the average of the overlap regions for the time being.

    for i=1,nrows do begin
      mlast = ord      ; save previous order number
      wlast = wtemp    ; save previous wavelengths
      iue3drd,filename,row(i),srec,extn,ord,n,wtemp,dw,ftemp,btemp,noise,$
         qtemp,ntemp,rtemp,spix,efld=abs(nf),/silent
      wtemp=wtemp+dindgen(n)*dw
      print,'Order: ',strtrim(fix(ord),2)+$
           '   '+strtrim(wtemp(0),2)+' - '+strtrim(wtemp(n-1),2)
      lastm = fix(ord)
      within,wave,wtemp,res
      temp=where(res eq 0,nres)
      if nres gt 0 then begin
         new_splice,main,mlast,wlast,ord,wtemp,wavg ; new splice point
;         owavg=total(wtemp(temp))/nres            ;Average wavelength (overlap)
;         print,'old wavg =',owavg,'new wavg =',wavg
         if min(wtemp) lt min(wave) then begin    ;Indices of kept points will
            indi=where(wtemp lt wavg)+spix-1      ;depend on whether the range
            indo=where(wave ge wavg)              ;of the current order is
         endif else begin                         ;shortward or longward of the
            indi=where(wtemp ge wavg)+spix-1      ;previous orders.  It is
            indo=where(wave lt wavg)              ;assumed that one does not
         endelse                                  ;entirely contain the other!
      endif else begin
         indi=indgen(n)+spix-1
         indo=indgen(n_elements(wave))
      endelse                                     ; ranges do not overlap
      wave=[wave(indo),wtemp(indi+1-spix)] 
      flux=[flux(indo),ftemp(indi)] 
      if nf(8) gt 0 then net=[net(indo),ntemp(indi)]     ;Don't bother with
      if nf(5) gt 0 then bkgrd=[bkgrd(indo),btemp(indi)] ;vectors not extracted
      if nf(6) gt 0 then begin
         if (keyword_set(noisecal)) then $
           sigma=[sigma(indo),noise(indi) * ftemp(indi)/ntemp(indi)] $
         else sigma=[sigma(indo),noise(indi)] 
      endif
      if nf(7) gt 0 then flags=[flags(indo),qtemp(indi)]
      if nf(9) gt 0 then ripple=[ripple(indo),rtemp(indi)]
    endfor ; i=1,nrows

;Sort the data so that wavelengths are monotonically increasing.

    sw=sort(wave) 
    wave=wave(sw) 
    flux=flux(sw)
    if nf(8) gt 0 then net=net(sw) else net=[0.]
    if nf(5) gt 0 then bkgrd=bkgrd(sw) else bkgrd=[0.]
    if nf(6) gt 0 then sigma=sigma(sw) else sigma=[0.]
    if nf(7) gt 0 then flags=flags(sw) else flags=[0]
    if nf(9) gt 0 then ripple=ripple(sw) else ripple=[0.]

 endif else begin   ;LOW dispersion

;If the user requested a specific aperture, figure out if that aperture
;is actually present. If so, note which row it is in.  Otherwise, assume
;(for single rows) that the row available is acceptable.

    if keyword_set(reqaper) then begin
       pos=strpos(ordap,strupcase(strtrim(reqaper,2)))
       row=where(pos eq 0)
       if (row(0) eq -1) then print,'Aperture "'+reqaper+'" not available.'
    endif else if naxis2 eq 1 then row=0

;If the row to extract has not yet been determined, and the user sets 
;NONINTER, assume LARGE is desired and note which row it's in. Otherwise,
;display the row(s) available and ask the user which they'd like.

    if (row(0) eq -1) then begin
       if n_elements(reqaper) eq 0 then ap='LARGE' else ap=reqaper
       if keyword_set(noninter) then row=where(ordap eq ap,ntemp) $
       else repeat begin
          for i=1,naxis2 do print,'        '+strtrim(i,2)+' - '+ordap(i-1)
          print,'        '+strtrim(i,2)+' - Exit'
          read,'Enter the number of the aperture to extract: ',tmp
          finished=(tmp ge 1) and (tmp le i)
          if not finished then print,'Invalid choice.'
          if (tmp eq i) then row=-1 else row=tmp-1
       endrep until finished ; interactive, double aperture
    endif ;row still not selected

;If the row has STILL not been selected, complain and return.

    if row(0) eq -1 then begin
       print,'Returning with no data extracted.'
       return
    endif ;row never selected

;Extract the desired data, making a note in the IUEDAC section of the header
;about which aperture was extracted.

    row=row(0)
    print,ordap(row)+' aperture data to be extracted.'
    stpar,main,'extaper',aper,err
    if string(aper) ne ordap(row) then addpar,main,'extaper',ordap(row),$
       ' aperture from which data is extracted','iuedac'
    iue3drd,filename,row+1,srec,extn,aper,n,wave,dw,flux,bkgrd,sigma,flags,$
       net,efld=abs(nf),/silent
    wave=wave+indgen(n)*dw
 endelse ;high or low dispersion
     
 decompose,filename,disk,path,name,extn,version
 addpar,main,'extfile',strupcase(strtrim(name + extn,2)),   $
    ' file from which data was extracted','iuedac'
;
;  trim the uncalibrated data points  (where nu flag equals -2)
;
 if (not keyword_set(uncalib)) then begin
   svar = size(flags)
   if (svar(0) eq 1) then keep = where((abs(flags) and 2) eq 0,count) $
   else begin
      print,''
      print,'Nu flags unavailable.'
      print,'Unable to trim the uncalibrated data points (where nu flags ' +  $
       'equal -2).'
      count = 0
   endelse  ; (svar(0) eq 1
   if (count ne 0) then begin
      oldsvar = svar(1)
      print,' '
      print,'Trimming uncalibrated data regions (nu flag = -2).'
      wave = wave(keep)
      flux = flux(keep)
      if nf(7) gt 0 then flags = flags(keep)
      if nf(6) gt 0 then sigma = sigma(keep)
      if nf(8) gt 0 then net = net(keep)
      if nf(5) gt 0 then bkgrd = bkgrd(keep)
      if (disp eq 'HIGH') then if nf(9) gt 0 then ripple=ripple(keep)
      print,strtrim(oldsvar-count,2) + ' points removed.    ' +    $
         strtrim(count,2) + ' data points remain.'
   endif  ; count ne 0
 endif
;
; check wrange
;
 if keyword_set(wrange) then begin
    within,wrange,wave,res
    temp=where(res eq 0,n)
    if (n gt 0) then begin
       flux=flux(temp) 
       wave=wave(temp)
       if nf(7) gt 0 then flags = flags(temp)
       if nf(6) gt 0 then sigma = sigma(temp)
       if nf(8) gt 0 then net = net(temp)
       if nf(5) gt 0 then bkgrd = bkgrd(temp)
       if (disp eq 'HIGH') then if nf(9) gt 0 then ripple=ripple(temp)
       print,strtrim(n,2)+' points lie within the wavelength range '+$
          strtrim(min(wrange),2)+' - '+strtrim(max(wrange),2)
    endif else begin
       print,' no data found within specified wavelength range. Returning'
       wave = 0
       flux = 0
       flags = 0
       sigma = 0
       net = 0
       bkgnd = 0
       ripple = 0
       return
    endelse
 endif ;wrange
 if (disp eq 'HIGH') then addpar,main,'COMMENT','Orders extracted: ' + $
         strtrim(string(firstm),2) +' to '+ $
         strtrim(string(lastm),2),' ','IUEDAC' 
 addpar,main,'COMMENT','Wavelengths originally extracted: '+ $
         strtrim(string(min(wave)),2) + ' to '+ $
         strtrim(string(max(wave)),2) + ' Angstroms',' ','IUEDAC' 
;
 return
 end  ; readmx
