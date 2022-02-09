;***************************************************************************
;+
;*NAME:
;	DATECONV
;
;*CLASS:
;       Data Conversion
;
;*CATEGORY:
;
;*PURPOSE:
;	Procedure to perform conversion of dates between any of 5 possible
;	formats:
;	format 1: real*8 scalar encoded as:
;		(year-1900)*1000 + day + hour/24. + min/24./60 + sec/24./60/60
;		where day is the day of year (1 to 366)
;	format 2: Vector encoded as:
;		date(0) = year (on output, year will be full year (e.g., 1987)
;                               but as input, can be specified as 1987 or 87)
;		date(1) = day of year (1 to 366)
;		date(2) = hour
;		date(3) = minute
;		date(4) = second
;	format 3: string (ascii text) encoded as one of the
;               following strings:
;		- DD-MON-YEAR HH:MM:SS.SS
;		   (eg.  14-JUL-1987 15:25:44.23)
;               - CCYY-MM-DD or CCYY-MM-DDThh:mm:ss
;                  (e.g., 1996-10-14 or 1996-10-14T10:14:36.123)
;                  the new FITS date convention
;               - DD/MM/YY (e.g., 22/07/93)
;                  the original FITS date convention
;               - DD/MM/YY HH:MM:SS.SS 
;                  (e.g., 19/10/98 10:30:20.50)
;                  only used as an input format
;	format 4: three element vector giving spacecraft time words
;               from ST telemetry packet.
;       format 5: integer number in form yymmdd or yyyymmdd (required for
;               dates after year 1999.
;
;*CALLING SEQUENCE:
;	DATECONV,DATE,TYPE,OUT
;
;*PARAMETERS:
;	DATE (REQ) (I) (0,1) (DSIL)
;            input date in one of the allowed formats.
;
;	TYPE (REQ) (I) (0) (S)
;            type of output format desired.  
;	     valid values:
;			'REAL'	- format 1
;			'VECTOR' - format 2
;			'STRING' - format 3
;                       'FITS'  - modified format 3 
;                       'NEWFITS' - modified format 3 
;                       'INTEGER' - format 5
;            TYPE can be abbreviated to the single character strings 'R',
;            'V', 'S', 'F', 'N', or 'I'.
;            Nobody wants to convert TO spacecraft time (I hope!)
;
;       OUT  (REQ) (O) (0,1) (BILFDS)
;            Converted date in specified format.
;
;*PROCEDURE:
;       Input date format is determined using SIZE information.
;       Year, day, hour, minute, and second are calculated and reformatted
;       as specified by TYPE. 
;
;*SYSTEM VARIABLES USED:
;
;*SUBROUTINES CALLED:
;       IUEGETTOK - Procedure version of gettok
;       PARCHECK
;
;*SIDE EFFECTS:
;
;*RESTRICTIONS:
;
;*PROCEDURE:
;
;*EXAMPLES:
;       To generate current date in a format suitable for a FITS keyword:
;          dateconv,!stime,'fits',out
;
;       To convert 1994 day 152 to yymmdd:
;          dateconv,94152.0,'I',out
;       out will equal 940601 (i.e., 1994, June 01)
;
;*NOTES:
;       Format 3 can be input in 1 of 4 string formats. Note the 1st
;          string format is that used by the IDL system variable !STIME.
;          The other 3 are FITS-related.
;       Format 5 actually produces a longword integer.
;       Assume all 2 digit year entries are 20th century dates.
;       Format 1 is not compatible with dates before 1900.
;
;*MODIFICATION HISTORY:
;	version  1  D. Lindler  July, 1987
;       adapted  for IDL version 2  J. Isensee  May, 1990
;       8/06/93  rwt convert function to a procedure, update prolog, use
;                IUEGETTOK instead of GETTOK, and add FITS option
;       8/03/94  rwt add format 5 option.
;      11/01/94  rwt fix bug when input date is in format 5
;      11/02/94  rwt allow FITS format strings as input, set internal year
;                parameter to be consistently 2 digits (i.e., 92 instead of 
;                1992),
;               and only use full year in output formats 'R' and 'S'.
;      02/10/95 tmw added check in the interger output block to allow for
;               leap years.
;      08/21/98 rwt make y2k compatible by making intrnal year
;               4 digits, and assume 2 digit input implies 20th century
;               dates.    
;      08/28/98 rwt add support for new FITS date convention.
;      10/14/98 rwt add 4th input string format
;       3/25/99 rwt keep 2 digits for hours in 'newsips' mode
;       4/01/99 rwt avoid 60 sec, 60 min, and 24 hr values
;-
;***************************************************************************
pro dateconv,date,type,out
;
npar = n_params()
if npar eq 0 then begin
   print,' DATECONV,DATE,TYPE,OUT'
   print," where type = 'REAL', 'VECTOR', 'STRING', 'FITS', 'NEWFITS', or 'INTEGER' "
   retall & end
parcheck,npar,[2,3],'DATECONV'
;
; data declaration
;
days = [0,31,28,31,30,31,30,31,31,30,31,30,31]
months = ['   ','JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG','SEP','OCT',$
	'NOV','DEC']
hour = 0     
minute = 0
sec = 0
;
; Determine type of input supplied
;
s = size(date) & ndim = s(0) & datatype = s(ndim+1)
if (ndim gt 0) then begin               ;vector?
	if ndim gt 1 then goto,notvalid
	if (s(1) ne 5) and (s(1) ne 3) then goto,notvalid
	if (s(1) eq 5) then form = 2 else form = 4
end else begin			;scalar input
        case datatype of
           0 : goto,notvalid    ; undefined
           2 : form = 5         ; integer
           3 : form = 5         ; integer
           7 : form = 3         ; string
           else: form = 1       ; all else assume numeric scalar
        endcase
end
;      -----------------------------------
;
;*** convert input to year,day,hour,minute,second
;
;      -----------------------------------
case form of

	1: begin					;real scalar
		idate = long(date)
                year = long(idate/1000)
		day = idate - year*1000
		fdate = date-idate
		fdate = fdate*24
		hour = fix(fdate)
		fdate = (fdate-hour)*60.0
		minute = fix(fdate)
		sec = float((fdate-minute)*60.0)
		year = year + 1900
	   end

	2: begin					;vector
		year = fix(date(0))
;                year = year - 1900 * (year gt 1900)
                year = year + 1900 * (year lt 100)
		day = fix(date(1))
		hour = fix(date(2))
		minute = fix(date(3))
		sec = float(date(4))
	   end

	3: begin					;string
		temp = date
                stest = strpos(temp,'/',1)
                if (stest gt 0) then begin      ; old FITS format input
                   sstest = strpos(temp,':',1)
                   if (sstest gt 0) then begin  ; old FITS with added hh:mm:ss.ss
                      ntemp = strmid(temp,9,11)
                      iuegettok,ntemp,':',hour
                      iuegettok,ntemp,':',minute
                      sec = float(strtrim(strmid(ntemp,0,5)))
                      temp = strmid(temp,0,9)
                   endif
                   iuegettok,temp,'/',day_of_month
                   iuegettok,temp,'/',mon
                   iuegettok,temp,'/',year
                   year = year + 1900 * (year lt 100) ; new line
                   if ((year mod 4 eq 0) and (year mod 100 ne 0)) or $
                      (year mod 400 eq 0) then days(2) = 29
                   day = fix(total(days(0:mon-1))+day_of_month)
                endif else begin               
                   stest = strpos(temp,'-',1)
                   if (stest le 2) then begin   ; IDL format input
                      iuegettok,temp,'-',day_of_month
                      iuegettok,temp,'-',month_name
                      iuegettok,temp,' ',year
;                      year = year - 1900 * (year gt 1900)
                      year = year + 1900 * (year lt 100)
                      iuegettok,temp,':',hour
                      iuegettok,temp,':',minute
                      sec = float(strtrim(strmid(temp,0,5)))
;
;	             convert to day of year from month/day_of_month
;                    correction for leap years
;                    if (fix(year) mod 4) eq 0 then days(2) = 29   ;add 1 to feb
;
                     if ((year mod 4 eq 0) and (year mod 100 ne 0)) or $
                        (year mod 400 eq 0) then days(2) = 29
;
; 	        determine month number
;
                     month_name = strupcase(month_name)
                     for mon = 1,12 do begin
			if month_name eq months(mon) then goto,found
                     endfor
                     print,'DATECONV -- invalid month name specified
                     retall
                     found:
	             day = fix(total(days(0:mon-1))+day_of_month) ; day of yr
                  endif else begin      ; new FITS format
                     iuegettok,temp,'-',year
                     iuegettok,temp,'-',mon
                     if ((year mod 4 eq 0) and (year mod 100 ne 0)) or $
                        (year mod 400 eq 0) then days(2) = 29
                     day_of_month = strmid(temp,0,2)
	             day = fix(total(days(0:mon-1))+day_of_month) ; day of yr
                     if (strlen(temp) gt 2) then begin
                        temp = strmid(temp,3,strlen(temp)-3)
                        iuegettok,temp,':',hour
                        iuegettok,temp,':',minute
                        sec = float(temp)
                     endif
                  endelse
                endelse
	   end

	4 : begin			;spacecraft time
		sc = double(date)
		sc = sc + (sc lt 0.0)*65536.	;get rid of neg. numbers 
;
;	     Determine total number of secs since midnight, JAN. 1, 1979
;
		secs = sc(2)/64 + sc(1)*1024 + sc(0)*1024*65536.
		secs = secs/8192.0d0		;convert from spacecraft units 
;
;	     Determine number of years 
;
		mins = secs/60.
		hours = mins/60.
		totdays = hours/24.
		years = totdays/365.
		years = fix(years)
;
;	     Compute number of leap years past 
;
		leapyears = (years+2)/4
;
; 	    Compute day of year 
;
		day = fix(totdays-years*365.-leapyears)
;
; 	    Correct for case of being right at end of leapyear
;
		if day lt 0 then begin
		  day = day+366
		  leapyears = leapyears-1
		  years = years-1
		end
;
;	     compute hour of day
;
		totdays = years*365.+day+leapyears
		hour = fix(hours - 24*totdays)
		tothours = totdays*24+hour
;
;	     compute minute
;
		minute = fix(mins-tothours*60)
		totmin = tothours*60+minute
;
;	     compute sec
;
		sec = secs-totmin*60
;
;	     compute actual year
;
;                year = years+79
		year = years+1979
;
;	     start day at one and not zero
;
		day=day+1
	   end
	5 : begin			;integer time
            leap = 0
            year = long(date/10000L)
            yrl = strlen(strtrim(string(date),2))
            case 1 of
               (yrl le 4): leap=0                       ; assume year is 1900
               (yrl eq 8): if ((year mod 4 eq 0) and (year mod 100 ne 0)) or $
               (year mod 400 eq 0) then leap=1 ; full year specified
               else: if (year mod 4 eq 0) then leap=1   ; assume 1900's
            endcase
            if (leap) then days(2) = 29 ; if year is a leap yr.
            mon = fix( (date - (year*10000L))/100 )
            tday = fix ( date - (year*10000L) - (mon*100) )
            day = fix(total(days(0:mon-1))) + tday
            year = year + 1900 * (year lt 100)
            end
endcase
;
; see if time entries need to be rounded up
; (hope day is lt 365!)
;
            if ((hour+minute+sec) gt 0.0) then begin
               if (sec ge 60.0) then begin
                  minute = minute + 1
                  sec = sec - 60.0
               endif
               if (minute ge 60.0) then begin
                  hour = hour + 1
                  minute = minute - 60.0
               endif
               if (hour ge 24.0) then begin
                  day = day + 1
                  hour = hour - 24.0
               endif 
           endif

;           ---------------------------------------
;
;   *****	Now convert to output format
;
;           ---------------------------------------
;
; is type a string
;
s = size(type)
if (s(0) ne 0) or (s(1) ne 7) then begin
	print,'DATECONV- Output type specification must be a string'
	retall
end
;
; check for leap year
;
if ((year mod 4 eq 0) and (year mod 100 ne 0)) or $
    (year mod 400 eq 0) then days(2) = 29

case strmid(strupcase(type),0,1) of

  	'V' : begin				;vector output
		out = fltarr(5)
;                year = year + 1900 * (year lt 1900)
		out(0) = year
		out(1) = day
		out(2) = hour
		out(3) = minute
		out(4) = sec
	     end
 
	'R' : begin				;floating point scalar
		out = sec/24.0d0/60./60. + minute/24.0d0/60. + hour/24.0d0 $
	   		+  day + (year-1900.)*1000d0
	      end

	'S' : begin				;string output 
;
;	     correction for leap years
;		if form ne 3 then $	;Was it already done?
;			if (fix(year) mod 4) eq 0 then days(2) = 29
;
;
;	     check for valid day
;
		if (day lt 1) or (day gt total(days)) then begin
		   print,'DATE1-- There are only',total(days),' in year ',year
		   retall
		end
;
;	     find month which day occurs
;
		day_of_month = day
		month_num = 1
		while day_of_month gt days(month_num) do begin
			day_of_month = day_of_month - days(month_num)
			month_num = month_num+1
		end

		month_name = months(month_num)
;
;	     encode into ascii_date
;
;		year = year+1900
                out = string(day_of_month,'(i2)') +'-'+ month_name +'-' + $
			string(year,'(i4)') + ' '+ $
			string(hour,'(i2)') +':'+ $
			strmid(string(minute+100,'(i3)'),1,2) + ':'+ $
			strmid(string(sec+100,'(f6.2)'),1,5)
  	   end

	'F' : begin				;FITS string output 
;
;	     correction for leap years
;
;		if form ne 3 then $	;Was it already done?
;			if (fix(year) mod 4) eq 0 then days(2) = 29
;
;	     check for valid day
;
		if (day lt 1) or (day gt total(days)) then begin
                  print,'DATECONV-- There are only',total(days),' in year ',year
                  retall
		end
;
;	     find month which day occurs
;
		day_of_month = day
		month_num = 1
		while day_of_month gt days(month_num) do begin
			day_of_month = day_of_month - days(month_num)
			month_num = month_num+1
		endwhile
		month_name = months(month_num)
;
;	     encode into ascii_date
;
                if (day_of_month lt 10) then $
                   sdom = '0'+string(day_of_month,'(i1)') else $
                   sdom = string(day_of_month,'(i2)')
                if (month_num lt 10) then $
                   mn = '0' + string(month_num,'(i1)') else $
                   mn = string(month_num,'(i2)')                 
;                out = sdom + '/' + mn + '/' + (string(year,'(i2)')
                if (year ge 1900) and (year lt 2000) then year = year - 1900
                out = sdom + '/' + mn + '/' + strtrim(string(year),2)
                if (year eq 0) then out = out + '0'  ; keep 2 zeroes or year
  	   end
	'N' : begin				;NEWFITS string output 
;
;	     check for valid day
;
		if (day lt 1) or (day gt total(days)) then begin
                  print,'DATECONV-- There are only',total(days),' in year ',year
                  retall
		end
;
;	     find month which day occurs
;
		day_of_month = day
		month_num = 1
		while day_of_month gt days(month_num) do begin
			day_of_month = day_of_month - days(month_num)
			month_num = month_num+1
		endwhile
;
;	     encode into ascii_date
;
                if (day_of_month lt 10) then $
                   sdom = '0'+string(day_of_month,'(i1)') else $
                   sdom = string(day_of_month,'(i2)')
                if (month_num lt 10) then $
                   mn = '0' + string(month_num,'(i1)') else $
                   mn = string(month_num,'(i2)')                 
                out = strtrim(string(year),2) + '-' + mn + '-' + sdom
                if ((hour+minute+sec) ne 0.0) then begin
                   out = out+'T'+ $
                   strmid(string(hour+100,'(i3)'),1,2) +':'+ $
                   strmid(string(minute+100,'(i3)'),1,2) + ':'+ $
                   strmid(string(sec+100,'(f6.2)'),1,5)
                endif
  	   end
	'I' : begin				;Integer output 
;
;	     correction for leap years
;
;		if form ne 3 then $	;Was it already done?
;			if (fix(year) mod 4) eq 0 then days(2) = 29
;
;	     check for valid day
;
		if (day lt 1) or (day gt total(days)) then begin
                  print,'DATECONV-- There are only',total(days),' in year ',year
	          retall
		end
;
;	     find month which day occurs
;
		day_of_month = day
		month_num = 1
		while day_of_month gt days(month_num) do begin
			day_of_month = day_of_month - days(month_num)
			month_num = month_num+1
		end
                if (year lt 2000) and (year ge 1900) then year = year - 1900
                out =  (year * 10000L) + (month_num * 100) + day_of_month
              end
	else: begin			;invalid type specified
		print,'DATECONV-- Invalid output type specified'
		print,'	It must be ''REAL'', ''STRING'', or ''VECTOR'''
		retall
	      end
endcase
return   ; dateconv
;
; invalid input date error section
;
notvalid:
print,'DATECONV -- invalid input date specified'
retall
end
