;****************************************************************************
;+
;*NAME:
;
;    FITSCON 
;
;*PURPOSE:
;
;    To convert FITS binary table data types (i.e. IEEE format)
;    to format used on current cpu.
;
;*CALLING SEQUENCE:
;
;    FITSCON,BYTE_EQ,TYPE,VARIABLE,TDIM=S
;
;*PARAMETERS:
;
;    BYTE_EQ (REQ) (I) (1) (B) 
;        The byte equivalent vector of the data variable to be converted.
;        (BYTE_EQ must be a vector.)
;
;    TYPE (REQ) (I) (0) (S)
;        Data type of VARIABLE described as a single character. Allowed
;        data types are: byte 'X', integer 'I', longword 'J', floating 
;        point 'E', double precision 'D', and string 'A'. 
;
;    VARIABLE (REQ) (O) (01) (BILFD)
;        Converted output vector.
;
;    TDIM (KEY) (I) (0) (S)
;        If specified, contains the value of the TDIMn keyword
;        describing the dimensions of a multi-dimensional array.
;
;*SIDE EFFECTS:
;
;*SYSTEM VARIABLES USED:
;
;*SUBROUTINES CALLED:
;
;	PARCHECK
;       FLAG_NAN
;       IUEGETTOK
;
;*EXAMPLE:
;
;*RESTRICTIONS:
;
;*NOTES:
;
;       Original version written for IDL version 1.
;       Data types are based on those allowed in the TFORM FITS keyword.
;       One element output vectors are converted to scalars.
;
;	tested with IDL Version 2.1.2 (sunos sparc)    11 Mar 93
;	tested with IDL Version 2.3.2 (vax vms)        11 Mar 93 
;       tested with IDL Version 2.1.2 (ultrix mipsel)  08 Nov 91
;       tested with IDL Version 2.2.0 (ultrix vax)     08 Nov 91
;
;*MODIFICATION HISTORY:
;
;	Version 1 Randy Thompson 2/8/91 (based on VTOS by John Hoegy)
;       3-13-91 RWT use algorithms derived by WTT in SUN2VAX procedure.
;       3-14-91 RWT replace SWAP with routine called SWAP_BYTES
;	3-27-91	PJL modified for unix/sun; renamed from FITSTOV to FITSSUN;
;		added PARCHECK; deleted SWAP_BYTES
;       4-23-91 GRA changed references to FITSTOV to FITSSUN (n_par eq 0)
;       7-25-91 RWT rename FITSCON, use TRANS_BYTES to allow conversion
;               to vms, Ultrix, unix, or DOS systems, and convert 1 element
;               arrays to scalars.
;      11/07/91 GRA defined cpupar = !version.arch for trans_bytes 
;               parameter; tested on sun, vax, and dec.
;       6/24/92 RWT change 383 to 386 for flagging IBM pcs
;       7/6/92  RWT properly handle B, P, C, & M formats
;      10/19/92 RWT add FLAGNAN to J,E,D,C, & M data types
;       3/11/93 RWT replace TRANS_BYTES with intrinsic byteorder command
;       9/22/93 RWT add multi-dimensional array support
;-
;****************************************************************************
 pro fitscon,byte_eq,type,variable,tdim=s
;
 if n_params(0) eq 0 then begin
    print,'PRO FITSCON,BYTE_EQ,TYPE,VARIABLE,tdim=s'
    retall
 endif  ; n_params(0)
 parcheck,n_params(0),3,'FITSCON'
;
 tmpv = byte_eq
 byte_elems = n_elements(tmpv)
 var_type = strupcase(type)
;
; check for multi-dimensional array 
;
 dimtest = 0
 if (keyword_set(s)) then begin
    dim = intarr(8)
    i = 0
    if (s ne '') then repeat begin
         iuegettok,s,',',tmp               ; search for "," delimiter
         dim(i) = fix(tmp)
         i = i + 1
    end until (s eq '') or (i gt 7)
    ind = where(dim,npts)                  ; remove 0's
    dim = dim(0:npts-1)
    if (npts ne 0) then dimtest = 1
    if (i eq 8) then print,'Warning: IDL limited to <= 8 dimensions'
endif
;    
 case var_type of
    'X': begin                      ; byte
         variable = byte_eq
         if (byte_elems eq 1) then variable = variable(0)
         if (dimtest) then begin
             tmp = make_array(dimension=dim,/byte)
             tmp(0) = variable
             variable = tmp
         endif
         return
         endcase  ; X
    'B': begin                      ; unsigned byte(?)
         variable = byte_eq
         if (byte_elems eq 1) then variable = variable(0)
         if (dimtest) then begin
             tmp = make_array(dimension=dim,/byte)
             tmp(0) = variable
             variable = tmp
         endif
         return
         endcase  ; B
    'I': begin                      ; integer
         var_elems = byte_elems / 2L
         variable = intarr(var_elems)
         variable = fix(tmpv,0,var_elems)
         byteorder,variable,/NTOHS
         if (var_elems eq 1) then variable = variable(0)
         if (dimtest) then begin
             tmp = make_array(dimension=dim,/int)
             tmp(0) = variable
             variable = tmp
         endif
         return
         endcase  ; I
    'J': begin                      ; longword
         var_elems = byte_elems / 4L
         variable = long(tmpv,0,var_elems)                ;convert for flagnan
         flagnan,variable,ind,count                       ;find NaNs
         variable = long(tmpv,0,var_elems)
         byteorder,variable,/NTOHL
         if (count ne 0) then begin
            print,strtrim(string(count),2),' special characters found.'
            print,' reassigned value = -99'
            variable(ind) = -99L        ;set NaN's to -99
         end
         if (var_elems eq 1) then variable = variable(0)
         if (dimtest) then begin
             tmp = make_array(dimension=dim,/long)
             tmp(0) = variable
             variable = tmp
         endif
         return
         endcase  ; J
    'P': begin                      ; var. length array(?)
         var_elems = byte_elems / 4L
         variable = long(tmpv,0,var_elems)                ;convert for flagnan
         flagnan,variable,ind,count                       ;find NaNs
         variable = long(tmpv,0,var_elems)
         byteorder,variable,/NTOHL
         if (count ne 0) then begin
            print,strtrim(string(count),2),' special characters found.'
            print,' reassigned value = -99'
            variable(ind) = -99L        ;set NaN's to -99
         end
         if (var_elems eq 1) then variable = variable(0)
         endcase  ; P
    'E': begin         		; floating point
         var_elems = byte_elems / 4L
         variable = float(tmpv,0,var_elems)               ;convert for flagnan
         flagnan,variable,ind,count                       ;find NaNs
         variable(0) = float(tmpv, 0, var_elems)
         byteorder,variable,/XDRTOF
         if (count ne 0) then begin
            print,strtrim(string(count),2),' special characters found.'
            print,' reassigned value = -99.9'
            variable(ind) = -99.9        ;set NaN's to -99
         end
         if (var_elems eq 1) then variable = variable(0)
         if (dimtest) then begin
             tmp = make_array(dimension=dim,/float)
             tmp(0) = variable
             variable = tmp
         endif
         return
         endcase  ; E
    'C': begin         		; complex(?)
         var_elems = byte_elems / 4L
         variable = float(tmpv,0,var_elems)               ;convert for flagnan
         flagnan,variable,ind,count                       ;find NaNs
         variable(0) = float(tmpv, 0, var_elems)
         byteorder,variable,/XDRTOF
         if (count ne 0) then begin
            print,strtrim(string(count),2),' special characters found.'
            print,' reassigned value = -99.9'
            variable(ind) = -99.9        ;set NaN's to -99
         end
         if (var_elems eq 1) then variable = variable(0)
         if (dimtest) then begin
             tmp = make_array(dimension=dim,/complex)
             tmp(0) = variable
             variable = tmp
         endif
         return
         endcase  ; C
    'D': begin    	         	; double precision
         var_elems = byte_elems / 8L 
         variable = double(tmpv,0,var_elems)              ;convert for flagnan
         flagnan,variable,ind,count                       ;find NaNs
         variable(0) = double(tmpv, 0, var_elems)
         byteorder,variable,/XDRTOD
         if (count ne 0) then begin
            print,strtrim(string(count),2),' special characters found.'
            print,' reassigned value = -99.9'
            variable(ind) = -99.9        ;set NaN's to -99
         end
         if (var_elems eq 1) then variable = variable(0)
         if (dimtest) then begin
             tmp = make_array(dimension=dim,/double)
             tmp(0) = variable
             variable = tmp
         endif
         return
         endcase  ; D
    'M': begin    	         	; complex double precision(?)
         var_elems = byte_elems / 8L 
         variable = double(tmpv,0,var_elems)              ;convert for flagnan
         flagnan,variable,ind,count                       ;find NaNs
         variable(0) = double(tmpv, 0, var_elems)
         byteorder,variable,/XDRTOD
         if (count ne 0) then begin
            print,strtrim(string(count),2),' special characters found.'
            print,' reassigned value = -99.9'
            variable(ind) = -99.9        ;set NaN's to -99
         end
         if (var_elems eq 1) then variable = variable(0)
         if (dimtest) then begin
             tmp = make_array(dimension=dim,/double)
             tmp(0) = variable
             variable = tmp
         endif
         return
         endcase  ; M
    'A': begin                      ; string
         variable = string(byte_eq)
         if (dimtest) then begin
             tmp = make_array(dimension=dim,/string)
             tmp(0) = variable
             variable = tmp
         endif
         return
         endcase  ; A
    'L': begin                      ; logical
         variable = string(byte_eq)
         if (dimtest) then begin
             tmp = make_array(dimension=dim,/string)
             tmp(0) = variable
             variable = tmp
         endif
         return
         endcase  ; L
   else: begin 			; unknown
         print,'Data type',var_type,' unknown, routine FITSCON'
         retall
         endelse
     endcase  ; var_type
 return
 end  ; fitscon
