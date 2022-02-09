PRO getTRIANGLES, dir, name_indx, FILE_EXIST, x1, y1, x2, y2, x3, y3, file_name=file_name, help=print_help
;+
; NAME:
;       getTRIANGLES
;
; PURPOSE:
;       This procedure reads in triangles
;
; CALLING SEQUENCE:

;
; INPUTS:
;
; KEYWORD PARAMETERS:
;      
; OUTPUTS:
;
; EXAMPLE:
;- 
;-------------------------OUTPUT IF HELP NEEDED-------------------------

if(keyword_set(print_help)) then begin 
 doc_library, 'getTRIANGLES'
 return
endif
;
;-----------------------DEFINE FILE-NAMES-------------------------------
;
IF(KEYWORD_SET (file_name)) THEN BEGIN
   fname=dir+'/'+file_name
ENDIF ELSE BEGIN
   dummy_name=STRING(name_indx, FORMAT='(I5.5)')
   fname=dir+'/spec_triangles_'+dummy_name+'.dat'
ENDELSE
;
;-----------------------CHECK IF FILES EXIST----------------------------
;
FILE_EXIST=FILE_TEST(fname)
IF(FILE_EXIST EQ 0) THEN BEGIN
   FILE_EXIST=0
   PRINT, 'FILE DOES NOT EXIST: ', fname
   RETURN
ENDIF
;
;------------------------READ IN DATA FROM FILES------------------------
;
readcol, fname, ip1, ip2, ip3, x1, y1, x2, y2, x3, y3
;
;
end
