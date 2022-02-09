PRO CONTOURPLOTS, ARR2D, ARRX, ARRY, RANGEX, RANGEY, $
                  titleStr, xtitleStr, ytitleStr, $
                  NCOLORS, NLEVELS_ISO, BOTTOM, $
                  LEVELS_CONTOUR, LEVELS_ISO, C_COLORS, $
                  CB_INDX, CB_TICKMARK_NAME, ctable=CTABLE, isoc=ISOC, ctbl_file=CTBL_FILE, $
                  isotropic=ISOTROPIC, irregular=IRREGULAR, ctitleStr=ctitleStr
;
;PREPARING DATA TO PLOT AXES
AXISX=[MIN(ARRX), MAX(ARRX)]
AXISY=[MIN(ARRY), MAX(ARRY)]
ORIGIN=[0., 0.]
;
;
IF(NOT KEYWORD_SET(ctable)) THEN BEGIN
;default:
   loadct, 13, ncolors=ncolors, bottom=bottom
ENDIF ELSE BEGIN
   IF(KEYWORD_SET(ctbl_file)) THEN BEGIN
      loadct, ctable, ncolors=ncolors, bottom=bottom, file=ctbl_file
   ENDIF ELSE BEGIN
      loadct, ctable, ncolors=ncolors, bottom=bottom
   ENDELSE
ENDELSE
;
;
CONTOUR, ARR2D, ARRX, ARRY, $
;POSITION=[x0, y0, x1, y1]
      position=[.1, .1, .8, .95], $
      levels=LEVELS_CONTOUR, $
      c_color=c_colors, $
      /cell_fill, $
      title=titleStr, $
      xtitle=xtitleStr, $
      ytitle=ytitleStr, $
      xrange=RANGEX, $
      yrange=RANGEY, $
      xstyle=1, ystyle=1 , $
      isotropic=ISOTROPIC, $
      irregular=IRREGULAR, $
      charsize=1.2
;
IF(KEYWORD_SET(ISOC)) THEN BEGIN
   CONTOUR, ARR2D, ARRX, ARRY, $
      levels=LEVELS_ISO, $
      color=0, $
      charsize=2., $
      c_labels=0, $
      /OVERPLOT
ENDIF
;
if(keyword_set(ctitleStr)) then begin
   COLORBAR, NCOLORS=NCOLORS, /VERTICAL, $
             divisions=nlevels_iso-1, $
             tickvalues=CB_INDX, $
             tickname=CB_TICKMARK_NAME, $
             POSITION=[.82, .1, .87, .95], $
             title=ctitleStr, /RIGHT, $
             charsize=1.2
endif else begin
   COLORBAR, NCOLORS=NCOLORS, /VERTICAL, $
             divisions=nlevels_iso-1, $
             tickvalues=CB_INDX, $
             tickname=CB_TICKMARK_NAME, $
             POSITION=[.82, .1, .87, .95], /RIGHT, $
             charsize=1.2
;             POSITION=[.9, .1, .95, .95], /RIGHT
endelse

;OPLOT, AXISX, ORIGIN, line=0, color=0
;OPLOT, ORIGIN, AXISY, line=0, color=0

END
