PRO CONTOURPLOTS_MAP, ARR2D, ARRX, ARRY, center_lat, center_lon, $
                      plimit, titleStr, $
                      NCOLORS, NLEVELS_ISO, BOTTOM, $
                      LEVELS_CONTOUR, LEVELS_ISO, C_COLORS, $
                      CB_INDX, CB_TICKMARK_NAME, ctable=CTABLE, isoc=ISOC, $
                      ctitleStr=ctitleStr

IF(NOT KEYWORD_SET(ctable)) THEN BEGIN
;default:
   loadct, 13, ncolors=ncolors, bottom=bottom
ENDIF ELSE BEGIN
   loadct, ctable, ncolors=ncolors, bottom=bottom
ENDELSE
;
;set up orthographic map projection
map_set, center_lat, center_lon, $
         position=[.1, .1, .8, .95], $
         title=titleStr, $
         limit=plimit, $
         /orthographic, $
         /horizon, $
         /isotropic
;         xstyle=1, ystyle=1 , $
;
contour, arr2d, arrx, arry, $
         levels=LEVELS_CONTOUR, $
         c_color=c_colors, $
         /overplot, $
         /cell_fill

IF(KEYWORD_SET(ISOC)) THEN BEGIN
   CONTOUR, arr2d, arrx, arry, $
         levels=LEVELS_ISO, $
         color=0, $
         charsize=2., $
         c_labels=0, $
         /OVERPLOT
ENDIF

if(keyword_set(ctitleStr)) then begin
   COLORBAR, NCOLORS=NCOLORS, /VERTICAL, $
             divisions=nlevels_iso-1, $
             tickvalues=CB_INDX, $
             tickname=CB_TICKMARK_NAME, $
             title=ctitleStr, $
             POSITION=[.82, .1, .87, .95], /RIGHT
endif else begin
   COLORBAR, NCOLORS=NCOLORS, /VERTICAL, $
             divisions=nlevels_iso-1, $
             tickvalues=CB_INDX, $
             tickname=CB_TICKMARK_NAME, $
             POSITION=[.82, .1, .87, .95], /RIGHT
endelse   

;
;overplot latitude and azimuth grid
;calculate number of grid points within plot
theta_dum=arry(where(arry ge plimit(0) and arry le plimit(2)))
nt_dum=n_elements(theta_dum)
phi_dum=arrx(where(arrx ge plimit(1) and arrx le plimit(3)))
np_dum=n_elements(phi_dum)
;
;maximum number of labels is 6
max_label_grid=6
;
;calculate increment of labels to be drawn
nlabel=max([nt_dum/max_label_grid, np_dum/max_label_grid, 1])
;
;map_grid, $
;   /box_axes, $       ;grid labels are at outer box
;  label=0, $;         ;=nlabel for each nlabel-th grid line is labeled
;   glinestyle=1, $    ;line style of grid
;;   latdel=4., $     ;latitude grid increment in degrees
;;   londel=4.,$      ;longitude grid increment in degrees
;   lons=arrx, $    ;user specified array of grid values
;   lats=arry
;;
END
