pro contourplots_map_2triple, arr2d_1, arr2d_2, arr2d_3, arr2d_4, arr2d_5, arr2d_6, $
                              arrphi_1, arrphi_2, arrphi_3, arrphi_4, arrphi_5, arrphi_6, $
                              arrth_1, arrth_2, arrth_3, arrth_4, arrth_5, arrth_6, $
                              ncolors, nlevels_iso, bottom, $
                              levels_contour, levels_iso, c_colors, $
                              cb_indx, cb_tickmark_name, $
                              plimit=plimit, $
                              center_lat=center_lat, $
                              center_lon=center_lon, $
                              titlestr1=titlestr1, $
                              titlestr2=titlestr2, $
                              titlestr3=titlestr3, $
                              titlestr4=titlestr4, $
                              titlestr5=titlestr5, $
                              titlestr6=titlestr6, $
                              xtitlestr=xtitlestr, $
                              ytitlestr=ytitlestr, $
                              ctitlestr=ctitlestr, $
                              isoc=isoc, $
                              ctable=ctable, $
                              isotropic=isotropic, $
                              cb_top=cb_top, $
                              help=print_help
;+
; NAME:
;       contourplots_map_2triple
; PURPOSE:
;	This procedure creates six coloured contour-plots on a 2rowsx3cols grid, along with
;	a colorbar on the top, optimized for charsize=1.2 (ps-output with
;	xsize=17.8 cm) (A&A paper style).
;       
; CALLING SEQUENCE:
;	contourplots_map_2triple, arr2d_1, arr2d_2, arr2d_3, arr2d_4, arr2d_5, arr2d_6, $
;                                 arrth_1, arrth_2, arrth_3, arrth_4, arrth_5, arrth_6, $
;                                 arrphi_1, arrphi_2, arrphi_3, arrphi_4, arrphi_5, arrphi_6, $
;	                          ncolors, nlevels_iso, bottom, levels_contour, $
;	                          levels_iso, c_colors, cb_inds, cb_tickmark_name
;
; INPUTS:
;       arr2d_1:  2-dimensional array to be plotted as contour (upper left)
;       arr2d_2:  2-dimensional array to be plotted as contour (upper right)
;       arr2d_3:  2-dimensional array to be plotted as contour (middle left)
;       arr2d_4:  2-dimensional array to be plotted as contour (middle right)
;       arr2d_5:  2-dimensional array to be plotted as contour (lower left)
;       arr2d_6:  2-dimensional array to be plotted as contour (lower right)
;       arrth_1:   theta-coordinate of array 1
;       arrth_2:   theta-coordinate of array 2
;       arrth_3:   theta-coordinate of array 3
;       arrth_4:   theta-coordinate of array 4
;       arrth_5:   theta-coordinate of array 5
;       arrth_6:   theta-coordinate of array 6
;       arrphi_1:   phi-coordinate of array 1
;       arrphi_2:   phi-coordinate of array 2
;       arrphi_3:   phi-coordinate of array 3
;       arrphi_4:   phi-coordinate of array 4
;       arrphi_5:   phi-coordinate of array 5
;       arrphi_6:   phi-coordinate of array 6
;       all others:   color-specifiers as output from get_contour_colors.pro
;                     or get_contour_colors2.pro 
;
; KEYWORD PARAMETERS:
;       help:   Set this keyword (flag) to show the documentation of this
;               procedure
;       plim:   Set this keyword to set the map-limits
;       center_lat: Set this keyword to set the center latitude of map
;       center_lon: Set this keyword to set the center longitude of map
;       titleStr(1,2,3,4):  Set this keyword to set the title of the individual plots
;       xtitleStr: Set this keyword to set the x-title of the plot
;       ytitleStr: Set this keyword to set the y-title of the plot
;       ctitleStr: Set this keyword to set the colorbar-title of the plot
;       isoc:   Set this keyword (flag) to overplot iso-contours as straight lines
;       ctable: Set this keyword to set the color-table
;       isotropic: Set this keyword (flag) to make an isotropic plot
;       cb_top:    Set this keyword (flag) to put the colorbar on top
;
; OUTPUTS:
;
; EXAMPLE:
;       contourplots_map_2triple, arr2d1, arr2d2, arr2d3, arr2d4, arr2d5, arr2d6, $
;                                 th1, th2, th3, th4, th5, th6, phi1, phi2,
;                                 phi3, phi4, phi5, phi6, ncolors, $
;                                 nlevels_iso, bottom, levels_contour, levels_iso, $
;                                 c_colors, cb_indx, cb_tickmark_name, $
;                                 plimit=[-90.d0,0.d0,90.d0,360.d0], $
;                                 center_lat=45., center_lon=45.
;                                 ctable=13, /isotropic
;
; NOTES:
;-
;
;-------------------------OUTPUT IF HELP NEEDED-------------------------
;
if(keyword_set(print_help)) then begin
   doc_library, 'contourplots_map_3double'
   return
endif
;
;-----------------------------------------------------------------------
;
;set colortable
if(not keyword_set(ctable)) then begin
;default:
   loadct, 13, ncolors=ncolors, bottom=bottom
endif else begin
   loadct, ctable, ncolors=ncolors, bottom=bottom
endelse
;
;-----------------------------first plot--------------------------------
;
;make smaller x-margin for multi plots
xmarg_default=!x.margin
ymarg_default=!y.margin
!x.margin=4.*!x.margin/5.
;!x.margin=5.*!x.margin/5.
;
;set the region where to set plot
if(keyword_set(cb_top)) then begin
  !P.REGION=[0.0, 0.415, 0.3333, 0.83]
endif else begin
   print, 'right color bar to be debugged'
   stop
endelse
;
;set up orthographic map projection
map_set, center_lat, center_lon, $
         title=titleStr1, $
         limit=plimit, $
         /orthographic, $
         /horizon, $
         charsize=1., $
         isotropic=isotropic
;
if(keyword_set(cb_top)) then begin
;set left position of colorbar
   xpos0_cb=!x.window(0)
endif
;
contour, arr2d_1, arrphi_1, arrth_1, $
         levels=LEVELS_CONTOUR, $
         c_color=c_colors, $
         /overplot, $
         /cell_fill

if(keyword_set(isoc)) then begin
   contour, arr2d_1, arrphi_1, arrth_1, $
         levels=levels_iso, $
         color=0, $
         charsize=2., $
         c_labels=0, $
         /overplot
endif
;
;---------------------------second plot---------------------------------
;
;set the region where to set plot
if(keyword_set(cb_top)) then begin
   !P.REGION=[0.3333, 0.415, 0.6666, 0.83]
endif else begin
   print, 'right color bar to be debugged'
   stop
endelse
;

;set up orthographic map projection
map_set, center_lat, center_lon, $
         title=titleStr2, $
         limit=plimit, $
         /orthographic, $
         /horizon, $
         charsize=1., $
         isotropic=isotropic, /noerase
;
contour, arr2d_2, arrphi_2, arrth_2, $
         levels=LEVELS_CONTOUR, $
         c_color=c_colors, $
         /overplot, $
         /cell_fill

if(keyword_set(isoc)) then begin
   contour, arr2d_2, arrphi_2, arrth_2, $
         levels=levels_iso, $
         color=0, $
         charsize=2., $
         c_labels=0, $
         /overplot
endif
;
;---------------------------third plot---------------------------------
;
;set the region where to set plot
if(keyword_set(cb_top)) then begin
   !P.REGION=[0.6666, .415, 1., 0.83]
endif else begin
   print, 'right color bar to be debugged'
   stop
endelse
;
;
;set up orthographic map projection
map_set, center_lat, center_lon, $
         title=titleStr3, $
         limit=plimit, $
         /orthographic, $
         /horizon, $
         charsize=1., $
         isotropic=isotropic, /noerase

if(keyword_set(cb_top)) then begin
;set right position of colorbar
   xpos1_cb=!x.window(1)
endif
;
contour, arr2d_3, arrphi_3, arrth_3, $
         levels=LEVELS_CONTOUR, $
         c_color=c_colors, $
         /overplot, $
         /cell_fill

if(keyword_set(isoc)) then begin
   contour, arr2d_3, arrphi_3, arrth_3, $
         levels=levels_iso, $
         color=0, $
         charsize=2., $
         c_labels=0, $
         /overplot
endif
;
;---------------------------foruth plot--------------------------------
;
;set the region where to set plot
if(keyword_set(cb_top)) then begin
   !P.REGION=[0., 0., .3333, 0.415]
endif else begin
   print, 'right color bar to be debugged'
   stop
endelse
;

;set up orthographic map projection
map_set, center_lat, center_lon, $
;         position=[.1, .1, .8, .95], $
         title=titleStr4, $
         limit=plimit, $
         /orthographic, $
         /horizon, $
         charsize=1., $
         isotropic=isotropic, /noerase
;
contour, arr2d_4, arrphi_4, arrth_4, $
         levels=LEVELS_CONTOUR, $
         c_color=c_colors, $
         /overplot, $
         /cell_fill

if(keyword_set(isoc)) then begin
   contour, arr2d_4, arrphi_4, arrth_4, $
         levels=levels_iso, $
         color=0, $
         charsize=2., $
         c_labels=0, $
         /overplot
endif
;
;---------------------------fifth plot---------------------------------
;
;set the region where to set plot
if(keyword_set(cb_top)) then begin
   !P.REGION=[0.3333, 0., .6666, 0.415]
endif else begin
   print, 'right color bar to be debugged'
   stop
endelse
;

;set up orthographic map projection
map_set, center_lat, center_lon, $
;         position=[.1, .1, .8, .95], $
         title=titleStr5, $
         limit=plimit, $
         /orthographic, $
         /horizon, $
         charsize=1., $
         isotropic=isotropic, /noerase
;
contour, arr2d_5, arrphi_5, arrth_5, $
         levels=levels_contour, $
         c_color=c_colors, $
         /overplot, $
         /cell_fill

if(keyword_set(isoc)) then begin
   contour, arr2d_5, arrphi_5, arrth_5, $
         levels=levels_iso, $
         color=0, $
         charsize=2., $
         c_labels=0, $
         /overplot
endif
;
;---------------------------sixth plot----------------------------------
;
;set the region where to set plot
if(keyword_set(cb_top)) then begin
   !P.REGION=[0.6666, 0., 1., 0.415]
endif else begin
   print, 'right color bar to be debugged'
   stop
endelse
;

;set up orthographic map projection
map_set, center_lat, center_lon, $
;         position=[.1, .1, .8, .95], $
         title=titleStr6, $
         limit=plimit, $
         /orthographic, $
         /horizon, $
         charsize=1., $
         isotropic=isotropic, /noerase
;
contour, arr2d_6, arrphi_6, arrth_6, $
         levels=levels_contour, $
         c_color=c_colors, $
         /overplot, $
         /cell_fill

if(keyword_set(isoc)) then begin
   contour, arr2d_6, arrphi_6, arrth_6, $
         levels=levels_iso, $
         color=0, $
         charsize=2., $
         c_labels=0, $
         /overplot
endif
;
;-----------------------------colorbar----------------------------------
;
if(keyword_set(cb_top)) then begin
   !P.REGION=[0.0, 0.83, 1., 1.]
endif else begin
   print, 'right color bar to be debugged'
   stop
endelse
;
region_width=!p.region(2)-!p.region(0)
region_height=!p.region(3)-!p.region(1)
pos=fltarr(4)
;
if(keyword_set(cb_top)) then begin
;x-arrangement depends on the position of the plot
   pos(0)=xpos0_cb
   pos(2)=xpos1_cb
;y-arrangement depends on region
   pos(1)=!p.region(1)+0.15*region_height
   pos(3)=!p.region(3)-0.55*region_height
endif else begin
;x-arrangement depends on region
   pos(0)=!p.region(0);+0.1*region_width
   pos(2)=!p.region(2)-0.8*region_width
;y-arrangement depends on the position of the plot
   pos(1)=!y.window(0)
   pos(3)=!y.window(1)
endelse
;
if(keyword_set(cb_top)) then begin
   vert=0
   horiz=1
   right=0
   top=1
endif else begin
   vert=1
   horiz=0
   right=1
   top=0
endelse
;
;
if(keyword_set(ctitleStr)) then begin
;/VERTICAL, /RIGHT

   COLORBAR, NCOLORS=NCOLORS, $
             divisions=nlevels_iso-1, $
             tickvalues=CB_INDX, $
             tickname=CB_TICKMARK_NAME, $
             POSITION=pos, $
             title=ctitleStr, $
             charsize=1.2, $
             VERTICAL=vert, $
             right=right, $
             top=top, $
             HORIZONTAL=horiz
endif else begin
   COLORBAR, NCOLORS=NCOLORS, $
             divisions=nlevels_iso-1, $
             tickvalues=CB_INDX, $
             tickname=CB_TICKMARK_NAME, $
             POSITION=pos, $
             charsize=1.2, $
             VERTICAL=vert, $
             right=right, $
             top=top, $
             HORIZONTAL=horiz
endelse
;
;
;set back to default
!x.margin=xmarg_default
!y.margin=ymarg_default
!p.region=0.
;
;
;
end
