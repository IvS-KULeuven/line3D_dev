PRO CONTOURPLOTS_TRIPLE, arr2d_1, arr2d_2, arr2d_3, $
                         arrx_1, arrx_2, arrx_3, $
                         arry_1, arry_2, arry_3, $
                         NCOLORS, NLEVELS_ISO, BOTTOM, $
                         LEVELS_CONTOUR, LEVELS_ISO, C_COLORS, $
                         CB_INDX, CB_TICKMARK_NAME, $
                         xlim=xlim, $
                         ylim=ylim, $
                         titleStr1=titleStr1, $
                         titleStr2=titleStr2, $
                         titleStr3=titleStr3, $
                         xtitleStr=xtitleStr, $
                         ytitleStr=ytitleStr, $
                         ctitleStr=ctitleStr, $
                         isoc=ISOC, $
                         ctable=CTABLE, $
                         ctbl_file=CTBL_FILE, $
                         isotropic=ISOTROPIC, $
                         irregular=IRREGULAR, $
                         cb_top=cb_top, $
                         ralfven1=ralfven1, $
                         ralfven2=ralfven2, $
                         ralfven3=ralfven3, $
                         oplot_grid=oplot_grid, $
                         oplot_circle=oplot_circle, $
                         oplot1_x1=oplot1_x1, oplot1_y1=oplot1_y1, $
                         oplot1_x2=oplot1_x2, oplot1_y2=oplot1_y2, $
                         oplot1_x3=oplot1_x3, oplot1_y3=oplot1_y3, $
                         oplot2_x1=oplot2_x1, oplot2_y1=oplot2_y1, $
                         oplot2_x2=oplot2_x2, oplot2_y2=oplot2_y2, $
                         oplot2_x3=oplot2_x3, oplot2_y3=oplot2_y3, $
                         oplot3_x1=oplot3_x1, oplot3_y1=oplot3_y1, $
                         oplot3_x2=oplot3_x2, oplot3_y2=oplot3_y2, $
                         oplot3_x3=oplot3_x3, oplot3_y3=oplot3_y3
;
;+
; NAME:
;       contourplots_triple
; PURPOSE:
;	This procedure creates three coloured contour-plots in a row, along with
;	a colorbar on the right (or on top), optimized for charsize=.8 (ps-output with
;	xsize=17.8 cm) for full page-width (A&A paper style).
;       
; CALLING SEQUENCE:
;	contourplots_triple, arr2d_1, arr2d_2, arr2d_3, arrx_1, arrx_2, $
;                            arrrx_3, arry_1, arry_2, arry_3, $
;	                     ncolors, nlevels_iso, bottom, levels_contour, $
;	                     levels_iso, c_colors, cb_inds, cb_tickmark_name
;
; INPUTS:
;       arr2d_1:  2-dimensional array to be plotted as contour (left)
;       arr2d_2:  2-dimensional array to be plotted as contour (middle)
;       arr2d_3:  2-dimensional array to be plotted as contour (right)
;       arrx_1:   x-coordinate of array 1, can be 2d-array when irregular is specified
;       arrx_2:   x-coordinate of array 2, can be 2d-array when irregular is specified
;       arrx_3:   x-coordinate of array 3, can be 2d-array when irregular is specified
;       arry_1:   y-coordinate of array 1, can be 2d-array when irregular is specified
;       arry_2:   y-coordinate of array 2, can be 2d-array when irregular is specified
;       arry_3:   y-coordinate of array 3, can be 2d-array when irregular is specified
;       all others:   color-specifiers as output from get_contour_colors.pro
;                     or get_contour_colors2.pro 
;
; KEYWORD PARAMETERS:
;       help:   Set this keyword (flag) to show the documentation of this
;               procedure
;       xlim:   Set this keyword to set the x-range of the plot
;       ylim:   Set this keyword to set the y-range of the plot
;       titleStr(1,2,3):  Set this keyword to set the title of the individual plots
;       xtitleStr: Set this keyword to set the x-title of the plot
;       ytitleStr: Set this keyword to set the y-title of the plot
;       ctitleStr: Set this keyword to set the colorbar-title of the plot
;       isoc:   Set this keyword (flag) to overplot iso-contours as straight lines
;       ctable: Set this keyword to set the color-table
;       ctbl_file: Set this keyword to read the color-table from a specified
;                  file
;       isotropic: Set this keyword (flag) to make an isotropic plot
;       irregular: Set this keyword (flag) if input data is irregulary spaced
;       cb_top:    Set this keyword (flag) to put the colorbar on top
;       oplot_grid:  Set this keyword (flag) to overplot the (regular) grid
;       oplot_circle: radius of a circle that shall be overplotted
;       oplot1_x1, oplot1_y1: Specifies any x,y-array that shall be
;                             overplotted for plot '1' (left plot)
;       oplot1_x2, oplot1_y2: Specifies any x,y-array that shall be
;                             overplotted for plot '1' (left plot)
;       oplot2_x1, oplot2_y1: Specifies any x,y-array that shall be
;                             overplotted for plot '2' (middle plot)
;       oplot2_x2, oplot2_y2: Specifies any x,y-array that shall be
;                             overplotted for plot '2' (middle plot)
;       oplot3_x1, oplot3_y1: Specifies any x,y-array that shall be
;                             overplotted for plot '3' (right plot)
;       oplot3_x2, oplot3_y2: Specifies any x,y-array that shall be
;                             overplotted for plot '3' (right plot)
;
; OUTPUTS:
;
; EXAMPLE:
;       CONTOUR_PLOTS_DOUBLE, arr2d1, arr2d2, arr2d3, x1, x2, x3, y1, y2, y3, ncolors, $
;                             nlevels_iso, bottom, levels_contour, levels_iso, $
;                             c_colors, cb_indx, cb_tickmark_name, $
;                             xlim=[0.,2.], ylim=[0.,2.], $ 
;                             ctable=13, /isotropic, /irregular
;
; NOTES:
;       so far, xlim, ylim, and all isocontours, etc. are the same for both plots
;-
;
;-------------------------OUTPUT IF HELP NEEDED-------------------------
;
if(keyword_set(print_help)) then begin
   doc_library, 'contourplots_triple'
   return
endif
;
;-----------------------------------------------------------------------
;
chsize_title=.8
chsize_axes=.8
;
;-----------------------------------------------------------------------
;
;set the colortable
if (not keyword_set(ctable)) then begin
;default:
   loadct, 13, ncolors=ncolors, bottom=bottom
endif else begin
   if (keyword_set(ctbl_file)) then begin
      loadct, ctable, ncolors=ncolors, bottom=bottom, file=ctbl_file
   endif else begin
      loadct, ctable, ncolors=ncolors, bottom=bottom
   endelse
endelse
;
;---------------------------first plot----------------------------------
;
;make smaller x-margin for multi plots
xmarg_default=!x.margin
ymarg_default=!y.margin
!x.margin=3.*!x.margin/4.
;!x.margin=5.*!x.margin/5.

;set the region where to set plot
if(keyword_set(cb_top)) then begin
   !P.REGION=[0.03, 0.0, 0.3433, 0.78]
endif else begin
   !P.REGION=[0.0, 0.0, .27333, 1.]
endelse
;
;create contour plot
CONTOUR, ARR2D_1, ARRX_1, ARRY_1, $
         levels=LEVELS_CONTOUR, $
         c_color=c_colors, $
         /cell_fill, $
         title=titleStr1, $
         xtitle=xtitleStr, $
         ytitle=ytitleStr, $
         xrange=xlim, $
         yrange=ylim, $
         xstyle=1, ystyle=1 , $
         isotropic=ISOTROPIC, $
         irregular=IRREGULAR, $
         charsize=chsize_title, $     ;overall charsize (to set smaller title)
         xcharsize=chsize_axes, $     ;x-axis charsize
         ycharsize=chsize_axes        ;y-axis charsize
;
if(keyword_set(cb_top)) then begin
;set left position of colorbar
   xpos0_cb=!x.window(0)
endif
;
;overplot the isocontours, if specified
if (keyword_set(isoc)) then begin
   CONTOUR, arr2d_1, arrx_1, arry_1, $
            levels=LEVELS_ISO, $
            color=0, $
            charsize=2., $
            c_labels=0, $
            /OVERPLOT
endif
;
;overplot photosphere
if(n_elements(oplot_circle) ne 0) then begin
   loadct, 0
   plotcircle, oplot_circle, /oplot
endif
;
;overplot magnetic field lines
if(keyword_set(ralfven1)) then begin
   loadct, 0
   plotbfield, nlines=8, ralfven=ralfven1, /oplot
   plotbfield_neg, nlines=8, ralfven=ralfven1, /oplot
endif
;
;overplot the grid
if(keyword_set(oplot_grid)) then begin
   loadct, 0   
   for i=0, n_elements(arrx_1)-1 do begin
      oplot, [arrx_1(i),arrx_1(i)], [ylim(0),ylim(1)]
   endfor
   for i=0, n_elements(arry_1)-1 do begin
      oplot, [xlim(0),xlim(1)], [arry_1(i),arry_1(i)]
   endfor
;
endif
;
;overplot any given input array
if(n_elements(oplot1_x1) ne n_elements(oplot1_y1)) then begin
   print, 'error in contourplots_triple: nd of oplot1_x1 and oplot1_y1 not the same'
   stop
endif
if(n_elements(oplot1_x2) ne n_elements(oplot1_y2)) then begin
   print, 'error in contourplots_triple: nd of oplot1_x2 and oplot1_y2 not the same'
   stop
endif
if(n_elements(oplot1_x3) ne n_elements(oplot1_y3)) then begin
   print, 'error in contourplots_triple: nd of oplot1_x3 and oplot1_y3 not the same'
   stop
endif
if(n_elements(oplot1_x1) gt 1) then begin
   loadct, 0
   oplot, oplot1_x1, oplot1_y1, line=0, thick=2.
endif
if(n_elements(oplot1_x2) gt 1) then begin
   loadct, 0
   oplot, oplot1_x2, oplot1_y2, line=0, thick=2.
endif
if(n_elements(oplot1_x3) gt 1) then begin
   loadct, 0
   oplot, oplot1_x3, oplot1_y3, line=0, thick=2.
endif
;
;---------------------------second plot---------------------------------
;
;set the colortable
if (not keyword_set(ctable)) then begin
;default:
   loadct, 13, ncolors=ncolors, bottom=bottom
endif else begin
   if (keyword_set(ctbl_file)) then begin
      loadct, ctable, ncolors=ncolors, bottom=bottom, file=ctbl_file
   endif else begin
      loadct, ctable, ncolors=ncolors, bottom=bottom
   endelse
endelse
;
;set the region where to set plot
if(keyword_set(cb_top)) then begin
   !P.REGION=[0.3333, 0.0, 0.6566, 0.78]
endif else begin
   !P.REGION=[0.27333, 0.0, 0.54666, 1.]
endelse
;
;create contour plot
CONTOUR, ARR2D_2, ARRX_2, ARRY_2, $
         levels=LEVELS_CONTOUR, $
         c_color=c_colors, $
         /cell_fill, $
         title=titleStr2, $
         xtitle=xtitleStr, $
         ytitle=ytitleStr, $
         xrange=xlim, $
         yrange=ylim, $
         xstyle=1, ystyle=1 , $
         isotropic=ISOTROPIC, $
         irregular=IRREGULAR, $
         /noerase, $
         charsize=chsize_title, $       ;overall charsize (to set smaller title)
         xcharsize=chsize_axes, $     ;x-axis charsize
         ycharsize=chsize_axes        ;y-axis charsize
;
;overplot the isocontours, if specified
if (keyword_set(isoc)) then begin
   CONTOUR, arr2d_2, arrx_2, arry_2, $
            levels=LEVELS_ISO, $
            color=0, $
            charsize=2., $
            c_labels=0, $
            /OVERPLOT
endif
;
;overplot photosphere
if(n_elements(oplot_circle) ne 0) then begin
   loadct, 0
   plotcircle, oplot_circle, /oplot
endif

;overplot magnetic field lines
if(keyword_set(ralfven2)) then begin
   loadct, 0
   plotbfield, nlines=8, ralfven=ralfven2, /oplot
   plotbfield_neg, nlines=8, ralfven=ralfven2, /oplot
endif
;
;overplot the grid
if(keyword_set(oplot_grid)) then begin
   loadct, 0   
   for i=0, n_elements(arrx_2)-1 do begin
      oplot, [arrx_2(i),arrx_2(i)], [ylim(0),ylim(1)]
   endfor
   for i=0, n_elements(arry_2)-1 do begin
      oplot, [xlim(0),xlim(1)], [arry_2(i),arry_2(i)]
   endfor
endif
;
;overplot any given input array
if(n_elements(oplot2_x1) ne n_elements(oplot2_y1)) then begin
   print, 'error in contourplots_triple: nd of oplot2_x1 and oplot2_y1 not the same'
   stop
endif
if(n_elements(oplot2_x2) ne n_elements(oplot2_y2)) then begin
   print, 'error in contourplots_triple: nd of oplot2_x2 and oplot2_y2 not the same'
   stop
endif
if(n_elements(oplot2_x3) ne n_elements(oplot2_y3)) then begin
   print, 'error in contourplots_triple: nd of oplot2_x3 and oplot2_y3 not the same'
   stop
endif
if(n_elements(oplot2_x1) gt 1) then begin
   loadct, 0
   oplot, oplot2_x1, oplot2_y1, line=0, thick=2.
endif
if(n_elements(oplot2_x2) gt 1) then begin
   loadct, 0
   oplot, oplot2_x2, oplot2_y2, line=0, thick=2.
endif
if(n_elements(oplot2_x3) gt 1) then begin
   loadct, 0
   oplot, oplot2_x3, oplot2_y3, line=0, thick=2.
endif

;
;---------------------------third plot----------------------------------
;
;set the colortable
if (not keyword_set(ctable)) then begin
;default:
   loadct, 13, ncolors=ncolors, bottom=bottom
endif else begin
   if (keyword_set(ctbl_file)) then begin
      loadct, ctable, ncolors=ncolors, bottom=bottom, file=ctbl_file
   endif else begin
      loadct, ctable, ncolors=ncolors, bottom=bottom
   endelse
endelse
;
;set the region where to set plot
if(keyword_set(cb_top)) then begin
   !P.REGION=[0.6666, 0.0, .97, 0.78]
endif else begin
   !P.REGION=[0.54666, 0.0, 0.82, 1.]
endelse
;
;create contour plot
CONTOUR, ARR2D_3, ARRX_3, ARRY_3, $
         levels=LEVELS_CONTOUR, $
         c_color=c_colors, $
         /cell_fill, $
         title=titleStr3, $
         xtitle=xtitleStr, $
         ytitle=ytitleStr, $
         xrange=xlim, $
         yrange=ylim, $
         xstyle=1, ystyle=1 , $
         isotropic=ISOTROPIC, $
         irregular=IRREGULAR, $
         /noerase, $
         charsize=chsize_title, $       ;overall charsize (to set smaller title)
         xcharsize=chsize_axes, $     ;x-axis charsize
         ycharsize=chsize_axes        ;y-axis charsize
;
if(keyword_set(cb_top)) then begin
;set right position of colorbar
   xpos1_cb=!x.window(1)
endif
;
;overplot the isocontours, if specified
if (keyword_set(isoc)) then begin
   CONTOUR, arr2d_3, arrx_3, arry_3, $
            levels=LEVELS_ISO, $
            color=0, $
            charsize=2., $
            c_labels=0, $
            /OVERPLOT
endif
;
;overplot photosphere
if(n_elements(oplot_circle) ne 0) then begin
   loadct, 0
   plotcircle, oplot_circle, /oplot
endif

;overplot magnetic field lines
if(keyword_set(ralfven3)) then begin
   loadct, 0
   plotbfield, nlines=8, ralfven=ralfven3, /oplot
   plotbfield_neg, nlines=8, ralfven=ralfven3, /oplot
endif
;
;overplot the grid
if(keyword_set(oplot_grid)) then begin
   loadct, 0   
   for i=0, n_elements(arrx_3)-1 do begin
      oplot, [arrx_3(i),arrx_3(i)], [ylim(0),ylim(1)]
   endfor
   for i=0, n_elements(arry_3)-1 do begin
      oplot, [xlim(0),xlim(1)], [arry_3(i),arry_3(i)]
   endfor
endif
;
;overplot any given input array
if(n_elements(oplot3_x1) ne n_elements(oplot3_y1)) then begin
   print, 'error in contourplots_triple: nd of oplot3_x1 and oplot3_y1 not the same'
   stop
endif
if(n_elements(oplot3_x2) ne n_elements(oplot3_y2)) then begin
   print, 'error in contourplots_triple: nd of oplot3_x2 and oplot3_y2 not the same'
   stop
endif
if(n_elements(oplot3_x3) ne n_elements(oplot3_y3)) then begin
   print, 'error in contourplots_triple: nd of oplot3_x3 and oplot3_y3 not the same'
   stop
endif
if(n_elements(oplot3_x1) gt 1) then begin
   loadct, 0
   oplot, oplot3_x1, oplot3_y1, line=0, thick=2.
endif
if(n_elements(oplot3_x2) gt 1) then begin
   loadct, 0
   oplot, oplot3_x2, oplot3_y2, line=0, thick=2.
endif
if(n_elements(oplot3_x3) gt 1) then begin
   loadct, 0
   oplot, oplot3_x3, oplot3_y3, line=0, thick=2.
endif

;
;-----------------------------colorbar----------------------------------
;
;set the colortable
if (not keyword_set(ctable)) then begin
;default:
   loadct, 13, ncolors=ncolors, bottom=bottom
endif else begin
   if (keyword_set(ctbl_file)) then begin
      loadct, ctable, ncolors=ncolors, bottom=bottom, file=ctbl_file
   endif else begin
      loadct, ctable, ncolors=ncolors, bottom=bottom
   endelse
endelse
;
if(keyword_set(cb_top)) then begin
   !P.REGION=[0.0, 0.78, 1., 1.]
endif else begin
   !P.REGION=[0.82, 0.0, 1., 1.]
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
             charsize=chsize_axes, $
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
             charsize=chsize_axes, $
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
end
