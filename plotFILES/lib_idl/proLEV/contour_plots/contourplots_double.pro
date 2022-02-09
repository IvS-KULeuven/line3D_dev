PRO CONTOURPLOTS_DOUBLE, arr2d_1, arr2d_2, arrx_1, arrx_2, arry_1, arry_2, $
                         NCOLORS, NLEVELS_ISO, BOTTOM, $
                         LEVELS_CONTOUR, LEVELS_ISO, C_COLORS, $
                         CB_INDX, CB_TICKMARK_NAME, $
                         xlim=xlim, $
                         ylim=ylim, $
                         titleStr1=titleStr1, $
                         titleStr2=titleStr2, $
                         xtitleStr=xtitleStr, $
                         ytitleStr=ytitleStr, $
                         ctitleStr=ctitleStr, $
                         isoc=ISOC, $
                         ctable=CTABLE, $
                         ctbl_file=CTBL_FILE, $
                         isotropic=ISOTROPIC, $
                         irregular=IRREGULAR, $
                         cb_top=cb_top, $
                         oplot_circle=oplot_circle, $
                         oplot_grid=oplot_grid, $
                         oplot1_x1=oplot1_x1, oplot1_y1=oplot1_y1, $
                         oplot1_x2=oplot1_x2, oplot1_y2=oplot1_y2, $
                         oplot2_x1=oplot2_x1, oplot2_y1=oplot2_y1, $
                         oplot2_x2=oplot2_x2, oplot2_y2=oplot2_y2, $
                         oplot1_style=oplot1_style, $
                         oplot2_style=oplot2_style, $ 
                         background1=background1
;
;+
; NAME:
;       contourplots_double
; PURPOSE:
;	This procedure creates two coloured contour-plot in a row, along with
;	a colorbar on the right, optimized for charsize=1.2 (ps-output with
;	xsize=17.8 cm) (A&A paper style).
;       
; CALLING SEQUENCE:
;	contourplots_double, arr2d_1, arr2d_2, arrx_1, arrx_2, arry_1, arry_2, $
;	                     ncolors, nlevels_iso, bottom, levels_contour, $
;	                     levels_iso, c_colors, cb_inds, cb_tickmark_name
;
; INPUTS:
;       arr2d_1:  2-dimensional array to be plotted as contour (left)
;       arr2d_2:  2-dimensional array to be plotted as contour (right)
;       arrx_1:   x-coordinate of array 1, can be 2d-array when irregular is specified
;       arrx_2:   x-coordinate of array 2, can be 2d-array when irregular is specified
;       arry_1:   y-coordinate of array 1, can be 2d-array when irregular is specified
;       arry_2:   y-coordinate of array 2, can be 2d-array when irregular is specified
;       all others:   color-specifiers as output from get_contour_colors.pro
;                     or get_contour_colors2.pro 
;
; KEYWORD PARAMETERS:
;       help:   Set this keyword (flag) to show the documentation of this
;               procedure
;       xlim:   Set this keyword to set the x-range of the plot
;       ylim:   Set this keyword to set the y-range of the plot
;       titleStr(1,2):  Set this keyword to set the title of the individual plots
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
;       oplot_circle: radius of a circle that shall be overplotted
;       oplot_grid:   Set this keyword (/flag) if x,y grid shall be
;                     overplotted
;       oplot1_x1, oplot1_y1: Specifies any x,y-array that shall be
;                             overplotted for plot '1' (left plot)
;       oplot1_x2, oplot1_y2: Specifies any x,y-array that shall be
;                             overplotted for plot '1' (left plot)
;       oplot2_x1, oplot2_y1: Specifies any x,y-array that shall be
;                             overplotted for plot '2' (right plot)
;       oplot2_x2, oplot2_y2: Specifies any x,y-array that shall be
;                             overplotted for plot '2' (right plot)
;       oplot1_style: Specifies the plot style of oplot1: 0 for black line, 1
;                     for red dots
;       oplot2_style: Specifies the plot style of oplot2: 0 for black line, 1
;                     for red dots
;       background1: Set this keyword (flag) to plot a background (only within
;                    plot region) in grey
;
; OUTPUTS:
;
; EXAMPLE:
;       CONTOUR_PLOTS_DOUBLE, arr2d1, arr2d2, x1, x2, y1, y2, ncolors, $
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
   doc_library, 'contourplots_double'
   return
endif
;
;---------------------------first plot----------------------------------
;
;make smaller x-margin for multi plots
xmarg_default=!x.margin
ymarg_default=!y.margin
!x.margin=4.*!x.margin/5.
;!x.margin=5.*!x.margin/5.

;set the region where to set plot
if(keyword_set(cb_top)) then begin
   !P.REGION=[0.03, 0.0, 0.5, 0.78]
endif else begin
   !P.REGION=[0.0, 0.0, .41, 1.]
endelse
;
;plot the background
if(keyword_set(background1)) then begin
   loadct, 0
   bg_arr=1.d0+fltarr(2,2)*0.d0
   contour, bg_arr, [xlim(0),xlim(1)], [ylim(0), ylim(1)], $
            levels=[1.], $
            c_colors=[150], $
            xstyle=1, ystyle=1, $
            charsize=1., $       ;overall charsize (to set smaller title)
            xcharsize=1.2, $     ;x-axis charsize
            ycharsize=1.2, $     ;y-axis charsize
            title=titleStr1, $
            xtitle=xtitleStr, $
            ytitle=ytitleStr, $         
            xrange=xlim, $
            yrange=ylim, $         
            isotropic=ISOTROPIC, $
            /cell_fill
endif
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
;create contour plot
if(keyword_set(background1)) then begin
   contour, arr2d_1, arrx_1, arry_1, $
            levels=levels_contour, $
            c_color=c_colors, $
            /cell_fill, $
            isotropic=isotrpoic, $
            irregular=irregular, $
            /overplot
endif else begin
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
         charsize=1., $       ;overall charsize (to set smaller title)
         xcharsize=1.2, $     ;x-axis charsize
         ycharsize=1.2        ;y-axis charsize
endelse
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
;overplot grid
if(keyword_set(oplot_grid)) then begin
   loadct, 0
   for i=0, n_elements(arrx_1)-1 do begin
      oplot, [arrx_1(i),arrx_1(i)], [min(arry_1), max(arry_1)], $
             color=0
   endfor
   for i=0, n_elements(arry_1)-1 do begin
      oplot, [min(arrx_1),max(arrx_1)], [arry_1(i), arry_1(i)], $
      color=0
   endfor
endif
;
;overplot any given input array
if(n_elements(oplot1_x1) ne n_elements(oplot1_y1)) then begin
   print, 'error in contourplots_double: nd of oplot1_x1 and oplot1_y1 not the same'
   stop
endif
if(n_elements(oplot1_x2) ne n_elements(oplot1_y2)) then begin
   print, 'error in contourplots_double: nd of oplot1_x2 and oplot1_y2 not the same'
   stop
endif
if(n_elements(oplot1_x1) gt 1) then begin
   loadct_lev, ci_black, ci_blue, ci_red, ci_cyan, ci_green, ci_magenta, ci_yellow, ci_white
   if(oplot1_style eq 1) then begin
      oplot, oplot1_x1, oplot1_y1, psym=1, thick=4., color=ci_red
   endif else begin
      oplot, oplot1_x1, oplot1_y1, line=0, color=ci_black, thick=2.
   endelse
endif
if(n_elements(oplot1_x2) gt 1) then begin
   loadct_lev, ci_black, ci_blue, ci_red, ci_cyan, ci_green, ci_magenta, ci_yellow, ci_white
   if(oplot1_style eq 1) then begin
      oplot, oplot1_x2, oplot1_y2, psym=1, thick=4., color=ci_red
   endif else begin
      oplot, oplot1_x2, oplot1_y2, line=0, color=ci_black, thick=2.
   endelse
endif
;
;re-plot the x,y-axes, if background-keyword is set (contour with overplot
;will overplot the tick-marks)
if(keyword_set(background1)) then begin
   loadct, 0
;plot lower x-axis
   axis, xlim(0), ylim(0), xaxis=0, xrange=[xlim(0),xlim(1)], /xs, charsize=1.2
;plot upper x-axis (without labels)
   axis, xlim(0), ylim(1), xaxis=1, xrange=[xlim(0),xlim(1)], /xs, xtickformat='(A1)'
;plot left y-axis
   axis, xlim(0), ylim(0), yaxis=0, yrange=[ylim(0),ylim(1)], /ys, charsize=1.2
;plot right y-axis (without labels)
   axis, xlim(1), ylim(0), yaxis=1, yrange=[ylim(0),ylim(1)], /ys, ytickformat='(A1)'
endif
;
;---------------------------second plot---------------------------------
;
;set the region where to set plot
if(keyword_set(cb_top)) then begin
   !P.REGION=[0.5, 0.0, .97, 0.78]
endif else begin
   !P.REGION=[0.41, 0.0, .82, 1.]
endelse
;
;plot the background
if(keyword_set(background1)) then begin
   loadct, 0
   bg_arr=1.d0+fltarr(2,2)*0.d0
   contour, bg_arr, [xlim(0),xlim(1)], [ylim(0), ylim(1)], $
            levels=[1.], $
            c_colors=[150], $
            xstyle=1, ystyle=1, $
            charsize=1., $       ;overall charsize (to set smaller title)
            xcharsize=1.2, $     ;x-axis charsize
            ycharsize=1.2, $     ;y-axis charsize
            title=titleStr2, $
            xtitle=xtitleStr, $
            ytitle=ytitleStr, $         
            xrange=xlim, $
            yrange=ylim, $         
            isotropic=ISOTROPIC, $
            /cell_fill, $
            /noerase
endif
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
;create contour plot
if(keyword_set(background1)) then begin
   contour, arr2d_2, arrx_2, arry_2, $
            levels=levels_contour, $
            c_color=c_colors, $
            /cell_fill, $
            isotropic=isotrpoic, $
            irregular=irregular, $
            /overplot
endif else begin
   CONTOUR, ARR2D_2, ARRX_2, ARRY_2, $
         levels=LEVELS_CONTOUR, $
         c_color=c_colors, $
         /cell_fill, $
         /noerase, $
         title=titleStr2, $
         xtitle=xtitleStr, $
         ytitle=ytitleStr, $
         xrange=xlim, $
         yrange=ylim, $
         xstyle=1, ystyle=1 , $
         isotropic=ISOTROPIC, $
         irregular=IRREGULAR, $
         charsize=1., $       ;overall charsize (to set smaller title)
         xcharsize=1.2, $     ;x-axis charsize
         ycharsize=1.2        ;y-axis charsize
endelse
;
if(keyword_set(cb_top)) then begin
;set right position of colorbar
   xpos1_cb=!x.window(1)
endif
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
;
;overplot grid
if(keyword_set(oplot_grid)) then begin
   loadct, 0
   for i=0, n_elements(arrx_2)-1 do begin
      oplot, [arrx_2(i),arrx_2(i)], [min(arry_2), max(arry_2)], $
             color=0
   endfor
   for i=0, n_elements(arry_2)-1 do begin
      oplot, [min(arrx_2),max(arrx_2)], [arry_2(i), arry_2(i)], $
      color=0
   endfor
endif
;
;overplot any given input array
if(n_elements(oplot2_x1) ne n_elements(oplot2_y1)) then begin
   print, 'error in contourplots_double: nd of oplot2_x1 and oplot2_y1 not the same'
   stop
endif
if(n_elements(oplot2_x2) ne n_elements(oplot2_y2)) then begin
   print, 'error in contourplots_double: nd of oplot2_x2 and oplot2_y2 not the same'
   stop
endif
if(n_elements(oplot2_x1) gt 1) then begin
   loadct_lev, ci_black, ci_blue, ci_red, ci_cyan, ci_green, ci_magenta, ci_yellow, ci_white
   if(oplot2_style eq 1) then begin
      oplot, oplot2_x1, oplot2_y1, psym=1, thick=4., color=ci_red
   endif else begin
      oplot, oplot2_x1, oplot2_y1, line=0, color=ci_black, thick=2.
   endelse
endif
if(n_elements(oplot2_x2) gt 1) then begin
   loadct_lev, ci_black, ci_blue, ci_red, ci_cyan, ci_green, ci_magenta, ci_yellow, ci_white
   if(oplot2_style eq 1) then begin
      oplot, oplot2_x2, oplot2_y2, psym=1, thick=4., color=ci_red
   endif else begin
      oplot, oplot2_x2, oplot2_y2, line=0, color=ci_black, thick=2.
   endelse
endif
;
;re-plot the x,y-axes, if background-keyword is set (contour with overplot
;will overplot the tick-marks)
if(keyword_set(background1)) then begin
   loadct, 0
;plot lower x-axis
   axis, xlim(0), ylim(0), xaxis=0, xrange=[xlim(0),xlim(1)], /xs, charsize=1.2
;plot upper x-axis (without labels)
   axis, xlim(0), ylim(1), xaxis=1, xrange=[xlim(0),xlim(1)], /xs, xtickformat='(A1)'
;plot left y-axis
   axis, xlim(0), ylim(0), yaxis=0, yrange=[ylim(0),ylim(1)], /ys, charsize=1.2
;plot right y-axis (without labels)
   axis, xlim(1), ylim(0), yaxis=1, yrange=[ylim(0),ylim(1)], /ys, ytickformat='(A1)'
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
end
