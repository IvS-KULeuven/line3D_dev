PRO CONTOURPLOTS_QUAD, arr2d_1, arr2d_2, arr2d_3, arr2d_4, $
                       arrx_1, arrx_2, arrx_3, arrx_4, $
                       arry_1, arry_2, arry_3, arry_4, $
                       NCOLORS, NLEVELS_ISO, BOTTOM, $
                       LEVELS_CONTOUR, LEVELS_ISO, C_COLORS, $
                       CB_INDX, CB_TICKMARK_NAME, $
                       xlim=xlim, $
                       ylim=ylim, $
                       titleStr1=titleStr1, $
                       titleStr2=titleStr2, $
                       titleStr3=titleStr3, $
                       titleStr4=titleStr4, $
                       xtitleStr=xtitleStr, $
                       ytitleStr=ytitleStr, $
                       ctitleStr=ctitleStr, $
                       isoc=ISOC, $
                       ctable=CTABLE, $
                       ctbl_file=CTBL_FILE, $
                       isotropic=ISOTROPIC, $
                       irregular=IRREGULAR, $
                       ralfven1=ralfven1, $
                       ralfven2=ralfven2, $
                       ralfven3=ralfven3, $
                       ralfven4=ralfven4
;
;+
; NAME:
;       contourplots_quad
; PURPOSE:
;	This procedure creates four coloured contour-plot in two rows and two
;	columns, additionally to a colorbar on the right, optimized for
;	charsize=1.2 (ps-output with xsize=17.8 cm) (A&A paper style).
;       
; CALLING SEQUENCE:
;	contourplots_quad, arr2d_1, arr2d_2, arr2d_3, arr2d_4, $
;                          arrx_1, arrx_2, arrx_3, arrx_4, $
;                          arry_1, arry_2, arry_3, arry_4, $
;	                   ncolors, nlevels_iso, bottom, levels_contour, $
;	                   levels_iso, c_colors, cb_inds, cb_tickmark_name
;
; INPUTS:
;       arr2d_1:  2-dimensional array to be plotted as contour (lower left)
;       arr2d_2:  2-dimensional array to be plotted as contour (lower right)
;       arr2d_3:  2-dimensional array to be plotted as contour (upper left)
;       arr2d_4:  2-dimensional array to be plotted as contour (upper right)
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
;       titleStr(1,2,3,4):  Set this keyword to set the title of the individual plots
;       xtitleStr: Set this keyword to set the x-title of the plot
;       ytitleStr: Set this keyword to set the y-title of the plot
;       ctitleStr: Set this keyword to set the colorbar-title of the plot
;       isoc:   Set this keyword (flag) to overplot iso-contours as straight lines
;       ctable: Set this keyword to set the color-table
;       ctbl_file: Set this keyword to read the color-table from a specified
;                  file
;       isotropic: Set this keyword (flag) to make an isotropic plot
;       irregular: Set this keyword (flag) if input data is irregulary spaced
;       ralfven(1,2,3,4): Set this keyword to overplot magnetic field lines
;
; OUTPUTS:
;
; EXAMPLE:
;       CONTOUR_PLOTS_QUAD, arr2d1, arr2d2, arr2d3, arr2d4, $
;                           x1, x2, x3, x4, y1, y2, y3, y4, ncolors, $
;                           nlevels_iso, bottom, levels_contour, levels_iso, $
;                           c_colors, cb_indx, cb_tickmark_name, $
;                           xlim=[0.,2.], ylim=[0.,2.], $ 
;                           ctable=13, /isotropic, /irregular
;
; NOTES:
;       so far, xlim, ylim, and all isocontours, etc. are the same for both plots
;-
;
;-------------------------OUTPUT IF HELP NEEDED-------------------------
;
if(keyword_set(print_help)) then begin
   doc_library, 'contourplots_quad'
   return
endif
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
!x.margin=4.*!x.margin/5.
;!x.margin=5.*!x.margin/5.
;
;set the region where to set plot
!P.REGION=[0.0, 0.0, .41, .5]
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
         charsize=1., $       ;overall charsize (to set smaller title)
         xcharsize=1.2, $     ;x-axis charsize
         ycharsize=1.2        ;y-axis charsize
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
;overplot the axes
xaxis=[MIN(ARRX_1), MAX(ARRX_1)]
yaxis=[MIN(ARRY_1), MAX(ARRY_1)]
origin=[0., 0.]
;OPLOT, xaxis, origin, line=0, color=0
;OPLOT, origin, yaxis, line=0, color=0
;
;overplot magnetic field lines
if(keyword_set(ralfven1)) then begin
   loadct, 0
   plotbfield, nlines=8, ralfven=ralfven1, /oplot
   plotbfield_neg, nlines=8, ralfven=ralfven1, /oplot
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
!P.REGION=[0.41, 0.0, .82, .5]
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
         charsize=1., $       ;overall charsize (to set smaller title)
         xcharsize=1.2, $     ;x-axis charsize
         ycharsize=1.2        ;y-axis charsize
;
;set upper position of colorbar
ypos0_cb=!y.window(0)
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
;overplot the axes
xaxis=[MIN(ARRX_2), MAX(ARRX_2)]
yaxis=[MIN(ARRY_2), MAX(ARRY_2)]
origin=[0., 0.]
;OPLOT, xaxis, origin, line=0, color=0
;OPLOT, origin, yaxis, line=0, color=0
;overplot magnetic field lines
;
if(keyword_set(ralfven2)) then begin
   loadct, 0
   plotbfield, nlines=8, ralfven=ralfven2, /oplot
   plotbfield_neg, nlines=8, ralfven=ralfven2, /oplot
endif
;
;----------------------------third plot---------------------------------
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
!P.REGION=[0., 0.5, .41, 1.]
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
         charsize=1., $       ;overall charsize (to set smaller title)
         xcharsize=1.2, $     ;x-axis charsize
         ycharsize=1.2        ;y-axis charsize
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
;overplot the axes
xaxis=[MIN(ARRX_3), MAX(ARRX_3)]
yaxis=[MIN(ARRY_3), MAX(ARRY_3)]
origin=[0., 0.]
;OPLOT, xaxis, origin, line=0, color=0
;OPLOT, origin, yaxis, line=0, color=0
;
;overplot magnetic field lines
if(keyword_set(ralfven3)) then begin
   loadct, 0
   plotbfield, nlines=8, ralfven=ralfven3, /oplot
   plotbfield_neg, nlines=8, ralfven=ralfven3, /oplot
endif
;
;---------------------------fourth plot---------------------------------
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
!P.REGION=[0.41, 0.5, .82, 1.]
;
;create contour plot
CONTOUR, ARR2D_4, ARRX_4, ARRY_4, $
         levels=LEVELS_CONTOUR, $
         c_color=c_colors, $
         /cell_fill, $
         title=titleStr4, $
         xtitle=xtitleStr, $
         ytitle=ytitleStr, $
         xrange=xlim, $
         yrange=ylim, $
         xstyle=1, ystyle=1 , $
         isotropic=ISOTROPIC, $
         irregular=IRREGULAR, $
         /noerase, $
         charsize=1., $       ;overall charsize (to set smaller title)
         xcharsize=1.2, $     ;x-axis charsize
         ycharsize=1.2        ;y-axis charsize
;
;set upper position of colorbar
ypos1_cb=!y.window(1)
;
;overplot the isocontours, if specified
if (keyword_set(isoc)) then begin
   CONTOUR, arr2d_4, arrx_4, arry_4, $
            levels=LEVELS_ISO, $
            color=0, $
            charsize=2., $
            c_labels=0, $
            /OVERPLOT
endif
;
;overplot the axes
xaxis=[MIN(ARRX_4), MAX(ARRX_4)]
yaxis=[MIN(ARRY_4), MAX(ARRY_4)]
origin=[0., 0.]
;OPLOT, xaxis, origin, line=0, color=0
;OPLOT, origin, yaxis, line=0, color=0

;overplot magnetic field lines
if(keyword_set(ralfven1)) then begin
   loadct, 0
   plotbfield, nlines=8, ralfven=ralfven4, /oplot
   plotbfield_neg, nlines=8, ralfven=ralfven4, /oplot
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
!P.REGION=[0.82, 0.0, 1., 1.]
;
region_width=!p.region(2)-!p.region(0)
region_height=!p.region(3)-!p.region(1)
pos=fltarr(4)
;
;x-arrangement depends on region
pos(0)=!p.region(0);+0.1*region_width
pos(2)=!p.region(2)-0.8*region_width
;y-arrangement depends on the position of the plot
pos(1)=ypos0_cb
pos(3)=ypos1_cb
;
if(keyword_set(ctitleStr)) then begin
   COLORBAR, NCOLORS=NCOLORS, $
             divisions=nlevels_iso-1, $
             tickvalues=CB_INDX, $
             tickname=CB_TICKMARK_NAME, $
             POSITION=pos, $
             title=ctitleStr, $
             charsize=1.2, $
             /VERTICAL, $
             /RIGHT
endif else begin
   COLORBAR, NCOLORS=NCOLORS, $
             divisions=nlevels_iso-1, $
             tickvalues=CB_INDX, $
             tickname=CB_TICKMARK_NAME, $
             POSITION=pos, $
             charsize=1.2, $
             /VERTICAL, $
             /RIGHT
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
