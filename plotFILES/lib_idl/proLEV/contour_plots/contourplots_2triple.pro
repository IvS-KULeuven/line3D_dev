PRO CONTOURPLOTS_2TRIPLE, arr2d_1, arr2d_2, arr2d_3, arr2d_4, arr2d_5, arr2d_6, $
                         arrx_1, arrx_2, arrx_3, arrx_4, arrx_5, arrx_6, $
                         arry_1, arry_2, arry_3, arry_4, arry_5, arry_6, $
                         NCOLORS, NLEVELS_ISO, BOTTOM, $
                         LEVELS_CONTOUR, LEVELS_ISO, C_COLORS, $
                         CB_INDX, CB_TICKMARK_NAME, $
                         xlim=xlim, $
                         ylim=ylim, $
                         titleStr1=titleStr1, $
                         titleStr2=titleStr2, $
                         titleStr3=titleStr3, $
                         titleStr4=titleStr4, $
                         titleStr5=titleStr5, $
                         titleStr6=titleStr6, $
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
                         ralfven4=ralfven4, $
                         ralfven5=ralfven5, $
                         ralfven6=ralfven6
;
;+
; NAME:
;       contourplots_2triple
; PURPOSE:
;	This procedure creates six coloured contour-plots with 2 rows, along with
;	a colorbar on the right (or on top), optimized for charsize=.8 (ps-output with
;	xsize=17.8 cm) for full page-width (A&A paper style).
;       
; CALLING SEQUENCE:
;	contourplots_triple, arr2d_1, arr2d_2, arr2d_3, arr2d_4, arr2d_5, $
;                            arr2d_6, $arrx_1, arrx_2, arrx_3, arrx_4, arrx_5,
;                            arrx_6, aarry_1, arry_2, arry_3, arry_4, arry_5,
;                            arry_6, ncolors, nlevels_iso, bottom, levels_contour, $
;	                     levels_iso, c_colors, cb_inds, cb_tickmark_name
;
; INPUTS:
;       layout:   arr2d_1-------arr2d_2-------arr2d_3
;                 arr2d_4-------arr2d_5-------arr2d_6
;            2-dimensional arrays to be plotted as contour
;       arrx_*:   x-coordinate of array 1-6, can be 2d-array when irregular is specified
;       arry_*:   y-coordinate of array 1-6, can be 2d-array when irregular is specified
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
;
; OUTPUTS:
;
; EXAMPLE:
;       CONTOUR_2TRIPLE, arr2d1, arr2d2, arr2d3, arr2d4, arr2d5, arr2d6, x1,
;                        x2, x3, x4, x5, x6, y1, y2, y3, y4, y5, y6,ncolors, $
;                        nlevels_iso, bottom, levels_contour, levels_iso, $
;                        c_colors, cb_indx, cb_tickmark_name, $
;                        xlim=[0.,2.], ylim=[0.,2.], $ 
;                        ctable=13, /isotropic, /irregular
;
; NOTES:
;       so far, xlim, ylim, and all isocontours, etc. are the same for both plots
;-
;
;-------------------------OUTPUT IF HELP NEEDED-------------------------
;
if(keyword_set(print_help)) then begin
   doc_library, 'contourplots_2triple'
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
   !P.REGION=[0.03, 0.39, 0.3233, 0.78]
endif else begin
   !P.REGION=[0.0, 0.5, .27333, 1.]
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
endif else begin
   ypos1_cb=!y.window(1)
endelse
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
if(keyword_set(cb_top)) then begin
   !P.REGION=[0.3533, 0.39, 0.6466, 0.78]
endif else begin
   !P.REGION=[0.27333, 0.5, 0.54666, 1.]
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

;overplot magnetic field lines
if(keyword_set(ralfven2)) then begin
   loadct, 0
   plotbfield, nlines=8, ralfven=ralfven2, /oplot
   plotbfield_neg, nlines=8, ralfven=ralfven2, /oplot
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
   !P.REGION=[0.6766, 0.39, .97, 0.78]
endif else begin
   !P.REGION=[0.54666, 0.5, 0.82, 1.]
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
if(keyword_set(cb_top)) then begin
   !P.REGION=[0.03, 0.0, 0.3233, 0.39]
endif else begin
   !P.REGION=[0.0, 0.0, .27333, .5]
endelse
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
         charsize=chsize_title, $     ;overall charsize (to set smaller title)
         xcharsize=chsize_axes, $     ;x-axis charsize
         ycharsize=chsize_axes        ;y-axis charsize
;
if(not keyword_set(cb_top)) then begin
;set right position of colorbar
   ypos0_cb=!y.window(0)
endif
;
;overplot the isocontours, if specified
if (keyword_set(isoc)) then begin
   CONTOUR, arr4d_1, arrx_4, arry_4, $
            levels=LEVELS_ISO, $
            color=0, $
            charsize=2., $
            c_labels=0, $
            /OVERPLOT
endif
;
;overplot magnetic field lines
if(keyword_set(ralfven4)) then begin
   loadct, 0
   plotbfield, nlines=8, ralfven=ralfven4, /oplot
   plotbfield_neg, nlines=8, ralfven=ralfven4, /oplot
endif
;
;--------------------------fifth plot-----------------------------------
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
   !P.REGION=[0.3533, 0.0, 0.6466, 0.39]
endif else begin
   !P.REGION=[0.27333, 0.0, 0.54666, .5]
endelse
;
;create contour plot
CONTOUR, ARR2D_5, ARRX_5, ARRY_5, $
         levels=LEVELS_CONTOUR, $
         c_color=c_colors, $
         /cell_fill, $
         title=titleStr5, $
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
   CONTOUR, arr2d_5, arrx_5, arry_5, $
            levels=LEVELS_ISO, $
            color=0, $
            charsize=2., $
            c_labels=0, $
            /OVERPLOT
endif

;overplot magnetic field lines
if(keyword_set(ralfven5)) then begin
   loadct, 0
   plotbfield, nlines=8, ralfven=ralfven5, /oplot
   plotbfield_neg, nlines=8, ralfven=ralfven5, /oplot
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
   !P.REGION=[0.6766, 0.0, .97, 0.39]
endif else begin
   !P.REGION=[0.54666, 0.0, 0.82, .5]
endelse
;
;create contour plot
CONTOUR, ARR2D_6, ARRX_6, ARRY_6, $
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
   CONTOUR, arr2d_6, arrx_6, arry_6, $
            levels=LEVELS_ISO, $
            color=0, $
            charsize=2., $
            c_labels=0, $
            /OVERPLOT
endif

;overplot magnetic field lines
if(keyword_set(ralfven6)) then begin
   loadct, 0
   plotbfield, nlines=8, ralfven=ralfven6, /oplot
   plotbfield_neg, nlines=8, ralfven=ralfven6, /oplot
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
   pos(1)=ypos0_cb
   pos(3)=ypos1_cb
endelse
print, !y.window
print, pos
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
