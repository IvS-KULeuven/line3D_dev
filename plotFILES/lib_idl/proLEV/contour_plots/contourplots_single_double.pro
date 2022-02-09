;
;-----------------------------------------------------------------------
;-----------------------------------------------------------------------
;-----------------------------------------------------------------------
;
pro set_xposition, marginl, marginr, eps, indx, x0, x1
;
;calculates x-position of a plot for 3 plots in a row
;
;input: marginr, marginl:  right and left margin in device coordinates
;       eps: distance between the three plots
;       indx: index for the plot (from 1 to 3)
;
;output: x0, x1: x-positions of the plot in device coordinates
;
;
;width of each plot:
delta = (1.d0-marginl-marginr-2.d0*eps)/3.d0
;
if(indx eq 1) then begin
   x0=marginl
   x1=x0+delta
endif
;
if(indx eq 2) then begin
   x0=marginl+delta+eps
   x1=x0+delta
endif
;
if(indx eq 3) then begin
   x0=marginl+delta+eps+delta+eps
   x1=x0+delta
endif
;
;print, delta
;stop

end
;
;-----------------------------------------------------------------------
;-----------------------------------------------------------------------
;-----------------------------------------------------------------------
;
pro contourplots_single_double, arr2d_1, arr2d_2, arr2d_3, $
                                arrx_1, arrx_2, arrx_3, $
                                arry_1, arry_2, arry_3, $
                                ncolors1, nlevels_iso1, bottom1, $
                                levels_contour1, levels_iso1, c_colors1, $
                                cb_indx1, cb_tickmark_name1, $
                                ncolors23, nlevels_iso23, bottom23, $
                                levels_contour23, levels_iso23, c_colors23, $
                                cb_indx23, cb_tickmark_name23, $
                                xlim=xlim, $
                                ylim=ylim, $
                                titlestr1=titlestr1, $
                                titlestr2=titlestr2, $
                                titlestr3=titlestr3, $
                                xtitlestr=xtitlestr, $
                                ytitlestr=ytitlestr, $
                                ctitlestr1=ctitlestr1, $
                                ctitlestr23=ctitlestr23, $
                                isoc=isoc, $
                                ctable=ctable, $
                                ctbl_file=ctbl_file, $
                                isotropic=isotropic, $
                                irregular=irregular, $
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
                                oplot3_x3=oplot3_x3, oplot3_y3=oplot3_y3, $
                                coplot1_arr=coplot1_arr, coplot1_x1=coplot1_x1, xoplot1_y1=coplot1_y1, $
                                coplot1_levels=coplot1_levels, coplot1_thick=coplot1_thick, coplot1_linestyle=coplot1_linestyle, $
                                coplot2_arr=coplot2_arr, coplot2_x1=coplot2_x1, coplot2_y1=coplot2_y1, $
                                coplot2_levels=coplot2_levels, coplot2_thick=coplot2_thick, coplot2_linestyle=coplot2_linestyle, $
                                coplot3_arr=coplot3_arr, coplot3_x1=coplot3_x1, coplot3_y1=coplot3_y1, $
                                coplot3_levels=coplot3_levels, coplot3_thick=coplot3_thick, coplot3_linestyle=coplot3_linestyle, $
                                background1_level3=background1_level3, $
                                background2_level3=background2_level3, $
                                background3_level3=background3_level3

;
;+
; NAME:
;       contourplots_single_double
; PURPOSE:
;	This procedure creates three coloured contour-plots in a row, along with
;	two colorbars top the top where the first colourbar holds for
;	contour-plot 1, and the second one for contour-plot 2 and 3
;       
; CALLING SEQUENCE:
;	contourplots_single_double, arr2d_1, arr2d_2, arr2d_3, arrx_1, arrx_2, $
;                                   arrrx_3, arry_1, arry_2, arry_3, $
;	                            ncolors1, nlevels_iso1, bottom1, levels_contour1, $
;	                            levels_iso1, c_colors1, cb_indx1, $
;	                            cb_tickmark_name1, $
;	                            ncolors23, nlevels_iso23, bottom23, levels_contour23, $
;	                            levels_iso23, c_colors23, cb_indx23, cb_tickmark_name23
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
;       ctitleStr1: Set this keyword to set the colorbar-title of plot 1
;       ctitleStr23: Set this keyword to set the colorbar-title of plot 2 & 3
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
;       coplot1, coplot1_x1, coplot1_y1: Specifies any 2d and x,y-array that shall be
;                             overplotted as contours for plot '1' (left plot)
;       coplot_levels1, coplot_thick1 coplot_linestyle1 style for overplotting contours
;       coplot2, coplot2_x1, coplot2_y1: Specifies any 2d and x,y-array that shall be
;                             overplotted as contours for plot '2' (middle plot)
;       coplot_levels2, coplot_thick2, coplot_linestyle2: style for overplotting contours
;       coplot3, coplot3_x1, coplot3_y1: Specifies any 2d and x,y-array that shall be
;                             overplotted as contours for plot '3' (right plot)
;       coplot_levels3, coplot_thick3, coplot_linestyle3: style for overplotting contours
;
;        background1_level3: background for plot 1 on level 3
;        background2_level3: background for plot 2 on level 3
;        background3_level3: background for plot 3 on level 3
;
; OUTPUTS:
;
; EXAMPLE:
;       contourplots_single_double, arr2d1, arr2d2, arr2d3, x1, x2, x3, y1, y2, y3, $
;                                   ncolors1, nlevels_iso1, bottom1, levels_contour1, $
;                                   levels_iso1, c_colors1, cb_indx1, cb_tickmark_name1, $
;                                   ncolors23, nlevels_iso23, bottom23, levels_contour23, $
;                                   levels_iso23, c_colors23, cb_indx23, cb_tickmark_name23, $
;                                   xlim=[0.,2.], ylim=[0.,2.], $ 
;                                   ctable=13, /isotropic, /irregular
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
if(not keyword_set(cb_top)) then begin
   print, 'colorbar to the right makes no sense in this routine'
   print, 'set keyword cb_top'
   stop
endif
;
;-----------------------------------------------------------------------
;
chsize_title=.8
chsize_axes=.8
;
;---------------------------first plot----------------------------------
;
;make smaller x-margin for multi plots
xmarg_default=!x.margin
ymarg_default=!y.margin
!x.margin=4.*!x.margin/5.
;
;make smaller x-margin for multi plots
marginl=0.06d0
marginr=0.04d0
eps=0.06d0
;
set_xposition, marginl, marginr, eps, 1, x0, x1
y0=0.1
y1=0.8
;
;create grey background with inner radius 1 and outer radius background1_level3
if(keyword_set(background1_level3)) then begin
   loadct, 0
   rmin=1.d0
   rmax=background1_level3
   ntheta=60
   bg_arr=1.d0+fltarr(2,ntheta)*0.d0
   r_arr=[rmin,rmax]
   bg_xarr=fltarr(2,ntheta)
   bg_yarr=fltarr(2,ntheta)
   for i=0, 1 do begin
      for j=0, ntheta-1 do begin
         bg_xarr(i,j)=r_arr(i)*sin(2.d0*!pi*j/(ntheta-1))
         bg_yarr(i,j)=r_arr(i)*cos(2.d0*!pi*j/(ntheta-1))
      endfor
   endfor
   contour, bg_arr, bg_xarr, bg_yarr, $
            levels=[1.], $
            c_colors=[200], $
            charsize=chsize_title, $       ;overall charsize (to set smaller title)
            xcharsize=chsize_axes, $     ;x-axis charsize
            ycharsize=chsize_axes, $     ;y-axis charsize
            title=titleStr1, $
            xtitle=xtitleStr, $
            ytitle=ytitleStr, $         
            xrange=xlim, $
            yrange=ylim, $
            xstyle=1, ystyle=1 , $
            isotropic=ISOTROPIC, $
            /cell_fill, $
            position=[x0,y0,x1,y1]
endif
;
;set the colortable
if (not keyword_set(ctable)) then begin
;default:
   loadct, 13, ncolors=ncolors1, bottom=bottom1
endif else begin
   if (keyword_set(ctbl_file)) then begin
      loadct, ctable, ncolors=ncolors1, bottom=bottom1, file=ctbl_file
   endif else begin
      loadct, ctable, ncolors=ncolors1, bottom=bottom1
   endelse
endelse
;
;-----------------------------------------------------------------------
;
;set the region where to set plot
;!p.region=[x0, y0, x1, y1]
;
;create contour plot
if(keyword_set(background1_level1)) then begin
   contour, arr2d_1, arrx_1, arry_1, $
            levels=levels_contour1, $
            c_color=c_colors1, $
            isotropic=ISOTROPIC, $
            irregular=IRREGULAR, $
            /overplot
endif else begin
   contour, arr2d_1, arrx_1, arry_1, $
            levels=levels_contour1, $
            c_color=c_colors1, $
            /cell_fill, $
            title=titlestr1, $
            xtitle=xtitlestr, $
            ytitle=ytitlestr, $
            xrange=xlim, $
            yrange=ylim, $
            xstyle=1, ystyle=1 , $
            position=[x0,y0,x1,y1], $
            isotropic=isotropic, $
            irregular=irregular, $
            charsize=chsize_title, $     ;overall charsize (to set smaller title)
            xcharsize=chsize_axes, $     ;x-axis charsize
            ycharsize=chsize_axes        ;y-axis charsize
endelse
;
;set left and right position of colorbar
xpos0_cb=!x.window(0)
xpos1_cb=!x.window(1)
;
;overplot the isocontours, if specified
if(keyword_set(isoc)) then begin
   contour, arr2d_1, arrx_1, arry_1, $
            levels=levels_iso1, $
            color=0, $
            charsize=2., $
            c_labels=0, $
            /overplot
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
;overplot contours
if(n_elements(coplot1_arr) gt 1) then begin
    contour, coplot1_arr, coplot1_x1, coplot1_y1, $
             levels=coplot1_levels, $
             c_linestyle=coplot1_linestyle, $
             c_thick=coplot1_thick, $
             /overplot
 endif
;
;create the colorbar
;set the colortable
if (not keyword_set(ctable)) then begin
;default:
   loadct, 13, ncolors=ncolors1, bottom=bottom1
endif else begin
   if(keyword_set(ctbl_file)) then begin
      loadct, ctable, ncolors=ncolors1, bottom=bottom1, file=ctbl_file
   endif else begin
      loadct, ctable, ncolors=ncolors1, bottom=bottom1
   endelse
endelse
;
region_width=x1-x0
region_height=1.-y1
pos=fltarr(4)
;
if(keyword_set(cb_top)) then begin
;x-arrangement depends on the position of the plot
   pos(0)=xpos0_cb
   pos(2)=xpos1_cb
;y-arrangement depends on region
   pos(1)=y1+0.15*region_height
   pos(3)=1.-0.55*region_height
endif
;
if(keyword_set(cb_top)) then begin
   vert=0
   horiz=1
   right=0
   top=1
endif
;
;
if(keyword_set(ctitleStr1)) then begin
   colorbar, ncolors=ncolors1, $
             divisions=nlevels_iso1-1, $
             tickvalues=cb_indx1, $
             tickname=cb_tickmark_name1, $
             position=pos, $
             title=ctitlestr1, $
             charsize=chsize_axes, $
             vertical=vert, $
             right=right, $
             top=top, $
             horizontal=horiz
endif else begin
   colorbar, ncolors=ncolors1, $
             divisions=nlevels_iso1-1, $
             tickvalues=cb_indx1, $
             tickname=cb_tickmark_name1, $
             position=pos, $
             charsize=chsize_axes, $
             vertical=vert, $
             right=right, $
             top=top, $
             horizontal=horiz
endelse
;
;---------------------------second plot---------------------------------
;
;set the region where to set plot
set_xposition, marginl, marginr, eps, 2, x0, x1
;
;create grey background with inner radius 1 and outer radius background2_level3
if (keyword_set(background2_level3)) then begin
   loadct, 0
   rmin=1.d0
   rmax=background2_level3
   ntheta=60
   bg_arr=1.d0+fltarr(2,ntheta)*0.d0
   r_arr=[rmin,rmax]
   bg_xarr=fltarr(2,ntheta)
   bg_yarr=fltarr(2,ntheta)
   for i=0, 1 do begin
      for j=0, ntheta-1 do begin
         bg_xarr(i,j)=r_arr(i)*sin(2.d0*!pi*j/(ntheta-1))
         bg_yarr(i,j)=r_arr(i)*cos(2.d0*!pi*j/(ntheta-1))
      endfor
   endfor
   contour, bg_arr, bg_xarr, bg_yarr, $
            levels=[1.], $
            c_colors=[200], $
            charsize=chsize_title, $       ;overall charsize (to set smaller title)
            xcharsize=chsize_axes, $     ;x-axis charsize
            ycharsize=chsize_axes, $     ;y-axis charsize
            title=titleStr2, $
            xtitle=xtitleStr, $
            ytitle=ytitleStr, $         
            xrange=xlim, $
            yrange=ylim, $
            xstyle=1, ystyle=1 , $
            isotropic=ISOTROPIC, $
            /cell_fill, $
            position=[x0,y0,x1,y1], $
            /noerase
endif
;
;set the colortable
if(not keyword_set(ctable)) then begin
;default:
   loadct, 13, ncolors=ncolors23, bottom=bottom23
endif else begin
   if (keyword_set(ctbl_file)) then begin
      loadct, ctable, ncolors=ncolors23, bottom=bottom23, file=ctbl_file
   endif else begin
      loadct, ctable, ncolors=ncolors23, bottom=bottom23
   endelse
endelse
;
;
;create contour plot
if(keyword_set(background2_level3)) then begin
   contour, arr2d_2, arrx_2, arry_2, $
            levels=levels_contour23, $
            c_color=c_colors23, $
            /cell_fill, $
            isotropic=ISOTROPIC, $
            irregular=IRREGULAR, $
            /overplot
endif else begin
   contour, arr2d_2, arrx_2, arry_2, $
            levels=levels_contour23, $
            c_color=c_colors23, $
            /cell_fill, $
            title=titlestr2, $
            xtitle=xtitlestr, $
            ytitle=ytitlestr, $
            xrange=xlim, $
            yrange=ylim, $
            xstyle=1, ystyle=1 , $
            position=[x0,y0,x1,y1], $
            isotropic=isotropic, $
            irregular=irregular, $
            /noerase, $
            charsize=chsize_title, $       ;overall charsize (to set smaller title)
            xcharsize=chsize_axes, $     ;x-axis charsize
            ycharsize=chsize_axes        ;y-axis charsize
endelse
;
;set left position of colorbar
xpos0_cb_test=x0
xpos0_cb=!x.window(0)
;
;overplot the isocontours, if specified
if(keyword_set(isoc)) then begin
   contour, arr2d_2, arrx_2, arry_2, $
            levels=levels_iso23, $
            color=0, $
            charsize=2., $
            c_labels=0, $
            /overplot
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
;overplot contours
if(n_elements(coplot2_arr) gt 1) then begin
    contour, coplot2_arr, coplot2_x1, coplot2_y1, $
             levels=coplot2_levels, $
             c_linestyle=coplot2_linestyle, $
             c_thick=coplot2_thick, $
            /overplot
endif
;
;---------------------------third plot----------------------------------
;
;set the region where to set plot
set_xposition, marginl, marginr, eps, 3, x0, x1
;
;create grey background with inner radius 1 and outer radius background3_level3
if (keyword_set(background3_level3)) then begin
   loadct, 0
   rmin=1.d0
   rmax=background3_level3
   ntheta=60
   bg_arr=1.d0+fltarr(2,ntheta)*0.d0
   r_arr=[rmin,rmax]
   bg_xarr=fltarr(2,ntheta)
   bg_yarr=fltarr(2,ntheta)
   for i=0, 1 do begin
      for j=0, ntheta-1 do begin
         bg_xarr(i,j)=r_arr(i)*sin(2.d0*!pi*j/(ntheta-1))
         bg_yarr(i,j)=r_arr(i)*cos(2.d0*!pi*j/(ntheta-1))
      endfor
   endfor
   contour, bg_arr, bg_xarr, bg_yarr, $
            levels=[1.], $
            c_colors=[200], $
            charsize=chsize_title, $       ;overall charsize (to set smaller title)
            xcharsize=chsize_axes, $     ;x-axis charsize
            ycharsize=chsize_axes, $     ;y-axis charsize
            title=titleStr3, $
            xtitle=xtitleStr, $
            ytitle=ytitleStr, $         
            xrange=xlim, $
            yrange=ylim, $
            xstyle=1, ystyle=1 , $
            isotropic=ISOTROPIC, $
            /cell_fill, $
            position=[x0,y0,x1,y1], $
            /noerase
endif
;
;set the colortable
if (not keyword_set(ctable)) then begin
;default:
   loadct, 13, ncolors=ncolors23, bottom=bottom23
endif else begin
   if (keyword_set(ctbl_file)) then begin
      loadct, ctable, ncolors=ncolors23, bottom=bottom23, file=ctbl_file
   endif else begin
      loadct, ctable, ncolors=ncolors23, bottom=bottom23
   endelse
endelse
;
;create contour plot
if(keyword_set(background3_level3)) then begin
   contour, arr2d_3, arrx_3, arry_3, $
            levels=levels_contour23, $
            c_color=c_colors23, $
            /cell_fill, $
            isotropic=isotropic, $
            irregular=irregular, $
            /overplot
endif else begin
   contour, arr2d_3, arrx_3, arry_3, $
            levels=levels_contour23, $
            c_color=c_colors23, $
            /cell_fill, $
            title=titlestr3, $
            xtitle=xtitlestr, $
            ytitle=ytitlestr, $
            xrange=xlim, $
            yrange=ylim, $
            xstyle=1, ystyle=1 , $
            isotropic=isotropic, $
            position=[x0,y0,x1,y1], $
            irregular=irregular, $
            /noerase, $
            charsize=chsize_title, $     ;overall charsize (to set smaller title)
            xcharsize=chsize_axes, $     ;x-axis charsize
            ycharsize=chsize_axes        ;y-axis charsize
endelse
;
if(keyword_set(cb_top)) then begin
;set right position of colorbar
   xpos1_cb=!x.window(1)
endif
;
;overplot the isocontours, if specified
if(keyword_set(isoc)) then begin
   contour, arr2d_3, arrx_3, arry_3, $
            levels=levels_iso23, $
            color=0, $
            charsize=2., $
            c_labels=0, $
            /overplot
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
;overplot contours
if(n_elements(coplot3_arr) gt 1) then begin
    contour, coplot3_arr, coplot3_x1, coplot3_y1, $
             levels=coplot3_levels, $
             c_linestyle=coplot3_linestyle, $
             c_thick=coplot3_thick, $
            /overplot
endif
;
;-----------------------------colorbar----------------------------------
;
;set the colortable
if (not keyword_set(ctable)) then begin
;default:
   loadct, 13, ncolors=ncolors23, bottom=bottom23
endif else begin
   if (keyword_set(ctbl_file)) then begin
      loadct, ctable, ncolors=ncolors23, bottom=bottom23, file=ctbl_file
   endif else begin
      loadct, ctable, ncolors=ncolors23, bottom=bottom23
   endelse
endelse
;


region_width=x1-xpos0_cb_test
region_height=1.-y1
pos=fltarr(4)
;
if(keyword_set(cb_top)) then begin
;x-arrangement depends on the position of the plot
   pos(0)=xpos0_cb
   pos(2)=xpos1_cb
;y-arrangement depends on region
   pos(1)=y1+0.15*region_height
   pos(3)=1.-0.55*region_height
endif
;
if(keyword_set(cb_top)) then begin
   vert=0
   horiz=1
   right=0
   top=1
endif
;
;
if(keyword_set(ctitleStr23)) then begin
;/VERTICAL, /RIGHT

   colorbar, ncolors=ncolors23, $
             divisions=nlevels_iso23-1, $
             tickvalues=cb_indx23, $
             tickname=cb_tickmark_name23, $
             position=pos, $
             title=ctitlestr23, $
             charsize=chsize_axes, $
             vertical=vert, $
             right=right, $
             top=top, $
             horizontal=horiz
endif else begin
   colorbar, ncolors=ncolors23, $
             divisions=nlevels_iso23-1, $
             tickvalues=cb_indx23, $
             tickname=cb_tickmark_name23, $
             position=pos, $
             charsize=chsize_axes, $
             vertical=vert, $
             right=right, $
             top=top, $
             horizontal=horiz
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
