PRO CONTOURPLOTS_SINGLE, ARR2D, ARRX, ARRY, $
                         NCOLORS, NLEVELS_ISO, BOTTOM, $
                         LEVELS_CONTOUR, LEVELS_ISO, C_COLORS, $
                         CB_INDX, CB_TICKMARK_NAME, $
                         xlim=xlim, $
                         ylim=ylim, $
                         titleStr=titleStr, $
                         xtitleStr=xtitleStr, $
                         ytitleStr=ytitleStr, $
                         ctitleStr=ctitleStr, $
                         isoc=ISOC, $
                         ctable=CTABLE, $
                         ctbl_file=CTBL_FILE, $
                         isotropic=ISOTROPIC, $
                         irregular=IRREGULAR, $
                         background1=BACKGROUND1, $
                         background2=BACKGROUND2, $
                         background3=BACKGROUND3, $
                         ogrid=ogrid, $
                         help=print_help
;
;+
; NAME:
;       contourplots_single
; PURPOSE:
;	This procedure creates a coloured contour-plot along with a colorbar
;	on the right, optimized for charsize=1.2 (ps-output with xsize=17.8
;	cm) (A&A paper style).
;       
; CALLING SEQUENCE:
;	contourplots_single, arr2d, arrx, arry, ncolors, nlevels_iso, $
;                            bottom, levels_contour, levels_iso, c_colors, $
;                            cb_inds, cb_tickmark_name
;
; INPUTS:
;       arr2d:  2-dimensional array to be plotted as contour
;       arrx:   x-coordinate of array, can be 2d-array when irregular is specified
;       arry:   y-coordinate of array, can be 2d-array when irregular is specified
;       all others:   color-specifiers as output from get_contour_colors.pro
;                     or get_contour_colors2.pro 
;
; KEYWORD PARAMETERS:
;       help:   Set this keyword (flag) to show the documentation of this
;               procedure
;       xlim:   Set this keyword to set the x-range of the plot
;       ylim:   Set this keyword to set the y-range of the plot
;       titleStr:  Set this keyword to set the title of the plot
;       xtitleStr: Set this keyword to set the x-title of the plot
;       ytitleStr: Set this keyword to set the y-title of the plot
;       ctitleStr: Set this keyword to set the colorbar-title of the plot
;       isoc:   Set this keyword (flag) to overplot iso-contours as straight lines
;       ctable: Set this keyword to set the color-table
;       ctbl_file: Set this keyword to read the color-table from a specified
;                  file
;       isotropic: Set this keyword (flag) to make an isotropic plot
;       irregular: Set this keyword (flag) if input data is irregulary spaced
;       background1: Set this keyword (flag) to plot a background (only within
;                    plot region) in grey
;       background2: Set this keyword (flag) to plot a background (only within
;                    plot region) in grey, but spare out the central star
;       background3: Set this keyword to a value, specifying the radius, up to
;                    which a background shall be plotted, and spare out the
;                    central star
;       ogrid:  Set this keyword to overplot the grid   
;
; OUTPUTS:
;
; EXAMPLE:
;       CONTOUR_PLOTS_SINGLE, arr2d, x, y, ncolors, nlevels_iso, bottom,
;                             levels_contour, levels_iso, c_colors, cb_indx,
;                             cb_tickmark_name, xlim=[0.,2.], ylim=[0.,2.], $
;                             ctable=13, /isotropic, /irregular
;-
;
;-------------------------OUTPUT IF HELP NEEDED-------------------------
;
if(keyword_set(print_help)) then begin
   doc_library, 'contourplots_single'
   return
endif
;
;-----------------------------------------------------------------------
;
;set x and y range
if (not keyword_set(xlim)) then xlim=[min(arrx),max(arrx)]
if (not keyword_set(ylim)) then ylim=[min(arry),max(arry)]
;
;make smaller x-margin for multi plots
xmarg_default=!x.margin
ymarg_default=!y.margin
;!x.margin=4.*!x.margin/5.
!x.margin=5.*!x.margin/5.
;
;
;set the region where to set plot
!P.REGION=[0.0, 0.0, .8, 1.]
;
;define background
if (keyword_set(background1)) then begin
   loadct, 0
   bg_arr=1.d0+fltarr(2,2)*0.d0
   contour, bg_arr, [xlim(0),xlim(1)], [ylim(0), ylim(1)], $
            levels=[1.], $
            c_colors=[150], $
            xstyle=1, ystyle=1, $
            charsize=1., $       ;overall charsize (to set smaller title)
            xcharsize=1.2, $     ;x-axis charsize
            ycharsize=1.2, $     ;y-axis charsize
            title=titleStr, $
            xtitle=xtitleStr, $
            ytitle=ytitleStr, $         
            xrange=xlim, $
            yrange=ylim, $         
            isotropic=ISOTROPIC, $
            /fill
endif      
;
if (keyword_set(background2)) then begin
   loadct, 0
   rmin=1.d0
   x1=max([xlim(0),xlim(1)])
   y1=max([ylim(0),ylim(1)])
   rmax=sqrt(x1^2 + y1^2)

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
            c_colors=[150], $
            xstyle=1, ystyle=1, $
            charsize=1., $       ;overall charsize (to set smaller title)
            xcharsize=1.2, $     ;x-axis charsize
            ycharsize=1.2, $     ;y-axis charsize
            title=titleStr, $
            xtitle=xtitleStr, $
            ytitle=ytitleStr, $         
            xrange=xlim, $
            yrange=ylim, $         
            isotropic=ISOTROPIC, $
            /fill
endif
;
if (keyword_set(background3)) then begin
   loadct, 0
   rmin=1.d0
   rmax=background3
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
            xstyle=1, ystyle=1, $
            charsize=1., $       ;overall charsize (to set smaller title)
            xcharsize=1.2, $     ;x-axis charsize
            ycharsize=1.2, $     ;y-axis charsize
            title=titleStr, $
            xtitle=xtitleStr, $
            ytitle=ytitleStr, $         
            xrange=xlim, $
            yrange=ylim, $         
            isotropic=ISOTROPIC, $
            /fill
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
;-----------------------------------------------------------------------
;
if(keyword_set(background1) or keyword_set(background2) or keyword_set(background3)) then begin
;create contour plot
   CONTOUR, ARR2D, ARRX, ARRY, $
            levels=LEVELS_CONTOUR, $
            c_color=c_colors, $
            /fill, $
            isotropic=ISOTROPIC, $
            irregular=IRREGULAR, $
            /overplot
endif else begin
   CONTOUR, ARR2D, ARRX, ARRY, $
            levels=LEVELS_CONTOUR, $
            c_color=c_colors, $
            /fill, $
            xstyle=1, ystyle=1 , $
            isotropic=ISOTROPIC, $
            irregular=IRREGULAR, $
            charsize=1., $       ;overall charsize (to set smaller title)
            xcharsize=1.2, $     ;x-axis charsize
            ycharsize=1.2, $     ;y-axis charsize
            title=titleStr, $
            xtitle=xtitleStr, $
            ytitle=ytitleStr, $         
            xrange=xlim, $
            yrange=ylim
endelse

;
;overplot the isocontours, if specified
if (keyword_set(isoc)) then begin
   CONTOUR, arr2d, arrx, arry, $
            levels=LEVELS_ISO, $
            color=0, $
            charsize=2., $
            c_labels=0, $
            /OVERPLOT
endif
;
;overplot the axes
xaxis=[MIN(ARRX), MAX(ARRX)]
yaxis=[MIN(ARRY), MAX(ARRY)]
origin=[0., 0.]
;OPLOT, xaxis, origin, line=0, color=0
;OPLOT, origin, yaxis, line=0, color=0

;
if(keyword_set(ogrid)) then begin
   nx=n_elements(arrx)
   nz=n_elements(arry)
   loadct, 0
   for i=0, nx-1 do begin
      oplot, [arrx(i),arrx(i)], $
             [arry(0),arry(nz-1)]
   endfor
   for i=0, nz-1 do begin
      oplot, [arrx(0),arrx(nx-1)], $
             [arry(i),arry(i)]
   endfor
endif
;
;-------------------------re-plot the axes------------------------------
;
if(keyword_set(background1) or keyword_set(background2) or keyword_set(background3)) then begin
;set zero color-table
   loadct, 0
;re-plot the x,y-axes, if background-keyword is set (contour with overplot
;will overplot the tick-marks)
;plot lower x-axis
   axis, xlim(0), ylim(0), xaxis=0, xrange=[xlim(0),xlim(1)], /xs, charsize=1.2
;plot upper x-axis (without labels)
   axis, xlim(0), ylim(1), xaxis=1, xrange=[xlim(0),xlim(1)], /xs, xtickformat='(A1)'
;plot left y-axis
   axis, xlim(0), ylim(0), yaxis=0, yrange=[ylim(0),ylim(1)], /ys, charsize=1.2
;plot right y-axis (without labels)
   axis, xlim(1), ylim(0), yaxis=1, yrange=[ylim(0),ylim(1)], /ys, ytickformat='(A1)'

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
endif
;
;-----------------------------colorbar----------------------------------
;
!P.REGION=[0.8, 0.0, 1., 1.]
;
region_width=!p.region(2)-!p.region(0)
pos=fltarr(4)
;x-arrangement depends on region
;pos(0)=!p.region(0)+0.1*region_width
;pos(2)=!p.region(2)-0.6*region_width
pos(0)=!p.region(0);+0.1*region_width
pos(2)=!p.region(2)-0.8*region_width
;
;y-arrangement depends on the position of the plot
pos(1)=!y.window(0)
pos(3)=!y.window(1)
;
;
;
if(keyword_set(ctitleStr)) then begin
   COLORBAR, NCOLORS=NCOLORS, /VERTICAL, $
             divisions=nlevels_iso-1, $
             tickvalues=CB_INDX, $
             tickname=CB_TICKMARK_NAME, $
             POSITION=pos, $
             title=ctitleStr, $
             /RIGHT, $
             charsize=1.2
endif else begin
   COLORBAR, NCOLORS=NCOLORS, /VERTICAL, $
             divisions=nlevels_iso-1, $
             tickvalues=CB_INDX, $
             tickname=CB_TICKMARK_NAME, $
             POSITION=pos, $
             /RIGHT, $
             charsize=1.2
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
