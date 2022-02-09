;
;-----------------------------------------------------------------------
;-----------------------------------------------------------------------
;-----------------------------------------------------------------------
;
pro set_position, x0p, x1p, y0p, y1p, x0c, x1c, y0c, y1c
;
nleft=!p.multi(0)
ncols=!p.multi(1)
nrows=!p.multi(2)
nplots=ncols*nrows

if(nleft eq 0.) then begin
;standard treatment if !p.multi=0
   x0p=0.
   x1p=0.8
   y0p=0.
   y1p=1.
   x0c=0.8
   x1c=1.
   y0c=0.
   y1c=1.
   return
endif else begin
;
;width and height of the plot area
   dx=1.d0/ncols
   dy=1.d0/nrows
;
;width and height of the plot itself
   dx_plot = dx*0.8d0
   dy_plot = dy*1.d0
;
;find row and column index of multiplot
   for i=1, ncols do begin
      for j=1, nrows do begin
         if(nleft eq 1-i+ncols*j) then begin
            indx_col = i
            indx_row = j
            goto, loopend
         endif
      endfor
   endfor
   loopend: 
;plot position
   x0p = (indx_col-1)*dx
   x1p = x0p+dx_plot
   y0p = (indx_row-1)*dy
   y1p = y0p+dy_plot
;
;colorbar position
   x0c = x1p
   y0c = y0p
   x1c = x0p+dx
   y1c = y0p+dy
;
endelse
;
end
;
;-----------------------------------------------------------------------
;-----------------------------------------------------------------------
;-----------------------------------------------------------------------
;
pro contourplots_single2, arr2d, arrx, arry, $
                         ncolors, nlevels_iso, bottom, $
                         levels_contour, levels_iso, c_colors, $
                         cb_indx, cb_tickmark_name, $
                         xlim=xlim, $
                         ylim=ylim, $
                         titlestr=titlestr, $
                         xtitlestr=xtitlestr, $
                         ytitlestr=ytitlestr, $
                         ctitlestr=ctitlestr, $
                         isoc=isoc, $
                         ctable=ctable, $
                         ctbl_file=ctbl_file, $
                         isotropic=isotropic, $
                         irregular=irregular, $
                         background1=background1, $
                         background2=background2, $
                         background3=background3
;
;+
; name:
;       contourplots_single2
; purpose:
;	this procedure creates a coloured contour-plot along with a colorbar
;	on the right, optimized for charsize=1.2 (ps-output with xsize=17.8
;	cm) (a&a paper style).
;       this procedure can be used if !p.multi environment variable is set
;       
; calling sequence:
;	contourplots_single2, arr2d, arrx, arry, ncolors, nlevels_iso, $
;                            bottom, levels_contour, levels_iso, c_colors, $
;                            cb_inds, cb_tickmark_name
;
; inputs:
;       arr2d:  2-dimensional array to be plotted as contour
;       arrx:   x-coordinate of array, can be 2d-array when irregular is specified
;       arry:   y-coordinate of array, can be 2d-array when irregular is specified
;       all others:   color-specifiers as output from get_contour_colors.pro
;                     or get_contour_colors2.pro 
;
; keyword parameters:
;       help:   set this keyword (flag) to show the documentation of this
;               procedure
;       xlim:   set this keyword to set the x-range of the plot
;       ylim:   set this keyword to set the y-range of the plot
;       titlestr:  set this keyword to set the title of the plot
;       xtitlestr: set this keyword to set the x-title of the plot
;       ytitlestr: set this keyword to set the y-title of the plot
;       ctitlestr: set this keyword to set the colorbar-title of the plot
;       isoc:   set this keyword (flag) to overplot iso-contours as straight lines
;       ctable: set this keyword to set the color-table
;       ctbl_file: set this keyword to read the color-table from a specified
;                  file
;       isotropic: set this keyword (flag) to make an isotropic plot
;       irregular: set this keyword (flag) if input data is irregulary spaced
;       background1: set this keyword (flag) to plot a background (only within
;                    plot region) in grey
;       background2: set this keyword (flag) to plot a background (only within
;                    plot region) in grey, but spare out the central star
;       background3: set this keyword to a value, specifying the radius, up to
;                    which a background shall be plotted, and spare out the
;                    central star
;
; outputs:
;
; example:
;       contour_plots_single2, arr2d, x, y, ncolors, nlevels_iso, bottom,
;                             levels_contour, levels_iso, c_colors, cb_indx,
;                             cb_tickmark_name, xlim=[0.,2.], ylim=[0.,2.], $
;                             ctable=13, /isotropic, /irregular
;-
;
;-------------------------output if help needed-------------------------
;
if(keyword_set(print_help)) then begin
   doc_library, 'contourplots_single2'
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
;set the region where to set plot
;!p.region=[0.0, 0.0, .8, 1.]
set_position, x0p, x1p, y0p, y1p, x0c, x1c, y0c, y1c
!p.region=[x0p,y0p,x1p,y1p]
;
;define background
if (keyword_set(background1)) then begin
   loadct, 0
   bg_arr=1.d0+fltarr(2,2)*0.d0
   contour, bg_arr, [xlim(0),xlim(1)], [ylim(0), ylim(1)], $
            levels=[1.], $
            c_colors=[150], $
            xstyle=1, ystyle=1, $
            charsize=1.2, $       ;overall charsize (to set smaller title)
            xcharsize=0.0001, $     ;x-axis charsize (will be set below)
            ycharsize=0.0001, $     ;y-axis charsize (will be set below)
            title=titlestr, $
            xtitle=xtitlestr, $
            ytitle=ytitlestr, $         
            xrange=xlim, $
            yrange=ylim, $         
            isotropic=isotropic, $
            /fill
;            /cell_fill
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
            charsize=1.2, $       ;overall charsize (to set smaller title)
            xcharsize=0.0001, $     ;x-axis charsize (will be set below)
            ycharsize=0.0001, $     ;y-axis charsize (will be set below)
            title=titlestr, $
            xtitle=xtitlestr, $
            ytitle=ytitlestr, $         
            xrange=xlim, $
            yrange=ylim, $         
            isotropic=isotropic, $
            /fill
;            /cell_fill
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
            charsize=1.2, $       ;overall charsize (to set smaller title)
            xcharsize=0.0001, $     ;x-axis charsize (will be set below)
            ycharsize=0.0001, $     ;y-axis charsize (will be set below)
            title=titlestr, $
            xtitle=xtitlestr, $
            ytitle=ytitlestr, $         
            xrange=xlim, $
            yrange=ylim, $         
            isotropic=isotropic, $
            /fill
;            /cell_fill
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
   contour, arr2d, arrx, arry, $
            levels=levels_contour, $
            c_color=c_colors, $
            /fill, $
;            /cell_fill, $
            isotropic=isotropic, $
            irregular=irregular, $
            /overplot
endif else begin
   contour, arr2d, arrx, arry, $
            levels=levels_contour, $
            c_color=c_colors, $
            /fill, $
;            /cell_fill, $
            xstyle=1, ystyle=1 , $
            isotropic=isotropic, $
            irregular=irregular, $
            charsize=1.2, $       ;overall charsize (to set smaller title)
            xcharsize=1.2, $     ;x-axis charsize
            ycharsize=1.2, $     ;y-axis charsize
            title=titlestr, $
            xtitle=xtitlestr, $
            ytitle=ytitlestr, $         
            xrange=xlim, $
            yrange=ylim
endelse

;
;overplot the isocontours, if specified
if (keyword_set(isoc)) then begin
   contour, arr2d, arrx, arry, $
            levels=levels_iso, $
            color=0, $
            charsize=2., $
            c_labels=0, $
            /overplot
endif
;
;overplot the axes
xaxis=[min(arrx), max(arrx)]
yaxis=[min(arry), max(arry)]
origin=[0., 0.]
;oplot, xaxis, origin, line=0, color=0
;oplot, origin, yaxis, line=0, color=0
;
;-------------------------re-plot the axes------------------------------
;
if(keyword_set(background1) or keyword_set(background2) or keyword_set(background3)) then begin
;set zero color-table
   loadct, 0
;re-plot the x,y-axes, if background-keyword is set (contour with overplot
;will overplot the tick-marks)
;plot lower x-axis
   axis, xlim(0), ylim(0), xaxis=0, xrange=[xlim(0),xlim(1)], /xs, charsize=1.2, xtitle=xtitlestr
;plot upper x-axis (without labels)
   axis, xlim(0), ylim(1), xaxis=1, xrange=[xlim(0),xlim(1)], /xs, xtickformat='(a1)'
;plot left y-axis
   axis, xlim(0), ylim(0), yaxis=0, yrange=[ylim(0),ylim(1)], /ys, charsize=1.2, ytitle=ytitlestr
;plot right y-axis (without labels)
   axis, xlim(1), ylim(0), yaxis=1, yrange=[ylim(0),ylim(1)], /ys, ytickformat='(a1)'

   
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
;!p.region=[0.8, 0.0, 1., 1.]
!p.region=[x0c,y0c,x1c,y1c]
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
if(keyword_set(ctitlestr)) then begin
   colorbar, ncolors=ncolors, /vertical, $
             divisions=nlevels_iso-1, $
             tickvalues=cb_indx, $
             tickname=cb_tickmark_name, $
             position=pos, $
             title=ctitlestr, $
             /right, $
             charsize=1.7
endif else begin
   colorbar, ncolors=ncolors, /vertical, $
             divisions=nlevels_iso-1, $
             tickvalues=cb_indx, $
             tickname=cb_tickmark_name, $
             position=pos, $
             /right, $
             charsize=1.7
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
