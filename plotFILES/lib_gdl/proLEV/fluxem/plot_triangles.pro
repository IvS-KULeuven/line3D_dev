pro plot_triangles, dir, fname=fname, obs_indx=obs_indx, oname=oname, ylim=ylim, xlim=xlim, oplot=oplot, cindx=cindx, help=print_help
;
;+
; NAME:
;	plot_triangles
;
; PURPOSE:
;	This procedure plots the tirangulation for the binary spec program
;
; CALLING SEQUENCE:
;
;	plot_triangles, dir
;
; INPUTS:
;	dir:	Directory, where appropriate files are stored.
;	
; KEYWORD PARAMETERS:
;       fname:  Set this keyword to the name of a file inside directory dir, 
;               from which the triangles shall be read in
;
;       obs_indx: Set this keyword to an integer, to describe from which file
;                 inside directory dir, the triangles are read in:
;                 e.g. obs_indx=1 reads from 'dir/spec_triangles_00001.dat'
;
;       oname:  Set this keyword to the output-name of a ps-file, if output
;               shall be plotted to that file
;
;       windx:  Set this keyword to an integer defining the window, to which
;               output is plotted
;
;       ylim:   Set this keyword to a 2-d array, which sets the yrange
;
;       xlim:   Set this keyword to a 2-d array, which sets the xrange
;
;       oplot:  Set this keyword (flag), to overplot in the current window
;               This keyword is not allowed, when oname is set
;
;       cindx:  Set this keyword to an integer, describing the color in which
;               the emergent profile is plotted (affects only overplotted
;               profiles)
;
;       isotropic: Set this keyword to make isotropic plot
;
;       help:   Set this keyword (flag) to show the documentation of the
;               procedure
;
; EXAMPLE:
;
;-
;
;-----------------------IF HELP IS NEEDED-------------------------------
;
IF(KEYWORD_SET(print_help)) THEN BEGIN
   doc_library, 'plot_triangles'
   return
ENDIF
;
IF (N_ELEMENTS(dir) EQ 0) THEN BEGIN
   PRINT, 'SYNTAX: '
   PRINT, 'plot_triangles, dir'
   PRINT, 'with dir=string , e.g. plot_triangles, "." '
   RETURN
ENDIF
;
IF (NOT KEYWORD_SET(obs_indx)) THEN obs_indx=0
;
;
;-----------------------------------------------------------------------
;
getTRIANGLES, dir, obs_indx, fexist, x1, y1, x2, y2, x3, y3, file_name=fname
;
;-----------------------------------------------------------------------
;
IF(fexist EQ 0) THEN BEGIN
   PRINT, 'ERROR IN OPENING FILES: RETURNING'
   RETURN
ENDIF
;
;---------------------CALCULATE PLOT RANGES-----------------------------
;
IF(NOT KEYWORD_SET(ylim)) THEN ylim=[min([y1,y2,y3]),max([y1,y2,y3])]     
IF(NOT KEYWORD_SET(xlim)) THEN xlim=[min([x1,x2,x3]),max([x1,x2,x3])]     
;
;------------------------TITLES-----------------------------------------
;
;on default
xtitleStr=textoidl('x_p')
xtitleStr=textoidl('y_p')
;
;----------------------PLOT---------------------------------------------
;
if(keyword_set(oname)) then begin
;set output file
   PRINT, "WRITING OUTPUT TO: ", ONAME
   set_plot,'ps'
   device,file=ONAME
   loadct, 0
endif

loadct, 0
if(keyword_set(cindx)) then loadct, 12
;
;-----------------------------------------------------------------------
;
plot, [xlim(0),xlim(0)],[ylim(0),ylim(0)], $
      xrange=xlim, $
      yrange=ylim, $
      xtitle=xtitlestr, $                   
      ytitle=ytitlestr

for i=0, n_elements(x1)-1 do begin
   oplot, [x1(i),x2(i)], [y1(i),y2(i)], color=cindx
   oplot, [x2(i),x3(i)], [y2(i),y3(i)], color=cindx
   oplot, [x3(i),x1(i)], [y3(i),y1(i)], color=cindx
endfor
;

IF KEYWORD_SET(ONAME) THEN BEGIN
   device, /close
   set_plot,'x'
ENDIF

;
;

END

