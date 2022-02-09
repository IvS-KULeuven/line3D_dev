pro get_contour_colors, arr2d, arrx, arry, ncolors, nlevels_iso, bottom, $
                      levels_final, levels_iso, c_colors_final,  $
                      cb_indx, cb_tickmark_name
;
;---------PREPARING INPUT-DATA AND STYLE-DATA FOR CONTOUR-PLOT----------
;
;ROUND ARR2D TO NEAREST INTEGER NUMBER (WILL GIVE ONLY ZEROS FOR ARR2D<1)
round_dez = 1.d0
arr2d_color = arr2d
arr2d_round = ROUND(arr2d_color * round_dez, /L64)/round_dez
;get the contour levels for arr2d_round:
;find unique values and sort them
contour_levels = arr2d_round[UNIQ(arr2d_round, sort(arr2d_round))]
;get number of levels
nlevels = N_ELEMENTS(contour_levels)
;
;ROUND ARR2D TO NEAREST INTEGER NUMBER, WHERE ROUNDING IS ONLY PERFORMED FOR
;ARR2D*I^10 UNTIL NLEVELS IS GE !D.TABLE_SIZE (MORE LEVELS THAN AVAILABLE
;COLORS ARE FOUND)
ncolors_max = !D.TABLE_SIZE
;
FOR I=1, 16 DO BEGIN
   IF(NLEVELS LT ncolors_max) THEN BEGIN
      PRINT, 'NLEVELS IS ZERO, NEW ROUNDING'
      round_dez = round_dez*10.d0
      arr2d_round = ROUND(arr2d_color * round_dez, /L64)/round_dez
      contour_levels = arr2d_round[UNIQ(arr2d_round, SORT(arr2d_round))]
      nlevels = N_ELEMENTS(contour_levels)
   ENDIF
ENDFOR
;
;----------------FOLLOWING ONLY IF nlevels EQ ncolors_max---------------

if(nlevels eq ncolors_max) then begin
;set ncolors to total table size
   ncolors=ncolors_max
   nlevels_final = nlevels
   levels_final=contour_levels
   c_colors_final=indgen(ncolors) 
endif
;
;-----------------FOLLOWING ONLY IF nlevels GT ncolors_max--------------
;
if(nlevels gt ncolors_max) then begin
;set the number of colors to total table size
   ncolors=ncolors_max
   nlevels_final=ncolors
   levels_final=fltarr(ncolors)*0.d0
   c_colors_final=indgen(ncolors)
;
;find those levels, which shall get an own color
;first level
   nlevels_rest=nlevels
   ncolors_rest=ncolors
   indx=0
   levels_final(0) = contour_levels(indx)
   nlevels_rest = nlevels_rest-1
   ncolors_rest = ncolors_rest-1
;
;all other levels: only delk'th level is taken
;note: delk is variable in order that complete range 
;      up to contour_levels(nlevels) is included
   for i=1, ncolors-1 do begin
      delk = nlevels_rest/ncolors_rest
      indx = indx+delk
      nlevels_rest=nlevels_rest-delk
      ncolors_rest=ncolors_rest-1
      levels_final(i) = contour_levels(indx)
   endfor
endif
;
;----------------FOLLOWING ONLY IF nlevels LT ncolors_max---------------
;
if(nlevels lt 2) then begin
   print, 'error: nlevels lt 2, not two unique values in arr2d'
   stop
endif
;
;
;
if(nlevels lt ncolors_max) then begin
;set the number of colors to number of available levels
   ncolors=nlevels
;set the number of levels to number of existing levels
   nlevels_final=nlevels
   levels_final = contour_levels
;set final c_colors to an array of length nlevels
   c_colors=indgen(ncolors_max)
   c_colors_final=indgen(nlevels_final)*0
;
;find those color-indices which shall get an own
;data level
;first index
   nlevels_rest=nlevels
   ncolors_rest=ncolors_max
   indx=0
   c_colors_final(0) = c_colors(indx)
   nlevels_rest = nlevels_rest-1
   ncolors_rest = ncolors_rest-1
;
;all other color-indices: only delk'th color-index is taken
   for i=1, nlevels-1 do begin
      delk = ncolors_rest/nlevels_rest
      indx = indx+delk
      nlevels_rest = nlevels_rest-1
      ncolors_rest = ncolors_rest-delk
      c_colors_final(i) = c_colors(indx)
   endfor
endif
;
;-----------------MAKE ISOCONTOURS AND INDEX FOR COLORBAR---------------
;
;----------------AT FIRST: IF NLEVELS_ISO GE NLEVELS_FINAL--------------
;
;make (at maximum) 6 isocontours
nlevels_iso = 6 
;
if(nlevels_iso ge nlevels_final) then begin
;reset nlevels_iso
   nlevels_iso=nlevels_final
;set isocontours to all levels availabe
   levels_iso=levels_final
;set corresponding color-indices
   cb_indx = c_colors_final
;rewrite c_colors_final since loadct automatically adapts that for
;   ncolors gt nlevels
   c_colors_final=indgen(nlevels_iso)
;set corresponding tickmarks
   cb_tickmark = levels_final
endif
;
;---NOW: NLEVELS_ISO LT NLEVELS_FINAL AND NLEVELS_FINAL LT NCOLORS_MAX--
;
if(nlevels_iso lt nlevels_final and nlevels_final lt ncolors_max) then begin
   print, nlevels_iso, nlevels_final, ncolors, ncolors_max
;
   levels_iso=fltarr(nlevels_iso)*0.d0
   cb_indx=indgen(nlevels_iso)*0
   level_indx=indgen(nlevels_iso)*0
;
   nlevels_iso_rest = nlevels_iso
   nlevels_rest = nlevels_final
   indx=0
   level_indx(0) = indx
   levels_iso(0) = levels_final(indx)
   cb_indx(0) = c_colors_final(indx)
   nlevels_rest = nlevels_rest-1
   nlevels_iso_rest = nlevels_iso_rest-1
;
   for i=1, nlevels_iso-1 do begin
      delk = nlevels_rest/nlevels_iso_rest
      indx = indx+delk
      nlevels_rest = nlevels_rest - delk
      nlevels_iso_rest = nlevels_iso_rest - 1
      levels_iso(i) = levels_final(indx)
      cb_indx(i) = c_colors_final(indx)
      level_indx(i) = indx
   endfor

   c_colors_final=indgen(nlevels_final)
   cb_tickmark = levels_final(level_indx)
;
endif
;
;---NOW: NLEVELS_ISO LT NLEVELS_FINAL AND NLEVELS_FINAL EQ NCOLORS_MAX--
;
if(nlevels_iso lt nlevels_final and nlevels_final eq ncolors_max) then begin
;
;define levels of isocontours and colorbar-index
;   (index of colors used for contours)
   levels_iso=fltarr(nlevels_iso)*0.d0
   cb_indx=indgen(nlevels_iso)*0
;
   nlevels_iso_rest = nlevels_iso
   nlevels_rest = nlevels_final
   indx=0
   levels_iso(0) = levels_final(indx)
   cb_indx(0) = c_colors_final(indx)
   nlevels_rest = nlevels_rest-1
   nlevels_iso_rest = nlevels_iso_rest-1
;
   for i=1, nlevels_iso-1 do begin
      delk = nlevels_rest/nlevels_iso_rest
      indx = indx+delk
      nlevels_rest = nlevels_rest - delk
      nlevels_iso_rest = nlevels_iso_rest - 1
      levels_iso(i) = levels_final(indx)
      cb_indx(i) = c_colors_final(indx)
   endfor
;
   cb_tickmark = levels_final(cb_indx)
endif
;
;cb_tickmark_name = STRING(cb_tickmark, FORMAT='(E10.3)')
cb_tickmark_name = STRING(cb_tickmark, FORMAT='(f6.2)')

;
;-----------------------MAKE THE REST-----------------------------------
;
bottom=0
;
;
END
