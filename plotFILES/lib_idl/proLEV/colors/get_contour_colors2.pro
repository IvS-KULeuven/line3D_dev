pro get_contour_colors2, clim, ncolors, nlevels_iso, bottom, $
                      levels_final, levels_iso, c_colors_final,  $
                      cb_indx, cb_tickmark_name

if(n_elements(clim) ne 2) then begin
   print, 'error in get_contour_colors2: n_elements(clim) ne 2'
endif
;
;define total number of colors
ncolors=!D.TABLE_SIZE
bottom=0
;
;define the levels from plot-limits
levels_final=clim(0)+indgen(ncolors)*(clim(1)-clim(0))/float(ncolors-1)
;
;define corresponding color-indices
c_colors_final=bottom+indgen(ncolors)*(ncolors-bottom)/ncolors
;
;define levels of iso-contours and colorbar-indices
nlevels_iso=6
cb_indx=bottom+indgen(nlevels_iso)*floor((ncolors-bottom)/float(nlevels_iso-1))
;
levels_iso=levels_final(cb_indx)
cb_tickmark=levels_final(cb_indx)

;cb_tickmark_name = STRING(cb_tickmark, FORMAT='(f6.2)')
tickmax=max(abs(cb_tickmark))
if(tickmax lt 10.) then begin
   cb_tickmark_name = STRING(cb_tickmark, FORMAT='(f5.2)')
endif else begin
   if(tickmax lt 100.) then begin
      cb_tickmark_name = STRING(cb_tickmark, FORMAT='(f6.2)')
   endif else begin
      if(tickmax lt 1000.) then begin
         cb_tickmark_name = STRING(cb_tickmark, FORMAT='(f7.2)')
      endif else begin
         if(tickmax lt 10000.) then begin
            cb_tickmark_name = STRING(cb_tickmark, FORMAT='(f8.2)')
         endif else begin
            cb_tickmark_name = STRING(cb_tickmark, FORMAT='(e9.2)')
         endelse
      endelse
   endelse
endelse
;
;
END
