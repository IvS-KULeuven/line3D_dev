pro get_linecenter, xobs, fnorm, xcenter, flim=flim, xlim=xlim, xobsl=xobsl, fnorml=fnorml, xobsr=xobsr, fnormr=fnormr
;
;get the center of an emission line with the barycenter method
;input: xobs in vth_fiducial
;       fnorm: normalized flux
;
;flim: only consider points left and right of a certain height
;      specified by flim (in units of max(fnorm)-1.)
;xlim: only consider points with xobs in [xobs(flim),xmax]
;
if(not keyword_set(flim)) then flim=1.d0
;
dfmax = max(fnorm)-1.d0
dfmin = 1.d0-min(fnorm)
;
if(dfmax ge dfmin) then begin
   labs=0 
   lemi=1  ;emission profile
endif else begin
   labs=1  ;absorption profile       
   lemi=0
endelse

if(lemi eq 1) then begin
   fp = 1.d0+(max(fnorm)-1.d0)*flim
endif else begin
   fp = 1.d0-(1.d0-min(fnorm))*flim
endelse
;
nd = n_elements(xobs)
;
;--------------find indices for left part of the profile----------------
;
;starting index
indxl2=0
for i=0, nd-1 do begin
;   print, i, fnorm(i), fmax, fmin
   if(fnorm(i) ge fp and lemi eq 1) then begin
      indxl2=i
      break
   endif
   if(fnorm(i) le fp and labs eq 1) then begin
      indxl2=i
      break
   endif   
endfor
;
if(indxl2 eq 0) then begin
   print, 'error in get_linecenter: indxl2 = 0'
   stop
endif
;
;interpolate xobs to find limit
xobsl2 = interpol_ypfct(fnorm(indxl2), fnorm(indxl2-1), xobs(indxl2), xobs(indxl2-1), fp)
;
if(not keyword_set(xlim)) then begin
   xobsl1 = min(xobs)
endif else begin
   xobsl1 = xobsl2 - xlim
endelse

;window, 0
;plot, xobs, fnorm

;print, indxl2, indxl2-1
;print, fnorm(indxl2), fp, fnorm(indxl2-1), fmin
;print, xobs(indxl2), xobs(indxl2-1)
;print, xobsl1, xobsl2, min(xobs)
;print, ''
;
if(xobsl1 ge xobsl2) then begin
   print, 'error in get linecenter: xobsl1 > xobsl2'
   stop
endif
;
if(xobsl1 lt min(xobs)) then begin
   print, 'error in get linecenter: xobsl1 < min(xobs)'
   stop
endif
;
;define flux profile with doubled resolution
indx = where(xobs ge xobsl1 and xobs le xobsl2)
ndl = n_elements(indx)*2
xobsl = xobsl1 + findgen(ndl)*(xobsl2-xobsl1)/(ndl-1)
fnorml = fltarr(ndl)*0.d0
for i=0, ndl-1 do begin
   find_indx, xobsl(i), xobs, nd, iim1, ii, help=print_help
   fnorml(i) = interpol_ypfct(xobs(iim1),xobs(ii),fnorm(iim1),fnorm(ii),xobsl(i))
endfor
;
;--------------find indices for right part of the profile---------------
;
if(lemi eq 1) then begin
   fp = 1.d0+(max(fnorm)-1.d0)*flim
endif else begin
   fp = 1.d0-(1.d0-min(fnorm))*flim
endelse
;
;starting index
indxr1=nd-1
for i=nd-1, 0, -1 do begin
   if(fnorm(i) ge fp and lemi eq 1) then begin
      indxr1=i
      break
   endif
   if(fnorm(i) le fp and labs eq 1) then begin
      indxr1=i
      break
   endif
endfor
;
if(indxr1 eq nd-1) then begin
   print, 'error in get_linecenter: indxr1 = nd-1'
   stop
endif
;
;interpolate xobs to find limit
xobsr1 = interpol_ypfct(fnorm(indxr1), fnorm(indxr1+1), xobs(indxr1), xobs(indxr1+1), fp)
;
if(not keyword_set(xlim)) then begin
   xobsr2 = xobsr1 + xobsl2-xobsl1
endif else begin
   xobsr2 = xobsr1 + xlim
endelse
;
if(xobsr2 le xobsr1) then begin
   print, 'error in get linecenter: xobsr2 < xobsr1'
   stop
endif
;
;print, max(xobs)
;print, xobsr2
;print, ''
if(xobsr2 gt max(xobs)) then begin
   print, 'error in get linecenter: xobsr2 > max(xobs)'
   stop
endif
;
;define flux profile with doubled resolution
indx = where(xobs ge xobsr1 and xobs le xobsr2)
ndr = n_elements(indx)*2
xobsr = xobsr1 + findgen(ndr)*(xobsr2-xobsr1)/(ndr-1)
fnormr = fltarr(ndr)*0.d0
for i=0, ndr-1 do begin
   find_indx, xobsr(i), xobs, nd, iim1, ii, help=print_help
   fnormr(i) = interpol_ypfct(xobs(iim1),xobs(ii),fnorm(iim1),fnorm(ii),xobsr(i))
endfor
;
;-------------------------get the barycenter----------------------------
;
;left part
suml1=0.d0
suml2=0.d0
for i=1, ndl-1 do begin
   suml1 = suml1 + 0.5d0*(xobsl(i-1)*fnorml(i-1)+xobsl(i)*fnorml(i))*(xobsl(i)-xobsl(i-1))
   suml2 = suml2 + 0.5d0*(fnorml(i-1)+fnorml(i))*(xobsl(i)-xobsl(i-1))
endfor
;
;right part
sumr1=0.d0
sumr2=0.d0
for i=1, ndr-1 do begin
   sumr1 = sumr1 + 0.5d0*(xobsr(i-1)*fnormr(i-1)+xobsr(i)*fnormr(i))*(xobsr(i)-xobsr(i-1))
   sumr2 = sumr2 + 0.5d0*(fnormr(i-1)+fnormr(i))*(xobsr(i)-xobsr(i-1))
endfor
;
;central part
sumc1=0.5d0*(fnorml(ndl-1)*xobsl(ndl-1)+fnormr(0)*xobsr(0))*(xobsr(0)-xobsl(ndl-1))
sumc2=0.5d0*(fnorml(ndl-1)+fnormr(0))*(xobsr(0)-xobsl(ndl-1))
;sumc1=0.d0
;sumc2=0.d0
;
;
;print, suml1, sumr1, sumc1
;print, suml2, sumr2, sumc2
xcenter=(suml1+sumr1+sumc1)/(suml2+sumr2+sumc2)




end
