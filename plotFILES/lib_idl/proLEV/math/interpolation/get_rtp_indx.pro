pro get_rtp_indx, rin, thetain, phiin, rarr, tharr, phiarr, nr, ntheta, nphi, indx_r1, indx_r2, indx_th1, indx_th2, indx_phi1, indx_phi2, expol
;
;-----------------------------------------------------------------------
;
;  finds indices of grid for a cube surrounding the point (rin,thetain,phiin)
;
;  input:  coordinates of point:             xin, yin, zin
;          dimension of grid:                ndxmax, ndymax, ndzmax
;          x,y,z-grid:                       xarr, yarr, zarr
;  output: indices of x-grid:                indx_x1, indx_x2
;          indices of y-grid:                indx_y1, indx_y2
;          indices of z-grid:                indx_z1, indx_z2
;          flag if extrapolation is needed:  expol
;
;-----------------------------------------------------------------------
;
expol=0
;
rmin=min(rarr)
rmax=max(rarr)
;
indx_r1=0
indx_r2=1
for i=1, nr-1 do begin
   if(rarr(i) ge rin) then begin
      indx_r1=i-1
      indx_r2=i
      break
   endif
endfor
;
if(rin gt rmax) then begin
   indx_r1=nr-2
   indx_r2=nr-1
   expol=1
endif
;
if(rin lt rmin) then begin
   indx_r1=0
   indx_r2=1
   expol=1
endif
;
;
;
thmin=min(tharr)
thmax=max(tharr)
indx_th1=0
indx_th2=1
for i=1, ntheta-1 do begin
   if(tharr(i) ge thetain) then begin
      indx_th1=i-1
      indx_th2=i
      break
   endif
endfor
;
if(thetain gt thmax) then begin
   indx_th1=ntheta-2
   indx_th2=ntheta-1
   expol=1
endif
;
if(thetain lt thmin) then begin
   indx_th1=0
   indx_th2=1
   expol=1
endif
;

;
phimin=min(phiarr)
phimax=max(phiarr)
indx_phi1=0
indx_phi2=1
for i=1, nphi-1 do begin
   if(phiarr(i) ge phiin) then begin
      indx_phi1=i-1
      indx_phi2=i
      break
   endif
endfor
;
if(phiin gt phimax) then begin
   indx_phi1=nphi-2
   indx_phi2=nphi-1
   expol=1
endif
;
if(phiin lt phimin) then begin
   indx_phi1=0
   indx_phi2=1
   expol=1
endif
   
;
;
;
end
