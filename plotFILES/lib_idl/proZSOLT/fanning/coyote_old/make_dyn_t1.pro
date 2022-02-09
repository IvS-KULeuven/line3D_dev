pro md1,bet=bet,inc=inc,col=col,ps=ps,con=con,phimi=phimi,phima=phima,ran=ran,rt=rt,file_sav=file_sav

;Make a dynamic Ha spetra to the paper, 
;HERE: NO MAX/MIN FLUX PROFILES..

;INPUT PARAMETERS 

;inclination: inc 
;obliquity:   bet  
 
!p.multi=[0,1,3]

;OBSERVATIONS FROM HOWARTH - as idl svae set
file = '/data-sims/Jon/formal_mhd/observations/HD191612/JonGrey/fdyn_obs_interp_phase.sav'
;restore,'../../observations/HD191612/JonGrey/fdyn_obs_interp.sav'
restore,file
;convert to velocities
zzo = fdyn_obs
xxo = xobs
l0 = 6562.84
cc = 2.998e5
xxo = (xxo-l0)*cc/l0
;---------------------
yyo = phase

;NOW SIMULATIONS 
;First interpolate background grid
   if not keyword_set(file_sav) then file_sav='flux_grid_paper_new.sav'
   dir = '/data-sims/Jon/formal_mhd/formal_3d_ha/grid_HD191612_2.5d_phot3/'
   file_sav = dir+file_sav
   restore,file_sav
   xobs = w
   alpha = alpha*!pi
   bet = bet*(!pi/180.)
   inc = inc*(!pi/180.)
   nel = 200
   ;MAX AND MIN PHASES 
   if not keyword_set(phima) then phima = 1.5
   phima = 2*!pi*phima
   if not keyword_set(phimi) then phimi = -0.5
   phimi = phimi*2.*!pi
   ;-----------------------------------------
   phi = findgen(nel)/(float(nel)-1.)*(phima-phimi)+phimi
   ;Linear phi-grid 
   cosa = cos(bet)*cos(inc)+cos(phi)*sin(bet)*sin(inc)
   alp = acos(cosa)
   alp_sym = !pi-alp 
   ;USE ASSUMED NORTH/SOUTH POLE SYMMETRY, see Sundqvist et al. 2011
   print,'Max/Min Angles/pi for given i,beta:'
   print,max(alp)/!pi,min(alp)/!pi
   print,'Min/Max SYMMETRY Angles/pi for given i,beta:'
   print,min(alp_sym)/!pi,max(alp_sym)/!pi
   print,'-------------------------------------------'
   print,'Now interpolate background grid'
   ;-----------------------------------------
   nint = n_elements(w)
   nf = n_elements(fdyn(0,*))
   flux_mat  = fltarr(nint,nel)
   fluxa_mat = fltarr(nint,nel)
   fluxe_mat = fltarr(nint,nel)
   print,'help: fdyn,flux_mat:'
   help,fdyn,flux_mat
   ;Now symmetrize north/south pole directions 
   for i=0,nint-1 do begin 
      fd1 = interpol(fdyn(i,*),alpha,alp)
      fd2 = interpol(fdyn(i,*),alpha,alp_sym)
      flux_mat(i,*)  = (fd1+fd2)*0.5

      fd1 = interpol(fadyn(i,*),alpha,alp)
      fd2 = interpol(fadyn(i,*),alpha,alp_sym)
      fluxa_mat(i,*)  = (fd1+fd2)*0.5

      fd1 = interpol(fedyn(i,*),alpha,alp)
      fd2 = interpol(fedyn(i,*),alpha,alp_sym)
      fluxe_mat(i,*)  = (fd1+fd2)*0.5

   endfor
   beta = phi/(2.*!pi)
   phase = beta
;DONE, NOW PLOT --------------------------

flux = flux_mat

zz = flux
xx = xobs
yy = phase
print,'help,xx,y,,zz:'
help,xx,yy,zz

if not keyword_set(ran) then ran = 600.
xr = [-ran,ran]
;Plot ranges! 

;convolve ? 
if keyword_set(con) then begin 
   np = n_elements(yy)
   xx = l0*(1.+xx/cc)
   ;wavelength required in convol..
   for i=0,np-1 do begin 
      zzd = zz(*,i)
      if not keyword_set(rt) then begin 
         convol,xx,zzd,xxc,zzc,vdop=con
      endif else begin 
         convol,xx,zzd,xxc,zzc,vmacro=con
      endelse
      ;RAD-TAN (a la Gray) or isotropic macro-turbulence ? 
      zz(*,i) = interpol(real_part(zzc),real_part(xxc),xx)
   endfor
   ;transform back
   xx = (xx-l0)*cc/l0
endif
;------------------------------------------


;DATA PREPARED, 
;NOW PLOTTING SPECIFCS 
;-----------------------------------------

if not keyword_set(col) then col=0
if (col eq 100) then cgloadct,0 else cgloadct,col
;colors 

;if keyword_set(ps) then begin 
;   ;Use an ansatz col=100 !
;   psfullps,col=col
   !p.font=3
   !p.charsize=1.3
   !x.charsize=1.3
   !y.charsize=1.3
   !p.thick=5.0
   !x.thick=3.0
   !y.thick=3.0
;endif
;print to ps-file? 

navailable=!d.table_size  
ncolors=navailable
;ASSUME that number of available colors is maximum of present color-table
colorstep = navailable / ncolors  
c_colors  = indgen(ncolors)*colorstep 
;Make manual stepsizes between levels 
ii = where(zzo eq max(zzo),count)
if (count eq 0) then message,'no max value ??'
z1 = zzo(ii(0))
ii = where(zz eq max(zz),count)
if (count eq 0) then message,'no max value ??'
zma = [z1,zz(ii(0))]
;MAX VALUES CHECKED FOR BOTH 
ii = where(zzo eq min(zzo),count)
if (count eq 0) then message,'no min value ??'
z1 = zzo(ii(0))
ii = where(zz eq min(zz),count)
if (count eq 0) then message,'no min value ??'
zmi = [z1,zz(ii(0))]
mini = min(zmi) 
maxi = max(zma)  

levels=float(ncolors) 
userlevels = findgen(levels) * (maxi - mini)/levels + mini
print,'minimum/maximum values/nb. levels:',mini,maxi,n_elements(userlevels)
print,'Min/Max for [obs,sim]:'
print,zmi
print,zma

;observations
pos = [0.15,0.2,0.45,0.95]
cgcontour,zzo,xxo,yyo,c_linestyle=0,/cell_fill,/follow,levels=userlevels,$
  /xstyle,/ystyle,c_colors=c_colors,xticks=3,$
  xcharsize=3.3,ycharsize=3.3,position=pos,$
  ytitle='Phase',xr=xr;, /traditional
;simulations
pos = [0.47,0.2,0.77,0.95]
cgcontour,zz,xx,yy,c_linestyle=0,/cell_fill,/follow,levels=userlevels,$
  /xstyle,/ystyle,c_colors=c_colors,xticks=3,$
  xcharsize=3.3,ycharsize=3.3,position=pos,$
  ytickname=replicate(' ',10),xtickname=[' ','-200','200',' '],$
  xr=xr;, /traditional
;Cell_fill or fill ?? 
cgtext,0.35,0.11,textoidl('Velocity [km/s]'),/normal,charsize=2.0
cgtext,0.2,0.05,'Observations',/normal,charsize=2.3
cgtext,0.53,0.05,'Simulations',/normal,charsize=2.3

;courtesy of D. Fanning, www.dfanning.com
pos = [0.85, 0.15, 0.9, 0.95] 
format='(F0.3)'
title ='Flux'
;if (col eq 100) then loadct,0 else loadct,col
cgcolorbar,range=[mini,maxi],ncolors=ncolors*colorstep,$
  format=format,charsize=2.0,position=pos,/vertical,title=title


;if keyword_set(ps) then begin 
;   psoff
   !p.font=3
   !p.charsize=1.0
   !x.charsize=1.0
   !y.charsize=1.0
;;   spawn,'ps2pdf idl.eps'
;endif

!p.multi=0
   

end

;@colorbar.pro
;@psfullps.pro
;@psoff.pro
;@integ.pro
;@legend.pro
