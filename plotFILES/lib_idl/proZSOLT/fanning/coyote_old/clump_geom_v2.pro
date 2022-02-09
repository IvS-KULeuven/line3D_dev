pro cg,ps=ps 

;makes a little venetian blind graphics... 

if keyword_set(ps) then begin 
   psfullps
   !p.font=3 
   !p.charsize=1.5
   !x.charsize=1.5
   !y.charsize=1.5
   !p.thick=6.0
   !x.thick=5.0
   !y.thick=5.0
endif

p1 = [0.15,0.15,0.95,0.95]
nb = 1000
tsh=14.0
ths2=2.0
xr=[-1.0,8]
yr=[-0.5,6.0]

;SPHERICAL CLUMP 

Xo = 1.5
Yo = 3.0
Roff = [Xo,Yo]
;Takes out center of epherical clump 
r0 = 0.5 
th0 = findgen(nb)/(nb-1.)*2.*!pi
x0 = r0*cos(th0)+Xo
y0 = r0*sin(th0)+Yo
;Equation of clump -- just vector addition
cgplot,x0,y0,/isotropic,/nodata,xr=xr,yr=yr,ys=1,xs=1,$
     xticks=1,yticks=1,ytickname=replicate(' ',10),$
     xtickname=replicate(' ',10)
cgplot,x0,y0,thick=tsh, /overplot
;Plotting clump

;------------------------------------

;n-hat direction photon -- should also penetrate radial fragment!  
xmax = xr(1)-1.0
xmin = xr(0)+0.5
p0 = Yo
;AT CENTER OF SPHERE!
yp0 = fltarr(nb) + p0
xp0 = findgen(nb)/(nb-1.)*(xmax-xmin) + xmin
dx=0.05
cgarrow,xp0(nb-1),yp0(nb-1),xp0(nb-1)+dx,yp0(nb-1),/data,thick=3.0
cgplot,xp0,yp0,linestyle=1, /overplot
xshf=0.2
yshf=0.0
cgtext,max(xp0)+xshf,max(yp0)+yshf,textoidl('n_u'),charsize=2.0
nhat = [max(xp0),0.0]

;NOW RADIAL VECTOR FOR SPHERICAL CLUMP 

dydx = Yo/(Xo-r0)  
;Gradient since vector cuts [0,0]
;Use many points, make arrow simpler
X2 = Xo 
Y2 = X2*dydx
xp0 = findgen(nb)/(nb-1.)*X2
yp0 = fltarr(nb)
yp0 = interpol([0.0,Y2],[0.0,X2],xp0)
cgplot,xp0,yp0,linestyle=2, /overplot
ddx=0.05
cgarrow,xp0(nb-1),yp0(nb-1),xp0(nb-1)+ddx,yp0(nb-1)+dydx*ddx,/data,thick=3.0

xshf =0.15
yshf =0.2 
cgtext,X2+xshf,Y2+yshf,textoidl('r_u'),charsize=2.0
xshf =-1.2
yshf =0.7
Xsp = X2+xshf
Ysp = Y2+yshf
cgtext,Xsp,Ysp,'Spherical clumps',charsize=2.0


;Ok, NOW radial fragment a la TILTED rectangle!  

;FIRST RADIUS VECTOR 
xmin = 5.5
xmax = 5.7
xmid = (xmax+xmin)*0.5
;Arbitray 
dydx3 = Yo/xmin
X3 = xmax+1.0
Y3 = X3*dydx3
xp0 = findgen(nb)/(nb-1.)*X3
yp0 = fltarr(nb)
yp0 = interpol([0.0,Y3],[0.0,X3],xp0)
cgplot,xp0,yp0,linestyle=2, /overplot
ddx=0.05
cgarrow,xp0(nb-1),yp0(nb-1),xp0(nb-1)+ddx,yp0(nb-1)+dydx3*ddx,/data,thick=3.0
xshf=0.2
yshf=0.1
cgtext,max(xp0)+xshf,max(yp0)+yshf,textoidl('r_u'),charsize=2.0
;RADIUS VECTOR DONE 

;DOT PROD. 
rvec = [X3,Y3]
nhat = nhat
costh = transpose(nhat)#rvec/(norm(rvec)*norm(nhat))
th = acos(costh)
print,'angle:',th

ymin = Yo-r0
ymax = Yo+r0
ymid = Yo

;Now apply rotation matrix to radial fragment, I do it manually, 
;since so simple - NEED 4 transforms, since 4 end-points
rot = [[cos(th),-sin(th)],[sin(th),cos(th)]]
;NOTE -- NEED TO PARALLELL MOVE VECTORS HERE! 
xminh = xmin
xmaxh = xmax
yminh = ymin
ymaxh = ymax
;HOLD THESE 
xmin = xmin-xmid
xmax = xmax-xmid
ymin = ymin-ymid
ymax = ymax-ymid
;MOVED TO CENTER
transl = [xmid,ymid]
v11 = [xmin,ymin] 
v1 = transpose(rot)#v11 + transl
v22 = [xmax,ymin]
v2 = transpose(rot)#v22 + transl
v33 = [xmax,ymax]
v3 = transpose(rot)#v33 + transl
v44 = [xmin,ymax]
v4 = transpose(rot)#v44 + transl
;ROTATED - NAD MOVED BACK 
 
cgplot,[v1(0),v2(0),v3(0),v4(0),v1(0)],[v1(1),v2(1),v3(1),v4(1),v1(1)],thick=tsh, /overplot
;plot,[v11(0),v22(0),v33(0),v44(0),v11(0)],[v11(1),v22(1),v33(1),v44(1),v11(1)],thick=tsh

xshf = 4.0
cgtext,Xsp+xshf,Ysp,'Radial fragments',charsize=2.0
xshf=0.4
yshf=0.07
cgtext,xmid+xshf,Yo+yshf,textoidl('\theta'),charsize=1.4
;Origo
Xo=0.0
Yo=0.0
r0 = 0.1 
th0 = findgen(nb)/(nb-1.)*2.*!pi
x0 = r0*cos(th0)+Xo
y0 = r0*sin(th0)+Yo
cgplot,x0,y0, /overplot

if keyword_set(ps) then begin 
   psoff
   !p.font=3 
   !p.charsize=1.0
   !x.charsize=1.0
   !y.charsize=1.0
   !p.thick=1.0
   spawn,'ps2pdf idl.ps'
endif


end

