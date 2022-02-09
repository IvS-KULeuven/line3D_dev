pro benchmark04, dir=DIR, windx=WINDX, oname=ONAME, xlim=XLIM, ylim=YLIM
;
;plots contours of searchlight along a direction
;
;------------------read all information from hdf5-file------------------
;
if(not keyword_set(dir)) then dir='.'
fname=dir+'/benchmark04.h5'
;
file_id = h5f_open(fname)
;
   group_id = h5g_open(file_id, 'dimensions')
      att_id=h5a_open_name(group_id, 'nr')
         nr=h5a_read(att_id)
         nr=nr(0)
      h5a_close, att_id
      att_id=h5a_open_name(group_id, 'dim_mu')
         dim_mu=h5a_read(att_id)
         dim_mu=dim_mu(0)
      h5a_close, att_id
      att_id=h5a_open_name(group_id, 'n1d_angdep')
         nr_angdep=h5a_read(att_id)
         nr_angdep=nr_angdep(0)
      h5a_close, att_id
   h5g_close, group_id
;
   r=fltarr(nr)
   r_angdep=fltarr(nr_angdep)
   mint_cr=fltarr(nr)
   mint_theo=fltarr(nr)
   group_id = h5g_open(file_id, 'coordinates')
      dset_id=h5d_open(group_id, 'r')
         r=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'r1d_angdep')
         r_angdep=h5d_read(dset_id)
      h5d_close, dset_id
   h5g_close, group_id
;
   nodes_mu = fltarr(dim_mu)
   group_id = h5g_open(file_id, 'angles')
      dset_id=h5d_open(group_id, 'nodes_mu')
         nodes_mu=h5d_read(dset_id)
      h5d_close, dset_id
   h5g_close, group_id
;
   group_id = h5g_open(file_id, 'solution_cr')
      dset_id=h5d_open(group_id, 'mint_sc')
         mint_sc=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'mint_fvm')
         mint_fvm=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'mint_theo')
         mint_theo=h5d_read(dset_id)
      h5d_close, dset_id
   h5g_close, group_id
;
   intsc_angdep=fltarr(nr_angdep,dim_mu)
   intfvm_angdep=fltarr(nr_angdep,dim_mu)
   group_id = h5g_open(file_id, 'solution_angdep')
      dset_id=h5d_open(group_id, 'intsc_angdep')
         intsc_angdep=h5d_read(dset_id)
      h5d_close, dset_id
      dset_id=h5d_open(group_id, 'intfvm_angdep')
         intfvm_angdep=h5d_read(dset_id)
      h5d_close, dset_id
   h5g_close, group_id
;
h5f_close, file_id
;
;-----------------------------------------------------------------------
;
;********************plot mean intensities******************************
;
;-----------------------define range------------------------------------
;
;
if(not keyword_set(xlim)) then begin
   xmin=min(r)
   xmax=max(r)
   xlim=[xmin, xmax]
endif
;
if(not keyword_set(ylim)) then begin
   ymin=min([min(mint_sc*r^2), min(mint_fvm*r^2), min(mint_theo*r^2)])
   ymax=max([max(mint_sc*r^2), max(mint_fvm*r^2), max(mint_theo*r^2)])
   ylim=[ymin, ymax]
endif
;
;----------------------------title-strings------------------------------
;
;titlestr=textoidl('i (\theta=') + string(acos(mu), format='(f9.5)') + $
;         textoidl(', \phi=') + string(phi, format='(f9.5)') + ')'
titlestr=''
xtitlestr=textoidl('r')
ytitlestr=textoidl('J*r^2')
;
;-----------------------------------------------------------------------
;
if(not keyword_set(windx)) then windx=0
;
if(keyword_set(oname)) then begin
   oname='ps_files/benchmark04_a.ps'
   print, "writing output to: ", oname
   set_plot,'ps'
   device,file=oname, $ ;xsize=19., ysize=26.7, xoffset=1., yoffset=1. , $
    decomposed=0, color=1, bits_per_pixel=8
endif else begin
   window, windx
   device, decomposed=0
   windx=windx+1
endelse
;
loadct_lev, ci_black, ci_blue, ci_red, ci_cyan, ci_green, ci_magenta, ci_yellow, ci_white
;
!p.multi=[2,1,2]
;
plot, r, mint_theo*r^2, $
      xrange=xlim, $
      yrange=ylim, $
      title=titleStr, $
      xtitle=xtitleStr, $
      ytitle=ytitleStr

oplot, r, mint_sc*r^2, $
      color=ci_blue, $
      line=2
oplot, r, mint_sc*r^2, color=ci_blue, psym=1

oplot, r, mint_fvm*r^2, $
      color=ci_red, $
      line=2
oplot, r, mint_fvm*r^2, color=ci_red, psym=1

for i=0, nr_angdep-1 do begin
   oplot, [r_angdep(i),r_angdep(i)],[0.,1.], line=1
endfor
;
lstr1='theoretical (from dilution)'
lstr2='2d short characteristics'
lstr3='2d finite volume method'
if(keyword_set(oname)) then begin
   legend, [lstr1, lstr2, lstr3], $
         psym=[0,0,0], $
         linestyle=[0,2,2], $
         color=[ci_black,ci_blue,ci_red], $
         textcolor=ci_black, $
         /right_legend
endif else begin
   legend, [lstr1, lstr2, lstr3], $
         psym=[0,0,0], $
         linestyle=[0,2,2], $
         color=[ci_white,ci_blue,ci_red], $
         textcolor=ci_white, $
         /right_legend
endelse
;
!p.multi=[1,1,2]
plot, [0.,max(r)],[1.,1.], yrange=[0.,2.],xrange=[0.,max(r)], $
      line=1, $
      xtitle='r', $
      ytitle='J(num)/J(theo)'
oplot, r, mint_sc/mint_theo, $
       color=ci_blue
oplot, r, mint_fvm/mint_theo, $
       color=ci_red
lstr1='2d short characteristics'
lstr2='2d finite volume method'
if(keyword_set(oname)) then begin
   legend, [lstr1, lstr2], $
         psym=[0,0], $
         linestyle=[0,0], $
         color=[ci_blue,ci_red], $
         textcolor=ci_black, $
         /right_legend
endif else begin
   legend, [lstr1, lstr2], $
         psym=[0,0], $
         linestyle=[0,0], $
         color=[ci_blue,ci_red], $
         textcolor=ci_white, $
         /right_legend
endelse
;
if keyword_set(oname) then begin
   device, /close
   set_plot,'x'
endif
;
!p.multi=0
;
;********************plot intensities as function of angle**************
;
;-----------------------define range------------------------------------
;
if(keyword_set(oname)) then begin
   oname='ps_files/benchmark04_b.ps'
   print, "writing output to: ", oname
   set_plot,'ps'
   device,file=oname, $ ;xsize=19., ysize=26.7, xoffset=1., yoffset=1. , $
    decomposed=0, color=1, bits_per_pixel=8
endif else begin
   window, windx, xsize=800, ysize=1000
   device, decomposed=0
   windx=windx+1
endelse
;
nline=ceil(sqrt(nr_angdep))
ncols=floor(sqrt(nr_angdep))
;
print, 'r   ', 'J_sc *r^2   ', 'J_theo * r^2   ', 'J_sc/J_theo   '
for i=0, nr_angdep-1 do begin
   !p.multi=[nline*ncols-i, ncols, nline]
   plot_intangdep, r_angdep(i), nodes_mu, intsc_angdep(i,*), intfvm_angdep(i,*), oname=oname
;
;perform integration
   jnum = 0.5d0*integ_trapez(nodes_mu, intsc_angdep(i,*))
   jtheo = 0.5d0*(1.d0-sqrt(1.d0-1.d0/r_angdep(i)^2))
   print, r_angdep(i), jnum*r_angdep(i)^2, jtheo*r_angdep(i)^2, jnum/jtheo
endfor
;
if keyword_set(oname) then begin
   device, /close
   set_plot,'x'
endif
;
!p.multi=0


end
;
;-----------------------------------------------------------------------
;-----------------------------------------------------------------------
;-----------------------------------------------------------------------
;
pro plot_intangdep, rp, nodes_mu, intsc_angdep, intfvm_angdep, oname=oname
;
;-----------------------define range------------------------------------
;
xlim=[min(nodes_mu), max(nodes_mu)]
;
ymin=min([min(intsc_angdep),min(intfvm_angdep)])
ymax=max([max(intsc_angdep),max(intfvm_angdep)])

dy=ymax-ymin
if(dy eq 0.) then begin
   ymin=ymin-0.1d0
   ymax=ymax+0.1d0
endif else begin
   ymin=ymin-0.1d0*dy
   ymax=ymax+0.1d0*dy
endelse
ylim=[ymin,ymax]
;
;----------------------------title-strings------------------------------
;
;titlestr=textoidl('i (\theta=') + string(acos(mu), format='(f9.5)') + $
;         textoidl(', \phi=') + string(phi, format='(f9.5)') + ')'
titlestr='at point r=' + string(rp, format='(f6.3)')
xtitlestr=textoidl('\mu')
ytitlestr=textoidl('I(\mu)')
;
;-----------------------------------------------------------------------
;
loadct_lev, ci_black, ci_blue, ci_red, ci_cyan, ci_green, ci_magenta, ci_yellow, ci_white
;

mu_star=sqrt(1.d0-1.d0/rp^2)
;xlim=[mu_star-0.1,mu_star+0.05]
xlim=[mu_star-0.2,mu_star+0.1]
plot, [mu_star,mu_star], [ylim(0),ylim(1)], $
      xrange=xlim, /xs, $
      yrange=ylim, /ys, $
      title=titleStr, $
      xtitle=xtitleStr, $
      ytitle=ytitleStr, $
      charsize=2., $
      line=1

oplot, nodes_mu, intsc_angdep, $
      color=ci_blue
oplot, nodes_mu, intsc_angdep, psym=1;, color=ci_blue

oplot, nodes_mu, intfvm_angdep, $
      color=ci_red
oplot, nodes_mu, intfvm_angdep, psym=1, color=ci_red
;
;
lstr1='2d short characteristics'
lstr2='2d finite volume method'
lstr3='mu_star'
;
if(keyword_set(oname)) then begin
   legend, [lstr1, lstr2, lstr3], $
         psym=[0,0,0], $
         linestyle=[1,0,0], $
         color=[ci_blue,ci_red,ci_black], $
         textcolor=ci_black
endif else begin
   legend, [lstr1, lstr2, lstr3], $
         psym=[0,0,0], $
         linestyle=[0,0,1], $
         color=[ci_blue,ci_red,ci_white], $
         textcolor=ci_white
endelse
;
end
