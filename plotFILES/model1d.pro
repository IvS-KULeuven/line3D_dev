pro model1d, dir=dir, windx=windx, oname=oname
;
;--------------------------read model-----------------------------------
;
if(not keyword_set(dir)) then dir='.'

fname=dir+'/model1d.h5'
;
read_model1d, fname, nr=nr1d, radius=r, rho1d=rho1d, velr1d=velr1d, t1d=t1d, vth1d=vth1d
;
   mdot1d = 4.d0*!pi*r^2*velr1d*rho1d*!yr/!msu
;
;-----------------------------------------------------------------------
;
if(not keyword_set(windx)) then windx=0
;
   window, windx, xsize=1000, ysize=1600
   device, decomposed=0
   windx=windx+1
;
   loadct_lev, ci_black, ci_blue, ci_red, ci_cyan, ci_green, ci_magenta, ci_yellow, ci_white
;
   !p.multi=[5,1,5]
   plot, r, rho1d, /ylog, $
         charsize=2.5, $
         xtitle=textoidl('r'), $
         ytitle=textoidl('\rho')

   !p.multi=[4,1,5]
   plot, r, velr1d, $
         charsize=2.5, $
         xtitle=textoidl('r'), $
         ytitle=textoidl('v_r')

   !p.multi=[3,1,5]
   plot, r, t1d, $
         charsize=2.5, $
         xtitle=textoidl('r'), $
         ytitle=textoidl('T')

   !p.multi=[2,1,5]
   plot, r, vth1d, $
         charsize=2.5, $
         xtitle=textoidl('r'), $
         ytitle=textoidl('v_{th}')

      !p.multi=[1,1,5]
   plot, r, mdot1d, /ylog, $
         charsize=2.5, $
         xtitle=textoidl('r'), $
         ytitle=textoidl('Mdot')
   


end
