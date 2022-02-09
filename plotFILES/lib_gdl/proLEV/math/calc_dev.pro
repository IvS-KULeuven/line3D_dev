pro calc_dev, x_old, x_new, nd, dev_max
;
;-----------------------------------------------------------------------
;---------calculates max deviation between two iterates/vectors---------
;
;                   dev_max=max((x_old-x_new)/x_new)
;
;   input: x_old, x_new: old/new iterates (vectors)
;          nd: dimension of vectors x_new, x_old
;   output: dev_max: maximum deviation between both vectors
;           x_old is overwritten with the deviation
;
;-----------------------------------------------------------------------
;
;
dev=fltarr(nd)*0.d0
;
dev=0.d0
;
indx = where(x_new ne 0.d0)
dev(indx) = (x_old(indx)-x_new(indx))/x_new(indx)
;
dev_max = max(abs(dev),indx_dev_max)
;
x_old=dev
;
print, 'max(dev) at 1d-grid-point:', indx_dev_max, dev_max
;
end
