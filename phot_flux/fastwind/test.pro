pro test
;
readcol, 'lb1_be/OUT.HALPHA_VT005', indx, fdum1, lambda, fdum2, profile1, profile2
;
;transform lambda to xobs=(nue-nue0)/vth*
cgs_clight=2.9979246d10
vth_fiducial=100.d5             ;100 km/s (in cm/s)
lambda0=6562.8399d0             ;halpha
;lambda0=4861.3330               ;hbeta
nue0=cgs_clight/lambda0/1.d-8
nue=cgs_clight/lambda/1.d-8
ddop_fiducial=nue0*vth_fiducial/cgs_clight
xobs = (nue-nue0)/ddop_fiducial
;
window, 0
plot, xobs, profile1

end
