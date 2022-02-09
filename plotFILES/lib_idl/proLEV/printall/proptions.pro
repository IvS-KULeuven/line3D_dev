;
;-----------------------------------------------------------------------
;-----------------------------------------------------------------------
;-----------------------------------------------------------------------
;
PRO PROPTIONS, OPT_AIT, OPT_NG, OPT_ALO, OPT_INVERT, SPATIAL_GRID, MU_GRID, PHI_GRID
;
PRINT, FORMAT='(A20, I4)', 'OPT_AIT = ', OPT_AIT
   IF(OPT_AIT EQ 0) THEN BEGIN
      PRINT, FORMAT='(A52)', '=> Aitkens-extrapolation is NOT included'
   ENDIF ELSE BEGIN
      PRINT, FORMAT='(A52)', '=> Aitkens-extrapolation is     included'
   ENDELSE
;
PRINT, FORMAT='(A20, I4)', 'OPT_NG = ', OPT_NG
   IF(OPT_NG EQ 0) THEN BEGIN
      PRINT, FORMAT='(A52)', '=> NG-extrapolation is NOT included'
   ENDIF ELSE BEGIN
      PRINT, FORMAT='(A52)', '=> NG-extrapolation is     included'
   ENDELSE
;
PRINT, FORMAT='(A20, I4)', 'OPT_ALO = ', OPT_ALO
   CASE OPT_ALO OF
      0: PRINT, FORMAT='(A52)', '=> Classical Lambda-Iteration is used'
      1: PRINT, FORMAT='(A52)', '=> ALI is used with diagonal ALO'
      2: PRINT, FORMAT='(A52)', '=> ALI is used with nearest neighbour ALO'
   ENDCASE
;
PRINT, FORMAT='(A20, I4)', 'OPT_INVERT = ', OPT_INVERT
   CASE OPT_INVERT OF
      0: PRINT, FORMAT='(A60)', '=> Inversion of ALO is performed by LU-decomposition (LAPACK)'
      1: PRINT, FORMAT='(A60)', '=> Inversion of ALO is performed by JOR-method (w=1.d0)'
      2: PRINT, FORMAT='(A60)', '=> Inversion of AlO is performed by SOR-method (w=1.d0)'
      3: PRINT, FORMAT='(A60)', '=> Inversion of ALO is performed by GMRES-method (SPARSKIT)'
   ENDCASE
;
PRINT, FORMAT='(A20, I4)', 'SPATIAL_GRID = ', SPATIAL_GRID
   CASE SPATIAL_GRID OF
      0: PRINT, FORMAT='(A52)', '=> Equidistant radial steps are used'
      1: PRINT, FORMAT='(A52)', '=> Equidistant velocity steps are used'
      2: PRINT, FORMAT='(A52)', '=> Equidistant tau_thomson steps are used'
      3: PRINT, FORMAT='(A52)', '=> See subroutine GRID1D_FINAL  for details'
      4: PRINT, FORMAT='(A52)', '=> See subroutine GRID1D_FINAL2 for details'
   ENDCASE
;
PRINT, FORMAT='(A20, I4)', 'MU_GRID = ', MU_GRID
   CASE MU_GRID OF
      0: PRINT, FORMAT='(A70)', '=> Integration with weights: equidistant theta steps'
      1: PRINT, FORMAT='(A70)', '=> Integration with weights: angles from 2D-procedure'
      2: PRINT, FORMAT='(A70)', '=> Integration with weights: mustar(r_i) are used as nodes'
      3: PRINT, FORMAT='(A70)', '=> Integration with weights: better resolution in ^q space'
      4: PRINT, FORMAT='(A70)', '=> Gauss-Legendre-Integration: over complete interval [-1,1]'
      5: PRINT, FORMAT='(A70)', '=> Gauss-Chebyshev-Integration: over complete interval [-1,1]'
      6: PRINT, FORMAT='(A70)', '=> Gauss-Legendre-Integration: nodes calculated in each octant'
      7: PRINT, FORMAT='(A70)', '=> Gauss-Chebyshev-Integration: nodes calculated in each octant'
   ENDCASE
;
PRINT, FORMAT='(A20, I4)', 'PHI_GRID = ', PHI_GRID
   CASE PHI_GRID OF
      0: PRINT, FORMAT='(A70)', '=> Integration with weights: equidistant phi steps'
      1: PRINT, FORMAT='(A70)', '=> Integration with weights: del_phi=del_theta'
      2: PRINT, FORMAT='(A70)', '=> Integration with weights: dOmega=const (Lobell & Blomme)'
      4: PRINT, FORMAT='(A70)', '=> Gauss-Legendre-Integration: over complete interval [0,2pi]'
      5: PRINT, FORMAT='(A70)', '=> Gauss-Chebyshev-Integration: over complete interval [0,2pi]'
      6: PRINT, FORMAT='(A70)', '=> Gauss-Legendre-Integration: del_phi=del_theta'
      7: PRINT, FORMAT='(A70)', '=> Gauss-Chebyshev-Integration: del_phi=del_theta'
   ENDCASE

END
