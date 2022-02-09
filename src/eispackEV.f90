subroutine balanc (nm, n, a, low, igh, scale)  
!
integer :: i, j, k, l, m, n, jj, nm, igh, low, iexc  
doubleprecision a (nm, n), scale (n)  
doubleprecision c, f, g, r, s, b2, radix  
logical :: noconv  
!
!     this subroutine is a translation of the algol procedure balance,
!     num. math. 13, 293-304(1969) by parlett and reinsch.
!     handbook for auto. comp., vol.ii-linear algebra, 315-326(1971).
!
!     this subroutine balances a real matrix and isolates
!     eigenvalues whenever possible.
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        n is the order of the matrix.
!
!        a contains the input matrix to be balanced.
!
!     on output
!
!        a contains the balanced matrix.
!
!        low and igh are two integers such that a(i,j)
!          is equal to zero if
!           (1) i is greater than j and
!           (2) j=1,...,low-1 or i=igh+1,...,n.
!
!        scale contains information determining the
!           permutations and scaling factors used.
!
!     suppose that the principal submatrix in rows low through igh
!     has been balanced, that p(j) denotes the index interchanged
!     with j during the permutation step, and that the elements
!     of the diagonal matrix used are denoted by d(i,j).  then
!        scale(j) = p(j),    for j = 1,...,low-1
!                 = d(j,j),      j = low,...,igh
!                 = p(j)         j = igh+1,...,n.
!     the order in which the interchanges are made is n to igh+1,
!     then 1 to low-1.
!
!     note that 1 is returned for igh if igh is zero formally.
!
!     the algol procedure exc contained in balance appears in
!     balanc  in line.  (note that the algol roles of identifiers
!     k,l have been reversed.)
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated august 1983.
!
!     ------------------------------------------------------------------
!
radix = 16.0d0  
!
b2 = radix * radix  
k = 1  
l = n  
goto 100  
!     .......... in-line procedure for row and
!                column exchange ..........
   20 scale (m) = j  
if (j.eq.m) goto 50  
!
do 30 i = 1, l  
   f = a (i, j)  
   a (i, j) = a (i, m)  
   a (i, m) = f  
   30 end do  
!
do 40 i = k, n  
   f = a (j, i)  
   a (j, i) = a (m, i)  
   a (m, i) = f  
   40 end do  
!
   50 goto (80, 130), iexc  
!     .......... search for rows isolating an eigenvalue
!                and push them down ..........
   80 if (l.eq.1) goto 280  
l = l - 1  
!     .......... for j=l step -1 until 1 do -- ..........
  100 do 120 jj = 1, l  
   j = l + 1 - jj  
!
   do 110 i = 1, l  
      if (i.eq.j) goto 110  
      if (a (j, i) .ne.0.0d0) goto 120  
  110    end do  
!
   m = l  
   iexc = 1  
   goto 20  
  120 end do  
!
goto 140  
!     .......... search for columns isolating an eigenvalue
!                and push them left ..........
  130 k = k + 1  
!
  140 do 170 j = k, l  
!
   do 150 i = k, l  
      if (i.eq.j) goto 150  
      if (a (i, j) .ne.0.0d0) goto 170  
  150    end do  
!
   m = k  
   iexc = 2  
   goto 20  
  170 end do  
!     .......... now balance the submatrix in rows k to l ..........
do 180 i = k, l  
  180 scale (i) = 1.0d0  
!     .......... iterative loop for norm reduction ..........
  190 noconv = .false.  
!
do 270 i = k, l  
   c = 0.0d0  
   r = 0.0d0  
!
   do 200 j = k, l  
      if (j.eq.i) goto 200  
      c = c + dabs (a (j, i) )  
      r = r + dabs (a (i, j) )  
  200    end do  
!     .......... guard against zero c or r due to underflow ..........
   if (c.eq.0.0d0.or.r.eq.0.0d0) goto 270  
   g = r / radix  
   f = 1.0d0  
   s = c + r  
  210    if (c.ge.g) goto 220  
   f = f * radix  
   c = c * b2  
   goto 210  
  220    g = r * radix  
  230    if (c.lt.g) goto 240  
   f = f / radix  
   c = c / b2  
   goto 230  
!     .......... now balance ..........
  240    if ( (c + r) / f.ge.0.95d0 * s) goto 270  
   g = 1.0d0 / f  
   scale (i) = scale (i) * f  
   noconv = .true.  
!
   do 250 j = k, n  
  250    a (i, j) = a (i, j) * g  
!
   do 260 j = 1, l  
  260    a (j, i) = a (j, i) * f  
!
  270 end do  
!
if (noconv) goto 190  
!
  280 low = k  
igh = l  
return  
end subroutine balanc
subroutine balbak (nm, n, low, igh, scale, m, z)  
!
integer :: i, j, k, m, n, ii, nm, igh, low  
doubleprecision scale (n), z (nm, m)  
doubleprecision s  
!
!     this subroutine is a translation of the algol procedure balbak,
!     num. math. 13, 293-304(1969) by parlett and reinsch.
!     handbook for auto. comp., vol.ii-linear algebra, 315-326(1971).
!
!     this subroutine forms the eigenvectors of a real general
!     matrix by back transforming those of the corresponding
!     balanced matrix determined by  balanc.
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        n is the order of the matrix.
!
!        low and igh are integers determined by  balanc.
!
!        scale contains information determining the permutations
!          and scaling factors used by  balanc.
!
!        m is the number of columns of z to be back transformed.
!
!        z contains the real and imaginary parts of the eigen-
!          vectors to be back transformed in its first m columns.
!
!     on output
!
!        z contains the real and imaginary parts of the
!          transformed eigenvectors in its first m columns.
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated august 1983.
!
!     ------------------------------------------------------------------
!
if (m.eq.0) goto 200  
if (igh.eq.low) goto 120  
!
do 110 i = low, igh  
   s = scale (i)  
!     .......... left hand eigenvectors are back transformed
!                if the foregoing statement is replaced by
!                s=1.0d0/scale(i). ..........
   do 100 j = 1, m  
  100    z (i, j) = z (i, j) * s  
!
  110 end do  
!     ......... for i=low-1 step -1 until 1,
!               igh+1 step 1 until n do -- ..........
  120 do 140 ii = 1, n  
   i = ii  
   if (i.ge.low.and.i.le.igh) goto 140  
   if (i.lt.low) i = low - ii  
   k = scale (i)  
   if (k.eq.i) goto 140  
!
   do 130 j = 1, m  
      s = z (i, j)  
      z (i, j) = z (k, j)  
      z (k, j) = s  
  130    end do  
!
  140 end do  
!
  200 return  
end subroutine balbak
subroutine cdiv (ar, ai, br, bi, cr, ci)  
doubleprecision ar, ai, br, bi, cr, ci  
!
!     complex division, (cr,ci) = (ar,ai)/(br,bi)
!
doubleprecision s, ars, ais, brs, bis  
s = dabs (br) + dabs (bi)  
ars = ar / s  
ais = ai / s  
brs = br / s  
bis = bi / s  
s = brs**2 + bis**2  
cr = (ars * brs + ais * bis) / s  
ci = (ais * brs - ars * bis) / s  
return  
end subroutine cdiv
subroutine elmhes (nm, n, low, igh, a, int)  
!
integer :: i, j, m, n, la, nm, igh, kp1, low, mm1, mp1  
doubleprecision a (nm, n)  
doubleprecision x, y  
integer :: int (igh)  
!
!     this subroutine is a translation of the algol procedure elmhes,
!     num. math. 12, 349-368(1968) by martin and wilkinson.
!     handbook for auto. comp., vol.ii-linear algebra, 339-358(1971).
!
!     given a real general matrix, this subroutine
!     reduces a submatrix situated in rows and columns
!     low through igh to upper hessenberg form by
!     stabilized elementary similarity transformations.
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        n is the order of the matrix.
!
!        low and igh are integers determined by the balancing
!          subroutine  balanc.  if  balanc  has not been used,
!          set low=1, igh=n.
!
!        a contains the input matrix.
!
!     on output
!
!        a contains the hessenberg matrix.  the multipliers
!          which were used in the reduction are stored in the
!          remaining triangle under the hessenberg matrix.
!
!        int contains information on the rows and columns
!          interchanged in the reduction.
!          only elements low through igh are used.
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated august 1983.
!
!     ------------------------------------------------------------------
!
la = igh - 1  
kp1 = low + 1  
if (la.lt.kp1) goto 200  
!
do 180 m = kp1, la  
   mm1 = m - 1  
   x = 0.0d0  
   i = m  
!
   do 100 j = m, igh  
      if (dabs (a (j, mm1) ) .le.dabs (x) ) goto 100  
      x = a (j, mm1)  
      i = j  
  100    end do  
!
   int (m) = i  
   if (i.eq.m) goto 130  
!     .......... interchange rows and columns of a ..........
   do 110 j = mm1, n  
      y = a (i, j)  
      a (i, j) = a (m, j)  
      a (m, j) = y  
  110    end do  
!
   do 120 j = 1, igh  
      y = a (j, i)  
      a (j, i) = a (j, m)  
      a (j, m) = y  
  120    end do  
!     .......... end interchange ..........
  130    if (x.eq.0.0d0) goto 180  
   mp1 = m + 1  
!
   do 160 i = mp1, igh  
      y = a (i, mm1)  
      if (y.eq.0.0d0) goto 160  
      y = y / x  
      a (i, mm1) = y  
!
      do 140 j = m, n  
  140       a (i, j) = a (i, j) - y * a (m, j)  
!
      do 150 j = 1, igh  
  150       a (j, m) = a (j, m) + y * a (j, i)  
!
  160    end do  
!
  180 end do  
!
  200 return  
end subroutine elmhes
subroutine eltran (nm, n, low, igh, a, int, z)  
!
integer :: i, j, n, kl, mm, mp, nm, igh, low, mp1  
doubleprecision a (nm, igh), z (nm, n)  
integer :: int (igh)  
!
!     this subroutine is a translation of the algol procedure elmtrans,
!     num. math. 16, 181-204(1970) by peters and wilkinson.
!     handbook for auto. comp., vol.ii-linear algebra, 372-395(1971).
!
!     this subroutine accumulates the stabilized elementary
!     similarity transformations used in the reduction of a
!     real general matrix to upper hessenberg form by  elmhes.
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        n is the order of the matrix.
!
!        low and igh are integers determined by the balancing
!          subroutine  balanc.  if  balanc  has not been used,
!          set low=1, igh=n.
!
!        a contains the multipliers which were used in the
!          reduction by  elmhes  in its lower triangle
!          below the subdiagonal.
!
!        int contains information on the rows and columns
!          interchanged in the reduction by  elmhes.
!          only elements low through igh are used.
!
!     on output
!
!        z contains the transformation matrix produced in the
!          reduction by  elmhes.
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated august 1983.
!
!     ------------------------------------------------------------------
!
!     .......... initialize z to identity matrix ..........
do 80 j = 1, n  
!
   do 60 i = 1, n  
   60    z (i, j) = 0.0d0  
!
   z (j, j) = 1.0d0  
   80 end do  
!
kl = igh - low - 1  
if (kl.lt.1) goto 200  
!     .......... for mp=igh-1 step -1 until low+1 do -- ..........
do 140 mm = 1, kl  
   mp = igh - mm  
   mp1 = mp + 1  
!
   do 100 i = mp1, igh  
  100    z (i, mp) = a (i, mp - 1)  
!
   i = int (mp)  
   if (i.eq.mp) goto 140  
!
   do 130 j = mp, igh  
      z (mp, j) = z (i, j)  
      z (i, j) = 0.0d0  
  130    end do  
!
   z (i, mp) = 1.0d0  
  140 end do  
!
  200 return  
end subroutine eltran
subroutine hqr2 (nm, n, low, igh, h, wr, wi, z, ierr)  
!
integer :: i, j, k, l, m, n, en, ii, jj, ll, mm, na, nm, nn, igh, &
 itn, its, low, mp2, enm2, ierr
doubleprecision h (nm, n), wr (n), wi (n), z (nm, n)  
doubleprecision p, q, r, s, t, w, x, y, ra, sa, vi, vr, zz, norm, &
 tst1, tst2
logical :: notlas  
!
!     this subroutine is a translation of the algol procedure hqr2,
!     num. math. 16, 181-204(1970) by peters and wilkinson.
!     handbook for auto. comp., vol.ii-linear algebra, 372-395(1971).
!
!     this subroutine finds the eigenvalues and eigenvectors
!     of a real upper hessenberg matrix by the qr method.  the
!     eigenvectors of a real general matrix can also be found
!     if  elmhes  and  eltran  or  orthes  and  ortran  have
!     been used to reduce this general matrix to hessenberg form
!     and to accumulate the similarity transformations.
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        n is the order of the matrix.
!
!        low and igh are integers determined by the balancing
!          subroutine  balanc.  if  balanc  has not been used,
!          set low=1, igh=n.
!
!        h contains the upper hessenberg matrix.
!
!        z contains the transformation matrix produced by  eltran
!          after the reduction by  elmhes, or by  ortran  after the
!          reduction by  orthes, if performed.  if the eigenvectors
!          of the hessenberg matrix are desired, z must contain the
!          identity matrix.
!
!     on output
!
!        h has been destroyed.
!
!        wr and wi contain the real and imaginary parts,
!          respectively, of the eigenvalues.  the eigenvalues
!          are unordered except that complex conjugate pairs
!          of values appear consecutively with the eigenvalue
!          having the positive imaginary part first.  if an
!          error exit is made, the eigenvalues should be correct
!          for indices ierr+1,...,n.
!
!        z contains the real and imaginary parts of the eigenvectors.
!          if the i-th eigenvalue is real, the i-th column of z
!          contains its eigenvector.  if the i-th eigenvalue is complex
!          with positive imaginary part, the i-th and (i+1)-th
!          columns of z contain the real and imaginary parts of its
!          eigenvector.  the eigenvectors are unnormalized.  if an
!          error exit is made, none of the eigenvectors has been found.
!
!        ierr is set to
!          zero       for normal return,
!          j          if the limit of 30*n iterations is exhausted
!                     while the j-th eigenvalue is being sought.
!
!     calls cdiv for complex division.
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated august 1983.
!
!     ------------------------------------------------------------------
!
ierr = 0  
norm = 0.0d0  
k = 1  
!     .......... store roots isolated by balanc
!                and compute matrix norm ..........
do 50 i = 1, n  
!
   do 40 j = k, n  
   40    norm = norm + dabs (h (i, j) )  
!
   k = i  
   if (i.ge.low.and.i.le.igh) goto 50  
   wr (i) = h (i, i)  
   wi (i) = 0.0d0  
   50 end do  
!
en = igh  
t = 0.0d0  
itn = 30 * n  
!     .......... search for next eigenvalues ..........
   60 if (en.lt.low) goto 340  
its = 0  
na = en - 1  
enm2 = na - 1  
!     .......... look for single small sub-diagonal element
!                for l=en step -1 until low do -- ..........
   70 do 80 ll = low, en  
   l = en + low - ll  
   if (l.eq.low) goto 100  
   s = dabs (h (l - 1, l - 1) ) + dabs (h (l, l) )  
   if (s.eq.0.0d0) s = norm  
   tst1 = s  
   tst2 = tst1 + dabs (h (l, l - 1) )  
   if (tst2.eq.tst1) goto 100  
   80 end do  
!     .......... form shift ..........
  100 x = h (en, en)  
if (l.eq.en) goto 270  
y = h (na, na)  
w = h (en, na) * h (na, en)  
if (l.eq.na) goto 280  
if (itn.eq.0) goto 1000  
if (its.ne.10.and.its.ne.20) goto 130  
!     .......... form exceptional shift ..........
t = t + x  
!
do 120 i = low, en  
  120 h (i, i) = h (i, i) - x  
!
s = dabs (h (en, na) ) + dabs (h (na, enm2) )  
x = 0.75d0 * s  
y = x  
w = - 0.4375d0 * s * s  
  130 its = its + 1  
itn = itn - 1  
!     .......... look for two consecutive small
!                sub-diagonal elements.
!                for m=en-2 step -1 until l do -- ..........
do 140 mm = l, enm2  
   m = enm2 + l - mm  
   zz = h (m, m)  
   r = x - zz  
   s = y - zz  
   p = (r * s - w) / h (m + 1, m) + h (m, m + 1)  
   q = h (m + 1, m + 1) - zz - r - s  
   r = h (m + 2, m + 1)  
   s = dabs (p) + dabs (q) + dabs (r)  
   p = p / s  
   q = q / s  
   r = r / s  
   if (m.eq.l) goto 150  
   tst1 = dabs (p) * (dabs (h (m - 1, m - 1) ) + dabs (zz) &
    + dabs (h (m + 1, m + 1) ) )
   tst2 = tst1 + dabs (h (m, m - 1) ) * (dabs (q) + dabs (r) )  
   if (tst2.eq.tst1) goto 150  
  140 end do  
!
  150 mp2 = m + 2  
!
do 160 i = mp2, en  
   h (i, i - 2) = 0.0d0  
   if (i.eq.mp2) goto 160  
   h (i, i - 3) = 0.0d0  
  160 end do  
!     .......... double qr step involving rows l to en and
!                columns m to en ..........
do 260 k = m, na  
   notlas = k.ne.na  
   if (k.eq.m) goto 170  
   p = h (k, k - 1)  
   q = h (k + 1, k - 1)  
   r = 0.0d0  
   if (notlas) r = h (k + 2, k - 1)  
   x = dabs (p) + dabs (q) + dabs (r)  
   if (x.eq.0.0d0) goto 260  
   p = p / x  
   q = q / x  
   r = r / x  
  170    s = dsign (dsqrt (p * p + q * q + r * r), p)  
   if (k.eq.m) goto 180  
   h (k, k - 1) = - s * x  
   goto 190  
  180    if (l.ne.m) h (k, k - 1) = - h (k, k - 1)  
  190    p = p + s  
   x = p / s  
   y = q / s  
   zz = r / s  
   q = q / p  
   r = r / p  
   if (notlas) goto 225  
!     .......... row modification ..........
   do 200 j = k, n  
      p = h (k, j) + q * h (k + 1, j)  
      h (k, j) = h (k, j) - p * x  
      h (k + 1, j) = h (k + 1, j) - p * y  
  200    end do  
!
   j = min0 (en, k + 3)  
!     .......... column modification ..........
   do 210 i = 1, j  
      p = x * h (i, k) + y * h (i, k + 1)  
      h (i, k) = h (i, k) - p  
      h (i, k + 1) = h (i, k + 1) - p * q  
  210    end do  
!     .......... accumulate transformations ..........
   do 220 i = low, igh  
      p = x * z (i, k) + y * z (i, k + 1)  
      z (i, k) = z (i, k) - p  
      z (i, k + 1) = z (i, k + 1) - p * q  
  220    end do  
   goto 255  
  225    continue  
!     .......... row modification ..........
   do 230 j = k, n  
      p = h (k, j) + q * h (k + 1, j) + r * h (k + 2, j)  
      h (k, j) = h (k, j) - p * x  
      h (k + 1, j) = h (k + 1, j) - p * y  
      h (k + 2, j) = h (k + 2, j) - p * zz  
  230    end do  
!
   j = min0 (en, k + 3)  
!     .......... column modification ..........
   do 240 i = 1, j  
      p = x * h (i, k) + y * h (i, k + 1) + zz * h (i, k + 2)  
      h (i, k) = h (i, k) - p  
      h (i, k + 1) = h (i, k + 1) - p * q  
      h (i, k + 2) = h (i, k + 2) - p * r  
  240    end do  
!     .......... accumulate transformations ..........
   do 250 i = low, igh  
      p = x * z (i, k) + y * z (i, k + 1) + zz * z (i, k + 2)  
      z (i, k) = z (i, k) - p  
      z (i, k + 1) = z (i, k + 1) - p * q  
      z (i, k + 2) = z (i, k + 2) - p * r  
  250    end do  
  255    continue  
!
  260 end do  
!
goto 70  
!     .......... one root found ..........
  270 h (en, en) = x + t  
wr (en) = h (en, en)  
wi (en) = 0.0d0  
en = na  
goto 60  
!     .......... two roots found ..........
  280 p = (y - x) / 2.0d0  
q = p * p + w  
zz = dsqrt (dabs (q) )  
h (en, en) = x + t  
x = h (en, en)  
h (na, na) = y + t  
if (q.lt.0.0d0) goto 320  
!     .......... real pair ..........
zz = p + dsign (zz, p)  
wr (na) = x + zz  
wr (en) = wr (na)  
if (zz.ne.0.0d0) wr (en) = x - w / zz  
wi (na) = 0.0d0  
wi (en) = 0.0d0  
x = h (en, na)  
s = dabs (x) + dabs (zz)  
p = x / s  
q = zz / s  
r = dsqrt (p * p + q * q)  
p = p / r  
q = q / r  
!     .......... row modification ..........
do 290 j = na, n  
   zz = h (na, j)  
   h (na, j) = q * zz + p * h (en, j)  
   h (en, j) = q * h (en, j) - p * zz  
  290 end do  
!     .......... column modification ..........
do 300 i = 1, en  
   zz = h (i, na)  
   h (i, na) = q * zz + p * h (i, en)  
   h (i, en) = q * h (i, en) - p * zz  
  300 end do  
!     .......... accumulate transformations ..........
do 310 i = low, igh  
   zz = z (i, na)  
   z (i, na) = q * zz + p * z (i, en)  
   z (i, en) = q * z (i, en) - p * zz  
  310 end do  
!
goto 330  
!     .......... complex pair ..........
  320 wr (na) = x + p  
wr (en) = x + p  
wi (na) = zz  
wi (en) = - zz  
  330 en = enm2  
goto 60  
!     .......... all roots found.  backsubstitute to find
!                vectors of upper triangular form ..........
  340 if (norm.eq.0.0d0) goto 1001  
!     .......... for en=n step -1 until 1 do -- ..........
do 800 nn = 1, n  
   en = n + 1 - nn  
   p = wr (en)  
   q = wi (en)  
   na = en - 1  
   if (q) 710, 600, 800  
!     .......... real vector ..........
  600    m = en  
   h (en, en) = 1.0d0  
   if (na.eq.0) goto 800  
!     .......... for i=en-1 step -1 until 1 do -- ..........
   do 700 ii = 1, na  
      i = en - ii  
      w = h (i, i) - p  
      r = 0.0d0  
!
      do 610 j = m, en  
  610       r = r + h (i, j) * h (j, en)  
!
      if (wi (i) .ge.0.0d0) goto 630  
      zz = w  
      s = r  
      goto 700  
  630       m = i  
      if (wi (i) .ne.0.0d0) goto 640  
      t = w  
      if (t.ne.0.0d0) goto 635  
      tst1 = norm  
      t = tst1  
  632       t = 0.01d0 * t  
      tst2 = norm + t  
      if (tst2.gt.tst1) goto 632  
  635       h (i, en) = - r / t  
      goto 680  
!     .......... solve real equations ..........
  640       x = h (i, i + 1)  
      y = h (i + 1, i)  
      q = (wr (i) - p) * (wr (i) - p) + wi (i) * wi (i)  
      t = (x * s - zz * r) / q  
      h (i, en) = t  
      if (dabs (x) .le.dabs (zz) ) goto 650  
      h (i + 1, en) = ( - r - w * t) / x  
      goto 680  
  650       h (i + 1, en) = ( - s - y * t) / zz  
!
!     .......... overflow control ..........
  680       t = dabs (h (i, en) )  
      if (t.eq.0.0d0) goto 700  
      tst1 = t  
      tst2 = tst1 + 1.0d0 / tst1  
      if (tst2.gt.tst1) goto 700  
      do 690 j = i, en  
         h (j, en) = h (j, en) / t  
  690       end do  
!
  700    end do  
!     .......... end real vector ..........
   goto 800  
!     .......... complex vector ..........
  710    m = na  
!     .......... last vector component chosen imaginary so that
!                eigenvector matrix is triangular ..........
   if (dabs (h (en, na) ) .le.dabs (h (na, en) ) ) goto 720  
   h (na, na) = q / h (en, na)  
   h (na, en) = - (h (en, en) - p) / h (en, na)  
   goto 730  
  720    call cdiv (0.0d0, - h (na, en), h (na, na) - p, q, h (na, na), &
    h (na, en) )
  730    h (en, na) = 0.0d0  
   h (en, en) = 1.0d0  
   enm2 = na - 1  
   if (enm2.eq.0) goto 800  
!     .......... for i=en-2 step -1 until 1 do -- ..........
   do 795 ii = 1, enm2  
      i = na - ii  
      w = h (i, i) - p  
      ra = 0.0d0  
      sa = 0.0d0  
!
      do 760 j = m, en  
         ra = ra + h (i, j) * h (j, na)  
         sa = sa + h (i, j) * h (j, en)  
  760       end do  
!
      if (wi (i) .ge.0.0d0) goto 770  
      zz = w  
      r = ra  
      s = sa  
      goto 795  
  770       m = i  
      if (wi (i) .ne.0.0d0) goto 780  
      call cdiv ( - ra, - sa, w, q, h (i, na), h (i, en) )  
      goto 790  
!     .......... solve complex equations ..........
  780       x = h (i, i + 1)  
      y = h (i + 1, i)  
      vr = (wr (i) - p) * (wr (i) - p) + wi (i) * wi (i) - q * q  
      vi = (wr (i) - p) * 2.0d0 * q  
      if (vr.ne.0.0d0.or.vi.ne.0.0d0) goto 784  
      tst1 = norm * (dabs (w) + dabs (q) + dabs (x) + dabs (y) &
       + dabs (zz) )
      vr = tst1  
  783       vr = 0.01d0 * vr  
      tst2 = tst1 + vr  
      if (tst2.gt.tst1) goto 783  
  784       call cdiv (x * r - zz * ra + q * sa, x * s - zz * sa - q * &
       ra, vr, vi, h (i, na), h (i, en) )
      if (dabs (x) .le.dabs (zz) + dabs (q) ) goto 785  
      h (i + 1, na) = ( - ra - w * h (i, na) + q * h (i, en) ) &
       / x
      h (i + 1, en) = ( - sa - w * h (i, en) - q * h (i, na) ) &
       / x
      goto 790  
  785       call cdiv ( - r - y * h (i, na), - s - y * h (i, en), &
       zz, q, h (i + 1, na), h (i + 1, en) )
!
!     .......... overflow control ..........
  790       t = dmax1 (dabs (h (i, na) ), dabs (h (i, en) ) )  
      if (t.eq.0.0d0) goto 795  
      tst1 = t  
      tst2 = tst1 + 1.0d0 / tst1  
      if (tst2.gt.tst1) goto 795  
      do 792 j = i, en  
         h (j, na) = h (j, na) / t  
         h (j, en) = h (j, en) / t  
  792       end do  
!
  795    end do  
!     .......... end complex vector ..........
  800 end do  
!     .......... end back substitution.
!                vectors of isolated roots ..........
do 840 i = 1, n  
   if (i.ge.low.and.i.le.igh) goto 840  
!
   do 820 j = i, n  
  820    z (i, j) = h (i, j)  
!
  840 end do  
!     .......... multiply by transformation matrix to give
!                vectors of original full matrix.
!                for j=n step -1 until low do -- ..........
do 880 jj = low, n  
   j = n + low - jj  
   m = min0 (j, igh)  
!
   do 880 i = low, igh  
      zz = 0.0d0  
!
      do 860 k = low, m  
  860       zz = zz + z (i, k) * h (k, j)  
!
      z (i, j) = zz  
  880 continue  
!
goto 1001  
!     .......... set error -- all eigenvalues have not
!                converged after 30*n iterations ..........
 1000 ierr = en  
 1001 return  
end subroutine hqr2
subroutine hqr (nm, n, low, igh, h, wr, wi, ierr)  
!  RESTORED CORRECT INDICES OF LOOPS (200,210,230,240). (9/29/89 BSG)
!
integer :: i, j, k, l, m, n, en, ll, mm, na, nm, igh, itn, its, &
 low, mp2, enm2, ierr
doubleprecision h (nm, n), wr (n), wi (n)  
doubleprecision p, q, r, s, t, w, x, y, zz, norm, tst1, tst2  
logical :: notlas  
!
!     this subroutine is a translation of the algol procedure hqr,
!     num. math. 14, 219-231(1970) by martin, peters, and wilkinson.
!     handbook for auto. comp., vol.ii-linear algebra, 359-371(1971).
!
!     this subroutine finds the eigenvalues of a real
!     upper hessenberg matrix by the qr method.
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        n is the order of the matrix.
!
!        low and igh are integers determined by the balancing
!          subroutine  balanc.  if  balanc  has not been used,
!          set low=1, igh=n.
!
!        h contains the upper hessenberg matrix.  information about
!          the transformations used in the reduction to hessenberg
!          form by  elmhes  or  orthes, if performed, is stored
!          in the remaining triangle under the hessenberg matrix.
!
!     on output
!
!        h has been destroyed.  therefore, it must be saved
!          before calling  hqr  if subsequent calculation and
!          back transformation of eigenvectors is to be performed.
!
!        wr and wi contain the real and imaginary parts,
!          respectively, of the eigenvalues.  the eigenvalues
!          are unordered except that complex conjugate pairs
!          of values appear consecutively with the eigenvalue
!          having the positive imaginary part first.  if an
!          error exit is made, the eigenvalues should be correct
!          for indices ierr+1,...,n.
!
!        ierr is set to
!          zero       for normal return,
!          j          if the limit of 30*n iterations is exhausted
!                     while the j-th eigenvalue is being sought.
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated september 1989.
!
!     ------------------------------------------------------------------
!
ierr = 0  
norm = 0.0d0  
k = 1  
!     .......... store roots isolated by balanc
!                and compute matrix norm ..........
do 50 i = 1, n  
!
   do 40 j = k, n  
   40    norm = norm + dabs (h (i, j) )  
!
   k = i  
   if (i.ge.low.and.i.le.igh) goto 50  
   wr (i) = h (i, i)  
   wi (i) = 0.0d0  
   50 end do  
!
en = igh  
t = 0.0d0  
itn = 30 * n  
!     .......... search for next eigenvalues ..........
   60 if (en.lt.low) goto 1001  
its = 0  
na = en - 1  
enm2 = na - 1  
!     .......... look for single small sub-diagonal element
!                for l=en step -1 until low do -- ..........
   70 do 80 ll = low, en  
   l = en + low - ll  
   if (l.eq.low) goto 100  
   s = dabs (h (l - 1, l - 1) ) + dabs (h (l, l) )  
   if (s.eq.0.0d0) s = norm  
   tst1 = s  
   tst2 = tst1 + dabs (h (l, l - 1) )  
   if (tst2.eq.tst1) goto 100  
   80 end do  
!     .......... form shift ..........
  100 x = h (en, en)  
if (l.eq.en) goto 270  
y = h (na, na)  
w = h (en, na) * h (na, en)  
if (l.eq.na) goto 280  
if (itn.eq.0) goto 1000  
if (its.ne.10.and.its.ne.20) goto 130  
!     .......... form exceptional shift ..........
t = t + x  
!
do 120 i = low, en  
  120 h (i, i) = h (i, i) - x  
!
s = dabs (h (en, na) ) + dabs (h (na, enm2) )  
x = 0.75d0 * s  
y = x  
w = - 0.4375d0 * s * s  
  130 its = its + 1  
itn = itn - 1  
!     .......... look for two consecutive small
!                sub-diagonal elements.
!                for m=en-2 step -1 until l do -- ..........
do 140 mm = l, enm2  
   m = enm2 + l - mm  
   zz = h (m, m)  
   r = x - zz  
   s = y - zz  
   p = (r * s - w) / h (m + 1, m) + h (m, m + 1)  
   q = h (m + 1, m + 1) - zz - r - s  
   r = h (m + 2, m + 1)  
   s = dabs (p) + dabs (q) + dabs (r)  
   p = p / s  
   q = q / s  
   r = r / s  
   if (m.eq.l) goto 150  
   tst1 = dabs (p) * (dabs (h (m - 1, m - 1) ) + dabs (zz) &
    + dabs (h (m + 1, m + 1) ) )
   tst2 = tst1 + dabs (h (m, m - 1) ) * (dabs (q) + dabs (r) )  
   if (tst2.eq.tst1) goto 150  
  140 end do  
!
  150 mp2 = m + 2  
!
do 160 i = mp2, en  
   h (i, i - 2) = 0.0d0  
   if (i.eq.mp2) goto 160  
   h (i, i - 3) = 0.0d0  
  160 end do  
!     .......... double qr step involving rows l to en and
!                columns m to en ..........
do 260 k = m, na  
   notlas = k.ne.na  
   if (k.eq.m) goto 170  
   p = h (k, k - 1)  
   q = h (k + 1, k - 1)  
   r = 0.0d0  
   if (notlas) r = h (k + 2, k - 1)  
   x = dabs (p) + dabs (q) + dabs (r)  
   if (x.eq.0.0d0) goto 260  
   p = p / x  
   q = q / x  
   r = r / x  
  170    s = dsign (dsqrt (p * p + q * q + r * r), p)  
   if (k.eq.m) goto 180  
   h (k, k - 1) = - s * x  
   goto 190  
  180    if (l.ne.m) h (k, k - 1) = - h (k, k - 1)  
  190    p = p + s  
   x = p / s  
   y = q / s  
   zz = r / s  
   q = q / p  
   r = r / p  
   if (notlas) goto 225  
!     .......... row modification ..........
   do 200 j = k, EN  
      p = h (k, j) + q * h (k + 1, j)  
      h (k, j) = h (k, j) - p * x  
      h (k + 1, j) = h (k + 1, j) - p * y  
  200    end do  
!
   j = min0 (en, k + 3)  
!     .......... column modification ..........
   do 210 i = L, j  
      p = x * h (i, k) + y * h (i, k + 1)  
      h (i, k) = h (i, k) - p  
      h (i, k + 1) = h (i, k + 1) - p * q  
  210    end do  
   goto 255  
  225    continue  
!     .......... row modification ..........
   do 230 j = k, EN  
      p = h (k, j) + q * h (k + 1, j) + r * h (k + 2, j)  
      h (k, j) = h (k, j) - p * x  
      h (k + 1, j) = h (k + 1, j) - p * y  
      h (k + 2, j) = h (k + 2, j) - p * zz  
  230    end do  
!
   j = min0 (en, k + 3)  
!     .......... column modification ..........
   do 240 i = L, j  
      p = x * h (i, k) + y * h (i, k + 1) + zz * h (i, k + 2)  
      h (i, k) = h (i, k) - p  
      h (i, k + 1) = h (i, k + 1) - p * q  
      h (i, k + 2) = h (i, k + 2) - p * r  
  240    end do  
  255    continue  
!
  260 end do  
!
goto 70  
!     .......... one root found ..........
  270 wr (en) = x + t  
wi (en) = 0.0d0  
en = na  
goto 60  
!     .......... two roots found ..........
  280 p = (y - x) / 2.0d0  
q = p * p + w  
zz = dsqrt (dabs (q) )  
x = x + t  
if (q.lt.0.0d0) goto 320  
!     .......... real pair ..........
zz = p + dsign (zz, p)  
wr (na) = x + zz  
wr (en) = wr (na)  
if (zz.ne.0.0d0) wr (en) = x - w / zz  
wi (na) = 0.0d0  
wi (en) = 0.0d0  
goto 330  
!     .......... complex pair ..........
  320 wr (na) = x + p  
wr (en) = x + p  
wi (na) = zz  
wi (en) = - zz  
  330 en = enm2  
goto 60  
!     .......... set error -- all eigenvalues have not
!                converged after 30*n iterations ..........
 1000 ierr = en  
 1001 return  
end subroutine hqr
subroutine rg (nm, n, a, wr, wi, matz, z, iv1, fv1, ierr)  
!
integer :: n, nm, is1, is2, ierr, matz  
doubleprecision a (nm, n), wr (n), wi (n), z (nm, n), fv1 (n)  
integer :: iv1 (n)  
!
!     this subroutine calls the recommended sequence of
!     subroutines from the eigensystem subroutine package (eispack)
!     to find the eigenvalues and eigenvectors (if desired)
!     of a real general matrix.
!
!     on input
!
!        nm  must be set to the row dimension of the two-dimensional
!        array parameters as declared in the calling program
!        dimension statement.
!
!        n  is the order of the matrix  a.
!
!        a  contains the real general matrix.
!
!        matz  is an integer variable set equal to zero if
!        only eigenvalues are desired.  otherwise it is set to
!        any non-zero integer for both eigenvalues and eigenvectors.
!
!     on output
!
!        wr  and  wi  contain the real and imaginary parts,
!        respectively, of the eigenvalues.  complex conjugate
!        pairs of eigenvalues appear consecutively with the
!        eigenvalue having the positive imaginary part first.
!
!        z  contains the real and imaginary parts of the eigenvectors
!        if matz is not zero.  if the j-th eigenvalue is real, the
!        j-th column of  z  contains its eigenvector.  if the j-th
!        eigenvalue is complex with positive imaginary part, the
!        j-th and (j+1)-th columns of  z  contain the real and
!        imaginary parts of its eigenvector.  the conjugate of this
!        vector is the eigenvector for the conjugate eigenvalue.
!
!        ierr  is an integer output variable set equal to an error
!           completion code described in the documentation for hqr
!           and hqr2.  the normal completion code is zero.
!
!        iv1  and  fv1  are temporary storage arrays.
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated august 1983.
!
!     ------------------------------------------------------------------
!
if (n.le.nm) goto 10  
ierr = 10 * n  
goto 50  
!
   10 call balanc (nm, n, a, is1, is2, fv1)  
call elmhes (nm, n, is1, is2, a, iv1)  
if (matz.ne.0) goto 20  
!     .......... find eigenvalues only ..........
call hqr (nm, n, is1, is2, a, wr, wi, ierr)  
goto 50  
!     .......... find both eigenvalues and eigenvectors ..........
   20 call eltran (nm, n, is1, is2, a, iv1, z)  
call hqr2 (nm, n, is1, is2, a, wr, wi, z, ierr)  
if (ierr.ne.0) goto 50  
call balbak (nm, n, is1, is2, fv1, n, z)  
   50 return  
end subroutine rg
