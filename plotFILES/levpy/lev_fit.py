import numpy as np
#
def fit_regrc(x,y,sigma):
#
#--------------fit constant function through data points----------------
#---------------------f(x) = a------------------------------------------
#
   nd=np.size(x)
   if nd != np.size(y):
      print('error in fit_regrl: np.size(y) not matching')
      exit()
   if nd != np.size(sigma):
      print('error in fit_regrl: np.size(sigma) not matching')
      exit()

   #calculate all required sums
   s=0.
   sy=0.
   for i in range(0,nd):
      yi=y[i]
      sigmai=sigma[i]
      sy = sy + yi/sigmai**2
      s = s + 1./sigmai**2

   #calculate a
   a = sy/s

   #calculate chi squared
   chi2=0.
   for i in range(0,nd):
      chi2 = chi2 + ((y[i]-a)/sigma[i])**2
      
   print('fit found with')
   print('a', a)
   print('chi2', chi2)
   print('')
   return a
#
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
#
def fit_regrl(x,y,sigma):
#
#--------------linear regression fit for straight line------------------
#---------------------f(x) = a + b*x------------------------------------
#
   nd=np.size(x)
   if nd != np.size(y):
      print('error in fit_regrl: np.size(y) not matching')
      exit()
   if nd != np.size(sigma):
      print('error in fit_regrl: np.size(sigma) not matching')
      exit()

   #calculate all required sums
   s=0.
   sx=0.
   sy=0.
   sxx=0.
   sxy=0.
   for i in range(0,nd):
      xi=x[i]
      yi=y[i]
      sigmai=sigma[i]
      sx = sx + xi/sigmai**2
      sy = sy + yi/sigmai**2
      sxx = sxx + xi**2/sigmai**2
      sxy = sxy + xi*yi/sigmai**2
      s = s + 1./sigmai**2

   #calculate a, b, sigma_a, sigma_b, cov_ab
   delta=sxx*s-sx**2
   a = (sxx*sy-sx*sxy)/delta
   b = (sxy*s-sx*sy)/delta
   sigma_a=sxx/delta
   sigma_b=s/delta
   cov_ab=-sx/delta

   #calculate chi squared
   chi2=0.
   for i in range(0,nd):
      chi2 = chi2 + ((y[i]-a-b*x[i])/sigma[i])**2
      
   print('fit found with')
   print('a', a)
   print('b', b)
   print('sigma_a', sigma_a)
   print('sigma_b', sigma_b)
   print('cov_ab', cov_ab)
   print('chi2', chi2)
   print('')
   return a,b
