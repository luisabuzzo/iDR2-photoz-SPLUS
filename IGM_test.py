from IGMabs import *
import numpy as np
import matplotlib.pyplot as plt

z1 = 1.0
z2 = 2.44
z3 = 4.0

# create a flat spectrum at z=0
lam = np.arange(1600)*1.0 + 500.
flx = lam*0.0 + 1.e-15

# redshift to z=z1 and apply IGM absorption
lam1 = lam * (z1 + 1)
flx1 = IGMabs(z1,lam1,flx)

# take the redshifted, absorbed spectrum at z=z1 and move it to z=z2, applying again IGM absorption
lam2 = lam1 * (z2+1) / (z1+1) 
flx2 = IGMabs(z2,lam2,flx1)

# for comparison, do it directly by shifting the z=0 spectrum to z=z2, applying IGM only once 
lam3 = lam * (z2 + 1) 
flx3 = IGMabs(z2,lam3,flx)

# take the redshifted, absorbed spectrum at z=z1 and move it to z=z3, applying again IGM absorption
lam4 = lam1 * (z3 + 1) / (z1+1) 
flx4 = IGMabs(z3,lam4,flx1)

# for comparison, do it directly by shifting the z=0 spectrum to z=z3, applying IGM only once 
lam5 = lam * (z3 + 1)  
flx5 = IGMabs(z3,lam5,flx)

plt.figure()

plt.plot(lam,flx,'b-',zorder=1)		# z=0 spectrum
plt.plot(lam1,flx1,'g-', zorder=2)	# z=z1 spectrum
plt.plot(lam2,flx2,'r-',zorder=3)	# z=z1 spectrum, re-redshifted to z=z2
plt.plot(lam3,flx3,'r:',zorder=4)	# z=0 spectrum, directly redshifted to z=z2
plt.plot(lam4,flx4,'k-',zorder=5)	# z=z1 spectrum, re-redshifted to z=z3
plt.plot(lam5,flx5,'k:',zorder=6)	# z=0 spectrum, directly redshifted to z=z3

plt.yscale('log')

plt.show()
