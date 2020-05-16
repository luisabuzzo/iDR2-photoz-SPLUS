import numpy as np

### make sure lam is the observed-frame wavelength ###
def IGMabs(Redshift, lam, flux):
	h = 6.625*pow(10,-27)   # Planck constant, ergs
	c = 2.999*pow(10,18)     # speed of light, Angstroms s-1

	ZabsA = (lam - 1216)/1216  # Ly-alpha
	ZabsB = (lam - 1026)/1026  # Ly-beta
	ZabsG = (lam - 973)/973    # Ly-gamma
	ZabsD = (lam - 950)/950    # Ly-delta
	ZabsL = (lam - 912)/912    # Ly continuum

	#Equation 12 in Madau (1995):
	tauA = 0.0036*np.power(lam/1216,3.46)

	#Equation 13 in Madau (1995):
	taumetals = 0.0017*np.power(lam/1216,1.68)

	#Equation 15 in Madau (1995):
	tauB = 1.7*pow(10,-3)*np.power(lam/1026,3.46)
	tauG = 1.2*pow(10,-3)*np.power(lam/973,3.46)
	tauD = 9.3*pow(10,-4)*np.power(lam/950,3.46)

	#find out at which observed wavelengths the absorption does not apply
	lam_range_A = np.where(ZabsA >= Redshift)
	lam_range_B = np.where(ZabsB >= Redshift)
	lam_range_G = np.where(ZabsG >= Redshift)
	lam_range_D = np.where(ZabsD >= Redshift)
	lam_range_L = np.where(ZabsL >= Redshift)
	tauA[lam_range_A[0]] = 0
	tauB[lam_range_B[0]] = 0
	tauG[lam_range_G[0]] = 0
	tauD[lam_range_D[0]] = 0
	taumetals[lam_range_A[0]] = 0

	#numeric approximation of Equation 16 in Madau (1995):
	Xc  = 1+ZabsL
	Xem = 1 + Redshift
	t1 = 0.25*np.power(Xc,3.0)*(pow(Xem,0.46)-np.power(Xc,0.46)) 
	t2 = 9.40*np.power(Xc,1.5)*(pow(Xem,0.18)-np.power(Xc,0.18))
	t3 = 0.70*np.power(Xc,3.0)*(np.power(Xc,-1.32)-pow(Xem,-1.32))
	t4 = 0.023*(pow(Xem,1.68)-np.power(Xc,1.68))
	tauL = t1 + t2 - t3 - t4
	tauL[lam_range_L[0]] = 0

	#sum up the contributions from the Lya line blanketing, metal lines, Lyman series line blanketing, and photo-electric absorption
	tau = tauA + taumetals + tauB + tauG + tauD + tauL

	#apply the attenuation to the input spectrum 
	absorbed_flux = flux * np.exp(-tau)
	return absorbed_flux
