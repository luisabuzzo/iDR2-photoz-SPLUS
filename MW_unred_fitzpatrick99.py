import numpy as np
from scipy.interpolate import CubicSpline
import matplotlib.pyplot as plt

def MW_unred_fitzpatrick99(wave,flux,ebv):
	verbose = 'no' #'yes'
	R_V = 3.1
	c2 = -0.824 + 4.171 * (R_V**(-1.0))
	c1 = 2.030 - 3.007 * c2
	x0 = 4.596 #um-1 (bump position)
	gamma = 0.99 #um-1 (bump width)
	c3 = 3.23 #(bump strength)
	c4 = 0.41 #(FUV curvature)
	x = 10000./ wave #Convert to inverse microns 
	curve = x*0. #Create a list of zeros with the same x size

	xcutuv = 10000.0/2700.0
	xspluv = np.zeros(2)
	xspluv[0] = 10000.0/2700.0
	xspluv[1] = 10000.0/2600.0

	#iuv => UV index:
	iuv = (x >= xcutuv)
	N_UV_where = np.where(iuv == True)
	N_UV = len(iuv[N_UV_where]) #number of elements 		
	if verbose == 'yes':
		print('iuv: ', iuv, 'N_UV: ', N_UV)

	#Optical and infrared index:
	iopir = x < xcutuv
	N_comp_where = np.where(iopir == True)
	iopir_xval = wave[N_comp_where]
	N_opir = len(iopir[N_comp_where]) #Number of elements, N_comp=N_opir

	if N_UV > 0.0:
		xuv = np.concatenate((xspluv,np.array(x[iuv])))
	else:
		xuv = xspluv

	if verbose == 'yes':
		print('xuv: ', xuv)

	xuv = np.array(xuv)
	yuv = c1  + c2 * xuv
	yuv = yuv + c3*xuv**2/((xuv**2-x0**2)**2 +(xuv*gamma)**2)
	yuv = yuv + c4*(0.5392*((xuv>5.9)-5.9)**2+0.05644*((xuv>5.9)-5.9)**3)
	yuv = yuv + R_V

	if verbose == 'yes':
		print('yuv: ', yuv)

	#save spline points in UV:		
	yspluv  = yuv[0:2] 
	if verbose == 'yes':
		print('yspluv: ', yspluv)
		
	if N_UV > 0.0:
		curve[iuv] = yuv[2:]
		
	values1 = [26500.0,12200.0,6000.0,5470.0,4670.0,4110.0]
	xsplopir = []
	initial_val = 0.
	xsplopir.append(initial_val)
	for i in values1:
		xsplopir_val = 10000.0/i
		xsplopir.append(xsplopir_val)

	if verbose == 'yes':
		print('xsplopir: ', xsplopir)
	
	#save spline points in IR:
	values2 = [0.0,0.26469,0.82925]
	ysplir = []
	for j in values2:
		ysplir_val = j*R_V/3.1
		ysplir.append(ysplir_val)

	if verbose == 'yes':
		print('ysplir: ', ysplir)

	#save spline points in VIS:
	ysplop = []
	ysplop1 = np.polynomial.polynomial.polyval(R_V, [-4.22809e-01, 1.00270, 2.13572e-04],tensor=True)
	ysplop.append(ysplop1)
	ysplop2 = np.polynomial.polynomial.polyval(R_V, [-5.13540e-02, 1.00216, -7.35778e-05],tensor=True)
	ysplop.append(ysplop2)
	ysplop3 = np.polynomial.polynomial.polyval(R_V, [ 7.00127e-01, 1.00184, -3.32598e-05],tensor=True)
	ysplop.append(ysplop3)
	ysplop4 = np.polynomial.polynomial.polyval(R_V, [ 1.19456, 1.01707, -5.46959e-03, 7.97809e-04, -4.45636e-05],tensor=True)
	ysplop.append(ysplop4)
	#print 'ysplop1: ', ysplop1, 'ysplop2: ', ysplop2, 'ysplop3: ', ysplop3, 'ysplop4: ', ysplop4, 'ysplop: ', ysplop

	#save spline points in VIS and IR together:
	ysplopir = np.concatenate((ysplir, ysplop), axis=0)

	if verbose == 'yes':
		print('ysplopir: ', ysplopir)

	#save spline points in UV, VIS and IR together:
	x_splval = np.concatenate((xsplopir,xspluv), axis=0)
	y_splval = np.concatenate((ysplopir,yspluv), axis=0)

	if verbose == 'yes':
		print('x: ', len(x_splval), 'y: ', len(y_splval))

	if N_opir > 0:
		cs = CubicSpline(x_splval,y_splval) #Create the function to use in the extindtion correction
		#--Plot the extiction law function (MW):

		curve[iopir] = cs(x[iopir])
		
		plotthis = 0
		if plotthis == 1:
			plt.figure(figsize=(6.5, 4))
			plt.plot(x_splval, y_splval, '--r')
			#plt.plot(x_splval, curve(x_splval), 'ok')
			plt.title('Extinction law for MW (Fitzpatrick et al. 1999)')
			plt.xlabel('(1/Wavelength)(um^-1)')
			plt.ylabel('A(Wavelength)/E(B-V)')
			plt.show()
			
	#Now apply extinction correction to input flux vector:			

	#Derive unreddened flux:
	funred = flux * 10.**(0.4*ebv*curve) 

	return funred

#def main():
#
#	lam = np.arange(10000) + 3000
#	flam = lam*0.0 + 1e-15
#	flam_unred01 = fm_unred(lam,flam,0.1)
#	flam_unred02 = fm_unred(lam,flam,0.2)
#	flam_unred03 = fm_unred(lam,flam,0.3)
#	flam_unred06 = fm_unred(lam,flam,0.6)
#	flam_unred10 = fm_unred(lam,flam,1.0)
#
#	plt.figure()
#	plt.xscale('log')
#	plt.yscale('log')
#	plt.plot(lam,flam,'k-')
#	plt.plot(lam,flam_unred01,'k:')
#	plt.plot(lam,flam_unred02,'k:')
#	plt.plot(lam,flam_unred03,'k:')
#	plt.plot(lam,flam_unred06,'k:')
#	plt.plot(lam,flam_unred10,'k:')
#	plt.show()
#
#if __name__ == '__main__':
#	main()



