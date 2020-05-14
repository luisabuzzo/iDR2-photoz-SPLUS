import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from astropy.io import fits,ascii
import pandas
import os
import os.path
from astropy.convolution import convolve, Gaussian1DKernel
from astropy.table import Table,vstack,hstack
from IGMabs import *

# add the filter transmission curves to an active plot
def PlotFilters(FilterNames, FilterCurves, FilterLen, Tmin, Tmax, col, zorder, linethick):

	for filter in range(0,len(FilterNames)):
		tmp_l = FilterCurves[filter,:,0]
		tmp_f = FilterCurves[filter,:,1]
		tmp_l = tmp_l[0:FilterLen[filter]]
		tmp_f = tmp_f[0:FilterLen[filter]]
		tmp_f = tmp_f / np.max(tmp_f)
		tmp_f = tmp_f * Tmax
		s = np.less(tmp_f, Tmin)
		tmp_f[s] = Tmin
		tmp_f[0] = Tmin
		tmp_f[-1] = Tmin
		plt.plot(tmp_l, tmp_f, '-', linewidth=linethick, color=col[filter], zorder=zorder)

	return 0

# read the filter transmission information given the location of a file that lists the filter files to be used
def ReadFilters(filterpath, filterlist):

	rowc = 0
	N_filters = 0
	FilterFile = pandas.read_table(filterpath + filterlist,names=['filters','type','filtershortname'],header=None,dtype={'filters': np.str, 'type': np.str, 'filtershortname': np.str},sep=' ',skipinitialspace=True,comment='#',skip_blank_lines=True)
	N_filters = FilterFile.shape[0]
	FilterName = np.array(FilterFile['filters'])
	FilterType = np.array(FilterFile['type'])
	FilterShortname = np.array(FilterFile['filtershortname'])	
	FilterLength = np.int_(np.zeros(N_filters))
	l_eff = np.zeros(N_filters) # effective filter wavelength using Gaussian approximation
	w_eff = np.zeros(N_filters) # effective filter width using Gaussian approximation
	FilterEdges = np.zeros((N_filters,2))

	max_FilterLength = 0
	print('Read ', N_filters, ' filter names from the filter file.')

	# find the longest filter. this will be used to define the 2D array that will hold all filters
	for filter in range(0,N_filters):
		filename = FilterName[filter]
		data = ascii.read(filterpath + filename,names=['wavelength','transmission'],converters={'wavelength': np.float, 'transmission': np.float},delimiter="\s")
		FilterLength[filter] = len(data) #data.shape[0]
	max_len_filters = max(FilterLength)	
	Filters = np.zeros((N_filters,max_len_filters,2)) 
	for filter in range(0,N_filters):
		filename = FilterName[filter]
		print('opening filter '+filename+' of type '+FilterType[filter])
		data = ascii.read(filterpath + filename,names=['wavelength','transmission'],converters={'wavelength': np.float, 'transmission': np.float},delimiter="\s")
		for i in range(0,FilterLength[filter]): 
			Filters[filter,i,0] = np.float_(data['wavelength'][i])
			Filters[filter,i,1] = np.float_(data['transmission'][i])
		lam = Filters[filter,0:FilterLength[filter],0]
		flx = Filters[filter,0:FilterLength[filter],1]

		# calculate effective wavelength and effective width
		l_eff[filter] = np.trapz(flx*lam,x=lam) / np.trapz(flx,x=lam)
		w_eff[filter] = 2. * (np.trapz((lam-l_eff[filter])**2 * flx,x=lam) / np.trapz(flx,x=lam))**0.5

		# interpolate the filter on a regular grid of 10,000 points between the 
		# minimum and maximum wavelength of the filter to select the range where the 
		# transmission is higher than some value (0.1%). Why? We want to make sure if the filter 
		# covers the spectrum. Suppose the filter file cuts off at some wavelength, and is filled with 
		# zeros beyond that, it would cover the object spectrum even though it is garbage...
		min_value = 0.001
		npoints = 10000
		lmin = np.min(lam)
		lmax = np.max(lam)
		dl = (lmax-lmin)/npoints
		tmplam = np.arange(npoints,dtype=float) * dl + lmin
		tmplam_Interp = np.interp(tmplam,lam,flx)
		good_range = np.where(np.greater(tmplam_Interp,min_value))[0]
		good_range1 = good_range[0]
		good_range2 = good_range[-1]
		FilterEdges[filter,0] = tmplam[good_range1]
		FilterEdges[filter,1] = tmplam[good_range2]

	return FilterName, Filters, FilterType, FilterLength, N_filters, l_eff, w_eff, FilterShortname, FilterEdges

# calculate the AB magnitudes given a spectrum and a set of filters	
def CalcMags(Wave, Spectrum, N_filters, FilterCurves, FilterLength, FilterType, l_eff, FilterEdges):

	# define output magnitude and flux density arrays	
	Mags = np.zeros(N_filters)
	F_nu = np.zeros(N_filters)

	# method of filter curve interpolation
	#intmode = 'filter' 	# interpolate filter curve on a regular grid of 10,000 points between minimum and max of the filter
	intmode = 'object' 	# interpolate filter curve on the wavelength axis of the object spectrum

	for filter in range(0,N_filters):
		int_upper = 0.
		int_lower = 0.
		TmpFilter = FilterCurves[filter,:,:]
		TmpFilter = TmpFilter[0:FilterLength[filter],:]
		lmin = TmpFilter[0,0]
		lmax = TmpFilter[-1,0]

		# minimum and maximum wavelength of the spectrum
		l1_sed = np.min(Wave)
		l2_sed = np.max(Wave)

		# the object spectrum needs to be defined from l1_filter_good and l2_filter_good, 
		# otherwise we will not be able to calculate the flux in this filter
		if l1_sed > FilterEdges[filter,0] or l2_sed < FilterEdges[filter,1]:
			print('WARNING filter ', filter, ' - The object spectrum does not fully cover \
the wavelength range of the filter (magnitude set to -99)')	
			Mags[filter] = -99.
			F_nu[filter] = -99

		# it is safe to calculate a magnitude for this spectrum and this filter
		else:

			if intmode == 'filter':		# interpolate the SED on the same wavelength grid as the filter
				npoints = 10000
				dl = (lmax-lmin)/npoints
				TmpWave = np.arange(npoints,dtype=float) * dl + lmin
				TmpFilterInterp = np.interp(TmpWave,TmpFilter[:,0],TmpFilter[:,1])
				TmpSpectrum = np.interp(TmpWave,Wave,Spectrum)

			# interpolate the filter on the wavelength grid of the SED and only for the wavelength range of the filter
			if intmode == 'object':		
				TmpFilterInterp = np.interp(Wave,TmpFilter[:,0],TmpFilter[:,1])

			# valid for a photon counting device (ccd) and an input spectrum in erg/s/cm2/A
			if FilterType[filter] == 'ccd':
				flux_integral_nominator = np.trapz(TmpFilterInterp*Spectrum*Wave, x=Wave)
				flux_integral_denominator = cLight * np.trapz(TmpFilterInterp/Wave, x=Wave)

			if flux_integral_denominator != 0 and flux_integral_nominator > 0:
				F_nu[filter] = flux_integral_nominator / flux_integral_denominator
				Mags[filter] = -2.5 * np.log10(F_nu[filter]) - 48.6
				
			else:
				print('WARNING filter ', filter, ' - flux integral contains negative or zero values \
(magnitude set to -99)')
				Mags[filter] = -99.
				F_nu[filter] = -99
	
	return Mags, F_nu

# read the filenames and object types from a list 
def ReadSpecFilenames(filename):

	print("Reading list of input spectra ", filename)
	Specfile = Table.read(filename, format='ascii', names=('Specfile','Spectype'))

	return Specfile

# apply a Gaussian kernel specified by the sigma to an input array and return the convolved spectrum
def ConvolveSpectrum(f, sigma):

	if sigma > 0:
		gauss_kernel = Gaussian1DKernel(stddev=sigma)
		f_smoothed = convolve(f, gauss_kernel, boundary='extend')
	else:
		f_smoothed = f

	return f_smoothed

# extend the wavelength range of a spectrum by padding it with some flux values
def pad_spectrum(lam,flux,lower_pad_lambda=0,lower_pad_flux=0.,upper_pad_lambda=0,upper_pad_flux=0.):

	# extend the spectrum to a lower wavelength range given by lower_padding and fill with zeros
	if lower_pad_lambda > 0:
		wbin = lam[1] - lam[0]
		dif = lam[0] - lower_padding
		nextrabins = np.int(np.ceil(dif / wbin))
		extrabins = np.arange(nextrabins) * wbin
		lam = np.concatenate((extrabins,lam))
		extraflux = np.zeros((nextrabins)) + lower_pad_flux
		flux = np.concatenate((extrabins,flux))

	# extend the spectrum to a higher wavelength range given by upper_padding and fill with zeros
	if upper_pad_lambda > 0:
		wbin = lam[-2] - lam[-1]
		dif = upper_pad_lambda - lam[-1]
		nextrabins = np.int(np.ceil(dif / wbin))
		extrabins = np.arange(nextrabins) * wbin
		lam = np.concatenate((lam,extrabins))
		extraflux = np.zeros((nextrabins)) + upper_pad_flux
		flux = np.concatenate((flux,extrabins))

	return lam, flux

# read an SDSS spectrum from fits file
def Read_SDSS_Spectrum(specpath,specname):

	df = fits.open(specpath + specname)
	spec = df[1].data
	info = df[2].data
	df.close()
	lam = np.array(10. ** spec.loglam)		# assumed to be Angstrom
	flux = np.array(spec.flux * 1.e-17) 	# assumed to be erg/s/cm2/A
	z = info['Z']

	#DataFrame = fits.open(SpecPath + specname)
	#Spec = DataFrame[0].data
	#SpecInfo = DataFrame[0].header
	#SpecInfo2 = DataFrame[1].header
	#Right_Ascension = SpecInfo['PLUG_RA']
	#Declination = SpecInfo['PLUG_DEC']
	
	return lam, flux, z #, Right_Ascension, Declination, Redshift

# TO BE DONE
def Read_Star_Spectrum(specname, SpecPath):

	SpecName = SpecPath + specname
	#print specname
	DataFrame = pandas.read_table(SpecName, sep="\s+",skiprows=3, usecols=[0,1], names = ['Wavelength','Flux'])
	#print DataFrame
	lam = DataFrame['Wavelength'].tolist()
	flux = DataFrame['Flux'].tolist()
	
	return lam, flux

def main():

	# current working directory
	path = os.getcwd() 

	# define other paths
	sed_path = path + '/templates/'
	filter_path = path + '/filters/'

	# define some constants
	global cLight
	GravCst = 4.3e-9 					# Mpc/Msun (km/s)^2		
	cLight = 2.9979e18 					# Angstrom/s
	h = 0.673							# h = H0/100
	Mpc2cm = 3.08568025 * pow(10,24) 	# cm per Mpc
	Lsol = 3.826e33 					# solar luminosity in erg/s
	omega_m = 0.315						# omega_matter
	omega_l = 1 - omega_m 				# omega_lambda in flat universe
	omega_k = 0.0						# curvature 

	# define a fixed set of 12 colors
	col = ['#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff'] #more colors here: [, '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000']

	# list of spectra for which to compute synthetic magnitudes and read the list
	SpecList = 'SDSS_targets.csv'
	SpecData = ReadSpecFilenames(sed_path + SpecList)

	# list of filters to be used and read the filter curves
	FilterList = 'S_PLUS_filterlist.csv'
	FilterNames, FilterCurves, FilterType, FilterLen, N_filters, l_eff, w_eff, FilterShortname, FilterEdges = ReadFilters(filter_path,FilterList)

	# read the VandenBerk2001 template:
	vdbFile = 'vandenberk.csv'
	vdb_template = pandas.read_table(sed_path + vdbFile, names=['lam','flam','flam_unc'],header=None,dtype={'lam': np.float, 'flam': np.float, 'flam_unc': np.float}, sep=',',skipinitialspace=True,comment='#',skip_blank_lines=True)
	vdb_lam = np.array(vdb_template['lam'])
	vdb_flam = np.array(vdb_template['flam'])
	# maybe apply some additional padding? --> use pad_spectrum()

	# read spectra from the input list and compute synthetic magnitudes
	for s in range(0, len(SpecData)):
		spec = SpecData[s]
		print("Processing spectrum {0:s} / {1:s} {2:s} {3:s}".format(str(s+1),str(len(SpecData)),sed_path+spec['Specfile'],spec['Spectype']))
		if spec['Spectype'] == 'QSO':
			lam, flux, z = Read_SDSS_Spectrum(sed_path,spec['Specfile'])
			lamspec_min = np.min(lam)
			lamspec_max = np.max(lam)

		# calculate the magnitudes 
		Mags, F_nu = CalcMags(lam, flux, N_filters, FilterCurves, FilterLen, FilterType, l_eff, FilterEdges)

		# check if it is necessary to extend the spectrum beyond its wavelengths
		shortest_wavelength_filters = np.min(FilterEdges[:,0])
		longest_wavelength_filters = np.max(FilterEdges[:,1])
		extend_lower = 'no'
		extend_higher = 'no'
		if lamspec_min > shortest_wavelength_filters:
			print("The spectrum does not cover the lower end of the filters.")
			extend_lower = 'yes'
		if lamspec_min > shortest_wavelength_filters:
			print("The spectrum does not cover the upper end of the filters.")
			extend_higher = 'yes'

		if extend_lower == 'yes' or extend_higher == 'yes':
		# if the template is a QSO and does not cover the full S-PLUS filter range, 
		# we use the vdBerk template to extend the spectrum in the following way:
		# 1) the vdBerk template is redshifted to the redshift of the QSO
		# 2) IGM absorption is applied to the vdBerk template using the Madau prescription and the redshift
		# 3) synthetic magnitudes of the redshifted, IGM-absorbed vdBerk template are calculated
		# 4) the vdBerk template is scaled to the i-band magnitude of the QSO (this scaling can be improved)
		# 5) we select the part of the scaled vdBerk spectrum that is below (and above?) the QSO spectral range
		# 6) we glue this part to the bottom and top ends of the QSO spectrum
	
			lam_orig = lam
			flux_orig = flux
			if spec['Spectype'] == 'QSO':

				vdb_lam_redshifted = (1.+z) * vdb_lam # redshift the template
				vdb_flam_abs = IGMabs(z, vdb_lam_redshifted, vdb_flam) # apply IGM absorption
				# calculate the filter magnitudes of the vdb01 template
				vdb_Mags, vdb_F_nu = CalcMags(vdb_lam_redshifted, vdb_flam_abs, N_filters, FilterCurves, FilterLen, FilterType, l_eff, FilterEdges)
				# scale the vdb01 template to match the current spectrum in the i-band:
				sel_scaling_band = (FilterNames == 'SPLUS/i.dat')
				magdif = vdb_Mags[sel_scaling_band]-Mags[sel_scaling_band]
				print("Magnitude difference: ", magdif)
				vdb_flam_scaled = vdb_flam_abs * 10.**(0.4*magdif)
				print('Recalculating vdBerk01 template magnitudes after scaling...')
				vdb_Mags, vdb_F_nu = CalcMags(vdb_lam_redshifted, vdb_flam_scaled, N_filters, FilterCurves, FilterLen, FilterType, l_eff, FilterEdges)
				# select the parts for gluing, with an additional buffer (to reduce the effects of bad parts of the spectrum at the extreme ends)
				lam_buffer = 300 #angstroms
				if extend_lower == 'yes':
					mask = np.where(np.greater(lam,lamspec_min+lam_buffer))
					lam = lam[mask]
					flux = flux[mask]
					lamspec_min_new = np.min(lam)
					sel = np.where(np.less(vdb_lam_redshifted,lamspec_min_new))
					lam = np.concatenate((vdb_lam_redshifted[sel],lam))
					flux = np.concatenate((vdb_flam_scaled[sel],flux))
				if extend_higher == 'yes':
					mask = np.where(np.less(lam,lamspec_max-lam_buffer))
					lam = lam[mask]
					flux = flux[mask]
					lamspec_max_new = np.max(lam)
					sel = np.where(np.greater(vdb_lam_redshifted,lamspec_max_new))
					lam = np.concatenate((lam,vdb_lam_redshifted[sel]))
					flux = np.concatenate((flux,vdb_flam_scaled[sel]))

			# our object spectrum should now cover the entire S-PLUS wavelength range and we recalculate the magnitudes
			Mags_corrected, F_nu_corrected = CalcMags(lam, flux, N_filters, FilterCurves, FilterLen, FilterType, l_eff, FilterEdges)

		# convert F_nu to F_lam for plotting purposes
		F_lam = np.zeros(len(Mags))
		for filter in range(len(Mags)):
			F_lam[filter] = cLight * F_nu[filter] / l_eff[filter]**2 
		if extend_lower == 'yes' or extend_higher == 'yes':
			F_lam_corrected = np.zeros(len(Mags_corrected))
			for filter in range(len(Mags_corrected)):
				F_lam_corrected[filter] = cLight * F_nu_corrected[filter] / l_eff[filter]**2 

		plt.figure(figsize=(12,6))
		LamMin = 2500. 
		LamMax = 11000.
		#sel = np.where(np.greater(flux,0.0))[0]
		Fmin = flux[0] # np.min(flux[sel])
		Fmax = 1.1 * np.max(flux)
		title = spec['Spectype'] + ' ' + spec['Specfile'] + ' z = ' + str(z[0])
		plt.title(title)
		plt.xlim([LamMin,LamMax])	
		plt.ylim([Fmin,Fmax])
		plt.yscale('log')	
		plt.xlabel('Wavelength ($\\mathrm{\\AA}$)',fontsize=15)
		plt.ylabel('Flux density (erg s$^{-1}$ cm$^{-2}$ $\\mathrm{\\AA}^{-1}$)',fontsize=15)

		# plot the filter curves
		PlotFilters(FilterNames, FilterCurves, FilterLen, Fmin, 0.8*Fmax, col, 1, 2)

		# plot the original object spectrum (red solid curve)
		plt.plot(lam_orig,ConvolveSpectrum(flux_orig,5),'r-',zorder=2,linewidth=1)

		# plot the synthetic magnitudes from the original spectrum (open circles)
		for p in range(0,len(l_eff)):
			plt.plot(l_eff[p],F_lam[p],marker='o',markerfacecolor='none',markersize=10,color=col[p], zorder=20)

		# plot the minimum and maximum wavelength of the original object spectrum (vertical dashed lines)
		plt.plot([lamspec_min,lamspec_min],[Fmin,Fmax],'r--')
		plt.plot([lamspec_max,lamspec_max],[Fmin,Fmax],'r--')

		if extend_lower == 'yes' or extend_higher == 'yes':

			# plot the minimum and maximum wavelength of the corrected spectrum (vertical dotted lines)
			plt.plot([lamspec_min+lam_buffer,lamspec_min+lam_buffer],[Fmin,Fmax],':',color='0.5')
			plt.plot([lamspec_max-lam_buffer,lamspec_max-lam_buffer],[Fmin,Fmax],':',color='0.5')

			# plot the corrected spectrum (black solid curve)
			plt.plot(lam,ConvolveSpectrum(flux,5),'k-',zorder=5)

			# plot the synthetic magnitudes from the corrected spectrum (filled circles)
			for p in range(0,len(l_eff)):
				plt.plot(l_eff[p],F_lam_corrected[p],marker='o',markersize=10, color=col[p],zorder=20)

			# plot the vdBerk template (grey dotted curve)
			plt.plot(vdb_lam_redshifted,ConvolveSpectrum(vdb_flam_scaled,5),':',color='0.5',zorder=10)

		plt.show()

if __name__ == '__main__':
     main()