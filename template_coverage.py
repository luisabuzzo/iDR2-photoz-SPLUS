# -*- coding: utf-8 -*-

import glob
import numpy as np
import os, sys
import pandas
import math
import getpass
import re
from random import randint
from astropy.io import fits
from astropy.io import ascii
import astropy.coordinates as coord
import astropy.units as u
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.lines as mlines
from matplotlib.legend import Legend
import extract_spectra

##################################################### DEFINE MAIN DIRECTORIES ######################################################################

if getpass.getuser() == "luisabuzzo":
    home = "/home/luisabuzzo/Work/Master/SPLUS"
elif getpass.getuser() == "mlgbuzzo":
    home = "/home/mlgbuzzo/LePhare/"

resulting_cat_dir = os.path.join(home, "results/catalogues/")
resulting_plot_dir = os.path.join(home,"results/plots/")

lephare_dev = os.path.join(home, "lephare_dev")
lepharework = os.path.join(home, "lepharework")

####################################################################################################################################################

def get_template_files(source,library):
	if library == "COSMOS":
		file = os.path.join(lephare_dev, "sed/{}/{}_SED/{}_MOD.list".format(source,library,library))
		with open(file) as f:
			line = f.readline()
			template = []
			for line in f:
				column = line.split()[0]
				files = os.path.join(lephare_dev,"sed/{}/{}".format(source,column))
				template.append(files)
			return template
	elif library == "BC03": 
		file = os.path.join(lephare_dev, "sed/{}/{}_CHAB/{}_MOD.list".format(source,library,library))
		with open(file) as f:
   			line = f.readline()
   			template = []
   			for line in f:
   				column = line.split()[0]
   				files = os.path.join(lephare_dev,"sed/{}/{}_CHAB/{}".format(source,library,column))
   				template.append(files)
   			return template
	else:
   		file = os.path.join(lephare_dev, "sed/{}/{}_MOD.list".format(source,source))
   		with open(file) as f:
   			line = f.readline()
   			template = []
   			for line in f:
   				column = line.split()[0]
   				files = os.path.join(lephare_dev, "sed/{}/{}".format(source,column))
   				template.append(files)
   			return template

def age():
	specs_age, specs_flux, specs_wavelength, nwavelengths, nages = extract_spectra.import_spectra(os.path.join(home, "lephare_dev/sed/GAL/BC03_CHAB/extracted_bc2003_lr_m42_chab_tau03_dust00.ised_ASCII")) # for 100 Myr
	a = specs_age
	return a

def get_bc03(template,age):
	specs_age, specs_flux, specs_wavelength, nwavelengths, nages = extract_spectra.import_spectra(template) # for 100 Myr
	a = specs_age
	age_bin = extract_spectra.find_closest_value(specs_age,age)
	age_myr = specs_age[age_bin][0]/1.e6
	lam = specs_wavelength
	flux = (specs_flux[:,age_bin]).flatten()
	return lam,flux

def apply_redshift(template,library,redshift,age):
	if library == "BC03":
		wavelength,flux = get_bc03(template,age)
	else:
		a = ascii.read(template)
		wavelength,flux = a['col1'],a['col2']
	q = wavelength*(1+redshift)
	return q,flux

def apply_extinction(template,redshift,library,extlaw,ebv):
	wavelength,flux = apply_redshift(template,library,redshift,age="all")
	if extlaw == "None":
		new_template = flux
		return wavelength,new_template
	else:
		law = ascii.read(os.path.join(lephare_dev, "ext/{}.dat".format(extlaw)))
		extinction = np.interp(wavelength,law['col1'],law['col2'])
		new_template = flux*10**(-0.4*extinction*ebv)
		return wavelength,new_template

def ReadFilters(File, FilterPath):

	#Imports all filter data
	# print("Importing Filter Data")
	rowc = 0
	N_filters = 0
	FilterFile = pandas.read_table(FilterPath+File,names=['filters','type'],header=None,dtype={'filters': np.str, 'type': np.str},sep=' ')
	N_filters = FilterFile.shape[0]
	FilterName = np.array(FilterFile['filters'])
	FilterType = np.array(FilterFile['type'])
	FilterLength = np.int_(np.zeros(N_filters))
	max_FilterLength = 0
	# print('There are ', N_filters, ' filters in the filter file.')
	# print('Finding longest filter...')
	for f in range(0,N_filters):
		filename = FilterPath+FilterName[f]
		data = pandas.read_table(filename,names=['wavelength','transmission'],header=None,skiprows=1,dtype={'wavelength': np.float, 'transmission': np.float},sep=' ')
		FilterLength[f] = data.shape[0]
	max_len_filters = max(FilterLength)	
	# print('Longest filter has ', max_len_filters, ' wavelength bins')

	Filters = np.empty((N_filters,max_len_filters,2))#Hard coded filter bin number max, max must be 1 < x < 30
	Filters.fill(0)
	
	for f in range(0,N_filters):
		filename = FilterPath+FilterName[f]
		#print 'opening filter '+filename+' of type '+FilterType[f] 
		data = pandas.read_table(filename,names=['wavelength','transmission'],header=None,skiprows=1,dtype={'wavelength': np.float, 'transmission': np.float},sep=' ')
		for i in range(0,FilterLength[f]): 
			Filters[f,i,0] = np.float_(data['wavelength'][i])
			Filters[f,i,1] = np.float_(data['transmission'][i])

	return FilterName,Filters,FilterType,FilterLength,N_filters

def getMags(Wave,Spectrum,N_filters,FilterCurves,FilterLength,FilterType):

	cLight = 2.9979e18

	Mags = np.empty(N_filters)
	Mags.fill(0)
	intmode = 0

	for filter in range(0,N_filters):
		int_upper = 0.
		int_lower = 0.
		TmpFilter = FilterCurves[filter,:,:]
		TmpFilter = TmpFilter[0:FilterLength[filter],:]
		
		lmin = TmpFilter[0,0]
		lmax = max(TmpFilter[:,0])
		l1 = np.where(np.greater(Wave,lmin))[0]
		l2 = np.where(np.less(Wave,lmax))[0]
		if len(l1) != 0 and len(l2) != 0:
	
			if intmode == 0:
		#interpolate filter on a grid of 10,000 points between the minimum and maximum wavelength of the filter
				npoints = 10000
				dl = (lmax-lmin)/npoints
				TmpWave = np.arange(npoints,dtype=float) * dl + lmin
				TmpFilterInterp = np.interp(TmpWave,TmpFilter[:,0],TmpFilter[:,1])
				#interpolate the SED on the same wavelength grid as the filter
				TmpSpectrum = np.interp(TmpWave,Wave,Spectrum)

			if intmode == 1:
		#interpolate filter on the wavelength grid, where the first point of the grid is the first point for which Lambda_SSP > Lambda_Filter_min 
		#and the last point the largest Lambda_SSP for which Lambda_SSP < Lambda_Filter_Max
				TmpWave = Wave[l1[0]:l2[-1]+1] 
				TmpSpectrum = Spectrum[l1[0]:l2[-1]+1]
				TmpFilterInterp = np.interp(TmpWave,TmpFilter[:,0],TmpFilter[:,1])

			if FilterType[filter] == 'ccd':
			#Lower and upper Integral of the magnitude determination function for a photon counting device and an input spctrum in f_lambda in erg/s/cm2/A
				int_upper = np.trapz(TmpFilterInterp*TmpSpectrum*TmpWave,x=TmpWave)
				int_lower = cLight * np.trapz(TmpFilterInterp/TmpWave,x=TmpWave)

			if int_lower != 0 and int_upper != 0:
				Mags[filter] = -2.5 * math.log10(int_upper/int_lower) - 48.6
			else:
				Mags[filter] = 99.
	
		else:
			Mags[filter] = 99.

	return Mags

def plot_objects(splus,catalog):
	if splus == 0:
		print("Do nothing")
	if splus == 1:
		data = fits.open(os.path.join(home, "data/iDR2/SPLUS_{}_master_catalogue_iDR2_december_2019.fits.fz".format(catalog)))
		cat = data[2].data
		data.close
		nrand = 10000
		sel = np.int_(np.zeros((nrand)))
		for p in range(0,nrand):
			num = randint(0, len(cat))
			sel[p] = num
		cat = cat[sel]
		print(len(cat))
		cat = cat[(cat['PhotoFlag']==0) & (cat['s2nDet'] > 10)]
		print(len(cat))
		mags_ug = cat['uJAVA_auto'] - cat['g_auto']
		mags_gr = cat['g_auto'] - cat['r_auto']
		mags_ri = cat['r_auto'] - cat['i_auto']
		mags_iz = cat['i_auto'] - cat['z_auto']
		return mags_ug,mags_gr,mags_ri,mags_iz
			

def plot(source,library,colour):
#read in the filters 
	template = get_template_files(source,library)
	FilterPath = os.path.join(lephare_dev, "filt/splus2/")
	FilterFile = 'filters.lst'
	FilterNames, FilterCurves, FilterType, FilterLen, N_filters = ReadFilters(FilterFile, FilterPath)
	Mags_ug = []
	Mags_gr = []
	Mags_ri = []
	Mags_iz = []
	if library == "COSMOS":
		fig = plt.figure()
		redshift = np.arange(0.0,1.5,0.1)
		extlaw = ["None","SB_calzetti","SMC_prevot","SB_calzetti_bump1","SB_calzetti_bump2"]
		ebv = [0.0,0.05,0.1,0.15,0.2,0.25,0.3,0.4,0.5]	
		for h in template:
			if 'Ell' in h:
				for i in redshift:
		 			print(h,i)
		 			w,f = apply_redshift(template=h,library=library,redshift=i,age="all")
		 			SED_WaveLength = w
		 			SED_Flux = f
	 				Mags = getMags(SED_WaveLength, SED_Flux, N_filters, FilterCurves, FilterLen, FilterType)
	 				Mags_ug.append(Mags[0] - Mags[5])
	 				Mags_gr.append(Mags[5] - Mags[7])
	 				Mags_ri.append(Mags[7] - Mags[9])
	 				Mags_iz.append(Mags[9] - Mags[11])
				if colour == "ug_gr":
					plt.plot(Mags_gr,Mags_ug,c='r',lw=0.8,alpha=0.3,zorder=100)
				elif colour == "gr_ri":
					plt.plot(Mags_ri,Mags_gr,c='r',lw=0.8,alpha=0.3,zorder=100)
				elif colour == "ri_iz":
					plt.plot(Mags_iz,Mags_ri,c='r',lw=0.8,alpha=0.3,zorder=100)
				Mags_ug = []
				Mags_gr = []
				Mags_ri = []
				Mags_iz = []

			else:
				for j in extlaw:
			 		for k in ebv:
			 			for i in redshift:
				 			print(h,i,j,k)
				 			w,f = apply_extinction(template=h,redshift=i,library="COSMOS",extlaw=j,ebv=k)
				 			SED_WaveLength = w
				 			SED_Flux = f
			 				Mags = getMags(SED_WaveLength, SED_Flux, N_filters, FilterCurves, FilterLen, FilterType)
			 				Mags_ug.append(Mags[0] - Mags[5])
			 				Mags_gr.append(Mags[5] - Mags[7])
			 				Mags_ri.append(Mags[7] - Mags[9])
			 				Mags_iz.append(Mags[9] - Mags[11])
			 			if 'SB' in h:
			 				if colour == "ug_gr":
			 					plt.plot(Mags_gr,Mags_ug,c='b',lw=0.8,alpha=0.3,zorder=90)
			 				elif colour == "gr_ri":
			 					plt.plot(Mags_ri,Mags_gr,c='b',lw=0.8,alpha=0.3,zorder=90)
			 				elif colour == "ri_iz":
			 					plt.plot(Mags_iz,Mags_ri,c='b',lw=0.8,alpha=0.3,zorder=90)
			 			else:
			 				if colour == "ug_gr":
			 					plt.plot(Mags_gr,Mags_ug,c='g',lw=0.8,alpha=0.3,zorder=90)
			 				elif colour == "gr_ri":
			 					plt.plot(Mags_ri,Mags_gr,c='g',lw=0.8,alpha=0.3,zorder=90)
			 				elif colour == "ri_iz":
			 					plt.plot(Mags_iz,Mags_ri,c='g',lw=0.8,alpha=0.3,zorder=90)
			 			Mags_ug = []
			 			Mags_gr = []
			 			Mags_ri = []
			 			Mags_iz = []
		line1 = mlines.Line2D([], [], color='r', linewidth=0.8)
		line2 = mlines.Line2D([], [], color='g', linewidth=0.8)
		line3 = mlines.Line2D([], [], color='b', linewidth=0.8)
		marker = plt.scatter([], [], color='orange', marker="*",s=15)
		leg = Legend(fig, handles=[line1,line2,line3,marker], labels=[r'Ell (COSMOS)', r'Disc (COSMOS)', r'SB (COSMOS)', r'STAR'], bbox_to_anchor=[0.4,0.89],frameon=True)
		fig.add_artist(leg)
	elif library == "BC03":
		age_wanted = age()
		fig = plt.figure()
		redshift = np.arange(0.0,1.5,0.1)	
		files = glob.glob(os.path.join(home, "lephare_dev/sed/GAL/BC03_CHAB/extracted_*"))
		for h in files:
			for a in age_wanted:
				for i in redshift:
					print(h,i)
					w,f = apply_redshift(template=h,library=library,redshift=i,age=a)
					SED_WaveLength = w
					SED_Flux = f
					Mags = getMags(SED_WaveLength, SED_Flux, N_filters, FilterCurves, FilterLen, FilterType)
					Mags_ug.append(Mags[0] - Mags[5])
					Mags_gr.append(Mags[5] - Mags[7])
					Mags_ri.append(Mags[7] - Mags[9])
					Mags_iz.append(Mags[9] - Mags[11])
				if colour == "ug_gr":
					plt.plot(Mags_gr,Mags_ug,c='magenta',lw=0.8,alpha=0.3,zorder=90)
				elif colour == "gr_ri":
					plt.plot(Mags_ri,Mags_gr,c='magenta',lw=0.8,alpha=0.3,zorder=90)
				elif colour == "ri_iz":
					plt.plot(Mags_iz,Mags_ri,c='magenta',lw=0.8,alpha=0.3,zorder=90)
				Mags_ug = []
				Mags_gr = []
				Mags_ri = []
				Mags_iz = []
			line1 = mlines.Line2D([], [], color='magenta', linewidth=0.8)
			marker = plt.scatter([], [], color='orange', marker="*",s=15)
			leg = Legend(fig, handles=[line1,marker], labels=[r'BC03', r'STAR'], bbox_to_anchor=[0.28,0.89],frameon=True)
			fig.add_artist(leg)
	# 
	elif source == "QSO":
		redshift = np.arange(0.0,4.0,0.1)	
		for h in template:
 			for i in redshift:
 				print(h,i)
	 			w,f = apply_redshift(template=h,library=library,redshift=i,age="all")
	 			SED_WaveLength = w
	 			SED_Flux = f
 				Mags = getMags(SED_WaveLength, SED_Flux, N_filters, FilterCurves, FilterLen, FilterType)
 				Mags_ug.append(Mags[0] - Mags[5])
 				Mags_gr.append(Mags[5] - Mags[7])
 				Mags_ri.append(Mags[7] - Mags[9])
 				Mags_iz.append(Mags[9] - Mags[11])
 			if colour == "ug_gr":
 				plt.plot(Mags_gr,Mags_ug,c='purple',lw=0.5,alpha=0.3,zorder=100)
 			elif colour == "gr_ri":
 				plt.plot(Mags_ri,Mags_gr,c='purple',lw=0.5,alpha=0.3,zorder=100)
 			elif colour == "ri_iz":
 				plt.plot(Mags_iz,Mags_ri,c='purple',lw=0.5,alpha=0.3,zorder=100)
 			Mags_ug = []
 			Mags_gr = []
 			Mags_ri = []
 			Mags_iz = []

	elif source == "STAR":
		for h in template:
			print(h)
			w,f = apply_redshift(template=h,library=library,redshift=0,age="all")
			SED_WaveLength = w
			SED_Flux = f
			Mags = getMags(SED_WaveLength, SED_Flux, N_filters, FilterCurves, FilterLen, FilterType)
			Mags_ug = Mags[0] - Mags[5]
			Mags_gr = Mags[5] - Mags[7]
			Mags_ri = Mags[7] - Mags[9]
			Mags_iz = Mags[9] - Mags[11]
			if colour == "ug_gr":
				plt.scatter(Mags_gr,Mags_ug,c='orange',s=20,marker='*',alpha=1,zorder=100)
			elif colour == "gr_ri":
				plt.scatter(Mags_ri,Mags_gr,c='orange',s=20,marker='*',alpha=1,zorder=100)
			elif colour == "ri_iz":
				plt.scatter(Mags_iz,Mags_ri,c='orange',s=20,marker='*',alpha=1,zorder=100)

def main(library,splus):
	if splus == 1:
		hydra = plot_objects(splus,"HYDRA")
		stripe82 = plot_objects(splus,"STRIPE82_SDSS")
		splus = plot_objects(splus,"SPLUS")
		colours = ["ug_gr","gr_ri","ri_iz"]
		for colour in colours:
			plot("GAL",library,colour)
			plot("STAR","all",colour)
			if colour == "ug_gr":
				ug = plt.scatter(hydra[1],hydra[0],c='gray',s=2,alpha=0.2,zorder=0)
				gr = plt.scatter(stripe82[1],stripe82[0],c='gray',s=2,alpha=0.2,zorder=0)
				ri = plt.scatter(splus[1],splus[0],c='gray',s=2,alpha=0.2,zorder=0)
				plt.xlabel(r'g-r',fontsize=15)
				plt.ylabel(r'u-g',fontsize=15)
				# plt.legend(handles=[ug,gr,ri], labels=[r'HYDRA', r'STRIPE82',r'MS'])
				plt.xlim(-1,3)
				plt.ylim(-1,3.5)
				plt.grid()
				plt.savefig(os.path.join(resulting_plot_dir,"template_coverage_all_{}_{}.png".format(library,colour)), bbox_inches="tight")
			elif colour == "gr_ri":
				ug = plt.scatter(hydra[2],hydra[1],c='gray',s=2,alpha=0.2,zorder=0)
				gr = plt.scatter(stripe82[2],stripe82[1],c='gray',s=2,alpha=0.2,zorder=0)
				ri = plt.scatter(splus[2],splus[1],c='gray',s=2,alpha=0.2,zorder=0)
				# plt.legend(handles=[ug,gr,ri], labels=[r'HYDRA', r'STRIPE82',r'MS'])
				plt.xlabel(r'r-i',fontsize=15)
				plt.ylabel(r'g-r',fontsize=15)
				plt.xlim(-1,3)
				plt.ylim(-1,3)
				plt.grid()
				plt.savefig(os.path.join(resulting_plot_dir,"template_coverage_all_{}_{}.png".format(library,colour)), bbox_inches="tight")
			elif colour == "ri_iz":
				ug = plt.scatter(hydra[3],hydra[2],c='gray',s=2,alpha=0.2,zorder=0)
				gr = plt.scatter(stripe82[3],stripe82[2],c='gray',s=2,alpha=0.2,zorder=0)
				ri = plt.scatter(splus[3],splus[2],c='gray',s=2,alpha=0.2,zorder=0)
				# plt.legend(handles=[ug,gr,ri], labels=[r'HYDRA', r'STRIPE82',r'MS'])
				plt.xlabel(r'i-z',fontsize=15)
				plt.ylabel(r'r-i',fontsize=15)
				plt.xlim(-1,3)
				plt.ylim(-1,3)
				plt.grid()
				plt.savefig(os.path.join(resulting_plot_dir,"template_coverage_all_{}_{}.png".format(library,colour)), bbox_inches="tight")
				plt.show()
	if splus == 0:
		hydra = plot_objects(splus,"HYDRA")
		stripe82 = plot_objects(splus,"STRIPE82_SDSS")
		splus = plot_objects(splus,"SPLUS")
		colours = ["ug_gr","gr_ri","ri_iz"]
		for colour in colours:
			plot("GAL",library,colour)
			plot("STAR","all",colour)
			if colour == "ug_gr":
				plt.xlabel(r'g-r',fontsize=15)
				plt.ylabel(r'u-g',fontsize=15)
				plt.xlim(-1,2.5)
				plt.ylim(-1,3.5)
				plt.grid()
				plt.savefig(os.path.join(resulting_plot_dir,"template_coverage_noobj_{}_{}.png".format(library,colour)), bbox_inches="tight")
			elif colour == "gr_ri":
				plt.xlabel(r'r-i',fontsize=15)
				plt.ylabel(r'g-r',fontsize=15)
				plt.xlim(-1,2.5)
				plt.ylim(-1,3)
				plt.grid()
				plt.savefig(os.path.join(resulting_plot_dir,"template_coverage_noobj_{}_{}.png".format(library,colour)), bbox_inches="tight")
				plt.show()
			elif colour == "ri_iz":
				plt.xlabel(r'i-z',fontsize=15)
				plt.ylabel(r'r-i',fontsize=15)
				plt.xlim(-1,2.5)
				plt.ylim(-1,3)
				plt.grid()
				plt.savefig(os.path.join(resulting_plot_dir,"template_coverage_noobj_{}_{}.png".format(library,colour)), bbox_inches="tight")
	
if __name__ == '__main__':
	main("COSMOS",1)
	main("BC03",1)
	

		