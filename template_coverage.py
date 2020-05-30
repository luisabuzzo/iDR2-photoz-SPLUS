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
from IGMabs import *
import synthetic_mags

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
	elif source == "QSO":
   		file = os.path.join(lephare_dev, "sed/{}/{}_SED_extended/{}_MOD.list".format(source,source,source))
   		with open(file) as f:
   			line = f.readline()
   			template = []
   			for line in f:
   				column = line.split()[0]
   				files = os.path.join(lephare_dev, "sed/{}/{}_SED_extended/{}".format(source,source,column))
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

def apply_redshift(source,template,library,redshift,age):
	if library == "BC03":
		wavelength,flux = get_bc03(template,age)
	elif source == "QSO_sdss":
		a = fits.open(template)
		a = a[1].data
		wavelength = 10**(a['loglam'])
		flux = a['flux']
	else:
		a = ascii.read(template)
		wavelength,flux = a['col1'],a['col2']
	q = wavelength*(1+redshift)
	return q,flux

def apply_extinction(source,template,redshift,library,extlaw,ebv):
	wavelength,flux = apply_redshift(source,template,library,redshift,age="all")
	if extlaw == "None":
		new_template = flux
		return wavelength,new_template
	else:
		law = ascii.read(os.path.join(lephare_dev, "ext/{}.dat".format(extlaw)))
		extinction = np.interp(wavelength,law['col1'],law['col2'])
		new_template = flux*10**(-0.4*extinction*ebv)
		return wavelength,new_template

def plot_objects(splus,catalog,photoflag,s2n,magcut):
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
		cat = cat[(cat['PhotoFlag'] == photoflag) & (cat['s2nDet'] > s2n) & (cat['r_auto'] < magcut)]
		print(len(cat))
		mags_ug = cat['uJAVA_auto'] - cat['g_auto']
		mags_gr = cat['g_auto'] - cat['r_auto']
		mags_ri = cat['r_auto'] - cat['i_auto']
		mags_iz = cat['i_auto'] - cat['z_auto']
		return mags_ug,mags_gr,mags_ri,mags_iz

def plot_qso14():	
	data = fits.open(os.path.join(home,'data/zspec/DR14Q_mags.fits'))
	cat = data[1].data
	data.close
	nrand = 10000
	sel = np.int_(np.zeros((nrand)))
	for p in range(0,nrand):
		num = randint(0, len(cat))
		sel[p] = num
	cat = cat[sel]
	cat = cat[cat['petroMag_u'] != 'null']
	redshift = cat['z']
	u = [float(i) for i in cat['petroMag_u']]
	g = [float(i) for i in cat['petroMag_g']]
	r = [float(i) for i in cat['petroMag_r']]
	i = [float(i) for i in cat['petroMag_i']]
	z = [float(i) for i in cat['petroMag_z']]
	mags_ug = [x-y for x,y in zip(u,g)]
	mags_gr = [x-y for x,y in zip(g,r)]
	mags_ri = [x-y for x,y in zip(r,i)]
	mags_iz = [x-y for x,y in zip(i,z)]
	return mags_ug,mags_gr,mags_ri,mags_iz,redshift
			

def plot_temp(source,library,colour,plot):
	template = get_template_files(source,library)

	FilterPath = os.path.join(lephare_dev, "filt/splus2/")
	FilterList = 'filters.lst'
	FilterNames, FilterCurves, FilterType, FilterLen, N_filters, l_eff, w_eff, FilterShortname, FilterEdges = synthetic_mags.ReadFilters(FilterPath,FilterList)
	Mags_ug = []
	Mags_gr = []
	Mags_ri = []
	Mags_iz = []
	if library == "COSMOS":
		# fig = plt.figure(figsize=(12,8))
		redshift = np.arange(0.0,1.5,0.1)
		extlaw = ["None"]#,"SB_calzetti","SMC_prevot","SB_calzetti_bump1","SB_calzetti_bump2"]
		ebv = [0.0]#,0.05,0.1,0.15,0.2,0.25,0.3,0.4,0.5]	
		for h in template:
			if 'Ell' in h:
				for i in redshift:
		 			print(h,i)
		 			w,f = apply_redshift(source=source,template=h,library=library,redshift=i,age="all")
		 			SED_WaveLength = w
		 			SED_Flux = f
	 				Mags, F_nu = synthetic_mags.CalcMags(SED_WaveLength, SED_Flux, N_filters, FilterCurves, FilterLen, FilterType, l_eff, FilterEdges, FilterNames)
	 				Mags_ug.append(Mags[0] - Mags[5])
	 				Mags_gr.append(Mags[5] - Mags[7])
	 				Mags_ri.append(Mags[7] - Mags[9])
	 				Mags_iz.append(Mags[9] - Mags[11])
				if colour == "ug_gr":
					if plot == "color-redshift":
						fig1 = plt.subplot(2,2,1)
						plt.plot(redshift,Mags_ug,c='r',lw=1,linestyle='-',alpha=0.5,zorder=80)
					elif plot == "color-color":
						fig1 = plt.subplot(1,3,1)
						plt.plot(Mags_gr,Mags_ug,c='r',lw=1,linestyle='-',alpha=0.5,zorder=80)
				elif colour == "gr_ri":
	 				if plot == "color-redshift":
		 				fig2 = plt.subplot(2,2,2)
		 				plt.plot(redshift,Mags_gr,c='r',lw=1,linestyle='-',alpha=0.5,zorder=80)
	 				elif plot == "color-color":
		 				fig2 = plt.subplot(1,3,3)
		 				plt.plot(Mags_ri,Mags_gr,c='r',lw=1,linestyle='-',alpha=0.5,zorder=80)

				elif colour == "ri_iz":
	 				if plot == "color-redshift":
		 				fig3 = plt.subplot(2,2,3)
		 				plt.plot(redshift,Mags_ri,c='r',lw=1,linestyle='-',alpha=0.5,zorder=80)
	 				elif plot == "color-color":
		 				fig3 = plt.subplot(1,3,3)
		 				plt.plot(Mags_iz,Mags_ri,c='r',lw=1,linestyle='-',alpha=0.5,zorder=80)

				elif colour == "iz":
	 				fig4 = plt.subplot(2,2,4)
	 				plt.plot(redshift,Mags_iz,c='r',lw=1,linestyle='-',alpha=0.5,zorder=80)
				Mags_ug = []
				Mags_gr = []
				Mags_ri = []
				Mags_iz = []
			else:
				for j in extlaw:
			 		for k in ebv:
			 			for i in redshift:
				 			print(h,i,j,k)
				 			w,f = apply_extinction(source=source,template=h,redshift=i,library="COSMOS",extlaw=j,ebv=k)
				 			SED_WaveLength = w
				 			SED_Flux = f
			 				Mags, F_nu = synthetic_mags.CalcMags(SED_WaveLength, SED_Flux, N_filters, FilterCurves, FilterLen, FilterType, l_eff, FilterEdges, FilterNames)
			 				Mags_ug.append(Mags[0] - Mags[5])
			 				Mags_gr.append(Mags[5] - Mags[7])
			 				Mags_ri.append(Mags[7] - Mags[9])
			 				Mags_iz.append(Mags[9] - Mags[11])
			 			if 'SB' in h:
			 				if colour == "ug_gr":
			 					if plot == "color-redshift":
					 				fig1 = plt.subplot(2,2,1)
					 				plt.plot(redshift,Mags_ug,c='b',lw=1,linestyle='-',alpha=0.5,zorder=80)
				 				elif plot == "color-color":
					 				fig1 = plt.subplot(1,3,1)
					 				plt.plot(Mags_gr,Mags_ug,c='b',lw=1,linestyle='-',alpha=0.5,zorder=80)
					 		elif colour == "gr_ri":
				 				if plot == "color-redshift":
					 				fig2 = plt.subplot(2,2,2)
					 				plt.plot(redshift,Mags_gr,c='b',lw=1,linestyle='-',alpha=0.5,zorder=80)
				 				elif plot == "color-color":
					 				fig2 = plt.subplot(1,3,3)
					 				plt.plot(Mags_ri,Mags_gr,c='b',lw=1,linestyle='-',alpha=0.5,zorder=80)
				 			elif colour == "ri_iz":
				 				if plot == "color-redshift":
					 				fig3 = plt.subplot(2,2,3)
					 				plt.plot(redshift,Mags_ri,c='b',lw=1,linestyle='-',alpha=0.5,zorder=80)
				 				elif plot == "color-color":
					 				fig3 = plt.subplot(1,3,3)
					 				plt.plot(Mags_iz,Mags_ri,c='b',lw=1,linestyle='-',alpha=0.5,zorder=80)
				 			elif colour == "iz":
				 				fig4 = plt.subplot(2,2,4)
				 				plt.plot(redshift,Mags_iz,c='b',lw=1,linestyle='-',alpha=0.5,zorder=80)
			 			else:
			 				if colour == "ug_gr":
			 					if plot == "color-redshift":
					 				fig1 = plt.subplot(2,2,1)
					 				plt.plot(redshift,Mags_ug,c='g',lw=1,linestyle='-',alpha=0.5,zorder=80)
				 				elif plot == "color-color":
					 				fig1 = plt.subplot(1,3,1)
					 				plt.plot(Mags_gr,Mags_ug,c='g',lw=1,linestyle='-',alpha=0.5,zorder=80)
					 		elif colour == "gr_ri":
				 				if plot == "color-redshift":
					 				fig2 = plt.subplot(2,2,2)
					 				plt.plot(redshift,Mags_gr,c='g',lw=1,linestyle='-',alpha=0.5,zorder=80)
				 				elif plot == "color-color":
					 				fig2 = plt.subplot(1,3,3)
					 				plt.plot(Mags_ri,Mags_gr,c='g',lw=1,linestyle='-',alpha=0.5,zorder=80)
				 			elif colour == "ri_iz":
				 				if plot == "color-redshift":
					 				fig3 = plt.subplot(2,2,3)
					 				plt.plot(redshift,Mags_ri,c='g',lw=1,linestyle='-',alpha=0.5,zorder=80)
				 				elif plot == "color-color":
					 				fig3 = plt.subplot(1,3,3)
					 				plt.plot(Mags_iz,Mags_ri,c='g',lw=1,linestyle='-',alpha=0.5,zorder=80)
				 			elif colour == "iz":
				 				fig4 = plt.subplot(2,2,4)
				 				plt.plot(redshift,Mags_iz,c='g',lw=1,linestyle='-',alpha=0.5,zorder=80)
			 			Mags_ug = []
			 			Mags_gr = []
			 			Mags_ri = []
			 			Mags_iz = []

	elif library == "BC03":
		age_wanted = age()
		fig = plt.figure(figsize=(12,8))
		redshift = np.arange(0.0,1.5,0.1)	
		files = glob.glob(os.path.join(home, "lephare_dev/sed/GAL/BC03_CHAB/extracted_*"))
		for h in files:
			for a in age_wanted:
				for i in redshift:
					print(h,i)
					w,f = apply_redshift(source=source,template=h,library=library,redshift=i,age=a)
					SED_WaveLength = w
					SED_Flux = f
					Mags, F_nu = synthetic_mags.CalcMags(SED_WaveLength, SED_Flux, N_filters, FilterCurves, FilterLen, FilterType, l_eff, FilterEdges, FilterNames)
					Mags_ug.append(Mags[0] - Mags[5])
					Mags_gr.append(Mags[5] - Mags[7])
					Mags_ri.append(Mags[7] - Mags[9])
					Mags_iz.append(Mags[9] - Mags[11])
				if colour == "ug_gr":
					if plot == "color-redshift":
		 				fig1 = plt.subplot(2,2,1)
		 				plt.plot(redshift,Mags_ug,c='magenta',lw=1,linestyle='-',alpha=0.5,zorder=80)
					elif plot == "color-color":
		 				fig1 = plt.subplot(1,3,1)
		 				plt.plot(Mags_gr,Mags_ug,c='magenta',lw=1,linestyle='-',alpha=0.5,zorder=80)
				elif colour == "gr_ri":
	 				if plot == "color-redshift":
		 				fig2 = plt.subplot(2,2,2)
		 				plt.plot(redshift,Mags_gr,c='magenta',lw=1,linestyle='-',alpha=0.5,zorder=80)
	 				elif plot == "color-color":
		 				fig2 = plt.subplot(1,3,3)
		 				plt.plot(Mags_ri,Mags_gr,c='magenta',lw=1,linestyle='-',alpha=0.5,zorder=80)
				elif colour == "ri_iz":
	 				if plot == "color-redshift":
		 				fig3 = plt.subplot(2,2,3)
		 				plt.plot(redshift,Mags_ri,c='magenta',lw=1,linestyle='-',alpha=0.5,zorder=80)
	 				elif plot == "color-color":
		 				fig3 = plt.subplot(1,3,3)
		 				plt.plot(Mags_iz,Mags_ri,c='magenta',lw=1,linestyle='-',alpha=0.5,zorder=80)
				elif colour == "iz":
	 				fig4 = plt.subplot(2,2,4)
	 				plt.plot(redshift,Mags_iz,c='magenta',lw=1,linestyle='-',alpha=0.5,zorder=80)
				Mags_ug = []
				Mags_gr = []
				Mags_ri = []
				Mags_iz = []
			line1 = mlines.Line2D([], [], color='pink', linewidth=0.8)
			line2 = mlines.Line2D([], [], color='purple',linestyle= '--', linewidth=1)
			marker = plt.scatter([], [], color='darkorange', marker="*",s=15)
			leg = Legend(fig, handles=[line1,line2,marker], labels=[r'BC03', r'QSO', r'STAR'], bbox_to_anchor=[0.28,0.89],frameon=True)
			fig.add_artist(leg)
	# 
	elif source == "QSO":
		redshift = np.arange(0.0,5.0,0.1)	
		for h in template:
 			for i in redshift:
 				print(h,i)
	 			w,f = apply_redshift(source=source,template=h,library=library,redshift=i,age="all")
	 			SED_WaveLength = w
	 			flux = f
	 			SED_Flux = IGMabs(i,w,flux)
 				Mags, F_nu = synthetic_mags.CalcMags(SED_WaveLength, SED_Flux, N_filters, FilterCurves, FilterLen, FilterType, l_eff, FilterEdges, FilterNames)
 				Mags_ug.append(Mags[0] - Mags[5])
 				Mags_gr.append(Mags[5] - Mags[7])
 				Mags_ri.append(Mags[7] - Mags[9])
 				Mags_iz.append(Mags[9] - Mags[11])
 			if colour == "ug_gr":
 				if plot == "color-redshift":
	 				fig1 = plt.subplot(2,2,1)
	 				plt.plot(redshift,Mags_ug,c='purple',lw=1,linestyle='--',alpha=0.3,zorder=70)
 				elif plot == "color-color":
	 				fig1 = plt.subplot(1,3,1)
	 				plt.plot(Mags_gr,Mags_ug,c='purple',lw=1,linestyle='--',alpha=0.3,zorder=70)
 			elif colour == "gr_ri":
 				if plot == "color-redshift":
	 				fig2 = plt.subplot(2,2,2)
	 				plt.plot(redshift,Mags_gr,c='purple',lw=1,linestyle='--',alpha=0.3,zorder=70)
 				elif plot == "color-color":
	 				fig2 = plt.subplot(1,3,3)
	 				plt.plot(Mags_ri,Mags_gr,c='purple',lw=1,linestyle='--',alpha=0.3,zorder=70)
 			elif colour == "ri_iz":
 				if plot == "color-redshift":
	 				fig3 = plt.subplot(2,2,3)
	 				plt.plot(redshift,Mags_ri,c='purple',lw=1,linestyle='--',alpha=0.3,zorder=70)
 				elif plot == "color-color":
	 				fig3 = plt.subplot(1,3,3)
	 				plt.plot(Mags_iz,Mags_ri,c='purple',lw=1,linestyle='--',alpha=0.3,zorder=70)
 			elif colour == "iz":
 				fig4 = plt.subplot(2,2,4)
 				plt.plot(redshift,Mags_iz,c='purple',lw=1,linestyle='--',alpha=0.3,zorder=70)
 			Mags_ug = []
 			Mags_gr = []
 			Mags_ri = []
 			Mags_iz = []


	elif source == "STAR":
		for h in template:
			print(h)
			w,f = apply_redshift(source=source,template=h,library=library,redshift=0,age="all")
			SED_WaveLength = w
			SED_Flux = f
			Mags, F_nu = synthetic_mags.CalcMags(SED_WaveLength, SED_Flux, N_filters, FilterCurves, FilterLen, FilterType, l_eff, FilterEdges, FilterNames)
			print(Mags)
			Mags_ug = Mags[0] - Mags[5]
			Mags_gr = Mags[5] - Mags[7]
			Mags_ri = Mags[7] - Mags[9]
			Mags_iz = Mags[9] - Mags[11]
			if colour == "ug_gr":
				plt.scatter(Mags_gr,Mags_ug,c='darkorange',s=20,marker='*',alpha=1,zorder=100)
			elif colour == "gr_ri":
				plt.scatter(Mags_ri,Mags_gr,c='darkorange',s=20,marker='*',alpha=1,zorder=100)
			elif colour == "ri_iz":
				plt.scatter(Mags_iz,Mags_ri,c='darkorange',s=20,marker='*',alpha=1,zorder=100)

def main():
	photoflag = 0
	s2n = 10
	magcut = 21
	# splus = 0 # No objects
	splus = 1 # Plot random sample of S-PLUS
	splus = "DR14Q"  # Plot random sample of DR14Q
	# plot = "color-color" # Make the plot of the templates in the color-color space
	plot = "color-redshift" # Make the plot of the templates in color-redshift space
	""" The keyword "library" only differentiates between COSMOS and BC03 templates. For STARs and QSO, use library = "all" """
	# library = "BC03"
	library = "COSMOS"

	if plot == "color-redshift":
		if splus == "DR14Q":
			# to be updated with spectroscopic redshifts of galaxies. For now, it plots the quasars from DR14
			fig = plt.figure(figsize=(12,8))
			DR14Q = plot_qso14()
			colours = ["ug_gr","gr_ri","ri_iz","iz"]
			for colour in colours:
				plot_temp("GAL",library,colour,plot)
				plot_temp("QSO","all",colour,plot)
				if colour == "ug_gr":
					line1 = mlines.Line2D([], [], color='r', linewidth=0.8)
					line2 = mlines.Line2D([], [], color='g', linewidth=0.8)
					line3 = mlines.Line2D([], [], color='b', linewidth=0.8)
					line4 = mlines.Line2D([], [], color='purple',linestyle= '--', linewidth=1)
					marker = plt.scatter([], [], s= 5,edgecolors='gray', facecolors='none')
					leg = Legend(fig, handles=[line1,line2,line3,line4,marker], labels=[r'Ell (COSMOS)', r'Disc (COSMOS)', r'SB (COSMOS)', r'QSO', r'DR14Q'], bbox_to_anchor=[0.4,0.89],frameon=True)
					fig.add_artist(leg)
					plt.xlabel(r'Redshift',fontsize=15)
					plt.scatter(DR14Q[4],DR14Q[0],s=0.3, edgecolors='gray', facecolors='none')
					plt.ylabel(r'u - g',fontsize=15)
					plt.ylim(-1,3)
					plt.xlim(0,5)
					plt.grid()
				if colour == "gr_ri":
					plt.xlabel(r'Redshift',fontsize=15)
					plt.scatter(DR14Q[4],DR14Q[1],s=0.3, edgecolors='gray', facecolors='none')
					plt.ylabel(r'g - r',fontsize=15)
					plt.ylim(-1,3)
					plt.xlim(0,5)
					plt.grid()
				elif colour == "ri_iz":
					plt.xlabel(r'Redshift',fontsize=15)
					plt.scatter(DR14Q[4],DR14Q[2],s=0.3, edgecolors='gray', facecolors='none')
					plt.ylabel(r'r - i',fontsize=15)
					plt.ylim(-1,3)
					plt.xlim(0,5)
					plt.grid()
				elif colour == "iz":
					plt.xlabel(r'Redshift',fontsize=15)
					plt.scatter(DR14Q[4],DR14Q[3],s=0.3, edgecolors='gray', facecolors='none')
					plt.ylabel(r'i - z',fontsize=15)
					plt.ylim(-1,3)
					plt.xlim(0,5)
					plt.grid()
					plt.savefig(os.path.join(resulting_plot_dir,"template_coverage_{}_{}_objects{}.png".format(plot,library,splus)), bbox_inches="tight")
					plt.show()
		if splus == 0:
			colours = ["ug_gr","gr_ri","ri_iz","iz"]
			fig = plt.figure(figsize=(12,8))
			for colour in colours:
				plot_temp("GAL",library,colour,plot)
				plot_temp("QSO","all",colour,plot)
				if colour == "ug_gr":
					plt.xlabel(r'redshift',fontsize=15)
					plt.ylabel(r'u - g',fontsize=15)
					plt.xlim(0,5)
					plt.ylim(-1,6)
					# line1 = mlines.Line2D([], [], color='purple',linestyle= '--', linewidth=1)
					# marker = plt.scatter([], [], edgecolors='gray',facecolors='none',s=1)
					# leg = Legend(fig, handles=[line1,marker], labels=[r'QSO - POLETTA', r'DR14Q'], bbox_to_anchor=[0.205,0.85],frameon=True)
					# fig.add_artist(leg)
					plt.grid()
				elif colour == "gr_ri":
					plt.xlabel(r'redshift',fontsize=15)
					plt.ylabel(r'g - r',fontsize=15)
					plt.xlim(0,5)
					plt.ylim(-1,4)
					plt.grid()
				elif colour == "ri_iz":
					plt.xlabel(r'redshift',fontsize=15)
					plt.ylabel(r'r - i',fontsize=15)
					plt.xlim(0,5)
					plt.ylim(-1,3)
					plt.grid()
				elif colour == "iz":
					plt.xlabel(r'redshift',fontsize=15)
					# plt.scatter(DR14Q[4],DR14Q[3],s=0.3, edgecolors='gray', facecolors='none')
					plt.ylabel(r'i - z',fontsize=15)
					plt.xlim(0,5)
					plt.ylim(-1,2)
					plt.grid()
					plt.savefig(os.path.join(resulting_plot_dir,"template_coverage_{}_{}_objects{}.png".format(plot,library,splus)), bbox_inches="tight")
					plt.show()

	elif plot == "color-color":
		if splus == 1:
			hydra = plot_objects(splus,"HYDRA",photoflag=photoflag,s2n=s2n,magcut=magcut)
			stripe82 = plot_objects(splus,"STRIPE82_SDSS",photoflag=photoflag,s2n=s2n,magcut=magcut)
			splus = plot_objects(splus,"SPLUS",photoflag=photoflag,s2n=s2n,magcut=magcut)
			colours = ["ug_gr","gr_ri","ri_iz"]
			for colour in colours:
				# plot("GAL",library,colour)
				plot("QSO","all",colour)
				plot("STAR","all",colour)
				if colour == "ug_gr":
					ug = plt.scatter(hydra[1],hydra[0],c='gray',s=2,alpha=0.2,zorder=0)
					gr = plt.scatter(stripe82[1],stripe82[0],c='gray',s=2,alpha=0.2,zorder=0)
					ri = plt.scatter(splus[1],splus[0],c='gray',s=2,alpha=0.2,zorder=0)
					plt.xlabel(r'g - r',fontsize=15)
					plt.ylabel(r'u - g',fontsize=15)
					plt.xlim(-1,3)
					plt.ylim(-1,3.5)
					plt.grid()
				if colour == "gr_ri":
					ug = plt.scatter(hydra[2],hydra[1],c='gray',s=2,alpha=0.2,zorder=0)
					gr = plt.scatter(stripe82[2],stripe82[1],c='gray',s=2,alpha=0.2,zorder=0)
					ri = plt.scatter(splus[2],splus[1],c='gray',s=2,alpha=0.2,zorder=0)
					plt.xlabel(r'r - i',fontsize=15)
					plt.ylabel(r'g - r',fontsize=15)
					plt.xlim(-1,3)
					plt.ylim(-1,3)
					plt.grid()
				elif colour == "ri_iz":
					ug = plt.scatter(hydra[3],hydra[2],c='gray',s=2,alpha=0.2,zorder=0)
					gr = plt.scatter(stripe82[3],stripe82[2],c='gray',s=2,alpha=0.2,zorder=0)
					ri = plt.scatter(splus[3],splus[2],c='gray',s=2,alpha=0.2,zorder=0)
					plt.xlabel(r'i - z',fontsize=15)
					plt.ylabel(r'r - i',fontsize=15)
					plt.xlim(-1,3)
					plt.ylim(-1,3)
					plt.grid()
					plt.savefig(os.path.join(resulting_plot_dir,"template_coverage_{}_{}_objects{}.png".format(plot,library,splus)), bbox_inches="tight")
					plt.show()
		if splus == 0:
			fig = plt.figure(figsize=(12,8))
			DR14Q = plot_qso14()
			colours = ["ug_gr","gr_ri","ri_iz"]
			for colour in colours:
				plot("GAL",library,colour)
				plot("QSO","all",colour)
				plot("STAR","all",colour)
				if colour == "ug_gr":
					plt.xlabel(r'g - r',fontsize=15)
					plt.ylabel(r'u - g',fontsize=15)
					plt.xlim(-1,3)
					plt.ylim(-1,3.5)
					plt.grid()
				if colour == "gr_ri":
					plt.xlabel(r'r - i',fontsize=15)
					plt.ylabel(r'g - r',fontsize=15)
					plt.xlim(-1,3)
					plt.ylim(-1,3)
					plt.grid()
				elif colour == "ri_iz":
					plt.xlabel(r'i - z',fontsize=15)
					plt.ylabel(r'r - i',fontsize=15)
					plt.xlim(-1,3)
					plt.ylim(-1,3)
					plt.grid()
					plt.savefig(os.path.join(resulting_plot_dir,"template_coverage_{}_{}_objects{}.png".format(plot,library,splus)), bbox_inches="tight")
					plt.show()
	
if __name__ == '__main__':
	main()
	

		