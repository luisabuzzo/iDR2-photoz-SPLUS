# -*- coding: utf-8 -*-
""" 

Created on 03/04/2020

Author : Maria Luísa Buzzo

"""

import glob
import numpy as np
import os, sys
import pandas as pd
import re
from astropy.io import fits
from astropy.io import ascii

##################################################### DEFINE MAIN DIRECTORIES ######################################################################

if getpass.getuser() == "luisabuzzo":
    home = "/home/luisabuzzo/Work/Master/SPLUS"
elif getpass.getuser() == "mlgbuzzo":
    home = "/home/mlgbuzzo/LePhare/"

resulting_cat_dir = os.path.join(home, "results/catalogues/")
resulting_plot_dir = os.path.join(home,"results/plots/")

os.environ['LEPHAREDIR'] = os.path.join(home, "lephare_dev")
os.environ['LEPHAREWORK'] = os.path.join(home, "lepharework")

####################################################################################################################################################

def get_master_catalogues(catalogue,dataset):
    """ Returns the names of the catalogues """
    if dataset == "iDR2":
    	data_dir = os.path.join(home, "data/iDR2")
    	if catalogue == "SPLUS":
        	final_cat = "SPLUS_{}_master_catalogue_iDR2_december_2019.fits.fz".format(catalogue)
    	elif catalogue == "HYDRA":
        	final_cat = "SPLUS_{}_master_catalogue_iDR2_december_2019.fits.fz".format(catalogue)
    	elif catalogue == "STRIPE82_SDSS":
        	final_cat = "SPLUS_{}_master_catalogue_iDR2_december_2019.fits.fz".format(catalogue)
    	return os.path.join(data_dir, final_cat)
    else:
        raise ValueError("Data set name not defined: {}".format(dataset))

def columns_master_catalogues(catalogue,dataset):
	""" Returns the columns in this catalog """
	if dataset == "iDR2":
		data_dir = os.path.join(home, "data/iDR2")
		f = get_master_catalogues(catalogue,dataset)
		master = pyfits.open(f)
		columns = master[2].columns
	return columns

def read_master_catalogues(catalogue,dataset):
	""" Read master catalogue """
	if dataset == "iDR2":
		f = get_master_catalogues(catalogue,dataset)
		master = fits.open(f)
		data = master[2].data
	return data

def corrected_catalogues(catalogue,mags,dataset):
	""" Replace undetected bands with respective errors and errors with -1 """
	bands = ['uJAVA', 'F378', 'F395', 'F410', 'F430', 'g', 'F515', 'r', 'F660', 'i', 'F861', 'z']
	ME = read_master_catalogues(catalogue,dataset)
	for f in bands:
		ME['{}_{}'.format(f,mags)][ME['{}_{}'.format(f,mags)] == 99.0] = ME['e{}_{}'.format(f,mags)][ME['{}_{}'.format(f,mags)] == 99.0]
		ME['e{}_{}'.format(f,mags)][ME['{}_{}'.format(f,mags)] == ME['e{}_{}'.format(f,mags)]] = -1.0
	return ME

def get_id(catalogue,mags,dataset):
	f = corrected_catalogues(catalogue,mags,dataset)
	ID = f['ID']
	ID_number = []
	if catalogue == "SPLUS":
		for k in ID:
			v = [int(s) for s in re.findall('\\d+',k)]
			ID_number.append(v[2])
	elif catalogue == "HYDRA":
		for k in ID:
			v = [int(s) for s in re.findall('\\d+',k)]
			ID_number.append(v[1])
	elif catalogue == "STRIPE82_SDSS":
		for k in ID:
			v = [int(s) for s in re.findall('\\d+',k)]
			ID_number.append(v[2])
	return ID_number

def calculate_context(catalogue,mags,dataset):
	""" Define missing bands for LePhare to ignore """
	f = corrected_catalogues(catalogue,mags,dataset)
	a = np.array([f['uJAVA_{}'.format(mags)],f['F378_{}'.format(mags)],f['F395_{}'.format(mags)],f['F410_{}'.format(mags)],f['F430_{}'.format(mags)],f['g_{}'.format(mags)],f['F515_{}'.format(mags)],
	f['r_{}'.format(mags)],f['F660_{}'.format(mags)],f['i_{}'.format(mags)],f['F861_{}'.format(mags)],f['z_{}'.format(mags)]])
	a = np.transpose(a)
	aux1 = 2*np.ones((1,12))
	aux2 = np.arange(0,12)
	d = np.where(abs(a)==99,0,np.power(aux1,aux2))
	c = np.ones((12,1))
	context = np.dot(d,c)
	context = [i[0] for i in context]
	return context

def lephare_catalogues(catalogue,mags,dataset):		
	""" Generate catalog in the correct format to run LePhare """
	ID = get_id(catalogue,mags,dataset)
	full = corrected_catalogues(catalogue,mags,dataset)
	context = calculate_context(catalogue,mags,dataset)
	LePhare = [ID, full['uJAVA_{}'.format(mags)], full['euJAVA_{}'.format(mags)], full['F378_{}'.format(mags)], full['eF378_{}'.format(mags)], full['F395_{}'.format(mags)], full['eF395_{}'.format(mags)], full['F410_{}'.format(mags)], full['eF410_{}'.format(mags)], 
	full['F430_{}'.format(mags)], full['eF430_{}'.format(mags)], full['g_{}'.format(mags)], full['eg_{}'.format(mags)], full['F515_{}'.format(mags)], full['eF515_{}'.format(mags)], full['r_{}'.format(mags)], full['er_{}'.format(mags)], full['F660_{}'.format(mags)], full['eF660_{}'.format(mags)],
	full['i_{}'.format(mags)], full['ei_{}'.format(mags)], full['F861_{}'.format(mags)], full['eF861_{}'.format(mags)], full['z_{}'.format(mags)], full['ez_{}'.format(mags)],context, full['FIELD']]
	return LePhare

def save_lephare_catalogue(catalogue,mags,dataset):
	""" Save new catalog"""
	data = lephare_catalogues(catalogue,mags,dataset)
	ascii.write(data, os.path.join(resulting_cat_dir,"LePhare_{}_{}_{}.dat".format(dataset,catalogue,mags)), names=['ID','uJAVA_{}'.format(mags),'euJAVA_{}'.format(mags), 'F378_{}'.format(mags), 'eF378_{}'.format(mags),
		'F395_{}'.format(mags), 'eF395_{}'.format(mags), 'F410_{}'.format(mags), 'eF410_{}'.format(mags), 'F430_{}'.format(mags), 'eF430_{}'.format(mags), 'g_{}'.format(mags), 'eg_{}'.format(mags),
		'F515_{}'.format(mags), 'eF515_{}'.format(mags), 'r_{}'.format(mags), 'er_{}'.format(mags), 'F660_{}'.format(mags), 'eF660_{}'.format(mags), 'i_{}'.format(mags), 'ei_{}'.format(mags),
		'F861_{}'.format(mags), 'eF861_{}'.format(mags), 'z_{}'.format(mags), 'ez_{}'.format(mags),'context','field'], overwrite=True, format="commented_header")
	return


def main():
	save_lephare_catalogue(catalogue="SPLUS", mags="auto",dataset="iDR2")
	save_lephare_catalogue(catalogue="HYDRA", mags="auto",dataset="iDR2")
	save_lephare_catalogue(catalogue="STRIPE82_SDSS", mags="auto",dataset="iDR2")
	return

if __name__ == "__main__":
	main()
		