import math
import matplotlib.pyplot as plt
import numpy as np
import getpass
import os
import astropy.units as u
from astropy.coordinates import ICRS, Galactic, FK4, FK5
from astropy.coordinates import SkyCoord
from astropy.coordinates import search_around_sky
from astropy.io import fits,ascii
from matplotlib import rc
import pandas as pd
from astropy.table import Table
from astropy.io import fits,ascii
import collections

if getpass.getuser() == "luisabuzzo":
    home = "/home/luisabuzzo/Work/Master/SPLUS"
elif getpass.getuser() == "mlgbuzzo":
    home = "/home/mlgbuzzo/LePhare/"

def sel_cat(rad,cat_stars,cat_splus,filt,photoflag):

	data = fits.open(cat_stars)
	catalog = data[1].data
	sel = (catalog['PhotoFlag'] == photoflag)
	starcat = catalog[sel]
	coords_stars = SkyCoord(ra=starcat['RA']*u.deg,dec=starcat['Dec']*u.deg, frame=FK5)

	data = fits.open(cat_splus)
	catalog = data[1].data
	sel = (catalog['PhotoFlag'] == photoflag)
	spluscat = catalog[sel]
	coords_splus = SkyCoord(ra=spluscat['RA']*u.deg,dec=spluscat['Dec']*u.deg, frame=FK5)

	idx_1, idx_2, d2d, d3d = coords_splus.search_around_sky(coords_stars,rad*u.arcsec)
	# print(idx_1, idx_2, d2d, d3d)
	mag_star = []
	mag_obj = []
	for i in idx_1:
		mag_star.append(starcat['Pmag'][i])
	for j in idx_2:
		mag_obj.append(spluscat['r_auto'][j])
	table = Table([idx_1,mag_star, idx_2,mag_obj, d2d, d3d],names=['index_starcoord','mag_star','index_spluscoord','mag_obj','distance_2d','distance_3d'])
	return table


def plot_radial_bins2(r,cat_stars,cat_splus,filt,photoflag,magstar_min,magstar_max):
	f = sel_cat(r,cat_stars,cat_splus,filt,photoflag)
	print(len(f))
	# f = ascii.read('numbercounts_astropy_{}arcsec.dat'.format(r))
	bin_mag_obj = np.arange(17,23,1)
	file = f[(f['mag_star'] > magstar_min) & (f['mag_star'] < magstar_max)]
	print(len(file))

	for j in range(6):
		file2 = file[file['mag_obj'] < bin_mag_obj[j]]
		globals()["Counts_magobj" + str(j)] = []
		a = collections.Counter(file2['index_starcoord'])
		for i in file2['index_starcoord']:
			globals()["Counts_magobj" + str(j)].append(a[i])
		d = {'col1': file2['index_spluscoord'], 'col2': globals()["Counts_magobj" + str(j)],'col3': file2['index_starcoord']}
		globals()["data" + str(j)] = pd.DataFrame(data=d)
		
	return globals()["data" + str(0)],globals()["data" + str(1)],globals()["data" + str(2)],globals()["data" + str(3)],globals()["data" + str(4)],globals()["data" + str(5)]

def main():

	catstars = os.path.join(home,'data/iDR2/match_STRIPE82_GSC.fits')
	catfile = os.path.join(home,'data/iDR2/STRIPE82_magtest_broad.fits')
	filt = 'r'
	photoflag = 0

	magstar_min = np.arange(10.0,15.0,1)

	radius = np.arange(4,42,2)
	label = ['mag_obj < 17','mag_obj < 18','mag_obj < 19','mag_obj < 20','mag_obj < 21','mag_obj < 22']
	marker= ['o', 's', '^', 'p', 'X','P']

	for i in magstar_min:
		print(i)
		plt.figure()
		for j in range(6):
			number_dens = []
			for r in radius:
				# plot_radial_bins2(r,catstars,catfile,filt,photoflag,i,i+1)
				df0,df1,df2,df3,df4,df5 = plot_radial_bins2(r,catstars,catfile,filt,photoflag,i,i+1)
				df0_in,df1_in,df2_in,df3_in,df4_in,df5_in = plot_radial_bins2(r-2,catstars,catfile,filt,photoflag,i,i+1)
				print(len(df4))
		
				area = np.pi*(r)**2 - np.pi*(r-2)**2
				df = pd.concat([locals()["df"+str(j)],locals()["df"+str(j)+"_in"]], keys=['col1', 'col2'])
				dropped = df.drop_duplicates(['col1'], keep=False)
				counts = np.median(dropped['col2'])
				number_density = [counts/area]
				number_dens.append(number_density)
				print(r,number_density)
			counts = np.where(np.isnan(number_dens), 0, number_dens)
			plt.plot(radius,counts,'o--',label=label[j],marker=marker[j],lw=0.5,markersize=3)

		plt.xlabel(r'$radius \,\,(arcsec)$')
		plt.ylabel(r'$ number\,\, density$')
		plt.title(r'{} < mag_star < {}'.format(i,i+1))
		plt.legend()
		plt.savefig(os.path.join(home,'results/plots/teststars_{}_photoflag{}_magmin{}.png'.format(filt,photoflag,i)))
	plt.show()

if __name__ == '__main__':
	main()
