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

# define path to the main directory
if getpass.getuser() == "luisabuzzo":
    home = "/home/luisabuzzo/Work/Master/SPLUS"
elif getpass.getuser() == "mlgbuzzo":
    home = "/home/mlgbuzzo/LePhare/"

def check_star_in_border(cat_stars,photoflag,survey,distance_to_border):
	# In this function, we check the star catalogue to see if there are 
	# stars close to other stars or to the borders of the tiles. If yes, we exclude such stars.
	# Important: this will be dependent on the maximum distance we are analysing from the star.

    splus_tiles_file = os.path.join('./data/all_pointings.csv')
    splus_tiles = Table.read(splus_tiles_file, format="csv")
    tilesize = 1.4 # degrees
    htilesize = tilesize / 2.
    
    data = fits.open(cat_stars)
    catalog = data[1].data
    sel = (catalog['PhotoFlag'] == photoflag)
    starcat = catalog[sel]

    coords_stars = SkyCoord(ra=starcat['RA']*u.deg,dec=starcat['Dec']*u.deg, frame=FK5)
    # First step: check if there are stars close to other stars, we do this with the function search_around_sky. 
    # The minimum distance allowed is defined by the parameter "distance_to_border".
    idx_stars, idx_stars2, d2d, d3d = coords_stars.search_around_sky(coords_stars,distance_to_border*u.arcsec)
    # in the output of search_around_sky, 
    # each star will appear at least once (the star is close to itself). 
    #If a star appears more than once, it is because it's close to other. 
    stars_near_stars = [x for x, y in collections.Counter(idx_stars).items() if y > 1] # with this function, we count how many stars appear more than once and save their index.
    
    # Below we check if the stars are close to any border of the SPLUS tiles.
    idx = np.where(splus_tiles['PID'] == survey)
    tiledata = splus_tiles[idx]
    linspace_npnts = 5
    coords = SkyCoord(ra=tiledata['RA'], dec=tiledata['DEC'], unit=(u.hourangle,u.degree))
    stars_near_borders = np.array([])
    for coord in coords:
    	# Define the coordinates of the borders of each tile.
    	l1 = coord.ra.value - htilesize
    	l2 = coord.ra.value + htilesize
    	l3 = l2
    	l4 = l1
    	b1 = coord.dec.value - htilesize
    	b2 = b1
    	b3 = coord.dec.value + htilesize
    	b4 = b3
    	seg1_l = np.linspace(l1,l2,linspace_npnts) 
    	seg1_b = np.linspace(b1,b2,linspace_npnts) 
    	seg2_l = np.linspace(l2,l3,linspace_npnts) 
    	seg2_b = np.linspace(b2,b3,linspace_npnts) 
    	seg3_l = np.linspace(l3,l4,linspace_npnts) 
    	seg3_b = np.linspace(b3,b4,linspace_npnts) 
    	seg4_l = np.linspace(l4,l1,linspace_npnts) 
    	seg4_b = np.linspace(b4,b1,linspace_npnts) 
    	coords_1 = SkyCoord(ra=seg1_l, dec=seg1_b, unit=(u.deg,u.deg))
    	coords_2 = SkyCoord(ra=seg2_l, dec=seg2_b, unit=(u.deg,u.deg))
    	coords_3 = SkyCoord(ra=seg3_l, dec=seg3_b, unit=(u.deg,u.deg))
    	coords_4 = SkyCoord(ra=seg4_l, dec=seg4_b, unit=(u.deg,u.deg))
    	# Check if there are stars close to any of the four borders.
    	idx_stars1, idx_coords1, d2d, d3d = coords_1.search_around_sky(coords_stars,distance_to_border*u.arcsec)
    	idx_stars2, idx_coords2, d2d, d3d = coords_2.search_around_sky(coords_stars,distance_to_border*u.arcsec)
    	idx_stars3, idx_coords3, d2d, d3d = coords_3.search_around_sky(coords_stars,distance_to_border*u.arcsec)
    	idx_stars4, idx_coords4, d2d, d3d = coords_4.search_around_sky(coords_stars,distance_to_border*u.arcsec)
    	# Below we append the stars that were identified to be close to any border.
    	stars_near_borders = np.concatenate([stars_near_borders, idx_stars1])
    	stars_near_borders = np.concatenate([stars_near_borders, idx_stars2])
    	stars_near_borders = np.concatenate([stars_near_borders, idx_stars3])
    	stars_near_borders = np.concatenate([stars_near_borders, idx_stars4])

    # Append the stars near stars and stars near borders.
    stars_near_something = np.concatenate([stars_near_stars,stars_near_borders])
    # Be sure that we are counting the star only once.
    stars_near_something = [x for x, y in collections.Counter(stars_near_something).items() if y > 1]

    # Exclude the stars that are near others/borders.
    coords_stars = np.delete(starcat,stars_near_something,axis=0)

    return coords_stars


def sel_cat(rad,cat_stars,cat_splus,filt,photoflag):
	# Input catalog is the one already excluding stars near others/borders.
	coords_stars = SkyCoord(ra=cat_stars['RA']*u.deg,dec=cat_stars['Dec']*u.deg, frame=FK5)

	data = fits.open(cat_splus)
	catalog = data[1].data
	sel = (catalog['PhotoFlag'] == photoflag)
	spluscat = catalog[sel]
	# Coordinates of the survey.
	coords_splus = SkyCoord(ra=spluscat['RA']*u.deg,dec=spluscat['Dec']*u.deg, frame=FK5)

	# Here we check how many objects from the SPLUS catalogs are around the stars in an increasing radius.
	idx_stars, idx_splus, d2d, d3d = coords_splus.search_around_sky(coords_stars,rad*u.arcsec)
	
	mag_star = []
	mag_obj = []

	for i in idx_stars:
		# Get the magnitudes of the stars.
		mag_star.append(cat_stars['Pmag'][i])
	for j in idx_splus:
		# Get the magnitudes of the SPLUS objects.
		mag_obj.append(spluscat['{}'.format(filt)][j])
	# Create a table including the index and magnitude of the stars and the index and magnitudes of all objects around each star.
	table = Table([idx_stars,mag_star, idx_splus,mag_obj, d2d, d3d],names=['index_starcoord','mag_star','index_spluscoord','mag_obj','distance_2d','distance_3d'])
	d = {'index_starcoord': idx_stars, 'mag_star': mag_star,'index_spluscoord': idx_splus, 'mag_obj':mag_obj}
	table = pd.DataFrame(data=d)

	return table


def plot_radial_bins2(f,magstar_min,magstar_max):
	# The analysis is done in bins of magnitude of the objects.
	bin_mag_obj = np.arange(17,23,1)
	# From the input file, we divide the file according to bins of magnitude of the stars.
	file = f[(f['mag_star'] > magstar_min) & (f['mag_star'] < magstar_max)]

	for j in range(6):
		# here we divide again the file, now in bins of magnitude of the objects.
		file2 = file[file['mag_obj'] < bin_mag_obj[j]]
		globals()["Counts_magobj" + str(j)] = []
		# Counts how many times each star appears. This represents the number of objects around the same star.
		a = collections.Counter(file2['index_starcoord'])
		for i in file2['index_starcoord']:
			# Here we save the number of counts respective to each star.
			globals()["Counts_magobj" + str(j)].append(a[i])
		# Create a table for each magnitude bin of the objects, this table includes the index of the object, the number of counts around the star and the index of the star.
		# Note: the intention of this table is to save the index of each object, so that we can control the objects that are only in the ring and not in the full radius,
		# so the number of counts and index of stars will appear more than once.
		# For example, if I have 2 objects (index 1 and 2) around a star, the number of counts will be 2 and the index of the star is "star1", the final table will be:
		# Idx_SPLUS Counts Idx_Star
		#    1         2     star1
		#    2         2     star1 
		d = {'col1': file2['index_spluscoord'], 'col2': globals()["Counts_magobj" + str(j)],'col3': file2['index_starcoord']}
		globals()["data" + str(j)] = pd.DataFrame(data=d)
		# here we return the table for each magnitude bin of the objects.
	return globals()["data" + str(0)],globals()["data" + str(1)],globals()["data" + str(2)],globals()["data" + str(3)],globals()["data" + str(4)],globals()["data" + str(5)]

def main():

	catstars = os.path.join(home,'data/iDR2/match_STRIPE82_GSC.fits')
	catfile = os.path.join(home,'data/iDR2/STRIPE82_magtest_broad.fits')
	filt = 'r_auto'
	photoflag = 0
	radius_steps = 10
	survey = 'STRIPE82'
	distance_to_border = 300 # arcsec

	magstar_min = np.arange(10.0,15.0,1)
	#minimum and maximum radius, I don't start in zero to don't count the star itself.
	radius = np.arange(2,314,radius_steps)
	label = ['mag_obj < 17','mag_obj < 18','mag_obj < 19','mag_obj < 20','mag_obj < 21','mag_obj < 22']
	marker= ['o', 's', '^', 'p', 'X','P']

	# Input the catalog of stars and check if the stars are close to others/borders of S-PLUS tiles.
	starcat = check_star_in_border(catstars,photoflag,survey,distance_to_border)

	for r in radius:
		# Calculate the number of objects around each star in the radius range defined above.
		print("Calculating objects around stars in a radius of: {} arcsec".format(r))
		# return the table defined in the function sel_cat
		globals()["Table_r" + str(r)] = sel_cat(r,starcat,catfile,filt,photoflag)

		# separate the data in bins of magnitude of the stars
	for i in magstar_min:
		plt.figure()
		#separate the data in bins of magnitude of the objects
		for j in range(6):
			number_dens = []
			for r in radius:
				if r <=302:
					# f_big and f_small are the objects found in the bigger radius (r+radius_steps) and in the smaller radius around the star.
					f_big = globals()["Table_r" + str(r+radius_steps)]
					f_small = globals()["Table_r" + str(r)]
					# Concatenate both tables (bigger and smaller)
					df = pd.concat([f_big,f_small])
					# Exclude objects that appear more than once, so that we end up only with the objects in the ring. 
					dropped = df.drop_duplicates(subset=['index_spluscoord'], keep=False)
					#Get the table defined in function plot_radial_bins2 for all magnitude bins of the objects.
					df0,df1,df2,df3,df4,df5 = plot_radial_bins2(dropped,i,i+1)
					# Calculate area of the ring.
					area = np.pi*((r+radius_steps)**2 - (r)**2)
					# Calculate median counts in each radius for each bin of magnitude of the stars.
					counts = np.median(locals()["df" + str(j)]['col2'])
					number_density = [counts/area]
					number_dens.append(number_density)
			# In some cases, the objects in the bigger radius are the exact same as the smaller, so when we exclude the duplicates, we end up with an empty array.
			# In that case, we substitute the resulting NaN with a zero.
			counts = np.where(np.isnan(number_dens), 0, number_dens)
			# plot radius vs. counts for each magnitude bin of the objects.
			plt.plot(radius[0:31],counts,'o--',label=label[j],marker=marker[j],lw=0.5,markersize=3)

		plt.xlabel(r'$radius \,\,(arcsec)$')
		plt.ylabel(r'$ number\,\, density$')
		plt.title(r'{} < mag_star < {}'.format(i,i+1))
		plt.legend()
		plt.savefig(os.path.join(home,'results/plots/stars_area_{}_photoflag{}_magmin{}_{}.png'.format(filt,photoflag,i,distance_to_border)))
	# plt.show()

if __name__ == '__main__':
	main()
