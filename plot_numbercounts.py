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

def check_star_in_border(cat_stars,photoflag,survey,distance_to_border):
    splus_tiles_file = os.path.join('./data/all_pointings.csv')
    splus_tiles = Table.read(splus_tiles_file, format="csv")
    tilesize = 1.4 # degrees
    htilesize = tilesize / 2.
    
    data = fits.open(cat_stars)
    catalog = data[1].data
    sel = (catalog['PhotoFlag'] == photoflag)
    starcat = catalog[sel]
    print(len(starcat))
    coords_stars = SkyCoord(ra=starcat['RA']*u.deg,dec=starcat['Dec']*u.deg, frame=FK5)
    idx_stars, idx_stars2, d2d, d3d = coords_stars.search_around_sky(coords_stars,distance_to_border*u.arcsec)
    stars_near_stars = [x for x, y in collections.Counter(idx_stars).items() if y > 1]
    idx = np.where(splus_tiles['PID'] == survey)
    tiledata = splus_tiles[idx]
    linspace_npnts = 5
    coords = SkyCoord(ra=tiledata['RA'], dec=tiledata['DEC'], unit=(u.hourangle,u.degree))
    stars_near_borders = np.array([])
    for coord in coords:
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
    	idx_stars1, idx_coords1, d2d, d3d = coords_1.search_around_sky(coords_stars,distance_to_border*u.arcsec)
    	idx_stars2, idx_coords2, d2d, d3d = coords_2.search_around_sky(coords_stars,distance_to_border*u.arcsec)
    	idx_stars3, idx_coords3, d2d, d3d = coords_3.search_around_sky(coords_stars,distance_to_border*u.arcsec)
    	idx_stars4, idx_coords4, d2d, d3d = coords_4.search_around_sky(coords_stars,distance_to_border*u.arcsec)
    	stars_near_borders = np.concatenate([stars_near_borders, idx_stars1])
    	stars_near_borders = np.concatenate([stars_near_borders, idx_stars2])
    	stars_near_borders = np.concatenate([stars_near_borders, idx_stars3])
    	stars_near_borders = np.concatenate([stars_near_borders, idx_stars4])

    stars_near_something = np.concatenate([stars_near_stars,stars_near_borders])
    stars_near_something = [x for x, y in collections.Counter(stars_near_something).items() if y > 1]
    print(len(stars_near_something))

    coords_stars = np.delete(starcat,stars_near_something,axis=0)
    print(len(coords_stars))
    return coords_stars


def sel_cat(rad,cat_stars,cat_splus,filt,photoflag):
	
	coords_stars = SkyCoord(ra=cat_stars['RA']*u.deg,dec=cat_stars['Dec']*u.deg, frame=FK5)

	data = fits.open(cat_splus)
	catalog = data[1].data
	sel = (catalog['PhotoFlag'] == photoflag)
	spluscat = catalog[sel]
	coords_splus = SkyCoord(ra=spluscat['RA']*u.deg,dec=spluscat['Dec']*u.deg, frame=FK5)

	idx_stars, idx_splus, d2d, d3d = coords_splus.search_around_sky(coords_stars,rad*u.arcsec)
	mag_star = []
	mag_obj = []
	for i in idx_stars:
		mag_star.append(cat_stars['Pmag'][i])
	for j in idx_splus:
		mag_obj.append(spluscat['{}'.format(filt)][j])
	table = Table([idx_stars,mag_star, idx_splus,mag_obj, d2d, d3d],names=['index_starcoord','mag_star','index_spluscoord','mag_obj','distance_2d','distance_3d'])
	return table


def plot_radial_bins2(f,magstar_min,magstar_max):

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
	filt = 'r_auto'
	photoflag = 0
	radius_steps = 10
	survey = 'STRIPE82'
	distance_to_border = 100 # arcsec

	magstar_min = np.arange(10.0,15.0,1)

	radius = np.arange(2,114,radius_steps)
	label = ['mag_obj < 17','mag_obj < 18','mag_obj < 19','mag_obj < 20','mag_obj < 21','mag_obj < 22']
	marker= ['o', 's', '^', 'p', 'X','P']

	starcat = check_star_in_border(catstars,photoflag,survey,distance_to_border)
	for r in radius:
		print(r)
		globals()["Table_r" + str(r)] = sel_cat(r,starcat,catfile,filt,photoflag)


	for i in magstar_min:
		print(i)
		plt.figure()
		for j in range(6): 
			number_dens = []
			for r in radius:
<<<<<<< HEAD
				if r <=102:
					f_big = globals()["Table_r" + str(r+radius_steps)]
					f_small = globals()["Table_r" + str(r)]
					df0,df1,df2,df3,df4,df5 = plot_radial_bins2(f_big,i,i+1)
					df0_in,df1_in,df2_in,df3_in,df4_in,df5_in = plot_radial_bins2(f_small,i,i+1)
					area = np.pi*((r+radius_steps)**2 - (r)**2)
					print(area)
					df = pd.concat([locals()["df"+str(j)],locals()["df"+str(j)+"_in"]], keys=['col1', 'col2'])
					dropped = df.drop_duplicates(['col1'], keep=False)
					print(len(df0),len(df0_in),len(dropped))
					counts = np.median(dropped['col2'])
					number_density = [counts/area]
					number_dens.append(number_density)
			counts = np.where(np.isnan(number_dens), 0, number_dens)
			plt.plot(radius[0:11],counts,'o--',label=label[j],marker=marker[j],lw=0.5,markersize=3)
=======
				# plot_radial_bins2(r,catstars,catfile,filt,photoflag,i,i+1)
				df0,df1,df2,df3,df4,df5 = plot_radial_bins2(r,catstars,catfile,filt,photoflag,i,i+1)
				df0_in,df1_in,df2_in,df3_in,df4_in,df5_in = plot_radial_bins2(r-2,catstars,catfile,filt,photoflag,i,i+1)
				print(len(df4))
		
				area = np.pi*(r)**2 - np.pi*(r-2)**2 #calculate the area of the ring
				df = pd.concat([locals()["df"+str(j)],locals()["df"+str(j)+"_in"]], keys=['col1', 'col2']) #compare the objects that are in both the bigger and smaller radius 
				dropped = df.drop_duplicates(['col1'], keep=False) #drop the duplicates, that way we are counting only the objects present in the ring.
				counts = np.median(dropped['col2'])
				number_density = [counts/area]
				number_dens.append(number_density)
				print(r,number_density)
			counts = np.where(np.isnan(number_dens), 0, number_dens) #replaces nan with 0, that way we can visualize in which radius there are no new objects appearing
			plt.plot(radius,counts,'o--',label=label[j],marker=marker[j],lw=0.5,markersize=3)
>>>>>>> 4f78d3be588344f85ba7342fc0cd3bf3375094dc

		plt.xlabel(r'$radius \,\,(arcsec)$')
		plt.ylabel(r'$ number\,\, density$')
		plt.title(r'{} < mag_star < {}'.format(i,i+1))
		plt.legend()
		plt.savefig(os.path.join(home,'results/plots/teststars_area_{}_photoflag{}_magmin{}_{}.png'.format(filt,photoflag,i,distance_to_border)))
	plt.show()

if __name__ == '__main__':
	main()
