from regions import read_ds9, write_ds9, CircleSkyRegion
from astropy.table import Table
from astropy.io import fits,ascii
import numpy as np
import getpass
from astropy.coordinates import SkyCoord, Angle
import os

def make_region(ra,dec,radius_arcsec,flag):
	strflag = str(flag)
	center = SkyCoord(ra, dec, unit='deg')
	radius = Angle(radius_arcsec, 'deg')
	meta = {'color': 'green','text': strflag}
	#meta = {'color': 'blue'}
	region = CircleSkyRegion(center, radius, meta)
	return region

# define path to the main directory
if getpass.getuser() == "luisabuzzo":
    home = "/home/luisabuzzo/Work/Master/SPLUS/"
elif getpass.getuser() == "mlgbuzzo":
    home = "/home/mlgbuzzo/LePhare/"
elif getpass.getuser() == "roderik":
    home = "/Users/roderik/Dropbox/splus/IDR2/"

catfile = os.path.join(home,'Catalogs/STRIPE82_all.fits')
data = fits.open(catfile)
cat = data[1].data

#sel = np.logical_and((cat['FIELD'] == 'STRIPE82-0160'),(cat['r_auto']/cat['er_auto'] < 5))
sel = np.logical_and((cat['FIELD'] == 'STRIPE82-0160'),(cat['r_auto'] < 21))

#sel = (cat['FIELD'] == 'STRIPE82-0160')
cat = cat[sel]

radius_arcsec = 2.0
regions = []
for i in range(0,len(cat)):
	regions.append(make_region(cat['RA'][i],cat['Dec'][i],radius_arcsec/3600.,np.int_(cat['PhotoFlag'][i]))) 
#filename = 'tile0160_sn_lt5.reg'
filename = 'tile0160_r_lt21.reg'

write_ds9(regions, filename)

