import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord,search_around_sky,match_coordinates_sky
from astropy.table import Table,vstack,hstack

# create a SkyCoord object
def create_sky_object(ra,dec):
	coords = SkyCoord(ra=ra, dec=dec, unit=(u.degree,u.degree))
	return coords

# checks if a given location is covered by the S-PLUS and IDR2 footprints within a search radius of tilesep_deg
# checks if an object exists in the iDR2 source catalog within a search radius of catsep_arcsec
# INPUT: RA (degrees), DEC (degrees), TILESEP_DEG (degrees), CATSEP_ARCSEC (ARCSEC)
def check_splus(ra,dec,tilesep_deg,catsep_arcsec):

	print()
	obj = create_sky_object(ra,dec)
	print("Input location (RA,DEC): ", obj.ra.value[0], 'deg', obj.dec.value[0], 'deg')

	# read SPLUS pointings file
	splus_tiles_file = 'data/all_pointings.csv'
	tiles = Table.read(splus_tiles_file, format="csv")
	tiles_coords = SkyCoord(ra=tiles['RA'], dec=tiles['DEC'], unit=(u.hourangle,u.degree))

	# read IDR2 pointings file
	idr2_tiles_file = 'data/idr2_pointings.csv'
	idr2_tiles = Table.read(idr2_tiles_file, format="csv")
	
	sky_radius_deg = tilesep_deg * u.degree

	print("Searching all S-PLUS pointings...")
	idx_data1, idx_data2, d2d, d3d = search_around_sky(obj,tiles_coords,sky_radius_deg) 

	if len(d2d) > 0:
		for i in range(0,len(d2d)):
			location_in_idr2 = '(NOT IN DR2)'
			sel = (idr2_tiles['FIELD'] == tiles['NAME'][idx_data2[i]])
			nfound = np.sum(sel)
			if nfound > 0:
				location_in_idr2 = '(IDR2)'
			print("Your location is within{0:4.1f} degrees of the center of pointings:".format(tilesep_deg))
			print("{0:s} {1:s}/{2:s} (separation:{3:5.2f})".format(location_in_idr2,tiles['PID'][idx_data2[i]],tiles['NAME'][idx_data2[i]],d2d[i].value))
	else:
		print("Nothing found.")

	print()
	print("Searching all iDR2 objects...")
	infile = 'Catalogs/merged_id.fits'
	poscat = Table.read(infile, format='fits')
	obscoords = SkyCoord(ra=poscat['RA'], dec=poscat['Dec'], unit=(u.degree,u.degree))
	sky_radius_deg = (catsep_arcsec/3600.) * u.degree
	idx_data1, idx_data2, d2d, d3d = search_around_sky(obj,obscoords,sky_radius_deg) 
	sor = np.argsort(d2d)
	idx_data1 = idx_data1[sor] 
	idx_data2 = idx_data2[sor] 
	d2d = d2d[sor] 
	d3d = d3d[sor] 
	if len(d2d) > 0:
		for i in range(0,len(d2d)):
			print("Object {0:s} in area/tile {1:7s}/{2:12s} with separation: {3:5.2f} arcsec)".format(str(poscat['ID'][idx_data2[i]]),poscat['SURVEY'][idx_data2[i]],poscat['FIELD'][idx_data2[i]],d2d[i].value*3600.))
	else:
		print("Nothing found.")

	return 
