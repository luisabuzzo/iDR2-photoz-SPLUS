import numpy as np
from astropy.io import fits,ascii
from astropy.table import Table,vstack,hstack

def get_tile_names(catalog_path, surveys, infiles):
# make csv files containing the names of each of the pointings contained in the iDR2 release catalogs

	for i in range(0,len(surveys)):
		base_infile = catalog_path + infiles[i]
		base_outfile = catalog_path + surveys[i]+'.csv'
		# open table as astropy table
		cat = Table.read(base_infile, format='fits')
		fields = np.array(cat['FIELD'])
		ufields = np.unique(fields)
		t = Table()
		t['FIELD'] = ufields
		print('Read catalog: ', surveys[i])
		print('Tiles: ', t)
		ascii.write(t, base_outfile, format='csv', overwrite=True)  

	return 0

def clean_IDR2files(infile, columnset, outfile):
# this routine starts from the original IDR2 catalogs and does some cleaning by  
# renaming columns and reducing the final output size by changing unnecessary long 
# data types and removing a number of unused columns
# main changes:
#	- uJAVA columns renamed to u (to be consistent with the other broad bands)
#	- dtype of FIELD changed to 15 character string
#	- dtype of ID changed from 34 character string (!) to an int32 by removing 
# 		the survey name ('SPLUS'), removing the field name (already in FIELD), 
#		removing the 'griz' label, and changing the object number into an integer
#		(each object can still be identified through the combination FIELD + ID)
#	- dtype of DEC changed to float64, just like RA
#	- all other floats changed to float32
#	- all petrosian magnitude columns were removed
#	- all aper magnitudes removed for narrowbands
#	- columns FWHM_n, MUMAX, A, B, THETA, FlRadDet, KrRadDet removed

	# open table as astropy table
	cat = Table.read(infile, format='fits')

	# change the IDs from a very long string to the simple object number
	id = cat['ID']
	new_ids = np.zeros((len(id)))
	for i in range(0,len(id)):
		idi = id[i] 
		new_ids[i] = np.int_(idi.split('.')[2])
	
	# show the columns
	colnames = cat.colnames
	print("This table contains the following columns:")
	print(cat.info)

	# copy the columns with new desired data types
	cat['FIELD2'] = np.array(cat['FIELD'],dtype='U15')
	cat['ID2'] = np.array(new_ids,dtype='int32')
	cat['Dec2'] = np.array(cat['Dec'],dtype='float64')
	cat['X2'] = np.array(cat['X'],dtype='float32')
	cat['Y2'] = np.array(cat['Y'],dtype='float32')
	cat['s2nDet2'] = np.array(cat['s2nDet'],dtype='float32')
	cat['FlRadDet2'] = np.array(cat['FlRadDet'],dtype='float32')
	cat['s2n_z_petro2'] = np.array(cat['s2n_z_petro'],dtype='float32')

	# remove the columns that have been copied
	cat.remove_columns(['FIELD','ID','Dec','X','Y','s2nDet','FlRadDet','s2n_z_petro'])

	# rename the new columns to the old columns that were removed
	cat.rename_column('FIELD2','FIELD')
	cat.rename_column('ID2','ID')
	cat.rename_column('Dec2','Dec')
	cat.rename_column('X2','X')
	cat.rename_column('Y2','Y')
	cat.rename_column('s2nDet2','s2nDet')
	cat.rename_column('FlRadDet2','FlRadDet')
	cat.rename_column('s2n_z_petro2','s2n_z_petro')

	# rename annoying uJAVA columns to u
	cat.rename_column('uJAVA_auto','u_auto')
	cat.rename_column('euJAVA_auto','eu_auto')
	cat.rename_column('s2n_uJAVA_auto','s2n_u_auto')
	cat.rename_column('uJAVA_petro','u_petro')
	cat.rename_column('euJAVA_petro','eu_petro')
	cat.rename_column('s2n_uJAVA_petro','s2n_u_petro')
	cat.rename_column('uJAVA_aper','u_aper')
	cat.rename_column('euJAVA_aper','eu_aper')
	cat.rename_column('s2n_uJAVA_aper','s2n_u_aper')

	if columnset == 'all':
		selected_columns = ['FIELD', 'ID', 'RA', 'Dec', 'X', 'Y', 'ISOarea', 's2nDet', 'PhotoFlag', \
							'FWHM', 'u_auto', 'eu_auto', 's2n_u_auto', 'F378_auto', 'eF378_auto', 's2n_F378_auto', \
		 					'F395_auto', 'eF395_auto', 's2n_F395_auto', 'F410_auto', 'eF410_auto', 's2n_F410_auto', \
		 					'F430_auto', 'eF430_auto', 's2n_F430_auto', 'g_auto', 'eg_auto', 's2n_g_auto', \
		    				'F515_auto', 'eF515_auto', 's2n_F515_auto',  'r_auto', 'er_auto', 's2n_r_auto', \
						    'F660_auto', 'eF660_auto', 's2n_F660_auto', 'i_auto', 'ei_auto', 's2n_i_auto', \
						    'F861_auto', 'eF861_auto', 's2n_F861_auto', 'z_auto', 'ez_auto', 's2n_z_auto', \
						    'u_aper', 'eu_aper', 's2n_u_aper', 'g_aper', 'eg_aper', 's2n_g_aper', \
						    'r_aper', 'er_aper', 's2n_r_aper', 'i_aper', 'ei_aper', 's2n_i_aper', \
						    'z_aper', 'ez_aper', 's2n_z_aper']


	### for a wider selection of columns, you can define a new column set here
	#if columnset == 'all0':
	#	selected_columns = ['FIELD', 'ID', 'RA', 'Dec', 'X', 'Y', 'ISOarea', 's2nDet', 'PhotoFlag', 'FWHM', 'MUMAX', \
	#	'A', 'B', 'THETA', 'FlRadDet', 'KrRadDet', 'u_auto', 'eu_auto', 's2n_u_auto', 'u_aper', 'eu_aper', 's2n_u_aper', \
	#	'F378_auto', 'eF378_auto', 's2n_F378_auto', 'F378_aper', 'eF378_aper', 's2n_F378_aper', 'F395_auto', 'eF395_auto',\
	#	 's2n_F395_auto','F395_aper', 'eF395_aper', 's2n_F395_aper', 'F410_auto', 'eF410_auto', 's2n_F410_auto', 'F410_aper',\
	#	  'eF410_aper', 's2n_F410_aper', 'F430_auto', 'eF430_auto', 's2n_F430_auto', 'F430_aper', 'eF430_aper', 's2n_F430_aper',\
	#	   'g_auto', 'eg_auto', 's2n_g_auto', 'g_aper', 'eg_aper', 's2n_g_aper', 'F515_auto', 'eF515_auto', 's2n_F515_auto',\
	#	    'F515_aper', 'eF515_aper', 's2n_F515_aper', 'r_auto', 'er_auto', 's2n_r_auto', 'r_aper', 'er_aper', 's2n_r_aper',\
	#	     'F660_auto', 'eF660_auto', 's2n_F660_auto', 'F660_aper', 'eF660_aper', 's2n_F660_aper', 'i_auto', 'ei_auto',\
	#	      's2n_i_auto', 'i_aper', 'ei_aper', 's2n_i_aper', 'F861_auto', 'eF861_auto', 's2n_F861_auto', 'F861_aper', \
	#	      'eF861_aper', 's2n_F861_aper', 'z_auto', 'ez_auto', 's2n_z_auto', 'z_aper', 'ez_aper', 's2n_z_aper']

	cat = cat[selected_columns] 

	# show updated column structure
	print("Updated table:")
	print(cat.info)

	# write to file
	print("Writing new table...")
	cat.write(outfile, format='fits', overwrite=True)

	return 0

def remove_columns(infile, columnset, outfile):
# this routine allows one to further down-select the number of columns in a given catalog
# when columnset = ID, only a very limited set of columns is maintained to allow some basic checks across the survey
# when columnset = ugriz, the columns are limited to auto and aper magnitudes for each broad band, plus sufficient additional columns related to IDs, positions and quality

	# open table as astropy table
	cat = Table.read(infile, format='fits')

	# show the columns
	colnames = cat.colnames
	print("This table contains the following columns:")
	print(cat.info)

	if columnset == 'id':
		selected_columns = ['FIELD', 'ID', 'RA', 'Dec', 'X', 'Y', 'ISOarea', 's2nDet', 'PhotoFlag', 'FWHM', 'r_auto', 'er_auto']

	if columnset == 'ugriz':
		selected_columns = ['FIELD', 'ID', 'RA', 'Dec', 'X', 'Y', 'ISOarea', 's2nDet', 'PhotoFlag', 'FWHM', 'u_auto', 'eu_auto', 's2n_u_auto', 'u_aper', 'eu_aper', 's2n_u_aper', 'g_auto', 'eg_auto', 's2n_g_auto', 'g_aper', 'eg_aper', 's2n_g_aper', 'r_auto', 'er_auto', 's2n_r_auto', 'r_aper', 'er_aper', 's2n_r_aper', 'i_auto', 'ei_auto', 's2n_i_auto', 'i_aper', 'ei_aper', 's2n_i_aper', 'z_auto', 'ez_auto', 's2n_z_auto', 'z_aper', 'ez_aper', 's2n_z_aper']

	cat = cat[selected_columns] 

	# show updated column structure
	print("Updated table:")
	print(cat.info)

	# write to file
	print("Writing new table...")
	cat.write(outfile, format='fits', overwrite=True)
	#hdu = bintablehdu(data=cat)
	#hdu.writeto(outfile)

	return 0

def make_smaller_catalogs(catalog_path, surveys, columnset, origfiles):
# cleaning and reducing catalog file size by limiting the number of columns
# calls clean_IDR2files() when the input files are the iDR2 release catalogs
# calls remove_columns() when the input files are already processed iDR2 catalogs

	for i in range(0,len(surveys)):

		base_name = catalog_path + surveys[i]
		if columnset == 'all' :
			clean_IDR2files(catalog_path + origfiles[i], columnset, base_name + '_' + columnset + '.fits')

		if columnset == 'id' or columnset == 'ugriz' :
			remove_columns(base_name + '_all.fits', columnset, base_name + '_' +columnset + '.fits')

	return 0

def merge_catalogs(catalog_path, surveys, outfile):
# routine to merge individual catalogs of the three sub-surveys into a single merged catalog
# the output is convenient for searching for objects in S-PLUS or plotting objects across 
# the sky and doing other statistics

	for i in range(0,len(surveys)):
		cat = Table.read(catalog_path + surveys[i]+'_id.fits', format='fits')
		cat['SURVEY'] = surveys[i]
		if i == 0:
			mcat = cat
			del cat
		if i > 0:
			mcat = vstack([mcat,cat])
			del cat

	# write to file
	print("Writing merged table...")
	mcat.write(catalog_path + outfile, format='fits', overwrite=True)

	return 0

def main():

	path = './'							# working directory
	catalog_path = path + 'Catalogs/'	# location of all input and output catalogs
	surveys = ['STRIPE82','HYDRA','MS']
	origfiles = 'SPLUS_' + surveys + '_master_catalogue_iDR2_december_2019.fits'	# original iDR2 release files
	origfiles[0] = 'SPLUS_STRIPE82_SDSS_master_catalogue_iDR2_december_2019.fits'	# fix because the STRIPE82 iDR2 filename has 'SDSS' in it

	### uncomment to make csv files containing the names of each of the pointings contained in the iDR2 release catalogs
	#get_tile_names(catalog_path, surveys, origfiles)

	### uncomment for initial cleaning, reducing and repairing of the input catalogs
	#make_smaller_catalogs(catalog_path, surveys, 'all', origfiles)

	### uncomment to make a broad-band catalog
	#make_smaller_catalogs(catalog_path, surveys, 'ugriz', origfiles)

	### uncomment to make an ID catalog
	#make_smaller_catalogs(catalog_path, surveys, 'id', origfiles)

	### uncomment to merge the three id catalogs
	#merge_catalogs(catalog_path, surveys, 'merged_id.fits')

	return

if __name__ == '__main__':
      main()

### CODE BY RODERIK OVERZIER, S-PLUS PROJECT###

