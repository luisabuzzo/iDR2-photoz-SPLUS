import math
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits,ascii

def make_plots(cat,starcat,filter,plotonscreen):

	fig, axs = plt.subplots(ncols=3, nrows=1, figsize=(12,6)) #, sharex='col', sharey='row', gridspec_kw={'hspace': 0, 'wspace': 0})

	axs[0].axis([13,25,-1.5,1.5])
	#axs[0].plot(cat['r_petro'],cat['r_auto']-cat['r_aper'],'k.',markersize=1,zorder=1)
	axs[0].scatter(starcat[filter+'_petro'],starcat[filter+'_auto']-starcat[filter+'_aper'],c=starcat['fwhm'],s=1,zorder=10,cmap='jet')
	axs[0].set_xlabel(filter+' (petro)')
	axs[0].set_ylabel(filter+' (auto) - '+filter+' (aper)')
	axs[0].plot([13,25],[0,0],'r:',zorder=50)

	axs[1].axis([13,25,-1.5,1.5])
	#axs[1].plot(cat['r_petro'],cat['r_petro']-cat['r_auto'],'k.',markersize=1,zorder=1)
	axs[1].scatter(starcat[filter+'_petro'],starcat[filter+'_petro']-starcat[filter+'_auto'],c=starcat['fwhm'],s=1,zorder=10,cmap='jet')
	axs[1].set_xlabel(filter+' (petro)')
	axs[1].set_ylabel(filter+' (petro) - '+filter+' (auto)')
	axs[1].plot([13,25],[0,0],'r:',zorder=50)

	axs[2].axis([13,25,-1.5,1.5])
	#axs[2].plot(cat['r_petro'],cat['r_petro']-cat['r_aper'],'k.',markersize=1,zorder=1)
	axs[2].scatter(starcat[filter+'_petro'],starcat[filter+'_petro']-starcat[filter+'_aper'],c=starcat['fwhm'], s=1,zorder=10,cmap='jet')
	axs[2].set_xlabel(filter+' (petro)')
	axs[2].set_ylabel(filter+' (petro) - '+filter+' (aper)')
	axs[2].plot([13,25],[0,0],'r:',zorder=50)

	plt.tight_layout()
	plt.savefig('results/STRIPE82_magtest_'+filter+'.png',dpi=150)

	if plotonscreen == 1:
		plt.show()


def mag_tests(cat,plotonscreen):

	sel = (cat['PhotoFlag'] == 0)
	cat = cat[sel]

	sel = (cat['fwhm'] < 2.0) & (cat['fwhm'] > 0.5)
	starcat = cat[sel]

	filters = ['u','g','r','i','z']
	for i in range(0,len(filters)):
		make_plots(cat,starcat,filters[i],plotonscreen)

	return 0

def main():

	path = './'
	catalog_path = path + 'Catalogs/'
	catfile = catalog_path + 'STRIPE82_magtest_broad.fits'
	data = fits.open(catfile)
	cat = data[1].data
	data.close
	cols = cat.columns
	colnames = np.array(cols.names)
	plotonscreen = 0
	mag_tests(cat,plotonscreen)

	return 0

if __name__ == '__main__':
	main()
