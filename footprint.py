#!/usr/bin/env python
import os
import time
from os import system
import numpy as np
import matplotlib.pyplot as plt
import healpy as H
import matplotlib.cm as cm
import matplotlib.patches as mpatches
from astropy.table import Table
from astropy.coordinates import SkyCoord,search_around_sky,match_coordinates_sky
import astropy.units as u
from astropy.coordinates.builtin_frames import GeocentricTrueEcliptic
from astropy.io import fits,ascii
from datetime import date
from random import randint

def plotEqTile(coords,tilesize,color,Tile,CenterPoints,CornerPoints,lw):
    cosdec = np.cos(np.radians(coords[1]))
    linspace_npnts = 3
    htilesize = tilesize / 2.

    l1 = coords[0] - htilesize / cosdec
    l2 = coords[0] + htilesize / cosdec
    l3 = l2
    l4 = l1
    b1 = coords[1] - htilesize
    b2 = b1
    b3 = coords[1] + htilesize
    b4 = b3
    seg1_l = np.linspace(l1,l2,linspace_npnts) 
    seg1_b = np.linspace(b1,b2,linspace_npnts) 
    seg2_l = np.linspace(l2,l3,linspace_npnts) 
    seg2_b = np.linspace(b2,b3,linspace_npnts) 
    seg3_l = np.linspace(l3,l4,linspace_npnts) 
    seg3_b = np.linspace(b3,b4,linspace_npnts) 
    seg4_l = np.linspace(l4,l1,linspace_npnts) 
    seg4_b = np.linspace(b4,b1,linspace_npnts) 

    if Tile:
        H.projplot(seg1_l,seg1_b,lonlat=True,color=color,lw=lw,coord=['C'])
        H.projplot(seg2_l,seg2_b,lonlat=True,color=color,lw=lw,coord=['C'])
        H.projplot(seg3_l,seg3_b,lonlat=True,color=color,lw=lw,coord=['C'])
        H.projplot(seg4_l,seg4_b,lonlat=True,color=color,lw=lw,coord=['C'])
        # NOTES: 
        # C = Celestial (Equatorial), G = Galactic, E = Ecliptic
        # if lonlat=True, then coords should be given in (lon,lat) with both in degrees 
    if CenterPoints:
        H.projplot(coords[0],coords[1],'2',lonlat=True,color=color,coord=['C'])
    if CornerPoints:
        H.projplot(l1,b1,'+',lonlat=True,color=color,lw=lw,coord=['C'])
        H.projplot(l2,b2,'+',lonlat=True,color=color,lw=lw,coord=['C'])
        H.projplot(l3,b3,'+',lonlat=True,color=color,lw=lw,coord=['C'])
        H.projplot(l4,b4,'+',lonlat=True,color=color,lw=lw,coord=['C'])
    return 0

def make_figure(moll=True, ecliptic=False, footprint=True, idr2=True, cmap='Greys', title='test', highlight=False, highlight_survey='HYDRA', highlight_tile='HYDRA-0011', path='./'):

    alpha = 0.3
    lw = 2

    dust_map = 'data/lambda_sfd_ebv.fits'
    bgmap = H.read_map(os.path.join(path, dust_map))

    if highlight:
        title = highlight_survey + '/' + highlight_tile

    # plot the map projection
    cmap = cm.get_cmap(cmap)
    vmin = np.log10(np.percentile(bgmap, 20))
    vmax = 1+np.log10(np.percentile(bgmap, 99))
    cmap.set_under('w') # Set the color for low out-of-range values when norm.clip = False to white
    if moll:
        cart = H.mollview(np.log10(bgmap), coord=['G', 'C'], cmap=cmap, cbar=False, notext=False, title=title, min=vmin, max=vmax)

        add_numbers_on_axes = 1
        if add_numbers_on_axes == 1:
            dra = -2
            ddec = +1
            H.projtext(  0+dra, 2+ddec,'0h', coord='C', lonlat=True)
            H.projtext( 30+dra, 2+ddec,'2h', coord='C', lonlat=True)
            H.projtext( 60+dra, 2+ddec,'4h', coord='C', lonlat=True)
            H.projtext( 90+dra, 2+ddec,'6h', coord='C', lonlat=True)
            H.projtext(120+dra, 2+ddec,'8h', coord='C', lonlat=True)
            H.projtext(150+dra, 2+ddec,'10h', coord='C', lonlat=True)
            H.projtext(210+dra, 2+ddec,'14h', coord='C', lonlat=True)
            H.projtext(240+dra, 2+ddec,'16h', coord='C', lonlat=True)
            H.projtext(270+dra, 2+ddec,'18h', coord='C', lonlat=True)
            H.projtext(300+dra, 2+ddec,'20h', coord='C', lonlat=True)
            H.projtext(330+dra, 2+ddec,'22h', coord='C', lonlat=True)
            H.projtext(  0+dra,30+ddec,'30', coord='C', lonlat=True)
            H.projtext(  0+dra,60+ddec,'60', coord='C', lonlat=True)
    else:
        cart = H.cartview(np.log10(bgmap), coord=['G', 'C'], cmap=cmap, cbar=False, notext=True, title=title, min=vmin, max=vmax)

    # draw the ecliptic
    if ecliptic:
        ra = np.linspace(-np.pi,np.pi,100)
        dec = np.zeros(100)+(90+23.5*np.sin(ra))*np.pi/180.
        H.projplot(dec,ra,'--', color='k', lw=2.)

    # draw the S-PLUS footprint
    if footprint:
        # read SPLUS pointings file
        splus_tiles_file = os.path.join(path,'data/all_pointings.csv')
        splus_tiles = Table.read(splus_tiles_file, format="csv")
        splus_tiles['obs'] = 0 # add an extra column to mark if the pointing has been observed or not
        tilesize = 1.4 # degrees
        tilecol = ['b','g','r','m']
        color_the_surveys = False
        surveys = ['SPLUS', 'STRIPE82', 'HYDRA', 'MC']
        for s in range(0,len(surveys)):
            if color_the_surveys:
                col = tilecol[s]
            else:
                col = '0.5'
            idx = np.where(splus_tiles['PID'] == surveys[s])
            tiledata = splus_tiles[idx]
            ntiles = len(tiledata)
            coords = SkyCoord(ra=tiledata['RA'], dec=tiledata['DEC'], unit=(u.hourangle,u.degree))
            for j in range(0,ntiles):
                plotEqTile([coords[j].ra.value,coords[j].dec.value],tilesize,col,Tile=True,CenterPoints=False,CornerPoints=False,lw=1)
 
        if highlight:
            idx = (splus_tiles['PID'] == highlight_survey) & (splus_tiles['NAME'] == highlight_tile)
            tiledata = splus_tiles[idx]
            coords = SkyCoord(ra=tiledata['RA'], dec=tiledata['DEC'], unit=(u.hourangle,u.degree))
            lon = coords[0].ra.value
            lat = coords[0].dec.value
            plotEqTile([lon,lat],tilesize,'r',Tile=True,CenterPoints=True,CornerPoints=False,lw=1)

    # draw the observed idr2 footprint
    if idr2:
        # read SPLUS pointings file
        splus_tiles_file = os.path.join(path,'data/all_pointings.csv')
        splus_tiles = Table.read(splus_tiles_file, format="csv")
        splus_tiles_names = np.array(splus_tiles['NAME'])
        # read IDR2 pointings file
        idr2_tiles_file = os.path.join(path,'data/idr2_pointings.csv')
        idr2_tiles = Table.read(idr2_tiles_file, format="csv")
        idr2_tiles_names = np.array(idr2_tiles['FIELD'])
        tilesize = 1.4 # degrees
        for tile in range(0,len(splus_tiles)):
            sel = (idr2_tiles_names == splus_tiles_names[tile])
            if np.sum(sel) > 0:
                tiledata = splus_tiles[tile]
                coords = SkyCoord(ra=tiledata['RA'], dec=tiledata['DEC'], unit=(u.hourangle,u.degree))
                plotEqTile([coords.ra.value,coords.dec.value],tilesize,'g',Tile=True,CenterPoints=False,CornerPoints=False,lw=1)

    H.graticule()   # draw the coordinate system axes

    return 0

def plot_objects(ra,dec,markersize,color,symbol):

    coords = SkyCoord(ra=ra, dec=dec, unit=(u.degree,u.degree))
    for k in range(0,len(coords)):
        H.projplot(coords[k].ra.value,coords[k].dec.value,symbol,lonlat=True,color=color,coord=['C'],markersize=markersize)

def main():

    path = './'
    catalog_path = 'Catalogs/'
    moll = True
    ecliptic = True
    footprint = True
    idr2 = True
    title = 'S-PLUS'
    highlight = False
    highlight_survey = 'HYDRA'
    highlight_tile = 'HYDRA-0011'

    # overplotting catalogs, make sure that coordinates are in degrees
    plot_catalog = 0    # 0 = do nothing
    # plot_catalog = 1    # overplot objects from S-PLUS catalogs
    # plot_catalog = 2    # overplot other things 
    # plot_catalog = 3    # overplot Julia spectroscopic catalog

    # make the sky projection with footprints
    make_figure(moll=moll, ecliptic=ecliptic, footprint=footprint, idr2=idr2, title=title, highlight=highlight, highlight_survey=highlight_survey, highlight_tile=highlight_tile, path=path)

    # overplot things
    if plot_catalog == 0:
        print('No objects overplotted.')

    if plot_catalog == 1:
        print('Overplotting S-PLUS catalogs...')

        # overplot S-PLUS sources
        catfile = 'merged_id.fits'
        data = fits.open(path + catalog_path + catfile)
        cat = data[1].data
        data.close
        # select 10000 random sources
        nrand = 10000
        sel = np.int_(np.zeros((nrand)))
        for p in range(0,nrand):
            num = randint(0, len(cat))
            sel[p] = num
        cat = cat[sel]
        plot_objects(cat['RA'],cat['Dec'],symbol='.',color='b',markersize=1) # plot them on the map

    if plot_catalog == 2:
        print('Overplotting some other catalog...')
        ra = [0,100,150]
        dec = [-30,0,30] 
        plot_objects(ra,dec,symbol='o',color='r',markersize=10) # plot them on the map

	if plot_catalog == 3:
        filename = 'SPECZ_SAMPLE_JULIA_UNIQUE_GALAXIES.fits'
		cat = fits.getdata(path + catalog_path + filename)
		# select 100000 random sources
		nrand = 100000
		sel = np.int_(np.zeros((nrand)))
		for p in range(0,nrand):
		    num = randint(0, len(cat))
		    sel[p] = num
		cat = cat[sel]
		plot_objects(cat['RA'],cat['DEC'],symbol='.',color='r',markersize=1) # plot them on the map		
	
	if highlight == True:
		plt.savefig(path + 'results/' + highlight_survey + '-' + highlight_tile + '_skyplot.png',transparent=True,dpi=600)
	else:	
		plt.savefig(path + 'skyplot.png',transparent=True,dpi=600)

	if plotonscreen == 1:	
		plt.show()
	else:
		plt.close()
        
    return
              
if __name__ == '__main__':
    main()
