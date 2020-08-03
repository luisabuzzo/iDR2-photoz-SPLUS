import healpy as H
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits,ascii
import os
import matplotlib.cm as cm
from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy import units as u

def add_skyplot_tile(coords,tilesize,color,Tile,CenterPoints,CornerPoints,lw,zorder=1):
    cosdec = np.cos(np.radians(coords[1]))
    linspace_npnts = 3
    htilesize = tilesize / 2.
    l1 = coords[0]-htilesize/cosdec
    l2 = coords[0]+htilesize/cosdec
    lon_pnts = [l1,l2,l2,l1]
    b1 = coords[1] - htilesize
    b3 = coords[1] + htilesize
    lat_pnts = [b1,b1,b3,b3]
    line_segments = [[np.linspace(lon_pnts[0],lon_pnts[1],linspace_npnts),np.linspace(lat_pnts[0],lat_pnts[1],linspace_npnts)],\
            [np.linspace(lon_pnts[1],lon_pnts[2],linspace_npnts),np.linspace(lat_pnts[1],lat_pnts[2],linspace_npnts)],\
            [np.linspace(lon_pnts[2],lon_pnts[3],linspace_npnts),np.linspace(lat_pnts[2],lat_pnts[3],linspace_npnts)],\
            [np.linspace(lon_pnts[3],lon_pnts[0],linspace_npnts),np.linspace(lat_pnts[3],lat_pnts[0],linspace_npnts)]]

    if Tile:
        for t in range(0,4):
            H.projplot(line_segments[t][0],line_segments[t][1],lonlat=True,color=color,lw=lw,coord=['C'],zorder=zorder)
    
    if CenterPoints:
        H.projplot(coords[0],coords[1],'2',lonlat=True,color=color,coord=['C'],zorder=zorder)
    
    if CornerPoints:
        for t in range(0,4):
            H.projplot(lon_pnts[t],lat_pnts[t],'+',lonlat=True,color=color,lw=lw,coord=['C'],zorder=zorder)
    
    return 0

def make_skyplot(moll=True, ecliptic=False, footprint=True, idr2=True, cmap='Greys', title='test', highlight=False, highlight_survey='HYDRA', highlight_tile='HYDRA_0011', path='./',bgmap='none'):
    # NOTES: 
    # C = Celestial (Equatorial), G = Galactic, E = Ecliptic
    # if lonlat=True, then coords should be given in (lon,lat) with both in degrees 
 
    alpha = 0.3
    lw = 2

    #dust_map = 'data/lambda_sfd_ebv.fits'
    #bgmap = H.read_map(os.path.join(path, dust_map))

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
            ddec = +3
            ra_positions = [0,30,60,90,120,150,210,240,270,300,330,0,0]
            dec_positions = [0,0,0,0,0,0,0,0,0,0,0,30,60]
            labels = ['0h','2h','4h','6h','8h','10h','14h','16h','18h','20h','22h','30','60']
            for l in range(0,len(labels)):
                H.projtext(ra_positions[l]+dra, dec_positions[l]+ddec, labels[l], coord='C', lonlat=True)
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
                add_skyplot_tile([coords[j].ra.value,coords[j].dec.value],tilesize,col,Tile=True,CenterPoints=False,CornerPoints=False,lw=1)
 
        if highlight:
            idx = (splus_tiles['PID'] == highlight_survey) & (splus_tiles['NAME'] == highlight_tile)
            tiledata = splus_tiles[idx]
            coords = SkyCoord(ra=tiledata['RA'], dec=tiledata['DEC'], unit=(u.hourangle,u.degree))
            lon = coords[0].ra.value
            lat = coords[0].dec.value
            add_skyplot_tile([lon,lat],tilesize,'r',Tile=True,CenterPoints=False,CornerPoints=False,lw=2,zorder=100)

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
                add_skyplot_tile([coords.ra.value,coords.dec.value],tilesize,'g',Tile=True,CenterPoints=False,CornerPoints=False,lw=1)

    H.graticule()   # draw the coordinate system axes

    return 0

def plot_objects(ra,dec,markersize,color,symbol):
    coords = SkyCoord(ra=ra, dec=dec, unit=(u.degree,u.degree))
    for k in range(0,len(coords)):
        H.projplot(coords[k].ra.value,coords[k].dec.value,symbol,lonlat=True,color=color,coord=['C'],markersize=markersize)

def skyplot(footprint=True, idr2=True, highlight=False, highlight_survey='HYDRA', highlight_tile='HYDRA-0011', plot_catalog=0, path='./', catalog_path='Catalogs/', plotonscreen=0, bgmap='none'):

    make_skyplot(moll=True, ecliptic=True, footprint=footprint, idr2=idr2, title='S-PLUS', highlight=highlight, highlight_survey=highlight_survey, highlight_tile=highlight_tile, path=path, bgmap=bgmap)

    # overplotting catalogs, make sure that coordinates are in degrees
    # plot_catalog = 0    # 0 = do nothing
    # plot_catalog = 1    # overplot objects from S-PLUS catalogs
    # plot_catalog = 2    # overplot other things 
    # plot_catalog = 3    # overplot Julia Gschwend spec-z catalog 


    if plot_catalog == 0:
        print('No objects overplotted.')

    if plot_catalog == 1:
        print('Overplotting S-PLUS catalogs...')

        # overplot S-PLUS sources
        catalog_path = path + catalog_path
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
        print('Overplotting some catalog...')
        ra = [0,100,150]
        dec = [-30,0,30] 
        plot_objects(ra,dec,symbol='o',color='r',markersize=10) # plot them on the map

    if plot_catalog == 3:
        cat = read_specz_catalog(path + catalog_path,'SPECZ_SAMPLE_JULIA_UNIQUE_GALAXIES.fits')
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

    return 0    
