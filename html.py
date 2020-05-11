#import math
#import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits,ascii
#from astropy.stats import sigma_clip
#from scipy.optimize import curve_fit
#from scipy.stats import gaussian_kde
from astropy import units as u
from astropy.coordinates import SkyCoord,search_around_sky,match_coordinates_sky
#import healpy as H
import os
import time
from os import system
#import matplotlib.cm as cm
#import matplotlib.patches as mpatches
from astropy.table import Table,vstack,hstack
from astropy.coordinates.builtin_frames import GeocentricTrueEcliptic
from datetime import date
#from random import randint

def make_html_tiles(area,tiles,splus_tiles):

    path = './html/' + area + '/'

    html_open = ['<!DOCTYPE html>','<html>','<head>','<style>img {display: block; margin-left: auto; margin-right: auto;}</style>','</head>','<body>']
    
    html_close = ['</body>','</html>']
    html_break = '<br>'

    for i in range(0,len(tiles)):

        html_code = []
        for k in range(0,len(html_open)):
            html_code.append(html_open[k]+'\n')

        tile = tiles[i]
        sel = (tile == splus_tiles['NAME'])
        seltile = splus_tiles[sel]
        racenter = str(np.array(seltile['RA']))
        deccenter = str(np.array(seltile['DEC']))

        idd = (tile['FIELD']).split('-')[1]
        tile_id = area + '-' + area + '-' + idd
        html_file = path + tile_id + '.html'

        html_title = '<h1>' + tile_id + '</h1>'
        html_info = ['Center position: RA = ', racenter, ' DEC = ', deccenter]
                    
        html_objects = ['Notable objects: ',\
                        'SDSS quasars: NN (cat)',\
                        'MW globular clusters: NN (cat)',\
                        'NGC galaxies: NN (cat)',\
                        'Planetary Nebulae: NN (cat)']

        html_code.append(html_title)
        for j in range(0,len(html_info)):
            html_code.append(html_info[j]+'\n')
        html_code.append(html_break)
        html_code.append(html_break)

        for j in range(0,len(html_objects)):
            html_code.append(html_objects[j]+'<br>\n')

        figrootname = tile_id + '_'
        figs = ['skyplot','trilogy0','skyradec','skyxy','skydens','numbercounts12','numbercounts12_depths','concentration','fwhm','snr','color']
        for j in range(0,len(figs)):
            if figs[j] == 'trilogy0':
                figname = figrootname+figs[j]+'.png'
            else:
                figname = figrootname+figs[j]+'.png'
            html_fig = '<img src="' + figname + '"' + ' width="800" onclick=' + '"' + "window.open('"+ figname + "'" + ", '_blank');" + '">'
            html_code.append(html_fig+'\n')
            #if j == 1 or j == 3:
            html_code.append('<br>'+'\n')

        for k in range(0,len(html_close)):
            html_code.append(html_close[k]+'\n')

        np.savetxt(html_file,[html_code],fmt="%s")
        print(html_file, ' saved')

    return

def make_html_global():

    areas = ['STRIPE82','SPLUS','HYDRA']
    splus_tiles_file = './data/all_pointings.csv'
    splus_tiles = Table.read(splus_tiles_file, format="csv")

# make the front page
    html_file = './html/splus.html'

    html_open = ['<!DOCTYPE html>','<head>','<style>body {color: green;}</style>',\
                '<style>a:active {color:#00FF00;}</style>','<style>a:visited {color: #FF0000;}</style>',\
                '<style>img {display: block; margin-left: auto; margin-right: auto;}</style>','</head>','<body>']
   
    html_title = '<p style="color:blue">S-PLUS internal Data Release 2</p>'
    html_skyplot = '<img style="width:50%" src="./skyplot.png" width="600" onclick="window.open(' + "'./skyplot.png', '_blank');" + '>'

    html_close = ['</body>','</html>']
    html_break = '<br>'

    html_code = []
    for k in range(0,len(html_open)): 
        html_code.append(html_open[k]+'\n')

    html_code.append(html_title+'\n')
    html_code.append(html_skyplot+'\n')

    for area in areas:

        html_code.append('<p style="color:blue">' + area + ' tiles</p>')
        html_code.append(html_break)

        tilefile = 'data/' + area + '.csv'
        tiles = Table.read(tilefile, format="csv")
        htmlstring = ''
        
        for i in range(0,len(tiles)):
            tile = tiles[i]
            idd = (tile['FIELD']).split('-')[1]
            ahref_open = '<a href="./' + area + '/' + area + '-' + area + '-' + idd + '.html">'
            ahref_close = '</a>'
            sep = ' | '
            if i < len(tiles)-1:
                htmlstring = htmlstring + ahref_open + idd + ahref_close + sep
            else:
                htmlstring = htmlstring + ahref_open + idd + ahref_close
        html_code.append(htmlstring)

    for k in range(0,len(html_close)):
        html_code.append(html_close[k]+'\n')

    np.savetxt(html_file,[html_code],fmt="%s")
    print("Main webpage created: ", html_file)

# make the individual pages

    for area in areas:
        tilefile = 'data/' + area + '.csv'
        tiles = Table.read(tilefile, format="csv")
        make_html_tiles(area,tiles,splus_tiles)


    return

def main():

    make_html_global()

if __name__ == '__main__':
      main()
