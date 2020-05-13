# iDR2-photoz-SPLUS
Contents:

1. create_lephare_catalogues.py - 
Description...
---
2. footprint.py - 
This routine makes a plot of the sky with the tiles from the S-PLUS footprint and the iDR2 overlaid. Other catalogs of sources can be added as an option. The program needs some files contained in the footprint_data.tar.gz package.
---
3. extract_spectra.py - 
This routine is used to manage the BC03 .ised_ASCII files and produce a more convenient output, which can then be read to quickly extract spectra of the SED at a given age from the library. 
---
4. html.py - 
This makes the iDR2 data analysis overview websites. It runs two routines directly after each other. The routine make_html_global() makes the main site (./html/splus.html) with the footprint plot and links to all the individual iDR2 tiles. The routine make_html_tiles() makes the websites for each individual tile. It assumes that all the plots are present in the location ./html/AREA/ (with AREA = STRIPE82,HYDRA,SPLUS) and that figures have the names ./html/AREA/AREA-TILENAME_PLOTTYPE.png. The PLOTTYPE currently produced by idr2.py are:
  
-skyplot = location of the tile on the sky

-trilogy0 = color image

-skyradec = object locations in RA,DEC, with magnitude indicated by symbol size

-skyxy = object locations in X,Y

-skydens = object density map of the tile

-numbercounts12 = Number counts in the 12 filters with fits of the 2,3,5,10 sigma turn-off points

-numbercounts12_depths = Results from the numbercounts plots

-concentration = Concentrations in gri

-fwhm = FWHM versus magnitude with Gaussian fits 

-snr = SNR versus magnitude for the broadbands

-color = broad-band color diagrams with point sources indicated

Wanted/would like:
- make the footprint plot clickable such that one can also go to a pointing directly by plotting on a location in the plot.
---

