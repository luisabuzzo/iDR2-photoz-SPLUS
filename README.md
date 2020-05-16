# iDR2-photoz-SPLUS
Contents:

1. synthetic_mags.py - 
Calculate synthetic magnitudes from spectra, performs checks and extends spectra using templates, scales everything to match additional broad-band constraints, ...
NEEDS:
- IGMabs.py
- filter_folder.tar.gz
- template_folder.tar.gz
---
2. create_lephare_catalogues.py - 
This routine is used to read the master catalogues, extract the desired columns for the photometric redshift estimation (magnitudes and erros), identify the missing bands and label them accordingly, in order for LePhare to ignore them/ use upper limit errors.
---
3. footprint.py - 
This routine makes a plot of the sky with the tiles from the S-PLUS footprint and the iDR2 overlaid. Other catalogs of sources can be added as an option. The program needs some files contained in the footprint_data.tar.gz package.
---
4. extract_spectra.py - 
This routine is used to manage the BC03 .ised_ASCII files and produce a more convenient output, which can then be read to quickly extract spectra of the SED at a given age from the library. 
---
5. check_splus.py - 
This routine checks if a given location is covered by the S-PLUS and IDR2 footprints within a given search radius, and checks whether there are any objects in the iDR2 source catalog that are within a given search radius.
---
6. html.py - 
This makes the iDR2 data analysis overview websites. It runs two routines directly after each other. The routine make_html_global() makes the main site (./html/splus.html) with the footprint plot and links to all the individual iDR2 tiles. The routine make_html_tiles() makes the websites for each individual tile. It assumes that all the plots are present in the location ./html/AREA/ (with AREA = STRIPE82,HYDRA,SPLUS) and that figures have the names ./html/AREA/AREA-TILENAME_PLOTTYPE.png. The PLOTTYPE currently produced by idr2.py are:
- skyplot = location of the tile on the sky
- trilogy0 = color image
- skyradec = object locations in RA,DEC, with magnitude indicated by symbol size
- skyxy = object locations in X,Y
- skydens = object density map of the tile
- numbercounts12 = Number counts in the 12 filters with fits of the 2,3,5,10 sigma turn-off points
- numbercounts12_depths = Results from the numbercounts plots
- concentration = Concentrations in gri
- fwhm = FWHM versus magnitude with Gaussian fits 
- snr = SNR versus magnitude for the broadbands
- color = broad-band color diagrams with point sources indicated

Wanted/would like:
- make the footprint plot clickable such that one can also go to a pointing directly by plotting on a location in the plot.
---
7. run_lephare.py -
This routine runs all the steps of LePhare for different types of objects: star, quasars and galaxies.
---
8. template_coverage.py -
This routine plots the star, quasar and galaxy templates (COSMOS and BC03) in the colour-colour space. The routine contains the option to plot the S-PLUS objects or not.
Working on: plot of the templates in the colour-redshift space.
---
