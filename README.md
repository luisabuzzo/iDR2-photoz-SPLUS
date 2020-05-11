# iDR2-photoz-SPLUS

========== TILES.PY =============

=================================


========== HTML.PY =============

This makes the iDR2 data analysis overview websites. It runs two routines. The first makes the main site (./html/splus.html) with links to all the individual iDR2 tiles, and the second makes the websites for each individual tile. It assumes that all the plots are present in the location ./html/<AREA>/ (with <AREA> = STRIPE82,HYDRA,SPLUS) and that figures have the names ./html/<AREA>/<AREA>-TILENAME_<plot_type>.png. The <plot_type> currently produced by idr2.py are:
  
  skyplot = location of the tile on the sky
  
  trilogy0 = color image
  
  skyradec = object locations in RA,DEC, with magnitude indicated by symbol size
  
  skyxy = object locations in X,Y
  
  skydens = object density map of the tile
  
  numbercounts12 = Number counts in the 12 filters with fits of the 2,3,5,10 sigma turn-off points
  
  numbercounts12_depths = Results from the numbercounts plots
  
  concentration = Concentrations in gri
  
  fwhm = FWHM versus magnitude with Gaussian fits 
  
  snr = SNR versus magnitude for the broadbands
  
  color = broad-band color diagrams with point sources indicated

To be done: 

- 

=================================

