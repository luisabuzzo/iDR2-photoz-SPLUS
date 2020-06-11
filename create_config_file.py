import numpy as np
import matplotlib.pyplot as plt
import glob
import os, sys
import time
from scipy import stats
from scipy.stats import gaussian_kde
import pandas as pd
from matplotlib import rc
import re
import getpass
import re
from astropy.io import fits
from astropy.io import ascii

##################################################### DEFINE MAIN DIRECTORIES #########################################################

if getpass.getuser() == "luisabuzzo":
    home = "/home/luisabuzzo/Work/Master/SPLUS"
elif getpass.getuser() == "mlgbuzzo":
    home = "/home/mlgbuzzo/LePhare/"

os.environ['LEPHAREDIR'] = os.path.join(home, "lephare_dev")
os.environ['LEPHAREWORK'] = os.path.join(home, "lepharework")

#######################################################################################################################################

def main():
	source = "GALAXY"
	gal_sed = "$LEPHAREDIR/sed/GAL/COSMOS_SED/COSMOS_MOD.list"
	cat_in = home +"/results/catalogues/LePhare_SPLUS_all_iDR2_matchsdss_auto.dat"
	cat_out = home +"/results/catalogues/SPLUS_all_LePhare.out"
	cat_source = "LONG"
	gallibin = "LIB_SPLUS_ALL"
	gallibout = "GAL_SPLUS_ALL"
	adapt_meth = "1" # 1 for mag, 2 for redshift, 3 for model
	filter_list = "splus/uJAVA.dat,splus/F378.dat,splus/F395.dat,splus/F410.dat,splus/F430.dat,splus/gSDSS.dat,splus/F515.dat,splus/rSDSS.dat,splus/F660.dat,splus/iSDSS.dat,splus/F861.dat,splus/zSDSS.dat"
	filter_file = "SPLUS_all.filt"
	z_step = "0.01,1.0,0.1"
	cosmology = "70,0.3,0.7"
	mod_extinc = "13,23,23,31,23,31,23,31"
	extinct_law = "SMC_prevot.dat,SB_calzetti.dat,SB_calzetti_bump1.dat,SB_calzetti_bump2.dat"
	ebv = "0.000,0.05,0.100,0.150,0.200,0.250,0.300,0.400,0.500"
	emlines = "YES"
	lib_ascii = "YES"
	para_out = "$LEPHAREDIR/config/zphot_output.para"
	bd_scale = "1"
	add_emlines = "YES"
	abs_mag = "-10.,-26."
	mag_ref = "8"
	z_range = "0.001,1.0"
	ebs_range = "0,0.5"
	prior = "8,6,8"  
	z_fix = "NO"
	z_interp = "YES"
	dz_win = "0.25"
	min_threshold = "0.1"
	abs_mag_method = "1"
	abs_mag_context = "-1"
	abs_mag_reference = "4"
	abs_mag_filt = "1,2,3,4"
	abs_mag_zbin = "0,0.5,1,1.5,2,3,3.5,4"
	spec_out = "NO"
	chi2_out = "NO"
	pdf_out = "NONE"
	col_num = "3"
	col_sigma = "3"
	autoadapt = "NO"
	adapt_bands = "8,6,8"
	adapt_maglim = "18.0,19.5"
	adapt_meth = "1"
	adapt_context = "4095.0"   
	adapt_zbin = "0.0001,0.1"
	adapt_modbin = "0,1000"
	error_adapt = "NO"

	if source == "GALAXY":
		with open(home + "/lephare_dev/config/SPLUS_"+str(source)+".para",'w') as file:
			file.write('##############################################################################\n'
	'#                CREATION OF LIBRARIES FROM SEDs List                        #\n'
	'# $LEPHAREDIR/source/sedtolib -t (S/Q/G) -c $LEPHAREDIR/config/zphot.para    #\n'
	'# help : $LEPHAREDIR/source/sedtolib -h (or -help)                           #\n'
	'##############################################################################\n'
	'#\n'
	'#------      GALAXY LIBRARY (ASCII or BINARY SEDs) \n'
	'GAL_SED	   '+str(gal_sed)+'  # GAL list (full path)\n'
	'GAL_FSCALE	1.					# Arbitrary Flux Scale\n'
	'GAL_LIB		'+str(gallibin)+'				# Bin. GAL LIBRARY ->\n'
	'							# $LEPHAREWORK/lib_bin\n'
	'#SEL_AGE   $LEPHAREDIR/sed/GAL/HYPERZ/AGE_GISSEL_ALL.dat # Age list(full path)\n'
	'							# (def=NONE)	\n'
	'AGE_RANGE  0.,13.e9                                     # Age Min-Max in yr\n'
	'FILTER_LIST	'+str(filter_list)+'\n'
	'TRANS_TYPE		0\n'
	'FILTER_CALIB	0\n'
	'FILTER_FILE	'+str(filter_file)+'\n'
	'GAL_LIB_IN		'+str(gallibin)+'\n'
	'GAL_LIB_OUT	'+str(gallibout)+'\n'
	'MAGTYPE		AB\n'
	'Z_STEP         '+str(z_step)+'\n'
	'COSMOLOGY		'+str(cosmology)+'\n'
	'MOD_EXTINC 	'+str(mod_extinc)+'\n'
	'EXTINC_LAW		'+str(extinct_law)+'\n'
	'EB_V           '+str(ebv)+' # E(B-V) (<50 values)\n'
	'EM_LINES 	    '+str(emlines)+'\n'
	'LIB_ASCII 		'+str(lib_ascii)+'\n'
	'CAT_IN   		'+str(cat_in)+'\n'
	'INP_TYPE     		M	\n'
	'CAT_MAG      		AB  \n'
	'CAT_FMT      		MEME  	\n'
	'CAT_LINES    		1,10000 \n'
	'CAT_TYPE     	    '+str(cat_source)+'	\n'
	'CAT_OUT	     	'+str(cat_out)+'\n'
	'PARA_OUT     		'+str(para_out)+' \n'
	'BD_SCALE        	'+str(bd_scale)+'	    \n'
	'GLB_CONTEXT     	-1	\n'
	'ERR_FACTOR      	1.0\n'
	'ERR_SCALE			-1.   \n'  
	'ZPHOTLIB        	'+str(gallibout)+'   \n'
	'ADD_EMLINES     	'+str(add_emlines)+'\n'
	'FIR_LIB         	NONE\n'
	'FIR_LMIN        	7.0 \n'
	'FIR_CONT        	-1\n'
	'FIR_SCALE       	-1\n'
	'FIR_FREESCALE   	NO  \n' 
	'FIR_SUBSTELLAR  	NO\n'
	'PHYS_LIB        	NONE \n'
	'PHYS_CONT       	-1\n'
	'PHYS_SCALE      	-1\n'
	'PHYS_NMAX       	100000\n'
	'MAG_ABS 		'+str(abs_mag)+'\n'
	'MAG_REF 		'+str(mag_ref)+'\n'
	'Z_RANGE         	'+str(z_range)+' \n'  
	'EBV_RANGE       	'+str(ebs_range)+' \n'
	'NZ_PRIOR     '+str(prior)+'       \n'    
	'ZFIX			'+str(z_fix)+'	\n'
	'Z_INTERP		'+str(z_interp)+'\n'
	'DZ_WIN          	'+str(dz_win)+' \n'    
	'MIN_THRES       	'+str(min_threshold)+'   \n'
	'MABS_METHOD		'+str(abs_mag_method)+'	\n'
	'MABS_CONTEXT    	'+str(abs_mag_context)+'  \n'
	'MABS_REF		'+str(abs_mag_reference)+'\n'
	'MABS_FILT       	'+str(abs_mag_filt)+'  \n'
	'MABS_ZBIN       	'+str(abs_mag_zbin)+'\n'
	'SPEC_OUT		'+str(spec_out)+'	\n'
	'CHI2_OUT        	'+str(chi2_out)+' \n'   
	'PDZ_OUT         	'+str(pdf_out)+'  \n'    
	'PDZ_MABS_FILT   	2,10,14\n'
	'FAST_E		NO 	\n'
	'COL_NUM			'+str(col_num)+' \n'	
	'COL_SIGMA		'+str(col_sigma)+'	\n'
	'COL_SEL			AND	\n'
	'AUTO_ADAPT		'+str(autoadapt)+'\n'
	'ADAPT_BAND 		'+str(adapt_bands)+'	\n'
	'ADAPT_LIM       	'+str(adapt_maglim)+'	\n'
	'ADAPT_POLY		1	\n'
	'ADAPT_METH      	'+str(adapt_meth)+'\n'
	'ADAPT_CONTEXT  		'+str(adapt_context)+'\n'    
	'ADAPT_ZBIN     		'+str(adapt_zbin)+'\n'      
	'ADAPT_MODBIN   		'+str(adapt_modbin)+'\n'
	'ERROR_ADAPT    		'+str(error_adapt)+'\n')

	if source == "QSO":
		with open(home + "/lephare_dev/config/SPLUS_"+str(source)+".para",'w') as file:
			file.write('##############################################################################\n'
	'#                CREATION OF LIBRARIES FROM SEDs List                        #\n'
	'# $LEPHAREDIR/source/sedtolib -t (S/Q/G) -c $LEPHAREDIR/config/zphot.para    #\n'
	'# help : $LEPHAREDIR/source/sedtolib -h (or -help)                           #\n'
	'##############################################################################\n'
	'#\n'
	'#------      QSO LIBRARY (ASCII SEDs) \n'
	'QSO_SED		'+str(qso_sed)+'       # QSO list (full path)\n'
	'QSO_FSCALE	1					# Arbitrary Flux Scale \n'
	'QSO_LIB		'+str(qsolibin)+'					# Bin. QSO LIBRARY ->\n'
	'							# $LEPHAREWORK/lib_bin\n'
	'#SEL_AGE   $LEPHAREDIR/sed/GAL/HYPERZ/AGE_GISSEL_ALL.dat # Age list(full path)\n'
	'							# (def=NONE)	\n'
	'AGE_RANGE  0.,13.e9                                     # Age Min-Max in yr\n'
	'FILTER_LIST	'+str(filter_list)+'\n'
	'TRANS_TYPE		0\n'
	'FILTER_CALIB	0\n'
	'FILTER_FILE	'+str(filter_file)+'\n'
	'GAL_LIB_IN		'+str(qsolibin)+'\n'
	'GAL_LIB_OUT	'+str(qsolibout)+'\n'
	'MAGTYPE		AB\n'
	'Z_STEP         '+str(z_step)+'\n'
	'COSMOLOGY		'+str(cosmology)+'\n'
	'MOD_EXTINC 	'+str(mod_extinc)+'\n'
	'EXTINC_LAW		'+str(extinct_law)+'\n'
	'EB_V           '+str(ebv)+' # E(B-V) (<50 values)\n'
	'EM_LINES 	    '+str(emlines)+'\n'
	'LIB_ASCII 		'+str(lib_ascii)+'\n'
	'CAT_IN   		'+str(cat_in)+'\n'
	'INP_TYPE     		M	\n'
	'CAT_MAG      		AB  \n'
	'CAT_FMT      		MEME  	\n'
	'CAT_LINES    		1,10000 \n'
	'CAT_TYPE     	    '+str(cat_source)+'	\n'
	'CAT_OUT	     	'+str(cat_out)+'\n'
	'PARA_OUT     		'+str(para_out)+' \n'
	'BD_SCALE        	'+str(bd_scale)+'	    \n'
	'GLB_CONTEXT     	-1	\n'
	'ERR_FACTOR      	1.0\n'
	'ERR_SCALE			-1.   \n'  
	'ZPHOTLIB        	'+str(qsolibout)+'   \n'
	'ADD_EMLINES     	'+str(add_emlines)+'\n'
	'FIR_LIB         	NONE\n'
	'FIR_LMIN        	7.0 \n'
	'FIR_CONT        	-1\n'
	'FIR_SCALE       	-1\n'
	'FIR_FREESCALE   	NO  \n' 
	'FIR_SUBSTELLAR  	NO\n'
	'PHYS_LIB        	NONE \n'
	'PHYS_CONT       	-1\n'
	'PHYS_SCALE      	-1\n'
	'PHYS_NMAX       	100000\n'
	'MAG_ABS 		'+str(abs_mag)+'\n'
	'MAG_REF 		'+str(mag_ref)+'\n'
	'Z_RANGE         	'+str(z_range)+' \n'  
	'EBV_RANGE       	'+str(ebs_range)+' \n'
	'NZ_PRIOR     '+str(prior)+'       \n'    
	'ZFIX			'+str(z_fix)+'	\n'
	'Z_INTERP		'+str(z_interp)+'\n'
	'DZ_WIN          	'+str(dz_win)+' \n'    
	'MIN_THRES       	'+str(min_threshold)+'   \n'
	'MABS_METHOD		'+str(abs_mag_method)+'	\n'
	'MABS_CONTEXT    	'+str(abs_mag_context)+'  \n'
	'MABS_REF		'+str(abs_mag_reference)+'\n'
	'MABS_FILT       	'+str(abs_mag_filt)+'  \n'
	'MABS_ZBIN       	'+str(abs_mag_zbin)+'\n'
	'SPEC_OUT		'+str(spec_out)+'	\n'
	'CHI2_OUT        	'+str(chi2_out)+' \n'   
	'PDZ_OUT         	'+str(pdf_out)+'  \n'    
	'PDZ_MABS_FILT   	2,10,14\n'
	'FAST_E		NO 	\n'
	'COL_NUM			'+str(col_num)+' \n'	
	'COL_SIGMA		'+str(col_sigma)+'	\n'
	'COL_SEL			AND	\n'
	'AUTO_ADAPT		'+str(autoadapt)+'\n'
	'ADAPT_BAND 		'+str(adapt_bands)+'	\n'
	'ADAPT_LIM       	'+str(adapt_maglim)+'	\n'
	'ADAPT_POLY		1	\n'
	'ADAPT_METH      	'+str(adapt_meth)+'\n'
	'ADAPT_CONTEXT  		'+str(adapt_context)+'\n'    
	'ADAPT_ZBIN     		'+str(adapt_zbin)+'\n'      
	'ADAPT_MODBIN   		'+str(adapt_modbin)+'\n'
	'ERROR_ADAPT    		'+str(error_adapt)+'\n')

	if source == "STAR":
		with open(home + "/lephare_dev/config/SPLUS_"+str(source)+".para",'w') as file:
			file.write('##############################################################################\n'
	'#                CREATION OF LIBRARIES FROM SEDs List                        #\n'
	'# $LEPHAREDIR/source/sedtolib -t (S/Q/G) -c $LEPHAREDIR/config/zphot.para    #\n'
	'# help : $LEPHAREDIR/source/sedtolib -h (or -help)                           #\n'
	'##############################################################################\n'
	'#\n'
	'#------      STELLAR LIBRARY (ASCII SEDs)\n'
	'STAR_SED	'+str(star_sed)+'	# STAR list (full path)\n'
	'STAR_FSCALE	1.				# Arbitrary Flux Scale \n'
	'STAR_LIB	'+str(starlibin)+'				# Bin. STAR LIBRARY ->\n'
	'							# $LEPHAREWORK/lib_bin\n'
	'#SEL_AGE   $LEPHAREDIR/sed/GAL/HYPERZ/AGE_GISSEL_ALL.dat # Age list(full path)\n'
	'							# (def=NONE)	\n'
	'AGE_RANGE  0.,13.e9                                     # Age Min-Max in yr\n'
	'FILTER_LIST	'+str(filter_list)+'\n'
	'TRANS_TYPE		0\n'
	'FILTER_CALIB	0\n'
	'FILTER_FILE	'+str(filter_file)+'\n'
	'GAL_LIB_IN		'+str(starlibin)+'\n'
	'GAL_LIB_OUT	'+str(starlibout)+'\n'
	'MAGTYPE 		AB\n'
	'Z_STEP         '+str(z_step)+'\n'
	'COSMOLOGY		'+str(cosmology)+'\n'
	'MOD_EXTINC 	'+str(mod_extinc)+'\n'
	'EXTINC_LAW		'+str(extinct_law)+'\n'
	'EB_V           '+str(ebv)+' # E(B-V) (<50 values)\n'
	'EM_LINES 	    '+str(emlines)+'\n'
	'LIB_ASCII 		'+str(lib_ascii)+'\n'
	'CAT_IN   		'+str(cat_in)+'\n'
	'INP_TYPE     		M	\n'
	'CAT_MAG      		AB  \n'
	'CAT_FMT      		MEME  	\n'
	'CAT_LINES    		1,10000 \n'
	'CAT_TYPE     	    '+str(cat_source)+'	\n'
	'CAT_OUT	     	'+str(cat_out)+'\n'
	'PARA_OUT     		'+str(para_out)+' \n'
	'BD_SCALE        	'+str(bd_scale)+'	    \n'
	'GLB_CONTEXT     	-1	\n'
	'ERR_FACTOR      	1.0\n'
	'ERR_SCALE			-1.   \n'  
	'ZPHOTLIB        	'+str(starlibout)+'   \n'
	'ADD_EMLINES     	'+str(add_emlines)+'\n'
	'FIR_LIB         	NONE\n'
	'FIR_LMIN        	7.0 \n'
	'FIR_CONT        	-1\n'
	'FIR_SCALE       	-1\n'
	'FIR_FREESCALE   	NO  \n' 
	'FIR_SUBSTELLAR  	NO\n'
	'PHYS_LIB        	NONE \n'
	'PHYS_CONT       	-1\n'
	'PHYS_SCALE      	-1\n'
	'PHYS_NMAX       	100000\n'
	'MAG_ABS 		'+str(abs_mag)+'\n'
	'MAG_REF 		'+str(mag_ref)+'\n'
	'Z_RANGE         	'+str(z_range)+' \n'  
	'EBV_RANGE       	'+str(ebs_range)+' \n'
	'NZ_PRIOR     '+str(prior)+'       \n'    
	'ZFIX			'+str(z_fix)+'	\n'
	'Z_INTERP		'+str(z_interp)+'\n'
	'DZ_WIN          	'+str(dz_win)+' \n'    
	'MIN_THRES       	'+str(min_threshold)+'   \n'
	'MABS_METHOD		'+str(abs_mag_method)+'	\n'
	'MABS_CONTEXT    	'+str(abs_mag_context)+'  \n'
	'MABS_REF		'+str(abs_mag_reference)+'\n'
	'MABS_FILT       	'+str(abs_mag_filt)+'  \n'
	'MABS_ZBIN       	'+str(abs_mag_zbin)+'\n'
	'SPEC_OUT		'+str(spec_out)+'	\n'
	'CHI2_OUT        	'+str(chi2_out)+' \n'   
	'PDZ_OUT         	'+str(pdf_out)+'  \n'    
	'PDZ_MABS_FILT   	2,10,14\n'
	'FAST_E		NO 	\n'
	'COL_NUM			'+str(col_num)+' \n'	
	'COL_SIGMA		'+str(col_sigma)+'	\n'
	'COL_SEL			AND	\n'
	'AUTO_ADAPT		'+str(autoadapt)+'\n'
	'ADAPT_BAND 		'+str(adapt_bands)+'	\n'
	'ADAPT_LIM       	'+str(adapt_maglim)+'	\n'
	'ADAPT_POLY		1	\n'
	'ADAPT_METH      	'+str(adapt_meth)+'\n'
	'ADAPT_CONTEXT  		'+str(adapt_context)+'\n'    
	'ADAPT_ZBIN     		'+str(adapt_zbin)+'\n'      
	'ADAPT_MODBIN   		'+str(adapt_modbin)+'\n'
	'ERROR_ADAPT    		'+str(error_adapt)+'\n')

if __name__ == "__main__":
	main()