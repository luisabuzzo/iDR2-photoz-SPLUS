# -*- coding: utf-8 -*-

import glob
import os, sys
from time import time, strftime

if getpass.getuser() == "luisabuzzo":
    home = "/home/luisabuzzo/Work/Master/SPLUS"
elif getpass.getuser() == "mlgbuzzo":
    home = "/home/mlgbuzzo/LePhare/"

os.environ['LEPHAREDIR'] = os.path.join(home,"lephare_dev")
os.environ['LEPHAREWORK'] = os.path.join(home,"lepharework")

out_path = os.path.join(home,"results/outputs/")

def main():
	config_files = glob.glob(os.path.join(home,"lephare_dev/config/*.para")) 
	for file in config_files:
		if "STAR" in file:
			sed_star = "$LEPHAREDIR/source/sedtolib -t S -c  " +str(file)
			os.system(sed_star)
			filters = "$LEPHAREDIR/source/filter -c  " + str(file)
			os.system(filters)
			mag_star = "$LEPHAREDIR/source/mag_star -c " + str(file)
			os.system(mag_star)
			zphota = "$LEPHAREDIR/source/zphota -c  " + str(file)
			os.system(zphota)

		elif "QSO" in file:
			sed_qso = "$LEPHAREDIR/source/sedtolib -t Q -c  " + str(file)
			os.system(sed_qso)
			filters = "$LEPHAREDIR/source/filter -c  " + str(file)
			os.system(filters)
			mag_qso = "$LEPHAREDIR/source/mag_gal -t Q -c  " + str(file) + " -EB_V 0"
			os.system(mag_qso)
			zphota = "$LEPHAREDIR/source/zphota -c  " + str(file)
			os.system(zphota)

		elif "GAL" in file:
			sed_gal = "$LEPHAREDIR/source/sedtolib -t G -c  " + str(file)
			os.system(sed_gal)
			filters = "$LEPHAREDIR/source/filter -c  " + str(file)
			os.system(filters)
			mag_gal = "$LEPHAREDIR/source/mag_gal -t G -c  " + str(file)
			os.system(mag_gal)
			zphota = "$LEPHAREDIR/source/zphota -c  " + str(file)
			os.system(zphota)

		else:
			sed_star = "$LEPHAREDIR/source/sedtolib -t S -c  " +str(file)
			os.system(sed_star)
			sed_qso = "$LEPHAREDIR/source/sedtolib -t Q -c  " + str(file)
			os.system(sed_qso)
			sed_gal = "$LEPHAREDIR/source/sedtolib -t G -c  " + str(file)
			os.system(sed_gal)
			filters = "$LEPHAREDIR/source/filter -c  " + str(file)
			os.system(filters)
			mag_star = "$LEPHAREDIR/source/mag_star -c " + str(file)
			os.system(mag_star)
			mag_qso = "$LEPHAREDIR/source/mag_gal -t Q -c  " + str(file) + " -EB_V 0"
			os.system(mag_qso)
			mag_gal = "$LEPHAREDIR/source/mag_gal -t G -c  " + str(file)
			os.system(mag_gal)
			zphota = "$LEPHAREDIR/source/zphota -c  " + str(file)
			os.system(zphota)

if __name__ == "__main__":
	main()
