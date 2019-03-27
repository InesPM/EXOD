#!/usr/bin/env python3
# coding=utf-8

'''
Getting rate from the events file header
'''

# Built-il imports
import sys
import os
from astropy.io import fits


# Extracting the data
if sys.argv[1] == "-fbktsr" :
	hdulist = fits.open(sys.argv[2])
	header  = hdulist['RATE'].header
	RATE    = header['FLCUTTHR']
	print(RATE)
	sys.exit(0)
	hdulist.close()

elif sys.argv[1] == "-mode" :
	hdulist  = fits.open(sys.argv[2])
	header   = hdulist[0].header
	DATAMODE = header['DATAMODE']
	SUBMODE  = header['SUBMODE']
	print(DATAMODE, SUBMODE)
	sys.exit(0)
	hdulist.close()
