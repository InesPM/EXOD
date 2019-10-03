#!/usr/bin/env python3
# coding=utf-8


########################################################################
#                                                                      #
# EXOD - EPIC-pn XMM-Newton Outburst Detector                          #
#                                                                      #
# Declaration of file names                                            #
#                                                                      #
# Inés Pastor Marazuela (2019) - ines.pastor.marazuela@gmail.com       #
#                                                                      #
########################################################################
"""
Declaration of the file names handled both by the detector and the renderer
"""

# Created files
LOG               = "log.txt"

VARIABILITY       = "variability_file.fits"
REGION            = "ds9_variable_sources.reg"

OUTPUT_IMAGE      = "variability.pdf"
OUTPUT_IMAGE_SRCS = "sources.pdf"
OUTPUT_IMAGE_ALL  = "variability_whole.pdf"

# Observation files

CLEAN_FILE        = "PN_clean.fits"
GTI_FILE          = "PN_gti.fits"
IMG_FILE          = "PN_image.fits"
RATE_FILE         = "PN_rate.fits"

# software installation paths

HEADAS = "/home/ines/astrosoft/heasoft-6.25/x86_64-pc-linux-gnu-libc2.27"
SAS    = "/home/ines/astrosoft/xmmsas_20180620_1732/setsas.sh"
