#!/usr/bin/env python3
# coding=utf-8


################################################################################
#                                                                              #
# EXOD - EPIC-pn XMM-Newton Outburst Detector                                  #
#                                                                              #
# DETECTOR from existing variability files                                     #
#                                                                              #
# Inés Pastor Marazuela (2019) - ines.pastor.marazuela@gmail.com               #
#                                                                              #
################################################################################
"""
Detecting variable sources for existing variability file by changing the detection level
"""

# Built-in imports

import sys
import os
import time
from functools import partial
from affine import Affine
from shutil import copy2

# Third-party imports

from math import *
from multiprocessing import Pool
from astropy.io import fits
from astropy import wcs
from astropy.coordinates import SkyCoord
from astropy import units as u
import numpy as np
import argparse

# Internal imports

from fits_extractor import *
from variability_utils import *
import file_names as FileNames
from file_utils import *
from detector import *

###
# Argument parser
###

parser = argparse.ArgumentParser()
parser.add_argument("evts", help="Path to the filtered events file", type=str)
parser.add_argument("init", help="Path to the initial input file", type=str)
parser.add_argument("fin", help="Path to the final file", type=str)
parser.add_argument("-obs", "--observation", dest="obs", help="Observation ID", default="", nargs='?', type=str)
parser.add_argument("-bs", "--box-size", dest="bs", help="Size of the detection box in pixel^2.", default=5, nargs='?', type=int)
parser.add_argument("-dl", "--detection-level", dest="dl", help="The number of times the median variability is required to trigger a detection.", default=10, nargs='?', type=float)
parser.add_argument("-tw", "--time-window", dest="tw", help="The duration of the time windows.", default=100.0, nargs='?', type=float)
parser.add_argument("-gtr", "--good-time-ratio", dest="gtr", help="Ratio of acceptability for a time window. Shall be between 0.0 and 1.0.", default=0.9, nargs='?', type=float)
parser.add_argument("-mta", "--max-threads-allowed", dest="mta", help="Maximal number of CPUs the program is allowed to use.", nargs='?', default=12, type=int)
parser.add_argument("-ol", "--output-log", dest="ol", help="tName of the general output file.", nargs='?', default="detected_sources", type=str)
args = parser.parse_args()

if args.init[-1] != '/' :
    args.init += '/'
if args.fin[-1] != '/' :
    args.fin += '/'

###
# Reading variability function
###

def read_file(file_path, counter=False, comment_token='#', separator=';') :
    """
    Function returning the content
    @param counter: True if it is the counters that are loaded, False if it is the variability
    @return: A list
    """

    data = np.zeros([12,64,200])
    ccd = 0; i = 0; j = 0
    with open(file_path) as f:
        for line in f:
            if comment_token not in line :
                data[ccd][i][j] = line
                if j < 199 :
                    j += 1
                elif j == 199 and i < 63 :
                    j = 0; i += 1
                elif j == 199 and i == 63 :
                    j = 0; i = 0; ccd += 1
    return data


###
# Starting
###

print('\n\t  RELEVEL Obs. {0}\n\t{1}'.format(args.obs,'-'*27))

original_time = time.time()
original_date = time.strftime("%d/%m/%Y %H:%M:%S", time.gmtime())

print('\n\n\tbox size = {0}\n\tdetection level = {1}\n\ttime window = {2}\n\tsignificance level = {3}\n'.format(args.bs, args.dl, args.tw, args.gtr))


###
# Opening existing files
###

log_f, var_f, var_per_tw_f, detected_var_areas_f, tws_f, detected_var_sources_f = open_files(sys.argv[3])
var_f.close()
output_log = open(args.ol, 'a')

copy2(args.init + FileNames.VARIABILITY, args.fin + FileNames.VARIABILITY)

v_matrix = read_file(args.fin + FileNames.VARIABILITY)

hdulist = fits.open(args.evts)
if len(hdulist) == 1 :
    os.system("bash /mnt/data/Ines/progs/filtering.sh /mnt/data/Ines/data/DR5 {0} /mnt/xmmcat/3xmm_pievli".format(args.obs))
header, dmin, dmax = extraction_info(args.evts)

log_f = open(sys.argv[2] + FileNames.LOG, "w+")
log_f.write('Command:\n\t')
log_f.write(' '.join([arg for arg in sys.argv]))
log_f.write('\n')
log_f.write("Creation of output files over.\n")


####
#
#   Detecting variable areas and sources
#
####

print('Detecting variable areas\t %.3f secondes' % (time.time() - original_time))
median = np.median([v_matrix[ccd][i][j] for ccd in range(12) for i in range(len(v_matrix[ccd])) for j in range(len(v_matrix[ccd][i]))])

# Avoiding a too small median value for detection
print('\nMedian\t\t',median)
if median < 0.75 :
    median = 0.75
    log_f.write('Median switched to 0.75. \n')

variable_areas = []

# Function for the pool of threads
variable_areas_detection_partial = partial(variable_areas_detection, median, args.bs, args.dl)
print('Box counts\t', args.dl * ((args.bs**2) * median))

# Performing parallel detection on each CCD
with Pool(args.mta) as p:
    variable_areas = p.map(variable_areas_detection_partial, v_matrix)

# Conversion pixels CCD en pixels ciel
w = wcs.WCS(header)
w.wcs.crpix = [header['REFXCRPX'], header['REFYCRPX']]
w.wcs.cdelt = [header['REFXCDLT']/15, header['REFYCDLT']]
w.wcs.crval = [header['REFXCRVL']/15, header['REFYCRVL']]
w.wcs.ctype = [header['REFXCTYP'], header['REFYCTYP']]
angle = header['PA_PNT']

# Writing sources to their files
cpt_source = 0

for ccd in range(12) :
    for source in variable_areas[ccd] :

        center_x = sum([p[0] for p in source]) / len(source)
        center_y = sum([p[1] for p in source]) / len(source)

        R = round(sqrt( (max([abs(p[0] - center_x) for p in source]))**2 + (max([abs(p[1] - center_y) for p in source]))**2 ))

        position = Source(cpt_source, ccd, center_x, center_y, R)
        (src_x, src_y) = transformation(position.x, position.y, dmax, dmin, angle)

        ra, dec = w.wcs_pix2world(src_x, src_y, 1)
        c   = SkyCoord(ra=ra*u.degree, dec=dec*u.degree)
        ra  = '{:.0f} {:.0f} {:.2f}'.format(c.ra.dms[0], c.ra.dms[1], c.ra.dms[2])
        dec = '{:.0f} {:.0f} {:.2f}'.format(c.dec.dms[0], c.dec.dms[1], c.dec.dms[2])

        # Avoiding bad pixels
        if [ccd, int(center_x)] not in [[4,11], [4,12], [4,13], [5,12], [10,28]] :
            cpt_source += 1
            detected_var_sources_f.write('{0};{1};{2};{3};{4};{5};{6}\n'.format(cpt_source, ccd + 1, center_x, center_y, R, ra, dec))

        for p in source :
            detected_var_areas_f.write('{0};{1};{2};{3}\n'.format(cpt_source, ccd + 1, p[0], p[1]))

print('Nb of sources\t',cpt_source)
log_f.write('Nb of sources\t{0}\n'.format(cpt_source))


####
#
#   End of the programme
#
####
output_log.write('{0} {1} {2} {3}\n'.format(args.obs, cpt_source, args.dl, args.tw))
log_f.write("# TOTAL EXECUTION TIME : %s seconds\n" % (time.time() - original_time))
close_files(log_f, var_f, var_per_tw_f, detected_var_areas_f, tws_f, detected_var_sources_f)
print(" # Total execution time Obs. {0} : {1} seconds\n".format(args.obs, (time.time() - original_time)))
