#!/usr/bin/env python3
# coding=utf-8

########################################################################
#                                                                      #
# EXOD - EPIC-pn XMM-Newton Outburst Detector                          #
#                                                                      #
# DETECTOR main programme                                              #
#                                                                      #
# Inés Pastor Marazuela (2019) - ines.pastor.marazuela@gmail.com       #
#                                                                      #
########################################################################
"""
Detector's main programme
"""


import sys
import os
import time
from functools import partial

# Third-party imports

from math import *
from multiprocessing import Pool
from astropy.io import fits
from astropy.table import Table
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

#############################################################################
#                                                                           #
# Parsing arguments                                                         #
#                                                                           #
#############################################################################

parser = argparse.ArgumentParser()
parser.add_argument("evts", help="Name of the clean observation file", type=str)
parser.add_argument("gti", help="Name of the GTI file", type=str)
parser.add_argument("path", help="Path to the folder containing the observation files", type=str)
parser.add_argument("out", help="Path to the folder where the output files will be stored", default=None, type=str)
parser.add_argument("-obs", "--observation", dest="obs", help="Observation ID", default=None, nargs='?', type=str)
parser.add_argument("-bs", "--box-size", dest="bs", help="Size of the detection box in pixel^2.", default=5, nargs='?', type=int)
parser.add_argument("-dl", "--detection-level", dest="dl", help="The number of times the median variability is required to trigger a detection.", default=10, nargs='?', type=float)
parser.add_argument("-tw", "--time-window", dest="tw", help="The duration of the time windows.", default=100.0, nargs='?', type=float)
parser.add_argument("-gtr", "--good-time-ratio", dest="gtr", help="Ratio of acceptability for a time window. Shall be between 0.0 and 1.0.", default=0.9, nargs='?', type=float)
parser.add_argument("-mta", "--max-threads-allowed", dest="mta", help="Maximal number of CPUs the program is allowed to use.", nargs='?', default=12, type=int)
parser.add_argument("-ol", "--output-log", dest="ol", help="tName of the general output file.", nargs='?', default="detected_sources", type=str)
args = parser.parse_args()

#############################################################################
#                                                                           #
# Functions                                                                 #
#                                                                           #
#############################################################################

def main_fct() :
    """
    Main function of the detector
    """
###
# Preliminaries
###

    # Counter for the overall execution time
    original_time = time.time()

    # Opening the output files
    log_f, var_f, var_per_tw_f, detected_var_areas_f, tws_f, detected_var_sources_f = open_files(args.out)
    output_log = open(args.ol, 'a')
    sys.stdout = Tee(sys.stdout, log_f)

# Replacing log_f.write by Tee class -> print

    print('Command:\n\t')
    print(' '.join([arg for arg in sys.argv]))
    print("\nCreation of output files over.\n")

    # Recovering the EVENTS list
    print('Recovering the EVENTS list\t %.3f seconds' % (time.time() - original_time))
    try :
        data = extraction_photons(args.evts)
        header, dmin, dmax = extraction_info(args.evts)

    except Exception as e:
        print("!!!!\nImpossible to extract photons. ABORTING.")
        close_files(log_f, var_f, var_per_tw_f, detected_var_areas_f, tws_f, detected_var_sources_f)
        exit(-2)

    # Defining header variability
    head_var_f = header
    head_var_f.append(card=('TW', args.tw, 'Time window'))
    head_var_f.append(card=('GTR', args.gtr, 'Good time ratio'))
    head_var_f.append(card=('DL', args.dl, 'Detection level'))
    head_var_f.append(card=('BS', args.bs, 'Box size'))
    data_var_f = Table(names=('VARIABILITY', 'RAWX', 'RAWY', 'CCDNR'), dtype=('f8', 'i2', 'i2', 'i2'))

    # Recovering GTI list
    try:
        print('Extracting data\t\t\t %.3f seconds' % (time.time() - original_time))
        gti_list = extraction_deleted_periods(args.gti)
        #print(len(gti_list))

    except Exception as e:
        print("!!!!\nImpossible to extract gti. ABORTING.")
        close_files(log_f, var_f, var_per_tw_f, detected_var_areas_f, tws_f, detected_var_sources_f)
        exit(-2)

    print("Extraction from FITS over.\n")

    time_windows = []
    t0_observation = min([evt['TIME'] for ccd in data for evt in ccd])
    tf_observation = max([evt['TIME'] for ccd in data for evt in ccd])

###
# Computing variability
###
    print('Computing variability\t\t %.3f seconds' % (time.time() - original_time))
    # Computing v_matrix
    v_matrix = []

    var_calc_partial = partial(variability_computation, gti_list, args.tw, args.gtr, t0_observation, tf_observation)

    with Pool(args.mta) as p:
        v_matrix = p.map(var_calc_partial, data)

    # Writing variability for each pixel
    for ccd in range(12) :
        for i in range(64):
            for j in range(200):
                #var_f.write(str(v_matrix[ccd][i][j]) + '\n')
                data_var_f.add_row([v_matrix[ccd][i][j], i+1, j+1, ccd+1])
    fits.writeto(var_f, data_var_f, header=head_var_f, overwrite=True)

###
# Detecting variable areas and sources
###

    print('Detecting variable areas\t %.3f seconds' % (time.time() - original_time))
    median = np.median([v_matrix[ccd][i][j] for ccd in range(12) for i in range(len(v_matrix[ccd])) for j in range(len(v_matrix[ccd][i]))])

    # Avoiding a too small median value for detection
    print('\nMedian\t\t{0}\n'.format(median))
    if median < 0.75 :
        median = 0.75
        print('Median switched to 0.75. \n')

    variable_areas = []

    # Currying the function for the pool of threads
    variable_areas_detection_partial = partial(variable_areas_detection, median, args.bs, args.dl)
    print('Box counts\t{0}'.format(args.dl * ((args.bs**2))))
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

    print('Nb of sources\t{0}\n'.format(cpt_source))


###
# End of the programme
###

    output_log.write('{0} {1} {2} {3}\n'.format(args.obs, cpt_source, args.dl, args.tw))
    output_log.close()
    #print("# TOTAL EXECUTION TIME : %s seconds\n" % (time.time() - original_time))
    close_files(log_f, var_f, var_per_tw_f, detected_var_areas_f, tws_f, detected_var_sources_f)

#############################################################################
#                                                                           #
# Main programme                                                            #
#                                                                           #
#############################################################################


if __name__ == '__main__':

    print('\n\t  DETECTOR Obs. {0}\n\t{1}'.format(args.obs,'-'*28))

    original_time = time.time()
    original_date = time.strftime("%d/%m/%Y %H:%M:%S", time.gmtime())

    print('\n\n\tbox size = {0}\n\tdetection level = {1}\n\ttime window = {2}\n\tsignificance level = {3}\n'.format(args.bs, args.dl, args.tw, args.gtr))

    main_fct()

    print(" # Total execution time OBS {0} : {1} seconds\n".format(args.obs, (time.time() - original_time)))
