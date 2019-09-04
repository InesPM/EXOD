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
import matplotlib.pyplot as plt

# Internal imports

from fits_extractor import *
from variability_utils import *
import file_names as FileNames
from file_utils import *
from renderer import *

########################################################################
#                                                                      #
# Parsing arguments                                                    #
#                                                                      #
########################################################################

parser = argparse.ArgumentParser()

# Path to files
parser.add_argument("-evts", help="Name of the clean observation file", type=str, nargs='?', default=FileNames.CLEAN_FILE)
parser.add_argument("-gti", help="Name of the GTI file", type=str, nargs='?', default=FileNames.GTI_FILE)
parser.add_argument("-img", help="Name of the image file", type=str, nargs='?', default=FileNames.IMG_FILE)
parser.add_argument("-path", help="Path to the folder containing the observation files", type=str)
parser.add_argument("-out", help="Path to the folder where the output files will be stored", default=None, type=str)

# Variability parameters
parser.add_argument("-bs", "--box-size", dest="bs", help="Size of the detection box in pixel^2.", default=5, nargs='?', type=int)
parser.add_argument("-dl", "--detection-level", dest="dl", help="The number of times the median variability is required to trigger a detection.", default=10, nargs='?', type=float)
parser.add_argument("-tw", "--time-window", dest="tw", help="The duration of the time windows.", default=100.0, nargs='?', type=float)
parser.add_argument("-gtr", "--good-time-ratio", dest="gtr", help="Ratio of acceptability for a time window. Shall be between 0.0 and 1.0.", default=0.9, nargs='?', type=float)
parser.add_argument("-mta", "--max-threads-allowed", dest="mta", help="Maximal number of CPUs the program is allowed to use.", nargs='?', default=12, type=int)

# Arguments set by default
parser.add_argument("-creator", dest="creator", help="User creating the variability files", nargs='?', default=os.environ['USER'], type=str)
parser.add_argument("-obs", "--observation", dest="obs", help="Observation ID", default=None, nargs='?', type=str)

# Boolean flags
parser.add_argument('--render', help='Plot variability output, produce pdf', action='store_true')
parser.add_argument('--ds9', help='Plot variability output in emerging ds9 window', action='store_true')

args = parser.parse_args()

# Modifying arguments
if args.path[-1] != '/' :
    args.path = args.path + '/'
if args.out[-1] != '/' :
    args.out = args.out + '/'
if args.out == None :
    args.out = args.path + '{}_{}_{}_{}'.format(args.dl, args.tw, args.bs, args.gtr)
args.evts = args.path + args.evts
args.gti  = args.path + args.gti
args.img  = args.path + args.img

########################################################################
#                                                                      #
# Functions                                                            #
#                                                                      #
########################################################################

def main_fct() :
    """
    Main function of the detector
    """
###
# Preliminaries
###
    #print(' '.join([arg for arg in sys.argv]))
    print(vars(args))

    print("""
        DETECTION LEVEL = {0}
        TIME WINDOW     = {1}
        BOX SIZE        = {2}
        GOOD TIME RATIO = {3}
        """.format(args.dl, args.tw, args.bs, args.gtr))

    # Counter for the overall execution time
    original_time = time.time()

    # Opening the output files
    log_f, var_f, reg_f = open_files(args.out)
    original = sys.stdout
    sys.stdout = Tee(sys.stdout, log_f)

# Replacing log_f.write by Tee class -> print

    # Recovering the EVENTS list
    print(' Recovering the events list\t {:7.2f} s'.format(time.time() - original_time))
    try :
        data, header = extraction_photons(args.evts)

        if args.obs == None :
            args.obs = header['OBS_ID']

    except Exception as e:
        print(" !!!!\nImpossible to extract photons. ABORTING.")
        close_files(log_f, var_f, var_per_tw_f, detected_var_areas_f, tws_f, detected_var_sources_f)
        exit(-2)

    # Parameters ready
    params = {
              "CREATOR" : args.creator,
              "DATE"    : time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime()),
              "OBS_ID"  : args.obs,
              "TW"      : args.tw,
              "GTR"     : args.gtr,
              "DL"      : args.dl,
              "BS"      : args.bs
             }

    # Recovering GTI list
    try:
        print(' Extracting data\t\t {:7.2f} s'.format(time.time() - original_time))
        gti_list = extraction_deleted_periods(args.gti)

    except Exception as e:
        print(" !!!!\nImpossible to extract gti. ABORTING.")
        close_files(log_f, var_f, var_per_tw_f, detected_var_areas_f, tws_f, detected_var_sources_f)
        exit(-2)

    time_windows = []
    t0_observation = min([evt['TIME'] for ccd in data for evt in ccd])
    tf_observation = max([evt['TIME'] for ccd in data for evt in ccd])

###
# Computing variability
###
    print(' Computing variability\t\t {:7.2f} s'.format(time.time() - original_time))
    # Computing v_matrix
    v_matrix = []

    var_calc_partial = partial(variability_computation, gti_list, args.tw, args.gtr, t0_observation, tf_observation)

    with Pool(args.mta) as p:
        v_matrix = p.map(var_calc_partial, data)

    # Aplying CCD configuration
    data_v = ccd_config(v_matrix)
    img_v  = data_transformation(data_v, header)

###
# Detecting variable areas and sources
###

    print(' Detecting variable sources\t {:7.2f} s'.format(time.time() - original_time))
    median = np.median([v_matrix[ccd][i][j] for ccd in range(12) for i in range(len(v_matrix[ccd])) for j in range(len(v_matrix[ccd][i]))])

    # Avoiding a too small median value for detection
    print('\n\tMedian\t\t{0}'.format(median))
    if median < 0.75 :
        median = 0.75
        print(' Median switched to 0.75. \n')

    variable_areas = []

    # Currying the function for the pool of threads
    variable_areas_detection_partial = partial(variable_areas_detection, median, args.bs, args.dl)
    print('\tBox counts\t{0}'.format(args.dl * ((args.bs**2))))
    # Performing parallel detection on each CCD
    with Pool(args.mta) as p:
        variable_areas = p.map(variable_areas_detection_partial, v_matrix)

    # Variable sources
    sources = variable_sources_position(variable_areas, args.obs, args.path, reg_f, log_f, args.img)

    print('\tNb of sources\t{0}\n'.format(len(sources)))

###
# Writing data to fits file
###

    fits_writer(img_v, sources, args.img, params, var_f)

###
# Plotting variability
###

    # Renderer
    if args.render :

        print(' Rendering variability image\t {:7.2f} s'.format(time.time() - original_time))

        render_variability(img_v, None, header, args.out + FileNames.OUTPUT_IMAGE, params, maximum_value=10)
        render_variability(img_v, sources, header, args.out + FileNames.OUTPUT_IMAGE_SRCS, params, maximum_value=10)

    # ds9
    if args.ds9 :
        ds9_renderer(var_f, reg_f)

###
# End of the programme
###

    #print('{0} {1} {2} {3}\n'.format(args.obs, len(sources), args.dl, args.tw))
    print(" # Total execution time OBS {0} : {1:.2f} s\n".format(args.obs, (time.time() - original_time)))
    log_f.close()

########################################################################
#                                                                      #
# Main programme                                                       #
#                                                                      #
########################################################################


if __name__ == '__main__':

    original_time = time.time()
    original_date = time.strftime("%d/%m/%Y %H:%M:%S", time.gmtime())

    main_fct()
