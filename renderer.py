#!/usr/bin/env python3
# coding=utf-8

################################################################################
#                                                                              #
# Variabilitectron - Searching variability into XMM-Newton data                #
#                                                                              #
# RENDERER main programme                                                      #
#                                                                              #
# Inés Pastor Marazuela (2018) - ines.pastor.marazuela@gmail.com               #
#                                                                              #
################################################################################
"""
Renderer's main programme
"""

# Built-in imports

from os.path import sys
import os
import shutil

# Third-party imports

from math import *
import scipy
from scipy import ndimage
import numpy as np
import matplotlib
matplotlib.use("Pdf")
from matplotlib import colors, image, transforms
import matplotlib.pyplot as plt
from pylab import figure, cm
from matplotlib.colors import LogNorm
from astropy.io import fits
from astropy import wcs
from astropy.wcs import WCS
import argparse

# Internal imports

import file_names as FileNames
from file_utils import *

################################################################################
#                                                                              #
# Utilities                                                                    #
#                                                                              #
################################################################################

# Argument parser

parser = argparse.ArgumentParser()
parser.add_argument("out", help="Path to the folder where the output files will be stored", type=str)
parser.add_argument("evts", help="Path to the clean observation file", type=str)
parser.add_argument("-obs", "--observation", dest="obs", help="Observation ID", default="", nargs='?', type=str)
parser.add_argument("-dl", "--detection-level", dest="dl", help="The number of times the median variability is required to trigger a detection.", default=10, nargs='?', type=float)
parser.add_argument("-tw", "--time-window", dest="tw", help="The duration of the time windows.", default=100.0, nargs='?', type=float)
parser.add_argument("-bs", "--box-size", dest="bs", help="The size of the detection box.", default=5, nargs='?', type=int)
args = parser.parse_args()

# Functions

def preprocess_data(input_data) :
    """
    Function building a matrix for the output image from a list of data built that way ;
    [dot for ccd in data for rawx in ccd for dot in rawx]
    @param input_data: A list
    @return: A 384*400 matrix
    """

    data = []

    # Generating the basis of the array
    for i in range(64*6) :
        data.append([])

    # Building the ccd data list
    # CCD order to follow :
    # 6     reversed(9)
    # 5     reversed(8)
    # 4     reversed(7)
    # 1     reversed(10)
    # 2     reversed(11)
    # 3     reversed(12)

    # Appending CCD 9 data

    i = 575
    while i >= 512 :
        data[i % 64].extend(input_data[int(512 + 31.5 + (512 + 31.5 - i))])
        i -= 1

    # Appending CCD 8 data
    i = 511
    while i >= 448 :
        data[(i % 64) + 64].extend(input_data[int(448 + 31.5 + (448 + 31.5 - i))])
        i -= 1

    # Appending CCD 7 data
    i = 447
    while i >= 384 :
        data[(i % 64) + 128].extend(input_data[int(384 + 31.5 + (384 + 31.5 - i))])
        i -= 1

    # Appending CCD 10 data
    i = 639
    while i >= 576 :
        data[(i % 64) + 192].extend(input_data[int(576 + 31.5 + (576 + 31.5 - i))])
        i -= 1

    # Appending CCD 11 data
    i = 703
    while i >= 640 :
        data[(i % 64) + 256].extend(input_data[int(640 + 31.5 + (640 + 31.5 - i))])
        i -= 1

    # Appending CCD 12 data
    i = 767
    while i >= 704 :
        #print(int(704 + 31.5 + (704 + 31.5 - i)))
        data[(i % 64) + 320].extend(input_data[int(704 + 31.5 + (704 + 31.5 - i))])
        i -= 1

    # Appending CCD 6 data
    for i in range(320, 384) :
        data[i % 64].extend(reversed(input_data[i]))

    # Appending CCD 5 data
    for i in range(256, 320) :
        data[(i % 64) + 64].extend(reversed(input_data[i]))

    # Appending CCD 4 data
    for i in range(192, 256) :
        data[(i % 64) + 128].extend(reversed(input_data[i]))

    # Appending CCD 1 data
    for i in range(0, 64) :
        data[(i % 64) + 192].extend(reversed(input_data[i]))

    # Appending CCD 2 data
    for i in range(64, 128) :
        data[(i % 64) + 256].extend(reversed(input_data[i]))

    # Appending CCD 3 data
    for i in range(128, 192) :
        data[(i % 64) + 320].extend(reversed(input_data[i]))

    return data


################################################################################

def render_variability_whole_image(data, sources, output_file, obs, maximum_value=None) :
    """
    Function rendering an from the matrix data.
    @param data: The matrix to render
    @param sources: The detected sources
    @param output_file: The path to the PDF file to be created
    @param maximum_value: The maximal value for the logarithmic scale
    """

    # Obtaining the WCS transformation parameters
    fileopen = args.evts 
    hdu = fits.open(fileopen)
    w = wcs.WCS(hdu[0].header)
    header  = hdu[0].header
    dat     = fits.getdata(fileopen)

    w.wcs.crpix = [header['REFXCRPX'], header['REFYCRPX']]
    w.wcs.cdelt = [header['REFXCDLT']/15, header['REFYCDLT']]
    w.wcs.crval = [header['REFXCRVL']/15, header['REFYCRVL']]
    w.wcs.ctype = [header['REFXCTYP'], header['REFYCTYP']]
    angle = header['PA_PNT']       # Degrees

    # Conversion pixels CCD en pixels ciel
    dmin  = [min(dat['X']), min(dat['Y'])]
    dmax  = [max(dat['X']), max(dat['Y'])]

    # Limite maximale de l'échelle des couleurs pour la normalisation par logarithme
    if maximum_value == None :
        maximum_value = max([max(tmp) for tmp in data])

    # Plotting the variability data
    plt.subplot(111, projection=w)
    ax = plt.gca()
    ax.set_facecolor('k')
    rotated_img = scipy.ndimage.rotate(data, angle, reshape = True)
    im = plt.imshow(rotated_img, cmap=cm.inferno, norm=LogNorm(vmin=0.9, vmax=maximum_value), extent=[dmin[0], dmax[0], dmin[1], dmax[1]])
    cbar = plt.colorbar(im)
    plt.suptitle('Obs. {0}'.format(args.obs), fontsize=14)
    plt.title("TW : {0} s    DL : {1}   BS : {2}".format(args.tw, args.dl, args.bs), fontsize=10)
    plt.xlabel('RA', fontsize=10)
    plt.ylabel('Dec', fontsize=10)
    cbar.ax.set_ylabel('Variability', fontsize=10)

    # Plotting the sources
    if sources != None :
        # Position of the sources
        src_x=[] ; src_y=[]
        for src in sources :
            (s_x, s_y) = transformation(int(src.x), int(src.y), dmax, dmin, angle)
            src_x.append(s_x); src_y.append(s_y)
        plt.plot(src_x, src_y, 'wo', alpha = 1, fillstyle='none')

    plt.savefig(output_file, pad_inches=0, bbox_inches='tight', dpi=500)



################################################################################


def print_help_r() :
    """
    Prints the help of the renderer
    """
    print("Use:")
    print("\tpython3 renderer.py <PATH_TO_DETECTOR_OUTPUT> [--gif]")
    print("Arguments:")
    print("\t--gif\tGenerates an animation of the counted events per pixel.\n\t\t\tWARN: Can be very memory consuming.")

################################################################################
#                                                                              #
# Main programme                                                               #
#                                                                              #
################################################################################


def main_r() :
    """
    Main function of the renderer
    """

    # Title
    print('\n\t  RENDERER Obs. {0}\n\t{1}'.format(args.obs,'-'*28))

    INPUT_FOLDER = args.out

    if INPUT_FOLDER[-1] != "/" :
        INPUT_FOLDER += "/"

    var_per_px = read_from_file(INPUT_FOLDER + FileNames.VARIABILITY)
    sources = read_sources_from_file(INPUT_FOLDER + FileNames.VARIABLE_SOURCES)
    data_matrix = preprocess_data(var_per_px)
    render_variability_whole_image(data_matrix, None, INPUT_FOLDER + FileNames.OUTPUT_IMAGE, args.obs)
    render_variability_whole_image(data_matrix, sources, INPUT_FOLDER + FileNames.OUTPUT_IMAGE_SRCS, args.obs)

################################################################################

if __name__ == '__main__':
    main_r()
