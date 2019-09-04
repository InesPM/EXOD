#!/usr/bin/env python3
# coding=utf-8

########################################################################
#                                                                      #
# Variabilitectron - Searching variability into XMM-Newton data        #
#                                                                      #
# RENDERER functions to be used with the DETECTOR                      #
#                                                                      #
# Inés Pastor Marazuela (2018) - ines.pastor.marazuela@gmail.com       #
#                                                                      #
########################################################################
"""
Renderer's main programme
"""

# Built-in imports

from os.path import sys
import os
import shutil

# Third-party imports

from math import *
#import scipy
#from scipy import ndimage
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

# Internal imports

import file_names as FileNames
from file_utils import *

def render_variability(data, sources, header, output_file, pars=None, maximum_value=None) :
    """
    Function rendering an from the matrix data.
    @param data: The matrix to render
    @param sources: The detected sources
    @param header: Header of the events file containing the information
    @param output_file: The path to the PDF file to be created
    @param pars: observation parameters
    @param maximum_value: The maximal value for the logarithmic scale
    """

    ###
    # Plot settings
    ###

    # Obtaining the WCS transformation parameters

    w = wcs.WCS(header)

    w.wcs.crpix = [header['REFXCRPX'], header['REFYCRPX']]
    w.wcs.cdelt = [header['REFXCDLT']/15, header['REFYCDLT']]
    w.wcs.crval = [header['REFXCRVL']/15, header['REFYCRVL']]
    w.wcs.ctype = [header['REFXCTYP'], header['REFYCTYP']]
    angle = header['PA_PNT']       # Degrees

    # Image limit
    dlim = [header['REFXLMIN'], header['REFXLMAX'], header['REFYLMIN'], header['REFYLMAX']]


    # Limite maximale de l'échelle des couleurs pour la normalisation par logarithme
    if maximum_value == None :
        maximum_value = max([max(tmp) for tmp in data])

    ###
    # Plotting
    ###

    # Plotting the variability data
    plt.subplot(111, projection=w)

    im = plt.imshow(data, cmap=cm.inferno, norm=LogNorm(vmin=1.0, vmax=maximum_value), extent=dlim)

    ax = plt.gca()
    ax.set_facecolor('k')
    cbar = plt.colorbar(im)

    # Plotting the sources
    if sources != None :
        if len(sources) != 0 :
            # Position of the sources
            plt.plot(sources['X'], sources['Y'], 'wo', alpha = 1, fillstyle='none')

    # Labels
    plt.xlabel('RA', fontsize=10)
    plt.ylabel('Dec', fontsize=10)
    cbar.ax.set_ylabel('Variability', fontsize=10)

    # Title
    if pars != None :
        plt.suptitle('OBS {0}'.format(pars['OBS_ID']), fontsize=14)
        plt.title("TW {0} s    DL {1}   BS {2}".format(pars['TW'], pars['DL'], pars['BS']), fontsize=10)

    plt.savefig(output_file, pad_inches=0, bbox_inches='tight', dpi=500)

########################################################################

def ds9_renderer(var_file, reg_file) :
    """
    Function rendering variability with ds9
    @param var_file: variability fits file
    @param reg_file: region file
    """
    command = "ds9 {0} -scale asinh -cmap bb -mode region -regionfile {1}".format(var_file, reg_file)
    process = subprocess.Popen(command.split(), stdout=subprocess.PIPE)
