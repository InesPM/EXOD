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
import numpy as np
import matplotlib
matplotlib.use("Pdf")
from matplotlib import colors, image, transforms, gridspec
import matplotlib.pyplot as plt
#import matplotlib.gridspec as gridspec
from pylab import figure, cm
#from matplotlib.colors import LogNorm
from astropy import wcs
#from astropy.wcs import WCS
from astropy.io import fits


# Internal imports

import file_names as FileNames
from file_utils import *

def render_variability(var_file, output_file, sources=True, pars=None,
        maximum_value=None) :
    """
    Function producing output pdf plots from the computed variability.
    @param var_file: fits file containing variability and sources data
    @param output_file: The path to the PDF file to be created
    @param sources: If the detected sources are plotted or not
    @param pars: observation parameters
    @param maximum_value: The maximal value for the logarithmic scale
    """

    ###
    # Plot settings
    ###

    # Opening variability file
    hdulist = fits.open(var_file)

    data   = hdulist[0].data
    src    = hdulist[1].data
    header = hdulist[0].header

    hdulist.close()

    # Obtaining the WCS transformation parameters

    w = wcs.WCS(header)

    w.wcs.crpix = [header['REFXCRPX'], header['REFYCRPX']]
    w.wcs.cdelt = [header['REFXCDLT']/15, header['REFYCDLT']]
    w.wcs.crval = [header['REFXCRVL']/15, header['REFYCRVL']]
    w.wcs.ctype = [header['REFXCTYP'], header['REFYCTYP']]
    angle = header['PA_PNT']       # Degrees

    # Image limit
    dlim = [header['REFXLMIN'], header['REFXLMAX'], header['REFYLMIN'],
            header['REFYLMAX']]


    # Limite maximale de l'échelle des couleurs pour la normalisation par logarithme
    if maximum_value == None :
        maximum_value = max([max(tmp) for tmp in data])

    ###
    # Plotting
    ###

    # Plotting the variability data
    plt.subplot(111, projection=w)

    im = plt.imshow(data, cmap=cm.inferno, norm=colors.LogNorm(vmin=1.0,
            vmax=maximum_value), extent=dlim)

    ax = plt.gca()
    ax.set_facecolor('k')
    cbar = plt.colorbar(im)

    # Plotting the sources
    if sources :
        if len(src) != 0 :
            # Position of the sources
            plt.plot(src['X'], src['Y'], 'wo', alpha = 1, fillstyle='none')

    ra  = ax.coords[0]
    dec = ax.coords[1]

    ra.display_minor_ticks(True)
    dec.display_minor_ticks(True)
    ax.tick_params(axis='both', which='both', direction='in', color='w',
            width=1)

    # Labels
    plt.xlabel('RA', fontsize=10)
    plt.ylabel('DEC', fontsize=10)
    cbar.ax.set_ylabel('Variability', fontsize=10)

    # Title
    #if pars != None :
    plt.title('OBS {0}'.format(header['OBS_ID']), fontsize=14)
    plt.text(0.5, 0.95, "TW {0} s    DL {1}   BS {2}".format(header['TW'],
            header['DL'], header['BS']), color='white', fontsize=10,
            horizontalalignment='center', transform = ax.transAxes)

    plt.savefig(output_file, pad_inches=0, bbox_inches='tight', dpi=500)

########################################################################

def render_variability_all(var_file0, var_file1, var_file2, var_file3,
        output_file, sources=True, pars=None, maximum_value=10) :

    var_files = [var_file0, var_file1, var_file2, var_file3]

    # Starting loop on the different parameters
    fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(9.5,8))
    gs1 = gridspec.GridSpec(2, 2, wspace=0.05, hspace=0.05)

    for i in range(len(var_files)) :

        hdulist = fits.open(var_files[i])

        data   = hdulist[0].data
        src    = hdulist[1].data
        header = hdulist[0].header

        hdulist.close()

        # Obtaining the WCS transformation parameters

        w = wcs.WCS(header)

        w.wcs.crpix = [header['REFXCRPX'], header['REFYCRPX']]
        w.wcs.cdelt = [header['REFXCDLT']/15, header['REFYCDLT']]
        w.wcs.crval = [header['REFXCRVL']/15, header['REFYCRVL']]
        w.wcs.ctype = [header['REFXCTYP'], header['REFYCTYP']]
        angle = header['PA_PNT']       # Degrees

        # Image limit
        dlim = [header['REFXLMIN'], header['REFXLMAX'], header['REFYLMIN'],
                header['REFYLMAX']]

        # Plotting the variability data
        plt.subplot(gs1[i], projection=w)
        ax = plt.gca()
        im = plt.imshow(data/header['DL'], cmap=cm.inferno,
                norm=colors.LogNorm(vmin=0.1, vmax=1.0), extent=dlim)

        # Plotting the sources
        if sources :
            if len(src) != 0 :
                # Position of the sources
                plt.plot(src['X'], src['Y'], 'wo', alpha = 1, fillstyle='none')

        plt.text(0.5, 0.92, "TW {0} s    DL {1}   BS {2}".format(header['TW'],
                header['DL'], header['BS']), color='white', fontsize=10,
                horizontalalignment='center', transform = ax.transAxes)

        ra  = ax.coords[0]
        dec = ax.coords[1]

        if i > 1 :
            ra.set_axislabel('RA', fontsize=12)
        if i % 2 == 0 :
            dec.set_axislabel('DEC', fontsize=12)
        if i < 2 :
            ra.set_ticklabel_visible(False)
        if i % 2 != 0 :
            dec.set_ticklabel_visible(False)

        ra.display_minor_ticks(True)
        dec.display_minor_ticks(True)
        ax.tick_params(axis='both', which='both', direction='in',
                color='w', width=1)

        ax.set_facecolor('k')

    fig.subplots_adjust(right=0.77)
    cbar_ax = fig.add_axes([0.8, 0.11, 0.02, 0.77])
    cbar    = fig.colorbar(im, cax=cbar_ax)
    cbar.ax.set_ylabel('$\mathcal{V}$ / DL', fontsize=12)
    fig.suptitle('OBS {0}'.format(header['OBS_ID']), x=0.5, y = 0.93,
            fontsize=18)

    plt.savefig(output_file, pad_inches=0, dpi=500, bbox_inches='tight')


########################################################################

def ds9_renderer(var_file, reg_file) :
    """
    Function rendering variability with ds9
    @param var_file: variability fits file
    @param reg_file: region file
    """
    command = "ds9 {0} -scale linear -cmap bb -mode region -regionfile {1}".format(
            var_file, reg_file)
    process = subprocess.Popen(command.split(), stdout=subprocess.PIPE)
