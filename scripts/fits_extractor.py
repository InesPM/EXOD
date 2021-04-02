#!/usr/bin/env python3
# coding=utf-8


########################################################################
#                                                                      #
# EXOD - EPIC-pn XMM-Newton Outburst Detector                          #                                                                       #
# DETECTOR utilities                                                   #
#                                                                      #
# Inés Pastor Marazuela (2019) - ines.pastor.marazuela@gmail.com       #
#                                                                      #
########################################################################
"""
Variability-related procedures specified into the documentation
"""

import time

# Built-in imports

from os.path import sys

# Third-party imports

from astropy.io import fits
from astropy import wcs
from astropy.table import Table


########################################################################


def extraction_photons(events_file):
    """
    Function extracting the E round list from its FITS events file.
    It alse returns the header information
    @param events_file: The events FITS file
    @return: The E round list
    @return: The events file header
    @return maximum, minimum
    @raise Exception: An exception from astropy if something went wrong
    """

    # Ouverture du fichier
    hdulist = fits.open(events_file)

    # Récupération des EVENTS
    events = hdulist[1].data
    header = hdulist[1].header
    events_filtered = []
    for i in range(12) :
        events_filtered.append([])


    for evt in events :
        events_filtered[int(evt['CCDNR'])-1].append(evt)


    hdulist.close()
    events_filtered_sorted = []
    for i in range(12) :
        events_filtered_sorted.append(sorted(events_filtered[i],
                key=lambda k: int(k['TIME'])))

    return events_filtered_sorted, header

########################################################################

# Deprecated
def extraction_info(events_file) :
    hdulist = fits.open(events_file)
    header  = hdulist[1].header
    data    = fits.getdata(events_file)
    dmin  = [min(data['X']), min(data['Y'])]
    dmax  = [max(data['X']), max(data['Y'])]
    return header, dmin, dmax

########################################################################


def extraction_deleted_periods(gti_file):
    """
    Function extracting the G round list from its Fits file.
    @param gti_file: The gti file
    @return: The G round list
    @raise Exception: An exception from astropy if something went wrong
    """
    hdulist = fits.open(gti_file)

    return hdulist[1].data

########################################################################

def fits_writer(data, sources, image, pars, file) :
    """
    Function writing the variability and sources to a fits file
    @param data: image of the variability data
    @param sources: detected variable sources
    @param image: image obtained with evselect, needed for the header
    @param pars: variability parameters used in the variability computation
    @param file: output file name
    """

    hdulist    = fits.open(image)
    header_img = hdulist[0].header
    hdulist.close()
    head_var_f = header_img
    head_var_f.append(card=('CREATOR', pars['CREATOR'], '[s] EXOD Time window'))
    head_var_f.append(card=('DATE', pars['DATE'], '[s] EXOD Time window'))
    head_var_f.append(card=('TW', pars['TW'], '[s] EXOD Time window'))
    head_var_f.append(card=('GTR', pars['GTR'], 'EXOD Good time ratio'))
    head_var_f.append(card=('DL', pars['DL'], 'EXOD Detection level'))
    head_var_f.append(card=('BS', pars['BS'], '[pix] EXOD Box size'))

    # data_var_f = Table(names=('VARIABILITY', 'RAWX', 'RAWY', 'CCDNR'),
    # dtype=('f8', 'i2', 'i2', 'i2'))

    # Creating fits file
    hdul_f   = fits.HDUList()
    hdul_var = fits.ImageHDU(data=data, header=head_var_f)
    hdul_src = fits.BinTableHDU(data=sources)
    hdul_f.append(hdul_var)
    hdul_f.append(hdul_src)

    # Writing to file
    hdul_f.writeto(file, overwrite=True)

    return True
