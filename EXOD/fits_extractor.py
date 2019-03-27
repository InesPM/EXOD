#!/usr/bin/env python3
# coding=utf-8


################################################################################
#                                                                              #
# EXOD - EPIC-pn XMM-Newton Outburst Detector                                  # 
#                                                                              #
# DETECTOR utilities                                                           #
#                                                                              #
# Inés Pastor Marazuela (2019) - ines.pastor.marazuela@gmail.com               #
#                                                                              #
################################################################################
"""
Variability-related procedures specified into the documentation
"""

# Built-in imports

from os.path import sys

# Third-party imports

from astropy.io import fits
from astropy import wcs


################################################################################


def extraction_photons(events_file):
    """
    Function extracting the E round list from its FITS events file.
    @param events_file: The events FITS file
    @return: The E round list
    @raise Exception: An exception from astropy if something went wrong
    """

    # Ouverture du fichier
    hdulist = fits.open(events_file)

    # Récupération des EVENTS
    events = hdulist[1].data
    events_filtres = []
    for i in range(12) :
        events_filtres.append([])


    for evt in events :
        events_filtres[int(evt['CCDNR'])-1].append(evt)


    hdulist.close()
    events_filtres_sorted = []
    for i in range(12) :
        events_filtres_sorted.append(sorted(events_filtres[i], key=lambda k: int(k['TIME'])))

    return events_filtres_sorted

################################################################################

def extraction_info(events_file) :
    hdulist = fits.open(events_file)
    header  = hdulist[0].header
    data    = fits.getdata(events_file)
    dmin  = [min(data['X']), min(data['Y'])]
    dmax  = [max(data['X']), max(data['Y'])]
    return header, dmin, dmax

################################################################################


def extraction_deleted_periods(gti_file):
    """
    Function extracting the G round list from its Fits file.
    @param gti_file: The gti file
    @return: The G round list
    @raise Exception: An exception from astropy if something went wrong
    """
    hdulist = fits.open(gti_file)

    return hdulist[1].data
