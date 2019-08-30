#!/usr/bin/env python3
# coding=utf-8


########################################################################
#                                                                      #
# EXOD - EPIC-pn XMM-Newton Outburst Detector                          # #                                                                      #
# Various utilities for both detector and renderer                     #
#                                                                      #
# Inés Pastor Marazuela (2019) - ines.pastor.marazuela@gmail.com       #
#                                                                      #
########################################################################
"""
Various resources for both detector and renderer
"""

# Built-in imports

import sys
import os
import time
from functools import partial

# Third-party imports

from math import *
from astropy.table import Table

# Internal imports

import file_names as FileNames


########################################################################
#                                                                      #
# Function and procedures                                              #
#                                                                      #
########################################################################


def open_files(folder_name) :
    """
    Function opening files and writing their legend.
    @param  folder_name:  The directory to create the files.
    @return: log_file, info_file, variability_file, counter_per_tw,
    detected_variable_areas_file, time_windows_file
    """

    # Fixing the name of the folder
    if folder_name[-1] != "/" :
        folder_name += "/"

    # Creating the folder if needed
    if not os.path.exists(folder_name):
        try :
            os.makedirs(folder_name)
        except :
            print("Error in creating output directory.\nABORTING", file=sys.stderr)
            exit(-1)

    # Declaring the files
    log_f = None
    info_f = None
    var_f = None
    var_per_tw_f = None
    detected_var_areas_f = None
    tws_f = None

    # Creating the log file
    try:
        log_file = open(folder_name + FileNames.LOG, "w+")

    except IOError as e:
        print("Error in creating log.txt.\nABORTING", file=sys.stderr)
        print(e, file=sys.stderr)
        close_files(log_f, info_f, var_f, var_per_tw_f, detected_var_areas_f, tws_f)
        print_help()
        exit(-1)

    # Creating the file to store variability per pixel
    try :
        # variability_file = open(folder_name + FileNames.VARIABILITY, "w+")
        # variability_file.write("# Variability for each pixel.\n")
        # variability_file = Table(names=('VARIABILITY', 'RAWX', 'RAWY', 'CCDNR'))
        variability_file = folder_name + FileNames.VARIABILITY

    except IOError as e:
        print("Error in creating variability_file.csv.\nABORTING", file=sys.stderr)
        print(e, file=sys.stderr)
        close_files(log_f, info_f, var_f, var_per_tw_f, detected_var_areas_f, tws_f)
        print_help()
        exit(-1)


    # Creating the file to store variability per pixel per time window
    try :
        variability_per_tw_file = open(folder_name + FileNames.EVTS_PX_TW, "w+")
        variability_per_tw_file.write("# Events count for each time window for each pixel.\n")

    except IOError as e:
        print("Error in creating events_count_per_px_per_tw.csv.\nABORTING", file=sys.stderr)
        print(e, file=sys.stderr)
        close_files(log_f, info_f, var_f, var_per_tw_f, detected_var_areas_f, tws_f)
        print_help()
        exit(-1)

    # Creating the file to store areas detected as variable
    try :
        detected_variable_areas_file = open(folder_name + FileNames.VARIABLE_AREAS, "w+")
        detected_variable_areas_file.write("# SOURCE_NUMBER;CCD;RAWX;RAWY\n")

    except IOError as e:
        print("Error in creating detected_variable_areas_file.csv.\nABORTING", file=sys.stderr)
        print(e, file=sys.stderr)
        close_files(log_f, info_f, var_f, var_per_tw_f, detected_var_areas_f, tws_f)
        print_help()
        exit(-1)

    # Creating the file to store acceptable time windows
    try :
        time_windows_file = open(folder_name + FileNames.TIME_WINDOWS, "w+")
        time_windows_file.write("# NUMBER;STARTING_TIME\n")

    except IOError as e:
        print("Error in creating time_windows_file.csv.\nABORTING", file=sys.stderr)
        print(e, file=sys.stderr)
        close_files(log_f, info_f, var_f, var_per_tw_f, detected_var_areas_f, tws_f)
        print_help()
        exit(-1)

    # Creating the file to store detected variable areas
    try :
        detected_variable_sources = open(folder_name + FileNames.VARIABLE_SOURCES, "w+")
        detected_variable_sources.write("# SOURCE_NUMBER;CCD;RAWX;RAWY;RADIUS;RA;DEC\n")

    except IOError as e:
        print("Error in creating detected_variable_sources.csv.\nABORTING", file=sys.stderr)
        print(e, file=sys.stderr)
        close_files(log_f, info_f, var_f, var_per_tw_f, detected_var_areas_f, tws_f)
        print_help()
        exit(-1)

    return log_file, variability_file, variability_per_tw_file, detected_variable_areas_file, time_windows_file, detected_variable_sources


########################################################################


def close_files(log_f, var_f, var_per_tw_f, detected_var_areas_f, tws_f, detected_var_sources_f) :
    """
    Function closing all files.
    """

    if log_f :
        log_f.close()

    if var_per_tw_f :
        var_per_tw_f.close()

    if detected_var_areas_f :
        detected_var_areas_f.close()

    if tws_f :
        tws_f.close()

    if detected_var_sources_f:
        detected_var_sources_f.close()

########################################################################

class Tee(object):
    """
    Class object that will print the output to the log_f file
    and to terminal
    """
    def __init__(self, *files):
        self.files = files
    def write(self, obj):
        for f in self.files:
            f.write(obj)
            f.flush() # If you want the output to be visible immediately
    def flush(self) :
        for f in self.files:
            f.flush()

########################################################################


def read_from_file(file_path, counter=False, comment_token='#', separator=';') :
    """
    Function returning the content
    @param counter: True if it is the counters that are loaded, False if it is the variability
    @return: A list
    """
    data = []
    i = 0
    nb_lignes = -1

    with open(file_path) as f:
        for line in f:
            #if len(line) > 2 and comment_token not in line :
            if comment_token not in line :
                if i % 200 == 0 :
                    data.append([])
                    nb_lignes += 1
                if not counter :
                    data[nb_lignes].append(float(line))
                else :
                    data[nb_lignes].append([float(tok) for tok in line.split(separator)])
                i += 1

    return data


########################################################################


def read_tws_from_file(file_path, comment_token='#', separator=';') :
    """
    Reads the list of time windows from its file.
    @return: A list of couples (ID of the TW, t0 of the TW)
    """
    tws = []

    with open(file_path) as f :
        for line in f :
            if len(line) > 2 and comment_token not in line :
                line_toks = line.split(separator)
                tws.append((int(line_toks[0]), float(line_toks[1])))

    return tws


########################################################################


class Source(object):
    """
    Datastructure providing easy storage for detected sources.\n

    Attributes:\n
    id_src:  The identifier number of the source\n
    ccd:     The CCD where the source was detected at\n
    rawx:    The x coordinate on the CCD\n
    rawy:    The y coordinate on the CCD\n
    r:       The radius of the variable area\n
    x:       The x coordinate on the output image\n
    y:       The y coordinate on the output image
    """

    def __init__(self, id_src, ccd, rawx, rawy, r):
        """
        Constructor for Source class. Computes the x and y attributes.
        @param id_src:  The identifier number of the source
        @param ccd:     The CCD where the source was detected at
        @param rawx:    The x coordinate on the CCD
        @param rawy:    The y coordinate on the CCD
        """
        super(Source, self).__init__()

        self.id_src = id_src
        self.ccd = ccd
        self.rawx = rawx
        self.rawy = rawy
        self.r = r

        self.x = self.rawy
        self.y = self.rawx


        self.x = 400 - self.rawy
        self.y = 64 - self.rawx

        # "Left" side CCDs
        if self.ccd >= 7 :
            self.x = 400 - self.x
            self.y = 64 - self.y

        if self.ccd in (2, 11) :
            self.y += 64
        elif self.ccd in (1, 10) :
            self.y += 128
        elif self.ccd in (4, 7) :
            self.y += 192
        elif self.ccd in (5, 8) :
            self.y += 256
        elif self.ccd in (6, 9) :
            self.y += 320


########################################################################


def read_sources_from_file(file_path, comment_token='#', separator=';') :
    """
    Reads the source from their file
    @return: A list of Source objects
    """

    sources = []

    with open(file_path) as f :
        for line in f :
            if len(line) > 2 and comment_token not in line :
                toks = line.split(separator)
                sources.append(Source(toks[0], int(toks[1]), float(toks[2]), float(toks[3]), float(toks[4])))

    return sources

########################################################################

def ccd_config(data_matrix) :
    """
    Arranges the variability data
    """
    data_v = []

    # XMM-Newton EPIC-pn CCD arrangement
    ccds = [[8,7,6,9,10,11],[5,4,3,0,1,2]]

    # Building data matrix
    for c in ccds[0] :
        data_v.extend(np.flipud(data_matrix[c]))
    i = 0
    for c in ccds[1] :
        m = np.flip(data_matrix[c])
        for j in range(64) :
            data_v[i] = np.append(data_v[i], m[63 - j])
            i += 1
    return data_v

########################################################################
#                                                                      #
# Geometrical transformations                                          #
#                                                                      #
########################################################################

def size(x, y, alpha) :
    x0 = min(abs((x * cos (radians(alpha)) - y * sin(radians(alpha)))/(cos(2 * radians(alpha)))), abs((y * cos (radians(alpha) - pi/2) - x * sin(radians(alpha) - pi/2))/(cos(2 * radians(alpha) - pi))))
    y0 = min(abs((y * cos (radians(alpha)) - x * sin(radians(alpha)))/(cos(2 * radians(alpha)))), abs((x * cos (radians(alpha) - pi/2) - y * sin(radians(alpha) - pi/2))/(cos(2 * radians(alpha) - pi))))
    return x0, y0

def rotation(x0, y0, alpha) :
    x = x0 * cos(radians(alpha)) - y0 * sin(radians(alpha))
    y = x0 * sin(radians(alpha)) + y0 * cos(radians(alpha))
    return x, y

def transformation(px, py, dmax, dmin, angle) :
    """
    Function Performing geometrical transformations on the pixels of the observation
    @param px    : rawx
    @param py    : rawy
    @param dmax  : extracted from the events files header: max x and y
    @param dmin  : extracted from the events file header: min x and y
    @param angle : rotation angle. header['PA_PNT']
    @return      : transformed x, y positions
    """

    dsize = [(dmax[0] - dmin[0]), (dmax[1] - dmin[1])]
    center = [dsize[0]/2 + dmin[0], dsize[1]/2 + dmin[1]]
    x_max, y_max = size(dsize[0], dsize[1], angle)

    (x, y) = rotation((px - 200)*x_max/400, (py - 192)*y_max/384, angle)
    x = x + center[0]
    y = y + center[1]
    return x, y
