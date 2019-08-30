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
import subprocess

# Third-party imports

from math import *
from astropy.table import Table
import numpy as np
import scipy.ndimage as nd

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

    def __init__(self, src):
        """
        Constructor for Source class. Computes the x and y attributes.
        @param src : source, output of variable_sources_position
            id_src:  The identifier number of the source
            ccd:     The CCD where the source was detected at
            rawx:    The x coordinate on the CCD
            rawy:    The y coordinate on the CCD
            r:       The raw pixel radius of the detected variable area
        """
        super(Source, self).__init__()

        self.id_src = src[0]
        self.ccd = src[1]
        self.rawx = src[2]
        self.rawy = src[3]
        self.rawr = src[4]
        self.x = None
        self.y = None
        self.skyr = self.rawr * 64
        self.ra = None
        self.dec = None
        self.r = self.skyr * 0.05 # arcseconds


    def sky_coord(self, path, out) :
        """
        Calculate sky coordinates with the sas task edet2sky.
        Return x, y, ra, dec
        """

        # Launching SAS commands, writing to output file
        in_file  = path + 'PN_image.fits'
        out_file = out + '/variable_sources.txt'

        if self.id_src == 0 : s = '>'
        else :      s = '>>'

        command = """
        export SAS_ODF=/home/ines/Documents/projects/EXOD/data/one_observation/0203542001;
        export SAS_CCF=/home/ines/Documents/projects/EXOD/data/one_observation/0203542001/ccf.cif;
        export HEADAS=/home/ines/astrosoft/heasoft-6.25/x86_64-pc-linux-gnu-libc2.27;
        . $HEADAS/headas-init.sh;
        . /home/ines/astrosoft/xmmsas_20180620_1732/setsas.sh;
        #echo "# Variable source -{0}-" {1} {6}/.txt;
        edet2sky datastyle=user inputunit=raw X={2} Y={3} ccd={4} calinfoset={5} -V 0 {1} {6}
        """.format(self.id_src, s, self.rawx, self.rawy, self.ccd, in_file, out_file)

        # Running command, writing to file
        process = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)
        time.sleep(0.5)

        # Reading
        det2sky = np.array([line.rstrip('\n') for line in open(out_file, 'r')])

        for i in range(len(det2sky)) :
            # Equatorial coordinates
            self.ra, self.dec = det2sky[np.where(det2sky == '# RA (deg)   DEC (deg)')[0][0] + 1].split()

            # Sky pixel coordinates
            self.x, self.y    = det2sky[np.where(det2sky == '# Sky X        Y pixel')[0][0] + 1].split()

        # Removing temporary output file
        os.remove(out_file)

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

def data_transformation(data, header) :
    """
    Performing geometrical transformations from raw coordinates to sky coordinates
    @param data: variability matrix
    @param header: header of the clean events file
    @return: transformed variability data
    """

    # Header information
    w = wcs.WCS(header)
    w.wcs.crpix = [header['REFXCRPX'], header['REFYCRPX']]
    w.wcs.cdelt = [header['REFXCDLT']/15, header['REFYCDLT']]
    w.wcs.crval = [header['REFXCRVL']/15, header['REFYCRVL']]
    w.wcs.ctype = [header['REFXCTYP'], header['REFYCTYP']]
    angle = header['PA_PNT']
    dlim = [header['REFXLMIN'], header['REFXLMAX'], header['REFYLMIN'], header['REFYLMAX']]

    # Transformations
    padX = (int((648-len(data))/2),)
    padY = (int((648-len(data[0]))/2),)
    dataR = np.flipud(nd.rotate(data, angle, reshape = True))
    dataT = np.pad(dataR, (padX, padY), 'constant', constant_values=0)

    return dataT
