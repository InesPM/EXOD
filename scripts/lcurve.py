#!/usr/bin/env python3
# coding=utf-8

########################################################################
#                                                                      #
# EXOD - Searching for fast transients into XMM-Newton data            #
#                                                                      #
# Generating lightcurve plots                                          #
#                                                                      #
# Inés Pastor Marazuela (2019) - ines.pastor.marazuela@gmail.com       #
#                                                                      #
########################################################################


# Built-in imports

from math import *
from os import path
from os.path import sys

# Third-party imports

import argparse
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.ticker import FormatStrFormatter
from scipy.stats import binned_statistic
from astropy.io import fits

###
# Parsing arguments
###

parser = argparse.ArgumentParser()
parser.add_argument("-path", dest="path", help="Path to the observation files", nargs='?', type=str)
parser.add_argument("-name", dest="name", help="Source name", nargs='?', type=str)
parser.add_argument("-obs", help="Observation identifier", nargs='?', type=str, default="")
parser.add_argument("-src", help="Path to the source's lightcurve fits file", nargs='?', type=str, default=None)
parser.add_argument("-bgd", help="Path to the background's lightcurve fits file", nargs='?', type=str, default=None)
parser.add_argument("-gti", help="Path to the GTI of the observation", nargs='?', type=str, default=None)
parser.add_argument("-tw", help="Time window", nargs='?', type=int, default=100)
parser.add_argument("-n", help="Lightcurve number", nargs='?', type=str, default="")
parser.add_argument("-pcs", dest="pcs", help="Chi-square probability of constancy", nargs='?', type=float, default=None)
parser.add_argument("-pks", dest="pks", help="Kolmogorov-Smirnov probability of constancy", nargs='?', type=float, default=None)
parser.add_argument("-mode", dest="mode", help="Plot style: monochrome / medium / color", nargs='?', type=str, default="medium")
args = parser.parse_args()

###
# Defining variables
###

# Path
if args.path[-1] == '/' :
    args.path = args.path[:-1]

# Source and background files
if args.src == None :
    lccorr = '{0}/{1}/lcurve_{2}/{3}_lccorr_{2}.lc'.format(args.path, args.obs, args.tw, args.name)
    print(args.name)
    if path.exists(lccorr) :
        args.src = lccorr
    else : 
        args.src = '{0}/{1}/lcurve_{2}/{3}_lc_{2}_src.lc'.format(args.path, args.obs, args.tw, args.name)
        args.bgd = '{0}/{1}/lcurve_{2}/{3}_lc_{2}_bgd.lc'.format(args.path, args.obs, args.tw, args.src)
        if not path.exists(args.name) :
            print('ERROR: Source File {0} does not exist'.format(args.src))
            sys.exit()
        if not path.exists(args.name) :
            print('ERROR: Background File {0} does not exist'.format(args.bgd))
            sys.exit()
        
# GTI file
if args.gti == None :
    args.gti = '{0}/{1}/PN_gti.fits'.format(args.path, args.obs)
    if not path.exists(args.gti) :
        print('ERROR: File {0} does not exist'.format(args.gti))
        sys.exit()
        
# Output file
src = (args.name).replace("_", "+")
out='{0}/{1}/lcurve_{2}/{3}_lc_{2}.pdf'.format(args.path, args.obs, args.tw, args.name)
print(out)

###
# Extracting information from fits files
###

# Source
hdu_src = fits.open(args.src)
data_src    = hdu_src[1].data
head_src    = hdu_src[1].header
hdu_src.close()

cts  = data_src[:]['RATE']
time = data_src[:]['TIME']
std  = data_src[:]['ERROR']

tstart = head_src['TSTART']
tstop  = head_src['TSTOP']
t0     = time[0]
tf     = time[-1]
time   = time - tstart

# Background
if args.bgd != None :
    hdu_bgd = fits.open(args.bgd)
    data_bgd    = hdu_bgd[1].data
    head_bgd    = hdu_bgd[1].header
    hdu_bgd.close()

    cts_bgd  = data_bgd[:]['RATE']
    time_bgd = data_bgd[:]['TIME']
    std_bgd  = data_bgd[:]['ERROR']

    for t_s in range(len(time)) :
        for t_b in range(len(time_bgd)) :
            if time[t_s] == time_bgd[t_s] :
                cts[t_s] == cts[t_s] - cts_bgd[t_b]

# GTI
if args.gti != None :
    hdulist_gti = fits.open(args.gti)
    data_gti    = hdulist_gti[1].data
    hdulist_gti.close()

    start = data_gti[:]['START']
    stop  = data_gti[:]['STOP']

    #time_bgd  = time_bgd - tstart
    if len(data_gti) > 1 :
        start_gti = np.insert(stop, 0, tstart) - tstart
        stop_gti  = np.insert(start, len(data_gti), tstop) - tstart
    elif len(data_gti) == 0 :
        start_gti = np.array([tstart])
        stop_gti  = np.array([tstop])
    else :
        start_gti = np.array([])
        stop_gti  = np.array([])

    # GTI inclusion
    if len(start_gti) > 0 :
        for i in range(len(start_gti)) :
            cdt  = np.where((time > start_gti[i]) & (time < stop_gti[i]))
            cts  = np.delete(cts, cdt)
            time = np.delete(time, cdt)
            std  = np.delete(std, cdt)
else : data_gti=[]

# Binning data

# if dt > tw :
#     time_bin = np.linspace()
#     cts_bin, cts_edges, binnumber = binned_statistic(cts)

cdt  = np.where(np.isfinite(cts) == True)
time = time[cdt]
cts  = cts[cdt]
std  = std[cdt]

for i in range(len(cts)) :
    if cts[i] < 0 :
        cts[i] = 0

time = time[cts != 0]
std  = std[cts != 0]
cts  = cts[cts != 0]

# Max, min, etc

xmin = time[0]
xmax = time[-1]
ymin = 0
ymax = np.max(cts + std)

med = np.median(cts)
max = np.argmax(cts)    # Argument of the max
min = np.argmin(cts)    # Argument of the min
med = np.median(cts)

if cts[max] - med > med - cts[min] :
    var = max - med
    index   = max
    var_min = med/(ymax-ymin) - ymin/(ymax-ymin)
    var_max = cts[index]/(ymax-ymin) - ymin/(ymax-ymin)
    label = "Maximum"
else :
    var = med - min
    index   = min
    var_min = cts[index]/(ymax-ymin) - ymin/(ymax-ymin)
    var_max = med/(ymax-ymin) - ymin/(ymax-ymin)
    label   = "Minimum"
###
# Plotting
###

#seaborn-colorblind
color = ['#01b4bc', '#009E73', '#D55E00', '#fa5457', '#f6d51f', '#56B4E9']
rcParams['font.family'] = 'serif'

fig, ax = plt.subplots(figsize=(7,5))

# Source

if "mono" in args.mode :
    # Data
    plt.errorbar(time, cts, yerr=std, fmt='o', color='k', markersize=4, elinewidth=1.0, zorder=6)
    plt.plot(time, cts, "-", linewidth=0.5, color='k', label="Source",zorder=6)
    # Max/min, median
    plt.axhline(med, xmin=0, xmax=1, linestyle='--', color='k', linewidth=1, zorder=4, label="Median")
    plt.axhline(cts[index], xmin=0, xmax=1, linestyle=':', color='k', linewidth=1, zorder=4, label=label)
    # GTI
    if len(data_gti) > 1 :
        for i in range(len(data_gti)) :
            mpl.rcParams['hatch.linewidth'] = 0.1
            ax.axvspan(start_gti[i], stop_gti[i], hatch='\\\\\\', facecolor='none', edgecolor='k', zorder=1)
            ax.axvline(start_gti[i], color='w', lw=3, zorder=2)
            ax.axvline(stop_gti[i], color='w', lw=3, zorder=3)
            #ax.Axes_fill_betweenx([-1,100], start_gti[i], stop_gti[i], )
if "med" in args.mode :
    # Data
    plt.plot(time, cts, "o-", linewidth=0.7, markersize=2, color='k', label="Source",zorder=2)
    plt.fill_between(time, cts - std, cts + std, alpha=0.3, color='c', zorder=2)
    # Max/min, median
    plt.axhline(med, xmin=0, xmax=1, linestyle='--', color='#54008c', linewidth=1, zorder=4, label="Median")
    plt.axhline(cts[index], xmin=0, xmax=1, linestyle=':', color='#54008c', linewidth=1, zorder=4, label=label)
    # GTI
    if len(data_gti) > 1 :
        for i in range(len(data_gti)) :
            mpl.rcParams['hatch.linewidth'] = 0.1
            ax.axvspan(start_gti[i], stop_gti[i], facecolor= 'k', alpha=0.2, edgecolor='None', zorder=1)
elif "colo" in args.mode :
    # Data
    plt.plot(time, cts, 'o', markersize=2, color=color[0], zorder=6)
    plt.plot(time, cts, "o-", linewidth=0.7, markersize=2, color=color[0], label="Source",zorder=2)
    plt.fill_between(time, cts - std, cts + std, alpha=0.2, color=color[0])
    # Lines
    plt.axhline(med, xmin=0, xmax=1, linestyle='--', color=color[3], linewidth=1, zorder=3, label="Median")
    plt.axhline(cts[index], xmin=0, xmax=1, linestyle=':', color=color[3], linewidth=1, zorder=4, label=label)
    # GTI
    if len(data_gti) > 1 :
        for i in range(len(data_gti)) :
            ax.axvspan(start_gti[i], stop_gti[i], alpha=0.2, color=color[4])
# Labels
#plt.legend(loc='upper right', fontsize=10)
plt.xlabel("Time (s)", fontsize=16)
plt.ylabel("counts s$^{-1}$", fontsize=16)
# Text
plt.text(0.03, 0.90, args.n, transform = ax.transAxes, fontsize=16)
plt.text(0.1, 0.90, "OBS {0}".format(args.obs), transform = ax.transAxes, fontsize=16)
plt.text(0.1, 0.80, src, transform = ax.transAxes, fontsize=16)
# Probabilities of constancy
plt.text(0.95, 0.90, r"P($\chi^2$) = {0:.2e} ".format(args.pcs), horizontalalignment='right', transform = ax.transAxes, fontsize=16)
plt.text(0.95, 0.80, r"P(KS) = {0:.2e} ".format(args.pks), horizontalalignment='right', transform = ax.transAxes, fontsize=16)
# Setup
plt.xlim(xmin, xmax)
plt.ylim(ymin, ymax)
ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
ax.xaxis.set_major_locator(plt.MaxNLocator(5))
plt.minorticks_on()
ax.yaxis.set_ticks_position('both')
ax.xaxis.set_ticks_position('both')
plt.tick_params(axis='both', which='both', direction='in', labelsize=14)
plt.savefig(out, pad_inches=0, bbox_inches='tight')
#plt.show()
