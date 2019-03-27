#!/usr/bin/env python3
# coding=utf-8

################################################################################
#                                                                              #
# EXOD - EPIC-pn XMM-Newton Outburst Detector                                  # #                                                                              #
# Generating lightcurve plots                                                  #
#                                                                              #
# Inés Pastor Marazuela (2019) - ines.pastor.marazuela@gmail.com               #
#                                                                              #
################################################################################


# Built-in imports

from math import *
from os.path import sys

# Third-party imports

import argparse
import numpy as np
import matplotlib as mpl
mpl.use("Pdf")
import matplotlib.pyplot as plt
from matplotlib import rcParams
from scipy.stats import binned_statistic
from astropy.io import fits
from astropy.time import Time

###
# Parsing arguments
###

parser = argparse.ArgumentParser()
parser.add_argument("-src", help="Path to the source's lightcurve fits file", nargs='?', type=str)
parser.add_argument("-bgd", help="Path to the background's lightcurve fits file", nargs='?', type=str)
parser.add_argument("-gti", help="Path to the GTI fits file", nargs='?', type=str)
parser.add_argument("-dtnb", help="Sampling time of the lightcurve", nargs='?', default=100, type=float)
parser.add_argument("-outdir", help="Path to the output directory", nargs='?', default=".", type=str)
parser.add_argument("-name", help="Name of the source", nargs='?', default="Source", type=str)
parser.add_argument("-text1", help="Text to add under the title", nargs='?', default="", type=str)
parser.add_argument("-obs", help="Observation number", nargs='?', default="", type=str)
parser.add_argument("-id", help="Id. of the variable source", nargs='?', default="", type=str)
args = parser.parse_args()

###
# Extracting information from fits files
###

hdulist_src = fits.open(args.src)
data_src    = hdulist_src[1].data
head_src    = hdulist_src[1].header
hdulist_src.close()

hdulist_bgd = fits.open(args.bgd)
data_bgd    = hdulist_bgd[1].data
hdulist_bgd.close()

hdulist_gti = fits.open(args.gti)
data_gti    = hdulist_gti[1].data
hdulist_gti.close()

cts_src  = data_src[:]['RATE']
time_src = data_src[:]['TIME']
std_src  = data_src[:]['ERROR']

cts_bgd  = data_bgd[:]['RATE']
time_bgd = data_bgd[:]['TIME']
std_bgd  = data_bgd[:]['ERROR']

start = data_gti[:]['START']
stop  = data_gti[:]['STOP']

tstart    = head_src['TSTART']
tstop     = head_src['TSTOP']
t0        = time_src[0]
tf        = time_src[-1]
time_src  = time_src - tstart
time_bgd  = time_bgd - tstart
if len(data_gti) > 1 :
    start_gti = np.insert(stop, 0, tstart) - tstart
    stop_gti  = np.insert(start, -1, tstop) - tstart
elif len(data_gti) == 0 :
    start_gti = np.array([tstart])
    stop_gti  = np.array([tstop])
else :
    start_gti = np.array([])
    stop_gti  = np.array([])

# GTI inclusion
if len(start_gti) > 0 :
    for i in range(len(start_gti)) :
        cdt_src  = np.where((time_src > start_gti[i]) & (time_src < stop_gti[i]))
        cts_src  = np.delete(cts_src, cdt_src)
        time_src = np.delete(time_src, cdt_src)
        std_src  = np.delete(std_src, cdt_src)

        cdt_bgd  = np.where((time_bgd > start_gti[i]) & (time_bgd < stop_gti[i]))
        cts_bgd  = np.delete(cts_bgd, cdt_bgd)
        time_bgd = np.delete(time_bgd, cdt_bgd)
        std_bgd  = np.delete(std_bgd, cdt_bgd)

# Removing null points
cdt_src  = np.where(cts_src == 0)
cts_src  = np.delete(cts_src, cdt_src)
time_src = np.delete(time_src, cdt_src)
std_src  = np.delete(std_src, cdt_src)

cdt_bgd  = np.where(cts_bgd == 0)
cts_bgd  = np.delete(cts_bgd, cdt_bgd)
time_bgd = np.delete(time_bgd, cdt_bgd)
std_bgd  = np.delete(std_bgd, cdt_bgd)

# Max, min, etc

xmin = time_src[0]
xmax = time_src[-1]
ymin = np.min((np.min(cts_bgd-std_bgd),np.min(cts_src-std_src)))
ymax = np.max((np.max(cts_bgd+std_bgd),np.max(cts_src+std_src)))

med_src = np.median(cts_src)
max_src = np.argmax(cts_src)    # Argument of the max
min_src = np.argmin(cts_src)    # Argument of the min
med_bgd = np.median(cts_src)

if cts_src[max_src] - med_src > med_src - cts_src[min_src] :
    var_src = max_src - med_src
    index   = max_src
    var_min = med_src/(ymax-ymin) - ymin/(ymax-ymin)
    var_max = cts_src[index]/(ymax-ymin) - ymin/(ymax-ymin)
    label = "Maximum"
else :
    var_src = med_src - min_src
    index   = min_src
    var_min = cts_src[index]/(ymax-ymin) - ymin/(ymax-ymin)
    var_max = med_src/(ymax-ymin) - ymin/(ymax-ymin)
    label   = "Minimum"
###
# Plotting
###

# Text
info1 = ''.join(["{0}\n".format(i) for i in args.text1.split(';')])
info2 = "Bin time    - {0} s\nStart time -  {1}\nStop time  - {2}".format(args.dtnb, Time(tstart, format='cxcsec').isot, Time(tstop, format='cxcsec').isot)
info3 = "Obs. {0} Source {1}".format(args.obs, args.id)

#seaborn-colorblind
color = ['#0072B2', '#009E73', '#D55E00', '#CC79A7', '#F0E442', '#56B4E9']
rcParams['font.family'] = 'serif'

fig = plt.figure(figsize=(5,21))
fig, (txt, ax) = plt.subplots(2,1, gridspec_kw = {'height_ratios':[1, 10]})
#ax = fig.add_subplot(111)
txt.axis("off")

# Source
ax.plot(time_src, cts_src, "o-", linewidth=0.7, markersize=2, color=color[0], label="Source",zorder=2)
ax.fill_between(time_src, cts_src - std_src, cts_src + std_src, alpha=0.3, color=color[0])
# Background
ax.plot(time_bgd, cts_bgd, "o-", linewidth=0.7, markersize=2, color=color[1], label="Background", zorder=1)
ax.fill_between(time_bgd, cts_bgd - std_bgd, cts_bgd + std_bgd, alpha=0.3, color=color[1])
# Lines
ax.axhline(med_src, xmin=0, xmax=1, linestyle='-', color=color[3], linewidth=1, zorder=3, label="Median")
ax.axhline(cts_src[index], xmin=0, xmax=1, linestyle='--', color=color[3], linewidth=1, zorder=4, label=label)
ax.plot(time_src[index], cts_src[index], 'x', color=color[3], zorder=5)
# GTI
if len(data_gti) > 1 :
    for i in range(len(data_gti)) :
        ax.axvspan(start_gti[i], stop_gti[i], alpha=0.1, color=color[4])
# Labels
txt.set_title(args.name.replace('_', '+'), fontsize=18)
txt.text(0.5, 0.9, info3, fontsize=10, verticalalignment='top', horizontalalignment='center')
txt.text(0.7, 0.5, info1, fontsize=10, verticalalignment='top', horizontalalignment='left')
txt.text(0.0, 0.5, info2, fontsize=10, verticalalignment='top', horizontalalignment='left')
ax.legend(loc='upper right', fontsize=8)
ax.set_xlabel("Time (s)", fontsize=14)
ax.set_ylabel("counts s$^{-1}$", fontsize=14)
ax.set_xlim(xmin, xmax)
ax.set_ylim(ymin, ymax)

#plt.show()
plt.minorticks_on()
plt.tick_params(axis='both', which='both', direction='in')
plt.savefig("{0}/{1}_lc_{2}.pdf".format(args.outdir, args.name, int(args.dtnb)))
