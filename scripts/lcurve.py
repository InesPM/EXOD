#!/usr/bin/env python3
# coding=utf-8

################################################################################
#                                                                              #
# Variabilitectron - Searching for fast transients into XMM-Newton data        #
#                                                                              #
# Generating lightcurve plots                                                  #
#                                                                              #
# Inés Pastor Marazuela (2018) - ines.pastor.marazuela@gmail.com               #
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
from matplotlib.ticker import FormatStrFormatter
from scipy.stats import binned_statistic
from astropy.io import fits
from astropy.time import Time

###
# Parsing arguments
###

parser = argparse.ArgumentParser()
parser.add_argument("-src", help="Path to the source's lightcurve fits file", nargs='?', type=str)
parser.add_argument("-bgd", help="Path to the source's lightcurve fits file", nargs='?', type=str)
parser.add_argument("-gti", help="Path to the GTI fits file", nargs='?', type=str)
parser.add_argument("-dtnb", help="Sampling time of the lightcurve", nargs='?', default=100, type=float)
parser.add_argument("-outdir", help="Path to the output directory", nargs='?', default=".", type=str)
parser.add_argument("-name", help="Name of the source", nargs='?', default="Source", type=str)
parser.add_argument("-Pcs", help="Chi square probability of constancy", nargs='?', default="", type=str)
parser.add_argument("-PKS", help="Kolmogorov-Smirnov probability of constancy", nargs='?', default="", type=str)
parser.add_argument("-obs", help="Observation number", nargs='?', default="", type=str)
parser.add_argument("-id", help="Id. of the variable source", nargs='?', default="", type=str)
args = parser.parse_args()

# Defining variables

input=(args.src).split("/")
src=input[-1][:14].replace("_", "+")
out=(input[-1][:14] + "_lc_{0}.pdf".format(int(args.dtnb)))

###
# Extracting information from fits files
###

hdu_src = fits.open(args.src)
data_src    = hdu_src[1].data
head_src    = hdu_src[1].header
hdu_src.close()

hdu_bgd = fits.open(args.bgd)
data_bgd    = hdu_bgd[1].data
head_bgd    = hdu_bgd[1].header
hdu_bgd.close()

hdulist_gti = fits.open(args.gti)
data_gti    = hdulist_gti[1].data
hdulist_gti.close()

cts  = data_src[:]['RATE']
time = data_src[:]['TIME']
std  = data_src[:]['ERROR']

cts_bgd  = data_bgd[:]['RATE']
time_bgd = data_bgd[:]['TIME']
std_bgd  = data_bgd[:]['ERROR']

for t_s in range(len(time)) :
    for t_b in range(len(time_bgd)) :
        if time[t_s] == time_bgd[t_s] :
            cts[t_s] == cts[t_s] - cts_bgd[t_b]

tstart    = head_src['TSTART']
tstop     = head_src['TSTOP']
#t0 = time[0]
time = time - tstart

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

start = data_gti[:]['START']
stop  = data_gti[:]['STOP']
if len(data_gti) > 1 :
    start_gti = np.insert(stop, 0, tstart) - tstart
    stop_gti  = np.insert(start, -1, tstop) - tstart
elif len(data_gti) == 0 :
    start_gti = np.array([tstart])
    stop_gti  = np.array([tstop])
else :
    start_gti = np.array([])
    stop_gti  = np.array([])

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
color = ['#0072B2', '#009E73', '#D55E00', '#CC79A7', '#F0E442', '#56B4E9']
rcParams['font.family'] = 'serif'

fig, ax = plt.subplots(figsize=(7,5))

# Source
plt.plot(time, cts, "o-", linewidth=0.7, markersize=2, color=color[0], label="Source",zorder=2)
plt.fill_between(time, cts - std, cts + std, alpha=0.2, color=color[0])
# Lines
plt.axhline(med, xmin=0, xmax=1, linestyle='-', color=color[3], linewidth=1, zorder=3, label="Median")
plt.axhline(cts[index], xmin=0, xmax=1, linestyle=':', color=color[3], linewidth=1, zorder=4, label=label)
# GTI
if len(data_gti) > 1 :
    for i in range(len(data_gti)) :
        ax.axvspan(start_gti[i], stop_gti[i], alpha=0.1, color=color[4])
# Labels
#plt.legend(loc='upper right', fontsize=10)
plt.xlabel("time (s)", fontsize=14)
plt.ylabel("counts s$^{-1}$", fontsize=14)
# Text
plt.text(0.02, 0.93, args.obs, transform = ax.transAxes, ha='left', fontsize=14)
plt.text(0.34, 0.93, args.id, transform = ax.transAxes, ha='right', fontsize=14)
plt.text(0.18, 0.86, src, transform = ax.transAxes, ha='center', fontsize=14)
plt.text(0.72, 0.93, "P$_{\chi^2} = $" + args.Pcs, transform = ax.transAxes, ha='left', fontsize=14)
plt.text(0.72, 0.86, "P$_{KS} = $" + args.PKS, transform = ax.transAxes, ha='left', fontsize=14)

plt.xlim(xmin, xmax)
plt.ylim(ymin, ymax)

#plt.show()
ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
ax.xaxis.set_major_locator(plt.MaxNLocator(5))
plt.minorticks_on()
ax.yaxis.set_ticks_position('both')
ax.xaxis.set_ticks_position('both')
plt.tick_params(axis='both', which='both', direction='in', labelsize=14)
plt.savefig(out, pad_inches=0, bbox_inches='tight')
