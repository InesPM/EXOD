#!/usr/bin/env python3
# coding=utf-8

import sys
import os
import time

# Third-party imports
import argparse

# Internal imports
import file_names as FileNames
from renderer import *

###
# Parsing arguments
###
parser = argparse.ArgumentParser()

# Path to files
parser.add_argument("-path", help="Path to the folder containing the observation files", type=str)
parser.add_argument("-out", help="Name of the output file", default=FileNames.OUTPUT_IMAGE_ALL, type=str)

args = parser.parse_args()

# Modifying arguments
if args.path[-1] != '/' :
    args.path = args.path + '/'
args.out = args.path + args.out

# Dd=efining paths to files

file0 = args.path + '/5_3_3_1.0/'   + FileNames.VARIABILITY
file1 = args.path + '/6_10_3_1.0/'  + FileNames.VARIABILITY
file2 = args.path + '/7_30_3_1.0/'  + FileNames.VARIABILITY
file3 = args.path + '/8_100_3_1.0/' + FileNames.VARIABILITY

###
# Applying renderer
###

render_variability_all(file0, file1, file2, file3, args.out)

print(args.out)
