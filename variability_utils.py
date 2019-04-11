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
Implementation of variability-related procedures specified into the documentation
"""

# Third-party imports

import numpy as np
from numpy import inf

################################################################################
#                                                                              #
# Variability computation: procedure count_events                              #
#                                                                              #
################################################################################

def variability_computation(gti, time_interval, acceptable_ratio, start_time, end_time, data) :
	"""
	Function implementing the variability calculation using average technique.
	@param  gti:     G round, the list of TW cut-off the observation
	@param  time_interval:   The duration of a time window
	@param  acceptable_ratio:  The acceptability ratio for a TW - good time ratio
	@param  start_time:  The t0 instant of the observation
	@param  end_time: THe tf instant of the observation
	@param  data:    E round, the list of events sorted by their TIME attribute
	@return: The matrix V_round
	"""

	# Defining the variables and matrices
	n_bins = int(np.ceil((end_time - start_time )/time_interval))
	stop_time = start_time + n_bins*time_interval
	if (stop_time - end_time)/time_interval > acceptable_ratio :
		n_bins = n_bins - 1
		stop_time = start_time + n_bins * time_interval

	V_mat = np.ones([64,200])
	counted_events = np.zeros([64,200,n_bins])
	time_windows = np.arange(start_time, stop_time, time_interval)
	projection_ratio = np.ones(n_bins)

	# GTI
	cdt_start = []
	cdt_stop  = []
	for l in range(len(gti['START'])):
		start = np.where((time_windows[:] < gti[l]['START']) & (time_windows[:] + time_interval > gti[l]['START']))[0]
		if len(start) != 0 :
			cdt_start.append(start[0])
		stop  = np.where((gti[l]['STOP'] > time_windows[:]) & (gti[l]['STOP'] < time_windows[:] + time_interval))[0]
		if len(stop) != 0 :
			cdt_stop.append(stop[0])

	# Counting events
	i = 0			# Data counts
	n_last = None	# Last time window with a stop on it
	for n in range(n_bins) :
		# Good time
		t0 = 0
		tf = 0
		if n in cdt_start :
			t0 = gti[cdt_start.index(n)]['START']
			n_last = None
		else :
			t0 = time_windows[n]
		if n in cdt_stop :
			tf = gti[cdt_stop.index(n)]['STOP']
			n_last = n
		else :
			tf = time_windows[n] + time_interval
		if n_last == None :
			good_time = tf - t0
			projection_ratio[n] = good_time / time_interval
		else :
			projection_ratio[n] = 0

		# Counting events
		while i < len(data) and data[i]['TIME'] <= time_windows[n] + time_interval :
			j = int(data[i]['RAWX'])-1
			k = int(data[i]['RAWY'])-1
			for x in range(j-1,j+2) :
				for y in range(k-1,k+2) :
					if 0<=x<64 and 0<=y<200 :
						counted_events[x][y][n] += 1
			i += 1
		i += 1

	# Correcting with projection ratio
	cdt = np.where(projection_ratio >= acceptable_ratio)[0]
	counted_events = counted_events[:,:,cdt] / projection_ratio[cdt]
	time_windows   = time_windows[cdt]

	# Computing variability
	if len(counted_events[0][0]) > 1 :
		for i in range(len(counted_events)) :
			for j in range(len(counted_events[i])) :
				max = np.amax(counted_events[i][j])
				min = np.amin(counted_events[i][j])
				med = np.median(counted_events[i][j])
				if med != 0 :
					V_mat[i][j] = np.amax([(max - med),np.absolute(min - med)])/med
				else :
					V_mat[i][j] = max

	elif len(counted_events[0][0]) == 1 :
		print("No data within the GTI")
	return V_mat


################################################################################
#                                                                              #
# Detecting variable areas                                                     #
#                                                                              #
################################################################################


def box_computations(variability_matrix, x, y, box_size) :
    """
    Function summing the variability values into a box.
    @param variability_matrix:  The V round matrix
    @param x:   The x coordinate of the top-left corner of the box
    @param y:   The y coordinate of the top-left corner of the box
    @param box_size:   The length of a side of a box
    @return: The sum of the variability for each pixel of the box
    """
    assert x <= 63 - box_size
    assert y <= 199 - box_size
    # Exception raised if box out of the limits of the CCD

    cpt = 0

    for i in range(x, x + box_size) :
        for j in range(y, y + box_size) :
            cpt += variability_matrix[x][y]

    return cpt


################################################################################


def __add_to_detected_areas(x, y, box_size, detected_areas) :
    """
    Function summing the variability values into a box.
    @param x:   The x coordinate of the top-left corner of the box
    @param y:   The y coordinate of the top-left corner of the box
    @param box_size:   The length of a side of a box
    @param detected_areas:  The A round set, containing already detected areas
    @return: The sum of the variability for each pixel of the box ??? Should be the A round set, containing the updated detected areas
    """
    box_set = {(a, b) for a in range(x, x + box_size) for b in range(y, y + box_size)}

    inserted = False
    i = 0
    while (not inserted) and i < len(detected_areas) :
        # If a part of the box has already been detected, the two sets of coordinates are merged...
        if len(box_set & detected_areas[i]) > 1 :
            detected_areas[i] |= box_set
            #Equivalent to detected_areas[i] = detected_areas[i] | box_set => add area of the box to the detected areas
            # ...and the loop is over
            inserted = True
        i += 1

    # If there has not been any merger :
    if not inserted :
        detected_areas.append(box_set)

    return detected_areas




################################################################################

def variable_areas_detection(lower_limit, box_size, detection_level, variability_matrix) :
	"""
	Function detecting variable areas into a variability_matrix.
	@param lower_limit:         The lower_limit value is the smallest variability value needed to consider a pixel variable
	@param box_size:            The size of the box (optional, default = 3)
	@param detection_level:     A factor for the limit of detection
	@param variability_matrix:  The matrix returned by variability_calculation
	@return: A list of sets of coordinates for each area detected as variable
	"""

	output = []

	box_count = 0

	#detection_level = 4.56602579 * log10(TW) + 0.09141909

	for i in range(len(variability_matrix) - box_size) :
		j = 0
		m = 0	# boxes above detection level counter
		while j < len(variability_matrix[i]) - box_size :
			box_count = box_computations(variability_matrix, i, j, box_size)

			# If there's nothing into the box, it is completely skipped
			if box_count == 0 :
				j += box_size

			else :
				if box_count > detection_level * ((box_size**2) * lower_limit) :
					output=__add_to_detected_areas(i, j, box_size, output)
					m += 1

				j += 1

	return output
