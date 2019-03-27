#!/bin/bash

################################################################################
#                                                                              #
# EXOD - EPIC-pn XMM-Newton Outburst Detector                                  #
#                                                                              #
# Events file generation and variability computation                           #
#                                                                              #
# Inés Pastor Marazuela (2019) - ines.pastor.marazuela@gmail.com               #
#                                                                              #
################################################################################

# bash script generating the filtered events file and the GTI file and computing the variability.

# bash variabilitectron.sh <FOLDER> <DL> <TW> <GTR>

################################################################################
#                                                                              #
# Defining the functions                                                       #
#                                                                              #
################################################################################

# Style functions
################################################################################

Title(){
  message=$1
	i=0; x='===='
	while [[ i -lt ${#message} ]]; do x='='$x; ((++i)); done
  echo -e "\n\t  $message \n\t$x"
}

title(){
  message=$1
	i=0; x=----
	while [[ i -lt ${#message} ]]; do x=-$x; ((++i)); done
  echo -e "\n\t  $message \n\t$x"
}

# Data analysis functions
################################################################################

preliminaries(){

	title "Preliminaries"
	path=$1/$2
	obs=$2

	# Setting up SAS
	cd $path
	export SAS_ODF=$path
	export SAS_CCF=$path/ccf.cif
	export HEADAS=/usr/local/heasoft-6.22.1/x86_64-unknown-linux-gnu-libc2.19/
	. $HEADAS/headas-init.sh
	. /usr/local/SAS/xmmsas_20170719_1539/setsas.sh

	if [ ! -f $path/ccf.cif ]; then cifbuild; else echo "CIF already built"; fi
	#if [ ! -f $path/*SUM.SAS ]; then odfingest; else echo "ODF already ingested"; fi
	#if [ ! -f $(ls $path/*ImagingEvts*) ]; then epproc; else echo "EP already processed"; fi

}

filtering(){

	title "Cleaning events file"
	path=$1/$2
	obs=$2
  #rate=0.5   #$3

	# File names
	org_file=$(ls /mnt/xmmcat/3xmm_pievli/*$obs*PIEVLI*)
	events_file=$path/PN_clean.fits
	gti_file=$path/PN_gti.fits
	rate_file=$path/PN_rate.fits

	# Creating GTI
	evselect table=$org_file withrateset=Y rateset=$rate_file maketimecolumn=Y timebinsize=100 makeratecolumn=Y expression='#XMMEA_EP && (PI in [10000:12000]) && (PATTERN==0)' -V 0

	fv $rate_file &
	read -p "Rate : " rate
	echo "Creating Good Time Intervals with threshold RATE=$rate"

	tabgtigen table=$rate_file expression="RATE<=$rate" gtiset=$gti_file -V 0

	# Cleaning events file
	evselect table=$org_file withfilteredset=Y filteredset=$events_file destruct=Y keepfilteroutput=T expression="#XMMEA_EP && gti($gti_file,TIME) && (PATTERN<=4) && (PI in [500:12000])" -V 0

	echo > $path/PN_processing.log "Rate: $rate"
}

# Main programme

FOLDER=$1
obs=$2

# Filtering observations
preliminaries $FOLDER $obs
filtering $FOLDER $obs
