#!/bin/bash

################################################################################
#                                                                              #
# EXOD - EPIC-pn XMM-Newton Outburst Detector                                  #
#                                                                              #
# Events file filtering                                                        #
#                                                                              #
# Inés Pastor Marazuela (2019) - ines.pastor.marazuela@gmail.com               #
#                                                                              #
################################################################################

###
# Parsing arguments                                                            #
###

# Default variables
RATE=0.5
FOLDER=/mnt/data/Ines/DR5

# Input variables
while [[ $# -gt 0 ]]; do
case "$1" in
  -f|--folder)           FOLDER=${2:-$FOLDER}
  shift; shift ;;
  -o|-obs|--observation) OBS=${2}
  shift; shift ;;
  -r|--rate)             RATE=${2:-$RATE}
  shift; shift ;;
esac
done

echo -e "\tFOLDER      = ${FOLDER}"
echo -e "\tOBSERVATION = ${OBS}"
echo -e "\tRATE        = ${RATE}"


path=$FOLDER/$OBS

###
# Defining functions
###

Title(){
  message=$1; i=0; x='===='
  while [[ i -lt ${#message} ]]; do x='='$x; ((++i)); done
  echo -e "\n\t  $message \n\t$x"
}

title(){
  message=$1; i=0; x=----
  while [[ i -lt ${#message} ]]; do x=-$x; ((++i)); done
  echo -e "\n\t  $message \n\t$x"
}

################################################################################
#                                                                              #
# Main programme                                                               #
#                                                                              #
################################################################################

Title "Filtering observation $OBS"

###
# Preliminaries
###

title "Preliminaries"

# Setting up SAS
cd $path
export SAS_ODF=$path
export SAS_CCF=$path/ccf.cif
export HEADAS=/home/ines/astrosoft/heasoft-6.25/x86_64-pc-linux-gnu-libc2.27
. $HEADAS/headas-init.sh
. /home/ines/astrosoft/xmmsas_20180620_1732/setsas.sh

cifbuild

###
# Filtering
###

title "Cleaning events file"

# File names
org_file=$(ls $path/*$OBS*PIEVLI*)
events_file=$path/PN_clean.fits
gti_file=$path/PN_gti.fits
rate_file=$path/PN_rate.fits
img_file=$path/PN_image.fits

# Creating GTI
evselect table=$org_file withrateset=Y rateset=$rate_file maketimecolumn=Y timebinsize=100 makeratecolumn=Y expression='#XMMEA_EP && (PI in [10000:12000]) && (PATTERN==0)' -V 0

if [[ $RATE != [0-9]* ]]; then
  echo "Opening PN_rate.fits" 
  fv $rate_file &
  read -p "Choose the GTI cut rate : " RATE
fi
echo "Creating Good Time Intervals with threshold RATE=$RATE"

tabgtigen table=$rate_file expression="RATE<=$RATE" gtiset=$gti_file -V 0

# Cleaning events file
evselect table=$org_file withfilteredset=Y filteredset=$events_file destruct=Y keepfilteroutput=T expression="#XMMEA_EP && gti($gti_file,TIME) && (PATTERN<=4) && (PI in [500:12000])" -V 0

ds9 $events_file -bin factor 64 -scale log -cmap bb -mode region &

# Creating image file
evselect table=$events_file imagebinning=binSize imageset=$img_file withimageset=yes xcolumn=X ycolumn=Y ximagebinsize=80 yimagebinsize=80 -V 0

echo > $path/PN_processing.log "Rate: $RATE"
echo "The end" 
date 
