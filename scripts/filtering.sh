#!/bin/bash

########################################################################
#                                                                      #
# EXOD - EPIC-pn XMM-Newton Outburst Detector                          #
#                                                                      #
# Events file filtering                                                #
#                                                                      #
# Inés Pastor Marazuela (2019) - ines.pastor.marazuela@gmail.com       #
#                                                                      #
########################################################################

###
# Parsing arguments
###

# Default variables
RATE=0.5
FOLDER=/home/ines/data
SCRIPTS=/home/ines/EXOD/scripts
INST=PN

# Input variables
while [[ $# -gt 0 ]]; do
case "$1" in
  -o|-obs|--observation) OBS=${2}
  shift; shift ;;
  -r|--rate)             RATE=${2:-$RATE}
  shift; shift ;;
  -i|--instrument)       INST=${2:-$INST}
  shift; shift ;;
  -f|--folder)           FOLDER=${2:-$FOLDER}
  shift; shift ;;
  -s|--scripts)          SCRIPTS=${2:-$SCRIPTS}
  shift; shift ;;
esac
done

path=$FOLDER/$OBS
log_file=$path/PN_processing.log

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

# Useful
########################################################################

var(){
  x=$1
  out=$(cat $SCRIPTS/file_names.py | grep ^$x | awk '{print $3}' | sed 's/"//g')
  echo $out
}

print_log(){
  log=$1; message=$2
  echo -e "$message"
  echo >> $log "$message"
}

########################################################################
#                                                                      #
# Main programme                                                       #
#                                                                      #
########################################################################

Title "Filtering observation $OBS"

start=$(date)
echo > $log_file ""
print_log $log_file "Start: $start"
print_log $log_file "FOLDER      = ${FOLDER}"
print_log $log_file "OBSERVATION = ${OBS}"
print_log $log_file "INSTRUMENT  = ${INST}"
print_log $log_file "RATE        = ${RATE}"

###
# Preliminaries
###

title "Preliminaries"

# Setting up SAS
cd $path
export SAS_ODF=$path
export SAS_CCF=$path/ccf.cif
export HEADAS=$(var HEADAS)
. $HEADAS/headas-init.sh
. $(var SAS)

if [ ! -f $path/ccf.cif ]; then cifbuild; fi

###
# Filtering
###

title "Cleaning events file"

# File names
if [ "$INST" == "PN" ]; then
  org_file=$(ls $path/*${OBS}PN*PIEVLI*)
  clean_file=$path/$(var CLEAN_FILE)
  gti_file=$path/$(var GTI_FILE)
  img_file=$path/$(var IMG_FILE)
  rate_file=$path/$(var RATE_FILE)
  l=P

elif [ "$INST" == "M1" ]; then
  org_file=$(ls $path/*${OBS}M1*MIEVLI*)
  clean_file=$path/M1_clean.fits
  gti_file=$path/$(var GTI_FILE)
  img_file=$path/M1_image.fits
  rate_file=$path/$(var RATE_FILE)
  l=M

elif [ "$INST" == "M2" ]; then
  org_file=$(ls $path/*${OBS}M2*MIEVLI*)
  clean_file=$path/M2_clean.fits
  gti_file=$path/$(var GTI_FILE)
  img_file=$path/M2_image.fits
  rate_file=$path/$(var RATE_FILE)
  l=M

else
  echo "ERROR: Instrument not recognized"
  exit
fi

print_log $log_file "RAW FILE   = ${org_file}"
print_log $log_file "CLEAN FILE = ${clean_file}"
print_log $log_file "GTI FILE   = ${gti_file}"
print_log $log_file "IMAGE FILE = ${img_file}"
print_log $log_file "RATE FILE  = ${rate_file}"

# Creating GTI
if [ "$INST" == "PN" ]; then
  title "Creating GTI"

  evselect table=$org_file withrateset=Y rateset=$rate_file maketimecolumn=Y \
    timebinsize=100 makeratecolumn=Y \
    expression='#XMMEA_EP && (PI in [10000:12000]) && (PATTERN==0)' -V 0

  if [[ $RATE != [0-9]* ]]; then
    echo "Opening PN_rate.fits"
    fv $rate_file &
    read -p "Choose the GTI cut rate : " RATE
  fi
  echo "Creating Good Time Intervals with threshold RATE=$RATE"

  tabgtigen table=$rate_file expression="RATE<=$RATE" gtiset=$gti_file -V 0
fi

# Cleaning events file
evselect table=$org_file withfilteredset=Y filteredset=$clean_file destruct=Y \
  keepfilteroutput=T \
  expression="#XMMEA_E$l && gti($gti_file,TIME) && (PATTERN<=4) && (PI in [500:12000])" \
  -V 0

#ds9 $events_file -bin factor 64 -scale log -cmap bb -mode region &

# Creating image file
evselect table=$clean_file imagebinning=binSize imageset=$img_file \
  withimageset=yes xcolumn=X ycolumn=Y ximagebinsize=80 yimagebinsize=80 -V 0

end=$(date)
print_log $log_file "End: $end"
