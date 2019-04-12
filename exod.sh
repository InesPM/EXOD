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

################################################################################
#                                                                              #
# Parsing arguments                                                            #
#                                                                              #
################################################################################

# Default variables
DL=8 ; TW=100 ; GTR=1.0 ; BS=3 ; CPUS=12
# Default folders
FOLDER=/home/ines/data
SCRIPTS=/home/ines/EXOD

# Input variables
while [[ $# -gt 0 ]]; do
case "$1" in
  # Variables
  -dl|--detection-level)  DL=${2:-$DL}
  shift; shift ;;
  -tw|--time-window)      TW=${2:-$TW}
  shift; shift ;;
  -gtr|--good-time-ratio) GTR=${2:-$GTR}
  shift; shift ;;
  -bs|--box-size)         BS=${2:-$BS}
  shift; shift ;;
  -cpus|--cpus)           CPUS=${2:-$CPUS}
  shift; shift ;;
  # Folders
  -f|--folder)            FOLDER=${2:-$FOLDER}
  shift; shift ;;
  -s|--scripts)           SCRIPTS=${2:-$SCRIPTS}
  shift; shift ;;
esac
done

echo -e "\tFOLDER          = ${FOLDER}"
echo -e "\tSCRIPTS         = ${SCRIPTS}\n"
echo -e "\tDETECTION LEVEL = ${DL}"
echo -e "\tTIME WINDOW     = ${TW}"
echo -e "\tGOOD TIME RATIO = ${GTR}" 
echo -e "\tBOX SIZE        = ${BS}"
echo -e "\tCPUS            = ${CPUS}"

################################################################################
#                                                                              #
# Defining functions                                                           #
#                                                                              #
################################################################################

# Style functions
################################################################################

Title(){
  message=$1; i=0; x='===='
  while [[ i -lt ${#message} ]]; do x='='$x; ((++i)); done
  echo -e "\n\t  $message \n\t$x"
}

title(){
  message=$1; i=0; x='----'
  while [[ i -lt ${#message} ]]; do x='-'$x; ((++i)); done
  echo -e "\n\t  $message \n\t$x"
}

waitForFinish()
{
  wait_file=$1
  # Waiting until file is created
  while [[ ! -f $wait_file ]]; do
    echo "waiting does not exist"
    sleep 10
  done
  # Waiting while waiting = True
  while [[ $(cat $wait_file) == True ]]; do
    echo "Waiting"
    sleep 10
  done
  # Erasing wait_file when waiting = False
  if [[ $(cat $wait_file) == False ]]; then rm $wait_file; fi
}


# Variability computation functions
################################################################################

filtering(){
  path_scripts=$1
  obs=$2

  if [ ! -f $path_obs/$obs/PN_clean.fits ] ; then
    if [[ $count -le $(( $nb_img )) ]]; then
      echo "bash $path_scripts/filtering.sh -f $FOLDER -o $obs" >> $path_obs/process_flt_${DL}_${TW}_${GTR}_${BS}
    else
      echo -n "bash $path_scripts/filtering.sh -f $FOLDER -o $obs; " >> $path_obs/process_flt_${DL}_${TW}_${GTR}_${BS}
      echo "echo 'False' > waiting_flt" >> $FOLDER/process_flt_${DL}_${TW}
    fi
  else
    if [[ $count -eq $(( $nb_img )) ]]; then echo "echo 'False' > waiting_flt" >> $FOLDER/process_flt_${DL}_${TW}_${GTR}_${BS}; fi
  fi
}

variabilitectron(){

  ###
  # Defining files and directories
  ###
  # Reading the arguments
  path_scripts=$1
  obs=$2

  # File names
  path=$FOLDER/$obs
  events_file=$path/PN_clean.fits
  gti_file=$path/PN_gti.fits

  ###
  # Variability computation
  ###
  # writing detector and renderer commands to file that will be run in parallel

  if [[ $count -lt $(( $nb_img )) ]]; then
    echo "python3 -W"ignore" $path_scripts/detector.py $events_file $gti_file $path/${DL}_${TW}_${GTR}_${BS} -obs $obs -bs $BS -dl $DL -tw $TW -gtr $GTR -mta 1 -ol $path_obs/variable_sources_${DL}_${TW}_${BS}" >> $path_obs/process_det_${DL}_${TW}_${GTR}_${BS}
    echo "python3 -W"ignore" $path_scripts/renderer.py $path/${DL}_${TW}_${GTR}_${BS} $events_file -obs $obs -tw $TW -dl $DL -bs $BS" >> $path_obs/process_ren_${DL}_${TW}_${GTR}_${BS}

  else
    echo -n "python3 -W"ignore" $path_scripts/detector.py $events_file $gti_file $path/${DL}_${TW}_${GTR}_${BS} -obs $obs -bs $BS -dl $DL -tw $TW -gtr $GTR -mta 1 -ol $path_obs/variable_sources_${DL}_${TW}_${BS} ; " >> $path_obs/process_det_${DL}_${TW}_${GTR}_${BS}
    echo "echo 'False' > $path_obs/waiting_det" >> $path_obs/process_det_${DL}_${TW}_${GTR}_${BS}
    echo -n "python3 -W"ignore" $path_scripts/renderer.py $path/${DL}_${TW}_${GTR}_${BS} $events_file -obs $obs -tw $TW -dl $DL -bs $BS ; " >> $path_obs/process_ren_${DL}_${TW}_${GTR}_${BS}
    echo "echo 'False' > $path_obs/waiting_ren" >> $path_obs/process_ren_${DL}_${TW}_${GTR}_${BS}
  fi
  ((++count))
}

lightcurves(){
  Title "Generating lightcurves"

  # Reading arguments
  path_obs=$1
  path_scripts=$2

  echo $path_obs

  # Output file
  echo "Observation Source DL TW P_chisq P_KS" >> $path_obs/sources_variability_${DL}_${TW}_${GTR}_${BS}

  # Reading file
  let i=0
  while IFS=$'\n' read -r line; do
  	line_data[i]="${line}"
  	((++i))
  done < $path_obs/variable_sources_${DL}_${TW}_${GTR}_${BS}

  # Starting loop
  for data in "${line_data[@]:1:i-1}"; do

    IFS=' ' read -r -a array <<< "$data"
    obs=${array[0]}
    src=${array[1]}
    dl=${array[2]}
    tw=${array[3]}

    count=1
    #if [[ $dl == $DL ]] && [[ $tw == $TW ]]; then
      while [[ $count -le $src ]] ; do
        echo "bash $path_scripts/lightcurve.sh $path_obs/$obs $path_scripts $count $DL $TW $path_obs/sources_variability_${DL}_${TW}" >> $path_obs/process_lc_${DL}_${TW}_${GTR}_${BS}
        ((++count))
      done
    #fi
  done
  echo "echo 'False' > $path_obs/waiting_lc" >> $path_obs/process_lc_${DL}_${TW}_${GTR}_${BS}
}

################################################################################
#                                                                              #
# Main programme                                                               #
#                                                                              #
################################################################################

# Retrieving a list of observations
cd $FOLDER
imgfull=(0*)
count=1

# Removing existing log files to avoid overwriting them
logs=(detected_sources_${DL}_${TW}_${GTR}_${BS} process_flt_${DL}_${TW}_${GTR}_${BS} process_det_${DL}_${TW}_${GTR}_${BS} process_ren_${DL}_${TW}_${GTR}_${BS} process_lc_${DL}_${TW}_${GTR}_${BS} sources_variability_${DL}_${TW}_${GTR}_${BS} variable_sources_${DL}_${TW}_${BS})
for l in ${logs[@]} ; do if [ -f $l ]; then rm $l; fi; done

echo "Observation sources DL TW" >> variable_sources_${DL}_${TW}_${BS}
echo "echo 'True' > $FOLDER/waiting_flt" >> process_flt_${DL}_${TW}_${GTR}_${BS}
echo "echo 'True' > $FOLDER/waiting_det" >> process_det_${DL}_${TW}_${GTR}_${BS}
echo "echo 'True' > $FOLDER/waiting_ren" >> process_ren_${DL}_${TW}_${GTR}_${BS}
echo "echo 'True' > $FOLDER/waiting_lc" >> process_lc_${DL}_${TW}_${GTR}_${BS}

###
# Launching the variability computation
###

start=$(date)
time {
  nb_img=${#imgfull[@]}
  echo $nb_img

  Title "Writing commands to files"
  # Writing commands to files to run them in parallel
  for obs in "${imgfull[@]}"; do
    # Filtering observations
    filtering $SCRIPTS $obs
    # Variability computation
    variabilitectron $SCRIPTS $obs
  done

  # Filtering observations
  Title "Filtering observations"
  bash $SCRIPTS/parallel.sh $FOLDER/process_flt_${DL}_${TW}_${GTR}_${BS} $CPUS
  waitForFinish $FOLDER/waiting_flt

  # Running detector
  Title "Applying detector"
  bash $SCRIPTS/parallel.sh $FOLDER/process_det_${DL}_${TW}_${GTR}_${BS} $CPUS
  waitForFinish $FOLDER/waiting_det

  # Running renderer
  Title "Applying renderer"
  bash $SCRIPTS/parallel.sh $FOLDER/process_ren_${DL}_${TW}_${GTR}_${BS} $CPUS
  waitForFinish $FOLDER/waiting_ren

  Title "Creating big pdf with observations"
  observations=($FOLDER/0*/${DL}_${TW}_${GTR}/sources.pdf)
  gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile=$FOLDER/variability_observations_${DL}_${TW}_${GTR}_${BS}.pdf ${observations[@]}

  # Generating lightcurves
  #Title "Generating lightcurves"
  #lightcurves $FOLDER $SCRIPTS $DL $TW $GTR $BS
  #bash $SCRIPTS/parallel.sh $FOLDER/process_lc_${DL}_${TW}_${GTR}_${BS} $CPUS
  #waitForFinish $FOLDER/waiting_lc

  #Title "Creating big pdf with sources"
  #sources=($FOLDER/0*/lcurve_${TW}/J*lc_${TW}.pdf)
  #gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile=$FOLDER/variable_sources_${DL}_${TW}_${GTR}_${BS}.pdf ${sources[@]}

  echo -e "\nTotal execution time for $nb_img obs. : "
}
end=$(date)
echo -e "Start : $start\nEnd  : $end"
