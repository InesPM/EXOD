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
FOLDER=/mnt/data/Ines/data
SCRIPTS=/mnt/data/Ines/EXOD/scripts

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

# Useful
########################################################################

var(){
  x=$1
  out=$(cat $SCRIPTS/file_names.py | grep ^$x | awk '{print $3}' | sed 's/"//g')
  echo $out
}

waitForFinish()
{
  STRING=$1;
  # wait until jobs have started
  sleep 1
  # check, twice, whether all done
  for i in 1 2 ; do
    job=99
    while [ $job -gt 0 ] ; do sleep 10; job=`top -b | head -n 40 | grep ${STRING} | wc -l`; done
  done
}


# Variability computation functions
################################################################################

filtering(){
  obs=$1

  if [ ! -f $FOLDER/$obs/PN_clean.fits ] ; then
    echo "bash $SCRIPTS/filtering.sh -f $FOLDER -o $obs" >> \
      $FOLDER/process_flt_${DL}_${TW}_${BS}_${GTR}
  fi
}

variability(){

  ###
  # Defining files and directories
  ###
  obs=$1

  # File names
  path=$FOLDER/$obs
  events_file=$path/$(var CLEAN_FILE)
  gti_file=$path/$(var GTI_FILE)
  img_file=$path/$(var IMG_FILE)

  ###
  # Variability computation
  ###
  # writing detector and renderer commands to file that will be run in parallel
  echo "python3 -W"ignore" $SCRIPTS/detector.py -path $path \
    -bs $BS -dl $DL -tw $TW -gtr $GTR -mta 1 --render" >> \
    $FOLDER/process_det_${DL}_${TW}_${BS}_${GTR}
  ((++count))
}

lightcurves(){

  # Output file
  echo "Observation Source DL TW P_chisq P_KS" >> \
    $FOLDER/sources_variability_${DL}_${TW}_${BS}_${GTR}

  for obs in ${observations[@]}; do
    var_src=$obs/${DL}_${TW}_${BS}_${GTR}/variable_sources.csv
    if [ -f $var_src ]; then
      if [ $(cat $var_src | wc -l) -gt 1 ]; then
        nsrc=$(($(cat $var_src | wc -l) - 1))
        echo "$obs $nsrc $DL $TW" >> \
          $FOLDER/variable_sources_${DL}_${TW}_${BS}_${GTR}
      fi
    fi
  done

  # Reading file
  let i=0
  while read line ; do
    line_data[$i]="${line}"
    ((++i))
  done < $FOLDER/variable_sources_${DL}_${TW}_${BS}_${GTR}

  # Starting loop
  for data in "${line_data[@]:1:i-1}"; do

    IFS=' ' read -r -a array <<< "$data"
    obs=${array[0]}
    src=${array[1]}
    dl=${array[2]}
    tw=${array[3]}
    echo "Number of sources $obs : $src"

    count=1
    while [[ $count -le $src ]] ; do
      echo "bash $SCRIPTS/lightcurve.sh -f $FOLDER -s $SCRIPTS -o $obs \
        -dl $DL -tw $TW -gtr $GTR -bs $BS -id $count" >> \
        $FOLDER/process_lc_${DL}_${TW}_${BS}_${GTR}
      ((++count))
    done
    #fi
  done
}

################################################################################
#                                                                              #
# Main programme                                                               #
#                                                                              #
################################################################################

# Retrieving a list of observations
cd $FOLDER
observations=(0*)
count=1

# Removing existing log files to avoid overwriting them
logs=(detected_sources_${DL}_${TW}_${BS}_${GTR} \
  process_flt_${DL}_${TW}_${BS}_${GTR} process_det_${DL}_${TW}_${BS}_${GTR} \
  process_ren_${DL}_${TW}_${BS}_${GTR} process_lc_${DL}_${TW}_${BS}_${GTR} \
  sources_variability_${DL}_${TW}_${BS}_${GTR} \
  variable_sources_${DL}_${TW}_${BS}_${GTR})
for l in ${logs[@]} ; do if [ -f $l ]; then rm $l; fi; done
#rm process_lc_${DL}_${TW}_${BS}_${GTR}

echo "Observation sources DL TW" >> variable_sources_${DL}_${TW}_${BS}_${GTR}

###
# Launching the variability computation
###

start=$(date)
time {
  nb_img=${#observations[@]}
  echo $nb_img

  Title "Writing commands to files"
  # Writing commands to files to run them in parallel
  for obs in "${observations[@]}"; do
    # Filtering observations
    filtering $obs
    # Variability computation
    variability $obs
  done

  # Filtering observations
  Title "Filtering observations"
  bash $SCRIPTS/parallel.sh $FOLDER/process_flt_${DL}_${TW}_${BS}_${GTR} $CPUS
  waitForFinish 'evselect'

  # Running detector
  Title "Applying detector"
  bash $SCRIPTS/parallel.sh $FOLDER/process_det_${DL}_${TW}_${BS}_${GTR} $CPUS
  waitForFinish 'python'

  Title "Creating big pdf with observations"
  files=$(ls */${DL}_${TW}_${BS}_${GTR}/sources.pdf)
  gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite \
    -sOutputFile=$FOLDER/variability_observations_${DL}_${TW}_${BS}_${GTR}.pdf \
    ${files[@]}

  # Generating lightcurves
  Title "Generating lightcurves"
  lightcurves
  bash $SCRIPTS/parallel.sh $FOLDER/process_lc_${DL}_${TW}_${BS}_${GTR} $CPUS
  sleep 10; waitForFinish 'evince'
  sleep 10; waitForFinish 'arfgen'
  sleep 10; waitForFinish 'epiclccorr'
  sleep 10; waitForFinish 'python'
  sleep 10

  title "Writing sources variability to a single file"

  for obs in ${observations[@]}; do
    if compgen -G "$obs/lcurve_${TW}/*xronos.log" > /dev/null; then
      var_files=$(ls $obs/lcurve_${TW}/*xronos.log)
      for var_src in ${var_files[@]}; do
        P_chisq=$(cat $var_src | grep "Chi-Square Prob of constancy" \
          | awk '{print $5}')
        P_KS=$(cat $var_src | grep "Kolm.-Smir. Prob of constancy" \
          |  awk '{print $5}')
        src=$(echo $var_src | sed 's/.*J/J/g' | sed 's/_.*//g')
        echo "$obs $src $DL $TW $P_chisq $P_KS" >> \
          sources_variability_${DL}_${TW}_${BS}_${GTR}
      done
    fi
  done

  Title "Creating big pdf with sources"
  files=$(ls */lcurve_${TW}/*.pdf)
  gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite \
    -sOutputFile=$FOLDER/lightcurves_${DL}_${TW}_${BS}_${GTR}.pdf ${files[@]}

  echo -e "\nTotal execution time for $nb_img obs. : "
}
end=$(date)
echo -e "Start : $start\nEnd  : $end"
