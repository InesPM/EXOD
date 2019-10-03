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

# bash variabilitectron.sh <FOLDER> <DL0> <DLf> <TW> <GTR> <BS0> <BSf>

# Defining the detection parameters
FOLDER=$1            	# FOLDER where obs are stored
DL0=$2                   # 10
DLf=$3
TW=$4                   # 100
GTR=$5                  # 1.0
BS0=$6                   # 5
BSf=$7
CPUS=10

# Folders where the different files necessary for the analysis are stored
SCRIPTS=/mnt/data/Ines/progs
FOLDER_CCF=/mnt/xmmcat/xmm-calib/ccf
FOLDER_EVTS=/mnt/xmmcat/3xmm_pievli
FOLDER_FBKT=/mnt/data/Ines/data/fbktsr_dr5

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
	i=0; x='----'
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
    sleep 100
  done
  # Erasing wait_file when waiting = False
  if [[ $(cat $wait_file) == False ]]; then rm $wait_file; fi
}


# Variability computation functions
################################################################################

variabilitectron(){

	###
	# Defining files and directories
	###
	# Reading the arguments
	obs=$1

	# File names
  path=$FOLDER/$obs
	events_file=$path/PN_clean.fits
	gti_file=$path/PN_gti.fits

	###
	# Variability calculation
	###
	# writing detector and renderer commands to file that will be run in parallel

	if [[ $count -lt $nb_img ]]; then
		echo "python3 -W"ignore" $SCRIPTS/relevel.py $events_file $path/${DL0}_${TW}_${GTR}_${BS0} $path/${DLf}_${TW}_${GTR}_${BSf} -obs $obs -bs $BSf -dl $DLf -tw $TW -gtr $GTR -mta 1 -ol $FOLDER/variable_sources_${DLf}_${TW}_${BSf}" >> $FOLDER/process_rel_${DLf}_${TW}_${BSf}
		echo "python3 -W"ignore" $SCRIPTS/renderer.py $path/${DLf}_${TW}_${GTR}_${BSf} $events_file -obs $obs -tw $TW -dl $DLf" >> $FOLDER/process_ren_${DLf}_${TW}_${BSf}

  else
    echo -n "python3 -W"ignore" $SCRIPTS/relevel.py $events_file $path/${DL0}_${TW}_${GTR}_${BS0} $path/${DLf}_${TW}_${GTR}_${BSf} -obs $obs -bs $BSf -dl $DLf -tw $TW -gtr $GTR -mta 1 -ol $FOLDER/variable_sources_${DLf}_${TW}_${BSf} ; " >> $FOLDER/process_rel_${DLf}_${TW}_${BSf}
    echo "echo 'False' > $FOLDER/waiting_rel" >> $FOLDER/process_rel_${DLf}_${TW}_${BSf}
    echo -n "python3 -W"ignore" $SCRIPTS/renderer.py $path/${DLf}_${TW}_${GTR}_${BSf} $events_file -obs $obs -tw $TW -dl $DLf ; " >> $FOLDER/process_ren_${DLf}_${TW}_${BSf}
    echo "echo 'False' > $FOLDER/waiting_ren" >> $FOLDER/process_ren_${DLf}_${TW}_${BSf}
	fi
	((++count))
}

lightcurves(){
  Title "Generating lightcurves"
  echo $FOLDER

  # Output file
  echo "Observation Source DL TW P_chisq P_KS" >> $FOLDER/sources_variability_${DLf}_${TW}_${BSf}

  # Reading file
  let i=0
  while IFS=$'\n' read -r line; do
  	line_data[i]="${line}"
  	((++i))
  done < $FOLDER/variable_sources_${DLf}_${TW}_${BSf}

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
        echo "bash $SCRIPTS/lightcurve.sh $FOLDER/$obs $SCRIPTS $count $DL $TW $FOLDER/sources_variability_${DLf}_${TW}_${BSf}" >> $FOLDER/process_lc_${DL}_${TW}_${BSf}
        ((++count))
      done
    #fi
  done
  echo "echo 'False' > $FOLDER/waiting_lc" >> $FOLDER/process_lc_${DL}_${TW}_${BS}
}

################################################################################
#                                                                              #
# Main programme                                                               #
#                                                                              #
################################################################################

# Retrieving a list of observations
cd $FOLDER
imgfull=(0*)
#files=(*PNS*PIEVLI*)
#nb_obs=${#files[@]}
count=1

# Removing existing log files to avoid overwriting them
logs=(detected_sources_${DLf}_${TW}_${BSf} process_rel_${DLf}_${TW}_${BSf} process_ren_${DLf}_${TW}_${BSf} process_lc_${DLf}_${TW}_${BSf} sources_variability_${DLf}_${TW}_${BSf} variable_sources_${DLf}_${TW}_${BSf})
for l in ${logs[@]} ; do if [ -f $l ]; then rm $l; fi; done

echo "Observation sources DL TW" >> variable_sources_${DLf}_${TW}_${BSf}
echo "echo 'True' > $FOLDER/waiting_rel" >> process_rel_${DLf}_${TW}_${BSf}
echo "echo 'True' > $FOLDER/waiting_ren" >> process_ren_${DLf}_${TW}_${BSf}
echo "echo 'True' > $FOLDER/waiting_lc" >> process_lc_${DLf}_${TW}_${BSf}

###
# Launching the variability computation
###

start=$(date)
time {
	nb_img=${#imgfull[@]}

  Title "Writing commands to files"
  # Writing commands to files to run them in parallel
	for obs in "${imgfull[@]}"; do
		# Variability computation
		variabilitectron $obs
	done

	# Running relevel
  Title "Applying detector"
  bash $SCRIPTS/parallel.sh $FOLDER/process_rel_${DLf}_${TW}_${BSf} $CPUS
  waitForFinish $FOLDER/waiting_rel

	# Running renderer
  Title "Applying renderer"
  bash $SCRIPTS/parallel.sh $FOLDER/process_ren_${DLf}_${TW}_${BSf} $CPUS
  waitForFinish $FOLDER/waiting_ren

  #Title "Creating big pdf with observations"
  #observations=($FOLDER/0*/${DLf}_${TW}_${GTR}_${BSf}/sources.pdf)
  #gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile=$FOLDER/variability_observations_${DLf}_${TW}_${BSf}.pdf ${observations[@]}

	# Generating lightcurves
  Title "Generating lightcurves"
  lightcurves
  #bash $SCRIPTS/parallel.sh $FOLDER/process_lc_${DLf}_${TW}_${BSf} $CPUS
  #waitForFinish $FOLDER/waiting_lc

  #Title "Creating big pdf with sources"
  #sources=($FOLDER/0*/lcurve_${TW}/J*lc_${TW}.pdf)
  #gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile=$FOLDER/variable_sources_${DLf}_${TW}_${BSf}.pdf ${sources[@]}

	echo -e "\nTotal execution time for $nb_img obs. : "
}
end=$(date)
echo -e "Start : $start\nEnd  : $end"
