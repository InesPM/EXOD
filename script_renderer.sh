#!/bin/bash

################################################################################
#                                                                              #
# Variabilitectron - Searching variability into XMM-Newton                     #
#                                                                              #
# Inés Pastor Marazuela (2018) - ines.pastor.marazuela@gmail.com               #
#                                                                              #
################################################################################

# bash /mnt/data/Ines/progs/script_renderer.sh

# Defining variables

FOLDER=/mnt/data/Ines/data/M31
SCRIPTS=/mnt/data/Ines/progs
CPUS=30

################################################################################
# Defining functions
function waitForFinish()
{
  local STRING;

  STRING=$3;

  # wait until jobs have started
  sleep 20

  # check, twice, whether all done
  for i in 1 2 ; do
    job=99
    while [ $job -gt 0 ] ; do sleep 10; top -b | head -n 40; job=`ps ax | grep ${STRING} | wc -l`; done
  done
  sleep 20
}

################################################################################
# Main programme

cd $FOLDER
observations=(0*)

# Removing existing log files to avoid overwriting them
logs=(process_det process_ren_all waiting_det waiting_ren_all)
for l in ${logs[@]} ; do if [ -f $l ]; then rm $l; fi; done

# Writing waiting files
echo "echo 'True' > $FOLDER/waiting_det" >> process_det
echo "echo 'True' > $FOLDER/waiting_ren_all" >> process_ren_all

# Loop on observation
for obs in ${observations[@]}; do

	# Missing parameters
	if [ ! -d $FOLDER/$obs/5_3_1.0_3 ] || [ $(cat $FOLDER/$obs/5_3_1.0_3/variability_file.csv | wc -l) -ne 153601 ]; then
		echo "python3 -W"ignore" $SCRIPTS/detector.py $FOLDER/$obs/PN_clean.fits $FOLDER/$obs/PN_gti.fits $FOLDER/$obs/5_3_1.0_3 -obs $obs -bs 3 -dl 5 -tw 3 -gtr 1.0 -mta 1 -ol $FOLDER/variable_sources_5_3_3; python3 -W"ignore" $SCRIPTS/renderer.py $FOLDER/$obs/5_3_1.0_3 $FOLDER/$obs/PN_clean.fits -obs $obs -tw 3 -dl 5" >> $FOLDER/process_det
	fi
	if [ ! -d $FOLDER/$obs/6_10_1.0_3 ] || [ $(cat $FOLDER/$obs/6_10_1.0_3/variability_file.csv | wc -l) -ne 153601 ]; then
		echo "python3 -W"ignore" $SCRIPTS/detector.py $FOLDER/$obs/PN_clean.fits $FOLDER/$obs/PN_gti.fits $FOLDER/$obs/6_10_1.0_3 -obs $obs -bs 3 -dl 6 -tw 10 -gtr 1.0 -mta 1 -ol $FOLDER/variable_sources_6_10_3; python3 -W"ignore" $SCRIPTS/renderer.py $FOLDER/$obs/6_10_1.0_3 $FOLDER/$obs/PN_clean.fits -obs $obs -tw 10 -dl 6" >> $FOLDER/process_det
	fi
	if [ ! -d $FOLDER/$obs/7_30_1.0_3 ] || [ $(cat $FOLDER/$obs/7_30_1.0_3/variability_file.csv | wc -l) -ne 153601 ]; then
		echo "python3 -W"ignore" $SCRIPTS/detector.py $FOLDER/$obs/PN_clean.fits $FOLDER/$obs/PN_gti.fits $FOLDER/$obs/7_30_1.0_3 -obs $obs -bs 3 -dl 7 -tw 30 -gtr 1.0 -mta 1 -ol $FOLDER/variable_sources_7_30_3; python3 -W"ignore" $SCRIPTS/renderer.py $FOLDER/$obs/7_30_1.0_3 $FOLDER/$obs/PN_clean.fits -obs $obs -tw 30 -dl 7" >> $FOLDER/process_det
	fi
	if [ ! -d $FOLDER/$obs/8_100_1.0_3 ] || [ $(cat $FOLDER/$obs/8_100_1.0_3/variability_file.csv | wc -l) -ne 153601 ]; then
		echo "python3 -W"ignore" $SCRIPTS/detector.py $FOLDER/$obs/PN_clean.fits $FOLDER/$obs/PN_gti.fits $FOLDER/$obs/8_100_1.0_3 -obs $obs -bs 3 -dl 8 -tw 100 -gtr 1.0 -mta 1 -ol $FOLDER/variable_sources_8_100_3; python3 -W"ignore" $SCRIPTS/renderer.py $FOLDER/$obs/8_100_1.0_3 $FOLDER/$obs/PN_clean.fits -obs $obs -tw 100 -dl 8" >> $FOLDER/process_det
	fi

	# Rendering all
	if [ ! -f $FOLDER/$obs/renderer/sources_all.pdf ]; then 
		echo "python3 -W"ignore" $SCRIPTS/renderer_all.py -in $FOLDER/$obs/5_3_1.0_3/ $FOLDER/$obs/6_10_1.0_3/ $FOLDER/$obs/7_30_1.0_3/ $FOLDER/$obs/8_100_1.0_3/ -out $FOLDER/$obs/renderer/ -evts $FOLDER/$obs/PN_clean.fits -obs $obs" >> $FOLDER/process_ren_all
	fi
	
done

# Finishing writing files
echo "echo 'False' > $FOLDER/waiting_det" >> process_det
echo "echo 'False' > $FOLDER/waiting_ren_all" >> process_ren_all

# Parallel

#bash $SCRIPTS/parallel.sh $FOLDER/process_det $CPUS
#waitForFinish '[p]'ython

bash $SCRIPTS/parallel.sh $FOLDER/process_ren_all $CPUS
#waitForFinish '[p]'ython

# Writing big pdf
#nb=(0 500 1000 1500 2000 2500 3000 3500 4000 4500 5000 5500)
#files=()

#for obs in ${observations[@]}; do if [ -f $FOLDER/$obs/renderer/sources_all.pdf ]; then files+=($FOLDER/$obs/renderer/sources_all.pdf); fi; done

#for n in ${nb[@]}; do gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile=$FOLDER/renderer/variability_observations_${n}-$(($n + 499)).pdf ${files[@]:$n:500}; done
