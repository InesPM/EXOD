#!/bin/bash

################################################################################
#                                                                              #
# EXOD - EPIC-pn XMM-Newton Outburst Detector                                  #
#                                                                              #
# Automatic lightcurve generation of the detected variable sources             #
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
DL=8 ; TW=100 ; GTR=1.0 ; BS=3; ID=1
# Default folders
FOLDER=/home/ines/data
SCRIPTS=/home/ines/EXOD

# Input variables
while [[ $# -gt 0 ]]; do
case "$1" in
  -o|-obs|--observation)  OBS=${2}
  shift; shift ;;
  # Variables
  -dl|--detection-level)  DL=${2:-$DL}
  shift; shift ;;
  -tw|--time-window)      TW=${2:-$TW}
  shift; shift ;;
  -gtr|--good-time-ratio) GTR=${2:-$GTR}
  shift; shift ;;
  -bs|--box-size)         BS=${2:-$BS}
  shift; shift ;;
  -id|--source-id)        ID=${2:-$ID}
  shift; shift ;;
  # Folders
  -f|--folder)            FOLDER=${2:-$FOLDER}
  shift; shift ;;
  -s|--scripts)           SCRIPTS=${2:-$SCRIPTS}
  shift; shift ;;
esac
done

output_log=$FOLDER/sources_variability_${DL}_${TW}_${GTR}_${BS}
path=$FOLDER/$OBS

# Style functions
########################################################################

title1(){
  message=$1; i=0; x='===='
  while [[ i -lt ${#message} ]]; do x='='$x; ((++i)); done
  echo -e "\n\t  $message \n\t$x\n"
}

title2(){
  message="$1 $OBS"; i=0; x='----'
  while [[ i -lt ${#message} ]]; do x='-'$x; ((++i)); done
  echo -e "\n\t  $message \n\t$x"
}

title3(){
  message="$1 $OBS"; echo -e "\n # $message"
}

# Useful
########################################################################

var(){
  x=$1
  out=$(cat scripts/file_names.py | grep ^$x | awk '{print $3}' | sed 's/"//g')
  echo $out
}

########################################################################
#                                                                      #
# Preliminaries                                                        #
#                                                                      #
########################################################################

#read arguments
start=`date +%s`

title1 "Lightcurve Obs. $OBS Src. $ID"

echo -e "\tFOLDER          = ${FOLDER}"
echo -e "\tSCRIPTS         = ${SCRIPTS}"
echo -e "\tOUTPUT LOG      = ${output_log}\n"
echo -e "\tDETECTION LEVEL = ${DL}"
echo -e "\tTIME WINDOW     = ${TW}"
echo -e "\tGOOD TIME RATIO = ${GTR}" 
echo -e "\tBOX SIZE        = ${BS}"

# Selecting the files and paths
clean_file=$path/$(var CLEAN_FILE)
gti_file=$path/$(var GTI_FILE)
img_file=$path/$(var IMG_FILE)
nosrc_file=$path/PN_sourceless.fits
sum_file=$(ls $path/*SUM.ASC)
fbk_file=$(ls $path/*$OBS*PNS*FBKTSR*)
path_out=$path/lcurve_${TW}
if [ ! -d $path_out ]; then mkdir $path_out; fi
cd $path

# Setting SAS tools
export SAS_ODF=$path
export SAS_CCF=$path/ccf.cif
export SAS_CCFPATH=/opt/ccfs
export HEADAS=$(var HEADAS)
. $HEADAS/headas-init.sh
. $(var SAS)

if [ ! -f $path/ccf.cif ]; then cifbuild; fi

################################################################################
#                                                                              #
# Source selection                                                             #
#                                                                              #
################################################################################

title2 "Preliminaries"

cd $path_out
if [ ! -f $path/${DL}_${TW}_${GTR}_${BS}/sources.pdf ]; then 
  python3 -W"ignore" $SCRIPTS/renderer.py $path/${DL}_${TW}_${GTR}_${BS} $clean_file -obs $OBS -tw $TW -dl $DL -bs $BS; fi
if [ ! -f $img_file ]; then
  evselect table=$clean_file imagebinning=binSize imageset=$img_file withimageset=yes xcolumn=X ycolumn=Y ximagebinsize=80 yimagebinsize=80 -V 0
fi

###
# Reading data from the detected_variable_sources file
###

data=$(cat $path/${DL}_${TW}_${GTR}_${BS}/detected_variable_sources.csv | grep "^${ID};")

###
# Defining source position
###

IFS=';' read -r -a array <<< "$data"
n=${array[0]}
ccd=${array[1]}
rawx=${array[2]}
rawy=${array[3]}
srcR=$((${array[4]} * 64))				# Pixel units

echo "ccd=$ccd, rawx=$rawx, rawy=$rawy, r=$srcR"

srcXY=$(ecoordconv imageset=$img_file coordtype=raw x=$rawx y=$rawy ccdno=$ccd | tee /dev/tty|grep 'X: Y:'|sed 's/X: Y: //'|sed 's/ /,/g'|sed 's/,//')

# Correcting source and background position
srcexp=$(eregionanalyse imageset=$img_file srcexp="(X,Y) in CIRCLE($srcXY,$srcR)" backval=0.1|tee /dev/tty|grep 'SASCIRCLE'|sed 's/SASCIRCLE: //g')
srcR=$(echo $srcexp | sed "s/(X,Y) in CIRCLE([0-9]*.[0-9]*,[0-9]*.[0-9]*,//" | sed "s/)//")
# arcsec
srcRas=$(bc <<< "scale=2; $srcR * 0.05")

###
# Source name and coordinates
###
srccoord=$(ecoordconv imageset=$img_file coordtype=raw x=$rawx y=$rawy ccdno=$ccd | tee /dev/tty|grep ' RA: DEC: ' | sed 's/ RA: DEC: //g')

# Right ascension
# Decimal degrees
RAd=$(echo $srccoord | awk '{print $1}')
# h
RA=$(bc <<< "scale=5; $RAd / 15")
if [[ $RA == .* ]]; then RA=0$RA; fi
h=$(echo $RA | sed "s/.[0-9]*$//g")
m=$(bc <<< "scale=5; ($RA - $h) * 60" | sed "s/.[0-9]*$//g")
s=$(bc <<< "scale=2; (($RA - $h) * 60 - $m) * 60" | sed "s/.[0-9]*$//")
if [ ${#h} == 1 ]; then h=0$h ; elif [ ${#h} == 0 ] ; then h=00 ; fi
if [ ${#m} == 1 ]; then m=0$m ; elif [ ${#m} == 0 ] ; then m=00 ; fi
if [ ${#s} == 1 ]; then s=0$s ; elif [ ${#s} == 0 ] ; then s=00; fi

# Declination
DEC=$(echo $srccoord | awk '{print $2}')
if [[ $DEC == -* ]]; then DEC=${DEC#"-"}; link="-"; else link="_"; fi
if [[ $DEC == .* ]]; then DEC=0$DEC; fi
dg=$(echo $DEC | sed "s/.[0-9]*$//g")
am=$(bc <<< "scale=5; ($DEC - $dg) * 60" | sed "s/.[0-9]*$//g")
as=$(bc <<< "scale=2; (($DEC - $dg) * 60 - $am) * 60" | sed "s/.[0-9]*$//")
if [ ${#dg} == 1 ]; then dg=0$dg ; elif [ ${#dg} == 0 ] ; then dg=00 ; fi
if [ ${#am} == 1 ]; then am=0$am ; elif [ ${#am} == 0 ] ; then am=00 ; fi
if [ ${#as} == 1 ]; then as=0$as ; elif [ ${#as} == 0 ] ; then as=00; fi
if [[ $link == "-" ]]; then DEC=$link$DEC; fi

sleep 1
# Source name
src=J$h$m$s$link$dg$am$as
echo -e "\n\t$src\n"

###
# Background extraction region
###

if [ ! -f $nosrc_file ]; then
evselect table=$clean_file withfilteredset=Y filteredset=$nosrc_file destruct=Y keepfilteroutput=T expression="region($fbk_file:REGION,X,Y)" -V 0
fi

bgdXY=$(ebkgreg withsrclist=no withcoords=yes imageset=$img_file x=$RAd y=$DEC r=$srcRas coordtype=EQPOS | grep 'X,Y Sky Coord.' | head -1 | awk '{print$5$6}')
bgdexp="(X,Y) in CIRCLE($bgdXY,$srcR)"

sleep 1
echo -e "\nExtracting data obs. $OBS source $src with the following coordinates: \n  Source     : $srcexp\n  Background : $bgdexp"


################################################################################
#                                                                              #
# Lightcurve generation                                                        #
#                                                                              #
################################################################################

title2 "Extracting lightcurve"
if [ ! -f $path/PN_gti.wi ]; then gti2xronwin -i $path/PN_gti.fits -o $path/PN_gti.wi; fi

title3 "evselect 0.0734 s"
evselect table=$clean_file energycolumn=PI expression="$srcexp" withrateset=yes rateset=$path_out/${src}_lc_0.0734_src.lc timebinsize=0.0734 maketimecolumn=yes makeratecolumn=yes -V 0
evselect table=$nosrc_file energycolumn=PI expression="$bgdexp" withrateset=yes rateset=$path_out/${src}_lc_0.0734_bgd.lc timebinsize=0.0734 maketimecolumn=yes makeratecolumn=yes -V 0

title3 "evselect $TW s"
evselect table=$clean_file energycolumn=PI expression="$srcexp" withrateset=yes rateset=$path_out/${src}_lc_${TW}_src.lc timebinsize=$TW maketimecolumn=yes makeratecolumn=yes -V 0
evselect table=$nosrc_file energycolumn=PI expression="$bgdexp" withrateset=yes rateset=$path_out/${src}_lc_${TW}_bgd.lc timebinsize=$TW maketimecolumn=yes makeratecolumn=yes -V 0

title3 "epiclccorr"
epiclccorr srctslist=$path_out/${src}_lc_${TW}_src.lc eventlist=$clean_file outset=$path_out/${src}_lccorr_${TW}.lc bkgtslist=$path_out/${src}_lc_${TW}_bgd.lc withbkgset=yes applyabsolutecorrections=yes -V 0
sleep 1

title3 "lcstats"
lcstats cfile1="$path_out/${src}_lccorr_${TW}.lc" window=$path/PN_gti.wi dtnb=${TW} nbint=1000000 tchat=2
P=$(lcstats cfile1="$path_out/${src}_lccorr_${TW}.lc" window=$path/PN_gti.wi dtnb=${TW} nbint=1000000 tchat=2 | grep "Prob of constancy")
sleep 1

P_chisq=$(echo $P | sed "s/Chi-Square Prob of constancy. //" | sed "s/ (0 means.*//")
P_KS=$(echo $P | sed "s/.*Kolm.-Smir. Prob of constancy //" |  sed "s/ (0 means.*//")

echo -e "Probabilities of constancy : \n\tP_chisq = $P_chisq\n\tP_KS    = $P_KS"

title3 "lcurve"
python3 $SCRIPTS/lcurve.py -src $path_out/${src}_lc_${TW}_src.lc -bgd $path_out/${src}_lc_${TW}_bgd.lc -gti $path/PN_gti.fits -dtnb $TW -outdir $path_out -name $src -Pcs "$P_chisq" -PKS "$P_KS" -obs $OBS -id $ID

echo -e " # Source $path_out/${src}_lc_${TW}.pdf"

end=`date +%s`
runtime=$((end-start))

###
# Writing output to files, ending script
###

echo -e > $path_out/${src}_region.txt "Source     = $srcexp\nBackground = $bgdexp\nTotal time = $runtime"
echo -e >> $output_log "$OBS $ID $src $DL $TW $P_chisq $P_KS"

echo -e " # Total time obs. $OBS : $runtime seconds"
echo -e "\nObservation $OBS ended\nTotal time = $runtime seconds"
