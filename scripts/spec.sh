#!/bin/bash

################################################################################
#                                                                              #
# Spectral fit                                                                 #
#                                                                              #
################################################################################


# Default variables
C=2 ; M=tbabs*pow ; GTR=1.0 ; BS=3; ID=1
# Default folders
FOLDER=/mnt/data/Ines/data/M31
SCRIPTS=/mnt/data/Ines/progs

# Input variables
while [[ $# -gt 0 ]]; do
case "$1" in
  -o|-obs|--observation)  OBS=${2}
  shift; shift ;;
  # Variables
  -c|--counts)  C=${2:-$C}
  shift; shift ;;
  -m|--model)      M=${2:-$M}
  shift; shift ;;
  # Folders
  -f|--folder)            FOLDER=${2:-$FOLDER}
  shift; shift ;;
  -s|--scripts)           SCRIPTS=${2:-$SCRIPTS}
  shift; shift ;;
esac
done

path=$FOLDER/$OBS

# Style functions
################################################################################

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

input(){
  message="$1"; var="$2"
  read -p "$(tput setaf 6)$message $(tput sgr 0)" out 
  printf -v $var $out
}

################################################################################
#                                                                              #
# Preliminaries                                                                #
#                                                                              #
################################################################################

#read arguments
start=`date +%s`

title1 "Spectral fit Obs. $OBS"

echo -e "\tFOLDER  = ${FOLDER}"
echo -e "\tSCRIPTS = ${SCRIPTS}"
echo -e "\tCOUNTS  = ${C}"
echo -e "\tMODEL   = ${M}"

# Selecting the files and paths
clean_file=$path/PN_clean.fits
img_file=$path/PN_image.fits
path_out=$path/spec
if [ ! -d $path_out ]; then mkdir $path_out; fi
cd $path

# Setting SAS tools
export SAS_ODF=$path
export SAS_CCF=$path/ccf.cif
export HEADAS=/usr/local/heasoft-6.22.1/x86_64-unknown-linux-gnu-libc2.19/
. $HEADAS/headas-init.sh
. /usr/local/SAS/xmmsas_20170719_1539/setsas.sh

if [ ! -f $path/ccf.cif ]; then cifbuild; fi

###
# Extraction
###
ds9 $img_file -scale log -cmap bb -mode region &

input "Proceed? [y/n] " reply
if [[ $reply = [n,N]* ]] ; then exit; fi

echo -e "\n"
input "Source position     [X] " srcX
input "Source position     [Y] " srcY 
input "Background position [X] " bgdX
input "Background position [Y] " bgdY
input "Radius              [R] " R
echo -e "\n"

title2 "evselect"

evselect table=$clean_file withspectrumset=yes spectrumset=$path_out/PN_spectrum_src.fits  energycolumn=PI spectralbinsize=5 withspecranges=yes specchannelmin=0 specchannelmax=20479 expression="(FLAG==0) && (PATTERN<=4) && ((X,Y) in CIRCLE($srcX,$srcY,$R))"

evselect table=$clean_file withspectrumset=yes spectrumset=$path_out/PN_spectrum_bgd.fits  energycolumn=PI spectralbinsize=5 withspecranges=yes specchannelmin=0 specchannelmax=20479 expression="(FLAG==0) && (PATTERN<=4) && ((X,Y) in CIRCLE($bgdX,$bgdY,$R))"

title2 "backscale"

backscale spectrumset=$path_out/PN_spectrum_src.fits badpixlocation=$clean_file
backscale spectrumset=$path_out/PN_spectrum_bgd.fits badpixlocation=$clean_file

title2 "gen"

rmfgen spectrumset=$path_out/PN_spectrum_src.fits rmfset=$path_out/PN.rmf

arfgen spectrumset=$path_out/PN_spectrum_src.fits arfset=$path_out/PN.arf withrmfset=yes rmfset=$path_out/PN.rmf badpixlocation=$clean_file detmaptype=psf

specgroup spectrumset=$path_out/PN_spectrum_src.fits mincounts=$C oversample=3 rmfset=$path_out/PN.rmf arfset=$path_out/PN.arf backgndset=$path_out/PN_spectrum_bgd.fits groupedset=$path_out/PN_spectrum_grp.fits

###
# Spectral fitting
###

#cd $path_out
#xspec $SCRIPTS/spec.tcl > $path_out/spec.txt






