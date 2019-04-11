#!/bin/bash

folder=$1
obs=$2
path=$folder/$obs
if [ ! -d $path ]; then mkdir $path; fi

cd $path

# Observation
curl -o $path/P${obs}PNS001PIEVLI.FTZ "http://nxsa.esac.esa.int/nxsa-sl/servlet/data-action-aio?obsno=${obs}&instname=PN&level=PPS&name=PIEVLI"

# Summary file
curl -o $path/sas.TAR "http://nxsa.esac.esa.int/nxsa-sl/servlet/data-action-aio?obsno=${obs}&level=ODF&extension=SAS"
tar -xf $path/sas.TAR -C $path
rm *ATS.FIT *TCS.FIT *RAS.ASC *ROS.ASC MANIFEST* *.TAR

# Fbktsr
curl -o $path/P${obs}PNS001FBKTSR0000.FTZ "http://nxsa.esac.esa.int/nxsa-sl/servlet/data-action-aio?obsno=${obs}&name=FBKTSR&instname=PN&level=PPS&extension=FTZ"

