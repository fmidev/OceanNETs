#!/bin/bash
set -x
exp='D30H'
which basename
basename

for exp in 'D30H' 'D40H' 'C40H' 'E40H'
do
    DATAPATH='/Volumes/ONETs-SSD/'${exp}'/'

    for file in ${DATAPATH}/*nc
    do  
        mkdir -p ${DATAPATH}/1x1
        echo ${file}
        ff="$(basename ${file})"
        outdata=${ff%.nc}.1x1.nc
        echo ${outdata}
        cdo remapbil,mygrid $file ${DATAPATH}/1x1/${outdata} 
    done  
done