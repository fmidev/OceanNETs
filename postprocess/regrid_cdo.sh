#!/bin/bash
set -x
exp='D30H'
which basename
basename

#for exp in 
#for exp in 'D30H' 'E40H'
for exp in 'E40H'
do
    DATAPATH='/Volumes/ONETs-SSD/'${exp}'/'

    for file in ${DATAPATH}/*nc
    do  
        mkdir -p ${DATAPATH}/1x1
        echo ${file}
        ff="$(basename ${file})"
        outdata=${ff%.nc}.1x1.nc
        echo ${outdata}
        if [ ! -f "${DATAPATH}/1x1/${outdata}" ];then
            if [[ ${file} = co2_AERmon.* ]]
            then
                cdo remapbil,mygrid -sellevel,1 $file ${DATAPATH}/1x1/${outdata}         
            else
                cdo remapbil,mygrid $file ${DATAPATH}/1x1/${outdata} 
            fi
        fi
    done  
done
