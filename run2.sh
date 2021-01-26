#!/bin/bash

source /opt/intel/composer_xe_2015/bin/compilervars.sh intel64
module load netcdf4.6.1-intel15.0

input_start=2017-06-01
input_end=2017-09-30
startdate=$(date -I -d "$input_start") || exit -1
enddate=$(date -I -d "$input_end")     || exit -1
lon=(97.5 95 97.5 100 92.5 95 97.5 100 95 97.5 100 102.5 97.5 100 102.5 105 100 102.5 105 102.5 105)
lat=(-10.4211 -8.52632 -8.52632 -8.52632 -6.63158 -6.63158 -6.63158 -6.63158 -4.73684 -4.73684 -4.73684 -4.73684 -2.84211 -2.84211 -2.84211 -2.84211 -0.947368 -0.947368 -0.947368 0.947368 0.947368)
d="$startdate"
while [[ "$(date -d "$d" +%Y%m%d)" -le "$(date -d "$enddate" +%Y%m%d)" ]]; do

    date=$(date -d "$d" +%Y%m%d)
    casename=tc$date
    path=/data/W.eddie/SPCAM/${casename}/atm/hist/
    for ((i = 0; i < ${#lon[@]}; ++i)); do
        target_lon=${lon[$i]}
        target_lat=${lat[$i]}

        sed -i "4s|\(casename=\"\).*\"$|\1${casename}\"|" getCRM2.f90
        sed -i "5s|\(path=\"\).*\"$|\1${path}\"|" getCRM2.f90
        sed -i "6s|\(target_lon=\).*\(, target_lat=\).* \! target grid|\1${target_lon}\2${target_lat} \! target grid|" getCRM2.f90

#ifort getCRM2.f90 -o getCRM2.exe -I/opt/netcdf-intel15/include -L/opt/netcdf-intel15/lib -lnetcdff -lnetcdf -mcmodel medium -shared-intel
        ifort getCRM2.f90 -o getCRM2.exe -I/opt/netcdf-intel15/include -L/opt/netcdf-intel15/lib -lnetcdff -lnetcdf -mcmodel medium -shared-intel -O3 -xHost -fp-model precise

        ./getCRM2.exe
    done
    d=$(date -I -d "$d + 1 day")

done
