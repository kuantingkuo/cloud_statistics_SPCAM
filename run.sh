#!/bin/bash

target_lon=120.
target_lat=-10.

source /opt/intel/composer_xe_2015/bin/compilervars.sh intel64
module load netcdf4.6.1-intel15.0

sed -i "6s|\(target_lon=\).*\(, target_lat=\)\S\+|\1${target_lon}\2${target_lat}|" getCRM.f90

ifort getCRM.f90 -o getCRM.exe -I/opt/netcdf-intel15/include -L/opt/netcdf-intel15/lib -lnetcdff -lnetcdf -mcmodel medium -shared-intel

./getCRM.exe
