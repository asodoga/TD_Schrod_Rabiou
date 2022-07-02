#!/bin/bash
QMLDIR=/home/elprof/QuantumModelLib
FPMDIR=build/gfortran*/TD_SCHROD
rm -r build
fpm build --flag  "-fopenmp"
gfortran -fopenmp -o TD_SCHROD.exe  $FPMDIR/app_TD_SCHROD.f90.o $FPMDIR/libTD_SCHROD.a $QMLDIR/libpot.a -lblas -llapack
./TD_SCHROD.exe <DAT_files/dat_dia> res


