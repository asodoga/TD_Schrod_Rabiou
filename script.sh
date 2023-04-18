#!/bin/bash
#
#FPMDIR=build/gfortran*/TD_SCHROD
#rm -r build
make
#fpm build --flag  "-fopenmp"
#gfortran -fopenmp -o TD_SCHRfourierMDIR/app_TD_SCHROD.f90.o $FPMDIR/libTD_SCHROD.a $QMLDIR/libpot.a $QMLDIR/libAD_dnSVM.a -lblas -llapack
#./TD_SCHROD.x <DAT_files/dat_fourier> resultat
./TD_SCHROD.x <DAT_files/dat_Hagedorn2> resultat
#mkdir -p ../results
#cp fort.* ../results
#rm fort.*
