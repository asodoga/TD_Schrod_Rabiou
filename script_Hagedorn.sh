#!/bin/bash
#
#FPMDIR=build/gfortran*/TD_SCHROD
#rm -r build
make
#fpm build --flag  "-fopenmp"
#gfortran -fopenmp -o TD_SCHRfourierMDIR/app_TD_SCHROD.f90.o $FPMDIR/libTD_SCHROD.a $QMLDIR/libpot.a $QMLDIR/libAD_dnSVM.a -lblas -llapack
#./TD_SCHROD.x <DAT_files/dat_fourier> resultat
mkdir -p ../results
./TD_SCHROD.x <DAT_files/dat_Hagedorn1> resultat
 rename -v 's/.dat/nb1.dat/' *.dat
cp *.dat ../results
rm *.dat
rm fort.*
rm resultat
#make clean
./TD_SCHROD.x <DAT_files/dat_Hagedorn2> resultat
 rename -v 's/.dat/nb2.dat/' *.dat
cp *.dat ../results
rm *.dat
rm fort.*
rm resultat
#make clean
./TD_SCHROD.x <DAT_files/dat_Hagedorn3> resultat
 rename -v 's/.dat/nb3.dat/' *.dat
cp *.dat ../results
rm *.dat
rm fort.*
rm resultat
#make clean
./TD_SCHROD.x <DAT_files/dat_Hagedorn4> resultat
 rename -v 's/.dat/nb4.dat/' *.dat
cp *.dat ../results
rm *.dat
rm fort.*
rm resultat
#make clean

