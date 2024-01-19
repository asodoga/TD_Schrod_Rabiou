#!/bin/bash
#
#FPMDIR=build/gfortran*/TD_SCHROD
#rm -r build
make clean
make
#fpm build --flag  "-fopenmp"
#gfortran -fopenmp -o TD_SCHRfourierMDIR/app_TD_SCHROD.f90.o $FPMDIR/libTD_SCHROD.a $QMLDIR/libpot.a $QMLDIR/libAD_dnSVM.a -lblas -llapack
mkdir -p ../results
#./TD_SCHROD.x <DAT_files/dat_vp.lua> resultat.lua
#./TD_SCHROD.x <DAT_files/dat_retinal_non_Hagedorn.lua> resultat.lua
./TD_SCHROD.x </home/rabiou_issa/rabiou.issa/TD_Schrod_Rabiou/DAT_files/data_std.lua> resultat_std.lua
#./TD_SCHROD.x <DAT_files/dat_fourier> resultat
#paste psi_Ha_hagedorn_taylor.dat psi_NHa_non_hagedorn_taylor.dat | awk '{print $2-$6,$3-$7 }'> psi_diff

 #rename -v 's/.dat/nb1.dat/' *.dat
#cp *.dat ../results
#rm *.dat
#rm fort.*
#rm resultat
#make clean
