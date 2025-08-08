#!/bin/bash

gfortran -o spectrum_std2D_t60.x spectre.f90

nb=$1
echo $nb

cd results_anhar-std2D-t60_nb$nb || { echo "Directory results_anhar-std2D-t60_nb$nb not found!"; exit 1; }

../spectrum_std2D_t60.x << ** > spect_std2D_t60.txt
&param Emin=0. Emax=8.0 file_auto='auto_cor_non_hagedorn_taylor.txt' option=1 dE=0.0001 /
**

mv file_spectrum.txt ../file_spectrum_anhar-std2D-t60-nb$nb.txt
cd ..