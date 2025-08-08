#!/bin/bash

gfortran -o spectrum_H2D_t60_ttt.x spectre.f90

nb=$1
echo $nb

cd results_anhar-h2D-t60-ttt_nb$nb || { echo "Directory results_anhar-h2D-t60-ttt_nb$nb not found!"; exit 1; }

../spectrum_H2D_t60_ttt.x << ** > spect_H2D_t60_ttt.txt
&param Emin=0. Emax=8.0 file_auto='auto_cor_hagedorn_taylor.txt' option=1 dE=0.0001 /
**

mv file_spectrum.txt ../file_spectrum_anhar-h2D-t60-ttt-nb$nb.txt
cd ..