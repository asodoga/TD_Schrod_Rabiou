gfortran -o spectrum_H.x spectre.f90

nb=$1
echo $nb
cd results_H_nb$nb

../spectrum_H.x << ** > spect_H.txt
&param Emin=0. Emax=5. file_auto='auto_cor_hagedorn_taylor.txt' option=2 dE=0.01 /

**

mv file_spectrum.txt ../file_spectrum_Hnb$nb.txt
cd ..
#gnuplot plot_spectrum_H.gp
